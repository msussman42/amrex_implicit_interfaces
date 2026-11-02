//#include <torch/torch.act.h>
#include <torch/torch.h>
#include <iostream>
#include <vector>

// 1. Define a Dynamic MLP that builds layers based on a dynamic vector configuration
struct DynamicMLP : torch::nn::Module {
    torch::nn::Sequential network;

    DynamicMLP(const std::vector<int64_t>& layers) {
        for (size_t i = 0; i < layers.size() - 1; ++i) {
            // Add Linear Layer dynamically
            network->push_back(torch::nn::Linear(layers[i], layers[i+1]));
            
            // Add Tanh activation for hidden layers (important for higher-order derivatives)
            if (i < layers.size() - 2) {
                network->push_back(torch::nn::Functional(torch::tanh));
            }
        }
        register_module("network", network);
    }

    torch::Tensor forward(torch::Tensor x) {
        return network->forward(x);
    }
};

int main() {
 
    //gradcheck diagnostic
    // 1. Inputs must be double precision (kDouble) and require gradients
    auto x_first = torch::randn({2, 3}, torch::dtype(torch::kDouble).requires_grad(true));
    auto y_first = torch::randn({2, 3}, torch::dtype(torch::kDouble).requires_grad(true));

    // 2. Define the functional operation to test
    auto func = [](torch::TensorList inputs) -> std::vector<torch::Tensor> {
        auto out = inputs[0] * inputs[1]; // e.g., element-wise multiplication
        return {out};
    };

  // 1. Define an input tensor with requires_grad = true (use float64/double for precision)
  torch::Tensor x = torch::tensor({2.0, 3.0}, torch::TensorOptions().dtype(torch::kFloat64).requires_grad(true));

  // 2. Forward pass via Autograd
  // Example function: y = sum(2 * x^3 + 5 * x)
  auto compute_func = [](torch::Tensor input) {
    return (2.0 * input.pow(3) + 5.0 * input).sum();
  };

  torch::Tensor y = compute_func(x);

  // 3. Backward pass (Backpropagation)
  y.backward();
  torch::Tensor autograd_grad = x.grad();

  std::cout << "Autograd Gradient:\n" << autograd_grad << std::endl;

  // 4. Compute Finite Difference Derivative manually
  double h = 1e-5;
  torch::Tensor finite_diff_grad = torch::zeros_like(x);

  for (int i = 0; i < x.numel(); ++i) {
    // Clone tensor data to perturb index i
    torch::Tensor x_plus = x.clone().detach();
    torch::Tensor x_minus = x.clone().detach();

 //   x_plus[i].item<double>() += h;
 //   x_minus[i].item<double>() -= h;
    x_plus[i] += h;
    x_minus[i]-= h;

    // Numerical approximation: f(x + h) - f(x - h) / (2 * h) (central difference)
    torch::Tensor y_plus = compute_func(x_plus);
    torch::Tensor y_minus = compute_func(x_minus);

    double diff = (y_plus.item<double>() - y_minus.item<double>()) / (2.0 * h);
    finite_diff_grad[i] = diff;
  }

  std::cout << "Finite Difference Gradient:\n" << finite_diff_grad << std::endl;

  // 5. Compare the two gradients
  double max_diff = torch::max(torch::abs(autograd_grad - finite_diff_grad)).item<double>();
  std::cout << "Max Absolute Difference: " << max_diff << std::endl;


    // Check for CUDA availability
    torch::Device device(torch::kCPU);
    if (torch::cuda::is_available()) {
        std::cout << "CUDA available! Training on GPU." << std::endl;
        device = torch::Device(torch::kCUDA);
    }

    // Define network topology dynamically: [Input, Hidden1, Hidden2, Output]
    std::vector<int64_t> layer_config = {1, 16, 32, 16, 1};
    auto model = std::make_shared<DynamicMLP>(layer_config);
    model->to(device);

    // Physical parameter: omega_0
    float omega_0 = 2.0f;

    // Generate Collocation Points (Interior domain points where physics is enforced)
    // We want to track gradients with respect to inputs to compute physical derivatives.
    auto t_collocation = torch::linspace(0.0, 4.0, 500, torch::TensorOptions().device(device))
                              .view({-1, 1})
                              .set_requires_grad(true);

    // Initial boundary condition tensors (t = 0)
    auto t_init = torch::zeros({1, 1}, torch::TensorOptions().device(device)).set_requires_grad(true);
    auto u_init_true = torch::ones({1, 1}, torch::TensorOptions().device(device)); // u(0) = 1
    auto du_init_true = torch::zeros({1, 1}, torch::TensorOptions().device(device)); // u'(0) = 0

    // Optimizer
    torch::optim::Adam optimizer(model->parameters(), torch::optim::AdamOptions(1e-3));

    // Training Loop
    for (int epoch = 1; epoch <= 3000; ++epoch) {
        optimizer.zero_grad();

        // --- 1. Compute Boundary / Initial Condition Loss ---
        auto u_init_pred = model->forward(t_init);
        
        // Compute first derivative at boundary using autograd
        auto du_init_pred = torch::autograd::grad(
            {u_init_pred}, {t_init}, 
            {torch::ones_like(u_init_pred)}, 
            /*create_graph=*/true, /*retain_graph=*/true, /*allow_unused=*/true
        )[0];

        auto loss_ic_u = torch::mse_loss(u_init_pred, u_init_true);
        auto loss_ic_du = torch::mse_loss(du_init_pred, du_init_true);

        // --- 2. Compute Physical Residual Loss (ODE) ---
        auto u_pred = model->forward(t_collocation);

        // First derivative: du/dt
        auto du_dt = torch::autograd::grad(
            {u_pred}, {t_collocation}, 
            {torch::ones_like(u_pred)}, 
            /*create_graph=*/true, /*retain_graph=*/true, /*allow_unused=*/true
        )[0];

        // Second derivative: d^2u/dt^2
        auto d2u_dt2 = torch::autograd::grad(
            {du_dt}, {t_collocation}, 
            {torch::ones_like(du_dt)}, 
            /*create_graph=*/true, /*retain_graph=*/true, /*allow_unused=*/true
        )[0];

        // Governing Equation Residual: u'' + (omega_0^2)*u = 0
        auto ode_residual = d2u_dt2 + (omega_0 * omega_0) * u_pred;
        auto loss_physics = torch::mse_loss(ode_residual, torch::zeros_like(ode_residual));

        // --- 3. Backpropagation ---
        // Summing the total loss components together
        auto total_loss = loss_ic_u + loss_ic_du + loss_physics;
        total_loss.backward();
        optimizer.step();

        // Print Progress
        if (epoch % 500 == 0) {
            std::cout << "Epoch: " << epoch 
                      << " | Total Loss: " << total_loss.item<float>() 
                      << " [IC u: " << loss_ic_u.item<float>() 
                      << ", Physics: " << loss_physics.item<float>() << "]" << std::endl;
        }
    }

    std::cout << "\nTraining Complete. Model is physics-informed!" << std::endl;
    return 0;
}
