//#include <torch/torch.act.h>
#include <torch/torch.h>
#include <iostream>
#include <vector>

// 1. Define a Dynamic MLP that builds layers based on a dynamic vector configuration
struct WavePINN : torch::nn::Module {
    torch::nn::Sequential network;

    WavePINN(const std::vector<int64_t>& layers) {
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


// Helper function to compute derivatives
std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor> 
  compute_derivatives(
    WavePINN& model, torch::Tensor x, torch::Tensor t) {
    
    // Concatenate inputs and enable gradient tracking
    auto xt = torch::cat({x, t}, /*dim=*/1);
    xt.requires_grad_(true);

    auto u = model->forward(xt);

    // First derivatives: du/d(xt)
    auto grads = torch::autograd::grad(
        /*outputs=*/{u},
        /*inputs=*/{xt},
        /*grad_outputs=*/{torch::ones_like(u)},
        /*create_graph=*/true,
        /*retain_graph=*/true
    );

    auto du_dxt = grads[0];
    auto du_dx = du_dxt.index({torch::indexing::Slice(), 0}).unsqueeze(1);
    auto du_dt = du_dxt.index({torch::indexing::Slice(), 1}).unsqueeze(1);

    // Second derivatives: d2u/dx2 and d2u/dt2
    auto du_dx_sum = du_dx.sum();
    auto du_dt_sum = du_dt.sum();

    auto d2u_dx2_grads = torch::autograd::grad({du_dx_sum}, {xt}, /*grad_outputs=*/{}, /*create_graph=*/true, /*retain_graph=*/true);
    auto d2u_dt2_grads = torch::autograd::grad({du_dt_sum}, {xt}, /*grad_outputs=*/{}, /*create_graph=*/true, /*retain_graph=*/true);

    auto d2u_dx2 = d2u_dx2_grads[0].index({torch::indexing::Slice(), 0}).unsqueeze(1);
    auto d2u_dt2 = d2u_dt2_grads[0].index({torch::indexing::Slice(), 1}).unsqueeze(1);

    return {u, du_dt, d2u_dx2, d2u_dt2};
}





int main() {

  torch::manual_seed(42);
 
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
    std::vector<int64_t> layer_config = {2, 16, 32, 16, 1};
    auto model = std::make_shared<WavePINN>(layer_config);
    model->to(device);

    // Optimizer
    torch::optim::Adam optimizer(model->parameters(), 
      torch::optim::AdamOptions(1e-3));

    float c=1.0f;
    int epochs=3000;

    // Training Loop
    for (int epoch = 1; epoch <= epochs; ++epoch) {
        optimizer.zero_grad();

// 1. Collocation points inside domain (x, t)
        auto x_pde = torch::rand({1000, 1}) * 2.0 - 1.0; // e.g., [-1, 1]
        auto t_pde = torch::rand({1000, 1});              // e.g., [0, 1]
        x_pde.requires_grad_(true);
        t_pde.requires_grad_(true);

        auto [u, du_dt, d2u_dx2, d2u_dt2] = 
          compute_derivatives(model, x_pde, t_pde);

// PDE residual: u_tt - c^2 * u_xx = 0
        auto pde_residual = d2u_dt2 - (c * c) * d2u_dx2;
        auto loss_pde = torch::mean(torch::square(pde_residual));

        // 2. Initial Conditions (t = 0)
        auto x_ic = torch::rand({200, 1}) * 2.0 - 1.0;
        auto t_ic = torch::zeros({200, 1});
        auto xt_ic = torch::cat({x_ic, t_ic},/*dim=*/ 1);
        xt_ic.requires_grad_(true);

        auto u_ic_pred = model->forward(xt_ic);

// Example IC: u(x, 0) = sin(pi * x), u_t(x, 0) = cos(pi*x)
        auto u_ic_true = torch::sin(M_PI * x_ic);
        auto loss_ic = torch::mean(torch::square(u_ic_pred - u_ic_true));

        auto u_ic2_true = torch::cos(M_PI * x_ic);

        // First derivatives: du/d(xt)
        auto grads = torch::autograd::grad(
         /*outputs=*/{u_ic_pred},
         /*inputs=*/{xt_ic},
         /*grad_outputs=*/{torch::ones_like(u_ic_pred)},
         /*create_graph=*/true,
         /*retain_graph=*/true
        );
        auto du_dxt_ic = grads[0];
        auto du_dx_ic = du_dxt_ic.index({torch::indexing::Slice(), 0}).
          unsqueeze(1);
        auto du_dt_ic = du_dxt_ic.index({torch::indexing::Slice(), 1}).
          unsqueeze(1);

        auto loss_ic2 = torch::mean(torch::square(du_dt_ic - u_ic2_true));

        // Total Loss
        auto loss = loss_pde + loss_ic + loss_ic2;

        loss.backward();
        optimizer.step();

        // Print Progress
        if (epoch % 200 == 0) {
            std::cout << "Epoch: " << epoch 
                      << " | Loss: " << loss.item<float>() 
                      << " | PDE: " << loss_pde.item<float>() 
                      << " [IC u: " << loss_ic.item<float>() 
                      << " [IC2 u: " << loss_ic2.item<float>() 
                      << "]" << std::endl;
        }
    }

    std::cout << "\nTraining Complete. Model is physics-informed!" << std::endl;
    return 0;
}
