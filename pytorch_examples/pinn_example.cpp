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
    auto x = torch::randn({2, 3}, torch::dtype(torch::kDouble).requires_grad(true));
    auto y = torch::randn({2, 3}, torch::dtype(torch::kDouble).requires_grad(true));

    // 2. Define the functional operation to test
    auto func = [](torch::TensorList inputs) -> std::vector<torch::Tensor> {
        auto out = inputs[0] * inputs[1]; // e.g., element-wise multiplication
        return {out};
    };

    // 3. Configure options (optional, defaults are usually sufficient)
    auto options = torch::autograd::GradcheckOptions()
        .eps(1e-6)      // Finite difference step size
        .atol(1e-5)     // Absolute tolerance
        .rtol(1e-3);    // Relative tolerance

    // 4. Run the check
    bool success = torch::autograd::gradcheck(func, {x, y}, options);

    std::cout << "Gradcheck passed: " << (success ? "YES" : "NO") << std::endl;

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
