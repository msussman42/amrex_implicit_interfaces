#include <torch/torch.h>
#include <iostream>
#include <vector>

// 1. Define a simple MLP Net for the PINN
struct PINNNet : torch::nn::Module {
    torch::nn::Linear fc1{nullptr}, fc2{nullptr}, fc3{nullptr};

    PINNNet(int input_dim, int hidden_dim, int output_dim) {
        fc1 = register_module("fc1", torch::nn::Linear(input_dim, hidden_dim));
        fc2 = register_module("fc2", torch::nn::Linear(hidden_dim, hidden_dim));
        fc3 = register_module("fc3", torch::nn::Linear(hidden_dim, output_dim));
    }

    torch::Tensor forward(torch::Tensor x) {
        x = torch::tanh(fc1->forward(x));
        x = torch::tanh(fc2->forward(x));
        x = fc3->forward(x);
        return x;
    }
};

int main() {
    // Set device
    torch::Device device(torch::kCPU);
    std::cout << "Training on CPU." << std::endl;

    // Initialize model
    PINNNet model(1, 32, 1);
    model.to(device);

    // Optimizer
    torch::optim::Adam optimizer(model.parameters(), torch::optim::AdamOptions(1e-3));

    // Generate collocation points in 1D (e.g., t or x in [0, 1])
    torch::Tensor x_collocation = torch::linspace(0.0, 1.0, 100, torch::TensorOptions().dtype(torch::kFloat32).requires_grad(true)).unsqueeze(1);
    x_collocation = x_collocation.to(device);

    for (int epoch = 0; epoch < 1000; ++epoch) {
        optimizer.zero_grad();

        // Forward pass
        torch::Tensor u = model.forward(x_collocation);

        // Compute first derivative: du/dx
        // Note: create_graph=True allows differentiating through the gradient for physics loss
        std::vector<torch::Tensor> grads = torch::autograd::grad(
            {u}, {x_collocation}, 
            /*grad_outputs=*/{torch::ones_like(u)}, 
            /*create_graph=*/true, 
            /*retain_graph=*/true
        );
        torch::Tensor du_dx = grads[0];

        // Example Physics Residual: e.g., du/dx + u = 0 (or a custom PDE/ODE residual)
        torch::Tensor physics_residual = du_dx + u; 
        torch::Tensor physics_loss = torch::mean(torch::square(physics_residual));

        // Boundary condition loss (e.g., at x = 0, u(0) = 1.0)
        torch::Tensor x_bc = torch::zeros({1, 1}, torch::TensorOptions().dtype(torch::kFloat32)).to(device);
        torch::Tensor u_bc_pred = model.forward(x_bc);
        torch::Tensor bc_loss = torch::mean(torch::square(u_bc_pred - 1.0));

        // Total loss
        torch::Tensor total_loss = physics_loss + bc_loss;

        // Backward and optimize
        total_loss.backward();
        optimizer.step();

        if (epoch % 100 == 0) {
            std::cout << "Epoch [" << epoch << "/1000] Loss: " 
                      << total_loss.item<float>() << std::endl;
        }
    }

    return 0;
}

