#include <torch/torch.h>
#include <iostream>
#include <vector>

struct DynamicMLP : torch::nn::Module {
    torch::nn::ModuleList layers;

    DynamicMLP(int64_t input_size, const std::vector<int64_t>& hidden_sizes, int64_t output_size) {
        int64_t in_features = input_size;

        // Dynamically add hidden layers based on the vector size
        for (const auto& next_features : hidden_sizes) {
            // Create the linear layer and register it via ModuleList
            layers->push_back(torch::nn::Linear(in_features, next_features));
            in_features = next_features;
        }

        // Add the final output layer
        layers->push_back(torch::nn::Linear(in_features, output_size));

        // Register the ModuleList so LibTorch tracks its parameters
        register_module("layers", layers);

    }

    torch::Tensor forward(torch::Tensor x) {
        // Iterate through all layers except the last one with ReLU or tanh activation
        for (size_t i = 0; i < layers->size() - 1; ++i) {
            std::cout << "forward with i= " << i << '\n';
            // Cast the abstract module to the concrete Linear module type
            //x = torch::relu(layers->ptr<torch::nn::LinearImpl>(i)->forward(x));
            x = torch::tanh(layers->ptr<torch::nn::LinearImpl>(i)->forward(x));
            std::cout << "after forward with i= " << i << '\n';
        }
        // Apply the final layer without activation
        x = layers->ptr<torch::nn::LinearImpl>(layers->size() - 1)->forward(x);
        std::cout << "after forward final layer " << '\n';
        return x;
    }

};

int main() {
    // Set device
    torch::Device device(torch::kCPU);
    std::cout << "Training on CPU." << std::endl;

    // Define architecture at runtime
    //int64_t input_dim = 10;
    int64_t input_dim = 1;
    std::vector<int64_t> hidden_dims = {64, 32, 16}; // 3 variable hidden layers
    //int64_t output_dim = 2;
    int64_t output_dim = 1;

    std::cout << "building the model \n";
    // Instantiate the model
    auto model = std::make_shared<DynamicMLP>(input_dim, hidden_dims, output_dim);
    std::cout << "model built \n";

    model->to(device);
    std::cout << "after model to device \n";

    // Optimizer
    torch::optim::Adam optimizer(model->parameters(), torch::optim::AdamOptions(1e-3));
    std::cout << "after optimizer creation \n";

    // Generate collocation points in 1D (e.g., t or x in [0, 1])
    torch::Tensor x_collocation = torch::linspace(0.0, 1.0, 100, torch::TensorOptions().dtype(torch::kFloat32).requires_grad(true)).unsqueeze(1);
    std::cout << "after linspace \n";
    x_collocation = x_collocation.to(device);
    std::cout << "after x_collocation to device  \n";

    for (int epoch = 0; epoch < 1000; ++epoch) {
        std::cout << "epoch= " << epoch << '\n';

        optimizer.zero_grad();
        std::cout << "after optimizer zero grad " << '\n';

        // Forward pass
        torch::Tensor u = model->forward(x_collocation);
        std::cout << "after forward pass " << '\n';

        // Compute first derivative: du/dx
        // Note: create_graph=True allows differentiating through the gradient for physics loss
        std::vector<torch::Tensor> grads = torch::autograd::grad(
            {u}, {x_collocation}, 
            /*grad_outputs=*/{torch::ones_like(u)}, 
            /*create_graph=*/true, 
            /*retain_graph=*/true
        );
        std::cout << "after grad " << '\n';
        torch::Tensor du_dx = grads[0];
        std::cout << "after du_dx " << '\n';

        // Example Physics Residual: e.g., du/dx + u = 0 (or a custom PDE/ODE residual)
        torch::Tensor physics_residual = du_dx + u; 
        std::cout << "after du_dx+u " << '\n';
        torch::Tensor physics_loss = torch::mean(torch::square(physics_residual));
        std::cout << "after physics_loss " << '\n';

        // Boundary condition loss (e.g., at x = 0, u(0) = 1.0)
        torch::Tensor x_bc = torch::zeros({1, 1}, torch::TensorOptions().dtype(torch::kFloat32)).to(device);
        std::cout << "after x_bc loss " << '\n';
        torch::Tensor u_bc_pred = model->forward(x_bc);
        std::cout << "after u_bc_pred " << '\n';
        torch::Tensor bc_loss = torch::mean(torch::square(u_bc_pred - 1.0));
        std::cout << "after bc_loss " << '\n';

        // Total loss
        torch::Tensor total_loss = physics_loss + bc_loss;
        std::cout << "after total_loss " << '\n';

        // Backward and optimize
        total_loss.backward();
        std::cout << "after total_loss backward" << '\n';
        optimizer.step();
        std::cout << "after optimizer.step()" << '\n';

        if (epoch % 100 == 0) {
            std::cout << "Epoch [" << epoch << "/1000] Loss: " 
                      << total_loss.item<float>() << std::endl;
        }
    }

    return 0;
}

