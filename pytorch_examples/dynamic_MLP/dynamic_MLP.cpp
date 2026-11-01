#include <torch/torch.h>
#include <iostream>
#include <vector>

// Define a module with a variable number of layers
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
        // Iterate through all layers except the last one with ReLU activation
        for (size_t i = 0; i < layers->size() - 1; ++i) {
            // Cast the abstract module to the concrete Linear module type
            x = torch::relu(layers->ptr<torch::nn::LinearImpl>(i)->forward(x));
        }
        // Apply the final layer without activation
        x = layers->ptr<torch::nn::LinearImpl>(layers->size() - 1)->forward(x);
        return x;
    }
};

int main() {
    // Define architecture at runtime
    int64_t input_dim = 10;
    std::vector<int64_t> hidden_dims = {64, 32, 16}; // 3 variable hidden layers
    int64_t output_dim = 2;

    // Instantiate the model
    auto model = std::make_shared<DynamicMLP>(input_dim, hidden_dims, output_dim);

    // Create a dummy input tensor (batch_size = 4, input_dim = 10)
    torch::Tensor input = torch::randn({4, input_dim});

    // Run forward pass
    torch::Tensor output = model->forward(input);

    // Print results and structure
    std::cout << "Input Shape: " << input.sizes() << std::endl;
    std::cout << "Output Shape: " << output.sizes() << std::endl;
    std::cout << "\nModel Architecture:\n" << *model << std::endl;

    return 0;
}

