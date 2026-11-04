#include <iomanip>
#include <iostream>
#include <numbers>
#include <torch/torch.h>

class FullyConnectedNet : public torch::nn::Module {
public:
    FullyConnectedNet(int input_size, int hidden_size, int output_size) {
        fc1 = register_module("fc1", torch::nn::Linear(input_size, hidden_size));
        fc2 = register_module("fc2", torch::nn::Linear(hidden_size, output_size));
    }

    torch::Tensor forward(torch::Tensor x) {
        x = torch::relu(fc1->forward(x));
        x = fc2->forward(x);
        return x;
    }

private:
    torch::nn::Linear fc1{nullptr};
    torch::nn::Linear fc2{nullptr};
};


int main() {
    FullyConnectedNet net(1, 128, 1);
    const double pi = std::numbers::pi;

    auto train_x = (torch::rand({512, 1}) * 2.0 - 1.0) * pi;
    auto train_y = torch::sin(train_x);

    torch::nn::MSELoss loss_function;
    torch::optim::Adam optimizer(
        net.parameters(),
        torch::optim::AdamOptions(0.01)
    );

    for (int epoch = 0; epoch < 1000; ++epoch) {
        auto predictions = net.forward(train_x);
        auto loss = loss_function(predictions, train_y);

        optimizer.zero_grad();
        loss.backward();
        optimizer.step();
    }

    auto final_loss = loss_function(net.forward(train_x), train_y);
    std::cout << std::fixed << std::setprecision(6)
              << "Final training loss: " << final_loss.item<float>() << '\n';

    auto test_x = torch::tensor({
        {-2.0 * pi},
        {-0.5 * pi},
        {0.0},
        {0.5 * pi},
        {2.0 * pi}
    });

    {
        torch::NoGradGuard no_grad;
        auto test_predictions = net.forward(test_x);

        std::cout << "x        predicted sin(x)   actual sin(x)\n";
        for (int i = 0; i < test_x.size(0); ++i) {
            auto x_value = test_x[i].item<float>();
            auto prediction = test_predictions[i].item<float>();
            auto actual = std::sin(x_value);
            std::cout << std::setw(7) << std::setprecision(3) << x_value
                      << "  " << std::setw(16) << std::setprecision(3) << prediction
                      << "  " << std::setw(14) << std::setprecision(3) << actual << '\n';
        }
    }

    return 0;
}
