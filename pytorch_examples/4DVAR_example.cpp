#include <torch/torch.h>
#include <iostream>
#include <vector>

// Define a simple dynamical model: x_{t+1} = alpha * x_t * (1 - x_t) [Logistic map]
torch::Tensor dynamical_model(const torch::Tensor& x, double alpha = 3.7) {
    return alpha * x * (1.0 - x);
}

// Define an observation operator: maps state space to observation space
torch::Tensor observation_operator(const torch::Tensor& x) {
    // For simplicity, let's assume we observe the state directly (Identity matrix H)
    return x;
}

int main() {
    // 1. Setup the Assimilation Window & Dimensions
    const int time_steps = 5;
    const int state_dim = 2;

    // True initial state (hidden from the assimilation system)
    torch::Tensor x_true = torch::tensor({0.4, 0.6});
    
    // Generate synthetic true trajectory and observations
    std::vector<torch::Tensor> observations;
    torch::Tensor x_curr = x_true.clone();
    
    // R: Observation error covariance matrix (Assuming variance = 0.01)
    torch::Tensor R_inv = torch::eye(state_dim) * (1.0 / 0.01); 

    std::cout << "--- Generating Synthetic Observations ---" << std::endl;
    for (int t = 0; t < time_steps; ++t) {
        x_curr = dynamical_model(x_curr);
        // Add random Gaussian noise to create observation y_t
        torch::Tensor noise = torch::randn({state_dim}) * 0.1;
        torch::Tensor y_t = observation_operator(x_curr) + noise;
        observations.push_back(y_t);
        std::cout << "Step " << t << " Obs: " << y_t << std::endl;
    }

    // 2. Background Information (Prior Guess)
    // x_b is a corrupted version of the true initial state
    torch::Tensor x_b = x_true + torch::tensor({0.15, -0.15});
    // B: Background error covariance matrix (Assuming variance = 0.04)
    torch::Tensor B_inv = torch::eye(state_dim) * (1.0 / 0.04);

    std::cout << "\nTrue x0:       " << x_true << std::endl;
    std::cout << "Background x0: " << x_b << std::endl;

    // 3. Define the Control Variable (Initial State x0 to optimize)
    // We start our optimization from the background guess
    torch::Tensor x0 = x_b.clone().set_requires_grad(true);

    // 4. Optimization Loop (Gradient Descent)
    // Using Adam optimizer to minimize the 4D-Var cost function
    torch::optim::Adam optimizer({x0}, torch::optim::AdamOptions(0.01));

    std::cout << "\n--- Starting 4D-Var Minimization ---" << std::endl;
    for (int epoch = 1; epoch <= 200; ++epoch) {
        optimizer.zero_grad();

        // Part 1 of Cost Function: Background penalty term
        torch::Tensor d_b = x0 - x_b;
        torch::Tensor J_b = 0.5 * torch::dot(d_b, torch::matmul(B_inv, d_b));

        // Part 2 of Cost Function: Observation penalty term aggregated over time
        torch::Tensor J_o = torch::zeros({});
        torch::Tensor x_traj = x0;

        for (int t = 0; t < time_steps; ++t) {
            // Propagate state forward in time
            x_traj = dynamical_model(x_traj);
            
            // Calculate departure from observation
            torch::Tensor d_o = observation_operator(x_traj) - observations[t];
            J_o += 0.5 * torch::dot(d_o, torch::matmul(R_inv, d_o));
        }

        // Total 4D-Var Cost Function: J = J_b + J_o
        torch::Tensor J = J_b + J_o;

        // Backward pass computes the exact gradient dJ/dx0 via automatic differentiation
        J.backward();
        optimizer.step();

        // Enforce physical constraints if necessary (e.g., bounds for logistic map)
        {
            torch::NoGradGuard no_grad;
            x0.clamp_(0.0, 1.0);
        }

        if (epoch % 20 == 0 || epoch == 1) {
            std::cout << "Epoch " << epoch << " | Total Cost J: " << J.item<double>() 
                      << " | x0 Guess: [" << x0[0].item<double>() << ", " << x0[1].item<double>() << "]" << std::endl;
        }
    }

    std::cout << "\n--- Assimilation Results ---" << std::endl;
    std::cout << "True Initial State:      " << x_true << std::endl;
    std::cout << "Background Initial State: " << x_b << std::endl;
    std::cout << "Estimated (Analysis) x0:  " << x0 << std::endl;

    return 0;
}

