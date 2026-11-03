#include <torch/torch.h>
#include <iostream>
#include <cmath>

using namespace torch::indexing;

int main() {
    // Check if CUDA is available for GPU acceleration
    torch::Device device(torch::kCPU);
    if (torch::cuda::is_available()) {
        device = torch::Device(torch::kCUDA);
        std::cout << "Using GPU: CUDA is available!" << std::endl;
    } else {
        std::cout << "Using CPU." << std::endl;
    }

    // --- Simulation Parameters ---
    int Nx = 500;                  // Number of spatial grid points
    float L = 1.0f;                // Domain length [0, L]
    float dx = L / (Nx - 1);       // Spatial step size
    
    float c = 1.0f;                // Wave propagation speed
    float T = 2.0f;                // Total simulation time
    float dt = 0.001f;             // Time step size
    int Nt = static_cast<int>(T / dt);

    // Verify CFL condition
    float cfl = c * (dt / dx);
    std::cout << "CFL Number: " << cfl << " (Must be <= 1.0 for stability)" << std::endl;
    if (cfl > 1.0f) {
        std::cerr << "Warning: CFL condition violated! Simulation may become unstable." << std::endl;
    }

    float alpha = std::pow(c * dt / dx, 2);

    // --- Grid and Initial Conditions ---
    auto x = torch::linspace(0, L, Nx, torch::kFloat32).to(device);
    
    // Initial displacement: u(x, 0) = sin(2 * pi * x / L)
    auto u_prev = torch::sin(2.0f * M_PI * x / L);
    
    // Initial velocity: u_t(x, 0) = 0 -> u^1 can be approximated from u^0
    auto u_curr = u_prev.clone();
    auto u_next = torch::zeros_like(x);

    // Apply first time step (assuming zero initial velocity)
    u_curr.slice(/*dim=*/0, /*start=*/1, /*end=*/Nx - 1) = 
        u_prev.slice(0, 1, Nx - 1) + 0.5f * alpha * (
            u_prev.slice(0, 2, Nx) - 
            2.0f * u_prev.slice(0, 1, Nx - 1) + 
            u_prev.slice(0, 0, Nx - 2)
        );

    // --- Time-Stepping Loop ---
    for (int n = 1; n < Nt; ++n) {
        // Vectorized interior point update using tensor slices
        auto u_center = u_curr.slice(0, 1, Nx - 1);
        auto u_right  = u_curr.slice(0, 2, Nx);
        auto u_left   = u_curr.slice(0, 0, Nx - 2);
        auto u_p      = u_prev.slice(0, 1, Nx - 1);

        // FDTD update formula
        u_next.slice(0, 1, Nx - 1) = 2.0f * u_center - u_p + alpha * (u_right - 2.0f * u_center + u_left);

        // Apply Dirichlet Boundary Conditions: u(0, t) = 0, u(L, t) = 0
        u_next[0] = 0.0f;
        u_next[Nx - 1] = 0.0f;

        // Shift time levels for the next iteration
        u_prev.copy_(u_curr);
        u_curr.copy_(u_next);

        // Optional: Print progress
        if (n % 200 == 0) {
            std::cout << "Step " << n << "/" << Nt 
                      << " | Max Amplitude: " << u_curr.max().item<float>() << std::endl;
        }
    }

    std::cout << "Wave simulation finished successfully!" << std::endl;
    return 0;
}
