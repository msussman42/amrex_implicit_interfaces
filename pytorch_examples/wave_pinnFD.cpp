#include <torch/torch.h>
#include <iostream>
#include <cmath>

using namespace torch::indexing;

int main() {

    // Set OpenMP and MKL environment variables programmatically (or set them in your shell before running)
    setenv("OMP_NUM_THREADS", "1", 1);
    setenv("MKL_NUM_THREADS", "1", 1);
    setenv("OPENBLAS_NUM_THREADS", "1", 1);
    setenv("VECLIB_MAXIMUM_THREADS", "1", 1);
    setenv("NUMEXPR_NUM_THREADS", "1", 1);

    // Set PyTorch / LibTorch thread counts
    at::set_num_threads(1);
    at::set_num_interop_threads(1); // Optional: controls inter-op parallelism


    // Check if CUDA is available for GPU acceleration
    torch::Device device(torch::kCPU);
    if (torch::cuda::is_available()) {
        device = torch::Device(torch::kCUDA);
        std::cout << "Using GPU: CUDA is available!" << std::endl;
    } else {
        std::cout << "Using CPU." << std::endl;
    }

    // --- Simulation Parameters ---
    int Nx = 1000;                  // Number of spatial grid points
    double L = 2.0;                // Domain length [0, L]
    double dx = L / (Nx - 1);       // Spatial step size (finite diff)
    
    double c = 1.0;                // Wave propagation speed
    double T = 1.0;                // Total simulation time
    double dt = 0.0005;             // Time step size
    int Nt = static_cast<int>(T / dt);

    // Verify CFL condition
    double cfl = c * (dt / dx);
    std::cout << "CFL Number: " << cfl << " (Must be <= 1.0 for stability)" << std::endl;
    if (cfl > 1.0) {
        std::cerr << "Warning: CFL condition violated! Simulation may become unstable." << std::endl;
    }

    double alpha = std::pow(c * dt / dx, 2);

    // --- Grid and Initial Conditions ---
    auto x = torch::linspace(-0.5*L, 0.5*L, Nx, torch::kFloat64).to(device);
    
    // Initial displacement: u(x, 0) = sin(2 * pi * x / L)
    auto u_prev = torch::sin(2.0 * M_PI * x / L);
    
    // Initial velocity: u_t(x, 0) = 0 -> u^1 can be approximated from u^0
    auto u_curr = u_prev.clone();
    auto u_next = torch::zeros_like(x);

    // Apply first time step (assuming zero initial velocity)
    //u_curr.slice(/*dim=*/0, /*start=*/1, /*end=*/Nx - 1) = 
    //    u_prev.slice(0, 1, Nx - 1) + 0.5f * alpha * (
    //        u_prev.slice(0, 2, Nx) - 
    //        2.0f * u_prev.slice(0, 1, Nx - 1) + 
    //        u_prev.slice(0, 0, Nx - 2)
    //    );
    u_curr = u_prev + 0.5 * alpha * (
	 torch::roll(u_prev,1,0)-  //shift down by 1 row
	 2.0*u_prev+ 
         torch::roll(u_prev,-1,0)); //shift up by 1
    double cur_time=dt;

    // --- Time-Stepping Loop ---
    for (int n = 1; n < Nt; ++n) {
        u_next = 2.0*u_curr-u_prev + alpha * (
	    torch::roll(u_curr,1,0)-  //shift down by 1 row
	    2.0*u_curr+ 
	    torch::roll(u_curr,-1,0)); //shift up by 1
        cur_time+=dt;
        // Vectorized interior point update using tensor slices
        //auto u_center = u_curr.slice(0, 1, Nx - 1);
        //auto u_right  = u_curr.slice(0, 2, Nx);
        //auto u_left   = u_curr.slice(0, 0, Nx - 2);
        //auto u_p      = u_prev.slice(0, 1, Nx - 1);

        // FDTD update formula
        //u_next.slice(0, 1, Nx - 1) = 2.0f * u_center - u_p + alpha * (u_right - 2.0f * u_center + u_left);

        // Apply Dirichlet Boundary Conditions: u(0, t) = 0, u(L, t) = 0
        //u_next[0] = 0.0f;
        //u_next[Nx - 1] = 0.0f;

        // Shift time levels for the next iteration
        u_prev.copy_(u_curr);
        u_curr.copy_(u_next);

        // Optional: Print progress
        if (n % 200 == 0) {
            std::cout << "Step " << n << "/" << Nt 
                      << " | Max Amplitude: " << u_curr.max().item<double>() << std::endl;
        }
	auto u_exact=torch::sin(2.0 * M_PI * x / L)*cos(M_PI*cur_time*c);
	auto abs_error=u_exact-u_curr;
	auto abs_error_norm = torch::mean(torch::square(abs_error));
        abs_error_norm=torch::sqrt(abs_error_norm);

        std::cout << "\n time,error " << cur_time << ' ' << abs_error_norm <<
		'\n';
    }

    std::cout << "Wave simulation finished successfully!" << std::endl;
    return 0;
}
