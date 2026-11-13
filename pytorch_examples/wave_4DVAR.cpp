#include <torch/torch.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace torch::indexing;

torch::Tensor dynamical_model_regular(const torch::Tensor& u,
	       const torch::Tensor& unm1,double alpha) {

   return 2.0*u-unm1 + alpha * (
     torch::roll(u,1,0)-  //shift down by 1 row
	    2.0*u+ 
	    torch::roll(u,-1,0)); //shift up by 1

}

// Define an observation operator: maps state space to observation space
torch::Tensor observation_operator(const torch::Tensor& u) {
    // For simplicity, let's assume we observe the state directly (Identity matrix H)
    return u;
}

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
    int Nx = 200;                  // Number of spatial grid points
    double L = 2.0;                // Domain length [0, L]
    double dx = L / (Nx - 1);       // Spatial step size (finite diff)
    
    double c = 1.0;                // Wave propagation speed
    double T = 1.0;                // Total simulation time
    double dt = 0.0025;             // Time step size
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
    
    // True initial state (hidden from the assimilation system)
    // Initial displacement: u(x, 0) = sin(2 * pi * x / L)
    auto u_true = torch::sin(2.0 * M_PI * x / L);
    auto observation_final = torch::sin(2.0 * M_PI * x / L)*
	    cos(M_PI*T*c);

    std::vector<torch::Tensor> observations;
    std::vector<torch::Tensor> R_inv;
    for (int n = 0; n < Nt; ++n) {
     auto u_zero=torch::zeros_like(u_true);

     auto sizes=observation_final.sizes();
     int64_t rows=sizes[0];

     if (n<Nt-1) {
      observations.push_back(u_zero);
      R_inv.push_back(torch::zeros({rows,rows},torch::dtype(torch::kFloat64))); 
     } else if (n==Nt-1) {
      observations.push_back(observation_final);
      std::cout << "sizes= " << sizes << '\n';
      R_inv.push_back(torch::eye(rows,rows,torch::dtype(torch::kFloat64))); 
     }
    }

    auto u_b=torch::zeros_like(u_true);

    auto u_true_sizes=u_true.sizes();
    std::cout << "u_true_sizes= " << u_true_sizes << '\n';
    int64_t u_true_rows=u_true_sizes[0];
    auto B_inv=torch::eye(u_true_rows,u_true_rows,
      torch::dtype(torch::kFloat64)); 

    torch::Tensor u0 = u_b.clone().set_requires_grad(true);
    torch::optim::Adam optimizer({u0}, torch::optim::AdamOptions(0.01));

    std::cout << "\n--- Starting 4D-Var Minimization ---" << std::endl;
    for (int epoch = 1; epoch <= 400; ++epoch) {
      optimizer.zero_grad();

      torch::Tensor d_b = u0 - u_b;
      torch::Tensor J_b = 
	   0.5 * 0.0 * torch::dot(d_b, torch::matmul(B_inv, d_b));
      torch::Tensor J_o = torch::zeros({});
      torch::Tensor u_traj = u0;
      auto u_curr = u_traj.clone();
      auto u_prev = u_traj.clone();

      for (int n = 0; n < Nt; ++n) {
    
       if (n==0) {
        u_traj = u_traj + 0.5 * alpha * (
	 torch::roll(u_traj,1,0)-  //shift down by 1 row
	 2.0*u_traj+ 
         torch::roll(u_traj,-1,0)); //shift up by 1
        u_curr.copy_(u_traj);
        u_prev.copy_(u0);
       } else {
        u_traj=dynamical_model_regular(u_curr,u_prev,alpha);
        u_prev.copy_(u_curr);
        u_curr.copy_(u_traj);
       }
       if ((n % 200 == 0)&&(1==0)) {
        std::cout << "epoch,Step " << epoch << ' ' << n << "/" << Nt 
         << " | Max Amplitude: " << u_curr.max().item<double>() << std::endl;
       }
       torch::Tensor d_o = observation_operator(u_traj) - observations[n];
       J_o += 0.5 * torch::dot(d_o, torch::matmul(R_inv[n], d_o));
      } // for (int n = 0; n < Nt; ++n) 

      // Total 4D-Var Cost Function: J = J_b + J_o
      torch::Tensor J = J_b + J_o;
      // Backward pass computes the exact gradient dJ/dx0 via 
      // automatic differentiation
      J.backward();
      optimizer.step();
      if (epoch % 20 == 0 || epoch == 1) {
       std::cout << "Epoch " << epoch << 
	       " | Total Cost J: " << J.item<double>() << '\n';
      }

      auto abs_error_initial=u0-u_true;
      auto abs_error_norm = torch::mean(torch::square(abs_error_initial));
      abs_error_norm=torch::sqrt(abs_error_norm);
      std::cout << "\n error_initial " << abs_error_norm << '\n';
      std::cout << "u0 when converged is the analysis (initial) state\n";
    } // for (int epoch = 1; epoch <= 200; ++epoch) 

    return 0;
}
