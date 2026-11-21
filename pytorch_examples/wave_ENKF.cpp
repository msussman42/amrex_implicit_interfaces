#include <torch/torch.h>
#include <iostream>
#include <functional>

class EnsembleKalmanFilter {
public:
    using DynamicsFn = std::function<torch::Tensor(
	torch::Tensor&, /* u */
        torch::Tensor&, /* unm1 */
	double,  /* alpha */
	int)>;   /* n */
    using ObsFn = std::function<torch::Tensor(torch::Tensor&)>;

    EnsembleKalmanFilter(
        int ensemble_size,
        int state_dim,
        int obs_dim,
        DynamicsFn f,
        ObsFn h,
        torch::Tensor Q,
        torch::Tensor R,
        torch::Device device = torch::kCPU
    ) : N_(ensemble_size),
        nx_(state_dim),
        nz_(obs_dim),
        f_(std::move(f)),
        h_(std::move(h)),
        Q_(Q.to(device)),
        R_(R.to(device)),
        device_(device) {
        
        // Precompute Cholesky factors for sampling Gaussian process & measurement noise
        L_Q_ = torch::linalg_cholesky(Q_);
        L_R_ = torch::linalg_cholesky(R_);
    }

    // Initialize ensemble around an initial state mean and covariance
    void initialize(
      torch::Tensor& x0, 
      torch::Tensor& x0_nm1,
      const torch::Tensor& P0) {

      //P0=L_P0 * L_P0^*
     auto L_P0 = torch::linalg_cholesky(P0.to(device_));
     auto noise = torch::randn({N_, nx_}, 
       torch::TensorOptions().device(device_));
      // Ensemble shape: [N, nx]
     ensemble_ = x0.to(device_).unsqueeze(0) + torch::matmul(noise, L_P0.t());
     ensemble_nm1_ = ensemble_.clone();
    }

    // Forecast (Predict) Step
    void predict(double alpha,int n) {
        // Propagate each ensemble member through the nonlinear dynamics
        // f_ expects [N, nx] and returns [N, nx]
        ensemble_ = f_(ensemble_,ensemble_nm1_,alpha,n);

        // Add process noise: w ~ N(0, Q)
        auto w = torch::matmul(
            torch::randn({N_, nx_}, torch::TensorOptions().device(device_)),
            L_Q_.t()
        );
        ensemble_ = ensemble_ + w;
    }

    // Analysis (Update) Step using Perturbed Observations
    void update(const torch::Tensor& z) {
        torch::Tensor z_dev = z.to(device_); // [nz]

        // 1. Map ensemble into observation space: [N, nz]
        torch::Tensor Z_ens = h_(ensemble_);

        // 2. Compute ensemble means
        torch::Tensor x_mean = ensemble_.mean(0, /*keepdim=*/true); // [1, nx]
        torch::Tensor z_mean = Z_ens.mean(0, /*keepdim=*/true);     // [1, nz]

        // 3. Compute ensemble anomalies (deviations from the mean)
        torch::Tensor A_x = ensemble_ - x_mean; // [N, nx]
        torch::Tensor A_z = Z_ens - z_mean;     // [N, nz]

        // 4. Compute sample covariances: P_xz and P_zz
        double norm = 1.0 / static_cast<double>(N_ - 1);
        torch::Tensor P_xz = norm * torch::matmul(A_x.t(), A_z); // [nx, nz]
        torch::Tensor P_zz = norm * torch::matmul(A_z.t(), A_z) + R_; // [nz, nz]

        // 5. Compute Kalman Gain: K = P_xz * P_zz^{-1}
        // Using torch::linalg_solve on transposes avoids explicit matrix inversion:
        // P_zz * K^T = P_xz^T  =>  K = (P_zz^{-1} * P_xz^T)^T
        torch::Tensor K = torch::linalg_solve(P_zz, P_xz.t(), /*left=*/true).t(); // [nx, nz]

        // 6. Generate perturbed observations for each ensemble member
        auto v = torch::matmul(
            torch::randn({N_, nz_}, torch::TensorOptions().device(device_)),
            L_R_.t()
        );
        torch::Tensor Z_perturbed = z_dev.unsqueeze(0) + v; // [N, nz]

        // 7. Update ensemble members
        torch::Tensor innovation = Z_perturbed - Z_ens; // [N, nz]
        ensemble_ = ensemble_ + torch::matmul(innovation, K.t()); // [N, nx]
    }

    torch::Tensor get_mean() const {
        return ensemble_.mean(0);
    }

    torch::Tensor get_covariance() const {
        auto A_x = ensemble_ - ensemble_.mean(0, true);
        return torch::matmul(A_x.t(), A_x) / static_cast<double>(N_ - 1);
    }

private:
    int N_, nx_, nz_;
    DynamicsFn f_;
    ObsFn h_;
    torch::Tensor Q_, R_, L_Q_, L_R_;
    torch::Tensor ensemble_;
    torch::Tensor ensemble_nm1_;
    torch::Device device_;
};

int main() {
    torch::manual_seed(42);
    torch::NoGradGuard no_grad; // Disable autograd overhead for standard filtering

    // Select GPU if available
    torch::Device device(torch::cuda::is_available() ? torch::kCUDA : torch::kCPU);
    std::cout << "Running EnKF on: " << (device.is_cuda() ? "CUDA" : "CPU") << "\n\n";

    const int N = 1000;   // Ensemble size
    const int Nx = 200;    // State dimension
    const int Nz = Nx;    // Observation dimension (position only)
    double L = 2.0;      // Domain length [0, L]
    double dx = L / (Nx - 1);       // Spatial step size (finite diff)
    double c = 1.0;                // Wave propagation speed
    double T = 1.0;                // Total simulation time
    const double dt = 0.002;
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

    // Nonlinear/Linear batch dynamics model: x_{k+1} = f(x_k, dt)
    auto dynamics_model = [](
      torch::Tensor& u, 
      torch::Tensor& unm1,
      double alpha,
      int n) {

     torch::Tensor u_next = u.clone();

      // Ensemble shape: [N, nx]
     if (n==0) {
      u_next = u + 0.5 * alpha * (
	 torch::roll(u,1,1)-  //a(i,j) <- a(i,j-1)
	 2.0*u+ 
         torch::roll(u,-1,1)); //a(i,j) <- a(i,j+1)
     } else {
      u_next=2.0*u-unm1 + alpha * (
       torch::roll(u,1,1)-  
       2.0*u+ 
       torch::roll(u,-1,1)); 
     }
     unm1.copy_(u);
     u.copy_(u_next);
     return u_next;
    };

    // Observation model: z_k = h(x_k) -> observes [x, y]
    auto obs_model = [](torch::Tensor& x) {
     torch::Tensor x_obs = x.clone();
     return x_obs;
    };

    double noise_variance=0.01;
    double noise_sd=std::sqrt(noise_variance);

    double process_noise_variance=0.5;
    double process_noise_sd=std::sqrt(process_noise_variance);

    // Noise covariances
    torch::Tensor Q = torch::eye(Nx) * process_noise_variance; // Process noise
    torch::Tensor R = torch::eye(Nz) * noise_variance; // Measurement noise

    // Instantiate filter
    EnsembleKalmanFilter enkf(N, Nx, Nz, dynamics_model, obs_model, 
		Q, R, device);

    auto u_true = torch::sin(2.0 * M_PI * x / L);
    auto unm1_true=u_true.clone();

    auto u0_guess = torch::sin(2.0 * M_PI * x / L)+
	    torch::randn({Nx},device)*process_noise_sd;  // u0 noise
    auto unm1_guess=u0_guess.clone();

    torch::Tensor P0=torch::eye(Nx, device)*process_noise_variance; //P0 noise

    enkf.initialize(u0_guess, unm1_guess, P0);

    // Run simulation loop
    for (int n = 0; n < Nt; ++n) {

     auto u_true_unsqueeze=u_true.unsqueeze(0);
     auto u_true_save=u_true.clone();
     auto unm1_true_unsqueeze=unm1_true.unsqueeze(0);

     u_true = dynamics_model(
	u_true_unsqueeze,
	unm1_true_unsqueeze,
	alpha,n).squeeze(0);

     unm1_true.copy_(u_true_save);
        
     auto u_true_unsqueeze_update=u_true.unsqueeze(0);

     torch::Tensor z = obs_model(u_true_unsqueeze_update).squeeze(0) + 
                          torch::randn({Nz}, device) * noise_sd;

     enkf.predict(alpha,n);
     enkf.update(z);

     torch::Tensor est = enkf.get_mean().cpu();
     torch::Tensor truth = u_true.cpu();
     auto abs_error=est-truth;
     auto abs_error_norm = torch::mean(torch::square(abs_error));
     abs_error_norm=torch::sqrt(abs_error_norm);

     std::cout << "n " << n 
            << " | abs_error_norm: " << abs_error_norm << '\n';

     double cur_time=(n+1)*dt;
     auto u_exact=torch::sin(2.0 * M_PI * x / L)*cos(M_PI*cur_time*c);
     auto abs_error_truth=u_exact-truth;
     auto abs_error_norm_truth = torch::mean(torch::square(abs_error_truth));
     abs_error_norm_truth=torch::sqrt(abs_error_norm_truth);
     std::cout << "n " << n 
            << " | abs_error_norm_truth: " << abs_error_norm_truth << '\n';
     
    }

    return 0;
}

