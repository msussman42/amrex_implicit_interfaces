#include <torch/torch.h>
#include <iostream>
#include <functional>

class KalmanFilter {
public:

	//Future: upwind method u_t + a(x,t)u_x=0
	//or conservation form upwind method u_t + (a(x,t)u)_x=0
    KalmanFilter(
        int state_dim,
        int obs_dim,
        torch::Tensor F, //state transition matrix
        torch::Tensor H, // Measurement matrix
        torch::Tensor Q, // process noise covariance
        torch::Tensor R, // measurement noise covariance
        torch::Tensor P, // Initial estimate error covariance
        torch::Tensor x0 // Initial state estimate
        torch::Device device = torch::kCPU
    ) : nx_(state_dim),
        nz_(obs_dim),
        F_(std::move(F)),
        H_(std::move(H)),
        Q_(Q.to(device)),
        R_(R.to(device)),
        P_(std::move(P)),
        x_(std::move(x0)),
        device_(device) {
        
        // Precompute Cholesky factors for sampling Gaussian process & measurement noise
        L_Q_ = torch::linalg_cholesky(Q_);
        L_R_ = torch::linalg_cholesky(R_);
        I_ = torch::eye(x_.size(0), x_.options());
    }

    // Initialize ensemble around an initial state mean and covariance
    void initialize(
      torch::Tensor& x0, 
      const torch::Tensor& P0) {

      //P0=L_P0 * L_P0^*
     auto L_P0 = torch::linalg_cholesky(P0.to(device_));
     auto noise = torch::randn({nx_}, 
       torch::TensorOptions().device(device_));
     x_ = x0.to(device_) + torch::matmul(noise, L_P0.t());
    }

    torch::Tensor predict(double alpha,int n) {

        x_ = torch::matmul(F_, x_);
        P_ = torch::matmul(torch::matmul(F_, P_), F_.transpose(0, 1)) + Q_;
        return x_;

    }

    torch::Tensor update(const torch::Tensor& z) {
        torch::Tensor z_dev = z.to(device_); // [nz]

        // Innovation (residual): y = z - H * x
        torch::Tensor y = z - torch::matmul(H_, x_);

        // Innovation covariance: S = H * P * H^T + R
        torch::Tensor PHt = torch::matmul(P_, H_.transpose(0, 1));
        torch::Tensor S = torch::matmul(H_, PHt) + R_;

        // Optimal Kalman Gain: K = P * H^T * S^-1
        // Using linalg_solve (S * K^T = (P * H^T)^T) is 
	// more numerically stable than matrix inverse
        torch::Tensor K = torch::linalg_solve(S, PHt.transpose(0, 1)).transpose(0, 1);

        x_ = x_ + torch::matmul(K, y);

        torch::Tensor I_KH = I_ - torch::matmul(K, H_);
        P_ = torch::matmul(torch::matmul(I_KH, P_), I_KH.transpose(0, 1))
           + torch::matmul(torch::matmul(K, R_), K.transpose(0, 1));

        return x_;

    }

private:
    int nx_, nz_;
    torch::Tensor F_,H_,Q_, R_, L_Q_, L_R_,P_,x_,I_;
    torch::Device device_;
};

int main() {
    torch::manual_seed(42);
    torch::NoGradGuard no_grad; // Disable autograd overhead for standard filtering

    // Select GPU if available
    torch::Device device(torch::cuda::is_available() ? torch::kCUDA : torch::kCPU);
    std::cout << "Running EnKF on: " << (device.is_cuda() ? "CUDA" : "CPU") << "\n\n";

    // Problem dimensions: 2D constant-velocity tracking [x, y, vx, vy]
    const int N = 1000;   // Ensemble size
    const int Nx = 400;    // State dimension
    const int Nz = Nx;    // Observation dimension (position only)
    double L = 2.0;      // Domain length [0, L]
    double dx = L / (Nx - 1);       // Spatial step size (finite diff)
    double c = 1.0;                // Wave propagation speed
    double T = 1.0;                // Total simulation time
    const double dt = 0.001;
    int Nt = static_cast<int>(T / dt);
    // Verify CFL condition
    double cfl = c * (dt / dx);
    std::cout << "CFL Number: " << cfl << " (Must be <= 1.0 for stability)" << std::endl;
    if (cfl > 1.0) {
        std::cerr << "Warning: CFL condition violated! Simulation may become unstable." << std::endl;
    }

    double alpha = c * dt / dx;

    // --- Grid and Initial Conditions ---
    auto x = torch::linspace(-0.5*L, 0.5*L, Nx, torch::kFloat64).to(device);
    auto x_shift = 
	torch::linspace(-0.5*L-c*T,0.5*L-c*T,Nx,torch::kFloat64).to(device);

    // Nonlinear/Linear batch dynamics model: x_{k+1} = f(x_k, dt)
    auto dynamics_model = [](
      torch::Tensor& u, 
      double alpha,
      int n) {

     torch::Tensor u_next = u.clone();

     // Ensemble shape: [N, nx]
     u_next=(1.0-alpha)*u+alpha*torch::roll(u,-1,1); //shift up by 1
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

    auto u0_guess = torch::sin(2.0 * M_PI * x / L)+
	    torch::randn({Nx},device)*process_noise_sd;  // u0 noise

    torch::Tensor P0=torch::eye(Nx, device)*process_noise_variance; //P0 noise

    enkf.initialize(u0_guess, P0);

    // Run simulation loop
    for (int n = 0; n < Nt; ++n) {

     auto u_true_unsqueeze=u_true.unsqueeze(0);
     auto u_true_save=u_true.clone();

     u_true = dynamics_model(
	u_true_unsqueeze,
	alpha,n).squeeze(0);

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
     auto x_shift_inner = 
	torch::linspace(-0.5*L-c*cur_time,0.5*L-c*cur_time, 
			Nx,torch::kFloat64).to(device);
     auto u_exact=torch::sin(2.0 * M_PI * x_shift_inner / L);
     auto abs_error_truth=u_exact-truth;
     auto abs_error_norm_truth = torch::mean(torch::square(abs_error_truth));
     abs_error_norm_truth=torch::sqrt(abs_error_norm_truth);
     std::cout << "n " << n 
            << " | abs_error_norm_truth: " << abs_error_norm_truth << '\n';
     
    }

    return 0;
}

