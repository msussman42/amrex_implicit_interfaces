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
        torch::Tensor x0, // Initial state estimate
        torch::Device device = torch::kCPU
    ) : nx_(state_dim),
        nz_(obs_dim),
        F_(F.to(device)),
        H_(H.to(device)),
        Q_(Q.to(device)),
        R_(R.to(device)),
        P_(P.to(device)),
        x_(x0.to(device)),
        device_(device) {
        
        I_ = torch::eye(x_.size(0), x_.options());
        L_Q_ = torch::linalg_cholesky(Q_);
    }

    // Initialize ensemble around an initial state mean and covariance
    void initialize(
      torch::Tensor& x0, 
      const torch::Tensor& P0) {

      //P0=L_P0 * L_P0^*
     auto L_P0 = torch::linalg_cholesky(P0.to(device_));
     auto opts=torch::TensorOptions().dtype(torch::kFloat64).device(device_);
     auto noise = torch::randn({nx_},opts); 
     x_ = x0.to(device_) + torch::matmul(noise, L_P0.t());
    }

    torch::Tensor predict(double alpha,int n) {

        x_ = torch::matmul(F_, x_);
        auto opts=torch::TensorOptions().dtype(torch::kFloat64).device(device_);
        auto noise = torch::randn({nx_},opts); 
        x_ = x_ + torch::matmul(noise, L_Q_.t());
        P_ = torch::matmul(torch::matmul(F_, P_), F_.transpose(0, 1)) + Q_;
        return x_;

    }

    torch::Tensor update(const torch::Tensor& z) {
        torch::Tensor z_dev = z.to(device_); // [nz]

        // Innovation (residual): y = z - H * x
        auto y = z - torch::matmul(H_, x_);

        // Innovation covariance: S = H * P * H^T + R
        auto PHt = torch::matmul(P_, H_.transpose(0, 1));
        auto S = torch::matmul(H_, PHt) + R_;

        // Optimal Kalman Gain: K = P * H^T * S^-1
        // Using linalg_solve (S * K^T = (P * H^T)^T) is 
	// more numerically stable than matrix inverse
        auto K = torch::linalg_solve(S, PHt.transpose(0, 1)).transpose(0, 1);

        x_ = x_ + torch::matmul(K, y);

        auto I_KH = I_ - torch::matmul(K, H_);
        P_ = torch::matmul(torch::matmul(I_KH, P_), I_KH.transpose(0, 1))
           + torch::matmul(torch::matmul(K, R_), K.transpose(0, 1));

        return x_;

    }

private:
    int nx_, nz_;
    torch::Tensor F_,H_,Q_, R_, P_,x_,I_,L_Q_;
    torch::Device device_;
};

int main() {
    torch::manual_seed(42);
    torch::NoGradGuard no_grad; // Disable autograd overhead for standard filtering

    // Select GPU if available
    torch::Device device(torch::cuda::is_available() ? torch::kCUDA : torch::kCPU);
    std::cout << "Running EnKF on: " << (device.is_cuda() ? "CUDA" : "CPU") << "\n\n";

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

    double noise_variance=0.01;
    double noise_sd=std::sqrt(noise_variance);

    double process_noise_variance=0.5;
    double process_noise_sd=std::sqrt(process_noise_variance);

    std::cout << "passed code at (1) \n";

    // Noise covariances
    auto opts=torch::TensorOptions().dtype(torch::kFloat64).device(device);
    torch::Tensor Q = torch::eye(Nx,opts) * process_noise_variance; // Process noise
    torch::Tensor R = torch::eye(Nz,opts) * noise_variance; // Measurement noise
    torch::Tensor H = torch::eye(Nz,opts);
    torch::Tensor F = torch::eye(Nx,opts);
    for (int i=0;i<Nx;i++) {
     F.index_put_({i,i},1.0-alpha);
     if (i==0) {
      F.index_put_({i,Nx-2},alpha);
     } else {
      F.index_put_({i,i-1},alpha);
     }
    }
    torch::Tensor P0=torch::eye(Nx, opts)*process_noise_variance; //P0 noise
    auto u0 = torch::sin(2.0 * M_PI * x / L);
    std::cout << "passed code at (2) \n";

    // Instantiate filter
    KalmanFilter kf(Nx, Nz, F,H, 
		Q, R, P0, u0, device);
    kf.initialize(u0, P0);

    std::cout << "passed code at (3) \n";

    auto u_true = torch::sin(2.0 * M_PI * x / L);

    // Run simulation loop
    for (int n = 0; n < Nt; ++n) {

     u_true = torch::matmul(F, u_true);

     auto opts=torch::TensorOptions().dtype(torch::kFloat64).device(device);
     auto z = u_true + torch::randn({Nz}, opts) * noise_sd;

     auto u_pred=kf.predict(alpha,n);
     auto est=kf.update(z);

     auto truth = u_true.cpu();
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

