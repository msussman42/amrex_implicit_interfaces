#include <torch/torch.h>
#include <iostream>
#include <functional>


class KalmanFilter {
public:

	//Future: upwind method u_t + a(x,t)u_x=0
	//or conservation form upwind method u_t + (a(x,t)u)_x=0
	//Also future detection of collapse for nonlinear problems:
	//compare "N samples" to "2N samples"; is the KL 
	//(Kullback-Leibler) divergence close to 0?
	//Q=approximation
	//P=truth
	//D_{KL}(P||Q)=sum_x P(x)log(P/Q)
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

    torch::Tensor get_covariance() const {
        auto P_output = P_;
        return P_output;
    }

private:
    int nx_, nz_;
    torch::Tensor F_,H_,Q_, R_, P_,x_,I_,L_Q_;
    torch::Device device_;
};



class EnsembleKalmanFilter {
public:
    using DynamicsFn = std::function<torch::Tensor(
	torch::Tensor&, /* u */
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
      const torch::Tensor& P0) {

      //P0=L_P0 * L_P0^*
     auto L_P0 = torch::linalg_cholesky(P0.to(device_));
     auto opts=torch::TensorOptions().dtype(torch::kFloat64).device(device_);
     auto noise = torch::randn({N_, nx_},opts); 
      // Ensemble shape: [N, nx]
     ensemble_ = x0.to(device_).unsqueeze(0) + torch::matmul(noise, L_P0.t());
    }

    // Forecast (Predict) Step
    void predict(double alpha,int n) {
        // Propagate each ensemble member through the nonlinear dynamics
        // f_ expects [N, nx] and returns [N, nx]
        ensemble_ = f_(ensemble_,alpha,n);

        // Add process noise: w ~ N(0, Q)
        auto opts=torch::TensorOptions().dtype(torch::kFloat64).device(device_);
        auto w = torch::matmul( torch::randn({N_, nx_},opts),L_Q_.t() );
        ensemble_ = ensemble_ + w;
    }

    // Analysis (Update) Step using Perturbed Observations
    void update(const torch::Tensor& z) {

        auto opts=torch::TensorOptions().dtype(torch::kFloat64).device(device_);
        auto z_dev = z.to(device_); // [nz]

        // 1. Map ensemble into observation space: [N, nz]
        auto Z_ens = h_(ensemble_);

        // 2. Compute ensemble means
        auto x_mean = ensemble_.mean(0, /*keepdim=*/true); // [1, nx]
        auto z_mean = Z_ens.mean(0, /*keepdim=*/true);     // [1, nz]

        // 3. Compute ensemble anomalies (deviations from the mean)
        auto A_x = ensemble_ - x_mean; // [N, nx]
        auto A_z = Z_ens - z_mean;     // [N, nz]

        // 4. Compute sample covariances: P_xz and P_zz
        double norm = 1.0 / static_cast<double>(N_ - 1);
        auto P_xz = norm * torch::matmul(A_x.t(), A_z);//[nx, nz]
        auto P_zz = norm * torch::matmul(A_z.t(), A_z) + R_;//[nz, nz]

        // 5. Compute Kalman Gain: K = P_xz * P_zz^{-1}
        // Using torch::linalg_solve on transposes avoids explicit matrix inversion:
        // P_zz * K^T = P_xz^T  =>  K = (P_zz^{-1} * P_xz^T)^T
        auto K = torch::linalg_solve(P_zz, P_xz.t(), /*left=*/true).t(); // [nx, nz]

        // 6. Generate perturbed observations for each ensemble member
        auto v = torch::matmul( torch::randn({N_, nz_},opts), L_R_.t());
        auto Z_perturbed = z_dev.unsqueeze(0) + v; // [N, nz]

        // 7. Update ensemble members
        auto innovation = Z_perturbed - Z_ens; // [N, nz]
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
    torch::Device device_;
};

// Computes KL(N(mu1, cov1) || N(mu2, cov2))
// mu1,cov1 is "truth"
// mu1, mu2: [batch_size, k]
// cov1, cov2: [batch_size, k, k] (must be symmetric positive-definite)
torch::Tensor kl_divergence_multivariate_gaussian(
    const torch::Tensor& mu1, const torch::Tensor& cov1,
    const torch::Tensor& mu2, const torch::Tensor& cov2) {
    
    int64_t k = mu1.size(-1);
    
    // Cholesky factors: cov = L * L^T
    auto L1 = torch::linalg_cholesky(cov1);
    auto L2 = torch::linalg_cholesky(cov2);
    
    // 1. Log determinant ratio: ln(|cov2| / |cov1|) = 2 * sum(log(diag(L2))) - 2 * sum(log(diag(L1)))
    auto logdet1 = 2.0 * L1.diagonal(0, -2, -1).log().sum(-1);
    auto logdet2 = 2.0 * L2.diagonal(0, -2, -1).log().sum(-1);
    auto logdet_ratio = logdet2 - logdet1;
    
    // 2. Trace term: tr(cov2^-1 * cov1) = || L2^-1 * L1 ||_F^2
    // torch::linalg_solve_triangular solves L2 * X=L1  for X
    auto L2_inv_L1 = torch::linalg_solve_triangular(L2, L1, /*upper=*/false);
    auto trace_term = L2_inv_L1.pow(2).sum({-2, -1});
    
    // 3. Quadratic term: (mu2 - mu1)^T * cov2^-1 * (mu2 - mu1) = || L2^-1 * (mu2 - mu1) ||_2^2
    auto diff = (mu2 - mu1).unsqueeze(-1); // [batch_size, k, 1]
    auto L2_inv_diff = torch::linalg_solve_triangular(L2, diff,/*upper=*/false);
    auto quad_term = L2_inv_diff.pow(2).sum({-2, -1});
    
    return 0.5 * (trace_term + quad_term - k + logdet_ratio);
}



int main() {
    torch::manual_seed(42);
    torch::NoGradGuard no_grad; // Disable autograd overhead for standard filtering

    // Select GPU if available
    torch::Device device(torch::cuda::is_available() ? torch::kCUDA : torch::kCPU);
    std::cout << "Running EnKF on: " << (device.is_cuda() ? "CUDA" : "CPU") << "\n\n";

    const int N = 2;   // Ensemble size
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

     auto u_next = u.clone();

     // Ensemble shape: [N, nx]
     u_next=(1.0-alpha)*u+alpha*torch::roll(u,-1,1); //shift up by 1
     u.copy_(u_next);
     return u_next;
    };

    // Observation model: z_k = h(x_k) -> observes [x, y]
    auto obs_model = [](torch::Tensor& x) {
     auto x_obs = x.clone();
     return x_obs;
    };

    double noise_variance=0.01;
    double noise_sd=std::sqrt(noise_variance);

    double process_noise_variance=0.5;
    double process_noise_sd=std::sqrt(process_noise_variance);

    auto opts=torch::TensorOptions().dtype(torch::kFloat64).device(device);

    // Noise covariances
    auto Q = torch::eye(Nx,opts) * process_noise_variance; // Process noise
    auto R = torch::eye(Nz,opts) * noise_variance; // Measurement noise

    auto H = torch::eye(Nz,opts);
    auto F = torch::eye(Nx,opts);
    for (int i=0;i<Nx;i++) {
     F.index_put_({i,i},1.0-alpha);
     if (i==0) {
      F.index_put_({i,Nx-2},alpha);
     } else {
      F.index_put_({i,i-1},alpha);
     }
    }

    auto P0=torch::eye(Nx, opts)*process_noise_variance; //P0 noise
    auto u0 = torch::sin(2.0 * M_PI * x / L);

    std::cout << "before kf declare\n";
    KalmanFilter kf(Nx, Nz, F,H, 
		Q, R, P0, u0, device);
    kf.initialize(u0, P0);

    std::cout << "before enkf declare\n";
    EnsembleKalmanFilter enkf(N, Nx, Nz, dynamics_model, obs_model, 
		Q, R, device);

    auto u_true = torch::sin(2.0 * M_PI * x / L);

    auto u0_guess = torch::sin(2.0 * M_PI * x / L)+
	    torch::randn({Nx},opts)*process_noise_sd;  // u0 noise

    std::cout << "before enkf.init\n";
    enkf.initialize(u0_guess, P0);

    std::cout << "before loop\n";

    // Run simulation loop
    for (int n = 0; n < Nt; ++n) {

     u_true = torch::matmul(F, u_true);

     auto opts=torch::TensorOptions().dtype(torch::kFloat64).device(device);
     auto z = u_true + torch::randn({Nz}, opts) * noise_sd;
     auto z_enkf = torch::matmul(H, z);

     enkf.predict(alpha,n);
     enkf.update(z_enkf);

     auto truth = u_true.cpu();
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


     auto u_pred_kf=kf.predict(alpha,n);
     auto est_kf=kf.update(z);
     auto est_cov_kf=kf.get_covariance().cpu();
     auto abs_error_kf=est_kf-truth;
     auto abs_error_norm_kf = torch::mean(torch::square(abs_error_kf));
     abs_error_norm_kf=torch::sqrt(abs_error_norm_kf);
     std::cout << "n " << n 
            << " | abs_error_norm_kf: " << abs_error_norm_kf << '\n';


     auto est = enkf.get_mean().cpu();
     auto est_cov=enkf.get_covariance().cpu();
     auto abs_error=est-truth;
     auto abs_error_norm = torch::mean(torch::square(abs_error));
     abs_error_norm=torch::sqrt(abs_error_norm);

     std::cout << "n " << n 
            << " | abs_error_norm: " << abs_error_norm << '\n';

     if (N>Nx) {
      auto kl_div=kl_divergence_multivariate_gaussian(
       est_kf,est_cov_kf,est,est_cov);

      std::cout << "n " << n 
            << " | kl_div: " << kl_div << '\n';
     }


     if (1==0) {
      auto kl_div_sanity1=kl_divergence_multivariate_gaussian(
       est_kf,est_cov_kf,est_kf,est_cov_kf);
      auto kl_div_sanity2=kl_divergence_multivariate_gaussian(
       est,est_cov,est,est_cov);
      auto kl_div_sanity3=kl_divergence_multivariate_gaussian(
       est,est_cov,est_kf,est_cov_kf);
      std::cout << "n " << n 
             << " | kl_div_sanity1: " << kl_div_sanity1 << '\n';
      std::cout << "n " << n 
             << " | kl_div_sanity2: " << kl_div_sanity2 << '\n';
      std::cout << "n " << n 
             << " | kl_div_sanity3: " << kl_div_sanity3 << '\n';
     }

    }

    return 0;
}

