#include <torch/torch.h>
#include <iostream>

class KalmanFilter {
public:
    KalmanFilter(
        torch::Tensor F, // State transition matrix
        torch::Tensor H, // Measurement matrix
        torch::Tensor Q, // Process noise covariance
        torch::Tensor R, // Measurement noise covariance
        torch::Tensor P, // Initial estimate error covariance
        torch::Tensor x0 // Initial state estimate
    ) : F_(std::move(F)), H_(std::move(H)), Q_(std::move(Q)),
        R_(std::move(R)), P_(std::move(P)), x_(std::move(x0)) {
        I_ = torch::eye(x_.size(0), x_.options());
    }

    // Step 1: Predict state and error covariance
    torch::Tensor predict() {
        // x = F * x
        x_ = torch::matmul(F_, x_);
        // P = F * P * F^T + Q
        P_ = torch::matmul(torch::matmul(F_, P_), F_.transpose(0, 1)) + Q_;
        return x_;
    }

    // Step 2: Update state with new measurement z
    torch::Tensor update(const torch::Tensor& z) {
        // Innovation (residual): y = z - H * x
        torch::Tensor y = z - torch::matmul(H_, x_);

        // Innovation covariance: S = H * P * H^T + R
        torch::Tensor PHt = torch::matmul(P_, H_.transpose(0, 1));
        torch::Tensor S = torch::matmul(H_, PHt) + R_;

        // Optimal Kalman Gain: K = P * H^T * S^-1
        // Using linalg_solve (S * K^T = (P * H^T)^T) is more numerically stable than matrix inverse
        torch::Tensor K = torch::linalg_solve(S, PHt.transpose(0, 1)).transpose(0, 1);

        // Updated state estimate: x = x + K * y
        x_ = x_ + torch::matmul(K, y);

        // Updated estimate covariance (Joseph form for numerical stability):
        // P = (I - K * H) * P * (I - K * H)^T + K * R * K^T
        torch::Tensor I_KH = I_ - torch::matmul(K, H_);
        P_ = torch::matmul(torch::matmul(I_KH, P_), I_KH.transpose(0, 1))
           + torch::matmul(torch::matmul(K, R_), K.transpose(0, 1));

        return x_;
    }

private:
    torch::Tensor F_, H_, Q_, R_, P_, x_, I_;
};

int main() {
    // Use double precision (float64) for better numerical stability in filtering
    auto opts = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);

    double dt = 0.1; // Time step

    // State vector: [position, velocity]^T
    torch::Tensor x0 = torch::tensor({{0.0}, {1.0}}, opts);

    // State transition matrix (constant velocity model)
    torch::Tensor F = torch::tensor({
        {1.0, dt},
        {0.0, 1.0}
    }, opts);

    // Measurement matrix (we only observe position)
    torch::Tensor H = torch::tensor({{1.0, 0.0}}, opts);

    // Process noise covariance
    torch::Tensor Q = torch::eye(2, opts) * 0.01;

    // Measurement noise covariance
    torch::Tensor R = torch::tensor({{0.5}}, opts);

    // Initial error covariance
    torch::Tensor P = torch::eye(2, opts) * 1.0;

    KalmanFilter kf(F, H, Q, R, P, x0);

    // Simulate noisy position measurements
    std::cout << "Step | Measured Pos | Estimated Pos | Estimated Vel\n";
    std::cout << "-----------------------------------------------------\n";

    for (int step = 1; step <= 10; ++step) {
        double true_pos = step * dt * 1.0;
        double noisy_pos = true_pos + torch::randn({1}, opts).item<double>() * 0.5;
        torch::Tensor z = torch::tensor({{noisy_pos}}, opts);

        kf.predict();
        torch::Tensor state = kf.update(z);

        printf("%4d | %12.4f | %13.4f | %13.4f\n",
               step,
               noisy_pos,
               state[0][0].item<double>(),
               state[1][0].item<double>());
    }

    return 0;
}

