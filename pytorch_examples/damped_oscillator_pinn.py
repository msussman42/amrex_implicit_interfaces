#pytorch.org
#PyTorch Build Stable (2.14.1)
#Your OS Linux (HPC)
#Package Pip
#Language Python
#Compute Platform CPU
#pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cpu
#python3 damped_oscillator_pinn.py
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 1. Define the Neural Network Architecture
class PINN(nn.Module):
    def __init__(self):
        super().__init__()
        # Use Tanh activation because we need continuous second derivatives
        self.net = nn.Sequential(
            nn.Linear(1, 32),
            nn.Tanh(),
            nn.Linear(32, 32),
            nn.Tanh(),
            nn.Linear(32, 1)
        )
        
    def forward(self, t):
        return self.net(t)

# 2. Define Physics Parameters (Damped Oscillator)
mu, k = 2.0, 20.0  # mu = damping coefficient, k = stiffness

# 3. Generate Training Data (Collocation Points)
# We don't need real data points, just places to evaluate the physics equation
t_physics = torch.linspace(0, 1, 500).view(-1, 1).requires_grad_(True)

# Boundary/Initial Conditions (at t = 0, u = 1 and du/dt = 0)
t_boundary = torch.tensor([[0.0]], requires_grad=True)

# 4. Initialize Network and Optimizer
model = PINN()
optimizer = optim.Adam(model.parameters(), lr=1e-3)

# 5. Training Loop
for epoch in range(3001):
    optimizer.zero_grad()
    
    # --- Compute Boundary/Initial Condition Loss ---
    u_b = model(t_boundary)
    # Derivative at boundary
    du_dt_b = torch.autograd.grad(u_b, t_boundary, torch.ones_like(u_b), create_graph=True)[0]
    
    loss_b1 = torch.mean((u_b - 1.0) ** 2)      # Initial displacement = 1
    loss_b2 = torch.mean((du_dt_b - 0.0) ** 2)  # Initial velocity = 0
    loss_boundary = loss_b1 + loss_b2

    # --- Compute Physics Residual Loss ---
    u = model(t_physics)
    
    # First derivative (du/dt)
    du_dt = torch.autograd.grad(u, t_physics, torch.ones_like(u), create_graph=True)[0]
    # Second derivative (d^2u/dt^2)
    d2u_dt2 = torch.autograd.grad(du_dt, t_physics, torch.ones_like(du_dt), create_graph=True)[0]
    
    # Governing ODE residual: must equal 0
    physics_residual = d2u_dt2 + mu * du_dt + k * u
    loss_physics = torch.mean(physics_residual ** 2)
    
    # --- Total Loss and Backpropagation ---
    # Scale boundary loss slightly higher to strictly enforce constraints
    total_loss = 2.0 * loss_boundary + loss_physics
    
    total_loss.backward()
    optimizer.step()
    
    if epoch % 500 == 0:
        print(f"Epoch {epoch:4d} | Total Loss: {total_loss.item():.6f} | "
              f"Boundary Loss: {loss_boundary.item():.6f} | Physics Loss: {loss_physics.item():.6f}")

# 6. Verify Results against Analytical Solution
def analytical_solution(t, mu, k):
    # Underdamped system exact solution
    w0 = np.sqrt(k)
    alpha = mu / 2
    w = np.sqrt(w0**2 - alpha**2)
    return np.exp(-alpha * t) * (np.cos(w * t) + (alpha / w) * np.sin(w * t))

# Plotting the evaluation
t_test = np.linspace(0, 1, 100)
t_test_tensor = torch.tensor(t_test, dtype=torch.float32).view(-1, 1)

with torch.no_grad():
    u_pred = model(t_test_tensor).numpy()
u_true = analytical_solution(t_test, mu, k)

plt.figure(figsize=(8, 4))
plt.plot(t_test, u_true, label='Exact Analytical', color='black', lw=2)
plt.plot(t_test, u_pred, '--', label='PINN Prediction', color='red', lw=2)
plt.title('PINN vs Analytical Solution (Damped Oscillator)')
plt.xlabel('Time (t)')
plt.ylabel('Displacement (u)')
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig('my_plot.png',bbox_inches='tight',dpi=300)
plt.close()
