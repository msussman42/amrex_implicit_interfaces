#pytorch.org
#PyTorch Build Stable (2.14.1)
#Your OS Linux (HPC)
#Package Pip
#Language Python
#Compute Platform CPU
#pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cpu
#python3 rand53.py
import torch
x = torch.rand(5, 3)
print(x)
