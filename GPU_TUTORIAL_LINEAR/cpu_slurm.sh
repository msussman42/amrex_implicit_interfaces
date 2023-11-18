#!/bin/bash
#SBATCH --job-name="cpu test"
#SBATCH --ntasks=1
#SBATCH --partition=quicktest
#SBATCH --mail-type="ALL"
#SBATCH --mail-user=msussman@fsu.edu
#SBATCH -t 00:04:00
#SBATCH --output=cpurun.out
#SBATCH --error=cpurun.err

# other partitions:
# backfill, backfill2, genacc_q, quicktest (10 minutes)
# mecfd18_q
# engineering_q
# engineering_long
# compiling options:
# module load gnu openmpi
# module load intel openmpi
# module load cuda
# sussman@hpc-login.rcc.fsu.edu

# To submit the job 
#   $ chmod +x NAME_Of_SCRIPT_FILE
#   $ sbatch NAME_Of_SCRIPT_FILE 

# To check the running jobs
#   $ squeue -u YOUR_USERNAME
#   $ squeue -u YOUR_USERNAME -o "%.18i %.14P %.20j %.2t %.10M %.4D %.4C %R" -S
#     "t,P,j"

# Cancel job
#   $ scancel JOB_ID

# Information about nodes on a partition
#   $ sinfo -p PARTTION_NAME
#   $ sinfo -p PARTTION_NAME -N -o "%.12N %.4D %.11T %.13C %.8z %.6m %.8f %10E"


# Check your partitions 
#   $ rcctool my:partitions

# More resourecs at: https://rcc.fsu.edu/doc/ [rcc.fsu.edu]
echo module load cuda-12.2
module load cuda-12.2
module list 
echo nvcc --version
nvcc --version
echo /usr/bin/nvidia-smi -L
#/usr/bin/nvidia-smi -L
echo /usr/bin/nvidia-smi --query-gpu=gpu_name,gpu_bus_id,vbios_version --format=csv
#/usr/bin/nvidia-smi --query-gpu=gpu_name,gpu_bus_id,vbios_version --format=csv
echo /usr/bin/nvidia-smi --query-gpu=timestamp,name,pci.bus_id,driver_version --format=csv
#/usr/bin/nvidia-smi --query-gpu=timestamp,name,pci.bus_id,driver_version --format=csv
gcc --version 
srun ~/CNSWAVE/CNS_CPU inputs

