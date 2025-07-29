#!/bin/bash
#SBATCH --job-name="1_150"
##SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH -C YEAR2022,amd
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --partition=engineering_q
#SBATCH --mail-type="ALL"
#SBATCH --mail-user=msussman@fsu.edu
#SBATCH -t 48:00:00
#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --open-mode=append
#SBATCH --exclusive



# Information about nodes on a partition
#   $ sinfo -p PARTTION_NAME
#   $ sinfo -p PARTTION_NAME -N -o "%.12N %.4D %.11T %.13C %.8z %.6m %.8f %10E"


# Check your partitions 
#   $ rcctool my:partitions

pwd;hostname;date
echo "running amrMPI (FABRIC) on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS tasks, each with $SLURM_CPUS_PER_TASK cores."

 module purge
 module load gnu/11.2.1
 module load openmpi
# # module load intel-openmpi

srun /gpfs/research/engineering/Kshoele/FabricDrop_New_20240120/amrex_implicit_interfaces/NS_fluids_lib/amr3d.gnu.FLOAT.MPI.ex inputs.FABRIC_DROP

