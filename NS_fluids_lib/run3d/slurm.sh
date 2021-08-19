#!/bin/bash
#SBATCH --job-name="fabric_drop_problem"
#SBATCH --ntasks=4
#SBATCH --partition=backfill
#SBATCH --mail-type="ALL"
#SBATCH --mail-user=msussman@fsu.edu
#SBATCH -t 03:59:00
#SBATCH --output=run.out
#SBATCH --error=run.err

# other partitions:
# backfill, backfill2, genacc_q
# mecfd18_q
# engineering_q
# engineering_long
# compiling options:
# module load gnu openmpi
# module load intel openmpi
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
pwd;hostname;date
echo "running amrMPI (FABRIC) on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS tasks, each with $SLURM_CPUS_PER_TASK cores."
module purge
module load gnu openmpi
#module load intel openmpi
srun ~/FABRIC_DROP/amrMPI inputs.FABRIC_DROP

