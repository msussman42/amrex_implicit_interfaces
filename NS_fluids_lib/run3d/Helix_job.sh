#!/bin/bash
#SBATCH --job-name="visc_F_1-10"
##SBATCH --ntasks=28                # Number of MPI tasks (i.e. processes)
##SBATCH --cpus-per-task=1            # Number of cores per MPI task 
#SBATCH  -N 1
#SBATCH  -n 28
##SBATCH --ntasks-per-node=28         # Maximum number of tasks on each node
##SBATCH --ntasks-per-socket=14        # Maximum number of tasks on each socket
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --partition=engineering_long
#SBATCH --mail-type="END,FAIL"
#SBATCH -t 1:00:00

echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

module purge
module load gnu openmpi
cd /gpfs/research/mecfd/Kshoele/AFTERMOVE_FILES/testvis/SystematicRunsX4/RUNFOLDER_1
srun --hint=nomultithread -n 28 --cpu_bind=cores --exclusive /gpfs/research/mecfd/Kshoele/AFTERMOVE_FILES/ViscoelasticCode/3_23_2021/amrex_implicit_interfaces-master/NS_fluids_lib/amr3d.gnu.MPI.ex inputs.viscoHelix 1>run1.out 2>run1.err  &
wait
#mpirun --mca btl_openib_allow_ib 1 -np 8 ../../../amr3d.gnu.MPI.ex inputs.viscoHelix 1>run1.out 2>run1.err
