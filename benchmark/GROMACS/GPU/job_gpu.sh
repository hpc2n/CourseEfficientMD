#!/bin/bash
# This is the job script to run a GROMACS example job 
# Remember to change the following to your own Project ID! 
#SBATCH -A Project_ID
#SBATCH -J Gromacs
#SBATCH -t 00:20:00
#SBATCH -n 4
# Asking for Type X of GPU
#SBATCH --gpus-per-node=Type:2

# It is a good idea to do a ml purge before loading other modules
ml purge > /dev/null 2>&1
ml GCC/12.3.0  OpenMPI/4.1.5
ml GROMACS/2023.3-CUDA-12.1.1-PLUMED-2.9.0

#Information about GPUs
nvidia-smi
echo $CUDA_VISIBLE_DEVICES 

gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_charmm2gmx.pdb -n index.ndx -p topol.top
sleep 10

#Three different ways to run this job:
#1. MPI version (Default)
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -dlb yes  -v -deffnm step4.1_equilibration
sleep 10

#2. MPI version (Offloading nb and pme to gpus)
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -nb gpu -pme gpu -npme 1 -dlb yes  -v -deffnm step4.1_equilibration
sleep 10

#3. Threaded-MPI version (Offloading nb and pme to gpus)
gmx mdrun -ntmpi $SLURM_NTASKS -nb gpu -pme gpu -npme 1  -dlb yes -v -deffnm step4.1_equilibration
