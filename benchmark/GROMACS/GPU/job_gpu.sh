#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -t 00:20:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 7
#K80 cards
#SBATCH --gres=gpu:k80:2
#V100 cards
##SBATCH --gres=gpu:v100:2
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --exclusive

ml purge  > /dev/null 2>&1 
ml GCC/8.3.0  CUDA/10.1.243  OpenMPI/3.1.4
ml GROMACS/2021-Python-3.7.4

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export reset_counters="-resetstep 10000 -nsteps 20000"

gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_charmm2gmx.pdb -n index.ndx -p topol.top

#Running with default options (MPI-Threaded version)
gmx mdrun -ntomp $SLURM_CPUS_PER_TASK -ntmpi $SLURM_NTASKS $reset_counters -dlb yes  -v -deffnm step4.1_equilibration

#Using 1 rank for PME and offloading NB and PME to the GPUs (MPI-Threaded version)
gmx mdrun -nb gpu -pme gpu -npme 1 -ntomp $SLURM_CPUS_PER_TASK -ntmpi $SLURM_NTASKS $reset_counters -dlb yes  -v -deffnm step4.1_equilibration

#Running with default options (MPI version)
srun gmx_mpi mdrun -ntomp $SLURM_CPUS_PER_TASK $reset_counters -dlb yes  -v -deffnm step4.1_equilibration

#Using 1 rank for PME and offloading NB and PME to the GPUs (MPI version)
srun gmx_mpi mdrun -nb gpu -pme gpu -npme 1 -ntomp $SLURM_CPUS_PER_TASK  $reset_counters -dlb yes  -v -deffnm step4.1_equilibration

