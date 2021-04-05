#!/bin/bash
#SBATCH -A SNICxxxx-yy-zz
#SBATCH -t 00:05:00
#SBATCH -n 4
#SBATCH -c 7
#For Skylake nodes uncomment the line below
###SBATCH --constraint=skylake
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err

ml purge  > /dev/null 2>&1 
ml GCC/8.3.0  CUDA/10.1.243  OpenMPI/3.1.4
ml GROMACS/2021-Python-3.7.4

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export reset_counters="-resetstep 1000 -nsteps 2000"

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    mdargs="-ntomp $SLURM_CPUS_PER_TASK"
else
    mdargs="-ntomp 1"
fi

gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_charmm2gmx.pdb -n index.ndx -p topol.top

#MPI-threaded version
gmx mdrun $reset_counters -ntmpi 4 -ntomp $SLURM_CPUS_PER_TASK -deffnm step4.1_equilibration

#MPI version
srun gmx_mpi mdrun $mdargs $reset_counters -npme 0 -dlb yes  -v -deffnm step4.1_equilibration


