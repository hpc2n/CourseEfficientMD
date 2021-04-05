#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -t 02:15:00
#SBATCH -N 2
#SBATCH -n 130
#SBATCH -p knl
#SBATCH --constraint=cache,quad
###SBATCH --ntasks-per-node=68
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --mail-type=END

#ml GCC/7.3.0-2.30 OpenMPI/3.1.1
#ml GROMACS/2018.3

#ml icc/2018.3.222-GCC-7.3.0-2.30 ifort/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
#ml GCC/7.3.0-2.30  OpenMPI/3.1.1
#ml GROMACS/2018.3

#ml GCC/7.3.0-2.30  OpenMPI/3.1.1
#ml GROMACS/2019

ml icc/2018.3.222-GCC-7.3.0-2.30 ifort/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
ml GROMACS/2019


# Automatic selection of options to mdrun_mpi depending on parameters given to SBATCH
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    md="-ntomp $SLURM_CPUS_PER_TASK"
else
    md="-ntomp 1"
fi


#tpr file taken from: www.mpibpc.mpg.de/grubmueller/bench 

numactl -H

#Threaded-MPI version
export OMP_NUM_THREADS=2
gmx mdrun -ntmpi 68 -npme 18 -ntomp 2 -pin on -pinoffset 0 -pinstride 2 -dlb yes -resetstep 5000 -nsteps 10000 -s benchRIB.tpr 

#MPI version
export OMP_NUM_THREADS=2
mpirun -np 130 gmx_mpi mdrun -ntomp 2 -pin on -pinoffset 0 -pinstride 2 -resetstep 5000 -nsteps 10000 -dlb auto -s benchRIB.tpr 



exit 0
