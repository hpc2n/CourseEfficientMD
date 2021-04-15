#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -J namd
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -n 28
#If you request a Skylake node uncomment the following two lines
##SBATCH --constraint=skylake
##SBATCH --exclusive

#Load modules necessary for running NAMD
ml purge  > /dev/null 2>&1 
ml GCC/9.3.0  OpenMPI/4.0.3
ml NAMD/2.14-mpi 

srun namd2 step4_equilibration.inp > output_mpi1.dat
srun namd2 step4_equilibration_mts.inp > output_mpi2.dat

