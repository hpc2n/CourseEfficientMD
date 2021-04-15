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

charmrun +p28 namd2 +setcpuaffinity step4_equilibration.inp > output_smp0.dat
charmrun +p28 namd2 step4_equilibration.inp > output_smp1.dat

