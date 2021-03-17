#!/bin/bash
#SBATCH -A staff
#SBATCH -J namd
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -n 28

#Load modules necessary for running NAMD
ml purge  > /dev/null 2>&1 
module add GCC/6.3.0-2.27  OpenMPI/2.0.2
module add NAMD/2.12-mpi

srun namd2 step4_equilibration.inp > output_mpi1.dat
srun namd2 step4_equilibration_mts.inp > output_mpi2.dat

#Load modules necessary for running NAMD
ml purge  > /dev/null 2>&1 
ml GCC/7.3.0-2.30  OpenMPI/3.1.1
ml NAMD/2.13-mpi

srun namd2 step4_equilibration.inp > output_mpi3.dat

