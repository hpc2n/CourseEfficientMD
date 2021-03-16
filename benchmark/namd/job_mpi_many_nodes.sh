#!/bin/bash
#SBATCH -A SNIC2020-9-25
#SBATCH -J namd
#SBATCH -t 00:30:00
#SBATCH -N 2
#SBATCH -n 56
#SBATCH -p batch

module add GCC/6.3.0-2.27  OpenMPI/2.0.2
module add NAMD/2.12-mpi

srun namd2 +setcpuaffinity step4_equilibration.inp > output_mpi.dat

exit 0
