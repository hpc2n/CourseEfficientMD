#!/bin/bash
#SBATCH -A staff
#Asking for 50 min.
#SBATCH -t 00:50:00
#Ask for 28 processes
#SBATCH -n 28
#SBATCH --output=job_o.out
#SBATCH --error=job_o.err

#ml gaussian/09.e.01-AVX
#ml GCC/7.3.0-2.30  OpenMPI/3.1.1

ml gaussian/16.C.01-AVX2
ml GCC/9.3.0  OpenMPI/4.0.3
ml SciPy-bundle/2020.03-Python-2.7.18
ml NAMD/2.14-mpi

#ml Python/2.7.15
#ml GCCcore/9.3.0
#ml Python/2.7.18
#ml NAMD/2.13-mpi

mpirun -np 28 namd2 4ake_eq_qmmm.conf > logfile_qmmm.txt
