#!/bin/bash
#SBATCH -A Project_ID
#Asking for 50 min.
#SBATCH -t 00:50:00
#Ask for 28 processes
#SBATCH -n 28
#SBATCH --output=job_o.out
#SBATCH --error=job_o.err

ml gaussian/16.C.01-AVX2
ml GCC/9.3.0  OpenMPI/4.0.3
ml SciPy-bundle/2020.03-Python-2.7.18
ml NAMD/2.14-mpi

mpirun -np 28 namd2 4ake_eq_qmmm.conf > logfile_qmmm.txt
