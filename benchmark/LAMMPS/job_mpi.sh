#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -t 02:10:00
#SBATCH -n 28


#Load modules necessary for running LAMMPS
ml GCC/8.3.0  OpenMPI/3.1.4
ml LAMMPS/3Mar2020-Python-3.7.4-kokkos

#Execute LAMMPS
#srun lmp -in step4.0_minimization.inp 
srun lmp -in step4.1_equilibration.inp
