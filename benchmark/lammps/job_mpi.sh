#!/bin/bash
#SBATCH -A staff
#SBATCH -t 02:10:00
#SBATCH -N 1
#SBATCH -n 28


#Load modules necessary for running LAMMPS
module load icc/2017.1.132-GCC-6.3.0-2.27 ifort/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load LAMMPS/31Mar17

#Execute LAMMPS
#srun lmp_intel_cpu -in step4.0_minimization.inp 
srun lmp_intel_cpu -in step4.1_equilibration.inp
