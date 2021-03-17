#!/bin/bash
#SBATCH -A staff
#SBATCH -J namd
#SBATCH -t 01:50:00
#SBATCH -N 1
#SBATCH -p knl
##SBATCH --constraint=cache,hemi
#SBATCH --constraint=quad,flat


#This one is working with multicore
#/home/p/pojedama/pfs/NAMD_Git-2017-11-16_Source/Linux-KNL-icc/namd2 +setcpuaffinity +ppn 126 step4_equilibration.inp > output_knl.dat
#/home/p/pojedama/pfs/NAMD_Git-2017-11-16_Source/Linux-KNL-icc/namd2 +ppn 126 step4_equilibration.inp > output_knl.dat


#I used it for  the course Feb19
/home/p/pojedama/pfs/NAMD_Git-2017-11-16_Source/Linux-KNL-icc/namd2 +setcpuaffinity +ppn 272 step4_equilibration.inp > output_knl.dat

numactl -H
numactl -m 1 /home/p/pojedama/pfs/NAMD_Git-2017-11-16_Source/Linux-KNL-icc/namd2 +setcpuaffinity +ppn 272 step4_equilibration.inp > output_knl2.dat 

ml purge
ml intel/2019a
/pfs/nobackup/home/p/pojedama/NAMD_2.12/NAMD_2.12_Source/Linux-x86_64-icc/namd2 +setcpuaffinity +ppn 272 step4_equilibration-tcl.inp > output_knl3.dat

numactl -m 1 /pfs/nobackup/home/p/pojedama/NAMD_2.12/NAMD_2.12_Source/Linux-x86_64-icc/namd2 +setcpuaffinity +ppn 272 step4_equilibration-tcl.inp > output_knl4.dat

