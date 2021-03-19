#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -J namd
#SBATCH -t 00:20:00
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --exclusive
#SBATCH --constraint=skylake
#Load modules necessary for running NAMD
module add icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132
module add NAMD/2.12-nompi

namd2 +p 28 +setcpuaffinity step4_equilibration.inp > output.dat


exit 0

#Running on scratch directory
rm -rf /scratch/test
export parent=/pfs/nobackup/home/p/pojedama/benchmarks/charmm-gui/namd
mkdir /scratch/test
rsync -avzh $parent/ /scratch/test/

cd /scratch/test 

namd2 +p 28 +setcpuaffinity step4_equilibration.inp > output_smp_scratch.dat

cd $parent
rsync -avzh /scratch/test/ $parent/

rm -rf /scratch/test

