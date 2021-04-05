#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -J namd
#SBATCH -t 00:08:00
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --gres=gpu:k80:2
#SBATCH --exclusive

echo $CUDA_VISIBLE_DEVICES


#Load modules necessary for running NAMD
ml purge  > /dev/null 2>&1 
ml GCC/5.4.0-2.26  CUDA/8.0.61_375.26  OpenMPI/2.0.2
ml NAMD/2.12-nompi

namd2 +p28 +setcpuaffinity +idlepoll +devices $CUDA_VISIBLE_DEVICES step4_equilibration.inp > output_gpu1.dat
namd2 +p28 +setcpuaffinity +idlepoll +devices $CUDA_VISIBLE_DEVICES step4_equilibration_mts.inp > output_gpu2.dat

exit 0
