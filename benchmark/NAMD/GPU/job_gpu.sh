#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -J namd
#SBATCH -t 00:08:00
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --gres=gpu:k80:2
#For V100 uncomment the following line and comment the line above
##SBATCH --gres=gpu:v100:2
#SBATCH --exclusive


#Load modules necessary for running NAMD
ml purge  > /dev/null 2>&1 
ml GCC/9.3.0  CUDA/11.0.2  OpenMPI/4.0.3
ml NAMD/2.14-nompi 

namd2 +p28 +setcpuaffinity +idlepoll +devices $CUDA_VISIBLE_DEVICES step4_equilibration.inp > output_gpu1.dat
namd2 +p28 +setcpuaffinity +idlepoll +devices $CUDA_VISIBLE_DEVICES step4_equilibration_mts.inp > output_gpu2.dat

