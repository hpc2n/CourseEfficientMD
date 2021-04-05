#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -J namd
#SBATCH -t 00:10:00
#SBATCH -n 14
#SBATCH --gres=gpu:k80:1
#SBATCH --exclusive

echo $CUDA_VISIBLE_DEVICES

#Load modules necessary for running NAMD
ml purge  > /dev/null 2>&1 
ml GCC/5.4.0-2.26  CUDA/8.0.61_375.26  OpenMPI/2.0.2
ml NAMD/2.12-nompi

#Using all GPUs
namd2 +p14 +setcpuaffinity +idlepoll +devices 0 step4_equilibration.inp > output_gpu_offload.dat

exit 0

#Running on scratch directory
rm -rf /scratch/test
export parent=/pfs/nobackup/home/p/pojedama/benchmarks/charmm-gui/namd
mkdir /scratch/test
rsync -avzh $parent/ /scratch/test/

cd /scratch/test 

namd2 +p28 +setcpuaffinity step4_equilibration.inp > output_gpu.dat

cd $parent
rsync -avzh /scratch/test/ $parent/

rm -rf /scratch/test

