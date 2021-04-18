#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -t 00:50:00
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --gres=gpu:k80:2
#For v100 uncomment the following line and comment out the previous one
##SBATCH --gres=gpu:v100:2
#SBATCH --exclusive
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --mail-type=END

ml purge  > /dev/null 2>&1 
ml GCC/9.3.0  CUDA/11.0.2  OpenMPI/4.0.3
ml OpenMM/7.5.0-Python-3.8.2 
          

# Equilibration
# step4

export init=step3_charmm2omm
export istep=step4_equilibration

python -u openmm_run.py -i ${istep}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -b ${init}.str -orst ${istep}.rst -odcd ${istep}.dcd > ${istep}.out

exit 0

