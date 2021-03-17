#!/bin/bash
#SBATCH -A staff
#SBATCH -t 00:50:00
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --gres=gpu:v100:2
#SBATCH --exclusive
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --mail-type=END

ml purge
ml icc/2017.1.132-GCC-5.4.0-2.26 ifort/2017.1.132-GCC-5.4.0-2.26 CUDA/8.0.44  impi/2017.1.132
ml OpenMM/7.1.1-Python-3.6.1
          

# Equilibration
# step4

export init=step3_charmm2omm
export istep=step4_equilibration

python -u openmm_run.py -i ${istep}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -b ${init}.str -orst ${istep}.rst -odcd ${istep}.dcd > ${istep}.out

exit 0

