#!/bin/bash
#SBATCH -A hpc2n2021-001
#SBATCH -t 00:25:00
#SBATCH -n 28
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --mail-type=END

ml purge  > /dev/null 2>&1 

ml icc/2017.1.132-GCC-5.4.0-2.26 ifort/2017.1.132-GCC-5.4.0-2.26  CUDA/8.0.44  impi/2017.1.132
#ml ifort/2017.1.132-GCC-5.4.0-2.26  CUDA/8.0.44  impi/2017.1.132

ml Amber/16-AmberTools-16-patchlevel-20-7-hpc2n 


####srun sander.MPI -ng 8 -i step4.0_minimization.mdin
export init="step3_charmm2amber"
#export istep="step4.0_minimization"
#srun sander.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${init}.rst7 -o ${istep}.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7

# Equilibration
# step4.1
export pstep="step4.0_minimization"
export istep="step4.1_equilibration"

srun sander.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${pstep}.rst7 -o ${istep}.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7 -x ${istep}.nc

