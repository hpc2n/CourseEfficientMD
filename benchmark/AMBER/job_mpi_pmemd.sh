#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -t 00:25:00
#SBATCH -n 28
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --mail-type=END

ml purge  > /dev/null 2>&1 

ml GCC/7.3.0-2.30  CUDA/9.2.88  OpenMPI/3.1.1
ml Amber/18-AmberTools-18-patchlevel-10-8 



####srun sander.MPI -ng 8 -i step4.0_minimization.mdin
export init="step3_charmm2amber"
#export istep="step4.0_minimization"
#srun sander.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${init}.rst7 -o ${istep}.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7

# Equilibration
# step4.1
export pstep="step4.0_minimization"
export istep="step4.1_equilibration"

srun pmemd.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${pstep}.rst7 -o ${istep}.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7 -x ${istep}.nc

