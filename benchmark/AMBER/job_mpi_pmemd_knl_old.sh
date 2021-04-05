#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -J amber
#SBATCH -t 00:50:00
#SBATCH -N 1
#SBATCH -n 68
#SBATCH -p knl
##SBATCH -c 1
##SBATCH --ntasks-per-core=4
#SBATCH --constraint=cache,quad

ml --force purge


module load icc/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196 ifort/2017.4.196-GCC-6.4.0-2.28  
module load Amber/16-AmberTools-17-patchlevel-0-7

#ml GCC/7.3.0-2.30  OpenMPI/3.1.1
#ml Amber/18-AmberTools-18-patchlevel-10-8-titr_res-100



####srun sander.MPI -ng 8 -i step4.0_minimization.mdin
export init="step3_charmm2amber"
#export istep="step4.0_minimization"
#srun sander.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${init}.rst7 -o ${istep}.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7

# Equilibration
# step4.1
export pstep="step4.0_minimization"
export istep="step4.1_equilibration"

#export I_MPI_PIN_MODE=pm
#export I_MPI_PIN_DOMAIN=auto 
export PHI_OMP_NUM_THREADS=2
#export KMP_AFFINITY=compact
#export KMP_STACKSIZE=10M

srun pmemd.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${pstep}.rst7 -o ${istep}.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7 -x ${istep}.nc
#srun pmemd.cuda.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${pstep}.rst7 -o ${istep}.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7 -x ${istep}.nc




exit 0
