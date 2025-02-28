#!/bin/bash
#SBATCH -A Project_ID
#SBATCH -t 00:20:00
#SATCH -n 2
# Asking for Type X of GPU and 2 cards
#SBATCH --gpus-per-node=Type:2
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --mail-type=END

ml purge  > /dev/null 2>&1 
ml GCC/11.3.0  OpenMPI/4.1.4
ml Amber/22.4-AmberTools-22.5-CUDA-11.7.0

#Device information
nvidia-smi

export num_dev=`echo $CUDA_VISIBLE_DEVICES | awk 'BEGIN{FS=","};{print NF}'`

# Minimization step
export init="step3_charmm2amber"
#export istep="step4.0_minimization"
#srun sander.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${init}.rst7 -o ${istep}.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7

# Equilibration
# step4.1
export pstep="step4.0_minimization"
export istep="step4.1_equilibration"

mpirun -np 2 pmemd.cuda.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${pstep}.rst7 -o ${istep}1.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7 -x ${istep}.nc

