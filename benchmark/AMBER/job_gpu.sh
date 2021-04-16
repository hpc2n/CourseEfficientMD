#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -t 00:20:00
#SBATCH -N 1
#SBATCH -n 28
#For k80 cards use the line below
#SBATCH --gres=gpu:k80:2
#For V100 cards use the line below
###SBATCH --gres=gpu:v100:2
###SBATCH --exclusive
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --mail-type=END

ml purge  > /dev/null 2>&1 
ml GCC/8.3.0  CUDA/10.1.243  OpenMPI/3.1.4 
ml Amber/18.17-AmberTools-19.12-Python-2.7.16 

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
mpirun -np $num_dev pmemd.cuda.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${pstep}.rst7 -o ${istep}2.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7 -x ${istep}.nc 
mpirun -np 6 pmemd.cuda.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${pstep}.rst7 -o ${istep}3.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7 -x ${istep}.nc 

