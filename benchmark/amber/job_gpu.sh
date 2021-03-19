#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -t 00:20:00
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --gres=gpu:k80:2
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --mail-type=END

ml purge  > /dev/null 2>&1 
ml GCC/7.3.0-2.30  CUDA/9.2.88  OpenMPI/3.1.1
ml Amber/18-AmberTools-18-patchlevel-10-8 


nvidia-smi

export num_dev=`echo $CUDA_VISIBLE_DEVICES | awk 'BEGIN{FS=","};{print NF}'`

#specify the cuda device
#export CUDA_VISIBLE_DEVICES=1

####srun sander.MPI -ng 8 -i step4.0_minimization.mdin
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

#Combining instances:
#mpirun -np 4 pmemd.cuda.MPI -O -i ${istep}.mdin -p ${init}.parm7 -c ${pstep}.rst7 -o ${istep}.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -ref ${init}.rst7 -x ${istep}.nc & 
#
#export istep2="step4.1_equilibration_parall"
#mpirun -np 24 pmemd.MPI -O -i ${istep2}.mdin -p ${init}.parm7 -c ${pstep}.rst7 -o ${istep2}.mdout -r ${istep2}.rst7 -inf ${istep2}.mdinfo -ref ${init}.rst7 -x ${istep2}.nc & 
#
#wait


exit 0
