#!/bin/bash
#SBATCH -A staff
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 7
#SBATCH --gres=gpu:k80:2
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --mail-type=END
#SBATCH --exclusive



ml purge  > /dev/null 2>&1 
ml GCC/8.3.0  CUDA/10.1.243  OpenMPI/3.1.4
ml GROMACS/2021-Python-3.7.4

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_charmm2gmx.pdb -n index.ndx -p topol.top
gmx mdrun -nb gpu -pme gpu -npme 1 -ntomp 7 -ntmpi $SLURM_NTASKS -resetstep 10000 -nsteps 20000 -dlb yes  -v -deffnm step4.1_equilibration
srun gmx_mpi mdrun -nb gpu -pme gpu -npme 1 -ntomp 7  -resetstep 10000 -nsteps 10000 -dlb yes  -v -deffnm step4.1_equilibration




exit 0
