#!/bin/bash
# This is the job script to run a GROMACS example job 
# Remember to change the following to your own Project ID! 
#SBATCH -A Project_ID
#SBATCH -J Gromacs
#SBATCH -t 00:22:00
#SBATCH -n *FIXME*

# It is a good idea to do a ml purge before loading other modules
ml purge > /dev/null 2>&1

ml GCC/13.2.0  OpenMPI/4.1.6
ml GROMACS/2024.1

export MDRUN='gmx_mpi mdrun'
gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_charmm2gmx.pdb -n index.ndx -p topol.top
mpirun -np $SLURM_NTASKS $MDRUN -dlb yes  -v -deffnm step4.1_equilibration

