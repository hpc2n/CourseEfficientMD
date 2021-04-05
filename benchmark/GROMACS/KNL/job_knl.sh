#!/bin/bash
#SBATCH -A SNICyyyy-xx-yy
#SBATCH -t 00:15:00
#SBATCH -N 1
#SBATCH -n 68
#SBATCH -p knl
#SBATCH --constraint=cache,quad
#SBATCH --output=job_str.out
#SBATCH --error=job_str.err
#SBATCH --mail-type=END

#ml icc/2017.4.196-GCC-6.4.0-2.28 ifort/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196
#ml GROMACS/2016.4

#ml GCC/7.3.0-2.30 OpenMPI/3.1.1
#ml GROMACS/2018.3

#ml icc/2018.3.222-GCC-7.3.0-2.30 ifort/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
#ml GCC/7.3.0-2.30  OpenMPI/3.1.1
#ml GROMACS/2018.3

#ml GCC/7.3.0-2.30  OpenMPI/3.1.1
#ml GROMACS/2019

ml icc/2018.3.222-GCC-7.3.0-2.30 ifort/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
ml GROMACS/2019


# Automatic selection of options to mdrun_mpi depending on parameters given to SBATCH
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    md="-ntomp $SLURM_CPUS_PER_TASK"
else
    md="-ntomp 1"
fi

gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_charmm2gmx.pdb -n index.ndx -p topol.top


numactl -H

export OMP_NUM_THREADS=2
mpirun -np 68 gmx_mpi mdrun -ntomp 2 -npme 18 -pin on -pinoffset 0 -pinstride 2 -dlb auto -v -deffnm step4.1_equilibration

exit 0
