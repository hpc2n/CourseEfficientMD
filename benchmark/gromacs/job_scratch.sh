#!/bin/bash
#SBATCH -A SNIC2020-9-25
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 7
#SBATCH --output=job_scratch.out
#SBATCH --error=job_scratch.err
#SBATCH --mail-type=END
#SBATCH --exclusive

ml purge

#ml GCC/6.4.0-2.28  CUDA/9.0.176  impi/2017.3.196
#ml GROMACS/2018
#.1
#ml GCC/7.3.0-2.30  CUDA/9.2.88  OpenMPI/3.1.1
#ml GROMACS/2018.3 


#ml GCC/5.4.0-2.26  CUDA/8.0.61_375.26  impi/2017.3.196
#ml GROMACS/2016.3-PLUMED 

ml GCC/8.2.0-2.31.1  CUDA/10.1.105  OpenMPI/3.1.3
ml GROMACS/2019.2

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    mdargs="-ntomp $SLURM_CPUS_PER_TASK"
else
    mdargs="-ntomp 1"
fi

export parent=/pfs/nobackup/home/p/pojedama/benchmarks/charmm-gui/gromacs
mkdir /scratch/test
rsync -avzh $parent/scratch/ /scratch/test/

cd /scratch/test 

gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_charmm2gmx.pdb -n index.ndx -p topol.top
gmx mdrun -ntmpi 4 $mdargs -deffnm step4.1_equilibration

cd $parent
rsync -avzh /scratch/test/ $parent/scratch/

rm -rf /scratch/test



exit 0
