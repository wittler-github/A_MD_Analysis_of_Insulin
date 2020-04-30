#!/bin/bash
#PBS -W group_list=ma4
#PBS -q longq
#PBS -N 3p2KJJpTc36m5		
#PBS -r n
#PBS -j oe
#PBS -o 3p2KJJpTc36m5.log

#PBS -l select=1:ncpus=16:mem=60G:mpiprocs=16,walltime=500:00:00

date -u '+%Y-%m-%d %TZ %a'
hostname
echo 'Run GROMACS'

cd /home/hpw572/MDRUNS/2KJJEXP/2KJJpTc36mno5/m5MD

source /usr/share/modules/init/bash
module load gromacs/5.0.4_intel_15.6

gmx_mpi grompp -f 3pmd.mdp -c 2pmd.gro -t 2pmd.cpt -p topol.top -po 3pmdout.mdp -o 3pmd.tpr
gmx_mpi mdrun -deffnm 3pmd 

date -u '+%Y-%m-%d %TZ %a'

