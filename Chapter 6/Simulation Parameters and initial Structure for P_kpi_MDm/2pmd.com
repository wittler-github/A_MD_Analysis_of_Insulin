#!/bin/bash
#PBS -W group_list=ma4
#PBS -q longq
#PBS -N 2p2KJJpTc36m5		
#PBS -r n
#PBS -j oe
#PBS -o 2p2KJJpTc36m5.log

#PBS -l select=1:ncpus=16:mem=60G:mpiprocs=16,walltime=500:00:00

date -u '+%Y-%m-%d %TZ %a'
hostname
echo 'Run GROMACS'

cd /home/hpw572/MDRUNS/2KJJEXP/2KJJpTc36mno5/m5MD

source /usr/share/modules/init/bash
module load gromacs/5.0.4_intel_15.6

gmx_mpi grompp -f 2pmd.mdp -c 1pmd.gro -t 1pmd.cpt -p topol.top -po 2pmdout.mdp -o 2pmd.tpr
gmx_mpi mdrun -deffnm 2pmd 

date -u '+%Y-%m-%d %TZ %a'

