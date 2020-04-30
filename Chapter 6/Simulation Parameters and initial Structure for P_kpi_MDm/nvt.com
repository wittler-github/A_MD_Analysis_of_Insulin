#!/bin/bash
#PBS -N nvt2KJJpTc36m5
#PBS -r n
#PBS -j oe
#PBS -o nvt2KJJpTc36m5.log
#PBS -l walltime=15:00:00,vmem=2GB
#PBS -l nodes=1:ppn=16

date -u '+%Y-%m-%d %TZ %a'
hostname
echo 'Run GROMACS'

cd /lustre/platr0008/wittler/MDRUNS/2KJJpTc36mno5/m5MD

module load gromacs/5.0.4-openmpi-gcc

gmx grompp -f nvt.mdp -c em.gro -p topol.top -po nvtout.mdp -o nvt.tpr
mdrun -deffnm nvt -ntomp 16 -ntmpi 1 

date -u '+%Y-%m-%d %TZ %a'

