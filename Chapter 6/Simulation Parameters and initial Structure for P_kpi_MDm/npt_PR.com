#!/bin/bash
#PBS -N nptpr2KJJpTc36m5
#PBS -r n
#PBS -j oe
#PBS -o nptpr2KJJpTc36m5.log
#PBS -l walltime=20:00:00,vmem=2GB
#PBS -l nodes=1:ppn=16

date -u '+%Y-%m-%d %TZ %a'
hostname
echo 'Run GROMACS'

cd /lustre/platr0008/wittler/MDRUNS/2KJJpTc36mno5/m5MD

module load gromacs/5.0.4-openmpi-gcc

gmx grompp -f npt_PR.mdp -c npt_B.gro -t npt_B.cpt -p topol.top -po npt_PRout.mdp -o npt_PR.tpr
gmx mdrun -deffnm npt_PR -ntomp 16 -ntmpi 1 

date -u '+%Y-%m-%d %TZ %a'

