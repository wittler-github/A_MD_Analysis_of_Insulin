#!/bin/bash
#PBS -N nptb2KJJpTc36m5
#PBS -r n
#PBS -j oe
#PBS -o nptb2KJJpTc36m5.log
#PBS -l walltime=20:00:00,vmem=2GB
#PBS -l nodes=1:ppn=16

date -u '+%Y-%m-%d %TZ %a'
hostname
echo 'Run GROMACS'

cd /lustre/platr0008/wittler/MDRUNS/2KJJpTc36mno5/m5MD

module load gromacs/5.0.4-openmpi-gcc

gmx grompp -f npt_B.mdp -c nvt.gro -t nvt.cpt -p topol.top -po npt_Bout.mdp -o npt_B.tpr
gmx mdrun -deffnm npt_B -ntomp 16 -ntmpi 1 

date -u '+%Y-%m-%d %TZ %a'

