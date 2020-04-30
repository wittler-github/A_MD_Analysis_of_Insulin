#!/bin/bash
#PBS -N MD2KJJpTc36m5
#PBS -j oe
#PBS -o MD2KJJpTc36m5.log
#PBS -r n 
#PBS -l walltime=05:00:00

date -u '+%Y-%m-%d %TZ %a'
hostname
echo 'Run GROMACS'

cd /lustre/platr0008/wittler/MDRUNS/2KJJpTc36mno5/m5MD

# Run all on trifid v. 5.0.4. trifid 
module load gromacs/5.0.4-openmpi-gcc

cp -rp ../startPDBmno5/2KJJpTc36m5.pdb .

# 1) Make Gromacs topology
gmx pdb2gmx -f 2KJJpTc36m5.pdb -o gmx.gro -ignh -merge all -water tip3p -ff charmm36-mar2014-NLG-NLE 

# 2) Define box
gmx editconf -f gmx.gro -o newbox.gro -c -d 1.0 -bt cubic -box 5.45 5.45 5.45   

# 3) Adding water solvent
gmx solvate -cp newbox.gro -o solv.gro -p topol.top -cs spc216.gro                                  
       
# 4) Replace water with NA and CL
gmx grompp -f em.mdp -c solv.gro -p topol.top -po em_ions_out.mdp -o ions.tpr	
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -pq 1 -np 10 -nname CL -nq -1 -nn 8

# 5) Energy Minimization           
m5EM=$(qsub em.com)
                         
# 6) Equilibration NVT 100ps
m5NVT=$(qsub -W depend=afterok:$m5EM nvt.com)	

# 7) Equilibration NPT (Berendsen) 400 [ps]
m5NPTB=$(qsub -W depend=afterok:$m5NVT npt_B.com)

# 8) Equilibration NPT (Parrinello-Rahman) 400 [ps]
m5NPTPR=$(qsub -W depend=afterok:$m5NPTB npt_PR.com)

# 9) 1'st Production MD Run 500ns
m5ONEPMD=$(qsub -W depend=afterok:$m5NPTPR 1pmd.com)

# 10) 2'nd Production MD Run 500ns
m5TWOPMD=$(qsub -W depend=afterok:$m5ONEPMD 2pmd.com)

# 11) 3'rd Production MD Run 500ns
m5THREEPMD=$(qsub -W depend=afterok:$m5TWOPMD 3pmd.com)


echo "$m5EM $m5NVT $m5NPTB $m5NPTPR $m5ONEPMD $m5TWOPMD $m5THREEPMD" > ../RESET/RESETm5.txt


date -u '+%Y-%m-%d %TZ %a'

