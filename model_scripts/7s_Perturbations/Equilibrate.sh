#!/bin/bash
#module load python3/recommended
Time_Eq=3e6
NumFrames=100
BP=0.01
BD=1.8
MP=1.0
MD=0.65
Nevery=1000
tstep=0.01
dens=0.05
for seed in 1 2 3
do
for BD in 1.8 #1.8
do
for BP in 0.005 #0.05 #0.01
do
for frac in 0.25 0.44 #0.56 #0.75
do
input='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Initial_DataFiles/dumbbells_7s_frac'${frac}'_n15680_bl3.00_Lx164.00_rho5.000e-02.lammpsdata'
foldernameadd='run_makebonds/run_dens'${dens}'_frac'${frac}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}
mkdir ${foldernameadd}
cp ${input} ${foldernameadd}'/data'
echo ${foldernameadd}
python3 build_Equilibrate.py ${foldernameadd} ${Time_Eq} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MP} ${BP} ${seed}
cd ${foldernameadd}
runscriptfile='runscript.sh'
sbatch $runscriptfile
cd .. 
cd ..
done
done
done
done