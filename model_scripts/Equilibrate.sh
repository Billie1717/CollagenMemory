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
#input='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Initial_DataFiles/dumbbells_lattice_slab_n48020_bl3.00_Lx164.00_rho1.500e-01.lammpsdata'
#input='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Initial_DataFiles/dumbbells_lattice_slab_n32000_bl3.00_Lx164.00_rho1.000e-01.lammpsdata'
input='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Initial_DataFiles/dumbbells_lattice_slab_n15680_bl3.00_Lx164.00_rho5.000e-02.lammpsdata'
#input='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Initial_DataFiles/dumbbells_lattice_slab_n8000_bl3.00_Lx164.00_rho2.500e-02.lammpsdata'
dens=0.05
for seed in 6 #7 #4 5 #1 2 3
do
for BD in 1.8 #1.8
do
for BP in 0.1 #0.05 #0.01
do
#input='runs_Equilib/run_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}'/restart/run_T1.3100000.restart'
input='run_makebonds_densvar/run_dens'${dens}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}'/restart/run_T1.775000.restart' #/run_T1.1085000.restart'
Time_Eq=6e5
NumFrames=20
foldernameadd='run_makebonds_densvar/runRS_dens'${dens}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}
mkdir ${foldernameadd}
cp ${input} ${foldernameadd}'/data'
echo ${foldernameadd}
#python3 build_EquilibrateRestart.py ${foldernameadd} ${Time_Eq} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MP} ${BP} ${seed}
python3 build_Equilibrate.py ${foldernameadd} ${Time_Eq} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MP} ${BP} ${seed}
cd ${foldernameadd}
runscriptfile='runscript.sh'
sbatch $runscriptfile
cd .. 
cd ..
done
done
done