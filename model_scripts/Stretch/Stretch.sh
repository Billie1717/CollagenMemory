#!/bin/bash
#module load python3/recommended
Time_Str=2e5
Time_Rel=2e7
Time_Eq=3e6
Time_InitRel=2e5
NumFrames=100
BP=0.01
BD=1.8
MP=1.0
MD=0.65
Nevery=1000
tstep=0.01
XStretch=0 
restartname='run_T1.1098000.restart' #'run_T1.1085000.restart'
dens=0.05
for XStretch in 1.0 #0.85 #1.0 #0.5 1
do
for seed in 6 #4 5 #2 3 #4 5
do
for BD in 1.8
do
for BP in 0.1 #0.01
do
tstep=0.01
#Scratch/Collagen/NargessProject/Alignment/run_makebonds_densvar/runRS_dens0.05_teq6e5_Frames20_tstep0.01_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed6/restart
#input='run_makebonds_densvar/runRS_dens'${dens}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}'/restart/'${restartname}
input='run_makebonds_densvar/runRS_dens0.05_teq6e5_Frames20_tstep0.01_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed6/restart/'${restartname}
BP=0.01
tstep=0.005
foldernameadd='run_Remodel/run_tStr'${Time_Str}'_tRel'${Time_Rel}'_tInitRel'${Time_InitRel}'_Xtretch'${XStretch}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}
mkdir ${foldernameadd}
cp ${input} ${foldernameadd}'/data'
echo ${foldernameadd}
python3 build_Stretch.py ${foldernameadd} ${Time_Str} ${Time_Rel} ${Time_InitRel} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MP} ${BP} ${XStretch} ${seed}
cd ${foldernameadd}
runscriptfile='runscript.sh'
sbatch $runscriptfile
cd .. 
cd ..
done
done
done
done