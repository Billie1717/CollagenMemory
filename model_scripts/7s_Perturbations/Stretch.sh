#!/bin/bash
#module load python3/recommended
Time_Str=1e5
Time_Rel=1e7
Time_Eq=3e6
Time_InitRel=2e5
NumFrames=100
BP=0.01
BD=1.8
MP=1.0
MD=0.65
Nevery=1000
tstep=0.01
XStretch=0.85
restartname='run_T1.2635000.restart'
dens=0.05
for seed in 1 2 3
do
for BD in 1.8 #1.8
do
for BP in 0.005 #0.05 #0.01
do
for frac in 0.25 0.44 #0.56 #0.75
do
input='run_makebonds/run_dens'${dens}'_frac'${frac}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}'/restart/'${restartname}
foldernameadd='runs_Stretch/run_frac'${frac}'_tStr'${Time_Str}'_tRel'${Time_Rel}'_tInitRel'${Time_InitRel}'_Xtretch'${XStretch}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}
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