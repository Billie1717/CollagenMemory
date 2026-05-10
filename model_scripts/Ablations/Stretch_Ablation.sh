#!/bin/bash
#module load python3/recommended
Time_Str=1e5
Time_Rel=3e6
Time_InitRel=1e5
NumFrames=100
BP=0.01
BD=1.8
MP=1.0
MD=0.65
Nevery=1000
tstep=0.005
XStretch=1 #In %
DataID=4 # Ablations 1 small : 2 bigger : 4 massive 
for XStretch in 0.0 # 0.85
do
for seed in 1 2
do
for BD in 1.8
do
for BP in 0.05 #0.01
do
input='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Ablations/datafiles/run_data'${DataID}'_T1.465000.restart'
foldernameadd='run_FirstAblations/runNoRemodAblate_tStr'${Time_Str}'_tRel'${Time_Rel}'_tInitRel'${Time_InitRel}'_DataID'${DataID}'_Xstretch'${XStretch}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}
mkdir ${foldernameadd}
cp ${input} ${foldernameadd}'/data'
echo ${foldernameadd}
#echo ${foldernameadd} ${Time_Str} ${Time_Rel} ${Time_Eq} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MP} ${BP} ${XStretch} ${seed}
python3 build_StretchAblate.py ${foldernameadd} ${Time_Str} ${Time_Rel} ${Time_InitRel} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MP} ${BP} ${XStretch} ${seed}
cd ${foldernameadd}
runscriptfile='runscript.sh'
sbatch $runscriptfile
cd .. 
cd ..
done
done
done
done