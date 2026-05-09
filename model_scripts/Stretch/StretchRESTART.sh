#!/bin/bash
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
tstep=0.005
XStretch=0 
restartname='run_T1.6060000.restart'
restartname='run_T1.9090000.restart'
restartname='run_T1.13130000.restart'
restartname='run_T1.12120000.restart'
restartname='run_T1.10100000.restart'
restartname='run_T1.23230000.restart'
restartname='run_T1.25250000.restart'
restartname='run_T1.22220000.restart'
restartname='run_T1.18180000.restart'
restartname='run_T1.20200000.restart'
restartname='run_T1.15150000.restart'
restartname='run_T1.14140000.restart'
restartname='run_T1.20200000.restart'
restartname='run_T1.28280000.restart'
restartname='run_T1.35350000.restart'
#restartname='run_T1.37370000.restart'
restartname='run_T1.41410000.restart'
dens=0.05
for XStretch in 1.0 #0.85 1.0 #0.85 #1.0 #0.5 1
do
for seed in 4 #2 3
do
for BD in 1.8
do
for BP in 0.01 #0.01
do
#input='run_Remodel/run_tStr'${Time_Str}'_tRel'${Time_Rel}'_tInitRel'${Time_InitRel}'_Xtretch'${XStretch}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}'/restart/'${restartname}
input='run_Remodel/runRESTARTRESTART_tStr'${Time_Str}'_tRel'${Time_Rel}'_tInitRel'${Time_InitRel}'_Xtretch'${XStretch}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}'/restart/'${restartname}
#input='run_Remodel/runRESTART_tStr'${Time_Str}'_tRel'${Time_Rel}'_tInitRel'${Time_InitRel}'_Xtretch'${XStretch}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}'/restart/'${restartname}
foldernameadd='run_Remodel/runRESTARTRESTARTRESTART_tStr'${Time_Str}'_tRel'${Time_Rel}'_tInitRel'${Time_InitRel}'_Xtretch'${XStretch}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}
mkdir ${foldernameadd}
cp ${input} ${foldernameadd}'/data'
echo ${foldernameadd}
python3 build_StretchRESTART.py ${foldernameadd} ${Time_Str} ${Time_Rel} ${Time_InitRel} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MP} ${BP} ${XStretch} ${seed}
cd ${foldernameadd}
runscriptfile='runscript.sh'
sbatch $runscriptfile
cd .. 
cd ..
done
done
done
done