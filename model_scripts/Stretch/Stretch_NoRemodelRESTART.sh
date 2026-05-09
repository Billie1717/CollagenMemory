#!/bin/bash
#module load python3/recommended
Time_Str=2e5
Time_Rel=2e7
Time_InitRel=0 #1e5
Time_Eq=3e6
NumFrames=100
BP=0.01
BD=1.8
MP=1.0
MD=0.65
Nevery=1000
tstep=0.005
restartname='run_T1.21008000.restart' #'run_T1.15352000.restart' #'run_T1.620000.restart' 
restartname='run_T1.35148000.restart'
#restartname='run_T1.40804000.restart'
#'run_T1.1085000.restart' #(for 0.025 0.05, and run_T1.620000.restart for 0.1 0.15)
XStretch=1.0
BD=1.8
for dens in 0.05 #0.025 #0.15 0.1 #0.15 #0.05 #0.1 0.15
do
for seed in 3 # 1 2 #2 3 4 5
do
for XStretch in 1.0 #1.0 #0.85 #1.0 1.5 #1.8
do
for BP in 0.1 #0.05 #0.01
do
input='run_NoRemodel/runRESTART_dens'${dens}'_Xstretch'${XStretch}'_trel'${Time_Rel}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}'/restart/'${restartname}
foldernameadd='run_NoRemodel/runRESTARTRESTART_dens'${dens}'_Xstretch'${XStretch}'_trel'${Time_Rel}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}
mkdir ${foldernameadd}
cp ${input} ${foldernameadd}'/data'
echo ${foldernameadd}
python3 build_StretchNoRemodelRESTART.py ${foldernameadd} ${Time_Str} ${Time_Rel} ${Time_InitRel} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MP} ${BP} ${XStretch} ${seed}
cd ${foldernameadd}
runscriptfile='runscript.sh'
sbatch $runscriptfile
cd .. 
cd ..
done
done
done
done