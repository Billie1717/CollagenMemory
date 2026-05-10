#!/bin/bash
#module load python3/recommended
#foldername='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/run_stretchNoRemodel_Xstretch1.5/run_dens'${dens}'_teq3e6_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1'
#foldername='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/run_stretchNoRemodel_Xstretch'${XStretch}'/run_dens'${dens}'_Xstretch'${XStretch}'_teq3e6_Frames100_tstep'${tstep}'_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed'${seed}
#foldername='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/run_stretchNoRemodel_FROMData2/run_dens'${dens}'_Xstretch'${XStretch}'_teq3e6_Frames100_tstep0.01_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed'${seed}
Time_Str=1e5
Time_Rel=3e6
Time_InitRel=0 #1e5
Time_Eq=3e6
NumFrames=100
BP=0.01
BD=1.8
MP=1.0
MD=0.65
Nevery=1000
tstep=0.01
XStretch=1.0
BD=1.8
BP=0.1
dens=0.05
for BP in 0.005 #0.15 0.1 #0.025 0.05 #0.15 #0.05 #0.1 0.15
do
for seed in 1 #2 3 #1 2 3 #4 5 6 #3 #2 3
do
for FracTagged in 0.44 #0.25 #0.44 #0.25 #0.56 0.75
do
for XStretch in 1.0 #0.85 #1.5 #0.05 #0.01
do
#foldername='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/7sPerturbation/runs_Stretch/run_frac'${FracTagged}'_tStr1e5_tRel1e7_tInitRel2e5_Xtretch'${XStretch}'_Frames100_tstep0.01_Nev1000_MP1.0_BP'${BP}'_MD0.65_BD1.8_seed'${seed}
foldername='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/7sPerturbation/runs_Stretch/runRESTART_frac0.44_tStr1e5_tRel1e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.005_MD0.65_BD1.8_seed1'
#foldername='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/7sPerturbation/run_makebonds/run_dens0.05_frac'${FracTagged}'_teq3e6_Frames100_tstep0.01_Nev1000_MP1.0_BP0.005_MD0.65_BD1.8_seed'${seed}
foldernameadd=${foldername}'/GaussianSpread/'
mkdir ${foldernameadd}
echo "added folder" ${foldernameadd}
python MetaScript_7s.py ${foldername} ${FracTagged}
done
done
done
done