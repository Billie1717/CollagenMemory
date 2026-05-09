#!/bin/bash
Time_Eq=3e6
NumFrames=100
BP=0.01
BD=1.8
MP=1.0
MD=0.65
Nevery=1000
tstep=0.01
for dens in 0.025 0.05 0.1 0.15
do
for seed in 1 2 3
do
for BD in 1.8 #1.8
do
for BP in 0.1 #0.05 #0.01
do
filepattern1='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/run_makebonds_densvar/run_dens'${dens}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}'/'
fileOut1='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/Data/HistLarge_dens'${dens}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_MP'${MP}'_BP'${BP}'_MD'${MD}'_BD'${BD}'_seed'${seed}'.txt'
echo ${filepattern1}
python /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/model2_dumbell/scripts/Analysis/Histogram_Types.py ${filepattern1} ${fileOut1}
done
done
done
done