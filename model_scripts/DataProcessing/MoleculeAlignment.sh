#!/bin/bash
#module load python3/recommended
NumMols=15680 #dens0.05 OG 15680
name='run_makebonds_densvar/run_dens0.05_teq3e6_Frames100_tstep0.01_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1/'
name2='run_makebonds_densvarrun_dens0.05_teq3e6_Frames100_tstep0.01_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1'
#name='run_stretchNoRemodel_Xstretch1/run_dens0.05_teq3e6_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1/'
#name2='run_stretchNoRemodel_Xstretch1run_dens0.05_teq3e6_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed'
name='runs_Stretch/run_tStr1e5_tRel1e6_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.02_MD0.65_BD1.8_seed1/'
name2='runs_Stretchrun_tStr1e5_tRel1e6_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.02_MD0.65_BD1.8_seed1'
name='runs_Stretch/run_tStr1e5_tRel1e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.001_MD0.65_BD1.8_seed2'
name2='runs_Stretchrun_tStr1e5_tRel1e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.001_MD0.65_BD1.8_seed2'
name='runs_Stretch/run_tStr1e5_tRel1e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.005_MD0.65_BD1.8_seed1'
name2='runs_Stretchrun_tStr1e5_tRel1e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.005_MD0.65_BD1.8_seed1'
name='runs_Stretch/run_tStr1e5_tRel1e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.001_MD0.65_BD1.8_seed4'
name2='runs_Stretchrun_tStr1e5_tRel1e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.001_MD0.65_BD1.8_seed4'
name='runs_Stretch/run_tStr1e5_tRel1e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.01_MD0.65_BD1.8_seed2'
name2='runs_Stretchrun_tStr1e5_tRel1e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.01_MD0.65_BD1.8_seed2'
name='run_Remodel/run_tStr2e5_tRel2e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.005_Nev1000_MP1.0_BP0.01_MD0.65_BD1.8_seed4'
name2='run_Remodelrun_tStr2e5_tRel2e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.005_Nev1000_MP1.0_BP0.01_MD0.65_BD1.8_seed4'
name='run_NoRemodel/run_dens0.05_Xstretch0.0_trel2e7_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed3'
name2='run_Remodelrun_dens0.05_Xstretch0.0_trel2e7_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed3'
#name='run_NoRemodel/runRESTART_dens0.05_Xstretch1.0_trel2e7_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1'
#name2='run_NoRemodelrunRESTART_dens0.05_Xstretch1.0_trel2e7_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1'
#name='run_Remodel/run_tStr2e5_tRel2e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.005_Nev1000_MP1.0_BP0.01_MD0.65_BD1.8_seed5'
#name2='run_Remodelrun_tStr2e5_tRel2e7_tInitRel2e5_Xtretch1.0_Frames100_tstep0.005_Nev1000_MP1.0_BP0.01_MD0.65_BD1.8_seed5'
# name='run_NoRemodel/run_dens0.05_Xstretch1.0_trel2e7_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed5'
# name2='run_NoRemodelrun_dens0.05_Xstretch1.0_trel2e7_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed5'
# name='runs_Stretch/run_tStr1e5_tRel1e7_tInitRel2e5_Xtretch0.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.005_MD0.65_BD1.8_seed1'
# name2='runs_Stretchrun_tStr1e5_tRel1e7_tInitRel2e5_Xtretch0.0_Frames100_tstep0.01_Nev1000_MP1.0_BP0.005_MD0.65_BD1.8_seed1'
# name='run_makebonds_densvar/run_dens0.05_teq3e6_Frames100_tstep0.01_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1'
# name2='run_makebonds_densvarrun_dens0.05_teq3e6_Frames100_tstep0.01_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1'
for seed in 1
do
filepattern0='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/'${name}'/dumplin/'
fileOut='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/Data/MolAlign'${name2}'.txt'
fileOut2='/nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/Data/MolAlignPBC'${name2}'.txt'
echo ${filepattern1}
for filename in ${filepattern0}*; do
echo ${filename}
python /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/model2_dumbell/scripts/Analysis/MoleculeAngle.py ${filename} ${fileOut} ${NumMols}
#python /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/model2_dumbell/scripts/Analysis/MoleculeAngle_PBC.py ${filename} ${fileOut2} ${NumMols}
done
done