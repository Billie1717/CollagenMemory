Workflow for creating gaussian spread images which can be directly analysed via the tool:

AFT alignment

# PointSpreadSingle.py

You could run this individually by e.g.:

mkdir /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/run_stretchNoRemodel_Xstretch1/run_dens0.025_teq3e6_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1//GaussianSpread/

python PointSpreadSingle.py /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/run_stretchNoRemodel_Xstretch1/run_dens0.025_teq3e6_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1/ dump.run_T1.0004185000.lammpstrj LASTFRAMEDEN0.025

python PointSpreadSingle.py ${trajdir} ${filename} ${dataname}

# MetaScript.py:
Or, you can create sbatch scripts for a selection of timepoints within that file using this script

mkdir /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/run_stretchNoRemodel_Xstretch1/run_dens0.025_teq3e6_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1//GaussianSpread/

python MetaScript.py /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/Alignment/run_stretchNoRemodel_Xstretch1/run_dens0.025_teq3e6_Frames100_tstep0.005_Nev1000_MP1.0_BP0.1_MD0.65_BD1.8_seed1/

python MetaScript.py: ${trajdir}

# GaussianSpreadImages.sh

This meta file can be run for multiple simulations with this bash script

TO NOTE/ EDIT:



Times sometimes need altering in MetaScript.py ::
IE dens 0.1 0.15 start form 620000 and go until 3720000
IE dens 0.025 0.05 start form 1085000  and go until 4185000 


# ImageCreation.py
mkdir /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/GaussianSpread/Example_data/Plots
mkdir /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/GaussianSpread/Example_data/Plots/vvmax_const12_sigma2.5/
mkdir /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/GaussianSpread/Example_data/Plots/vvmax_const12_sigma1.5/
mkdir /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/GaussianSpread/Example_data/Plots/vvmax_const12_sigma2.0/

$ python ImageCreation.py {GaussianDir}
(first time):
$ python3 -m venv path/to/venv
$ pip install pandas matplotlib seaborn

(other times):
$ source path/to/venv/bin/activate

$ python ImageCreation.py /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/GaussianSpread/Example_data/



