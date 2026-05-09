The lammps_src code contains custom lammps reaction package source codes. The lammps version used for this project was lammps-15Jun2023. Lammps must be built with these replaces reaction packages in the src/REACTION folder.

Each of the folders in model_scripts has the lammps inputs and instructions for running the simulation experiments: Stretch (Fig3), 7s_Perturbations (Fig4), and Ablations (FigS3). The GaussianSpread data processing was used for the analysis in Figs 3&4. The folder model_example has a basic annotated script of the collagen model. Other data processing (like finding the stress or molecule alignment) are detailed in DataProcessing folder. 

Figures has plotting scripts and data for the various figures in the paper.

To create initial config files follow instructions in model_scripts/Initial_config/Readme.md

Documenting where the final version of the simulations are in ISTA cluster (For myself):


Stretch:
1.) To equilibrate network bach Equilibrate.sh and change path to initial configuration in Initial_config/dumbbells*
2.) Run equilibrate script with mpirun -np 1 /nfs/scistore26/saricgrp/bmeadowc/Scratch/lammps-15Jun2023/src/lmp_mpi -in collagen.in
3.) The networks in the paper were equilibrated for 3e6 timesteps
4.) Run stretch and relax with and without remodelling with Stretch.sh which calls build_Stretch.py and Stretch_NoRemodel.sh which calls build_StretchNoRemodel.py. The RESTART files are for restarting the relaxation after stretch, which in the paper was utilised multiple times to have long enough simulations to see networks relax alignment via remodelling.
--> these simulations are now ready for data processing using GaussianSpread.