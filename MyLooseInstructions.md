The lammps_src code contains custom lammps reaction package source codes. The lammps version used for this project was lammps-15Jun2023. Lammps must be built with these replaces reaction packages in the src/REACTION folder.

Each of the folders in model_scripts has the lammps inputs and instructions for running the simulation experiments: Stretch (Fig3), 7s_Perturbations (Fig4), and Ablations (FigS3). These folders contain annotated and simplified scripts for understanding how the simulations are working. The GaussianSpread data processing was used for the analysis in Figs 3&4.  Other data processing (like finding the stress or molecule alignment or counting the molecule binding types) are detailed in DataProcessing folder. 

Figures has plotting scripts and data for the various figures in the paper.

To create initial config files follow instructions in model_scripts/Initial_config/Readme.md

Documenting where the final version of the simulations are in ISTA cluster (For myself):


Stretch:
1.) To equilibrate network bach Equilibrate.sh and change path to initial configuration in Initial_config/dumbbells*
2.) Run equilibrate script with mpirun -np 1 /nfs/scistore26/saricgrp/bmeadowc/Scratch/lammps-15Jun2023/src/lmp_mpi -in collagen.in
3.) The networks in the paper were equilibrated for 3e6 timesteps
4.) Run stretch and relax with and without remodelling with Stretch.sh which calls build_Stretch.py and Stretch_NoRemodel.sh which calls build_StretchNoRemodel.py. The RESTART files are for restarting the relaxation after stretch, which in the paper was utilised multiple times to have long enough simulations to see networks relax alignment via remodelling.
--> these simulations are now ready for data processing using GaussianSpread.

To understand the stretch-remodel code, the script build_Stretch_SimpleAnnotated.py is a commented and simplified version of the code for readability. It is not a working version. 

7s_Perturbations:
The 7s perturbation simulations have different templates and molecule types. We introduce a new bead type which represents 7s ends which cannot bind. The new type is type 9 and a different initial config must be used which creates a number of molecules where some fraction have 'perturbed' type 9 ends. Type 9s do not participate in any binding reactions. This is opposed to type 2 which are the 'normal' 7s binding ends which have not yet bound. The templates are in the templates folder. 

Ablations:
These simulations are run on very large networks. The equilibrated input script is in Ablations/SimulationsStartFromRestartFile/. The networks are stretched, relaxed for a short time then ablated. Some simulations are also run for pure stretched networks with no ablations for benchmarking stress (Figure S3a). The templates are the same for those use in Stretch simulations.

Understanding template structure:
The bonding/unbonding of the collagen IV molecules happen via fix react and the templates. One bonding pattern will be specified by a pre, post and map template file. The pre and post specify the structure of a bonding pattern before and after the reaction occurs and the map instructs how to translate the molecule identifiers/properties from the pre to post molecule. Due to the many combinations of types of molecules, there are 100s of template patterns. More information about templates can be found: https://docs.lammps.org/fix_bond_react.html. 


TO do (myself):
Put templates in relevant folder(s) (ssh from cluster)
Gaussian spread data analysis documentation
Sort data (cluster) > delete not used, move used to home/fs