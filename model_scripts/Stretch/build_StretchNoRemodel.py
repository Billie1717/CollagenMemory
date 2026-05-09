from __future__ import division, print_function
import numpy as np
#import matplotlib.pyplot as plt
import math
import time
import random
import sys, getopt
import os


def main():

    ################ Parameters ################
    workingdir = str(sys.argv[1])
    #${foldernameadd} ${Time_Eq} ${NumFrames} ${Nevery} ${MD} ${BD} ${MP} ${BP} ${seed}
    tStr = str(sys.argv[2])
    tRel = str(sys.argv[3])
    tEq = str(sys.argv[4])
    Frames =  str(sys.argv[5])
    Nevery = str(sys.argv[6])
    tstep = str(sys.argv[7])
    MakeDist = str(sys.argv[8])
    BreakDist = str(sys.argv[9])
    MakeProb = str(sys.argv[10])
    BreakProb = str(sys.argv[11])
    Xstretch = str(sys.argv[12])
    seed  = str(sys.argv[13])
    
    #cutoff = np.sqrt((ChemBond-np.log(0.01))/6)+Blen #Cutoff the making probabbilty at Prob = 0.01, with bond constant 6 
    #print("cutoff"+str(cutoff))
    ################ write in.local ################
    print('Writing collagen.in\n')
    outfile2 = workingdir + '/collagen.in'
    fo=open(outfile2,'w')
    fo.write(
        '''#####################################################
#PARAMETERS

#NUMBER OF STEPS OF RUNS (stages)
''')
    #fo.write('variable tmix			equal 	'+str(tmix)+'	#NVT steps for molecule mixing\n')
    fo.write('variable tstretch			equal 	'+str(tStr)+'	#NVT steps for bond equilibration\n')
    fo.write('variable trelax			equal 	'+str(tRel)+'	#NVT steps for bond equilibration\n')
    fo.write('variable tinitrel			equal 	'+str(tEq)+'	#NVT steps for bond equilibration\n')
    fo.write('variable Frames			equal 	'+str(Frames)+'	#NVT steps for bond equilibration\n')
    fo.write(
        '''

#THERMODYNAMIC PARAMETERS
variable mytemp         equal 	1.0	    	#Temperature

#INFORMATION PRINTING PERIODS
variable ttotal     equal   $(v_trelax+v_tstretch)		#total timesteps
variable thermodump     equal   $(floor(v_ttotal/5e4))		#Thermo info dump period
variable trajdump       equal   $(floor(v_ttotal/v_Frames))    	#Config info dump period
variable trestart       equal   $(floor(v_ttotal/50))    	#Restart configuration saving period\n''')

    fo.write(
        '''
#OTHER SIMULATION PARAMETERS
''')
    fo.write('variable tstep        	equal   '+str(tstep)+'   	#Integration time step for first part of simulation')
    fo.write(
        '''
variable damp   		equal   0.1	   	#Langevin thermostat damping coefficient [F_friction = -(m/damp)*v]
variable neigh_cutoff	equal	0.4  		#Neighbor list cutoff
variable lj_cutoff      equal   1.122462	#For WCA: 1.122462048 

#new variables Fernanda
variable RminF  equal  0.
variable RmaxFB  equal  30

variable NeveryFast equal 5

''')
    fo.write('variable RmaxF  equal '+str(MakeDist)+' #Distance to make a bond\n')
    fo.write('variable probF  equal '+str(MakeProb)+' #Distance to make a bond\n')
    fo.write('variable RminFB  equal '+str(BreakDist)+' #Distance to make a bond\n')
    fo.write('variable probFB  equal '+str(BreakProb)+' #Distance to make a bond\n')
    fo.write(
        '''

variable lseed_base equal 13579
variable lseed1 equal ${lseed_base}+1 #NC1 creation
variable lseed2 equal ${lseed_base}+2 #dimer creation
variable lseed3 equal ${lseed_base}+2 #line-trimer creation
variable lseed4 equal ${lseed_base}+2 #triangle-trimer creation
variable lseed5 equal ${lseed_base}+2 #open-tetramer creation
variable lseed6 equal ${lseed_base}+2 #rhomboid-tetramer creation
variable lseed7 equal ${lseed_base}+2 #squared-tetramer creation
variable lseed8 equal ${lseed_base}+2 #squared-tetramer break
variable lseed9 equal ${lseed_base}+2 #rhomboid-tetramer break
variable lseed10 equal ${lseed_base}+2 #open-tetramer break
variable lseed11 equal ${lseed_base}+2 #triangle-trimer break
variable lseed12 equal ${lseed_base}+2 #line-trimer break
variable lseed13 equal ${lseed_base}+2 #dimer break
variable lseed14 equal ${lseed_base}+2 #NC1 break
''')
    fo.write('variable NeverySlow equal '+str(Nevery)+' #how often we attempt reactions\n')
    fo.write('variable k_nc1			equal	6.0 		#Elastic constant A type bond\n')
    fo.write('variable r0_nc1			equal	0.5 		#Rest length A type bond\n\n')

    fo.write('variable k_7s			equal	6.0 		#Elastic constant B type bond\n')
    fo.write('variable r0_7s			equal	0.5 		#Rest length B type bond\n\n')

    fo.write('variable max_bonds_7s	equal   3				#Max. number of bonds of type A that each particle can form\n')
    fo.write(
        '''variable zero_cutoff    equal   6.0				#Cutoff for zero potential; must be larger that max bonding distance 

#SEEDS
variable lseed			file 	lseed.dat	#Seed for fix langevin
variable vseed			file 	vseed.dat	#Seed for creation of velocities

#STRINGS
variable fname 			string 	data		#name of the data file to be read to start the simulation
variable thermofile     string  thermo.dat	#name of the file with thermodynamic infos
variable simname 		string 	run_T${mytemp}	#name of this simulation

#####################################################
#INITALIZATION
units           lj          	#use lj units
boundary        p p f       	#periodic boundary conditions
atom_style     	molecular

#BOND STYLE: 
bond_style      harmonic

#ANGLE STYLE: 
angle_style      harmonic

#TEMPORARY PAIR STYLE TO START SIMULATION
pair_style		lj/cut 1.0

#Starting from an initial configuration WITHOUT BONDS (format lammpsdata):
#read_data       ${fname} extra/bond/types 2 extra/bond/per/atom 10 extra/special/per/atom 100 extra/angle/types 3 extra/angle/per/atom 20
read_restart ${fname}
#reset_timestep  0
variable box_size equal xhi-xlo

#GROUPS
#group  			nc1 type 1
#group  			7s  type 2 3 4 5 6 7
#group  			7sBreaking  type 2 3 4 5
#group  			7s_trimer  type 4

#BONDS COEFFICIENTS

bond_coeff     	1 	100.0	3.0
bond_coeff     	2 ${k_nc1} ${r0_nc1}
bond_coeff		3 ${k_7s} ${r0_7s}

#ANGLE COEFFICIENTS\n
''')
    fo.write('angle_coeff     1 	1.0	180.0\n')
    fo.write('angle_coeff     2 	1.0	155.0\n')
    fo.write('angle_coeff     3 	1.0	60.0\n')
    fo.write('variable x_stretch_Frac equal     '+str(Xstretch)+'\n')
    fo.write('variable x_stretch equal     $(v_x_stretch_Frac*v_box_size/2) \n')
    fo.write(
        '''

#NEIGHBOR LISTS
neighbor		${neigh_cutoff} bin 		#value = skin = extra distance beyond cutoff
neigh_modify	every 10 delay 0 check yes
variable box_size equal xhi-xlo

#PAIR STYLE (WCA)
pair_style     	hybrid/overlay lj/cut 1.0 zero 1.0
pair_coeff      * * lj/cut 0.0 0.0 ${lj_cutoff}
pair_coeff      * * zero ${zero_cutoff}

#SET TIME STEP
timestep	${tstep}

####-PRINTING STUFF-##############################

#COMPUTE BOND PROPERTIES
compute 		1 all property/local batom1 batom2 btype
compute 		4 all property/local atype aatom1 aatom2 aatom3
compute 		2 all pe bond
compute 		3 all pe angle
compute perAtomStress all stress/atom NULL
compute totalStressX all reduce sum c_perAtomStress[1]
compute totalStressY all reduce sum c_perAtomStress[2]
compute totalStressZ all reduce sum c_perAtomStress[3]
compute forceBond all bond/local force
compute forceBondAll all reduce ave c_forceBond

#DEFINITION OF VARIABLES FOR USE IN FIX PRINT
variable step	equal step
variable temp	equal temp
variable etot	equal etotal
variable ke		equal ke
variable peV	equal pe
variable press	equal press
variable peBond	equal c_2
variable peAngle	equal c_3
variable StressX equal c_totalStressX
variable StressY equal c_totalStressY
variable StressZ equal c_totalStressZ
variable Fbond equal c_forceBondAll

#THERMODYNAMIC INFO
thermo_style	custom step etotal ke pe temp press c_totalStressX c_totalStressY c_totalStressZ
thermo		${thermodump}
thermo_modify	norm no

##################################################
#NVT (NVE+LANGEVIN) RUN

fix 		fix_wall_lo 	all wall/lj126 zlo EDGE 1.0 1.0 ${lj_cutoff} #Lower wall
fix 		fix_wall_hi 	all wall/lj126 zhi EDGE 1.0 1.0 ${lj_cutoff} #Upper wall
fix 		fix_lengevin	all langevin ${mytemp} ${mytemp} ${damp} ${lseed}
fix	fNVE	all	nve 



#CREATE VELOCITIES USING MAXWELL-BOLTZMANN DISTRIBUTION
velocity 	all create ${mytemp} ${vseed} dist gaussian mom yes rot yes 
run 0
velocity	all scale ${mytemp}

fix		fix_print all print ${thermodump} "${step} ${etot} ${ke} ${peBond} ${peAngle} ${temp} ${press} ${StressX} ${StressY} ${StressZ} ${Fbond}" file ${thermofile} screen no title "step etot ke peBond peAngle temp press stressX stressY stressZ AvBondForce"

#Save configurations in linear time
dump		2 all custom ${trajdump} dumplin/dump.${simname}.*.lammpstrj type id xu yu zu
dump_modify	2 sort id
dump_modify	2 pad 10
dump_modify	2 format line "%d %d %.8f %.8f %.8f"

#Save bond information
dump 		3 all local ${trajdump} dumplin_bonds/bonds.${simname}.* index c_1[1] c_1[2] c_1[3] #bond index, id of atom 1 and 2 in the bond, bond type
dump_modify	3 format line "%d %0.0f %0.0f %0.0f"
dump_modify	3 pad 10


restart         ${trestart} restart/${simname}.*.restart

run 0

special_bonds lj 0.0 1.0 1.0

write_restart	restart.postmakebonds
fix                    fDeform all deform 1 x delta -$(v_x_stretch) $(v_x_stretch) remap x
run                    ${tstretch}
unfix                  fDeform
write_restart	restart.poststretch
                
run		${trelax}

''')
    
    fo.close()

    ################ write runscript ################
    runscriptfilename =  workingdir +"/runscript.sh"
    print('Writing runscript.sh\n')
    f = open(runscriptfilename, "w")
    f.write(
        '''#!/bin/bash
#
#SBATCH -J MakingBreakingOutputs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1000M 
#SBATCH --time=100:00:00 
#SBATCH --mail-user=ucapbbm@ucl.ac.uk
#SBATCH --mail-type=END
#SBATCH --export=NONE
#SBATCH --exclude=eta332
module load gcc
module load openmpi

mkdir dumplin dumplin_bonds restart
printf $RANDOM > lseed.dat
printf $RANDOM > vseed.dat

mpirun -np 1 /nfs/scistore15/saricgrp/bmeadowc/Scratch/lammps-15Jun2023/src/lmp_mpi -in collagen.in
\n''')
    f.close()

##SBATCH --constraint="epsilon|delta|beta|leonid|serbyn|gamma"     
if __name__ == "__main__":
    main()
