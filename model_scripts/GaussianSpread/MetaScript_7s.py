from __future__ import division, print_function
import numpy as np
#import matplotlib.pyplot as plt
import math
import time
import random
import sys, getopt
import os
import subprocess

def main():

    ################ Parameters ################
    trajdir = str(sys.argv[1]) 
    FracTagged = str(sys.argv[2]) 
    workingdir = str(sys.argv[1]) +'/GaussianSpread/'
    Frame = 101000 #31000 #101000
    Last =  5858000 #3100000 #12928000
    First=  4040000 #0000000 #2727000 #1414000 
    time = np.arange(First,Last+Frame*2,Frame*2)

    time_strList=[]
    for i in range(len(time)):
        time_str = str(time[i])
        if len(time_str)<5:
            time_strList.append('000000'+time_str)
        elif len(time_str)<6:
            time_strList.append('00000'+time_str)
        elif len(time_str)<7:
            time_strList.append('0000'+time_str)
        elif len(time_str)<8:
            time_strList.append('000'+time_str)
        else:
            time_strList.append('00'+time_str)
            

    ################ write in.local ################
    for t in range(len(time_strList)):
    ################ write runscript ################
        runscriptfilename =  workingdir +"/runscript"+time_strList[t]+".sh"
        print('Writing runscript.sh\n')
        f = open(runscriptfilename, "w")
        f.write(
            '''#!/bin/bash
#
#SBATCH -J Analysis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=500M 
#SBATCH --time=10:00:00 
#SBATCH --mail-user=ucapbbm@ucl.ac.uk
#SBATCH --mail-type=END
#SBATCH --export=NONE 
''')
        f.write('#SBATCH --chdir='+workingdir+'\n')
        f.write(
            '''

module load lammps/20220623b
module load python
pip install ovito
''')
        f.write('''#!/bin/bash\n
#
''')
        f.write('trajdir="'+trajdir+'"\n')
        f.write('filename="dump.run_T1.'+time_strList[t]+'.lammpstrj"\n')
        f.write('dataname=Gaussian'+time_strList[t]+'\n')
        f.write('fractagged='+FracTagged+'\n')
        f.write('python /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/GaussianSpread/PointSpreadSingle_7s.py ${trajdir} ${filename} ${dataname} ${fractagged}')
        
        f.close()
        command = 'sbatch '+runscriptfilename
        print(command)
        # Execute the command
        subprocess.run(command, shell=True)

##SBATCH --constraint="epsilon|delta|beta|leonid|serbyn|gamma"     
if __name__ == "__main__":
    main()
