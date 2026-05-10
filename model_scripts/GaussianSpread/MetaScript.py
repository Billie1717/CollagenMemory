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
    workingdir = str(sys.argv[1]) +'/GaussianSpreadF'+str(FracTagged)+'/'
    Frame = 202000 #31000 #202000
    Last = 60802000 #41410000  #3100000 #20806000  #28482000 #15352000 #12928000 #38178000 #21210000 #41006000 #25250000 #16968000 #19190000 #21614000 #9292000 #4242000 #14140000 #11312000 #7474000 #6060000 #10807000 #6363000 #4185000 #3565000 #3720000 #4185000  #and go until 3720000
    First= 35148000 #28280000 #14140000 #15150000 #1212000 #21008000 #1212000 #20200000 #18180000 #22220000 #25250000 #23230000 #1212000 #21008000 #13736000 #1212000 #9090000 #10100000 #12120000 #4343000 #4242000 #1111000 #4242000 #1111000 #6060000 #1111000 #6060000 #3030000 #6060000 #4185000 #1085000 #465000 #620000 #1085000 #620000
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
#SBATCH --exclude=eta332
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
        f.write('python /nfs/scistore26/saricgrp/bmeadowc/Scratch/Collagen/NargessProject/GaussianSpread/PointSpreadSingle.py ${trajdir} ${filename} ${dataname} ${fractagged}')
        
        f.close()
        command = 'sbatch '+runscriptfilename
        print(command)
        # Execute the command
        subprocess.run(command, shell=True)

##SBATCH --constraint="epsilon|delta|beta|leonid|serbyn|gamma"     
if __name__ == "__main__":
    main()
