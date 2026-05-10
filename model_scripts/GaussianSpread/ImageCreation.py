import sys
sys.path.append("..")
import pandas
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image
import io
import os

clrs = ['#a6611a','#dfc27d','#80cdc1','#018571','#f5f5f5']
datadir  = str(sys.argv[1]) +'/'
plotsdir  = datadir+'/Plots/'

sigmas = "2.5 3.0 3.5".split()
Frame = 202000 #202000 #202000 #101000 #31000 #101000 #31000 #101000 #31000 #101000 #31000
Last =   60802000 #20806000 #15352000 #3100000 #5858000 #37370000 #5858000 #3100000 #40804000 #12928000 #3100000 #12928000 #3100000 #38178000 #41006000 #25250000 #12120000 #21614000 #13938000 #21210000 #21614000 #21210000 #12928000 #21210000 #6666000 #747400014140000 # #6161000 #4061000 #3720000 #4185000  #3720000 # 4185000 
First= 35148000 #14140000 #100000-Frame*20 #2790000 #000000 #4040000 #21008000 #4040000 #000000 #20200000 #9090000 #1212000 # 465000 #899000 #1212000 #0000000 #23230000 #20200000 #18180000 #22220000 #25250000 #23230000 #21008000 #13736000 #9090000 ##10100000 #12120000 #1212000 #2727000 #10100000 #1212000 #4242000 #1111000 #10850004242000 # #620000 #465000 #1085000 #465000 # 1085000  #620000  # 1085000 
time = np.arange(First,Last+Frame*2,Frame*2)
#time = np.arange(First,Last,Frame*2)
print('Time',time)
seed = "1".split()
vvmax = 30 #13 #30 #16 #30 #12 #16 [16 for dens0.05 stretched 100% networks, 12 for 150% networks, 30 for not stretched]
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
fig2,ax2 = plt.subplots(1,len(sigmas),figsize = (8*len(sigmas),6))

for ss in range(len(seed)):
    for t in range(len(time_strList)):
        
        fig1,ax1 = plt.subplots(1,figsize = (4.22,4.25))
        
        sigmachosen = "3.5 3.0".split() #1.5 2.0 
        for s in range(1):
            fname1 = 'Gaussian'+time_strList[t]+"sigma"+str(sigmachosen[s])+".txt"
            if os.path.exists(datadir+fname1):
                #df = pandas.read_csv(filepath, sep=' ', header=None)
            
                df = pandas.read_csv(datadir+fname1, sep = ' ',header =None)
                print(fname1)
                 #12 
                sns.heatmap(df,ax=ax1, cmap="Greys_r",cbar=False,vmin = 0,vmax = vvmax,alpha=(1)) 
                ax1.invert_yaxis()
                plt.axis('off')
                plt.gca().set_frame_on(False)
                plt.gcf().set_facecolor('white')
                # Create your plot
                # Convert the Matplotlib figure to a PIL image
                buf = io.BytesIO()
                plt.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
                buf.seek(0)
                img = Image.open(buf)
    
                # Convert to RGB mode (optional, only if the image is in RGBA mode)
                if img.mode == 'RGBA':
                    img = img.convert('RGB')
                plotname = 'Gaussian'+time_strList[t]+"vvmax_const"+str(vvmax)+"sigma"+str(sigmachosen[s])+".tif"
                # Save the image without transparency
                print("saving image for time"+str(time_strList[t]))
                img.save(plotsdir+plotname,bbox_inches='tight', pad_inches=0)
            else:
                print("Didn't find file"+datadir+fname1)



