import random
import sys
sys.path.append("..")
#import pandas
import numpy as np
from ovito.io import *
from ovito.modifiers import *
from ovito.io import export_file
from ovito.data import *



def spreadzedNoStretch(x1,y,means,sigma1,sigma2,sigma3):
    z1 = np.zeros((len(y),len(x1)))
    z2 = np.zeros((len(y),len(x1)))
    z3 = np.zeros((len(y),len(x1)))
    print(y[0],y[-1],abs((y[-1]+y[0])))
    for i in range(len(x1)):
        #if i%(int(len(x1)/20))==0:
        #    print(i*100/len(x1),"% done")
        for j in range(len(y)):
            for m in range(len(means[0,:])):
                z1[j,i] += np.exp(-(((x1[i]-means[0,m])/sigma1)**2+((y[j]-means[1,m])/sigma1)**2))
                z2[j,i] += np.exp(-(((x1[i]-means[0,m])/sigma2)**2+((y[j]-means[1,m])/sigma2)**2))
                z3[j,i] += np.exp(-(((x1[i]-means[0,m])/sigma3)**2+((y[j]-means[1,m])/sigma3)**2))
                if means[1,m] > abs((y[-1]+y[0])):
                    print("got here")
                    meansT = means[1,m]-y[-1]
                    z1[j,i] += np.exp(-(((x1[i]-means[0,m])/sigma1)**2+((y[j]-meansT)/sigma1)**2))
                    z2[j,i] += np.exp(-(((x1[i]-means[0,m])/sigma2)**2+((y[j]-meansT)/sigma2)**2))
                    z3[j,i] += np.exp(-(((x1[i]-means[0,m])/sigma3)**2+((y[j]-meansT)/sigma3)**2))
    return z1,z2,z3


def main():
    
    xhi = 1.6454482700000000*100 #331.976 #
    xlo = 0 
    yhi = 1.6454482700000000*100 #331.976 #
    ylo = 0
    xtot = xhi-xlo
    
    
    NumPixels = 200
    sigma1 = 2.5
    sigma2 = 3.0
    sigma3 = 3.5
    trajdir = str(sys.argv[1]) +'/dumplin/'
    filename = str(sys.argv[2])
    FracTagged = 1.0 #float(sys.argv[4])
    datadir  = str(sys.argv[1]) +'/GaussianSpread/'
    dataname= str(sys.argv[3]) 
    
    print(trajdir+filename)
    pipeline = import_file(trajdir+filename)
    pipeline.modifiers.append(WrapPeriodicImagesModifier())
    data = pipeline.compute()
    positions = data.particles['Position']
    types = data.particles['Particle Type']
    xcoords = positions[:, 0]
    ycoords = positions[:, 1]
    z_coordinates = positions[:, 2]
    print(len(xcoords))
    NumCoords = len(xcoords) #3215 #len(xcoords)
    
    # Choose subset of indices j where type > 1
    eligible_indices = [j for j in range(NumCoords) if types[j] == 9] #only choosing broken type 7s to tag 
    num_to_tag = int(FracTagged * len(eligible_indices))
    tagged_indices = random.sample(eligible_indices, num_to_tag)  # No repeats
    
    # Allocate coords array
    coords = np.zeros((2, num_to_tag))
    
    for i, j in enumerate(tagged_indices):
        coords[0, i] = xcoords[j]
        coords[1, i] = ycoords[j]
        
    MeshX = np.linspace(xlo,xhi,NumPixels) #np.arange(xlo,xhi,0.5)
    MeshY = np.linspace(ylo,yhi,NumPixels)
    zspread1 = spreadzedNoStretch(MeshX,MeshY,coords,sigma1,sigma2,sigma3)
    
    outfile1 = datadir+ dataname
    f1 = open(outfile1+"sigma"+str(sigma1)+".txt",'w')    
    f2 = open(outfile1+"sigma"+str(sigma2)+".txt",'w')   
    f3 = open(outfile1+"sigma"+str(sigma3)+".txt",'w')  
    
    for j in range(len(MeshY)):
        if j!=0:
            f1.write('\n')
            f2.write('\n')
            f3.write('\n')
        for i in range(len(MeshX)):
            f1.write(str(zspread1[0][j][i])+" ")
            f2.write(str(zspread1[1][j][i])+" ")
            f3.write(str(zspread1[2][j][i])+" ")
    f1.close()
    f2.close()
    f3.close()
        
if __name__ == "__main__":
    main()
