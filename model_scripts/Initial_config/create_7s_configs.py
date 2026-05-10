#The unit of length is the particle diameter.

from numpy import *
from numpy.random import *
import random

box_x=float(input("Enter box length (X = Y): "))
bond_length=float(input("Enter bond length: "))
density=float(input("Enter density: "))  # molecules per unit volume
frac=float(input("Enter fraction of 7s: "))  # molecules per unit volume

lattice_z = 2.
nmol_side_z = 5

box_y = box_x
box_z=lattice_z*(nmol_side_z+1)  # fixed z box


volume = box_x * box_y * box_z
nmol = int(volume * density)

nmol_side_z = 5
nmol_side_x = int((nmol / nmol_side_z) ** 0.5)
nmol_side_y = nmol_side_x

nmol = nmol_side_x * nmol_side_y * nmol_side_z  # recompute to ensure it’s a perfect grid


nmol=nmol_side_x*nmol_side_y*nmol_side_z
natoms=2*nmol


lattice_x = box_x / nmol_side_x
lattice_y = box_y / nmol_side_y

with open('dumbbells_7s_frac%.2f_n%d_bl%.2f_Lx%.2f_rho%.3e.lammpsdata'%(frac,nmol,bond_length,box_x,density),'w') as fout:
    fout.write("LAMMPS data file\n\n%d atoms\n%d bonds\n\n9 atom types\n1 bond types\n\n0.0 %f xlo xhi\n0.0 %f ylo yhi\n0.0 %f zlo zhi\n\nMasses\n\n1 1\n2 1\n3 1\n4 1\n5 1\n6 1\n7 1\n8 1\n9 1\n\nBond Coeffs # harmonic\n\n1 10 3\n\nAtoms\n\n"%(natoms,nmol,box_x,box_y,box_z))
    atom_id=0
    mol_id=0
    for ix in range(nmol_side_x):
        for iy in range(nmol_side_y):
            for iz in range(nmol_side_z):
                # make randomly a fraction of them 
                r = random.random()  # generates float in [0, 1)
                if r < frac:
                    EndType = 2  
                else:
                    EndType = 9
                mol_id+=1
                r1=array([ix*lattice_x,iy*lattice_y,(iz+1)*lattice_z])
                atom_id+=1
                fout.write("%d %d 1 %.15f %.15f %.15f\n"%(atom_id,mol_id,r1[0],r1[1],r1[2]))
                r2=r1+bond_length*array([1,1,0.5])/sqrt(1*1+1*1+0.5*0.5)
                atom_id+=1
                fout.write("%d %d %d %.15f %.15f %.15f\n"%(atom_id,mol_id,EndType,r2[0],r2[1],r2[2]))
    mol_id=0
    atom_id=0
    fout.write("\nBonds\n\n")    
    for ix in range(nmol_side_x):
        for iy in range(nmol_side_y):
            for iz in range(nmol_side_z):
                mol_id+=1
                atom_id+=1
                fout.write("%d 1 %d %d\n"%(mol_id,atom_id,atom_id+1))
                atom_id+=1
        
    


