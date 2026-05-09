Normal ablation size:

## For larger sims for ablations:

Number of molecules = 40000
Bond lenggth = 3.0
Density = 0.05

## For massive ablation sims:

Number of molecules = 80000
Bond lenggth = 3.0
Density = 0.05


Then need to change file to have header:


...
7 atom types
...

Masses

1 1
2 1
3 1
4 1
5 1
6 1
7 1

Bond Coeffs # harmonic

1 10 3

## For Large alignment sims:

Number of molecules = 20000 Bond lenggth = 3.0 Density = 0.05

## For Denser alignment sims:

Number of molecules = 40000 Bond lenggth = 3.0 Density = 0.1

## for 7S pertubrations:

Enter box length (X = Y): 164
Enter bond length: 3.0
Enter density: 0.05
Enter fraction of 7s: 0.44

OR 
Enter box length (X = Y): 164
Enter bond length: 3.0
Enter density: 0.05
Enter fraction of 7s: 0.25


