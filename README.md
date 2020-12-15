# ensemble-average-q

#### Create NWChem input files for the calculation of restrained electrostatic potential (RESP) atomic partial charges according to GLYCAMs ensemble-averaging procedure.

Assuming you have run a GROMACS simulation of your molecule, you can run 

    gmx trjconv -f md.xtc -s md.tpr -b 980 -dt 1000 -center -pbc mol
    
with the desired spacing of intervals and select your interesting molecule for cenetering and output. The outfile (trajout.xtc) and a topology file of your molecule, e.g. an .itp file is used as input for ensemble-average-q.py to create the NWChem files with all dihedral angles frozen to the values in the trajectory. 

#### NOTE!
The number of NWChem input files is determined by the number of snapshots in the trajectory. NWChem options should be changed in the python script according to your molecule of interest. 



Required python packages: \
Numpy \
MDAnalysis 

###### The following are the results from a 200 ns GROMACS simulation of alpha-L-GulA-OMe using ACPYPE to make GLYCAM compatible with GROMACS. 200 equally spaced snapshots were extracted from the trajectory and their geometry optimized the HF/6-31G* level of theory with all dihedral angles frozen to their MD values using NWChem. 

| Atom type | Mean charge | Std. dev. | Comment |
| --- | --- | --- | --- |
H1  |   0.000 |  0.000 |                            
CH3 |   0.264 |  0.000 |                            
H2  |   0.000 |  0.000 |                            
H3  |   0.000 |  0.000 |                            
 O  |  -0.458 |  0.000 |                            
C1  |   0.333 |  0.049 |                            
H1  |   0.000 |  0.000 |                            
O5  |  -0.499 |  0.043 |                            
C5  |   0.191 |  0.069 |                            
H5  |   0.000 |  0.000 |                            
C6  |   0.906 |  0.053 |                            
O6B |  -0.847 |  0.018 |                            
O6A |  -0.847 |  0.018 |                            
C4  |   0.284 |  0.094 |                            
H4  |   0.000 |  0.000 |                            
O4  |  -0.786 |  0.049 |terminal unit               
H4O |   0.430 |  0.044 |                            
C3  |   0.298 |  0.083 |                            
H3  |   0.000 |  0.000 |                            
O3  |  -0.747 |  0.040 |                            
H3O |   0.403 |  0.025 |                            
C2  |   0.445 |  0.097 |                            
H2  |   0.000 |  0.000 |                            
O2  |  -0.783 |  0.041 |                            
H2O |   0.413 |  0.022 |                            
O4  |  -0.550 |   -    |= O4 + H4O + O + CH3           
H1O |   0.445 |   -    |taken from glycam default 
O1  |  -0.639 |   -    |taken from glycam default 
