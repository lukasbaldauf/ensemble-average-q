# ensemble-average-q

#### Create NWChem input files for the calculation of ensemble-average RESP atomic partial charges according to GLYCAMs ensemble-averaging procedure.

#### NOTE!
The number of NWChem input files is determined by the number of snapshots in the trajectory. NWChem options should be changed in the python script according to your molecule of interest. 

Required python packages: \
Numpy \
MDAnalysis 

###### The following are the results from a 200 ns GROMACS simulation of alpha-L-GulA-OMe in the 1C4 ring conformation and beta-D-manA-OMe in the 4C1 ring conformation, using ACPYPE to make GLYCAM compatible with GROMACS. 200 equally spaced snapshots were extracted from the trajectory and their geometry optimized the HF/6-31G* level of theory using NWChem. All dihedral angles were frozen to their MD snapshot values.
| Atom type |av. q gulA |std. q gulA |av. q manA |std. q manA| Comment |
| --- | --- | --- | --- | --- | --- |
H1  |   0.000 |  0.000 |  0.000 | 0.000 |            
CH3 |   0.264 |  0.000 |  0.264 | 0.000 |            
H2  |   0.000 |  0.000 |  0.000 | 0.000 |            
H3  |   0.000 |  0.000 |  0.000 | 0.000 |            
 O  |  -0.458 |  0.000 | -0.458 | 0.000 |            
C1  |   0.333 |  0.049 |  0.251 | 0.044 |            
H1  |   0.000 |  0.000 |  0.000 | 0.000 |            
O5  |  -0.499 |  0.043 | -0.369 | 0.042 |            
C5  |   0.191 |  0.069 |  0.017 | 0.076 |            
H5  |   0.000 |  0.000 |  0.000 | 0.000 |            
C6  |   0.906 |  0.053 |  0.928 | 0.051 |            
O6B |  -0.847 |  0.018 | -0.842 | 0.022 |            
O6A |  -0.847 |  0.018 | -0.842 | 0.022 |            
C4  |   0.284 |  0.094 |  0.479 | 0.110 |            
H4  |   0.000 |  0.000 |  0.000 | 0.000 |            
O4  |  -0.786 |  0.049 | -0.766 | 0.039 | terminal unit                              
H4O |   0.430 |  0.044 |  0.412 | 0.022 |            
C3  |   0.298 |  0.083 |  0.217 | 0.100 |            
H3  |   0.000 |  0.000 |  0.000 | 0.000 |            
O3  |  -0.747 |  0.040 | -0.711 | 0.035 |            
H3O |   0.403 |  0.025 |  0.406 | 0.019 |            
C2  |   0.445 |  0.097 |  0.264 | 0.070 |            
H2  |   0.000 |  0.000 |  0.000 | 0.000 |            
O2  |  -0.783 |  0.041 | -0.654 | 0.032 |            
H2O |   0.413 |  0.022 |  0.404 | 0.017 |            
O4  |  -0.550 |   -    | -0.548 |   -   | = O4 + H4O + O + CH3                              
H1O |   0.445 |   -    |  0.445 |   -   | glycam default for H1O 
O1  |  -0.639 |   -    | -0.639 |   -   | glycam default for O1  



Assuming you have run a GROMACS simulation of your molecule, you can run 

    gmx trjconv -f md.xtc -s md.tpr -b 980 -dt 1000 -center -pbc mol
    
with the desired spacing of intervals and select your interesting molecule for cenetering and output. The outfile (trajout.xtc) and a topology file of your molecule, e.g. an .itp file is used as input for ensemble-average-q.py to create the NWChem files with all dihedral angles frozen to the values in the trajectory. 
