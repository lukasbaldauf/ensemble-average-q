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

###### The following are the results from 200 snapshots of a 200 ns GROMACS simulation using ACPYPE of alpha-L-GulA with a methyl group attached to O1 at the HF/6-31G* level of theory with all dihedral angles frozen to their MD values.  

| Atom type | Mean charge | Std. dev. |
| --- | --- | --- |
   H3  |    0.000   |   0.000
  CH3  |    0.203   |   0.015
   H2  |    0.000   |   0.000
   H3  |    0.000   |   0.000
   O   |   -0.431   |   0.032
   C1  |    0.436   |   0.045
   H1  |    0.000   |   0.000
   O5  |   -0.526   |   0.042
   C5  |    0.221   |   0.076
   H5  |    0.000   |   0.000
   C6  |    0.891   |   0.061
  O6B  |   -0.840   |   0.021
  O6A  |   -0.840   |   0.021
   C4  |    0.238   |   0.096
   H4  |    0.000   |   0.000
terminal O4  |   -0.769   |   0.046
  H4O  |    0.426   |   0.043
   C3  |    0.347   |   0.079
   H3  |    0.000   |   0.000
   O3  |   -0.765   |   0.035
  H3O  |    0.410   |   0.022
   C2  |    0.353   |   0.088
   H2  |    0.000   |   0.000
   O2  |   -0.764   |   0.042
  H2O  |    0.410   |   0.025
  O4*  |   
  H1O^ |
  O1^  |
