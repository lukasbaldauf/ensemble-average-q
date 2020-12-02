import MDAnalysis as mda
from MDAnalysis.coordinates import XYZ
from MDAnalysis.analysis.dihedrals import Dihedral
import numpy as np
start = """start 1st-job # Job name

charge -1
geometry units angstroms print xyz autosym # Cartesian coordinates of atoms
"""
end = """basis # Name of basis set                                   
  O  library STO-2G                                     
  H  library STO-2G                                     
  C  library STO-2G                                        
end                                                         
                                                            
scf                                                         
 maxiter 50       
 direct                                          
end
                                                     
esp
 spacing 0.02   # around 30 000 points
 constrain 0 1  # constrain zero charge on atom 1
 constrain equal 1 3 4 7 10 15 19 23    # constrain equal charges (zero)
 constrain equal 12 13      # constrain equal charges carboxyl
 restrain hyperbolic 0.01
end

                                             
driver                                                      
 maxiter 300                                                
end                                                         
                                                            
task scf optimize 
task esp"""

mol = mda.Universe('test-xyz1.xyz',topology_format='XYZ',guess_bonds=True)

# Create geometry input 
xyz = mol.coord.positions.astype(str) # coordiantes in Å
names = mol.atoms.names.astype(str).reshape(-1,1) # atom names
xyz_save = np.append(names,xyz,axis=1) # input data .xyz file format

di = mol.dihedrals 

# Create dihedral input to specify constant
dihedral_constant = 'All' # angles to hold constant during geometry optimization
di_indices = (di.indices+1).astype(str) # dihedral indices
n_dih = len(di_indices) # number of dihedral angles in molecule
di_indices = np.append(np.array([['Torsion']]*n_dih),di_indices,axis=1) # append 'Torsion' for nwchem format it is a dihedral angle
di_angles = np.append(di_indices,np.round_(Dihedral(di).run().angles,5).astype(str).T,axis=1) # append angles 
if dihedral_constant == 'All':
    di_save = np.append(di_angles,np.array([['constant']]*n_dih),axis=1) # append 'constant' to specify constant values
else:
    None



            
for i in range(1):
    
    path = 'program-input'+str(i)+'.nw'

                
    writefile = open(path,'w')
    writefile.write(start+'\n')
    
    
    writefile = open(path,'ab')
    np.savetxt(writefile,xyz_save,fmt=['%8s','%16s','%16s','%16s'])
    writefile.close()
    
    writefile = open(path,'a')
    writefile.write(' '*4+'zcoord\n')
    writefile.close()
    
    writefile = open(path,'ab')
    np.savetxt(writefile,di_save,fmt=['%15s','%4s','%4s','%4s','%4s','%12s','%12s'])
    writefile.close()
    
    writefile = open(path,'a')
    writefile.write(' '*4+'end\n'+'end\n')
    writefile.write(end)
    writefile.close()
    
    
    


