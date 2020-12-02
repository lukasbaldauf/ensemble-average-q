import MDAnalysis as mda
from MDAnalysis.coordinates import XYZ
from MDAnalysis.analysis.dihedrals import Dihedral
import numpy as np

# For Gromacs:
    #$ gmx trjconv -f md.xtc -s md.tpr -b 980 -dt 1000 -center -pbc mol
    # select interesting molecule for cenetering and output
    # outfile is trajout.xtc
    
# NWCHEM input options
start = """
start 1st-job # job name

charge -1 # total molecular charge
geometry units angstroms print xyz autosym # extracted xyz coordinates are put below this line
"""
end = """
basis                                  
  O  library 6-31G*                        
  H  library 6-31G*                                     
  C  library 6-31G*                                        
end
                                                            
scf                                                         
 maxiter 50     # maximum number of scf iterations
 direct         # do not store integrals                 
end
                                                     
esp
 spacing 0.02   # around 30 000 points for a methyl glycoside
 constrain 0 1  # constrain zero charge on atom 1 (haliphatic hydrogen)
 constrain equal 1 3 4 7 10 15 19 23    # constrain equal charges (zero charge on all hydrogens)
 constrain equal 12 13      # constrain equal charges on carboxyl groups
 restrain hyperbolic 0.01   # restraint used in GLYCAM
end

                                             
driver                                                      
 maxiter 50     # maximum number of geometry optimization steps                                           
end                                                         
                                                            
task scf optimize  
task esp #ESP with restraints (RESP)
"""

mol = mda.Universe('gulA.itp','980ps-1ns-spacing.xtc')
di = mol.dihedrals 
all_di_angles = np.round_(Dihedral(di).run().angles,5).astype(str)
names = mol.atoms.names.astype(str).reshape(-1,1) # atom names
for i,name in enumerate(names):
    names[i] = name[0][0]
    
for i in range(mol.trajectory.n_frames):
    # Create geometry input 
    xyz = mol.coord.positions.astype(str) # coordiantes in Å
    
    xyz_save = np.append(names,xyz,axis=1) # input data .xyz file format
    
    
    
    
    # Create dihedral input to specify constant
    dihedral_constant = 'All' # angles to hold constant during geometry optimization
    di_indices = (di.indices+1).astype(str) # dihedral indices
    n_dih = len(di_indices) # number of dihedral angles in molecule
    di_indices = np.append(np.array([['Torsion']]*n_dih),di_indices,axis=1) # append 'Torsion' for nwchem format it is a dihedral angle
    di_angles = np.append(di_indices,all_di_angles[i].reshape(-1,1),axis=1) # append angles 
    if dihedral_constant == 'All':
        di_save = np.append(di_angles,np.array([['constant']]*n_dih),axis=1) # append 'constant' to specify constant values
    else:
        None



            

    
    path = 'nwchem_ens-av-q_input'+str(i+1)+'.nw'

                
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
    mol.trajectory.next()   
    
    
    


