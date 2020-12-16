import MDAnalysis as mda
from MDAnalysis.coordinates import XYZ
from MDAnalysis.analysis.dihedrals import Dihedral
import numpy as np

# For Gromacs:
    #$ gmx trjconv -f md.xtc -s md.tpr -b 980 -dt 100000 -center -pbc mol
    # select interesting molecule for cenetering and output
    # outfile is trajout.xtc
    


TOPOLOGY_FILE = 'gulA.itp'
TRAJECTORY_FILE = 'traj.xtc'
OUTPUT_FILE = 'gulA_ens_av_q' # +'frame_nr.nw'
ATOM_SELECTION = 'not name C1 C2 C3 C4 C5 O5' # atoms that dont get dihedral angles fixed

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
end
                                                     
esp
 spacing 0.02   # around 30 000 points for a methyl glycoside
 constrain 0 1  # constrain zero charge on atom 1 (haliphatic hydrogen)
 constrain equal 1 3 4 7 10 15 19 23    # constrain equal charges (zero charge on all hydrogens)
 constrain 0.264 2 # -O-CH3 carbon
 constrain -0.458 5 # -O-CH3 oxygen
 constrain equal 12 13 # carboxyl oxygens
 restrain hyperbolic 0.01   # restraint used in GLYCAM
 
end

                                             
driver                                                      
 maxiter 50     # maximum number of geometry optimization steps
end            

task scf optimize
task esp #ESP with restraints (RESP)
"""

mol = mda.Universe(TOPOLOGY_FILE,TRAJECTORY_FILE)
mol_deselect = mol.select_atoms(ATOM_SELECTION)

di = mol_deselect.dihedrals
all_di_angles = np.round_(Dihedral(di).run().angles,5).astype(str)
di_indices = (di.indices+1).astype(str) # dihedral indices
n_dih = len(di_indices) # number of dihedral angles in molecule
di_indices = np.append(np.array([['Torsion']]*n_dih),di_indices,axis=1)


names = mol.atoms.names.astype(str).reshape(-1,1) # atom names

print('#'*52)
print('Writing',mol.atoms.n_atoms,'atom coordinates in',mol.trajectory.n_frames,'frames')
print('Constraining',n_dih,'out of',len(mol.dihedrals.indices),'dihedral angles')
print('#'*52)
for i,name in enumerate(names):
    names[i] = name[0][0]


for i in range(mol.trajectory.n_frames):
    # Create geometry input 
    xyz = mol.coord.positions.astype(str) # coordiantes in Å
    xyz_save = np.append(names,xyz,axis=1) # input data .xyz file format
    

    
    # Create dihedral input to specify constant
    di_angles = np.append(di_indices,all_di_angles[i].reshape(-1,1),axis=1) # append angles 
    di_save = np.append(di_angles,np.array([['constant']]*n_dih),axis=1) # append 'constant' to specify constant values

        

    
    path = OUTPUT_FILE+str(i+101)+'.nw'

                
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
    if i != mol.trajectory.n_frames-1:
        mol.trajectory.next()
    else:
        n = i
        break
 
print('\n'  +str(n+1)+ ' NWChem input files written to ' + OUTPUT_FILE + 'i.nw')
    


