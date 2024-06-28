import sys, warnings
import numpy as np
from ase.io import iread
from hbond import *

def get_OH_dict(struct):
    """
    Calculates the phi and theta angles for the water molecules

    rmax:        maximum radial distance
    nbins:       number of bins to use
    fend:        how many frames to use from the end of the trajectory
    mode:        if this is "save", then this will save the results

    """

    warnings.filterwarnings('ignore')

    path = "{}/simulation.pos_0.pdb"
    Ts = [str(20*(i+1)) for i in range(30)]
    index = "50:850"

    OH_dict = {}
    cell_dict = {}
    H_bond_dict = {}

    for T in Ts:

        print("Starting all frames for {} K for {} phase".format(T, struct))

        try:
        
            Os, Hs, cells = ase_get_OH_pairs(path.format(T), index=index, fail_safely=True)
            traj, O_bonds, H_bonds = ase_H_bond_atoms(path.format(T), index=index)
        
            OH_dict[T] = (Os, Hs)
            cell_dict[T] = cells
            H_bond_dict[T] = (traj, O_bonds, H_bonds)

        except:
           
            print("Failed {} K for {} phase".format(T, struct)) 
            continue;
        
    np.save("{}_OH_dict.npy".format(struct), OH_dict)
    np.save("{}_cell_dict.npy".format(struct), cell_dict)
    np.save("{}_H_bond_dict.npy".format(struct), H_bond_dict)

struct = sys.argv[1]

print("Starting {} phase".format(struct))
get_OH_dict(struct)
print("--------------------\n")
