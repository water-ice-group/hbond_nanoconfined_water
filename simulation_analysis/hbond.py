import ase, warnings
from ase.io import read, write, iread
import numpy as np

warnings.filterwarnings('ignore')

def ase_get_OH_pairs(path, index="-1", OH_cutoff=1.2, fail_safely=False, FRAME_THRESHOLD=0.99):
    
    """
    This method takes the trajectory frames at `path` in the provided `index` range and returns numpy
    arrays with the oxygen and hydrogen coordinates arranged such that the oxygen with coordinates at
    index (t,O) of the `Os` array is covalently bonded (NOT Hydrogen bonded) to the hydrogens with 
    coordinates at index (t,O) of the `Hs` array (refer to the shapes below if this is confusing).

    `OH_cutoff` is the cutoff distance in Angstroms for an OH bond.

    This is a modified version that allows the program to throw away problematic frames if the
    `fail_safely` tag is set to True. The minimum acceptable proportion of accepted frames can be
    specified via the `FRAME_THRESHOLD` agument, which I've conservatively set to 99%.
    
    Important note: this uses ase's `get_distances` method, which is a bit slower than manually
    shifting the cell to account for PBCs of the hydrogen atoms.
    
    The shape of the resulting vectors are:
    
    `Os`: (num_frames, num_Os, 3)
    `Hs`: (num_frames, num_Os, 2, 3)
    `cells`: (num_frames, 3, 3)
    
    """
    
    Os = []
    Hs = []
    cells = []

    all_frame_count = 0
    
    # This iterates through all of the specified frames in the provided path
    for (i, struct) in enumerate(iread(path, index=index)):
        
        all_frame_count += 1

        atoms = np.array(struct.get_chemical_symbols())
        positions = struct.get_positions()
        Os_idxs = np.where(atoms == 'O')[0]
        Hs_idxs = np.where(atoms == 'H')[0]
        
        local_Hs = []

        frame_success = True
        
        for O_idx in Os_idxs:
            
            OH_dists = struct.get_distances(O_idx, Hs_idxs, mic=True)
            bonded_Hs_idxs = np.where(OH_dists < OH_cutoff)[0]
            
            # This is triggered if there is an oxygen without exactly two hydrogen atoms bonded to it
            if (np.where(OH_dists < OH_cutoff)[0].shape[0] != 2):
                
                frame_success = False

                if fail_safely == True:
                    break;
                else:
                    print(np.sort(OH_dists)[:5])
                    raise Exception(
                        "PR: Something wrong with detecting OH bonds (oxygen without 2 hydrogens at frame {})".format(i)
                    )
                
            curr_Hs = np.zeros((2,3))
            OH_vecs = struct.get_distances(O_idx, Hs_idxs[bonded_Hs_idxs], mic=True, vector=True)
            curr_Hs[0,:] = positions[O_idx] + OH_vecs[0]
            curr_Hs[1,:] = positions[O_idx] + OH_vecs[1]
            local_Hs.append(curr_Hs)

        if frame_success == False:
            continue;

        Os.append(positions[np.where(atoms == 'O')])
        Hs.append(local_Hs)
        cells.append(struct.cell)
        
    accept_ratio = len(Os) / all_frame_count
    if (accept_ratio < FRAME_THRESHOLD):
        raise Exception(
            "PR: only {} percent of frames were kept... not enough to proceed without concern".format(accept_ratio)
        )

    return np.array(Os), np.array(Hs), np.array(cells)

def ase_H_bond_atoms(path, index="-1", OH_cutoff=1.2, OO_CUTOFF=3.5, OHO_CUTOFF=30):
    
    """
    This method takes as arguments the numpy arrays `Os` and `Hs` with dimensions (num_frames, num_Os, 3) and
    (num_frames, num_Os, 2, 3), which contain the Cartesian coordinates for the oxygens and hydrogens in the
    system such that the oxygens are lined up with the hydrogens that they are bonded to.
    
    The returned value `traj` is the trajectory for the ratio of hydrogen bonds as determined by the `OO_CUTOFF`
    (in Angstroms) and the `OHO_CUTOFF` (in degrees) provided as arguments. The ratio is the ratio of the number
    or hydrogen bonds to the number of oxygens in the frame.
    
    `traj`:          (num_frames, )
    `O_O_bonds`:     (num_frames, * ):
        for each frame, there is a *Python* list of tuples (i,j) where oxygens i and j are H-bonded
    `H_bonds`:       (num_frames, * ):
        for each frame, there is a *Python* list of indexes of the hydrogens that are H-bonded
    
    """
    
    traj = []
    O_O_bonds = []
    H_bonds = []
    
    for (frame_idx, struct) in enumerate(iread(path, index=index)):
        
        #if frame_idx % 100 == 0:
        #    print("Starting frame {}".format(i))
        #print(frame_idx)
        
        atoms = np.array(struct.get_chemical_symbols())
        Os_idxs = np.where(atoms == 'O')[0]
        Hs_idxs = np.where(atoms == 'H')[0]
        
        num_Os = Os_idxs.shape[0]
        num_Hs = Hs_idxs.shape[0]
        
        num_bonds = 0
        curr_O_O_bonds = []
        curr_H_bonds = []
        
        for O_idx in Os_idxs:
            
            OO_dists = struct.get_distances(O_idx, Os_idxs, mic=True)
            OH_dists = struct.get_distances(O_idx, Hs_idxs, mic=True)
            
            close_Os = Os_idxs[np.where((OO_dists < OO_CUTOFF) & (OO_dists > 1e-2))[0]]
            close_Hs = Hs_idxs[np.where(OH_dists < OH_cutoff)[0]]
            
            if close_Os.shape[0] == 0:
                continue;
            
            OO_vecs = struct.get_distances(O_idx, close_Os, mic=True, vector=True)
            OH_vecs = struct.get_distances(O_idx, close_Hs, mic=True, vector=True)
            
            for i in range(OO_vecs.shape[0]):
                for j in range(OH_vecs.shape[0]):

                    OO = OO_vecs[i,:]
                    OH = OH_vecs[j,:]
                    angle = np.arccos(np.dot(OO, OH) / (np.linalg.norm(OO) * np.linalg.norm(OH))) * 180 / np.pi

                    if angle < OHO_CUTOFF:
                        num_bonds += 1
                        curr_O_O_bonds.append((O_idx, close_Os[i]))
                        curr_H_bonds.append(close_Hs[j])
            
        traj.append(num_bonds / num_Os)
        O_O_bonds.append(curr_O_O_bonds)
        H_bonds.append(curr_H_bonds)
        
    return np.array(traj), O_O_bonds, H_bonds
