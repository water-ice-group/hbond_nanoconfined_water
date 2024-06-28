import numpy as np

# Input paths
H_bond_path = "../raw_thermo/{}_H_bond_dict.npy"
OH_path = "../raw_thermo/{}_OH_dict.npy"

# Output path
probs_path = "{}_{}K_2D_probs.txt"
nH_path = "{}_{}K_2D_nH.txt"

NUM_BINS = 24
OH_cutoff=1.2

phi_range = (-np.pi/2.0, np.pi/2.0)
theta_range = (0.0, np.pi)

def shift(idx):
    return (np.pi*idx) / NUM_BINS

"""
The below is written to be flexible wrt the # of O's in anticipation of larger unit cells.
Format for each element in the `structs` array is (`struct_name`, `num_Os`, `Ts[3]`).
"""

structs = [
    ("r2", 144, [100,280,360,380]),
    ("h", 144, [20, 100, 160]),
    ("p", 144, [20, 100, 160]),
    ("big_r2", 576, [100,280,360]),
]

def convert(num, big):
    if big:
        return int((num-2)/3)
    else:
        return int(num/3)

for (struct, num_Os, Ts) in structs:

    if struct == "big_r2":
        big = True
    else:
        big = False
    
    H_bond_dict = np.load(H_bond_path.format(struct), allow_pickle=True).item()
    OH_dict = np.load(OH_path.format(struct), allow_pickle=True).item()

    for T in Ts:

        """
        This chunk computes `H_coord`, which is a ( traj_len, num_Os ) shape matrix
        that shows how many H bonds each water participates in at a given point
        in the trajectory. Note that this chunk doesn't care whether these are
        accepted or donated H bonds. Hence, it is just the H bond coordination.
        """

        _, donor_list, _ = H_bond_dict[str(T)]
        traj_len = len(donor_list)

        H_coord = np.zeros((traj_len, num_Os))

        for (i,step) in enumerate(donor_list):
            for (d,a) in step:
                O_d, O_a = convert(d,big) , convert(a,big)
                H_coord[i,O_d] += 1
                H_coord[i,O_a] += 1

        """
        We first compute `thetas` which is a ( traj_len, num_Os ) shape matrix
        that shows the theta angle of each water molecule at a given time point,
        and then we compute a similar `phis` matrix.
        """

        Os, Hs = OH_dict[str(T)]

        traj_len = Os.shape[0]
        thetas = np.zeros((traj_len, num_Os))

        for t in range(traj_len):
            for O_i in range(num_Os):

                OH0 = Hs[t,O_i,0,:] - Os[t,O_i,:]
                OH1 = Hs[t,O_i,1,:] - Os[t,O_i,:]

                if Hs[0,O_i,0,0] > Hs[0,O_i,1,0]:
                    norm = np.cross(OH1, OH0)
                    norm = norm / np.linalg.norm(norm)
                else:
                    norm = np.cross(OH0, OH1)
                    norm = norm / np.linalg.norm(norm)

                if norm[0] >= 0:
                    theta = np.arctan(norm[2] / norm[0])
                else:
                    theta = np.arctan(norm[2] / norm[0]) + np.pi

                thetas[t,O_i] = theta
    
        reshaped_Os = Os.reshape((Os.shape[0], Os.shape[1], 1, Os.shape[2]))

        # The oxygens and hydrogens provided as arguments should have PBC issues taken care
        # of already, i.e., the Cartesian distance between each oxygen and its corresponding
        # hydrogens should already be close enough to be seen as bonded.
        OH_dists = np.linalg.norm(reshaped_Os - Hs, axis=3)
        if np.max(OH_dists) > OH_cutoff:
            raise Exception("PR: Os and Hs provided to get_PR_OP are not properly aligned to OH bonds")

        OH_vectors = (Hs - reshaped_Os) / np.linalg.norm(Hs - reshaped_Os, axis=3, keepdims=True)
        rot_vectors = np.average(OH_vectors, axis=2) / np.linalg.norm(np.average(OH_vectors, axis=2), axis=2, keepdims=True)
        z_unit = np.array([0.0, 0.0, 1.0])
        phis = np.pi / 2.0 - np.arccos(np.dot(rot_vectors, z_unit))

        """
        Now we compute `num_pts` and `avgs`, which are ( NUM_BINS , NUM_BINS ) shape
        matrices that contain the number of water molecules in each 2D bin along
        phi and theta as computed above. The first index corresponds to phi, and the
        second index corresponds to theta.
        """

        num_pts = np.zeros((NUM_BINS,NUM_BINS))
        avgs = np.zeros((NUM_BINS,NUM_BINS))

        H_coord = H_coord.flatten()
        phis = phis.flatten()
        thetas = (thetas.flatten() + 2.0 * np.pi) % np.pi

        for i in range(NUM_BINS):
            for j in range(NUM_BINS):

                filtered_data1 = H_coord[np.where(
                    (phis > phi_range[0]+shift(i)) & (phis < phi_range[0] + shift(i+1)) &
                    (thetas > theta_range[0]+shift(j)) & (thetas < theta_range[0] + shift(j+1))
                )]
                filtered_data2 = H_coord[np.where(
                    (phis > phi_range[0]+shift(NUM_BINS-i-1)) & (phis < phi_range[0] + shift(NUM_BINS-i)) &
                    (thetas > theta_range[0]+shift(j)) & (thetas < theta_range[0] + shift(j+1))
                )]
                filtered_data = np.concatenate([filtered_data1, filtered_data2])

                num_pts[i,j] = filtered_data.shape[0]
                avgs[i,j] = np.average(filtered_data)

        num_pts = num_pts / np.sum(num_pts) # normalize

        HEADER = "first index corresponds to phi; second index corresponds to theta"
        np.savetxt(probs_path.format(struct,T), num_pts, header=HEADER)
        np.savetxt(nH_path.format(struct, T), avgs, header=HEADER)
