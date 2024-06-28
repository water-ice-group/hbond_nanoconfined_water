import numpy as np
from sklearn.cluster import KMeans

# Input paths
OH_path = "../raw_data/{}_OH_dict.npy"

# Output path
traj_path = "{}_{}K_sigma_traj.txt"
fep_path = "{}_{}K_sigma_FEP.txt"

NUM_BINS=16
OH_cutoff = 1.2
sigma_range = (-np.pi/4.0, np.pi/4.0)

def shift(idx):
    return (np.pi*idx) / NUM_BINS

"""
The below is written to be flexible wrt the # of O's in anticipation of larger unit cells.
Format for each element in the `structs` array is (`struct_name`, `num_Os`, `Ts[3]`).
"""

structs = [
    ("r2", 144, 12, [100,280,380]),
    ("big_r2", 576, 24, [100,280,380]),
]

def convert(num, big):
    if big:
        return int((num-2)/3)
    else:
        return int(num/3)

for (struct, num_Os, num_rows, Ts) in structs:

    if struct == "big_r2":
        big = True
    else:
        big = False

    OH_dict = np.load(OH_path.format(struct), allow_pickle=True).item()

    for T in Ts:
        
        """
        This code uses the same `phi` angle calculation as the rest of the codebase. Then we use a
        somewhat hacky way to get the row identities, but we check pretty carefully that the hacky
        way succeeds. If you're getting an exception here, it's worth looking into the details of
        this part to check what's going on.
        """

        Os, Hs = OH_dict[str(T)]
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
        This part is supposed to identify the different rows of water molecules... this is the
        hacky part that you have to be careful about. At the end, we just easily compute sigma
        values from the resulting assignments.
        """
        
        y_coords = Os[:,:,1]
        y_coords_shape = y_coords.shape
        y_coords = y_coords.reshape((-1,1))
        
        kmeans = KMeans(n_clusters=num_rows)
        row_labels = kmeans.fit(y_coords)
        row_labels = kmeans.predict(y_coords)
        row_labels = np.reshape(row_labels, y_coords_shape)
        y_coords = np.reshape(y_coords, y_coords_shape)
        
        centroids = (kmeans.cluster_centers_).flatten()
        ordering = np.argsort(centroids)
        
        row_nums = np.ones(y_coords_shape) * np.inf
        
        for r in range(y_coords_shape[0]):
            for c in range(y_coords_shape[1]):
                row_nums[r,c] = np.where(ordering == row_labels[r,c])[0]
        
        for i in range(num_rows):
            row_counts = np.count_nonzero(row_nums == i, axis=1)
            expected_count = int(num_Os / num_rows)
            if not (np.count_nonzero(row_counts == expected_count) == y_coords.shape[0]):
                raise Exception("Missing water molecules in row {}.".format(i))
        
        switch = (row_nums % 2) * 2.0 - 1.0
        sigmas = np.average(np.multiply(phis, switch), axis=1)
        
        symmetrized_sigmas = np.concatenate([sigmas, -1.0*sigmas], axis=None)
        hist, bin_edges = np.histogram(symmetrized_sigmas, bins=NUM_BINS, range=sigma_range)
        fep = -np.log(hist)
        fep = fep - np.min(fep)
        bin_edges = bin_edges[:-1] + (bin_edges[1] - bin_edges[0]) / 2.0
        
        # save sigmas trajectory
        
        HEADER = "trajectory of sigma values"
        np.savetxt(traj_path.format(struct,T), sigmas, header=HEADER)
        
        # save fep
        
        if not (np.shape(fep) == np.shape(bin_edges)):
            raise Exception("Something very strange has gone wrong with the numpy histogramming...")
        
        fep_save_data = np.zeros((np.shape(fep)[0],2))
        fep_save_data[:,0] = bin_edges
        fep_save_data[:,1] = fep
        HEADER = "(centered) sigma ; free energy"
        np.savetxt(fep_path.format(struct,T), fep_save_data, header=HEADER)
