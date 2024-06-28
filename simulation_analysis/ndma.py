import numpy as np

path = "../raw_thermo/{}_H_bond_dict.npy"

structs = [
    ("h", "h", 0, 8, "Hexagonal", "0.1", 144),
    ("p", "p", 0, 9, "Pentagonal", "0.3", 144),
    ("r2", "r2", 0, 19, "Flat-Rhombic", "2.0", 144),
    ("liq_h", "liq_h", 8, 30, "Liquid", "0.1", 144),
    ("r2", "hexatic", 19, 30, "Hexatic", "2.0", 144),
    ("big_r2", "big_r2", 0, 19, "Flat-Rhombic", "2.0", 576),
]

MAX_BOND = 3
NDNA_LEN = (MAX_BOND + 1) * (MAX_BOND + 1)

def convert(num, big):
    if big:
        return int((num-2)/3)
    else:
        return int(num/3)

for (s,save_name,idx0,idx1,name,P,num_Os) in structs:

    if s == "big_r2":
        big = True
    else:
        big = False

    H_bond_dict = np.load(path.format(s), allow_pickle=True).item()

    Ts = [str(20*(i+1)) for i in range(idx0,idx1)]

    struct_ndma = np.zeros(((idx1-idx0),(NDNA_LEN+2)))
    struct_ndma[:,0] = Ts

    for (i_T,T) in enumerate(Ts):

        try:

            num_samples = 0
            temp_ndma = np.zeros(((NDNA_LEN+1),)) # allow room for an "other" column

            _, O_bonds, _ = H_bond_dict[T]

            for frame in O_bonds:

                tracker = np.zeros((num_Os,2)) # first column is donor, second column is acceptor

                for (d,a) in frame:
                    tracker[convert(d,big),0] += 1
                    tracker[convert(a,big),1] += 1

                for i in range(num_Os):

                    d, a = tracker[i,0], tracker[i,1]
                    ndma_idx = int(d*(MAX_BOND+1) + a)
                    num_samples += 1

                    if (ndma_idx < NDNA_LEN):
                        temp_ndma[ndma_idx] += 1.0
                    else:
                        temp_ndma[-1] += 1.0

            temp_ndma = temp_ndma / num_samples
            struct_ndma[i_T,1:] = temp_ndma

        except:

            print("{} structure at {} K failed".format(save_name, T))
            struct_ndma[i_T,1:] = np.inf

    HEADER = "T\t 0D0A\t 0D1A\t 0D2A\t ...\t 1D0A\t 1D1A\t ...\t OTHER"
    np.savetxt("{}_ndma.txt".format(save_name), struct_ndma, header=HEADER)
