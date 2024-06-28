import numpy as np
import scipy.stats

def convert(num, big):
    if big:
        return int((num-2)/3)
    else:
        return int(num/3)

s = None
save_name = None
idx0, idx1 = 0, 19
name = None
num_Os = 144
big = False

TIMESTEP = 0.5 # how long is one timestep in the numpy file (in femtoseconds)?

path = "r2_wt_H_bond_dict.npy"

H_bond_dict = np.load(path.format(s), allow_pickle=True).item()

Ts = [str(20*(i+1)) for i in range(idx0,idx1)]
Ts = ["220", "280", "320", "380"]
lifetime_Ts = np.zeros((len(Ts),3))
lifetime_Ts[:,0] = Ts

for (i_T,T) in enumerate(Ts):

    _, O_bonds, _ = H_bond_dict[T]

    num_frames = len(O_bonds)
    time_series = np.zeros((num_frames,num_Os))

    for (i_frame,frame) in enumerate(O_bonds):

        tracker = np.zeros((num_Os,2)) # first column is donor, second column is acceptor

        for (d,a) in frame:
            tracker[convert(d,big),0] += 1
            tracker[convert(a,big),1] += 1

        for i in range(num_Os):

            d, a = tracker[i,0], tracker[i,1]

            if (d >= 2) and (a >= 1):
                time_series[i_frame,i] = 1

    lifetimes = []
    flag, lifetime = False, 0
    for i_O in range(num_Os):
        for i_frame in range(num_frames):
            if (flag == True) and (time_series[i_frame,i_O] == 1): # already in state, and *still* in state
                lifetime = lifetime + 1
            elif (flag == True) and (time_series[i_frame,i_O] == 0): # alread in state, but *leaving* state
                flag = False
                lifetimes.append(lifetime)
            elif (time_series[i_frame,i_O] == 1): # wasn't in state before, but *entering* state
                flag = True
                lifetime = 1
        flag = False

    lifetimes = np.array(lifetimes)
    weights = np.ones_like(lifetimes) / len(lifetimes)
    
    lifetime_Ts[i_T,1] = np.average(lifetimes * TIMESTEP)
    lifetime_Ts[i_T,2] = scipy.stats.sem(lifetimes * TIMESTEP)

np.savetxt("lifetime_2D1A.txt", lifetime_Ts)
