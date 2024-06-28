# Simulation Analysis

The `hbond.py` file contains two functions that process the `pdb` files produced by the i-PI simulations in `../simulation_scripts`. The `get_OH_dict.py` file shows how to use these functions to produce the `npy` files that will be used by the other analysis scripts. These `npy` files are all that are needed to run the below analysis scripts.

Here is a list of the analysis scripts provided, along with a description of what they do:
- `process_2D.py`: Computes the theta-phi 2D histograms. This includes the log-probability histogram and the number of hydrogen bonds histogram.
- `get_nH.py`: Calculates the statistics for the number of putative/geometric hydrogen bonds in a molecular dynamics trajectory.
- `ndma.py`: Computes the NDMA hydrogen bonding probability distributions (the NDMA notation is defined in the main text).
- `lifetime_2D1A.py`: Computes the statistics for the lifetime of the 2D1A state.
- `sigma.py`: Calculates the trajectory of $\sigma$ order parameter values as a function of simulation time, as well as the free energy profile along $\sigma$.
