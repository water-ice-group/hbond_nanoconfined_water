import os
import sys
import numpy as np

from ase.calculators.socketio import SocketClient
from ase.io import read
from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators.mixing import LinearCombinationCalculator

from confining_potential_calculator import ConfiningPotentialMorseCalculator

# Define atoms object
print ("Reading atoms object.")
atoms = read("init-ase.xyz", 0)
#atoms = read("debug.xyz", 0)

# Set CP2K calculator #################
calcs = []
n_committee = 0

atom_types = {'O' : 2, 'H' : 1}

print ("Setting up potential.")
for i in range(n_committee):
    cmd_nn = [
            "#!/bin/bash",
            "variable runnerDir       string \"/work/e05/e05/xr223/potentials/revPBE0D3-CS/nnp-" + str(i+1) + "/\"",
            "variable runnerCutoff    equal  12.",
            "pair_style nnp dir ${runnerDir} showew no resetew yes maxew 1000000  cflength 1.889726 cfenergy 0.036749",
            "pair_coeff * * ${runnerCutoff}",
            ]

    calcs.append(LAMMPSlib(atom_types=atom_types, lmpcmds=cmd_nn, tmp_dir="./", log_file='./log.lammps', keep_alive=True))

calcs.append(ConfiningPotentialMorseCalculator(57.8 * 1e-3, 3.85, 0.92, 5))

weights = np.ones(n_committee + 1)
weights[0:n_committee] /=  n_committee 
#weights[0:n_committee] =  0.0 
#weights[-1] =  0.0
#weights = [1]

LinearCombinationCalculator(calcs, weights, atoms)
# Create Client
# inet
port = 11111
host = "r2_nqes_4000_400_cp"
print ("Setting up socket.")
client = SocketClient(unixsocket=host)
print ("Running socket.")
client.run(atoms)
print ("setting up calculator")
