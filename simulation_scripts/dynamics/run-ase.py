import os
import sys
import numpy as np

from ase.io import read
from ase.calculators.socketio import SocketClient
from confining_potential_calculator import ConfiningPotentialMorseCalculator

# Define atoms object
atoms = read("init-ase.xyz", 0)

# Set calculator #################
width = 5.0
D0 = np.array([ 90.1069 , 92.3613 , 98.9979 ])
x0 = np.array([ 3.10529 , 3.45671 , 3.37265 ])
a  = np.array([ 1.30431 , 1.43449 , 1.34725 ])
calc = ConfiningPotentialMorseCalculator(D0 * 1e-3, x0, a, width)
#calc = ConfiningPotentialMorseCalculator(57.8 * 1e-3, 3.85, 0.92, width)
atoms.set_calculator(calc)

# Create Client
host = "wt_r2_2000_220_stat_cl_05.0_cp"
client = SocketClient(unixsocket=host)
client.run(atoms)
