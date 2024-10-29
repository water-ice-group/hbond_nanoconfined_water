# Simulation Scripts

This directory contains examples of the i-PI input files for running molecular dynamics simulations of nanoconfined water.

The initial configurations for each phase are in `../structs`.

The three sub-directories here correspond to:

- `MLP`: the files for the committee of 8 machine learning interatomic potentials used in this work
  NOTE: at simulation time, only the first model in `MLP/train_001` is used to generate forces
- `classical`: classical molecular dynamics of nanoconfined water with a thermostat
- `dynamical`: classical molecular dynamics of nanoconfined water without a thermostat (for dynamics)
- `quantum`: path integral molecular dynamics of nanoconfined water that account for nuclear quantum effects

