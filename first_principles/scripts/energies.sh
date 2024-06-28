#!/bin/bash
green='\033[0;32m'
nocolor='\033[0m'

#Hartree to mEv
#conversion_factor=2625.4996
conversion_factor=627.5096080305927
total_wat_molec=144

xyz_chain='chain_r2/r2_chain.xyz'
num_wat_chain=$(head -n 1 "$xyz_chain")
wat_molec_per_chain=$(echo "scale=6; $num_wat_chain / 3" | bc)



#What is the cohesive / binding energy of a water molecule in a chain?
echo -e "${green}"
echo "# TOTAL COHESIVE ENERGIES"
filename_chain='chain_r2/r2_chain.out'
chain_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_chain")
filename_wat='wat_r2/r2_wat.out'
wat_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_wat")
total_energy_mev=$(echo "scale=6; $chain_energy * $conversion_factor / $wat_molec_per_chain  - $wat_energy * $conversion_factor  " | bc)
total_energy_mev_wat=$(echo "scale=6; $total_energy_mev " | bc)
echo "cohesive energy of a water molecule in a chain (in kcal/mol): $total_energy_mev_wat"

#What is the cohesive energy of a water molecule in r2?
filename_total='total_r2/r2_total.out'
total_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_total")
filename_wat='wat_r2/r2_wat.out'
wat_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_wat")
total_energy_mev=$(echo "scale=6; $total_energy * $conversion_factor / $total_wat_molec  - $wat_energy * $conversion_factor  " | bc)
echo "cohesive energy of a water molecule in r2 (in kcal/mol): $total_energy_mev"

#What is the cohesive energy of the chains in r2?
filename_total='total_r2/r2_total.out'
total_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_total")
filename_chain='chain_r2/r2_chain.out'
chain_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_chain")
total_energy_mev=$(echo "scale=6; $total_energy * $conversion_factor / $total_wat_molec  - $chain_energy * $conversion_factor / $wat_molec_per_chain " | bc)
#total_energy_mev=$(echo "scale=6; $total_energy_mev / $wat_molec_per_chain " | bc)
echo "cohesive energy of the chains in r2 (in kcal/mol): $total_energy_mev"

#What is the role on dispersion in stabilizing a water molecule in a chain?
echo ""
echo "# COHESIVE ENERGIES WITHOUT DISPERSION"
filename_chain='chain_r2/r2_chain.out'
chain_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_chain")
dispersion_energy=$(awk '/Dispersion energy:/ {print $3; exit}' "$filename_chain")
chain_energy=$(echo "scale=6; $chain_energy - $dispersion_energy  " | bc)
filename_wat='wat_r2/r2_wat.out'
wat_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_wat")
dispersion_energy=$(awk '/Dispersion energy:/ {print $3; exit}' "$filename_wat")
wat_energy=$(echo "scale=6; $wat_energy - $dispersion_energy  " | bc)
total_energy_mev=$(echo "scale=6; $chain_energy * $conversion_factor / $wat_molec_per_chain  - $wat_energy * $conversion_factor  " | bc)
total_energy_mev_wat=$(echo "scale=6; $total_energy_mev " | bc)
echo "cohesive energy of a water molecule in a chain (in kcal/mol): $total_energy_mev_wat"

#What is the role on dispersion in stabilizing a water molecule in r2?
filename_total='total_r2/r2_total.out'
total_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_total")
dispersion_energy=$(awk '/Dispersion energy:/ {print $3; exit}' "$filename_total")
total_energy=$(echo "scale=6; $total_energy - $dispersion_energy  " | bc)
filename_wat='wat_r2/r2_wat.out'
wat_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_wat")
dispersion_energy=$(awk '/Dispersion energy:/ {print $3; exit}' "$filename_wat")
wat_energy=$(echo "scale=6; $wat_energy - $dispersion_energy  " | bc)
total_energy_mev=$(echo "scale=6; $total_energy * $conversion_factor / $total_wat_molec  - $wat_energy * $conversion_factor  " | bc)
echo "cohesive energy of a water molecule in r2 (in kcal/mol): $total_energy_mev"

#What is the role on dispersion in stabilizing a chain in r2?
filename_total='total_r2/r2_total.out'
total_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_total")
dispersion_energy=$(awk '/Dispersion energy:/ {print $3; exit}' "$filename_total")
total_energy=$(echo "scale=6; $total_energy - $dispersion_energy  " | bc)
filename_chain='chain_r2/r2_chain.out'
chain_energy=$(awk '/Total energy:/ {print $3; exit}' "$filename_chain")
dispersion_energy=$(awk '/Dispersion energy:/ {print $3; exit}' "$filename_chain")
chain_energy=$(echo "scale=6; $chain_energy - $dispersion_energy  " | bc)
total_energy_mev=$(echo "scale=6; $total_energy * $conversion_factor / $total_wat_molec  - $chain_energy * $conversion_factor / $wat_molec_per_chain " | bc)
#total_energy_mev=$(echo "scale=6; $total_energy_mev / $wat_molec_per_chain " | bc)
echo "cohesive energy of the chains in r2 (in kcal/mol): $total_energy_mev"


echo -e "${nocolor}"
