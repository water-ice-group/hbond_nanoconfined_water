#!/bin/bash
#SBATCH -J opt
#SBATCH --time=24:00:00
#SBATCH -n1 -N1
#SBATCH --job-name=r2_2000_20_stat_cl_05.0

module load gcc mpi/openmpi/gnu7 mkl

source ~/.bashrc
conda activate py3.8

ipi=~/source/i-pi-dev/bin/i-pi
lmp=~/source/lammps/src/lmp_mpi

rm /tmp/ipi_r2_2000_20_stat_cl_05.0*

if [ -f "RESTART" ]; then
${ipi} RESTART > log.i-pi &
else
${ipi} input.xml > log.i-pi &
fi

sleep 30

${lmp} < in1.lmp > log.lmp &
python run-ase.py > log.cp &

wait

