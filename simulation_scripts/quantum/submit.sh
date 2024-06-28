#!/bin/bash
#SBATCH --job-name=run-disordered-b32-t400-p4000
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=e05-surfin-mic
#SBATCH --partition=standard
#SBATCH --qos=standard

module load PrgEnv-cray
module load cray-python

# Avoid any unintentional OpenMP threading by setting OMP_NUM_THREADS
export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/work/e05/e05/xr223/.local
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.8/site-packages:$PYTHONPATH
source /work/e05/e05/xr223/source/n2p2/env.sh

ipi=/work/e05/e05/xr223/source/i-pi/bin/i-pi
lmp=/work/e05/e05/xr223/source/n2p2/bin/lmp_mpi

rm /tmp/ipi_*

if [ -f "simulation.restart" ]; then
    python3 ${ipi} simulation.restart > log.i-pi &
else
    python3 ${ipi} input.xml &
fi

sleep 30

for x in {1..32}
do
    ${lmp} < in1.lmp > log.lmp &
    python3 run-ase.py > log.cp &
done
        
wait
