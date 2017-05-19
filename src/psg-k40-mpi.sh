#!/bin/bash
##SBATCH --job-name="hoomd-k40-mpi"
##SBATCH --partition=ivb_k40
##SBATCH --ntasks-per-node=6
##SBATCH --nodes=4
##SBATCH -t 02:00:00

LIMIT="################################################################"
function new() {
    for i in `seq 10`; do
        echo ""
    done
}

function rexp () {
    echo $LIMIT
    echo Running: mpirun -n $1 python polymer_brush.py --mode=$2 --user=$3
    mpirun -n $1 python polymer_brush.py --mode=$2 --user=$3
    new;
}


echo $LIMIT
echo GPU n=8
echo $LIMIT

rexp 1 gpu 8;
rexp 2 gpu 8;
rexp 4 gpu 8;
rexp 8 gpu 8;
rexp 16 gpu 8;

echo $LIMIT
echo GPU n=16
echo $LIMIT
rexp 1 gpu 16;
rexp 2 gpu 16;
rexp 4 gpu 16;
rexp 8 gpu 16;
rexp 16 gpu 16;

echo $LIMIT
echo GPU n=32
echo $LIMIT
rexp 1 gpu 32;
rexp 2 gpu 32;
rexp 4 gpu 32;
rexp 8 gpu 32;
rexp 16 gpu 32;


echo $LIMIT
echo CPU n=8
echo $LIMIT
rexp 16 cpu 8;
rexp 8 cpu 8;
rexp 4 cpu 8;
rexp 2 cpu 8;
rexp 1 cpu 8;

echo $LIMIT
echo CPU n=16
echo $LIMIT
rexp 16 cpu 16;
rexp 8 cpu 16;
rexp 4 cpu 16;
rexp 2 cpu 16;
rexp 1 cpu 16;


echo $LIMIT
echo CPU n=32
echo $LIMIT
rexp 16 cpu 32;
rexp 8 cpu 32;
rexp 4 cpu 32;
rexp 2 cpu 32;
rexp 1 cpu 32;
