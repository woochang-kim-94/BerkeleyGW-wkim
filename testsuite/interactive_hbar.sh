#!/bin/bash -l

# To run interactively:
# qsub -I hbar.scr
# cd $PBS_O_WORKDIR
# sh interactive_hbar.sh

export TEMPDIRPATH=$SCRATCH/tmp

export MPIEXEC="`which mpiexec` --mca btl gm,self"
make check-parallel | tee test.out

echo "Exit status = $?"
