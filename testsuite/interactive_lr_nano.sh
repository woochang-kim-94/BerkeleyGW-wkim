#!/bin/bash -l

# To run interactively:
# qsub -I lawrencium.scr [or nano.scr]
# cd $PBS_O_WORKDIR
# sh interactive_lr_nano.sh

export TEMPDIRPATH=$SCRATCH/tmp

export MPIEXEC="`which mpiexec` --mca btl openib,self"
make check-parallel | tee test.out

echo "Exit status = $?"
