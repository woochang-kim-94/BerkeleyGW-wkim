#!/bin/bash -l

export OMP_NUM_THREADS=2
export BGW_NUM_PROCS=2
#./run_testsuite.sh -b ../../bin 2>&1 | tee test.out
./run_testsuite.sh -b ../../bin
