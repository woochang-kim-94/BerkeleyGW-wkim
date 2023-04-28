#!/bin/bash -l

# Script to automatically generate complex version of Si-EPM test file
# Input is Si.test, output is Si_cplx.test. All matches should be identical
# between the two versions, ensured by dynamic generation by this script.
# David Strubbe, 16 Oct 2010

OUTPUT=`echo $1 | sed 's/\./_cplx\./'`
sed 's/(real)/(complex)/;s/.real.x/.cplx.x/' $1 > $OUTPUT
