#!/bin/bash -l
# gen_doxygen.sh
#
# Create an input file for doxygen (Doxyfile). We have to use this script
# (instead of simple wildcards) because doxygen gets confused with Fortran
# preprocessing macros.
#
# Felipe H. da Jornada, Berkeley, October 2011

path='../../../Common ../../../Epsilon ../../../Sigma ../../../BSE ../../../PlotXct ../../../MeanField ../../../NonLinearOptics ../../../Visual ../../../GW+C'
exts='*.f90 *.F90 *.cpp *.c *.cc'
cond='-false'
for ext in $exts; do
	cond="$cond -or ( -name $ext ) "
done
cond="( ( $cond ) -and ( -not -name *_inc.f90 ) )"
input=`find $path $cond | tr "\n" ' '`

base=${PWD%/doxygen}

cat Doxyfile.template | sed -e "s!%BASE%!$base!" -e "s!%INPUT%!$input!" -e "s!%INCLUDE_PATH%!$path!" > Doxyfile
