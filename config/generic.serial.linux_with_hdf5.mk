# arch.mk for BerkeleyGW codes
#
# Suitable for Generic single processor linux 
# Tested on laptop running Ubuntu 12.04.  You will need to have yum 
# or apt-get installed liblapack-dev, fftw-dev, gfortran, g++, 
# libhdf5-serial-dev packages
#
# J. Deslippe
# Jun, 2013, NERSC

COMPFLAG  = -DGNU
MATHFLAG  = -DHDF5
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

#########################################################################
#   NOTE: This arch.mk is used by a buildslave. The compiler flags are  #
#   optimized to debug the code and not for code performance.           #
#########################################################################
FCPP    = cpp -C
#F90free = gfortran -ffree-form -ffree-line-length-none
F90free = gfortran -ffree-form -ffree-line-length-none -fbounds-check -Wall -pedantic-errors -std=gnu
LINK    = gfortran
FOPTS   = -O3
FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

CC_COMP = g++ -Wall -pedantic-errors -std=c++0x
C_COMP  = gcc -Wall -pedantic-errors -std=c99
C_LINK  = g++
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/lib/ -lfftw
FFTWINCLUDE  = /usr/share/doc/fftw-dev/fortran
LAPACKLIB    = -L/usr/lib/ -llapack -lblas
HDF5LIB      = -L/usr/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
HDF5INCLUDE  = /usr/include
