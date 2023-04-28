# arch.mk for BerkeleyGW codes
#
# Suitable for generic single-processor MacOS 
# Use macports to install gcc11 (or latest gcc), fftw-3 +gfortran, veclibfort
# Add location of the binaries (/opt/local/bin by default) to PATH
#
# D. Strubbe, 2021

COMPFLAG  = -DGNU
MATHFLAG  = -DUSEFFTW3
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

VERSION = -mp-11
FCPP    = cpp$(VERSION) -C
F90free = gfortran$(VERSION) -ffree-form -ffree-line-length-none -fno-second-underscore
LINK    = gfortran$(VERSION)
FOPTS   = -O3
#FOPTS   = -g -O0 -fbounds-check -fbacktrace
FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

CC_COMP = g++$(VERSION)
C_COMP  = gcc$(VERSION)
C_LINK  = g++$(VERSION)
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/opt/local/lib/ -lfftw3
#FFTWLIB      = -L/opt/local/lib/ -ldfftw  # for FFTW2
FFTWINCLUDE  = /opt/local/include
LAPACKLIB    = -L/opt/local/lib/ -lveclibfort

#To build with HDF5, install port hdf5 +hl +fortran +gfortran +mpich and set:
#MATHFLAG += -DHDF5
#HDF5INCLUDE  = /opt/local/include
#HDF5LIB      = -L/opt/local/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
