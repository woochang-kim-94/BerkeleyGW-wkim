# arch.mk for BerkeleyGW codes
#
# Suitable for generic MacOS with multiple processors (e.g. any modern Mac)
# Use macports to install mpich-default +gfortran, fftw-3 +gfortran, veclibfort,
# scalapack +mpich
# Add location of the binaries (/opt/local/bin by default) to PATH
#
# With OpenMPI rather than MPICH, ScaLAPACK in diagonalization for
# BSE with complex wavefunctions may hang.
#
# D. Strubbe, 2021

COMPFLAG  = -DGNU
PARAFLAG  = -DMPI
MATHFLAG  = -DUSEFFTW3 -DUSESCALAPACK
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
# Some error will occur compiling with DEBUG in current version.
#DEBUGFLAG = -DDEBUG

VERSION = -mp-11
MPI     = mpich
FCPP    = cpp$(VERSION) -C
F90free = $(LINK) -ffree-form -ffree-line-length-none -fno-second-underscore
LINK    = mpif90-$(MPI)-mp
#FOPTS   = -g -C -O0 -fbounds-check -fbacktrace -fallow-argument-mismatch
FOPTS   = -O3 -fallow-argument-mismatch
FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

# NOTE: this is not gcc, but clang.
CC_COMP = mpicxx-$(MPI)-mp
C_COMP  = mpicc-$(MPI)-mp
C_LINK  = $(CC_COMP)
C_OPTS  = -O3
C_PARAFLAG  = -DPARA
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/opt/local/lib/ -lfftw3
#FFTWLIB      = -L/opt/local/lib/ -ldfftw  # for FFTW2
FFTWINCLUDE  = /opt/local/include
LAPACKLIB    = -L/opt/local/lib/ -lveclibfort
SCALAPACKLIB = -L/opt/local/lib/ -lscalapack

#To build with HDF5, install port hdf5 -cxx +hl +fortran +gfortran +mpich
# and set:
#MATHFLAG += -DHDF5
#HDF5INCLUDE  = /opt/local/include
#HDF5LIB      = -L/opt/local/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz

TESTSCRIPT = MPIEXEC=`which mpiexec-$(MPI)-mp` make check-parallel
