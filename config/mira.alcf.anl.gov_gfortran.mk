# arch.mk for BerkeleyGW codes
#
# almost-suitable for Mira at ALCF with gfortran compiler
#
# Run the following command before compiling:
# soft add +mpiwrapper-gnu
#
#########################################
# WARNING: COMPILES BUT FAILS TESTSUITE #
#########################################
#
# Felipe H. da Jornada
# Aug 2013, UCB

COMPFLAG  = -DGNU
PARAFLAG  = -DMPI #-DOMP
MATHFLAG  = -DUSESCALAPACK -DUSEFFTW3 #-DHDF5
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG


# We use the -frecord-marker=4 and -fconvert=little-endian to add
# compatibility with binary WFNs generated on x86 architectures.
FCPP    = /usr/bin/cpp -C
F90free = mpif90 -ffree-form -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fconvert=little-endian -frecord-marker=4
LINK    = mpif90 -fdefault-real-8 -fdefault-double-8 -fconvert=little-endian -frecord-marker=4
FOPTS   = -O0 -g 
FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

# Use the below for debugging (catching segfaults, unused variables, etc)
FOPTS   = -g3 -ggdb -fbounds-check -Wall 
FNOOPTS = $(FOPTS)

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpicxx
C_COMP  = mpicc
C_LINK  = mpicxx
C_OPTS  = -O0 -g
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
BASEPATH     = /soft/libraries/alcf/current/gcc
FFTWPATH     = $(BASEPATH)/FFTW3
FFTWLIB      = -L$(FFTWPATH)/lib -lfftw3
FFTWINCLUDE  = $(FFTWPATH)/include
BLASPATH     = $(BASEPATH)/BLAS
BLASLIB      = -L$(BLASPATH)/lib -lblas
LAPACKPATH   = $(BASEPATH)/LAPACK
LAPACKLIB    = -L$(LAPACKPATH)/lib -llapack $(BLASLIB)
SCALAPACKPATH= $(BASEPATH)/SCALAPACK
SCALAPACKLIB = -L$(SCALAPACKPATH)/lib -lscalapack $(LAPACKLIB)
# to use HDF5, you'll need to `module load hdf5-parallel` and include -DHDF5 in MATHFLAG
#HDF5LIB      = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
#HDF5INCLUDE  = $(HDF5_INCLUDE_OPTS)
PERFORMANCE  = 

TESTSCRIPT = qsub -t 30 -n 1 --mode script mira.scr
