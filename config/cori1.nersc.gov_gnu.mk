# arch.mk for BerkeleyGW codes using gnu compiler
# Features: MPI, LibSci, HDF5, dynamic linking. ** No threading! **
#
# Suitable for Cori Phase 1 at NERSC
#
# *** THIS IS *NOT* THE RECOMMENDED ARCH.MK FOR CORI PHASE 1 ***
# *** We recommend you use cori1.nersc.gov_intel.mk instead. ***
#
# Load the following modules to compile **AND TO RUN** BerkeleyGW
# module swap PrgEnv-intel PrgEnv-gnu
# module load cray-fftw cray-hdf5-parallel
#
# Felipe H. da Jornada (2015)
#
# This version *almost* passes the testsuite, r6922. The graphene absp test fails
# because of the degeneracy of the coarse k-points, and the Si/haydock test
# fails because of the optimization flags. This is a known issue with the
# haydock code. Threading does not currently work.

# Precompiler options
#
COMPFLAG  = -DGNU
PARAFLAG  = -DMPI #-DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = /usr/bin/cpp -ansi 
F90free = ftn -ffree-form -ffree-line-length-none -fno-second-underscore
# If you don`t use the -dynamic flag, your binary files will take a total of ~3GB!
LINK    = ftn -dynamic
FOPTS   = -O3 -funroll-loops -funsafe-math-optimizations
FNOOPTS = $(FOPTS)
MOD_OPT = -J 
INCFLAG = -I

# Debugging options:
#F90free = ftn -ffree-form -ffree-line-length-none -fno-second-underscore -fbounds-check -fbacktrace -Wall #-finit-real=nan -finit-integer=-100 -ffpe-trap=invalid
#LINK    = ftn -fbounds-check -fbacktrace -Wall -dynamic #-ffpe-trap=invalid
#FOPTS   = -g -g3 -C
#FNOOPTS = $(FOPTS)

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = CC
C_COMP  = cc
C_LINK  = CC -dynamic
C_OPTS  = -O3 -ffast-math
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
# FHJ: LibSci already provides FFT routines. But are they threaded?!
# FFTWLIB      = -L$(FFTW_DIR)/lib  -lfftw3  # -lfftw3_omp
FFTWINCLUDE  = $(FFTW_INC)
PERFORMANCE  = 

LAPACKLIB = 

HDF5_LDIR    =  $(HDF5_DIR)/lib
HDF5LIB      =  $(HDF5_LDIR)/libhdf5hl_fortran.a \
                $(HDF5_LDIR)/libhdf5_hl.a \
                $(HDF5_LDIR)/libhdf5_fortran.a \
                $(HDF5_LDIR)/libhdf5.a -lz -ldl
HDF5INCLUDE  = $(HDF5_DIR)/include

TESTSCRIPT = sbatch cori1.scr
