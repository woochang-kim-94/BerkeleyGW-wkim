# arch.mk for BerkeleyGW codes
#
# Recommended arch.mk for Babbage at NERSC running natively on the MIC
# Supports OMP threading + HDF5.
#
# Before you compile, do:
# module load hdf5-parallel/mic-1.8.10-p1 fftw/3.3.4-mic
#
# Note: testsuite doesn't pass b/c of a bug in NERSC's file system with babbage.
#
# Felipe H. da Jornada
# 2015, UC Berkeley
#
# Precompiler options
#
COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = /usr/bin/cpp -C
F90free = mpiifort -free -openmp -mmic -static-intel
LINK    = mpiifort -openmp -mmic -static-intel
FOPTS   = -ip -no-ipo -O3 -no-prec-div -fp-model fast=2 -align array64byte
FNOOPTS = -O2 -no-prec-div -fp-model fast=2 -align array64byte 
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpiicpc -openmp -mmic -static-intel
C_COMP  = mpiicc -openmp -mmic -static-intel
C_LINK  = mpiicpc -openmp -mmic -static-intel
C_OPTS   = -ip -no-ipo -O3 -no-prec-div -fp-model fast=2 -static-intel
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
MKLPATH      = $(MKLROOT)/lib/mic
FFTWINCLUDE  = $(FFTW_INC)
FFTWLIB      = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group
LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a

HDF5PATH     =  $(HDF5_DIR)
HDF5LIB      =  $(HDF5PATH)/lib/libhdf5hl_fortran.a \
                $(HDF5PATH)/lib/libhdf5_hl.a \
                $(HDF5PATH)/lib/libhdf5_fortran.a \
                $(HDF5PATH)/lib/libhdf5.a $(ZLIB_DIR)/lib/libz.a -ldl
HDF5INCLUDE  = $(HDF5_DIR)/include

TESTSCRIPT = sbatch babbage.scr
