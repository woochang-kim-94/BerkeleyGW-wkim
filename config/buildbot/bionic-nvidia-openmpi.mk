# arch.mk for BerkeleyGW codes
#
# Configuration file for a buildbot worker running on Amazon EC2.
# Ubuntu 18.04 "Bionic", Nvidia HPC SDK 20.11 + FFTW3 + HDF5, threaded
#
# All tests passed for commit 54fc2545
#
# William Huhn
# Feb 2021

COMPFLAG  = -DPGI
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DUSEMR3 #-DOPENACC -DOMP_TARGET -DUSEELPA -DUSEPRIMME
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

#########################################################################
#   NOTE: This arch.mk is used by the testsuite. The compiler flags are #
#   optimized to debug the code and not for code performance.           #
#########################################################################
FCPP    = cpp -C -nostdinc
#F90free = mpif90 -Mfree -mp -C -Mchkptr -Mchkstk -traceback
F90free = mpif90 -Mfree -mp=multicore
LINK    = mpif90 -mp=multicore
# -O2 and above causes the Benzene-SAPO and Si-EPM tests to fail
FOPTS   = -O1
FNOOPTS = -O0
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA
CC_COMP  = mpicxx -mp
C_COMP  = mpicc -mp
C_LINK  = mpicxx -mp
C_OPTS  = -O0
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWROOT     = /opt/fftw/
FFTWLIB      = $(FFTWROOT)/lib/libfftw3_omp.a \
               $(FFTWROOT)/lib/libfftw3.a
FFTWINCLUDE  = $(FFTWROOT)/include

LAPACKLIB    = -llapack -lblas
SCALAPACKLIB = -lscalapack

HDF5ROOT     = /opt/hdf5
HDF5LIB      = $(HDF5ROOT)/lib/libhdf5hl_fortran.a \
               $(HDF5ROOT)/lib/libhdf5_hl.a \
               $(HDF5ROOT)/lib/libhdf5_fortran.a \
               $(HDF5ROOT)/lib/libhdf5.a -ldl
HDF5INCLUDE  = $(HDF5ROOT)/include

TESTSCRIPT = bash singularity.scr
