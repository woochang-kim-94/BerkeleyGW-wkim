# arch.mk for BerkeleyGW codes
#
# Configuration file for a buildbot worker running on Amazon EC2.
# Centos7, pgi compiler + LAPACK + ScaLAPACK + openmpi + HDF5, threaded
#
# All tests passed, r7296
#
# Felipe H da Jornada
# Jul 2018

COMPFLAG  = -DPGI
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DUSEMR3 #-DUSEELPA -DUSEPRIMME
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

#########################################################################
#   NOTE: This arch.mk is used by the testsuite. The compiler flags are #
#   optimized to debug the code and not for code performance.           #
#########################################################################
FCPP    = cpp -C -nostdinc
#F90free = mpif90 -Mfree -mp -C -Mchkptr -Mchkstk -traceback -Mallocatable=03
F90free = mpif90 -Mfree -mp -Mallocatable=03
LINK    = mpif90 -mp
# -fast, and even -O2, give me an out-of-memory error
FOPTS   = -O0
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA
#CC_COMP  = mpiCC -mp #-traceback
#C_COMP  = mpicc -mp #-traceback
#C_LINK  = mpiCC -mp #-traceback
CC_COMP  = mpiCC -mp
C_COMP  = mpicc -mp
C_LINK  = mpiCC -mp
C_OPTS  = -O0
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = /opt/pgi/fftw-3.3.5/lib/libfftw3_omp.a /opt/pgi/fftw-3.3.5/lib/libfftw3.a
FFTWINCLUDE  = /opt/pgi/fftw-3.3.5/include

LDIR         = /opt/pgi/linux86-64/18.4/lib
LAPACKLIB    = $(LDIR)/liblapack.a $(LDIR)/libblas.a
SCALAPACKLIB = $(LDIR)/scalapack/scalapack-2.0.2/openmpi-2.1.2/lib/libscalapack.a

HDF5_LDIR    =  /opt/pgi/hdf5-1.10.1/lib
HDF5LIB      =  $(HDF5_LDIR)/libhdf5hl_fortran.a \
                $(HDF5_LDIR)/libhdf5_hl.a \
                $(HDF5_LDIR)/libhdf5_fortran.a \
                $(HDF5_LDIR)/libhdf5.a -lz -ldl
HDF5INCLUDE  = /opt/pgi/hdf5-1.10.1/include

TESTSCRIPT = bash singularity.scr
