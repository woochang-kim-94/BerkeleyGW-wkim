# arch.mk for BerkeleyGW codes
#
# Configuration file for a buildbot worker running on Amazon EC2.
# Ubuntu xenial, nag compiler + LAPACK + ScaLAPACK + openmpi + HDF5, not threaded
#
# All tests passed, r7296
#
# Felipe H da Jornada
# Jul 2018

# FHJ: We disable OMP threading because it is too slow with hyperthreads @ EC2.
COMPFLAG  = -DNAG
PARAFLAG  = -DMPI #-DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DUSEMR3 #-DHDF5 -DUSEELPA -DUSEPRIMME
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG -DVERBOSE

#########################################################################
#   NOTE: This arch.mk is used by the testsuite. The compiler flags are #
#   optimized to debug the code and not for code performance.           #
#########################################################################
FCPP    = cpp -C -nostdinc
# FHJ: Note that we disabled threading.
F90free = mpif90 -colour -free -f2008 -ieee=full -C -mismatch_all -I $(PREFIX)/BSE/full_solver -dcfuns -u #-openmp
LINK    = mpif90 -colour -f2008 -ieee=full -C -Wl,-fopenmp #-openmp
#F90free = mpif90 -colour -free -f2008 -ieee=full -C -openmp -mismatch_all -I $(PREFIX)/BSE/full_solver -dcfuns -u
#LINK    = mpif90 -colour -f2008 -ieee=full -C -openmp -Wl,-fopenmp
#F90free = mpif90 -colour -free -f2008 -ieee=full -C -gline -mismatch_all -I $(PREFIX)/BSE/full_solver -dcfuns #-u
#LINK    = mpif90 -colour -f2008 -ieee=full -C -gline -Wl,-fopenmp
# -fast, and even -O2, give me an out-of-memory error
FOPTS   = -O3
FNOOPTS = $(FOPTS)
MOD_OPT = -mdir 
INCFLAG = -I

C_PARAFLAG  = -DPARA
CC_COMP  = mpiCC -g -fopenmp
C_COMP  = mpicc -g -fopenmp
C_LINK  = mpiCC -g -fopenmp
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = /opt/nag/lib/libfftw3_omp.a /opt/nag/lib/libfftw3.a
FFTWINCLUDE  = /opt/nag/include

LDIR         = /opt/nag/lib
LAPACKLIB    = $(LDIR)/liblapack.a /opt/nag/OpenBLAS-0.2.20/lib/libopenblas.a
SCALAPACKLIB = $(LDIR)/libscalapack.a

#HDF5_LDIR    =  /opt/nag/lib
#HDF5LIB      =  $(HDF5_LDIR)/libhdf5hl_fortran.a \
#                $(HDF5_LDIR)/libhdf5_hl.a \
#                $(HDF5_LDIR)/libhdf5_fortran.a \
#                $(HDF5_LDIR)/libhdf5.a -ldl
#HDF5INCLUDE  = /opt/nag/include

TESTSCRIPT = OMP_NUM_THREADS=1 BGW_CHECK_FLAGS='-g binary' bash singularity.scr
