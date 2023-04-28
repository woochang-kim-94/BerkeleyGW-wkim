# arch.mk for BerkeleyGW codes
#
# almost-suitable for Mira at ALCF with ibm compiler
# In order to use this arch.mk, edit ~/.soft and make +mpiwrapper-xl your
# default compiler. Then, type "resoft".
#
#########################################
# WARNING: COMPILES BUT FAILS TESTSUITE #
#########################################
#
# Jack Deslippe
# Sep 2013, NERSC

COMPFLAG  = -DXLF
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSEFFTW3 -DUSESCALAPACK -DHDF5
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = /usr/bin/cpp -C -P
F90free = mpixlf90_r -qsmp=omp -qthreaded
LINK    = mpixlf90_r -qsmp=omp -qthreaded
FOPTS   = -qp -qalias=noaryovrlp:nointptr -O2
FNOOPTS = $(FOPTS)
MOD_OPT = -qmoddir=
INCFLAG = -I

# Use the below for debugging (catching segfaults, unused variables, etc)
#FOPTS   = -O9 -g9 -qcheck=all
#FNOOPTS = $(FOPTS)

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpixlcxx_r -qsmp=omp -qthreaded
C_COMP  = mpixlc_r -qsmp=omp -qthreaded
C_LINK  = mpixlcxx_r -qsmp=omp -qthreaded
C_OPTS  = -O2
C_DEBUGFLAG = -O0 -g9

REMOVE  = /bin/rm -f

# Math Libraries
#
BASEPATH     = /soft/libraries/alcf/current/xl
FFTWPATH     = $(BASEPATH)/FFTW3
FFTWLIB      = -L$(FFTWPATH)/lib -lfftw3_omp -lfftw3
FFTWINCLUDE  = $(FFTWPATH)/include

#Using ESSL for BLAS
BLASPATH     = $(BASEPATH)/BLAS
BLASLIB      = -L$(BLASPATH)/lib -lblas
LAPACKPATH   = $(BASEPATH)/LAPACK
LAPACKLIB    = -L$(LAPACKPATH)/lib -llapack $(BLASLIB)
SCALAPACKPATH= $(BASEPATH)/SCALAPACK
SCALAPACKLIB = -L$(SCALAPACKPATH)/lib -lscalapack $(LAPACKLIB)


# to use HDF5, you'll need to `module load hdf5-parallel` and include -DHDF5 in MATHFLAG
HDF5LIB      = -L/soft/libraries/hdf5/1.8.10/cnk-xl/current/lib/ -lhdf5hl_fortran -lhdf5_hl \
               -lhdf5_fortran -lhdf5 -L/soft/libraries/alcf/current/xl/ZLIB/lib/ -lz
HDF5INCLUDE  = /soft/libraries/hdf5/1.8.10/cnk-xl/current/include/
PERFORMANCE  = 

TESTSCRIPT = qsub -t 30 -n 1 --mode script mira.scr
