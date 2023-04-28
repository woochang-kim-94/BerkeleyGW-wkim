# arch.mk for BerkeleyGW codes
#
# suitable for cruller.mit.edu, with Intel compilers and MVAPICH
# Status: works only for 4 procs
#
# Run the following command before compiling:
# module load mvapich2-2.0/intel
#
# D. Strubbe
# Oct 2012, MIT

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI
MATHFLAG  = -DUSESCALAPACK -DUSEFFTW3
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C
F90free = mpif90 -free
LINK    = mpif90 -i-static -openmp
FOPTS   = -O3 -no-prec-div
FNOOPTS = -O2 -no-prec-div
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG = -DPARA
CC_COMP = mpicxx
C_COMP  = mpicc
C_LINK  = mpicxx
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
# Using MKL for FFTW3 gives wrong answers!!!
FFTWPATH     = /opt/fftw3/mvapich2-2.0-intel
FFTWLIB      = $(FFTWPATH)/lib/libfftw3.a
FFTWINCLUDE  = $(FFTWPATH)/include
MKLPATH      = /opt/intel/cmkl/10.2.5.035/lib/em64t
LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a

TESTSCRIPT = qsub cruller.scr
