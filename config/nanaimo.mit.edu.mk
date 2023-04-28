# arch.mk for BerkeleyGW codes
#
# suitable for nanaimo.mit.edu
#
# Run the following command before compiling:
# module load mvapich2-2.1/intel
#
# Status: all tests passed on 16 procs, r6941
#
# D. Strubbe
# Dec 2015, MIT

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI
MATHFLAG  = -DUSESCALAPACK -DUSEFFTW3
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C
F90free = mpif90 -free
LINK    = mpif90 -openmp
FOPTS   = -O3 -no-prec-div
FNOOPTS = -O2 -no-prec-div
#FOPTS   = -g -O0 -traceback -check all -Warn all
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
# Using MKL for FFTW3
FFTWINCLUDE  = /opt/intel/mkl/include/fftw
MKLPATH      = /opt/intel/mkl/lib/intel64
LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a

TESTSCRIPT = qsub nanaimo.scr
