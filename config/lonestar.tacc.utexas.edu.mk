# arch.mk for BerkeleyGW codes
#
# suitable for Lonestar at TACC
#
# D. Prendergast
# May, 2007, UCB

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI
MATHFLAG  = -DUSESCALAPACK
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = /usr/bin/cpp -C
F90free = mpif90 -free
LINK    = mpif90
FOPTS   = -O3 -xT
FNOOPTS = -O2 -xT
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA
CC_COMP = mpicxx
C_COMP  = mpicc
C_LINK  = mpicxx
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWPATH     = /opt/apps/intel11_1/mvapich2_1_6/fftw2/2.1.5
FFTWLIB      = $(FFTWPATH)/lib/libfftw.a
FFTWINCLUDE  = $(FFTWPATH)/include
MKLPATH      = /opt/apps/intel/11.1/mkl/lib/em64t
LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a \
               $(MKLPATH)/libmkl_sequential.a $(MKLPATH)/libmkl_core.a \
               -Wl,--end-group -lpthread
SLKPATH = /opt/apps/intel11_1/mvapich2_1_6/scalapack/1.8.0/lib
SCALAPACKLIB = $(SLKPATH)/libscalapack.a $(SLKPATH)/blacs_MPI-LINUX-0.a \
               $(SLKPATH)/blacsF77init_MPI-LINUX-0.a $(SLKPATH)/blacs_MPI-LINUX-0.a

TESTSCRIPT = qsub lonestar.scr
