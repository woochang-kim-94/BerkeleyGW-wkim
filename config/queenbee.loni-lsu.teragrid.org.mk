# arch.mk for BerkeleyGW codes
#
# suitable for QueenBee at LONI
# testsuite passed 2/4/2011 r2260
#
# D. Strubbe
# November 2010, UC Berkeley

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI
MATHFLAG  = -DUSESCALAPACK
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C
F90free = mpif90 -free
LINK    = mpif90 -i-static
FOPTS   = -O3 -no-prec-div
FNOOPTS = -O2 -no-prec-div
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpicxx
C_COMP  = mpicc
C_LINK  = mpicxx
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWPATH     = /usr/local/packages/fftw/2.1.5/intel-11.1
FFTWLIB      = $(FFTWPATH)/lib/libfftw.a
FFTWINCLUDE  = $(FFTWPATH)/include
MKLPATH      = /usr/local/compilers/Intel/mkl-10.2/lib/em64t
LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a

TESTSCRIPT = qsub queenbee.scr
