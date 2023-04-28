# arch.mk for BerkeleyGW codes
#
# suitable for  PARAM Yuva at CDAC, Pune
# all tests passed 12 October 2013, r5710
#
# M. Jain
# Oct 2013, IISc

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI
MATHFLAG  = -DUSESCALAPACK
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C
F90free = mpiifort -free
LINK    = mpiifort
#FOPTS   = -g -O0
FOPTS   = -O2 -assume byterecl -g -traceback -par-report0 -vec-report0
FNOOPTS = -g -O0
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpiicpc
C_COMP  = mpiicc
C_LINK  = mpiicpc
C_OPTS  = -O0 
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWPATH     = /home/external/iisc/manish
FFTWLIB      = -L$(FFTWPATH)/lib -lfftw
FFTWINCLUDE  = $(FFTWPATH)/include
MKLPATH      = /opt/intel/mkl/lib/intel64
LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group
# -lpthread
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a

TESTSCRIPT = qsub yuva.scr
