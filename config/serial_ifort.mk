# arch.mk for BerkeleyGW codes
#
# single processor with ifort, MKL, ia32
# note that Ubuntu 10.10 is 'unsupported' by Intel
# but it seems to work fine nonetheless
#
# D. Strubbe
# Nov 2010, UCB

COMPFLAG  = -DINTEL
MATHFLAG  = 
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C
F90free = ifort -free
LINK    = ifort
FOPTS   = -O3
# -O3 generates correct code but takes forever
FNOOPTS = -O2
MOD_OPT = -module 
INCFLAG = -I

CC_COMP = icc
C_COMP  = icc
C_LINK  = icc
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/lib -lfftw
FFTWINCLUDE  = /usr/share/doc/fftw-dev/fortran
MKLPATH = /auto/opt/intel/mkl/lib/ia32
LAPACKLIB = -Wl,--start-group $(MKLPATH)/libmkl_intel.a $(MKLPATH)/libmkl_sequential.a \
            $(MKLPATH)/libmkl_core.a -Wl,--end-group -lpthread
