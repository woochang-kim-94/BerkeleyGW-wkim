# arch.mk for BerkeleyGW codes
#
# single processor with ifort, MKL, 32-bit on x86_64, suitable for Tigger
# note that Ubuntu 10.10 is 'unsupported' by Intel
# but it seems to work fine nonetheless
#
# D. Strubbe
# Jan 2011, UCB

COMPFLAG  = -DINTEL
MATHFLAG  = 
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C
F90free = ifort -free
LINK    = ifort
# -warn all can often lead to compiler crashes!
FOPTS   = -O3
#FOPTS   = -O3 -g -check all -Warn all -traceback
# -O3 generates correct code, but especially with -warn all takes forever and can crash the compiler
FNOOPTS = -O3
#FNOOPTS = -O3 -g -check all -Warn all -traceback
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
MKLPATH = /auto/opt/intel/mkl/lib/intel64
LAPACKLIB = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
            $(MKLPATH)/libmkl_core.a -Wl,--end-group -lpthread
