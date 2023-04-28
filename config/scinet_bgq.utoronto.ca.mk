# arch.mk for BerkeleyGW codes
#
# suitable for BGQ machines at SciNet
# https://support.scinet.utoronto.ca/wiki/index.php/BGQ
#
# D.Gruner
# Jan 2013, SciNet
#
# Run the following command before compiling:
# module load   xlf/14.1       vacpp/12.1     mpich2/xl      essl/5.1       lapack/3.4.2   fftw/2.1.5

# Precompiler options
#
COMPFLAG  = -DXLF
PARAFLAG  = -DMPI
MATHFLAG  = -DUSESCALAPACK
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = $(SCINET_XLF_BASE)/exe/cpp
F90free = mpif90
LINK    = mpif90
FOPTS   = -O3 -qhot -qarch=qp -qtune=qp -qalias=noaryovrlp
FNOOPTS = -O2
MOD_OPT = -qmoddir=
INCFLAG = -I

C_PARAFLAG = -DPARA
CC_COMP = mpicxx
C_COMP  = mpicc
C_LINK  = mpicxx
C_OPTS  = -O3 -qarch=qp -qtune=qp
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWPATH     = $(SCINET_FFTW_BASE)
FFTWLIB      = -L$(FFTWPATH)/lib -ldfftw
FFTWINCLUDE  = $(FFTWPATH)/include
LAPACKLIB    = $(SCINET_ESSL_LIB)/libesslbg.a \
		$(SCINET_LAPACK_LIB)/liblapack.a  \
		$(SCINET_LAPACK_LIB)/liblapacke.a  \
		$(SCINET_LAPACK_LIB)/librefblas.a  \
		$(SCINET_LAPACK_LIB)/libtmglib.a
SCALAPACKLIB = $(SCINET_LAPACK_LIB)/libscalapack.a 

TESTSCRIPT =
