# arch.mk for BerkeleyGW codes
#
# Suitable for tigger.berkeley.edu with MPICH1
# Needs fftw-dev and installation of PGI compiler suite
#
# D. Strubbe
# January 2011, UCB

COMPFLAG  = -DPGI
PARAFLAG  = -DMPI
MATHFLAG  = -DUSESCALAPACK -DUNPACKED
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C
F90free = mpif90 -Mfree
LINK    = mpif90
FOPTS   = -fast
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpiCC
C_COMP  = mpicc
C_LINK  = mpiCC
C_OPTS  = -fast
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/lib/ -lfftw
FFTWINCLUDE  = /usr/share/doc/fftw-dev/fortran
LAPACKLIB    = -lacml
SCALAPACKPATH = /opt/pgi/linux86-64/10.9/mpi/mpich/lib
SCALAPACKLIB = $(SCALAPACKPATH)/libscalapack.a $(SCALAPACKPATH)/blacs_MPI-LINUX-0.a \
$(SCALAPACKPATH)/blacsF77init_MPI-LINUX-0.a $(SCALAPACKPATH)/blacs_MPI-LINUX-0.a

TESTSCRIPT = MPIEXEC=/opt/pgi/linux86-64/10.9/mpi/mpich/bin/mpiexec make check-parallel
