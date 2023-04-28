# arch.mk for BerkeleyGW codes
#
# Suitable for Ubuntu 10.10 in parallel.
# Install packages: liblapack-dev, fftw-dev, gfortran, g++, mpi-default-dev
# Build BLACS according to http://www.open-mpi.org/faq/?category=mpi-apps#blacs
# Build ScaLAPACK according to http://www.open-mpi.org/faq/?category=mpi-apps#scalapack,
#
# D. Strubbe
# January 2011, UCB

COMPFLAG  = -DGNU
PARAFLAG  = -DMPI
MATHFLAG  = -DUSESCALAPACK
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

#########################################################################
#   NOTE: This arch.mk is used by a buildslave. The compiler flags are  #
#   optimized to debug the code and not for code performance.           #
#########################################################################
FCPP    = cpp -C
#F90free = /usr/bin/mpif90 -ffree-form -ffree-line-length-none -Wall
F90free = /usr/bin/mpif90 -ffree-form -ffree-line-length-none -fbounds-check -Wall -pedantic-errors -std=gnu
#F90free = /usr/bin/mpif90 -ffree-form -ffree-line-length-none -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
# denormal not supported, according to runtime error
# precision is flagged by cputime: http://gcc.gnu.org/onlinedocs/gcc-4.5.3/gfortran/Debugging-Options.html
# underflow is flagged by dlamch from lapack
LINK    = /usr/bin/mpif90
FOPTS   = -O3
#FOPTS   = -g -C
FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

C_PARAFLAG  = -DPARA
CC_COMP  = /usr/bin/mpiCC -Wall -pedantic-errors -std=c++0x
C_COMP  = /usr/bin/mpicc -Wall -pedantic-errors -std=c99
C_LINK  = /usr/bin/mpiCC
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/lib/ -lfftw
FFTWINCLUDE  = /usr/share/doc/fftw-dev/fortran
LAPACKLIB    = -L/usr/lib/ -llapack
BLACSDIR     = /usr/lib
BLACS        = $(BLACSDIR)/blacs_MPI-LINUX-1.a $(BLACSDIR)/blacsF77init_MPI-LINUX-1.a $(BLACSDIR)/blacs_MPI-LINUX-1.a
SCALAPACKLIB = /usr/lib/libscalapack-1.a $(BLACS)

#Ubuntu BLACS/SCALAPACK packages don't work properly sometimes
#BLACS        = -lblacs-openmpi -lblacsF77init-openmpi -lblacs-openmpi
#SCALAPACK    = -lscalapack-openmpi $(BLACS)

#MKL ScaLAPACK fails completely in complex version
#MKLPATH      = /auto/opt/intel/mkl/lib/intel64
#LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_gf_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread
#SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a 

#need to export MPIEXEC=/usr/bin/mpirun if this is not default in `which mpiexec`
TESTSCRIPT = make check-parallel
