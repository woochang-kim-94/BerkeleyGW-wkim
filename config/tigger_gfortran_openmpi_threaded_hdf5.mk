# arch.mk for BerkeleyGW codes
#
# Suitable for Ubuntu 10.10 in parallel.
# Install packages: liblapack-dev, fftw-dev, gfortran, g++, mpi-default-dev
# Build BLACS according to http://www.open-mpi.org/faq/?category=mpi-apps#blacs
# Build ScaLAPACK according to http://www.open-mpi.org/faq/?category=mpi-apps#scalapack,
#
# Jamal I. Mustafa
# July 2014, UC Berkeley

COMPFLAG  = -DGNU
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DHDF5 -DUSEFFTW3
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

#########################################################################
#   NOTE: This arch.mk is used by a buildslave. The compiler flags are  #
#   optimized to debug the code and not for code performance.           #
#########################################################################
FCPP    = cpp -C
#F90free = /usr/bin/mpif90 -ffree-form -ffree-line-length-none -Wall
F90free = /usr/bin/mpif90 -ffree-form -ffree-line-length-none -fbounds-check -Wall -pedantic-errors -std=gnu -fopenmp
#F90free = /usr/bin/mpif90 -ffree-form -ffree-line-length-none -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
# denormal not supported, according to runtime error
# precision is flagged by cputime: http://gcc.gnu.org/onlinedocs/gcc-4.5.3/gfortran/Debugging-Options.html
# underflow is flagged by dlamch from lapack
LINK    = /usr/bin/mpif90 -fopenmp
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
FFTWLIB      = -L/usr/lib/ -lfftw3 -lfftw3_threads
FFTWINCLUDE  = /usr/include -I /home/buildbot/fftw-3.3.3/include
# FHJ: Note: manually compiled threaded LAPACK and ScaLAPACK on tigger
LAPACKLIB    = /opt/lib/liblapack.a -L/usr/lib -lblas
SCALAPACKLIB = /opt/lib/libscalapack.a

HDF5LIB      = -L/usr/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
HDF5INCLUDE  = /usr/include

#need to export MPIEXEC=/usr/bin/mpirun if this is not default in `which mpiexec`
TESTSCRIPT = make check-parallel
