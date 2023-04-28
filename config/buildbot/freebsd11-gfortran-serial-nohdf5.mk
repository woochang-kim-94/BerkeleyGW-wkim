# arch.mk for BerkeleyGW codes
#
# Suitable for generic single-processor FreeBSD 11.2
#
# You'll need to instal the following packages/ports:
# gmake gcc6 blas lapack fftw3
#
# and also manually compile HDF5 with Fortran bindings.
#
# Compile with (GNU) gmake, not (BSD) make
#
# All tests passed Jul 2018, r7296
# FHJ

COMPFLAG  = -DGNU
MATHFLAG  = -DUSEFFTW3 -DUNPACKED #-DHDF5
# FHJ: Note: we are not using HDF5 right now because the library
# that comes with the binary package does not include Fortran bindings,

# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

# See: https://lists.freebsd.org/pipermail/freebsd-ports/2016-August/104469.html
_GCC_VER=6
LOCALBASE=/usr/local

MAKE    = gmake
FCPP    = cpp$(_GCC_VER) -C -nostdinc
F90free = gfortran$(_GCC_VER) -ffree-form -ffree-line-length-none
LINK    = gfortran$(_GCC_VER) -Wl,-rpath=$(LOCALBASE)/lib/gcc$(_GCC_VER)
FOPTS   = -O3
FNOOPTS = -O1
MOD_OPT = -J
INCFLAG = -I

CC_COMP = g++$(_GCC_VER)
C_COMP  = gcc$(_GCC_VER)
C_LINK  = g++$(_GCC_VER) -Wl,-rpath=$(LOCALBASE)/lib/gcc$(_GCC_VER)
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/local/lib/ -lfftw3
FFTWINCLUDE  = /usr/local/include
LAPACKLIB    = -L/usr/local/lib/ -llapack -lblas


#HDF5_DIR     =  /opt/hdf5-1.10.1
#HDF5_LDIR    =  $(HDF5_DIR)/lib
#HDF5LIB      =  $(HDF5_LDIR)/libhdf5hl_fortran.a \
#                $(HDF5_LDIR)/libhdf5_hl.a \
#                $(HDF5_LDIR)/libhdf5_fortran.a \
#                $(HDF5_LDIR)/libhdf5.a -lz -ldl
#HDF5INCLUDE  = $(HDF5_DIR)/include

# FHJ: Note that we only run tests that have the "binary" test group, i.e.,
# tests that can work without HDF5 support.
TESTSCRIPT = BGW_CHECK_FLAGS='-g binary' $(MAKE) check
