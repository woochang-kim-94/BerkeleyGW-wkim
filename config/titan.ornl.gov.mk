# arch.mk for BerkeleyGW codes
#
# suitable for Titan at ORNL with PGI
# all tests passed in r5035
#
# D. Strubbe
# Jan 2013, MIT

COMPFLAG  = -DPGI
PARAFLAG  = -DMPI
#PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED
#MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = /usr/bin/cpp -C
F90free = ftn -Mfree -mp
LINK    = ftn -mp
#F90free = ftn -Mfree
#LINK    = ftn
FOPTS   = -fast
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = CC
C_COMP  = cc
C_LINK  = CC
C_OPTS  = -fast
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWPATH     = /opt/fftw/2.1.5.4
FFTWLIB      = -L$(FFTWPATH)/lib -ldfftw
FFTWINCLUDE  = $(FFTWPATH)/include
#HDF5PATH     = /opt/cray/hdf5-parallel/1.8.8/pgi/119
#HDF5LIB      = -L$(HDF5PATH)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
#HDF5INCLUDE  = $(HDF5PATH)/include

PERFORMANCE  = 

TESTSCRIPT = qsub titan.scr
