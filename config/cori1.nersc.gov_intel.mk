# arch.mk for BerkeleyGW codes
#
# suitable for Cori Phase 1 at NERSC
#
# JRD
# 2016, NERSC
#
# Run the following command before compiling:
# module load cray-hdf5-parallel && module unload darshan && module swap intel intel/17.0.2.174
#
# Status: most tests pass (r7121). The failing tests have only a small
# deviation, and are due the different k-point picked for the interpolation.
# This should cause no harm. --FHJ (02/25/2017)
#
# Precompiler options

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DUSEMR3
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = /usr/bin/cpp -C -nostdinc
F90free = ftn -free -qopenmp
LINK    = ftn -qopenmp
FOPTS   = -fast -no-ip -no-ipo -align array64byte
#FOPTS   = -fast -no-ip -no-ipo -align array64byte -g -debug inline-debug-info -traceback -check all -ftrapuv -init=snan
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = CC -qopenmp
C_COMP  = cc -qopenmp
C_LINK  = CC -qopenmp
C_OPTS  = -fast -no-ip -no-ipo -align #-g -traceback
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
#FFTWPATH     = 
FFTWLIB      = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a \
               $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl 
FFTWINCLUDE  = $(MKLROOT)/include/fftw/

HDF5_LDIR    =  $(HDF5_DIR)/lib
HDF5LIB      =  $(HDF5_LDIR)/libhdf5hl_fortran.a \
                $(HDF5_LDIR)/libhdf5_hl.a \
                $(HDF5_LDIR)/libhdf5_fortran.a \
                $(HDF5_LDIR)/libhdf5.a -lz -ldl
HDF5INCLUDE  = $(HDF5_DIR)/include

PERFORMANCE  = 

LAPACKLIB = $(FFTWLIB)

TESTSCRIPT = sbatch cori1.scr
