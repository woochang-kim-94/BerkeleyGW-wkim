# arch.mk for BerkeleyGW codes
#
# Configuration file for a buildbot worker running on Amazon EC2.
# Ubuntu 20.04, intel compiler + MKL + HDF5, threaded
#
# All tests passed for commit d608418a4a53990ed7b652b225120aebbcd8df05
#
# Mauro Del Ben
# Aug 2019
# (updated by Will Huhn, March 2021)

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DUSEMR3 # -DUSEPRIMME
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

#########################################################################
#   NOTE: This arch.mk is used by the testsuite. The compiler flags are #
#   optimized to debug the code and not for code performance.           #
#########################################################################
FCPP    = /usr/bin/cpp -C  -nostdinc #  -ansi
# F90free = mpiifort -free -O0 -qopenmp  -g -traceback  -debug inline-debug-info -traceback -check all
# LINK    = mpiifort -free -O0 -qopenmp  -g -traceback  -debug inline-debug-info -traceback -check all
# FOPTS   = -O0
#
# F90free = mpiifort -free -qopenmp -O0 -fp-model precise -no-ip -no-ipo -align array64byte -g -debug inline-debug-info -traceback -check all,noarg_temp_created -ftrapuv -init=snan -fpe0
# LINK    = mpiifort -free -qopenmp -O0 -fp-model precise -no-ip -no-ipo -align array64byte -g -debug inline-debug-info -traceback -check all,noarg_temp_created -ftrapuv -init=snan -fpe0
# FOPTS   = -O0
#
F90free = mpiifort  -free -O3 -qopenmp  -g  -traceback
LINK    = mpiifort  -free -O3 -qopenmp  -g  -traceback
FOPTS   = -O3
FNOOPTS = $(FOPTS)
MOD_OPT = -module ./
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP  = mpiicpc -qopenmp
C_COMP   = mpiicc  -qopenmp
C_LINK   = mpiicc  -qopenmp
C_OPTS   = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
               $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
               -Wl,--end-group -lpthread -lm  # /global/homes/m/mdelben/cori/BGW_subspace/dlopen_lib/dlopen.o
FFTWINCLUDE  = $(MKLROOT)/include/fftw/

LAPACKLIB = $(FFTWLIB)

HDF5_LDIR    =  /opt/lib/hdf5-1.10.5/lib
HDF5LIB      =  $(HDF5_LDIR)/libhdf5hl_fortran.a \
                $(HDF5_LDIR)/libhdf5_hl.a \
                $(HDF5_LDIR)/libhdf5_fortran.a \
                $(HDF5_LDIR)/libhdf5.a -lz -ldl
HDF5INCLUDE  = /opt/lib/hdf5-1.10.5/include

TESTSCRIPT = bash singularity.scr
