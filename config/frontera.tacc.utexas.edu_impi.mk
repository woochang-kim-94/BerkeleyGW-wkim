# arch.mk for BerkeleyGW codes
#
# suitable for Frontera at TACC, Cascade Lake (CLX) architecture
# Uses intel compilers + impi
#
# BERKELEYGW DOES NOT WORK WITH INTEL/18 DUE TO A COMPILER BUG!
# ------------------------------------------------------------
#
# You'll need to run:
# module load arpack impi intel/19.0.5 phdf5
#
# Use 'make -j 8' for parallel build on Frontera
#
# Zhenglu Li
# July 2020, Berkeley

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI -DOMP 
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C -nostdinc
F90free = mpiifort -xCORE-AVX512 -free -qopenmp -ip -no-ipo
LINK    = mpiifort -xCORE-AVX512 -qopenmp -ip -no-ipo
# We need the -fp-model procise to pass the testsuite.
FOPTS   = -O3 -fp-model source
FNOOPTS = -O2 -fp-model source -no-ip
#FOPTS   = -g -O0 -check all -Warn all -traceback
#FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpiicpc -xCORE-AVX512
C_COMP  = mpiicc -xCORE-AVX512
C_LINK  = mpiicpc -xCORE-AVX512
C_OPTS  = -O3 -ip -no-ipo -qopenmp
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
MKLPATH      = $(MKLROOT)/lib/intel64

FFTWLIB      =	-Wl,--start-group \
		$(MKLPATH)/libmkl_intel_lp64.a \
		$(MKLPATH)/libmkl_intel_thread.a \
		$(MKLPATH)/libmkl_core.a \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl
FFTWINCLUDE  = $(MKLROOT)/include/fftw


LAPACKLIB    = -Wl,--start-group \
		$(MKLPATH)/libmkl_intel_lp64.a \
		$(MKLPATH)/libmkl_intel_thread.a \
		$(MKLPATH)/libmkl_core.a \
		$(MKLPATH)/libmkl_blacs_intelmpi_lp64.a \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a

HDF5PATH     = $(TACC_HDF5_LIB)
HDF5LIB      =	$(HDF5PATH)/libhdf5hl_fortran.a \
		$(HDF5PATH)/libhdf5_hl.a \
		$(HDF5PATH)/libhdf5_fortran.a \
		$(HDF5PATH)/libhdf5.a \
		$(HDF5PATH)/libsz.a \
		-lz
HDF5INCLUDE  = $(HDF5PATH)/../include

TESTSCRIPT = 
