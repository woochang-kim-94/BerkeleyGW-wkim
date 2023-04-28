# arch.mk for BerkeleyGW codes
#
# suitable for Cori Phse 2 (KNL) at NERSC
#
# JRD
# 2016, NERSC
#
# MW
# 2020, use shared libs
# Run the following command before compiling:
# module swap craype-haswell craype-mic-knl && module unload darshan && module load cray-hdf5-parallel && module swap intel intel/19.0.0.117
#
# Precompiler options

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DUSEMR3 #-DUSEELPA -DUSEPRIMME
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = /usr/bin/cpp -C -nostdinc
F90free = ftn -shared -free -qopenmp
LINK    = ftn -shared -qopenmp
FOPTS   = -fast -no-ip -no-ipo -align array64byte
#FOPTS   = -fast -no-ip -no-ipo -align array64byte -traceback
#FOPTS   = -fast -no-ip -no-ipo -align array64byte -g -debug inline-debug-info -traceback -check all -ftrapuv -init=snan
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = CC -qopenmp -shared
C_COMP  = cc -qopenmp -shared
C_LINK  = CC -qopenmp -shared
C_OPTS  = -fast -no-ip -no-ipo -align #-g -traceback
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWPATH     = 
FFTWLIB      = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a \
               $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl 

#FFTWLIB      = 
FFTWINCLUDE  = $(MKLROOT)/include/fftw/

HDF5_LDIR    =  $(HDF5_DIR)/lib
HDF5LIB      =  $(HDF5_LDIR)/libhdf5hl_fortran.a \
                $(HDF5_LDIR)/libhdf5_hl.a \
                $(HDF5_LDIR)/libhdf5_fortran.a \
                $(HDF5_LDIR)/libhdf5.a -lz -ldl
HDF5INCLUDE  = $(HDF5_DIR)/include

PERFORMANCE  = 

LAPACKLIB = $(FFTWLIB)
#PRIMMELIB = ${HOME}/codes/primme-2.1/lib/libprimme.a

#ELPA_DIR = ${HOME}/codes/elpa/2017.05.002-libsci-knl-omp
#ELPAINCLUDE = ${ELPA_DIR}/include/elpa/modules
#ELPALIB = ${ELPA_DIR}/lib/libelpa_openmp.a

TESTSCRIPT = sbatch cori2.scr
