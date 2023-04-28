# arch.mk for BerkeleyGW codes
#
# suitable for Lawrencium at LBNL with HDF5
#
# Run the following command before compiling:
# module load intel/2013_sp1.4.211 openmpi hdf5/1.8.13-intel-p
#
# Last uptdate: July 2016, T. Rangel
# D. Strubbe
# Jan 2013, MIT

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI
MATHFLAG  = -DUSESCALAPACK -DHDF5 -DUNPACKED
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C
F90free = mpif90 -free
LINK    = mpif90 -i-static -openmp
FOPTS   = -O3 -no-prec-div
FNOOPTS = -O2 -no-prec-div
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG = -DPARA
CC_COMP = mpiCC
C_COMP  = mpicc
C_LINK  = mpiCC
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWPATH     = /global/software/sl-6.x86_64/modules/intel/2013.5.192/fftw/2.1.5-intel
FFTWLIB      = $(FFTWPATH)/lib/libdfftw.a
FFTWINCLUDE  = $(FFTWPATH)/include
MKLPATH      = /global/software/sl-6.x86_64/modules/langs/intel/2013.5.192/mkl/lib/intel64
LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a
HDF5DIR      = /global/software/sl-6.x86_64/modules/intel/2013.5.192/hdf5/1.8.11-intel-p
HDF5LIB      = -L$(HDF5DIR)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
HDF5INCLUDE  = $(HDF5DIR)/include

TESTSCRIPT = sbatch hbar.scr
