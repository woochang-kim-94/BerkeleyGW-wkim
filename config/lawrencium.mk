# arch.mk for BerkeleyGW codes
#
# MDB
# 2018, LBL
# module load intel/2016.4.072
# module load openmpi/2.0.2-intel
# module load mkl/2016.4.072
# module load fftw/3.3.6-intel
# module load hdf5/1.8.18-intel-p
#
# The same as above
# module load intel/2016.4.072 ; module load openmpi/2.0.2-intel ; module load mkl/2016.4.072 ; module load fftw/3.3.6-intel ; module load hdf5/1.8.18-intel-p

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI  -DOMP
MATHFLAG  = -DUSESCALAPACK  -DUSEFFTW3  -DUNPACKED  -DHDF5 #  -DUSEELPA
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
# DEBUGFLAG = -DDEBUG -DVERBOSE
# NOTE: DEBUGFLAG is uncommented because we use this arch.mk file for a buildslave.

FCPP    = cpp -ansi
F90free = mpif90 -free -qopenmp
LINK    = mpif90 -qopenmp  #  -static
FOPTS   = -O3 -no-prec-div
FNOOPTS = -O2 -no-prec-div
MOD_OPT = -module ./
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
FFTWPATH     = /global/software/sl-7.x86_64/modules/intel/2016.4.072/fftw/3.3.6-intel
FFTWLIB      = $(FFTWPATH)/lib/libfftw3_omp.a  $(FFTWPATH)/lib/libfftw3.a
FFTWINCLUDE  = $(FFTWPATH)/include
MKLPATH      = /global/software/sl-7.x86_64/modules/langs/intel/2016.4.072/mkl/lib/intel64
LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_openmpi_lp64.a -Wl,--end-group  -qopenmp  -lpthread
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a
HDF5DIR      = /global/software/sl-7.x86_64/modules/intel/2016.4.072/hdf5/1.8.18-intel-p
HDF5LIB      = -L$(HDF5DIR)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
HDF5INCLUDE  = $(HDF5DIR)/include

# ELPALIB = /global/home/users/mdelben/ELPA/elpa-2017.05.002/lib/libelpa.a -lstdc++ $(LAPACKLIB) $(SCALAPACKLIB)
# ELPAINCLUDE = /global/home/users/mdelben/ELPA/elpa-2017.05.002/modules

TESTSCRIPT = sbatch hbar.scr
