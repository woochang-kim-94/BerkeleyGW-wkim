# arch.mk for BerkeleyGW codes
#
# Suitable for Comet at SDSC with gnu + openmpi_ib + atlas environment.
# Supports MPI, OpenMP and HDF5.
#
# FHJ: We are definitely not using the most optimized linear algebra libraries
# and MPI implementation available, but this arch.mk in one of the few stable
# ones for this machine.
# 
# To compile the code, you'll need to:
# module purge && module load gnutools gnu openmpi_ib hdf5 fftw lapack atlas
#
# To run the code, you'll need to include the following lines in your script:
# module purge && module load gnu openmpi_ib
#
# All tests pass @ r6954.
#
# Felipe H. da Jornada

COMPFLAG  = -DGNU
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -ansi
F90free = mpif90 -ffree-form -ffree-line-length-none -fopenmp -fno-second-underscore
LINK    = mpif90 -fopenmp -ldl
FOPTS   = -O3 -funroll-loops -funsafe-math-optimizations
FNOOPTS = -O3 -funroll-loops -funsafe-math-optimizations
#FOPTS   = -g -O0 -ffpe-trap=invalid -fbounds-check -fbacktrace -Wall
#FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpicxx
C_COMP  = mpicc
C_LINK  = mpicxx
C_OPTS  = -O3 -fopenmp
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
FFTWLIB      =	$(FFTWHOME)/lib/libfftw3_omp.a \
		$(FFTWHOME)/lib/libfftw3.a
FFTWINCLUDE  = $(FFTWHOME)/include

SCALAPACKHOME = /opt/scalapack/gnu/openmpi_ib
LAHOME       = $(ATLASHOME)/lib
LAPACKLIB    = $(LAHOME)/liblapack.a $(LAHOME)/libf77blas.a $(LAHOME)/libcblas.a $(LAHOME)/libatlas.a
SCALAPACKLIB = $(SCALAPACKHOME)/lib/libscalapack.a

HDF5PATH     = $(HDF5HOME)/lib
HDF5LIB      =	$(HDF5PATH)/libhdf5hl_fortran.a \
		$(HDF5PATH)/libhdf5_hl.a \
		$(HDF5PATH)/libhdf5_fortran.a \
		$(HDF5PATH)/libhdf5.a \
		-lz
HDF5INCLUDE  = $(HDF5HOME)/include

TESTSCRIPT = sbatch comet.scr
