# arch.mk for BerkeleyGW codes
#
# Recommended arch.mk for Stampede at TACC.
# Uses intel compilers + mvapich2
#
# You'll need to run
# module load fftw3/3.3.4 phdf5/1.8.16
#
# Make sure mvapich2 is loaded, and not impi
#
# All tests passed 5 May 2014, r6033
#
# Felipe H. da Jornada
# May 2014, Berkeley

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI -DOMP 
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -ansi
F90free = mpif90 -free -openmp -no-ipo -ip
LINK    = mpif90 -openmp -no-ipo -ip
# We need the -fp-model source to pass the testsuite.
FOPTS   = -O3 -xhost -fp-model source
FNOOPTS = -O2 -xhost -fp-model source -no-ip

#F90free = mpif90 -free -openmp
#LINK    = mpif90 -openmp
#FOPTS   = -g -O0 -check all -Warn all -traceback -fp-stack-check -fstack-protector-all -ftrapuv
#FNOOPTS = $(FOPTS)

MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpicxx -xhost
C_COMP  = mpicc -xhost
C_LINK  = mpicxx -xhost
C_OPTS  = -O3 -xhost -no-ipo -ip -openmp
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      =	$(TACC_FFTW3_LIB)/libfftw3_omp.a \
		$(TACC_FFTW3_LIB)/libfftw3.a
FFTWINCLUDE  = $(TACC_FFTW3_INC)

MKLPATH      = $(MKLROOT)/lib/intel64
LAPACKLIB    = -Wl,--start-group \
		$(MKLPATH)/libmkl_intel_lp64.a \
		$(MKLPATH)/libmkl_intel_thread.a \
		$(MKLPATH)/libmkl_core.a \
		$(MKLPATH)/libmkl_blacs_intelmpi_lp64.a \
		-Wl,--end-group -lpthread -lm
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a

HDF5PATH     = $(TACC_HDF5_LIB)
HDF5LIB      =	$(HDF5PATH)/libhdf5hl_fortran.a \
		$(HDF5PATH)/libhdf5_hl.a \
		$(HDF5PATH)/libhdf5_fortran.a \
		$(HDF5PATH)/libhdf5.a \
		$(HDF5PATH)/libsz.a \
		-lz
HDF5INCLUDE  = $(HDF5PATH)/../include

TESTSCRIPT = sbatch stampede.scr
