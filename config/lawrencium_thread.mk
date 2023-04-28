# arch.mk for BerkeleyGW codes
#
# suitable for Lawrencium at LBNL after SL6 upgrade
# need modules intel and openmpi
# uses threaded FFTW3 and threaded MKL with OpenMP
#
# D. Strubbe
# April 2012, UC Berkeley
# all tests passed r4660 3 Aug 2012

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUSEFFTW3 -DUNPACKED
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
DEBUGFLAG = -DDEBUG -DVERBOSE
# NOTE: DEBUGFLAG is uncommented because we use this arch.mk file for a buildslave.

FCPP    = cpp -C
F90free = mpif90 -free -openmp
LINK    = mpif90 -i-static -openmp
FOPTS   = -O3 -no-prec-div
#FNOOPTS = $(FOPTS)
# -O3 gives correct code but takes forever
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
FFTWPATH     = /global/software/sl-6.x86_64/modules/intel/2013.5.192/fftw/3.3.3-intel
FFTWLIB      = $(FFTWPATH)/lib/libfftw3_omp.a $(FFTWPATH)/lib/libfftw3.a
FFTWINCLUDE  = $(FFTWPATH)/include

MKLPATH      = /global/software/sl-6.x86_64/modules/langs/intel/2013.5.192/mkl/lib/intel64
LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -openmp -lpthread
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a 

#Uncomment these lines there to use threads:
#export OMP_NUM_THREADS=2
#export BGW_TEST_MPI_NPROCS=2
TESTSCRIPT = sbatch hbar.scr
