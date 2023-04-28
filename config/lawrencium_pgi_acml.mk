# arch.mk for BerkeleyGW codes
#
# suitable for Lawrencium at LBNL after SL6 upgrade
# PGI compiler, ACML LAPACK, MKL ScaLAPACK
#
# Run the following command before compiling:
# module load pgi/13.10, openmpi

# do not use PGI 12.8, which would have the compiler error in Sigma:
#Lowering Error: symbol components$sd is an inconsistent array descriptor
#PGF90-F-0000-Internal compiler error. Errors in Lowering       1 (input.f90: 813)
#PGF90/x86-64 Linux 12.8-0: compilation aborted

#
# D. Strubbe
# Aug 2012, MIT

COMPFLAG  = -DPGI
PARAFLAG  = -DMPI
MATHFLAG  = -DUNPACKED -DUSESCALAPACK
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
DEBUGFLAG = -DDEBUG -DVERBOSE
# NOTE: DEBUGFLAG is uncommented because we use this arch.mk file for a buildslave.

FCPP    = cpp -C
F90free = mpif90 -Mfree
LINK    = mpif90
FOPTS   = -fast
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA
CC_COMP = mpiCC
C_COMP  = mpicc
C_LINK  = mpiCC
C_OPTS  = -fast
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWPATH     = /global/software/sl-6.x86_64/modules/pgi/13.10/fftw/2.1.5-pgi
FFTWLIB      = $(FFTWPATH)/lib/libdfftw.a
FFTWINCLUDE  = $(FFTWPATH)/include
LAPACKLIB    =  /global/software/sl-6.x86_64/modules/pgi/13.10/acml/5.3.1-pgi/pgi64/lib/libacml.a
MKLPATH      = /global/software/sl-6.x86_64/modules/langs/intel/2013.5.192/mkl/lib/intel64
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread

TESTSCRIPT = sbatch hbar.scr
