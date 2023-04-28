# arch.mk for BerkeleyGW codes
#
# suitable for Summit ORNL
#
# MDB
# 2021, ORNL
# 
# Do:
# module swap xl  pgi/20.4;  module load netlib-lapack ; module load netlib-scalapack ; module load fftw ; module load hdf5 ; module load essl ; module load cuda
#
#
COMPFLAG  = -DPGI
PARAFLAG  = -DMPI  -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DOPENACC 
#

FCPP    = /usr/bin/cpp  -P -ansi  -nostdinc   #  -C  -P  -E  -nostdinc
F90free = mpif90 -Mfree -mp -acc -ta=tesla -Mcuda -Mcudalib=cublas,cufft -Mcuda=lineinfo -Minfo=mp # -g -traceback
LINK    = mpif90 -Mfree -mp -acc -ta=tesla -Mcuda -Mcudalib=cublas,cufft -Mcuda=lineinfo -Minfo=mp # -g -traceback
FOPTS   = -O1 # -fast 
FNOOPTS = $(FOPTS)
#MOD_OPT =  -J
MOD_OPT =  -module  
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpiCC 
C_COMP  = mpicc 
C_LINK  = mpiCC 
C_OPTS  = -O1 -mp # -fast
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# this must be linked if non ESSL blas library is missing are missing 
#               -L$(OLCF_NETLIB_LAPACK_ROOT)/lib64/ -llapack -lblas \

FFTWLIB      = \
               -L$(OLCF_ESSL_ROOT)/lib64/ -lessl \
               -L$(OLCF_FFTW_ROOT)/lib/ -lfftw3 -lfftw3_threads -lfftw3_omp \
               -L$(OLCF_NETLIB_SCALAPACK_ROOT)/lib/ -lscalapack  \
                 ${CUDALIB}  -lstdc++
FFTWINCLUDE  = $(OLCF_FFTW_ROOT)/include/  
PERFORMANCE  =

HDF5_LDIR    = $(OLCF_HDF5_ROOT)/lib/
HDF5LIB      = $(HDF5_LDIR)/libhdf5hl_fortran.a \
               $(HDF5_LDIR)/libhdf5_hl.a \
               $(HDF5_LDIR)/libhdf5_fortran.a \
               $(HDF5_LDIR)/libhdf5.a -lm -lz -ldl  -lstdc++ 
HDF5INCLUDE  = $(HDF5_LDIR)/../include


LAPACKLIB = -L$(OLCF_ESSL_ROOT)/lib64/ -lessl  -L$(OLCF_NETLIB_LAPACK_ROOT)/lib64/ -llapack 
SCALAPACKLIB = -L$(OLCF_NETLIB_SCALAPACK_ROOT)/lib/ -lscalapack 
#
TESTSCRIPT = 
#
