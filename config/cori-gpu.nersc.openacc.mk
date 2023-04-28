# arch.mk for BerkeleyGW codes
#
# suitable for Cori-GPU
#
# MDB
# 2021, NERSC
# PGI compiler
#
# Do:
# module purge ; module load cgpu ;  module load nvhpc/21.3 ; module load cuda/11.2.2 ; module load openmpi/4.0.3
#
COMPFLAG  = -DPGI
PARAFLAG  = -DMPI  -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DOPENACC

FCPP    = /usr/bin/cpp  -P -ansi  -nostdinc  -C  -E  -std=c11
F90free = mpif90 -Mfree -acc -mp -gpu=cc70 -Mcuda=cuda11.2 -Mcudalib=cublas,cufft -Mcuda=lineinfo -traceback -Minfo=mp
LINK    = mpif90        -acc -mp -gpu=cc70 -Mcuda=cuda11.2 -Mcudalib=cublas,cufft -Minfo=mp
FOPTS   = -fast -Mfree -Mlarge_arrays
FNOOPTS = $(FOPTS)
#MOD_OPT =  -J
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpiCC
C_COMP  = mpicc
C_LINK  = mpiCC
C_OPTS  = -fast -mp
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

FFTW_DIR=/project/projectdirs/nesap/pgi-19.7/fftw-3.3.8/
FFTWLIB      = $(FFTW_DIR)/lib/libfftw3.a \
               $(FFTW_DIR)/lib/libfftw3_omp.a \
               /usr/common/software/pgi/19.7/linux86-64/19.7/lib/libblas.a \
               ${CUDALIB}  -lstdc++
FFTWINCLUDE  =$(FFTW_DIR)/include/
PERFORMANCE  =

SCALAPACKLIB = /project/projectdirs/nesap/pgi-19.7/scalapack-2.0.2/libscalapack.a  /usr/common/software/pgi/19.7/linux86-64/19.7/lib/liblapack.a  /usr/common/software/pgi/19.7/linux86-64/19.7/lib/libblas.a
LAPACKLIB = /usr/common/software/pgi/19.7/linux86-64/19.7/lib/libblas.a /usr/common/software/pgi/19.7/linux86-64/19.7/lib/liblapack.a /usr/common/software/pgi/19.7/linux86-64/19.7/lib/libblas.a

HDF5_LDIR    =  /project/projectdirs/nesap/pgi-19.7/hdf5-1.8.17/lib/
HDF5LIB      =  $(HDF5_LDIR)/libhdf5hl_fortran.a \
                $(HDF5_LDIR)/libhdf5_hl.a \
                $(HDF5_LDIR)/libhdf5_fortran.a \
                $(HDF5_LDIR)/libhdf5.a -lz -ldl
HDF5INCLUDE  =/project/projectdirs/nesap/pgi-19.7/hdf5-1.8.17/include
