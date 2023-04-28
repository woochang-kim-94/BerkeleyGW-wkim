# arch.mk for BerkeleyGW codes
#
# Configuration file for a singularity (http://singularity.lbl.gov/) image
# running ubuntu 16.04 or 17.10, compiled with gfortran, openmpi, and hdf5.
# Note that BerkeleyGW should be compiled (and run) inside the singularity
# image. Eg:
# singularity exec -B /run:/run <image.img> make all-flavors -j 4
#
# All tests passed, r7196
#
# Felipe H da Jornada
# Nov 2017

COMPFLAG  = -DGNU
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DUSEMR3 #-DUSEELPA -DUSEPRIMME
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

#########################################################################
#   NOTE: This arch.mk is used by a buildslave. The compiler flags are  #
#   optimized to debug the code and not for code performance.           #
#########################################################################
FCPP    = cpp -C -nostdinc
F90free = mpif90 -ffree-form -ffree-line-length-none -Wall -pedantic-errors -std=gnu -fopenmp
LINK    = mpif90 -fopenmp
FOPTS   = -g -O3 -fbacktrace -fcheck=bounds,do,mem,pointer,recursion
FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

C_PARAFLAG  = -DPARA
CC_COMP  = mpiCC -Wall -pedantic-errors -std=c++0x -fopenmp
C_COMP  = mpicc -Wall -pedantic-errors -std=c99 -fopenmp
C_LINK  = mpiCC -fopenmp
C_OPTS  = -g -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = /opt/fftw-3.3.5/lib/libfftw3_omp.a /opt/fftw-3.3.5/lib/libfftw3.a
FFTWINCLUDE  = /opt/fftw-3.3.5/include

LAPACKLIB    = /lib64/liblapack.a /opt/OpenBLAS-0.2.20/lib/libopenblas.a
SCALAPACKLIB = /lib64/libscalapack.a

HDF5_LDIR    =  /usr/lib/x86_64-linux-gnu
HDF5LIB      =  $(HDF5_LDIR)/libhdf5_openmpihl_fortran.a \
                $(HDF5_LDIR)/libhdf5_openmpi_hl.a \
                $(HDF5_LDIR)/libhdf5_openmpi_fortran.a \
                $(HDF5_LDIR)/libhdf5_openmpi.a -lz -lsz -ldl
HDF5INCLUDE  = /usr/include/hdf5/openmpi

TESTSCRIPT = bash singularity.scr
