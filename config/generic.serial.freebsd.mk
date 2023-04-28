# arch.mk for BerkeleyGW codes
#
# Suitable for generic single-processor FreeBSD 9.0
# 
# You'll need to instal the following packages/ports:
# gmake, blas, lapack, fftw3, gcc47
#
# Compile with (GNU) gmake, not (BSD) make
#
# All tests passed 14 Dec 2012, r4956
#
# F. H. da Jornada
# Dec, 2012, UCB

COMPFLAG  = -DGNU
MATHFLAG  = -DUSEFFTW3 -DUNPACKED
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

MAKE    = gmake
FCPP    = cpp47 -C
F90free = gfortran47 -ffree-form -ffree-line-length-none
LINK    = gfortran47
FOPTS   = -O3
FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

CC_COMP = g++47
C_COMP  = gcc47
C_LINK  = g++47
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/local/lib/ -lfftw3
FFTWINCLUDE  = /usr/local/include
LAPACKLIB    = -L/usr/local/lib/ -llapack -lblas
