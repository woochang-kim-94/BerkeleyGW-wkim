# arch.mk for BerkeleyGW codes
#
# Suitable for Generic single processor linux 
# Tested on netbook running Ubuntu 10.04.  You will need to have yum 
# or apt-get installed liblapack-dev, fftw-dev, gfortran, g++ packages 
#
# J. Deslippe
# Jun, 2008, UCB

COMPFLAG  = -DGNU
MATHFLAG  = 
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

#########################################################################
#   NOTE: This arch.mk is used by a buildslave. The compiler flags are  #
#   optimized to debug the code and not for code performance.           #
#########################################################################
FCPP    = cpp -C
#F90free = gfortran -ffree-form -ffree-line-length-none -Wall
F90free = gfortran -ffree-form -ffree-line-length-none -fbounds-check -Wall -pedantic-errors -std=gnu
LINK    = gfortran
FOPTS   = -O3
FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

CC_COMP = g++ -Wall -pedantic-errors -std=c++0x
C_COMP  = gcc -Wall -pedantic-errors -std=c99
C_LINK  = g++
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/lib/ -lfftw
FFTWINCLUDE  = /usr/share/doc/fftw-dev/fortran
LAPACKLIB    = -L/usr/lib/ -llapack -lblas
