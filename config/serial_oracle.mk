# arch.mk for BerkeleyGW codes
#
# FHJ: FIXME - no longer supported
#
# SolarisStudio12.2-linux-x86 serial
# needs packages liblapack-dev, fftw-dev
#
# D. Strubbe
# December 2010, UCB

COMPFLAG  = -DORACLE
MATHFLAG  = 

#########################################################################
#   NOTE: This arch.mk is used by a buildslave. The compiler flags are  #
#   optimized to debug the code and not for code performance.           #
#########################################################################
FCPP    = cpp -P -C
#F90free = sunf90 -free
F90free = sunf90 -free -C -O0 -g -xcheck
LINK    = sunf90
FOPTS   = #-O3
# sigma_main may cause segfault or system crash at compile time if -O3!
FNOOPTS = #-O2
MOD_OPT = -moddir=
INCFLAG = -I

# sunCC segfaults mysteriously at the very beginning of surface.cpp and icm.cpp
CC_COMP = g++ -g -g3
C_COMP  = gcc -g -g3
C_LINK  = g++
C_OPTS  = -O0
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/lib/ -lfftw
FFTWINCLUDE  = /usr/share/doc/fftw-dev/fortran
LAPACKLIB    = -L/usr/lib/ -llapack -lblas
