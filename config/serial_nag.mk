# arch.mk for BerkeleyGW codes
#
# NAG serial
# needs packages liblapack-dev, fftw-dev, g++ packages 
# "If you have a dead NAG, you need to buy a new donkey."
# Status: all tests passed -O0. 19 Oct. 2011, r3449
# absorption is broken at -O1 and higher due to something in Common.
#
# D. Strubbe
# September 2011, UCB

COMPFLAG  = -DNAG
MATHFLAG  = 
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C
F90free = nagfor -free
LINK    = nagfor
FOPTS   = -O0
FNOOPTS = $(FOPTS)
MOD_OPT = -mdir 
INCFLAG = -I

CC_COMP = g++
C_COMP  = gcc
C_LINK  = g++
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/lib/ -lfftw
FFTWINCLUDE  = /usr/share/doc/fftw-dev/fortran
LAPACKLIB    = -L/usr/lib/ -llapack
