# arch.mk for BerkeleyGW codes
#
# open64 serial
# needs packages liblapack-dev, fftw-dev, g++ packages 
#
# D. Strubbe
# November 2010, UCB

COMPFLAG  = -DOPEN64
MATHFLAG  = -DUNPACKED
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -P -C
F90free = openf90 -freeform -fno-second-underscore -intrinsic=EVERY
LINK    = openf90
FOPTS   = -O3
# -O3 generates correct code but takes forever
FNOOPTS = -O2
MOD_OPT = -module 
INCFLAG = -I

CC_COMP = openCC
C_COMP  = opencc
C_LINK  = openCC
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/lib/ -lfftw
FFTWINCLUDE  = /usr/share/doc/fftw-dev/fortran
LAPACKLIB    = -L/usr/lib/ -llapack
