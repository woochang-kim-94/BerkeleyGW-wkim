# arch.mk for BerkeleyGW codes
#
# open64 serial, 64-bit
# needs packages liblapack-dev, fftw-dev, g++ packages 
# tested with open64-4.2.1

# Status:
# openCC doesn't work with gcc 4.4, only with 4.2 or 4.1
# openmp test compilation fails saying
# /opt/open64-4.2.1//lib/gcc-lib/x86_64-open64-linux/4.2.1/libopenmp.a(omp_init.o):(.eh_frame+0x12): undefined reference to `__gxx_personality_v0'

# D. Strubbe
# January 2011, UCB

COMPFLAG  = -DOPEN64
MATHFLAG  = -DUNPACKED
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -P -C
F90free = openf90 -m64 -freeform -fno-second-underscore -intrinsic=EVERY
LINK    = openf90 -m64
FOPTS   = -O3
# -O3 generates correct code but takes forever
FNOOPTS = -O2
MOD_OPT = -module 
INCFLAG = -I

CC_COMP = g++
C_COMP  = gcc
C_LINK  = g++
#C_COMP  = openCC -m64
#C_LINK  = openCC -m64
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/lib/ -lfftw
FFTWINCLUDE  = /usr/share/doc/fftw-dev/fortran
LAPACKLIB    = -L/usr/lib/ -llapack
