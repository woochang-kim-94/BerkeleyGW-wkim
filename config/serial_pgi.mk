# arch.mk for BerkeleyGW codes
#
# Suitable for generic single-processor Linux 
# Needs fftw-dev and installation of PGI compiler suite
#
# D. Strubbe
# September 2010, UCB

COMPFLAG  = -DPGI
MATHFLAG  = -DUNPACKED
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

FCPP    = cpp -C
F90free = pgf90 -Mfree
LINK    = pgf90
FOPTS   = -fast
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

CC_COMP = pgCC
C_COMP  = pgcc
C_LINK  = pgCC
C_OPTS  = -fast
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L/usr/lib/ -lfftw
FFTWINCLUDE  = /usr/share/doc/fftw-dev/fortran
LAPACKLIB    = -lacml
#LAPACKLIB    = -L/usr/lib/ -llapack
