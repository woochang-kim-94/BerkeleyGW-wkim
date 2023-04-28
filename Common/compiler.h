! Created Sept 2011 by DAS. 
! Define characteristics of various compilers, via compiler symbols (e.g. -DGNU) 
! to be used directly from the arch.mk files, and then defining what we need to do 
! for that compiler via the symbols for various properties (e.g. NOSIZEOF). 
! Ideally, to support a new compiler, one need only change this file, adding a 
! new block to define what -DNEWCOMPILER would mean. 
! NOTE: of course, Makefile-level issues still need to be handled in common-rules.mk 

#ifdef PGI
#  define COMPILER_STR "PGI"
#endif

#ifdef INTEL
#  define SIZEOF64
#  define COMPILER_STR "INTEL"
#endif

! very ancient version may require NOSIZEOF 
#ifdef GNU
#  define SIZEOF64
#  define COMPILER_STR "GNU"
#endif

! FHJ: Support for Open64 will be removed shortly in favor of OpenUH
! open64 is very similar to path, it is an open-sourced version of it
! omp_lib.f90 needed to do OpenMP, see common-rules.mk.
#ifdef OPEN64
#  define NOSIZEOF
#  define NO_OPENMP_MOD
#  define COMPILER_STR "OPEN64"
#endif
#ifdef XLF
#  define COMPILER_STR "XLF" 
#endif 
 
#ifdef ORACLE
#  define NOSIZEOF 
#  define COMPILER_STR "ORACLE" 
#endif 
 
#ifdef NAG 
#  define NOSIZEOF 
#  define COMPILER_STR "NAG" 
#endif 
 
#ifdef ABSOFT 
#  define COMPILER_STR "ABSOFT"
#endif

! cce 7.4.4 and before support sizeof for intrinsic types, but need NOSIZEOF_TYPE
! cce 8.0.0 and later do not allow sizeof for multidimensional arrays, requiring us
! to turn sizeof off everywhere. Why would Cray do this?
#ifdef CRAY
#  define NOSIZEOF
#  define COMPILER_STR "CRAY"
#endif

#ifndef COMPILER_STR
#  define COMPILER_STR "UNKNOWN"
#endif

! It is considered a bug in OPEN64 that sizeof will not work in our code. 
#if defined NOSIZEOF
#  define SIZEOF(x) 1
#elif defined SIZEOFTRICK
#  define SIZEOF(x) int(size(x),8)*int(sizeof(reshape(x,(/ 1 /))),8)
#else
#  define SIZEOF(x) sizeof(x)
#endif

! on some platforms there is a different return value for sizeof if build is 64-bit 
#if defined SIZEOF64 && defined _LP64 
#  define INTSIZEOF integer*8
#elif defined SIZEOFTRICK
#  define INTSIZEOF integer*8
#else
#  define INTSIZEOF integer
#endif

! Intrinsic module for OpenMP. Almost all compilers that support OpenMP provide
! a "omp_lib.mod" module, though the OpenMP standard allow them to only ship a
! "omp_lib.h" Fortran header.
#ifdef OMP
#  ifdef NO_OPENMP_MOD
#    define USEOMPLIB
#  else
#    define USEOMPLIB use omp_lib
#  endif
#else
#  define USEOMPLIB
#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
