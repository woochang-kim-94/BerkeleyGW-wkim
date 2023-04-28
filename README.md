BerkeleyGW
==========
[![Buildbot status](https://buildbot.berkeleygw.org:884/status/master)](https://buildbot.berkeleygw.org)

BerkeleyGW is a GW / Bethe-Salpeter equation computer package.


Documentation
-------------
The user documentation is gathered in the directory `documentation/users/`.


Compilation
-----------
We have tested the code extensively with various configurations, and support the following compilers and libraries:

  * Operating systems: Linux, AIX, MacOS
  * Fortran compilers (required): pgf90, ifort, gfortran, g95, openf90, sunf90, pathf90, crayftn, af90 (Absoft), nagfor, xlf90 (experimental)
  * C compilers (required): pgcc, icc, gcc, opencc, pathcc, craycc, clang
  * C++ compilers (required): pgCC, icpc, g++, openCC, pathCC, crayCC, clang++
  * MPI implementation (optional): OpenMPI, MPICH1, MPICH2, MVAPICH2, Intel MPI
  * LAPACK/BLAS implementation (required): NetLib, ATLAS, Intel MKL, ACML, Cray LibSci
  * ScaLAPACK/BLACS implementation (required if MPI is used): NetLib, Cray LibSci, Intel MKL, AMD
  * FFTW (required): versions 3.3.x

  1. Architecture-specific Makefile-include files appropriate for various supercomputers as well as
     for using standard Ubuntu or Macports packages are provided in the config directory. Copy or
     link one to the top directory. Example:

     cp config/stampede.tacc.utexas.edu_threaded_hdf5.mk arch.mk

  2. Edit it to fit your needs. Refer to config/README for details.

  3. Copy flavor_real.mk or flavor_cplx.mk to flavor.mk to set flavor. 
     Complex may always be used. Real may be used for systems that
     have both inversion (about the origin) and time-reversal symmetry, and will be
     faster and use less memory.

  4. Stay in the root directory and run the following commands 
     to compile the various codes:
 
     MAIN:
     make epsilon
     make sigma
     make bse

     MISC:
     make plotxct
     make nonlinearoptics
     make meanfield
     make visual

     EVERYTHING:
     make all                       All codes/packages (for selected flavor)
     make -j all                    parallel build (for selected flavor)
     make all-flavors               everything, for both flavors
     make -j all-flavors            parallel build of everything, for both flavors
     make install INSTDIR=          install binaries, library, examples, testsuite into specified prefix

     These commands will generate executables with extension .x in 
     the source directory and symbolic links with the same name in 
     the bin directory.

  5. Test your build. See the user documentation for more info.


License
-------
BerkeleyGW is distributed under the Berkeley Software Distribution (BSD) license.
For license information see "license.txt" as well as "Copyright.txt"


Contributors
------------
See the documentation for the complete list of authors.
