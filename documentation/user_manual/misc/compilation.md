# Compiling BerkeleyGW


## Prerequisites

We have tested BerkeleyGW extensively with various configurations, and support
the following compilers and libraries:

* Operating systems: Linux, AIX, MacOS
* Fortran compilers (required, and must support Fortran 2003): `pgf90`,
  `ifort`, `gfortran`, `g95`, `openf90`, `sunf90`, `pathf90`, `crayftn`, `af90`
(Absoft), `nagfor`, `xlf90`
* C compilers (required): `pgcc`, `icc`, `gcc`, `opencc`, `pathcc`, `craycc`,
  `clang`
* C++ compilers (required): `pgCC`, `icpc`, `g++`, `openCC`, `pathCC`,
  `crayCC`, `clang++`
* MPI implementation (optional): OpenMPI, MPICH1, MPICH2, MVAPICH2, Intel MPI
* OpenMP threading (optional)
* LAPACK/BLAS implementation (required): NetLib, ATLAS, Intel MKL, ACML, Cray
  LibSci
* ScaLAPACK/BLACS implementation (required if MPI is used): NetLib, Cray
  LibSci, Intel MKL, AMD
* FFTW (required): versions 3.3.x (preferred) or 2.x.


## Compilation instructions

The instructions below represent a typical case for compiling BerkeleyGW:

1. Architecture-specific Makefile-include files appropriate for various
   supercomputers as well as for using standard Ubuntu or Macports packages are
   provided in the `config/` directory. Copy the configuration file that most
   resembles your architecture and name it `arch.mk`. Example:
   ```shell
   cp config/stampede.tacc.utexas.edu_threaded_hdf5.mk arch.mk
   ```

2. Edit the `arch.mk` file to fit your needs. Refer to the [compilation flags
   documentation](compilation-flags.md) for details.

3. Type `make` to see the available compilation targets. If the `arch.mk` is
   properly configured, you only need to type:
   ```shell
   make all-flavors -j 16
   ```
   where `16` can be substituted by the number of concurrent build processes.

3. If there is any error in the make processes, it may be beneficial to first
   build the real version of the code. You can repeatedly type `make real` and
   fix any issues in the `arch.mk`.  Once the code builds successfully, clean
   the code and make the complex version of the code:
   ```shell
   make clean-flavored
   make cplx -j 16
   ```

   Note that the complex version may always be used. Real may be used for
   systems that have both inversion (about the origin) and time-reversal
   symmetry, and will be faster and use less memory.

5. Test your build. See testsuite/README for more info.


## BerkeleyGW library

BerkeleyGW also provides a library that can be used to facilitate writing
mean-field quantities. Some external DFT codes, such as PARATEC, use this
library to directly write wavefunction files in BerkeleyGW format.

To build the BerkeleyGW library, type `make library` to create `libBGW_wfn.a`
and `wfn_rho_vxc_io_m.mod` (and dependent modules), and then compile with `-I
library/` and link `library/libBGW_wfn.a` with the other code.

An `m4` macro for configure scripts is provided in the `library` directory, for
use in linking to this library. Codes linking the library should `use` the
module `wfn_rho_vxc_io_m`.

To generate real wavefunctions, a Gram-Schmidt process should be used.
The appropriate parallelization scheme etc. will be dependent on the code,
and cannot be easily handled in a general manner here, but examples can
be found in `MeanField/SIESTA/real.f90` and `MeanField/ESPRESSO/pw2bgw.f90`
(routine `real_wfng`).
