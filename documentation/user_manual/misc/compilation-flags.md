# Compilation flags in `arch.mk`

Below, we describe the flags that can be defined in the machine-dependent
`arch.mk` file necessary to compile BerkeleyGW. Also, refer to the examples in
the `config/` directory.

!!! Note "Reuse first!"
    It is likely that you can simply copy an `arch.mk` in the `config/` folder
    that will work on your platform without any modification. In the worst case
    scenario, pick a template `arch.mk` from the `config/` folder that most
    closely matches your architecture, and customize it as necessary instead of
    creating a new `arch.mk` from scratch!


## General options

These flags control which compiler, parallelization scheme, and external
libraries are available, and hence dramatically affect the performance of the
code and which components of the code are enabled.

#### `COMPFLAG`: Fortran compiler selector
Use this to select the Fortran compiler. Select only one of the following
options:
```makefile
COMPFLAG = -D{INTEL, PGI, GNU, XLF, ABSOFT, NAG, OPEN64, ORACLE, CRAY}
```

If you have a different one than these, modifications in `Common/compiler.h`
and `Common/common-rules.mk` may be be necessary. The BerkeleyGW developers'
team can assist you in this case.

----

#### `PARAFLAG`: Parallelization flags
More than one option is possible. Possible flags are:

- `-DMPI`: compiles BerkeleyGW with support for MPI distributed-memory
  parallelization.
- `-DOMP`: compiles BerkeleyGW with support for OpenMP shared-memory
  parallelization. In this case, add appropriate compiler-dependent flag to
  `F90free` and `LINK` (documented below):
    - `ifort`, `sunf90`, `absoft`, `pathscale`, `open64`: `-openmp`.
    - `pgi`: `-mp=nonuma`.
    - `gfortran`: `-fopenmp`.
    - `xlf`: `-qsmp=omp`.
    - `crayftn`: on by default (to turn off use `-h noomp`).
    - `nag`: add `-openmp`. You may also need to add `-Wl,-fopenmp` to `LINK`.
    - `g95` does not support OpenMP.

Typical compilations will include both MPI and OpenMP parallelization, with
```makefile
PARAFLAG = -DMPI -DOMP
```

----

#### `MATHFLAG`: Mathematical libraries
More than open option is possible. Possible flags are:

- `-DUSESCALAPACK`: enables usage of ScaLAPACK (required for MPI parallel
  builds).
- `-DUSEESSL`: uses IBM's ESSL library instead of LAPACK in some parts of the
  code.
- `-DUNPACKED`: uses unpacked rather than packed representation of the
  Hamiltonian in EPM. Packed LAPACK operations give bad results eventually in
  Si-EPM kernel 'max x' test for ACML, Cray LibSci, and MKL sometimes.
- `-DUSEFFTW3`: uses FFTW3 instead of FFTW2, which is recommended for better
  performance, and enables threading.
- `-DHDF5`: enables usage of HDF5 for writing `epsmat` and `bsemat` files. The
  produced files will be `eps0mat.h5`, `epsmat.h5` and `bsemat.h5` and are used
  in place of `eps0mat`, `epsmat`, `bsedmat`, `bsexmat` (the latter two are
  combined in one `.h5` file).  Using `-DHDF5` gives you substantially better
  IO performance in many cases by taking advantage of HDF5's parallel IO
  options and other options. However, you must have the HDF5 libraries
  installed on your system and define `HDF5LIB` and `HDF5INCLUDE` below.
- `-DUSEPRIMME`: enables support for the PRIMME library for iterative
  diagonalization. Currently only used in `ParaBands`.
- `-DUSEMR3`: enables support for the MR3 solver in ScaLAPACK. Used in BSE for
  solving the full BSE without the Tamm-Dancoff approximation. Most moderns
  ScaLAPACK implementations, including the Netlib reference implementation,
  support MR3.
- `-DUSEVORO`: enable support for voro++ program to compute Voronoi regions.

Typical compilations will use
```makefile
MATHFLAG = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5
```

----

#### `DEBUGFLAG` and `C_DEBUGFLAG`: Debugging flags
More than open option is possible. Possible flags are:

- `-DDEBUG`: enables extra memory checking in the code. This typically slows down the code by up to ~20%.
- `-DVERBOSE`: equivalent to running all codes with the maximum verbosity level.

`DEBUGFLAG` and `C_DEBUGFLAG` are used separately by Fortran and C/C++ programs, respectively.
Don't set any of those options unless a developer told you to, as they will slow down
the code and make the output much more verbose (almost unreadable).

Typical compilations will not set `DEBUGFLAG` nor `C_DEBUGFLAG`.

----

## Fortran 2003 compiler options


#### `FCPP`: Fortran preprocessor
The best option is to use gcc's c preprocessor (often in `/bin/cpp`), or a C
compiler (even `mpicc`) with the `-E` and/or `-x c` option. The options `-C`
(and not `-ansi`) should generally be used.  Add `-P` if the Fortran compiler
will not accept indication of line numbers. You may also need to add
`-nostdinc` and `-nostdinc++` if your Fortran compiler complains about lines
such as `#include ...`.

A typical usage on Linux systems is:
```makefile
FCPP = cpp -C -nostdinc -nostdinc++
```

----

#### `F90free`: Fortran compiler
Fortran 2003 compiler and its flags for F90 and f90 files (free source form).
Typical syntax is:
```makefile
F90free = {mpif90, mpiifort, ftn, ifort, pgf90, gfortran, ...} [flags]
```

See note in the section on `PARAFLAG` regarding flags for OpenMP.  Sometimes,
flags such as `-fno-second-underscore` or `-assume 2underscores` may be
necessary to make linking succeed against C libraries such as FFTW, depending
on what compiler built them.

A typical usage with the `gfortran` compiler is:
```makefile
F90free = mpif90 -ffree-form -ffree-line-length-none -fno-second-underscore -fopenmp
```

----

#### `LINK`: Fortran linker
Compiler and any flag(s) required for *linking*. Typical syntax is:
```makefile
LINK = {mpif90, mpiifort, ftn, ifort, pgf90, gfortran, ...} [flags]
```

A typical usage with the `gfortran` compiler is:
```makefile
LINK = mpif90 -fopenmp
```

----

#### `FOPTS`: Fortran optimization flags
Optimization flags depend on the compiler. A typical usage with the `gfortran`
compiler is:
```makefile
FOPTS = -O3
```

A typical usage with the `gfortran` compiler to debug the code is:
```makefile
FOPTS = -g -C -O0 -fbounds-check -fbacktrace -finit-real=snan -Wall -ffpe-trap=invalid,zero,overflow
```

----

#### `FNOOPTS`: Fortran conservative optimization flags
These are used for compiling the `epsilon.f90` and `sigma.f90` codes, which can
have trouble in high optimization due to their long length.

If below does not work,
```makefile
FNOOPTS = $(FOPTS)
```

try:
```makefile
FNOOPTS = -O2
```

----

#### `MOD_OPT`: Output path of Fortran modules.
Fortran compiler flag to specify location where module should be created.  This
is different for each compiler, so please refer to the examples. A typical
example for `gfortran` is:
```makefile
MOD_OPT = -J
```

----

#### `INCFLAG`: Search path for Fortran modules.
Compiler flag to specify where to look for modules. This flag is compiler
dependent. Example for `gfortran`:
```makefile
INCFLAG = -I
```

----

## C/C++ compiler options

  
#### `C_PARAFLAG`: Parallel flag for C++ programs.
Use this to compile C++ in parallel with MPI with the option:
```makefile
C_PARAFLAG = -DPARA
```

Leave blank for serial. MPICH, MVAPICH, Intel MPI require also
`-DMPICH_IGNORE_CXX_SEEK`, or an error may occur such as: `SEEK_SET is #defined
but must not be for the C++ binding of MPI. Include mpi.h before stdio.h`

----

#### `C_COMP`: C compiler
Include any needed flags. Typical syntax is:
```makefile
C_COMP = {mpicc, mpiicc, cc, icc, pgcc, gcc, ...} [options]
```

A typical usage with the `gcc` compiler suite is:
```makefile
C_COMP = mpiccc -fopenmp
```

----

#### `CC_COMP`: C++ compiler
Include any needed flags. Typical syntax is:
```makefile
CC_COMP = {mpiCC, mpicxx, mpiicpc, CC, icpc, pgCC, g++, ...} [options]
```

A typical usage with the `gcc` compiler suite is:
```makefile
CC_COMP = mpicxx -fopenmp
```

----

#### `C_LINK`: C++ linker
Include any needed flags. Typical syntax is:
```makefile
CC_COMP = {mpiCC, mpicxx, mpiicpc, CC, icpc, pgCC, g++, ...} [options]
```

A typical usage with the `gcc` compiler suite is:
```makefile
CC_COMP = mpicxx -fopenmp
```

----

#### `C_OPTS`: C/C++ optimization flags
Typical usage with the `gcc` compiler suite:
```makefile
C_OPTS = -O3
```

----

## Math Libraries

Static linking is highly recommended when running large-scale parallel
calculations. Note that sometime we define variables such as `FFTWPATH` and
`HDF5PATH`, which are not required by BerkeleyGW but which make it easier to
unify some paths, such as `FFTWLIB` and `FFTWINCLUDE`.


#### `FFTWLIB`: FFTW linking
Link line for FFTW2 or FFTW3 library. Sometimes needs to be linked with `-ldfftw`, and sometimes with `-lfftw`.

Example for FFTW2 or non-threaded FFTW3:
```makefile
FFTWLIB = $(FFTWPATH)/lib/libfftw.a # static linking, recommended
FFTWLIB = -L$(FFTWPATH)/lib -lfftw # dynamic linking

```

Example for threaded FFTW3 (i.e., if `-DUSEFFTW3` is used in `MATHFLAGS`).
```makefile
FFTWLIB = $(FFTWPATH)/lib/libfftw3_omp.a $(FFTWPATH)/lib/libfftw3.a # static linking, recommended
FFTWLIB = $(FFTWPATH)/lib -lfftw3_omp -lfftw3 # dynamic linking
```

Examples for MKL:
```makefile
FFTWLIB = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
               -Wl,--start-group \
               $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
               $(MKLROOT)/lib/intel64/libmkl_core.a \
               $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
               $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
               -Wl,--end-group -lpthread -lm -ldl
```
and:
```makefile
FFTWLIB = -Wl,--start-group \
                $(MKLPATH)/libmkl_intel_lp64.a \
                $(MKLPATH)/libmkl_intel_thread.a \
                $(MKLPATH)/libmkl_core.a \
                -Wl,--end-group -liomp5 -lpthread -lm -ldl
```

!!! Note
    To use FFTW3 from MKL, do not link the `fftw3xf` interfaces.

----

#### `FFTWINCLUDE`: FFTW headers

FFTW2 requires include file `fftw_f77.i`; sometimes supercomputer installations
will not have it. If not, [find it
online](http://www.fifi.org/doc/fftw-dev/fortran/fftw_f77.i) and copy it. FFTW3
requires include file `fftw3.f` instead.

Typical syntax:
```makefile
FFTWINCLUDE = $(FFTWPATH)/include
```

Example for MKL:
```makefile
FFTWINCLUDE = $(MKLROOT)/include/fftw
```

----

#### `LAPACKLIB`: LAPACK linking

Different styles will apply depending on the package used.
Must provide BLAS as well as LAPACK: some packages will do this in one library file,
others will not. See examples in this directory for how to link with MKL, ACML,
ATLAS, Ubuntu packages, etc. For MKL, if you will use ScaLAPACK, it is most
convenient to include BLACS here too. For further info on linking with MKL, see
the [link-line advisor](http://software.intel.com/sites/products/mkl/).

Typical generic syntax on Ubuntu:
```makefile
LAPACKLIB = -L/usr/lib/ -llapack -lblas
```

Another syntax on Ubuntu, statically linking against Netlib LAPACK and
OpenBLAS:
```makefile
LAPACKLIB = /lib64/liblapack.a /opt/OpenBLAS-0.2.20/lib/libopenblas.a
```

Typical syntax on MacOS, statically linking against Veclib:
```makefile
LAPACKLIB = /opt/local/lib/libvecLibFort.a
```

Typical example for MKL:
```makefile
LAPACKLIB = -Wl,--start-group \
                $(MKLPATH)/libmkl_intel_lp64.a \
                $(MKLPATH)/libmkl_intel_thread.a \
                $(MKLPATH)/libmkl_core.a \
                $(MKLPATH)/libmkl_blacs_intelmpi_lp64.a \
                -Wl,--end-group -liomp5 -lpthread -lm -ldl
```

----

#### `SCALAPACKLIB`: ScaLAPACK linking
Specify if you have `-DUSESCALAPACK` in `MATHFLAG`. BLACS must be provided here
too, if it is not provided in LAPACKLIB above. See examples in this directory
for different packages, as with LAPACKLIB.


Typical example for Ubuntu:
```makefile
SCALAPACKLIB = /lib64/libscalapack.a
```

Typical example for MacOS:
```makefile
SCALAPACKLIB = -L/opt/local/lib/ -lscalapack
```

Typical example for MKL:
```makefile
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a
```

----

#### `HDF5LIB`: HDF5 linking

Only necessary if you have `-DHDF5` in `MATHFLAG`. If compiling in parallel,
you must have HDF5 compiled with parallel support or linking will fail. Path as
below should generally work. Sometimes `-lsz` too, or static linkage, might be
required. Note that `-lz` is not from HDF5, but its dependency zlib, and
`libdl` may be required on your system even when linking static HDF5 due to
plugin support in 1.8.11 and higher.

Example usage:
```makefile
HDF5LIB = $(HDF5PATH)/libhdf5hl_fortran.a \
          $(HDF5PATH)/libhdf5_hl.a \
          $(HDF5PATH)/libhdf5_fortran.a \
          $(HDF5PATH)/libhdf5.a \
          $(HDF5PATH)/libsz.a \
          -lz
```
and:
```makefile
HDF5LIB = -L$(HDF5PATH)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
```

----

#### `HDF5INCLUDE`: HDF5 headers

See notes on `HDF5LIB`. Typical usage:
```makefile
HDF5INCLUDE = $(HDF5PATH)/include
```

----


## Advanced options
These flags don't need to be typically set, and the default values work on most
systems.


#### `CC_LINK_LIBS`: Additional C++ libraries
Any flag or library that you need to pass to your Fortran compiler to link
object files compiled from C++. If this field is blank, we use a set of
compiler-dependent default values, hard coded in `Common/common-rules.mk`.
This flag is currently only being used if you compile your code with the
`-DUSEVORO` flag.

Normally you don't need to set `CC_LINK_LIBS`. In some rare cases, you may need
to use:
```makefile
CC_LINK_LIBS = -lstdc++
```

----

#### `REMOVE`: Remove files
Default is:
```makefile
REMOVE = /bin/rm -f
```

----

#### `AR`: Make archives
Command for making archives `(.a)`. Default is:
```makefile
AR = /usr/bin/ar
```


Usually that is fine and `AR` does not need to be specified, but
to do interprocedural optimizations (IPO) with ifort, you need `xiar` instead:
```makefile
AR = /usr/bin/env xiar
```

----

#### `PERFORMANCE`: Performance/profiling

To be used with an external tool such as IPM on NERSC. Default is empty:
```makefile
PERFORMANCE = 
```

----

#### `TESTSCRIPT`: Testsuite job script

Command to submit testsuite job script from testsuite directory.
Usage of `qsub` is for PBS. If no scheduler, put `make check-parallel` here.
In serial, delete this line.

Example:
```makefile
TESTSCRIPT = qsub {architecture}.scr
```
