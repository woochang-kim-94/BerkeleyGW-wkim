ABOUT PARABANDS
===============

ParaBands is a code that builds and diagonalizes the DFT Hamiltonian. Instead
of using iterative algorithms like most DFT programs, ParaBands builds the
dense DFT Hamiltonian in plane wave basis set and calls a direct solver to
obtain the requested number of empty bands. This approach is suitable when the
desired number of bands is a large fraction (more than 5%) of the size of the
Hamiltonian.

ParaBands has the following features:

  - Support for 3 different classes of solvers using:
    - LAPACK/ScaLAPACK (Bisection, Divide and Conquer, and MR3)
    - ELPA
    - PRIMME

  - Parallelization over k-point pools

  - Fully parallel HDF5 output for the wave function file

  - Hybrid MPI and OpenMP parallelism

Note that you MUST compile BerkeleyGW with HDF5 support, otherwise ParaBands
will not be able to write any wave function file. Also, ParaBands currently
only works on top of calculations performed with Quantum ESPRESSO.

ParaBands was written by Felipe H. da Jornada, and it is inspired and motivated
by the SAPO code written by Georgy Samsonidze.


BASIC USAGE
===========

In order to run ParaBands, you will first need to perform a mean-field
calculation with Quantum ESPRESSO with the desired k-points, and with a small
number of bands (containing at least one fully empty band) and export the
following quantities:
- WFN_in: the wave function file (we need this to get the el. occupations);
- VKB: the self-consistent potential in G-space; and
- VSC: the Kleinman-Bylander non-local projectors in G-space.
Note that only QE can export the VSC and VKB files right now.

If you plan on doing self-energy calculations, you should also export the
exchange-correlation potential in G-space (VXC), because degenerate wave
functions generated by ParaBands will likely differ from those generated by the
mean-field program, so vxc.dat = <nk|VXC|nk> will not be consistent.

The recommended input for pw2bgw.inp has the following structure:

```
&input_pw2bgw
  prefix = 'SYSTEM_NAME'
  wfng_flag = .true.
  wfng_file = 'WFN_in'
  real_or_complex = 2
  wfng_kgrid = .true.
  wfng_nk1 = 6
  wfng_nk2 = 6
  wfng_nk3 = 6
  wfng_dk1 = 0.0
  wfng_dk2 = 0.0
  wfng_dk3 = 0.0
  rhog_flag = .false.
  vxc_flag = .false.
  vxcg_flag = .true.
  vxcg_file = 'VXC'
  vscg_flag = .true.
  vscg_file = 'VSC'
  vkbg_flag = .true.
  vkbg_file = 'VKB'
/
```

The documentation for the ParaBands input file is in parabands.inp. Note that
only the input and output file names are actually required to run ParaBands.

Note that the output wavefunction file is written in HDF5 format. However, only
Fortran binary I/O of wavefunctions is supported in the rest of BerkeleyGW. You
will need to use the utility `hdf2wfn.x` to convert to Fortran binary format.


FAQs
====

Q) How many atoms can ParaBands handle?
A) ParaBands does not care (much) about the number of atoms. The bottleneck is
the size of the Hamiltonian in G-space.

Q) Ok, what is the maximum size of the Hamiltonian ParaBands can handle?
A) This is equivalent to asking "how long does it take to diagonalize a matrix
of size NxN". Here are some benchmarks on the Cori2 supercomputer at NERSC,
using the ELPA solver, for a single k-point, and to obtain all bands:

```
  |      N | Nodes |  Time  |
  |--------|-------|--------|
  |  17000 |     4 |   10 s |
  |  60000 |    16 | 10 min |
  | 100000 |    64 | 15 min |
  | 150000 |   256 | 17 min |
  | 200000 |   512 | 25 min |
  | 250000 |   512 | 38 min |
  | 320000 |   512 | 68 min |
  | 416000 |  2560 | 80 min |
```

Note that I/O time is not included.

Q) How do I determine the size of the Hamiltonian before I run ParaBands?
A) The easiest way is to run an scf calculation with QE, and look for the part
on "Parallelization Info" at the very top. Example:

```
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         282     282     70                19098    19098    2384
     Max         284     284     74                19112    19112    2396
     Sum       18167   18167   4543              1222827  1222827  152845
```

If your calculation has a single k-point, the row "Sum", column "PW" has the
number of G-vectors to represent the wavefunction (152845), which is the size
of the Hamiltonian constructed by ParaBands. If you have more k-points in your
calculation, that number contains all possible G-vectors needed to represent
any wavefunction for all k-points, which is typically ~10% larger than the
number of G-vectors at any
given k-point.

Q) Are there any other caveats?
A) ParaBands sometimes doesn't work with the real version. If you see a warning
that the difference between input and output energies is large, try to run
ParaBands again with the complex version of the code and complex input
wavefunctions.

PERFORMANCE TIPS
================

  - The best solver algorithm for ParaBands is ELPA. Make sure you compile
    BerkeleyGW with support for ELPA (see section below).

  - If you are using a parallel file system, make sure to stripe the output
    directory before running ParaBands.

  - Out of the ScaLAPACK solvers, Bisection is typically the best. But you
    might experience some performance gain with the Divide and Conquer
    algorithm if you ask for all bands from your DFT Hamiltonian.

  - When using a direct solver (LAPACK or ScaLAPCK, which is the default), you
    won't save much run time on the diagonalization if you ask for a smaller
    number of bands. But you may save some run time overall because of the I/O.

  - The performance of each solver depends a lot on the LAPACK/ScaLAPACK
    implementation, MPI stack, and interconnect. You are encouraged to test
    all solvers on a small system when you run ParaBands on a new architecture.

  - Use the default block_size of 32 when running ParaBands without threads.
    If you increase the number of threads, you'll need to increase the block
    size to 64 or 128 to get any performance gain.

  - Note that the memory required to run ScaLAPACK increases with the block_size!


MR3 SOLVER
==========

Not all ScaLAPACK implementations support Relatively Robust Representations (MR3).
You'll need to include the flag -DUSEMR3 in the MATHFLAGS variable in your
arch.mk file to enable this solver.


ELPA SOLVER
===========

Using the ELPA solver is **HIGHLY RECOMMENDED**. To enable it, you should:

  1) Download and compile ELPA **VERSION 2017.05.001 OR ABOVE**.

  2) Add the flag -DUSEELPA in the MATHFLAGS variable in your arch.mk file.

  3) Add a variable ELPALIB to your arch.mk which points to libelpa.a

  4) Add a variable ELPAINCLUDE to your arch.mk which points to the directory
     that contains the modules generated by ELPA. This has typically the
     structure: ELPA_PREFIX/include/elpa-VERSION/modules

You can control the kernel that ELPA uses by setting the environment variables
COMPLEX_ELPA_KERNEL and REAL_ELPA_KERNEL. Take a look at
ELPA/src/elpa2_utilities.F90 for the available kernels for your architecture. A
somewhat safe and good option is to include the following lines in your
submission script:

```
export REAL_ELPA_KERNEL=REAL_ELPA_KERNEL_SSE export
COMPLEX_ELPA_KERNEL=COMPLEX_ELPA_KERNEL_SSE
```


PRIMME SOLVER
=============

The PRIMME solver is only useful if you want to obtain a very small number of
bands, though it is a very robust solver.  If you wish to run ParaBands with
PRIMME, you will need to:

  1) Download PRIMME version 2.1.

  2) Compile PRIMME.

  3) Add the flag -DUSEPRIMME in the MATHFLAGS variable in your arch.mk file.

  4) Add a variable PRIMMELIB to your arch.mk which points to libprimme.a.


TODO
====

FHJ's wish list:

  - Memory report
  - Restart feature
  - Build VKB from PP file (almost there!)
