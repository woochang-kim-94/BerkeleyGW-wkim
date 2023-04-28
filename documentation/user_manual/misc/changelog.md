# Changelog


## Berkeley 2.1 (Jul/2019)
BerkeleyGW 2.1 is the first version based on a new and more modern coding
infrastructure, with better support for new compilers and improved consistency
checks. For the end users, the main noticeable features include increased
performance and bug fixes for I/O operations involving the new HDF5 file format
and for the subspace code. We also improved considerably the documentation of
the code, with a new and expanded user manual.

Some highlights and features for the end user:

- New user manual for the code, which comes bundled with the code, and which is
  also available online: <http://manual.berkeleygw.org>

- New wrapper for the [StochasticGW](http://stochasticgw.com) code.

- Bug fixes when writing HDF5 wavefunctions in parallel which could cause the
  code to hang (relevant to the ParaBands code).

- Improved error checking for operations involving HDF5 files.

- Improved support for compilers, including PGI and NAG.

- Improved performance and stability for building BerkeleyGW in parallel with
  new dependency system.

- Improved performance of HDF5 routines for the subspace code and for reading
  `chimat.h5` files.


## Berkeley 2.0 (May/2018)
BerkeleyGW 2.0 represents the culmination of nearly two years of development
effort, and this release contains a number of important new features and
capabilities including:

Main new features:

1. Initial release of ParaBands: a new tool for efficiently generating
wave-function files including many empty orbitals required for BerkeleyGW
calculations.

2. Full BSE calculations that do not employ the Tamm-Dancoff approximation.

3. Improved algorithms for k-point sampling in 2D, which include the newly
proposed nonuniform neck subsampling (NNS) and the cluster sampling
interpolation (CSI) algorithms.

4. Accelerated full-frequency GW calculations through the use of a low-rank
subspace approximation for expressing the dielectric matrix. In fact,
large-scale full-frequency GW calculations are now faster than calculations
using plasmon-pole models!

5. Significant performance improvements throughout, but particularly in the
calculation of the full-frequency dielectric matrix and evaluation of the
full-frequency Sigma operator. Continued optimizations were made throughout the
package for multi- and many-core architectures including Intel Xeon-Phi, which
allows BerkeleyGW to scale half a million cores on Cori 2 for large-scale
calculations!

6. Improved user and developer documentation, as well as a new quick reference
guide (see the link on the top of the page).


## Berkeley 1.2 (Aug/2016)
F. H. da Jornada, J. Deslippe, D. Vigil-Fowler, J. I. Mustafa, T. Rangel,
F. Bruneval, F. Liu, D. Y. Qiu, D. A. Strubbe, G. Samsonidze, J. Lischner.

Features marked with [*] change the default behavior of the code relative to
BGW-1.1 and may cause a small change in the numbers produced by the code.


### New features

1. Added new methods to deal with the frequency-dependence of the
   polarizability:

   1. Added Godby-Needs (GN) plasmon-pole model.
   2. Added spectral method for real-axis (RA) full-frequency (FF)
      calculations.  The spectral method allows one to compute the dielectric
      matrix much faster for many frequency points on the real axis. The old
      method for RA-FF calculations is now refered to as the Adler-Wiser
      method. Still, RA is no longer the default FF method (see next item).
   3. Added Contour-deformation (CD) full-frequency (FF) formalism.  The CD is
      the recommended (and default) scheme for calculations performing an
      explicit evaluation of the dynamical effects of the dielectric matrix
      without plasmon-pole models. It requires the evaluation of the dielectric
      matrix on both real and imaginary frequencies, but with typically a much
      smaller number of frequencies. The previous formalism is now refered to
      as real-axis (RA) full-frequency (FF) formalism. CD is now the default
      method for FF calculations.

2. Released non-linear optics post-processing utility: This codes allows one to
   compute the inter-exciton transitions in non-linear-optics experiments.

3. Improved usability in the code output and added dynamic verbosity switch:
   BerkeleyGW has now a much simpler and neater output. The code gives time
   estimates for most tasks it performs -- at least for those that are more
   time-consuming.  There is also a run-time flag supported by all codes to
   switch the amount of verbosity the code produces. The old compilation flag
   `-DVERBOSE` is now deprecated.

4. Improved usability in the code input: Several input parameters, such as band
   occupations, are automatically detected. If not set, parameters such as the
   cutoff for the bare and screened Coulomb potential are also read from the
   wave function and dielectric matrix, respectively.

5. Automatic solution of Dyson's equation in FF calculations:
   BerkeleyGW now automatically solves Dyson's equation for FF calculations. It uses
   a varying number of frequencies to find the graphical solution for the
   quasiparticle energies, and perform extrapolation (with appropriate warning
   messages) if no intersection could be found. The code will also output the files
   `eqp0.dat` and `eqp1.dat` directly, containing the off-shell and on-shell
   quasiparticle energies.

6. Added ABINIT wrapper: It is possible to interface the ABINIT code with
   BerkeleyGW. For now, the wrapper can only output complex-valued wavefunctions,
   even for systems with inversion symmetry.


### Improvements

1. Added support for Xeon-Phi Knight's Landing (KNL)-based systems.

2. Improved OpenMP scaling throughout.  Using OpenMP is now recommended for
   large-scale computations.

3. Added new and more performant HDF5 file formats for `epsmat` and `bsemat`
   matrices.  This is recommend as default for all builds. Note that files are
   generally not compatible with BerkeleyGW 1.1 and earlier. Binary `epsmat`
   files can be converted to the HDF5 format with the `epsmat_old2hdf5` utility.
   Older `epsmat.h5` files from BerkeleyGW-1.1-beta2 can be converted to the
   current format with the `epsmat_hdf5_upgrade` utility.

4. Added parallelization over frequencies in the epsilon code.  This allows for
   better scalability in full-frequency calculations. This scheme is supported
   in the old Adler-Wiser and in the new Contour-Deformation formalisms.

5. Improved the Monte-Carlo average scheme used to compute the Coulomb
   potential. [*] For bulk systems, the new default scheme is to compute the
   average of Coulomb potential v(q+G) for all q-points and G-vectors. We use a
   hybrid scheme to make this evaluation faster. While the code is slightly
   slower for small calculations, results should converge much faster with
   k-point sampling. One can still use the old defaults with the input option
   `cell_average_cutoff 1.0d-12`. Note that no change was performed for 2D, 1D
   and 0D systems, or 3D metals.

6. Rewrote k-point interpolation engine. [*] The previous k-point interpolation
   scheme used in the `haydock`, `absorption`, and `inteqp` codes used a greedy
   algorithm to search for the closest coarse-grid points around each fine-grid
   points, which lead to discontinuities of the interpolands.  The new
   interpolation engine uses QHull, which is bundled with BerkeleyGW, to
   perform a Delaunay tessellation of the coarse-grid k-points. The new
   interpolation scheme is much more robust. You can use the old interpolation
   scheme with the `greedy_interpolation` flag.

7. Generalized scheme for k-point interpolation.  The k-point interpolation is
   used in the haydock, absorption, and inteqp codes.  The new scheme is useful
   for metals and systems with valence-conduction band character mixing in the
   BSE.


### Misc

1. Changed logic to deal with invalid GPP frequency modes. [*] When the code
   performs a HL-GPP or a GN-GPP calculation and finds an invalid mode with
   frequency with $\omega_{G,G'}^2 < 0$, the code will now treat that mode
   within the static COHSEX approximation (i.e., move frequency to infinity).
   The previous behavior was to "find" a purely complex mode frequency and
   relax the causality constraint. One can switch which strategy to use with
   the `invalid_gpp_mode` input flag.

2. Many bug fixes and performance improvement.
