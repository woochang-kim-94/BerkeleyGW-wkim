# `Absorption` code overview

This code is the second half of the BSE code. It interpolates the coarse grid
electron-hole kernel onto the fine grid and then diagonalized the BSE equation.
The output is the electron-hole eigenfunctions and eigenvalues as well the
absorption spectra.

## Summary of input and output files

### Required input files

- [`absorption.inp`](absorption-keywords.md) : Input parameters.

- `WFN_fi[.h5]` : Wavefunctions in unshifted fine grid (conduction and valence for
  momentum operator, conduction for velocity operator). Information of spinor wavefunctions
  is initialized automatically.

- `WFNq_fi[.h5]` : Wavefunctions in shifted fine grid (not needed for momentum
  operator, valence for velocity operator).

- `WFN_co[.h5]` : Wavefunctions on unshifted coarse grid.  Must be the same as used
  for Kernel.

- `eps0mat[.h5]` : Must be same as used in Kernel.  The file has a `.h5`
  extension if the code is compiled with HDF5 support (for specification, see
  `epsmat.h5.spec`).

- `epsmat[.h5]` : Must be same as used in Kernel.  The file has a `.h5`
  extension if the code is compiled with HDF5 support (for specification, see
  `epsmat.h5.spec`).

- `bsedmat` : BSE matrix elements in coarse grid, direct part. This should be
  generated with Kernel code using same `WFN_co[.h5]`.

- `bsexmat` : BSE exchange matrix elements.  This should be generated with
  Kernel code using same `WFN_co[.h5]`.

- `bsemat.h5` : Includes data from both of above if compiled with HDF5 support.
  For specification, see `bsemat.h5.spec`. Note that you can use the
  `bsemat_hdf5_upgrade.py` utility to upgrade a `bsemat.h5` file generated with
  BerkeleyGW prior to 1.2 (r6974).

### Additional input

- `WFNq_co[.h5]`: Coarse grid wavefunctions for Finite Q calcualtion, shifted relative 
  `WFN_co[.h5]` by an amount [`exciton_Q_shift`](kernel-keywords.md#exciton_Q_shift).

- `eqp.dat`: A list of quasiparticle energy corrections for the bands in
  `WFN_fi[.h5]`.  Used if [`eqp_corrections`](absorption-keywords.md#eqp_corrections)
  is set in [`absorption.inp`](absorption-keywords.md).

- `eqp_q.dat`: A list of quasiparticle energy corrections for the bands in
  `WFNq_fi[.h5]`.  Used if [`eqp_corrections`](absorption-keywords.md#eqp_corrections)
  is set in [`absorption.inp`](absorption-keywords.md).

- `eqp_co.dat` : A list of quasiparticle energy corrections for the bands in
  `WFN_co[.h5]`.  Used if
  [`eqp_co_corrections`](absorption-keywords.md#eqp_co_corrections) is set in
  [`absorption.inp`](absorption-keywords.md).

### WFN files in HDF5 (.h5) format

The wavefunction files will be read in HDF5 format (i.e. from `WFN_fi.h5`, `WFNq_fi.h5`, `WFN_co.h5` and `WFNq_co.h5` files)
if the code is compiled with HDF5 support and the  [`use_wfn_hdf5`](absorption-keywords.md#use_wfn_hdf5) is set in input.

### Auxiliary files

The files below are output files from previous runs, and may be used as input
to speed up calculations:

- `dtmat` : Transformation matrices, dcc/dvv use for interpolation between
  coarse and fine grid.  This file must be consistent with your `bsedmat` and
  `bsexmat` files and corresponding coarse and fine wavefunctions. Note: the
  file format for `dtmat` was changed in BerkeleyGW 1.1.0 (r5961).

- `vmtxel` : Optical matrix elements (velocity or momentum) between single
  particle states.

- `epsdiag.dat` :  Diagonal elements of dielectric matrix on the q-grid.  Must
  be consistent with `epsmat` and `eps0mat`.

- `eigenvalues.dat` : Contains electron-hole eigenvalues and transition matrix
  elements.

### Output files

- `eigenvalues.dat`: Has eigenvalues/transition matrix elements of e-h states,
  eigenvalues in eV, matrix elements in atomic units.

- `eigenvalues_noeh.dat`: As above but for non-interacting e-h transitions,
  showing k-point, bands, and spin for each.

- `eigenvectors`: Has the excitonic wavefunctions in Bloch space: $A^S_{vck}$.

- `absorption_eh.dat`: Dielectric function and density of excitonic states.
  There are four Columns: energy (in eV), $\varepsilon_2$, $\varepsilon_1$, and
  the DOS. The DOS is normalized, $\int DOS(\omega)\,d\omega = 1$.

- `absorption_noeh.dat`: Non-interacting (RPA) dielectric function and joint
  density of states. There are four Columns: energy (in eV) , $\varepsilon_2$,
  $\varepsilon_1$, and the JDOS. The JDOS is normalized, $\int
  JDOS(\omega)\,d\omega = 1$.

- `dvmat_norm.dat`: The norms of the `dvv` overlap matrices between the valence
  band k on the fine grid and the closest k-point on the coarse grid

- `dcmat_norm.dat`: The norms of the `dcc` overlap matrices between the
  conduction band k on the fine grid and the closest k-point on the coarse grid

- `eqp.dat`, `eqp_q.dat`: Quasiparticle corrections for `WFN_fi[.h5]` and `WFNq_fi[.h5]`
  interpolated from the coarse grid if
  [`eqp_co_corrections`](absorption-keywords.md#eqp_co_corrections) is used.

- `bandstructure.dat`: Same as `eqp.dat` and `eqp_q.dat` but in a format
  suitable for plotting as a bandstructure.

## Details: Finite Q conventions

As a direct result of the convention dicussed in the kernel overview,
running a finite Q calculation with the input:

[`exciton_Q_shift`](absorption-keywords.md#exciton_Q_shift) `Qflag` `Q_x`, `Q_y`, `Q_z`

will yield excitons with:

$$ 
\Psi_{S,-Q}(r_e, r_h) = \sum_{cvk} A^{S,-Q}_{cvk} \psi_{ck}(r_e) \psi^\star_{vk+Q}(r_h)  
$$

where $Q = (Q_x, Q_y, Q_z)$. In other words, the shift provided in 
[`exciton_Q_shift`](kernel-keywords.md#exciton_Q_shift) is negative the excitons
center of mass momentum.

## Tricks and hints

- To optimize distribution of work among PEs, choose the number of PEs so that
  $n_k n_c n_v$ in the fine grid is a multiple of the number of PEs. The
  parallelization is first done over over k-points.

- Check if the transformation matrices have norm close to 1! The norms are
  written to the files `dvmat_norm.dat` and `dcmat_norm.dat`.

- Unfolding of irreducible BZ: if you want to skip the unfolding and use the
  set of k-points in `WFN_fi[.h5]` as a sampling of the whole BZ, specify the
  `no_symmetries_*` options in [`absorption.inp`](absorption-keywords.md).

- The [`number_eigenvalues`](absorption-keywords.md#number_eigenvalues) keyword:
  using this keyword tells the code to store only the first
  [`number_eigenvalues`](absorption-keywords.md#number_eigenvalues)
  eigenvalues/eigenvectors. So far, this option is implemented only with the
  ScaLAPACK diagonalization routines.

- Analyzing eigenvectors: we provide a tool called `summarize_eigenvectors` to
  read in and analyze eigenvectors for a group of exciton states. For specific
  states specified, it sums the total contribution from each k-point. Please
  see the example input `summarize_eigenvectors.inp`.
