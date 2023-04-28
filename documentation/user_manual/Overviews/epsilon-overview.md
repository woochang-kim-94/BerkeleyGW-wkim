# `Epsilon` code overview

`Epsilon` is the name of the code that generates the polarizability matrix and
the inverse dielectric matrix for a bulk or nanoscale system. The main result
of the `epsilon` code is the generation of the `epsmat[.h5]` and `eps0mat[.h5]`
files which can be used in a Sigma or BSE calculation.

## Summary of input and output files

### Required input files

- [`epsilon.inp`](epsilon-keywords.md): Input parameters.

- `WFN`: This is linked to the unshifted grid. Information of spinor wavefunctions
  is initialized automatically.

- `WFNq`: This is linked to the shifted grid.  When a shift is not required,
  this can be the same file as `WFN`. This happens for calculations on:
    - semiconductors with spherical or box truncation;
    - metals / graphene-like screening with the `is_q0` flag set to 2.

### Auxiliary input files

These files are output files from previous runs, and may be used as input to speed up calculation:

- `chimat[.h5]`: Polarizability matrix for $q\ne0$. No need to recalculate
  matrix elements.  The file has a `.h5` extension if the code is compiled with
  HDF5 support (for specification, see `epsmat.h5.spec`).

- `chi0mat[.h5]`: Polarizability matrix for $q\rightarrow0$. No need to
  recalculate matrix elements.  The file has a `.h5` extension if the code is
  compiled with HDF5 support (for specification, see `epsmat.h5.spec`).

The files below are used if
[`eqp_corrections`](epsilon-keywords.md#eqp_corrections) is set in
[`epsilon.inp`](epsilon-keywords.md). The corrected eigenvalues are used for
constructing the polarizability matrix.

- `eqp.dat`: A list of quasiparticle energy corrections for the bands in `WFN`.
- `eqp_q.dat`: A list of quasiparticle energy corrections for the bands in
  `WFNq`.

### Output Files

- `espmat[.h5]`: Inverse dielectric matrix for $q\ne0$.  The file has a `.h5`
  extension if the code is compiled with HDF5 support (for specification, see
  `epsmat.h5.spec`).

- `eps0mat[.h5]`: Inverse dielectric matrix for $q\rightarrow0$.  The file has
  a `.h5` extension if the code is compiled with HDF5 support (for
  specification, see `epsmat.h5.spec`).

- `epsilon.log`: The log file containing values of chimat and epsmat.

- `chi_converge.dat`: Convergence chi with respect to empty orbitals. Columns:

    - number of conduction bands;

    - $\mathrm{Re}\, \chi(G=G'=0,q)$, extrapolated
      $\mathrm{Re}\,\chi(G=G'=0,q)$;

    - $\mathrm{Re}\, \chi(G=G'=G_\mathrm{max},q)$, extrapolated
      $\mathrm{Re}\,\chi(G=G'=G_\mathrm{max},q)$.
