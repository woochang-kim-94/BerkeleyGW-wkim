# `Kernel` code overview

This code constructs the direct and exchange Kernel matrix on the coarse grid.
This is done essentially by computing Eqs 34, 35 and 42-46 of [Rohlfing and
Louie, PRB 62, 4927 (2000)](https://doi.org/10.1103/PhysRevB.62.4927).

## Summary of input and output files

### Required input files

- [`kernel.inp`](kernel-keywords.md): Input parameters.

- `WFN_co[.h5]`: Wavefunctions on coarse grid. Recommended: use an unshifted grid of
  the same size as the q-grid in `epsmat[.h5]`.  Shift will increase number of
  q-vectors needed in `epsmat[.h5]`. The file will be read in HDF5 format (i.e. from file
  with `.h5` extension) if the code is compiled with HDF5 support and the
  [`use_wfn_hdf5`](kernel-keywords.md#use_wfn_hdf5) is set in input. Information of spinor
  wavefunctions is initialized automatically.

- `epsmat[.h5]`: Inverse dielectric matrix ($q\ne 0$), created using
  [`epsilon`](epsilon-overview.md).  Must contain the all $q=k-k'$ generated from
  `WFN_co[.h5]`, including with symmetry if
  [`use_symmetries_coarse_grid`](kernel-keywords.md#use_symmetries_corase_grid)
  is set. The file has a `.h5` extension if the code was built with HDF5
  support.

- `eps0mat[.h5]`: Inverse dielectric matrix ($q\rightarrow 0$), created using
  [`epsilon`](epsilon-overview.md) The file has a `.h5` extension if the code was
  built with HDF5 support. Note Kernel does not require quasiparticle
  eigenvalues.  It may be run in parallel with [`sigma`](sigma-overview.md)

### Additional input

- `WFNq_co[.h5]`: Coarse-grid wavefunctions for finite-Q calculations. 
   The file will be read in HDF5 format (i.e. from file with `.h5` extension) 
   if the code is compiled with HDF5 support and the
  [`use_wfn_hdf5`](kernel-keywords.md#use_wfn_hdf5) is set in input.


### Output files

- `bsedmat`: Direct kernel matrix elements on unshifted coarse grid.

- `bsexmat`: Exchange kernel matrix elements on unshifted coarse grid.

- `bsemat.h5`: Includes data from both of above if compiled with HDF5 support.
  For specification, see `bsemat.h5.spec`.


## Details: wings of epsilon

For semiconductors, the wings of $\chi$ have terms of the following form:

$$
\langle vk | e^{i(G+q)r} | ck+q \rangle\langle ck+q | e^{-iqr} | vk \rangle
$$

The matrix element on the right is $\langle u_{ck+q} | u_{vk} \rangle$ where
$u$ is the periodic part of the Bloch function. From k.p perturbation theory,
this matrix element is proportional to q. The matrix element on the left with a
non-zero G is typically roughly a constant as a function of q for small q (q
being a small addition to G).

Thus for a general G-vector, $\chi_\mathrm{wing}(G,q) \propto q$. This directly
leads to the wings of the screened untruncated Coulomb interaction being
proportional to 1/q. Note that this function changes sign as $q \rightarrow
-q$. Thus, when treating the $q=0$ point, we set the value of the wing to zero
(the average of its value in the mini-Brillouin zone (mBZ).

For G-vectors on high-symmetry lines, some of the matrix elements on the left
of the equation above will be zero for $q=0$, and therefore proportional to q.
For such cases, $\chi_\mathrm{wing}(G,q) \propto q^2$, and the wings of the
screened Coulomb interaction are constant as a function of q. However, setting
the $q\rightarrow 0$ wings to zero still gives us, at worst, linear convergence
to the final answer with increased k-point sampling, because the  $q\rightarrow
0$ point represents an increasingly smaller area in the BZ. Thus, we still zero
out the $q\rightarrow 0$ wings, as discussed in [A Baldereschi and E Tosatti,
Phys. Rev.
B 17, 4710 (1978)](https://doi.org/10.1103/PhysRevB.17.4710).

In the future, it may be worthwhile to have the user calculate $\chi$ /
$\epsilon$ at two q-points (a third q-point at $q=0$ is known) in order to
compute the linear and quadratic coefficients of each $\chi_\mathrm{wing}(G,q)$
so that all the correct analytic limits can be taken. This requires a lot of
messy code and more work for the user for only marginally better convergence
with respect to k-points (the wings tend to make a small difference, and this
procedure would matter for only a small set of the G-vectors).

It is important, as always, for the user to converge their calculation with
respect to both the coarse k-point grid used in sigma and kernel as well as
with the fine k-point grid in absorption.

## Details: Finite Q kernel conventions

The internal convention for setting up the finite Q kernel differ from those 
discussed in [Phys. Rev. Lett. 115, 176801
(2015)](https://doi.org/10.1103/PhysRevLett.115.176801).

Specifically when the [`exciton_Q_shift`](kernel-keywords.md#exciton_Q_shift) flag 
is used, a finite Q shift is applied to the valence state, not the conduction state. 
Explicity, with the convention $| cvkQ \rangle \equiv | ck \rangle |vk+Q \rangle$:

$$
\langle cvkQ | K^x | c'v'k'Q \rangle = \sum_G M^*_{vc}(k,Q,G) v(Q+G) M_{v'c'}(k',Q,G)
$$

$$
\langle cvkQ | K^d | c'v'k'Q \rangle = -\sum_{GG'} M^*_{cc'}(k,Q,G) W_{GG'}(Q) M_{v'v}(k+Q,q,G)
$$

where $M_{nn'}(k,Q,G) = \langle nk+Q |e^{i(Q+G) \cdot r} |n'k \rangle$

## Tricks and hints

-   To optimize distribution of work among MPI processors, choose the number of
    processors to divide:

    - $n_k^2$, if $n_\mathrm{pes} \le n_k^2$; or
    - $n_k^2 n_c^2$, if $n_\mathrm{pes} \le n_k^2 n_c^2$; or
    - $n_k^2 n_c^2 n_v^2$, if $n_\mathrm{pes} \le n_k^2 n_c^2 n_v^2$,

    where $n_\mathrm{pes}$ is the number of MPI ranks, $n_k$ is the number of
    symmetry-unfolded k-points, and $n_v$ ($n_c$) is the number of valence
    (conduction) bands included in the `kernel` calculation.

-   Converters from old versions of file formats to current version are available
    in version 2.4.
