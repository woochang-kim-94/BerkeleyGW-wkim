# `SAPO` code overview

!!! warning
    The SAPO documentation is under development!

The SAPO (Simple Approximate Physical Orbitals) code reads DFT
wavefunctions from WFN file, generates additional wavefunctions on top
of them, and writes both to another WFN file. It can generate plane waves,
decompose them into irreducible representations of the space group of the
Bravais lattice or the crystal, and orthonormalize them with respect to
DFT wavefunctions. It can insert SIESTA wavefunctions read from auxiliary
WFN file in between plane waves and orthonormalize them altogether.
It can apply a random variation to the plane waves and SIESTA wavefunctions
or generate completely random wavefunctions and orthonormalize them.
It can correct eigenvalues during orthonormalization or perform subspace
diagonalization. It can drop the plane waves and SIESTA wavefunctions which
have a large overlap with DFT wavefunctions. It can keep wavefunctions with
the eigenvalues in a given energy range. It can extract eigenvalues and
plane wave coefficients from WFN file for plotting purposes.

The SAPO method is described in [Phys. Rev. Lett. 107, 186404
(2011)](https://doi.org/10.1103/PhysRevLett.107.186404).

A few notes:

- When using iterative Davidson diagonalization, it is recommended to generate
  `VXC` rather than `vxc.dat` in `pw2bgw` run, then compute matrix elements of `VXC`
  between `SAPO` wavefunctions. This is more consistent than using `vxc.dat`
  because `vxc.dat` is computed between ESPRESSO wavefunctions, and those might
  be slightly different from `SAPO` wavefunctions.

- `wfng_input_file` should contain all occupied orbitals (and may or may not
  contain some unoccupied orbitals), otherwise occupations written to
  `wfng_output_file` will be wrong.

- Symmetrization of planewaves (`sapo_symmetry .gt. 0`) has no effect on a SAPO
  calculation since it is just a linear combination in a degenerate subspace.
  It is currently disabled because it requires symmetry subroutines from
  Quantum ESPRESSO which are not compatible with our BSD-type license.  To
  enable it: open `MeanField/SAPO/pw.f90`; comment out `die`; uncomment calls
  to `s_axis_to_cart`, `find_group`, `set_irr_rap`, and `divide_class`; and
  link the corresponding routines from Quantum ESPRESSO.

## `SAPO` input file (`sapo.inp`)

Below is a typical example of a `SAPO` input file, `sapo.inp`:

```Fortran
&input_sapo
   wfng_input_file = 'wfng.in'
   wfng_aux_file = ''
   vlr_input_file = ''
   vnlg_input_file = ''
   wfng_output_file = 'wfng.out'
   sapo_band_number = 0
   sapo_planewave_min = 1
   sapo_planewave_max = 580
   sapo_energy_shift = 0.0
   sapo_energy_match = .true.
   sapo_symmetry = 2
   sapo_print_ir = .true.
   aux_flag = .false.
   aux_band_min = 0
   aux_band_max = 0
   aux_energy_shift = 0.0
   sapo_random = 1
   sapo_random_ampl = 1.0d-3
   sapo_overlap_flag = .false.
   sapo_overlap_max = 0.0
   sapo_orthonormal = .true.
   sapo_ortho_block = 0
   sapo_ortho_order = 0
   sapo_ortho_energy = .false.
   sapo_energy_sort = .false.
   sapo_hamiltonian = .false.
   sapo_energy_range = .false.
   sapo_energy_min = 0.0
   sapo_energy_max = 0.0
   sapo_check_norm = .false.
   sapo_plot_kpoint = 0
   sapo_plot_spin = 0
   sapo_plot_bandmin = 0
   sapo_plot_bandmax = 0
   sapo_plot_pwmin = 0
   sapo_plot_pwmax = 0
   sapo_eigenvalue = .false.
   sapo_projection = 0
   sapo_amplitude = .false.
   sapo_ampl_num = 0
   sapo_ampl_del = 0.0
   sapo_ampl_brd = 0.0
/
```

In this example, the wavefunctions are read from `wfng.in` file and written to
`wfng.out` file. The parameter `sapo_band_number` specifies the number of
wavefunctions to be read from the input file (if set to 0 all the wavefunctions
will be read). The PW wavefunctions will be constructed from plane waves
ranging from 1 (`sapo_planewave_min`) to 580 (`sapo_planewave_max`). The
kinetic energies of the plane waves will be shifted to match the eigenvalues of
the input wavefunctions (`sapo_energy_match`). The plane waves will be
decomposed into irreducible representations of the space group of the Bravais
lattice (`sapo_symmetry = 2`), and the decomposition will be printed to the
standard output (`sapo_print_ir`). The auxiliary wavefunctions will not be
inserted in between the PW wavefunctions (`aux_flag`). A random variation
(`sapo_random = 1`) with a small amplitude (`sapo_random_ampl = 1.0d-3`) will
be added to the PW wavefunctions. The PW wavefunctions will be orthonormalized
with respect to the valence and conduction bands (`sapo_orthonormal`) in the
ascending order with respect to energy (`sapo_ortho_order = 0`). The
orthonormality check will not be performed (`sapo_check_norm`). The energy
eigenvalues, the projections of wavefunctions onto plane waves, and the squared
absolute values of amplitudes of the wavefunctions with respect to kinetic
energies of plane waves will not be plotted (`sapo_eigenvalue`,
`sapo_projection`, and `sapo_amplitude`).

The auxiliary wavefunctions can be used in two different ways:

1. Use hybrid valence wavefunctions as `wfng_input_file` and LDA conduction
   wavefunctions as `wfng_aux_file`, set `sapo_planewave_min` and
   `sapo_planewave_max` to 0, `sapo_orthonormal` to `.true`. and
   `sapo_ortho_block` to 1. This way you will get a good starting point for
   nscf hybrid calculations in PARATEC or ESPRESSO.

2. Use Siesta wrapper in `MeanField/SIESTA` directory to construct the resonant
   states from SIESTA wavefunctions (look into `MeanField/SIESTA/README` for
   details on how to use `siesta2bgw`), feed them to the SAPO code through
   `wfng_aux_file` parameter, and generate continuum states from plane waves by
   setting `sapo_planewave_max` to a large number. Also set `sapo_orthonormal`,
   `sapo_ortho_block` and `sapo_ortho_order` to orthonormalize the resonant and
   continuum states in the desired order with respect to the bound states from
   PARATEC or ESPRESSO. You may want to set `sapo_ortho_energy` to `.true.` for
   correcting the eigenvalues during the orthonormalization, or to set
   `sapo_hamiltonian` to `.true.` for correcting both the wavefunctions and the
   eigenvalues by diagonalizing the Hamiltonian. The latter requires the local
   potential file in Gaussian Cube format (`vlr_input_file`), and optionally
   the non-local pseudopotential file (`vnlg_input_file`).  Finally, set
   `sapo_energy_sort` to `.true.` for sorting the resonant and continuum states
   by their eigenvalues before writing them to the output file. You may also
   set `sapo_energy_range` to `.true.` for keeping the eigenvalues in the
   energy range of `sapo_energy_min` to `sapo_energy_max`.

!!! hint
    If you don't know or don't want to manually set `sapo_planewave_min` and
    `aux_band_min`, set both of them to 1, `sapo_overlap_flag` to `.true.`, and
    `sapo_overlap_max` to 0.3.  This will automatically throw away
    plane-waves and auxiliary states that have large overlaps (scalar
    products > 0.3) with DFT states and among themselves. The scalar
    products of PW and AUX states with DFT states will be written to files
    `overlap_dft_k#_s#.dat`, and the scalar products of PW states with AUX
    states will be written to files `overlap_aux_k#_s#.dat`, so you can
    inspect what states were thrown away after the run. <!---This useful
    feature is informally known as boverlap, in honor of Brad Malone who
    first proposed this idea. This confirms once again that evolution is a
    byproduct of laziness.--->
