# Clustered Sampling Interpolation (CSI)

The CSI allows you to obtain converged optical absorption spectra of quasi-2D
materials with a relatively coarse k-point grid. Before you proceed with the
CSI calculation, make sure you understand the formalism and the first
mean-field calculation that has to be performed, as [documented in the overview
of the subsampling approach](subsample-overview.md). You may also want to first
compute the electronic self-energy using the [NNS method](NNS.md).


## Pre-setup: `kernel` on a coarse grid
First, run `kernel.x` as usual on a coarser, uniform grid. Note that the
converged grid for the `kernel` calculation is typically finer than the one
used to converge `epsilon` and `sigma` calculations as in a NNS calculation.
Use `epsmat` and `esp0mat` files from a usual, uniform `epsilon` calculation.


## Setup: `setup_subsampling_csi.x`
Analogous to the NNS, run `setup_subsampling_csi.x format WFN_co WFN nk_fi_x
nk_fi_y nk_fi_z`, where `format` is either `ASCII`, `BIN`, or `HDF5`, `WFN_co`
is the coarse wave function used to calculate the BSE kernel, and `nk_fi_x
nk_fi_y nk_fi_z` are the find kgrid to be used in BSE. This must be an integer
multiple of the coarse grid.

The output files are:

- `epsilon_q0s.inp`: contains part of the `epsilon.inp` file needed to generate
  `epsmat` for the clustered points.
- `kpoints_wfnq.dat`: contains a list of kpoints needed for the `WFNq` used to
  calculate `epsmat`.
- `kpoints_sub_*.dat`: contains a list of subsampled k-points surrounding each
  coarse point. There is one file per coarse point.
- `subsample.inp`: contains the header for an input file needed during
  absorption.

Use the k-points in `kpoints_wfnq.dat` to generate a `WFNq` file with the
mean-field code.


## Mean-field and `epsilon` calculations

- Run epsilon using `epsilon_q0s.inp` as the basis for `epsilon.inp` and link
  to the `WFNq` file generated above. The output file is called `eps0mat.h5`, but
  it is referred to as `epsmat_sub.h5` for clarity in the rest of this
  document.

- Generate a `WFN` file containing the k-points in each `kpoints_sub_*.dat`
  file.  For example, for a 2x2 coarse grid, you will have 4 `kpoints_sub_*.dat`
  files, so you will need to generate 4 wavefunctions, `WFN_1`, `WFN_2`, `WFN_3`,
  `WFN_4`.

- Run kernel for every `WFN_1`, `WFN_2`, etc. The recommended way to do this is to
  have a separate directory for each `WFN_*` file and then set up symbolic links to
  the correct input files. For example: In `directory_1/` link `eps0mat.h5` to
  the same `eps0mat.h5` used in the coarse kernel calculation. Link the 
  subsampled dielectric function into `epsmat.h5` with `ln -s epsmat_sub.h5 epsmat.h5`,
  (bear in mind our naming convention here), and the wavefunction as `ln -s WFN_1 WFN_co`.
  In `directory_2/`, link both `eps0mat.h5` and `epsmat.h` as in `directory_1/`, but
  use the corresponding wavefunction with `ln -s WFN_2 WFN_co`. Repeat for
  other directorie.


## `Kernel` calculation
In `kernel.inp`, make sure to add the following flags:
`no_symmetries_coarse_grid`, and `patched_sampling_co`.

!!! Warning 
    When you run the `kernel` calculation, make sure that $n_k^2$ > number of
    processors, where $n_k$ is the number of k-points. This is a current
    limitation of how the CSI is implemented.


## `Absorption` calculation

- Move the generated `subsample.inp` file to the directory where you will run
  your absorption calculation.

- Open and edit `subsample.inp`. You can find a sample `subsample.inp` file
  that explains the file format in `BSE/subsample.inp` in your BerkeleyGW
  source directory.

- Run absorption. Link the `bsemat` coarse file, `WFN_co` is same as in the
  coarse kernel calculation, and `WFN_fi` is as usual, with a size that has to
  be converged. In `absorption.inp` include the following flag: `subsample_line
  cutoff`, where `cutoff` is your coarse grid spacing. For instance, if your
  coarse grid is `30x30x1`, the cutoff is 0.0333. Everything else in
  `absorption.inp` should be the same as your previous calculation. Make sure
  that `subsample.inp` is in the same directory as `absorption.inp`.
