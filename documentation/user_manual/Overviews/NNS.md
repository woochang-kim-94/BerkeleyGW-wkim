# Nonuniform Neck Subsampling (NNS)

The NNS allows you to obtain converged electronic self-energy on quasi-2D
materials with a relatively coarse k-point grid.  Before you proceed with the
NNS calculation, make sure you understand the formalism and the first
mean-field calculation that has to be performed, as [documented in the overview
of the subsampling approach](subsample-overview.md).


## Setup: `setup_subsampling_nns.x`
Run `setup_subsampling_nns.x`, which is located in the `Epsilon` subdirectory
of BerkeleyGW. The syntax is: `setup_subsampling_nns.x format WFN`, where
`format` specifies the format of the input WFN file, which can be either
`ASCII`, `BIN`, or`HDF5`. If you generated the file `WFN` with Quantum
ESPRESSO, you will want to select `BIN`; if you generated a `WFN.h5` file with
[ParaBands](parabands-overview.md), then you want to select `HDF5`. The script
also supports other optional parameters, which are documented by running the
script without any argument. The default options work well for many quasi-2D
semiconductors, including monolayer MoS2.

The output of `setup_subsampling_nns.x` is a list of k-points to generate
q-shifted a `WFNq` file. There are two ways you can generate this: you can
generate a single `WFNq` file which includes all possible shifted k points and
which allow you to compute $\varepsilon(q)$ for all possible wavevectors $q$.
This is the simpler way to use the script, but the generated `WFNq` may too
large and you may run out of memory depending on your calculation. The other
option is to split the calculation, and generate a series of `WFNq` files, one
for each calculation of $\varepsilon(q)$ for different $q$. This requires a bit
more extra work, such as running several instances of Quantum ESPRESSO, and
separately using each generated `WFNq` file to compute $\varepsilon(q)$ for a
different $q$ point, and merging the generated `eps0mat.h5` files. However,
this will use significantly less memory in BerkeleyGW.

The output files produced by `setup_subsampling_nns.x` are:

- `kpoint_all.dat`: list of all the shifted k points needed to generate a
  single `WFNq` file.
- `kpoint_xxx.dat`: list of k points needed to generate a `WFNq` file of the
  specific `xxx` qpoint by the `epsilon` code.
- `epsilon_q0s.inp`: list of qpoints for the `epsilon` input file.
- `subweights.dat`: file to be used as input in the `sigma` calculation.


## Mean-field calculations

Choose one of the two following ways:

- Use the list of all k-points from `kpoint_all.dat`, then use `pw2bgw` to
  generate one (often very large) `WFNq` file. This is the easiest way to
  operate, but you may run out of memory.

- Use a separate list of k-points for each qpoint separately. In this case
  compute each different `WFNq` file in a separate directory, use
  `kpoint_xxx.dat` for the list of kpoints of each one, and generate `WFNq` for
  each one separately


## Dielectric matrix: `epsilon`
Use the `WFN` file that you provided to `setup_subsampling_nns.x` for the
`epsilon` calculation, which should include many empty bands. In `epsilon.inp`,
replace the usual gamma qpoint of `[0 0 0]` with the q-list from
`epsilon_q0s.inp`. The uniform grid kpoint list is added as usual. Use the
keyword `subsample` that should also be present in `epsilon_q0s.inp`. Also,
remember to use the needed truncation scheme for 2D systems, such as
`cell_slab_truncation`. See an example `epsilon.inp` file below. BerkeleyGW
will produce an `eps0mat.h5` file which now contains a number of
$q\rightarrow0$ points required for the NNS calculation!

!!! Note
    You can perform epsilon for each qpoint separately. In the end, you will
    need to combine all `eps0mat.h5` files together and all `eps0mat.h5` files
    together. To do that you can use the BGW script `epsmat_hdf5_merge.py`.


## Self-energy: `sigma`
Run sigma with `eps0mat` and `epsmat` as usual, and in `sigma.inp` file:

- Use the keyword `subsample` in the input file, as well as the truncation you
  used in `epsilon.inp`, e.g., `cell_slab_truncation`.
- Copy/link the file `subweights.dat` produced by `setup_subsampling_nns.x` to
  the directory in which you run sigma.


If you are interested in computing optical absorption spectra, you may also
want to read about the [CSI method](CSI.md).
