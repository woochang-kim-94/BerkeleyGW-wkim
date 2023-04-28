# MeanField utilities

Here are a number of miscellaneous utilities related to mean-field files and calculations:

### `kgrid.x`

Utility to generate k-points for BerkeleyGW calculations, originally written
for Quantum ESPRESSO DFT calculations. More information [here](espresso-overview.md#kgridx).
The format of the input file for `kgrid.x` is described [here](kgrid-input.md).


### `degeneracy_check.x`

Determine which numbers of bands do not break a degenerate subspace, to avoid warnings or
errors in Epsilon, Sigma, and BSE codes. Supply any number of binary WFN files as command-
line arguments, and a list of numbers that are acceptable for all k-points in all files
will be written to standard output.


### `mf_convert_wrapper.sh`

Converts WFN/RHO/VXC files between binary and ASCII. This can be useful for moving files
correctly between big- and little-endian machines, or for examining the contents of a
file (though wfn_rho_vxc_info.x is more user-friendly). Using this wrapper script, the
flavor (real/complex), type of file, and binary/ASCII is determined automatically.


### `wfn_rho_vxc_info.x`

Prints the contents of the header of a (binary) WFN, RHO, or VXC file in a clearly labeled
format, for human inspection.


### `wfnmerge.x`

Merges many WFN files into one. It assumes that all input files have the 
same number and same ordering of G-vectors for the charge density.
The number and name of input files is read from `wfnmerge.inp`, as well
is the kgrid and the kshift.

Q: Why would someone want to use this?

A: For example, maybe to do a converged calculation one needs to include a
    large number of unoccupied states. There may not be the resources (either
    CPUs or wallclock) to do this all in one shot, so it may be beneficial to
    split up a calculation of kpoints (e.g., 4x4x4 MP grid) into smaller
    pieces and then add up the final wavefunctions.

Q: What is the deal with kpoint weights?

A: Quantum Espresso renormalizes the kpoint weights that you use in a 
   calculation. If you are splitting up a MP grid into different groups of 
   kpoints, and each group is normalized, the relative weights between
   kpoints in different groups is no longer correct. BerkeleyGW does not 
   make use of the kpoint weights (except for in determining whether the 
   Fermi level is reasonable for metallic calculations). So for most uses
   the fact that these weights are incorrect does not matter. If it matters
   to you, you can simply modify wfnmerge to read in the kpoint weights of
   your choosing.

### `scissors2eqp.x`

Write an `eqp.dat` file based on a WFN file and scissors parameters.
For testing equivalence of `eqp_corrections` and scissors.


### `fix_occ.x`

Fix the occupations of a WFN file. This is useful if you have split a
large nscf calculations into smaller blocks, and you want the occupations
of the merged WFN to be consistent. It can also fix inconsistencies in the
occupations and in the Fermi energy for metals.


### `wfn_dotproduct.x`

Form the overlap between the bands in two wavefunction files.
If they are copies of the same file, you can check orthonormality.
Only bands at corresponding k-points are considered, since
the overlap is zero by Bloch's theorem if the k-points differ.
Warning: a known issue is possible errors from gmap when relating
k-points by symmetry operations.


### `wfn_time_reversal.x`

Give a WFN file with a k-grid reduced by time-reversal symmetry, and get a new WFN file
with the k-points unfolded by time-reversal symmetry. This can save up to half the time
in the non-self-consistent DFT calculations for input to BerkeleyGW. It can be used
after calculating k-points from `kgrid.x` with time-reversal symmetry, or from 'automatic'
k-points from Quantum ESPRESSO. Either use VXC, or add the new k-points to vxc.dat by
hand according to the rule `<n,-k|Vxc|m,-k> = <n,k|Vxc|m,k>*`.

### `hdf2wfn.x`

Converts a WFN from binary Fortran format to the HDF5 format.


### `wfn2hdf.x`

Converts a WFN fromthe HDF5 format to the binary Fortran format.


#BEGIN_INTERNAL_ONLY
The files below are unreleased because they are still under development and/or not fully
validated.

### `analyzebz.x`

Creates fullbz maps from WFN file, Real/Complex flavor.


### `wfnreduce.x`

Reduces the cutoff of WFN files, Real/Complex flavor. Not properly working.

#END_INTERNAL_ONLY
