# Paratec code

Paratec is a simple DFT code optimized for small- and medium-sized systems, and
with a tight integration with BerkeleyGW. Please, refer to the
[full documentation for PARATEC 5.1.12](https://drive.google.com/open?id=0BxTQVSgElL8nSnZDbEdhRVBoclE).


## Wrapper

PARATEC output for BerkeleyGW is controlled by flags 
`gw_shift`, `gwc`, `gwr`, `gwscreening`, `gwcscreening`, and `vxc_matrix_elements`.
The flags can be combined with an underscore: e.g. output_flags `gwr_gwscreening`
`gwr` and `gwc` are incompatible; `gwscreening` and `gwcscreening` are incompatible.

Main flags:

-   `gw_shift q1 q2 q3`

    Generates q-shifted grid, q-vector is in crystal coordinates
    in units of reciprocal lattice vectors (for `WFNq`, `WFNq_fi`)
    This variable does the same job as the `kgrid.x` utility.

-   `output_flags gwc`

    Writes complex wavefunctions in G-space, for systems without
    inversion symmetry about the origin, to file `WFN` (for all codes).

-   `output_flags gwr`

    Writes real wavefunctions in G-space, for systems with
    inversion symmetry about the origin, to file `WFN` (for all codes)

-   `output_flags gwscreening`

    Writes charge density in G-space to file `RHO`,
    exchange-correlation potential in G-space to file `VXC`,
    and matrix elements of exchange-correlation potential to file
    `vxc.dat` (for Sigma). Real if possible and `gwc` not set, else complex.

-   `output_flags gwcscreening`

    Like `gwscreening`, except forces complex even if real is possible.

-   `vxc_matrix_elements diagmin diagmax offdiagmin offdiagmax`

    Specifies the range of bands for which to compute diagonal and
    off-diagonal matrix elements of exchange-correlation potential
    (for Sigma, in conjunction with output_flags gwscreening)

Other key input flags:

-   `pw_job {scf, nonselfcon}`

    Use scf for initial calculation, nonselfcon for generating
    BerkeleyGW outputs. (bandstructure does not seem to work)

-   `energy_cutoff ecut`

    Plane-wave cutoff for wavefunctions, in Ry.

-   `number_kpoints`

    - set to 0  to use k_grid and reduce with symmetries
    - set to -1 to use k_grid and do not reduce with symmetries
    - set to any other number to read from file KPOINTS

-   `k_grid nx ny nz`

    3 integers specifying Monkhorst-Pack k-grid dimensions

-   `k_grid_shift dx dy dz`

    Monkhorst-Pack k-grid shifts (typically 0.0 or 0.5)

-   `number_bands nb`

    Number of bands to use in calculation. Fraction actually useful
    or written to BerkeleyGW output determined by next variable.

-   `eigspacefrac frac`

    Fraction of bands to converge. Setting a higher number_bands
    and lower eigspacefrac can make the calculation more efficient
    depending on the diagonalization scheme. $0 < frac <= 1.0$.

You can find the actual input files for PARATEC and BerkeleyGW in 
`examples/DFT`, in PARATEC subdirectories for each example.

There are also `bgw2para` and `rho2cd` utilities that convert 
BerkeleyGW files `WFN` and `RHO` to PARATEC format. This may be 
useful, for example, if one generates the plane waves on top of the 
valence and conduction bands (look into [SAPO](sapo-overview.md) for details) 
and wants to diagonalize them further in PARATEC. There are no input 
files; `bgw2para` takes as argument the wfn filename,
and it creates files `WFN$n.$s` and `BAND` needed for `PARATEC`. 
Similarly, `rho2cd` requires file `RHO` and it creates file `CD`.


## Utilities


### `kptlist.pl`

Extracts a formatted list of k-points from PARATEC file for use in the Sigma code


### `qptlist.pl`

Extracts a formatted list of q-points from PARATEC file for use in the Epsilon code

## Additional information

To build with support for BerkeleyGW output, in `arch.mk` set the
line `GWWFNPATH` to `BerkeleyGW/library`, and add `-DBGW` to `M4OPTLIBS`.


Literature:

-   B. G. Pfrommer, J. Demmel, and H. Simon, "Unconstrained Energy Functionals
    for Electronic Structure Calculations," J. Comp. Phys. 150, 287 (1999).
-   B. G. Pfrommer, M. Cote, S. G. Louie, and M. L. Cohen, "Relaxation of
    Crystals with the Quasi-Newton Method," J. Comp. Phys. 131, 233 (1997).
-   Mathieu Taillefumier, Delphine Cabaret, Anne-Marie Flank, and Francesco
    Mauri, "X-ray absorption near-edge structure calculations with
    pseudopotentials: Application to the K-edge in diamond and alpha-quartz,"
    Phys. Rev. B 66, 195107 (2002).
-   <http://cmsn.drupalgardens.com/sites/cmsn.drupalgardens.com/files/CMSN_Newsletter_Vol2No2.pdf>

The pseudopotentials for PARATEC can be generated with the `fhi98PP` program
which is available for download at
<http://www.fhi-berlin.mpg.de/th/fhi98md/fhi98PP/>
