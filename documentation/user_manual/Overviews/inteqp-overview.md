# `Inteqp` code overview

This code takes the `eqp_co.dat` file from Sigma and
interpolates it from the coarse to the fine grid using wavefunction
projections (to resolve band crossings) and linear interpolation
in the Brillouin zone.

## Summary of input and output files

### Required input files:

- [`inteqp.inp`](inteqp-keywords.md) : Input parameters.

- `WFN_fi` : Wavefunctions in unshifted fine grid.

- `WFN_co` : Wavefunctions on coarse grid.

- `eqp_co.dat` : A list of quasiparticle energy corrections for the bands in WFN_co.

### Additional input

- `WFNq_fi` : Wavefunctions in shifted fine grid (used for valence bands if use_velocity is selected)

### Auxiliary Files 

The files are output files from previous runs - used as input to speed up calculation.

- `dtmat` : Transformation matrices, dcc/dvv use for interpolation
           between coarse and fine grid.  This file must be consistent
           with your bsedmat and bsexmat files and corresponding
           coarse and fine wavefunctions.
           NOTE: the file format for dtmat was changed in BerkeleyGW 1.1.0 (r5961)

### Output Files

- `bandstructure.dat` : The GW bandstructure on the fine grid.

- `eqp.dat` : Quasiparticle energy corrections for the bands in WFN_fi.

- `eqp_q.dat` : Quasiparticle energy corrections for the bands in WFNq_fi (if use_velocity).

- `dvmat_norm.dat` : The norms of the dvv overlap matrices between the valence
                     band k on the fine grid and the closest k-point on the coarse grid.
                     
- `dcmat_norm.dat` : The norms of the dcc overlap matrices between the conduction
                     band k on the fine grid and the closest k-point on the coarse grid
