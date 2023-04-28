# `NonLinearOptics` code overview

This code is intended to be run after a full GW/BSE calculation and
generates the UltraFast or 2-photon spectrum of a material by calculating
the optical transition strength between the lowest excited exciton state
(user inputed) and higher excited exciton states.

!!! Danger
    The `NonLinearOptics` ode should be considered alpha quality. It has been
    tested on only a few systems. Use at your own risk and validate and, if
    necessary, contact the developers.

### Required input files

- [`nonlinearoptics.inp`](nonlinearoptics-keywords.md) : Input parameters.  

- `eigenvectors` : The exciton wavefunctions in Bloch space.

- `WFN_fi` : Bloch wavefunctions on the unshifted fine grid.

- `WFNq_fi` : Bloch wavefunctions on the shifted fine grid. Only needed if velocity operator is used.

<!--
### Additional input

- `eqp.dat` : This is not yet implemented (energies are read from A_svck)
              If it were to work, it would probably be:
              A list of quasiparticle energy corrections for the bands
              in WFN_fi and WFNq_fi. If absent the scissors operator from
              nonlinearoptics.inp is used.
-->

### Auxiliary Files

The files are output files from previous runs - used as input to speed up calculation

- `vmtxel_nl`: The interband optical matrix elements.  This is different
               from the file `vmtxel` created by absorption in that it contains
               valence-valence and conduction-conduction matrix elements.
               Specify `read_vmtxel_uf` in the `nonlinearoptics.inp` file and files
               `WFN_fi` and `WFNq_fi` are not needed.
