# `icm` (image charge model) code

`icm.x` reads the wavefunction file in Gaussian Cube or 
XCrySDen XSF format and computes the Coulomb integral within 
an image charge model.

It is based on [surface code](surface.md), so the input parameter 
file formats for the two codes are quite similar.

## Input file structure (`icm.inp`)

```
inputfilename C6H6.b_15.cube
inputfileformat cube
threshold 0.99
threshold_power 1
coulomb_power 1
mirrorplaneorigin
0.0 0.0 -2.0
mirrorplanenormal
0.0 0.0 1.0
mirrorplaneunit angstrom
uc F
uco
0.0 0.0 0.0
ucv
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
ucu latvec
sc T
sco
-0.5 -0.5 -0.5
scv
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
scu latvec
```

## Documentation

In the above example, the HOMO wavefunction of benzene is read from Gaussian
Cube file. The wavefunction is placed in the center of the supercell (see the
meaning of `uc`, `uco`, `ucv`, `ucu`, `sc`, `sco`, `scv`, `scu` parameters in
the [surface code](surface.md)). The parts of the wavefunction outside an
isosurface that contains 99% of the charge density are dropped (parameters
`threshold` and `threshold_power` have the same meaning as `isovalue` and
`power` in the [surface code](surface.md)). Parameter `coulomb_power` tells the
code whether the wavefunction in the Coulomb integral needs to be squared. Set
both powers to 1 if the wavefunction file contains the squared amplitude as
produced by ESPRESSO, and to 2 for the linear amplitude as in PARATEC or
SIESTA. The mirror plane is defined by parameters `mirrorplaneorigin` and
`mirrorplanenormal`, in the above example it is parallel to the xy plane
crossing the z axis at -2 Angstrom. Problems may occur with non-orthogonal unit
cells; use of cubic cells is recommended.
