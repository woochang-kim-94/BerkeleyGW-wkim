# `PlotXct` code

## Overview

PlotXct plots the exciton wavefunction in real space based on output of
BSE/absorption (with full diagonalization).

Reads files `WFN_fi`/`WFNq_fi` (from mean field) and eigenvectors (output from
BSE/absorption code when using diagonalization, not Haydock) and plots selected
exciton state in real space according to:

$$
\Psi(r,r',s) = \sum_{cvk} A_{svck} \; u_{sck}(r) u_{svk}^*(r')
$$

where $r'$, the hole coordinate, is fixed at some point and $r$, the electron
coordinate, runs over a supercell of given size (typically, equivalent to the
fine k-grid), with a real-space mesh defined by the FFT grid of the
wavefunction files. A separate plot can be generated for each spin s.

<!---
An alternate approach is to project onto a Wannier function $|W_n\rangle$:

   $\Psi(r,n,s) = \sum_{cvk} A_{svck} \; u_{sck}(r) \langle u_{svk}| W_n\rangle$

Output is written in format suitable for Data Explorer (DX).
--->

The input option `restrict_kpoints` reduces the sum over k-points
above to a sum over the ones that give most of the contribution
to the norm of eigenvectors. This is handy if there are many k-points
but only a few of them give sizable contribution.

Input files:

- `WFN_fi`, `WFNq_fi`: wavefunction files
- `eigenvectors`:  output from BSE/absorption code
- `plotxct.inp`: `plotxct` input parameters
<!---
       (seedname.chk        checkpoint from Wannier90 code)
--->

## `PlotXct` Instructions:

1. Copy/edit `plotxct.inp`. You'll need the:
    - index of the exciton state you want to plot
    - hole position in lattice coordinates
    - supercell dimensions (or k-grid)
    - spin component you want to plot (if spin-polarized)

2. Run `plotxct.x`.  The ASCII file `xct.[state]_s[spin].a3Dr` will be created.
    `a3Dr` means that the file is in ASCII and contains 3D data in $r$ real
    space.  The header of this file contains information on:
    - state index
    - state energy in eV
    - hole position in atomic units (a.u.)
    - spin component plotted
    - supercell lattice vectors in a.u.
    - number of discretization points

    followed by three columns of data corresponding to the complex values of the
    electronic part of the excitonic wavefunction, and its magnitude squared.
    The values are written with very low precision due to file-size concerns.


3. Convert to a format readable by the plotting utility of your choice.  Use
     the `volume.py` utility:  `volume.py imfn imff ivfn ivff ovfn ovff phas
     cplx [hole]`. Parameters:
       - `imfn`: input matter file name
       - `imff`: input matter file format:
               `{mat|paratec|vasp|espresso|siesta|tbpw|xyz|xsf}`
       - `ivfn`: input volumetric file name
       - `ivff`: input volumetric file format `{a3dr}`
       - `ovfn`: output volumetric file name
       - `ovff`: output volumetric file format `{cube|xsf}`
       - `phas`: remove wavefunction phase `{true|false}`
       - `cplx`: complex wavefunction `{re|im|abs|abs2}`
       - `hole`: plot hole `{true|false}`

     You can specify whether you want a `.xsf` (XCrysDen) or `.cube` (Gaussian Cube)
     file format. You can remove or keep the arbitrary phase of the wavefunction
     by setting parameter "phas" to true or false. You can plot the re, im, abs,
     or abs2 part of the wavefunction.

     You must specify a path to the input file (example: `../01-scf/input`)
     in one of the supported formats (DFT and others) so that the script can
     obtain the atomic positions from the file.

     You can choose whether to display or not the position of the hole by setting
     parameter "hole" to true or false.

     The resulting `.xsf` or `.cube` file can be viewed in XCrysDen or converted to
     POV-Ray script format using [`Surface` code](surface.md).

## Example
```shell
~/BerkeleyGW/Visual/volume.py ../01-scf/input paratec xct.000001_s1.a3Dr a3dr xct.000001.abs2.xsf xsf false abs2 false
xcrysden --xsf xct.000001.abs2.xsf
```
