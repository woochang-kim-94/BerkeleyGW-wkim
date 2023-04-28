# `Visual` package

BerkeleyGW ships with a number of visualization utilities to simplify your work
with the DFT codes and the BerkeleyGW package. Here you will find the following
scripts and codes from this `Visual` package.

## Surface

Surface is a C++ code for generating an isosurface of a volumetric scalar field
(such as the wave function, charge density, or local potential). The scalar
field is read from Gaussian Cube or XCrySDen XSF file, the surface
triangulation is performed using marching cubes or marching tetrahedra
algorithm, and the isosurface is written in the POV-Ray scripting language.
The final image is rendered using the ray-tracing program POV-Ray.  Running the
code requires a fairly complicated input parameter file, described with an
example in [surface.inp](surface.md).


## Matter

Matter is a python library for manipulating atomic structures with periodic
boundary conditions. It can translate and rotate the atomic coordinates,
generate supercells, assemble atomic systems from fragments, and convert
between different file formats. The supported file formats are `mat`, paratec,
vasp, espresso, siesta, `xyz`, `xsf`, and povray.

`mat` is the native file format of the library. paratec, vasp, espresso, and
siesta represent the formats used by different plane-wave and local-orbital DFT
codes. `xyz` is a simple format supported by many molecular viewers, and `xsf`
is the internal format of XCrySDen viewer.  povray stands for a scripting
language used by the ray-tracing program POV-Ray.

<!-- The core of the library consists of files common.py, matrix.py, and
matter.py. Script convert.py is a command-line driver that performs the basic
operations supported by the library. Script link.py is a molecular assembler
that can be used to rotate and link two molecules together. Script gsphere.py
generates a real-space grid and a sphere of G-vectors given lattice parameters
and a kinetic energy cutoff.  This helps to estimate the number of unoccupied
states needed in GW calculations. Script average.py takes an average of the
scalar field on the faces or in the volume of the unit cell. This is used to
determine the vacuum level in DFT calculations. Script volume.py converts the
a3Dr file produced by PARATEC or BerkeleyGW to Gaussian Cube or XCrySDen XSF
format.  --> The main end-user scripts from this library are:

- `convert.py`: performs basic conversion operations supported by the library.
- `link.py`: molecular assembler that can be used to rotate and link two
  molecules together.
- `gsphere.py`: generates a real-space grid and a sphere of G-vectors given
  lattice parameters and a kinetic energy cutoff.  This helps to estimate the
  number of unoccupied states needed in GW calculations.
- `average.py`: takes an average of the scalar field on the faces or in the
  volume of the unit cell. This is used to determine the vacuum level in DFT
  calculations.
- `volume.py`: converts an a3Dr file produced by PARATEC or BerkeleyGW to
  Gaussian Cube or XCrySDen XSF format.

Each script requires a different set of command-line arguments.  Running
individual scripts without arguments displays a list of all possible
command-line arguments and a short description of each argument.


### Examples

#### Benzene Charge Density

This example demonstrates how to plot isosurfaces. Go to directory
`BGW/Visual/benzene` and run ESPRESSO:

```shell
sh link
```

Customize and submit `script`. (suggested ncpu = 36, walltime = 0:30:00)

This will produce the Gaussian Cube file `rho.cube` containing the charge
density of the benzene molecule. For your convenience, the compressed file
`rho.cube.tgz` is included in the directory.

We can generate an isosurface that contains 90% of the charge density using the
marching cubes algorithm with the smooth triangles and render it in POV-Ray:

```shell
../surface.x rho_uc.inp
povray -A0.3 +W640 +H480 +Irho_uc.pov
```

Examine the `rho_uc.gif` file to find that the pieces of the benzene molecule
are placed in the corners of the unit cell. To assemble the pieces together we
define a supercell in the `rho_sc.inp` file. This is done by setting the
parameter `sc` to `T` and by placing the supercell origin (`sco`) at position
(-0.5 -0.5 -0.5) in crystal coordinates (`scu` which stands for supercell units
is set to latvec). The lattice vectors of the supercell are not changed (`scv`
or supercell vectors in the parameter file). One can render a new image:

```shell
../surface.x rho_sc.inp
povray -A0.3 +W640 +H480 +Irho_sc.pov
```

Now the charge density comes up nice and centered in the middle of the
supercell.


#### Nanotube Exciton Wavefunction

In this example we will plot an isosurface of the exciton wavefunction. Run
`BGW/examples/DFT/swcnt_8-0/` fully. Then go to directory
`BGW/examples/DFT/swcnt_8-0/8-absorption/`.

PlotXct will produce the a3Dr file `xct.020_s1.a3Dr` containing a wavefunction
of the 20th exciton in the (8,0) nanotube.  Then `volume.py` will convert the
a3Dr file to Gaussian Cube:

```shell
volume.py ../ESPRESSO/1-scf/in espresso xct.020_s1.a3Dr a3dr xct020.cube cube false abs2 true
```

The command will produce the Gaussian Cube file `xct020.cube` that contains the
squared absolute value of the wavefunction (`cplx = abs2`) without the sign
(`phas = false`) and the position of the hole (hole = true). Then surface will
generate an isosurface that contains 90% of the charge density:

```shell
../surface.x xct020.inp
```

You can then render it in POV-Ray (which you must install separately):

```shell
povray -A0.3 +W2560 +H1920 +Ixct020.pov
```

The resulting image `xct020.gif` is included in the example directory.
Alternately, you can visualize the Cube file directly in XCrySDen.


#### Rock-Salt Lattice

Here you will learn how to make large supercells of bulk crystals.  Go to
directory `BGW/Visual/rocksalt` where you will find the `nacl.mat` file.
Render it in POV-Ray:

```shell
../convert.py nacl.mat mat nacl.pov povray
povray -A0.3 +W640 +H480 +Inacl.pov
```

The unit cell contains only two atoms. Let us make a supercell that contains a
fractional number of unit cells. The supercell is described by the `sc1` file
which defines the origin of the supercell, the lattice vectors of the
supercell, and the units in which these quantities are given. These are
equivalent to `sco`, `scv`, and `scu` in the surface parameter file in the
first example. The last line in the `sc1` file corresponds to `sct` (supercell
translation) in the surface parameter file. We will get back to it later on.
Let us render the supercell:

```shell
../convert.py nacl.mat mat nacl_sc1.pov povray supercell sc1
povray -A0.3 +W640 +H480 +Inacl_sc1.pov
```

Note the broken bonds on the faces of the supercell. That is because the
supercell contains a fractional number of the unit cells, so the translational
symmetry is broken which results in the Na-Na and Cl-Cl bonds. We disable the
translational symmetry of the supercell by setting `sct` to false (the last
line in the `sc2` file). This enforces the translational symmetry of the
underlying unit cell structure. Render the new supercell in POV-Ray:

```shell
../convert.py nacl.mat mat nacl_sc2.pov povray supercell sc2
povray -A0.3 +W640 +H480 +Inacl_sc2.pov
```

Now the bonds on the faces of the supercell come out right.


#### Organic Molecule Synthesis

In this example you will assemble the bipolar molecule, bithiophene naphthalene
diimide, from the donor and acceptor units, bithiophene and naphthalene
diimide. Go to directory `BGW/Visual/btnd` and open files `nd.xyz` and `bt.xyz`
in any molecular viewer, for example Jmol.  Identify the atoms you want to link
together, these are the two carbon atoms, # 15 in `nd.xyz` and # 8 in `bt.xyz`.
You need to remove the two hydrogen atoms, # 21 in `nd.xyz` and # 14 in
`bt.xyz`, and place # 15 and # 8 at the distance of 1.49 Angstrom from each
other along the 15-21-14-8 line. This is done by invoking the following
command:

```shell
../link.py nd.xyz xyz bt.xyz xyz btnd.xyz xyz 15 21 8 14 0.0 degree 1.49 angstrom
```

You can also rotate `bt.xyz` around the 15-21-14-8 axis by an arbitrary angle,
although the rotation angle is set to zero in the above command.

#### Estimate the number of unoccupied states for `epsilon` and `sigma`

Go to directory `BGW/examples/DFT/benzene/0-gsphere/` and run:

```shell
gsphere.py g_eps.in g_eps.out
```

Examine input file `g_eps.in`. It contains the lattice vectors in bohr, the
cutoff energy in Ry (epsilon_cutoff from `../5-gw/epsilon.inp` or
`screened_coulomb_cutoff` from `../5-gw/sigma.inp`), `0 0 0` for the FFT grid
to be determined automatically, and 'true' to sort the G-vectors by their
kinetic energies. Examine output file `g_eps.out`. Look for the number of
G-vectors, `ng`. You will find `ng = 2637`, meaning that you will need at least
that many unoccupied states. Of course the number of unoccupied states is a
convergence parameter, but `gsphere.py` gives you an idea in which range of
values should you check the convergence.

Now run the following:

```shell
gsphere.py g_rho.in g_rho.out
```

Input file `g_rho.in` is similar to `g_eps.in` but it has a different cutoff
energy which is 4 times `ecutwfc` from `../ESPRESSO/1-scf/in`.  This is the
charge density cutoff or ecutrho in espresso documentation.  Note that the last
line in `g_rho.in` is set to `false`. This is because using built-in Python
sort function on that many G-vectors (up to `ecutrho`) will freeze your python
interpreter. Examine output file `g_rho.out`. You should see the following:

```
 grid = ( 128 72 140 ) -- minimal
        ( 128 72 140 ) -- factors 2, 3, 5, 7, 1*11, 1*13
        ( 128 72 144 ) -- factors 2, 3, 5
```

The third line is the size of FFT grid in espresso calculation.  (Well, not
necessarily, because espresso algorithm starting with version 5 is smarter than
what is implemented in `gsphere.py`.) The z-component is equal to 144 meaning
that for an optimal load-balancing you should run espresso on a number of cores
which is a divisor of 144 (that is, 144, 72, 36, 18, etc. cores).


### Notes:

Matter library may lack support for some advanced features of PARATEC, VASP,
ESPRESSO and SIESTA formats. For example, `LatticeParameters` and `ZMatrix` are
not implemented for SIESTA.  This can be added to functions `paratec_read`,
`paratec_write`, `vasp_read`, `vasp_write`, `espresso_read`, `espresso_write`,
`siesta_read` and `siesta_write` in file `matter.py`.

Also note that all the additional information not related to the crystal
structure (k-points, energy cutoffs, etc.) will be lost and replaced with the
default parameters when converting between PARATEC, VASP, ESPRESSO and SIESTA
formats.
