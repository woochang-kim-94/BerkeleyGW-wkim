# BerkeleyGW overview

BerkeleyGW is a [free, open source](license.md), and massively parallel
computational package for electron excited-state properties that is based on
the many-body perturbation theory employing the *ab initio* GW and GW plus
Bethe-Salpeter equation methodology.

It is able to calculate accurate electronic and optical properties in materials
of different dimensionalities and complexity, from bulk semiconductors and
metals to nanostructured materials and molecules.  

It can be used in conjunction with many external and well-established
density-functional theory codes for ground-state properties, including PARATEC,
Abinit, PARSEC, Quantum ESPRESSO, OCTOPUS and SIESTA. These codes are used to
generate initial files, containing the ground-state density and wavefunctions
from density-functional theory. In addition, BerkeleyGW also ships with two
codes to generate a large number of empty states for GW calculations:
[SAPO](sapo-overview.md) and [ParaBands](parabands-overview.md). See the page on
[mean-field](meanfield.md) calculations for further information.

After you [compile](compilation.md) and test BerkeleyGW, we suggest you follow the
following tutorials on how to run calculations with BerkeleyGW:

1. GW calculation: 
    1. [`epsilon`](https://berkeleygw.org/tutorial-epsilon/): evaluating the
       dielectric screening 
    2. [`sigma`](https://berkeleygw.org/tutorial-sigma/): calculating the
       electronic self-energy
2. Bethe-Salpeter equation (BSE) calculation:
    1. [`kernel`](https://berkeleygw.org/tutorial/tutorial-bethe-salpeter-equation/):
       calculating the electron-hole interaction kernel
    2. [`absorption`](https://berkeleygw.org/tutorial/tutorial-bethe-salpeter-equation/):
       computing neutral optical excitation properties, such as optical
       absorption spectrum.

For practical details of how to perform each part of the calculation,
as well as an overview of input keywords, see the
[Typical workflow](overview-workflow.md) section.

![General BerkeleyGW workflow](Overview_scheme.svg "General BerkeleyGW workflow")

