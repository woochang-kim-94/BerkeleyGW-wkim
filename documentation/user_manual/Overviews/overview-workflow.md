# Typical workflow

In order to to compute GW quasiparticle energies and optical absorption spectra
including excitonic effects with BerkeleyGW, one is required to run the
following steps:

 - A [`mean-field`](meanfield.md) program of their choice, such as Quantum ESPRESSO, PARATEC,
   Abinit, PARSEC, OCTOPUS, and SIESTA, and generate a set of files;
 - The [`epsilon`](epsilon-overview.md) executable from the BerkeleyGW package, to
   generate the dielectric matrix;
 - The [`sigma`](sigma-overview.md) executable from the BerkeleyGW package, to
   compute the GW quasiparticle energies;
 - The [`kernel`](kernel-overview.md) executable from the BerkeleyGW package, to
   compute the BSE kernel matrix elements on a coarse k-point grid; and
 - The [`absorption`](absorption-overview.md) executable from the BerkeleyGW
   package, to interpolate the BSE kernel matrix elements to a fine k-point
   grid, diagonalize the BSE Hamiltonian, and compute the optical absorption
   spectrum.

A typical workflow for this set of calculations with BerkeleyGW, including
input and output files, is summarized in the diagram below:

![Typical BerkeleyGW workflow](BerkeleyGW-workflow.svg "Typical BerkeleyGW workflow")

In addition, there is a package to automate these workflows and simply
calculations using the high-level Python language,
[`BGWpy`](https://github.com/BerkeleyGW/BGWpy).
