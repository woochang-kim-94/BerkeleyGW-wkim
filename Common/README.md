================================================================================
BerkeleyGW: Common directory (developer info)
================================================================================

  Version 2.0   (May, 2018)
  To be announced.

  Version 1.2   (Aug, 2016)
  F. H. da Jornada, J. Deslippe, D. Vigil-Fowler, J. I. Mustafa, T. Rangel,
  F. Bruneval, F. Liu, D. Y. Qiu, D. A. Strubbe, G. Samsonidze, J. Lischner.

--------------------------------------------------------------------------------

Description:

The ~BGW/Common directory contains common files shared by different components of 
the BerkeleyGW package. The most important files are described below:

1. typedefs.f90

This file contains the derived data types that are used throughout 
the code.

2. nrtype.f90

This file contains some general constants and the definitions of 
single- and double-precision real and complex types following the 
convention of the Numerical Recipes book. The definitions of DP 
and DPC are used throughout the code.

It also contains the fundamental physical constants. The units used 
throughout the code are Rydberg atomic units with a few exceptions, 
such as eV in Sigma and Hartree atomic units in EPM. The difference 
between Hartree and Rydberg atomic units is sketched below.

Hartree atomic units:
hbar=1  e=1        m=1    a0=hbar^2/(m*e^2)=1  E=hbar^2/(2*m*a0^2)=0.5

Rydberg atomic units:
hbar=1  e=sqrt(2)  m=0.5  a0=hbar^2/(m*e^2)=1  E=hbar^2/(2*m*a0^2)=1

3. f_defs.h

This file contains preprocessor definitions for the built-in Fortran 
functions. This helps to avoid accidental conversion from double to 
single precision.
