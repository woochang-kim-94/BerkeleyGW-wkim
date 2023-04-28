# Mean-field calculations

The BerkeleyGW software package uses many-body perturbation-theory formalisms;
therefore, one needs to provide a reasonable mean-field starting point for the
perturbation-theory calculations. For most application, density-functional
theory (DFT) based on semilocal functionals provide a good starting point for
GW and GW-BSE calculations.

The following mean-field codes are currently supported by BerkeleyGW:

- [Quantum ESPRESSO](espresso-overview.md)
- [Abinit](abi2bgw-input.md)
- [Octopus](octopus-overview.md)
- [PARATEC](paratec-overview.md)
- [PARSEC](http://parsec.ices.utexas.edu/){:target="_blank"}
- [EPM (Empirical Pseudopotential Method)](epm-overview.md)
- [RMGDFT](rmgdft-overview.md)
- [Siesta](siesta-overview.md)
- [JDFTx](jdftx-overview.md)

In addition, we include [a wrapper](bgw2sgw-overview.md) to convert mean-field
quantities to the [StochasticGW](http://stochasticgw.com) code.

We provide a library to easily write mean-field-related quantities in the
format used by BerkeleyGW.

BerkeleyGW requires the following quantities from mean-field codes:

 - Mean-field eigenvalues and eigenvectors, stored in `WFN` files, for
   arbitrary k-point grids.
 - The mean-field exchange-correlation matrix elements, `vxc.dat`, or
   exchange-correlation matrix in reciprocal space, `VXC`.
 - Ground-state charge density $\rho(G)$, stored in `RHO` -- only for
   calculations based on the generalized plasmon-pole (GPP) model.

For further information on how to use each mean-field code with the appropriate
wrapper for BerkeleyGW, we refer to the BerkeleyGW tutorials. We also
present some notes on potential issues arising when running GW calculations
starting from different mean-field codes in a [separate page](meanfield-details.md).
