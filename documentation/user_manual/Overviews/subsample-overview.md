# Subsampling (NNS and CSI)

GW and BSE calculations are considerably more computationally demanding when
applied to systems with reduced dimensionality, since the electronic
confinement leads to a slower convergence of sums over the Brillouin zone due
to a much more complicated screening environment that manifests in the 'head'
and 'neck' elements of the dielectric matrix. BerkeleyGW implements two schemes
to speedup GW and BSE calculations for systems with reduced dimensionality: the
nonuniform neck subsampling method (NNS) and the clustered sampling
interpolation (CSI) method, which can lead to orders-of-magnitude savings in
computer time.  Please see [PRB 95, 035109
(2017)](https://doi.org/10.1103/PhysRevB.95.035109) for further information,
and please cite the paper if use this functionality!


## Mean-field

The first step before running NNS and CSI is to perform a ground-state DFT
calculation on a relatively coarse k-point grid that is enough to sample
different transitions in the Brillouin zone. This is roughly the same k-grid
that one would use to converge a DFT calculation. For instance, for monolayer
MoS2, a uniform, Gamma-centered `6x6x1` k-point grid is enough. Make sure you
generate a number of unoccupied states necessary to converge your GW
calculation. We will refer to the generated wavefunction file as `WFN`. The
subsampling scheme both NNS and CSI) will then automatically generate a set of
shifted wavefunction files which will allow you to compute the dielectric
matrix at arbitrarily small wavevectors $q$. Note that these shifted calculations
will only involve occupied states, so they are fairly inexpensive.

!!! Note
    For now, the script only automatically generates k-points for Quantum
    ESPRESSO. However, it is fairly straightforward to convert the list of
    k-points to be used with another mean-field code.


## NNS and CSI

After you perform the mean-field calculation, follow the documentation to
perform an [NNS calculation](NNS.md) to obtain the electronic self-energy of a
quasi-2D system, or an [CSI calculation](CSI.md) to obtain the optical absorption
spectrum.
