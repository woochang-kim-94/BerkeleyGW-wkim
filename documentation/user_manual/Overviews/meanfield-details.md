# Mean-field calculations: tricks and hints


## Tricks and hints

### Vxc0

`Vxc0` is the $G=0$ component of the exchange-correlation potential $V_{xc}$. 
`Vxc0` is added to the eigenvalues in the wavefunction files produced 
by PARATEC and ESPRESSO. `Vxc0` is also included in the `VXC` and `vxc.dat` 
files produced by PARATEC and ESPRESSO (the `VXC` file contains $Vxc(G)$ 
and the `vxc.dat` file contains matrix elements of $V_{xc}$). The two `Vxc0` 
terms cancel out when calculating the quasiparticle corrections 
in the Sigma code.

### Vacuum level

To correct the DFT eigenvalues for the vacuum level take the average 
of the electrostatic potential on the faces of the unit cell in the 
non-periodic directions and subtract it from the DFT eigenvalues. 
The electrostatic potential is defined as $V_{elec} = V_{bare} + V_{hartree}$, 
while the total potential is $V_{tot} = V_{bare} + V_{hartree} + V_{xc}$, hence 
$V_{tot}$ contains `Vxc0` and $V_{elec}$ does not. The average of $V_{elec}$ is 
fed into the Sigma code using keywords `avgpot` and `avgpot_outer`. 
The potentials can be generated with PARATEC and ESPRESSO. 
In PARATEC, use `elecplot` for $V_{elec}$ or `potplot` for $V_{tot}$, then 
convert from `a3dr` to cube using the `Visual/volume.py` script. 
In ESPRESSO, use `iflag=3`, `output_format=6`, and `plot_num=11` 
for $V_{elec}$ or `plot_num=1` for $V_{tot}$. The averaging is done with 
the `Visual/average.py` script.

### Unit cell size

If you truncate the Coulomb interaction, make sure that the size 
of the unit cell in non-periodic directions is at least two times 
larger than the size of the charge distribution. This is needed 
to avoid spurious interactions between periodic replicas but at 
the same time not to alter interactions within the same unit cell. 
Run `Visual/surface.x` to plot an isosurface that contains 99% of the 
charge density (see Visual/README for instructions on how to do this). 
The code will print the size of the box that contains the isosurface 
to stdout. Multiply the box dimensions in non-periodic directions 
by two to get the minimum size of the unit cell.

### Inversion symmetry

When using the real flavor of the code, make sure the inversion 
symmetry has no associated fractional translation (if it does, 
shift the coordinate origin). Otherwise `WFN`, `VXC`, `RHO` (and `VSC` 
in SAPO) have non-vanishing imaginary parts which are simply 
dropped in the real flavor of the code. This won't do any 
good to your calculation.


## Pitfalls

There are some notorious cases where typical DFT calculations might not provide
a good starting point for one-shot perturbation-theory calculations. Examples
include:

 - Germanium crystal, which is often predicted to be metallic at DFT. This
   issue can be remedied by either using another mean-field starting point, or
   performing some sort of self-consistent iteration, for instance, based on
   the static COHSEX approximation.

 - Molecules such as silane (SiH4), in which semilocal DFT often yields a LUMO
   orbital that bound, where in reality it is unbound. These systems can also
   be remedied by using a different starting mean-field point, or performing
   self-consistent GW calculations. For molecules, another common approach is
   known as best $G$, best $W$, where one picks one mean-field starting point
   to write the Green's function $G$ in $\Sigma=iGW$, such as Hartree-Fock, and
   another one to compute the polarizability $\chi^0$ used to construct $W$,
   such as LDA.

 - Some strongly correlated systems, especially those with partially filled $f$
   orbitals. Systems such as transition-metal oxides are often tackled with a
   Hubbard-type of correction scheme, such as DFT+U.
