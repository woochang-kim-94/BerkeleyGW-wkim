# SIESTA wrapper

SIESTA is a localized-orbital DFT code, available at:
<https://departments.icmab.es/leem/siesta/>
It can be used in conjunction with BerkeleyGW, especially to generate unoccupied states
with the SAPO utility.


!!! Warning
    Make sure you understand when you can use localized basis sets for GW
    calculations. When if doubt, use a more well-tested and validated wrapper
    based on Quantum ESPRESSO, PARATEC, or Abinit.


### `siesta2bgw.x`

An example `siesta2bgw.x` calculation is provided in `examples/DFT/benzene`. 
Please see this example for usage and example siesta and denchar input
files.

The SIESTA wrapper reads lattice vectors, atomic species and
positions from `SystemLabel.XV` file, k-points from `SystemLabel.KP`
file, eigenvalues from `SystemLabel.EIG` file, wavefunctions in R-space
from `SystemLabel{.K$k}.WF$n{.UP|.DOWN}{.REAL|.IMAG}.cube` file, symmetry
group, kinetic energy cutoffs for wavefunctions and charge density,
k-grid dimensions and offsets from the input file, generates reciprocal
lattice vectors, rotation matrices, fractional translations, FFT grid and
G-vectors, performs FFT from R-space to G-space, and writes the selected
range of bands to the output WFN file. Make sure to set up the real
space grid in `SIESTA/DENCHAR` utility same as FFT grid. If the two grids
differ the SIESTA wrapper will return an error. The output WFN file
can be used as an auxiliary wavefunction for the SAPO code.

A full GW/BSE calculation with SIESTA wavefunctions also requires RHO and 
`VXC` files. These files can be generated from `SystemLabel.RHO{.UP|.DN}.cube`, 
`SystemLabel.VT{.UP|.DN}.cube` and `SystemLabel.VH.cube` files obtained with 
`SIESTA/GRID2CUBE` utility on SIESTA real space grid. Make sure SIESTA real 
space grid is equivalent to FFT grid by setting parameter `MeshCutoff` 
in SIESTA input file same as kinetic energy cutoff for charge density.

Input is read from file `siesta2bgw.inp`:

```
&input_siesta2bgw
   systemlabel      = 'ammonia'
   wfng_output_file = 'wfng.lo'
   rhog_output_file = ''
   vxcg_output_file = ''
   ecutwfn          = 60.0
   ecutrho          = 240.0
   wfng_nk1         = 1
   wfng_nk2         = 1
   wfng_nk3         = 1
   wfng_dk1         = 0.0
   wfng_dk2         = 0.0
   wfng_dk3         = 0.0
   wfng_ref_flag    = .true.
   wfng_ref_kpoint  = 1
   wfng_ref_spin    = 1
   wfng_ref_band    = 4
   wfng_ref_energy  = -5.87
   wfng_band_flag   = .true.
   wfng_band_min    = 5
   wfng_band_max    = 28
   wfng_energy_flag = .false.
   wfng_energy_min  = 0.0
   wfng_energy_max  = 0.0
   wfng_gamma_real  = .true.
/
```

Here, lattice vectors, k-points, eigenvalues and wavefunctions in R-space 
are read from `ammonia.XV`, `ammonia.KP`, `ammonia.EIG` and `ammonia.WF$n.cube `
files, respectively. The kinetic-energy cutoffs for 
wavefunction and charge density are set to 60 and 240 Ry, respectively. 
The k-grid is set to the Gamma-point. The SIESTA eigenvalues are shifted 
in energy so that the HOMO state (1st k-point, 1st spin, 4th band) appears 
at -5.87 eV. The real parts of resonant SIESTA wavefunctions (from 5th to 
28th band) are read from Gaussian Cube files and the imaginary parts are 
set to zero (`wfng_gamma_real`), the wavefunctions are brought from R-space 
to G-space and written to `wfng.lo` file. The SIESTA charge density and 
exchange-correlation potential files are not generated.
