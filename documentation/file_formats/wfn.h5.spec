Group: mf_header
Content: Information from the WFN file.
#======================================================================

	Dataset: versionnumber
	Type: integer
	Rank: 0
	Value: Version number of this particular file.

	Dataset: flavor
	Type: integer
	Rank: 0
	Value: 1 for Real, 2 for CPLX.


	Group: kpoints
	Content: Variables related to the k-point grid.
	#======================================================================

		Dataset: nspin
		Type: integer
		Rank: 0
		Value: number of spins.

		Dataset: nspinor
		Type: integer
		Rank: 0
		Value: number of spinors, 2 for fully relativistic calculations, 1 if not.

		Dataset: nrk
		Type: integer
		Rank: 0
		Value: number of k-points in irreducible Brillouin zone (IBZ).

		Dataset: mnband
		Type: integer
		Rank: 0
		Value: number of bands.

		Dataset: ngkmax
		Type: integer
		Rank: 0
		Value: maximum number of G-vectors over all-kpoints.

		Dataset: ecutwfc
		Type: double
		Rank: 0
		Value: wavefunction cutoff, in Ry.

		Dataset: kgrid
		Type: integer
		Rank: 1
		Dims(1): 3
		Value: k-grid used in mean-field calculation.

		Dataset: shift
		Type: double
		Rank: 1
		Dims(1): 3
		Value: k-grid shift used in mean-field calculation, in units of
		1/kgrid. For instance, if the kgrid is (2,2,2) and kshift is
		(0.5,0.5,0.5), the kpoints should be (0.25, 0.25, 0.25), (0.75,
		0.25, 0.25), etc.

		Dataset: ngk
		Type: integer
		Rank: 1
		Dims(1): nrk
		Value: number of G-vectors at each k-point. 

		Dataset: ifmin
		Type: integer
		Rank: 2
		Dims(1): nrk
		Dims(2): nspin
		Value: lowest occupied band at each k-point. Numbering starts with 1 (not 0).

		Dataset: ifmax
		Type: integer
		Rank: 2
		Dims(1): nrk
		Dims(2): nspin
		Value: highest occupied band at each k-point. Numbering starts with 1 (not 0).

		Dataset: w
		Type: double
		Rank: 1
		Dims(1): nrk
		Value: weight of each k-point. 

		Dataset: rk
		Type: double
		Rank: 2
		Dims(1): 3
		Dims(2): nrk
		Value: k-points.

		Dataset: el
		Type: double
		Rank: 3
		Dims(1): mnband
		Dims(2): nrk
		Dims(3): nspin
		Value: mean-field energies. (In Ry)

		Dataset: occ
		Type: double
		Rank: 3
		Dims(1): mnband
		Dims(2): nrk
		Dims(3): nspin
		Value: occupations, between 0 and 1.


	Group: gspace
	Content: Variables related to the g-space. 
	#======================================================================

		Dataset: ng
		Type: integer
		Rank: 0
		Value: number of G-vectors in the "full" G-space. Note that we
		typically use a much smaller number of G-vectors to represent
		the wave functions and even the charge density. Note that ng
		should be large enough to represent the charge density, i.e.,
		`ng >= product(FFTgrid)`.

		Dataset: ecutrho
		Type: double
		Rank: 0
		Value: charge density cutoff, in Ry.

		Dataset: FFTgrid
		Type: integer
		Rank: 1
		Dims(1): 3
		Value: FFT grid. 

		Dataset: components
		Type: integer
		Rank: 2
		Dims(1): 3
		Dims(2): ng
		Value: G-vectors in the "full" G-space. See comments on ng.


	Group: symmetry
	Content: Variables related to symmetry operations. We use spglib, so
	convention is rotate, then translate.
	#======================================================================

		Dataset: ntran
		Type: integer
		Rank: 0
		Value: number of symmetries.

		Dataset: cell_symmetry
		Type: integer
		Rank: 0
		Value: cell type (0 = cubic, 1 = hexagonal). 

		Dataset: mtrx
		Type: integer
		Rank: 3
		Dims(1): 3
		Dims(2): 3
		Dims(3): ntran
		Value: symmetry matrices. 

		Dataset: tnp
		Type: integer
		Rank: 2
		Dims(1): 3
		Dims(2): ntran
		Value: fractional translations.


	Group: crystal
	Content: Variables related to crystal structure. 
	#======================================================================

		Dataset: celvol
		Type: double
		Rank: 0
		Value: cell volume (a.u.).

		Dataset: recvol
		Type: double
		Rank: 0
		Value: reciprocal lattice volume (a.u.). 

		Dataset: alat
		Type: double
		Rank: 0
		Value: lattice constant (a.u.).

		Dataset: blat
		Type: double
		Rank: 0
		Value: reciprocal lattice constant (a.u.). 

		Dataset: nat
		Type: integer
		Rank: 0
		Value: number of atoms.

		Dataset: avec
		Type: double
		Rank: 2
		Dims(1): 3
		Dims(2): 3
		Value: lattice vectors (unit of alat). In FORTRAN index
		convention, a vector $[xr]$ in cartesian coordinates can be
		obtained from a vector $[xc]$ in crystal coordinates via $[xr]
		= alat [avec] [xc]$.

		Dataset: bvec
		Type: double
		Rank: 2
		Dims(1): 3
		Dims(2): 3
		Value: reciprocal lattice vectors (units of blat). In FORTRAN
		index convention, a k-point $[kr]$ in cartesian coordinates can
		be obtained from a k-point $[kc]$ in crystal coordinates via
		$[kr] = blat [bvec] [kc]$.

		Dataset: adot
		Type: double
		Rank: 2
		Dims(1): 3
		Dims(2): 3
		Value: metric tensor in real space (a.u.). In FORTRAN index
		convention, a the length $l$ of a vector $[xc]$ in crystal
		coordinates is given by $l = [xc]^T [adot] [xc]$.

		Dataset: bdot
		Type: double
		Rank: 2
		Dims(1): 3
		Dims(2): 3
		Value: metric tensor in reciprocal space (a.u.). In FORTRAN
		index convention, a the length $l$ of a k-point $[kc]$ in crystal
		coordinates is given by $l = [kc]^T [bdot] [kc]$.

		Dataset: atyp
		Type: integer
		Rank: 1
		Dims(1): nat
		Value: atomic species (atomic number, e.g. for carbon, `atyp=6`).

		Dataset: apos
		Type: double
		Rank: 2
		Dims(1): 3
		Dims(2): nat
		Value: atomic positions (in cartesian coordinates, in units of alat).


Group: wfns
Content: Wavefunction G-vectors and coefficients. 
#======================================================================

	Dataset: gvecs
	Type: integer
	Rank: 2
	Dims(1): 3
	Dims(2): ngktot (sum of ngk over all k-points)
	Value: G-vectors for each k-point, listed consecutively, i.e. G-vectors
	for kpt 1, G-vectors for k-point 2, etc

	Dataset: coeffs
	Type: double
	Rank: 4
	Dims(1): 1 for REAL, 2 for CPLX
	Dims(2): ngktot 
	Dims(3): nspin*nspinor
	Dims(4): mnband
	Value: Wavefunction coefficients. We have all G-vector components for
	all k-points consecutively for a given spin and band, i.e. all
	G-vectors for k-point 1, all G-vectors for k-point 2, etc.

