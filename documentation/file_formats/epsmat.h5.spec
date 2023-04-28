Group: mf_header
Content: Header from the [WFN file](../wfn_h5_spec).
#===============================================================================

Group: eps_header
Content: Everything but the dielectric matrices.
#===============================================================================

	Dataset: versionnumber
	Type: integer
	Rank: 0
	Value: Version number of this particular file.

	Dataset: flavor
	Type: integer
	Rank: 0
	Value: 1 for Real, 2 for CPLX.


	Group: params
	Content: General parameters.
	#===============================================================================

		Dataset: matrix_type
		Type: integer
		Rank: 0
		Value: Controls what matrix we compute: 0=epsilon^{-1} (epsmat.h5), 
		1=epsilon (not implemented), 2=chi0 (chimat.h5).

		Dataset: has_advanced
		Type: integer/logical
		Rank: 0
		Value: Whether we store the advanced matrix. This is the default in
		BerkeleyGW 1.1.x, but only for FF calculations with CPLX flavor.

		Dataset: nmatrix
		Type: integer
		Rank: 0
		Value: Number of matrices we store. If we store dielectric matrices,
		this is (has_advanced+1). For polarizability calculations, this is
		(has_advanced+1)*nspin.

		Dataset: matrix_flavor
		Type: integer
		Rank: 0
		Value: Whether we represent the dielectric matrix with real or complex numbers.
		This will be 2 (i.e., we store complex numbers), unless we have a static
		calculation run with the real code.

		Dataset: icutv
		Type: integer
		Rank: 0
		Value: Truncation flag. See Epsilon/inread.f90 for the meaning of each icutv.

		Dataset: ecuts
		Type: double
		Rank: 0
		Value: Cutoff for screened Coulomb interaction.

		Dataset: nband
		Type: integer
		Rank: 0
		Value: Total number of bands included in the epsilon calculation.

		Dataset: efermi
		Type: double
		Rank: 0
		Value: Fermi level determined by the code, in ryd.

#BEGIN_INTERNAL_ONLY
		Dataset: intraband_flag
		Type: integer
		Rank: 0
		Value: See Epsilon/epsilon.inp.

		Dataset: intraband_overlap_min
		Type: double
		Rank: 0
		Value: See Epsilon/epsilon.inp.
#END_INTERNAL_ONLY

		Dataset: subsampling
		Type: integer/logical
		Rank: 0
		Value: See Epsilon/epsilon.inp.

		Dataset: subspace
		Type: logical
		Rank: 0
		Value: true if the epsilon matrices have been generated using
		the static subspace approximation (only works in combination
		with full-frequency contour-deformation approach)


	Group: qpoints
	Content: Q-points-related datasets
	#===============================================================================

		Dataset: nq
		Type: integer
		Rank: 0
		Value: Number of q-points.

		Dataset: qpts
		Rank: 2
		Dims(1): 3 #(the three crystal coordinates)
		Dims(2): nq
		Value: Q-points.

		Dataset: qgrid
		Type: integer
		Rank: 1
		Dims(1): 3
		Value: Q-grid used in epsilon calculation.

		Dataset: qpt_done
		Type: integer/logical
		Rank: 1
		Dims(1): nq
		Value: 1 is the calculation for a particular q-point is done, 0 if not.


	Group: freqs
	Content: Frequency-related datasets
	#===============================================================================

		Dataset: freq_dep
		Type: integer
		Rank: 0
		Value: Frequency dependency (same as in epsilon.inp and sigma.inp).

		Dataset: nfreq
		Type: integer
		Rank: 0
		Value: Number of frequency points.

		Dataset: nfreq_imag
		Type: integer
		Rank: 0
		Value: Number of imaginary frequency points. Used for CD calculations.

		Dataset: freqs
		Type: double
		Rank: 2
		Dims(1): 2 #(Real/complex part)
		Dims(2): nfreq
		Value: Frequencies, including broadening.


	Group: gspace
	Content: G-vectors-related datasets
	#===============================================================================

		Dataset: nmtx
		Type: integer
		Rank: 1
		Dims(1): nq
		Value: Number of matrix elements we actually compute for each q-point.

		Dataset: nmtx_max
		Type: integer
		Rank: 0
		Value: Same as `maxval(nmtx(1:nq))`.

		Dataset: ekin
		Type: double
		Rank: 2
		Dims(1): ng
		Dims(2): nq
		Value: Kinetic energies $|G+q|^2$, in Rydberg. Note: the
		G-vectors are sorted wrt the RHO G-space $|G|^2$, and not wrt
		epsilon G-space $|G+q|^2$!

		Dataset: gind_eps2rho
		Type: integer
		Rank: 2
		Dims(1): ng
		Dims(2): nq
		Value: Given a row/column ig_eps from the epsilon matrix (in
		the epsilon G-space), map_eps2rho(ig_eps,iq) is the index of
		the corresponding G-vector in the RHO G-space.

		Dataset: gind_rho2eps
		Type: integer
		Rank: 2
		Dims(1): ng
		Dims(2): nq
		Value: Given a G-vector ig from the RHO G-space,
		map_rho2eps(ig,iq) is the corresponding row/column in the
		epsilon matrix (i.e., the index in the epsilon G-space).

	Group: subspace
	Content: Static subspace approximation (SSA) related parameters, it is
	created only if the calcuation is performed within the static subspace
	approximation (this forces to check the existence of the group and thus
	help to keep portability for .h5 matrices generated with older versions
	of the code). Again this works (for now) only in combination with
	full-frequency contour-deformation approach
	#===============================================================================

		Dataset: keep_full_eps_static
		Type: logical
		Rank: 0
		Value: Set to true if the inverse static dielectric matrix has
		to be retained in the eps0mat.h5 and epsmat.h5

		Dataset: matrix_in_subspace_basis
		Type: logical
		Rank: 0
		Value: Set to true if the inverse dielectric matrices are stored in the subspace basis

		Dataset: eps_eigenvalue_cutoff
		Type: double
		Rank: 0
		Value: Set the cutoff for the selection of the eigenvectors
		which will build the static subspace (namely all eigenvectors
		with eigenvalue larger than eps_eigenvalue_cutoff)
  
		Dataset: neig_max
		Type: integer
		Rank: 0
		Value: Define the maximum size (number of eigenvectors) that
		will be used for the generation of the static subspace, can be
		given in input (see nbasis_subspace in epsilon.inp) otherwise
		it will be set automatically to 20% of the max nmtx.

		Dataset: neig
		Rank: 1
		Dims(1): nq
		Value: give the actual number of eigenvectors employed for the
		generation of the static subspace for each q-point



Group: mats
Content: Matrix-elements-related datasets
#===============================================================================

	Dataset: matrix
	Type: double
	Rank: 6
	Dims(1): matrix_flavor
	Dims(2): nmtx_max #rows of the matrix
	Dims(3): nmtx_max #columns of the matrix
	Dims(4): nfreq
	Dims(5): nmatrix
	Dims(6): nq
	Value: Matrix elements (see matrix_type dataset). Rows/columns are
	sorted wrt epsilon(q) G-space, $|q+G|^2$. Note: for each q-point, we only
	really compute values up to nmtx(q).

	Dataset: matrix-diagonal
	Type: double
	Rank: 3
	Dims(1): matrix_flavor
	Dims(2): nmtx_max #diagonal elements of the matrix
	Dims(3): nq
	Value: Static diagonal elements from "matrix" dataset. Not used when
	matrix_type==2.

	Dataset: matrix_subspace
	Type: double
	Rank: 6
	Dims(1): matrix_flavor
	Dims(2): neig_max #rows of the matrix
	Dims(3): neig_max #columns of the matrix
	Dims(4): nfreq
	Dims(5): nmatrix
	Dims(6): nq
	Value: Matrix elements of the inverse dielectric matrix in the subspace
	basis, for each q-point, we only really compute and store values up to
	neig(q)

	Dataset: matrix_eigenvec
	Type: double
	Rank: 6
	Dims(1): matrix_flavor
	Dims(2): nmtx_max #rows of the matrix
	Dims(3): neig_max #columns of the matrix
	Dims(4): 1	#frequency independent 
	Dims(5): nmatrix
	Dims(6): nq
	Value: Eigenvectors of the static dielectric matrix for each q-point,
	we only really compute and store the rows/cols to nmtx(q)/neig(q)
	(defined with rank 6 so we can reuse the same routines used for the
	previous arrays)

	Dataset: matrix_fulleps0
	Type: double
	Rank: 6
	Dims(1): matrix_flavor
	Dims(2): nmtx_max #rows of the matrix
	Dims(3): nmtx_max #columns of the matrix
	Dims(4): 1	#frequency independent 
	Dims(5): nmatrix
	Dims(6): nq
	Value: Static dielectric matrix for each q-point (if
	keep_full_eps_static is set to true), we only really compute and store
	the rows/cols to nmtx(q) (defined with rank 6 so we can reuse the same
	routines used for the previous arrays)
