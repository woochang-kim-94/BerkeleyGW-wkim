Group: mf_header
Content: Header from the [WFN file](../wfn_h5_spec).
#===============================================================================

Group: eps_header
Content: Header from the [epsmat file](../epsmat_h5_spec).
#===============================================================================

Group: bse_header
Content: Header from the [bsemat file](../bsemat_h5_spec).
#===============================================================================

Group: exciton_header
Content: All the dimensions and scalar quantities
#===============================================================================

	Dataset: version
	Type: integer
	Rank: 0
	Value: Version number of this particular file.

	Dataset: flavor
	Type: integer
	Rank: 0
	Value: 1 if arrays are real, 2 if arrays are complex


	Group: params
	Content: General parameters
	#=============================================================================
	
		Dataset: bse_hamiltonian_size
		Type: integer
		Rank: 0
		Value: Size of the BSE Hamiltonian, or the A sub-block case of
		non-TDA calculation, that is, ns * nk * nv * nc.
		
		Dataset: evec_sz
		Type: integer
		Rank: 0
		Value: Size of the flat eigenvectors. Depends on whether
		Tamm-Dancoff approximation is used. If TDA is used then it is
		equal to bse_hamiltonian_size, otherwise it is twice that value
		
		Dataset: spin_kernel
		Type: integer
		Rank: 0
		Value: 0 for spin-triplet, 1 for spin-singlet, 2 for
		local-field, 3 for spinor calculations, respectively.
		
		Dataset: nevecs
		Type: integer
		Rank: 0
		Value: The number of eigenvectors present in the file
		
		Dataset: ns
		Type: integer
		Rank: 0
		Value: The number of spins.
		
		Dataset: nc
		Type: integer
		Rank: 0
		Value: The number of conduction bands
		
		Dataset: nv
		Type: integer
		Rank: 0
		Value: The number of valence bands
		
		Dataset: use_tda
		Type: integer/logical
		Rank: 0
		Value: 1 TDA is used, 0 otherwise
	
	
	Group: kpoints
	Content: K-points-related datasets.
	#=============================================================================
	
		Dataset: nk
		Type: integer
		Rank: 0
		Value: The number of k points.
		
		Dataset: kpts
		Type: double
		Rank: 2
		Dims(1): 3
		Dims(2): nk
		Value: The list of k-points.
		
		Dataset: nQ
		Type: integer
		Rank: 0
		Value: The number of Q points.
		
		Dataset: exciton_Q_shifts
		Type: double
		Rank: 2
		Dims(1): 3
		Dims(2): nQ
		Value: Shifted wavevectors used for finite Q momentum. Note the shifts
		are applied to the valence state, and the resulting exciton has center 
		mass momentum MINUS Q.


Group: exciton_data
Content: All the arrays
#===============================================================================

	Dataset: eigenvalues
	Type: double
	Rank: 1
	Dims(1): number_or_eigenvectors
	Value: The BSE eigenvalues
	
	Dataset: eigenvectors
	Type: double
	Rank: 5 or 6
	Dims(1): 2 (if flavor==2). Otherwise, this dimension is suppressed.
	Dims(2): ns
	Dims(3): nv
	Dims(4): nc
	Dims(5): nk
	Dims(6): nevecs
	Dims(7): nQ
	Value: The electron-hole coefficients, that is, the right BSE
	eigenvectors Avc. This array is used when TDA is true or false.
	
	Dataset: eigenvectors_left
	Type: double
	Rank: 5 or 6
	Dims(1): 2 (if flavor==2). Otherwise, this dimension is suppressed.
	Dims(2): ns
	Dims(3): nv
	Dims(4): nc
	Dims(5): nk
	Dims(6): nevecs
	Dims(7): nQ
	Value: The electron-hole coefficients, that is, the left BSE
	eigenvectors Avc. This array is used when TDA is false.
	
	Dataset: eigenvectors_deexcitation
	Type: double
	Rank: 5 or 6
	Dims(1): 2 (if flavor==2). Otherwise, this dimension is suppressed.
	Dims(2): ns
	Dims(3): nv
	Dims(4): nc
	Dims(5): nk
	Dims(6): nevecs
	Dims(7): nQ
	Value: The electron-hole coefficients, that is, the right BSE
	eigenvectors  Bvc with negative eigenvalues. This array is used when
	TDA is false.
	
	Dataset: eigenvectors_deexcitation_left
	Type: double
	Rank: 5 or 6
	Dims(1): 2 (if flavor==2). Otherwise, this dimension is suppressed.
	Dims(2): ns
	Dims(3): nv
	Dims(4): nc
	Dims(5): nk
	Dims(6): nevecs
	Dims(7): nQ
	Value: The electron-hole coefficients, that is, the left BSE
	eigenvectors Bvc with negative eigenvalues. This array is used when TDA
	is false.
