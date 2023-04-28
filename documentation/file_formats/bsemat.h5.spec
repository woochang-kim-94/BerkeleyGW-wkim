Group: mf_header
Content: Header from the [WFN file](../wfn_h5_spec).
#===============================================================================

Group: eps_header
Content: Header from the [epsmat file](../epsmat_h5_spec).
#===============================================================================

Group: bse_header
Content: Everything but the kernel matrix elements.
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

		Dataset: screening
		Type: integer
		Rank: 0
		Value: What type of screening was used.

		Dataset: icutv
		Type: integer
		Rank: 0
		Value: Truncation flag. See Epsilon/inread.f90 for the meaning of each icutv.

		Dataset: ecuts
		Type: double
		Rank: 0
		Value: Cutoff for screened Coulomb interaction.

		Dataset: ecutg
		Type: double
		Rank: 0
		Value: Cutoff for the bare Coulomb interaction.

		Dataset: efermi
		Type: double
		Rank: 0
		Value: Fermi level determined by the code, in eV.

		Dataset: theory
		Type: integer
		Rank: 0
		Value: 0 for BSE, 1 for TDDFT.

		Dataset: nblocks
		Type: integer
		Rank: 0
		Value: How many transitions blocks are there in the kernel matrix?
		Possible values are:
		- 1 for restricted TDA kernel: vc -> v`c`
		- 2 for restricted non-TDA kernel: {vc,cv} -> {v`c`,c`v`}  [not implemented]
		- 4 for extended kernel: {n1,n2} -> {n1`,n2`}

		Dataset: storage
		Type: integer
		Rank: 0
		Value: How we store the kernel matrix. Possible values are:
		0 for no symmetries, i.e., {n1,n2} -> {n1`,n2`}

		Dataset: nmat
		Type: integer
		Rank: 0
		Value: Number of matrices stored in /mats

		Dataset: energy_loss
		Type: integer
		Rank: 0
		Value: 1 if this an energy-loss calculation. Otherwise, 0.


	Group: bands
	Content: Bands-related datasets.
	#===============================================================================

		Dataset: nvb
		Type: integer
		Rank: 0
		Value: Number of valence bands included in the kernel computation.

		Dataset: ncb
		Type: integer
		Rank: 0
		Value: Number of conduction bands included in the kernel computation.

		Dataset: n1b
		Type: integer
		Rank: 0
		Value: Number of bands of first kind we loop over in the kernel matrix.
		This is nvb for restricted kernels, and nvb+ncb for extended kernels.

		Dataset: n2b
		Type: integer
		Rank: 0
		Value: Number of bands of second kind we loop over in the kernel matrix.
		This is ncb for restricted kernels, and nvb+ncb for extended kernels.

		Dataset: ns
		Type: integer
		Rank: 0
		Value: Number of spins. This is 1 unless this is a spin-polarized
		calculation with a single spinor component, on which case this is 2.

		Dataset: nspinor
		Type: integer
		Rank: 0
		Value: Number of spinor components. Either 1 or 2.


	Group: kpoints
	Content: K-points-related datasets.
	#===============================================================================

		Dataset: nk
		Type: integer
		Rank: 0
		Value: Number of k-points in the full zone.

		Dataset: kpts
		Rank: 2
		Dims(1): 3 #(the three crystal coordinates)
		Dims(2): nk
		Value: K-points.

		Dataset: kgrid
		Type: integer
		Rank: 1
		Dims(1): 3
		Value: K-grid used in epsilon calculation.

		Dataset: qflag
		Type: integer
		Rank: 0
		Value: 0 for finite Q calculation with arb. Q (deprecated)
		       1 for Q=0 calculation (default)
		       2 for Q commensuare with WFN_co (under development)

		Dataset: exciton_Q_shift
		Type: double
		Rank: 1
		Dims(1): 3
		Value: Q-shift for finite Q calculation, note the shift is applied to 
		the valence states, and the excitons computed with this kernel have
		center of mass moment equal to MINUS Q.

		Dataset: patched_sampling
		Type: logical
		Rank: 0
		Value: Was this bsemat calculated on a patch?

Group: mats
Content: Matrix-elements-related datasets.
#===============================================================================

	Dataset: head, wing, body, exchange, fxc
	Type: double
	Rank: 6
	Dims(1): flavor
	Dims(2:3): n1b
	Dims(4:5): n2b
	Dims(6:7): nk*ns
	Value: Kernel matrix elements times V/8*Pi * Ry. The indices are
	stored as follows:
	\n\n
	- For restricted kernels (default) n1b = nv and n2b = nc. The
	  *valence* bands (dims 2 and 3) are counted *downwards* from the
	  Fermi level, and the *conduction* bands (dims 4 and 5) are counted
	  *upwards* from the Fermi level. So, restricted kernels separate
	  valence-band from conduction-band indices. Without spin:
	  \n
	  $K^d(v,v',c,c',k,k') =$ $-\int d^3r d^3r' \;
	      \psi_{ck}(r)^* \psi_{c'k'}(r) W(r,r')
	      \psi_{vk}(r') \psi_{v'k'}(r')^*$
	  \n
	  $K^x(v,v',c,c',k,k') =$ $\int d^3r d^3r' \;
	     \psi_{ck}(r)^* \psi_{vk}(r) v(r-r')
	     \psi_{c'k'}(r') \psi_{v'k'}(r')^*$
	\n\n
	- For extended kernels, n1b = n2b = nv + nc, so all bands (dims 2 to
	  5) are treated on the same footing. The first nv indices refer to
	  valence states (counted downwards from the Fermi level), and the
	  following nc indices refer to conduction states (counted upwards
	  from the Fermi level). Without spin:
	  \n
	  $K^d(n,n',m,m',k,k') =$ $- \int d^3r d^3r' \;
	     \psi_{mk}(r)^* \psi_{m'k'}(r) W(r,r')
	     \psi_{nk}(r') \psi_{n'k'}(r')^*$
	  \n
	  $K^x(v,v',c,c',k,k') =$ $\int d^3r d^3r' \;
	     \psi_{mk}(r)^* \psi_{nk}(r) v(r-r')
	     \psi_{m'k'}(r') \psi_{n'k'}(r')^*$
	  \n
