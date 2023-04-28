
Mandatory keywords are uncommented and preceded by documentation,
while optional keywords are commented:
    

    # Energy cutoff for the dielectric matrix, in Ry. The dielectric matrix
    # $\varepsilon_{GG'}$ will contain all G-vectors with kinetic energy $|q+G|^2$
    # up to this cutoff.
    #[unit=Ry]
    epsilon_cutoff           35.0

    # Total number of bands (valence+conduction) to sum over. Defaults to the
    # number of bands in the WFN file minus 1.
    #number_bands            1000


An empty commented line will allow you to write a more structured block:

    # This flags specifies the frequency dependence of the inverse dielectric matrix:
    #
    # - Set to 0 to compute the static inverse dielectric matrix (default).
    # - Set to 2 to compute the full frequency dependent inverse dielectric matrix.
    # - Set to 3 to compute the two frequencies needed for Godby-Needs GPP model.
    #frequency_dependence 0


Input variables that are closely related should be linked by commented lines:

    # Increase in the frequency step for the non-uniform frequency grid.
    #delta_frequency_step 1.0
    #
    # Frequency step for the linear grid for the spectral function method.
    # Defaults to [[delta_frequency]].
    #delta_sfrequency 0.1
    #
    # Increase in frequency step for the non-uniform grid for the spectral function method.
    # Defaults to [[delta_frequency_step]]
    #delta_sfrequency_step 1.0


Two empty commented lines will allow you to write a block for several input variables.

    # Symmetry specification                                                                                                                                         
    #
    #
    #no_symmetries_coarse_grid
    #use_symmetries_coarse_grid                                                     
    #                                                                               
    #no_symmetries_fine_grid
    #use_symmetries_fine_grid


Enclose internal features:

    #BEGIN_INTERNAL_ONLY

    # This keyword will not appear in the assembled manual.
    #secret_keyword

    #END_INTERNAL_ONLY


Use include statements for variables that are common to several codes:

    %include common.inp                                                             
    %include scissors.inp


You can reference other variables:

    # This variable is used in conjunction with [[two_photon_job]].
    #number_of_final_states 500
