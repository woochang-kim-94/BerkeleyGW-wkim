# This file was generated by make_deps.py. Do not edit this directly.


#Compile dependencies:
algos_epsilon.o algos_epsilon_m.mod: \
	$(COMMON)/algos_common_m.mod \
	$(COMMON)/message_m.mod \
	$(COMMON)/peinfo_m.mod \
	$(COMMON)/push_pop_m.mod
chi_convergence.o chi_convergence_m.mod: \
	$(COMMON)/global_m.mod
chi_summation.o chi_summation_m.mod: \
	$(COMMON)/accel_linalg_m.mod \
	$(COMMON)/accel_memory_m.mod \
	$(COMMON)/blas_m.mod \
	$(COMMON)/elpa_interface_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/inversion_m.mod \
	$(COMMON)/io_utils_m.mod \
	$(COMMON)/misc_m.mod \
	$(COMMON)/scalapack_m.mod \
	$(COMMON)/timing_m.mod \
	$(COMMON)/vcoul_generator_m.mod \
	algos_epsilon_m.mod \
	lin_denominator_m.mod \
	mtxelmultiply_m.mod
eps0sym.o : \
	$(COMMON)/epsread_hdf5_m.mod \
	$(COMMON)/epswrite_hdf5_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/write_matrix_m.mod
epsascbin.o : \
	$(COMMON)/global_m.mod
epsbinasc.o : \
	$(COMMON)/global_m.mod
epsilon_main.o : \
	$(COMMON)/accel_fft_m.mod \
	$(COMMON)/blas_m.mod \
	$(COMMON)/epsread_hdf5_m.mod \
	$(COMMON)/epswrite_hdf5_m.mod \
	$(COMMON)/fftw_m.mod \
	$(COMMON)/fullbz_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/gmap_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/io_utils_m.mod \
	$(COMMON)/irrbz_m.mod \
	$(COMMON)/read_matrix_m.mod \
	$(COMMON)/references_m.mod \
	$(COMMON)/scalapack_m.mod \
	$(COMMON)/sort_m.mod \
	$(COMMON)/subgrp_m.mod \
	$(COMMON)/tile_m.mod \
	$(COMMON)/timing_m.mod \
	$(COMMON)/vcoul_generator_m.mod \
	$(COMMON)/write_matrix_m.mod \
	$(COMMON)/write_program_header_m.mod \
	algos_epsilon_m.mod \
	chi_convergence_m.mod \
	chi_summation_m.mod \
	epsinv_m.mod \
	genwf_eps_m.mod \
	input_m.mod \
	input_q_m.mod \
	mtxel_m.mod \
	mtxelmultiply_m.mod \
	rqstar_m.mod
epsinv.o epsinv_m.mod: \
	$(COMMON)/global_m.mod \
	$(COMMON)/inversion_m.mod \
	$(COMMON)/io_utils_m.mod \
	$(COMMON)/misc_m.mod \
	$(COMMON)/scalapack_m.mod \
	$(COMMON)/timing_m.mod \
	$(COMMON)/vcoul_generator_m.mod \
	$(COMMON)/write_matrix_m.mod
epsinvomega.o : \
	$(COMMON)/epsread_hdf5_m.mod \
	$(COMMON)/global_m.mod
epsmat_hdf5_upgrade.o : \
	$(COMMON)/epswrite_hdf5_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/hdf5_io_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/wfn_io_hdf5_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod \
	$(COMMON)/write_matrix_m.mod
epsmat_intp.o : \
	$(COMMON)/fullbz_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/sort_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod \
	intp_utils_m.mod
epsmat_merge.o : \
	$(COMMON)/global_m.mod
epsmat_old2hdf5.o : \
	$(COMMON)/epswrite_hdf5_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/vcoul_generator_m.mod \
	$(COMMON)/wfn_io_hdf5_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod \
	$(COMMON)/write_matrix_m.mod
epsomega.o : \
	$(COMMON)/epsread_hdf5_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/inversion_m.mod \
	$(COMMON)/read_matrix_m.mod \
	$(COMMON)/scalapack_m.mod
genwf_eps.o genwf_eps_m.mod: \
	$(COMMON)/fftw_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/gmap_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/sort_m.mod \
	$(COMMON)/susymmetries_m.mod \
	$(COMMON)/timing_m.mod \
	genwf_mpi_m.mod
genwf_mpi.o genwf_mpi_m.mod: \
	$(COMMON)/find_kpt_match_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/gmap_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/misc_m.mod \
	$(COMMON)/sort_m.mod \
	$(COMMON)/susymmetries_m.mod \
	$(COMMON)/timing_m.mod
input.o input_m.mod: \
	$(COMMON)/checkbz_m.mod \
	$(COMMON)/createpools_m.mod \
	$(COMMON)/epswrite_hdf5_m.mod \
	$(COMMON)/eqpcor_m.mod \
	$(COMMON)/fftw_m.mod \
	$(COMMON)/fullbz_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/hdf5_io_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/io_utils_m.mod \
	$(COMMON)/irrbz_m.mod \
	$(COMMON)/misc_m.mod \
	$(COMMON)/scissors_m.mod \
	$(COMMON)/sort_m.mod \
	$(COMMON)/subgrp_m.mod \
	$(COMMON)/timing_m.mod \
	$(COMMON)/wfn_io_hdf5_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod \
	inread_m.mod
input_q.o input_q_m.mod: \
	$(COMMON)/eqpcor_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/hdf5_io_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/io_utils_m.mod \
	$(COMMON)/misc_m.mod \
	$(COMMON)/scissors_m.mod \
	$(COMMON)/timing_m.mod \
	$(COMMON)/wfn_io_hdf5_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod
inread.o inread_m.mod: \
	$(COMMON)/global_m.mod \
	$(COMMON)/inread_common_m.mod \
	$(COMMON)/references_m.mod \
	$(COMMON)/scissors_m.mod \
	algos_epsilon_m.mod
intp_utils.o intp_utils_m.mod: \
	$(COMMON)/global_m.mod
lin_denominator.o lin_denominator_m.mod: \
	$(COMMON)/global_m.mod \
	$(COMMON)/sort_m.mod
mtxel.o mtxel_m.mod: \
	$(COMMON)/accel_fft_m.mod \
	$(COMMON)/accel_memory_m.mod \
	$(COMMON)/fftw_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/misc_m.mod \
	$(COMMON)/timing_m.mod \
	algos_epsilon_m.mod \
	lin_denominator_m.mod
mtxelmultiply.o mtxelmultiply_m.mod: \
	$(COMMON)/accel_linalg_m.mod \
	$(COMMON)/accel_memory_m.mod \
	$(COMMON)/blas_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/scalapack_m.mod \
	$(COMMON)/timing_m.mod \
	algos_epsilon_m.mod
ploteps.o : \
	$(COMMON)/global_m.mod
printchi.o : \
	$(COMMON)/global_m.mod
rqstar.o rqstar_m.mod: \
	$(COMMON)/global_m.mod \
	$(COMMON)/misc_m.mod
setup_subsampling_nns.o : \
	$(COMMON)/checkbz_m.mod \
	$(COMMON)/fullbz_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/irrbz_m.mod \
	$(COMMON)/random_m.mod \
	$(COMMON)/subgrp_m.mod \
	$(COMMON)/wfn_io_hdf5_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod

#Link dependencies:
eps0sym.x eps0sym$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/epsread_hdf5.o \
	$(COMMON)/epswrite_hdf5.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/io_utils.o \
	$(COMMON)/message.o \
	$(COMMON)/nrtype.o \
	$(COMMON)/os.o \
	$(COMMON)/peinfo.o \
	$(COMMON)/push_pop.o \
	$(COMMON)/scalapack.o \
	$(COMMON)/timing.o \
	$(COMMON)/wfn_io_hdf5.o \
	$(COMMON)/write_matrix.o \
	eps0sym.o
epsascbin.x epsascbin$(FLAVOR).x: \
	$(GLOBALOBJS) \
	epsascbin.o
epsbinasc.x epsbinasc$(FLAVOR).x: \
	$(GLOBALOBJS) \
	epsbinasc.o
epsilon.x epsilon$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/accel_fft.o \
	$(COMMON)/accel_linalg.o \
	$(COMMON)/accel_memory.o \
	$(COMMON)/algos_common.o \
	$(COMMON)/bessel.o \
	$(COMMON)/blas.o \
	$(COMMON)/check_inversion.o \
	$(COMMON)/checkbz.o \
	$(COMMON)/createpools.o \
	$(COMMON)/elpa_interface.o \
	$(COMMON)/epsread_hdf5.o \
	$(COMMON)/epswrite_hdf5.o \
	$(COMMON)/eqpcor.o \
	$(COMMON)/fft_parallel.o \
	$(COMMON)/fftw.o \
	$(COMMON)/find_kpt_match.o \
	$(COMMON)/fullbz.o \
	$(COMMON)/gmap.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/input_utils.o \
	$(COMMON)/inread_common.o \
	$(COMMON)/inversion.o \
	$(COMMON)/io_utils.o \
	$(COMMON)/irrbz.o \
	$(COMMON)/lapack.o \
	$(COMMON)/message.o \
	$(COMMON)/minibzaverage.o \
	$(COMMON)/misc.o \
	$(COMMON)/nrtype.o \
	$(COMMON)/os.o \
	$(COMMON)/peinfo.o \
	$(COMMON)/push_pop.o \
	$(COMMON)/random.o \
	$(COMMON)/read_matrix.o \
	$(COMMON)/references.o \
	$(COMMON)/scalapack.o \
	$(COMMON)/scissors.o \
	$(COMMON)/sort.o \
	$(COMMON)/splines.o \
	$(COMMON)/subgrp.o \
	$(COMMON)/susymmetries.o \
	$(COMMON)/tile.o \
	$(COMMON)/timing.o \
	$(COMMON)/trunc_cell_box.o \
	$(COMMON)/trunc_cell_box_d.o \
	$(COMMON)/trunc_cell_wire.o \
	$(COMMON)/trunc_scell_box_d.o \
	$(COMMON)/typedefs.o \
	$(COMMON)/vcoul_generator.o \
	$(COMMON)/version.o \
	$(COMMON)/wfn_io_hdf5.o \
	$(COMMON)/wfn_rho_vxc_io.o \
	$(COMMON)/write_matrix.o \
	$(COMMON)/write_program_header.o \
	algos_epsilon.o \
	chi_convergence.o \
	chi_summation.o \
	epsilon_main.o \
	epsinv.o \
	genwf_eps.o \
	genwf_mpi.o \
	input.o \
	input_q.o \
	inread.o \
	lin_denominator.o \
	mtxel.o \
	mtxelmultiply.o \
	rqstar.o
epsinvomega.x epsinvomega$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/epsread_hdf5.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/message.o \
	$(COMMON)/nrtype.o \
	$(COMMON)/os.o \
	$(COMMON)/peinfo.o \
	$(COMMON)/push_pop.o \
	$(COMMON)/scalapack.o \
	$(COMMON)/timing.o \
	epsinvomega.o
epsmat_hdf5_upgrade.x epsmat_hdf5_upgrade$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/blas.o \
	$(COMMON)/check_inversion.o \
	$(COMMON)/epswrite_hdf5.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/input_utils.o \
	$(COMMON)/io_utils.o \
	$(COMMON)/message.o \
	$(COMMON)/nrtype.o \
	$(COMMON)/os.o \
	$(COMMON)/peinfo.o \
	$(COMMON)/push_pop.o \
	$(COMMON)/scalapack.o \
	$(COMMON)/sort.o \
	$(COMMON)/timing.o \
	$(COMMON)/wfn_io_hdf5.o \
	$(COMMON)/wfn_rho_vxc_io.o \
	$(COMMON)/write_matrix.o \
	epsmat_hdf5_upgrade.o
epsmat_intp.x epsmat_intp$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/blas.o \
	$(COMMON)/check_inversion.o \
	$(COMMON)/fullbz.o \
	$(COMMON)/misc.o \
	$(COMMON)/sort.o \
	$(COMMON)/wfn_rho_vxc_io.o \
	epsmat_intp.o \
	intp_utils.o
epsmat_merge.x epsmat_merge$(FLAVOR).x: \
	$(GLOBALOBJS) \
	epsmat_merge.o
epsmat_old2hdf5.x epsmat_old2hdf5$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/bessel.o \
	$(COMMON)/blas.o \
	$(COMMON)/check_inversion.o \
	$(COMMON)/epswrite_hdf5.o \
	$(COMMON)/fft_parallel.o \
	$(COMMON)/fftw.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/input_utils.o \
	$(COMMON)/io_utils.o \
	$(COMMON)/message.o \
	$(COMMON)/minibzaverage.o \
	$(COMMON)/misc.o \
	$(COMMON)/nrtype.o \
	$(COMMON)/os.o \
	$(COMMON)/peinfo.o \
	$(COMMON)/push_pop.o \
	$(COMMON)/random.o \
	$(COMMON)/scalapack.o \
	$(COMMON)/sort.o \
	$(COMMON)/timing.o \
	$(COMMON)/trunc_cell_box.o \
	$(COMMON)/trunc_cell_box_d.o \
	$(COMMON)/trunc_cell_wire.o \
	$(COMMON)/trunc_scell_box_d.o \
	$(COMMON)/vcoul_generator.o \
	$(COMMON)/wfn_io_hdf5.o \
	$(COMMON)/wfn_rho_vxc_io.o \
	$(COMMON)/write_matrix.o \
	epsmat_old2hdf5.o
epsomega.x epsomega$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/blas.o \
	$(COMMON)/epsread_hdf5.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/input_utils.o \
	$(COMMON)/inversion.o \
	$(COMMON)/lapack.o \
	$(COMMON)/message.o \
	$(COMMON)/nrtype.o \
	$(COMMON)/os.o \
	$(COMMON)/peinfo.o \
	$(COMMON)/push_pop.o \
	$(COMMON)/read_matrix.o \
	$(COMMON)/scalapack.o \
	$(COMMON)/timing.o \
	epsomega.o
ploteps.x ploteps$(FLAVOR).x: \
	$(GLOBALOBJS) \
	ploteps.o
printchi.x printchi$(FLAVOR).x: \
	$(GLOBALOBJS) \
	printchi.o
setup_subsampling_nns.x setup_subsampling_nns$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/blas.o \
	$(COMMON)/check_inversion.o \
	$(COMMON)/checkbz.o \
	$(COMMON)/fullbz.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/irrbz.o \
	$(COMMON)/misc.o \
	$(COMMON)/random.o \
	$(COMMON)/sort.o \
	$(COMMON)/subgrp.o \
	$(COMMON)/wfn_io_hdf5.o \
	$(COMMON)/wfn_rho_vxc_io.o \
	setup_subsampling_nns.o
