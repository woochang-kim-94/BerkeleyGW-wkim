# This file was generated by make_deps.py. Do not edit this directly.


#Compile dependencies:
bgw_mpi.o bgw_mpi_m.mod: \
	$(COMMON)/global_m.mod
compare_wfns.o : \
	$(COMMON)/hdf5_io_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/misc_m.mod
diag_driver.o diag_driver_m.mod: \
	$(COMMON)/global_m.mod \
	diag_elpa_m.mod \
	diag_primme_m.mod \
	diag_scalapack_m.mod \
	distribution_m.mod \
	inread_m.mod
diag_elpa.o diag_elpa_m.mod: \
	$(COMMON)/global_m.mod \
	$(COMMON)/scalapack_m.mod \
	distribution_m.mod \
	inread_m.mod
diag_primme.o diag_primme_m.mod: \
	$(COMMON)/blas_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/lapack_m.mod \
	$(COMMON)/scalapack_m.mod \
	diag_elpa_m.mod \
	diag_scalapack_m.mod \
	distribution_m.mod \
	inread_m.mod \
	primme_m.mod
diag_scalapack.o diag_scalapack_m.mod: \
	$(COMMON)/blas_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/lapack_m.mod \
	$(COMMON)/scalapack_m.mod \
	distribution_m.mod \
	inread_m.mod
distribution.o distribution_m.mod: \
	$(COMMON)/blas_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/lapack_m.mod \
	$(COMMON)/scalapack_m.mod \
	kpoint_pool_m.mod
hamiltonian.o hamiltonian_m.mod: \
	$(COMMON)/global_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/io_utils_m.mod \
	$(COMMON)/misc_m.mod \
	diag_driver_m.mod \
	distribution_m.mod \
	inread_m.mod \
	pseudopot_m.mod
inread.o inread_m.mod: \
	$(COMMON)/global_m.mod \
	kpoint_pool_m.mod
iteration_data.o iteration_data_m.mod: \
	$(COMMON)/global_m.mod
kpoint_pool.o kpoint_pool_m.mod: \
	$(COMMON)/global_m.mod
parabands.o : \
	$(COMMON)/global_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/io_utils_m.mod \
	$(COMMON)/misc_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod \
	$(COMMON)/write_program_header_m.mod \
	diag_driver_m.mod \
	distribution_m.mod \
	hamiltonian_m.mod \
	inread_m.mod \
	iteration_data_m.mod \
	pseudopot_m.mod \
	pseudopot_vkb_m.mod \
	wfn_io_m.mod
pseudopot.o pseudopot_m.mod: \
	$(COMMON)/global_m.mod \
	bgw_mpi_m.mod \
	kpoint_pool_m.mod
pseudopot_vkb.o pseudopot_vkb_m.mod: \
	$(COMMON)/global_m.mod \
	$(COMMON)/misc_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod \
	bgw_mpi_m.mod \
	kpoint_pool_m.mod \
	pseudopot_m.mod
split_spin.o : \
	$(COMMON)/global_m.mod \
	$(COMMON)/io_utils_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod \
	inread_m.mod \
	pseudopot_vkb_m.mod \
	wfn_io_m.mod
wfn_io.o wfn_io_m.mod: \
	$(COMMON)/check_inversion_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/hdf5_io_m.mod \
	$(COMMON)/input_utils_m.mod \
	$(COMMON)/misc_m.mod \
	$(COMMON)/wfn_io_hdf5_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod \
	distribution_m.mod \
	inread_m.mod \
	iteration_data_m.mod

#Link dependencies:
compare_wfns.x compare_wfns$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/blas.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/input_utils.o \
	$(COMMON)/misc.o \
	compare_wfns.o
parabands.x parabands$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/blas.o \
	$(COMMON)/check_inversion.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/input_utils.o \
	$(COMMON)/io_utils.o \
	$(COMMON)/lapack.o \
	$(COMMON)/misc.o \
	$(COMMON)/scalapack.o \
	$(COMMON)/sort.o \
	$(COMMON)/version.o \
	$(COMMON)/wfn_io_hdf5.o \
	$(COMMON)/wfn_rho_vxc_io.o \
	$(COMMON)/write_program_header.o \
	bgw_mpi.o \
	diag_driver.o \
	diag_elpa.o \
	diag_primme.o \
	diag_scalapack.o \
	distribution.o \
	hamiltonian.o \
	inread.o \
	iteration_data.o \
	kpoint_pool.o \
	parabands.o \
	primme.o \
	pseudopot.o \
	pseudopot_vkb.o \
	wfn_io.o
split_spin.x split_spin$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/blas.o \
	$(COMMON)/check_inversion.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/input_utils.o \
	$(COMMON)/io_utils.o \
	$(COMMON)/lapack.o \
	$(COMMON)/misc.o \
	$(COMMON)/scalapack.o \
	$(COMMON)/sort.o \
	$(COMMON)/wfn_io_hdf5.o \
	$(COMMON)/wfn_rho_vxc_io.o \
	bgw_mpi.o \
	distribution.o \
	inread.o \
	iteration_data.o \
	kpoint_pool.o \
	pseudopot.o \
	pseudopot_vkb.o \
	split_spin.o \
	wfn_io.o
