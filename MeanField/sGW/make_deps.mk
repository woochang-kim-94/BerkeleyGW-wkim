# This file was generated by make_deps.py. Do not edit this directly.


#Compile dependencies:
bgw2sgw.o : \
	$(COMMON)/fftw_m.mod \
	$(COMMON)/global_m.mod \
	$(COMMON)/io_utils_m.mod \
	$(COMMON)/wfn_io_hdf5_m.mod \
	$(COMMON)/wfn_rho_vxc_io_m.mod

#Link dependencies:
bgw2sgw.x bgw2sgw$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/check_inversion.o \
	$(COMMON)/fftw.o \
	$(COMMON)/hdf5_io.o \
	$(COMMON)/hdf5_io_data.o \
	$(COMMON)/hdf5_io_safe.o \
	$(COMMON)/io_utils.o \
	$(COMMON)/message.o \
	$(COMMON)/nrtype.o \
	$(COMMON)/os.o \
	$(COMMON)/peinfo.o \
	$(COMMON)/push_pop.o \
	$(COMMON)/sort.o \
	$(COMMON)/timing.o \
	$(COMMON)/wfn_io_hdf5.o \
	$(COMMON)/wfn_rho_vxc_io.o \
	bgw2sgw.o
