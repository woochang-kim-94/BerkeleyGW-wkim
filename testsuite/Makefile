include ../arch.mk

list:
	@echo
	@echo "testsuite Makefile"
	@echo
	@echo "Available make targets:"
	@echo
	@echo "make check"
	@echo "make check-save"
ifdef TESTSCRIPT
	@echo "make check-jobscript"
	@echo "make check-jobscript-save"
endif

CHECK=uptime; sh make_complex_test.sh Si-EPM/Si.test && sh make_complex_test.sh Si-EPM/Si_hdf5.test && sh make_complex_test.sh Si-EPM_subspace/Si_subspace.test && ./run_testsuite.sh ${BGW_CHECK_FLAGS}
SETVAR=SAVETESTDIRS=yes

ifdef TESTSCRIPT
check-jobscript:
	$(TESTSCRIPT)
# keep working directories
check-jobscript-save:
	$(SETVAR) make check-jobscript
endif

check:
	$(CHECK) -w

# Use version below instead to keep work directories,
# for later examination and debugging.
check-save:
	$(SETVAR) make check

# To be called from testsuite job scripts
# Only call directly if MPI is used without scheduler
check-parallel:
	$(CHECK)
