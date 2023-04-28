# BerkleyGW: BSE Makefile
#
# Most rules are defined in Common/common-rules.mk, and dependencies in
# make_deps.mk. Rerun dev-scripts/make_deps.py if you change any dependency!

PREFIX=..
include $(PREFIX)/Common/common-rules.mk

ifeq ($(findstring -DUSESCALAPACK,$(MATHFLAG)),-DUSESCALAPACK)
# SSEIG library
SSEIG_ROOT = $(PREFIX)/BSE/SSEIG
SSEIG_LIB = $(SSEIG_ROOT)/libsseig.a
$(SSEIG_LIB):
	cd $(SSEIG_ROOT) && $(MAKE) lib
clean-bse-sseig:
	cd $(SSEIG_ROOT) && $(MAKE) clean

# Full solver library
FULL_SOLVER_ROOT = $(PREFIX)/BSE/full_solver
FULL_SOLVER_LIB = $(FULL_SOLVER_ROOT)/libfull_solver.a
$(FULL_SOLVER_LIB):
	cd $(FULL_SOLVER_ROOT) && $(MAKE) lib
clean-bse-full-solver:
	cd full_solver && $(MAKE) clean

clean: clean-bse-sseig clean-bse-full-solver
.PHONY: clean-bse-sseig clean-bse-full-solver
endif
# FHJ: flavored clean is clean without clean-bse-sseig and clean-bse-full-solver
clean-flavored: do-clean

# put truly external libs here, i.e. not part of our buildsystem at all
LIBS = $(PRIMMELIB) $(TILE_LIBS) $(TILE_LIBS_EXT) $(ELPALIB) $(SCALAPACKLIB) $(LAPACKLIB) $(FFTWLIB) $(HDF5LIB) $(PERFORMANCE)

default: all

all: bse tools

bse: kernel$(FLAVOR).x absorption$(FLAVOR).x inteqp$(FLAVOR).x

tools: python_scripts bseascbin$(FLAVOR).x bsebinasc$(FLAVOR).x bsemat_check$(FLAVOR).x
tools: summarize_eigenvectors$(FLAVOR).x setup_subsampling_csi.x
#BEGIN_INTERNAL_ONLY
tools: forces$(FLAVOR).x vmtxel$(FLAVOR).x
#END_INTERNAL_ONLY


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk

# Extra dependencies for absorption code.
# Must appear after including make_deps.mk!
absorption$(FLAVOR).x: $(FULL_SOLVER_LIB) $(SSEIG_LIB)
