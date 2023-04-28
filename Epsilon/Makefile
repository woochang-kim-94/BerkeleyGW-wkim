# BerkleyGW: Epsilon Makefile
#
# Most rules are defined in Common/common-rules.mk, and dependencies in
# make_deps.mk. Rerun dev-scripts/make_deps.py if you change any dependency!

PREFIX=..
include $(PREFIX)/Common/common-rules.mk

LIBS = $(TILE_LIBS) $(TILE_LIBS_EXT) $(ELPALIB) $(SCALAPACKLIB) $(LAPACKLIB) $(FFTWLIB) $(HDF5LIB) $(PERFORMANCE)

default: all

all: epsilon tools

epsilon: epsilon$(FLAVOR).x

TOOLS = epsmat_merge epsascbin epsbinasc ploteps epsinvomega epsomega printchi eps0sym
ifeq ($(findstring -DHDF5,$(MATHFLAG)),-DHDF5)
TOOLS += epsmat_old2hdf5 epsmat_hdf5_upgrade
endif
tools: python_scripts $(addsuffix $(FLAVOR).x, $(TOOLS)) setup_subsampling_nns.x
#BEGIN_INTERNAL_ONLY
tools: epsmat_intp.x
#END_INTERNAL_ONLY

# FHJ: no difference between flavored and full clean
clean-flavored: clean


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk
