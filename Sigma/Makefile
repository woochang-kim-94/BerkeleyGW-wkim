# BerkleyGW: Sigma Makefile
#
# Most rules are defined in Common/common-rules.mk, and dependencies in
# make_deps.mk. Rerun dev-scripts/make_deps.py if you change any dependency!

PREFIX=..
include $(PREFIX)/Common/common-rules.mk

LIBS = $(TILE_LIBS) $(TILE_LIBS_EXT) $(ELPALIB) $(SCALAPACKLIB) $(LAPACKLIB) $(FFTWLIB) $(HDF5LIB) $(PERFORMANCE)

default: all

all: sigma tools

sigma: sigma$(FLAVOR).x

tools: offdiag$(FLAVOR).x sig2wan.x python_scripts 
#BEGIN_INTERNAL_ONLY
tools: offdiag_QSGW.x sympert.x
#END_INTERNAL_ONLY

# FHJ: no difference between flavored and full clean
clean-flavored: clean


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk
