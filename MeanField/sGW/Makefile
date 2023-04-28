# BerkleyGW: Meanfield/Utilities Makefile
#
# Most rules are defined in Common/common-rules.mk, and dependencies in
# make_deps.mk. Rerun dev-scripts/make_deps.py if you change any dependency!

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

LIBS = $(TILE_LIBS_EXT) $(SCALAPACKLIB) $(LAPACKLIB) $(FFTWLIB) $(HDF5LIB) $(PERFORMANCE)

default: all

all: bgw2sgw.x

# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk

# External dependencies:
bgw2sgw.x: $(COMMON)/periodic_table.o
