# BerkleyGW: NonLinearOptics Makefile
#
# Most rules are defined in Common/common-rules.mk, and dependencies in
# make_deps.mk. Rerun dev-scripts/make_deps.py if you change any dependency!

PREFIX=..
include $(PREFIX)/Common/common-rules.mk

LIBS = $(SCALAPACKLIB) $(LAPACKLIB) $(FFTWLIB) $(HDF5LIB)

default: all

all: nonlinearoptics$(FLAVOR).x

# FHJ: no difference between flavored and full clean
clean-flavored: clean


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk
