# BerkleyGW: Meanfield/ESPRESSO Makefile
#
# Most rules are defined in Common/common-rules.mk, and dependencies in
# make_deps.mk. Rerun dev-scripts/make_deps.py if you change any dependency!

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

LIBS = $(SPGLIB)/libsymspg.a $(LAPACKLIB)

default: all

all: kgrid.x python_scripts

# FHJ: Extra dependencies
kgrid.x: $(SPGLIB)/libsymspg.a


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk
