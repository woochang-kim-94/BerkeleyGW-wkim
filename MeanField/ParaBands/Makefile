# BerkleyGW: Meanfield/ParaBands Makefile
#
# Most rules are defined in Common/common-rules.mk, and dependencies in
# make_deps.mk. Rerun dev-scripts/make_deps.py if you change any dependency!

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

LIBS =
ifeq ($(findstring -DUSEPRIMME,$(MATHFLAG)),-DUSEPRIMME)
LIBS += $(PRIMMELIB)
endif
ifeq ($(findstring -DUSEELPA,$(MATHFLAG)),-DUSEELPA)
LIBS += $(ELPALIB)
endif
LIBS += $(SCALAPACKLIB) $(LAPACKLIB) $(HDF5LIB)

default: all

all: parabands$(FLAVOR).x split_spin$(FLAVOR).x python_scripts tests
ifeq ($(findstring -DHDF5,$(MATHFLAG)),-DHDF5)
all: compare_wfns$(FLAVOR).x
endif

TEST_SRCS = $(wildcard test/*/*.mako)
TEST_OBJS = $(patsubst %.mako,%,$(TEST_SRCS))
%.test: %.test.mako
	./test/gen_test.py $< $@
tests: $(TEST_OBJS)


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk

diag_primme.o diag_primme_m.mod: primme_aux.o
parabands$(FLAVOR).x: primme_aux.o
