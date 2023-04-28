# BerkleyGW: Meanfied/ABINIT Makefile
#
# Most rules are defined in Common/common-rules.mk, and dependencies in
# make_deps.mk. Rerun dev-scripts/make_deps.py if you change any dependency!

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

LIBS = $(LAPACKLIB) $(FFTWLIB)

EXECS = abi2bgw.x check_eig.x getbox.x eqp_convert.x

default: all

all: $(EXECS) python_scripts

# FHJ: no difference between flavored and full clean
clean-flavored: clean

# FHJ: Extra dependencies
$(EXECS): $(COMMON)/periodic_table.o


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk
