## Makefile (D. Strubbe May 2010)
##

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

OBJS = $(GLOBALOBJS) $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/check_inversion.o $(COMMON)/sort.o

default: bgw2para perls

PERLS = kptlist.pl qptlist.pl
perls: $(addprefix $(PREFIX)/bin/,$(PERLS))

$(PREFIX)/bin/%.pl: %.pl
	ln -sf $(PWD)/$^ $(PREFIX)/bin

bgw2para: bgw2para$(FLAVOR).x rho2cd.x

bgw2para$(FLAVOR).x: bgw2para.o iolib.o band.o $(COMMON)/input_utils.o $(COMMON)/blas.o $(OBJS)
	$(LINK) $(FOPTS) -o $@ $^ $(LAPACKLIB)
	$(INSTALL_CMD)

rho2cd.x: rho2cd.o $(OBJS)
	$(LINK) $(FOPTS) -o $@ $^
	$(INSTALL_CMD)


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk
