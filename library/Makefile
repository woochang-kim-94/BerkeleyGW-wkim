## Makefile (D. Strubbe Dec 2011)
##

PREFIX=..
include $(PREFIX)/Common/common-rules.mk

default: library
all: library
# FHJ: no difference between flavored and full clean
clean-flavored: clean

# NOTE: this target is 'flavorless'; the result is the same whether -DCPLX is used or not.
LIBRARY_MODS = $(GLOBMOD) wfn_rho_vxc_io_m.mod check_inversion_m.mod sort_m.mod
library: libBGW_wfn.a $(LIBRARY_MODS)

# Other codes are not expected to use any of the modules below besides
# wfn_rho_vxc_io_m.mod. However, some compilers require them to be present in
# order to use wfn_rho_vxc_io_m.mod since it depends on them.
$(LIBRARY_MODS): %: $(COMMON)/%
	ln -sf $^ .

libBGW_wfn.a : $(GLOBALOBJS) $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/check_inversion.o $(COMMON)/sort.o
	$(AR) ru $@ $^
