## Makefile (D. Prendergast Jul 2008)
##

PREFIX=..
include $(PREFIX)/Common/common-rules.mk

default: common
all: common
clean: clean-qhull clean-phonon
# FHJ: flavored clean is clean without clean-qhull
clean-flavored: do-clean

# The real content is in common-rules.mk

# FHJ: all rules for files in Common/ in common-rules.mk include a leading
# Common/ folder. This is not convenient if you just want to go to Common/
# and build a specific target. The rule below allows one to directly build a
# file without having to specify the Common/ subdirectory, i.e., you can type
# make kernel_io.o
# instead of
# make ../Common/kernel_io.o
OBJS = $(ALL_COMOBJ) $(GLOBOBJ)
$(OBJS): %: ../Common/%
$(OBJS:.o=_m.mod): %: ../Common/%
