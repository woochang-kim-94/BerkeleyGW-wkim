## Makefile (G. Samsonidze Jan 2010)
##

PREFIX=..
include $(PREFIX)/Common/common-rules.mk

default: surface python_scripts

surface: surface.x

surface.x: surface.o $(COMMON)/wfn_utils.o $(COMMON)/periodic_table.o
	$(C_LINK) $(C_OPTS) -o $@ $^
	$(INSTALL_CMD)

surface.o : $(COMMON)/wfn_utils.h
