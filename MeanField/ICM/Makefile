## Makefile (G. Samsonidze Jan 2010)
##

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

default: icm

icm: icm.x

icm.x: icm.o $(COMMON)/wfn_utils.o $(COMMON)/periodic_table.o
	$(C_LINK) $(C_OPTS) -o $@ $^
	$(INSTALL_CMD)

icm.o : $(COMMON)/wfn_utils.h
