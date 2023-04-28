## Makefile (Felipe H. da Jornada -  Aug 2014)
## This is not the Makefile that comes with voro++.
## It has been simplified and made compatible with the arch.mk of BerkeleyGW.

PREFIX:=../..
include $(PREFIX)/Common/common-rules.mk

BGW_HEADERS := libtile_voro.h
BGW_SRCS := libtile_voro.cc test_voro.f90
BGW_OBJS := $(subst .f90,.o,$(BGW_SRCS:.cc=.o))

# List of the common source files
objs_ := cell.o common.o container.o unitcell.o v_compute.o c_loops.o \
     v_base.o wall.o pre_container.o container_prd.o
objs := $(addprefix src/,$(objs_))
src := $(patsubst %.o,%.cc,$(objs))

include Makefile.dep

#Don`t use $(TILE_LIBS) here, since we might have to use multiple libraries in
# the future to implement different tessellation routines.
LIB := libtile_voro.a

default : all
lib : $(LIB)
all :
	$(MAKE) $(LIB) && $(MAKE) test_voro.x
$(LIB) : $(objs)

.PHONY: lib all clean-src clean-local

clean: clean-src clean-local

clean-src:
	cd src && $(MAKE) clean

clean-local:
	$(MAKE_CLEAN) *.mod

$(BGW_OBJS) : $(BGW_HEADERS)

$(LIB) : $(objs) libtile_voro.o 
	$(AR) ru $@ $^

test_voro.o: $(COMMON)/tile_m.mod

test_voro.x : test_voro.o $(GLOBALOBJS) $(COMMON)/tile.o $(TILE_LIBS)
	$(LINK) $(FOPTS) -o $@ $^ $(TILE_LIBS_EXT) $(LAPACKLIB)
