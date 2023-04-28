## Makefile (Felipe H. da Jornada -  Dec 2012)
## This is not the Makefile that comes with Qhull.
## It has been simplified and made compatible with the arch.mk of BerkeleyGW.

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

QHULL_HEADERS = geom.h io.h libqhull.h mem.h merge.h poly.h qhull_a.h qset.h \
random.h stat.h user.h
QHULL_SRCS = geom2.c geom.c global.c io.c libqhull.c mem.c merge.c poly2.c \
poly.c qset.c random.c rboxlib.c stat.c user.c usermem.c userprintf.c \
userprintf_rbox.c
QHULL_OBJS = $(QHULL_SRCS:.c=.o)
BGW_HEADERS = libtile_qhull.h libtile_qhull_find.h qhull_blas.h
BGW_SRCS = libtile_qhull.c libtile_qhull_find.c test_delaunay_c.c test_delaunay_f.f90
BGW_OBJS = $(subst .f90,.o,$(BGW_SRCS:.c=.o))
ALL_OBJS = $(QHULL_OBJS) $(BGW_OBJS)

#Don`t use $(TILE_LIBS) here, since we might have to use multiple libraries in
# the future to implement different tessellation routines.
LIB = libtile_qhull.a

default : all
lib : $(LIB)
all :
	$(MAKE) $(LIB) && $(MAKE) test_delaunay_f.x test_delaunay_c.x

.PHONY: lib all

$(ALL_OBJS) : $(QHULL_HEADERS)
$(BGW_OBJS) : $(BGW_HEADERS)

$(LIB) : libtile_qhull.o libtile_qhull_find.o $(QHULL_OBJS)
	$(AR) ru $@ $^

test_delaunay_f.o: $(COMMON)/tile_m.mod

LIBS = $(TILE_LIBS_EXT) $(LAPACKLIB)

test_delaunay_f.x : test_delaunay_f.o $(GLOBALOBJS) $(COMMON)/tile.o $(TILE_LIBS)
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)

test_delaunay_c.x : test_delaunay_c.o $(LIB)
	$(C_LINK) $(C_OPTS) -o $@ $^ $(LAPACKLIB)
