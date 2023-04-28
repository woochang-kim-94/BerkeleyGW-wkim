# BerkleyGW: TDGW Makefile

PREFIX=..
include $(PREFIX)/Common/common-rules.mk

# add debug flags only for this folder
FOPTS += -g -debug inline-debug-info -traceback -check all -ftrapuv -init=snan

LIBS = $(TILE_LIBS) $(TILE_LIBS_EXT) $(ELPALIB) $(SCALAPACKLIB) $(LAPACKLIB) $(FFTWLIB) $(HDF5LIB) $(PERFORMANCE)

default: all

all: dynamics

dynamics: dynamics.x

dynamics.x: \
	$(GLOBALOBJS) \
	../Common/write_program_header.o \
	../Common/version.o \
  ../Common/io_utils.o \
	dynamics_main.o \
	input_dynamics.o \
	equation_of_motion.o \
	utility.o \
	propagation.o \
	propagation_q.o \
	measurement.o

dynamics_main.o: \
	$(GLOBALOBJS) \
	input_dynamics.o \
	equation_of_motion.o \
	propagation.o \
	propagation_q.o \
	measurement.o

input_dynamics.o: \
	$(GLOBALOBJS) \
	equation_of_motion.o \
	utility.o

equation_of_motion.o: \
	$(GLOBALOBJS) \
	utility.o

utility.o: \
	$(GLOBALOBJS)

propagation.o: \
	$(GLOBALOBJS) \
	utility.o \
	equation_of_motion.o \
	measurement.o

propagation_q.o: \
	$(GLOBALOBJS) \
	utility.o \
	equation_of_motion.o \
	measurement.o

measurement.o: \
	$(GLOBALOBJS) \
	equation_of_motion.o
