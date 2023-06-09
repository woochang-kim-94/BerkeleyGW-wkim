# This file was generated by make_deps.py. Do not edit this directly.


#Compile dependencies:
extrapolate_dynmat.o : \
	$(COMMON)/global_m.mod \
	$(COMMON)/lapack_m.mod
kgrid.o : \
	$(COMMON)/global_m.mod \
	$(COMMON)/misc_m.mod \
	$(COMMON)/symmetries_m.mod \
	kgrid_routines_m.mod
kgrid_routines.o kgrid_routines_m.mod: \
	$(COMMON)/global_m.mod \
	$(COMMON)/misc_m.mod

#Link dependencies:
extrapolate_dynmat.x extrapolate_dynmat$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/lapack.o \
	extrapolate_dynmat.o
kgrid.x kgrid$(FLAVOR).x: \
	$(GLOBALOBJS) \
	$(COMMON)/blas.o \
	$(COMMON)/misc.o \
	$(COMMON)/sort.o \
	$(COMMON)/symmetries.o \
	kgrid.o \
	kgrid_routines.o
