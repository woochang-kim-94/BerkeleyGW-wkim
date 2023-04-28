$(if $(wildcard $(PREFIX)/arch.mk),,$(error Error: Please create arch.mk from config/ directory for machine-dependent configuration))
include $(PREFIX)/arch.mk

ifneq ($(wildcard $(PREFIX)/flavor.mk),)
  include $(PREFIX)/flavor.mk
endif
# it is ok not to be present for some targets such as 'clean' and 'all-flavors'
# if it is needed, it will cause an error according to the flavor.mk target below.

ifeq ($(TYPEFLAG),-DCPLX)
  FLAVOR  = .cplx
else
  FLAVOR  = .real
endif

ifeq ($(findstring -DFFTW3,$(MATHFLAG)),-DFFTW3)
  $(warning -DFFTW3 is deprecated; please use -DUSEFFTW3 in MATHFLAG in arch.mk file)
#  MATHFLAG := $(subst -DFFTW3,-DUSEFFTW3,$(MATHFLAG))
  MATHFLAG += -DUSEFFTW3
endif

ifeq ($(findstring -DUSESCALAPACK,$(MATHFLAG)),-DUSESCALAPACK)
  ifneq ($(findstring -DMPI,$(PARAFLAG)),-DMPI)
    $(error -DUSESCALAPACK flag requires -DMPI; please check your arch.mk file)
  endif
endif

#Default values for the C++ library, stored in $(CC_LINK_LIBS), when linking
#C++ code with Fortran linker. NOTE: we don`t have a default value for xlf!
ifeq ($(CC_LINK_LIBS),)
  ifeq ($(findstring -DINTEL,$(COMPFLAG)),-DINTEL)
    CC_LINK_LIBS = -cxxlib
  else ifeq ($(findstring -DPGI,$(COMPFLAG)),-DPGI)
    CC_LINK_LIBS = -pgcpplibs
  else ifeq ($(findstring -DORACLE,$(COMPFLAG)),-DORACLE)
    CC_LINK_LIBS = -lCstd -lCrun -lc
  else
#Default case: works with GNU,G95,CRAY,OPEN64 plus many combinations of
# compilers, such as NAG+g++, ABSOFT+g++, etc.
    CC_LINK_LIBS = -lstdc++
  endif
endif

ifeq ($(findstring USESCALAPACK,$(MATHFLAG)),)
  SCALAPACKLIB =
endif

# Fortran preprocessing
CPPOPT  = $(TYPEFLAG) $(COMPFLAG) $(PARAFLAG) $(MATHFLAG) $(DEBUGFLAG)

# C/C++
C_CPPOPT  = $(TYPEFLAG) $(COMPFLAG) $(C_PARAFLAG) $(C_DEBUGFLAG)

COMMON = $(PREFIX)/Common
SPGLIB = $(PREFIX)/MeanField/spglib-1.0.9

PHONON_ROOT = $(PREFIX)/Common/Phonon
QHULL_ROOT = $(PREFIX)/Common/qhull
VORO_ROOT = $(PREFIX)/Common/voro++
TILE_LIBS = $(QHULL_ROOT)/libtile_qhull.a
ifeq ($(findstring -DUSEVORO,$(MATHFLAG)),-DUSEVORO)
  TILE_LIBS += $(VORO_ROOT)/libtile_voro.a
  TILE_LIBS_EXT = $(CC_LINK_LIBS)
endif

# this one is for C, C++, and cpp step for Fortran
INCLUDE = -I$(COMMON)
# this one is for Fortran compilation
FTNINC = $(INCFLAG) $(COMMON) $(INCFLAG) $(FFTWINCLUDE)
ifeq ($(findstring -DHDF5,$(MATHFLAG)),-DHDF5)
  FTNINC += $(INCFLAG) $(HDF5INCLUDE)
endif
ifeq ($(findstring -DUSEELPA,$(MATHFLAG)),-DUSEELPA)
  FTNINC += $(INCFLAG) $(ELPAINCLUDE)
endif
#FHJ: We are bundling PRIMME`s include files in Common/ for now.
#ifeq ($(findstring -DUSEPRIMME,$(MATHFLAG)),-DUSEPRIMME)
#  FTNINC += $(INCFLAG) $(PRIMMEINCLUDE)
#endif

MAKE_CLEAN = -$(REMOVE) *.o *.p.f *~ core *__genmod.f90 *.a
# ifort -warn all creates __genmod.f90 files as pseudo-modules
INSTALL_CMD = ln -sf $(PWD)/$@ $(PREFIX)/bin
F90_CMD = $(F90free) $(FTNINC) -c $(FOPTS) $(basename $<).p.f -o $(basename $<).o $(MOD_OPT)$(dir $<)
F90_CMD_NOOPT = $(F90free) $(FTNINC) -c $(FNOOPTS) $(basename $<).p.f -o $(basename $<).o $(MOD_OPT)$(dir $<)
F90_ASM_CMD = $(F90free) $(FTNINC) -S $(FOPTS) $(basename $@).p.f -o $(basename $@).s $(MOD_OPT)$(dir $<)
# $(MOD_OPT) directs where to put the resulting *.mod file

#clang C-preprocessing treats files incorrectly if they have .F90 extension
ifeq ($(findstring clang,$(FCPP)),clang)
  f90_CPP = cp $(basename $<).f90 $(basename $<)_cp.F90; $(FCPP) $(INCLUDE) $(CPPOPT) $(basename $<)_cp.F90 > $(basename $<).p.f; $(REMOVE) $(basename $<)_cp.F90
else
  f90_CPP = $(FCPP) $(INCLUDE) $(CPPOPT) $< > $(basename $<).p.f
endif
F90_CPP = $(FCPP) -P $(INCLUDE) $(CPPOPT) $< > $(basename $<).p.f
ifneq (,$(filter $(COMPFLAG),-DOPEN64 -DABSOFT -DCRAY))
# these compilers name all modules uppercase
MODLINK = @NAME=$(dir $*)`echo $(notdir $*) | tr '[:lower:]' '[:upper:]'`_M.mod; test ! -e $$NAME || ln -sf $(PWD)/$$NAME $(PWD)/$(basename $<)_m.mod
endif

default_goal: default
# GNU make 3.80 and earlier don't have .DEFAULT_GOAL

# all objects to be made from Common directory must appear here
ALL_COMOBJ = bessel.o fullbz.o subgrp.o gmap.o find_kpt_match.o \
      scalapack.o inversion.o mtxel_optical.o \
      write_matrix.o read_matrix.o minibzaverage.o vcoul_generator.o \
      trunc_cell_wire.o trunc_cell_box.o \
      trunc_cell_box_d.o trunc_scell_box_d.o fixwings.o check_screening.o \
      checkgriduniformity.o checkbz.o eqpcor.o createpools.o \
      read_cube.o fft_parallel.o norm.o genwf.o irrbz.o \
      check_inversion.o write_program_header.o version.o wfn_rho_vxc_io.o read_rho_vxc.o \
      kernel_io.o random.o sort.o blas.o scalapack.o lapack.o fftw.o misc.o essl.o input_utils.o \
      splines.o groupk.o symmetries.o scissors.o hdf5_io.o wfn_io_hdf5.o epsread_hdf5.o epswrite_hdf5.o \
      io_utils.o cells.o tile.o inread_common.o print_version_info.o \
      bgw_mpi.o elpa_interface.o evecs.o test_evecs.o scalapack_hl.o references.o \
      algos_common.o accel_fft.o accel_linalg.o accel_memory.o

ALL_COMOBJ += susymmetries.o
ALL_COMMON_OBJ = $(addprefix $(COMMON)/,$(ALL_COMOBJ))

ALL_COMOBJ_C = periodic_table.o wfn_utils.o
ALL_COMMON_OBJ_C = $(addprefix $(COMMON)/,$(ALL_COMOBJ_C))

# These are involved in the modules which all routines must use.
# Note that object files written in C produce no modules.
GLOBOBJ_F := global.o typedefs.o nrtype.o push_pop.o message.o peinfo.o timing.o scalapack_aux.o os.o
GLOBOBJ_C := os_c.o
GLOBOBJ := $(GLOBOBJ_F) $(GLOBOBJ_C)
GLOBMOD := $(GLOBOBJ_F:.o=_m.mod)
GLOBALOBJS := $(addprefix $(COMMON)/,$(GLOBOBJ))
GLOBALMODS := $(addprefix $(COMMON)/,$(GLOBMOD))


$(PREFIX)/flavor.mk :
	@if [ ! -f $(PREFIX)/flavor.mk ]; then echo "Error: Please create flavor.mk from flavor_real.mk or flavor_cplx.mk."; false ; fi


# FHJ: Manual rules that cannot the solved: they involve languages other than Fortran:
$(COMMON)/wfn_utils.o: $(COMMON)/wfn_utils.h $(COMMON)/periodic_table.h $(COMMON)/periodic_table.o
$(COMMON)/periodic_table.o: $(COMMON)/periodic_table.h
$(COMMON)/f_defs.h: $(COMMON)/compiler.h $(PREFIX)/flavor.mk
$(COMMON)/global.o $(COMMON)/global_m.mod: $(COMMON)/f_defs.h
$(COMMON)/os.o $(COMMON)/os_m.mod: $(COMMON)/os_c.o

# FHJ: Automatically populate variables associated with the last git commit and
# corresponding author & date in version_git.h. We also store the version of
# the public release, or the last public release on which this git commit is
# based on, in version_base.h.

# FHJ: Determine if this is a development or public release.
IS_GIT := $(shell git status -u no 1>/dev/null 2>&1 && echo '.true.' || echo '.false.')
ifeq ($(IS_GIT),.true.)
# FHJ: If this is a git version, get the following variables from the last commit.
GIT_AUTHOR := $(shell git log -1 --pretty=format:'%an <%ae>')
GIT_COMMIT := $(shell git --no-pager describe --tags --always --dirty)
GIT_DATE := $(shell git log -1 --pretty=format:'%ai')
else
GIT_AUTHOR :=
GIT_COMMIT :=
GIT_DATE :=
endif
# FHJ: We put these define statements, read by version.f90, into version_git.h
define VERSION_GIT_H
#define IS_GIT $(IS_GIT)
#define GIT_AUTHOR "$(GIT_AUTHOR)"
#define GIT_COMMIT "$(GIT_COMMIT)"
#define GIT_DATE "$(GIT_DATE)"
endef
export VERSION_GIT_H

# FHJ: version_git.h has a dependency against the phony version_force target,
# but it will only be changed if VERSION_GIT_H changed. Note that, if this is
# not a git release, IS_GIT will be .false. and we will only refer to
# VERSION_BASE inside version.f90.
.PHONY: version_force
$(COMMON)/version_git.h: version_force
	@echo "$$VERSION_GIT_H" | cmp -s - $@ || echo "$$VERSION_GIT_H" > $@

# FHJ: verision.f90 depends on the following header files.
# Note that only version_git.h is automatically generated, version_base.h
# must be manually changed every time there is a release!
$(COMMON)/version.o $(COMMON)/version_m.mod: $(COMMON)/version_base.h $(COMMON)/version_git.h


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include $(COMMON)/make_deps.mk


ifeq ($(COMPFLAG),-DOPEN64)
# other compilers provide this module, so only Open64 needs its own here
# To use it, copy in the file http://www.open64.net/doc/db/d6b/omp__lib_8f-source.html to this directory
# and uncomment the dependencies below,  needed to deal with the fact that the module name doesn't end on _m.
# $(COMMON)/peinfo_m.mod : $(COMMON)/OMP_LIB.mod
# $(COMMON)/OMP_LIB.mod : $(COMMON)/omp_lib.o
endif



common: $(ALL_COMMON_OBJ) $(ALL_COMMON_OBJ_C) $(GLOBALOBJS) qhull print_version_info test_evecs
ifeq ($(findstring -DUSEVORO,$(MATHFLAG)),-DUSEVORO)
  common: voro
endif
#BEGIN_INTERNAL_ONLY
# FHJ: Targets that can only be built after the previous conditions were met
common: phonon
#END_INTERNAL_ONLY

print_version_info : $(COMMON)/print_version_info.x

test_evecs : $(COMMON)/test_evecs$(FLAVOR).x

$(COMMON)/test_evecs$(FLAVOR).x : $(COMMON)/test_evecs.o $(GLOBALOBJS) $(COMMON)/evecs.o $(COMMON)/hdf5_io.o $(COMMON)/wfn_io_hdf5.o
	$(LINK) $(FOPTS) -o $@ $^ $(HDF5LIB)
	$(INSTALL_CMD)


$(TILE_LIBS): qhull
ifeq ($(findstring -DUSEVORO,$(MATHFLAG)),-DUSEVORO)
  $(TILE_LIBS): voro
endif

ifeq ($(findstring -DUSEPRIMME,$(MATHFLAG)),-DUSEPRIMME)
common: $(COMMON)/test_primme$(FLAVOR).x

$(COMMON)/test_primme$(FLAVOR).x:
	$(LINK) $(FOPTS) -o $@ $^ $(PRIMMELIB) $(LAPACKLIB)
	$(INSTALL_CMD)
endif


#BEGIN_INTERNAL_ONLY
phonon: $(ALL_COMMON_OBJ) $(ALL_COMMON_OBJ_C) $(GLOBALOBJS)
	cd $(PHONON_ROOT) && $(MAKE)
clean-phonon:
	cd $(PHONON_ROOT) && $(MAKE) clean
cleanall-phonon:
	cd $(PHONON_ROOT) && $(MAKE) cleanall
#END_INTERNAL_ONLY

qhull: $(GLOBALOBJS)
	cd $(QHULL_ROOT) && $(MAKE) lib
clean-qhull:
	cd $(QHULL_ROOT) && $(MAKE) clean
cleanall-qhull:
	cd $(QHULL_ROOT) && $(MAKE) cleanall

voro: $(GLOBALOBJS)
	cd $(VORO_ROOT) && $(MAKE) lib

spglib: $(SPGLIB)/libsymspg.a
clean-spglib:
	cd $(SPGLIB) && $(MAKE) clean
cleanall-spglib:
	cd $(SPGLIB) && $(MAKE) cleanall

SPGLIB_SRC = $(addprefix $(SPGLIB)/, cell.c debug.c hall_symbol.c kpoint.c lattice.c mathfunc.c pointgroup.c primitive.c \
refinement.c site_symmetry.c sitesym_database.c spacegroup.c spg_database.c \
spglib.c symmetry.c spglib_f.c spglib_f_meta.c)
SPGLIB_OBJ = $(SPGLIB_SRC:.c=.o)

SPGLIB_HEADERS = $(addprefix $(SPGLIB)/, cell.h debug.h hall_symbol.h kpoint.h lattice.h mathfunc.h \
pointgroup.h primitive.h refinement.h site_symmetry.h sitesym_database.h \
spacegroup.h spg_database.h spglib.h symmetry.h)

SPGLIB_OBJ : $(SPGLIB_HEADERS)

ifndef AR
AR = /usr/bin/ar
endif
$(SPGLIB)/libsymspg.a : $(SPGLIB_OBJ)
	$(AR) ru $@ $^


# remove all implicit suffix rules
# FHJ: Note: uncommenting the line below breaks the $(AR) variable
#MAKEFLAGS += --no-builtin-variables --no-builtin-rules
.SUFFIXES:

#GNU and gfortran don't update mod files unless interfaces have changed, frustrating make's dependency-checking
ifneq ($(findstring -DG,$(COMPFLAG)),)
RM_MOD_CMD = @$(REMOVE) $(basename $<)_m.mod
endif
ifeq ($(findstring -P,$(FCPP)),)
# if not running with cpp -P, keep .p.f
RM_P_F_CMD = @$(REMOVE) $(basename $<).p.f
endif
# both files are made at once by this rule
%.o %_m.mod : %.f90
	$(RM_MOD_CMD)
	$(f90_CPP)
	$(F90_CMD)
	$(RM_P_F_CMD)
	$(MODLINK)

# Files including _inc.f90 by cpp have their own rules to keep .p.f, and not insert line markers,
# so you can tell which of the various preprocessed versions caused any runtime error.
%.o %_m.mod : %.F90 %_inc.f90
	$(RM_MOD_CMD)
	$(F90_CPP)
	$(F90_CMD)
	$(MODLINK)

# FHJ: Use this rule if you want to generate an assembly file.
# You should only do this if you care about fine optimization.
%.s : %.p.f
	$(F90_ASM_CMD)

# rules to make .p.f by hand for debugging purposes
%.p.f : %.f90
	$(f90_CPP)
%.p.f : %.F90 %_inc.f90
	$(F90_CPP)

%.o : %.cpp
	$(CC_COMP) $(INCLUDE) $(C_CPPOPT) -c $(C_OPTS) $< -o $@

%.o : %.cc
	$(CC_COMP) $(INCLUDE) $(C_CPPOPT) -c $(C_OPTS) $< -o $@

# Alas, it seems that if you compile a .c with gcc or .cpp with g++ you get a leading underscore in the symbol
# table, but it you compile a .c with g++ you do not. Therefore we must make a distinction between C and C++ (for spglib).
%.o : %.c
	$(C_COMP) $(INCLUDE) $(C_CPPOPT) -c $(C_OPTS) $< -o $@


# FHJ: Default rule for executables. Note that we rely on the variables $(LIBS)
# to link against libraries. This variable may need to be ajusted in each
# individual Makefile as necessary.
%.x :
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)
	$(INSTALL_CMD)

# FHJ: automatic generate list of python files to install under bin/
PYTHON_SCRIPTS = $(wildcard *.py)
python_scripts: $(addprefix $(PREFIX)/bin/,$(PYTHON_SCRIPTS))

SHELL_SCRIPTS = $(wildcard *.sh)
shell_scripts: $(addprefix $(PREFIX)/bin/,$(SHELL_SCRIPTS))

# FHJ: install scripts
$(PREFIX)/bin/%.py: %.py
	ln -sf $(PWD)/$^ $(PREFIX)/bin

$(PREFIX)/bin/%.sh: %.sh
	ln -sf $(PWD)/$^ $(PREFIX)/bin

clean: do-clean
do-clean:
	$(MAKE_CLEAN) *.mod

clean-keepmod:
	$(MAKE_CLEAN)

cleanall:
	$(MAKE_CLEAN) *.mod *.x
	-$(REMOVE) -r *.x.dSYM
# note: .x.dSYM is a debug file for MacOS.

donkey:
	@$(PREFIX)/LOGO/print_logo.sh $(FLAVOR) && sleep 1

# all targets which are not the name of a file should be listed here
.PHONY: donkey list default all_ all print_version_info \
check check-save manual voro pre install default_goal \
do-clean clean clean-flavored clean-common \
clean-spglib clean-manual \
cleanall cleanall-bin \
utilities clean-utilities cleanall-utilities cleanall-spglib \
clean-keepmod clean-keepmod-common \
clean-keepmod-utilities clean-cpp clean-wfnutils \
qhull clean-qhull cleanall-qhull clean-keepmod-qhull \
phonon clean-phonon cleanall-phonon clean-keepmod-phonon \
abi2bgw clean-abi2bgw cleanall-abi2bgw clean-keepmod-abi2bgw \
epm siesta sapo icm parabands
