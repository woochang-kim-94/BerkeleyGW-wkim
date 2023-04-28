## Makefile (G. Samsonidze Aug 2008)
##

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

SRC_EPM = atomBasisMod.f90 atomPotentialMod.f90 cgParamsMod.f90 conjGradMod.f90 densityArrayMod.f90 eigenStatesMod.f90 graphMod.f90 gridMod.f90 gVectorsMod.f90 hamiltonianMod.f90 initializeMod.f90 kPointsMod.f90 latticeMod.f90 minimizeMod.f90 normalizeMod.f90 pwHamMod.f90 structureMod.f90 sysParams.f90 tagHandlerMod.f90 tiseOpsMod.f90 typeMod.f90
OBJ_EPM = $(SRC_EPM:.f90=.o)
OBJ_FFT = $(COMMON)/fftw.o $(GLOBALOBJS)

# object files from Common directory
LIBOBJS = $(COMMON)/blas.o $(COMMON)/lapack.o
COMMONOBJS = $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/check_inversion.o $(COMMON)/misc.o $(COMMON)/sort.o $(COMMON)/symmetries.o $(LIBOBJS)

LIBS = $(LAPACKLIB) $(FFTWLIB)

default: all

all: epm epm2bgw

list:
	@echo
	@echo "EPM makefile"
	@echo
	@echo "make epm      - compile the EPM code"
	@echo "make epm2bgw  - compile the BGW wrapper"
	@echo "make all      - compile all the codes"
	@echo "make clean    - remove all the object files"
	@echo "make cleanall - remove all the object files and the executables"
	@echo

epm: epm.x

epm.x: $(OBJ_EPM) $(OBJ_FFT) $(LIBOBJS) epm.o
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)
	$(INSTALL_CMD)

epm2bgw: epm2bgw.x

epm2bgw.x: epm2bgw.o $(OBJ_EPM) $(OBJ_FFT) $(COMMONOBJS) $(SPGLIB)/libsymspg.a
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)
	$(INSTALL_CMD)

# dependencies

epm2bgw.o : $(COMMON)/message_m.mod $(COMMON)/wfn_rho_vxc_io_m.mod $(COMMON)/fftw_m.mod $(COMMON)/check_inversion_m.mod $(COMMON)/symmetries_m.mod $(COMMON)/nrtype_m.mod
tagHandlerMod.o : $(COMMON)/message_m.mod
epm.o epm2bgw.o : $(OBJ_EPM)
atomBasisMod.o : latticeMod.o sysParams.o tagHandlerMod.o $(COMMON)/message_m.mod
atomPotentialMod.o : sysParams.o
cgParamsMod.o : sysParams.o tagHandlerMod.o
conjGradMod.o : eigenStatesMod.o hamiltonianMod.o initializeMod.o minimizeMod.o normalizeMod.o sysParams.o tiseOpsMod.o typeMod.o
densityArrayMod.o : sysParams.o
eigenStatesMod.o : sysParams.o
graphMod.o : sysParams.o eigenStatesMod.o $(COMMON)/message_m.mod
gridMod.o : sysParams.o
gVectorsMod.o : structureMod.o sysParams.o tagHandlerMod.o eigenStatesMod.o densityArrayMod.o
hamiltonianMod.o : eigenStatesMod.o sysParams.o typeMod.o $(COMMON)/blas_m.mod $(COMMON)/lapack_m.mod $(COMMON)/message_m.mod
initializeMod.o : eigenStatesMod.o sysParams.o typeMod.o
kPointsMod.o : graphMod.o latticeMod.o sysParams.o tagHandlerMod.o
latticeMod.o : sysParams.o tagHandlerMod.o
minimizeMod.o : eigenStatesMod.o normalizeMod.o sysParams.o typeMod.o
normalizeMod.o : eigenStatesMod.o sysParams.o typeMod.o
pwHamMod.o : $(COMMON)/fftw_m.mod atomPotentialMod.o cgParamsMod.o gVectorsMod.o eigenStatesMod.o sysParams.o typeMod.o structureMod.o kPointsMod.o densityArrayMod.o tagHandlerMod.o
structureMod.o : atomBasisMod.o latticeMod.o sysParams.o tagHandlerMod.o
tiseOpsMod.o : sysParams.o eigenStatesMod.o typeMod.o $(COMMON)/blas_m.mod
typeMod.o : sysParams.o
$(OBJ_EPM) : $(COMMON)/push_pop_m.mod
