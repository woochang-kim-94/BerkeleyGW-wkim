## Makefile
##

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

##
## Common files
##

# object files from $(COMMON) directory
COMOBJ = blas.o hdf5_io.o write_program_header.o version.o io_utils.o evecs.o wfn_io_hdf5.o
COMMONOBJS = $(addprefix $(COMMON)/,$(COMOBJ)) $(GLOBALOBJS) $(TILE_LIBS)

SRC = gkq.f90 compute_gss.f90 xct_ph_maps.f90

OBJ = $(SRC:.f90=.o)

LIBS = $(HDF5LIB) $(LAPACKLIB) $(SCALAPACKLIB)


default: test_gkq print_gkq exciton_phonon python_scripts
test_gkq: test_gkq$(FLAVOR).x
print_gkq: print_gkq$(FLAVOR).x
exciton_phonon: exciton_phonon$(FLAVOR).x
clean-flavored: clean

$(PREFIX)/bin/%$(FLAVOR).x:
	ln -sf $(PWD)/$^ $(PREFIX)/bin

test_gkq$(FLAVOR).x: test_gkq.o $(OBJ) $(COMMONOBJS)
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)
	$(INSTALL_CMD)

print_gkq$(FLAVOR).x: print_gkq.o $(OBJ) $(COMMONOBJS)
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)	
	$(INSTALL_CMD)

exciton_phonon$(FLAVOR).x: exciton_phonon.o $(OBJ) $(COMMONOBJS)
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)
	$(INSTALL_CMD)


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk
