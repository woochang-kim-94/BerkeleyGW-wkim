## Makefile
##

PREFIX=..
include $(PREFIX)/Common/common-rules.mk

##
## Common files
##

# object files from $(COMMON) directory
COMOBJ = version.o
COMMONOBJS = $(addprefix $(COMMON)/,$(COMOBJ)) $(GLOBALOBJS)

SRC = CumExp.f90

OBJ = $(SRC:.f90=.o)


default: cumexp
# FHJ: GW+C is flavorless => clean-flavored doesn`t do anthing
clean-flavored:

# at least currently GW+C is flavorless
cumexp: cumexp.x

cumexp.x: $(OBJ) $(COMMONOBJS)
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)	
	$(INSTALL_CMD)


# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk
