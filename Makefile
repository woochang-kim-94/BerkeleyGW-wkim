## Makefile (D. Prendergast Jul 2008)
##

PREFIX=.
include $(PREFIX)/Common/common-rules.mk

pre:
	$(MAKE) donkey && $(MAKE) common

default: list

list:
	@echo
	@echo "BerkeleyGW Makefile"
	@echo "==================="
	@echo
	@echo "Below are the available make targets for the BerkeleyGW software package."
	@echo "The recommended target is 'make all-flavors'."
	@echo
	@echo "If you are only interested in building a specific flavor of the code and"
	@echo "if you do not wish to execute the testsuite, you can also just issue"
	@echo "'make real' or 'make cplx' to build just the real or complex version of the"
	@echo "code, respectively."
	@echo
	@echo "Likewise, if the make processes fails, you might benefit from first building the"
	@echo "real and then the complex version of the code separately. Just make sure you type"
	@echo "'make clean-flavored' before you make the other flavor."
	@echo
	@echo "RECOMMENDED TARGET:"
	@echo "> make all-flavors           Builds the real and complex flavors of the code. Same as:"
	@echo "                             make clean && make real && make clean-flavored && make cplx"
	@echo
	@echo "Other targets:"
	@echo "  make real                  Builds just the real flavor of the code. Same as:"
	@echo "                             cp flavor_real.mk flavor.mk && make all"
	@echo "  make cplx                  Builds just the complex flavor of the code. Same as:"
	@echo "                             cp flavor_cplx.mk flavor.mk && make all"
	@echo "  make all                   Builds the current flavor (in flavor.mk) of the code"
	@echo "  make install               Install BerkeleyGW into INSTDIR"
	@echo "  make clean                 Clean all object, module, and preprocessed files"
	@echo "  make clean-flavored        Clean targets that depend on flavor.mk"
	@echo "  make cleanall              Clean everything, including binary (*.x) file"
	@echo
	@echo "Running the testsuite:"
	@echo "  make check                 Runs the testsuite of BerkeleyGW"
	@echo "  make check-save            Same as 'check', but always saves the result"
ifdef TESTSCRIPT
	@echo "  make check-jobscript       Runs the testsuite via a submission script"
	@echo "  make check-jobscript-save  Same as 'check-jobscript', but always saves the result"
endif
	@echo
	@echo "Targets for specific subpackages and using the current flavor (in flavor.mk):"
	@echo "  make common                Common components of the code"
	@echo "  make library               BerkeleyGW I/O library"
	@echo "  make epsilon               Epsilon code"
	@echo "  make sigma                 Sigma code"
	@echo "  make bse                   BSE code"
	@echo "  make plotxct               PlotXct utility"
	@echo "  make meanfield             Mean-field codes and utilities"
	@echo "  make nonlinearoptics       Nonlinear optics code"
#BEGIN_INTERNAL_ONLY
	@echo "  make gw+c                  GW + Cumulant"
#END_INTERNAL_ONLY
	@echo "  make visual                Visualization tools"
	@echo

help: list


# FHJ: This is the master code which defines the rules to make and clean the targets.
# The rules we define are make-*, clean-*, clean-flavored-*, and cleanall-*.
DIRS = Common bin library Epsilon Sigma BSE PlotXct MeanField NonLinearOptics Visual
#BEGIN_INTERNAL_ONLY
DIRS += GW+C
#END_INTERNAL_ONLY

# FHJ: Defines targets such as make-Epsilon.
MAKE_DIRS = $(addprefix make-,$(DIRS))
# FHJ: Full clean: clean everything, flavorfull and flavorless.
CLEAN_DIRS = $(addprefix clean-,$(DIRS))
# FHJ: Flavored clean: clean everything that depends on the flavor, but saves flavorless builds.
CLEAN_FLAVORED_DIRS = $(addprefix clean-flavored-,$(DIRS))
CLEANALL_DIRS = $(addprefix cleanall-,$(DIRS))
.PHONY: $(MAKE_DIRS) $(CLEAN_DIRS) $(CLEAN_FLAVORED_DIRS) $(CLEANALL_DIRS)
# FHJ: These fancy rules are called static pattern rules. We used them instead
# of the implicit pattern rules because they play nice with .PHONY.
$(MAKE_DIRS): make-%:
	cd $* && $(MAKE)
$(CLEAN_DIRS): clean-%:
	cd $* && $(MAKE) clean
$(CLEAN_FLAVORED_DIRS): clean-flavored-%:
	cd $* && $(MAKE) clean-flavored
$(CLEANALL_DIRS): cleanall-%:
	cd $* && $(MAKE) cleanall

# FHJ: These are just shortcuts if the user wants to build
# only a specific subdirectory
# Note that the rule common is already define in Common/common-rules.mk
library: common
	$(MAKE) make-library
epsilon: common
	$(MAKE) make-Epsilon
sigma: common
	$(MAKE) make-Sigma
bse: common
	$(MAKE) make-BSE
plotxct: common
	$(MAKE) make-PlotXct
meanfield: common
	$(MAKE) make-MeanField
nonlinearoptics: common
	$(MAKE) make-NonLinearOptics
.PHONY: common library epsilon sigma bse plotxct meanfield nonlinearoptics
.PHONY: all donkey all-flavors real cplx

#
# All
#
all_: $(MAKE_DIRS)
all:
	$(MAKE) pre && $(MAKE) all_

all-flavors:
	$(MAKE) clean
	$(MAKE) real
	$(MAKE) clean-flavored
	$(MAKE) cplx

real:
	rm -f flavor.mk
	cp flavor_real.mk flavor.mk
	$(MAKE) all

cplx:
	rm -f flavor.mk
	cp flavor_cplx.mk flavor.mk
	$(MAKE) all

install: all
ifdef INSTDIR
	mkdir -p $(INSTDIR)/bin
	install bin/*.x $(INSTDIR)/bin/
	install bin/*.sh $(INSTDIR)/bin/
	install bin/*.py $(INSTDIR)/bin/
	install bin/*.pl $(INSTDIR)/bin/

	mkdir -p $(INSTDIR)/lib
	install library/*.a $(INSTDIR)/lib/
	mkdir -p $(INSTDIR)/include
	install library/*.mod $(INSTDIR)/include/

	mkdir -p $(INSTDIR)/share
	mkdir -p $(INSTDIR)/share/BerkeleyGW
# install cannot work on a whole directory
	cp -rf examples  $(INSTDIR)/share/BerkeleyGW/
	cp -rf testsuite $(INSTDIR)/share/BerkeleyGW/
else
	$(error Error: Please define installation prefix INSTDIR via 'make install INSTDIR='.)
endif

check: all
	cd testsuite && $(MAKE) check

check-save: all
	cd testsuite && $(MAKE) check-save

ifdef TESTSCRIPT
check-jobscript: all
	cd testsuite && $(MAKE) check-jobscript
check-jobscript-save: all
	cd testsuite && $(MAKE) check-jobscript-save
endif

#
# Clean
#
clean: $(CLEAN_DIRS)
clean-flavored: $(CLEAN_FLAVORED_DIRS)
cleanall: cleanall-bin $(CLEANALL_DIRS)

