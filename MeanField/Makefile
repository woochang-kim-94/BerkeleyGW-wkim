## Makefile (D. Strubbe Feb 2010)
##

PREFIX=..
include $(PREFIX)/Common/common-rules.mk

default: all

list:
	@echo
	@echo "MeanField Makefile"
	@echo
	@echo "make default"
	@echo "make all"
	@echo "make parabands"
	@echo "make sapo"
	@echo "make epm"
	@echo "make abi2bgw"
	@echo "make siesta2bgw"
	@echo "make bgw2para"
	@echo "make sgw"
	@echo "make kgrid"
	@echo "make icm"
	@echo "make spglib"
	@echo "make utilities"
	@echo
	@echo "make clean"
	@echo "make clean-keepmod"
	@echo "make cleanall"
	@echo

pre:
	$(MAKE) common && $(MAKE) spglib

# All
#
all_: epm siesta2bgw bgw2para kgrid icm utilities sapo abi2bgw parabands sgw
all:
	$(MAKE) pre && $(MAKE) all_

parabands:
	cd ParaBands && $(MAKE)

sgw:
	cd sGW && $(MAKE)

sapo:
	cd SAPO && $(MAKE)

epm:
	cd EPM && $(MAKE)

abi2bgw:
	cd ABINIT && $(MAKE)

siesta2bgw:
	cd SIESTA && $(MAKE)

bgw2para:
	cd PARATEC && $(MAKE)

kgrid:
	cd ESPRESSO && $(MAKE)

icm:
	cd ICM && $(MAKE)

utilities:
	cd Utilities && $(MAKE)

# Clean
#
# FHJ: Targets that depend on the flavor go below:
clean-flavored: clean-sapo clean-siesta2bgw clean-abi2bgw
# FHJ: Include flavorless target below:
clean: clean-flavored clean-epm clean-bgw2para clean-kgrid clean-icm clean-spglib clean-utilities clean-parabands clean-sgw
#

clean-parabands:
	cd ParaBands && $(MAKE) clean
#
clean-sgw:
	cd sGW && $(MAKE) clean
#
clean-sapo:
	cd SAPO && $(MAKE) clean
#
clean-epm:
	cd EPM && $(MAKE) clean
#
clean-abi2bgw:
	cd ABINIT && $(MAKE) clean
#
clean-siesta2bgw:
	cd SIESTA && $(MAKE) clean
#
clean-bgw2para:
	cd PARATEC && $(MAKE) clean
#
clean-kgrid:
	cd ESPRESSO && $(MAKE) clean
#
clean-icm:
	cd ICM && $(MAKE) clean
#
clean-utilities:
	cd Utilities && $(MAKE) clean

# Cleanall
#
cleanall: cleanall-sapo cleanall-epm cleanall-siesta2bgw cleanall-bgw2para cleanall-parabands \
	cleanall-kgrid cleanall-icm cleanall-utilities cleanall-spglib cleanall-abi2bgw cleanall-sgw
#
cleanall-parabands:
	cd ParaBands && $(MAKE) cleanall
#
cleanall-sgw:
	cd sGW && $(MAKE) cleanall
#
cleanall-sapo:
	cd SAPO && $(MAKE) cleanall
#
cleanall-epm:
	cd EPM && $(MAKE) cleanall
#
cleanall-abi2bgw:
	cd ABINIT && $(MAKE) cleanall
#
cleanall-siesta2bgw:
	cd SIESTA && $(MAKE) cleanall
#
cleanall-bgw2para:
	cd PARATEC && $(MAKE) cleanall
#
cleanall-kgrid:
	cd ESPRESSO && $(MAKE) cleanall
#
cleanall-icm:
	cd ICM && $(MAKE) cleanall
#
cleanall-utilities:
	cd Utilities && $(MAKE) cleanall

# Clean-keepmod
#
clean-keepmod: clean-keepmod-sapo clean-keepmod-epm clean-keepmod-siesta2bgw clean-keepmod-parabands \
	clean-keepmod-bgw2para clean-keepmod-kgrid clean-keepmod-utilities clean-keepmod-abi2bgw clean-keepmod-sgw
#
clean-keepmod-parabands:
	cd ParaBands && $(MAKE) clean-keepmod
#
clean-keepmod-sgw:
	cd sGW && $(MAKE) clean-keepmod
#
clean-keepmod-sapo:
	cd SAPO && $(MAKE) clean-keepmod
#
clean-keepmod-epm:
	cd EPM && $(MAKE) clean-keepmod
#
clean-keepmod-abi2bgw:
	cd ABINIT && $(MAKE) clean-keepmod
#
clean-keepmod-siesta2bgw:
	cd SIESTA && $(MAKE) clean-keepmod
#
clean-keepmod-bgw2para:
	cd PARATEC && $(MAKE) clean-keepmod
#
clean-keepmod-kgrid:
	cd ESPRESSO && $(MAKE) clean-keepmod
#
clean-keepmod-utilities:
	cd Utilities && $(MAKE) clean-keepmod
