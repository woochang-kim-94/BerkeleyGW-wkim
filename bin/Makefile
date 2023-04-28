PREFIX=..
include $(PREFIX)/Common/common-rules.mk

cleanall: cleanall-bin

# Get rid of the symbolic links and the executables they point to
# For scripts, only get rid of the symlinks
cleanall-bin:
	-for FILE in *.x; do $(REMOVE) $${FILE} `readlink $${FILE}`; done
	-find . \( \( -name \*.py -or -name \*.sh -or -name \*.pl \) -and -type l \) -delete

.PHONY: cleanall-bin

#
uninstall: cleanall-bin
