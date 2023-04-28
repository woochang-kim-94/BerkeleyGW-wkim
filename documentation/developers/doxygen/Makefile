.PHONY: default clean

default: doxygen

doxygen: Doxyfile
	doxygen Doxyfile

Doxyfile: Doxyfile.template
	./gen_doxyfile.sh

clean:
	-rm -rf doxygen Doxyfile
