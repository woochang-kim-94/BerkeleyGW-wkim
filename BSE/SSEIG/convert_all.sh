#!/bin/bash

# Script to import the SSEIG library from Meiyue Shao to the BerkeleyGW repository.
# This script essentially converts from Fortran77 to Fortran90 and fixed single
# quotation marks.
# 
# Felipe H. da Jornada (2015)

#The following line is just to test the find statement
##find . -name \*.f -exec sh -c 'echo {} `dirname {}`/`basename {} .f`.f90' \;

nf=$(find . -name \*.f | wc -l)
echo "Converting $nf files"

#Convert from fixed form to free form. Assumes that you installed fixed2free2 with:
#pip install --user fixed2free2
find . -name \*.f -exec sh -c '~/.local/bin/fixed2free2.py {} > `dirname {}`/`basename {} .f`.f90' \;

#Convert single quotation marks to backticks
find . -name \*.f90 -exec sed -i.bak $'s/^\\([^\']*\\)\'\\([^\']*\\)$/\\1`\\2/' '{}' \;

echo "All done! Don't forget to: rm *.f"

#If you want to compare the new .f90 files against the old ones in directory
#OLD_DIR, use the following code:
#OLD_DIR=../SSEIG
#for f in *.f90; do cmp $f $OLD_DIR/$f || vimdiff $f $OLD_DIR/$f; done
