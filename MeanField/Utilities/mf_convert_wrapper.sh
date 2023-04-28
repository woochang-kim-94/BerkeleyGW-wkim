#!/bin/bash -l

# David Strubbe, October 2010
# Wrapper for mf_convert.x
# Automatically detects whether input file is ascii or binary
# and passes on to the executable. Fortran cannot do this.

if [ $# -ne 2 ]; then
    echo "Usage: $0 infile outfile" >&2
    exit 999
fi

file_type=`file -b $1`

if [[ $file_type == *text* ]]; then
    FORM=A2B
else
    FORM=B2A
fi

echo bash -c \"$SHEXEC `dirname $0`/mf_convert.x $FORM $1 $2\" >&2
bash -c "$SHEXEC `dirname $0`/mf_convert.x $FORM $1 $2"
