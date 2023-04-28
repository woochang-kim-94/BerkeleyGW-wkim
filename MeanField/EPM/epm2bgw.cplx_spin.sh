#!/bin/bash -l

# David Strubbe, April 2011
# Wrapper to run epm2bgw.x with complex wavefunctions and spin-polarized,
# using same input file as for real unpolarized wavefunctions.
# To reduce duplication in testsuite.

if [ -n "$SHEXEC" ]; then
  echo bash -c \"$SHEXEC `dirname $0`/epm2bgw.x\" >&2
fi
( echo "real_or_complex 2"
echo "nspin 2"
cat - ) | bash -c "$SHEXEC `dirname $0`/epm2bgw.x"
