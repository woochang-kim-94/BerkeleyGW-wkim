#!/bin/bash -l

# David Strubbe, October 2010
# Wrapper to run epm2bgw.x with complex wavefunctions, using same input
# file as for real wavefunctions. To reduce duplication in testsuite.

if [ -n "$SHEXEC" ]; then
  echo bash -c \"$SHEXEC `dirname $0`/epm2bgw.x\" >&2
fi
( echo "real_or_complex 2"
cat - ) | bash -c "$SHEXEC `dirname $0`/epm2bgw.x"
