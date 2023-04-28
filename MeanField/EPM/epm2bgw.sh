#!/bin/bash -l

# David Strubbe, Aug 2012
# Dummy wrapper that just runs epm2bgw.x.
# To provide same interface as with epm2bgw.cplx.sh, for use in testsuite.

if [ -n "$SHEXEC" ]; then
  echo bash -c \"$SHEXEC `dirname $0`/epm2bgw.x\" >&2
fi
(cat - ) | bash -c "$SHEXEC `dirname $0`/epm2bgw.x"
