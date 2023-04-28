#!/bin/bash -l

# Script to patch ESPRESSO and install pw2bgw and bgw2pw in an ESPRESSO build,
# as described in the headers of those files.
# Usage: Go to espresso directory (e.g. cd espresso-4.3.2).
# Then run this script (e.g. sh ../BerkeleyGW/MeanField/ESPRESSO/version-4.3.2/install.sh)
# D. Strubbe Nov 2011

DIR=$(cd $(dirname "$0"); pwd)
QEDIR=$PWD

cd ..
patch -p0 -i $DIR/espresso-4.3.2-patch

cd $QEDIR
ln -s $DIR/pw2bgw.f90 PP/
ln -s $DIR/bgw2pw.f90 PP/
ln -s $DIR/make.bgw PP/

echo "include make.bgw" >> PP/Makefile

make pw
#BEGIN_INTERNAL_ONLY
make ph
#END_INTERNAL_ONLY
cd PP && make pw2bgw.x bgw2pw.x
