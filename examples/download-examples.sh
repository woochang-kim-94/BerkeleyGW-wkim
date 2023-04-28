#!/bin/bash -l

mkdir -p downloaded
curl -L https://api.github.com/repos/BerkeleyGW/BerkeleyGW-examples/tarball | tar -xz --strip-components 1 -f - -C downloaded
