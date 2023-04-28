#!/usr/bin/env python

# Update bsemat.h5 files generated prior to October 2020
# Renames center_of_mass to _exciton_Q_shift
# Diana Y. Qiu (March 2021)

import sys
import h5py
if len(sys.argv)!=2:
    print('Usage: %s bsemat_in.h5'%(sys.argv[0]))
    sys.exit(1)

fname_in = sys.argv[1]
f=h5py.File(fname_in,'r+')

f['bse_header/kpoints/exciton_Q_shift']=f['bse_header/kpoints/center_mass_q']

del f['bse_header/kpoints/center_mass_q']
