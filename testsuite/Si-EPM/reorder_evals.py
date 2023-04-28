#!/usr/bin/env python

# Parses the file eigenvalues.dat, and output the positive eigenvalues sorted
# by energy. We need this script because the order of the eigenvalues changes
# depending on the solver we use. The situation might change in the future.

# Felipe H. da Jornada
# Jan 2015, Berkeley

import sys
import numpy as np

fname_in, fname_out = sys.argv[1:3]
with open(fname_in) as f_in:
    with open(fname_out,'w') as f_out:
        # Copy header
        for i in range(4):
            f_out.write(f_in.readline())

        # Reorder positive evals
        data = np.loadtxt(f_in)
        cond = data[:,0] > 0
        data = data[cond]
        order = np.argsort(data[:,0])
        np.savetxt(f_out, data[order])
