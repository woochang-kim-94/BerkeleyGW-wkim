#!/usr/bin/env python

# Parses the output of wfn_dotproduct.x, calculates the sum of the overlaps
# within a degenerate subset, and summarizes the worst error in calculated
# overlaps.

# Felipe H. da Jornada
# May 2014, Berkeley

from __future__ import print_function
import itertools as it
import numpy as np


TOL_DEG = 1e-5

# This is the same as np.unique or np.unique1d with return_index=True.
# I copied it here to make it compatible with older versions of numpy.
def unique_indexed(ar1):
    ar = np.asanyarray(ar1).flatten()
    perm = ar.argsort()
    aux = ar[perm]
    flag = np.concatenate(([True], aux[1:] != aux[:-1]))
    return aux[flag], perm[flag]


# For each block corresponding to a given k-point, this function sums the
# overlap of a band with its degenerate subspace. It also prints the worst
# overlap, and makes sure it doesn't include bands that are degenerate with
# the highest band in the WFN file.
def process_kpt(ikpt, iter_lines):
    kpt = next(iter_lines).split(':')[1][:-1]
    print('*** k-point %d: %s ***'%(ikpt, kpt))

    line = next(iter_lines)
    while not line.startswith('    band'):
        line = next(iter_lines)
    lines = ''.join(iter_lines)
    data = np.fromstring(lines, dtype=np.float64, sep=' ').reshape(-1,7)

    bands1,idx_unique1 = unique_indexed(data[:,0])
    bands1 = bands1.astype(np.int)
    bands2,idx_unique2 = unique_indexed(data[:,2])
    bands2 = bands2.astype(np.int)
    nb = min(len(bands1), len(bands2))
    bands1 = bands1[:nb]
    bands2 = bands2[:nb]
    # FIXME: temporarily disabling this check
    #assert np.all(bands1==bands2), "Number of bands doesn't match between two WFN files."

    en1 = data[idx_unique1,1][:nb]
    en2 = data[idx_unique2,3][:nb]
    assert np.allclose(en1, en2), "Energy mismatch between two WFN files."

    last_degens = []
    overlaps = []

    print('%6s %12s %12s %20s'%('band', 'first deg.', 'last deg.', 'subspace overlap'))
    # For each band ib1, search for the degenerate subspace and sum the overlap
    for ib1 in bands1:
        cond_degen = (data[:,0]==ib1) & (np.fabs(data[:,3]-en1[ib1-1]) < TOL_DEG)
        deg_bands = data[cond_degen,2]
        deg_first = int(deg_bands[0])
        deg_last = int(deg_bands[-1])
        overlap = data[cond_degen,6].sum()
        overlaps.append(overlap)
        last_degens.append(deg_last)
        print('%6d %12d %12d %20.16f'%(ib1, deg_first, deg_last, overlap))

    # Remove bands that might be degenerate to bands not included in the WFN
    last_degens = np.array(last_degens)
    cond_valid = np.array(last_degens)<bands1[-1]
    overlaps = np.array(overlaps)[cond_valid]
    errors = np.fabs(overlaps-1)
    worst_error = np.max(errors)
    print('')
    print('> Worst overlap error for k-point %d: %.3e'%(ikpt,worst_error))
    print('')
    print('')


# Loads `fname` and separates the file into blocks corresponding to different 
# k-points.
def main(fname):
    def new_kpt(line):
        return line.startswith('k-point')

    with open(fname,'r') as f:
        ikpt = 0
        for is_delimiter, iter_lines in it.groupby(f, new_kpt):
            if is_delimiter:
                ikpt += 1
            elif ikpt>0:
                process_kpt(ikpt, iter_lines)


if __name__=='__main__':
    import sys

    if len(sys.argv)<1:
        sys.exit('usage: %s wfn_dotproduct.log'%(sys.argv[0]))

    main(sys.argv[1])
