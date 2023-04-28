#!/usr/bin/env python

import h5py
import numpy as np
from progress.bar import Bar


def merge_spins(fname_spin1, fname_spin2, fname_out):
    f1 = h5py.File(fname_spin1, 'r')
    f2 = h5py.File(fname_spin2, 'r')
    fout = h5py.File(fname_out, 'w')
    assert f1['mf_header/kpoints/nspin'][()] == 1
    assert f2['mf_header/kpoints/nspin'][()] == 1

    print 'Copying mean-field header'
    fout.copy(f1['mf_header'], 'mf_header')
    fout['mf_header/kpoints/nspin'][()] = 2

    def merge_dsets(dset, axis=0):
        if dset in fout:
            del fout[dset]
        if axis!=0:
            print 'Reading dataset {} from file 1'.format(dset)
        mat1 = f1[dset][()]
        if axis!=0:
            print 'Reading dataset {} from file 2'.format(dset)
        mat2 = f2[dset][()]
        mat12 = np.concatenate((mat1,mat2), axis=axis)
        assert mat12.dtype==mat1.dtype
        assert mat12.dtype==mat2.dtype
        if axis!=0:
            print 'Writing dataset {}'.format(dset)
        fout[dset] = mat12

    print 'Copying G-vectors'
    fout.create_group('/wfns')
    fout.copy(f1['/wfns/gvecs'], '/wfns/gvecs')

    print 'Merging spin-related datasets in mean-field header'
    for dset in ('ifmin', 'ifmax', 'el', 'occ'):
        merge_dsets('mf_header/kpoints/{}'.format(dset))

    print 'Merging WFN coefficients'
    assert np.allclose(f1['/wfns/coeffs'].shape, f2['/wfns/coeffs'].shape)
    shape = f1['/wfns/coeffs'].shape
    nb = shape[0]
    shape = list(shape)
    shape[1] = 2
    fout.create_dataset('/wfns/coeffs', shape=shape, dtype=f1['/wfns/coeffs'].dtype)
    bar = Bar('Copying WFN coefficients', max=nb, bar_prefix=' [', bar_suffix='] ',
        fill='#', suffix='%(percent)d%% - Remaining: %(eta_td)s')
    for ib in range(nb):
        bar.next()
        fout['/wfns/coeffs'][ib,0] = f1['/wfns/coeffs'][ib,0]
        fout['/wfns/coeffs'][ib,1] = f2['/wfns/coeffs'][ib,0]
    bar.finish()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fname_spin1', help=
        'Input WFN.h5 file for the 1st spin channel')
    parser.add_argument('fname_spin2', help=
        'Input WFN.h5 file for the 2nd spin channel')
    parser.add_argument('fname_out', help=
        'Output WFN.h5 file with the merged spin channels')
    args = parser.parse_args()

    merge_spins(args.fname_spin1, args.fname_spin2, args.fname_out)
