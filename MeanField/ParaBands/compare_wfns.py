#!/usr/bin/env python

import numpy as np
from scipy.spatial import cKDTree
import h5py
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('fname_ref', help='Reference WFN.h5 file.')
parser.add_argument('fname_new', help='New WFN.h5 file to check against reference.')
parser.add_argument('--terse', default=False, action='store_true',
                    help='Return raw numbers without any additional explanation.')
parser.add_argument('--nb', type=int, help='Number of bands to compare. Defaults to all.')
parser.add_argument('--tol_deg', type=float, default=1e-10,
                    help='Tolerance for degenerate states, in Ry.')
parser.add_argument('--no_deg', default=False, action='store_true',
                    help=('Do not try to find a non-degenerate subspace. Only '
                    'use this if the WFN does not slice a degenerate subspace'))
args = parser.parse_args()

f1 = h5py.File(args.fname_ref, 'r')
f2 = h5py.File(args.fname_new, 'r')

def get_wfns(f, nb=None):
    if nb is None:
        wfns = f['wfns/coeffs'][()]
    else:
        wfns = f['wfns/coeffs'][:nb]
    if wfns.shape[-1]==2:
        wfns = wfns.view(np.complex128)
    return wfns[...,0]


def get_and_check(dset, err_msg):
    data = f1[dset][()]
    assert np.all(data==f2[dset][()]), err_msg
    return data




nk = get_and_check('mf_header/kpoints/nrk', 'number of k-point mismatch')
ns = get_and_check('mf_header/kpoints/nspin', 'number of spin components mismatch')
ngk = get_and_check('mf_header/kpoints/ngk', 'number of G-vectors per k-point mismatch')
nb1 = f1['mf_header/kpoints/mnband'][()]
nb2 = f2['mf_header/kpoints/mnband'][()]
nb = min(nb1, nb2)
if args.nb is not None:
    assert args.nb<=nb
    nb = args.nb

# 1) Check that the eigenvalues match
en1 = f1['mf_header/kpoints/el'][:,:,:nb]
en2 = f2['mf_header/kpoints/el'][:,:,:nb]
assert en1.shape==en2.shape, 'eigenvalue shape mismatch'
err_evals = np.amax(np.abs(en1 - en2))

wfns1 = get_wfns(f1, nb=nb)
wfns2 = get_wfns(f2, nb=nb)
assert wfns1.shape==wfns2.shape, 'wavefunction shape mismatch'

# We make a few consistency checks, where C = C[ib,ig], i.e. Python convention:
# 2-3) Check that each WFN file is orthonormal to within itself, C C* ~ 1,
# 4) Check that both WFN file spawn the same identity, C2 C2* - C1 C1*.
#    This is useful for pseudobands.
# 5) Check that the reassembled Hamiltonian is the same, where the
#    reassembled Hamiltonian is constructed as H = C* E C, on a non-degenerate
#    subspace.

err_ortho1 = 0.
err_ortho2 = 0.
err_cross = 0.
err_ham = 0.

gstart = 0
for ik in range(nk):
    gend = gstart + ngk[ik]
    # For each k-point, find mapping array between G-vectors in original file,
    # f1, and the new file, f2.
    gvecs1 = f1['wfns/gvecs'][gstart:gend]
    gvecs2 = f2['wfns/gvecs'][gstart:gend]
    tree1 = cKDTree(gvecs1)
    dist, ind = tree1.query(gvecs2)
    assert np.allclose(dist,0), (
        'G-vectors mismatch for k-point #{}'.format(ik+1))
    assert np.all(gvecs1[ind]==gvecs2), (
        'Error mapping G-vectors for k-point #{}'.format(ik+1))

    for ispin in range(ns):
        wfns1_ks = wfns1[:nb,ispin,gstart:gend][:,ind]
        wfns2_ks = wfns2[:nb,ispin,gstart:gend]

        # Check 2-3: WFNs orthonormal to within thenselves
        err1 = wfns1_ks.dot(wfns1_ks.conj().T) - np.eye(nb)
        err_ortho1 = max(err_ortho1, np.amax(np.abs(err1)))
        err2 = wfns2_ks.dot(wfns2_ks.conj().T) - np.eye(nb)
        err_ortho2 = max(err_ortho2, np.amax(np.abs(err2)))

        # Check 4: cross WFN orthonormality
        err12 = wfns2_ks.dot(wfns2_ks.conj().T) - wfns1_ks.dot(wfns1_ks.conj().T)

        # Check 5: reassembled Hamiltonian
        # First, find the last non-degenerate band. This test is very picky
        # with degenerancies!
        e1 = en1[ispin,ik]
        e2 = en2[ispin,ik]
        if args.no_deg:
            nb_ = nb
        else:
            for ib in range(nb-2,-1,-1):
                if ((abs(e1[ib]-e1[ib+1]) > args.tol_deg) and 
                    (abs(e2[ib]-e2[ib+1]) > args.tol_deg)):
                    break
            else:
                raise Exception('All bands are degenerate')
            nb_ = ib+1
        ham1 = wfns1_ks[:nb_].conj().T.dot(e1[:nb_,None]*wfns1_ks[:nb_])
        ham2 = wfns2_ks[:nb_].conj().T.dot(e2[:nb_,None]*wfns2_ks[:nb_])
        ham_diff = ham1 - ham2
        norm = np.linalg.norm
        err_ham = max(err_ham, norm(ham_diff)/(norm(ham1)*norm(ham2)))

    gstart = gend

if args.terse:
    print('{:0.3e}'.format(err_evals))
    print('{:0.3e}'.format(err_ortho1))
    print('{:0.3e}'.format(err_ortho2))
    print('{:0.3e}'.format(err_cross))
    print('{:0.3e}'.format(err_ham))
else:
    print('Maximum error in:')
    print('- eigenvalues |E_new - E_ref|: {:0.3e} Ry'.format(err_evals))
    print('- ref. wavefunctions orthonormality |C_ref^H C_ref - 1|: {:0.3e}'.format(err_ortho1))
    print('- new wavefunctions orthonormality |C_new^H C_new - 1|: {:0.3e}'.format(err_ortho2))
    print('- cross wavefunctions orthonormality |C_new^H C_new - C_ref^H C_ref|: {:0.3e}'.format(err_cross))
    print('- reassembled Hamiltonian |H_ref - H_new|/(|H_ref||H_new|): {:0.3e}'.format(err_ham))
