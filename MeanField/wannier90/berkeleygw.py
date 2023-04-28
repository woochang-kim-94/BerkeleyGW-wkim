# berkeleygw.py --- Author: Zhenglu Li, lzlwell@berkeleygw.edu
# handles berkeleygw wavefunctions I/O

import numpy as np
import h5py
from copy import copy

import common
from periodic import periodic_system


def create_wfn(wfn_in, wfn_out, wavef):
  ''' Create wannier-converted wfn_out using template wfn_in and wavef data '''

  print('Creating Wannier-converted wavefunction in BerkeleyGW hdf5 format')

  f_in = h5py.File(wfn_in, 'r')
  f_out = h5py.File(wfn_out, 'w')

  print()
  print('reading from', wfn_in)
  print('writing', wfn_out)
  print()

  # check consistency
  wfn_consistency_check(f_in, wavef)

  ##
  # copy header
  f_out.copy(f_in['mf_header'], 'mf_header')
  f_in.close()

  # modify gspace
  f_out['mf_header/gspace/components'][()] = wavef.gspace[:,:]

  # modify kpoints
  #
  # band energies, from wannier (eV) to berkeleygw (rydberg)
  del f_out['mf_header/kpoints/el']

  # add spin dimension, and only consider nspin=1 case
  nspin = 1
  eigvals = np.zeros((nspin, wavef.nkpt, wavef.nband), dtype=np.float)
  eigvals[0,:,:] = copy(wavef.eigenvalues[:,:]) / common._rydberg2eV
  f_out.create_dataset('mf_header/kpoints/el', data=eigvals)

  # ifmax
  del f_out['mf_header/kpoints/ifmax']

  ifmax = np.zeros(wavef.nkpt, dtype=np.int)
  for ik in range(wavef.nkpt):
    for ib in range(wavef.nband):
      if wavef.eigenvalues[ik, ib] < common.efermi + common._tiny:
        ifmax[ik] += 1
  
  f_out.create_dataset('mf_header/kpoints/ifmax', data=ifmax)

  # ifmin
  del f_out['mf_header/kpoints/ifmin']

  ifmin = np.ones(wavef.nkpt, dtype=np.int)
  f_out.create_dataset('mf_header/kpoints/ifmin', data=ifmin)

  # kgrid, set kgrid to zero for general patch/fine/high-sym-line sampling
  f_out['mf_header/kpoints/kgrid'][:] = 0

  # mnband
  f_out['mf_header/kpoints/mnband'][()] = wavef.nband

  # ngk
  del f_out['mf_header/kpoints/ngk']
  f_out.create_dataset('mf_header/kpoints/ngk', data=wavef.ngk)

  # ngkmax
  f_out['mf_header/kpoints/ngkmax'][()] = max(wavef.ngk)

  # nrk, number of kpoints
  f_out['mf_header/kpoints/nrk'][()] = wavef.nkpt

  # occ
  del f_out['mf_header/kpoints/occ']

  # consider only nspin = 1 case
  nspin = 1
  occ = np.zeros((nspin, wavef.nkpt, wavef.nband), dtype=np.float)
  for ik in range(wavef.nkpt):
    for ib in range(wavef.nband):
      if wavef.eigenvalues[ik, ib] < common.efermi - common._small:
        occ[0, ik, ib] = 1.0

      elif wavef.eigenvalues[ik, ib] > common.efermi + common._small:
        occ[0, ik, ib] = 0.0
      
      else:
        occ[0, ik, ib] = 0.5  # for graphene-like case

  f_out.create_dataset('mf_header/kpoints/occ', data=occ)

  # rk, kpoints coordinates in fractional coordinates
  del f_out['mf_header/kpoints/rk']
  f_out.create_dataset('mf_header/kpoints/rk', data=wavef.kpoints)

  # w, weights of kpoints, always set to 1/N_k, not used by BerkeleyGW
  del f_out['mf_header/kpoints/w']

  weights = np.zeros(wavef.nkpt, dtype=np.float)
  weights[:] = 1.0 / wavef.nkpt

  f_out.create_dataset('mf_header/kpoints/w', data=weights)

  ##
  # create wfns
  f_out.create_group('wfns')

  # coeffs and gvecs, consider only nspinor = 1 case
  nspinor = 1
  ngk = copy(wavef.ngk)

  # last index is size 2 for real and imaginary parts respectively
  coeffs = np.zeros((wavef.nband, nspinor, sum(ngk), 2), dtype=np.float)
  gvecs = np.zeros((sum(ngk), 3), dtype=np.int)

  for ik in range(wavef.nkpt):
    # get cg for ik
    if common.wvfn_disk_io:
      cg = wavef.load_wvfn_from_disk('g', ik)
    else:
      cg = copy(wavef.wavefuns_g[ik])

    coeffs[:, 0, sum(ngk[:ik]) : sum(ngk[:ik+1]), 0] = cg[:,:].real
    coeffs[:, 0, sum(ngk[:ik]) : sum(ngk[:ik+1]), 1] = cg[:,:].imag

    gvecs[sum(ngk[:ik]) : sum(ngk[:ik+1]), :] = wavef.gvecs[ik]

  f_out.create_dataset('wfns/coeffs', data=coeffs)
  f_out.create_dataset('wfns/gvecs', data=gvecs)

  ##
  f_out.close()

  return


def wfn_consistency_check(fwfn, wavef):
  ''' Check consistency between wfn and wavef and other aspects'''

  # lattice vectors
  if diff_arrays(fwfn['mf_header/crystal/avec'][()], wavef.latvec) > common._tiny:
    raise ValueError('bgw: inconsistent lattice vectors between fwfn and wavef.')

  # flavor
  if fwfn['mf_header/flavor'][()] != 2:
    raise ValueError('bgw: wavefunction template must be in complex version with flavor = 2.')

  # fft grid
  if diff_arrays(fwfn['mf_header/gspace/FFTgrid'][()], wavef.fftgrid) > common._tiny:
    raise ValueError('bgw: inconsistent FFTgrid between fwfn and wavef.')

  # gspace size
  if fwfn['mf_header/gspace/ng'][()] != len(wavef.gspace):
    raise ValueError('bgw: inconsistent size of gspace between fwfn and wavef.')

  # cutoff
  if fwfn['mf_header/kpoints/ecutwfc'][()] != wavef.ecut:
    raise ValueError('bgw: inconsistent wavefunction energy cutoff between fwfn and wavef.')

  # nspin and nspinor
  if (fwfn['mf_header/kpoints/nspin'][()] != 1) or (fwfn['mf_header/kpoints/nspinor'][()] != 1):
    raise ValueError('bgw: currently only supports nspin = 1 and nspinor = 1 case.')

  return


def diff_arrays(a, b):
  ''' Compute the absolute difference between all elements of two arrays '''

  diffab= a - b
  diff = np.sum(np.absolute(diffab))

  if diff > common._tiny:
    print('Warning: difference between arrays too large')
    print(a)
    print(b)

  return diff
