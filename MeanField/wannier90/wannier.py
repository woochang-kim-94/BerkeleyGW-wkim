# wannier.py  ---  Author: Zhenglu Li, lzlwell@berkeley.edu
# handles wannier90 preprocessing UNK files and 
# wannier90-interpolated UWK files, with special modifications

import numpy as np


def load_u_file(fname, volume=1.0):
  ''' Read UNK and UWK file'''

  with open(fname) as ufile:
    lines = []
    for line in ufile:
      lines.append(line)
    
    ln = lines[0].split() # first line
    ngx = int(ln[0])
    ngy = int(ln[1])
    ngz = int(ln[2])
    ikpoint = int(ln[3])
    nstates = int(ln[4]) # number of bloch or wannier states

    # now we declare udata on FFT grid
    udata = np.zeros((nstates, ngz*ngy*ngx), dtype=np.complex)

    iline = 0
    # read data from following lines
    for ist in range(nstates):
      count = 0
      for igz in range(ngz):
        for igy in range(ngy):
          for igx in range(ngx):
            iline += 1
            ln = lines[iline].split()

            udata[ist, count] = float(ln[0]) + 1.0j*float(ln[1])
            count += 1

  # some DFT codes generated (e.g., Siesta) UNK file is directly normalized to 1/V, where V = cell volume in Ang**3
  # such that \int dV u^* u = 1.0
  # we recover this factor if volume is speficied, consistent with wannier convention
  udata[:,:] = udata[:,:] * np.sqrt(volume)

  fftgrid = np.zeros(3, dtype=np.int)
  fftgrid[0] = ngx
  fftgrid[1] = ngy
  fftgrid[2] = ngz

  #print('File', fname, 'has been successfully loaded. Averaged norm is ', check_norm(udata))

  return udata, fftgrid


def check_norm(udata):
  '''Check the normalizaiton of the unk or uwk data'''

  norm = np.average(udata.conjugate()*udata, axis=None)  # average over all elements
  if abs(norm - 1.0) > 0.1:
    print('Warning: norm of u(r) wavefunction deviates from unity. norm =', norm.real)

  return norm.real


def unk_overlap(unk_bra, unk_ket, fftgrid_mask=None):
  '''Compute <unk_bra|unk_ket>, both unk data contain only one band'''

  nfft = len(unk_bra)

  if np.shape(unk_bra) != np.shape(unk_ket):
    raise ValueError('unk_bra and unk_ket have different dimensions.')

  if fftgrid_mask is None:
    ovlp = np.sum(unk_bra.conjugate()*unk_ket, axis=None)  # sum over all axes
  else:
    if np.shape(fftgrid_mask) != np.shape(unk_ket):
      raise ValueError('fftgrid_mask and unk_ket have different dimensions.')
    ovlp = np.sum(unk_bra.conjugate()*unk_ket*fftgrid_mask, axis=None)  # sum over all axes

  ovlp /= nfft

  return ovlp


def uwk_to_unk_winterp(prefix, kidx, uwk_data):
  '''Rotate uwk in wannier basis to unk_winterp in bloch basis
     read prefix_interp_rotmat.dat for rotation coefficients'''

  nwann = len(uwk_data)
  nband_wann = nwann # in the wannier subspace, number of bands equals number of wannier functions

  unk_winterp = np.zeros(np.shape(uwk_data), dtype=np.complex)
  eig_winterp = np.zeros(nband_wann, dtype=np.float)

  with open(prefix + '_interp_rotmat.dat') as finterp:
    lines = []
    for line in finterp:
      lines.append(line)

    for iline in range(len(lines)):
      ln = lines[iline].split()

      if len(ln) == 2:
        if ln[0] == 'ikp:' and int(ln[1]) == kidx:
          for ib in range(nband_wann):
            iline += 4

            # read eigenvalue
            lneig = lines[iline].split()
            eig_winterp[ib] = float(lneig[0])

            iline += 1
            for iw in range(nwann):
              iline += 1
  
              lncoef = lines[iline].split()
              coeff = float(lncoef[0]) + 1.0j * float(lncoef[1])
  
              unk_winterp[ib,:] += coeff * uwk_data[iw,:]

  return unk_winterp, eig_winterp

