# periodic.py --- Author: Zhenglu Li, lzlwell@berkeley.edu
# periodic_system class, central object containing wavefunctions and system info and operations

import numpy as np
from copy import copy
import os
import pickle

import common
import wannier as wan
import planewaves as pw

class periodic_system:
  '''
  Define a periodic system and its wavefunctions
  '''

  def __init__(self, prefix, latvec, ecut):
    ''' Initialization '''

    self.prefix = prefix

    # lattice vector in bohr
    self.latvec = np.zeros((3,3), dtype=np.float)
    self.latvec[:,:] = copy(latvec[:,:])

    # reciprocal lattice vector
    self.recipvec = np.zeros((3,3), dtype=np.float)
    self._get_recipvec()

    # energy cutoff for wavefunctions
    self.ecut = ecut

    self.nkpt = 0 # number of kpoints
    self.fftgrid = np.zeros(3, dtype=np.int)  # FFT grid

    return


## public functions
  def print_info(self):
    ''' Print information '''

    print()
    print('=== Information ===')
    print('system name:', self.prefix)
    print('lattice vectors:')
    print(self.latvec)
    print('reciprocal lattice vectors:')
    print(self.recipvec)
    print('energy cutoff:', self.ecut, 'Ry')
    print('number of kpoints:', self.nkpt)
    print('number of bands:', self.nband)
    print('FFT grid:', self.fftgrid)
    print('===================')
    print()

    return


  def read_kpoints_list(self):
    name = 'kpoints.list'

    kpoints = self._read_points_list(name)

    self.nkpt = len(kpoints)
    self.kpoints = np.zeros((self.nkpt, 3), dtype=float)
    self.kpoints[:,:] = copy(kpoints[:,:])

    return


  def write_kpoints_list(self, nsplit=0):

    kpath = 'kpts_split/'

    try:
      os.stat(kpath)
    except:
      os.mkdir(kpath)

    name = kpath + self.prefix + '_kpoints'
    self._write_points_list(self.kpoints, name)

    if nsplit > 0:
      # split kpoints to write to different files
      # mainly used for wannier-interpolation manual parallelization
      nktot = len(self.kpoints)
      if nsplit > nktot:
        raise ValueError('write_kpoints_list: nsplit too large.')

      if nktot%nsplit == 0:
        nk = np.int(nktot/nsplit)

        for isplit in range(nsplit):
          name = kpath + self.prefix + '_kpoints' + '.split_' + str(isplit+1)
          self._write_points_list(self.kpoints[nk*isplit:nk*(isplit+1), :], name)
      else:
        nk = np.int(nktot/nsplit)

        for isplit in range(nsplit - 1):
          name = kpath + self.prefix + '_kpoints' + '.split_' + str(isplit+1)
          self._write_points_list(self.kpoints[nk*isplit:nk*(isplit+1), :], name)

        # write last isplit, all the remaining kpionts
        isplit = nsplit -1
        name = kpath + self.prefix + '_kpoints' + '.split_' + str(isplit+1)
        self._write_points_list(self.kpoints[nk*isplit:, :], name)

    return


  def read_dft_unk_wavefuns(self):
    '''
    Read DFT-UNK wavefunctions
    '''

    wpath = common.wpath
 
    wvfns = []
 
    if common.wvfn_disk_io:
      try:
        os.stat('wvfn_save/')
      except:
        os.mkdir('wvfn_save/')

    for ik in range(self.nkpt):
      kidx = ik + 1  # k index for file labels
 
      print('Reading DFT-UNK wavefunctions k =', kidx, '/', self.nkpt)

      f_unk = wpath + 'UNK' + '{:05}'.format(kidx) + '.1'
 
      # different DFT codes normalize to different volume
      # here we use default volume 1.0 and then normalize it
      unk_data, fftgrid = wan.load_u_file(f_unk)

      nband = len(unk_data)

      for ib in range(nband):
        # normalize to exact 1.0
        norm = wan.check_norm(unk_data[ib,:])
        unk_data[ib,:] /= np.sqrt(norm)

      if common.wvfn_disk_io:
        with open('wvfn_save/wvfn_r_ik_' + str(ik) + '.pkl', 'wb') as output:
          pickle.dump(unk_data, output, pickle.HIGHEST_PROTOCOL)
      else:
        wvfns.append(unk_data)

    # store unk self.wavefuns_r array
    if common.wvfn_disk_io:
      self.wavefuns_r = None  # stored on disk
    else:
      self.wavefuns_r = np.stack(wvfns)
      wvfns = None

    # save nband
    self.nband = nband

    # save fftgrid
    self.fftgrid = copy(fftgrid)

    # load eigenvalues
    eigvals = np.loadtxt(wpath + self.prefix + '.eig')
    if len(eigvals) != self.nkpt * self.nband:
      raise ValueError('read_dft_unk_wavefuns: .eig file dimension mismatch.')

    # store eigenvalues, without Fermi shift
    self.eigenvalues = np.zeros((self.nkpt, self.nband), dtype=np.float)

    count = 0
    for ik in range(self.nkpt):
      for ib in range(self.nband):
        # store eigenvalues, without Fermi shift
        self.eigenvalues[ik, ib] = eigvals[count, 2]
        count += 1

    return


  def read_wannier_wavefuns(self, nsplit=0):
    ''' 
    Read wannier interpolated wavefunctions

    Load UWK files and _interp_rotmat.dat files, and rotate
    uwk to unk
    '''

    wpath = common.wpath

    # wavefunctions in real-space FFT grid basis
    # basis fftgrid points stored in self.fftgrid
    wvfns = []
    eigvals = []

    # set up parameters for split case
    if nsplit > self.nkpt:
      raise ValueError('read_wannier_wavefuns: nsplit too large.')

    if common.wvfn_disk_io:
      try:
        os.stat('wvfn_save/')
      except:
        os.mkdir('wvfn_save/')

    nksplit = 0
    if nsplit > 0:
      nksplit = np.int(self.nkpt/nsplit)

    isplit = 0
    iksplit = 0
    for ik in range(self.nkpt):
      kidx = ik + 1  # k index for file labels

      progress = np.int(100 * kidx / self.nkpt)
      if (progress % 10 == 0) and (progress > 0):
        print('Reading Wannier-interpolated wavefunctions', progress, '%')

      # read wannier-interp files, check whether klist was split
      if nsplit == 0:
        f_uwk = wpath + 'wannier-interp/' + 'UWK' + '{:05}'.format(kidx) + '.interp'

        # read uwk data and coefficients, and then rotate to unk
        uwk_data, fftgrid = wan.load_u_file(f_uwk)  # contains all wannier basis
        unk_winterp, eig_winterp = wan.uwk_to_unk_winterp(wpath + 'wannier-interp/' + self.prefix, kidx, uwk_data)

      elif nsplit > 0:
        f_uwk = wpath + 'wannier-interp-split_' + str(isplit+1) \
                        + '/UWK' + '{:05}'.format(iksplit+1) + '.interp'

        # read uwk data and coefficients, and then rotate to unk
        uwk_data, fftgrid = wan.load_u_file(f_uwk)  # contains all wannier basis
        unk_winterp, eig_winterp = wan.uwk_to_unk_winterp(wpath + 'wannier-interp-split_' + str(isplit+1) \
                                                      + '/' + self.prefix, iksplit+1, uwk_data)

        # update the looper
        iksplit += 1
        # for first nsplit-1 folders, we check when we finish reading the folder
        if isplit < nsplit - 1:
          # done reading this folder
          if iksplit%nksplit == 0:
            isplit += 1
            iksplit = 0

      nband = len(unk_winterp)

      # normalize to exact 1.0
      for ib in range(nband):
        norm = wan.check_norm(unk_winterp[ib,:])
        unk_winterp[ib,:] /= np.sqrt(norm)

      if common.wvfn_disk_io:
        with open('wvfn_save/wvfn_r_ik_' + str(ik) + '.pkl', 'wb') as output:
          pickle.dump(unk_winterp, output, pickle.HIGHEST_PROTOCOL)
      else:
        wvfns.append(unk_winterp)

      eigvals.append(eig_winterp)

    # store unk self.wavefuns_r array
    if common.wvfn_disk_io:
      self.wavefuns_r = None  # stored on disk
    else:
      self.wavefuns_r = np.stack(wvfns)
      wvfns = None

    # save nband
    self.nband = nband

    # save fftgrid
    self.fftgrid = copy(fftgrid)

    # store eigenvalues, without Fermi shift
    self.eigenvalues = np.stack(eigvals)
    eigvals = None

    return


  def unk_to_planewaves(self):
    ''' Do Fourier transform to convert u(r) to c(G) '''

    # get ngk and all gvecs, and the needed gvecs' indices, based on energy cutoff
    ngk, gspace, gvecs_idx_k = pw.get_ngk_gvecs(self.fftgrid, self.kpoints, self.ecut, self.recipvec)

    # save ngk and gspace
    self.ngk = copy(ngk)
    self.gspace = copy(gspace)

    # prepare self.wavefuns_g, this is a list since different k-points have different ngk
    self.wavefuns_g = []

    # prepare self.gvecs, this is a list since different k-points have different ngk
    self.gvecs = []

    # double check
    if common.wvfn_disk_io:
      try:
        os.stat('wvfn_save/')
      except:
        os.mkdir('wvfn_save/')

    # ik loop
    for ik in range(self.nkpt):
      print('Converting wavefunctions k =', ik + 1, '/', self.nkpt)

      # get unk(r) for ik
      if common.wvfn_disk_io:
        unkr = self.load_wvfn_from_disk('r', ik)
      else:
        unkr = copy(self.wavefuns_r[ik])

      # get gvecs_k for ik
      gkidx = gvecs_idx_k[sum(ngk[:ik]) : sum(ngk[:ik+1])]
      gvecs_k = gspace[gkidx]

      if len(gvecs_k) != ngk[ik]:
        raise ValueError('unk_to_planewaves: gvecs_k length and ngk mismatch.')

      cg = pw.unkr_fourier_transform_cg(unkr, self.fftgrid, gvecs_k, self.latvec, self.recipvec)

      if common.wvfn_disk_io:
        with open('wvfn_save/wvfn_g_ik_' + str(ik) + '.pkl', 'wb') as output:
          pickle.dump(cg, output, pickle.HIGHEST_PROTOCOL)
      else:
        self.wavefuns_g.append(cg)

      self.gvecs.append(gvecs_k)

    return


  def load_wvfn_from_disk(self, mode, ik):
    ''' Load wavefunction from disk '''

    if not common.wvfn_disk_io:
      raise ValueError('load_wvfn_from_disk: this function cannot be called if common.wvfn_disk_io is False.')

    if (mode != 'r') and (mode != 'g'):
      raise ValueError("load_wvfn_from_disk only takes mode value of 'r' or 'g'.")

    # mode = 'r' or 'g'
    with open('wvfn_save/wvfn_' + mode + '_ik_' + str(ik) + '.pkl', 'rb') as inp:
      wvfn_ik = pickle.load(inp)

    return wvfn_ik

 
## private functions

  def _get_recipvec(self):
    self.recipvec = np.linalg.inv(self.latvec)
    self.recipvec = self.recipvec.transpose()

    # note that with the 2\pi factor it is the standard definition,
    # however different from some dft codes
    # here we use the standard definition
    self.recipvec *= 2.0*np.pi

    return


  def _write_points_list(self, points, name):
    ''' Write data points list into file '''

    with open(name, 'w') as ptf:
      npts = len(points)

      ptf.write(str(npts) + '\n')
      for i in range(npts):
        ptf.write(f'{points[i,0]:18.12f}' + f'{points[i,1]:18.12f}' + f'{points[i,2]:18.12f}' + '\n')

    return


  def _read_points_list(self, name):
    ''' Read data points list from file '''

    with open(name, 'r') as ptf:
      lines = []
      for line in ptf:
        lines.append(line)

      ln = lines[0].split()
      npts = int(ln[0])

      points = np.zeros((npts, 3), dtype=np.float)

      for ip in range(npts):
        ln = lines[ip+1].split()

        points[ip,0] = float(ln[0])
        points[ip,1] = float(ln[1])
        points[ip,2] = float(ln[2])

    return points


  def _crystal_volume(self):
    ''' Crystal colume in Ang^3 '''
    a = self.latvec[0,:]
    b = self.latvec[1,:]
    c = self.latvec[2,:]

    volume = np.dot( np.cross(a, b), c )

    if volume < 0.0:
      raise ValueError('_crystal_volume: negative crystal volume.')

    return volume

# end of class periodic_system

