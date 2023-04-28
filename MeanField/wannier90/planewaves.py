
import numpy as np
from copy import copy
import wannier as wan


def get_ngk_gvecs(fftgrid, klist, ecut, recipvec):
  ''' Get the array ngk(nk) '''

  # in rydberg atomic unit
  # length is in bohr, energy is in rydberg
  # hbar = 2m = 1
  # therefore kinetic energy = (k+G)^2 in Ry, where (k+G) is in reciprocal cartesian coordinate

  nk = len(klist)
  ngk = np.zeros(nk, dtype=np.int)

  print('pw: generating all planewaves')

  ngx = fftgrid[0]
  ngy = fftgrid[1]
  ngz = fftgrid[2]

  # find global G space
  #
  gvs = []

  for igx in range(-int(ngx/2), int(ngx/2)+1):
    for igy in range(-int(ngy/2), int(ngy/2)+1):
      for igz in range(-int(ngz/2), int(ngz/2)+1):
        gv = np.array([igx, igy, igz], dtype=np.int)

        gvs.append(gv)

  gvs_pool = np.stack(gvs)

  # convert G from fraction to cartesian
  gv_c = frac_to_cart(recipvec, gvs_pool)

  # compute G^2
  ekin = (gv_c * gv_c).sum(axis=1)

  # sort kinetic energy from small to large
  isort = np.argsort(ekin)

  # find cutoff index in ordered array
  icut = len(ekin)
  for ig in range(len(ekin)):
    # keep gvecs with kinetic energy < ecut_rho = ecut * 4.0
    if ekin[isort[ig]] > ecut * 4.0:
      icut = ig
      break

  gspace = copy(gvs_pool[isort[:icut]])  # all G-vectors needed in calculation

  # find k-dependent G sphere
  #
  gkidx = []

  for ik in range(nk):
    print('pw: generating ngk and gkvecs for k =', ik + 1, '/', nk)

    # k+G
    kgvecs = klist[ik,:] + gspace[:,:]

    # convert k+G from fraction to cartesian
    kgv_c = frac_to_cart(recipvec, kgvecs)

    # computed all (k+G)^2
    ekin = (kgv_c * kgv_c).sum(axis=1)

    # indices of which the kinetic energy from (k+G) is smaller than ecut
    idx = np.where(ekin < ecut)[0]

    # number of G-vectors for every k point
    ngk[ik] = len(idx)

    # save idx
    gkidx.append(idx)

  # flatten the gkidx
  gvecs_idx_k = np.concatenate(gkidx, axis=None)

  if len(gvecs_idx_k) != sum(ngk):
    raise ValueError('Length of gkvecs_idx and sum of ngk mismatch.')

  print('pw: total number of g-vectors', len(gspace))
  print('pw: ngk array')
  print(ngk)
  print('pw: gspace')
  print(gspace)
  print()

  return ngk, gspace, gvecs_idx_k


def unkr_fourier_transform_cg(unkr, fftgrid, gvecs_k, latvec, recipvec):
  ''' Fourier transform of u_nk(r) to c(G) for single k-point '''

  # c_g = \sum_r e^{-i g.r} u(r) / N_fft

  nband = len(unkr)  # unk(r) here contains only one k-point

  ng = len(gvecs_k)  # ngk for one k-point
  cg = np.zeros((nband, ng), dtype=np.complex)

  #print('pw-ft: generating all real-space grid points')

  ngx = fftgrid[0]
  ngy = fftgrid[1]
  ngz = fftgrid[2]
  nfft = ngx * ngy * ngz

  rvs = []

  # fastest axis is x, then y, and z is the slowest axis
  # wannier90 UNK/UWK convention
  for igz in range(ngz):
    for igy in range(ngy):
      for igx in range(ngx):
        rv = np.array([igx/ngx, igy/ngy, igz/ngz], dtype=np.float)
        rvs.append(rv)

  rvecs = np.stack(rvs)  # all r-vectors considered

  # convert r gird points from fraction to cartesian
  rvecs_c = frac_to_cart(latvec, rvecs)

  # convert gkvecs from fraction to cartesian
  gvecs_k_c = frac_to_cart(recipvec, gvecs_k)

  #print('pw-ft: performing fourier transform')
  for ig in range(ng):
    # get phase factors exp(-i g.r)
    g_dot_r = np.dot(rvecs_c, gvecs_k_c[ig])
    phase = np.exp(-1.0j * g_dot_r)

    for ib in range(nband):
      cg[ib, ig] = np.dot(phase[:], unkr[ib,:]) / nfft

  # normalize
  for ib in range(nband):
    # planewave normalization sum(cg*cg) = 1
    norm = np.sum(cg[ib,:].conjugate() * cg[ib,:])
    if abs(norm - 1.0) > 0.1:
      print('Warning: norm of planewave wavefunction deviates from unity. norm =', norm.real)

    cg[ib,:] /= np.sqrt(norm.real)

  return cg


def frac_to_cart(cellvec, frac_coord):
  ''' fractional to cartesian coordinate conversion '''

  cart_coord = np.zeros(np.shape(frac_coord), dtype=float)

  if len(np.shape(frac_coord)) == 2:
    for i in range(len(frac_coord)):
      cart_coord[i,0] = frac_coord[i,0]*cellvec[0,0] + frac_coord[i,1]*cellvec[1,0] + frac_coord[i,2]*cellvec[2,0]
      cart_coord[i,1] = frac_coord[i,0]*cellvec[0,1] + frac_coord[i,1]*cellvec[1,1] + frac_coord[i,2]*cellvec[2,1]
      cart_coord[i,2] = frac_coord[i,0]*cellvec[0,2] + frac_coord[i,1]*cellvec[1,2] + frac_coord[i,2]*cellvec[2,2]
  elif len(np.shape(frac_coord)) == 1:
    cart_coord[0] = frac_coord[0]*cellvec[0,0] + frac_coord[1]*cellvec[1,0] + frac_coord[2]*cellvec[2,0]
    cart_coord[1] = frac_coord[0]*cellvec[0,1] + frac_coord[1]*cellvec[1,1] + frac_coord[2]*cellvec[2,1]
    cart_coord[2] = frac_coord[0]*cellvec[0,2] + frac_coord[1]*cellvec[1,2] + frac_coord[2]*cellvec[2,2]
  else:
    raise ValueError('frac/cart conversion: coord array dimension is incorrect.')

  return cart_coord


def cart_to_frac(cellvec, cart_coord):
  ''' cartesian to fractional coordinate conversion'''

  frac_coord = np.zeros(np.shape(cart_coord), dtype=float)

  invcellvec = np.linalg.inv(cellvec)
  frac_coord = frac_to_cart(invcellvec, cart_coord)

  return frac_coord

