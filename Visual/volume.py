#!/usr/bin/env python

#-------------------------------------------------------------------------------
#
#   volume.py
#   converts volumetric data from a3Dr to Gaussian Cube or XCrySDen XSF
#   written by Georgy Samsonidze (November 2008)
#
#-------------------------------------------------------------------------------

from common import format_index2mat, boolean_input
from matter import mat_check, mat_species, mat_format, mat_replicate, mat_read, paratec_read, vasp_read, espresso_read, siesta_read, tbpw_read, xyz_read, xsf_read

def main(argv = None):
   if argv is None:
      argv = sys.argv
   if not (len(argv) == 9 or len(argv) == 10):
      sys.exit("\n   Usage: volume.py imfn imff ivfn ivff ovfn ovff phas cplx [hole]\n" +
               "   imfn = input matter file name\n" +
               "   imff = input matter file format\n" +
               "          (mat|paratec|vasp|espresso|siesta|tbpw|xyz|xsf)\n" +
               "   ivfn = input volumetric file name\n" +
               "   ivff = input volumetric file format (a3dr)\n" +
               "   ovfn = output volumetric file name\n" +
               "   ovff = output volumetric file format (cube|xsf)\n" +
               "   phas = remove wavefunction phase (true|false)\n" +
               "   cplx = complex wavefunction (re|im|abs|abs2)\n" +
               "   hole = plot hole (true|false)\n")
   imfn = argv[1]
   imff = argv[2].lower()
   ivfn = argv[3]
   ivff = argv[4].lower()
   ovfn = argv[5]
   ovff = argv[6].lower()
   if imff == 'mat':
      ierr, mat = mat_read(imfn)
   elif imff == 'paratec':
      ierr, mat = paratec_read(imfn)
   elif imff == 'vasp':
      ierr, mat = vasp_read(imfn)
   elif imff == 'espresso':
      ierr, mat = espresso_read(imfn)
   elif imff == 'siesta':
      ierr, mat = siesta_read(imfn)
   elif imff == 'tbpw':
      ierr, mat = tbpw_read(imfn)
   elif imff == 'xyz':
      ierr, mat = xyz_read(imfn)
   elif imff == 'xsf':
      ierr, mat = xsf_read(imfn)
   else:
      sys.exit("\n   Error: invalid input format %s\n" % imff)
   if ierr != 0:
      sys.exit("\n   Error: unable to read file %s\n" % imfn)
   ierr = mat_check(mat)
   if ierr != 0:
      sys.exit("\n   Error: invalid content in file %s\n" % imfn)
   if ivff == 'a3dr':
      ierr, state, energy, holpos, latvec, grid, psi = a3dr_read(ivfn)
   else:
      sys.exit("\n   Error: invalid input format %s\n" % ivff)
   if ierr != 0:
      sys.exit("\n   Error: unable to read file %s\n" % ivfn)
   format = format_index2mat[0]
   so = [0.0, 0.0, 0.0]
# temporary hack for a3dr format
# lattice vectors in a3dr file are written in %.6f format
# in paratec benzene example, lattice vectors read from a3dr file
# differ from lattice vectors read from paratec input file
# this results in identical atoms written to rho.cube
# (total of 16 atoms instead of 12 atoms)
# this hack computes the supercell dimensions and
# replaces the lattice vectors read from a3dr file
# with the lattice vectors read from the input file
# multiplied by the supercell dimensions
   mat = mat_format(mat, mat['fc'], mat['fo'], format, mat['fp'])
   a1 = mat['lv'][0]
   a2 = mat['lv'][1]
   a3 = mat['lv'][2]
   b1 = latvec[0]
   b2 = latvec[1]
   b3 = latvec[2]
   n1 = int(round((a1[0] * b1[0] + a1[1] * b1[1] + a1[2] * b1[2]) / (a1[0]**2 + a1[1]**2 + a1[2]**2)))
   n2 = int(round((a2[0] * b2[0] + a2[1] * b2[1] + a2[2] * b2[2]) / (a2[0]**2 + a2[1]**2 + a2[2]**2)))
   n3 = int(round((a3[0] * b3[0] + a3[1] * b3[1] + a3[2] * b3[2]) / (a3[0]**2 + a3[1]**2 + a3[2]**2)))
   b1 = [n1 * a1[0], n1 * a1[1], n1 * a1[2]]
   b2 = [n2 * a2[0], n2 * a2[1], n2 * a2[2]]
   b3 = [n3 * a3[0], n3 * a3[1], n3 * a3[2]]
   latvec = [b1, b2, b3]
# end of hack
   sv = latvec
   smin, smax, si, mat = mat_replicate(mat, so, sv, format)
   if mat['fo'] != format and mat['fv'] != format and mat['fp'] != format:
      sys.exit("\n   Error: unable to make supercell\n")
   if len(argv) > 9 and holpos != []:
      ierr, hole = boolean_input(argv[9])
      if ierr != 0:
         sys.exit("\n   Error: invalid hole %s\n" % argv[9])
      if hole != 0:
         mat['as'].append(0)
         mat['ap'].append(holpos)
         mat['na'] += 1
         mat['ns'] = len(mat_species(mat['na'], mat['as']))
   ierr, phas = boolean_input(argv[7])
   if ierr != 0:
      sys.exit("\n   Error: invalid phase %s\n" % argv[7])
   if phas != 0:
      rmax = 0.0
      for ii in range(grid[0]):
         for jj in range(grid[1]):
            for kk in range(grid[2]):
               rr = math.sqrt(math.pow(psi[kk][jj][ii][0], 2) + math.pow(psi[kk][jj][ii][1], 2))
               if rr > rmax:
                  rmax = rr
                  phase = math.atan2(psi[kk][jj][ii][1], psi[kk][jj][ii][0])
      cs = math.cos(phase)
      sn = math.sin(phase)
      for ii in range(grid[0]):
         for jj in range(grid[1]):
            for kk in range(grid[2]):
               xx = psi[kk][jj][ii][0]
               yy = psi[kk][jj][ii][1]
               psi[kk][jj][ii][0] = cs * xx + sn * yy
               psi[kk][jj][ii][1] = -sn * xx + cs * yy
   ierr, cplx = cplx_input(argv[8])
   if ierr != 0:
      sys.exit("\n   Error: invalid cplx %s\n" % argv[8])
   phi = []
   for ii in range(grid[0]):
      phi.append([])
      for jj in range(grid[1]):
         phi[ii].append([])
         for kk in range(grid[2]):
            if cplx == 0:
               chi = psi[kk][jj][ii][0]
            elif cplx == 1:
               chi = psi[kk][jj][ii][1]
            elif cplx == 2:
               chi = math.sqrt(math.pow(psi[kk][jj][ii][0], 2) + math.pow(psi[kk][jj][ii][1], 2))
            elif cplx == 3:
               chi = math.pow(psi[kk][jj][ii][0], 2) + math.pow(psi[kk][jj][ii][1], 2)
            else:
               chi = 0.0
            if phas != 0 and (cplx == 2 or cplx == 3) and psi[kk][jj][ii][0] < 0.0:
               chi = -chi
            phi[ii][jj].append(chi)
   if ovff == 'cube':
      ii = 0
   elif ovff == 'xsf':
      ii = 1
   else:
      sys.exit("\n   Error: invalid output format %s\n" % ovff)
   fc = format_index2mat[ii]
   fo = format_index2mat[ii]
   fv = format_index2mat[ii]
   fp = format_index2mat[ii]
   mat = mat_format(mat, fc, fo, fv, fp)
   if mat['fc'] != fc or mat['fo'] != fo or mat['fv'] != fv or mat['fp'] != fp:
      sys.exit("\n   Error: invalid content in file %s\n" % imfn)
   na = mat['na']
   aspecies = []
   for ii in range(na):
      aspecies.append(mat['as'][ii])
   acharge = []
   for ii in range(na):
      acharge.append(float(mat['as'][ii]))
   aposition = []
   for ii in range(na):
      aposition.append([])
      for jj in range(3):
         aposition[ii].append(mat['ap'][ii][jj])
   origin = []
   for jj in range(3):
      origin.append(mat['lo'][jj])
   latvec = []
   for ii in range(3):
      latvec.append([])
      for jj in range(3):
         latvec[ii].append(mat['lv'][ii][jj])
   step = []
   for ii in range(3):
      step.append([])
      for jj in range(3):
         step[ii].append(latvec[ii][jj] / float(grid[ii]))
   if ovff == 'cube':
      ierr = cub_write(ovfn, state, energy, na, aspecies, acharge, aposition, grid, origin, step, phi)
   elif ovff == 'xsf':
      ierr = xsf_write(ovfn, state, energy, na, aspecies, aposition, grid, origin, latvec, phi)
   else:
      sys.exit("\n   Error: invalid output format %s\n" % ovff)
   if ierr != 0:
      sys.exit("\n   Error: unable to write file %s\n" % ovfn)
   return 0

def a3dr_read(file):
   ierr = 0
   exciton_state = 0
   exciton_energy = 0.0
   hole_position = []
   lattice_vectors = []
   grid_size = []
   wavefunction = []
   try:
      hh = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      while True:
         ss = hh.readline()
         tt = ss.split()
         tt.append('')
         if tt[1] == 'ie':
            exciton_state = int(tt[3])
         if tt[1] == 'e':
            exciton_energy = float(tt[3])
         if tt[1] == 'rh':
            for jj in range(3):
               hole_position.append(float(tt[jj + 3]))
         if tt[1] == 'unit':
            for ii in range(3):
               lattice_vectors.append([])
               ss = hh.readline()
               tt = ss.split()
               for jj in range(3):
                  lattice_vectors[ii].append(float(tt[jj + 3]))
         if tt[1] == 'ni':
            for ii in range(3):
               grid_size.append(int(tt[ii + 3]))
         if tt[1] == 'real':
            break
      for kk in range(grid_size[2]):
         wavefunction.append([])
         for jj in range(grid_size[1]):
            wavefunction[kk].append([])
            for ii in range(grid_size[0]):
               while True:
                  ss = hh.readline()
                  tt = ss.split()
                  if len(tt) != 0:
                     break
               if len(tt) > 1:
                  wavefunction[kk][jj].append([])
                  for ll in range(2):
                     wavefunction[kk][jj][ii].append(float(tt[ll]))
# temporary hack for a3dr format
# in case of real wavefunctions, paratec does not write imaginary column
# to a3dr file, while plotxct does write zero imaginary column
               else:
                  wavefunction[kk][jj].append([float(tt[0]), 0.0])
# end of hack
      hh.close()
      for ii in range(3):
         if grid_size[ii] < 1:
            ierr = 1
   return ierr, exciton_state, exciton_energy, hole_position, lattice_vectors, grid_size, wavefunction

def cub_write(file, state, energy, na, aspecies, acharge, aposition, grid, origin, step, vol):
   ierr = 0
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      if state == 0 and abs(energy) < 1.0e-6:
         ss = 'This file was produced with PARATEC\n'
         ss += '\n'
      else:
         ss = 'exciton state = %i\n' % state
         ss += 'exciton energy = %.4f eV\n' % energy
      hh.write(ss)
      ss = ' %4i' % na
      for jj in range(3):
         ss += ' %11.6f' % origin[jj]
      ss += '\n'
      hh.write(ss)
      for ii in range(3):
         ss = ' %4i' % grid[ii]
         for jj in range(3):
            ss += ' %11.6f' % step[ii][jj]
         ss += '\n'
         hh.write(ss)
      for ii in range(na):
         ss = ' %4i' % aspecies[ii]
         ss += ' %11.6f' % acharge[ii]
         for jj in range(3):
            ss += ' %11.6f' % aposition[ii][jj]
         ss += '\n'
         hh.write(ss)
      ss = ''
      for ii in range(grid[0]):
         for jj in range(grid[1]):
            for kk in range(grid[2]):
               ss += ' %12.5E' % vol[ii][jj][kk]
               if (kk + 1) % 6 == 0 or kk + 1 == grid[2]:
                  ss += '\n'
                  hh.write(ss)
                  ss = ''
      hh.close()
   return ierr

def xsf_write(file, state, energy, na, aspecies, aposition, grid, origin, latvec, vol):
   ierr = 0
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      hh.write('CRYSTAL\n')
      hh.write('PRIMVEC\n')
      for ii in range(3):
         ss = ''
         for jj in range(3):
            ss += ' %11.6f' % latvec[ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.write('CONVVEC\n')
      for ii in range(3):
         ss = ''
         for jj in range(3):
            ss += ' %11.6f' % latvec[ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.write('PRIMCOORD\n')
      ss = ' %4i %i\n' % (na, 1)
      hh.write(ss)
      for ii in range(na):
         ss = ' %4i' % aspecies[ii]
         for jj in range(3):
            ss += ' %11.6f' % aposition[ii][jj]
         for jj in range(3):
            ss += ' %11.6f' % 0.0
         ss += '\n'
         hh.write(ss)
      hh.write('BEGIN_BLOCK_DATAGRID_3D\n')
      tt = 'exciton_state_%i_energy_%.4f_eV' % (state, energy)
      ss = '%s\n' % tt
      hh.write(ss)
      ss = 'BEGIN_DATAGRID_3D_%s\n' % tt
      hh.write(ss)
      ss = ''
      for ii in range(3):
         jj = grid[ii] + 1
         ss += ' %11i' % jj
      ss += '\n'
      hh.write(ss)
      ss = ''
      for ii in range(3):
         ss += ' %11.6f' % origin[ii]
      ss += '\n'
      hh.write(ss)
      for ii in range(3):
         ss = ''
         for jj in range(3):
            ss += ' %11.6f' % latvec[ii][jj]
         ss += '\n'
         hh.write(ss)
      ss = ''
      for kk in range(grid[2] + 1):
         if kk < grid[2]:
            pk = kk
         else:
            pk = 0
         for jj in range(grid[1] + 1):
            if jj < grid[1]:
               pj = jj
            else:
               pj = 0
            for ii in range(grid[0] + 1):
               if ii < grid[0]:
                  pi = ii
               else:
                  pi = 0
               ss += ' %12.5E' % vol[pi][pj][pk]
               if (ii + 1) % 6 == 0 or ii == grid[0]:
                  ss += '\n'
                  hh.write(ss)
                  ss = ''
      hh.write('END_DATAGRID_3D\n')
      hh.write('END_BLOCK_DATAGRID_3D\n')
      hh.close()
   return ierr

def cplx_input(cplx):
   ierr = 0
   index = -1
   cplx = cplx.lower()
   if cplx == 'r' or cplx == 're' or cplx == 'real':
      index = 0
   elif cplx == 'i' or cplx == 'im' or cplx == 'imag' or cplx == 'imaginary':
      index = 1
   elif cplx == 'a' or cplx == 'abs' or cplx == 'absolute':
      index = 2
   elif cplx == 'a2' or cplx == 'abs2' or cplx == 'absolute2':
      index = 3
   else:
      ierr = 1
   return ierr, index

if __name__ == "__main__":
   import sys
   import math
   sys.exit(main())

