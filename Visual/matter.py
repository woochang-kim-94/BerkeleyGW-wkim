#!/usr/bin/env python

#-------------------------------------------------------------------------------
#
#   matter.py
#   library of operations on dictionary mat
#   written by Georgy Samsonidze (October 2008)
#
#-------------------------------------------------------------------------------

# FIXME: module, classes, numpy, string formatting, 2to3

import math
from common import eps9, inf9, bohr, format_index2mat, format_espresso2mat, format_mat2espresso, format_siesta2mat, factor_siesta2mat, format_mat2siesta, format_tbpw2mat, format_mat2tbpw, atomsize, bondradius, bondtolerance, minimumbondingdistance, strictvalence, vertex, edge, face, periodic_table, index_by_number, index_by_symbol
from matrix import vector_dot_matrix, matrix_dot_vector, matrix_dot_matrix, determinant_matrix, invert_matrix

def main(argv = None):
   if argv is None:
      argv = sys.argv
   sys.exit("\n   Module matter.py contains the following names:\n" +
            "   mat_create     - create dictionary mat\n" +
            "   mat_check      - check dictionary mat\n" +
            "   mat_species    - get list of species from dictionary mat\n" +
            "   mat_format     - change format of dictionary mat\n" +
            "   mat_lattice    - set lattive parameters for dictionary mat\n" +
            "   mat_remove     - remove atom from dictionary mat\n" +
            "   mat_merge      - merge two dictionaries mat\n" +
            "   mat_replicate  - replicate dictionary mat into supercell\n" +
            "   mat_translate  - translate dictionary mat\n" +
            "   mat_rotate     - rotate dictionary mat\n" +
            "   mat_read       - read dictionary mat from file mat\n" +
            "   mat_write      - write dictionary mat to file mat\n" +
            "   paratec_read   - read dictionary mat from PARATEC input file\n" +
            "   paratec_write  - write dictionary mat to PARATEC input file\n" +
            "   vasp_read      - read dictionary mat from VASP POSCAR file\n" +
            "   vasp_write     - write dictionary mat to VASP POSCAR file\n" +
            "   espresso_read  - read dictionary mat from QE input file\n" +
            "   espresso_write - write dictionary mat to QE input file\n" +
            "   siesta_read    - read dictionary mat from SIESTA input file\n" +
            "   siesta_write   - write dictionary mat to SIESTA input file\n" +
            "   tbpw_read      - read dictionary mat from TBPW input file\n" +
            "   tbpw_write     - write dictionary mat to TBPW input file\n" +
            "   xyz_read       - read dictionary mat from XYZ file\n" +
            "   xyz_write      - write dictionary mat to XYZ file\n" +
            "   xsf_read       - read dictionary mat from XCrysDen XSF file\n" +
            "   xsf_write      - write dictionary mat to XCrysDen XSF file\n" +
            "   wien_write     - write dictionary mat to WIEN2k XYZ file\n" +
            "   povray_write   - write dictionary mat to POV-Ray file\n")

def mat_create():
   mat = {}
   mat['nm'] = ''
   mat['ns'] = 0
   mat['na'] = 0
   mat['fc'] = format_index2mat[0]
   mat['lc'] = 0.0
   mat['fo'] = format_index2mat[0]
   mat['lo'] = [0.0, 0.0, 0.0]
   mat['fv'] = format_index2mat[0]
   mat['lv'] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
   mat['fp'] = format_index2mat[0]
   mat['as'] = []
   mat['ap'] = []
   return mat

def mat_check(mat):
   ierr = 0
   if 'nm' in mat and 'ns' in mat and 'na' in mat and 'fc' in mat and 'lc' in mat and 'fo' in mat and 'lo' in mat and 'fv' in mat and 'lv' in mat and 'fp' in mat and 'as' in mat and 'ap' in mat:
      if type(mat['nm']) is not str:
         ierr = 1
      if type(mat['ns']) is not int:
         ierr = 1
      if type(mat['na']) is not int:
         ierr = 1
      if mat['fc'] != format_index2mat[0] and mat['fc'] != format_index2mat[1]:
         ierr = 1
      if type(mat['lc']) is not float:
         ierr = 1
      if mat['fo'] != format_index2mat[0] and mat['fo'] != format_index2mat[1] and mat['fo'] != format_index2mat[2] and mat['fo'] != format_index2mat[3]:
         ierr = 1
      if len(mat['lo']) != 3:
         ierr = 1
      for jj in range(3):
         if type(mat['lo'][jj]) is not float:
            ierr = 1
      if mat['fv'] != format_index2mat[0] and mat['fv'] != format_index2mat[1] and mat['fv'] != format_index2mat[2]:
         ierr = 1
      if len(mat['lv']) != 3:
         ierr = 1
      for ii in range(3):
         if len(mat['lv'][ii]) != 3:
            ierr = 1
      for ii in range(3):
         for jj in range(3):
            if type(mat['lv'][ii][jj]) is not float:
               ierr = 1
      if mat['fp'] != format_index2mat[0] and mat['fp'] != format_index2mat[1] and mat['fp'] != format_index2mat[2] and mat['fp'] != format_index2mat[3]:
         ierr = 1
      if len(mat['as']) != mat['na']:
         ierr = 1
      for ii in range(mat['na']):
         if type(mat['as'][ii]) is not int:
            ierr = 1
      if len(mat_species(mat['na'], mat['as'])) != mat['ns']:
         ierr = 1
      if len(mat['ap']) != mat['na']:
         ierr = 1
      for ii in range(mat['na']):
         if len(mat['ap'][ii]) != 3:
            ierr = 1
      for ii in range(mat['na']):
         for jj in range(3):
            if type(mat['ap'][ii][jj]) is not float:
               ierr = 1
   else:
      ierr = 1
   return ierr

def mat_species(mat_na, mat_as):
   cs = []
   for ii in range(mat_na):
      if mat_as[ii] not in cs:
         cs.append(mat_as[ii])
   return cs

def mat_format(mat, fc, fo, fv, fp):
   mat['fc'], mat['lc'] = mat_fc_lc(mat['fc'], mat['lc'], fc)
   mat['fv'], mat['lv'] = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], fv)
   mat['fo'], mat['lo'] = mat_fo_lo(mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fo'], mat['lo'], fo)
   mat['fp'], mat['ap'] = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], fp)
   return mat

def mat_lattice(mat, axes, distance, format):
   if format != format_index2mat[3]:
      if mat['fo'] == format_index2mat[3]:
         mat['fo'], mat['lo'] = mat_fo_lo(mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fo'], mat['lo'], mat['fv'])
      if mat['fp'] == format_index2mat[3]:
         mat['fp'], mat['ap'] = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], mat['fv'])
      fr, rmin, rmax = mat_range(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
      if fr == format:
         fo = fr
         lo = []
         for jj in range(3):
            lo.append(rmin[jj] - distance / 2)
         fo, lo = mat_fo_lo(mat['fc'], mat['lc'], mat['fv'], mat['lv'], fo, lo, mat['fo'])
         mat['fo'] = fo
         mat['lo'] = lo
         fv = fr
         lv = []
         for ii in range(3):
            lv.append([])
            for jj in range(3):
               if jj == ii:
                  lv[ii].append(rmax[jj] - rmin[jj] + distance)
               else:
                  lv[ii].append(0.0)
         fv, lv = mat_fv_lv(mat['fc'], mat['lc'], fv, lv, mat['fv'])
         if axes[0] != 0 and axes[1] != 0 and axes[2] != 0:
            mat['fv'] = fv
            mat['lv'] = lv
         elif fv == mat['fv']:
            for ii in range(3):
               if axes[ii] != 0:
                  for jj in range(3):
                     mat['lv'][ii][jj] = lv[ii][jj]
   return mat

def mat_remove(mat, index):
   if index >= 0 and index < mat['na']:
      mat['na'] -= 1
      del mat['as'][index]
      del mat['ap'][index]
      mat['ns'] = len(mat_species(mat['na'], mat['as']))
   return mat

def mat_merge(matc, mats):
   mat = {}
   if matc['fc'] == format_index2mat[0] and mats['fc'] == format_index2mat[0] and matc['fo'] == format_index2mat[0] and mats['fo'] == format_index2mat[0] and matc['fv'] == format_index2mat[0] and mats['fv'] == format_index2mat[0] and matc['fp'] == format_index2mat[0] and mats['fp'] == format_index2mat[0]:
      format = format_index2mat[0]
   if matc['fc'] == format_index2mat[1] and mats['fc'] == format_index2mat[1] and matc['fo'] == format_index2mat[1] and mats['fo'] == format_index2mat[1] and matc['fv'] == format_index2mat[1] and mats['fv'] == format_index2mat[1] and matc['fp'] == format_index2mat[1] and mats['fp'] == format_index2mat[1]:
      format = format_index2mat[1]
   else:
      format = format_index2mat[0]
   fcc, lcc = mat_fc_lc(matc['fc'], matc['lc'], format)
   fcs, lcs = mat_fc_lc(mats['fc'], mats['lc'], format)
   fvc, lvc = mat_fv_lv(matc['fc'], matc['lc'], matc['fv'], matc['lv'], format)
   fvs, lvs = mat_fv_lv(mats['fc'], mats['lc'], mats['fv'], mats['lv'], format)
   foc, loc = mat_fo_lo(matc['fc'], matc['lc'], matc['fv'], matc['lv'], matc['fo'], matc['lo'], format)
   fos, los = mat_fo_lo(mats['fc'], mats['lc'], mats['fv'], mats['lv'], mats['fo'], mats['lo'], format)
   fpc, apc = mat_fp_ap(matc['na'], matc['fc'], matc['lc'], matc['fv'], matc['lv'], matc['fp'], matc['ap'], format)
   fps, aps = mat_fp_ap(mats['na'], mats['fc'], mats['lc'], mats['fv'], mats['lv'], mats['fp'], mats['ap'], format)
   if fcc == format and fcs == format and fvc == format and fvs == format and foc == format and fos == format and fpc == format and fps == format:
      mat['nm'] = '%s %s' % (matc['nm'], mats['nm'])
      mat['na'] = matc['na'] + mats['na']
      mat['fc'] = matc['fc']
      mat['lc'] = matc['lc']
      mat['fo'] = format
      mat['lo'] = loc
      mat['fv'] = format
      mat['lv'] = []
      for ii in range(3):
         mat['lv'].append([])
         for jj in range(3):
            vc = lvc[ii][jj]
            vs = lvs[ii][jj] + los[jj] - loc[jj]
            if abs(vc) > abs(vs):
               lv = vc
            else:
               lv = vs
            mat['lv'][ii].append(lv)
      mat['fp'] = format
      mat['as'] = []
      mat['ap'] = []
      for ii in range(matc['na']):
         mat['as'].append(matc['as'][ii])
         mat['ap'].append(apc[ii])
      for ii in range(mats['na']):
         mat['as'].append(mats['as'][ii])
         mat['ap'].append(aps[ii])
      mat['ns'] = len(mat_species(mat['na'], mat['as']))
   fcc = matc['fc']
   fvc = matc['fv']
   foc = matc['fo']
   fpc = matc['fp']
   mat['fc'], mat['lc'] = mat_fc_lc(mat['fc'], mat['lc'], fcc)
   mat['fv'], mat['lv'] = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], fvc)
   mat['fo'], mat['lo'] = mat_fo_lo(mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fo'], mat['lo'], foc)
   mat['fp'], mat['ap'] = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], fpc)
   return mat

def mat_replicate(mat, so, sv, format):
   smin = []
   smax = []
   si = []
   imax = 1000
   if format == format_index2mat[3]:
      format = mat['fv']
      so = vector_dot_matrix(so, mat['lv'])
      sv = matrix_dot_matrix(sv, mat['lv'])
   fv, lv = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], format)
   fo, lo = mat_fo_lo(mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fo'], mat['lo'], format)
   fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fo == format and fv == format and fp == format:
      sn, sw = normal_make(sv)
      ln, lw = normal_make(lv)
      scenter = []
      for jj in range(3):
         scenter.append(so[jj])
      for ii in range(3):
         for jj in range(3):
            scenter[jj] += sv[ii][jj] / 2
      scorner = []
      for ii in range(8):
         scorner.append([])
         for jj in range(3):
            scorner[ii].append(scenter[jj])
            for kk in range(3):
               scorner[ii][jj] += sv[kk][jj] * (vertex[ii][kk] - 0.5)
      for ii in range(3):
         smin.append(imax)
         smax.append(-imax)
         for jj in range(8):
            pp = 0.0
            for kk in range(3):
               pp += scorner[jj][kk] * ln[ii][kk]
            pp /= lw[ii]
            kk = int(pp)
            if pp < 0.0:
               kk -= 1
            else:
               kk += 1
            if kk < smin[ii]:
               smin[ii] = kk
            if kk > smax[ii]:
               smax[ii] = kk
      amin = []
      amax = []
      for ii in range(3):
         amin.append(imax)
         amax.append(-imax)
         for jj in range(mat['na']):
            pp = 0.0
            for kk in range(3):
               pp += ap[jj][kk] * ln[ii][kk]
            pp /= lw[ii]
            kk = int(pp)
            if pp < 0.0:
               kk -= 1
            else:
               kk += 1
            if kk < amin[ii]:
               amin[ii] = kk
            if kk > amax[ii]:
               amax[ii] = kk
      jj = 0
      for ii in range(3):
         if abs(amin[ii]) > jj:
            jj = abs(amin[ii])
         if abs(amax[ii]) > jj:
            jj = abs(amax[ii])
      for ii in range(3):
         smin[ii] -= jj
         smax[ii] += jj
      pscenter = []
      for ii in range(3):
         pscenter.append(0.0)
         for jj in range(3):
            pscenter[ii] += scenter[jj] * sn[ii][jj]
      psmin = []
      psmax = []
      for ii in range(3):
         psmin.append(pscenter[ii] - 0.5 * sw[ii] - eps9)
         psmax.append(pscenter[ii] + 0.5 * sw[ii] - eps9)
      sna = 0
      sas = []
      sap = []
      for ii in range(smin[0], smax[0] + 1):
         for jj in range(smin[1], smax[1] + 1):
            for kk in range(smin[2], smax[2] + 1):
               for ll in range(mat['na']):
                  sa = []
                  for mm in range(3):
                     sa.append(lv[0][mm] * ii + lv[1][mm] * jj + lv[2][mm] * kk + ap[ll][mm])
                  pa = []
                  for mm in range(3):
                     pa.append(0.0)
                     for nn in range(3):
                        pa[mm] += sa[nn] * sn[mm][nn]
                  flag = 1
                  for mm in range(3):
                     if pa[mm] < psmin[mm] or pa[mm] > psmax[mm]:
                        flag = 0
                  if flag != 0:
                     sna += 1
                     sas.append(mat['as'][ll])
                     sap.append(sa)
                     si.append([ll, ii, jj, kk])
      mat['na'] = sna
      mat['ns'] = len(mat_species(sna, sas))
      mat['fo'] = format
      mat['lo'] = so
      mat['fv'] = format
      mat['lv'] = sv
      mat['fp'] = format
      mat['as'] = sas
      mat['ap'] = sap
   return smin, smax, si, mat

def mat_translate(mat, translation, format):
   fo, lo = mat_fo_lo(mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fo'], mat['lo'], format)
   fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fo == format and fp == format:
      for jj in range(3):
         lo[jj] += translation[jj]
      for ii in range(mat['na']):
         for jj in range(3):
            ap[ii][jj] += translation[jj]
      mat['fo'], mat['lo'] = mat_fo_lo(mat['fc'], mat['lc'], mat['fv'], mat['lv'], fo, lo, mat['fo'])
      mat['fp'], mat['ap'] = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], fp, ap, mat['fp'])
   return mat

def mat_rotate(mat, rotation):
   if mat['fo'] == format_index2mat[3] or mat['fp'] == format_index2mat[3]:
      format = format_index2mat[0]
   else:
      format = mat['fp']
   fv, lv = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], format)
   fo, lo = mat_fo_lo(mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fo'], mat['lo'], format)
   fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fo == format and fv == format and fp == format:
      lo = matrix_dot_vector(rotation, lo)
      for ii in range(3):
         lv[ii] = matrix_dot_vector(rotation, lv[ii])
      for ii in range(mat['na']):
         ap[ii] = matrix_dot_vector(rotation, ap[ii])
      mat['fv'], mat['lv'] = mat_fv_lv(mat['fc'], mat['lc'], fv, lv, mat['fv'])
      mat['fo'], mat['lo'] = mat_fo_lo(mat['fc'], mat['lc'], mat['fv'], mat['lv'], fo, lo, mat['fo'])
      mat['fp'], mat['ap'] = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], fp, ap, mat['fp'])
   return mat

def mat_fc_lc(mat_fc, mat_lc, fc):
   if fc == format_index2mat[0] and mat_fc == format_index2mat[1]:
      lc = mat_lc / bohr
   elif fc == format_index2mat[1] and mat_fc == format_index2mat[0]:
      lc = mat_lc * bohr
   else:
      fc = mat_fc
      lc = mat_lc
   return fc, lc

def mat_fv_lv(mat_fc, mat_lc, mat_fv, mat_lv, fv):
   flag = 1
   if fv == format_index2mat[0] and mat_fv == format_index2mat[1]:
      factor = 1.0 / bohr
   elif fv == format_index2mat[1] and mat_fv == format_index2mat[0]:
      factor = bohr
   elif (fv == format_index2mat[0] or fv == format_index2mat[1]) and mat_fv == format_index2mat[2]:
      fc, lc = mat_fc_lc(mat_fc, mat_lc, fv)
      factor = lc
   elif fv == format_index2mat[2] and (mat_fv == format_index2mat[0] or mat_fv == format_index2mat[1]):
      fc, lc = mat_fc_lc(mat_fc, mat_lc, mat_fv)
      if abs(lc) > eps9:
         factor = 1.0 / lc
      else:
         flag = 0
   else:
      flag = 0
   if flag == 0:
      fv = mat_fv
      factor = 1.0
   lv = []
   for ii in range(3):
      lv.append([])
      for jj in range(3):
         lv[ii].append(mat_lv[ii][jj] * factor)
   return fv, lv

def mat_fo_lo(mat_fc, mat_lc, mat_fv, mat_lv, mat_fo, mat_lo, fo):
   flag = 1
   if fo == format_index2mat[0] and mat_fo == format_index2mat[1]:
      factor = 1.0 / bohr
   elif fo == format_index2mat[1] and mat_fo == format_index2mat[0]:
      factor = bohr
   elif (fo == format_index2mat[0] or fo == format_index2mat[1]) and mat_fo == format_index2mat[2]:
      fc, lc = mat_fc_lc(mat_fc, mat_lc, fo)
      factor = lc
   elif fo == format_index2mat[2] and (mat_fo == format_index2mat[0] or mat_fo == format_index2mat[1]):
      fc, lc = mat_fc_lc(mat_fc, mat_lc, mat_fo)
      if abs(lc) > eps9:
         factor = 1.0 / lc
      else:
         flag = 0
   elif (fo == format_index2mat[0] or fo == format_index2mat[1] or fo == format_index2mat[2]) and mat_fo == format_index2mat[3]:
      fv, lv = mat_fv_lv(mat_fc, mat_lc, mat_fv, mat_lv, fo)
      if fv == fo:
         matrix = lv
         flag = 2
      else:
         flag = 0
   elif fo == format_index2mat[3] and (mat_fo == format_index2mat[0] or mat_fo == format_index2mat[1] or mat_fo == format_index2mat[2]):
      fv, lv = mat_fv_lv(mat_fc, mat_lc, mat_fv, mat_lv, mat_fo)
      if fv == mat_fo:
         ierr, matrix = invert_matrix(lv)
         if ierr == 0:
            flag = 2
         else:
            flag = 0
      else:
         flag = 0
   else:
      flag = 0
   if flag == 0:
      fo = mat_fo
      factor = 1.0
   if flag != 2:
      lo = []
      for jj in range(3):
         lo.append(mat_lo[jj] * factor)
   else:
      lo = vector_dot_matrix(mat_lo, matrix)
   return fo, lo

def mat_fp_ap(mat_na, mat_fc, mat_lc, mat_fv, mat_lv, mat_fp, mat_ap, fp):
   flag = 1
   if fp == format_index2mat[0] and mat_fp == format_index2mat[1]:
      factor = 1.0 / bohr
   elif fp == format_index2mat[1] and mat_fp == format_index2mat[0]:
      factor = bohr
   elif (fp == format_index2mat[0] or fp == format_index2mat[1]) and mat_fp == format_index2mat[2]:
      fc, lc = mat_fc_lc(mat_fc, mat_lc, fp)
      factor = lc
   elif fp == format_index2mat[2] and (mat_fp == format_index2mat[0] or mat_fp == format_index2mat[1]):
      fc, lc = mat_fc_lc(mat_fc, mat_lc, mat_fp)
      if abs(lc) > eps9:
         factor = 1.0 / lc
      else:
         flag = 0
   elif (fp == format_index2mat[0] or fp == format_index2mat[1] or fp == format_index2mat[2]) and mat_fp == format_index2mat[3]:
      fv, lv = mat_fv_lv(mat_fc, mat_lc, mat_fv, mat_lv, fp)
      if fv == fp:
         matrix = lv
         flag = 2
      else:
         flag = 0
   elif fp == format_index2mat[3] and (mat_fp == format_index2mat[0] or mat_fp == format_index2mat[1] or mat_fp == format_index2mat[2]):
      fv, lv = mat_fv_lv(mat_fc, mat_lc, mat_fv, mat_lv, mat_fp)
      if fv == mat_fp:
         ierr, matrix = invert_matrix(lv)
         if ierr == 0:
            flag = 2
         else:
            flag = 0
      else:
         flag = 0
   else:
      flag = 0
   if flag == 0:
      fp = mat_fp
      factor = 1.0
   ap = []
   if flag != 2:
      for ii in range(mat_na):
         ap.append([])
         for jj in range(3):
            ap[ii].append(mat_ap[ii][jj] * factor)
   else:
      for ii in range(mat_na):
         ap.append(vector_dot_matrix(mat_ap[ii], matrix))
   return fp, ap

def mat_range(mat_na, mat_fc, mat_lc, mat_fv, mat_lv, mat_fp, mat_ap, fr):
   fp, ap = mat_fp_ap(mat_na, mat_fc, mat_lc, mat_fv, mat_lv, mat_fp, mat_ap, fr)
   if fp != fr:
      fr = fp
   if mat_na > 0:
      rmin = [inf9, inf9, inf9]
      rmax = [-inf9, -inf9, -inf9]
      for jj in range(3):
         for ii in range(mat_na):
            if ap[ii][jj] < rmin[jj]:
               rmin[jj] = ap[ii][jj]
            if ap[ii][jj] > rmax[jj]:
               rmax[jj] = ap[ii][jj]
   else:
      rmin = [0.0, 0.0, 0.0]
      rmax = [0.0, 0.0, 0.0]
   return fr, rmin, rmax

def normal_make(vector):
   normal = []
   weight = []
   normal.append([])
   normal[0].append(vector[1][1] * vector[2][2] - vector[1][2] * vector[2][1])
   normal[0].append(vector[1][2] * vector[2][0] - vector[1][0] * vector[2][2])
   normal[0].append(vector[1][0] * vector[2][1] - vector[1][1] * vector[2][0])
   normal.append([])
   normal[1].append(vector[2][1] * vector[0][2] - vector[2][2] * vector[0][1])
   normal[1].append(vector[2][2] * vector[0][0] - vector[2][0] * vector[0][2])
   normal[1].append(vector[2][0] * vector[0][1] - vector[2][1] * vector[0][0])
   normal.append([])
   normal[2].append(vector[0][1] * vector[1][2] - vector[0][2] * vector[1][1])
   normal[2].append(vector[0][2] * vector[1][0] - vector[0][0] * vector[1][2])
   normal[2].append(vector[0][0] * vector[1][1] - vector[0][1] * vector[1][0])
   for ii in range(3):
      weight.append(0.0)
      for jj in range(3):
         weight[ii] += math.pow(normal[ii][jj], 2)
      weight[ii] = math.sqrt(weight[ii])
   for ii in range(3):
      if abs(weight[ii]) < eps9:
         weight[ii] = 1.0
   for ii in range(3):
      for jj in range(3):
         normal[ii][jj] /= weight[ii]
   for ii in range(3):
      weight[ii] = 0.0
      for jj in range(3):
         weight[ii] += vector[ii][jj] * normal[ii][jj]
   for ii in range(3):
      for jj in range(3):
         normal[ii][jj] *= weight[ii]
   for ii in range(3):
      weight[ii] = 0.0
      for jj in range(3):
         weight[ii] += math.pow(normal[ii][jj], 2)
   for ii in range(3):
      if abs(weight[ii]) < eps9:
         weight[ii] = 1.0
   return normal, weight

def mat_read(file):
   ierr = 0
   mat = mat_create()
   mat['fc'] = format_index2mat[0]
   mat['lc'] = 1.0
   mat['fo'] = format_index2mat[0]
   mat['fv'] = format_index2mat[0]
   mat['fp'] = format_index2mat[0]
   try:
      hh = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      rr = hh.readlines()
      hh.close()
      ii = 1
      ic = 0
      io = 0
      iv = 0
      ip = 0
      for ss in rr:
         tt = ss.split()
         if len(tt) > 0:
            if tt[0].lower() == 'name_of_matter':
               for jj in range(1, len(tt)):
                  if jj > 1:
                     mat['nm'] += ' '
                  mat['nm'] += tt[jj]
            elif tt[0].lower() == 'number_of_species':
               mat['ns'] = int(tt[1])
            elif tt[0].lower() == 'number_of_atoms':
               mat['na'] = int(tt[1])
            elif tt[0].lower() == 'lattice_constant':
               mat['fc'] = tt[1].lower()
               ic = ii
            elif tt[0].lower() == 'lattice_origin':
               mat['fo'] = tt[1].lower()
               io = ii
            elif tt[0].lower() == 'lattice_vectors':
               mat['fv'] = tt[1].lower()
               iv = ii
            elif tt[0].lower() == 'atomic_species_and_atomic_positions':
               mat['fp'] = tt[1].lower()
               ip = ii
         ii += 1
      if ic > 0:
         ss = rr[ic]
         tt = ss.split()
         mat['lc'] = float(tt[0])
      if io > 0:
         ss = rr[io]
         tt = ss.split()
         for jj in range(3):
            mat['lo'][jj] = float(tt[jj])
      if iv > 0:
         for ii in range(3):
            ss = rr[iv + ii]
            tt = ss.split()
            for jj in range(3):
               mat['lv'][ii][jj] = float(tt[jj])
      if ip > 0:
         for ii in range(mat['na']):
            ss = rr[ip + ii]
            tt = ss.split()
            mat['as'].append(int(tt[0]))
            mat['ap'].append([])
            for jj in range(3):
               mat['ap'][ii].append(float(tt[jj + 1]))
   return ierr, mat

def mat_write(file, mat):
   ierr = 0
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      hh.write('\n')
      ss = 'name_of_matter %s\n' % mat['nm']
      hh.write(ss)
      ss = 'number_of_species %i\n' % mat['ns']
      hh.write(ss)
      ss = 'number_of_atoms %i\n' % mat['na']
      hh.write(ss)
      hh.write('\n')
      ss = 'lattice_constant %s\n' % mat['fc']
      hh.write(ss)
      ss = '     %13.9f\n' % mat['lc']
      hh.write(ss)
      hh.write('\n')
      ss = 'lattice_origin %s\n' % mat['fo']
      hh.write(ss)
      ss = '    '
      for jj in range(3):
         ss += ' %13.9f' % mat['lo'][jj]
      ss += '\n'
      hh.write(ss)
      hh.write('\n')
      ss = 'lattice_vectors %s\n' % mat['fv']
      hh.write(ss)
      for ii in range(3):
         ss = '    '
         for jj in range(3):
            ss += ' %13.9f' % mat['lv'][ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.write('\n')
      ss = 'atomic_species_and_atomic_positions %s\n' % mat['fp']
      hh.write(ss)
      for ii in range(mat['na']):
         ss = ' %3i' % mat['as'][ii]
         for jj in range(3):
            ss += ' %13.9f' % mat['ap'][ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.write('\n')
      hh.close()
   return ierr

def paratec_read(file):
   ierr = 0
   mat = mat_create()
   mat['nm'] = 'paratec'
   mat['fc'] = format_index2mat[0]
   mat['lc'] = 1.0
   mat['fv'] = format_index2mat[0]
   mat['fp'] = format_index2mat[3]
   try:
      hh = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      rr = hh.readlines()
      hh.close()
      ii = 1
      iv = 0
      ip = 0
      flag = 0
      for ss in rr:
         tt = ss.split()
         if len(tt) > 0:
            if tt[0].lower() == 'volume':
               volume = float(tt[1])
               flag = 1
            elif tt[0].lower() == 'coordinates_absolute':
               mat['fp'] = format_index2mat[0]
            elif tt[0].lower() == 'begin':
               if tt[1].lower() == 'latticevecs':
                  iv = ii
               elif tt[1].lower() == 'coordinates':
                  ip = ii
         ii += 1
      if iv > 0:
         for ii in range(3):
            ss = rr[iv + ii]
            tt = ss.split()
            for jj in range(3):
               mat['lv'][ii][jj] = float(tt[jj + 1])
      if ip > 0:
         ii = 0
         na = 0
         ns = 0
         number = 0
         while True:
            ss = rr[ip + ii]         
            tt = ss.split()
            if tt[0].lower() == 'newtype':
               index = index_by_symbol(tt[1])
               number = periodic_table[index]['number']
               ns += 1
            elif tt[0].lower() == 'coord':
               mat['as'].append(number)
               mat['ap'].append([])
               for jj in range(3):
                  mat['ap'][na].append(float(tt[jj + 1]))
               na += 1
            else:
               break
            ii += 1
         mat['na'] = na
         mat['ns'] = ns
      if flag != 0:
         ierr, det = determinant_matrix(mat['lv'])
         if abs(det) > eps9:
            mat['fc'] = format_index2mat[0]
            mat['lc'] = math.pow(volume / det, 1.0 / 3.0)
            mat['fv'] = format_index2mat[2]
         else:
            ierr = 1
   format = mat['fv']
   mat['fo'] = format
   fr, rmin, rmax = mat_range(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fr == format:
      for jj in range(3):
         mat['lo'][jj] = (rmin[jj] + rmax[jj]) / 2
         for ii in range(3):
            mat['lo'][jj] -= mat['lv'][ii][jj] / 2
   return ierr, mat

def paratec_write(file, mat):
   ierr = 0
   nsafety = 5
   eigspacefrac = 0.5
   cs = mat_species(mat['na'], mat['as'])
   nelec = 0
   for ii in range(mat['na']):
      nelec += periodic_table[index_by_number(mat['as'][ii])]['nvelec']
   nbnd = nelec / 2
   if nelec % 2 != 0:
      nbnd += 1
   ratio = nbnd / eigspacefrac
   nbnd = int(ratio)
   if nbnd < ratio - eps9:
      nbnd += 1
   nbndfft = nbnd + nsafety - 1
   if mat['fv'] == format_index2mat[1]:
      fv, lv = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], format_index2mat[0])
      if fv != format_index2mat[0]:
         ierr = 1
   else:
      fv = mat['fv']
      lv = mat['lv']
   if fv == format_index2mat[2]:
      fc, lc = mat_fc_lc(mat['fc'], mat['lc'], format_index2mat[0])
      if fc == format_index2mat[0]:
         ierr, det = determinant_matrix(lv)
         volume = det * math.pow(lc, 3)
      else:
         ierr = 1
   if mat['fp'] == format_index2mat[1] or mat['fp'] == format_index2mat[2]:
      fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format_index2mat[0])
      if fp != format_index2mat[0]:
         ierr = 1
   else:
      fp = mat['fp']
      ap = mat['ap']
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      hh.write('\n')
      hh.write('begin pw_jobs\n')
      hh.write('pw_job scf\n')
      hh.write('end pw_jobs\n')
      hh.write('\n')
      hh.write('begin latticevecs\n')
      for ii in range(3):
         ss = 'coord'
         for jj in range(3):
            ss += ' %13.9f' % lv[ii][jj]
         ss += '\n'
         hh.write(ss)
      if fv == format_index2mat[2]:
         ss = 'volume %.9f\n' % volume
         hh.write(ss)
      hh.write('end latticevecs\n')
      hh.write('\n')
      if fp == format_index2mat[0]:
         hh.write('coordinates_absolute\n')
      hh.write('begin coordinates\n')
      for ii in range(mat['ns']):
         index = index_by_number(cs[ii])
         symbol = periodic_table[index]['symbol']
         ss = 'newtype %s\n' % symbol
         hh.write(ss)
         for jj in range(mat['na']):
            if mat['as'][jj] == cs[ii]:
               ss = 'coord'
               for kk in range(3):
                  ss += ' %13.9f' % ap[jj][kk]
               ss += '\n'
               hh.write(ss)
      hh.write('end coordinates\n')
      hh.write('\n')
      hh.write('pp_format 3\n')
      hh.write('begin pseudopotential\n')
      for ii in range(mat['ns']):
         number = periodic_table[index_by_number(cs[ii])]['number']
         nvelec = periodic_table[index_by_number(cs[ii])]['nvelec']
         period = periodic_table[index_by_number(cs[ii])]['period']
         group = periodic_table[index_by_number(cs[ii])]['group']
         #
         # find occupations of s,p,d,f shells
         #
         if period < 2 or group < 3:
            occupation = [nvelec, 0, 0, 0]
         elif period < 8 and group > 12:
            occupation = [2, nvelec - 2, 0, 0]
         elif period < 8:
            occupation = [2, 0, nvelec - 2, 0]
         else:
            occupation = [2, 0, 0, nvelec - 2]
         #
         # exceptions are hard-coded below
         #
         if number == 57:
            occupation = [2, 0, 1, 0]
         if number == 58:
            occupation = [2, 0, 1, 1]
         if number == 89:
            occupation = [2, 0, 1, 0]
         if number == 90:
            occupation = [2, 0, 2, 0]
         if number == 91:
            occupation = [2, 0, 1, 2]
         if number == 92:
            occupation = [2, 0, 1, 3]
         #
         # set local channel
         #
         if occupation[3] != 0:
            lloc = 3
         elif occupation[2] != 0:
            lloc = 2
         elif occupation[1] != 0:
            lloc = 1
         else:
            lloc = 0
         #
         # write pseudopotential information
         #
         ss = 'pp_data %i %i.0 %i.0 %i.0 %i.0\n' % (tuple([lloc + 1]) + tuple(occupation))
         hh.write(ss)
      hh.write('end pseudopotential\n')
      hh.write('\n')
      hh.write('exchange_corelation ceperley_alder\n')
      hh.write('number_of_spins 1\n')
      hh.write('diagonalization_method Grassmann_metal\n')
      hh.write('energy_cutoff 60.0\n')
      hh.write('submatrix_energy_cutoff 3.0\n')
      hh.write('accuracy_diag 1.0e-10\n')
      hh.write('diagsafety 2.0\n')
      hh.write('number_kpoints 0\n')
      hh.write('k_grid 1 1 1\n')
      hh.write('k_grid_shift 0.0 0.0 0.0\n')
      ss = 'number_bands %i\n' % nbnd
      hh.write(ss)
      ss = 'number_bands_fft %i\n' % nbndfft
      hh.write(ss)
      ss = 'eigspacefrac %.1f\n' % eigspacefrac
      hh.write(ss)
      hh.write('occupy_levels normal\n')
      hh.write('smearing_energy 0.001\n')
      hh.write('\n')
      hh.write('mix_method pulay_kerker\n')
      hh.write('energy_convergence_criterion 1.0e-10\n')
      hh.write('max_iter_scfloop 50\n')
      hh.write('screening_type previous\n')
      hh.write('\n')
      hh.close()
   return ierr

def vasp_read(file):
   ierr = 0
   mat = mat_create()
   mat['fc'] = format_index2mat[1]
   mat['fv'] = format_index2mat[2]
   mat['fp'] = format_index2mat[3]
   try:
      hh = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      rr = hh.readlines()
      hh.close()
      #
      # TEMP ASE FIX BEGIN
      #
      #ss = rr[0]
      #tt = ss.split()
      #mat['nm'] = tt[0]
      #
      mat['nm'] = 'vasp'
      #
      # TEMP ASE FIX END
      #
      ss = rr[1]
      tt = ss.split()
      mat['lc'] = float(tt[0])
      for ii in range(3):
         ss = rr[2 + ii]
         tt = ss.split()
         for jj in range(3):
            mat['lv'][ii][jj] = float(tt[jj])
      if mat['lc'] < 0.0:
         volume = -mat['lc']
         ierr, det = determinant_matrix(mat['lv'])
         if abs(det) > eps9:
            mat['lc'] = math.pow(volume / det, 1.0 / 3.0)
      ss = rr[5]
      species = ss.split()
      ss = rr[6]
      tt = ss.split()
      mat['ns'] = len(tt)
      for ii in range(mat['ns']):
         kk = int(tt[ii])
         index = index_by_symbol(species[ii])
         number = periodic_table[index]['number']
         for jj in range(kk):
            mat['as'].append(number)
         mat['na'] += kk
      if rr[7][0].lower() == 's':
         kk = 1
      else:
         kk = 0
      if rr[7 + kk][0].lower() == 'c' or rr[7 + kk][0].lower() == 'k':
         mat['fp'] = format_index2mat[2]
      for ii in range(mat['na']):
         ss = rr[8 + kk + ii]
         tt = ss.split()
         mat['ap'].append([])
         for jj in range(3):
            mat['ap'][ii].append(float(tt[jj]))
   format = mat['fv']
   mat['fo'] = format
   fr, rmin, rmax = mat_range(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fr == format:
      for jj in range(3):
         mat['lo'][jj] = (rmin[jj] + rmax[jj]) / 2
         for ii in range(3):
            mat['lo'][jj] -= mat['lv'][ii][jj] / 2
   return ierr, mat

def vasp_write(file, mat):
   ierr = 0
   cs = mat_species(mat['na'], mat['as'])
   if mat['fc'] != format_index2mat[1]:
      fc, lc = mat_fc_lc(mat['fc'], mat['lc'], format_index2mat[1])
      if fc != format_index2mat[1]:
         ierr = 1
   else:
      fc = mat['fc']
      lc = mat['lc']
   if mat['fv'] != format_index2mat[2]:
      fv, lv = mat_fv_lv(fc, lc, mat['fv'], mat['lv'], format_index2mat[2])
      if fv != format_index2mat[2]:
         ierr = 1
   else:
      fv = mat['fv']
      lv = mat['lv']
   if mat['fp'] == format_index2mat[0] or mat['fp'] == format_index2mat[1]:
      fp, ap = mat_fp_ap(mat['na'], fc, lc, fv, lv, mat['fp'], mat['ap'], format_index2mat[3])
      if fp != format_index2mat[3]:
         ierr = 1
   else:
      fp = mat['fp']
      ap = mat['ap']
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      #
      # TEMP ASE FIX BEGIN
      #
      # replace the name of the system with the list of atomic species
      # for compatibility with the Atomic Simulation Environment (ASE)
      #
      #ss = mat['nm'] + '\n'
      #
      ss = ''
      for ii in range(mat['ns']):
         index = index_by_number(cs[ii])
         symbol = periodic_table[index]['symbol']
         ss += ' %s' % symbol
      ss += '\n'
      #
      # TEMP ASE FIX END
      #
      hh.write(ss)
      ss = ' %13.9f\n' % lc
      hh.write(ss)
      for ii in range(3):
         ss = ''
         for jj in range(3):
            ss += ' %13.9f' % lv[ii][jj]
         ss += '\n'
         hh.write(ss)
      ss = ''
      for ii in range(mat['ns']):
         index = index_by_number(cs[ii])
         symbol = periodic_table[index]['symbol']
         ss += ' %s' % symbol
      ss += '\n'
      hh.write(ss)
      ss = ''
      for ii in range(mat['ns']):
         number = 0
         for jj in range(mat['na']):
            if mat['as'][jj] == cs[ii]:
               number += 1
         ss += ' %i' % number
      ss += '\n'
      hh.write(ss)
      if fp == format_index2mat[2]:
         ss = 'Cartesian\n'
      else:
         ss = 'Direct\n'
      hh.write(ss)
      for ii in range(mat['ns']):
         for jj in range(mat['na']):
            if mat['as'][jj] == cs[ii]:
               ss = ''
               for kk in range(3):
                  ss += ' %13.9f' % ap[jj][kk]
               #
               # last column is atomic species, inspired by pymatgen
               #
               ss += '  %2s' % periodic_table[index_by_number(mat['as'][jj])]['symbol']
               ss += '\n'
               hh.write(ss)
      hh.close()
   return ierr

def espresso_read(file):
   ierr = 0
   mat = mat_create()
   mat['fc'] = format_index2mat[0]
   mat['lc'] = 1.0
   mat['fv'] = format_index2mat[0]
   mat['fp'] = format_index2mat[0]
   try:
      hh = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      uu = hh.readlines()
      hh.close()
      rr = []
      for vv in uu:
         if len(vv.replace(' ', '').replace('\t', '').replace('\n', '')) > 0:
            ww = vv.split(',')
            for xx in ww:
               rr.append(xx)
      crystallographic_constants = 0
      ii = 1
      iv = 0
      ip = 0
      for ss in rr:
         tt = ss.replace('=', ' = ').replace('\'', '').replace('\"', '').replace('{', '').replace('}', '').replace('(', '').replace(')', '').split()
         if len(tt) > 0:
            if tt[0].lower() == 'prefix':
               for jj in range(2, len(tt)):
                  if jj > 2:
                     mat['nm'] += ' '
                  mat['nm'] += tt[jj]
            elif tt[0].lower() == 'ibrav':
               ibrav = int(tt[2])
            elif tt[0].lower() == 'celldm1':
               crystallographic_constants = 1
               mat['fc'] = format_index2mat[0]
               celldm1 = float(tt[2])
               mat['lc'] = celldm1
               mat['fv'] = format_index2mat[2]
               mat['fp'] = format_index2mat[2]
            elif tt[0].lower() == 'celldm2':
               celldm2 = float(tt[2])
            elif tt[0].lower() == 'celldm3':
               celldm3 = float(tt[2])
            elif tt[0].lower() == 'celldm4':
               celldm4 = float(tt[2])
            elif tt[0].lower() == 'celldm5':
               celldm5 = float(tt[2])
            elif tt[0].lower() == 'celldm6':
               celldm6 = float(tt[2])
            elif tt[0].lower() == 'a':
               crystallographic_constants = 2
               mat['fc'] = format_index2mat[1]
               aa = float(tt[2])
               mat['lc'] = aa
               mat['fv'] = format_index2mat[2]
               mat['fp'] = format_index2mat[2]
            elif tt[0].lower() == 'b':
               # to avoid confusion with boron
               if tt[1] == '=':
                  bb = float(tt[2])
            elif tt[0].lower() == 'c':
               # to avoid confusion with carbon
               if tt[1] == '=':
                  cc = float(tt[2])
            elif tt[0].lower() == 'cos_ab':
               cos_gamma = float(tt[2])
            elif tt[0].lower() == 'cos_ac':
               cos_beta = float(tt[2])
            elif tt[0].lower() == 'cos_bc':
               cos_alpha = float(tt[2])
            elif tt[0].lower() == 'ntyp':
               mat['ns'] = int(tt[2])
            elif tt[0].lower() == 'nat':
               mat['na'] = int(tt[2])
            elif tt[0].lower() == 'cell_parameters':
               if len(tt) > 1:
                  fv = tt[1].lower()
                  if fv == format_mat2espresso[format_index2mat[0]]:
                     mat['fc'] = format_index2mat[0]
                     mat['lc'] = 1.0
                     mat['fv'] = format_index2mat[2]
                  elif fv == format_mat2espresso[format_index2mat[1]]:
                     mat['fc'] = format_index2mat[1]
                     mat['lc'] = 1.0
                     mat['fv'] = format_index2mat[2]
               iv = ii
            elif tt[0].lower() == 'atomic_positions':
               if len(tt) > 1:
                  fp = tt[1].lower()
                  if fp in format_espresso2mat:
                     mat['fp'] = format_espresso2mat[fp]
                  else:
                     mat['fp'] = 'latconst'
               else:
                  mat['fp'] = 'latconst'
               ip = ii
         ii += 1
      if iv > 0:
         for ii in range(3):
            ss = rr[iv + ii]
            tt = ss.split()
            for jj in range(3):
               mat['lv'][ii][jj] = float(tt[jj])
      if ip > 0:
         for ii in range(mat['na']):
            ss = rr[ip + ii]
            tt = ss.split()
            mat['as'].append(periodic_table[index_by_symbol(tt[0])]['number'])
            mat['ap'].append([])
            for jj in range(3):
               mat['ap'][ii].append(float(tt[jj + 1]))
      if ibrav == 1:
         mat['lv'] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
      elif ibrav == 2:
         mat['lv'] = [[-0.5, 0.0, 0.5], [0.0, 0.5, 0.5], [-0.5, 0.5, 0.0]]
      elif ibrav == 3:
         mat['lv'] = [[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5]]
      elif ibrav == 4:
         if crystallographic_constants == 2:
            celldm3 = cc / aa
         mat['lv'] = [[1.0, 0.0, 0.0], [-0.5, math.sqrt(3.0) / 2.0, 0.0], [0.0, 0.0, celldm3]]
      elif ibrav == 5:
         if crystallographic_constants == 1:
            cos_alpha = celldm4
         tx = math.sqrt((1.0 - cos_alpha) / 2.0)
         ty = math.sqrt((1.0 - cos_alpha) / 6.0)
         tz = math.sqrt((1.0 + 2.0 * cos_alpha) / 3.0)
         mat['lv'] = [[tx, -ty, tz], [0.0, 2.0 * ty, tz], [-tx, -ty, tz]]
      elif ibrav == -5:
         if crystallographic_constants == 1:
            cos_alpha = celldm4
         tx = math.sqrt((1.0 - cos_alpha) / 2.0)
         ty = math.sqrt((1.0 - cos_alpha) / 6.0)
         tz = math.sqrt((1.0 + 2.0 * cos_alpha) / 3.0)
         uu = tz - 2.0 * math.sqrt(2.0) * ty
         vv = tz + math.sqrt(2.0) * ty
         mat['lv'] = [[uu / math.sqrt(3.0), vv / math.sqrt(3.0), vv / math.sqrt(3.0)], [vv / math.sqrt(3.0), uu / math.sqrt(3.0), vv / math.sqrt(3.0)], [vv / math.sqrt(3.0), vv / math.sqrt(3.0), uu / math.sqrt(3.0)]]
      elif ibrav == 6:
         if crystallographic_constants == 2:
            celldm3 = cc / aa
         mat['lv'] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, celldm3]]
      elif ibrav == 7:
         if crystallographic_constants == 2:
            celldm3 = cc / aa
         mat['lv'] = [[0.5, -0.5, celldm3 / 2.0], [0.5, 0.5, celldm3 / 2.0], [-0.5, -0.5, celldm3 / 2.0]]
      elif ibrav == 8:
         if crystallographic_constants == 2:
            celldm2 = bb / aa
            celldm3 = cc / aa
         mat['lv'] = [[1.0, 0.0, 0.0], [0.0, celldm2, 0.0], [0.0, 0.0, celldm3]]
      elif ibrav == 9:
         if crystallographic_constants == 2:
            celldm2 = bb / aa
            celldm3 = cc / aa
         mat['lv'] = [[0.5, celldm2 / 2.0, 0.0], [-0.5, celldm2 / 2.0, 0.0], [0.0, 0.0, celldm3]]
      elif ibrav == 10:
         if crystallographic_constants == 2:
            celldm2 = bb / aa
            celldm3 = cc / aa
         mat['lv'] = [[0.5, 0.0, celldm3 / 2.0], [0.5, celldm2 / 2.0, 0.0], [0.0, celldm2 / 2.0, celldm3 / 2.0]]
      elif ibrav == 11:
         if crystallographic_constants == 2:
            celldm2 = bb / aa
            celldm3 = cc / aa
         mat['lv'] = [[0.5, celldm2 / 2.0, celldm3 / 2.0], [-0.5, celldm2 / 2.0, celldm3 / 2.0], [-0.5, -celldm2 / 2.0, celldm3 / 2.0]]
      elif ibrav == 12:
         if crystallographic_constants == 1:
            cos_gamma = celldm4
         else:
            celldm2 = bb / aa
            celldm3 = cc / aa
         sin_gamma = math.sqrt(1.0 - math.pow(cos_gamma, 2))
         mat['lv'] = [[1.0, 0.0, 0.0], [celldm2 * cos_gamma, celldm2 * sin_gamma, 0.0], [0.0, 0.0, celldm3]]
      elif ibrav == -12:
         if crystallographic_constants == 1:
            cos_beta = celldm5
         else:
            celldm2 = bb / aa
            celldm3 = cc / aa
         sin_beta = math.sqrt(1.0 - math.pow(cos_beta, 2))
         mat['lv'] = [[1.0, 0.0, 0.0], [0.0, celldm2, 0.0], [sin_beta, 0.0, celldm3 * cos_beta]]
      elif ibrav == 13:
         if crystallographic_constants == 1:
            cos_gamma = celldm4
         else:
            celldm2 = bb / aa
            celldm3 = cc / aa
         sin_gamma = math.sqrt(1.0 - math.pow(cos_gamma, 2))
         mat['lv'] = [[0.5, 0.0, -celldm3 / 2.0], [celldm2 * cos_gamma, celldm2 * sin_gamma, 0.0], [0.5, 0.0, celldm3 / 2.0]]
      elif ibrav == 14:
         if crystallographic_constants == 1:
            cos_alpha = celldm4
            cos_beta = celldm5
            cos_gamma = celldm6
         else:
            celldm2 = bb / aa
            celldm3 = cc / aa
         sin_gamma = math.sqrt(1.0 - math.pow(cos_gamma, 2))
         mat['lv'] = [[1.0, 0.0, 0.0], [celldm2 * cos_gamma, celldm2 * sin_gamma, 0.0], [celldm3 * cos_beta, celldm3 * (cos_alpha - cos_beta * cos_gamma) / sin_gamma, celldm3 * math.sqrt(1.0 + 2.0 * cos_alpha * cos_beta * cos_gamma - math.pow(cos_alpha, 2) - math.pow(cos_beta, 2) - math.pow(cos_gamma, 2)) / sin_gamma]]
   format = mat['fv']
   mat['fo'] = format
   fr, rmin, rmax = mat_range(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fr == format:
      for jj in range(3):
         mat['lo'][jj] = (rmin[jj] + rmax[jj]) / 2
         for ii in range(3):
            mat['lo'][jj] -= mat['lv'][ii][jj] / 2
   return ierr, mat

def espresso_write(file, mat):
   ierr = 0
   cs = mat_species(mat['na'], mat['as'])
   nelec = 0
   for ii in range(mat['na']):
      nelec += periodic_table[index_by_number(mat['as'][ii])]['nvelec']
   nbnd = nelec / 2
   if nelec % 2 != 0:
      nbnd += 1
   if mat['fv'] == format_index2mat[2] or mat['fp'] == format_index2mat[2]:
      celldm = 1
   else:
      celldm = 0
   if mat['fp'] in format_mat2espresso:
      fp = format_mat2espresso[mat['fp']]
   else:
      ierr = 1
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      hh.write('&CONTROL\n')
      ss = '   prefix = \'%s\'\n' % mat['nm']
      hh.write(ss)
      hh.write('   calculation = \'scf\'\n')
      hh.write('   verbosity = \'high\'\n')
      hh.write('   wf_collect = .false.\n')
      hh.write('   tstress = .true.\n')
      hh.write('   tprnfor = .true.\n')
      hh.write('   outdir = \'.\'\n')
      hh.write('   wfcdir = \'.\'\n')
      hh.write('   pseudo_dir = \'.\'\n')
      hh.write('/\n')
      hh.write('&SYSTEM\n')
      hh.write('   ibrav = 0\n')
      if celldm != 0:
         if mat['fc'] == format_index2mat[0]:
            ss = '   celldm(1) = %.9f\n' % mat['lc']
         else:
            ss = '   a = %.9f\n' % mat['lc']
         hh.write(ss)
      ss = '   nat = %i\n' % mat['na']
      hh.write(ss)
      ss = '   ntyp = %i\n' % mat['ns']
      hh.write(ss)
      ss = '   nbnd = %i\n' % nbnd
      hh.write(ss)
      hh.write('   ecutwfc = 60.0\n')
      hh.write('/\n')
      hh.write('&ELECTRONS\n')
      hh.write('   electron_maxstep = 100\n')
      hh.write('   conv_thr = 1.0d-10\n')
      hh.write('   mixing_mode = \'plain\'\n')
      hh.write('   mixing_beta = 0.7\n')
      hh.write('   mixing_ndim = 8\n')
      hh.write('   diagonalization = \'david\'\n')
      hh.write('   diago_david_ndim = 4\n')
      hh.write('   diago_full_acc = .true.\n')
      hh.write('/\n')
      ss = 'CELL_PARAMETERS'
      if mat['fv'] == format_index2mat[0]:
          ss += ' bohr'
      elif mat['fv'] == format_index2mat[1]:
          ss += ' angstrom'
      else:
          ss += ' alat'
      ss += '\n'
      hh.write(ss)
      for ii in range(3):
         ss = '    '
         for jj in range(3):
            ss += ' %13.9f' % mat['lv'][ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.write('ATOMIC_SPECIES\n')
      for ii in range(mat['ns']):
         index = index_by_number(cs[ii])
         ss = '  %2s %7.3f  %s%s\n' % (periodic_table[index]['symbol'], periodic_table[index]['mass'], periodic_table[index]['symbol'], '.UPF')
         hh.write(ss)
      ss = 'ATOMIC_POSITIONS %s\n' % fp
      hh.write(ss)
      for ii in range(mat['na']):
         ss = '  %2s' % periodic_table[index_by_number(mat['as'][ii])]['symbol']
         for jj in range(3):
            ss += ' %13.9f' % mat['ap'][ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.write('K_POINTS gamma\n')
      hh.close()
   return ierr

def siesta_read(file):
   ierr = 0
   mat = mat_create()
   mat['fc'] = format_index2mat[0]
   mat['lc'] = 0.0
   mat['fv'] = format_index2mat[2]
   mat['fp'] = format_index2mat[0]
   try:
      hh = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      rr = hh.readlines()
      hh.close()
      ii = 1
      iv = 0
      il = 0
      ip = 0
      for ss in rr:
         tt = ss.split()
         if len(tt) > 0:
            tt[0] = tt[0].replace('-', '').replace('_', '').replace('.', '')
            if tt[0].lower() == 'systemlabel':
               for jj in range(1, len(tt)):
                  if jj > 1:
                     mat['nm'] += ' '
                  mat['nm'] += tt[jj]
            elif tt[0].lower() == 'numberofspecies':
               mat['ns'] = int(tt[1])
            elif tt[0].lower() == 'numberofatoms':
               mat['na'] = int(tt[1])
            elif tt[0].lower() == 'latticeconstant':
               mat['lc'] = float(tt[1])
               fc = tt[2].lower()
               if fc in format_siesta2mat:
                  mat['lc'] *= factor_siesta2mat[fc]
                  mat['fc'] = format_siesta2mat[fc]
               else:
                  ierr = 1
            elif tt[0].lower() == 'atomiccoordinatesformat':
               fp = tt[1].lower()
               if fp in format_siesta2mat:
                  mat['fp'] = format_siesta2mat[fp]
               else:
                  ierr = 1
            elif tt[0].lower() == '%block':
               if tt[1].lower() == 'latticevectors':
                  iv = ii
               elif tt[1].lower() == 'chemicalspecieslabel':
                  il = ii
               elif tt[1].lower() == 'atomiccoordinatesandatomicspecies':
                  ip = ii
         ii += 1
      if iv > 0:
         for ii in range(3):
            ss = rr[iv + ii]
            tt = ss.split()
            for jj in range(3):
               mat['lv'][ii][jj] = float(tt[jj])
      if il > 0:
         index = []
         number = []
         for ii in range(mat['ns']):
            ss = rr[il + ii]
            tt = ss.split()
            index.append(int(tt[0]))
            number.append(int(tt[1]))
      if ip > 0:
         for ii in range(mat['na']):
            ss = rr[ip + ii]
            tt = ss.split()
            mat['ap'].append([])
            for jj in range(3):
               mat['ap'][ii].append(float(tt[jj]))
            mat['as'].append(number[index.index(int(tt[3]))])
   format = mat['fv']
   mat['fo'] = format
   fr, rmin, rmax = mat_range(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fr == format:
      for jj in range(3):
         mat['lo'][jj] = (rmin[jj] + rmax[jj]) / 2
         for ii in range(3):
            mat['lo'][jj] -= mat['lv'][ii][jj] / 2
   return ierr, mat

def siesta_write(file, mat):
   ierr = 0
   cs = mat_species(mat['na'], mat['as'])
   if mat['fc'] in format_mat2siesta:
      fc = format_mat2siesta[mat['fc']]
   else:
      ierr = 1
   fv, lv = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], format_index2mat[2])
   if fv != format_index2mat[2]:
      ierr = 1
   if mat['fp'] in format_mat2siesta:
      fp = format_mat2siesta[mat['fp']]
   else:
      ierr = 1
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      hh.write('\n')
      ss = 'SystemName      %s\n' % mat['nm']
      hh.write(ss)
      ss = 'SystemLabel     %s\n' % mat['nm']
      hh.write(ss)
      ss = 'NumberOfSpecies %i\n' % mat['ns']
      hh.write(ss)
      ss = 'NumberOfAtoms   %i\n' % mat['na']
      hh.write(ss)
      hh.write('\n')
      ss = 'LatticeConstant %.9f %s\n' % (mat['lc'], fc)
      hh.write(ss)
      hh.write('\n')
      ss = 'AtomicCoordinatesFormat %s\n' % fp
      hh.write(ss)
      hh.write('\n')
      hh.write('%block LatticeVectors\n')
      for ii in range(3):
         ss = ''
         for jj in range(3):
            ss += ' %13.9f' % lv[ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.write('%endblock LatticeVectors\n')
      hh.write('\n')
      hh.write('%block ChemicalSpeciesLabel\n')
      for ii in range(mat['ns']):
         ss = '%i %i %s\n' % (ii + 1, cs[ii], periodic_table[index_by_number(cs[ii])]['symbol'])
         hh.write(ss)
      hh.write('%endblock ChemicalSpeciesLabel\n')
      hh.write('\n')
      hh.write('%block AtomicCoordinatesAndAtomicSpecies\n')
      for ii in range(mat['na']):
         ss = ''
         for jj in range(3):
            ss += ' %13.9f' % mat['ap'][ii][jj]
         an = cs.index(mat['as'][ii]) + 1
         ss += '  %i' % an
         ss += '\n'
         hh.write(ss)
      hh.write('%endblock AtomicCoordinatesAndAtomicSpecies\n')
      hh.write('\n')
      hh.write('UserBasis T\n')
      hh.write('\n')
      hh.write('SolutionMethod diagon\n')
      hh.write('\n')
      hh.write('MD.TypeOfRun CG\n')
      hh.write('MD.VariableCell F\n')
      hh.write('MD.NumCGsteps 0\n')
      hh.write('\n')
      hh.write('%block kgrid_Monkhorst_Pack\n')
      hh.write('1 0 0 +0.0\n')
      hh.write('0 1 0 +0.0\n')
      hh.write('0 0 1 +0.0\n')
      hh.write('%endblock kgrid_Monkhorst_Pack\n')
      hh.write('\n')
      hh.close()
   return ierr

def tbpw_read(file):
   ierr = 0
   mat = mat_create()
   mat['nm'] = 'tbpw'
   mat['fc'] = format_index2mat[0]
   mat['lc'] = 0.0
   mat['fv'] = format_index2mat[2]
   mat['fp'] = format_index2mat[3]
   try:
      hh = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      rr = hh.readlines()
      hh.close()
      ii = 1
      iv = 0
      il = 0
      ip = 0
      for ss in rr:
         tt = ss.split()
         if len(tt) > 0:
            if tt[0].lower() == 'numberofspecies':
               mat['ns'] = int(tt[1])
            elif tt[0].lower() == 'numberofatoms':
               mat['na'] = int(tt[1])
            elif tt[0].lower() == 'latticeconstant':
               mat['lc'] = float(tt[1])
            elif tt[0].lower() == 'latticeconstantformat':
               fc = tt[1].lower()
               if fc in format_tbpw2mat:
                  mat['fc'] = format_tbpw2mat[fc]
               else:
                  ierr = 1
            elif tt[0].lower() == 'atomiccoordinatesformat':
               fp = tt[1].lower()
               if fp in format_tbpw2mat:
                  mat['fp'] = format_tbpw2mat[fp]
               else:
                  ierr = 1
            if tt[0].lower() == 'latticevectors':
               iv = ii
            elif tt[0].lower() == 'chemicalspecieslabel':
               il = ii
            elif tt[0].lower() == 'atomiccoordinatesandatomicspecies':
               ip = ii
         ii += 1
      if iv > 0:
         for ii in range(3):
            ss = rr[iv + ii]
            tt = ss.split()
            for jj in range(3):
               mat['lv'][ii][jj] = float(tt[jj])
      if il > 0:
         index = []
         number = []
         for ii in range(mat['ns']):
            ss = rr[il + ii]
            tt = ss.split()
            index.append(int(tt[0]))
            number.append(int(tt[1]))
      if ip > 0:
         for ii in range(mat['na']):
            ss = rr[ip + ii]
            tt = ss.split()
            mat['ap'].append([])
            for jj in range(3):
               mat['ap'][ii].append(float(tt[jj]))
            mat['as'].append(number[index.index(int(tt[3]))])
   format = mat['fv']
   mat['fo'] = format
   fr, rmin, rmax = mat_range(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fr == format:
      for jj in range(3):
         mat['lo'][jj] = (rmin[jj] + rmax[jj]) / 2
         for ii in range(3):
            mat['lo'][jj] -= mat['lv'][ii][jj] / 2
   return ierr, mat

def tbpw_write(file, mat):
   ierr = 0
   cs = mat_species(mat['na'], mat['as'])
   nelec = 0
   for ii in range(mat['na']):
      nelec += periodic_table[index_by_number(mat['as'][ii])]['nvelec']
   nbnd = nelec / 2
   if mat['fc'] in format_mat2tbpw:
      fc = format_mat2tbpw[mat['fc']]
   else:
      ierr = 1
   fv, lv = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], format_index2mat[2])
   if fv != format_index2mat[2]:
      ierr = 1
   if mat['fp'] in format_mat2tbpw:
      if mat['fp'] == format_index2mat[0] or mat['fp'] == format_index2mat[1]:
         fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format_index2mat[2])
         if fp != format_index2mat[2]:
            ierr = 1
      else:
         fp = mat['fp']
         ap = mat['ap']
      fp = format_mat2tbpw[fp]
   else:
      ierr = 1
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      hh.write('\n')
      hh.write('DiagonalizationSwitch 0\n')
      hh.write('#AbsoluteTolerance -1.0d0\n')
      hh.write('\n')
      hh.write('#InputEnergiesInEV\n')
      hh.write('EnergyCutoff 11.0\n')
      hh.write('\n')
      hh.write('NumberOfDimensions 3\n')
      ss = 'LatticeConstantFormat %s\n' % fc
      hh.write(ss)
      ss = 'LatticeConstant %.9f\n' % mat['lc']
      hh.write(ss)
      hh.write('LatticeVectors\n')
      for ii in range(3):
         ss = ''
         for jj in range(3):
            ss += ' %13.9f' % lv[ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.write('\n')
      ss = 'NumberOfAtoms %i\n' % mat['na']
      hh.write(ss)
      ss = 'NumberOfSpecies %i\n' % mat['ns']
      hh.write(ss)
      hh.write('ChemicalSpeciesLabel\n')
      for ii in range(mat['ns']):
         ss = '%i %i %s\n' % (ii + 1, cs[ii], periodic_table[index_by_number(cs[ii])]['symbol'])
         hh.write(ss)
      ss = 'AtomicCoordinatesFormat %s\n' % fp
      hh.write(ss)
      hh.write('AtomicCoordinatesAndAtomicSpecies\n')
      for ii in range(mat['na']):
         ss = ''
         for jj in range(3):
            ss += ' %13.9f' % ap[ii][jj]
         an = cs.index(mat['as'][ii]) + 1
         ss += '  %i' % an
         ss += '\n'
         hh.write(ss)
      hh.write('\n')
      ss = 'NumberOfBands %i\n' % (2 * nbnd)
      hh.write(ss)
      ss = 'NumberOfOccupiedBands %i\n' % nbnd
      hh.write(ss)
      hh.write('\n')
      hh.write('KPointsScale reciprocallatticevectors\n')
      hh.write('KPointsList 1\n')
      hh.write('0.000000000  0.000000000  0.000000000  1.0\n')
      hh.write('\n')
      hh.close()
   return ierr

def xyz_read(file):
   ierr = 0
   mat = mat_create()
   mat['nm'] = 'xyz'
   mat['fc'] = format_index2mat[1]
   mat['lc'] = 1.0
   mat['fv'] = format_index2mat[1]
   mat['fp'] = format_index2mat[1]
   try:
      hh = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      rr = hh.readlines()
      hh.close()
      ss = rr[0]
      tt = ss.split()
      mat['na'] = int(tt[0])
      ss = rr[1]
      tt = ss.split()
      mat['nm'] = ''
      for jj in range(len(tt)):
         if jj > 0:
            mat['nm'] += ' '
         mat['nm'] += tt[jj]
      for ii in range(mat['na']):
         ss = rr[ii + 2]
         tt = ss.split()
         mat['as'].append(periodic_table[index_by_symbol(tt[0])]['number'])
         mat['ap'].append([])
         for jj in range(3):
            mat['ap'][ii].append(float(tt[jj + 1]))
   mat['ns'] = len(mat_species(mat['na'], mat['as']))
   format = mat['fv']
   mat['fo'] = format
   axes = [1, 1, 1]
   distance = 10.0
   mat = mat_lattice(mat, axes, distance, format)
   return ierr, mat

def xyz_write(file, mat):
   ierr = 0
   format = format_index2mat[1]
   fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fp != format:
      ierr = 1
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      ss = '%i\n' % mat['na']
      hh.write(ss)
      ss = '%s\n' % mat['nm']
      hh.write(ss)
      for ii in range(mat['na']):
         ss = '  %2s' % periodic_table[index_by_number(mat['as'][ii])]['symbol']
         for jj in range(3):
            ss += ' %13.9f' % ap[ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.close()
   return ierr

def xsf_read(file):
   ierr = 0
   mat = mat_create()
   mat['nm'] = 'xsf'
   mat['fc'] = format_index2mat[1]
   mat['lc'] = 1.0
   mat['fv'] = format_index2mat[1]
   mat['fp'] = format_index2mat[1]
   try:
      hh = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      rr = hh.readlines()
      hh.close()
      ii = 1
      iv = 0
      ip = 0
      for ss in rr:
         tt = ss.split()
         if len(tt) > 0:
            if tt[0].lower() == 'primvec':
               iv = ii
            elif tt[0].lower() == 'primcoord':
               ip = ii
         ii += 1
      if iv > 0:
         for ii in range(3):
            ss = rr[iv + ii]
            tt = ss.split()
            for jj in range(3):
               mat['lv'][ii][jj] = float(tt[jj])
      if ip > 0:
         ss = rr[ip]
         tt = ss.split()
         mat['na'] = int(tt[0])
         for ii in range(mat['na']):
            ss = rr[ip + ii + 1]
            tt = ss.split()
            if tt[0].isdigit():
               mat['as'].append(int(tt[0]))
            else:
               index = index_by_symbol(tt[0])
               number = periodic_table[index]['number']
               mat['as'].append(number)
            mat['ap'].append([])
            for jj in range(3):
               mat['ap'][ii].append(float(tt[jj + 1]))
   mat['ns'] = len(mat_species(mat['na'], mat['as']))
   format = mat['fv']
   mat['fo'] = format
   fr, rmin, rmax = mat_range(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fr == format:
      for jj in range(3):
         mat['lo'][jj] = (rmin[jj] + rmax[jj]) / 2
         for ii in range(3):
            mat['lo'][jj] -= mat['lv'][ii][jj] / 2
   return ierr, mat

def xsf_write(file, mat):
   ierr = 0
   format = format_index2mat[1]
   fv, lv = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], format)
   fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fv != format or fp != format:
      ierr = 1
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
            ss += ' %11.6f' % lv[ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.write('CONVVEC\n')
      for ii in range(3):
         ss = ''
         for jj in range(3):
            ss += ' %11.6f' % lv[ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.write('PRIMCOORD\n')
      ss = ' %4i %i\n' % (mat['na'], 1)
      hh.write(ss)
      for ii in range(mat['na']):
         ss = ' %4i' % mat['as'][ii]
         for jj in range(3):
            ss += ' %11.6f' % ap[ii][jj]
         for jj in range(3):
            ss += ' %11.6f' % 0.0
         ss += '\n'
         hh.write(ss)
      hh.close()
   return ierr

def wien_write(file, mat):
   ierr = 0
   formatv = format_index2mat[2]
   if mat['fp'] != format_index2mat[3]:
      formatp = format_index2mat[2]
   else:
      formatp = format_index2mat[3]
   fv, lv = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], formatv)
   fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], formatp)
   if fv != formatv or fp != formatp:
      ierr = 1
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      if mat['fc'] == format_index2mat[0]:
         ss = 'BOHR'
      else:
         ss = 'ANG'
      if fp == format_index2mat[2]:
         ss += ' C L\n'
      else:
         ss += ' F L\n'
      hh.write(ss)
      ss = '     %13.9f\n' % mat['lc']
      hh.write(ss)
      for ii in range(3):
         ss = '    '
         for jj in range(3):
            ss += ' %13.9f' % lv[ii][jj]
         ss += '\n'
         hh.write(ss)
      for ii in range(mat['na']):
         ss = '  %2s' % periodic_table[index_by_number(mat['as'][ii])]['symbol']
         for jj in range(3):
            ss += ' %13.9f' % ap[ii][jj]
         ss += '\n'
         hh.write(ss)
      hh.close()
   return ierr

def povray_bonds_st(mat):
   ierr = 0
   bonds = []
   format = format_index2mat[1]
   fv, lv = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], format)
   fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fv == format and fp == format:
      for i1 in range(mat['na']):
         bonds.append([])
         n1 = index_by_number(mat['as'][i1])
         v1 = periodic_table[n1]['valence']
         r1 = periodic_table[n1]['rcov']
         ii = 0
         for i2 in range(mat['na']):
            n2 = index_by_number(mat['as'][i2])
            v2 = periodic_table[n2]['valence']
            r2 = periodic_table[n2]['rcov']
            for j1 in range(-1, 2):
               for j2 in range(-1, 2):
                  for j3 in range(-1, 2):
                     if (j3 != 0 or j2 != 0 or j1 != 0 or i2 != i1) and n2 != 0 and n1 != 0:
                        pp = []
                        for jj in range(3):
                           pp.append(lv[2][jj] * j3 + lv[1][jj] * j2 + lv[0][jj] * j1 + ap[i2][jj] - ap[i1][jj])
                        rr = 0.0
                        for jj in range(3):
                           rr += math.pow(pp[jj], 2)
                        rr = math.sqrt(rr)
                        if rr < minimumbondingdistance:
                           ierr = 1
                        if rr < r1 + r2 + bondtolerance:
                           if strictvalence == 1 and ii >= v1:
                              ierr = 1
                           bonds[i1].append([])
                           for jj in range(3):
                              bonds[i1][ii].append(ap[i1][jj] + pp[jj] * r1 / (r1 + r2))
                           ii += 1
   else:
      ierr = 1
   return ierr, bonds

def povray_bonds_ut(mat, umat, smin, smax, si):
   ierr = 0
   bonds = []
   format = format_index2mat[1]
   fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   ufv, ulv = mat_fv_lv(umat['fc'], umat['lc'], umat['fv'], umat['lv'], format)
   ufp, uap = mat_fp_ap(umat['na'], umat['fc'], umat['lc'], umat['fv'], umat['lv'], umat['fp'], umat['ap'], format)
   if fp == format and ufv == format and ufp == format:
      for i1 in range(mat['na']):
         bonds.append([])
         n1 = index_by_number(mat['as'][i1])
         v1 = periodic_table[n1]['valence']
         r1 = periodic_table[n1]['rcov']
         ii = 0
         for i2 in range(umat['na']):
            n2 = index_by_number(umat['as'][i2])
            v2 = periodic_table[n2]['valence']
            r2 = periodic_table[n2]['rcov']
            for j1 in range(smin[0] - 1, smax[0] + 2):
               for j2 in range(smin[1] - 1, smax[1] + 2):
                  for j3 in range(smin[2] - 1, smax[2] + 2):
                     if (j3 != si[i1][3] or j2 != si[i1][2] or j1 != si[i1][1] or i2 != si[i1][0]) and n2 != 0 and n1 != 0:
                        pp = []
                        for jj in range(3):
                           pp.append(ulv[2][jj] * j3 + ulv[1][jj] * j2 + ulv[0][jj] * j1 + uap[i2][jj] - ap[i1][jj])
                        rr = 0.0
                        for jj in range(3):
                           rr += math.pow(pp[jj], 2)
                        rr = math.sqrt(rr)
                        if rr < minimumbondingdistance:
                           ierr = 1
                        if rr < r1 + r2 + bondtolerance:
                           if strictvalence == 1 and ii >= v1:
                              ierr = 1
                           bonds[i1].append([])
                           for jj in range(3):
                              bonds[i1][ii].append(ap[i1][jj] + pp[jj] * r1 / (r1 + r2))
                           ii += 1
   else:
      ierr = 1
   return ierr, bonds

def povray_write(file, mat, box, basis, bonds):
   cameralocation = [-0.1, -0.5, 0.9]
   cameralookat = [0.0, 0.0, 0.0]
   camerafactor = 3.0
   ierr = 0
   format = format_index2mat[1]
   fv, lv = mat_fv_lv(mat['fc'], mat['lc'], mat['fv'], mat['lv'], format)
   fo, lo = mat_fo_lo(mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fo'], mat['lo'], format)
   fp, ap = mat_fp_ap(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
   if fo != format or fv != format or fp != format:
      ierr = 1
   lvcenter = []
   for jj in range(3):
      lvcenter.append(lo[jj])
   for ii in range(3):
      for jj in range(3):
         lvcenter[jj] += lv[ii][jj] / 2
   if box != 0:
      lvcorner = []
      for ii in range(8):
         lvcorner.append([])
         for jj in range(3):
            lvcorner[ii].append(lvcenter[jj])
            for kk in range(3):
               lvcorner[ii][jj] += lv[kk][jj] * (vertex[ii][kk] - 0.5)
   lvlength = []
   for ii in range(3):
      lvlength.append(0.0)
      for jj in range(3):
         lvlength[ii] += math.pow(lv[ii][jj], 2)
      lvlength[ii] = math.sqrt(lvlength[ii])
   lvscale = 0.0
   for ii in range(3):
      lvscale += lvlength[ii]
   lvscale /= 3
   cameradistance = 0.0
   for ii in range(3):
      cameradistance += math.pow(cameralocation[ii] - cameralookat[ii], 2)
   cameradistance = math.sqrt(cameradistance)
   camerascale = camerafactor * lvscale / cameradistance
   cameralocationlv = []
   for jj in range(3):
      cameralocationlv.append(lvcenter[jj] / camerascale + cameralocation[jj])
   cameralookatlv = []
   for jj in range(3):
      cameralookatlv.append(lvcenter[jj] / camerascale + cameralookat[jj])
   cs = mat_species(mat['na'], mat['as'])
   symbol = []
   ratom = []
   rbond = []
   color = []
   for ii in range(mat['ns']):
      index = index_by_number(cs[ii])
      symbol.append(periodic_table[index]['symbol'].lower())
      ratom.append(atomsize * periodic_table[index]['rvdw'])
      rbond.append(bondradius)
      color.append(periodic_table[index]['color'])
   atomicsymbol = []
   for ii in range(mat['na']):
      index = index_by_number(mat['as'][ii])
      atomicsymbol.append(periodic_table[index]['symbol'].lower())
   try:
      hh = open(file, 'w')
   except:
      ierr = 1
   if ierr == 0:
      hh.write('\n')
      ss = '// name of matter: %s\n' % mat['nm']
      hh.write(ss)
      hh.write('\n')
      ss = '#declare camera_location = <%.2f, %.2f, %.2f>;\n' % tuple(cameralocationlv)
      hh.write(ss)
      ss = '#declare camera_look_at = <%.2f, %.2f, %.2f>;\n' % tuple(cameralookatlv)
      hh.write(ss)
      ss = '#declare camera_scale = %.2f;\n' % camerascale
      hh.write(ss)
      hh.write('#declare light_location = camera_location - camera_look_at;\n')
      hh.write('#declare light_scale = 1e6;\n')
      hh.write('#declare color_light = rgb <2.00, 2.00, 2.00>;\n')
      hh.write('#declare color_background = rgb <0.00, 0.00, 0.00>;\n')
      if box != 0:
         hh.write('#declare radius_frame = 0.01;\n')
         hh.write('#declare color_frame = rgb <0.75, 0.75, 0.75>;\n')
         hh.write('#declare color_box = rgbf <1.00, 1.00, 1.00, 0.75>;\n')
      if basis != 0:
         axislength = lvscale / 3
         axisradius = axislength / 50
         ss = '#declare length_axis = %.2f;\n' % axislength
         hh.write(ss)
         ss = '#declare radius_axis = %.2f;\n' % axisradius
         hh.write(ss)
         hh.write('#declare color_axis_x = rgb <1.00, 0.00, 0.00>;\n')
         hh.write('#declare color_axis_y = rgb <0.00, 1.00, 0.00>;\n')
         hh.write('#declare color_axis_z = rgb <0.00, 0.00, 1.00>;\n')
         hh.write('#declare length_arow = 0.2 * length_axis;\n')
         hh.write('#declare radius_arow = 2.0 * radius_axis;\n')
         hh.write('#declare color_arow_x = color_axis_x;\n')
         hh.write('#declare color_arow_y = color_axis_y;\n')
         hh.write('#declare color_arow_z = color_axis_z;\n')
      for ii in range(mat['ns']):
         ss = '#declare radius_atom_%s = %.2f;\n' % (symbol[ii], ratom[ii])
         hh.write(ss)
         if cs[ii] != 0:
            ss = '#declare radius_bond_%s = %.2f;\n' % (symbol[ii], rbond[ii])
            hh.write(ss)
         ss = '#declare color_atom_%s = rgb <%.2f, %.2f, %.2f>;\n' % (tuple([symbol[ii]]) + tuple(color[ii]))
         hh.write(ss)
         if cs[ii] != 0:
            ss = '#declare color_bond_%s = rgb <%.2f, %.2f, %.2f>;\n' % (tuple([symbol[ii]]) + tuple(color[ii]))
            hh.write(ss)
      hh.write('\n')
      hh.write('camera { location camera_location sky <0.00, 0.00, 1.00> up <0.00, 0.00, 1.00> right <-1.33, 0.00, 0.00> direction <0.00, -1.00, 0.00> look_at camera_look_at scale camera_scale }\n')
      hh.write('light_source { light_location color color_light shadowless scale light_scale }\n')
      hh.write('background { color color_background }\n')
      hh.write('\n')
      if box != 0:
         hh.write('union {\n')
         for ii in range(8):
            ss = 'sphere { <%.9f, %.9f, %.9f>, radius_frame }\n' % tuple(lvcorner[ii])
            hh.write(ss)
         for ii in range(12):
            ss = 'cylinder { <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, radius_frame }\n' % (tuple(lvcorner[edge[ii][0]]) + tuple(lvcorner[edge[ii][1]]))
            hh.write(ss)
         hh.write('texture { pigment { color color_frame } }\n')
         hh.write('}\n')
         hh.write('\n')
         hh.write('union {\n')
         for ii in range(12):
            ss = 'triangle { <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f> }\n' % (tuple(lvcorner[face[ii][0]]) + tuple(lvcorner[face[ii][1]]) + tuple(lvcorner[face[ii][2]]))
            hh.write(ss)
         hh.write('texture { pigment { color color_box } }\n')
         hh.write('}\n')
         hh.write('\n')
      if basis != 0:
         hh.write('union {\n')
         hh.write('cylinder { <0.00, 0.00, 0.00>, <length_axis, 0.00, 0.00>, radius_axis texture { pigment { color color_axis_x } } }\n')
         hh.write('cylinder { <0.00, 0.00, 0.00>, <0.00, length_axis, 0.00>, radius_axis texture { pigment { color color_axis_y } } }\n')
         hh.write('cylinder { <0.00, 0.00, 0.00>, <0.00, 0.00, length_axis>, radius_axis texture { pigment { color color_axis_z } } }\n')
         hh.write('cone { <length_axis, 0.00, 0.00>, radius_arow <length_axis + length_arow, 0.00, 0.00>, 0.00 texture { pigment { color color_arow_x } } }\n')
         hh.write('cone { <0.00, length_axis, 0.00>, radius_arow <0.00, length_axis + length_arow, 0.00>, 0.00 texture { pigment { color color_arow_y } } }\n')
         hh.write('cone { <0.00, 0.00, length_axis>, radius_arow <0.00, 0.00, length_axis + length_arow>, 0.00 texture { pigment { color color_arow_z } } }\n')
         hh.write('}\n')
         hh.write('\n')
      hh.write('union {\n')
      for ii in range(mat['na']):
         ss = 'sphere { <%.9f, %.9f, %.9f>, radius_atom_%s texture { pigment { color color_atom_%s } } }\n' % (tuple(ap[ii]) + tuple([atomicsymbol[ii]]) + tuple([atomicsymbol[ii]]))
         hh.write(ss)
      hh.write('}\n')
      hh.write('\n')
      hh.write('union {\n')
      for ii in range(mat['na']):
         for jj in range(len(bonds[ii])):
            ss = 'cylinder { <%.9f, %.9f, %.9f>, <%.9f, %.9f, %.9f>, radius_bond_%s texture { pigment { color color_bond_%s } } }\n' % (tuple(ap[ii]) + tuple(bonds[ii][jj]) + tuple([atomicsymbol[ii]]) + tuple([atomicsymbol[ii]]))
            hh.write(ss)
      hh.write('}\n')
      hh.write('\n')
      hh.close()
   return ierr

if __name__ == "__main__":
   import sys
   sys.exit(main())

