#!/usr/bin/env python

#-------------------------------------------------------------------------------
#
#   convert.py
#   command-line driver for library matter
#   written by Georgy Samsonidze (October 2008)
#
#-------------------------------------------------------------------------------

from common import eps9, inf9, format_mat2index, format_index2mat, periodic_table, index_by_number, boolean_input
from matrix import rotation_matrix
from matter import mat_check, mat_format, mat_lattice, mat_replicate, mat_translate, mat_rotate, mat_fp_ap, mat_range, mat_read, mat_write, paratec_read, paratec_write, vasp_read, vasp_write, espresso_read, espresso_write, siesta_read, siesta_write, tbpw_read, tbpw_write, xyz_read, xyz_write, xsf_read, xsf_write, wien_write, povray_bonds_st, povray_bonds_ut, povray_write

def main(argv = None):
   if argv is None:
      argv = sys.argv
   if not (len(argv) == 5 or (len(argv) > 5 and ((argv[5].lower() == 'units' and len(argv) == 9) or (argv[5].lower() == 'vacuum' and len(argv) == 9) or (argv[5].lower() == 'supercell' and (len(argv) == 7 or len(argv) == 20)) or (argv[5].lower() == 'translate' and len(argv) == 10) or (argv[5].lower() == 'rotate' and (len(argv) == 8 or len(argv) == 9)) or (argv[5].lower() == 'center' and len(argv) == 8) or (argv[5].lower() == 'minbox' and len(argv) == 7)))):
      sys.exit("\n   Usage: convert.py ifn iff ofn off [operation parameters]\n" +
               "      ifn = input file name         ofn = output file name\n" +
               "      iff = input file format       off = output file format\n" +
               "   format = mat|paratec|vasp|espresso|siesta|tbpw|xyz|xsf|wien|povray\n" +
               "                                   (wien & povray are only for output)\n\n" +
               "  operation  parameters\n" +
               "  ---------------------\n" +
               "      units  lattice_constant_units lattice_vector_units atomic_position_units\n" +
               "     vacuum  axes separation_distance separation_distance_units\n" +
               "  supercell  ( sc_file | sc_origin sc_vectors sc_units sc_translation )\n" +
               "  translate  translation_vector translation_vector_units\n" +
               "     rotate  axis ( rotation_angle | adjacent_side opposite_side )\n" +
               "     center  axes center_type\n" +
               "     minbox  axes\n\n" +
               "             axis = x | y | z             rotation_angle in degrees\n" +
               "             axes = x & y & z             *_side in arbitrary cartesian units\n" +
               "   sc_translation = true | false          center_type = geometric | mass\n" +
               "          *_units = bohr | angstrom | latconst | latvec\n")
   ifn = argv[1]
   ift = argv[2].lower()
   ofn = argv[3]
   oft = argv[4].lower()
   if ift == 'mat':
      ierr, mat = mat_read(ifn)
   elif ift == 'paratec':
      ierr, mat = paratec_read(ifn)
   elif ift == 'vasp':
      ierr, mat = vasp_read(ifn)
   elif ift == 'espresso':
      ierr, mat = espresso_read(ifn)
   elif ift == 'siesta':
      ierr, mat = siesta_read(ifn)
   elif ift == 'tbpw':
      ierr, mat = tbpw_read(ifn)
   elif ift == 'xyz':
      ierr, mat = xyz_read(ifn)
   elif ift == 'xsf':
      ierr, mat = xsf_read(ifn)
   else:
      sys.exit("\n   Error: invalid input format %s\n" % ift)
   if ierr != 0:
      sys.exit("\n   Error: unable to read file %s\n" % ifn)
   ierr = mat_check(mat)
   if ierr != 0:
      sys.exit("\n   Error: invalid content in file %s\n" % ifn)
   if oft == 'povray':
      st = 1
      umat = {}
      smin = []
      smax = []
      si = []
   if len(argv) != 5:
      if argv[5].lower() == 'units':
         ierr, fc = format_input(argv[6], 2)
         if ierr != 0:
            sys.exit("\n   Error: invalid lattice_constant_units %s\n" % argv[6])
         ierr, fo = format_input(argv[7], 3)
         if ierr != 0:
            sys.exit("\n   Error: invalid lattice_vector_units %s\n" % argv[7])
         ierr, fp = format_input(argv[8], 4)
         if ierr != 0:
            sys.exit("\n   Error: invalid atomic_position_units %s\n" % argv[8])
         if fo == format_index2mat[3]:
            fv = format_index2mat[2]
         else:
            fv = fo
         mat = mat_format(mat, fc, fo, fv, fp)
      elif argv[5].lower() == 'vacuum':
         axes = axes_input(argv[6])
         try:
            distance = float(argv[7])
         except:
            sys.exit("\n   Error: invalid separation_distance %s\n" % argv[7])
         ierr, format = format_input(argv[8], 3)
         if ierr != 0:
            sys.exit("\n   Error: invalid separation_distance_units %s\n" % argv[8])
         mat = mat_lattice(mat, axes, distance, format)
      elif argv[5].lower() == 'supercell':
         if len(argv) == 7:
            ierr = 0
            try:
               hh = open(argv[6], 'r')
            except:
               ierr = 1
            if ierr == 0:
               rr = hh.readlines()
               hh.close()
               ss = rr[0]
               tt = ss.split()
               so = [float(tt[0]), float(tt[1]), float(tt[2])]
               sv = []
               for ii in range(3):
                  ss = rr[ii + 1]
                  tt = ss.split()
                  sv.append([float(tt[0]), float(tt[1]), float(tt[2])])
               ss = rr[4]
               tt = ss.split()
               jer, format = format_input(tt[0], 4)
               ierr += jer
               ss = rr[5]
               tt = ss.split()
               jer, st = boolean_input(tt[0])
               ierr += jer
            if ierr != 0:
               sys.exit("\n   Error: unable to read supercell_file %s\n" % argv[6])
         else:
            try:
               so = [float(argv[6]), float(argv[7]), float(argv[8])]
            except:
               sys.exit("\n   Error: invalid supercell_origin %s %s %s\n" % (argv[6], argv[7], argv[8]))
            try:
               sv = [[float(argv[9]), float(argv[10]), float(argv[11])], [float(argv[12]), float(argv[13]), float(argv[14])], [float(argv[15]), float(argv[16]), float(argv[17])]]
            except:
               sys.exit("\n   Error: invalid supercell_vectors %s %s %s %s %s %s %s %s %s\n" % (argv[9], argv[10], argv[11], argv[12], argv[13], argv[14], argv[15], argv[16], argv[17]))
            ierr, format = format_input(argv[18], 4)
            if ierr != 0:
               sys.exit("\n   Error: invalid supercell_units %s\n" % argv[18])
            ierr, st = boolean_input(argv[19])
            if ierr != 0:
               sys.exit("\n   Error: invalid supercell_translation %s\n" % argv[19])
         if st == 0:
            umat = copy.deepcopy(mat)
         smin, smax, si, mat = mat_replicate(mat, so, sv, format)
      elif argv[5].lower() == 'translate':
         try:
            translation = [float(argv[6]), float(argv[7]), float(argv[8])]
         except:
            sys.exit("\n   Error: invalid translation_vector %s %s %s\n" % (argv[6], argv[7], argv[8]))
         ierr, format = format_input(argv[9], 4)
         if ierr != 0:
            sys.exit("\n   Error: invalid translation_vector_units %s\n" % argv[9])
         mat = mat_translate(mat, translation, format)
      elif argv[5].lower() == 'rotate':
         ierr, axis = axis_input(argv[6])
         if len(argv) == 8:
            try:
               angle = float(argv[7]) * math.pi / 180.0
            except:
               sys.exit("\n   Error: invalid rotation_angle %s\n" % argv[7])
         else:
            try:
               adjacent = float(argv[7])
            except:
               sys.exit("\n   Error: invalid adjacent_side %s\n" % argv[7])
            try:
               opposite = float(argv[8])
            except:
               sys.exit("\n   Error: invalid opposite_side %s\n" % argv[8])
            angle = math.atan2(opposite, adjacent)
         rotation = rotation_matrix(axis, angle)
         mat = mat_rotate(mat, rotation)
      elif argv[5].lower() == 'center':
         axes = axes_input(argv[6])
         ierr, center_type = center_input(argv[7])
         if ierr != 0:
            sys.exit("\n   Error: invalid center_type %s\n" % argv[7])
         format = mat['fp']
         center = []
         if center_type == 0:
            fr, rmin, rmax = mat_range(mat['na'], mat['fc'], mat['lc'], mat['fv'], mat['lv'], mat['fp'], mat['ap'], format)
            if fr != format:
               sys.exit("\n   Error: failed center_type geometric\n")
            for jj in range(3):
               center.append((rmin[jj] + rmax[jj]) / 2)
         else:
            if format == format_index2mat[3]:
               sys.exit("\n   Error: center_type mass incompatible with atomic_position_units %s\n" % format)
            total_mass = 0.0
            for jj in range(3):
               center.append(0.0)
            for ii in range(mat['na']):
               atomic_mass = periodic_table[index_by_number(mat['as'][ii])]['mass']
               total_mass += atomic_mass
               for jj in range(3):
                  center[jj] += mat['ap'][ii][jj] * atomic_mass
            if abs(total_mass) > eps9:
               for jj in range(3):
                  center[jj] /= total_mass
            else:
               sys.exit("\n   Error: failed center_type mass\n")
         translation = []
         for jj in range(3):
            translation.append(-center[jj] * axes[jj])
         mat = mat_translate(mat, translation, format)
      elif argv[5].lower() == 'minbox':
         axes = axes_input(argv[6])
         for ii in range(3):
            if axes[ii] == 1:
               for jj in range(3):
                  if jj != ii:
                     minwidth = inf9
                     minangle = 0.0
                     for kk in range(180):
                        angle = float(kk) * math.pi / 180.0
                        rotation = rotation_matrix(jj, angle)
                        cmat = copy.deepcopy(mat)
                        rmat = mat_rotate(cmat, rotation)
                        format, rmin, rmax = mat_range(rmat['na'], rmat['fc'], rmat['lc'], rmat['fv'], rmat['lv'], rmat['fp'], rmat['ap'], 'bohr')
                        cmat = []
                        rmat = []
                        width = rmax[ii] - rmin[ii]
                        if width < minwidth - eps9:
                           minwidth = width
                           minangle = angle
                     rotation = rotation_matrix(jj, minangle)
                     cmat = copy.deepcopy(mat)
                     mat = mat_rotate(cmat, rotation)
                     cmat = []
   if oft == 'povray':
      box = 1
      basis = 1
      if st != 0:
         ierr, bonds = povray_bonds_st(mat)
      else:
         ierr, bonds = povray_bonds_ut(mat, umat, smin, smax, si)
      if ierr != 0:
         sys.exit("\n   Error: failed povray_bonds\n")
   if oft == 'mat':
      ierr = mat_write(ofn, mat)
   elif oft == 'paratec':
      ierr = paratec_write(ofn, mat)
   elif oft == 'vasp':
      ierr = vasp_write(ofn, mat)
   elif oft == 'espresso':
      ierr = espresso_write(ofn, mat)
   elif oft == 'siesta':
      ierr = siesta_write(ofn, mat)
   elif oft == 'tbpw':
      ierr = tbpw_write(ofn, mat)
   elif oft == 'xyz':
      ierr = xyz_write(ofn, mat)
   elif oft == 'xsf':
      ierr = xsf_write(ofn, mat)
   elif oft == 'wien':
      ierr = wien_write(ofn, mat)
   elif oft == 'povray':
      ierr = povray_write(ofn, mat, box, basis, bonds)
   else:
      sys.exit("\n   Error: invalid output format %s\n" % oft)
   if ierr != 0:
      sys.exit("\n   Error: unable to write file %s\n" % ofn)
   return 0

def format_input(format, range):
   ierr = 0
   index = -1
   format = format.lower()
   if format in format_mat2index:
      index = format_mat2index[format]
   else:
      ierr = 1
   if index >= range:
      ierr = 1
   return ierr, format

def axis_input(axis):
   ierr = 0
   index = -1
   axis = axis.lower()
   if axis == 'x':
      index = 0
   elif axis == 'y':
      index = 1
   elif axis == 'z':
      index = 2
   else:
      ierr = 1
   return ierr, index

def axes_input(axes):
   flags = [0, 0, 0]
   axes = axes.lower()
   if axes.find('x') > -1:
      flags[0] = 1
   if axes.find('y') > -1:
      flags[1] = 1
   if axes.find('z') > -1:
      flags[2] = 1
   return flags

def center_input(center):
   ierr = 0
   index = -1
   center = center.lower()
   if center == 'geometric':
      index = 0
   elif center == 'mass':
      index = 1
   else:
      ierr = 1
   return ierr, index

if __name__ == "__main__":
   import sys
   import math
   import copy
   sys.exit(main())

