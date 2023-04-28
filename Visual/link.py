#!/usr/bin/env python

#-------------------------------------------------------------------------------
#
#   link.py
#   assembles a molecule from two groups of atoms
#   written by Georgy Samsonidze (October 2008)
#
#-------------------------------------------------------------------------------

from common import eps9, inf9, format_index2mat
from matrix import align_vector, matrix_dot_matrix, rotation_matrix
from matter import mat_check, mat_format, mat_remove, mat_merge, mat_translate, mat_rotate, mat_read, mat_write, paratec_read, paratec_write, vasp_read, vasp_write, espresso_read, espresso_write, siesta_read, siesta_write, tbpw_read, tbpw_write, xyz_read, xyz_write, xsf_read, xsf_write

def main(argv = None):
   if argv is None:
      argv = sys.argv
   if len(argv) != 15:
      sys.exit("\n   Usage: link.py fnc ffc fns ffs fna ffa icb icd isb isd acs facs rcs frcs\n" +
               "   Description: links the side chain read from file fns to the core structure\n" +
               "   read from file fnc and writes the obtained adduct to file fna. Variables\n" +
               "   ffc, ffs, ffa set to mat|paratec|vasp|espresso|siesta|tbpw|xyz|xsf specify\n" +
               "   the file formats. The core and side atoms with indices icb and isb are\n" +
               "   bonded together, while the core and side atoms with indices icd and isd\n" +
               "   are detached. The side chain is oriented in such a way that bonds icb-icd,\n" +
               "   isd-isb, and icb-isb are all collinear. Variable acs gives the angle by\n" +
               "   which the side chain is rotated around the icb-isb axis. Variable rcs\n" +
               "   gives the distance between icb and isb atoms. Variables facs and frcs set\n" +
               "   to radian or degree and bohr or angstrom define the units of acs and rcs.\n")
   fnc = argv[1]
   ffc = argv[2].lower()
   fns = argv[3]
   ffs = argv[4].lower()
   fna = argv[5]
   ffa = argv[6].lower()
   try:
      icb = int(argv[7]) - 1
   except:
      icb = -1
   try:
      icd = int(argv[8]) - 1
   except:
      icd = -1
   try:
      isb = int(argv[9]) - 1
   except:
      isb = -1
   try:
      isd = int(argv[10]) - 1
   except:
      isd = -1
   try:
      acs = float(argv[11])
   except:
      sys.exit("\n   Error: acs %s\n" % argv[11])
   facs = argv[12].lower()
   if facs == 'r' or facs == 'rad' or facs == 'radian':
      pass
   elif facs == 'd' or facs == 'deg' or facs == 'degree':
      acs *= math.pi / 180.0
   else:
      sys.exit("\n   Error: facs %s\n" % argv[12])
   try:
      rcs = float(argv[13])
   except:
      rsc = -inf9
   if rcs < eps9:
      sys.exit("\n   Error: rcs %s\n" % argv[13])
   frcs = argv[14].lower()
   if frcs == 'b' or frcs == 'bohr':
      format = format_index2mat[0]
   elif frcs == 'a' or frcs == 'ang' or frcs == 'angstrom':
      format = format_index2mat[1]
   else:
      sys.exit("\n   Error: frcs %s\n" % argv[14])
   if ffc == 'mat':
      ierr, matc = mat_read(fnc)
   elif ffc == 'paratec':
      ierr, matc = paratec_read(fnc)
   elif ffc == 'vasp':
      ierr, matc = vasp_read(fnc)
   elif ffc == 'espresso':
      ierr, matc = espresso_read(fnc)
   elif ffc == 'siesta':
      ierr, matc = siesta_read(fnc)
   elif ffc == 'tbpw':
      ierr, matc = tbpw_read(fnc)
   elif ffc == 'xyz':
      ierr, matc = xyz_read(fnc)
   elif ffc == 'xsf':
      ierr, matc = xsf_read(fnc)
   else:
      sys.exit("\n   Error: invalid input format %s\n" % ffc)
   if ierr != 0:
      sys.exit("\n   Error: unable to read file %s\n" % fnc)
   ierr = mat_check(matc)
   if ierr != 0:
      sys.exit("\n   Error: invalid content in file %s\n" % fnc)
   matc = mat_format(matc, format, format, format, format)
   if matc['fc'] != format or matc['fo'] != format or matc['fv'] != format or matc['fp'] != format:
      sys.exit("\n   Error: unable to change format in file %s\n" % fnc)
   if icb < 0 or icb > matc['na'] - 1:
      sys.exit("\n   Error: icb %s\n" % argv[7])
   if icd < 0 or icd > matc['na'] - 1 or icd == icb:
      sys.exit("\n   Error: icd %s\n" % argv[8])
   vc = []
   for jj in range(3):
      vc.append(matc['ap'][icd][jj] - matc['ap'][icb][jj])
   rc = 0.0
   for jj in range(3):
      rc += math.pow(vc[jj], 2)
   rc = math.sqrt(rc)
   if rc < eps9:
      sys.exit("\n   Error: atomic positions in file %s\n" % fnc)
   if ffs == 'mat':
      ierr, mats = mat_read(fns)
   elif ffs == 'paratec':
      ierr, mats = paratec_read(fns)
   elif ffs == 'vasp':
      ierr, mats = vasp_read(fns)
   elif ffs == 'espresso':
      ierr, mats = espresso_read(fns)
   elif ffs == 'siesta':
      ierr, mats = siesta_read(fns)
   elif ffs == 'tbpw':
      ierr, mats = tbpw_read(fns)
   elif ffs == 'xyz':
      ierr, mats = xyz_read(fns)
   elif ffs == 'xsf':
      ierr, mats = xsf_read(fns)
   else:
      sys.exit("\n   Error: invalid input format %s\n" % ffs)
   if ierr != 0:
      sys.exit("\n   Error: unable to read file %s\n" % fns)
   ierr = mat_check(matc)
   if ierr != 0:
      sys.exit("\n   Error: invalid content in file %s\n" % fns)
   mats = mat_format(mats, format, format, format, format)
   if mats['fc'] != format or mats['fo'] != format or mats['fv'] != format or mats['fp'] != format:
      sys.exit("\n   Error: unable to change format in file %s\n" % fns)
   if isb < 0 or isb > mats['na'] - 1:
      sys.exit("\n   Error: isb %s\n" % argv[9])
   if isd < 0 or isd > mats['na'] - 1 or isd == isb:
      sys.exit("\n   Error: isd %s\n" % argv[10])
   vs = []
   for jj in range(3):
      vs.append(mats['ap'][isb][jj] - mats['ap'][isd][jj])
   rs = 0.0
   for jj in range(3):
      rs += math.pow(vs[jj], 2)
   rs = math.sqrt(rs)
   if rs < eps9:
      sys.exit("\n   Error: atomic positions in file %s\n" % fns)
   numberc, axesc, anglesc = align_vector(vc)
   numbers, axess, angless = align_vector(vs)
   number = numbers + 1 + numberc
   axes = []
   angles = []
   for ii in range(numbers):
      axes.append(axess[ii])
      angles.append(angless[ii])
   axes.append(2)
   angles.append(acs)
   for ii in range(numberc - 1, -1, -1):
      axes.append(axesc[ii])
      angles.append(-anglesc[ii])
   rotation = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
   for ii in range(number):
      rotation = matrix_dot_matrix(rotation_matrix(axes[ii], angles[ii]), rotation)
   mats = mat_rotate(mats, rotation)
   translation = [0.0, 0.0, 0.0]
   for jj in range(3):
      translation[jj] += matc['ap'][icb][jj] - mats['ap'][isb][jj] + vc[jj] * rcs / rc
   mats = mat_translate(mats, translation, format)
   matc = mat_remove(matc, icd)
   mats = mat_remove(mats, isd)
   mata = mat_merge(matc, mats)
   if ffa == 'mat':
      ierr = mat_write(fna, mata)
   elif ffa == 'paratec':
      ierr = paratec_write(fna, mata)
   elif ffa == 'vasp':
      ierr = vasp_write(fna, mata)
   elif ffa == 'espresso':
      ierr = espresso_write(fna, mata)
   elif ffa == 'siesta':
      ierr = siesta_write(fna, mata)
   elif ffa == 'tbpw':
      ierr = tbpw_write(fna, mata)
   elif ffa == 'xyz':
      ierr = xyz_write(fna, mata)
   elif ffa == 'xsf':
      ierr = xsf_write(fna, mata)
   else:
      sys.exit("\n   Error: invalid output format %s\n" % ffa)
   if ierr != 0:
      sys.exit("\n   Error: unable to write file %s\n" % fna)
   return 0

if __name__ == "__main__":
   import sys
   import math
   sys.exit(main())

