#!/usr/bin/env python

#-------------------------------------------------------------------------------
#
#   average.py
#   averages potential V(x, y, z) over the volume of the unit cell; on the
#   faces x = xpos | y = ypos | z = zpos; over (x, y) planes for the slab
#   in z direction; over the bulk-like region of the slab zmin < z < zmax
#   written by Georgy Samsonidze (April 2009)
#
#-------------------------------------------------------------------------------

from common import bohr, rydberg, inf9

def main(argv = None):
   if argv is None:
      argv = sys.argv
   if not ((len(argv) == 4 and argv[3].lower() == 'volume') or len(argv) == 7):
      sys.exit("\n   Usage: average.py ivfn ivff mode [xpos ypos zpos | zout nsa ofn]\n" +
               "   ivfn = input volumetric file name\n" +
               "   ivff = input volumetric file format (cube)\n" +
               "   mode = volume | faces | planes | slab\n" +
               "   xpos, ypoz, zpos = positions of yz, zx, xy faces in crystal coordinates,\n" +
               "       if outside [0, 1) then the face is not included in the averaging\n" +
               "   zout = position anywhere outside the slab in Angstrom, used to place\n" +
               "       the slab in the center of the unit cell\n" +
               "   nsa  = number of surface atoms on each side of the slab\n" +
               "   ofn  = output file name\n\n" +
               "   mode = volume : potential is averaged over the volume of the unit cell\n" +
               "   mode = faces  : potential is averaged on xpos, ypoz, zpos faces\n" +
               "   mode = planes : slab is centered using zout, V is averaged over (x, y)\n" +
               "       planes and written to ofn (two columns, z in Angstrom & V in eV)\n" +
               "   mode = slab   : slab is centered using zout, boundaries of bulk-like\n" +
               "       region (zmin, zmax) are determined using nsa (requires nsa > 0),\n" +
               "       V is averaged over zmin < z < zmax, zmin & zmax in Angstrom and\n" +
               "       V in eV are written to ofn\n")
   ivfn = argv[1]
   ivff = argv[2].lower()
   mode = argv[3].lower()
   if mode == 'volume':
      pass
   elif mode == 'faces':
      try:
         xpos = float(argv[4])
      except:
         sys.exit("\n   Error: invalid xpos %s\n" % argv[4])
      try:
         ypos = float(argv[5])
      except:
         sys.exit("\n   Error: invalid ypos %s\n" % argv[5])
      try:
         zpos = float(argv[6])
      except:
         sys.exit("\n   Error: invalid zpos %s\n" % argv[6])
   elif mode == 'planes' or mode == 'slab':
      try:
         zout = float(argv[4])
      except:
         sys.exit("\n   Error: invalid zout %s\n" % argv[4])
      try:
         nsa = int(argv[5])
         if mode == 'slab' and nsa < 1:
            raise Exception
      except:
         sys.exit("\n   Error: invalid nsa %s\n" % argv[5])
      ofn = argv[6]
   else:
      sys.exit("\n   Error: invalid mode %s\n" % mode)

   if ivff == 'cube':
      ierr, na, aspecies, acharge, aposition, grid, origin, step, vol = cub_read(ivfn)
   else:
      sys.exit("\n   Error: invalid input format %s\n" % ivff)
   if ierr != 0:
      sys.exit("\n   Error: unable to read file %s\n" % ivfn)

   if mode == 'volume' or mode == 'faces':
      if mode == 'volume':
         title = 'over the volume of the unit cell'
      if mode == 'faces':
         ipos = int(round(xpos * float(grid[0])))
         jpos = int(round(ypos * float(grid[1])))
         kpos = int(round(zpos * float(grid[2])))
         title = 'on the faces'
         nfaces = 0
         if ipos >=0 and ipos < grid[0]:
            title += ' x = %f' % xpos
            nfaces += 1
         if jpos >=0 and jpos < grid[1]:
            title += ' y = %f' % ypos
            nfaces += 1
         if kpos >=0 and kpos < grid[2]:
            title += ' z = %f' % zpos
            nfaces += 1
         if nfaces < 1:
            sys.exit("\n   Error: no faces found %s %s %s\n" % (argv[4], argv[5], argv[6]))
      print('\n   ' + ivfn + ' : averaging ' + title)

      navg = 0
      vavg = 0.0
      vmin = inf9
      vmax = -inf9

      if mode == 'volume':
         for ii in range(grid[0]):
            for jj in range(grid[1]):
               for kk in range(grid[2]):
                  vcur = vol[ii][jj][kk]
                  vavg += vcur
                  if vcur < vmin:
                     vmin = vcur
                  if vcur > vmax:
                     vmax = vcur
         navg += grid[0]*grid[1]*grid[2]

      if mode == 'faces':
         if ipos >=0 and ipos < grid[0]:
            ii = ipos
            for jj in range(grid[1]):
               for kk in range(grid[2]):
                  vcur = vol[ii][jj][kk]
                  vavg += vcur
                  if vcur < vmin:
                     vmin = vcur
                  if vcur > vmax:
                     vmax = vcur
            navg += grid[1]*grid[2]

         if jpos >=0 and jpos < grid[1]:
            jj = jpos
            for ii in range(grid[0]):
               for kk in range(grid[2]):
                  vcur = vol[ii][jj][kk]
                  vavg += vcur
                  if vcur < vmin:
                     vmin = vcur
                  if vcur > vmax:
                     vmax = vcur
            navg += grid[0]*grid[2]

         if kpos >=0 and kpos < grid[2]:
            kk = kpos
            for ii in range(grid[0]):
               for jj in range(grid[1]):
                  vcur = vol[ii][jj][kk]
                  vavg += vcur
                  if vcur < vmin:
                     vmin = vcur
                  if vcur > vmax:
                     vmax = vcur
            navg += grid[0]*grid[1]

      vavg /= float(navg)
      vavg *= rydberg
      vmin *= rydberg
      vmax *= rydberg
      print("\n   Average potential = %.6f eV" % vavg)
      print("   Minimum potential = %.6f eV" % vmin)
      print("   Maximum potential = %.6f eV\n" % vmax)

   if mode == 'planes' or mode == 'slab':
      zstep = step[2][2]
      zlen = zstep * float(grid[2])
      zlist = []
      for ia in range(na):
         zlist.append(aposition[ia][2])
      for ia in range(na):
         fdummy = (zlist[ia] - zout) / zlen
         idummy = int(fdummy)
         if fdummy < 0.0:
            idummy -= 1
         zlist[ia] -= zlen * float(idummy)
      zlist.sort()
      zcenter = (zlist[0] + zlist[na - 1]) / 2.0
      n0 = int((zcenter - zlen / 2.0) / zstep)
      s1 = ''
      s2 = ''

      if mode == 'planes':
         vavg = []
         for kk in range(grid[2]):
            xx = 0.0
            for ii in range(grid[0]):
               for jj in range(grid[1]):
                  xx += vol[ii][jj][kk]
            xx /= float(grid[0] * grid[1])
            xx *= rydberg
            vavg.append(xx)

         zcoord = 0.0
         for kk in range(grid[2] + 1):
            k0 = kk + n0
            fdummy = float(k0) / float(grid[2])
            idummy = int(fdummy)
            if fdummy < 0.0:
               idummy -= 1
            k0 -= grid[2] * idummy
            if kk == 0:
               s1 += ivfn + ' : vacuum level = %f eV' % vavg[k0]
            s2 += '%f %f\n' % (zcoord, vavg[k0])
            zcoord += step[2][2] * bohr

      if mode == 'slab':
         zstart = (zlist[nsa - 1] + zlist[nsa]) / 2.0
         n1 = int(zstart / zstep)
         n2 = n1 + 1
         z1 = zstep * float(n1)
         z2 = zstep * float(n2)
         if zstart - z1 < z2 - zstart:
            w1 = ((z1 + z2) / 2.0 - zstart) / zstep
            w2 = 1.0
         else:
            w1 = 0.0
            w2 = 1.0 - (zstart - (z1 + z2) / 2.0) / zstep

         zend = (zlist[na - 1 - nsa] + zlist[na - nsa]) / 2.0
         n3 = int(zend / zstep)
         n4 = n3 + 1
         z3 = zstep * float(n3)
         z4 = zstep * float(n4)
         if zend - z3 < z4 - zend:
            w3 = 1.0 - ((z3 + z4) / 2.0 - zend) / zstep
            w4 = 0.0
         else:
            w3 = 1.0
            w4 = (zend - (z3 + z4) / 2.0) / zstep

         wlist = []
         wlist.append(w1)
         wlist.append(w2)
         for kk in range(n2 + 1, n3):
            wlist.append(1.0)
         wlist.append(w3)
         wlist.append(w4)

         vavg = 0.0
         norm = 0.0
         for ii in range(grid[0]):
            for jj in range(grid[1]):
               for kk in range(n1, n4 + 1):
                  vavg += vol[ii][jj][kk] * wlist[kk - n1]
                  norm += wlist[kk - n1]
         vavg /= norm
         vavg *= rydberg

         zstart -= zstep * float(n0)
         zend -= zstep * float(n0)
         zstart *= bohr
         zend *= bohr

         s1 += ivfn + ' : average in bulk-like region = %f eV' % vavg

         s2 += '%f %f\n' % (zstart, -1.0e3)
         s2 += '%f %f\n\n' % (zstart, 1.0e3)
         s2 += '%f %f\n' % (zend, -1.0e3)
         s2 += '%f %f\n\n' % (zend, 1.0e3)
         s2 += '%f %f\n' % (-1.0e3, vavg)
         s2 += '%f %f\n' % (1.0e3, vavg)

      zshift = zcenter / zlen - 0.5
      s1 += ' zshift = %f a3' % zshift
      print(s1)

      hh = open(ofn, 'w')
      hh.write(s2)
      hh.close()

   return 0

def cub_read(file):
   ierr = 0
   na = 0
   aspecies = []
   acharge = []
   aposition = []
   grid = []
   origin = []
   step = []
   vol = [[[]]]
   try:
      hh = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      for ii in range(2):
         ss = hh.readline()
      ss = hh.readline()
      tt = ss.split()
      na = int(tt[0])
      origin = []
      for jj in range(3):
         origin.append(float(tt[jj + 1]))
      grid = []
      step = []
      for ii in range(3):
         ss = hh.readline()
         tt = ss.split()
         grid.append(int(tt[0]))
         step.append([])
         for jj in range(3):
            step[ii].append(float(tt[jj + 1]))
      for ii in range(na):
         ss = hh.readline()
         tt = ss.split()
         aspecies.append(int(tt[0]))
         acharge.append(float(tt[1]))
         aposition.append([])
         for jj in range(3):
            aposition[ii].append(float(tt[jj + 2]))
      nn = grid[0] * grid[1] * ((grid[2] - 1) // 6 + 1)
      ii = 0
      jj = 0
      kk = 0
      for mm in range(nn):
         ss = hh.readline()
         tt = ss.split()
         for ll in range(6):
            if kk < grid[2]:
               vol[ii][jj].append(float(tt[ll]))
               kk += 1
         if kk == grid[2]:
            kk = 0
            jj += 1
            if jj < grid[1]:
               vol[ii].append([])
            else:
               kk = 0
               jj = 0
               ii += 1
               if ii < grid[0]:
                  vol.append([[]])
      hh.close()
   return ierr, na, aspecies, acharge, aposition, grid, origin, step, vol

if __name__ == "__main__":
   import sys
   import math
   sys.exit(main())

