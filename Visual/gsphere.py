#!/usr/bin/env python

#-------------------------------------------------------------------------------
#
#   gsphere.py
#   generates a sphere of G-vectors given lattice parameters and energy cutoff
#   written by Georgy Samsonidze (September 2009)
#
#-------------------------------------------------------------------------------

from common import eps9, inf9, boolean_input

def main(argv = None):
   if argv is None:
      argv = sys.argv
   if len(argv) != 3:
      sys.exit("\n   Usage: gsphere.py ifn ofn\n" +
                "   Description: reads system info from file ifn, generates the G-space\n" +
               "   and writes it to file ofn. ifn contains:\n" +
               "   * lattice vectors (in bohr),\n" +
               "   * energy cutoff (in Ry),\n" +
               "   * FFT grid (set to zero to determine automatically)\n" +
               "   * flag (true or false) for sorting G-vectors, and\n"+
               "   * a k-point at which to generate G-vectors (crystal coordinates).\n" +
               "   The number of G-vectors is roughly equal to\n" +
               "   N = V * E^(3/2) / (6 * pi^2),\n" +
               "   where V is the volume of the unit cell (in bohr^3)\n" +
               "   and E is the energy cutoff (in Ry).\n")
   ifn = argv[1]
   ofn = argv[2]
   ierr = 0
   try:
      hh = open(ifn, 'r')
   except:
      ierr = 1
   if ierr == 0:
      aa = []
      for ii in range(3):
         aa.append([])
         ss = hh.readline()
         tt = ss.split()
         for jj in range(3):
            aa[ii].append(float(tt[jj]))
      ss = hh.readline()
      tt = ss.split()
      ecut = float(tt[0])
      ss = hh.readline()
      tt = ss.split()
      nr = [int(tt[0]), int(tt[1]), int(tt[2])]
      ss = hh.readline()
      tt = ss.split()
      ierr, fsort = boolean_input(tt[0])
      ss = hh.readline()
      tt = ss.split()
      kpt = [float(tt[0]), float(tt[1]), float(tt[2])]
      hh.close()
   if ierr != 0:
      sys.exit("\n   Error: unable to read file %s\n" % ifn)

   al = math.sqrt(aa[0][0] * aa[0][0] + aa[0][1] * aa[0][1] + aa[0][2] * aa[0][2])
   for ii in range(3):
      for jj in range(3):
         aa[ii][jj] /= al
   vol = 0.0
   vol += aa[0][0] * (aa[1][1] * aa[2][2] - aa[2][1] * aa[1][2])
   vol -= aa[0][1] * (aa[1][0] * aa[2][2] - aa[2][0] * aa[1][2])
   vol += aa[0][2] * (aa[1][0] * aa[2][1] - aa[2][0] * aa[1][1])
   bl = 2.0 * math.pi / al
   bb = [[], [], []]
   bb[0].append((aa[1][1] * aa[2][2] - aa[1][2] * aa[2][1]) / vol)
   bb[0].append((aa[1][2] * aa[2][0] - aa[1][0] * aa[2][2]) / vol)
   bb[0].append((aa[1][0] * aa[2][1] - aa[1][1] * aa[2][0]) / vol)
   bb[1].append((aa[2][1] * aa[0][2] - aa[2][2] * aa[0][1]) / vol)
   bb[1].append((aa[2][2] * aa[0][0] - aa[2][0] * aa[0][2]) / vol)
   bb[1].append((aa[2][0] * aa[0][1] - aa[2][1] * aa[0][0]) / vol)
   bb[2].append((aa[0][1] * aa[1][2] - aa[0][2] * aa[1][1]) / vol)
   bb[2].append((aa[0][2] * aa[1][0] - aa[0][0] * aa[1][2]) / vol)
   bb[2].append((aa[0][0] * aa[1][1] - aa[0][1] * aa[1][0]) / vol)
   vol = abs(vol) * al * al * al
   gcutm = ecut / (bl * bl)

   try:
      hh = open(ofn, 'w')
   except:
      ierr = 1
   if ierr == 0:
      hh.write('\n')
      for ii in range(3):
         ss = ' a%1i =' % (ii + 1)
         for jj in range(3):
            ss += ' %13.9f' % (al * aa[ii][jj])
         ss += '\n'
         hh.write(ss)
      hh.write('\n')
      for ii in range(3):
         ss = ' b%1i =' % (ii + 1)
         for jj in range(3):
            ss += ' %13.9f' % (bl * bb[ii][jj])
         ss += '\n'
         hh.write(ss)
      hh.write('\n')
      ss = ' volume = %.9f\n' % vol
      hh.write(ss)
      hh.write('\n')
      ss = ' cutoff = %.9f\n' % ecut
      hh.write(ss)
      hh.write('\n')
   if ierr != 0:
      sys.exit("\n   Error: unable to write file %s\n" % ofn)

   if nr[0] < 1 or nr[1] < 1 or nr[2] < 1:
      for ii in range(3):
         nr[ii] = 2 * int(math.sqrt(gcutm) * math.sqrt(aa[ii][0] * aa[ii][0] + aa[ii][1] * aa[ii][1] + aa[ii][2] * aa[ii][2])) + 1
   nr6 = []
   for ii in range(3):
      nr6.append(nr[ii])
      while not check_nr_size(nr6[ii], 6):
         nr6[ii] += 1
   nr3 = []
   for ii in range(3):
      nr3.append(nr[ii])
      while not check_nr_size(nr3[ii], 3):
         nr3[ii] += 1

   if ierr == 0:
      ss = ' grid = ('
      for ii in range(3):
         ss += ' %i' % nr[ii]
      ss += ' ) -- minimal\n'
      hh.write(ss)
      ss = '        ('
      for ii in range(3):
         ss += ' %i' % nr6[ii]
      ss += ' ) -- factors 2, 3, 5, 7, 1*11, 1*13\n'
      hh.write(ss)
      ss = '        ('
      for ii in range(3):
         ss += ' %i' % nr3[ii]
      ss += ' ) -- factors 2, 3, 5\n'
      hh.write(ss)
      hh.write('\n')

   grange = []
   for idir in range(3):
      grange.append((nr[idir] - 1) // 2)

   if ierr == 0:
      ss = ' max g-vector indices = ('
      for idir in range(3):
         ss += " %i" % grange[idir]
      ss += ' )\n'
      hh.write(ss)
      hh.write('\n')
      ss = ' k-point ='
      for ii in range(3):
         ss += ' %13.9f' % kpt[ii]
      ss += '\n'
      hh.write(ss)
      hh.write('\n')

   ng = 0
   g_cryst = []
   g_cart = []
   e_kin = []
   for ii in range(-grange[0], grange[0] + 1):
      for jj in range(-grange[1], grange[1] + 1):
         for kk in range(-grange[2], grange[2] + 1):
            ifl = float(ii) + kpt[0]
            jfl = float(jj) + kpt[1]
            kfl = float(kk) + kpt[2]
            gx = ifl * bb[0][0] + jfl * bb[1][0] + kfl * bb[2][0]
            gy = ifl * bb[0][1] + jfl * bb[1][1] + kfl * bb[2][1]
            gz = ifl * bb[0][2] + jfl * bb[1][2] + kfl * bb[2][2]
            g2 = gx * gx + gy * gy + gz * gz
            if bl * bl * g2 < ecut:
               ng += 1
               g_cryst.append([ii, jj, kk])
               gx = float(ii) * bb[0][0] + float(jj) * bb[1][0] + float(kk) * bb[2][0]
               gy = float(ii) * bb[0][1] + float(jj) * bb[1][1] + float(kk) * bb[2][1]
               gz = float(ii) * bb[0][2] + float(jj) * bb[1][2] + float(kk) * bb[2][2]
               g_cart.append([bl * gx, bl * gy, bl * gz])
               e_kin.append(bl * bl * g2)
   g_sort = []
   for ii in range(ng):
      g_sort.append(ii)
   if fsort:
      e_sort = []
      for ii in range(ng):
         e_sort.append(e_kin[ii])
      e_sort.sort()
      for ii in range(ng):
         jj = e_sort.index(e_kin[ii])
         e_sort[jj] = inf9
         g_sort[jj] = ii

   if ierr == 0:
      ss = ' ng = %i\n' % ng
      hh.write(ss)
      hh.write('\n')
      hh.write('     n      g1    g2    g3          gx        gy        gz              ekin\n')
      hh.write(' ---------------------------------------------------------------------------\n')
      for ii in range(ng):
         ss = ' %5i  ' % (ii + 1)
         for jj in range(3):
            ss += ' %5i' % g_cryst[g_sort[ii]][jj]
         ss += '  '
         for jj in range(3):
            ss += ' %9.3f' % g_cart[g_sort[ii]][jj]
         ss += ' %17.9f' % e_kin[g_sort[ii]]
         ss += '\n'
         hh.write(ss)
      hh.write('\n')
      hh.close()

   return 0

def check_nr_size(nr, nf):
   factor = [2, 3, 5, 7, 11, 13]
   power = [0, 0, 0, 0, 0, 0]
   remainder = nr
   for ii in range(nf):
      jj = int(math.log(remainder) / math.log(factor[ii])) + 1
      for kk in range(jj):
         if remainder % factor[ii] == 0:
            remainder //= factor[ii]
            power[ii] += 1
   return remainder == 1 and power[4] <= 1 and power[5] <= 1

if __name__ == "__main__":
   import sys
   import math
   sys.exit(main())

