#!/usr/bin/env python
# irrbz.py
#
# Usage: irrbz.py irrbz prefix
#
# Extracts irreducible wedge from siesta output files prefix.KP 
# and prefix.EIG and generates prefix.KP.irrbz and prefix.EIG.irrbz.
# Reads in lattice vectors from prefix.XV and k-points from irrbz
# generated by kgrid.x.
#
# Georgy Samsonidze, UCB, December 2010

eps = 1.0e-5

def main(argv = None):
   if argv is None:
      argv = sys.argv
   argc = len(argv)
   self = "irrbz.py"
   if argv[0][-len(self):] != self:
      sys.exit("\n   Rename script to %s\n" % self)
   if argc != 3:
      sys.exit("\n   Usage: %s irrbz prefix\n" % self)
   irrbz = argv[1]
   prefix = argv[2]
   kpin = prefix + '.KP'
   kpout = prefix + '.KP.irrbz'
   eigin = prefix + '.EIG'
   eigout = prefix + '.EIG.irrbz'

   h = open(irrbz, 'r')
   s = h.readline()
   if s[0:16] != 'K_POINTS crystal':
      sys.exit("\n   Error: irrbz file %s\n" % irrbz)
   s = h.readline()
   nirr = int(s)
   kirr = []
   for i in xrange(nirr):
      s = h.readline()
      t = s.split()
      k1 = float(t[0])
      k2 = float(t[1])
      k3 = float(t[2])
      kirr.append([k1,k2,k3])
   h.close()

   h = open(prefix + '.XV', 'r')
   avec = []
   for i in xrange(3):
      s = h.readline()
      t = s.split()
      ax = float(t[0])
      ay = float(t[1])
      az = float(t[2])
      avec.append([ax,ay,az])
   h.close()
   bvec = [[], [], []]
   celvol = 0.0
   celvol += avec[0][0] * (avec[1][1] * avec[2][2] - avec[2][1] * avec[1][2])
   celvol -= avec[0][1] * (avec[1][0] * avec[2][2] - avec[2][0] * avec[1][2])
   celvol += avec[0][2] * (avec[1][0] * avec[2][1] - avec[2][0] * avec[1][1])
   bvec[0].append(2.0 * math.pi * (avec[1][1] * avec[2][2] - avec[1][2] * avec[2][1]) / celvol)
   bvec[0].append(2.0 * math.pi * (avec[1][2] * avec[2][0] - avec[1][0] * avec[2][2]) / celvol)
   bvec[0].append(2.0 * math.pi * (avec[1][0] * avec[2][1] - avec[1][1] * avec[2][0]) / celvol)
   bvec[1].append(2.0 * math.pi * (avec[2][1] * avec[0][2] - avec[2][2] * avec[0][1]) / celvol)
   bvec[1].append(2.0 * math.pi * (avec[2][2] * avec[0][0] - avec[2][0] * avec[0][2]) / celvol)
   bvec[1].append(2.0 * math.pi * (avec[2][0] * avec[0][1] - avec[2][1] * avec[0][0]) / celvol)
   bvec[2].append(2.0 * math.pi * (avec[0][1] * avec[1][2] - avec[0][2] * avec[1][1]) / celvol)
   bvec[2].append(2.0 * math.pi * (avec[0][2] * avec[1][0] - avec[0][0] * avec[1][2]) / celvol)
   bvec[2].append(2.0 * math.pi * (avec[0][0] * avec[1][1] - avec[0][1] * avec[1][0]) / celvol)

   h = open(kpin, 'r')
   s = h.readline()
   nful = int(s)
   kful = []
   for i in xrange(nful):
      s = h.readline()
      t = s.split()
      j = int(t[0])
      kx = float(t[1])
      ky = float(t[2])
      kz = float(t[3])
      wk = float(t[4])
      if j != i + 1:
         sys.exit("\n   Error: KP file %s\n" % kpin)
      kful.append([kx,ky,kz,wk])
   h.close()

   kidx = []
   for i in xrange(nirr):
      for j in xrange(nful):
         kx = kful[j][0]
         ky = kful[j][1]
         kz = kful[j][2]
         k1 = (kx * avec[0][0] + ky * avec[0][1] + kz * avec[0][2]) / (2.0 * math.pi)
         k2 = (kx * avec[1][0] + ky * avec[1][1] + kz * avec[1][2]) / (2.0 * math.pi)
         k3 = (kx * avec[2][0] + ky * avec[2][1] + kz * avec[2][2]) / (2.0 * math.pi)
         d1 = kirr[i][0] - k1
         d2 = kirr[i][1] - k2
         d3 = kirr[i][2] - k3
         l1 = abs(d1 - round(d1)) + abs(d2 - round(d2)) + abs(d3 - round(d3))
         d1 = kirr[i][0] + k1
         d2 = kirr[i][1] + k2
         d3 = kirr[i][2] + k3
         l2 = abs(d1 - round(d1)) + abs(d2 - round(d2)) + abs(d3 - round(d3))
         if l1 < eps or l2 < eps:
            kidx.append(j)
            break
   if len(kidx) != nirr:
      sys.exit("\n   Error: cannot match k-points\n")

   kout = []
   for i in xrange(nirr):
      k1 = kirr[i][0]
      k2 = kirr[i][1]
      k3 = kirr[i][2]
      kx = (k1 * bvec[0][0] + k2 * bvec[1][0] + k3 * bvec[2][0])
      ky = (k1 * bvec[0][1] + k2 * bvec[1][1] + k3 * bvec[2][1])
      kz = (k1 * bvec[0][2] + k2 * bvec[1][2] + k3 * bvec[2][2])
      wk = kful[kidx[i]][3]
      kout.append([kx,ky,kz,wk])

   h = open(kpout, 'w')
   s = '%6i\n' % nirr
   h.write(s)
   for i in xrange(nirr):
      s = '%6i%12.6f%12.6f%12.6f%15.6f\n' % (i + 1, kout[i][0], kout[i][1], kout[i][2], kout[i][3])
      h.write(s)
   h.close()

   h = open(eigin, 'r')
   s = h.readline()
   efermi = float(s)
   s = h.readline()
   t = s.split()
   nb = int(t[0])
   ns = int(t[1])
   nk = int(t[2])
   if nk != nful:
      sys.exit("\n   Error: EIG file %s\n" % eigin)
   nline = (nb * ns) / 10
   if 10 * nline < nb * ns:
      nline += 1
   en = []
   for i in xrange(nful):
      en.append([])
      for j in xrange(nline):
         s = h.readline()
         if j == 0:
            t = s.split()
            if int(t[0]) != i + 1:
               sys.exit("\n   Error: EIG file %s\n" % eigin)
         en[i].append(s)
   h.close()

   h = open(eigout, 'w')
   s = '%14.4f\n' % efermi
   h.write(s)
   s = '%6i%6i%6i\n' % (nb, ns, nirr)
   h.write(s)
   for i in xrange(nirr):
      for j in xrange(nline):
         if j == 0:
            s = '%5i' % (i + 1)
            s += en[kidx[i]][j][5:]
            h.write(s)
         else:
            h.write(en[kidx[i]][j])
   h.close()

   return 0

if __name__ == "__main__":
   import sys
   import math
   sys.exit(main())
