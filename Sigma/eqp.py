#!/usr/bin/env python
# eqp.py
#
# Usage: eqp.py edft|ecor|eqp0|eqp1|eqp0p|eqp1p sigma_hp.log eqp.dat [edft.dat]
#
# Extracts quasiparticle energies from sigma_hp.log file produced by Sigma
# Writes eqp.dat file that can be read in Epsilon/Sigma/BSE(absorption/inteqp)
#
# 1st parameter specifies whether to write edft|ecor|eqp0|eqp1|eqp0p|eqp1p energies
#
# OPTIONAL: if edft.dat is provided, k-dependent scissors operators are applied to
# the bands found in edft.dat but not present in sigma_hp.log. These are typically
# low-energy occupied and high-energy unoccupied bands far from the Fermi level for
# which quasiparticle corrections are not computed. edft.dat is extracted from WFN
# using MeanField/Utilities/scissors2eqp.x. THIS FEATURE IS HIGHLY EXPERIMENTAL!!!
#
# Georgy Samsonidze, UCB, August 2009

def main(argv = None):
   if argv is None:
      argv = sys.argv
   argc = len(argv)
   if argc != 4 and argc != 5:
      self = '/' + argv[0]
      self = self[len(self)-self[::-1].find('/'):]
      sys.exit("\n   Usage: %s edft|ecor|eqp0|eqp1|eqp0p|eqp1p sigma_hp.log eqp.dat\n                 [edft.dat]\n" % self)
   seqp = argv[1].lower()
   fni = argv[2]
   fno = argv[3]
   is_scissor = (argc == 5)
   if is_scissor:
      fni2 = argv[4]

   if seqp == 'edft':
      ieqp = 1
   elif seqp == 'ecor':
      ieqp = 2
   elif seqp == 'eqp0':
      ieqp = 8
   elif seqp == 'eqp1':
      ieqp = 9
   elif seqp == 'eqp0p':
      ieqp = 12
   elif seqp == 'eqp1p':
      ieqp = 13
   else:
      sys.exit("\n   Error: eqp type %s\n" % seqp)
 
   ierr = 0
   ifrq = 1
   try:
      h = open(fni, 'r')
      r = h.readlines()
      h.close()
   except:
      ierr = 1
   if ierr == 0:
      kpt = []
      kidx = []
      spin = []
      pos = []
      i = 0
      for s in r:
         if s.find("frequency_dependence") != -1:
            t = s.split()
            freq_dep = int(t[1])
            print("Frequency dependence = %s" % freq_dep)
            if freq_dep > 1:
               ifrq = 2
         if s.find("ik =") != -1:
            t = s.split()
            kpt.append([float(t[2]),float(t[3]),float(t[4])])
            kidx.append(int(t[7]))
            spin.append(int(t[10]))
            pos.append(i+2+ifrq)
         i += 1
      nks = len(kidx)
      ns = max(spin)
      nb = []
      band = []
      edft = []
      eqp = []
      minZnk = 1.0e6
      maxZnk = -1.0e6
      for i in range(nks):
         band.append([])
         edft.append([])
         eqp.append([])
         j = 0
         while 1:
            t = r[pos[i]+j*ifrq].split()
            if len(t) < 11:
               break
            band[i].append(int(t[0]))
            edft[i].append(float(t[1]))
            eqp[i].append(float(t[ieqp]))
            Znk = float(t[len(t) - 1])
            if Znk < minZnk:
                minZnk = Znk
            if Znk > maxZnk:
                maxZnk = Znk
            j += 1
         nb.append(j)
   else:
      sys.exit("\n   Error: could not open input file %s\n" % fni)
   ss = 'min Znk = %.6f\nmax Znk = %.6f' % (minZnk, maxZnk)
   print(ss)
   eps = 1.0e-3
   if minZnk < 0.0 - eps or maxZnk > 1.0 + eps:
       print('Check linearity of Sigma(E) using finite_difference_form = -3')

   if is_scissor:
      tol = 1.0e-5
      ierr = 0
      try:
         h = open(fni2, 'r')
         r = h.readlines()
         h.close()
      except:
         ierr = 1
      if ierr == 0:
         edft2 = []
         nbtot = int(r[0].split()[3]) / ns
         for i in range(nks / ns):
            t = r[i * (nbtot * ns + 1)].split()
            kpt2 = [float(t[0]),float(t[1]),float(t[2])]
            if abs(kpt[i][0] - kpt2[0]) > tol or abs(kpt[i][1] - kpt2[1]) > tol or abs(kpt[i][2] - kpt2[2]) > tol:
               sys.exit("\n   Error: k-point mismatch in %s\n" % fni2)
            for k in range(ns):
               edft2.append([])
               for j in range(nbtot):
                  t = r[i * (nbtot * ns + 1) + k * nbtot + j + 1].split()
                  edft2[i * ns + k].append(float(t[2]))
         eqp2 = []
         for i in range(nks):
            eqp2.append([])
            nbmin = band[i][0]
            nbmax = band[i][nb[i] - 1]
            demin = eqp[i][0] - edft[i][0]
            demax = eqp[i][nb[i] - 1] - edft[i][nb[i] - 1]
            for j in range(nbtot):
                if j < nbmin - 1:
                   eqp2[i].append(edft2[i][j] + demin)
                elif j < nbmax:
                   if abs(edft[i][j - nbmin + 1] - edft2[i][j]) > tol:
                      sys.exit("\n   Error: energy mismatch in %s\n" % fni2)
                   eqp2[i].append(eqp[i][j - nbmin + 1])
                else:
                   eqp2[i].append(edft2[i][j] + demax)
         nb2 = []
         for i in range(nks):
            nb2.append(nbtot)
         band2 = []
         for i in range(nks):
            band2.append([])
            for j in range(nbtot):
               band2[i].append(j + 1)
      else:
         sys.exit("\n   Error: could not open input file %s\n" % fni2)
   else:
      edft2 = edft
      eqp2 = eqp
      nb2 = nb
      band2 = band

   ierr = 0
   try:
      h = open(fno, 'w')
   except:
      ierr = 1
   if ierr == 0:
      for i in range(nks):
         if spin[i] == 1:
            s = '%13.9f%13.9f%13.9f%8i\n' % (kpt[i][0], kpt[i][1], kpt[i][2], ns * nb2[i])
            h.write(s)
         else:
            if kpt[i] != kpt[i-1]:
               sys.exit("\n   Error: k-point/spin order\n")
         for j in range(nb2[i]):
            s = '%8i%8i%15.9f%15.9f\n' % (spin[i], band2[i][j], edft2[i][j], eqp2[i][j])
            h.write(s)
      h.close()
   else:
      sys.exit("\n   Error: output file %s\n" % fno)

   return 0

if __name__ == "__main__":
   import sys
   import math
   import string
   sys.exit(main())
