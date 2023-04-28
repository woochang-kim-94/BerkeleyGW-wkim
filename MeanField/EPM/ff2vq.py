#!/usr/bin/env python
# ff2vq.py
#
# Reads EPM form factors from file form_factors.dat,
# fits them to the chosen functional form of the V(q) potential,
# writes potential coefficients to file v_of_q.dat,
# writes new form factors computed from V(q) to file vq2ff.dat.
#
# Georgy Samsonidze, UCB, 2011-02-04

from math import pi, exp, pow
from numpy import *
from scipy.optimize import leastsq

tol = 1.0e-6
bohr = 0.52917721092

def fitvq(material, a, V3S, V8S, V11S, V3A, V4A, V11A, hovq, hoff):
  #
  b = pow(2.0 * pi * bohr / a, 2)
  #
  v = lambda c, x: c[0] * (x - c[1]) / (c[2] * exp(c[3] * x) - 1.0)
  c0 = array([0.54, 2.19, 2.06, 0.48])
  e = lambda c, x, y: (v(c, x) - y)
  x = array([3.0 * b, 8.0 * b, 11.0 * b, 1.0e3 * b])
  y = array([V3S, V8S, V11S, 0.0])
  c, i = leastsq(e, c0, args = (x, y), maxfev = 10000)
  #
  if abs(V3A) + abs(V4A) + abs(V11A) < tol:
    #
    ca = c
    cb = array([0.0, 0.0, 0.0, 0.0])
    #
    s = '%-4s  %4.2f ' % (material, a)
    s += ' %6.3f' % v(c, 3.0 * b)
    s += ' %6.3f' % v(c, 8.0 * b)
    s += ' %6.3f' % v(c, 11.0 * b)
    s += ' %6.3f' % 0.0
    s += ' %6.3f' % 0.0
    s += ' %6.3f' % 0.0
    s += '\n'
    hoff.write(s)
    #
  else:
    #
    V4S = v(c, 4.0 * b)
    #
    va = lambda c, x: c[0] * (x - c[1]) / (exp(c[2] * (x + c[3])) + 1.0)
    c0 = array([1.0, 3.0, 1.0, 3.0])
    ea = lambda c, x, y: (va(c, x) - y)
    x = array([3.0 * b, 4.0 * b, 11.0 * b, 1.0e6 * b])
    y = array([V3S + V3A, V4S + V4A, V11S + V11A, 0.0])
    ca, i = leastsq(ea, c0, args = (x, y), maxfev = 10000)
    #
    vb = lambda c, x: c[0] * (x - c[1]) / (exp(c[2] * (x + c[3])) + 1.0)
    c0 = array([1.0, 3.0, 1.0, 3.0])
    eb = lambda c, x, y: (vb(c, x) - y)
    x = array([3.0 * b, 4.0 * b, 11.0 * b, 1.0e6 * b])
    y = array([V3S - V3A, V4S - V4A, V11S - V11A, 0.0])
    cb, i = leastsq(eb, c0, args = (x, y), maxfev = 10000)
    #
    s = '%-4s  %4.2f ' % (material, a)
    s += ' %6.3f' % ((va(ca, 3.0 * b) + vb(cb, 3.0 * b)) / 2.0)
    s += ' %6.3f' % ((va(ca, 8.0 * b) + vb(cb, 8.0 * b)) / 2.0)
    s += ' %6.3f' % ((va(ca, 11.0 * b) + vb(cb, 11.0 * b)) / 2.0)
    s += ' %6.3f' % ((va(ca, 3.0 * b) - vb(cb, 3.0 * b)) / 2.0)
    s += ' %6.3f' % ((va(ca, 4.0 * b) - vb(cb, 4.0 * b)) / 2.0)
    s += ' %6.3f' % ((va(ca, 11.0 * b) - vb(cb, 11.0 * b)) / 2.0)
    s += '\n'
    hoff.write(s)
    #
  s = '%-4s  %4.2f ' % (material, a)
  s += ' %9.6f %9.6f %9.6f %9.6f' % (ca[0], ca[1], ca[2], ca[3])
  s += ' %9.6f %9.6f %9.6f %9.6f' % (cb[0], cb[1], cb[2], cb[3])
  s += '\n'
  hovq.write(s)
  #

hiff = open('form_factors.dat', 'r')
hovq = open('v_of_q.dat', 'w')
hoff = open('vq2ff.dat', 'w')

s = hiff.readline()
hovq.write('      a      C_1^A     C_2^A     C_3^A     C_4^A     C_1^B     C_2^B     C_3^B     C_4^B\n')
hoff.write('      a      V3S    V8S    V11S   V3A    V4A    V11A \n')

while 1:
  s = hiff.readline()
  if len(s) < 2:
    break
  t = s.split()
  material = t[0]
  a = float(t[1])
  V3S = float(t[2])
  V8S = float(t[3])
  V11S = float(t[4])
  V3A = float(t[5])
  V4A = float(t[6])
  V11A  = float(t[7])
  fitvq(material, a, V3S, V8S, V11S, V3A, V4A, V11A, hovq, hoff)

hovq.write('\n')
hovq.write('elemental compound: V(q) = C_1 * (q^2 - C_2) / (C_3 * exp(C_4 * q^2) - 1)\n')
hovq.write('binary compound: V(q) = C_1 * (q^2 - C_2) / (exp(C_3 * (q^2 + C_4)) + 1)\n')
hovq.write('\n')
hovq.write('a in Angstrom, V in Rydberg, q in bohr^-1\n')
hoff.write('\n')
hoff.write('a in Angstrom, VG^2S & VG^2A in Rydberg, G^2 in 2*pi/a\n')

hiff.close()
hovq.close()
hoff.close()
