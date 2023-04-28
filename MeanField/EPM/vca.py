#!/usr/bin/env python
# vca.py
#
# Reads V(q) potential coefficients from file v_of_q.dat,
# employs the virtual crystal approximation to compute
# hybrid form factors, writes them to file vca_ff.dat.
# The potential mixing is controlled by identifiers
# host_material and doping_level.
#
# Georgy Samsonidze, UCB, 2011-02-06

from sys import exit
from math import pi, exp, pow

tol = 1.0e-6
bohr = 0.52917721092

host_material = 'Si'
doping_level = 0.05

def v_elemental_of_q(q, c):
  return c[0] * (q - c[1]) / (c[2] * exp(c[3] * q) - 1.0)

def v_binary_of_q(q, c):
  return c[0] * (q - c[1]) / (exp(c[2] * (q + c[3])) + 1.0)

def v_of_b(b, ca, cb):
  if abs(cb[0]) + abs(cb[1]) + abs(cb[2]) + abs(cb[3]) < tol:
    V3S = v_elemental_of_q(3.0 * b, ca)
    V8S = v_elemental_of_q(8.0 * b, ca)
    V11S = v_elemental_of_q(11.0 * b, ca)
    V3A = 0.0
    V4A = 0.0
    V11A = 0.0
  else:
    V3S = (v_binary_of_q(3.0 * b, ca) + v_binary_of_q(3.0 * b, cb)) / 2.0
    V8S = (v_binary_of_q(8.0 * b, ca) + v_binary_of_q(8.0 * b, cb)) / 2.0
    V11S = (v_binary_of_q(11.0 * b, ca) + v_binary_of_q(11.0 * b, cb)) / 2.0
    V3A = (v_binary_of_q(3.0 * b, ca) - v_binary_of_q(3.0 * b, cb)) / 2.0
    V4A = (v_binary_of_q(4.0 * b, ca) - v_binary_of_q(4.0 * b, cb)) / 2.0
    V11A = (v_binary_of_q(11.0 * b, ca) - v_binary_of_q(11.0 * b, cb)) / 2.0
  return [V3S, V8S, V11S], [V3A, V4A, V11A]

def mixff(hybrid_material, doping_level, host_a, host_ca, host_cb, dopant_a, dopant_ca, dopant_cb, da, dVS, dVA, hoff):
  #
  hybrid_a = pow((1 - doping_level) * pow(host_a, 3) + doping_level * pow(dopant_a, 3), 1.0/3.0)
  hybrid_b = pow(2.0 * pi * bohr / hybrid_a, 2)
  host_VS, host_VA = v_of_b(hybrid_b, host_ca, host_cb)
  dopant_VS, dopant_VA = v_of_b(hybrid_b, dopant_ca, dopant_cb)
  hybrid_VS = []
  hybrid_VA = []
  for i in xrange(3):
    hybrid_VS.append((1 - doping_level) * host_VS[i] + doping_level * dopant_VS[i])
    hybrid_VA.append((1 - doping_level) * host_VA[i] + doping_level * dopant_VA[i])
    #
  s = '%-9s  %6.4f ' % (hybrid_material, hybrid_a)
  s += ' %8.5f %8.5f %8.5f' % (hybrid_VS[0], hybrid_VS[1], hybrid_VS[2])
  s += ' %8.5f %8.5f %8.5f' % (hybrid_VA[0], hybrid_VA[1], hybrid_VA[2])
  s += '\n'
  hoff.write(s)
  #
  host_b = pow(2.0 * pi * bohr / host_a, 2)
  pure_VS, pure_VA = v_of_b(host_b, host_ca, host_cb)
  da.append(hybrid_a - host_a)
  for i in xrange(3):
    dVS[i].append(hybrid_VS[i] - pure_VS[i])
    dVA[i].append(hybrid_VA[i] - pure_VA[i])
    #

hivq = open('v_of_q.dat', 'r')
hoff = open('vca_ff.dat', 'w')

s = hivq.readline()
hoff.write('           a        V3S      V8S      V11S     V3A      V4A      V11A\n')

i = 0
host_index = -1
material = []
a = []
ca = []
cb = []
while 1:
  s = hivq.readline()
  if len(s) < 2:
    break
  t = s.split()
  if t[0] == host_material:
    host_index = i
  material.append(t[0])
  a.append(float(t[1]))
  ca.append([])
  ca[i].append(float(t[2]))
  ca[i].append(float(t[3]))
  ca[i].append(float(t[4]))
  ca[i].append(float(t[5]))
  cb.append([])
  cb[i].append(float(t[6]))
  cb[i].append(float(t[7]))
  cb[i].append(float(t[8]))
  cb[i].append(float(t[9]))
  i += 1
if host_index == -1:
  sys.exit("\n   Error: invalid host material\n")

da = []
dVS = [[], [], []]
dVA = [[], [], []]
for dopant_index in xrange(len(material)):
  if dopant_index == host_index:
    continue
  hybrid_material = '%s-%s' % (material[host_index], material[dopant_index])
  mixff(hybrid_material, doping_level, a[host_index], ca[host_index], cb[host_index], a[dopant_index], ca[dopant_index], cb[dopant_index], da, dVS, dVA, hoff)
s = 'minimum   %7.4f ' % min(da)
s += ' %8.5f %8.5f %8.5f' % (min(dVS[0]), min(dVS[1]), min(dVS[2]))
s += ' %8.5f %8.5f %8.5f' % (min(dVA[0]), min(dVA[1]), min(dVA[2]))
s += '\n'
s += 'maximum   %7.4f ' % max(da)
s += ' %8.5f %8.5f %8.5f' % (max(dVS[0]), max(dVS[1]), max(dVS[2]))
s += ' %8.5f %8.5f %8.5f' % (max(dVA[0]), max(dVA[1]), max(dVA[2]))
s += '\n\n'
s += 'doping level = %.4f, minimum/maximum = deviation from pure material' % doping_level
s += '\n'
hoff.write(s)

hoff.write('\n')
hoff.write('a in Angstrom, VG^2S & VG^2A in Rydberg, G^2 in 2*pi/a\n')

hivq.close()
hoff.close()
