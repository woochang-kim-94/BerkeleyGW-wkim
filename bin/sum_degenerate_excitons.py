#!/usr/bin/env python

# DAS 3 April 2012
# Average the transition matrix elements of degenerate excitons.
# For use in testsuite, but possibly useful in other contexts too.
# Reads from eigenvalues.dat. No arguments. Writes to standard out.

from __future__ import print_function

f_in = open('eigenvalues.dat', 'r')

dipole_sum = 0.0
energy_prev = -1.0
multiplicity = 0

print("      energy    multiplicity    sum(|dipole|^2)")
for line in f_in:
    if line[:1] == '#':
        continue
    linesplit = line.split()
    energy = float(linesplit[0])
    dipole = float(linesplit[1])
    if abs(energy - energy_prev) < 1e-6:
        dipole_sum += dipole
        multiplicity += 1
    else:
        if multiplicity > 0:
            print("%12.6f  %14i  %17.6f"%(energy_prev, multiplicity, dipole_sum))
        dipole_sum = dipole
        energy_prev = energy
        multiplicity = 1

f_in.close()
