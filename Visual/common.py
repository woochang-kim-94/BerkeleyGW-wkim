#!/usr/bin/env python

#-------------------------------------------------------------------------------
#
#   common.py
#   library of common data structures
#   written by Georgy Samsonidze (October 2008)
#
#-------------------------------------------------------------------------------

def main(argv = None):
   if argv is None:
      argv = sys.argv
   sys.exit("\n   Module common.py contains common data structures\n")

eps9 = 1.0e-9
inf9 = 1.0e+9

bohr = 0.52917721092
rydberg = 13.60569253
hartree = 27.21138505

format_mat2index = {'bohr': 0, 'angstrom': 1, 'latconst': 2, 'latvec': 3}
format_index2mat = {0: 'bohr', 1: 'angstrom', 2: 'latconst', 3: 'latvec'}
format_espresso2mat = {'bohr': 'bohr', 'angstrom': 'angstrom', 'alat': 'latconst', 'crystal': 'latvec'}
format_mat2espresso = {'bohr': 'bohr', 'angstrom': 'angstrom', 'latconst': 'alat', 'latvec': 'crystal'}
format_siesta2mat = {'bohr': 'bohr', 'notscaledcartesianbohr': 'bohr', 'ang': 'angstrom', 'notscaledcartesianang': 'angstrom', 'scaledcartesian': 'latconst', 'fractional': 'latvec', 'scaledbylatticevectors': 'latvec', 'nm': 'angstrom', 'cm': 'angstrom', 'm': 'angstrom'}
factor_siesta2mat = {'bohr': 1.0, 'ang': 1.0, 'nm': 1.0e1, 'cm': 1.0e8, 'm': 1.0e10}
format_mat2siesta = {'bohr': 'bohr', 'angstrom': 'ang', 'latconst': 'scaledcartesian', 'latvec': 'fractional'}
format_tbpw2mat = {'bohr': 'bohr', 'angstrom': 'angstrom', 'scaledcartesian': 'latconst', 'scaledbylatticevectors': 'latvec'}
format_mat2tbpw = {'bohr': 'bohr', 'angstrom': 'angstrom', 'latconst': 'scaledcartesian', 'latvec': 'scaledbylatticevectors'}

atomsize = 0.20
bondradius = 0.15
bondtolerance = 0.45
minimumbondingdistance = 0.40
strictvalence = 0

vertex = [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
edge = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]]
face = [[0, 1, 2], [0, 3, 2], [4, 5, 6], [4, 7, 6], [0, 1, 5], [0, 4, 5], [1, 2, 6], [1, 5, 6], [2, 3, 7], [2, 6, 7], [3, 0, 4], [3, 7, 4]]

periodic_table = [
{'number':   0, 'symbol': 'X',   'name': 'Dummy',         'nvelec':  0, 'valence':  0, 'period':  0, 'group':  0, 'rcov': 0.00, 'rvdw': 0.00, 'mass':   0.000, 'color': [0.07, 0.50, 0.70]},
{'number':   1, 'symbol': 'H',   'name': 'Hydrogen',      'nvelec':  1, 'valence':  1, 'period':  1, 'group':  1, 'rcov': 0.31, 'rvdw': 1.10, 'mass':   1.008, 'color': [0.75, 0.75, 0.75]},
{'number':   2, 'symbol': 'He',  'name': 'Helium',        'nvelec':  2, 'valence':  0, 'period':  1, 'group': 18, 'rcov': 0.28, 'rvdw': 1.40, 'mass':   4.003, 'color': [0.85, 1.00, 1.00]},
{'number':   3, 'symbol': 'Li',  'name': 'Lithium',       'nvelec':  1, 'valence':  1, 'period':  2, 'group':  1, 'rcov': 1.28, 'rvdw': 1.81, 'mass':   6.941, 'color': [0.80, 0.50, 1.00]},
{'number':   4, 'symbol': 'Be',  'name': 'Beryllium',     'nvelec':  2, 'valence':  2, 'period':  2, 'group':  2, 'rcov': 0.96, 'rvdw': 1.53, 'mass':   9.012, 'color': [0.76, 1.00, 0.00]},
{'number':   5, 'symbol': 'B',   'name': 'Boron',         'nvelec':  3, 'valence':  4, 'period':  2, 'group': 13, 'rcov': 0.84, 'rvdw': 1.92, 'mass':  10.811, 'color': [1.00, 0.71, 0.71]},
{'number':   6, 'symbol': 'C',   'name': 'Carbon',        'nvelec':  4, 'valence':  4, 'period':  2, 'group': 14, 'rcov': 0.76, 'rvdw': 1.70, 'mass':  12.011, 'color': [0.40, 0.40, 0.40]},
{'number':   7, 'symbol': 'N',   'name': 'Nitrogen',      'nvelec':  5, 'valence':  4, 'period':  2, 'group': 15, 'rcov': 0.71, 'rvdw': 1.55, 'mass':  14.007, 'color': [0.05, 0.05, 1.00]},
{'number':   8, 'symbol': 'O',   'name': 'Oxygen',        'nvelec':  6, 'valence':  2, 'period':  2, 'group': 16, 'rcov': 0.66, 'rvdw': 1.52, 'mass':  15.999, 'color': [1.00, 0.05, 0.05]},
{'number':   9, 'symbol': 'F',   'name': 'Fluorine',      'nvelec':  7, 'valence':  1, 'period':  2, 'group': 17, 'rcov': 0.57, 'rvdw': 1.47, 'mass':  18.998, 'color': [0.50, 0.70, 1.00]},
{'number':  10, 'symbol': 'Ne',  'name': 'Neon',          'nvelec':  8, 'valence':  0, 'period':  2, 'group': 18, 'rcov': 0.58, 'rvdw': 1.54, 'mass':  20.180, 'color': [0.70, 0.89, 0.96]},
{'number':  11, 'symbol': 'Na',  'name': 'Sodium',        'nvelec':  1, 'valence':  1, 'period':  3, 'group':  1, 'rcov': 1.66, 'rvdw': 2.27, 'mass':  22.990, 'color': [0.67, 0.36, 0.95]},
{'number':  12, 'symbol': 'Mg',  'name': 'Magnesium',     'nvelec':  2, 'valence':  2, 'period':  3, 'group':  2, 'rcov': 1.41, 'rvdw': 1.73, 'mass':  24.305, 'color': [0.54, 1.00, 0.00]},
{'number':  13, 'symbol': 'Al',  'name': 'Aluminium',     'nvelec':  3, 'valence':  6, 'period':  3, 'group': 13, 'rcov': 1.21, 'rvdw': 1.84, 'mass':  26.982, 'color': [0.75, 0.65, 0.65]},
{'number':  14, 'symbol': 'Si',  'name': 'Silicon',       'nvelec':  4, 'valence':  6, 'period':  3, 'group': 14, 'rcov': 1.11, 'rvdw': 2.10, 'mass':  28.085, 'color': [0.50, 0.60, 0.60]},
{'number':  15, 'symbol': 'P',   'name': 'Phosphorus',    'nvelec':  5, 'valence':  6, 'period':  3, 'group': 15, 'rcov': 1.07, 'rvdw': 1.80, 'mass':  30.974, 'color': [1.00, 0.50, 0.00]},
{'number':  16, 'symbol': 'S',   'name': 'Sulfur',        'nvelec':  6, 'valence':  6, 'period':  3, 'group': 16, 'rcov': 1.05, 'rvdw': 1.80, 'mass':  32.065, 'color': [0.70, 0.70, 0.00]},
{'number':  17, 'symbol': 'Cl',  'name': 'Chlorine',      'nvelec':  7, 'valence':  1, 'period':  3, 'group': 17, 'rcov': 1.02, 'rvdw': 1.75, 'mass':  35.453, 'color': [0.12, 0.94, 0.12]},
{'number':  18, 'symbol': 'Ar',  'name': 'Argon',         'nvelec':  8, 'valence':  0, 'period':  3, 'group': 18, 'rcov': 1.06, 'rvdw': 1.88, 'mass':  39.948, 'color': [0.50, 0.82, 0.89]},
{'number':  19, 'symbol': 'K',   'name': 'Potassium',     'nvelec':  1, 'valence':  1, 'period':  4, 'group':  1, 'rcov': 2.03, 'rvdw': 2.75, 'mass':  39.098, 'color': [0.56, 0.25, 0.83]},
{'number':  20, 'symbol': 'Ca',  'name': 'Calcium',       'nvelec':  2, 'valence':  2, 'period':  4, 'group':  2, 'rcov': 1.76, 'rvdw': 2.31, 'mass':  40.078, 'color': [0.24, 1.00, 0.00]},
{'number':  21, 'symbol': 'Sc',  'name': 'Scandium',      'nvelec':  3, 'valence':  6, 'period':  4, 'group':  3, 'rcov': 1.70, 'rvdw': 2.30, 'mass':  44.956, 'color': [0.90, 0.90, 0.90]},
{'number':  22, 'symbol': 'Ti',  'name': 'Titanium',      'nvelec':  4, 'valence':  6, 'period':  4, 'group':  4, 'rcov': 1.60, 'rvdw': 2.15, 'mass':  47.867, 'color': [0.75, 0.76, 0.78]},
{'number':  23, 'symbol': 'V',   'name': 'Vanadium',      'nvelec':  5, 'valence':  6, 'period':  4, 'group':  5, 'rcov': 1.53, 'rvdw': 2.05, 'mass':  50.941, 'color': [0.65, 0.65, 0.67]},
{'number':  24, 'symbol': 'Cr',  'name': 'Chromium',      'nvelec':  6, 'valence':  6, 'period':  4, 'group':  6, 'rcov': 1.39, 'rvdw': 2.05, 'mass':  51.996, 'color': [0.54, 0.60, 0.78]},
{'number':  25, 'symbol': 'Mn',  'name': 'Manganese',     'nvelec':  7, 'valence':  8, 'period':  4, 'group':  7, 'rcov': 1.39, 'rvdw': 2.05, 'mass':  54.938, 'color': [0.61, 0.48, 0.78]},
{'number':  26, 'symbol': 'Fe',  'name': 'Iron',          'nvelec':  8, 'valence':  6, 'period':  4, 'group':  8, 'rcov': 1.32, 'rvdw': 2.05, 'mass':  55.845, 'color': [0.88, 0.40, 0.20]},
{'number':  27, 'symbol': 'Co',  'name': 'Cobalt',        'nvelec':  9, 'valence':  6, 'period':  4, 'group':  9, 'rcov': 1.26, 'rvdw': 2.00, 'mass':  58.933, 'color': [0.94, 0.56, 0.63]},
{'number':  28, 'symbol': 'Ni',  'name': 'Nickel',        'nvelec': 10, 'valence':  6, 'period':  4, 'group': 10, 'rcov': 1.24, 'rvdw': 2.00, 'mass':  58.693, 'color': [0.31, 0.82, 0.31]},
{'number':  29, 'symbol': 'Cu',  'name': 'Copper',        'nvelec': 11, 'valence':  6, 'period':  4, 'group': 11, 'rcov': 1.32, 'rvdw': 2.00, 'mass':  63.546, 'color': [0.78, 0.50, 0.20]},
{'number':  30, 'symbol': 'Zn',  'name': 'Zinc',          'nvelec': 12, 'valence':  6, 'period':  4, 'group': 12, 'rcov': 1.22, 'rvdw': 2.10, 'mass':  65.380, 'color': [0.49, 0.50, 0.69]},
{'number':  31, 'symbol': 'Ga',  'name': 'Gallium',       'nvelec':  3, 'valence':  3, 'period':  4, 'group': 13, 'rcov': 1.22, 'rvdw': 1.87, 'mass':  69.723, 'color': [0.76, 0.56, 0.56]},
{'number':  32, 'symbol': 'Ge',  'name': 'Germanium',     'nvelec':  4, 'valence':  4, 'period':  4, 'group': 14, 'rcov': 1.20, 'rvdw': 2.11, 'mass':  72.640, 'color': [0.40, 0.56, 0.56]},
{'number':  33, 'symbol': 'As',  'name': 'Arsenic',       'nvelec':  5, 'valence':  3, 'period':  4, 'group': 15, 'rcov': 1.19, 'rvdw': 1.85, 'mass':  74.922, 'color': [0.74, 0.50, 0.89]},
{'number':  34, 'symbol': 'Se',  'name': 'Selenium',      'nvelec':  6, 'valence':  2, 'period':  4, 'group': 16, 'rcov': 1.20, 'rvdw': 1.90, 'mass':  78.960, 'color': [1.00, 0.63, 0.00]},
{'number':  35, 'symbol': 'Br',  'name': 'Bromine',       'nvelec':  7, 'valence':  1, 'period':  4, 'group': 17, 'rcov': 1.20, 'rvdw': 1.83, 'mass':  79.904, 'color': [0.65, 0.16, 0.16]},
{'number':  36, 'symbol': 'Kr',  'name': 'Krypton',       'nvelec':  8, 'valence':  0, 'period':  4, 'group': 18, 'rcov': 1.16, 'rvdw': 2.02, 'mass':  83.798, 'color': [0.36, 0.72, 0.82]},
{'number':  37, 'symbol': 'Rb',  'name': 'Rubidium',      'nvelec':  1, 'valence':  1, 'period':  5, 'group':  1, 'rcov': 2.20, 'rvdw': 3.03, 'mass':  85.468, 'color': [0.44, 0.18, 0.69]},
{'number':  38, 'symbol': 'Sr',  'name': 'Strontium',     'nvelec':  2, 'valence':  2, 'period':  5, 'group':  2, 'rcov': 1.95, 'rvdw': 2.49, 'mass':  87.620, 'color': [0.00, 1.00, 0.00]},
{'number':  39, 'symbol': 'Y',   'name': 'Yttrium',       'nvelec':  3, 'valence':  6, 'period':  5, 'group':  3, 'rcov': 1.90, 'rvdw': 2.40, 'mass':  88.906, 'color': [0.58, 1.00, 1.00]},
{'number':  40, 'symbol': 'Zr',  'name': 'Zirconium',     'nvelec':  4, 'valence':  6, 'period':  5, 'group':  4, 'rcov': 1.75, 'rvdw': 2.30, 'mass':  91.224, 'color': [0.58, 0.88, 0.88]},
{'number':  41, 'symbol': 'Nb',  'name': 'Niobium',       'nvelec':  5, 'valence':  6, 'period':  5, 'group':  5, 'rcov': 1.64, 'rvdw': 2.15, 'mass':  92.906, 'color': [0.45, 0.76, 0.79]},
{'number':  42, 'symbol': 'Mo',  'name': 'Molybdenum',    'nvelec':  6, 'valence':  6, 'period':  5, 'group':  6, 'rcov': 1.54, 'rvdw': 2.10, 'mass':  95.960, 'color': [0.33, 0.71, 0.71]},
{'number':  43, 'symbol': 'Tc',  'name': 'Technetium',    'nvelec':  7, 'valence':  6, 'period':  5, 'group':  7, 'rcov': 1.47, 'rvdw': 2.05, 'mass':  98.000, 'color': [0.23, 0.62, 0.62]},
{'number':  44, 'symbol': 'Ru',  'name': 'Ruthenium',     'nvelec':  8, 'valence':  6, 'period':  5, 'group':  8, 'rcov': 1.46, 'rvdw': 2.05, 'mass': 101.070, 'color': [0.14, 0.56, 0.56]},
{'number':  45, 'symbol': 'Rh',  'name': 'Rhodium',       'nvelec':  9, 'valence':  6, 'period':  5, 'group':  9, 'rcov': 1.42, 'rvdw': 2.00, 'mass': 102.906, 'color': [0.04, 0.49, 0.55]},
{'number':  46, 'symbol': 'Pd',  'name': 'Palladium',     'nvelec': 10, 'valence':  6, 'period':  5, 'group': 10, 'rcov': 1.39, 'rvdw': 2.05, 'mass': 106.420, 'color': [0.00, 0.41, 0.52]},
{'number':  47, 'symbol': 'Ag',  'name': 'Silver',        'nvelec': 11, 'valence':  6, 'period':  5, 'group': 11, 'rcov': 1.45, 'rvdw': 2.10, 'mass': 107.868, 'color': [0.88, 0.88, 1.00]},
{'number':  48, 'symbol': 'Cd',  'name': 'Cadmium',       'nvelec': 12, 'valence':  6, 'period':  5, 'group': 12, 'rcov': 1.44, 'rvdw': 2.20, 'mass': 112.411, 'color': [1.00, 0.85, 0.56]},
{'number':  49, 'symbol': 'In',  'name': 'Indium',        'nvelec':  3, 'valence':  3, 'period':  5, 'group': 13, 'rcov': 1.42, 'rvdw': 2.20, 'mass': 114.818, 'color': [0.65, 0.46, 0.45]},
{'number':  50, 'symbol': 'Sn',  'name': 'Tin',           'nvelec':  4, 'valence':  4, 'period':  5, 'group': 14, 'rcov': 1.39, 'rvdw': 1.93, 'mass': 118.701, 'color': [0.40, 0.50, 0.50]},
{'number':  51, 'symbol': 'Sb',  'name': 'Antimony',      'nvelec':  5, 'valence':  3, 'period':  5, 'group': 15, 'rcov': 1.39, 'rvdw': 2.17, 'mass': 121.760, 'color': [0.62, 0.39, 0.71]},
{'number':  52, 'symbol': 'Te',  'name': 'Tellurium',     'nvelec':  6, 'valence':  2, 'period':  5, 'group': 16, 'rcov': 1.38, 'rvdw': 2.06, 'mass': 127.600, 'color': [0.83, 0.48, 0.00]},
{'number':  53, 'symbol': 'I',   'name': 'Iodine',        'nvelec':  7, 'valence':  1, 'period':  5, 'group': 17, 'rcov': 1.39, 'rvdw': 1.98, 'mass': 126.904, 'color': [0.58, 0.00, 0.58]},
{'number':  54, 'symbol': 'Xe',  'name': 'Xenon',         'nvelec':  8, 'valence':  0, 'period':  5, 'group': 18, 'rcov': 1.40, 'rvdw': 2.16, 'mass': 131.293, 'color': [0.26, 0.62, 0.69]},
{'number':  55, 'symbol': 'Cs',  'name': 'Caesium',       'nvelec':  1, 'valence':  1, 'period':  6, 'group':  1, 'rcov': 2.44, 'rvdw': 3.43, 'mass': 132.905, 'color': [0.34, 0.09, 0.56]},
{'number':  56, 'symbol': 'Ba',  'name': 'Barium',        'nvelec':  2, 'valence':  2, 'period':  6, 'group':  2, 'rcov': 2.15, 'rvdw': 2.68, 'mass': 137.327, 'color': [0.00, 0.79, 0.00]},
{'number':  57, 'symbol': 'La',  'name': 'Lanthanum',     'nvelec':  3, 'valence': 12, 'period':  9, 'group':  3, 'rcov': 2.07, 'rvdw': 2.50, 'mass': 138.905, 'color': [0.44, 0.83, 1.00]},
{'number':  58, 'symbol': 'Ce',  'name': 'Cerium',        'nvelec':  4, 'valence':  6, 'period':  9, 'group':  4, 'rcov': 2.04, 'rvdw': 2.48, 'mass': 140.116, 'color': [1.00, 1.00, 0.78]},
{'number':  59, 'symbol': 'Pr',  'name': 'Praseodymium',  'nvelec':  5, 'valence':  6, 'period':  9, 'group':  5, 'rcov': 2.03, 'rvdw': 2.47, 'mass': 140.908, 'color': [0.85, 1.00, 0.78]},
{'number':  60, 'symbol': 'Nd',  'name': 'Neodymium',     'nvelec':  6, 'valence':  6, 'period':  9, 'group':  6, 'rcov': 2.01, 'rvdw': 2.45, 'mass': 144.240, 'color': [0.78, 1.00, 0.78]},
{'number':  61, 'symbol': 'Pm',  'name': 'Promethium',    'nvelec':  7, 'valence':  6, 'period':  9, 'group':  7, 'rcov': 1.99, 'rvdw': 2.43, 'mass': 145.000, 'color': [0.64, 1.00, 0.78]},
{'number':  62, 'symbol': 'Sm',  'name': 'Samarium',      'nvelec':  8, 'valence':  6, 'period':  9, 'group':  8, 'rcov': 1.98, 'rvdw': 2.42, 'mass': 150.360, 'color': [0.56, 1.00, 0.78]},
{'number':  63, 'symbol': 'Eu',  'name': 'Europium',      'nvelec':  9, 'valence':  6, 'period':  9, 'group':  9, 'rcov': 1.98, 'rvdw': 2.40, 'mass': 151.964, 'color': [0.38, 1.00, 0.78]},
{'number':  64, 'symbol': 'Gd',  'name': 'Gadolinium',    'nvelec': 10, 'valence':  6, 'period':  9, 'group': 10, 'rcov': 1.96, 'rvdw': 2.38, 'mass': 157.250, 'color': [0.27, 1.00, 0.78]},
{'number':  65, 'symbol': 'Tb',  'name': 'Terbium',       'nvelec': 11, 'valence':  6, 'period':  9, 'group': 11, 'rcov': 1.94, 'rvdw': 2.37, 'mass': 158.925, 'color': [0.19, 1.00, 0.78]},
{'number':  66, 'symbol': 'Dy',  'name': 'Dysprosium',    'nvelec': 12, 'valence':  6, 'period':  9, 'group': 12, 'rcov': 1.92, 'rvdw': 2.35, 'mass': 162.500, 'color': [0.12, 1.00, 0.78]},
{'number':  67, 'symbol': 'Ho',  'name': 'Holmium',       'nvelec': 13, 'valence':  6, 'period':  9, 'group': 13, 'rcov': 1.92, 'rvdw': 2.33, 'mass': 164.930, 'color': [0.00, 1.00, 0.61]},
{'number':  68, 'symbol': 'Er',  'name': 'Erbium',        'nvelec': 14, 'valence':  6, 'period':  9, 'group': 14, 'rcov': 1.89, 'rvdw': 2.32, 'mass': 167.259, 'color': [0.00, 0.90, 0.46]},
{'number':  69, 'symbol': 'Tm',  'name': 'Thulium',       'nvelec': 15, 'valence':  6, 'period':  9, 'group': 15, 'rcov': 1.90, 'rvdw': 2.30, 'mass': 168.934, 'color': [0.00, 0.83, 0.32]},
{'number':  70, 'symbol': 'Yb',  'name': 'Ytterbium',     'nvelec': 16, 'valence':  6, 'period':  9, 'group': 16, 'rcov': 1.87, 'rvdw': 2.28, 'mass': 173.054, 'color': [0.00, 0.75, 0.22]},
{'number':  71, 'symbol': 'Lu',  'name': 'Lutetium',      'nvelec':  3, 'valence':  6, 'period':  6, 'group':  3, 'rcov': 1.87, 'rvdw': 2.27, 'mass': 174.967, 'color': [0.00, 0.67, 0.14]},
{'number':  72, 'symbol': 'Hf',  'name': 'Hafnium',       'nvelec':  4, 'valence':  6, 'period':  6, 'group':  4, 'rcov': 1.75, 'rvdw': 2.25, 'mass': 178.490, 'color': [0.30, 0.76, 1.00]},
{'number':  73, 'symbol': 'Ta',  'name': 'Tantalum',      'nvelec':  5, 'valence':  6, 'period':  6, 'group':  5, 'rcov': 1.70, 'rvdw': 2.20, 'mass': 180.948, 'color': [0.30, 0.65, 1.00]},
{'number':  74, 'symbol': 'W',   'name': 'Tungsten',      'nvelec':  6, 'valence':  6, 'period':  6, 'group':  6, 'rcov': 1.62, 'rvdw': 2.10, 'mass': 183.840, 'color': [0.13, 0.58, 0.84]},
{'number':  75, 'symbol': 'Re',  'name': 'Rhenium',       'nvelec':  7, 'valence':  6, 'period':  6, 'group':  7, 'rcov': 1.51, 'rvdw': 2.05, 'mass': 186.207, 'color': [0.15, 0.49, 0.67]},
{'number':  76, 'symbol': 'Os',  'name': 'Osmium',        'nvelec':  8, 'valence':  6, 'period':  6, 'group':  8, 'rcov': 1.44, 'rvdw': 2.00, 'mass': 190.230, 'color': [0.15, 0.40, 0.59]},
{'number':  77, 'symbol': 'Ir',  'name': 'Iridium',       'nvelec':  9, 'valence':  6, 'period':  6, 'group':  9, 'rcov': 1.41, 'rvdw': 2.00, 'mass': 192.217, 'color': [0.09, 0.33, 0.53]},
{'number':  78, 'symbol': 'Pt',  'name': 'Platinum',      'nvelec': 10, 'valence':  6, 'period':  6, 'group': 10, 'rcov': 1.36, 'rvdw': 2.05, 'mass': 195.078, 'color': [0.90, 0.85, 0.68]},
{'number':  79, 'symbol': 'Au',  'name': 'Gold',          'nvelec': 11, 'valence':  6, 'period':  6, 'group': 11, 'rcov': 1.36, 'rvdw': 2.10, 'mass': 196.967, 'color': [0.80, 0.82, 0.12]},
{'number':  80, 'symbol': 'Hg',  'name': 'Mercury',       'nvelec': 12, 'valence':  6, 'period':  6, 'group': 12, 'rcov': 1.32, 'rvdw': 2.05, 'mass': 200.590, 'color': [0.71, 0.71, 0.76]},
{'number':  81, 'symbol': 'Tl',  'name': 'Thallium',      'nvelec':  3, 'valence':  3, 'period':  6, 'group': 13, 'rcov': 1.45, 'rvdw': 1.96, 'mass': 204.383, 'color': [0.65, 0.33, 0.30]},
{'number':  82, 'symbol': 'Pb',  'name': 'Lead',          'nvelec':  4, 'valence':  4, 'period':  6, 'group': 14, 'rcov': 1.46, 'rvdw': 2.02, 'mass': 207.200, 'color': [0.34, 0.35, 0.38]},
{'number':  83, 'symbol': 'Bi',  'name': 'Bismuth',       'nvelec':  5, 'valence':  3, 'period':  6, 'group': 15, 'rcov': 1.48, 'rvdw': 2.07, 'mass': 208.980, 'color': [0.62, 0.31, 0.71]},
{'number':  84, 'symbol': 'Po',  'name': 'Polonium',      'nvelec':  6, 'valence':  2, 'period':  6, 'group': 16, 'rcov': 1.40, 'rvdw': 1.97, 'mass': 209.000, 'color': [0.67, 0.36, 0.00]},
{'number':  85, 'symbol': 'At',  'name': 'Astatine',      'nvelec':  7, 'valence':  1, 'period':  6, 'group': 17, 'rcov': 1.50, 'rvdw': 2.02, 'mass': 210.000, 'color': [0.46, 0.31, 0.27]},
{'number':  86, 'symbol': 'Rn',  'name': 'Radon',         'nvelec':  8, 'valence':  0, 'period':  6, 'group': 18, 'rcov': 1.50, 'rvdw': 2.20, 'mass': 222.000, 'color': [0.26, 0.51, 0.59]},
{'number':  87, 'symbol': 'Fr',  'name': 'Francium',      'nvelec':  1, 'valence':  1, 'period':  7, 'group':  1, 'rcov': 2.60, 'rvdw': 3.48, 'mass': 223.000, 'color': [0.26, 0.00, 0.40]},
{'number':  88, 'symbol': 'Ra',  'name': 'Radium',        'nvelec':  2, 'valence':  2, 'period':  7, 'group':  2, 'rcov': 2.21, 'rvdw': 2.83, 'mass': 226.000, 'color': [0.00, 0.49, 0.00]},
{'number':  89, 'symbol': 'Ac',  'name': 'Actinium',      'nvelec':  3, 'valence':  6, 'period': 10, 'group':  3, 'rcov': 2.15, 'rvdw': 2.00, 'mass': 227.000, 'color': [0.44, 0.67, 0.98]},
{'number':  90, 'symbol': 'Th',  'name': 'Thorium',       'nvelec':  4, 'valence':  6, 'period': 10, 'group':  4, 'rcov': 2.06, 'rvdw': 2.40, 'mass': 232.038, 'color': [0.00, 0.73, 1.00]},
{'number':  91, 'symbol': 'Pa',  'name': 'Protactinium',  'nvelec':  5, 'valence':  6, 'period': 10, 'group':  5, 'rcov': 2.00, 'rvdw': 2.00, 'mass': 231.036, 'color': [0.00, 0.63, 1.00]},
{'number':  92, 'symbol': 'U',   'name': 'Uranium',       'nvelec':  6, 'valence':  6, 'period': 10, 'group':  6, 'rcov': 1.96, 'rvdw': 2.30, 'mass': 238.029, 'color': [0.00, 0.56, 1.00]},
{'number':  93, 'symbol': 'Np',  'name': 'Neptunium',     'nvelec':  7, 'valence':  6, 'period': 10, 'group':  7, 'rcov': 1.90, 'rvdw': 2.00, 'mass': 237.050, 'color': [0.00, 0.50, 1.00]},
{'number':  94, 'symbol': 'Pu',  'name': 'Plutonium',     'nvelec':  8, 'valence':  6, 'period': 10, 'group':  8, 'rcov': 1.87, 'rvdw': 2.00, 'mass': 244.060, 'color': [0.00, 0.42, 1.00]},
{'number':  95, 'symbol': 'Am',  'name': 'Americium',     'nvelec':  9, 'valence':  6, 'period': 10, 'group':  9, 'rcov': 1.80, 'rvdw': 2.00, 'mass': 243.060, 'color': [0.33, 0.36, 0.95]},
{'number':  96, 'symbol': 'Cm',  'name': 'Curium',        'nvelec': 10, 'valence':  6, 'period': 10, 'group': 10, 'rcov': 1.69, 'rvdw': 2.00, 'mass': 247.070, 'color': [0.47, 0.36, 0.89]},
{'number':  97, 'symbol': 'Bk',  'name': 'Berkelium',     'nvelec': 11, 'valence':  6, 'period': 10, 'group': 11, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 247.070, 'color': [0.54, 0.31, 0.89]},
{'number':  98, 'symbol': 'Cf',  'name': 'Californium',   'nvelec': 12, 'valence':  6, 'period': 10, 'group': 12, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 251.080, 'color': [0.63, 0.21, 0.83]},
{'number':  99, 'symbol': 'Es',  'name': 'Einsteinium',   'nvelec': 13, 'valence':  6, 'period': 10, 'group': 13, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 252.080, 'color': [0.70, 0.12, 0.83]},
{'number': 100, 'symbol': 'Fm',  'name': 'Fermium',       'nvelec': 14, 'valence':  6, 'period': 10, 'group': 14, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 257.100, 'color': [0.70, 0.12, 0.73]},
{'number': 101, 'symbol': 'Md',  'name': 'Mendelevium',   'nvelec': 15, 'valence':  6, 'period': 10, 'group': 15, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 258.100, 'color': [0.70, 0.05, 0.65]},
{'number': 102, 'symbol': 'No',  'name': 'Nobelium',      'nvelec': 16, 'valence':  6, 'period': 10, 'group': 16, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 259.100, 'color': [0.74, 0.05, 0.53]},
{'number': 103, 'symbol': 'Lr',  'name': 'Lawrencium',    'nvelec':  3, 'valence':  6, 'period':  7, 'group':  3, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 262.110, 'color': [0.78, 0.00, 0.40]},
{'number': 104, 'symbol': 'Rf',  'name': 'Rutherfordium', 'nvelec':  4, 'valence':  6, 'period':  7, 'group':  4, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 265.120, 'color': [0.80, 0.00, 0.35]},
{'number': 105, 'symbol': 'Db',  'name': 'Dubnium',       'nvelec':  5, 'valence':  6, 'period':  7, 'group':  5, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 268.130, 'color': [0.82, 0.00, 0.31]},
{'number': 106, 'symbol': 'Sg',  'name': 'Seaborgium',    'nvelec':  6, 'valence':  6, 'period':  7, 'group':  6, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 271.130, 'color': [0.85, 0.00, 0.27]},
{'number': 107, 'symbol': 'Bh',  'name': 'Bohrium',       'nvelec':  7, 'valence':  6, 'period':  7, 'group':  7, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 270.000, 'color': [0.88, 0.00, 0.22]},
{'number': 108, 'symbol': 'Hs',  'name': 'Hassium',       'nvelec':  8, 'valence':  6, 'period':  7, 'group':  8, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 277.150, 'color': [0.90, 0.00, 0.18]},
{'number': 109, 'symbol': 'Mt',  'name': 'Meitnerium',    'nvelec':  9, 'valence':  6, 'period':  7, 'group':  9, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 276.150, 'color': [0.92, 0.00, 0.15]},
{'number': 110, 'symbol': 'Ds',  'name': 'Darmstadtium',  'nvelec': 10, 'valence':  6, 'period':  7, 'group': 10, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 281.160, 'color': [0.93, 0.00, 0.14]},
{'number': 111, 'symbol': 'Rg',  'name': 'Roentgenium',   'nvelec': 11, 'valence':  6, 'period':  7, 'group': 11, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 280.160, 'color': [0.94, 0.00, 0.13]},
{'number': 112, 'symbol': 'Cn',  'name': 'Copernicium',   'nvelec': 12, 'valence':  6, 'period':  7, 'group': 12, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 285.170, 'color': [0.95, 0.00, 0.12]},
{'number': 113, 'symbol': 'Uut', 'name': 'Ununtrium',     'nvelec':  3, 'valence':  6, 'period':  7, 'group': 13, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 284.180, 'color': [0.96, 0.00, 0.11]},
{'number': 114, 'symbol': 'Fl',  'name': 'Flerovium',     'nvelec':  4, 'valence':  6, 'period':  7, 'group': 14, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 289.190, 'color': [0.97, 0.00, 0.10]},
{'number': 115, 'symbol': 'Uup', 'name': 'Ununpentium',   'nvelec':  5, 'valence':  6, 'period':  7, 'group': 15, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 288.190, 'color': [0.98, 0.00, 0.09]},
{'number': 116, 'symbol': 'Lv',  'name': 'Livermorium',   'nvelec':  6, 'valence':  6, 'period':  7, 'group': 16, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 293.000, 'color': [0.99, 0.00, 0.08]},
{'number': 117, 'symbol': 'Uus', 'name': 'Ununseptium',   'nvelec':  7, 'valence':  6, 'period':  7, 'group': 17, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 294.000, 'color': [0.99, 0.00, 0.07]},
{'number': 118, 'symbol': 'Uuo', 'name': 'Ununoctium',    'nvelec':  8, 'valence':  6, 'period':  7, 'group': 18, 'rcov': 1.60, 'rvdw': 2.00, 'mass': 294.000, 'color': [0.99, 0.00, 0.06]}]

def index_by_number(number):
   for element in periodic_table:
      if element['number'] == number:
         return periodic_table.index(element)

def index_by_symbol(symbol):
   for element in periodic_table:
      if element['symbol'].lower() == symbol.lower():
         return periodic_table.index(element)

def boolean_input(boolean):
   ierr = 0
   index = -1
   boolean = boolean.lower()
   if boolean == 'f' or boolean == 'false' or boolean == 'n' or boolean == 'no':
      index = 0
   elif boolean == 't' or boolean == 'true' or boolean == 'y' or boolean == 'yes':
      index = 1
   else:
      ierr = 1
   return ierr, index

if __name__ == "__main__":
   import sys
   sys.exit(main())

