# define input parameters
# define constants
# always imported as a whole module

import numpy as np
from datetime import datetime


## input parameters (mandatory)
#convert = 'split-wannier-kpoints'   # preprocessing option
#convert = 'dft-unk'   # converting DFT-UNK wavefunctions
convert = 'wannier-interp-uwk'   # converting Wannier-interpolated UWK with rotmat

prefix = 'graphene'  # wannier90 prefix
wpath = 'unk/qe/'  # wannier run folders containing UNK or UWK files
nsplit = 0  # number of split of kpoints for wannier interpolation, for no-split, set to 0
wvfn_disk_io = False  # save intermediate wavefunctions on disk or in memory
restart = False  # restart from main object wavef saved in main folder

ecut = 60.0  # rydberg
efermi = 0.577598909  # eV, for insulators, set a value in the gap
latvec = np.array([[ 2.323902005,   4.025116356,   0.000000000], \
                   [-2.323902005,   4.025116356,   0.000000000], \
                   [ 0.000000000,   0.000000000,  15.000000000]])  # bohr

wfn_in = 'WFN_template.h5'  # input template wavefunction in hdf5
wfn_out = 'WFN_wannier.h5'  # output wannier-converted wavefunction in hdf5


## constants, start with '_', do not change
_tiny = 1.0e-8
_small = 1.0e-3
_rydberg2eV = 13.6056980659


## functions
def print_time(info_note):
  now = datetime.now()

  current_time = now.strftime("%H:%M:%S")
  print(info_note, current_time, '\n')

  return

