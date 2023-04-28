#!/usr/bin/env python

# This script merges an arbitrary number of epsmat.h5 files (in HDF5 format).
# The default output is epsmat_merged.h5, and it can be configured.
# The script requires h5py and numpy to work. On hopper, you should type:
# module load python h5py

# Felipe H. da Jornada <jornada@berkeley.edu> (2013)

import h5py
import sys
import os
from optparse import OptionParser
import numpy as np

usage = '%prog: [options] epsmat_1 [epsmat_2] [...]'
parser = OptionParser(usage=usage)
parser.add_option('--take-header-from', metavar='F_NUM', type='int',
                  help="don't check for consistency of the headers."
                       "Instead, use the header from the argument number F_NUM.")
parser.add_option('-o', '--output', default='epsmat_merged.h5', type='str',
                  help="output file. Defaults to 'epsmat_merged.h5'.")
parser.add_option('-f', '--force', default=False, action='store_true',
                  help="don't ask any questions, and overwrite output file if"
                       "it exists.")
(options, args) = parser.parse_args()
if len(args)<2:
    parser.error("wrong number of arguments")

print('Output file:\n> %s\n'%(options.output))
print('Input files:')
for fname in args:
    print('< %s'%(fname))
    if not os.path.exists(fname):
        raise IOError("File %s doesn't exist."%(fname))

fs_in = [h5py.File(fname, 'r') for fname in args]
if options.take_header_from is None:
    ind = 0
else:
    assert options.take_header_from >=1 and options.take_header_from <= len(args), \
        'invalid value (%d) for option --take-header-from.\n'%(options.take_header_from) + \
        'Expecting and integer between %d and %d'%(1,len(args))
    ind = options.take_header_from-1
f1 = fs_in[ind]

if (not options.force) and os.path.exists(options.output):
    # Fix Python 2.x.
    try: input = raw_input
    except NameError: pass
    ans = input('\nOutput file %s already exists. Would you '
                    'like to overwrite it? (y,N) '%(options.output))
    if len(ans)==0 or ans[0].lower()!='y':
        exit(0)
    print("Overwriting output file %s"%(options.output))

###############################################################################
print("\nAnalyzing files and checking consistency against %s\n"%(args[ind]))
###############################################################################

#FHJ: these dsets *must* match between different files
dsets = ['eps_header/versionnumber', 'eps_header/flavor', \
    'mf_header/gspace/ng', 'mf_header/gspace/components', \
    'eps_header/freqs/nfreq', 'eps_header/freqs/freqs']
dsets += [v.name for v in f1['eps_header/params'].values()
    if v.name.split('/')[-1] not in ('ecuts','nband','efermi',)]

#FHJ: match these dsets unless we use the --take-header-from option
if options.take_header_from is None:
    dsets += ['eps_header/qpoints/qgrid', 'mf_header/flavor']
    def parse_optional(name, obj):
        if isinstance(obj, h5py.Dataset):
            dsets.append(obj.name)
    f1['mf_header'].visititems(parse_optional)
else:
    print('  Note: skipping extra consistency checks because of --take-header-from')

#Checking consistency among files wrt f1
for f2 in fs_in:
    if f2==f1: continue
    print('  Checking %s'%(f2.filename))
    for dname in dsets:
        if not (f1[dname].shape==f2[dname].shape and
            np.allclose(f1[dname][()], f2[dname][()])):
            raise AttributeError('Incompatible values for dataset %s'%(dname))
print('  All files are consistent!')

flavor = f1['eps_header/flavor'][()]
ng = f1['mf_header/gspace/ng'][()]
nfreq = f1['eps_header/freqs/nfreq'][()]
freq_dep = f1['eps_header/freqs/freq_dep'][()]
ecuts = f1['eps_header/params/ecuts'][()]
print("\nSummary:\n")
print("  Flavor: %s"%((flavor==2 and 'Complex') or 'Real'))
print("  Frequency dependency: %d"%(freq_dep))
print("  Number of frequencies: %d"%(nfreq))
print("  Number of G-vectors: %d"%(ng))
print("  Epsilon cutoff: %f"%(ecuts))
nmtx_max = np.amax([f['eps_header/gspace/nmtx_max'][()] for f in fs_in])
print("  Maximum rank of the merged epsilon matrix (nmtx): %d"%(nmtx_max))
nq_tot = np.sum([f['eps_header/qpoints/nq'][()] for f in fs_in])
print("  Total number of q-points after merging: %d"%(nq_tot))

###############################################################################
print("\nWriting output parameters\n")
###############################################################################

f_out = h5py.File(options.output, 'w')
f_out.copy(f1['mf_header'], 'mf_header')
f_out.copy(f1['eps_header'], 'eps_header')
f_out['eps_header/qpoints/nq'][()] = nq_tot
f_out['eps_header/gspace/nmtx_max'][()] = nmtx_max

#Delete and merge all datasets related to nq
dsets_merge = {
    'eps_header/qpoints/qpts': (nq_tot,3),
    'eps_header/qpoints/qpt_done': (nq_tot,),
    'eps_header/gspace/nmtx': (nq_tot,),
    'eps_header/gspace/ekin': (nq_tot,ng),
    'eps_header/gspace/gind_eps2rho': (nq_tot,ng),
    'eps_header/gspace/gind_rho2eps': (nq_tot,ng)}
for name, shape in dsets_merge.items():
    del f_out[name]
    f_out.create_dataset(name, shape, dtype=f1[name].dtype)
    f_out[name][()] = np.concatenate([f[name][()] for f in fs_in], axis=0)

print("  Ok!")

###############################################################################
print("\nWriting output matrices\n")
###############################################################################

f_out.create_group('mats')
name = 'mats/matrix-diagonal'
shape = np.array(f1[name].shape)
shape[[0,1]] = [nq_tot, nmtx_max] #shape = (nq, nmtx_max, 1|2)
f_out.create_dataset(name, shape, f1[name].dtype)
name = 'mats/matrix'
shape = np.array(f1[name].shape)
shape[[0,3,4]] = [nq_tot, nmtx_max, nmtx_max] #shape = (nq, 1|2|3|4, nfreq, nmtx_max, nmtx_max, 1|2)
f_out.create_dataset(name, shape, f1[name].dtype)
name = 'eps_header/gspace/vcoul'
del f_out[name]
shape = [nq_tot, nmtx_max] #shape = (nq, nmtx_max)
f_out.create_dataset(name, shape, f1[name].dtype)

iq_glob = 0
for f in fs_in:
    nq = f['eps_header/qpoints/nq'][()]
    for iq in range(nq):
        nmtx_q = f['eps_header/gspace/nmtx'][iq]
        print("  Dealing with q=%d/%d (rank=%d)"%(iq_glob+1, nq_tot, nmtx_q))
        name = 'mats/matrix-diagonal'
        f_out[name][iq_glob,:nmtx_q,:] = (
            f[name][iq,:nmtx_q,:] )
        name = 'mats/matrix'
        f_out[name][iq_glob,:,:,:nmtx_q,:nmtx_q,:] = (
            f[name][iq,:,:,:nmtx_q,:nmtx_q,:] )
        name = 'eps_header/gspace/vcoul'
        f_out[name][iq_glob,:nmtx_q] = (
            f[name][iq,:nmtx_q] )
        iq_glob += 1
f_out.close()

print('\nAll done!\n')
