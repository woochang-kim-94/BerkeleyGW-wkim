#!/usr/bin/env python

# This utility will upgrade a bsemat.h5 file generated with BerkeleyGW-1.1-beta
# (version 1) to the version used by BerkeleyGW-1.2.0 (version 2).
#
# Run this script with the --help option for information about input and 
# output arguments.
#
# Felipe H. da Jornada (2016)


from __future__ import print_function
import numpy as np
import h5py
import re


VER_BSE_HDF5 = 2
MATRICES = ('head', 'wing', 'body', 'exchange', 'fxc')


class KeywordsParser(object):


    def __init__(self, input):
        self.input = input

    def _find_keyword(self, kw):
        '''Returns True is the keyword kw is present in the input'''
        ret = re.search(' *{}'.format(kw), self.input)
        return ret is not None

    def _read_keyword_scalar(self, kw, fcast, default=None):
        '''Returns the scalar argument following the keyword kw, or None if
        it's not found.'''
        ret = re.search(' *{} *([^ #\n]*)'.format(kw), self.input)
        return fcast(ret.group(1)) if ret is not None else default

    def _read_keyword_vector(self, kw, fcast, default=None, vec_sz=3):
        '''Returns the vector argument following the keyword kw, or None if
        it's not found.'''
        ret = re.search(' *{}' + vec_sz*' *([^ #\n]*)'.format(kw), self.input)
        return np.array(map(fcast, ret.groups())) if ret is not None else default

    def map_keywords_toggle(self, map):
        for var, _dict in map.iteritems():
            # Get default value for the variable
            val = _dict.get(None, None)
            for kw, val_ in _dict.iteritems():
                if self._find_keyword(kw):
                    val = val_
            if val is None:
                raise ValueError('Variable {} was not set and has no default value.'.format(var))
            setattr(self, var, val)

    def map_keywords_scalar(self, map, fcast, default=None):
        for var, kw in map.iteritems():
            val = self._read_keyword_scalar(kw, fcast, default)
            if val is None:
                raise ValueError('Required input keyword {} not found.'.format(kw))
            setattr(self, var, val)

    def map_keywords_vector(self, map, fcast, default=None):
        for var, kw in map.iteritems():
            val = self._read_keyword_vector(kw, fcast, default)
            if val is None:
                raise ValueError('Required input keyword {} not found.'.format(kw))
            setattr(self, var, val)

    def write_datasets(self, g, kws):
        '''Write datasets specified b kws into group g, where kws is a list of str.'''

        for kw in kws:
            value = getattr(self, kw)
            shape = ()
            if isinstance(value, np.ndarray):
                shape = value.shape
                dtype = value.dtype
            else:
                if isinstance(value, float):
                    dtype = 'f8'
                elif isinstance(value, int):
                    dtype = 'i4'
                else:
                    raise ValueError('Unknown type of {}={} with type={}.'.
                                     format(kw, value, type(value)))

            g.create_dataset(kw, shape, dtype, value)


def file_needs_upgrade(f):
    print('Checking version of file {}.'.format(f.filename))
    if not 'versionnumber' in f:
        if 'bse_header' in f and 'versionnumber' in f['bse_header']:
            version = f['bse_header/versionnumber'][()]
            if version==VER_BSE_HDF5:
                print('File is already in version 3!')
                return False
        raise IOError('Unknown file version: {}'.format(version))
    version = f['versionnumber'][()]
    if version!=1:
        raise IOError('Unknown file version: {}'.format(version))
    return True


def upgrade(f, kernel_inp):
    print('Upgrading file {}.'.format(f.filename))
    print('Reading previous parameters from {}.'.format(f.filename))
    with open(kernel_inp) as f_inp:
        kp = KeywordsParser(f_inp.read())

    kp.flavor, kp.nk, kp.n2b, kp.n1b, kp.ns = map(int, f['info'])
    kp.flavor += 1
    assert kp.flavor==1 or kp.flavor==2, 'Wrong flavor in input file'
    kp.kpts = f['kpoints'][()]

    print('Parsing kernel input file {}.'.format(kernel_inp))
    # Parse toggle argument
    kp.map_keywords_toggle({
        'screening': {None: 0,
                      'screening_semiconductor': 0,
                      'screening_graphene': 1,
                      'screening_metal': 2},
        'icutv': {None: 0,
                  'spherical_truncation': 2,
                  'cell_wire_truncation': 4,
                  'cell_box_truncation': 5,
                  'cell_slab_truncation': 6},
        'nblocks': {None: 1,
                   'extended_kernel': 4},
        'theory': {None: 0,
                   'bse': 0,
                   'tddft': 1},
        'energy_loss': {None: 0,    
                        'energy_loss': 1},
        'patched_sampling': {None: 0,
                             'patched_sampling_co': 1}
    })

    # Parse required integer argument
    kp.map_keywords_scalar({
        'nvb': 'number_val_bands',
        'ncb': 'number_cond_bands',
    }, int)

    # Parse optional float argument
    kp.map_keywords_scalar({
        'short_range_frac_fock': 'short_range_frac_fock',
        'long_range_frac_fock': 'long_range_frac_fock',
        'screening_length': 'screening_length',
        'ecuts': 'screened_coulomb_cutoff',
        'ecutg': 'bare_coulomb_cutoff'
    }, float, default=0.)

    # Parse optional vector argument
    z3 = np.zeros((3,), dtype=np.float64)

    # Parse stuff related to center of mass q. This is slightly different..
    kp.qflag = 1
    kp.center_mass_q = z3.copy()
    tmp = kp._read_keyword_vector('qflag', float, vec_sz=4)
    if tmp is not None:
        kp.qflag = int(tmp[0])
        kp.center_mass_q = map(float, tmp[1:4])

    # Input fermi_level might be relative or absolute. This is too complicated to handle!
    kp.efermi = 0.
    kp.storage = 0
    kp.nspinor = 1 # Hard coded!
    kp.kgrid = np.array([0,0,0])
    kp.nmat = sum(map(lambda mat: 1 if mat in f else 0, MATRICES))

    # Note: we won`t create a mf_header nor an eps_header group.
    print('Writing new parameters to {}.'.format(f.filename))
    g = f.create_group('bse_header')
    g.create_dataset('flavor', (), 'i4', kp.flavor)
    g.create_dataset('versionnumber', (), 'i4', VER_BSE_HDF5)

    g = f.create_group('bse_header/params')
    fields = (
        'screening', 'icutv', 'ecuts', 'ecutg', 'efermi',
        'theory', 'nblocks', 'storage', 'nmat', 'energy_loss'
        )
    kp.write_datasets(g, fields)

    g = f.create_group('bse_header/bands')
    fields = (
        'nvb', 'ncb', 'n1b', 'n2b', 'ns', 'nspinor'
        )
    kp.write_datasets(g, fields)

    g = f.create_group('bse_header/kpoints')
    fields = (
        'nk', 'kpts', 'kgrid', 'qflag',
        'center_mass_q', 'patched_sampling'
        )
    kp.write_datasets(g, fields)

    print('Reorganizing file structure')
    f.create_group('mats')
    for matrix in MATRICES:
        if matrix in f:
            f['/mats/'+matrix] = f[matrix]
            del f[matrix]

    print('Cleaning up old datasets')
    del f['versionnumber']
    del f['info']
    del f['kpoints']

    print('All done!\n')


if __name__=="__main__":
    import argparse
    import shutil
    from os.path import exists

    desc = ('Given a bsemat.h5 version 1, generated with BerkeleyGW-1.1-beta '
            '(version 1) and the kernel.inp file used in the kernel '
            'calculation, this script will upgrade the file to the current '
            ' version (=2) used in BerkeleyGW 1.2.0.')

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('bsemat',
        help=('Input/output bsemat.h5 file. The script will overwrite the '
        'input file, but by default it will also create a backup of the '
        'bsemat.h5 file with a ".bck" extension.'))
    parser.add_argument('kernel_inp',
        help='Input kernel.inp file used in the kernel calculation.')
    parser.add_argument('-n', dest='backup', default=True, action='store_false',
        help='Don\'t create a backup of the input file. Default is to create.')
    parser.add_argument('-f', dest='force', default=False, action='store_true',
        help=('Force creation of backup file even if it already exists. '
        'Default is to halt the program if a backup exists and a new backup '
        'is requested.'))
    args = parser.parse_args()

    f = h5py.File(args.bsemat, 'r')
    try:
        do_upgrade = file_needs_upgrade(f)
    finally:
        f.close()

    if do_upgrade:
        if args.backup:
            fname_backup = args.bsemat + '.bck'
            if exists(fname_backup) and not args.force:
                raise OSError('Backup file {} already exists!'.format(fname_backup))
            shutil.copyfile(args.bsemat, fname_backup)

        f = h5py.File(args.bsemat, 'r+')
        upgrade(f, args.kernel_inp)
        f.close()
