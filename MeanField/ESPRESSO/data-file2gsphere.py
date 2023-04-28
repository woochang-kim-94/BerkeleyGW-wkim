#!/usr/bin/env python

# Script to create a gsphere.inp file automatically from a Quantum Espresso
# data-file.xml file.
#
# For use with Visual/gsphere.py.
#
# On command line, specify:
# (1) path to data-file.xml from espresso calculation
# (2) desired name for gsphere.py input file
# (3) g-sphere cutoff (in Ry)
# (4, optional) logical for g-vector ordering (default is false)
# (5, optional) k-point for use with g-sphere (default is origin)
# 
# e.g., data-file2gsphere.py gsphere.inp --gcut 25
#
# Confirmed to work with Python version 2.7.9.
# Incompatible with Python version 2.6.6.
#
# Bradford A. Barker (Sept, 2015), after FHJ data-file2kgrid.py

import xml.etree.cElementTree as ET

def create_gsphere_inp(datafile, f_gsphere, gcut, order_gvecs, kpoint):

    tree = ET.parse(datafile)

    for paths in ('output/atomic_structure/cell', 'CELL/DIRECT_LATTICE_VECTORS/a1'):
        node = tree.find('output/atomic_structure/cell')
        if node is not None:
            a1 = node.find('a1').text.strip()
            a2 = node.find('a2').text.strip()
            a3 = node.find('a3').text.strip()
            break
    else:
        raise Exception('Could read lattice vector from XML file. Has the XML file changed?')
    f_gsphere.write(' {}\n'.format(a1))
    f_gsphere.write(' {}\n'.format(a2))
    f_gsphere.write(' {}\n'.format(a3))

    f_gsphere.write(' {}\n'.format(gcut))

    for path in ('output/basis_set/fft_grid', 'PLANE_WAVES/FFT_GRID'):
        node = tree.find(path)
        if node is not None:
            fftgrid  = [node.get(nr) for nr in ('nr1', 'nr2', 'nr3')]
            break
    else:
        raise Exception('Could read FFT grid vector from XML file. Has the XML file changed?')
    f_gsphere.write('{} {} {}\n'.format(*fftgrid))

    def write_bool(b):
        if b:
            f_gsphere.write('true\n')
        else:
            f_gsphere.write('false\n')

    write_bool(order_gvecs)

    f_gsphere.write((3*'{:.12f} '+'\n').format(*kpoint))


if __name__=="__main__":
    from argparse import ArgumentParser

    desc = ('Creates a gsphere.inp file given a data-file.xml from a Quantum'
    ' Espresso save directory.')
    parser = ArgumentParser(description=desc)

    group = parser.add_argument_group('required input and output files')
    group.add_argument('datafile', help='input data-file.xml from a QE save directory')
    group.add_argument('gsphere_inp', help='output kgrid.inp file to generate')

    group = parser.add_argument_group('required value of gvector cutoff')
    group.add_argument('--gcut', help='value of gvector cutoff (units of Ry)')

    group = parser.add_argument_group('order gvectors (optional)')
    group.add_argument('--order_gvecs', default=False, action='store_true',
        help='Whether to order gvectors. Defaults to false.')

    group = parser.add_argument_group('kpoint specification (optional)')
    group.add_argument('--kpoint', type=float, default=[0.0,0.0,0.0], nargs=3,
    metavar=('k1', 'k2', 'k3'),
        help='k-point about which to calculate the G-sphere. Defaults to origin.')

    args = parser.parse_args()

    with open(args.gsphere_inp, 'w') as f_gsphere:
        create_gsphere_inp(args.datafile, f_gsphere, args.gcut, args.order_gvecs, args.kpoint)
