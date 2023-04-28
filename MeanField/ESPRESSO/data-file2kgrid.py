#!/usr/bin/env python

# Script to create a kgrid.inp file automatically from a Quantum Espresso
# data-file.xml file.
#
# Felipe H. da Jornada (May, 2015)


import xml.etree.cElementTree as ET


def create_kgrid_inp(datafile, f_kgrid, kgrid, kshift, qshift, use_trs,
    out_cart, out_oct):

    f_kgrid.write((3*'{} '+'\n').format(*kgrid))
    f_kgrid.write((3*'{:.12f} '+'\n').format(*kshift))
    f_kgrid.write((3*'{:.12f} '+'\n').format(*qshift))

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
    f_kgrid.write(' {}\n'.format(a1))
    f_kgrid.write(' {}\n'.format(a2))
    f_kgrid.write(' {}\n'.format(a3))

    try:
        nat = int(tree.find('output/atomic_structure').get('nat'))
        alat = float(tree.find('output/atomic_structure').get('alat')) #in bohr
    except:
        nat = int(tree.find('IONS/NUMBER_OF_ATOMS').text.strip())
    f_kgrid.write('{}\n'.format(nat))

    root = tree.find('output/atomic_structure/atomic_positions')
    if root is not None:
        # New XML reports positions in bohrs
        symbols = []
        positions = []
        for atom in root.iter('atom'):
            positions.append(atom.text.strip())
            symbols.append(atom.get('name'))
        assert len(symbols)==nat
        unique_symbols = list(set(symbols))
        for iat in range(nat):
            symbol = symbols[iat]
            ityp = unique_symbols.index(symbol)
            tau = positions[iat]
            f_kgrid.write((' {} {}\n').format(ityp, tau))
    else:
        for iat in range(nat):
            node = tree.find('IONS/ATOM.{}'.format(iat+1))
            ityp = node.get('INDEX')
            tau = node.get('tau')
            f_kgrid.write((' {} {}\n').format(ityp, tau))

    for path in ('output/basis_set/fft_grid', 'PLANE_WAVES/FFT_GRID'):
        node = tree.find(path)
        if node is not None:
            fftgrid  = [node.get(nr) for nr in ('nr1', 'nr2', 'nr3')]
            break
    else:
        raise Exception('Could read FFT grid vector from XML file. Has the XML file changed?')
    f_kgrid.write('{} {} {}\n'.format(*fftgrid))


    def write_bool(b):
        if b:
            f_kgrid.write('.true.\n')
        else:
            f_kgrid.write('.false.\n')

    write_bool(use_trs)
    write_bool(out_cart)
    write_bool(out_oct)
    write_bool(False)


if __name__=="__main__":
    from argparse import ArgumentParser

    desc = ('Creates a kgrid.inp file given a data-file.xml from a Quantum'
    ' Espresso save directory.')
    parser = ArgumentParser(description=desc)

    group = parser.add_argument_group('required input and output files')
    group.add_argument('datafile', help='input data-file.xml from a QE save directory')
    group.add_argument('kgrid_inp', help='output kgrid.inp file to generate')

    group = parser.add_argument_group('kgrid specification')
    group.add_argument('--kgrid', type=int, default=[1,1,1], nargs=3,
    metavar=('nk1', 'nk2', 'nk3'),
        help='k-point density in each reciprocal-lattice direction. Defaults to 1.')
    group.add_argument('--kshift', type=float, default=[0.,0.,0.], nargs=3,
        metavar=('ks1', 'ks2', 'ks3'),
        help='k-shift in each reciprocal-lattice direction. Defaults to 0.0.')
    group.add_argument('--qshift', type=float, default=[0.,0.,0.], nargs=3,
        metavar=('qs1', 'qs2', 'qs3'),
        help='q0-shift in each reciprocal-lattice direction. Defaults to 0.0.')
    group.add_argument('--use_trs', default=False, action='store_true',
        help='Whether to use time-reversal symmetries. Defaults to false.')

    group = parser.add_argument_group('output control')
    group.add_argument('--output_cartesian', dest='out_cart', default=False,
        action='store_true', help=('Write output in Cartesian coordinates.'
        ' Default is false, which writes in crystal coordinates.'))
    group.add_argument('--output_octopus', dest='out_oct', default=False,
        action='store_true',  help=('Write output in Octopus format.'
        ' Default is false, which writes in Quantum Espresso format.'))

    args = parser.parse_args()

    with open(args.kgrid_inp, 'w') as f_kgrid:
        create_kgrid_inp(args.datafile, f_kgrid, args.kgrid, args.kshift,
            args.qshift, args.use_trs, args.out_cart, args.out_oct)
