#!/usr/bin/env python 
"""
Script to convert an Abinit wavefunction in NetCDF format into HDF5 format
for BerkeleyGW.

The best way to use this script is to import the function 'abi2bgw_nc_h5'
and use it in a python script.

Alternatively, the script can be used interactively by feeding the arguments
through the standard input. e.g.:

    $ abi2bgw_nc_h5 < stdin

with stdin containing the lines

    abi_wfn_fname
    bgw_wfn_fname
    k1 k2 k3
    sk1 sk2 sk3
    cell_symmetry

where
    
    abi_wfn_fname =
        Name of the NetCDF wavefunction file produced by Abinit

    bgw_wfn_fname =
        Name of the output HDF5 wavefunction for BerkeleyGW

    k1 k2 k3 =
        Number of divisions in the k-point grid

    sk1 sk2 sk3 = 
        Shift for the k-point grid sk1 sk2 sk3

    cell_symmetry =
        Symmetry of the unit cell (0=cubic, 1=hexagonal)


Author: Gabriel Antonius
"""

import numpy as np
import netCDF4 as nc
import h5py


def main():
    arguments = get_arguments_from_user()
    abi2bgw_nc_h5(**arguments)


# =========================================================================== #


def abi2bgw_nc_h5(abi_wfn_fname, bgw_wfn_fname,
                  kgrid,
                  shiftk=[0.,0.,0.],
                  cell_symmetry=0,
                  ):
    """
    Read the Abinit wavefunction in NetCDF format,
    convert it for BerkeleyGW in HDF5 format.

    Arguments
    ---------

    abi_wfn_fname: str
        Name of the NetCDF wavefunction file produced by Abinit

    bgw_wfn_fname: str
        Name of the output HDF5 wavefunction for BerkeleyGW

    kgrid: int(3)
        Number of divisions in the k-point grid k1 k2 k3

    shiftk: float(3)
        Shift for the k-point grid sk1 sk2 sk3

    cell_symmetry: int
        Symmetry of the unit cell (0=cubic, 1=hexagonal)
    """

    Ha_to_Ry = 2

    wfn = WavefunctionFileH5()

    with nc.Dataset(abi_wfn_fname, 'r') as ds:

        # Read the dimensions
        wfn.flavor = len(ds.dimensions['real_or_complex_coefficients'])
        wfn.nspin = len(ds.dimensions['number_of_spins'])
        wfn.nspinor = len(ds.dimensions['number_of_spinor_components'])
        wfn.mnband = len(ds.dimensions['max_number_of_states'])
        wfn.nrk = len(ds.dimensions['number_of_kpoints'])
        wfn.ntran = len(ds.dimensions['number_of_symmetry_operations'])
        wfn.nat = len(ds.dimensions['number_of_atoms'])
        wfn.ngkmax = len(ds.dimensions['max_number_of_coefficients'])

        n1 = len(ds.dimensions['number_of_grid_points_vector1'])
        n2 = len(ds.dimensions['number_of_grid_points_vector2'])
        n3 = len(ds.dimensions['number_of_grid_points_vector3'])

        # ==== crystal ==== #
        wfn.alat = 1.
        wfn.avec = ds.variables['primitive_vectors'][...]
        wfn.adot = np.matmul(wfn.avec, wfn.avec.transpose()) * wfn.alat ** 2

        wfn.blat = 2 * np.pi
        wfn.bvec = np.linalg.inv(wfn.avec * wfn.alat).transpose()
        wfn.bdot = np.matmul(wfn.bvec, wfn.bvec.transpose()) * wfn.blat ** 2

        wfn.celvol = np.linalg.det(wfn.avec) * wfn.alat ** 3
        wfn.recvol = np.linalg.det(wfn.bvec) * wfn.blat ** 3

        wfn.atyp = np.zeros((wfn.nat), dtype=np.int)
        znucl = ds.variables['atomic_numbers'][...]
        typat =  ds.variables['atom_species'][...]
        for i, j in enumerate(typat):
            wfn.atyp[i] = znucl[j-1]

        wfn.apos = np.zeros((wfn.nat,3), dtype=np.float)
        xred = ds.variables['reduced_atom_positions'][...]
        for iat in range(wfn.nat):
            x = xred[iat,:]
            wfn.apos[iat] = np.dot(wfn.avec.transpose(), x) * wfn.alat


        # ==== kpoints ==== #
        wfn.rk = ds.variables['reduced_coordinates_of_kpoints'][...]
        wfn.w = ds.variables['kpoint_weights'][...]

        assert len(ds.dimensions['number_of_spins']) == 1, 'nshiftk must be 1'
        #wfn.shiftk = ds.variables['shiftk'][0,:]
        wfn.kgrid = np.array(kgrid)
        wfn.shift = np.array(shiftk)
        # This is not always present. It will have to be provided on input.
        #kptrlatt = ds.variables['kptrlatt'][...]

        wfn.ecutwfc = float(ds.variables['ecut_eff'][...]) * Ha_to_Ry

        wfn.el = ds.variables['eigenvalues'][...] * Ha_to_Ry
        wfn.occ = ds.variables['occupations'][...] * wfn.nspin / 2.

        wfn.ifmin = np.ones((wfn.nspin,wfn.nrk), dtype=np.int)
        wfn.ifmax = np.ones((wfn.nspin,wfn.nrk), dtype=np.int)

        for ispin in range(wfn.nspin):
            for ikpt in range(wfn.nrk):
                for iband in range(wfn.mnband):
                    occ = wfn.occ[ispin,ikpt,iband]
                    if occ <= wfn.nspin / 2.:
                        break
                wfn.ifmax[ispin,ikpt] = iband  # Fortran convention

        wfn.ngk = ds.variables['number_of_coefficients'][...]
        wfn.ngksum = np.sum(wfn.ngk)

        # ==== symmetry ==== #

        wfn.cell_symmetry = cell_symmetry  # 0 = cubic, 1 = hexagonal.

        wfn.mtrx = np.zeros((wfn.ntran,3,3), dtype=np.int)
        wfn.tnp = ds.variables['reduced_symmetry_translations'][...] * wfn.blat
        for i in range(wfn.ntran):
            M =  ds.variables['reduced_symmetry_matrices'][i,:,:]
            wfn.mtrx[i,:,:] = np.linalg.inv(np.transpose(M))


        # ==== gspace ==== #
        wfn.ecutrho = 4 * wfn.ecutwfc

        gmax = np.sqrt(2 * wfn.ecutrho / Ha_to_Ry)
        gnorms = [np.linalg.norm(wfn.bvec[i]) * wfn.blat for i in range(3)]
        igmax = [ int(np.ceil(gmax / n)) for n in gnorms ]

        #wfn.FFTgrid = np.array([2 * i + 1 for i in igmax], dtype=np.int)
        wfn.FFTgrid = np.array([n1,n2,n3], dtype=np.int)

        components = list()
        for ig1 in range(-igmax[0], igmax[0]+1):
            for ig2 in range(-igmax[1], igmax[1]+1):
                for ig3 in range(-igmax[2], igmax[2]+1):

                    Gred = np.array([ig1,ig2,ig3])
                    Gcart = np.dot(wfn.bvec.transpose(), Gred) * wfn.blat

                    ecut = (sum(Gcart ** 2) / 2.) * Ha_to_Ry
                    if ecut <= wfn.ecutrho:
                        components.append(Gred)

        wfn.ng = len(components)
        wfn.components = np.array(components, dtype=np.int)
        

        # ==== wfns ==== #

        wfn.gvecs = np.zeros((wfn.ngksum,3), dtype=np.int)

        wfn.coeffs = np.zeros((wfn.mnband,wfn.nspin,wfn.ngksum,2),
                               dtype=np.float)

        igk_start = 0
        for ik, ng in enumerate(wfn.ngk):
            igk_end = igk_start + ng

            wfn.gvecs[igk_start:igk_end,:] = (
                ds.variables['reduced_coordinates_of_plane_waves'][ik,:ng,:])

            for ispin in range(wfn.nspin):
              for ispinor in range(wfn.nspinor):
                istot = ispin * wfn.nspinor + ispinor
                for iband in range(wfn.mnband):
                  for icplx in range(wfn.flavor):
                    wfn.coeffs[iband,istot,igk_start:igk_end,icplx] = (
                    ds.variables['coefficients_of_wavefunctions']
                                [ispin,ik,iband,ispinor,:ng,icplx])

            igk_start = igk_end

    # Write the file
    wfn.write_h5(bgw_wfn_fname)


# =========================================================================== #


class WavefunctionFileH5(object):

    def __init__(self):

        # flavor
        self.flavor = np.int(2)

        # ==== crystal ==== #
        self.nat = np.int()  # Number of atoms

        self.adot = np.zeros((3,3), dtype=np.float)
        self.alat = np.float()
        self.avec = np.zeros((3,3), dtype=np.float)

        self.apos = np.zeros((self.nat,3), dtype=np.float)
        self.atyp = np.zeros((self.nat), dtype=np.int)

        self.bdot = np.zeros((3,3), dtype=np.float)
        self.blat = np.float()
        self.bvec = np.zeros((3,3), dtype=np.float)

        self.celvol = np.float()
        self.recvol = np.float()

        # ==== gspace ==== #
        self.ng = np.int()  # Number of G vectors for density
        self.FFTgrid = np.zeros((3), dtype=np.int)
        self.components = np.zeros((self.ng,3), dtype=np.int)  # G vectors for density
        self.ecutrho = np.float()

        # ==== kpoints ===== #
        self.nspin = np.int()    # Number of spin polarizations
        self.nspinor = np.int()  # Number of spinor components
        self.nrk = np.int()      # Number of kpoints
        self.mnband = np.int()   # Number of bands

        self.ecutwfc = np.float()
        self.ifmin = np.ones((self.nspin,self.nrk), dtype=np.int)
        self.ifmax = np.ones((self.nspin,self.nrk), dtype=np.int)
        self.el = np.zeros((self.nspin,self.nrk,self.mnband), dtype=np.float)
        self.occ = np.zeros((self.nspin,self.nrk,self.mnband), dtype=np.float)

        # Number of G vectors at each k-point
        self.ngk = np.zeros((self.nrk), dtype=np.int)
        self.ngkmax = np.int()
        self.ngksum = np.int()

        self.rk = np.zeros((self.nrk,3), dtype=np.float)
        self.shift = np.zeros((3), dtype=np.float)
        self.w = np.zeros((self.nrk), dtype=np.float)

        self.kgrid = np.zeros((3), dtype=np.int)

        # ==== symmetry ==== #
        self.cell_symmetry = np.int()
        self.ntran = np.int()     # Number of symmetry operations
        self.versionnumber = np.int(1)

        self.mtrx = np.zeros((self.ntran,3,3), dtype=np.int)
        self.tnp = np.zeros((self.ntran,3), dtype=np.float)

        # ==== wfns ==== #
        self.coeffs = np.zeros((self.mnband,self.nspin,self.ngksum,2),
                               dtype=np.float)

        # G vectors for wfn
        self.gvecs = np.zeros((self.ngksum,3), dtype=np.int)

    def write_h5(self, fname):

        f = h5py.File(name=fname, mode='w')
        
        f.create_group('mf_header')
        f.create_group('mf_header/crystal')
        f.create_group('mf_header/gspace')
        f.create_group('mf_header/kpoints')
        f.create_group('mf_header/symmetry')
        f.create_group('wfns')

        sorting_dict = {
            '/mf_header/crystal/adot' : self.adot,
            '/mf_header/crystal/alat' : self.alat,
            '/mf_header/crystal/apos' : self.apos,
            '/mf_header/crystal/atyp' : self.atyp,
            '/mf_header/crystal/avec' : self.avec,
            '/mf_header/crystal/bdot' : self.bdot,
            '/mf_header/crystal/blat' : self.blat,
            '/mf_header/crystal/bvec' : self.bvec,
            '/mf_header/crystal/celvol' : self.celvol,
            '/mf_header/crystal/nat' : self.nat,
            '/mf_header/crystal/recvol' : self.recvol,
            '/mf_header/flavor' : self.flavor,
            '/mf_header/gspace/FFTgrid' : self.FFTgrid,
            '/mf_header/gspace/components' : self.components,
            '/mf_header/gspace/ecutrho' : self.ecutrho,
            '/mf_header/gspace/ng' : self.ng,
            '/mf_header/kpoints/ecutwfc' : self.ecutwfc,
            '/mf_header/kpoints/el' : self.el,
            '/mf_header/kpoints/ifmax' : self.ifmax,
            '/mf_header/kpoints/ifmin' : self.ifmin,
            '/mf_header/kpoints/kgrid' : self.kgrid,
            '/mf_header/kpoints/mnband' : self.mnband,
            '/mf_header/kpoints/ngk' : self.ngk,
            '/mf_header/kpoints/ngkmax' : self.ngkmax,
            '/mf_header/kpoints/nrk' : self.nrk,
            '/mf_header/kpoints/nspin' : self.nspin,
            '/mf_header/kpoints/nspinor' : self.nspinor,
            '/mf_header/kpoints/occ' : self.occ,
            '/mf_header/kpoints/rk' : self.rk,
            '/mf_header/kpoints/shift' : self.shift,
            '/mf_header/kpoints/w' : self.w,
            '/mf_header/symmetry/cell_symmetry' : self.cell_symmetry,
            '/mf_header/symmetry/mtrx' : self.mtrx,
            '/mf_header/symmetry/ntran' : self.ntran,
            '/mf_header/symmetry/tnp' : self.tnp,
            '/mf_header/versionnumber' : self.versionnumber,
            '/wfns/coeffs' : self.coeffs,
            '/wfns/gvecs' : self.gvecs,
            }

        for name, data in sorting_dict.items():

            if 'shape' in dir(data):
                shape = data.shape
            else:
                shape = ()

            if 'dtype' in dir(data):
                typ = data.dtype
            else:
                typ = type(data)

            if typ in (float, np.float, np.float32, np.float64):
                dtype = h5py.h5t.IEEE_F64LE
            elif typ in (int, np.int, np.int32, np.int64):
                dtype = h5py.h5t.STD_I32LE
            else:
                raise Exception('Data type not recognized for variable '
                                '{} : {}'.format(name, typ))


            f.create_dataset(name,shape=shape,data=data,dtype=dtype)

        f.close()


# =========================================================================== #


def get_arguments_from_user():
    """
    Prompt the user for the file names and the paramenters
    needed to connvert the wavefunction file.
    """

    arguments = dict()

    arguments['abi_wfn_fname'] = raw_input(
        'Enter the name of the NetCDF wavefunction file produced by Abinit: ')

    arguments['bgw_wfn_fname'] = raw_input(
        'Enter the name of the output HDF5 wavefunction for BerkeleyGW: ')

    kgrid = raw_input(
        'Enter the number of divisions in the k-point grid k1 k2 k3 : ')
    kgrid = [int(ki) for ki in kgrid.split()]
    assert len(kgrid) == 3, 'You must enter 3 integers for the kpoint grid: k1 k2 k3'
    arguments['kgrid'] = kgrid

    shiftk = raw_input(
        'Enter the shift for the k-point grid sk1 sk2 sk3 : ')
    shiftk = [float(ki) for ki in shiftk.split()]
    assert len(shiftk) == 3, 'You must enter 3 floats for the shift of the kpoint grid: sk1 sk2 sk3'
    arguments['shiftk'] = shiftk

    cell_symmetry = raw_input(
        'Enter the symmetry of the unit cell (0=cubic, 1=hexagonal) : ')
    arguments['cell_symmetry'] = int(cell_symmetry)

    return arguments


# =========================================================================== #


if __name__ == '__main__':
    main()
