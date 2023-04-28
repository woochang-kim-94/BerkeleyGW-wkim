#!/usr/bin/env python 
"""
Author: Gabriel Antonius
"""
from __future__ import print_function
from builtins import input
import numpy as np
from numpy import zeros
import netCDF4 as nc
import h5py

# Electron mass over atomical mass unit
me_amu = 5.4857990965007152E-4

def main():
    arguments = get_arguments_from_user()
    arguments['verbose'] = True  # Used for tests.
    abi2bgw_gkq(**arguments)

# =========================================================================== #


def abi2bgw_gkq(abi_gkk_fname, abi_ddb_fname, bgw_gkq_fname, verbose=False):
    """
    Read the Abinit GKK.nc file (in NetCDF format),
    convert it for BerkeleyGW in HDF5 format.

    Arguments
    ---------

    abi_gkk_fname: str
        Name of the NetCDF electron-phonon coupling matrix (GKK.nc) file
        produced by Abinit.

    abi_ddb_fname: str
        Name of the NetCDF dynamical matrix (DDB.nc) files
         produced by Abinit.

    bgw_gkq_fname: str
        Name of the output HDF5 gkq file for BerkeleyGW.

    Keyword arguments
    -----------------

    verbose:
        Print more information 

    """

    gkq = GkqFileH5()

    gkq.read_abinit_gkk(abi_gkk_fname, abi_ddb_fname)

    gkq.write_h5(bgw_gkq_fname)

    if verbose:
        print('abi2bgw_gkq: Success. 0')


# =========================================================================== #


class GkqFileH5(object):

    def __init__(self):

        self.ns = int()
        self.nk = int()
        self.nband = int()
        self.mband = int()
        self.nq = 1
        self.nat = int()
        self.ndir = 3
        self.cplx = 2

        # Arrays
        self.amu = None
        self.kpts = None
        self.qpts = None
        self.g_nu = None
        #self.g_at = None

    @property
    def nmode(self):
        return self.ndir * self.nat

    @property
    def natom(self):
        return self.nat

    @property
    def is_gamma(self):
        return np.allclose(self.qpts[0,:], 0.)

    def allocate(self):

        self.amu = np.zeros(self.nat, dtype=np.float)
        self.kpts = np.zeros((self.nk,3), dtype=np.float)
        self.qpts = np.zeros((self.nq,3), dtype=np.float)
        #self.g_at = np.zeros((self.nq, self.nmode, self.ns, self.nk,
        #                        self.nband, self.mband, self.cplx),
        #                        dtype=np.float)

        self.g_nu = np.zeros((self.nq, self.nmode, self.ns, self.nk,
                                self.nband, self.mband, self.cplx),
                                dtype=np.float)

    def write_h5(self, fname):

        f = h5py.File(name=fname, mode='w')
        
        f.create_group('gkq_header')
        f.create_group('gkq_data')

        sorting_dict = {
            '/gkq_header/ns' : self.ns,
            '/gkq_header/nk' : self.nk,
            '/gkq_header/nband' : self.nband,
            '/gkq_header/mband' : self.mband,
            '/gkq_header/nq' : self.nq,
            '/gkq_header/nat' : self.nat,
            '/gkq_header/ndir' : self.ndir,
            '/gkq_header/nmode' : self.nmode,
            '/gkq_data/amu' : self.amu,
            '/gkq_data/kpts' : self.kpts,
            '/gkq_data/qpts' : self.qpts,
            #'/gkq_data/g_at' : self.g_at,
            '/gkq_data/g_nu' : self.g_nu,
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


            f.create_dataset(name, shape=shape, data=data, dtype=dtype)

        f.close()

    def read_abinit_gkk(self, gkk_fname, ddb_fname):

        with nc.Dataset(gkk_fname, 'r') as ds:

            # Read the dimensions
            self.ns = len(ds.dimensions['number_of_spins'])
            self.nk = len(ds.dimensions['number_of_kpoints'])
            self.nband = len(ds.dimensions['product_mband_nsppol']) / self.ns
            self.mband = self.nband
            self.nat = len(ds.dimensions['number_of_atoms'])

            # Read the arrays
            amu = ds.variables['atomic_mass_units'][...]
            typat = ds.variables['atom_species'][...]

            kpts = ds.variables['reduced_coordinates_of_kpoints'][...]
            qpt = ds.variables['current_q_point'][...]

            # nband, nat, ndir, nk, nband * nspin * cplx,
            gkk = ds.variables['second_derivative_eigenenergies_actif'][...]

        # Reshape gkk in a more convenient form

        with nc.Dataset(ddb_fname, 'r') as ds:

            #self.nat = len(ds.dimensions['number_of_atoms'])
            #self.ndir = len(ds.dimensions['number_of_cartesian_directions'])  # 3
            self.ntypat = len(ds.dimensions['number_of_atom_species'])

            self.typat = ds.variables['atom_species'][:self.nat]
            self.amu = ds.variables['atomic_masses_amu'][:self.ntypat]
            self.rprim = ds.variables['primitive_vectors'][:self.ndir,:self.ndir]
            self.gprimd = np.linalg.inv(np.matrix(self.rprim))
            self.xred = ds.variables['reduced_atom_positions'][:self.nat,:self.ndir]
            self.qred = ds.variables['q_point_reduced_coord'][:]

            # The d2E/dRdR' matrix
            self.E2D = np.zeros((self.nat, self.ndir, self.nat, self.ndir), dtype=np.complex)
            self.E2D.real = ds.variables['second_derivative_of_energy'][:,:,:,:,0]
            self.E2D.imag = ds.variables['second_derivative_of_energy'][:,:,:,:,1]
            self.E2D = np.einsum('aibj->bjai', self.E2D)  # Indicies are reversed when writing them from Fortran.

        # Initialize the arrays
        self.allocate()

        self.kpts[...] = kpts
        self.qpts[0,:] = qpt

        for iat in range(self.nat):
            self.amu[iat] = amu[typat[iat]-1]

        # Reorder indices
        gkk = gkk.reshape(self.mband, self.nat, self.ndir, self.nk,
                          self.nband, self.ns, self.cplx)
        gkk = np.einsum('maiknsc->aisknmc', gkk)
        gkk_cplx = gkk[...,0] + 1j * gkk[...,1]
        #gkk = gkk.reshape(self.nmode, self.ns, self.nk, self.nband, self.mband, self.cplx)
        #self.g_at[0,...] = gkk[...]

        # Transform into mode basis
        polvec = self.get_reduced_displ()
        gkk_mode_cplx = np.einsum('aisknm,oia->osknm', gkk_cplx, polvec)
        self.g_nu[0,...,0] = np.real(gkk_mode_cplx)
        self.g_nu[0,...,1] = np.imag(gkk_mode_cplx)


    def get_mass_scaled_dynmat_cart(self):
        """
        Format the dynamical matrix in a 3Nx3N matrix,
        scale with masses, and transform into Cartesian coordinates.
        """
        # Retrive the amu for each atom
        amu = zeros(self.natom)
        for ii in np.arange(self.natom):
          jj = self.typat[ii]
          amu[ii] = self.amu[jj-1]
    
        # Transform from 2nd-order matrix (non-cartesian coordinates, 
        # masses not included, asr not included ) from self to
        # dynamical matrix, in cartesian coordinates, asr not imposed.
        E2D_cart = zeros((3,self.natom,3,self.natom),dtype=complex)
        for ii in np.arange(self.natom):
          for jj in np.arange(self.natom):
            for dir1 in np.arange(3):
              for dir2 in np.arange(3):
                for dir3 in np.arange(3):
                  for dir4 in np.arange(3):
                    E2D_cart[dir1,ii,dir2,jj] += (self.E2D[ii,dir3,jj,dir4] *
                            self.gprimd[dir1,dir3] * self.gprimd[dir2,dir4])

        # Reduce the 4 dimensional E2D_cart matrice to 2 dimensional Dynamical matrice
        # with scaled masses.
        Dyn_mat = zeros((3*self.natom,3*self.natom),dtype=complex)
        for ii in np.arange(self.natom):
          for dir1 in np.arange(3):
            ipert1 = ii * 3 + dir1
            for jj in np.arange(self.natom):
              for dir2 in np.arange(3):
                ipert2 = jj * 3 + dir2

                Dyn_mat[ipert1,ipert2] = (E2D_cart[dir1,ii,dir2,jj] *
                                          me_amu / np.sqrt(amu[ii]*amu[jj]))
    
        # Hermitianize the dynamical matrix
        #dynmat = np.matrix(Dyn_mat)
        dynmat = Dyn_mat
        dynmat = 0.5 * (dynmat + dynmat.transpose().conjugate())

        return dynmat


    def compute_dynmat(self, asr=None, zero_negative=True):
        """
        Diagonalize the dynamical matrix.
    
        Returns:
          omega: the frequencies, in Ha
          eigvect: the eigenvectors, in reduced coord
        """
    
        # Retrive the amu for each atom
        amu = zeros(self.natom)
        for ii in np.arange(self.natom):
          jj = self.typat[ii]
          amu[ii] = self.amu[jj-1]
    
        dynmat = self.get_mass_scaled_dynmat_cart()
    
        # Diagonalize the matrix
        eigval, eigvect = np.linalg.eigh(dynmat)
    
        # Scale the eigenvectors 
        for ii in np.arange(self.natom):
          for dir1 in np.arange(3):
            ipert = ii * 3 + dir1
            eigvect[ipert] = eigvect[ipert] * np.sqrt(me_amu / amu[ii])

        # Nullify imaginary frequencies
        if zero_negative:
            for i, eig in enumerate(eigval):
              if eig < 0.0:
                warnings.warn("An eigenvalue is negative with value: {} ... but proceed with value 0.0".format(eig))
                eigval[i] = 0.0

        # Impose the accoustic sum rule
        if asr and self.is_gamma:
          eigval[0] = 0.0
          eigval[1] = 0.0
          eigval[2] = 0.0

        # Frequencies
        self.omega = np.sqrt(np.abs(eigval)) * np.sign(eigval)
        self.eigvect = eigvect
    
        return self.omega, self.eigvect

    def get_reduced_displ(self):
        """
        Compute the mode eigenvectors, scaled by the mode displacements
        Also transform from cartesian to reduced coordinates.

        Returns: polvec[nmode,3,natom]
        """

        # Minimal value for omega (Ha)
        omega_tolerance = 1e-5

        self.polvec = zeros((self.nmode,3,self.natom), dtype=complex)
        xi_at = zeros(3, dtype=complex)

        omega, eigvect = self.compute_dynmat()

        for imode in range(self.nmode):
            # Skip mode with zero frequency (leave displacements null)
            if omega[imode].real < omega_tolerance:
              continue

            z0 = 1. / np.sqrt(2.0 * omega[imode].real)

            for iatom in np.arange(self.natom):
                for idir in range(3):
                    xi_at[idir] = eigvect[3*iatom+idir,imode] * z0

                for idir in range(3):
                    for jdir in range(3):
                        self.polvec[imode,idir,iatom] += xi_at[jdir] * self.gprimd[jdir,idir]

        return self.polvec
    

# =========================================================================== #


def get_arguments_from_user():
    """
    Prompt the user for the file names and the paramenters
    needed to connvert the wavefunction file.
    """

    arguments = dict()

    arguments['abi_gkk_fname'] = input(
        'Enter the name of the NetCDF gkk file produced by Abinit (GKK.nc):\n')

    arguments['abi_ddb_fname'] = input(
        'Enter the name of the NetCDF dynamical matrix file produced by Abinit (DDB.nc):\n')

    arguments['bgw_gkq_fname'] = input(
        'Enter the name of the output HDF5 gkq file for BerkeleyGW:\n')

    return arguments


# =========================================================================== #


if __name__ == '__main__':
    main()
