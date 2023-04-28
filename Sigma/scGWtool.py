#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# This python script is used to do subspace self-consistent GW with BerkeleyGW package.   #
# This is used between GW steps: reading off-diagonal matrix elements, diagonalize H^{GW} #
# matrix (subspace), update wavefunctions, rotate energy and self-energy matrices from    #
# last step, and generate related files to be used for next iteration.                    #
#                                                                                         #
# File needed: sigma_hp_col.log sigma_hp_row.log                                          #
#              (rename sigma_hp.log from two separate calculations:                       #
#               one evaluate with column states, the other row states.                    #
#               see manual for off-diagolan sigma calculation for more clarifications.)   #
#              wfn_old.h5 (last step WFN, binary can be converted to h5, or vice versa    #
# File output: wfn_new.h5 vxc_new.dat eqp_rotate.dat                                      #
#              (use wfn_new for next step Sigma, use vxc_new.dat as vxc.dat.              #
#               eqp_rotate.dat is the rotation matrix kept for analysis.)                 #
#                                                                                         #
#                                                                                         #
# Note: scGW using this script updates both eigvalues and wavefunctions. However,         #
#       W remains the same by default (still GW_0 level), unless one recomputes W (but    #
#       may have vertex correction self-consistency issue.)                               #
#       Certain approximations are still involved since we did not explicitly construct   #
#       other mean-field terms with new wavefunctions (but just using basis rotation).    #
#                                                                                         #       
#       Spin unpolarized case and full spinor case are supproted for now.                 #
#                                                                                         #
# Author: Zhenglu Li                                                                      #
# Date: 08/10/2015                                                                        #
#                                                                                         #
# Modified from Zhenglu Li original script by Bradford A. Barker                          #
# Modified to allow spinor wavefunctions,                                                 #
# and Upper/Lower triangular matrices...                                                  #
# Date Modified: April 2018.                                                              #
#                                                                                         #
# Last modified by Jiawei Ruan, 2020.                                                     #
#                                                                                         #
# Thanks to Diana Qiu, Ting Cao, Felipe da Jornada, Meng Wu                               #
#                                                                                         #
###########################################################################################

#
# Set the following line in sigma.inp:
#      sigma_matrix   -1/-2   0
# we are averaging Sigma evaluated at two energies, from colomn and row
# However they are usually very close to each other, so using one might be enough
#
#
# To keep consistency with BerkeleyGW, we now force WFN_outer and WFN_iner use same number
# of kpoints in the same order, although not necessarily.
# kmap is defined under this assumption, or errors will occur.
#


import sys
import h5py
import numpy as np


def read_sigma_hp():
    r"""
    This function reads sigma_hp_col.log sigma_hp_row file, construct three matrices:
        hammf[k,i,j], hamqp[k,i,j], sigma[k,i,j]
        kmap, bmap, kpts

    NOTE: kmap and bmap both map to the real indices starting from 1, when using as
          python indices, ALWAYS -1

    return [hammf, hamqp, sigma, kmap, bmap]
    """

    fsig = open("sigma_hp_col.log", "r")
    print()
    print("Reading sigma_hp_col.log file ...")

    ln = fsig.readlines()
    fsig.close()

    # First time reading, to get general information

    nk = 0
    nb = 0
    # 0 for full matrix (default); 1 for lower triangular; 2 for upper
    mattype = 0

    for il in range(len(ln)):
        sp = ln[il].split()

        if len(sp) >= 3:
            if sp[0] == "band_index":
                bandmin = int(sp[1])
                bandmax = int(sp[2])

            if sp[0] == "sigma_matrix":
                if sp[1] == "-1" and sp[2] == "0":
                    print("Using full matrix")
                elif sp[1] == "-1" and sp[2] == "-1":
                    print("Using lower triangular matrix")
                    mattype = 1
                elif sp[1] == "-1" and sp[2] == "1":
                    print("Using upper triangular matrix")
                    mattype = 2
                else:
                    print("Now only supports Sigma with")
                    print("sigma_matrix   -1   -1/0/+1")
                    sys.exit(0)

            if sp[0] == "k" and sp[1] == "=":
                nk = nk + 1

    nb = bandmax - bandmin + 1

    print()
    print("Info from sigma_hp_col.log:")
    print("number of kpoints calculated:", nk)
    print("subspace size is", nb, "x", nb,
          "within bandmin, bandmax:", bandmin, bandmax)
    print()

    # Now define the matrices

    hamqp = np.zeros((nk, nb, nb), dtype=np.complex)
    hammf = np.zeros((nk, nb, nb), dtype=np.complex)
    sigma = np.zeros((nk, nb, nb), dtype=np.complex)
    vxc = np.zeros((nk, nb, nb), dtype=np.complex)

    kmap = np.zeros(nk, dtype=np.int)
    bmap = np.zeros(nb, dtype=np.int)
    kpts = np.zeros((nk, 3), dtype=np.float)

    # generate bmap

    for ib in range(nb):
        bmap[ib] = bandmin + ib

    # Second time reading to get hammf, kmap

    ik = -1

    for il in range(len(ln)):
        sp = ln[il].split()

        if len(sp) >= 2:
            if sp[0] == "k" and sp[1] == "=":
                ik = ik + 1
                kmap[ik] = int(sp[7])
                kpts[ik] = [float(sp[2]), float(sp[3]), float(sp[4])]

                for ib in range(nb):
                    spb = ln[il+3+ib].split()

                    if int(spb[0]) != bmap[ib]:
                        print("Band index not matching in reading sigma_hp_col.log")
                        sys.exit(0)

                    hammf[ik][ib][ib] = float(spb[1])

    # Third time reading to get sigma, vxc
    #read_sigma_vxc(ln, nb, sigma, vxc, bmap)
    read_sigma_vxc(mattype, ln, nb, sigma, vxc, bmap)

    # Now we open sigma_hp_row, and do the same step as the third step above, add them up
    # and then divided by 2, we then get the correct sigma and vxc matrices

    # BAB note: Can insert a consistency check for sigma_hp_row, but we are assuming same parameters as _col.log

    fsig = open("sigma_hp_row.log", "r")
    print()
    print("Reading sigma_hp_row.log file ...")

    ln = fsig.readlines()
    fsig.close()
    #read_sigma_vxc(ln, nb, sigma, vxc, bmap)
    read_sigma_vxc(mattype, ln, nb, sigma, vxc, bmap)

    sigma = sigma / 2.0
    vxc = vxc / 2.0
    # Now we have sigma and vxc matrices

    # H^{GW} = H^{DFT} - V_{xc} + \Sigma
    # under the same basis

    hamqp = hammf - vxc + sigma

    print()
    print("Quasi-particle energies from sigma_hp.log")
    for ik in range(nk):
        print("kpoint", ik, ", kmap", kmap[ik])
        for ib in range(nb):
            print("   Eqp0 (", bmap[ib], ") = ", hamqp[ik][ib][ib].real)

    return (hammf, hamqp, sigma, bmap, kmap, kpts)


# def read_sigma_vxc(ln, nb, sigma, vxc, bmap):
def read_sigma_vxc(mattype, ln, nb, sigma, vxc, bmap):
    r"""
    Called by read_sigma_hp, read vxc and sigma matrices
    """

    ik = -1

    for il in range(len(ln)):
        sp = ln[il].split()

        if len(sp) >= 3:
            if sp[0] == "n" and sp[1] == "m" and sp[2] == "l":
                ik = ik + 1

                ilb = 0

                # Create arrays for ib1 and ib2 based on full, lower tri, upper tri
                # instead of putting if statements over the choices of for-loops
                if mattype == 0:
                    ib1array = list(range(nb))
                    ib2array = list(range(nb))
                elif mattype == 1:
                    ib1array = []
                    ib2array = list(range(nb))
                    for ib2 in ib2array:
                        ib1array.append(list(range(ib2, nb)))
                elif mattype == 2:
                    ib1array = list(range(nb))
                    ib2array = []
                    for ib1 in ib1array:
                        ib2array.append(list(range(ib1, nb)))

                # Cycle through calculated matrix elements of Sig and VXC matrices
                for ib1 in ib1array:
                    # for ib2 in ib2array[ib1]:
                    for ib2 in ib2array:
                        spb = ln[il + ilb + 1].split()
                        if int(spb[0]) != bmap[ib1] or int(spb[1]) != bmap[ib2] or spb[3] != "real":
                            print(
                                "Info not matching in off-diag reading in sigmg_hp_col.log")
                            sys.exit(0)
                        # Initialization MUST be ZEROS
                        sigma[ik][ib1][ib2] += float(spb[7])
                        # Initialization MUST be ZEROS
                        vxc[ik][ib1][ib2] += float(spb[8])

                        spb = ln[il + ilb + 2].split()
                        if int(spb[0]) != bmap[ib1] or int(spb[1]) != bmap[ib2] or spb[3] != "imag":
                            print(
                                "Info not matching in off-diag reading in sigmg_hp_col.log")
                            sys.exit(0)
                        sigma[ik][ib1][ib2] += 1.0j * float(spb[7])
                        vxc[ik][ib1][ib2] += 1.0j * float(spb[8])

                        ilb = ilb + 2

                # Reconstruct full matrix if using upper or lower triangular matrices...
                if mattype == 1:
                    for ib1 in range(nb):
                        for ib2 in range(ib1+1, nb):
                            sigma[ik][ib1][ib2] = sigma[ik][ib2][ib1]
                            vxc[ik][ib1][ib2] = vxc[ik][ib2][ib1]

                elif mattype == 2:
                    for ib1 in range(nb):
                        for ib2 in range(ib1+1, nb):
                            sigma[ik][ib2][ib1] = sigma[ik][ib1][ib2]
                            vxc[ik][ib2][ib1] = vxc[ik][ib1][ib2]

    return


def diag(mat, deg=False, herm=True):  # RJW: change herm to true.
    r"""
    Diagonalize the matrix on one kpoint, and return eigenvalues and eigenvectors
    Tricks related to degeneracy and hermicity can be played here

    NOTE: we keep the order of eigenvector matrix becuase we can directly use it later
          However in rotating wavefunctions, we need to be very careful to the
          correspondence of eigenvector

    We sort according to eigenvalus, so that we have well-defined valence/conduction  bands

    deg: if do degeneracy check and wipe out off-diags
    herm: if force hermicity
    """

    if herm == True:
        mat = 0.5 * (mat + mat.transpose().conjugate())

    #eval, evec = np.linalg.eig(mat)
    eval, evec = np.linalg.eigh(mat)

    # Do sorting
    idx = eval.argsort()
    eval = eval[idx]
    evec = evec[:, idx]

    return (eval, evec)


def matrot(mat, vec):
    return np.dot(np.dot(vec.transpose().conjugate(), mat), vec)


def extrap(hammfk, hamqpk, bmap, bandmfk, scissor=False):
    r"""
    Do extrapolation based on the rotated mf and qp energies in subspace
    Return extrapolated bandqp for all bands at this kpoint

    NOTE: bmap maps to the band index starting from 1, however in python
          the index starts from 0

    This function is modified from Diana's extrapolate_g0w0.py script
    """

    nb = len(bmap)
    nbtot = len(bandmfk)

    bandqpk = np.zeros(nbtot, dtype=np.float)

    # index of one band in python is 1 smaller than the band number starting from 1
    idxmin = bmap[0] - 1
    idxmax = bmap[nb-1] - 1

    # When doing scissor, one more info is needed, the vbm/cbm at this kpoint
    if scissor == True:
        print("Scissor oprator not implemented yet. Only rigid shift method.")
        sys.exit(0)

    if scissor == False:
        dlow = hamqpk[0][0].real - hammfk[0][0].real
        dhigh = hamqpk[nb-1][nb-1].real - hammfk[nb-1][nb-1].real

        for ib in range(nbtot):
            if ib < idxmin:
                bandqpk[ib] = bandmfk[ib] + dlow
            elif ib > idxmax:
                bandqpk[ib] = bandmfk[ib] + dhigh

        for ib in range(nb):
            bandqpk[bmap[ib]-1] = hamqpk[ib][ib].real

    return bandqpk


def wfn_update(hammf, hamqp, bmap, kmap, evecs):
    r"""
    Read old wavefunction and copy to a new wavefunction
    Update new wavefunction with the rotation matrix
    Update energies in new wavefunctions
    New wavefunction file will be written

    Since inner and outer wavefunctions are forced (somehow) to be the same one,
    we don't need to find out kmap because we already got it from sigma_hp.log

    We will read out total number of bands in wavefunctions
    We will read out band energies at all available kpoints
    We will read out ngk, which includes number of G-vectors at each kpoint
    We will get nktot from ngk

    This function is modified based on Felipe's pseudobands.py script
    """

    print()
    print("Reading and copying wavefunctions ...")

    wfnold = "wfn_old.h5"
    wfnnew = "wfn_new.h5"

    f_in = h5py.File(wfnold)
    f_out = h5py.File(wfnnew, "w")

    print()
    print("Reading from", wfnold)
    print("Writing to", wfnnew)
    print()

    nbtot = f_in['mf_header/kpoints/mnband'][()]
    print("Number of bands in wavefunctions is", nbtot)

    f_out.copy(f_in['mf_header'], 'mf_header')
    f_out.create_group('wfns')
    f_out.copy(f_in['wfns/gvecs'], 'wfns/gvecs')

    f_out['mf_header/kpoints/mnband'][()] = nbtot

    print()
    print("Copying wavefunctions ...")

    # [nbtot, nspin*nspinor, ngktot, {1 if REAL, 2 if CPLX}}]
    shape = list(f_in['wfns/coeffs'].shape)
    # BAB
    print(shape)
    ngk = f_in['mf_header']['kpoints']['ngk'].value
    print((ngk[0]))
    #shape[2] = ngk[0]
    # end bAB
    shape[0] = nbtot
    nspinor = shape[1]  # assuming nspin=1, nspinor=1 or 2
    f_out.create_dataset('wfns/coeffs', shape, 'd')
    #f_out.create_dataset('wfns/coeffs', shape)
    # does this have four indices? yes-ish... looks like real/cplx is treated differently now....
    f_out['wfns/coeffs'][:nbtot, :, :, :] = f_in['wfns/coeffs'][:nbtot, :, :, :]
    # f_out['wfns/coeffs'][:nbtot,:,:] = f_in['wfns/coeffs'][:nbtot,:,:ngk[0]]

    print("Finish copying wavefunctions")
    print()

    # We now update the new wavefunction

    # BAB
    #ngk = f_in['mf_header/kpoints/ngk'][()]
    #ntotk = len(ngk)
    # end BAB

    print("Now rotating the small block of wfns and overwriting")

    nb = len(bmap)
    nk = len(kmap)

    bandmin = bmap[0]
    bandmax = bmap[nb-1]

    # BAB
    # for ik in range(1,nk):
    for ik in range(0, nk):    # RJW
        # end BAB
        for isp in range(nspinor):
            # BAB note: probably add in a loop over spinor index
            print()

            wfncoef_ik = np.zeros((nb, ngk[kmap[ik]-1]), dtype=np.complex)
            wfn_rotate = np.zeros((nb, ngk[kmap[ik]-1]), dtype=np.complex)
            numabv = get_ngk_above(ngk, kmap[ik]-1)

            for ib in range(bandmin-1, bandmax):   # python starts from 0
                wfncoef_ik[ib-(bandmin-1), :] = f_in['wfns/coeffs'][ib, isp, numabv: numabv+ngk[kmap[ik]-1], 0] \
                    + 1.0j * f_in['wfns/coeffs'][ib, isp,
                                                 numabv: numabv+ngk[kmap[ik]-1], 1]

            for ib1 in range(nb):
                for ib2 in range(nb):
                    # Pay attention to evec index order
                    wfn_rotate[ib1, :] += evecs[ik][ib2][ib1] * \
                        wfncoef_ik[ib2, :]

            for ib in range(bandmin-1, bandmax):
                f_out['wfns/coeffs'][ib, isp, numabv: numabv +
                                     ngk[kmap[ik]-1], 0] = wfn_rotate[ib-(bandmin-1), :].real
                f_out['wfns/coeffs'][ib, isp, numabv: numabv +
                                     ngk[kmap[ik]-1], 1] = wfn_rotate[ib-(bandmin-1), :].imag

    print()
    print("Extrapolating quasi-particle band energies, writing new energies into new wavefunctions")

    bandmf = np.zeros((nk, nbtot), dtype=np.float)
    bandqp = np.zeros((nk, nbtot), dtype=np.float)

    # BAB note: not anymore!
    # ispin = 0  # nspin = 1, ispin = 0 is the only index
    ryd = 13.60569193

    # BAB note: no need to loop over spinor index since this is energy  #RJW: WARNNING, spin-index dispear
    for ik in range(nk):

        for ib in range(nbtot):
            bandmf[ik][ib] = f_in['mf_header/kpoints/el'][:, kmap[ik]-1, ib] * ryd

        bandqp[ik, :] = extrap(hammf[ik], hamqp[ik], bmap, bandmf[ik])[:]

        for ib in range(nbtot):
            f_out['mf_header/kpoints/el'][:,
                                          kmap[ik]-1, ib] = bandqp[ik][ib] / ryd

    f_in.close()
    f_out.close()

    print()
    print("Updating wavefunctions completed!")
    print()

    return bandqp


def get_ngk_above(ngk, kidx):
    r"""
    kidx is the index of kpoints that is being rotated,
    and this kidx MUST BE python index, i.e. kmap - 1
    """
    num = 0
    for i in range(kidx):
        num = num + ngk[i]

    return num


def write_eqp_vxc(bandqp, hamqp, sigma, bmap, kpts):
    r"""
    This function writes eqp_outer.dat and vxc_new.dat

    In eqp_outer.dat, MF and QP energies are both set to QP energies
    In vxc_new.dat, it is actually not vxc, but rotated Sigma matrix from last GW step
    """

    nb = len(bmap)
    nk = len(kpts)
    nbtot = len(bandqp[0])

    feqp = open("eqp_outer.dat", "w")

    for ik in range(nk):
        ln = "   " + str(kpts[ik][0]) + "   " + str(kpts[ik][1]) + "   " + str(kpts[ik][2]) \
            + "      " + str(nbtot) + "\n"
        feqp.write(ln)

        for ib in range(nbtot):
            ln = "        1    " + str(ib+1) + "   " \
                + str(bandqp[ik][ib]) + "   " + str(bandqp[ik][ib]) + "\n"
            feqp.write(ln)

    feqp.close()

    fvxc = open("vxc_new.dat", "w")

    vxc = sigma - hamqp
    for ik in range(nk):
        for ib in range(nb):
            vxc[ik][ib][ib] += hamqp[ik][ib][ib].real

    for ik in range(nk):
        ln = "   " + str(kpts[ik][0]) + "   " + str(kpts[ik][1]) + "   " + str(kpts[ik][2]) \
            + "      " + str(nb) + "      " + str(nb*nb) + "\n"
        fvxc.write(ln)

        for ib in range(nb):
            ln = "        1    " + str(bmap[ib]) + "   " \
                + str(vxc[ik][ib][ib].real) + "   " + \
                str(vxc[ik][ib][ib].imag) + "\n"
            fvxc.write(ln)

        for ib1 in range(nb):
            for ib2 in range(nb):
                ln = "        1    " + str(bmap[ib1]) + "   " + str(bmap[ib2]) + "   " \
                    + str(vxc[ik][ib1][ib2].real) + "   " + \
                    str(vxc[ik][ib1][ib2].imag) + "\n"
                fvxc.write(ln)

    fvxc.close()

    return


def scgwtool():
    r"""
    This is the major function that controls the flow
    """

    print("Start scGWtool.py")

    # Read sigma_hp.log, and initialize matrices
    print()
    print()
    print("============================== Step 1/5 ==============================")
    print("==                                                                  ==")
    print("==  Reading sigma_hp files and constructing Hqp sub-space matrices  ==")
    print("==                                                                  ==")
    print("======================================================================")
    print()
    (hammf0, hamqp0, sigma0, bmap, kmap, kpts) = read_sigma_hp()

    print("Reading sigma0 is done ")

    nk = len(kmap)
    nb = len(bmap)

    evals = np.zeros((nk, nb), dtype=np.complex)
    evecs = np.zeros((nk, nb, nb), dtype=np.complex)

    hamqp = np.zeros((nk, nb, nb), dtype=np.complex)
    hammf = np.zeros((nk, nb, nb), dtype=np.complex)
    sigma = np.zeros((nk, nb, nb), dtype=np.complex)

    # diagonalize at each kpoint
    print()
    print("============================== Step 2/5 ==============================")
    print("==                                                                  ==")
    print("==                      Diagonalizing matrices                      ==")
    print("==                                                                  ==")
    print("======================================================================")
    print()

    for ik in range(nk):
        (evals[ik], evecs[ik]) = diag(hamqp0[ik])

    # rotate matrices
    print()
    print("============================== Step 3/5 ==============================")
    print("==                                                                  ==")
    print("==              Rotating Hmf, Hqp, and Sigma matrices               ==")
    print("==                                                                  ==")
    print("======================================================================")
    print()

    for ik in range(nk):
        hammf[ik] = matrot(hammf0[ik], evecs[ik])
        hamqp[ik] = matrot(hamqp0[ik], evecs[ik])
        sigma[ik] = matrot(sigma0[ik], evecs[ik])

    # wavefunctions
    print()
    print("============================== Step 4/5 ==============================")
    print("==                                                                  ==")
    print("==                     Updating wavefunctions                       ==")
    print("==                                                                  ==")
    print("======================================================================")
    print()

    bandqp = wfn_update(hammf, hamqp, bmap, kmap, evecs)

    # write eqp_outer.dat and vxc_new.dat
    print()
    print("============================== Step 5/5 ==============================")
    print("==                                                                  ==")
    print("==              Writing eqp_outer.dat and vxc_new.dat               ==")
    print("==                                                                  ==")
    print("======================================================================")
    print()

    write_eqp_vxc(bandqp, hamqp, sigma, bmap, kpts)

    print("scGWtool.py finished.")

    return


def main():

    scgwtool()

    return


if __name__ == "__main__":
    main()
