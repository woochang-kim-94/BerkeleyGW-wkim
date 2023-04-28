/******************************************************************************
*
*  libtile_voro           Originally by FHJ     Last Modified 08/18/2014 (FHJ)
*
*    BerkeleyGW wrapper for the voro++ tessellation routines.
*
*    This file is part of the BerkeleyGW package.
*
******************************************************************************/

#ifndef libtile_voro_h
#define libtile_voro_h
extern "C" {
int voro_pts2bz(double *M, double *pts_cart, int *npts, int *umks);
int voro_pts2bz_(double *M, double *pts_cart, int *npts, int *umks);
int voro_get_kpts_volumes(double *M, double *pts_cart, int *npts, double *vols);
int voro_get_kpts_volumes_(double *M, double *pts_cart, int *npts, double *vols);
}
#endif
