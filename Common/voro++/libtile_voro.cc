/******************************************************************************
*
*  libtile_voro           Originally by FHJ     Last Modified 08/18/2014 (FHJ)
*
*    BerkeleyGW wrapper for the voro++ tessellation routines.
*
*    This file is part of the BerkeleyGW package.
*
******************************************************************************/

#include <cstdio>
#include "src/voro++.hh"
#include "libtile_voro.h"

extern "C" {
int voro_pts2bz(double *M, double *pts_cart, int *npts, int *umks) {
#if defined(DEBUG) && defined(VERBOSE)
	printf("Lattice Vectors:\n");
	printf("v1: %f %f %f\n", M[0], 0., 0.);
	printf("v2: %f %f %f\n", M[1], M[2], 0.);
	printf("v3: %f %f %f\n", M[3], M[4], M[5]);
	printf("\n");
#endif

	/* M = {xx, xy, yy, xz, yz, zz} */
	voro::container_periodic con(M[0], M[1], M[2], M[3], M[4], M[5], 1,1,1,8);
	con.put(0,0,0,0);

	double rx, ry, rz, *pp;
	int gx, gy, gz, i, ipt, *umk;
	for (ipt=0, pp=pts_cart, umk=umks; ipt<*npts; ipt++, pp+=3, umk+=3) {
		if(!con.find_voronoi_cell(pp[0],pp[1],pp[2],rx,ry,rz,i)) {
			fprintf(stderr,"ERROR: find_voronoi_cell error for point # %d\n", ipt);
			return -1;
		}
		gz = (int) round(rz/M[5]);
		gy = (int) round((ry-gz*M[4])/M[2]);
		gx = (int) round((rx-gy*M[1]-gz*M[3])/M[0]);
		umk[0] = -gx;
		umk[1] = -gy;
		umk[2] = -gz;
	}

	return 0;
}

int voro_pts2bz_(double *M, double *pts_cart, int *npts, int *umks) {
	return voro_pts2bz(M, pts_cart, npts, umks);
}


int voro_get_kpts_volumes(double *M, double *pts_cart, int *npts, double *vols) {
#if defined(DEBUG) && defined(VERBOSE)
	printf("Lattice Vectors:\n");
	printf("v1: %f %f %f\n", M[0], 0., 0.);
	printf("v2: %f %f %f\n", M[1], M[2], 0.);
	printf("v3: %f %f %f\n", M[3], M[4], M[5]);
	printf("\n");
#endif

	/* M = {xx, xy, yy, xz, yz, zz} */
	voro::container_periodic con(M[0], M[1], M[2], M[3], M[4], M[5], 1,1,1,8);
	int ipt;
	double *pp;
	for (ipt=0, pp=pts_cart; ipt<*npts; ipt++, pp+=3){
		con.put(ipt, pp[0], pp[1], pp[2]);
	}

	for (ipt=0; ipt<*npts; ipt++) vols[ipt] = 0.;
	voro::c_loop_all_periodic clap(con);
	voro::voronoicell c;
	if (clap.start()) {
		do if (con.compute_cell(c,clap)){
			vols[clap.pid()] = c.volume();
		} while (clap.inc());
	} else {
		fprintf(stderr,"ERROR: Could not loop over voronoi cells!\n");
		return -1;
	}
	return 0;
}

int voro_get_kpts_volumes_(double *M, double *pts_cart, int *npts, double *vols) {
	return voro_get_kpts_volumes(M, pts_cart, npts, vols);
}
}
