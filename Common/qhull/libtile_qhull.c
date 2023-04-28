/******************************************************************************
*
*  libtile_qhull           Originally by FHJ     Last Modified 12/04/2012 (FHJ)
*
*    BerkeleyGW wrapper for the Qhull tessellation routines.
*    Right now, only the Delaunay triangulation routines is implemented.
*
*    This file is part of the BerkeleyGW package.
*
******************************************************************************/

#include <stdarg.h>
#include "qhull_a.h"
#include "libtile_qhull.h"
#include "libtile_qhull_find.h"

int _inode;

void qhull_info(const char *message, ...)
{
	va_list args;

	va_start(args, message);
	fprintf(stdout, "\nINFO [libtile_qhull @ %d]: ", _inode);
	vfprintf(stdout, message, args);
	fprintf(stdout, "\n\n");
	va_end (args);
}
void qhull_warn(const char *message, ...)
{
	va_list args;

	va_start(args, message);
	fprintf(stderr, "\nWARNING [libtile_qhull @ %d]: ", _inode);
	vfprintf(stderr, message, args);
	fprintf(stderr, "\n\n");
	va_end (args);
}
void qhull_error(const char *message, ...)
{
	va_list args;

	va_start(args, message);
	fprintf(stderr, "\nERROR [libtile_qhull @ %d]: ", _inode);
	vfprintf(stderr, message, args);
	fprintf(stderr, "\n\n");
	va_end (args);
}


/**
 * Calculates the Delaunay triangulation/tetrahedralization of a set of points.
 * Note: right now, there can be only one active Delaunay triangulation at
 * a given time.
 * 
 * @param points [inout] Fortran array (ndim, npts)
 * @param npts [in] number of points.
 * @param ndim [in] number of dimensions
 */
int qhull_init_delaunay_(double points[], const int *npts, const int *ndim, const int *inode)
{
	boolT ismalloc= False;	 /* don`t free points in qh_freeqhull() */
	char flags[250];
	FILE *outf=NULL;
	FILE *errf=stderr;
	int ierr;

	_inode = *inode;
	sprintf(flags, QHULL_STR);
	ierr = qh_new_qhull(*ndim, *npts, points, ismalloc, flags, outf, errf);
	if (ierr) {
		qhull_error("could not create Delaunay triangulation (code %d)", ierr);
		return ierr;
	}
	ierr = init_simplex_facet_array();
	if (ierr) {
		qhull_error("could not get simplicial facets (code %d)", ierr);
		return ierr;
	}
	ierr = init_barycentric_transforms();
	if (ierr) {
		qhull_error("could not get barycentric transformations (code %d)", ierr);
		return ierr;
	}
	return ierr;
}
int qhull_init_delaunay(double points[], const int *npts, const int *ndim, const int *inode)
{
	return qhull_init_delaunay_(points, npts, ndim, inode);
}


/** 
 * Frees buffers associated to the Delaunay triangulation.
 */
int qhull_free_delaunay_()
{
	int curlong, totlong;      /* mem remaining after qh_memfreeshort */

	free_second_neighbors();
	qh_freeqhull(!qh_ALL);                 /* long mem */
	qh_memfreeshort (&curlong, &totlong);  /* short mem/mem allocator */
	if (curlong || totlong) {
		qhull_warn("did not free %d bytes of long memory (%d pieces)\n", 
							totlong, curlong);
	}
	free_simplex_facet_array();
	free_barycentric_transforms();
	return curlong || totlong;
}
int qhull_free_delaunay()
{
	return qhull_free_delaunay_();
}


/**
 * Finds the Delaunay triangle/tetrahedron that encloses `point`.
 *
 * @param point [in] array (ndim).
 * @param indices [out] array (ndim+1) of the (Fortran) indices of the vertices.
 * @param coefs [out] coefficients of point in barycentric coordinates.
 */
int qhull_find_delaunay_simplex_(const double point[], int indices[], double coefs[])
{
	int isimplex, i;
	int *verts;
	const int ndim = qh input_dim;

	isimplex = find_simplex(point, coefs);
	if (isimplex<0) {
		for (i=0; i<ndim+1; i++) indices[i] = 0;
		qhull_warn("could not find simplex");
		/* XXX: should we return the closest simplex instead? */
		return -1;
	} else {
		verts = simplices + (ndim+1)*isimplex;
		for (i=0; i<ndim+1; i++) indices[i] = verts[i]+1;
		return isimplex;
	}
}

int qhull_find_delaunay_simplex(const double point[], int indices[], double coefs[])
{
	return qhull_find_delaunay_simplex_(point, indices, coefs);
}

/**
 * Returns the total number of simplices obtained from the Delaunay
 * triangulation.
 */
int qhull_get_num_simplices_(int *num_simplices)
{
	*num_simplices = nsimplex;
	return 0;
}

int qhull_get_num_simplices(int *num_simplices)
{
	return qhull_get_num_simplices_(num_simplices);
}

/**
 * Returns the indices of the points that define each simplex.
 */
int qhull_get_simplices_(int indices[])
{
	const int ndim = qh input_dim;
	int ii;
	for (ii=0; ii<nsimplex*(ndim+1); ii++) indices[ii] = simplices[ii]+1;
	return 0;
}

int qhull_get_simplices(int indices[])
{
	return qhull_get_simplices_(indices);
}

/**
 * Returns the indices of the neighbors for each simplex.
 */
int qhull_get_neighbors_(int nei[])
{
	const int ndim = qh input_dim;
	int ii;
	for (ii=0; ii<nsimplex*(ndim+1); ii++) nei[ii] = neighbors[ii]+1;
	return 0;
}

int qhull_get_neighbors(int nei[])
{
	return qhull_get_neighbors_(nei);
}

/*
 * Routines to compute second nearest neighbors
 */

/*
 * Computes second-nearest neighbors and barycentric coefficients for point `point`.
 * The buffers indices and coefs must be ndim+1 + `max_num_second_neighbors` big.
 * The first ndim+1 points will correspond to the first neighbors.
 * nnei is set to the total number of neighbors found, with
 * nnei <= ndim + 1 + max_num_second_neighbors
 */
int qhull_find_delaunay_simplex_with_second_neighbors_(const double point[], int indices[], double coefs[], int *nnei)
{
	const int ndim = qh input_dim;
	int ipt, isimplex, ierr;

	if (max_num_second_neighbors==-1) {
		qhull_error("second neighbors not initialized");
		return -1;
	}
	for (ipt=0; ipt < ndim + 1 + max_num_second_neighbors; ipt++) {
		indices[ipt] = 0;
		coefs[ipt] = 0.;
	}
	isimplex = qhull_find_delaunay_simplex_(point, indices, coefs);
	if (isimplex<0) return isimplex;
	// Only call get_second_neighbors for the indices associated with the second neighbors.
	ierr = get_second_neighbors(isimplex, indices+ndim+1, nnei);
	if (ierr<0) return ierr;
	// nnei is the *total* number of neighbors!
	*nnei = *nnei + ndim + 1;
	// Convert from C to Fortran indicing:
	for (ipt=ndim+1; ipt<(*nnei); ipt++) indices[ipt]++;
	return isimplex;
}

int qhull_find_delaunay_simplex_with_second_neighbors(const double point[], int indices[], double coefs[], int *nnei)
{
	return qhull_find_delaunay_simplex_with_second_neighbors_(point, indices, coefs, nnei);
}

int qhull_init_second_neighbors()
{
	return init_second_neighbors();
}

int qhull_init_second_neighbors_()
{
	return qhull_init_second_neighbors();
}

int qhull_free_second_neighbors()
{
	return free_second_neighbors();
}

int qhull_free_second_neighbors_()
{
	return qhull_free_second_neighbors();
}

int qhull_get_max_num_second_neighbors(int *max_nei)
{
	*max_nei = max_num_second_neighbors;
	return 0;
}

int qhull_get_max_num_second_neighbors_(int *max_nei)
{
	return qhull_get_max_num_second_neighbors(max_nei);
}
