/******************************************************************************
*
*  libtile_qhull_find_h    Originally by FHJ     Last Modified 04/18/2013 (FHJ)
*
*    BerkeleyGW wrapper for the Qhull tessellation routines.
*    This file implements routines to find the Delaunay simplex that encloses
*    a particular point. The code is based on scipy`s wrapper for qhull:
*    https://github.com/scipy/scipy/blob/master/scipy/spatial/qhull.pyx
*
*    This file is part of the BerkeleyGW package.
*
******************************************************************************/

#ifndef libtile_qhull_find_h
#define libtile_qhull_find_h

#define TOL_INSIDE 1e-10
#define MAX_VERTEX_NEIGHBORS 32

/* These two arrays are initialized by get_simplex_facet_array() */

/** (nsimplex, ndim+1) in C order
 *
 * Indices of the points forming the simplices in the triangulation. */
extern int *simplices;

/** (nsimplex, ndim+1) in C order
 *
 * Indices of neighbor simplices for each simplex.
 * The kth neighbor is opposite to the kth vertex.
 * For simplices at the boundary, -1 denotes no neighbor. */
extern int *neighbors;

extern int nsimplex; /* Number of simplicial facets = nsimplex */

/** (nsimplex, ndim+1, ndim) in C order
 *
 * Affine transform from ``x`` to the barycentric coordinates ``c``.
 *      This is defined by::
 *
 *          T c = x - r
 *
 *      At vertex ``j``, ``c_j = 1`` and the other coordinates zero.
 *
 *      For simplex ``i``, ``transform[i,:ndim,:ndim]`` contains
 *      inverse of the matrix ``T``, and ``transform[i,ndim,:]``
 *      contains the vector ``r``.
 *
 *      If the simplex is degenerate or nearly degenerate, its
 *      barycentric transform contains NaNs. */
extern double *transform;

/* Array with sets containing the second-nearest-neighbor vertices for a given
 * lower Delaunay simplex. The vertices in the set first_neighbors are not
 * included. */
extern setT **second_neighbors;
/* Maximum number of second neighbors; -1 if second neighbors not initialized
 * (via init_second_neighbors).
 */
extern int max_num_second_neighbors;


int init_simplex_facet_array();
int free_simplex_facet_array();
int init_barycentric_transforms();
int free_barycentric_transforms();
/*
int find_simplex_bruteforce(const coordT *x, double *coef);
int find_simplex_directed(const coordT *x, double *coef);
*/
int find_simplex(const coordT *x, double *coef);
int init_second_neighbors();
int get_second_neighbors(int isimplex, int *indices, int *num_second);
int free_second_neighbors();
#endif
