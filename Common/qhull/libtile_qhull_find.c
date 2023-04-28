/******************************************************************************
*
*  libtile_qhull_find      Originally by FHJ     Last Modified 04/18/2013 (FHJ)
*
*    BerkeleyGW wrapper for the Qhull tessellation routines.
*    This file implements routines to find the Delaunay simplex that encloses
*    a particular point. The code is based on scipy`s wrapper for qhull:
*    https://github.com/scipy/scipy/blob/master/scipy/spatial/qhull.pyx
*
*    This file is part of the BerkeleyGW package.
*
******************************************************************************/

#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "qhull_a.h"
#include "libtile_qhull.h"
#include "libtile_qhull_find.h"
#include "qhull_blas.h"

#if defined(DEBUG) && defined(VERBOSE)
#define DEBUG_VERBOSE
#endif

int nsimplex;
int *simplices;
int *neighbors;
double *transform;
static unsigned char *transf_valid;
static int last_simplex;
static int first_valid;
/* Maps global list of facets (qh facet->id) to lower Delaunay simplices */
static int *id_map;
/* Array with sets containing the vertices defining a given lower Delaunay
 * simplex. */
//setT **first_neighbors;
/* Array with sets containing only the second-nearest-neighbor vertices for a
 * given lower Delaunay simplex. */
setT **second_neighbors;
int max_num_second_neighbors = -1;

/**
  * Dumps the facet and neighbor list from Qhull into arrays simplices and
  * neighbors. Note that some simplices are degenerate, which are taken care
  * by init_barycentric_transforms and fix_neighbors.
  */
int init_simplex_facet_array()
{
	facetT* facet;
	facetT* neighbor;
	vertexT *vertex;
	int i, j;
	int facet_ndim;
	const int ndim = qh input_dim;

	facet_ndim = ndim + 1;
	last_simplex = -1;

	/* Note: qh facet_id = id of the next facet to be created = number of facets */
	id_map = (int*) malloc((qh facet_id)*sizeof(int));
	for (i=0; i<qh facet_id; i++) id_map[i] = -1;

	j = 0;
	FORALLfacets {
		/* Only include lower Delaunay facets in our list */
		if (facet->upperdelaunay == qh UPPERdelaunay) {
			/* Make sure all lower Delaunay facets are simplicial */
			if ( (!facet->simplicial) && \
			     ((qh_setsize(facet->vertices) != facet_ndim) || \
			      (qh_setsize(facet->neighbors) != facet_ndim)) ) {
				qhull_error("non-simplical facet encountered: %d vertices", \
					qh_setsize(facet->vertices));
				return -1;
			}
			id_map[facet->id] = j;
			j++;
		}
	}
	nsimplex = j;

#ifdef VERBOSE
	if (_inode==0) {
		fflush(stdout);
		qhull_info("convex hull has %d facets, %d of which are lower Delaunay.", \
			qh num_facets, nsimplex);
		fflush(stdout);
	}
#endif

	simplices = (int*) malloc(nsimplex*facet_ndim*sizeof(int));
	neighbors = (int*) malloc(nsimplex*facet_ndim*sizeof(int));

	j = 0;
	FORALLfacets {
		/* Make sure the facet is lower Delaunay */
		if (facet->upperdelaunay != qh UPPERdelaunay) continue;

		/* Save vertex info */
		for (i=0; i<facet_ndim; i++) {
			vertex = facet->vertices->e[i].p;
			simplices[j*facet_ndim + i] = qh_pointid(vertex->point);
		}

		/* Save neighbor info */
		for (i=0; i<facet_ndim; i++) {
			neighbor = facet->neighbors->e[i].p;
			neighbors[j*facet_ndim + i] = id_map[neighbor->id];
		}

		j++;
	}
	return 0;
}

int free_simplex_facet_array()
{
	free(id_map);
	free(simplices);
	free(neighbors);
	return 0;
}

/**
  * Qhull often outputs simplices which are degenerate (i.e., have zero or
  * almost zero volume). This function makes sure each simplex only has
  * non-degenerate neighbors, which decreases the need for brute force searches
  * later on.
  */
static void fix_neighbors()
{
	int isimplex, jsimplex, *cur_vertices;
	const int ndim = qh input_dim;
	int facet_ndim = ndim + 1;
	int i, j, ir, jr, iv, tmp;
	int *sorted_vertices, match;
	int *problematic, nprob, iprob, jprob, nfixed;

#ifdef VERBOSE
	if (_inode==0) {
		fflush(stdout);
		qhull_info("fixing the neighbors list.");
	}
#endif
	nfixed = 0;
	/* List of problematic simplices */
	problematic = (int*) malloc(nsimplex*sizeof(int));
	nprob = 0;
	/* Loop over all ridges, and create a list of sorted vertices for each ridge.
	We will use that to simplify the search for neighbors. */
	sorted_vertices = (int*) malloc(nsimplex*facet_ndim*ndim*sizeof(int));
	cur_vertices = sorted_vertices;
	for (isimplex=0; isimplex<nsimplex; isimplex++) {
		if (!transf_valid[isimplex]) continue;
		/* Figure out if this simplex is problematic. */
		for (ir=0; ir<facet_ndim; ir++) {
			jsimplex = neighbors[isimplex*facet_ndim + ir];
			if (jsimplex!=-1) {
				if (!transf_valid[jsimplex]) {
					problematic[nprob++] = isimplex;
					break;
				}
			}
		}
		if (ir==facet_ndim) continue;

		/* Loop over all ridges from simplex isimplex */
		cur_vertices = sorted_vertices + isimplex*ndim*facet_ndim;
		for (ir=0; ir<facet_ndim; ir++) {
			/* Construct the vertex list for ridge ir, which
			contains ndim vertices */
			tmp = 0;
			for (iv=0; iv<facet_ndim; iv++) {
				if (iv==ir) continue;
				cur_vertices[tmp] = simplices[isimplex*facet_ndim + iv];
				tmp++;
			}
			/* Sort that list */
			for (i=1; i<ndim; i++) {
				j = i;
				while (j>0 && cur_vertices[j] < cur_vertices[j-1]) {
					tmp = cur_vertices[j];
					cur_vertices[j] = cur_vertices[j-1];
					cur_vertices[j-1] = tmp;
					j--;
				}
			}
			cur_vertices += ndim;
		}
	}

	/* Loop over all problematic simplices */
	for (iprob=0; iprob<nprob; iprob++){
		isimplex = problematic[iprob];
		/* Loop over all ridges from simplex isimplex */
		for (ir=0; ir<facet_ndim; ir++) {
			/* Only consider ridges that point to invalid simplices */
			tmp = neighbors[isimplex*facet_ndim + ir];
			if (tmp==-1) continue;
			if (transf_valid[tmp]) continue;
			/* Search all other problematic simplices that share the same ridge */
			for (jprob=iprob+1; jprob<nprob; jprob++){
				jsimplex = problematic[jprob];
				/* Loop over all ridges from simplex jsimplex */
				for (jr=0; jr<facet_ndim; jr++) {
					/* See if both ridges match */
					match = 1;
					for (iv=0; iv<ndim; iv++) {
						if (sorted_vertices[isimplex*facet_ndim*ndim + ir*ndim + iv]\
						!= sorted_vertices[jsimplex*facet_ndim*ndim + jr*ndim + iv]) {
							match = 0;
							break;
						}
					}
					if (match) {
						 neighbors[isimplex*facet_ndim + ir] = jsimplex;
						 neighbors[jsimplex*facet_ndim + jr] = isimplex;
						 nfixed++;
						 break;
					}
				}
				if (match) break;
			}
		}
	}

	free(sorted_vertices);
	free(problematic);
#ifdef VERBOSE
	if (_inode==0) {
		qhull_info("fixed %d neighbors.", nfixed);
		fflush(stdout);
	}
#endif
}

static double matrix_norm1(int n, double *a)
{
	double colsum, maxsum = 0.0;
	int i,j;

	for (i=0; i<n; i++) {
		colsum = 0.0;
		for (j=0; j<n; j++) {
			colsum += fabs(*(a++));
		}
		maxsum = (maxsum>colsum)? maxsum : colsum;
	}
	return maxsum;
}

/**
 * Calculates all barycentric transforms for all simplices. This also finds out
 * which simplices are degenerate (i.e., have zero or almost zero volume),
 * calls fix_neighbors to fix neighbors list.
 */
int init_barycentric_transforms()
{
	coordT *all_coords;
	double *T, *cur_transf;
	int i, j, isimplex, last_vertex, jvertex, *cur_vertices;
	double *last_coord;
	const int ndim = qh input_dim;
	int n, nrhs, lda, ldb, *ipiv, *iwork, info;
	double *work, anorm, rcond;
	int ii, jj;
	/* This guarantees that ~three digits must be correct. */
	double rcond_limit = 1000.0 * TOL_INSIDE;
	int non_sing_transf = 0;

	T = (double*) calloc(ndim*ndim, sizeof(double));
	transform = (double*) calloc(nsimplex*(ndim+1)*(ndim), sizeof(double));
	transf_valid = (unsigned char*) calloc(nsimplex, sizeof(unsigned char));
	ipiv = (int*) malloc(ndim*sizeof(int));
	iwork = (int*) malloc(ndim*sizeof(int));
	work = (double*) malloc(((ndim<4)?4:ndim)*ndim*sizeof(double));
	all_coords = qh first_point; /* (npts,ndim) */
	n = nrhs = lda = ldb = ndim;
	first_valid = -1;

#ifdef VERBOSE
	if (_inode==0) {
		fflush(stdout);
		qhull_info("calculating barycentric transformations X via A*X=B.");
	}
#endif

	cur_transf = transform; /* T^{-1} for current isimplex */
	cur_vertices = simplices; /* list of vertices id for current isimplex */
	for (isimplex=0; isimplex<nsimplex; isimplex++) {
		last_vertex = cur_vertices[ndim];
		last_coord = all_coords + last_vertex*(ndim+1);
		for (i=0; i<ndim; i++) { /* loop over dimensions */
			cur_transf[i*ndim + i] = 1.0;
			cur_transf[ndim*ndim + i] = *last_coord;
			for (j=0; j<ndim; j++) { /* loop over vertices (except last one) */
				jvertex = cur_vertices[j];
				T[i*ndim + j] = all_coords[jvertex*(ndim+1) + i] - (*last_coord);
			}
			last_coord++;
		}

#ifdef DEBUG_VERBOSE
		if (_inode==0) {
			printf("\nSIMPLEX: %d\n", isimplex);
			printf("Original vectors =\n");
			for (ii=0; ii<ndim; ii++) {
				printf("  ");
				for (jj=0; jj<ndim+1; jj++) {
					jvertex = cur_vertices[jj];
					printf(" %11.6f",all_coords[jvertex*(ndim+1) + ii]);
				}
				printf("\n");
			}
			printf("A =\n");
			for (ii=0; ii<ndim; ii++) {
				printf("  ");
				for (jj=0; jj<ndim; jj++) {
					printf(" %11.6f", T[jj*ndim + ii]);
				}
				printf("\n");
			}
		}
#endif

		/* LU decomposition */
		qh_dgetrf(&n, &n, T, &lda, ipiv, &info);

		/* Check condition number */
		if (!info) {
			anorm = matrix_norm1(ndim, T);
			qh_dgecon("1", &n, T, &lda, &anorm, &rcond,
				work, iwork, &info);
			if (rcond < rcond_limit) info = 1;
		}

		if (!info) {
			/* Compute transform */
			qh_dgetrs("N", &n, &nrhs, T, &lda, ipiv, cur_transf, &ldb, &info);
		}

		if (!info) {
#ifdef DEBUG_VERBOSE
			if (_inode==0) {
				printf("X =\n");
				for (ii=0; ii<ndim; ii++) {
					printf("  ");
					for (jj=0; jj<ndim+1; jj++) {
						printf(" %11.6f", cur_transf[jj*ndim + ii]);
					}
					printf("\n");
				}
				printf("1/(condition number) = %g\n", rcond);
			}
#endif
			transf_valid[isimplex] = 1;
			non_sing_transf++;
			if (first_valid==-1) first_valid = isimplex;
#ifdef DEBUG_VERBOSE
		} else if (_inode==0) {
			printf("Transformation is singular for simplex %d.\n", isimplex);
#endif
		}
		cur_transf += ndim*(ndim+1);
		cur_vertices += ndim+1;
	}

#ifdef VERBOSE
	if (_inode==0) {
		printf("\nThere are %d non-singular Delaunay simplices.\n", non_sing_transf);
		qhull_info("finished calculating barycentric transformations.");
		fflush(stdout);
	}
#endif

	free(work);
	free(iwork);
	free(ipiv);
	free(T);
	fix_neighbors();
	if (first_valid==-1) {
		qhull_error("all simplices are degenerate.");
		return 1;
	}
	return 0;
}

int free_barycentric_transforms()
{
	free(transf_valid);
	free(transform);
	return 0;
}

/**
 * Finds the barycentric coordinates `c` for point `x` using the affine 
 * transformation `cur_transf`. Only the component `c[i]` is calculated.
 *
 * @param ndim [in] number of input dimensions (=qh input_dim).
 * @param cur_transf [in] affine transformation for the current simplex.
 * @param x [in] coordinate of the fine point.
 * @param c [out] barycentric coordinates for point `x`.
 * @param i [in] which component of `c` should be computed.
 */
static void get_barycentric_coord_single(const int ndim, \
	const double *cur_transf, const coordT *x, double *c, const int i)
{
	int j;

	if (i==ndim) {
		c[ndim] = 1.0;
		for (j=0; j<ndim; j++) c[ndim] -= c[j];
    	} else {
		c[i] = 0.0;
		for (j=0; j<ndim; j++) {
			c[i] += cur_transf[ndim*i + j] * (x[j] - cur_transf[ndim*ndim + j]);
		}
	}
}

/**
 * Determines if point `x` is inside the simplex determined by the affine
 * transformation `cur_transf`. In case of success, all barycentric coordinates
 * are stored in `c`.
 *
 * @param ndim [in] number of input dimensions (=qh input_dim).
 * @param cur_transf [in] affine transformation for the current simplex.
 * @param x [in] coordinate of the fine point.
 * @param c [out] barycentric coordinates for point `x`.
 * @return 1 if the point is inside, 0 if it`s outside.
 */
static int is_barycentric_inside(const int ndim, const double *cur_transf, const coordT *x, double *c)
{
	int i, j;

	c[ndim] = 1.0;
	for (i=0; i<ndim; i++) {
		c[i] = 0.0;
		for (j=0; j<ndim; j++) {
			c[i] += cur_transf[ndim*i + j] * (x[j] - cur_transf[ndim*ndim + j]);
		}
		if ( (c[i]<-TOL_INSIDE) || (c[i]>1+TOL_INSIDE) ) return 0;
		c[ndim] -= c[i];
	}
	if ( (c[ndim]<-TOL_INSIDE) || (c[ndim]>1+TOL_INSIDE) ) return 0;
	return 1;
}

/**
 * Finds the simplex enclosing point `x` using a brute force method.
 *
 * @param x [in] coordinate of the fine point.
 * @param c [out] barycentric coordinates for point `x`.
 * @return local index of lower Delaunay simplex, -1 for failure.
 */
static int find_simplex_bruteforce(const coordT *x, double *c)
{
	int isimplex;
	const double* cur_transf;
	const int ndim = qh input_dim;
	const unsigned char *cur_valid;

#ifdef VERBOSE
	qhull_warn("switching to brute force method.");
#endif
	cur_transf = transform;
	cur_valid = transf_valid;
	for (isimplex=0; isimplex<nsimplex; isimplex++) {
		if ( *cur_valid ) {
			if ( is_barycentric_inside(ndim, cur_transf, x, c) ) return isimplex;
		}
		cur_transf += ndim*(ndim+1);
		cur_valid++;
	}
	return -1;
}

/**
 * Finds the simplex enclosing point `x` using a downhill algorithm.
 * If everything goes wrong, a brute force method is employed.
 *
 * @param x [in] coordinate of the fine point.
 * @param c [out] barycentric coordinates for point `x`.
 * @return local index of lower Delaunay simplex, -1 for failure.
 */
static int find_simplex_directed(const coordT *x, double *c)
{
	int k, simplex_next, inside, isimplex, cycle_k;
	int ok;
	const int ndim = qh input_dim;
	double *cur_transf;

	isimplex = last_simplex;
	last_simplex = -1;

	if ((isimplex<0)||(isimplex>=nsimplex)) {
		isimplex = first_valid;
	}

	/*The maximum iteration count: it should be large enough so that
	the algorithm usually succeeds, but smaller than nsimplex so
	that for the cases where the algorithm fails, the main cost
	still comes from the brute force search.*/
	for (cycle_k=0; (cycle_k<1+nsimplex/4)&&(isimplex!=-1); cycle_k++) {
#ifdef DEBUG_VERBOSE
		qhull_info("cycle_k=%d, isimplex=%d, valid=%d", \
			cycle_k, isimplex, transf_valid[isimplex]);
#endif
		/* If the current simplex is degenerate, hop to a neighbor */
		if (!transf_valid[isimplex]) {
			ok = 0;
			for (k=0; k<ndim+1; k++) {
				simplex_next = neighbors[(ndim+1)*isimplex + k];
				if (simplex_next!=last_simplex) {
					last_simplex = isimplex;
					isimplex = simplex_next;
					ok = 1;
					break;
				}
			}
			if (!ok) {
				/* If we got here, there is no valid neighbor.
				Let`s just try to go somewhere else then... */
#ifdef DEBUG_VERBOSE
				qhull_info("could`t find any neighbor for simplex %d.",	isimplex);
#endif
				last_simplex = isimplex;
				isimplex++;
			}
			continue;
		}

		cur_transf = transform + isimplex*ndim*(ndim+1);
		inside = 1;
		for (k=0; k<ndim+1; k++){
			get_barycentric_coord_single(ndim, cur_transf, x, c, k);
			if (c[k] < -TOL_INSIDE) {
				simplex_next = neighbors[(ndim+1)*isimplex + k];
				if (simplex_next==-1) {
					/* The point is probably outside the Delaunay convex hull */
#ifdef DEBUG_VERBOSE
					qhull_info("cycle_k=%d, isimplex=%d, valid=%d\
	It seems like the point is outside the Delaunay convex hull.\n\
	Falling back to the brute force method to be sure.",\
	cycle_k, isimplex, transf_valid[isimplex]);
#endif
					inside = 0;
					break;
				}
				last_simplex = isimplex;
				isimplex = simplex_next;
				inside = -1;
				break;
			} else if (c[k] > 1.0 + TOL_INSIDE ) {
				/* There must be at least one negative coefficient.
				Let`s try to find the corresponding neighbor.*/
				inside = 0;
				continue;
			}
		}

		if (inside==1) {
			/* Great success! */
#ifdef DEBUG_VERBOSE
			qhull_info("found the right simplex.");
#endif
			last_simplex = isimplex;
			return isimplex;
		} else if (inside==0) {
			/* We've failed utterly. Fall back to brute force. */
			break;
		}
		/* if we got here, we hopped to another simplex */
	}
	/* If we got here, we couldn`t find the right simplex. Switch to
	brute force method. */

	last_simplex = find_simplex_bruteforce(x, c);
	return last_simplex;
}


/**
 * Finds the simplex enclosing point `x` and its barycentric coordinates.
 *
 * @param x [in] coordinate of the fine point.
 * @param coefs [out] barycentric coordinates for point `x`.
 * @return local index of lower Delaunay simplex, -1 for failure.
 */
int find_simplex(const coordT *x, double *coefs)
{
#if 1
	return find_simplex_directed(x, coefs);
#else
	return find_simplex_bruteforce(x, coefs);
#endif
}

/**
  * Precompute vertex second neighbor list.
  */
int init_second_neighbors()
{
	int npts = qh num_points;
	int isimplex, iv, jv, ivert, jvert;
	int facet_ndim;
	const int ndim = qh input_dim;
	/* Number of nverts in a simplex given a dimensionality:
	 * 1 -> 1
	 * 2 -> 3
	 * 3 -> 6
	 * ndim -> nverts = ndim*(ndim+1)/2
	 */
	const int nverts = (ndim*(ndim+1))/2;
	const int buf_vert = 9, buf_simpl = 16;
	facet_ndim = ndim + 1;
	int cur_vertex, last_vertex, sz;
	facetT *facet;
	intptr_t iaddr, jaddr;
	setT *set;
	setT *first_neighbors;
	setT **vertex_edges;

	/* Loop over all simplices isimplex. Double loop over all vertices iv
	 * and jv that make up each simplex, iv < jv. For each vertex iv, add vertex
	 * jv to the list of connected vertices, and vice-versa. Use qhull's setT
	 * structure to avoid double counting. */
	// WARNING: need to add 1 to any address iaddr or jaddr, otherwise
	// iaddr==0 will be interpreted as the null pointer.

	vertex_edges = malloc(npts*sizeof(setT));
	for (ivert=0; ivert<npts; ivert++){
		vertex_edges[ivert] = qh_setnew(buf_vert);
	}

	/* Loop over all simplices */
	for (isimplex=0; isimplex<nsimplex; isimplex++) {
		/* Don't need to consider invalid simplices */
		if (!transf_valid[isimplex]) continue;
		for (iv=0; iv<facet_ndim-1; iv++) {
			ivert = simplices[isimplex*facet_ndim + iv];
			iaddr = ivert + 1;
			for (jv=iv+1; jv<facet_ndim; jv++) {
				jvert = simplices[isimplex*facet_ndim + jv];
				jaddr = jvert + 1;
				//qh_setaddsorted(&vertex_edges[ivert], (void*) jaddr);
				//qh_setaddsorted(&vertex_edges[jvert], (void*) iaddr);
				qh_setunique(&vertex_edges[ivert], (void*) jaddr);
				qh_setunique(&vertex_edges[jvert], (void*) iaddr);
			}
		}
	}
#if 0
	printf("Number of neighbors per vertex:");
	for (ivert=0; ivert<npts; ivert++){
		if (!vertex_edges[ivert]) continue;
		printf("%03d - %03d\n", ivert, qh_setsize(vertex_edges[ivert]));
	}
	fflush(stdout);
	fflush(stderr);
#endif

	/* Now, loop over all simplices and figure out the number of first and
	 * second neighbors */
	first_neighbors = qh_setnew(facet_ndim);
	second_neighbors = malloc(nsimplex*sizeof(setT));
	for (isimplex=0; isimplex<nsimplex; isimplex++) {
		qh_setzero(first_neighbors, 0, facet_ndim);
		/* Don't need to consider invalid simplices */
		if (!transf_valid[isimplex]) continue;
		first_neighbors = qh_setnew(facet_ndim);
		/* Build list of first neighbors, i.e., vertices of simplex */
		for (iv=0; iv<facet_ndim; iv++) {
			ivert = simplices[isimplex*facet_ndim + iv];
			iaddr = ivert + 1;
			qh_setaddsorted(&first_neighbors, (void*) iaddr);
		}
		/* Build list of second-nearest neighbors */
		second_neighbors[isimplex] = qh_setnew(buf_simpl);
		for (iv=0; iv<facet_ndim; iv++) {
			ivert = simplices[isimplex*facet_ndim + iv];
			set = vertex_edges[ivert];
			sz = qh_setsize(set);
			// Loop over all connected vertices associated with ivert.
			for (jv=0; jv<sz; jv++){
				jaddr = (intptr_t) set->e[jv].p;
				/* Only include explicit second neighbors here */
				if (!qh_setin(first_neighbors, (void*) jaddr)) {
					qh_setaddsorted(&second_neighbors[isimplex], (void*) jaddr);
				}
			}
		}
	}
	qh_setfree(&first_neighbors);
#if 0
	printf("Number of second neighbors per simplex:\n");
	for (isimplex=0; isimplex<nsimplex; isimplex++) {
		/* Don't need to consider invalid simplices */
		if (!transf_valid[isimplex]) continue;
		printf("%03d - %03d\n", isimplex,
				qh_setsize(second_neighbors[isimplex]));
	}
	fflush(stdout);
	fflush(stderr);
#endif
	for (isimplex=0; isimplex<nsimplex; isimplex++) {
		/* Don't need to consider invalid simplices */
		if (!transf_valid[isimplex]) continue;
		sz = qh_setsize(second_neighbors[isimplex]);
		max_num_second_neighbors = sz>max_num_second_neighbors? sz: max_num_second_neighbors;
	}

	/* Free edges set and array of sets */
	for (ivert=0; ivert<npts; ivert++){
		qh_setfree(&vertex_edges[ivert]);
	}
	free(vertex_edges);

	return 0;
}

/* Remeber: everything here is in C order/starting indexing! */
int get_second_neighbors(int isimplex, int *indices, int *num_second)
{
	int iv;
	intptr_t iaddr;
	setT *set;

	set = second_neighbors[isimplex];
	*num_second = qh_setsize(set);
	for (iv=0; iv<(*num_second); iv++){
		iaddr = (intptr_t) set->e[iv].p;
		// Subtract 1 because we had to add 1 before to avoid
		// confusions with the NULL pointers;
		indices[iv] = iaddr - 1;
	}
	return 0;
}

int free_second_neighbors()
{
	int npts = qh num_points;
	int ivert, isimplex;

	if (max_num_second_neighbors==-1) return 0;

	/* Free all sets with neighbors and array of sets */
	for (isimplex=0; isimplex<nsimplex; isimplex++){
		if (!transf_valid[isimplex]) continue;
		qh_setfree(&second_neighbors[isimplex]);
	}
	free(second_neighbors);
	max_num_second_neighbors = -1;
	return 0;
}
