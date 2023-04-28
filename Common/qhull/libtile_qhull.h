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

#ifndef libtile_qhull_h
#define libtile_qhull_h

extern int _inode;

/* Exported functions. We add flavors with and without leading underscore to
   make this wrapper compatible with Fortran and C. */
int qhull_init_delaunay_(double points[], const int *npts, const int *dim, const int *inode);
int qhull_init_delaunay(double points[], const int *npts, const int *dim, const int *inode);
int qhull_find_delaunay_simplex_(const double point[], int indices[], double coefs[]);
int qhull_find_delaunay_simplex(const double point[], int indices[], double coefs[]);
int qhull_find_delaunay_simplex_with_second_neighbors_(const double point[], int indices[], double coefs[], int *nnei);
int qhull_find_delaunay_simplex_with_second_neighbors(const double point[], int indices[], double coefs[], int *nnei);
int qhull_get_num_simplices_(int *num_simplices);
int qhull_get_num_simplices(int *num_simplices);
int qhull_get_simplices_(int indices[]);
int qhull_get_simplices(int indices[]);
int qhull_get_neighbors_(int neighbors[]);
int qhull_get_neighbors(int neighbors[]);
int qhull_free_delaunay_();
int qhull_free_delaunay();

int qhull_init_second_neighbors();
int qhull_init_second_neighbors_();
int qhull_free_second_neighbors();
int qhull_free_second_neighbors_();
int qhull_get_max_num_second_neighbors(int *max_nei);
int qhull_get_max_num_second_neighbors_(int *max_nei);

void qhull_info(const char *message, ...);
void qhull_warn(const char *message, ...);
void qhull_error(const char *message, ...);

#ifdef DEBUG
#define QHULL_STR "qhull s d Tcv Qt Qbb Qc Qz"
#else
#define QHULL_STR "qhull s d Qt Qbb Qc Qz"
#endif

#endif
