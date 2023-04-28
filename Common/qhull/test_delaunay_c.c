/******************************************************************************
*
*  test_delaunay_c        Originally by FHJ     Last Modified 12/04/2012 (FHJ)
*
*    A simple C program that tests the BerkeleyGW C bindings for Qhull.
*
*    This file is part of the BerkeleyGW package.
*
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "libtile_qhull.h"

void print_vect(void *fp, int id, double *vec, int dims)
{
	int jdim_;

	if (id>=0) fprintf(fp, " %4d ", id);
	for (jdim_=0; jdim_<dims; jdim_++)
		fprintf(fp, "%6.3f ", vec[jdim_]);
	fprintf(fp, "\n");
}

int main(int argc, char **argv)
{
	int dim_, npts, nrep;
	int ip, idim_, div, rem;
	int *indices, *base;
	const int inode = 0;
	double *points, *point, *pt, *coefs;
	double *coefs2;
	int *indices2;
	int nnei, nnei_max;
	char show_coords = 1;
	void *fp;
	int npts_orig;

	nrep = 1000;
	printf("Testing libtile_qhull - C bindings.\n");
	fp = fopen("points.dat", "r");
	fscanf(fp, "%d %d %d", &dim_, &npts, &npts_orig);
	points = (double*) malloc(sizeof(double)*dim_*npts);
	point = points;
	for (ip=0; ip<npts; ip++) {
		if (dim_==1) {
			fscanf(fp, "%lf", point);
		} else if (dim_==2) {
			fscanf(fp, "%lf %lf", point, point+1);
		} else if (dim_==3) {
			fscanf(fp, "%lf %lf %lf", point, point+1, point+2);
		}
		point += dim_;
	}
        fclose(fp);

	printf("Dimensions: %d\n", dim_);
	printf("Number of points: %d\n", npts);
	printf("Number of points w/o ghost points: %d\n", npts_orig);

	if (show_coords) {
		printf("\nInput coordinates:\n");
		for (ip=0; ip<npts; ip++) {
			point = points + ip*dim_;
			print_vect(stdout, ip+1, point, dim_);
		}
	}

	indices = (int*) malloc(sizeof(int) * (dim_+1));
	coefs = (double*) malloc(sizeof(double) * (dim_+1));
	pt = (double*) calloc(dim_, sizeof(double));
	pt[0] = 0.125;
	pt[1] = 0.75;
	printf("\nTest coordinate:\n");
	print_vect(stdout, 0, pt, dim_);
	printf("\n");

	printf("Calling init_delaunay.\n");
	fflush(stdout);
	fflush(stderr);
	qhull_init_delaunay(points, &npts, &dim_, &inode);
	fflush(stdout);
	fflush(stderr);

	printf("\nFinding max number of second neighbors:\n");
	fflush(stdout);
	fflush(stderr);
	qhull_init_second_neighbors();
	qhull_get_max_num_second_neighbors(&nnei_max);
	printf("%d\n\n", nnei_max);

	indices2 = (int*) malloc(sizeof(int) * (dim_ + 1 + nnei_max));
	coefs2 = (double*) malloc(sizeof(double) * (dim_ + 1 + nnei_max));
	printf("Calling find_delaunay_simplex_with_second_neighbors.\n");

	qhull_find_delaunay_simplex_with_second_neighbors(pt, indices2, coefs2, &nnei);
	printf("\nFound vertices:\n");
	for (ip=0; ip<nnei; ip++) {
		point = points + (indices2[ip]-1)*dim_;
		print_vect(stdout, indices2[ip], point, dim_);
	}
	printf("\nCoefficients: ");
	print_vect(stdout, -1, coefs2, nnei);

        printf("\nWriting output data to file 'points_out_c.dat'.\n");
        // Write header
        fp = fopen("points_out_c.dat", "w");
        fprintf(fp, "%d %d %d %d\n", dim_, npts, npts_orig, nnei);
        // Input coordinates
        for (ip=0; ip<npts; ip++) {
                point = points + ip*dim_;
                print_vect(fp, ip+1, point, dim_);
        }
        fprintf(fp, "\n");
        // First + second neighbors
	for (ip=0; ip<nnei; ip++) {
		point = points + (indices2[ip]-1)*dim_;
		print_vect(fp, indices2[ip], point, dim_);
	}
        fprintf(fp, "\n");
        // Original point
        print_vect(fp, 0, pt, dim_);
        fclose(fp);

	printf("\n");
	printf("Calling find_delaunay_simplex %8d times.\n", nrep);
	for (ip=0; ip<nrep; ip++)
		qhull_find_delaunay_simplex(pt, indices, coefs);
	printf("Calling free_delaunay.\n");
	qhull_free_delaunay();

	printf("\nFound vertices:\n");
	for (ip=0; ip<dim_+1; ip++) {
		point = points + (indices[ip]-1)*dim_;
		print_vect(stdout, indices[ip], point, dim_);
	}
	printf("\nCoefficients: ");
	print_vect(stdout, -1, coefs, dim_+1);

	free(indices);
	free(coefs);
	free(pt);
	free(points);
	printf("\nAll Done!\n");
	return 0;
}
