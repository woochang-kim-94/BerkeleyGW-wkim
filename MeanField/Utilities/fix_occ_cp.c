/*
* 
* fix_occ_cp     Originally By FHJ      Last Modified 10/17/2010 (FHJ)
*
* This file contains the routine copy_tail, which copies the rest of the
*  WFN produced by fix_occ. Using this C function can be way faster than
*  using the Fortran read_gvecs subroutine, since it doesn`t have all
*  the consistency checks, and it is also way more fun! 
*
* (Originally by FHJ)
*
*/

#include <stdlib.h>
#include <stdio.h>

// This is the buffer size for copying the file. You can make this as big
// as you want, as long as the memory fits in the HEAP. In particular, Lustre
// file systems should benefit considerably from larger buffer sizes.
#define BUF_SIZE 65536

// Copies the rest of the wave function with fixed header
// It is way faster (and more fun) to do this in C than in Fortran,
//  since we are avoiding all overhead of doing consistency checks
void copy_tail_(unsigned char *ierror){
	FILE *fp_in, *fp_out;
	char fname_in[256], fname_out[256], *ptr;
	int sz, cnt;
	char buf[BUF_SIZE];

	*ierror=0;

	fp_in = fopen("fix_occ.inp","r");

	if ( !fp_in ){
		printf("Error opening input file '%s'\n", "fix_occ.inp");
		*ierror=1;
		return;
	}

	if (!fgets(fname_in, 255, fp_in)){
		printf("Error reading input file '%s'\n", fname_in);
		*ierror=2;
		return;
	}
	if (!fgets(fname_out, 255, fp_in)){
		printf("Error reading output file '%s'\n", fname_out);
		*ierror=2;
		return;
	}

	fclose(fp_in);

	//remove \n from file names	
	ptr=fname_in;
	while(*ptr!='\n'||!*ptr) ptr++; *ptr=0;
	ptr=fname_out;
	while(*ptr!='\n'||!*ptr) ptr++; *ptr=0;

	fp_in = fopen(fname_in, "rb");
	fp_out = fopen(fname_out, "rb+");

	if ( (!fp_in) || (!fp_out) ){
		if (!fp_in){
			printf("Error opening input file '%s'\n", fname_in);
		}else{
			printf("Error opening output file '%s'\n", fname_out);
		}
		*ierror=3;
		return;
	}

	//go to the end of the files
	fseek(fp_out, 0L, SEEK_END);
	sz = ftell(fp_out);	
	printf(" Header size: %d bytes\n", sz);
	fseek(fp_in, sz, SEEK_SET);

	//copy
	while (!feof(fp_in)){
		if ((cnt=fread(buf, 1, BUF_SIZE, fp_in))!=BUF_SIZE){
			if (!feof(fp_in)){
				printf("Error reading input file (read %d bytes)\n", cnt);
				*ierror=ferror(fp_in);
				return;
			}
		}
		if ((int)fwrite(buf, 1, cnt, fp_out)!=cnt){
			printf("Error writing output WFN file\n");
			*ierror=5;
			return;
		}
	}	

}
void copy_tail(unsigned char *ierror){
	copy_tail_(ierror);
}
