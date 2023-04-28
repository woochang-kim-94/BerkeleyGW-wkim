#include <stdio.h>

extern void open_primme_file(FILE **fp, int *ik, int *ip);
extern void open_primme_file_(FILE **fp, int *ik, int *ip);
extern void close_primme_file(FILE **fp, int *ik);
extern void close_primme_file_(FILE **fp, int *ik);


void open_primme_file(FILE **fp, int *ik, int *ip) {
	char fname[128];

	sprintf(fname, "primme_%d.log", *ip);
	*fp = fopen(fname, "a");
	fprintf(*fp, "\n");
	//fprintf(*fp, "------------------------------------------------------------------------------\n");
	fprintf(*fp, "PRIMME output for k-point %d\n", *ik);
	//fprintf(*fp, "------------------------------------------------------------------------------\n");
	fprintf(*fp, "\n");
	fflush(*fp);
}
void open_primme_file_(FILE **fp, int *ik, int *ip) {
	open_primme_file(fp, ik, ip);
}


void close_primme_file(FILE **fp, int *ik) {
	fprintf(*fp, "\n");
	//fprintf(*fp, "------------------------------------------------------------------------------\n");
	fprintf(*fp, "Done with k-point %d\n", *ik);
	//fprintf(*fp, "------------------------------------------------------------------------------\n");
	fprintf(*fp, "\n");
	fprintf(*fp, "==============================================================================\n");
	fflush(*fp);
	fclose(*fp);
}
void close_primme_file_(FILE **fp, int *ik) {
	close_primme_file(fp, ik);
}
