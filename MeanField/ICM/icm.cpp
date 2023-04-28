/*------------------------------------------------------------------------------

   icm.cpp
   computes the Coulomb integral within an image-charge model
   written by Georgy Samsonidze (January 2010)

   Based on surface.cpp for generating an isosurface of the scalar field.
   Objects used by ICM refactored out into Common/wfn_utils.cpp by DAS

   The following common objects are used from Common/wfn_utils.cpp
   macros: MAX
   structures: CARTESIAN, ELEMENT
   constants: MAXCHAR, EPS9, INF9, BOHR, vertex, periodic_table
   functions: lowercase, erasechar, parse, cub_read, xsf_read, cell_set,
   normal_make, box_check, scalar_clone, inversion, isovalue_scale

   The Coulomb integral is from J. B. Neaton, M. S. Hybertsen, and S. G. Louie,
   Renormalization of Molecular Electronic Levels at Metal-Molecule Interfaces,
   Phys. Rev. Lett. 97, 216405 (2006).

   If the wavefunction overlaps with the image plane, 1 / | r - r' | is averaged
   over r' in the grid cell by Monte Carlo algorithm. Make sure two of the
   lattice vectors lie in the image plane and the third one is normal to it,
   otherwise the number of (r,r') pairs to be averaged increases significantly
   and the code will refuse to run.

   Exactly one of the two blocks in brackets below should appear. Supply the
   mirror plane lines to do the image-charge model; supply a second filename
   to compute the Coulomb integral of two different wavefunctions or densities.

--------------------------------------------------------------------------------

   Input is read from file icm.inp   |   example (HOMO of benzene)

   inputfilename file                |   inputfilename C6H6.b_15.cube
   inputfileformat cube|xsf          |   inputfileformat cube
   threshold (0.0,1.0]               |   threshold 0.99
   threshold_power 0|1|2             |   threshold_power 1
   coulomb_power 1|2                 |   coulomb_power 1
[
   mirrorplaneorigin                 |   mirrorplaneorigin
   <mpox> <mpoy> <mpoz>              |   0.0 0.0 -2.0
   mirrorplanenormal                 |   mirrorplanenormal
   <mpnx> <mpny> <mpnz>              |   0.0 0.0 1.0
   mirrorplaneunit bohr|angstrom     |   mirrorplaneunit angstrom
][
   inputfilename2 file               |   inputfilename2 C6H6.b_15.cube
   inputfileformat2 cube|xsf         |   inputfileformat2 cube
]
   uc T|F                            |   uc F
   uco                               |   uco
   <ucox> <ucoy> <ucoz>              |   0.0 0.0 0.0
   ucv                               |   ucv
   <ucv1x> <ucv1y> <ucv1z>           |   1.0 0.0 0.0
   <ucv2x> <ucv2y> <ucv2z>           |   0.0 1.0 0.0
   <ucv3x> <ucv3y> <ucv3z>           |   0.0 0.0 1.0
   ucu bohr|angstrom|latvec          |   ucu latvec
   sc T|F                            |   sc T
   sco                               |   sco
   <scox> <scoy> <scoz>              |   -0.5 -0.5 -0.5
   scv                               |   scv
   <scv1x> <scv1y> <scv1z>           |   1.0 0.0 0.0
   <scv2x> <scv2y> <scv2z>           |   0.0 1.0 0.0
   <scv3x> <scv3y> <scv3z>           |   0.0 0.0 1.0
   scu bohr|angstrom|latvec          |   scu latvec
   renormalize                       |   F

--------------------------------------------------------------------------------

In the above example, the HOMO wavefunction of benzene is read from a
Gaussian Cube file. The wavefunction is placed in the center of the
supercell (see the meaning of uc, uco, ucv, ucu, sc, sco, scv, scu
parameters in Visual/surface.cpp). The parts of the wavefunction
outside an isosurface that contains 99% of the charge density are
dropped (parameters threshold and threshold_power have the same
meaning as isovalue and power in Visual/surface.cpp). Parameter
coulomb_power tells the code whether the wavefunction in the Coulomb
integral needs to be squared. Set both powers to 1 if the wavefunction
file contains the squared amplitude as produced by ESPRESSO, and to 2
for the linear amplitude as in PARATEC or SIESTA. The mirror plane
is defined by parameters mirrorplaneorigin and mirrorplanenormal;
in the above example it is parallel to the xy plane crossing the
z axis at -2 Angstrom.

------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
//#include <time.h>
#include "wfn_utils.h"

#ifdef PARA
#include "mpi.h"
#endif

const int NCELL = 2;
const int NRAND = 2500000;
const double HARTREE = 27.21138505;

// wtf, these are global variables!!
double ravg[2*NCELL+1][2*NCELL+1][2*NCELL+1];

struct cell_info
{
   bool uc, sc;
   int ucf, scf;
   CARTESIAN sfv[3], uco, ucv[3], sco, scv[3];
   bool renorm;
};

double volume(const CARTESIAN *vec)
{
   return fabs(vec[0].xx * (vec[1].yy * vec[2].zz - vec[2].yy * vec[1].zz) -
	       vec[0].yy * (vec[1].xx * vec[2].zz - vec[2].xx * vec[1].zz) +
	       vec[0].zz * (vec[1].xx * vec[2].yy - vec[2].xx * vec[1].yy));
}

int par_read(const char *pfn, char *ifn, char *iff, char *ifn2, char *iff2, double *threshold, int *threshold_power, int *coulomb_power, CARTESIAN *mpo, CARTESIAN *mpn)
{
   int ierr, icount = 0, icheck = 0, icount2 = 0, icount_mirror = 0;
   char s1[MAXCHAR], s2[MAXCHAR], s3[MAXCHAR];
   char* trash;
   FILE *hh;
   bool is_bohr;

   ifn2[0] = '\0';
   iff2[0] = '\0';

   hh = fopen(pfn, "r");
   if (hh == NULL)
      return -1;

   while (!feof(hh))
   {
      strncpy(s1, "\0", 1);
      trash = fgets(s1, MAXCHAR, hh);
      parse(s1, s2, s3);
      lowercase(s2);
      if (strcmp(s2, "inputfilename") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         strncpy(ifn, s3, MAXCHAR);
         ifn[MAXCHAR - 1] = '\0';
      }
      else if (strcmp(s2, "inputfileformat") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         strncpy(iff, s3, MAXCHAR);
         iff[MAXCHAR - 1] = '\0';
      }
      if (strcmp(s2, "inputfilename2") == 0)
      {
         icount2++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         strncpy(ifn2, s3, MAXCHAR);
         ifn2[MAXCHAR - 1] = '\0';
      }
      else if (strcmp(s2, "inputfileformat2") == 0)
      {
         icount2++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         strncpy(iff2, s3, MAXCHAR);
         iff2[MAXCHAR - 1] = '\0';
      }
      else if (strcmp(s2, "threshold") == 0)
      {
         icount++;
         *threshold = atof(s3);
         if (*threshold < EPS9 || *threshold > 1.0 + EPS9)
	 {
	    fprintf(stderr, "Value for 'threshold' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "threshold_power") == 0)
      {
         icount++;
         *threshold_power = atoi(s3);
         if (*threshold_power < 0 || *threshold_power > 2)
	 {
	    fprintf(stderr, "Value for 'threshold_power' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "coulomb_power") == 0)
      {
         icount++;
         *coulomb_power = atoi(s3);
         if (*coulomb_power < 1 || *coulomb_power > 2)
	 {
	    fprintf(stderr, "Value for 'coulomb_power' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "mirrorplaneorigin") == 0)
      {
         icount_mirror++;
         ierr = fscanf(hh, "%le%le%le\n", &mpo->xx, &mpo->yy, &mpo->zz);
      }
      else if (strcmp(s2, "mirrorplanenormal") == 0)
      {
         icount_mirror++;
         ierr = fscanf(hh, "%le%le%le\n", &mpn->xx, &mpn->yy, &mpn->zz);
      }
      else if (strcmp(s2, "mirrorplaneunit") == 0)
      {
         icount_mirror++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "bohr") == 0)
            is_bohr = true;
         else if (strcmp(s3, "angstrom") == 0)
            is_bohr = false;
         else
	 {
	    fprintf(stderr, "Value for 'mirrorplaneunit' is not legal.\n");
            icheck--;
	 }
      }
   }

   if(icount2 > 0 && icount_mirror > 0) {
      fprintf(stderr, "Set parameters about only one of second filename or mirror plane.\n");
      return -1;
   }

   if(icount2 <= 0 && icount_mirror <= 0) {
      fprintf(stderr, "Must set parameters about one of second filename or mirror plane.\n");
      return -1;
   }

   if (icount_mirror > 0 && !is_bohr)
   {
      mpo->xx /= BOHR;
      mpo->yy /= BOHR;
      mpo->zz /= BOHR;
   }

   ierr = fclose(hh);
   if (ierr != 0)
      return -1;

   if (icount != 5)
   {
      fprintf(stderr, "Not all non-cell parameters were set.\n");
      return -1;
   }

   if (icount2 > 0 && icount2 != 2)
   {
      fprintf(stderr, "Not all parameters about second filename were set.\n");
      return -1;
   }

   if (icount_mirror > 0 && icount_mirror != 3)
   {
      fprintf(stderr, "Not all parameters about mirror plane were set.\n");
      return -1;
   }

   if (icheck < 0)
      return -1;

   return 0;
}

int par_read_cell(const char *pfn, cell_info & ci, char suffix)
{
   int ii, ierr, icount = 0, icheck = 0;
   char s1[MAXCHAR], s2[MAXCHAR], s3[MAXCHAR];
   char* trash;
   FILE *hh;

   hh = fopen(pfn, "r");
   if (hh == NULL)
      return -1;

   while (!feof(hh))
   {
      strncpy(s1, "\0", 1);
      trash = fgets(s1, MAXCHAR, hh);
      parse(s1, s2, s3);
      lowercase(s2);

      if(suffix != '\0')
      {
	 if(s2[strlen(s2)-1] == suffix)
	 {
	    // remove suffix
	    s2[strlen(s2)-1] = '\0';
	 }
	 else
	 {
	    continue;
	 }
      }

      if (strcmp(s2, "uc") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "f") == 0 || strcmp(s3, "false") == 0 || strcmp(s3, "n") == 0 || strcmp(s3, "no") == 0)
            ci.uc = false;
         else if (strcmp(s3, "t") == 0 || strcmp(s3, "true") == 0 || strcmp(s3, "y") == 0 || strcmp(s3, "yes") == 0)
            ci.uc = true;
         else
	 {
	    fprintf(stderr, "Value for 'uc' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "uco") == 0)
      {
         icount++;
         ierr = fscanf(hh, "%le%le%le\n", &ci.uco.xx, &ci.uco.yy, &ci.uco.zz);
      }
      else if (strcmp(s2, "ucv") == 0)
      {
         icount++;
         for (ii = 0; ii < 3; ii++)
            ierr = fscanf(hh, "%le%le%le\n", &ci.ucv[ii].xx, &ci.ucv[ii].yy, &ci.ucv[ii].zz);
      }
      else if (strcmp(s2, "ucu") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "bohr") == 0)
            ci.ucf = 0;
         else if (strcmp(s3, "angstrom") == 0)
            ci.ucf = 1;
         else if (strcmp(s3, "latvec") == 0)
            ci.ucf = 2;
         else
	 {
	    fprintf(stderr, "Value for 'ucu' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "sc") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "f") == 0 || strcmp(s3, "false") == 0 || strcmp(s3, "n") == 0 || strcmp(s3, "no") == 0)
            ci.sc = false;
         else if (strcmp(s3, "t") == 0 || strcmp(s3, "true") == 0 || strcmp(s3, "y") == 0 || strcmp(s3, "yes") == 0)
            ci.sc = true;
         else
	 {
	    fprintf(stderr, "Value for 'sc' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "sco") == 0)
      {
         icount++;
         ierr = fscanf(hh, "%le%le%le\n", &ci.sco.xx, &ci.sco.yy, &ci.sco.zz);
      }
      else if (strcmp(s2, "scv") == 0)
      {
         icount++;
         for (ii = 0; ii < 3; ii++)
            ierr = fscanf(hh, "%le%le%le\n", &ci.scv[ii].xx, &ci.scv[ii].yy, &ci.scv[ii].zz);
      }
      else if (strcmp(s2, "scu") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "bohr") == 0)
            ci.scf = 0;
         else if (strcmp(s3, "angstrom") == 0)
            ci.scf = 1;
         else if (strcmp(s3, "latvec") == 0)
            ci.scf = 2;
         else
	 {
	    fprintf(stderr, "Value for 'scu' is not legal.\n");
            icheck--;
	 }
      }
      else if (strcmp(s2, "renormalize") == 0)
      {
         icount++;
         erasechar(s3, ' ');
         erasechar(s3, '\t');
         lowercase(s3);
         if (strcmp(s3, "f") == 0 || strcmp(s3, "false") == 0 || strcmp(s3, "n") == 0 || strcmp(s3, "no") == 0)
            ci.renorm = false;
         else if (strcmp(s3, "t") == 0 || strcmp(s3, "true") == 0 || strcmp(s3, "y") == 0 || strcmp(s3, "yes") == 0)
            ci.renorm = true;
         else
	 {
	    fprintf(stderr, "Value for 'renorm' is not legal.\n");
            icheck--;
	 }
      }
   }

   ierr = fclose(hh);
   if (ierr != 0)
      return -1;

   if (icount != 9)
   {
      if(suffix == '\0')
      {
	 fprintf(stderr, "Not all cell parameters were set.\n");
      }
      else
      {
	 fprintf(stderr, "Not all cell parameters were set, with suffix %c.\n", suffix);
      }
   }

   if (icount != 9 || icheck < 0)
      return -1;

   return 0;
}

int scalar_trunc(double threshold, scalar_field & sf)
{
   int ii, jj, kk, imin, imax, jmin, jmax, kmin, kmax, ni0, nj0, nk0, ni1, nj1, nk1;
   double ***scalar0 = NULL;

   ni0 = sf.ni;
   nj0 = sf.nj;
   nk0 = sf.nk;

   imin = ni0;
   imax = -1;
   jmin = nj0;
   jmax = -1;
   kmin = nk0;
   kmax = -1;

   for (ii = 0; ii < ni0; ii++)
   for (jj = 0; jj < nj0; jj++)
   for (kk = 0; kk < nk0; kk++)
      if (fabs(sf.scalar[ii][jj][kk]) > threshold)
      {
         if (ii < imin)
            imin = ii;
         if (ii > imax)
            imax = ii;
         if (jj < jmin)
            jmin = jj;
         if (jj > jmax)
            jmax = jj;
         if (kk < kmin)
            kmin = kk;
         if (kk > kmax)
            kmax = kk;
      }

   ni1 = imax - imin + 1;
   nj1 = jmax - jmin + 1;
   nk1 = kmax - kmin + 1;

   if (ni1 < 1 || nj1 < 1 || nk1 < 1)
      return -1;

   if (ni1 != ni0 || nj1 != nj0 || nk1 != nk0)
   {
      sf.ni = ni1;
      sf.nj = nj1;
      sf.nk = nk1;

      scalar0 = sf.scalar;

      sf.scalar = allocate_scalar(ni1, nj1, nk1);

      for (ii = 0; ii < ni1; ii++)
      for (jj = 0; jj < nj1; jj++)
      for (kk = 0; kk < nk1; kk++)
         sf.scalar[ii][jj][kk] = scalar0[ii + imin][jj + jmin][kk + kmin];

      deallocate_scalar(scalar0, ni0, nj0);

      sf.origin.xx += sf.stepv[0].xx * double(imin) + sf.stepv[1].xx * double(jmin) + sf.stepv[2].xx * double(kmin);
      sf.origin.yy += sf.stepv[0].yy * double(imin) + sf.stepv[1].yy * double(jmin) + sf.stepv[2].yy * double(kmin);
      sf.origin.zz += sf.stepv[0].zz * double(imin) + sf.stepv[1].zz * double(jmin) + sf.stepv[2].zz * double(kmin);
   }

   if (ni1 == ni0 || nj1 == nj0 || nk1 == nk0)
   {
      fprintf(stderr, "WARNING: The selected isovalue overlaps with the edge of the cell.\n");
      fprintf(stderr, "Beware, the Coulomb integral does not use periodic boundary conditions.\n");
   }

   return 0;
}

int scalar_norm(int coulomb_power, double *wfnorm, scalar_field & sf, bool renorm)
{
   int ii, jj, kk, ierr = 0;
   double vv, ww = 0.0;

   for (ii = 0; ii < sf.ni; ii++)
   for (jj = 0; jj < sf.nj; jj++)
   for (kk = 0; kk < sf.nk; kk++)
      if (coulomb_power == 1)
         sf.scalar[ii][jj][kk] = fabs(sf.scalar[ii][jj][kk]);
      else
         sf.scalar[ii][jj][kk] = sf.scalar[ii][jj][kk] * sf.scalar[ii][jj][kk];

   vv = volume(sf.stepv) / (BOHR * BOHR * BOHR);

   if (vv < EPS9)
      ierr = -1;

   for (ii = 0; ii < sf.ni; ii++)
   for (jj = 0; jj < sf.nj; jj++)
   for (kk = 0; kk < sf.nk; kk++)
         ww += sf.scalar[ii][jj][kk];

   if (ww < EPS9)
      ierr = -1;

   ww *= vv;
   *wfnorm = ww;

   if(renorm)
   {
      for (ii = 0; ii < sf.ni; ii++)
      for (jj = 0; jj < sf.nj; jj++)
      for (kk = 0; kk < sf.nk; kk++)
         sf.scalar[ii][jj][kk] /= *wfnorm;
   }

   return ierr;
}

CARTESIAN mirror_transform(const CARTESIAN *in, const CARTESIAN *mpn)
{
  double pp;
  CARTESIAN dd;

  pp = in->xx * mpn->xx + in->yy * mpn->yy + in->zz * mpn->zz;
  dd.xx = in->xx - 2 * pp * mpn->xx;
  dd.yy = in->yy - 2 * pp * mpn->yy;
  dd.zz = in->zz - 2 * pp * mpn->zz;

  return dd;
}

// convert origin and step vectors from Angstrom to Bohr
void convert_to_bohr(scalar_field & sf1)
{
   int ii;

   sf1.origin.xx /= BOHR;
   sf1.origin.yy /= BOHR;
   sf1.origin.zz /= BOHR;

   for (ii = 0; ii < 3; ii++)
   {
      sf1.stepv[ii].xx /= BOHR;
      sf1.stepv[ii].yy /= BOHR;
      sf1.stepv[ii].zz /= BOHR;
   }
}

int mirror_plane(CARTESIAN *mpo, CARTESIAN *mpn, scalar_field & sf1, scalar_field & sf2)
{
   int ii;
   double pp;
   CARTESIAN aa, dd;

   pp = mpn->xx * mpn->xx + mpn->yy * mpn->yy + mpn->zz * mpn->zz;
   if (pp < EPS9)
      return -1;
   pp = sqrt(pp);
   if (fabs(1.0 - pp) > EPS9)
   {
      mpn->xx /= pp;
      mpn->yy /= pp;
      mpn->zz /= pp;
   }

   for (ii = 0; ii < 4; ii++)
   {
      if (ii == 0)
      {
         aa.xx = sf1.origin.xx - mpo->xx;
         aa.yy = sf1.origin.yy - mpo->yy;
         aa.zz = sf1.origin.zz - mpo->zz;
      }
      else
      {
         aa.xx = sf1.stepv[ii - 1].xx;
         aa.yy = sf1.stepv[ii - 1].yy;
         aa.zz = sf1.stepv[ii - 1].zz;
      }

      dd = mirror_transform(&aa, mpn);

      if (ii == 0)
      {
         sf2.origin.xx = mpo->xx + dd.xx;
         sf2.origin.yy = mpo->yy + dd.yy;
         sf2.origin.zz = mpo->zz + dd.zz;
      }
      else
      {
         sf2.stepv[ii - 1].xx = dd.xx;
         sf2.stepv[ii - 1].yy = dd.yy;
         sf2.stepv[ii - 1].zz = dd.zz;
      }
   }

   return 0;
}

// shortest step vector of v1 or v2, times NCELL. possibly it doesn't need to be v2 also?
void set_cutoff(const CARTESIAN *v1, const CARTESIAN *v2, double *rcutoff, double *r2cutoff)
{
   int ii;
   double l2, l2min = INF9;

   for (ii = 0; ii < 3; ii++)
   {
      l2 = v1[ii].xx * v1[ii].xx + v1[ii].yy * v1[ii].yy + v1[ii].zz * v1[ii].zz;
      if (l2min > l2)
         l2min = l2;

      l2 = v2[ii].xx * v2[ii].xx + v2[ii].yy * v2[ii].yy + v2[ii].zz * v2[ii].zz;
      if (l2min > l2)
         l2min = l2;
   }
   *r2cutoff = l2min * double(NCELL * NCELL);
   *rcutoff = sqrt(*r2cutoff);
}

// Do the two scalar fields overlap?
bool check_overlap(int rank, const scalar_field & sf1, const scalar_field & sf2, double rcutoff)
{
   bool flag_overlap = false;
   int ii, nplus = 0, nminus = 0, iminus, ngrid1, ngrid2, dplus = 0, dminus = 0;
   double vl1[3], vl2[3], vdotv, vcosv, vdotr, pp[4], dd[4];
   CARTESIAN rr[4];

   for (ii = 0; ii < 3; ii++)
   {
      vl1[ii] = sqrt(sf1.stepv[ii].xx * sf1.stepv[ii].xx + sf1.stepv[ii].yy * sf1.stepv[ii].yy + sf1.stepv[ii].zz * sf1.stepv[ii].zz);
      vl2[ii] = sqrt(sf2.stepv[ii].xx * sf2.stepv[ii].xx + sf2.stepv[ii].yy * sf2.stepv[ii].yy + sf2.stepv[ii].zz * sf2.stepv[ii].zz);
      vdotv = sf1.stepv[ii].xx * sf2.stepv[ii].xx + sf1.stepv[ii].yy * sf2.stepv[ii].yy + sf1.stepv[ii].zz * sf2.stepv[ii].zz;
      vcosv = vdotv / (vl1[ii] * vl2[ii]);
      if (fabs(vcosv - 1.0) < EPS9)
         nplus++;
      if (fabs(vcosv + 1.0) < EPS9)
      {
         nminus++;
         iminus = ii;
      }
   }

   if (nplus == 2 && nminus == 1)
   {
      switch (iminus)
      {
         case 0:
            ngrid1 = sf1.ni;
	    ngrid2 = sf2.ni;
         case 1:
            ngrid1 = sf1.nj;
	    ngrid2 = sf2.nj;
         case 2:
            ngrid1 = sf1.nk;
	    ngrid2 = sf2.nk;
      }

      rr[0].xx = sf1.origin.xx;
      rr[0].yy = sf1.origin.yy;
      rr[0].zz = sf1.origin.zz;
      rr[1].xx = sf1.origin.xx + sf1.stepv[iminus].xx * double(ngrid1 - 1);
      rr[1].yy = sf1.origin.yy + sf1.stepv[iminus].yy * double(ngrid1 - 1);
      rr[1].zz = sf1.origin.zz + sf1.stepv[iminus].zz * double(ngrid1 - 1);
      rr[2].xx = sf2.origin.xx;
      rr[2].yy = sf2.origin.yy;
      rr[2].zz = sf2.origin.zz;
      rr[3].xx = sf2.origin.xx + sf2.stepv[iminus].xx * double(ngrid2 - 1);
      rr[3].yy = sf2.origin.yy + sf2.stepv[iminus].yy * double(ngrid2 - 1);
      rr[3].zz = sf2.origin.zz + sf2.stepv[iminus].zz * double(ngrid2 - 1);

      for (ii = 0; ii < 4; ii++)
      {
         vdotr = sf1.stepv[iminus].xx * rr[ii].xx + sf1.stepv[iminus].yy * rr[ii].yy + sf1.stepv[iminus].zz * rr[ii].zz;
         pp[ii] = vdotr / vl1[iminus];
      }

      dd[0] = pp[0] - pp[2];
      dd[1] = pp[0] - pp[3];
      dd[2] = pp[1] - pp[2];
      dd[3] = pp[1] - pp[3];

      for (ii = 0; ii < 4; ii++)
         if (dd[ii] > 0.0)
            dplus++;
         else
            dminus++;

      if (dplus == 4 || dminus == 4)
      {
         for (ii = 0; ii < 4; ii++)
            if (fabs(dd[ii]) < rcutoff)
               flag_overlap = true;
      }
      else
        flag_overlap = true;
   }
   else
   {
      if (rank == 0)
      {
	 fprintf(stderr, "    The image plane is not parallel/perpendicular to the lattice vectors.\n");
	 fprintf(stderr, "    Skipping wavefunction overlap check and Monte Carlo averaging.\n");
	 fprintf(stderr, "    Your job may fail if the wavefunction overlaps with its image.\n\n");
      }
   }

   return flag_overlap;
}

CARTESIAN grid_offset(const scalar_field & sf1, const scalar_field & sf2)
{
   int ii, nn;
   double lv2, vdotd;
   CARTESIAN dd;

   dd.xx = sf1.origin.xx - sf2.origin.xx;
   dd.yy = sf1.origin.yy - sf2.origin.yy;
   dd.zz = sf1.origin.zz - sf2.origin.zz;

   for (ii = 0; ii < 3; ii++)
   {
      lv2 = sf1.stepv[ii].xx * sf1.stepv[ii].xx + sf1.stepv[ii].yy * sf1.stepv[ii].yy + sf1.stepv[ii].zz * sf1.stepv[ii].zz;
      vdotd = sf1.stepv[ii].xx * dd.xx + sf1.stepv[ii].yy * dd.yy + sf1.stepv[ii].zz * dd.zz;
      nn = double_to_int(vdotd / lv2);

      dd.xx = dd.xx - sf1.stepv[ii].xx * double(nn);
      dd.yy = dd.yy - sf1.stepv[ii].yy * double(nn);
      dd.zz = dd.zz - sf1.stepv[ii].zz * double(nn);
   }

   return dd;
}

void rand_init()
{
   int seed;
   //time_t seconds;

   //   seconds = time(NULL);
   //   seed = (int)seconds;
   seed = 5000;
   srand(seed);
}

void mc_average(int rank, int size, const CARTESIAN *offset, const CARTESIAN *v1, const CARTESIAN *v2)
{
   int ii, jj, kk, ll, mm, lmin, lmax, navg;
   int nlocal[2*NCELL+1][2*NCELL+1][2*NCELL+1];
#ifdef PARA
   int ndummy[2*NCELL+1][2*NCELL+1][2*NCELL+1];
#endif
   double rinv, r2, rr[3];
   double rlocal[2*NCELL+1][2*NCELL+1][2*NCELL+1];
#ifdef PARA
   double rdummy[2*NCELL+1][2*NCELL+1][2*NCELL+1];
#endif
   CARTESIAN dd, crand;

#ifdef VERBOSE
   if (rank == 0)
      printf("    averaging 1 / | r - r' | over r' in the grid cell\n\n");
#endif

   lmin = double_to_int(double(NRAND) * double(rank) / double(size));
   lmax = double_to_int(double(NRAND) * double(rank + 1) / double(size));

   for (ii = -NCELL; ii <= NCELL; ii++)
   for (jj = -NCELL; jj <= NCELL; jj++)
   for (kk = -NCELL; kk <= NCELL; kk++)
   {
      nlocal[ii+NCELL][jj+NCELL][kk+NCELL] = lmax - lmin;
      rlocal[ii+NCELL][jj+NCELL][kk+NCELL] = 0.0;
   }

   for (ll = lmin; ll < lmax; ll++)
   {
      for (mm = 0; mm < 3; mm++)
         rr[mm] = double(rand()) / double(RAND_MAX) - 0.5;

      crand.xx = v2[0].xx * rr[0] + v2[1].xx * rr[1] + v2[2].xx * rr[2];
      crand.yy = v2[0].yy * rr[0] + v2[1].yy * rr[1] + v2[2].yy * rr[2];
      crand.zz = v2[0].zz * rr[0] + v2[1].zz * rr[1] + v2[2].zz * rr[2];

      for (ii = -NCELL; ii <= NCELL; ii++)
      for (jj = -NCELL; jj <= NCELL; jj++)
      for (kk = -NCELL; kk <= NCELL; kk++)
      {
	 dd.xx = offset->xx + v1[0].xx * double(ii) + v1[1].xx * double(jj) + v1[2].xx * double(kk);
	 dd.yy = offset->yy + v1[0].yy * double(ii) + v1[1].yy * double(jj) + v1[2].yy * double(kk);
	 dd.zz = offset->zz + v1[0].zz * double(ii) + v1[1].zz * double(jj) + v1[2].zz * double(kk);

         r2 = (dd.xx - crand.xx) * (dd.xx - crand.xx) + (dd.yy - crand.yy) * (dd.yy - crand.yy) + (dd.zz - crand.zz) * (dd.zz - crand.zz);
         if (r2 > EPS9)
	    rlocal[ii+NCELL][jj+NCELL][kk+NCELL] += 1.0 / sqrt(r2);
         else
	    nlocal[ii+NCELL][jj+NCELL][kk+NCELL]--;
      }
   }

#ifdef PARA
   MPI_Allreduce(nlocal, ndummy, pow(2*NCELL+1, 3), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

   for (ii = -NCELL; ii <= NCELL; ii++)
   for (jj = -NCELL; jj <= NCELL; jj++)
   for (kk = -NCELL; kk <= NCELL; kk++)
      nlocal[ii+NCELL][jj+NCELL][kk+NCELL] = ndummy[ii+NCELL][jj+NCELL][kk+NCELL];

   MPI_Allreduce(rlocal, rdummy, pow(2*NCELL+1, 3), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);

   for (ii = -NCELL; ii <= NCELL; ii++)
   for (jj = -NCELL; jj <= NCELL; jj++)
   for (kk = -NCELL; kk <= NCELL; kk++)
      rlocal[ii+NCELL][jj+NCELL][kk+NCELL] = rdummy[ii+NCELL][jj+NCELL][kk+NCELL];
#endif

   for (ii = -NCELL; ii <= NCELL; ii++)
   for (jj = -NCELL; jj <= NCELL; jj++)
   for (kk = -NCELL; kk <= NCELL; kk++)
      ravg[ii+NCELL][jj+NCELL][kk+NCELL] = rlocal[ii+NCELL][jj+NCELL][kk+NCELL] / double(nlocal[ii+NCELL][jj+NCELL][kk+NCELL]);
}

double tab_average(const CARTESIAN *v1, const CARTESIAN *p1, const CARTESIAN *p2, const CARTESIAN *offset)
{
   bool f_offset, f_range;
   int ii, nn[3];
   double lv2, vdotd, rinv;
   CARTESIAN dd;

   dd.xx = p1->xx - p2->xx - offset->xx;
   dd.yy = p1->yy - p2->yy - offset->yy;
   dd.zz = p1->zz - p2->zz - offset->zz;

   for (ii = 0; ii < 3; ii++)
   {
      lv2 = v1[ii].xx * v1[ii].xx + v1[ii].yy * v1[ii].yy + v1[ii].zz * v1[ii].zz;
      vdotd = v1[ii].xx * dd.xx + v1[ii].yy * dd.yy + v1[ii].zz * dd.zz;
      nn[ii] = double_to_int(vdotd / lv2);

      dd.xx = dd.xx - v1[ii].xx * double(nn[ii]);
      dd.yy = dd.yy - v1[ii].yy * double(nn[ii]);
      dd.zz = dd.zz - v1[ii].zz * double(nn[ii]);
   }

   f_offset = fabs(dd.xx) < EPS9 && fabs(dd.yy) < EPS9 && fabs(dd.zz) < EPS9;
   f_range = abs(nn[0]) <= NCELL && abs(nn[1]) <= NCELL && abs(nn[2]) <= NCELL;
   if (f_offset && f_range)
      rinv = ravg[nn[0]+NCELL][nn[1]+NCELL][nn[2]+NCELL];
   else
   {
     if(~f_offset) fprintf(stderr, "\nError: f_offset failed in tab_average.\n");
     if(~f_range)  fprintf(stderr, "\nError: f_range failed in tab_average.\n");
     rinv = -1.0;
   }

   return rinv;
}

int coulomb_integral(int rank, int size, const CARTESIAN *offset, double r2cutoff, double *energy,
		     const scalar_field & sf1, const scalar_field & sf2)
{
  int i1, j1, k1, i2, j2, k2;
#ifdef VERBOSE
   int pnew, pold = -1;
#endif
   double r2, rinv, vol1, vol2, ww = 0.0;
#ifdef PARA
   double rdummy;
#endif
   CARTESIAN p1, p2;

   vol1 = volume(sf1.stepv);
   if (vol1 < EPS9)
   {
      fprintf(stderr, "\nError: cell volume 1 = 0 in coulomb_integral\n\n");
      return -1;
   }

   vol2 = volume(sf2.stepv);
   if (vol2 < EPS9)
   {
      fprintf(stderr, "\nError: cell volume 2 = 0 in coulomb_integral\n\n");
      return -1;
   }

   for (i1 = 0; i1 < sf1.ni; i1++)
   for (j1 = 0; j1 < sf1.nj; j1++)
   for (k1 = 0; k1 < sf1.nk; k1++)
   {

#ifdef PARA
      if ((i1 * sf1.nj * sf1.nk + j1 * sf1.nk + k1) % size != rank)
         continue;
#endif

#ifdef VERBOSE
      if (rank == 0)
      {
         pnew = double_to_int(100.0 * double(i1 * sf1.nj * sf1.nk + j1 * sf1.nk + k1) / double(sf1.ni * sf1.nj * sf1.nk));
         if (pnew != pold)
         {
            printf("    completed %3i%%\n", pnew);
            pold = pnew;
         }
      }
#endif

      p1.xx = sf1.origin.xx + sf1.stepv[0].xx * double(i1) + sf1.stepv[1].xx * double(j1) + sf1.stepv[2].xx * double(k1);
      p1.yy = sf1.origin.yy + sf1.stepv[0].yy * double(i1) + sf1.stepv[1].yy * double(j1) + sf1.stepv[2].yy * double(k1);
      p1.zz = sf1.origin.zz + sf1.stepv[0].zz * double(i1) + sf1.stepv[1].zz * double(j1) + sf1.stepv[2].zz * double(k1);

      for (i2 = 0; i2 < sf2.ni; i2++)
      for (j2 = 0; j2 < sf2.nj; j2++)
      for (k2 = 0; k2 < sf2.nk; k2++)
      {
         p2.xx = sf2.origin.xx + sf2.stepv[0].xx * double(i2) + sf2.stepv[1].xx * double(j2) + sf2.stepv[2].xx * double(k2);
         p2.yy = sf2.origin.yy + sf2.stepv[0].yy * double(i2) + sf2.stepv[1].yy * double(j2) + sf2.stepv[2].yy * double(k2);
         p2.zz = sf2.origin.zz + sf2.stepv[0].zz * double(i2) + sf2.stepv[1].zz * double(j2) + sf2.stepv[2].zz * double(k2);

         r2 = (p1.xx - p2.xx) * (p1.xx - p2.xx) + (p1.yy - p2.yy) * (p1.yy - p2.yy) + (p1.zz - p2.zz) * (p1.zz - p2.zz);

         if (r2 > r2cutoff)
            rinv = 1.0 / sqrt(r2);
         else
            rinv = tab_average(sf1.stepv, &p1, &p2, &(*offset));

         if (rinv < -EPS9)
	   {
	     fprintf(stderr, "\nError: rinv < 0 in coulomb_integral\n\n");
	     return -1;
	   }

         ww += sf1.scalar[i1][j1][k1] * sf2.scalar[i2][j2][k2] * rinv;
      }
   }

#ifdef PARA
   MPI_Barrier(MPI_COMM_WORLD);
#endif

   ww *= 0.5 * HARTREE * vol1 * vol2;
#ifdef PARA
   MPI_Reduce(&ww, &rdummy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
   ww = rdummy;
#endif
   *energy = ww;

#ifdef VERBOSE
   if (rank == 0)
      printf("\n");
#endif

   return 0;
}

int terminate(int ierr)
{
#ifdef PARA
   if (ierr != 0)
      MPI_Abort(MPI_COMM_WORLD, ierr);
   MPI_Finalize();
#endif

   return ierr;
}

int main(int argc, char *argv[])
{
   bool flag_overlap;
   int threshold_power, coulomb_power, rank, size, ii, jj, kk, ifile = 0, nfiles;
#ifdef PARA
   int idummy, info;
#endif
   int ierr = 0, na = 0;
   double threshold, wfnorm[2], energy, isovalue[2], rcutoff, r2cutoff;
   char pfn[MAXCHAR] = "icm.inp";
   char ifn[2][MAXCHAR], iff[2][MAXCHAR];
   CARTESIAN mpo, mpn, offset;
   cell_info ci;
   scalar_field sf[2];

#ifdef PARA
   info = MPI_Init(&argc, &argv);
   if (info != MPI_SUCCESS)
   {
      fprintf(stderr, "\nMPI initialization failed\n\n");
      return terminate(-1);
   }
   info = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   info = MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
   rank = 0;
   size = 1;
#endif

   if (rank == 0)
      ierr = par_read(pfn, ifn[0], iff[0], ifn[1], iff[1], &threshold, &threshold_power, &coulomb_power, &mpo, &mpn);
#ifdef PARA
   info = MPI_Bcast(&ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif
   if (ierr != 0)
   {
      if (rank == 0)
	fprintf(stderr, "\nError: failed to read %s\n\n", pfn);
      return terminate(-1);
   }

   if(strcmp(ifn[1], "") == 0)
   {
      nfiles = 1;
   }
   else
   {
      nfiles = 2;
   }
#ifdef PARA
   info = MPI_Bcast(&nfiles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif

   for(ifile = 0; ifile < nfiles; ifile++) {
      if (rank == 0)
      {
	 if(ifile == 0)
	 {
	    ierr = par_read_cell(pfn, ci, '\0');
	 }
	 else
	 {
	    ierr = par_read_cell(pfn, ci, '2');
	 }
      }
#ifdef PARA
      info = MPI_Bcast(&ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif
      if (ierr != 0)
      {
	 if (rank == 0)
	    fprintf(stderr, "\nError: failed to read %s\n\n", pfn);
	 return terminate(-1);
      }

      if (rank == 0)
      {
	 if (strcmp(iff[ifile], "cube") == 0)
	 {
	    ierr = cub_read(ifn[ifile], &na, ci.sfv, sf[ifile]);
	 }
	 else if (strcmp(iff[ifile], "xsf") == 0)
	 {
	    ierr = xsf_read(ifn[ifile], &na, ci.sfv, sf[ifile]);
	 }
	 else
	 {
	    ierr = -2;
	 }

	 // we do not need atom info
	 if (as != NULL) delete [] as;
	 if (ap != NULL) delete [] ap;
      }
#ifdef PARA
      info = MPI_Bcast(&ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif
      if (ierr == -2)
      {
	 if (rank == 0)
	    fprintf(stderr, "\nError: unrecognized input format %s\n\n", iff[ifile]);
	 return terminate(-1);
      }
      else if (ierr == -1)
      {
	 if (rank == 0)
	    fprintf(stderr, "\nError: failed to read %s\n\n", ifn[ifile]);
	 return terminate(-1);
      }

      // FIXME: moving file 2 via uco or sco does not seem to work?

//      printf("  Lattice vectors from file %s:\n", ifn[ifile]);
//      for (ii = 0; ii < 3; ii++)
//	 printf("    v[%i] = %.6f %.6f %.6f (Ang)\n", ii, ci.sfv[ii].xx, ci.sfv[ii].yy, ci.sfv[ii].zz);

      if (rank == 0)
	 cell_set(sf[ifile].origin, ci.sfv, ci.uc, &ci.uco, ci.ucv, ci.ucf, ci.sc, &ci.sco, ci.scv, ci.scf);

      if (rank == 0)
	 if (ci.sc)
	    ierr = scalar_clone(ci.uco, ci.ucv, ci.sco, ci.scv, sf[ifile]);
#ifdef PARA
      info = MPI_Bcast(&ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif
      if (ierr != 0)
      {
	 if (rank == 0)
	    fprintf(stderr, "\nError: failed to build supercell\n\n");
	 return terminate(-1);
      }

      if (rank == 0)
      {
	 if (threshold < 1.0 - EPS9)
	 {
	    isovalue[ifile] = isovalue_scale(threshold_power, threshold, sf[ifile]);
	    ierr = scalar_trunc(isovalue[ifile], sf[ifile]);
	 }
	 else
	 {
	    isovalue[ifile] = 0.0;
	    ierr = 0;
	 }
      }

#ifdef PARA
      info = MPI_Bcast(&ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif

      if (ierr != 0)
      {
	 if (rank == 0)
	    fprintf(stderr, "\nError: failed in scalar_trunc\n\n");
	 return terminate(-1);
      }

      if (rank == 0)
	 ierr = scalar_norm(coulomb_power, &wfnorm[ifile], sf[ifile], ci.renorm);
#ifdef PARA
      info = MPI_Bcast(&ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif
      if (ierr != 0)
      {
	 if (rank == 0)
	    fprintf(stderr, "\nError: failed in scalar_norm\n\n");
	 return terminate(-1);
      }

      if (rank == 0)
      {
	 convert_to_bohr(sf[ifile]);
	 if(nfiles == 1)
	 {
	    ierr = mirror_plane(&mpo, &mpn, sf[0], sf[1]);
	 }
	 else
	 {
	    ierr = 0;
	 }
      }
#ifdef PARA
      info = MPI_Bcast(&ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#endif
      if (ierr != 0)
      {
	 if (rank == 0)
	    fprintf(stderr, "\nError: failed in mirror_plane\n\n");
	 return terminate(-1);
      }

#ifdef PARA
      info = MPI_Bcast(&sf[ifile].ni, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      info = MPI_Bcast(&sf[ifile].nj, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      info = MPI_Bcast(&sf[ifile].nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

      info = MPI_Bcast(&sf[ifile].origin, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
      info = MPI_Bcast(sf[ifile].stepv, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);

      if (rank != 0)
      {
	 sf[ifile].scalar = allocate_scalar(sf[ifile].ni, sf[ifile].nj, sf[ifile].nk);
      }

      for (ii = 0; ii < sf[ifile].ni; ii++)
	 for (jj = 0; jj < sf[ifile].nj; jj++)
	    info = MPI_Bcast(sf[ifile].scalar[ii][jj], sf[ifile].nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
#endif
   }

   if(nfiles == 1)
   {
#ifdef PARA
      info = MPI_Bcast(&sf[1].origin, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
      info = MPI_Bcast(sf[1].stepv, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
#endif

      sf[1].ni = sf[0].ni;
      sf[1].nj = sf[0].nj;
      sf[1].nk = sf[0].nk;
      sf[1].scalar = sf[0].scalar;
   }

   if (rank == 0)
   {
      printf("\n");
      for(ifile = 0; ifile < 2; ifile++)
      {
	 printf("  Parameters for grid %i:\n", ifile+1);
	 printf("    grid  = %i %i %i\n", sf[ifile].ni, sf[ifile].nj, sf[ifile].nk);
	 printf("    Origin and step vectors\n");
	 printf("    o%i    = %.6f %.6f %.6f (Bohr)\n", ifile+1, sf[ifile].origin.xx, sf[ifile].origin.yy, sf[ifile].origin.zz);
	 for (ii = 0; ii < 3; ii++)
	    printf("    v%i[%i] = %.6f %.6f %.6f (Bohr)\n", ifile+1, ii, sf[ifile].stepv[ii].xx, sf[ifile].stepv[ii].yy, sf[ifile].stepv[ii].zz);
	 if(ifile < nfiles)
	 {
	    printf("    wfn isovalue  = %.15f\n", isovalue[ifile]);
	    printf("    wfn norm      = %.15f\n", wfnorm[ifile]);
	 }
	 printf("\n");
      }
   }

   set_cutoff(sf[0].stepv, sf[1].stepv, &rcutoff, &r2cutoff);

   if(nfiles == 1)
   {
      flag_overlap = check_overlap(rank, sf[0], sf[1], rcutoff);
   }
   else
   {
      // TODO: new algorithm to check this.
      flag_overlap = true;
   }

   if (flag_overlap)
   {
      if(rank == 0)
      {
	 printf("Wavefunction overlaps with its image. Using Monte Carlo averaging.\n");
      }
      offset = grid_offset(sf[0], sf[1]);
      rand_init();
      mc_average(rank, size, &offset, sf[0].stepv, sf[1].stepv);
   }
   else
   {
      for (ii = -NCELL; ii <= NCELL; ii++)
      for (jj = -NCELL; jj <= NCELL; jj++)
      for (kk = -NCELL; kk <= NCELL; kk++)
         ravg[ii+NCELL][jj+NCELL][kk+NCELL] = -1.0;
   }

   ierr = coulomb_integral(rank, size, &offset, r2cutoff, &energy, sf[0], sf[1]);
#ifdef PARA
   info = MPI_Allreduce(&ierr, &idummy, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
   ierr = idummy;
#endif
   if (ierr != 0)
   {
      if (rank == 0)
	fprintf(stderr, "\nError: failed in coulomb_integral with ierr = %i\n\n", ierr);
      return terminate(-1);
   }
   if (rank == 0)
      printf(" coulomb integral = %.15f eV\n\n", energy);

   deallocate_scalar(sf[0].scalar, sf[0].ni, sf[0].nj);

   if (nfiles == 2)
   {
      deallocate_scalar(sf[1].scalar, sf[1].ni, sf[1].nj);
   }

   return terminate(0);
}

