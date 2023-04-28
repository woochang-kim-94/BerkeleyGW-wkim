!============================================================================
!
! Routines:
!
! (1) checkgriduniformity()  Originally by JRD  Last Modified: 3/16/2011 (das)
!
!   To obtain the correct result for integrals in reciprocal space, the k+G
!   sampling should be similar in each direction. This routine determines the
!   sampling and writes a warning if it is too non-uniform.
!
!   gsm: xgrid = sig%qgrid in Sigma/main.f90, kp%kgrid in BSE/intkernel.f90
!
!   FHJ: We now simply compute the SVD decomposition of bvec * kgrid and
!   compare the largest/smallest singular vectors.
!
!   FHJ: Temporarily disabling this check. I am not sure what we do here is
!   meaningful at all, since we never construct the Voronoi cell.
!
!   Truncation:
!   icutv = 0, none        : 0D truncation, compare x, y, and z
!   icutv = 2, spherical   : 3D truncation, no relevant directions
!   icutv = 4, cell_wire   : 2D truncation, z is uniform by itself automatically
!   icutv = 5, cell_box    : 3D truncation, no relevant directions
!   icutv = 6, cell_slab   : 1D truncation, compare x and y
!
!============================================================================

#include "f_defs.h"

module checkgriduniformity_m

  use global_m
  implicit none

  private

  public :: checkgriduniformity

contains

subroutine checkgriduniformity(xgrid, crys, icutv)
  integer, intent(in) :: xgrid(3) !< k-/q-point grid
  type(crystal), intent(in) :: crys
  integer, intent(in) :: icutv !< code for type of truncation

  real(DP), parameter :: TOL_Ratio = 10.0d0
  integer, parameter :: lwork=20
  integer :: maxdir, info, jj
  real(DP) :: bb(3,3), singvals(3), dummy(1), work(lwork)
  character(len=32) :: sfmt

  PUSH_SUB(checkgriduniformity)

  ! FHJ: not sure this routine makes any sense.
  POP_SUB(checkgriduniformity)
  return

  if(icutv==TRUNC_NONE) then
    maxdir = 3
  else if(icutv==TRUNC_SLAB) then
    maxdir = 2
  else
    write(6,'(1x,a)') "No k+G sampling uniformity to check for the selected truncation scheme."
    POP_SUB(checkgriduniformity)
    return
  endif

  ! Compute miniBZ lattice vectors
  do jj = 1, 3
    bb(:,jj) = crys%blat * crys%bvec(:,jj) / dble(xgrid(jj))
  enddo

  ! Compute SVD of miniBZ lattice vectors
  singvals(:) = 0d0
  call dgesvd('N', 'N', maxdir, maxdir, bb(1:maxdir,1:maxdir), maxdir, &
    singvals, dummy, 1, dummy, 1, work, lwork, info)
  if (info/=0) then
    write(0,*) 'ERROR: call to dgesvd returned info=', info
    call die('Call to dgesvd failed')
  endif

  ! Singular vectors should give us the rough bounding box of the miniBZ.
  ! Normalize them relative to the largest one.
  singvals(1:maxdir) = singvals(1:maxdir) / singvals(1)
  write(sfmt,'(a,i1,a)') "(1x,a,", maxdir,"(1x,f8.6))"
  !This is in Cartesian coordinates
  write(6,trim(sfmt)) "Lenghts of the semi-axes of the Voronoi cell around each k+G point:", &
    singvals(1:maxdir)

  if ((singvals(1)/singvals(maxdir)) > TOL_RATIO) then
    write(0,'(/a)') "WARNING: detected non-uniform bounding box for the k+G Voronoi cells."
    write(0,'(a/)') "Check that the unit cell definition and k-/q-point samplings are reasonable."
  endif

  ! Note: for a cell with a small angle between two equivalent lattice vectors, it will
  ! be impossible to satisfy this criterion without using a different k-grid in those
  ! two equivalent directions, breaking crystal symmetry and possibly causing other problems.
  ! It is unclear what should be done in such a case. --DAS

  ! Do nothing at all. Our MC algorithm should capture nearly extreme cases --FHJ

  POP_SUB(checkgriduniformity)

  return

end subroutine checkgriduniformity

end module checkgriduniformity_m
