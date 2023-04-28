!==============================================================================
!
!  test_voro             Originally by FHJ     Last Modified 08/18/2014 (FHJ)
!
!    A simple Fortran program that tests the BerkeleyGW bindings for voro++.
!
!    This file is part of the BerkeleyGW package.
!
!==============================================================================

#include "f_defs.h"

program test_voro

  use global_m
  use tile_m
  implicit none

  integer, parameter :: nk = 3
  integer, parameter :: npts = nk**3 + 1
  real(DP) :: bdot(3,3), pts_crys(3,npts), vols(npts)
  integer :: umks(3,npts)
  integer :: ix, iy, iz, ipt
  
  PUSH_SUB(test_voro)

  !TODO: make this configurable!
  write(6,'(a)') "Testing libtile_voro - Fortran bindings."
 
 ipt = 0
  do ix=1,nk
    do iy=1,nk
      do iz=1,nk
        ipt = ipt + 1
        pts_crys(:,ipt) = (/ix,iy,iz/)/dble(nk)
        !pts_crys(:,ipt) = ((/ix,iy,iz/)/dble(nk))**2
        !pts_crys(:,ipt) = (/(ix-nk/2d0)**2,(iy-nk/2d0)**2,(iz-nk/2d0)**2/)/dble(nk)
      enddo
    enddo
  enddo
  ipt = ipt + 1
  pts_crys(:,ipt) = (/0.5d0,0.5d0,0.5d0/)

  bdot(:,:) = 0d0
  bdot(1,1) = 2d0
  bdot(2,2) = 1d0
  bdot(3,3) = 0.5d0
  bdot(1,1) = 1d0
  bdot(2,2) = 1d0
  bdot(3,3) = 1d0

  write(6,*)
  write(6,'(a)') 'Original points:'
  write(6,'(3(f12.9,1x))') pts_crys
  call pts2bz(bdot, pts_crys, umks)
  call get_kpts_volumes(bdot, pts_crys, vols)
  write(6,*)
  write(6,'(a)') 'Mapped points, volume, and umklapp vectors:'
  do ipt=1,npts
    write(6,'(3(f12.9,1x),2x,f12.9,2x,3(i6,1x))') pts_crys(:,ipt)+umks(:,ipt), vols(ipt), umks(:,ipt)
  enddo


  write(6,'(/,a)') "All Done!"

  POP_SUB(test_voro)

end program test_voro
