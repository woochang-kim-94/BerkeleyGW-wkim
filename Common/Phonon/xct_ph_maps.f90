!==============================================================================
!
! module xct_ph_maps_m
!
! Various mapping functions between exciton eigenvectors indices
! and electron-phonon coupling matrix indices.
!
! Originally by GKA (2018)
!
!==============================================================================

#include "f_defs.h"

module xct_ph_maps_m

use global_m
use gkq_m
use evecs_m

implicit none

public

contains


!> Get the indices mapping k to k+q within a list of k-point.
!! The k-points must be a regular grid and qpt must be part of that grid.
subroutine get_kq_map(kq_map, kpts, qpt, nk)
  integer, allocatable, intent(out) :: kq_map(:)    !< kq_map(nk) such that
                                !! kpts(:,kq_map(ii)) = kpts(:,ii) + qpt(:)
  integer, intent(in) :: nk           !< number of k-points
  real(dp), intent(in) :: kpts(3,nk)  !< number of k-points
  real(dp), intent(in) :: qpt(3)      !< q-points

  real(DP) :: delta,tol, Gi
  real(DP) :: kk(3), kq(3), dk(3)
  integer :: ik, ikq, ii

  PUSH_SUB(get_kq_map)

  ! GKA: This routine could be optimized for large lists of k-points
  !      by arranging the k-points by shells

  tol = TOL_Small

  SAFE_ALLOCATE(kq_map, (nk))

  do ik=1,nk

    kq = kpts(:,ik) + qpt
    do ii=1,3
      if (abs(kq(ii)) .ge. 1.0_dp) then
        Gi = floor(abs(kq(ii)))
        kq(ii) = kq(ii) - sign(Gi, kq(ii))
      end if
    enddo

    ikq=0
    delta=0.1d0
    do while((delta.gt.tol).and.(ikq.lt.nk))
      ikq=ikq+1      
      dk = kq - kpts(:,ikq)

      do ii=1,3
        if (abs(dk(ii)) .ge. 1.0_dp) then
          Gi = floor(abs(dk(ii)))
          dk(ii) = dk(ii) - sign(Gi, dk(ii))
        end if
      enddo

      delta=sqrt(sum(dk(1:3) ** 2))
    enddo

    if (delta.lt.tol) then
      kq_map(ik) = ikq
    else
      if(peinf%inode.eq.0) then
        write(0,*) 'Could not find point equivalent to ', (kq(ii),ii=1,3)
      endif
      call die('k-point not found.', only_root_writes = .true.)
    endif

  enddo

  POP_SUB(get_kq_map)

end subroutine get_kq_map


!> Compute the mapping between bands indices in gkq and bands in Avc
!! such that, for example, the index iv in evecs%Avc and v_map(iv) in gkq%g_nu
!! refer to the same band
subroutine get_band_maps(v_map, c_map, nv, nc, nocc)
  integer, allocatable, intent(out) :: v_map(:)   !< v_map(nv)
  integer, allocatable, intent(out) :: c_map(:)   !< c_map(nc)
  integer, intent(in) :: nv                       !< Number of valence bands in Avc
  integer, intent(in) :: nc                       !< Number of conduction bands in Avc
  integer, intent(in) :: nocc                     !< Number of occupied bands

  integer :: ic, iv

  PUSH_SUB(get_band_maps)

  SAFE_ALLOCATE(v_map, (nv))
  SAFE_ALLOCATE(c_map, (nc))

  do iv=1,nv
    v_map(iv) = nocc + 1 - iv
  end do

  do ic=1,nc
    c_map(ic) = nocc + ic
  end do

  POP_SUB(get_band_maps)

end subroutine get_band_maps

end module xct_ph_maps_m
