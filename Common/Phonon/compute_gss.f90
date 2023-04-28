!==============================================================================
!
! subroutine compute_gss
!
! Compute the exciton-phonon coupling matrix elements.
!
! Originally by GKA (2018)
!
!==============================================================================

#include "f_defs.h"

module compute_gss_m

  use blas_m
  use global_m
  implicit none

  private

  public :: &
    compute_gss

contains

!> Compute the exciton-phonon coupling matrix elements for one q-point
!! and a single pair of excitons.
subroutine compute_gss(gss_nu, gkq_nu, Avc, Avc_q, kq_map, v_map, c_map, &
                       ns, nk, nc, nv, mband, nband, nmode)
  integer, intent(in) :: ns     !< Number of spins
  integer, intent(in) :: nk     !< Number of k-points
  integer, intent(in) :: nv     !< Number of valence bands
  integer, intent(in) :: nc     !< Number of conduction bands
  integer, intent(in) :: mband  !< Number of bands at k in gkq
  integer, intent(in) :: nband  !< Number of bands at k+q in gkq
  integer, intent(in) :: nmode  !< Number of phonon modes
  SCALAR, intent(out) :: gss_nu(ns,nmode)                !< exciton-phonon couplig matrix
  SCALAR, intent(in) :: gkq_nu(mband,nband,nk,ns,nmode)  !< electron-phonon couplig matrix
  SCALAR, intent(in) :: Avc(ns,nv,nc,nk)    !< exciton eigenvectors at Q
  SCALAR, intent(in) :: Avc_q(ns,nv,nc,nk)  !< exciton eigenvectors at Q+q
  integer, intent(in) :: kq_map(nk) !< mapping between k and k+q
  integer, intent(in) :: v_map(nv)  !< mapping between valence bands in Avc and bands in gkq
  integer, intent(in) :: c_map(nc)  !< mapping between conduction bands in Avc and bands in gkq


  SCALAR, allocatable :: u_Q(:), u_Qq(:) ! Flat eigenvectors
  SCALAR, allocatable :: u_Q_el(:), u_Qq_el(:) ! Flat eigenvectors (nc) 
  SCALAR, allocatable :: u_Q_ho(:), u_Qq_ho(:) ! Flat eigenvectors (nv) 
  SCALAR, allocatable :: gkq_el(:,:,:) ! (nmode,nc,nc)
  SCALAR, allocatable :: gkq_ho(:,:,:) ! (nmode,nv,nv)
  SCALAR, allocatable :: gA_el(:), gA_ho(:)
  SCALAR :: fac_el, fac_ho
  integer :: iv1_g,iv2_g,ic1_g,ic2_g,imode
  integer :: is, ik, ikq, ic, icp, iv, ivp


  PUSH_SUB(compute_gss)

  SAFE_ALLOCATE(u_Q_el, (nc))
  SAFE_ALLOCATE(u_Qq_el, (nc))
  SAFE_ALLOCATE(gA_el, (nc))
  SAFE_ALLOCATE(gkq_el, (nmode,nc,nc))
  SAFE_ALLOCATE(u_Q_ho, (nv))
  SAFE_ALLOCATE(u_Qq_ho, (nv))
  SAFE_ALLOCATE(gA_ho, (nc))
  SAFE_ALLOCATE(gkq_ho, (nmode,nv,nv))

  gss_nu = ZERO

  fac_el =   ONE / nk
  fac_ho = - ONE / nk

  ! Outer matrix multiplication
  do is=1,ns  ! FIXME: This is bad because spin is the innermost index
    do ik=1,nk

      ikq = kq_map(ik)

      ! electron coupling term
      do iv=1,nv

        iv1_g = v_map(iv)

        do ic=1,nc

          ! Align flat vectors
          u_Q_el(ic) = Avc(is,iv,ic,ik)
          u_Qq_el(ic) = Avc_q(is,iv,ic,ik)

          ! Reorder matrix
          ic1_g = c_map(ic)
          do icp=1,nc
            ic2_g = c_map(icp)
            gkq_el(:,ic,icp) = gkq_nu(ic1_g,ic2_g,ik,is,:)  ! FIXME: This is bad because mode is the outermost index
          end do
        end do

        ! Perform vector-matrix-vector contraction
        ! by computing  g | A >  then < Aq | g | A >
        do imode=1,nmode
#ifdef CPLX
          call ZGEMM('n', 'n', nc, 1, nc, ONE, gkq_el(imode,:,:), nc, &
                     u_Q_el(:), nc, ZERO, gA_el(:), nc)
          call ZGEMM('t', 'n', 1, 1, nc, fac_el, u_Qq_el, nc, &
                     gA_el, nc, ONE, gss_nu(is,imode), 1)
#else
          call DGEMM('n', 'n', nc, 1, nc, ONE, gkq_el(imode,:,:), nc, &
                     u_Q_el(:), nc, ZERO, gA_el(:), nc)
          call DGEMM('t', 'n', 1, 1, nc, fac_el, u_Qq_el, nc, &
                     gA_el, nc, ONE, gss_nu(is,imode), 1)
#endif
        end do
      end do  ! iv

      ! hole coupling term
      do ic=1,nc

        ic1_g = c_map(ic)

        do iv=1,nv

          ! Align flat vectors
          u_Q_ho(iv) = Avc(is,iv,ic,ikq)
          u_Qq_ho(iv) = Avc_q(is,iv,ic,ik)

          ! Reorder matrix
          iv1_g = v_map(iv)
          do ivp=1,nv
            iv2_g = v_map(ivp)
            ! Note the inversion of iv1_g, iv2_g, compared to ic1_g, ic2_g
            gkq_el(:,iv,ivp) = gkq_nu(iv2_g,iv1_g,ik,is,:)  ! FIXME: This is bad because mode is the outermost index
          end do

        end do

        ! Perform vector-matrix-vector contraction
        ! by computing  g | A >  then < Aq | g | A >
        do imode=1,nmode
#ifdef CPLX
          call ZGEMM('n', 'n', nc, 1, nc, ONE, gkq_ho(imode,:,:), nc, &
                     u_Q_ho(:), nc, ZERO, gA_ho(:), nc)
          call ZGEMM('t', 'n', 1, 1, nc, fac_ho, u_Qq_ho, nc, &
                     gA_ho, nc, ONE, gss_nu(is,imode),1)
#else
          call DGEMM('n', 'n', nc, 1, nc, ONE, gkq_ho(imode,:,:), nc, &
                     u_Q_ho(:), nc, ZERO, gA_ho(:), nc)
          call DGEMM('t', 'n', 1, 1, nc, fac_ho, u_Qq_ho, nc, &
                     gA_ho, nc, ONE, gss_nu(is,imode),1)
#endif
        end do  ! imode
      end do  ! ic
    end do  ! ik
  end do  ! is


  SAFE_DEALLOCATE(u_Q_el)
  SAFE_DEALLOCATE(u_Qq_el)
  SAFE_DEALLOCATE(gA_el)
  SAFE_DEALLOCATE(gkq_el)
  SAFE_DEALLOCATE(u_Q_ho)
  SAFE_DEALLOCATE(u_Qq_ho)
  SAFE_DEALLOCATE(gA_ho)
  SAFE_DEALLOCATE(gkq_ho)

  POP_SUB(compute_gss)

end subroutine compute_gss


end module compute_gss_m
