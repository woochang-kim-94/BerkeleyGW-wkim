!============================================================================
!
! Routines:
!
! (1) mtxel_occ()        Originally By ?         Last Modified 8/15/2015 (FHJ)
!
!     Adapted from subroutine mtxel. Subroutine computes required matrix
!     elements that involve only occupied states, of the form <nk|exp{i g.r}|nk>
!     They are required for COHSEX and exact CH calculations.
!
!     input   n,m                band indices
!     input   ncoul              number of matrix elements required
!     input   isrtrq             index array for g-vectors in
!                                <nk|exp{i g.r}|nk>
!     output  aqs                matrix elements required
!
!============================================================================

#include "f_defs.h"

module mtxel_occ_m

  use fftw_m
  use global_m
  use misc_m
  implicit none

  private

  public :: &
    mtxel_occ

contains

subroutine mtxel_occ(n,m,gvec,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: n, m
  type (gspace), intent(in) :: gvec
  type (wfnkstates), intent(in) :: wfnk
  type (kpoints), intent(inout) :: kp
  integer, intent(in) :: ncoul
  integer, intent(in) :: isrtrq(gvec%ng)
  SCALAR, intent(out) :: aqs(ncoul)
  integer, intent(in) :: ispin

!-------------------------
! If we are using FFT to calculate matrix elements...

! We use FFT to compute <u_nk|e^(iG.r)|u_nk> elements where
! u_nk is the periodic part of the wave function.
! The calculation is done in real space, and integration over
! the grid is replaced by the sum over the grid points p:
!
! <u_nk|e^(iG.r)|u_nk>  =
!     Volume/Np * sum_p { conj(u_nk(p))*e^(iG.p)*u_nk(p) }
!
! Since u_nk(p) = Volume^-0.5 * sum_G { cnk(G)*e^(iG.p) },
! and FFT is defined as FFT(cnk,+,p) = sum_G { cnk(G)*e^{+iG.p} },
! we must compute
!
! <u_nk|e^(iG.r)|u_nk>
!   = 1/Np * sum_p { conj(FFT(cnk,+,p))*e^(iG.p)*FFT(cnk,+,p) }
!   = 1/Np * FFT(conj(FFT(cnk,+,:)).*FFT(cnk,+,:),+,G)
!
! where .* is a point by point multiplication on the grid

  complex(DPC), dimension(:,:,:), allocatable :: fftbox1,fftbox2
  integer :: jsp
  integer, dimension(3) :: Nfft
  real(DP) :: scale
  SCALAR, dimension(:), allocatable :: tmparray

  PUSH_SUB(mtxel_occ)

! Compute size of FFT box we need and scale factor

  call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)

! Allocate FFT boxes

  SAFE_ALLOCATE(fftbox1, (Nfft(1),Nfft(2),Nfft(3)))
  SAFE_ALLOCATE(fftbox2, (Nfft(1),Nfft(2),Nfft(3)))

! Put the data for band n into FFT box 1 and do the FFT,zk(:,1)

  SAFE_ALLOCATE(tmparray, (ncoul))

  do jsp = ispin,ispin*kp%nspinor

    call put_into_fftbox(wfnk%nkpt,wfnk%zk((n-1)*wfnk%nkpt+1:,jsp),gvec%components,wfnk%isrtk,fftbox1,Nfft)
    call do_FFT(fftbox1,Nfft,1)

! We need the complex conjugate of the |nk> band actually

    call conjg_fftbox(fftbox1,Nfft)

! Now we get the matrix elements:
!  Get n wave function and put it into box 2,
!  do FFT,
!  multiply by box1 contents,
!  do FFT again,
!  and extract the resulting matrix elements

    call put_into_fftbox(wfnk%nkpt,wfnk%zk((m-1)*wfnk%nkpt+1:,jsp),gvec%components,wfnk%isrtk,fftbox2,Nfft)
    call do_FFT(fftbox2,Nfft,1)
    call multiply_fftboxes(fftbox1,fftbox2,Nfft)
    call do_FFT(fftbox2,Nfft,1)
    call get_from_fftbox(ncoul,tmparray,gvec%components,isrtrq,fftbox2,Nfft,scale)
    if (kp%nspinor.eq.1 .or. jsp.eq. 1) then
      aqs(:) = tmparray(:)
    else
      aqs(:) = aqs(:) + tmparray(:)
    endif

  enddo

  SAFE_DEALLOCATE(tmparray)

! We are done, so deallocate FFT boxes

  SAFE_DEALLOCATE(fftbox1)
  SAFE_DEALLOCATE(fftbox2)
  
  POP_SUB(mtxel_occ)
  
  return
end subroutine mtxel_occ

end module mtxel_occ_m
