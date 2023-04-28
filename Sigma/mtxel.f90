!============================================================================
!
! Routines:
!
! (1) mtxel()   Originally By ?         Last Modified 7/8/2008 (JRD)
!
!     Subroutine computes required matrix elements
!     of the form <nn,k|exp{i(q+G).r}|n1,k-q> = M_{nn,n1}(k-q,q,G)
!
!     FHJ: Note that eqns. (20-24) of the BerkleyGW arxiv paper use the
!     quantity [M_{n``,n}(k,-q,-G)]^*, which is the same as what we compute
!     here, M_{nn,n1}(k-q,q,G) = aqs(G,n1). The indices in the arxiv paper and
!     in the code are related, respectively, by:
!       - n and n` <-> nn
!       - n`` <-> n1
!
!     input   nn                 band index for "outer" band
!     input   ncoul              number of matrix elements required
!     input   isrtrq             index array for g-vectors in
!                                <nk|exp(i(q+g).r)|n1k-q>
!     output  aqs                matrix elements required
!
!============================================================================

#include "f_defs.h"

module mtxel_m

  use accel_memory_m
  use accel_fft_m, only: accel_mtxel_sig, &
                                 accel_zero_box, accel_put_into_fftbox, accel_box_multiply, accel_get_from_fftbox, &
                                 accel_run_fft
  use algos_sigma_m
  use fftw_m
  use global_m
  use misc_m

  implicit none

  private

  public :: &
    mtxel

contains

subroutine mtxel(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: nn
  type (gspace), intent(in) :: gvec
  type (wfnkqstates), intent(in) :: wfnkq
  type (wfnkstates), intent(in) :: wfnk
  integer, intent(in) :: ncoul
  integer, intent(in) :: isrtrq(gvec%ng)
  SCALAR, intent(out) :: aqs(ncoul,peinf%ntband_max)
  integer, intent(in) :: ispin
  type (kpoints), intent(inout) :: kp

  PUSH_SUB(mtxel)

  select case (mtxel_algo)
  case (CPU_ALGO)
    call mtxel_cpu(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  case (OPENACC_ALGO)
    call mtxel_openacc(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  case (OMP_TARGET_ALGO)
    call mtxel_omp_target(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  case default
    call die("Invald algorithm for mtxel", only_root_writes = .true.)
  end select

  POP_SUB(mtxel)

  return
end subroutine mtxel

subroutine mtxel_cpu(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: nn
  type (gspace), intent(in), target :: gvec
  type (wfnkqstates), intent(in), target :: wfnkq
  type (wfnkstates), intent(in), target :: wfnk
  integer, intent(in), target :: isrtrq(gvec%ng)
  type (kpoints), intent(inout) :: kp
  integer, intent(in) :: ncoul
  SCALAR, intent(out) :: aqs(ncoul,peinf%ntband_max)
  integer, intent(in) :: ispin

  integer :: n1,jsp

! We use FFT to compute <u_nn,k|e^(iG.r)|u_n1,k-q> elements where
! u_nk is the periodic part of the wave function.
! The calculation is done in real space, and integration over
! the grid is replaced by the sum over the grid points p:
!
! <u_nn,k|e^(iG.r)|u_n1,k-q>  =
!     Volume/Np * sum_p { conj(u_nn,k(p))*e^(iG.p)*u_n1k-q(p) }
!
! Since u_nk(p) = Volume^-0.5 * sum_G { cnk(G)*e^(iG.p) },
! and FFT is defined as FFT(cnk,+,p) = sum_G { cnk(G)*e^{+iG.p} },
! we must compute
!
! <u_nn,k|e^(iG.r)|u_n1,k-q>
!   = 1/Np * sum_p { conj(FFT(c_nn k,+,p))*e^(iG.p)*FFT(c_n1 k-q,+,p) }
!   = 1/Np * FFT(conj(FFT(c_nn k,+,:)).*FFT(c_n1 k-q,+,:),+,G)
!
! where .* is a point by point multiplication on the grid

  complex(DPC), dimension(:,:,:), allocatable :: fftbox1,fftbox2
  integer, dimension(3) :: Nfft
  real(DP) :: scale
  SCALAR, dimension(:), allocatable :: tmparray

! Compute size of FFT box we need and scale factor

  PUSH_SUB(mtxel_cpu)

  call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)

! Allocate FFT boxes

  SAFE_ALLOCATE(fftbox1, (Nfft(1),Nfft(2),Nfft(3)))
  SAFE_ALLOCATE(fftbox2, (Nfft(1),Nfft(2),Nfft(3)))

! Put the data for band nn into FFT box 1 and do the FFT,zk(:,1)

  SAFE_ALLOCATE(tmparray, (ncoul))

  do jsp = ispin,ispin*kp%nspinor

    call put_into_fftbox(wfnk%nkpt,wfnk%zk((nn-1)*wfnk%nkpt+1:,jsp),gvec%components,wfnk%isrtk,fftbox1,Nfft)
    call do_FFT(fftbox1,Nfft,1)
    ! We need the complex conjugate of u_{nn,k)(r) for the cross correlation
    call conjg_fftbox(fftbox1,Nfft)

! Now we loop over the n1 states and get the matrix elements:
!  Get n1 wave function and put it into box 2,
!  do FFT, get u_{n1k-q}(r)
!  multiply by box1 contents, get 
!  do FFT again,
!  and extract the resulting matrix elements
! Now we loop over the n1 states and get the matrix elements:
! 1. Get conduction wave function and put it into box 2,
! 2. do FFT, get u_{n1,k-q}(r)
! 3. multiply by box1 contents, get F(r) = [u_{nn,k)(r)]^* u_{n1,k-q}(r)
! 4. do FFT again, and extract the resulting matrix elements,
!    <nn,k|exp{i(q+G).r}|n1,k-q>

    do n1=1,peinf%ntband_node
      call put_into_fftbox(wfnkq%nkpt,wfnkq%zkq((n1-1)*wfnkq%nkpt+1:,jsp),gvec%components,wfnkq%isrtkq,fftbox2,Nfft)
      call do_FFT(fftbox2,Nfft,1)
      call multiply_fftboxes(fftbox1,fftbox2,Nfft)
      call do_FFT(fftbox2,Nfft,1)
      call get_from_fftbox(ncoul,tmparray,gvec%components,isrtrq,fftbox2,Nfft,scale)
      if (kp%nspinor.eq.1 .or. jsp.eq. 1) then
        aqs(:,n1) = tmparray(:)
      else
        aqs(:,n1) = aqs(:,n1) + tmparray(:)
      endif
    enddo
  enddo

  SAFE_DEALLOCATE(tmparray)

! We are done, so deallocate FFT boxes

  SAFE_DEALLOCATE(fftbox1)
  SAFE_DEALLOCATE(fftbox2)

  POP_SUB(mtxel_cpu)

  return
end subroutine mtxel_cpu

subroutine mtxel_openacc(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: nn
  type (gspace), intent(in):: gvec
  type (wfnkqstates), intent(in) :: wfnkq
  type (wfnkstates), intent(in) :: wfnk
  integer, intent(in) :: isrtrq(gvec%ng)
  type (kpoints), intent(inout) :: kp
  integer, intent(in) :: ncoul
  SCALAR, intent(out) :: aqs(ncoul,peinf%ntband_max)
  integer, intent(in) :: ispin
#ifdef OPENACC
  integer, dimension(3) :: Nfft
  real(DP)   :: scale
  integer :: inx_start, inx_end
  integer :: n1_global, iblock, n_start, n_end
  integer :: n1, jsp

  complex(DPC), dimension(:,:,:), allocatable :: fftbox1

  PUSH_SUB(mtxel_openacc)

! Compute size of FFT box we need and scale factor
  call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)

  !$acc enter data copyin(Nfft)
  ! copy in arrays
  accel_mtxel_sig%gvec_comp        = gvec%components
  call accel_update_to(accel_mtxel_sig%gvec_comp, mtxel_algo)
  accel_mtxel_sig%gvec_isrtx_k     = wfnk%isrtk
  call accel_update_to(accel_mtxel_sig%gvec_isrtx_k, mtxel_algo)
  accel_mtxel_sig%gvec_isrtx_q     = wfnkq%isrtkq
  call accel_update_to(accel_mtxel_sig%gvec_isrtx_q, mtxel_algo)
  accel_mtxel_sig%gvec_isrtx_mtxel = isrtrq
  call accel_update_to(accel_mtxel_sig%gvec_isrtx_mtxel, mtxel_algo)

  do jsp = ispin,ispin*kp%nspinor
    inx_start = (nn-1)*wfnk%nkpt+1
    inx_end   = inx_start + wfnk%nkpt - 1
    accel_mtxel_sig%eqp_vec(:) = CMPLX(0.0D+00,0.0D+00)
    accel_mtxel_sig%eqp_vec(1:wfnk%nkpt) = wfnk%zk(inx_start:inx_end,jsp)
    !$acc update device( accel_mtxel_sig%eqp_vec ) async(1)

    ! zerofy FFT box
    call accel_zero_box(accel_mtxel_sig%eqp_box, Nfft, mtxel_algo, queue=1)
    ! put
    call accel_put_into_fftbox( wfnk%nkpt, accel_mtxel_sig%eqp_vec, accel_mtxel_sig%gvec_comp, &
                              accel_mtxel_sig%gvec_isrtx_k, accel_mtxel_sig%eqp_box, Nfft, 1.0D+00, mtxel_algo, queue=1 )
    ! run FFT
    call accel_run_fft( accel_mtxel_sig%eqp_box, accel_mtxel_sig%eqp_plan, 1, mtxel_algo )
    !$acc wait(1)

    ! Inner loop over blocks
    do iblock = 1, peinf%ntband_node, accel_mtxel_sig%mtxel_band_block_size
      n_start = iblock
      n_end   = MIN(n_start + accel_mtxel_sig%mtxel_band_block_size - 1 , peinf%ntband_node)

      n1 = 0
      do n1_global = n_start, n_end
        n1 = n1 + 1

        inx_start = (n1_global-1)*wfnkq%nkpt+1
        inx_end   = inx_start + wfnkq%nkpt - 1
        accel_mtxel_sig%bands_vec(n1)%vec(:) = CMPLX(0.0D+00,0.0D+00)
        accel_mtxel_sig%bands_vec(n1)%vec(1:wfnkq%nkpt) = wfnkq%zkq(inx_start:inx_end,jsp)
        !$acc update device( accel_mtxel_sig%bands_vec(n1)%vec ) async(n1+1)
        accel_mtxel_sig%mtxel_vec(n1)%vec(:) = CMPLX(0.0D+00,0.0D+00)
        !$acc update device( accel_mtxel_sig%mtxel_vec(n1)%vec ) async(n1+1)

        ! zerofy FFT box
        call accel_zero_box(accel_mtxel_sig%mtxel_box(n1)%box, Nfft, mtxel_algo, queue=(n1+1))

        ! put
        call accel_put_into_fftbox(wfnkq%nkpt, accel_mtxel_sig%bands_vec(n1)%vec, accel_mtxel_sig%gvec_comp, accel_mtxel_sig%gvec_isrtx_q, &
                                 accel_mtxel_sig%mtxel_box(n1)%box, Nfft, 1.0D+00, mtxel_algo, queue=(n1+1))

        ! FFT
        call accel_run_fft(accel_mtxel_sig%mtxel_box(n1)%box, accel_mtxel_sig%mtxel_plan(n1), 1, mtxel_algo)

        ! multiply
        call accel_box_multiply(accel_mtxel_sig%mtxel_box(n1)%box, accel_mtxel_sig%eqp_box, Nfft, 1, mtxel_algo, queue=(n1+1) )

        ! FFT
        call accel_run_fft( accel_mtxel_sig%mtxel_box(n1)%box, accel_mtxel_sig%mtxel_plan(n1), 1, mtxel_algo)

        ! get
        call accel_get_from_fftbox(ncoul, accel_mtxel_sig%mtxel_vec(n1)%vec, accel_mtxel_sig%gvec_comp, accel_mtxel_sig%gvec_isrtx_mtxel, &
                                 accel_mtxel_sig%mtxel_box(n1)%box, Nfft, scale, mtxel_algo, queue=(n1+1) )

        ! copy to host
        !$acc update self( accel_mtxel_sig%mtxel_vec(n1)%vec ) async(n1+1)

      end do ! n1

      n1 = 0
      do n1_global = n_start, n_end
        n1 = n1 + 1
        !$acc wait(n1+1)

        if (kp%nspinor.eq.1 .or. jsp.eq. 1) then
          aqs(:,n1_global) = accel_mtxel_sig%mtxel_vec(n1)%vec(:)
        else
          aqs(:,n1_global) = aqs(:,n1_global) + accel_mtxel_sig%mtxel_vec(n1)%vec(:)
        endif
      end do ! n1_global
    end do ! iblock
  end do ! jsp
  !$acc exit data delete(Nfft)

#else
  PUSH_SUB(mtxel_openacc)

  call die("OpenACC version of mtxel requested, but OpenACC not compiled&
           & into this executable", only_root_writes = .true.)
#endif

  POP_SUB(mtxel_openacc)

  return
end subroutine mtxel_openacc

subroutine mtxel_omp_target(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: nn
  type (gspace), intent(in):: gvec
  type (wfnkqstates), intent(in) :: wfnkq
  type (wfnkstates), intent(in) :: wfnk
  integer, intent(in) :: isrtrq(gvec%ng)
  type (kpoints), intent(inout) :: kp
  integer, intent(in) :: ncoul
  SCALAR, intent(out) :: aqs(ncoul,peinf%ntband_max)
  integer, intent(in) :: ispin
#ifdef OMP_TARGET
  integer, dimension(3) :: Nfft
  real(DP)   :: scale
  integer :: inx_start, inx_end
  integer :: n1_global, iblock, n_start, n_end
  integer :: n1, jsp

  complex(DPC), dimension(:,:,:), allocatable :: fftbox1

  PUSH_SUB(mtxel_omp_target)

! Compute size of FFT box we need and scale factor
  call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)

  ! copy in arrays
  accel_mtxel_sig%gvec_comp        = gvec%components
  call accel_update_to(accel_mtxel_sig%gvec_comp, mtxel_algo)
  accel_mtxel_sig%gvec_isrtx_k     = wfnk%isrtk
  call accel_update_to(accel_mtxel_sig%gvec_isrtx_k, mtxel_algo)
  accel_mtxel_sig%gvec_isrtx_q     = wfnkq%isrtkq
  call accel_update_to(accel_mtxel_sig%gvec_isrtx_q, mtxel_algo)
  accel_mtxel_sig%gvec_isrtx_mtxel = isrtrq
  call accel_update_to(accel_mtxel_sig%gvec_isrtx_mtxel, mtxel_algo)

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
  !XXX Need to fix the nowait and stream stuff XXX!
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!

  do jsp = ispin,ispin*kp%nspinor
    inx_start = (nn-1)*wfnk%nkpt+1
    inx_end   = inx_start + wfnk%nkpt - 1
    accel_mtxel_sig%eqp_vec(:) = CMPLX(0.0D+00,0.0D+00)
    accel_mtxel_sig%eqp_vec(1:wfnk%nkpt) = wfnk%zk(inx_start:inx_end,jsp)
    call accel_update_to(accel_mtxel_sig%eqp_vec, mtxel_algo)

    ! zerofy FFT box
    call accel_zero_box(accel_mtxel_sig%eqp_box, Nfft, mtxel_algo, queue=1)
    ! put
    call accel_put_into_fftbox( wfnk%nkpt, accel_mtxel_sig%eqp_vec, accel_mtxel_sig%gvec_comp, &
                              accel_mtxel_sig%gvec_isrtx_k, accel_mtxel_sig%eqp_box, Nfft, 1.0D+00, mtxel_algo, queue=1 )
    ! run FFT
    call accel_run_fft( accel_mtxel_sig%eqp_box, accel_mtxel_sig%eqp_plan, 1, mtxel_algo )
    !XX !need to be fixed in OMP_TARGET $acc wait(1)

    ! Inner loop over blocks
    do iblock = 1, peinf%ntband_node, accel_mtxel_sig%mtxel_band_block_size
      n_start = iblock
      n_end   = MIN(n_start + accel_mtxel_sig%mtxel_band_block_size - 1 , peinf%ntband_node)

      n1 = 0
      do n1_global = n_start, n_end
        n1 = n1 + 1

        inx_start = (n1_global-1)*wfnkq%nkpt+1
        inx_end   = inx_start + wfnkq%nkpt - 1
        accel_mtxel_sig%bands_vec(n1)%vec(:) = CMPLX(0.0D+00,0.0D+00)
        accel_mtxel_sig%bands_vec(n1)%vec(1:wfnkq%nkpt) = wfnkq%zkq(inx_start:inx_end,jsp)
        call accel_update_to(accel_mtxel_sig%bands_vec(n1)%vec, mtxel_algo)
        accel_mtxel_sig%mtxel_vec(n1)%vec(:) = CMPLX(0.0D+00,0.0D+00)
        call accel_update_to(accel_mtxel_sig%mtxel_vec(n1)%vec, mtxel_algo)

        ! zerofy FFT box
        call accel_zero_box(accel_mtxel_sig%mtxel_box(n1)%box, Nfft, mtxel_algo, queue=(n1+1))

        ! put
        call accel_put_into_fftbox(wfnkq%nkpt, accel_mtxel_sig%bands_vec(n1)%vec, accel_mtxel_sig%gvec_comp, accel_mtxel_sig%gvec_isrtx_q, &
                                 accel_mtxel_sig%mtxel_box(n1)%box, Nfft, 1.0D+00, mtxel_algo, queue=(n1+1))

        ! FFT
        call accel_run_fft(accel_mtxel_sig%mtxel_box(n1)%box, accel_mtxel_sig%mtxel_plan(n1), 1, mtxel_algo)

        ! multiply
        call accel_box_multiply(accel_mtxel_sig%mtxel_box(n1)%box, accel_mtxel_sig%eqp_box, Nfft, 1, mtxel_algo, queue=(n1+1) )

        ! FFT
        call accel_run_fft( accel_mtxel_sig%mtxel_box(n1)%box, accel_mtxel_sig%mtxel_plan(n1), 1, mtxel_algo)

        ! get
        call accel_get_from_fftbox(ncoul, accel_mtxel_sig%mtxel_vec(n1)%vec, accel_mtxel_sig%gvec_comp, accel_mtxel_sig%gvec_isrtx_mtxel, &
                                 accel_mtxel_sig%mtxel_box(n1)%box, Nfft, scale, mtxel_algo, queue=(n1+1) )

        ! copy to host
        call accel_update_from(accel_mtxel_sig%mtxel_vec(n1)%vec, mtxel_algo)
      end do ! n1

      n1 = 0
      do n1_global = n_start, n_end
        n1 = n1 + 1
        !XXX need to be fixed !$acc wait(n1+1)

        if (kp%nspinor.eq.1 .or. jsp.eq. 1) then
          aqs(:,n1_global) = accel_mtxel_sig%mtxel_vec(n1)%vec(:)
        else
          aqs(:,n1_global) = aqs(:,n1_global) + accel_mtxel_sig%mtxel_vec(n1)%vec(:)
        endif
      end do ! n1_global
    end do ! iblock
  end do ! jsp
#else
  PUSH_SUB(mtxel_omp_target)

  call die("OpenMP Target version of mtxel requested, but OpenMP Target&
           & not compiled into this executable", only_root_writes = .true.)
#endif

  POP_SUB(mtxel_omp_target)
  return
end subroutine mtxel_omp_target

end module mtxel_m
