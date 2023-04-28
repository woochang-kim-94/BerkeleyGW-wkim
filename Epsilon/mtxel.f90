!===========================================================================
!
! Routines()
!
! (1) mtxel()   Originally By (?)               Last Modified 5/6/2012 (FHJ)
!
!     Compute matrix elements (gme) for valence state iv with all
!     conduction bands and for all G-vectors.
!
!                <c, k, ispin|exp{-i(q+G).r}|v, k+q, ispin> = [M_{vc}(k, q, G)]^*
!
!     On exit, 
!       pol%eden(band, spin) = 1/(e_val-e_cond) = energy denominators
!       pol%gme(band, g-vector, spin) = plane wave matrix elements
!       pol%isrtx   orders the |G(i)|^2   i = 1, pol%nmtx
!       vwfn%isort  orders |qk+g|^2    (in vwfn type)
!
!       energies are apparently assumed in Rydbergs.
!
!===========================================================================

#include "f_defs.h"

module mtxel_m

  use algos_epsilon_m
  use accel_memory_m
  use accel_fft_m, only: accel_zero_box, accel_box_multiply, &
      accel_put_into_fftbox, accel_run_fft, accel_get_from_fftbox, accel_mtxel_eps, &
      accel_run_fft_batched, accel_box_multiply_batched, accel_zero_box_batched, &
      accel_put_into_fftbox_batched, accel_get_from_fftbox_batched
  use global_m
  use fftw_m
  use misc_m
  use lin_denominator_m
  use timing_m, only: timing => epsilon_timing

#if defined (OPENACC) || defined (OMP_TARGET)
  use openacc
  use cufft
#endif

  implicit none

  private

  public:: mtxel_init_FFT_cond, mtxel_free_FFT_cond, mtxel

contains

!> Precalculates the FFTs for all conduction bands
subroutine mtxel_init_FFT_cond(gvec, pol, cwfn, kp)
  type (gspace), intent(in):: gvec
  type (polarizability), intent(inout):: pol
  type (conduction_wfns), intent(inout):: cwfn
  type (kpoints), intent(inout):: kp

  real(DP):: scale
  integer, dimension(3):: Nfft
  integer:: j, offset_g, jsp, ict

  PUSH_SUB(mtxel_init_FFT_cond)

  if (kp%nspinor > 1) then
      call die('Epsilon on Steroids only supports one spin at the moment.')
  endif

  call timing%start(timing%opt_fft)
  call timing%start(timing%opt_fft_fft)
  call setup_FFT_sizes(pol%FFTgrid, Nfft, scale)
  SAFE_ALLOCATE(cwfn%wfn_fft, (Nfft(1), Nfft(2), Nfft(3), peinf%ncownactual))

  jsp = 1
! JRD XXX We should we be doing a many_fft call here or what is the point of allocating a bigger array
  do j = 1, peinf%ncownactual
    offset_g = (j-1)*cwfn%ngc

    call put_into_fftbox(cwfn%ngc, cwfn%zc(offset_g+1:,jsp), gvec%components, cwfn%isort, cwfn%wfn_fft(:,:,:,j), Nfft)

    call do_FFT(cwfn%wfn_fft(:,:,:,j), Nfft, 1)
  enddo

  call timing%stop(timing%opt_fft)
  call timing%stop(timing%opt_fft_fft)

  POP_SUB(mtxel_init_FFT_cond)
  return

end subroutine mtxel_init_FFT_cond

!> Frees ffts_cond buffer
subroutine mtxel_free_FFT_cond(cwfn)
  type (conduction_wfns), intent(inout):: cwfn

  PUSH_SUB(mtxel_free_FFT_cond)

  SAFE_DEALLOCATE_P(cwfn%wfn_fft)

  POP_SUB(mtxel_free_FFT_cond)
  return

end subroutine mtxel_free_FFT_cond

subroutine mtxel(iv, gvec, vwfn, cwfn, pol, ispin, irk, kp, kpq, rank_mtxel, kfact)
  integer, intent(in):: iv
  type (gspace), intent(in):: gvec
  type (valence_wfns), intent(in):: vwfn
  type (conduction_wfns), intent(in):: cwfn
  type (polarizability), intent(inout):: pol
  integer, intent(in):: ispin, irk
  type (kpoints), intent(inout):: kp, kpq
  integer, intent(in):: rank_mtxel
  real(DP), intent(in):: kfact

  PUSH_SUB(mtxel)

  select case (mtxel_algo)
  case (CPU_ALGO)
    call mtxel_cpu(iv, gvec, vwfn, cwfn, pol, ispin, irk, kp, kpq, &
                   rank_mtxel, kfact)
  case (OPENACC_ALGO)
    if (accel_mtxel_eps%n_ffts_per_batch .gt. 1) then
      call mtxel_accel_batched(iv, gvec, vwfn, cwfn, pol, ispin, irk, kp, kpq, &
                               rank_mtxel, kfact)
    else
      call mtxel_openacc(iv, gvec, vwfn, cwfn, pol, ispin, irk, kp, kpq, &
                         rank_mtxel, kfact)
    endif
  case (OMP_TARGET_ALGO)
    if (accel_mtxel_eps%n_ffts_per_batch .gt. 1) then
      call mtxel_accel_batched(iv, gvec, vwfn, cwfn, pol, ispin, irk, kp, kpq, &
                               rank_mtxel, kfact)
    else
      call mtxel_omp_target(iv, gvec, vwfn, cwfn, pol, ispin, irk, kp, kpq, &
                            rank_mtxel, kfact)
    end if
  case default
    call die("Invald algorithm for mtxel", only_root_writes = .true.)
  end select

  POP_SUB(mtxel)
end subroutine mtxel

subroutine mtxel_cpu(iv, gvec, vwfn, cwfn, pol, ispin, irk, kp, kpq, rank_mtxel, kfact)
  integer, intent(in):: iv
  type (gspace), intent(in):: gvec
  type (valence_wfns), intent(in):: vwfn
  type (conduction_wfns), intent(in):: cwfn
  type (polarizability), intent(inout):: pol
  type (kpoints), intent(inout):: kp, kpq
  integer, intent(in):: ispin, irk
  integer, intent(in):: rank_mtxel
  real(DP), intent(in):: kfact

  integer:: ic_loc, ic, ic_FE, offset_g, jsp, ig, iv_loc
  integer:: freq_idx
  integer:: ia, ib, ict

  integer, dimension(3):: Nfft
  real(DP):: scale, eden
  complex(DPC), dimension(:,:,:), allocatable:: fftbox1, fftbox2
  SCALAR, dimension(:), allocatable:: tmparray
  logical:: keep_transition

  PUSH_SUB(mtxel_cpu)

  call timing%start(timing%mtxel_fft)

  ! FHJ: notation:
  ! iv_loc -> local valence band index
  ! iv -> global valence band index (iv == 1 := lowest-energy band)
  ! ic_loc -> local conduction band index
  ! ic -> global conduction band (ic == nv+1 := lowest-energy unocc. band)
  ! ic_FE -> global conduction band, but starting at FE
  !          (ic_FE == 1 := first cond band). Rarely used.

  iv_loc = peinf%indexv(iv)

  if(pol%nfreq_group .gt. 1 .and. pol%gcomm .eq. 0) then
    freq_idx = 1
  else
    freq_idx = rank_mtxel+1
  endif

!--------------------------
! Use FFTs to calculate matrix elements

! Compute size of FFT box we need

  call setup_FFT_sizes(pol%FFTgrid, Nfft, scale)
! Allocate FFT boxes
  SAFE_ALLOCATE(fftbox2, (Nfft(1), Nfft(2), Nfft(3)))

! Put the data for valence band iv into FFT box 1 and do the FFT

  if (pol%os_opt_ffts /= 2) then
    SAFE_ALLOCATE(fftbox1, (Nfft(1), Nfft(2), Nfft(3)))
  endif

! JRD XXX This needs to be threaded
!$OMP PARALLEL DO collapse(2)
  do ic_loc = 1, peinf%ncownactual
    do ig = 1, pol%nmtx
      pol%gme(ig, ic_loc, iv_loc, ispin, irk, freq_idx) = ZERO
    enddo
  enddo

  SAFE_ALLOCATE(tmparray, (pol%nmtx))

  do jsp = ispin, ispin*kp%nspinor

    if (pol%os_opt_ffts /= 2) then
      call put_into_fftbox(vwfn%ngv, vwfn%zv(:,jsp), gvec%components, vwfn%isort, fftbox1, Nfft)
      call do_FFT(fftbox1, Nfft, 1)
      ! We need the complex conjugate of u_{vk+q)(r) for the cross correlation
      call conjg_fftbox(fftbox1, Nfft)
    endif

! Now we loop over the conduction states and get the matrix elements:
! 1. Get conduction wave function and put it into box 2, 
! 2. do FFT, get u_{ck}(r)
! 3. multiply by box1 contents, get F(r) = [u_{vk+q)(r)]^* u_{ck}(r)
! 4. do FFT again, and extract the resulting matrix elements and put the into pol
! We conjugate the final result since we really want < ck|e^{-i(q+G).r}|vk+q>
! but we have calculated < vk+q|e^{i(q+G).r}|ck>.

    do ic_loc = 1, peinf%ncownactual
      ic_FE = peinf%invindexc(ic_loc)
      ic = vwfn%nband+ic_FE
      offset_g = (ic_loc-1)*cwfn%ngc

      if (pol%os_opt_ffts == 2) then
        ! FHJ: optimization level 2 precomputed all the FFTs

        call timing%start(timing%fft_put)

!$OMP PARALLEL DO collapse(2) PRIVATE(ia, ib, ict)
        do ia = 1, Nfft(3)
        do ib = 1, Nfft(2)
          do ict = 1, Nfft(1)
            fftbox2(ict, ib, ia) = vwfn%wfn_fft(ict, ib, ia, iv_loc) * cwfn%wfn_fft(ict, ib, ia, ic_loc)
          enddo
        enddo
        enddo
!$OMP END PARALLEL DO

        call timing%stop(timing%fft_put)
     else
        if (pol%os_opt_ffts == 1) then

          call timing%start(timing%fft_put)

          ! FHJ: Optimization level 1 precalculated at least these cond. FFTs
!$OMP PARALLEL DO collapse(2) PRIVATE (ia, ib, ict)
          do ia = 1, Nfft(3)
          do ib = 1, Nfft(2)
            do ict = 1, Nfft(1)
              fftbox2(ict, ib, ia) = fftbox1(ict, ib, ia) * cwfn%wfn_fft(ict, ib, ia, ic_loc)
            enddo
          enddo
          enddo
!$OMP END PARALLEL DO

          call timing%stop(timing%fft_put)
        else
          call put_into_fftbox(cwfn%ngc, cwfn%zc(offset_g+1:,jsp), gvec%components, cwfn%isort, fftbox2, Nfft)
          call do_FFT(fftbox2, Nfft, 1)
          call multiply_fftboxes(fftbox1, fftbox2, Nfft)
        endif
      endif

      call do_FFT(fftbox2, Nfft, 1)
      call get_from_fftbox(pol%nmtx, tmparray, gvec%components, pol%isrtx, fftbox2, Nfft, scale)

      ! Get energy denominator (static), or transition energy (FF)
      call get_eden(gvec, vwfn, cwfn, kp, kpq, iv, ic, offset_g, iv_loc, &
                    ic_loc, ispin, irk, freq_idx, pol, eden)
      
      call update_gme(tmparray, kp, eden, kfact, iv, ic, jsp, ic_loc, &
                      iv_loc, ispin, irk, freq_idx, pol) 

      call timing%stop(timing%fft_mltply)
    end do  ! ic_loc
  end do  ! jsp

  SAFE_DEALLOCATE(tmparray)

! We are done, so deallocate FFT boxes
  if (pol%os_opt_ffts .ne. 2) then
    SAFE_DEALLOCATE(fftbox1)
  endif
  SAFE_DEALLOCATE(fftbox2)


! End FFT Case
!---------------------------
  call timing%stop(timing%mtxel_fft)

  POP_SUB(mtxel_cpu)
end subroutine mtxel_cpu

subroutine mtxel_openacc(iv, gvec, vwfn, cwfn, pol, ispin, irk, kp, kpq, rank_mtxel, kfact)
  integer, intent(in):: iv
  type (gspace), intent(in):: gvec
  type (valence_wfns), intent(in):: vwfn
  type (conduction_wfns), intent(in):: cwfn
  type (polarizability), intent(inout):: pol
  type (kpoints), intent(inout):: kp, kpq
  integer, intent(in):: ispin, irk
  integer, intent(in):: rank_mtxel
  real(DP), intent(in):: kfact
#ifdef OPENACC
  integer:: ic_loc, ic, ic_FE, offset_g, jsp, ig, iv_loc
  integer:: freq_idx
  integer:: ib, ict

  integer:: i_block, i_block_loc, block_size, i_band_start, i_band_end, n1_loc
  integer:: iii

  integer, dimension(3):: Nfft
  real(DP):: scale, eden

  PUSH_SUB(mtxel_openacc)

  call timing%start(timing%mtxel_fft)

  ! FHJ: notation:
  ! iv_loc -> local valence band index
  ! iv -> global valence band index (iv == 1 := lowest-energy band)
  ! ic_loc -> local conduction band index
  ! ic -> global conduction band (ic == nv+1 := lowest-energy unocc. band)
  ! ic_FE -> global conduction band, but starting at FE
  !          (ic_FE == 1 := first cond band). Rarely used.

  iv_loc = peinf%indexv(iv)

  if(pol%nfreq_group .gt. 1 .and. pol%gcomm .eq. 0) then
    freq_idx = 1
  else
    freq_idx = rank_mtxel+1
  endif

  if (pol%os_opt_ffts .ne. 0) then
    call die("os_opt_ffts not implemented for accelerated MTXEL", only_root_writes = .true.)
  end if

!--------------------------
! Use FFTs to calculate matrix elements

  ! Only for scale
  call setup_FFT_sizes(pol%FFTgrid, Nfft, scale)
  !$acc enter data copyin(Nfft)
  block_size = accel_mtxel_eps%mtxel_band_block_size

! Put the data for valence band iv into FFT box 1 and do the FFT

! JRD XXX This needs to be threaded
!$OMP PARALLEL DO collapse(2)
  do ic_loc = 1, peinf%ncownactual
    do ig = 1, pol%nmtx
      pol%gme(ig, ic_loc, iv_loc, ispin, irk, freq_idx) = ZERO
    enddo
  enddo

  if ( iv_loc == 1 ) then
    !$OMP parallel do collapse(2) default(shared) private(jsp, ic_loc, ic_FE, ic, offset_g)
    do jsp = ispin, ispin*kp%nspinor
      do ic_loc = 1, peinf%ncownactual
        ic_FE = peinf%invindexc(ic_loc)
        ic = vwfn%nband+ic_FE
        offset_g = (ic_loc-1)*cwfn%ngc
        accel_mtxel_eps%inner_bands_2d(1:peinf%ncownactual*cwfn%ngc, jsp) = cwfn%zc(1:peinf%ncownactual*cwfn%ngc, jsp)
      end do
    end do
    !$OMP end parallel do
    !$acc update device(accel_mtxel_eps%inner_bands_2d)
  end if

  accel_mtxel_eps%gvec_comp   = 0
  accel_mtxel_eps%outer_sort  = 0
  accel_mtxel_eps%inner_sort  = 0
  accel_mtxel_eps%result_sort = 0
  ! copy in arrays
  accel_mtxel_eps%gvec_comp = gvec%components
  call accel_update_to(accel_mtxel_eps%gvec_comp, mtxel_algo)
  accel_mtxel_eps%outer_sort = vwfn%isort
  call accel_update_to(accel_mtxel_eps%outer_sort, mtxel_algo)
  accel_mtxel_eps%inner_sort = cwfn%isort
  call accel_update_to(accel_mtxel_eps%inner_sort, mtxel_algo)
  accel_mtxel_eps%result_sort = pol%isrtx
  call accel_update_to(accel_mtxel_eps%result_sort, mtxel_algo)

  do jsp = ispin, ispin*kp%nspinor
    accel_mtxel_eps%outer_bands = ZERO
    accel_mtxel_eps%outer_bands(1:vwfn%ngv) = vwfn%zv(1:vwfn%ngv, jsp)
    !$acc update device(accel_mtxel_eps%outer_bands) async(1)
    
    ! FFT valence band
    call accel_zero_box(accel_mtxel_eps%outer_fftbox, Nfft, mtxel_algo, queue = 1)
    call accel_put_into_fftbox(vwfn%ngv, accel_mtxel_eps%outer_bands, accel_mtxel_eps%gvec_comp, &
                             accel_mtxel_eps%outer_sort, accel_mtxel_eps%outer_fftbox, &
                             Nfft, 1.0D+00, mtxel_algo, queue = 1)
    call accel_run_fft(accel_mtxel_eps%outer_fftbox, accel_mtxel_eps%outer_plan, 1, mtxel_algo)

    i_block_loc = 0
    do i_block = 1, peinf%ncownactual, block_size
      i_block_loc  = i_block_loc+1
      i_band_start = i_block 
      i_band_end   = MIN( (i_band_start+block_size-1), peinf%ncownactual )

      n1_loc = 0
      do ic_loc = i_band_start, i_band_end
        n1_loc = n1_loc+1
        ic_FE = peinf%invindexc(ic_loc)
        ic = vwfn%nband+ic_FE
        offset_g = (ic_loc-1)*cwfn%ngc

        call accel_zero_box(accel_mtxel_eps%inner_box(n1_loc)%box, Nfft, mtxel_algo, queue=(n1_loc+1))
        call accel_put_into_fftbox(cwfn%ngc, accel_mtxel_eps%inner_bands_2d(offset_g+1:offset_g+cwfn%ngc, jsp), &
                                   accel_mtxel_eps%gvec_comp, &
                                   accel_mtxel_eps%inner_sort, accel_mtxel_eps%inner_box(n1_loc)%box, Nfft, &
                                   1.0D+00, mtxel_algo, queue=(n1_loc+1))
        call accel_run_fft(accel_mtxel_eps%inner_box(n1_loc)%box, accel_mtxel_eps%inner_plan(n1_loc), 1, mtxel_algo)

        ! sync the outer stream only first cycle
        if (i_block_loc == 1 .and. n1_loc == 1) then
          !$acc wait(1)
        end if

        call accel_box_multiply(accel_mtxel_eps%inner_box(n1_loc)%box, &
                                accel_mtxel_eps%outer_fftbox, Nfft, 1, mtxel_algo, queue=(n1_loc+1))

        call accel_run_fft(accel_mtxel_eps%inner_box(n1_loc)%box, accel_mtxel_eps%inner_plan(n1_loc), 1, mtxel_algo)

        call accel_get_from_fftbox(pol%nmtx, accel_mtxel_eps%result_a(n1_loc)%vec, &
                                   accel_mtxel_eps%gvec_comp, accel_mtxel_eps%result_sort, &
                                   accel_mtxel_eps%inner_box(n1_loc)%box, Nfft, scale, mtxel_algo, queue=(n1_loc+1))

        !$acc update self(accel_mtxel_eps%result_a(n1_loc)%vec(1:pol%nmtx)) async(n1_loc+1)
      end do

      ! sync and update host
      n1_loc = 0
      do ic_loc = i_band_start, i_band_end
        n1_loc = n1_loc+1
        ic_FE = peinf%invindexc(ic_loc)
        ic = vwfn%nband+ic_FE
        offset_g = (ic_loc-1)*cwfn%ngc

        ! syncronize stream here
        !$acc wait(n1_loc+1)

        ! Get energy denominator (static), or transition energy (FF)
        call get_eden(gvec, vwfn, cwfn, kp, kpq, iv, ic, offset_g, iv_loc, &
                      ic_loc, ispin, irk, freq_idx, pol, eden)

        call update_gme(accel_mtxel_eps%result_a(n1_loc)%vec, kp, eden, kfact, iv, ic, jsp, &
                        ic_loc, iv_loc, ispin, irk, freq_idx, pol)
      end do

    end do  ! i_block
  end do  ! jsp
  !$acc exit data delete(Nfft)

! End FFT Case
!---------------------------
  call timing%stop(timing%mtxel_fft)
#else
  PUSH_SUB(mtxel_openacc)

  call die("OpenACC version of mtxel requested, but OpenACC not compiled&
           & into this executable", only_root_writes = .true.)
#endif

  POP_SUB(mtxel_openacc)
end subroutine mtxel_openacc

subroutine mtxel_accel_batched(iv, gvec, vwfn, cwfn, pol, ispin, irk, kp, kpq, rank_mtxel, kfact)
  integer, intent(in):: iv
  type (gspace), intent(in):: gvec
  type (valence_wfns), intent(in):: vwfn
  type (conduction_wfns), intent(in):: cwfn
  type (polarizability), intent(inout):: pol
  type (kpoints), intent(inout):: kp, kpq
  integer, intent(in):: ispin, irk
  integer, intent(in):: rank_mtxel
  real(DP), intent(in):: kfact

  integer:: ic_loc_temp, ic_FE, jsp, ig, iv_loc
  integer, allocatable:: ic(:), offset_g(:), ic_loc(:)
  SCALAR, dimension(:,:), allocatable:: tmparray
  integer:: i_batch, i_in_batch, n_ffts_in_batch
  integer:: freq_idx
  integer:: ia, ib, ict

  integer, dimension(3):: Nfft
  real(DP):: scale, eden
  logical:: keep_transition

  PUSH_SUB(mtxel_cpu)

  call timing%start(timing%mtxel_fft)

  if (pol%os_opt_ffts .ne. 0) then
    call die("os_opt_ffts not implemented for accelerated MTXEL", only_root_writes = .true.)
  end if

  ! FHJ: notation:
  ! iv_loc -> local valence band index
  ! iv -> global valence band index (iv == 1 := lowest-energy band)
  ! ic_loc -> local conduction band index
  ! ic -> global conduction band (ic == nv+1 := lowest-energy unocc. band)
  ! ic_FE -> global conduction band, but starting at FE
  !          (ic_FE == 1 := first cond band). Rarely used.

  iv_loc = peinf%indexv(iv)

  if(pol%nfreq_group .gt. 1 .and. pol%gcomm .eq. 0) then
    freq_idx = 1
  else
    freq_idx = rank_mtxel+1
  endif

!--------------------------
! Use FFTs to calculate matrix elements

! Compute size of FFT box we need

  call setup_FFT_sizes(pol%FFTgrid, Nfft, scale)

! Allocate batched inner loop quantities
  SAFE_ALLOCATE(ic_loc, (accel_mtxel_eps%n_ffts_per_batch))
  SAFE_ALLOCATE(ic, (accel_mtxel_eps%n_ffts_per_batch))
  SAFE_ALLOCATE(offset_g, (accel_mtxel_eps%n_ffts_per_batch))

! JRD XXX This needs to be threaded
!$OMP PARALLEL DO collapse(2)
  do ic_loc_temp = 1, peinf%ncownactual
    do ig = 1, pol%nmtx
      pol%gme(ig, ic_loc_temp, iv_loc, ispin, irk, freq_idx) = ZERO
    enddo
  enddo

  if ( iv_loc == 1 ) then
    !$OMP parallel do collapse(2) default(shared) private(jsp, ic_loc_temp, ic_FE, ic, offset_g)
    do jsp = ispin, ispin*kp%nspinor
      do ic_loc_temp = 1, peinf%ncownactual
        ic_FE = peinf%invindexc(ic_loc_temp)
        ic = vwfn%nband+ic_FE
        offset_g = (ic_loc_temp-1)*cwfn%ngc
        accel_mtxel_eps%inner_bands_2d(1:peinf%ncownactual*cwfn%ngc, jsp) = cwfn%zc(1:peinf%ncownactual*cwfn%ngc, jsp)
      end do
    end do
    !$OMP end parallel do
    call accel_update_to(accel_mtxel_eps%inner_bands_2d, mtxel_algo)
  end if

  accel_mtxel_eps%gvec_comp   = 0
  accel_mtxel_eps%outer_sort  = 0
  accel_mtxel_eps%inner_sort  = 0
  accel_mtxel_eps%result_sort = 0
  ! copy in arrays
  accel_mtxel_eps%gvec_comp = gvec%components
  call accel_update_to(accel_mtxel_eps%gvec_comp, mtxel_algo)
  accel_mtxel_eps%outer_sort = vwfn%isort
  call accel_update_to(accel_mtxel_eps%outer_sort, mtxel_algo)
  accel_mtxel_eps%inner_sort = cwfn%isort
  call accel_update_to(accel_mtxel_eps%inner_sort, mtxel_algo)
  accel_mtxel_eps%result_sort = pol%isrtx
  call accel_update_to(accel_mtxel_eps%result_sort, mtxel_algo)

  SAFE_ALLOCATE(tmparray, (pol%nmtx, accel_mtxel_eps%n_ffts_per_batch))

  call accel_enter_data_map_to(Nfft, mtxel_algo)
  call accel_enter_data_map_alloc(offset_g, mtxel_algo)
  call accel_enter_data_map_alloc(tmparray, mtxel_algo)
  do jsp = ispin, ispin*kp%nspinor
    accel_mtxel_eps%outer_bands = ZERO
    accel_mtxel_eps%outer_bands(1:vwfn%ngv) = vwfn%zv(1:vwfn%ngv, jsp)
    call accel_update_to(accel_mtxel_eps%outer_bands, mtxel_algo)

    call accel_zero_box(accel_mtxel_eps%outer_fftbox, Nfft, mtxel_algo)
    call accel_put_into_fftbox(vwfn%ngv, accel_mtxel_eps%outer_bands, accel_mtxel_eps%gvec_comp, &
                               accel_mtxel_eps%outer_sort, accel_mtxel_eps%outer_fftbox, Nfft, &
                               1.0D+00, mtxel_algo)
    call accel_run_fft(accel_mtxel_eps%outer_fftbox, accel_mtxel_eps%outer_plan, 1, mtxel_algo)

! Now we loop over the conduction states and get the matrix elements:
! 1. Get conduction wave function and put it into box 2, 
! 2. do FFT, get u_{ck}(r)
! 3. multiply by box1 contents, get F(r) = [u_{vk+q)(r)]^* u_{ck}(r)
! 4. do FFT again, and extract the resulting matrix elements and put the into pol
! We conjugate the final result since we really want < ck|e^{-i(q+G).r}|vk+q>
! but we have calculated < vk+q|e^{i(q+G).r}|ck>.

    ! Do a round-robin distribution of the conduction band into batches
    do i_batch = 1, (peinf%ncownactual-1) / accel_mtxel_eps%n_ffts_per_batch+1
      ! Create indexing arrays for batch and create batches of initial
      ! wavefunctions
      n_ffts_in_batch = 0
      call accel_zero_box_batched(accel_mtxel_eps%inner_fftbox_batched, Nfft, accel_mtxel_eps%n_ffts_per_batch, mtxel_algo)
      do i_in_batch = 1, accel_mtxel_eps%n_ffts_per_batch
        ic_loc(i_in_batch) = (i_batch-1) * accel_mtxel_eps%n_ffts_per_batch+i_in_batch
        if (ic_loc(i_in_batch) .gt. peinf%ncownactual) cycle

        n_ffts_in_batch = n_ffts_in_batch+1
        ic_FE = peinf%invindexc(ic_loc(i_in_batch))
        ic(i_in_batch) = vwfn%nband+ic_FE
        offset_g(i_in_batch) = (ic_loc(i_in_batch)- 1) * cwfn%ngc
      end do
      call accel_update_to(offset_g, mtxel_algo)
      call accel_put_into_fftbox_batched(cwfn%ngc, accel_mtxel_eps%inner_bands_2d(:,jsp), accel_mtxel_eps%gvec_comp, &
                                 accel_mtxel_eps%inner_sort, accel_mtxel_eps%inner_fftbox_batched, Nfft, &
                                 1.0D+00, offset_g, n_ffts_in_batch, mtxel_algo)

      ! Do the hard work of FFTs across batch
      call accel_run_fft_batched(accel_mtxel_eps%inner_fftbox_batched, accel_mtxel_eps%inner_plan_batched, 1, mtxel_algo)
      call accel_box_multiply_batched(accel_mtxel_eps%inner_fftbox_batched, accel_mtxel_eps%outer_fftbox, Nfft, 1, n_ffts_in_batch, mtxel_algo)
      call accel_run_fft_batched(accel_mtxel_eps%inner_fftbox_batched, accel_mtxel_eps%inner_plan_batched, 1, mtxel_algo)

      ! Extract finished FFTs and update results
      call accel_get_from_fftbox_batched(pol%nmtx, tmparray, accel_mtxel_eps%gvec_comp, accel_mtxel_eps%result_sort, &
                                   accel_mtxel_eps%inner_fftbox_batched, Nfft, scale, n_ffts_in_batch, mtxel_algo)
      call accel_update_from(tmparray, mtxel_algo)

      ! Get energy denominator (static), or transition energy (FF)
      do i_in_batch = 1, n_ffts_in_batch
        call get_eden(gvec, vwfn, cwfn, kp, kpq, iv, ic(i_in_batch), offset_g(i_in_batch), iv_loc, &
                    ic_loc(i_in_batch), ispin, irk, freq_idx, pol, eden)
        call update_gme(tmparray(:, i_in_batch), kp, eden, kfact, iv, ic(i_in_batch), jsp, ic_loc(i_in_batch), &
                        iv_loc, ispin, irk, freq_idx, pol)
        call timing%stop(timing%fft_mltply)
      end do  ! i_in_batch
    end do  ! i_batch
  end do  ! jsp
  call accel_exit_data_map_delete(Nfft, mtxel_algo)
  call accel_exit_data_map_delete(offset_g, mtxel_algo)
  call accel_exit_data_map_delete(tmparray, mtxel_algo)

! End FFT Case
!---------------------------
  call timing%stop(timing%mtxel_fft)

  SAFE_DEALLOCATE(tmparray)

  POP_SUB(mtxel_cpu)
end subroutine mtxel_accel_batched

subroutine mtxel_omp_target(iv, gvec, vwfn, cwfn, pol, ispin, irk, kp, kpq, rank_mtxel, kfact)
  integer, intent(in):: iv
  type (gspace), intent(in):: gvec
  type (valence_wfns), intent(in):: vwfn
  type (conduction_wfns), intent(in):: cwfn
  type (polarizability), intent(inout):: pol
  type (kpoints), intent(inout):: kp, kpq
  integer, intent(in):: ispin, irk
  integer, intent(in):: rank_mtxel
  real(DP), intent(in):: kfact
#ifdef OMP_TARGET
  integer:: ic_loc, ic, ic_FE, offset_g, jsp, ig, iv_loc
  integer:: freq_idx
  integer:: ib, ict

  integer, dimension(3):: Nfft
  real(DP):: scale, eden

  PUSH_SUB(mtxel_omp_target)

  call timing%start(timing%mtxel_fft)

  if (pol%os_opt_ffts .ne. 0) then
    call die("os_opt_ffts not implemented for accelerated MTXEL", only_root_writes = .true.)
  end if

  ! FHJ: notation:
  ! iv_loc -> local valence band index
  ! iv -> global valence band index (iv == 1 := lowest-energy band)
  ! ic_loc -> local conduction band index
  ! ic -> global conduction band (ic == nv+1 := lowest-energy unocc. band)
  ! ic_FE -> global conduction band, but starting at FE
  !          (ic_FE == 1 := first cond band). Rarely used.

  iv_loc = peinf%indexv(iv)

  if(pol%nfreq_group .gt. 1 .and. pol%gcomm .eq. 0) then
    freq_idx = 1
  else
    freq_idx = rank_mtxel+1
  endif

!--------------------------
! Use FFTs to calculate matrix elements

  ! Only for scale
  call setup_FFT_sizes(pol%FFTgrid, Nfft, scale)

! Put the data for valence band iv into FFT box 1 and do the FFT

! JRD XXX This needs to be threaded
!$OMP PARALLEL DO collapse(2)
  do ic_loc = 1, peinf%ncownactual
    do ig = 1, pol%nmtx
      pol%gme(ig, ic_loc, iv_loc, ispin, irk, freq_idx) = ZERO
    enddo
  enddo

  ! copy in arrays
  accel_mtxel_eps%gvec_comp = gvec%components
  call accel_update_to(accel_mtxel_eps%gvec_comp, mtxel_algo)
  accel_mtxel_eps%outer_sort = vwfn%isort
  call accel_update_to(accel_mtxel_eps%outer_sort, mtxel_algo)
  accel_mtxel_eps%inner_sort = cwfn%isort
  call accel_update_to(accel_mtxel_eps%inner_sort, mtxel_algo)
  accel_mtxel_eps%result_sort = pol%isrtx
  call accel_update_to(accel_mtxel_eps%result_sort, mtxel_algo)

  do jsp = ispin, ispin*kp%nspinor
    accel_mtxel_eps%outer_bands(1:vwfn%ngv) = vwfn%zv(1:vwfn%ngv, jsp)
    call accel_update_to(accel_mtxel_eps%outer_bands, 1, vwfn%ngv, mtxel_algo)
    accel_mtxel_eps%inner_bands(1:peinf%ncownactual*cwfn%ngc) = cwfn%zc(1:peinf%ncownactual*cwfn%ngc, jsp)
    call accel_update_to(accel_mtxel_eps%inner_bands, 1, &
                       peinf%ncownactual*cwfn%ngc, mtxel_algo)

    call accel_zero_box(accel_mtxel_eps%outer_fftbox, Nfft, mtxel_algo, queue = 1)
    call accel_put_into_fftbox(vwfn%ngv, accel_mtxel_eps%outer_bands, accel_mtxel_eps%gvec_comp, &
                               accel_mtxel_eps%outer_sort, accel_mtxel_eps%outer_fftbox, &
                               Nfft, 1.0D+00, mtxel_algo, queue = 1)
    call accel_run_fft(accel_mtxel_eps%outer_fftbox, accel_mtxel_eps%outer_plan, 1, mtxel_algo)

! Now we loop over the conduction states and get the matrix elements:
! 1. Get conduction wave function and put it into box 2, 
! 2. do FFT, get u_{ck}(r)
! 3. multiply by box1 contents, get F(r) = [u_{vk+q)(r)]^* u_{ck}(r)
! 4. do FFT again, and extract the resulting matrix elements and put the into pol
! We conjugate the final result since we really want < ck|e^{-i(q+G).r}|vk+q>
! but we have calculated < vk+q|e^{i(q+G).r}|ck>.

    do ic_loc = 1, peinf%ncownactual
      ic_FE = peinf%invindexc(ic_loc)
      ic = vwfn%nband+ic_FE
      offset_g = (ic_loc-1)*cwfn%ngc

      call accel_zero_box(accel_mtxel_eps%inner_fftbox, Nfft, mtxel_algo, queue = 1)
      call accel_put_into_fftbox(cwfn%ngc, accel_mtxel_eps%inner_bands(offset_g+1:), accel_mtxel_eps%gvec_comp, &
                                 accel_mtxel_eps%inner_sort, accel_mtxel_eps%inner_fftbox, Nfft, &
                                 1.0D+00, mtxel_algo, queue = 1)
      call accel_run_fft(accel_mtxel_eps%inner_fftbox, accel_mtxel_eps%outer_plan, 1, mtxel_algo)
      call accel_box_multiply(accel_mtxel_eps%inner_fftbox, &
                              accel_mtxel_eps%outer_fftbox, Nfft, 1, mtxel_algo, queue = 1)

      call accel_run_fft(accel_mtxel_eps%inner_fftbox, accel_mtxel_eps%outer_plan, 1, mtxel_algo)
      call accel_get_from_fftbox(pol%nmtx, accel_mtxel_eps%result_array, &
                               accel_mtxel_eps%gvec_comp, accel_mtxel_eps%result_sort, &
                               accel_mtxel_eps%inner_fftbox, Nfft, scale, mtxel_algo, queue = 1)
      call accel_update_from(accel_mtxel_eps%result_array, 1, pol%nmtx, mtxel_algo)

      ! Get energy denominator (static), or transition energy (FF)
      call get_eden(gvec, vwfn, cwfn, kp, kpq, iv, ic, offset_g, iv_loc, &
                    ic_loc, ispin, irk, freq_idx, pol, eden)

      call update_gme(accel_mtxel_eps%result_array, kp, eden, kfact, iv, ic, jsp, &
                      ic_loc, iv_loc, ispin, irk, freq_idx, pol)

      call timing%stop(timing%fft_mltply)
    end do  ! ic_loc
  end do  ! jsp

! End FFT Case
!---------------------------
  call timing%stop(timing%mtxel_fft)
#else
  PUSH_SUB(mtxel_omp_target)

  call die("OpenMP Target version of mtxel requested, but OpenMP Target&
           & not compiled into this executable", only_root_writes = .true.)
#endif

  POP_SUB(mtxel_omp_target)
end subroutine mtxel_omp_target

subroutine update_gme(tmparray, kp, eden, kfact, iv, ic, jsp, ic_loc, &
                      iv_loc, ispin, irk, freq_idx, pol)
  implicit none

  SCALAR, dimension(:), intent(in):: tmparray
  type (kpoints), intent(in):: kp
  real(DP), intent(in):: eden
  real(DP), intent(in):: kfact
  integer, intent(in):: iv
  integer, intent(in):: ic
  integer, intent(in):: jsp
  integer, intent(in):: ic_loc
  integer, intent(in):: iv_loc
  integer, intent(in):: ispin
  integer, intent(in):: irk
  integer, intent(in):: freq_idx
  type (polarizability), intent(inout):: pol

  logical:: keep_transition
  integer:: ia
  

  keep_transition = .true.

!#BEGIN_INTERNAL_ONLY
  select case (pol%intraband_flag)
    case (1)  ! FHJ: keep only intraband transitions
      keep_transition = abs(tmparray(1))**2 >= pol%intraband_overlap_min
    case (2)  ! FHJ: keep only interband transitions
      keep_transition = abs(tmparray(1))**2 < pol%intraband_overlap_min
    case (3)  ! FHJ: only include transitions inside the protection window
      keep_transition = all([iv, ic]>=pol%protection_window(1)).and.&
        all([iv, ic]<=pol%protection_window(2))
    case (4)  ! FHJ: block transitions inside the protection window
      keep_transition = .not.(all([iv, ic]>=pol%protection_window(1)).and.&
        all([iv, ic]<=pol%protection_window(2)))
    case (5)  ! DYQ: block transitions where the occ. state is outside & unocc. is inside the window
      keep_transition = .not. (iv < pol%protection_window(1).and.&
        ic <= pol%protection_window(2))
  end select
  keep_transition = keep_transition .and. (ic > pol%num_cond_bands_ignore)
!#END_INTERNAL_ONLY

  call timing%start(timing%fft_mltply)

  if (keep_transition) then
!$OMP PARALLEL DO
    do ia = 1, pol%nmtx
      pol%gme(ia, ic_loc, iv_loc, ispin, irk, freq_idx) = &
        pol%gme(ia, ic_loc, iv_loc, ispin, irk, freq_idx) + MYCONJG(tmparray(ia))*kfact
    end do
!$OMP END PARALLEL DO
  end if

  if (kp%nspinor .eq. 1 .or. jsp .eq. 2) then
    if (pol%freq_dep .eq. 0) then
!$OMP PARALLEL DO
      do ia = 1, pol%nmtx
        pol%gme(ia, ic_loc, iv_loc, ispin, irk, 1) = &
          pol%gme(ia, ic_loc, iv_loc, ispin, irk, 1) * &
          sqrt(-1D0*eden)
      end do
!$OMP END PARALLEL DO
    end if
  end if

end subroutine update_gme

!> FHJ: Get energy denominator (static), or transition energy (FF)
subroutine get_eden(gvec, vwfn, cwfn, kp, kpq, iv, ic, offset_g, iv_loc, &
                    ic_loc, ispin, irk, freq_idx, pol, eden)
  implicit none

  type (gspace), intent(in):: gvec
  type (valence_wfns), intent(in):: vwfn
  type (conduction_wfns), intent(in):: cwfn
  type (kpoints), intent(in):: kp
  type (kpoints), intent(in):: kpq
  integer, intent(in):: iv
  integer, intent(in):: ic
  integer, intent(in):: offset_g
  integer, intent(in):: iv_loc
  integer, intent(in):: ic_loc
  integer, intent(in):: ispin
  integer, intent(in):: irk
  integer, intent(in):: freq_idx
  type (polarizability), intent(inout):: pol
  real(DP), intent(out):: eden

  integer:: ig
  type(cvpair_info):: lin_edenTemp
  real(DP):: eval, econd, occ_v, occ_c, occ_diff
  real(DP):: vk(2), vkq(2)
  real(DP)  :: ewin_min, ewin_max  ! WK


  PUSH_SUB(mtxel.get_eden)

  call timing%start(timing%mtxel_denom)

  ! FHJ: See convention for conduction/valence band indices above.
  eval = vwfn%ev(iv, ispin) 
  econd = cwfn%ec(ic, ispin)
  eden = 0d0

  ! WK: Energy window
  ewin_min = pol%ewin_min/RYD  !pol%ewin_min is in eV so converting now
  ewin_max = pol%ewin_max/RYD  !pol%ewin_max is in eV so converting now

  if( allocated(pol%occq)) then
    occ_v = pol%occq(iv+vwfn%ncore_excl, vwfn%idx_kp, ispin, 1) 
  else
    occ_v = pol%occ(iv+vwfn%ncore_excl, vwfn%idx_kp, ispin, 1)
  endif
  occ_c = pol%occ(ic+vwfn%ncore_excl, cwfn%idx_kp, ispin, 1)

  occ_diff = occ_v-occ_c

  ! FHJ: Note that eden means different things depending on pol%freq_dep
  ! static: eden := 1/(transition energy).
  ! FF: eden := (transition energy). (I know...)
  ! In the static case, we lump sqrt(eden) into the matrix elements gme.
  ! In the FF case, we have to put it by hand for each frequency we evaluate chi.
  if (pol%freq_dep == 0) then

    if(eval-econd < TOL_Degeneracy .and. occ_diff > TOL_Zero) then
      ! avoid dividing by zero or making eden > 0
      eden = occ_diff / (eval-econd)
    else
      eden = 0d0  ! in this case, occ_diff = 0 too
    endif

    !> WK check energy if eval, econd in (ewin_min, ewin_max): 
    if ((pol%energy_window) .and. (pol%energy_window_mode == 1) .and. &
        (eval > ewin_min) .and. (eval < ewin_max) .and. & 
        (econd > ewin_min) .and. (econd < ewin_max)) then
      eden = 0d0  ! WK: in this case we put zero in eden 
    endif
    !> WK check energy
    
    !> WK check energy if eval or econd out of (ewin_min, ewin_max): 
    if ((pol%energy_window) .and. (pol%energy_window_mode == 2) .and. &
        ((eval <= ewin_min) .or. (eval >= ewin_max) .or. & 
         (econd <= ewin_min) .or. (econd >= ewin_max))) then
      eden = 0d0  ! WK:in this case we put zero in eden 
    endif
    !> WK check energy
    

!#BEGIN_INTERNAL_ONLY
    ! LZT: do linearized energy denominator
    if(eval-econd < TOL_Degeneracy .and. abs(eval-econd)<pol%lin_denominator)then
      !calculate velocities from k.p perturbation theory
      !valence bands velocity
      vk(1:2) = kpq%rk(1:2, vwfn%idx_kp)
      do ig = 1, vwfn%ngv
        vk = vk+gvec%components(1:2, vwfn%isort(ig)) * abs(vwfn%zv(ig, ispin))**2
      enddo

      !conduction bands velocity
      vkq(1:2) = kp%rk(1:2, cwfn%idx_kp)
      do ig = 1, cwfn%ngc
        vkq = vkq+gvec%components(1:2, cwfn%isort(ig)) * abs(cwfn%zc(offset_g+ig, ispin))**2
      enddo

      !calculate integral
      eden = integrate_mbz(kp, kpq%rk(1:3, vwfn%idx_kp), &
        kpq%el(iv, vwfn%idx_kp, ispin), kp%el(ic, cwfn%idx_kp, ispin), &
        vk, vkq, 0d0, pol%efermi/ryd)

    endif
    !end LZT linearized denominator
!#END_INTERNAL_ONLY

  elseif (pol%freq_dep == 2 .or. pol%freq_dep == 3) then

    ! FHJ: In chi_summation, we explicitly neglect transitions if eden < TOL_ZERO.
    ! That`s the way we keep track of forbidden transitions.
    if(eval-econd < TOL_Degeneracy .and. occ_diff > TOL_Zero) then
!      eden = (eval-econd) / occ_diff
      eden = (eval-econd)
    else
      eden = 0.0d0
    endif

!#BEGIN_INTERNAL_ONLY
    ! LZT: store energies and velocities for linearized energy denominator
    if(pol%lin_denominator > TOL_Zero) then
      lin_edenTemp%ec = econd
      lin_edenTemp%ev = eval
      lin_edenTemp%vltc = eval-econd < TOL_Degeneracy

      !calculate velocities from k.p perturbation theory
      !valence bands velocity
      vk(1:2) = kpq%rk(1:2, vwfn%idx_kp)
      do ig = 1, vwfn%ngv
        vk = vk+gvec%components(1:2, vwfn%isort(ig)) * abs(vwfn%zv(ig, ispin))**2
      enddo
      lin_edenTemp%vv = vk
      lin_edenTemp%idx_kp = vwfn%idx_kp

      !conduction bands velocity
      vkq(1:2) = kp%rk(1:2, cwfn%idx_kp)
      do ig = 1, cwfn%ngc
        vkq = vkq+gvec%components(1:2, cwfn%isort(ig)) * abs(cwfn%zc(offset_g+ig, ispin))**2
      enddo
      lin_edenTemp%vc = vkq
    endif
    !finished storing energies and velocities for linearized energy denom
!#END_INTERNAL_ONLY

    pol%edenDyn(iv_loc, ic_loc, ispin, irk, freq_idx) = eden
    pol%occ_diff(iv_loc, ic_loc, ispin, irk, freq_idx) = occ_diff
!#BEGIN_INTERNAL_ONLY
    if (pol%lin_denominator > TOL_Zero) then
      pol%lin_edenDyn(iv_loc, ic_loc, ispin, irk, freq_idx) = lin_edenTemp
    endif
!#END_INTERNAL_ONLY

  endif  ! pol%freq_dep

  call timing%stop(timing%mtxel_denom)

  POP_SUB(mtxel.get_eden)

end subroutine get_eden

end module mtxel_m
