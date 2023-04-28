!===========================================================================
!
! Module mtxelmultiply_m
!
! (1) mtxelMultiply()   Refactored from main.f90 By (PWD)     
!                                           Last Modified 10/21/2010 (PWD)
!
!     Combine the <c,k|e^(-i(q+G).r)|v,k+q><v,k+q|e^(i(q+G`).r)|c,k>
!                 -------------------------------------------------
!                                  energy denominator
!     Input:
!       pol%eden(band,spin) = 1/(e_val-e_cond) = energy denominators
!       pol%gme(band,g-vector,spin) = plane wave matrix elements
!       pol%isrtx   orders the |G(i)|^2   i=1,pol%nmtx
!       vwfn%isort  orders |qk+g|^2    (in vwfn type)
!
!       energies are apparently assumed in Rydbergs.
!
!===========================================================================

#include "f_defs.h"

module mtxelmultiply_m

  use algos_epsilon_m
  use accel_linalg_m
  use accel_memory_m
  use blas_m
  use global_m
  use scalapack_m
  use timing_m, only: timing => epsilon_timing

  implicit none

  private

  public :: &
    mtxelMultiply

contains

  subroutine mtxelmultiply(scal,ntot,nrk,nst,fact,vwfn,gmetempr,gmetempc,chilocal, &
    polgme,pol,indt,pht,ipe,ispin,nrow_max,ncol_max)
    type (scalapack), intent(in) :: scal
    integer, intent(in) :: ntot
    integer, intent(in) :: nrk
    integer, intent(in) :: nst(:) !< (nrk)
    real(DP), intent(in) :: fact
    type (valence_wfns), intent(in) :: vwfn
    SCALAR, dimension(:,:), intent(inout) :: gmetempr(:,:),gmetempc(:,:)
    SCALAR, dimension(:,:), intent(inout) :: chilocal(:,:)
    SCALAR, intent(in) :: polgme(:,:,:,:,:,:)
    type (polarizability), intent(inout) :: pol
    integer, dimension(:,:,:), intent(in) :: indt
    SCALAR, dimension(:,:,:), intent(in) :: pht
    integer, intent(in) :: ipe
    integer, intent(in) :: ispin
    integer, intent(in) :: nrow_max,ncol_max

    SCALAR :: negfact
    
    PUSH_SUB(mtxelMultiply)
    
    negfact = -1D0*fact
    
    call timing%start(timing%chi_sum_prep)
    call chi_sum_prep(scal,nrk,nst,gmetempr,gmetempc,chilocal, &
                      polgme,pol,indt,pht,ipe,ispin)
    call timing%stop(timing%chi_sum_prep)
    
    call timing%start(timing%chi_sum_gemm)
    
    if (ntot.ne.0) then
      if (.not. pol%accel_chisum_fulloffload) then
       if ( chi_summation_algo /= CPU_ALGO ) then
         select case(chi_summation_algo)
           case (OPENACC_ALGO)
#ifdef OPENACC
             !$acc update device(gmetempr, gmetempc, chilocal) async(1)
#else
             call die_algos("OpenACC")
#endif
           case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
             !$omp target update to(gmetempr, gmetempc, chilocal)
#else
             call die_algos("OpenMP Target")
#endif
           case default
             call die("Invald algorithm for chi_summation.", only_root_writes = .true.)
         end select
       end if
      end if
      call accel_xgemm('n', 't', &
                     nrow_max, ncol_max, ntot, &
                     negfact, &
                     gmetempr, nrow_max, &
                     gmetempc, ncol_max, &
                     ZERO, &
                     chilocal, nrow_max, &
                     chi_summation_algo)

      if (chi_summation_algo /= CPU_ALGO) then
        select case(chi_summation_algo)
          case (OPENACC_ALGO)
#ifdef OPENACC
            !$acc update self(chilocal) async(1)
#else
      call die_algos("OpenACC")
#endif
          case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
            !$omp target update from(chilocal)
#else
            call die_algos("OpenMP Target")
#endif
          case default
            call die("Invald algorithm for chi_summation.", only_root_writes = .true.)
        end select
      end if

    end if

    call timing%stop(timing%chi_sum_gemm)

    POP_SUB(mtxelMultiply)
    
    return
    
  end subroutine mtxelMultiply
  
  subroutine chi_sum_prep(scal,nrk,nst,gmetempr,gmetempc,chilocal, &
    polgme,pol,indt,pht,ipe,ispin)
    type (scalapack), intent(in) :: scal
    integer, intent(in) :: nrk
    integer, intent(in) :: nst(:) !< (nrk)
    SCALAR, dimension(:,:), intent(inout) :: gmetempr(:,:),gmetempc(:,:)
    SCALAR, dimension(:,:), intent(inout) :: chilocal(:,:)
    SCALAR, intent(in) :: polgme(:,:,:,:,:,:)
    type (polarizability), intent(inout) :: pol
    integer, dimension(:,:,:), intent(in) :: indt
    SCALAR, dimension(:,:,:), intent(in) :: pht
    integer, intent(in) :: ipe
    integer, intent(in) :: ispin

    integer, allocatable :: tmprowindex(:),tmpcolindex(:)
    SCALAR, allocatable :: tmprowph(:),tmpcolph(:)
    integer :: irk, iv, j, it, icurr, itot, mytot

    SAFE_ALLOCATE(tmprowindex,(scal%nprd(ipe)))
    SAFE_ALLOCATE(tmpcolindex,(scal%npcd(ipe)))
    SAFE_ALLOCATE(tmprowph,(scal%nprd(ipe)))
    SAFE_ALLOCATE(tmpcolph,(scal%npcd(ipe)))

    if (chi_summation_algo /= CPU_ALGO .and. pol%accel_chisum_fulloffload) then
#ifdef OPENACC
      !$acc enter data create(tmprowindex,tmpcolindex,tmprowph,tmpcolph) async(1)
#else
      call die_algos("OpenACC")
#endif
    end if

    itot = 0

    do irk = 1, nrk
      do it = 1, nst(irk)

   if ( chi_summation_algo == OPENACC_ALGO .and. pol%accel_chisum_fulloffload) then
#ifdef OPENACC
        !$acc parallel private(icurr) present(scal%nprd, indt, pht, scal%imyrowd, scal%imyrowd, tmprowindex, tmprowph) &
        !$acc& async(1) vector_length(512)
        !$acc loop gang vector
        do icurr=1,scal%nprd(ipe)
          tmprowindex(icurr) = indt(scal%imyrowd(icurr,ipe),it,irk)
          tmprowph(icurr) = pht(scal%imyrowd(icurr,ipe),it,irk)
        enddo 
        !$acc end parallel

        !$acc parallel private(icurr) present(scal%npcd, indt, pht, scal%imycold, scal%imycold, tmpcolindex, tmpcolph) &
        !$acc& async(1) vector_length(512)
        !$acc loop gang vector
        do icurr=1,scal%npcd(ipe)
          tmpcolindex(icurr) = indt(scal%imycold(icurr,ipe),it,irk)
          tmpcolph(icurr) = pht(scal%imycold(icurr,ipe),it,irk)
        enddo
        !$acc end parallel

        !$acc parallel private(iv, j, mytot, icurr) present(peinf, scal%nprd, scal%npcd, gmetempr, gmetempc, polgme, &
        !$acc& tmprowindex, tmpcolindex, tmprowph, tmpcolph) &
        !$acc& async(1) vector_length(512)
        !$acc loop gang vector collapse(2)
        do iv = 1,peinf%nvownactual
          do j = 1, peinf%ncownactual

            mytot = itot + (iv-1)*peinf%ncownactual + j

            gmetempr(:,mytot) = ZERO
            gmetempc(:,mytot) = ZERO

! JRD XXX the index here probably generates gather instructions
! May want to also consider streaming stores

            do icurr=1,scal%nprd(ipe)
              gmetempr(icurr,mytot)=polgme( &
                tmprowindex(icurr),j,iv, &
                ispin,irk,pol%nfreq_group)* &
                tmprowph(icurr)
            enddo

            do icurr=1,scal%npcd(ipe)
              gmetempc(icurr,mytot)= &
                MYCONJG(polgme(tmpcolindex(icurr),j,iv,ispin,irk,pol%nfreq_group)*tmpcolph(icurr))
            enddo

          enddo ! j
        enddo ! iv
        !$acc end parallel

        !$acc wait(1)
#else
      call die_algos("OpenACC")
#endif
   else ! .not. chi_summation_algo == OPENACC_ALGO .or. ...
        do icurr=1,scal%nprd(ipe)
          tmprowindex(icurr) = indt(scal%imyrowd(icurr,ipe),it,irk)
          tmprowph(icurr) = pht(scal%imyrowd(icurr,ipe),it,irk)
        enddo
        do icurr=1,scal%npcd(ipe)
          tmpcolindex(icurr) = indt(scal%imycold(icurr,ipe),it,irk)
          tmpcolph(icurr) = pht(scal%imycold(icurr,ipe),it,irk)
        enddo

!$OMP PARALLEL DO collapse(2) private (mytot,iv,j,icurr)
        do iv = 1,peinf%nvownactual
          do j = 1, peinf%ncownactual
            mytot = itot + (iv-1)*peinf%ncownactual + j

            gmetempr(:,mytot) = ZERO
            gmetempc(:,mytot) = ZERO

! JRD XXX the index here probably generates gather instructions
! May want to also consider streaming stores
            do icurr=1,scal%nprd(ipe)
              gmetempr(icurr,mytot)=polgme( &
                tmprowindex(icurr),j,iv, &
                ispin,irk,pol%nfreq_group)* &
                tmprowph(icurr)
            enddo

            do icurr=1,scal%npcd(ipe)
              gmetempc(icurr,mytot)= &
                MYCONJG(polgme(tmpcolindex(icurr),j,iv,ispin,irk,pol%nfreq_group)*tmpcolph(icurr))
            enddo ! icurr
          enddo ! j
        enddo ! iv
!$OMP END PARALLEL DO

   end if ! chi_summation_algo

        itot = itot + peinf%nvownactual*peinf%ncownactual

      enddo ! it
    enddo ! irk

    if (chi_summation_algo /= CPU_ALGO .and. pol%accel_chisum_fulloffload) then
#ifdef OPENACC
      !$acc exit data delete(tmprowindex,tmpcolindex,tmprowph,tmpcolph) async(1)
#else
      call die_algos("OpenACC")
#endif
    end if

    SAFE_DEALLOCATE(tmprowindex)
    SAFE_DEALLOCATE(tmpcolindex)    
    SAFE_DEALLOCATE(tmprowph)
    SAFE_DEALLOCATE(tmpcolph)    

  end subroutine

end module mtxelmultiply_m
