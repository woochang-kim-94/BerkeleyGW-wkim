!===============================================================================
!
! Routines:
!
! (1) genwf_mpi()                             Last Modified 6/09/2010 (gsm)
!
! Generates valence and conduction wavefunctions at k-point rkq = rkn - rq
! from wavefunctions at k-point rkn read from structures distributed in memory.
!
!===============================================================================

#include "f_defs.h"

module genwf_mpi_m

  use global_m
  use gmap_m
  use input_utils_m
  use misc_m
  use sort_m
  use susymmetries_m
  implicit none

  private

  public :: &
    genwf_mpi, &
    time_reversal, &
    shift_phonq_g0

contains

subroutine genwf_mpi(rkq,syms,gvec,crys,kp,sig,wfnkq,wfnkqmpi)
!#BEGIN_INTERNAL_ONLY
!#END_INTERNAL_ONLY
  real(DP), intent(in) :: rkq(3)
  type (symmetry), intent(in) :: syms
  type (gspace), intent(in) :: gvec
  type (crystal), intent(in) :: crys
  type (kpoints), intent(in) :: kp
  type (siginfo), intent(in) :: sig
  type (wfnkqstates), intent(inout) :: wfnkq !< pointers are allocated outside
  type (wfnkqmpiinfo), intent(in) :: wfnkqmpi
!  logical, intent(in) :: lrcalc !< whether this is a linear-response phonon calc.

  integer :: nkpti
  integer :: nbandi
  integer, allocatable :: isrtk(:)
  SCALAR, allocatable :: zin(:,:)
#ifdef CPLX
  SCALAR, allocatable :: zkqtemp(:,:)
  complex(DP) :: umtrx(2,2)
#endif

  integer :: ii,iin,ik,irk,ikn1,it,iband,ispin, igk
  integer :: naddv
  integer :: nkpt2,kg0(3)
  integer, allocatable :: isorti(:),ind(:)
  real(DP) :: del,qk(3)
  real(DP), allocatable :: xnorm(:,:)
  real(DP), allocatable :: ekin(:)
  SCALAR, allocatable :: ph(:)
  logical :: found

!---------------------------------------
! find rotation (R), translation (g0), and reduced k-point (rk)
! such that krq = rkn - rq = R(rk) + g0
! it = rotation, kg0 = translation, irk = reduced k-point

  PUSH_SUB(genwf_mpi)

  found = .false.
  irk_loop: do irk=1,kp%nrk
    it_loop: do it=1,syms%ntran
      do ii=1,3
        qk(ii) = DOT_PRODUCT(dble(syms%mtrx(ii,:,it)),kp%rk(:,irk))
        del = rkq(ii) - qk(ii)
        if (del .ge. 0.0d0) kg0(ii) = del + TOL_Small
        if (del .lt. 0.0d0) kg0(ii) = del - TOL_Small
        if (abs(del-kg0(ii)) .gt. TOL_Small) cycle it_loop
      enddo
      found = .true.
      exit irk_loop
    enddo it_loop
  enddo irk_loop
  if(.not. found) call die('genwf: rkq mismatch')
  ikn1 = irk 
  wfnkq%kmq_pt = ikn1 ! OAH index of kpoint in original kpoint list that corresponds to k-q

!---------------------------------------
! Read in (rk-q) wavefunction

  SAFE_ALLOCATE(ekin, (gvec%ng))
  SAFE_ALLOCATE(isrtk, (gvec%ng))
  SAFE_ALLOCATE(isorti, (gvec%ng))
  SAFE_ALLOCATE(ind, (gvec%ng))
  SAFE_ALLOCATE(ph, (gvec%ng))
  nkpti = wfnkqmpi%nkptotal(ikn1)
  nbandi = sig%ntband
  isrtk(1:nkpti) = wfnkqmpi%isort(1:nkpti,ikn1)
  do ik = 1, sig%nspin
    wfnkq%ekq(1:nbandi,ik) = wfnkqmpi%el(1:nbandi,ik,ikn1)
  enddo
  qk(1:3) = wfnkqmpi%qk(1:3,ikn1)
  wfnkq%nkpt=nkpti
  if (any(abs(kp%rk(1:3,ikn1)-qk(1:3)) .gt. TOL_Small)) then
    call die('genwf: kp mismatch')
  endif

!---------------------------------------
! Compute inverse to array isort
! isrtk orders |rk + G|^2
! isrtkq  orders |rkn - rq + G|^2
! with krq = rkn - rq = R(rk) + g0

  isorti(:)=0
  do ii=1,nkpti
    isorti(isrtk(ii))=ii
  enddo
  SAFE_DEALLOCATE(isrtk)

!---------------------------------------
! compute kinetic energy for rkq = rk-q

  call kinetic_energies(gvec, crys%bdot, ekin, qvec = rkq)

!---------------------------------------
! sort array ekin to ascending order
! store indices in array isort

  call sortrx(gvec%ng, ekin, wfnkq%isrtkq, gvec = gvec%components)

!---------------------------------------
! map indices for rk(ikn1) to those for rkq

  call gmap(gvec,syms,nkpti,it,kg0, &
    wfnkq%isrtkq,isorti,ind,ph,.true.)

!---------------------------------------
! loop over wfnkq%zkq wavefunctions
! map zin wfns onto wfnkq%zkq

  SAFE_ALLOCATE(xnorm, (sig%ntband,sig%nspin))
  SAFE_ALLOCATE(zin, (wfnkq%nkpt,sig%nspin*kp%nspinor))
  SAFE_ALLOCATE(wfnkq%zkq, (peinf%ntband_node*wfnkq%nkpt,sig%nspin*kp%nspinor))
  do iin=1,peinf%ntband_node
    iband=wfnkqmpi%band_index(iin,ikn1)
    nkpt2=wfnkqmpi%nkptotal(ikn1)
    do ispin = 1, sig%nspin*kp%nspinor
      zin(1:nkpt2,ispin)=wfnkqmpi%cg(1:nkpt2,iin,ispin,irk)
    enddo
    if (iband.ne.peinf%indext(iin)) then
!      if(peinf%inode.eq.0) then
        write(0,*) 'iband=',iband,' <> ',peinf%indext(iin)
!      endif
      call die('genwf: iband error')
    endif
    
    naddv = (iin-1)*wfnkq%nkpt

!---------------------------------------
! loop over components of zv

    do ispin=1,sig%nspin*kp%nspinor
      do igk=1,wfnkq%nkpt
        ! ZL: if ind(igk) is zero (which happens for dWFN 
        !     the new gvec list for a kponit is different from 
        !     that in WFN, and also that generated from gsphere
        !     so we must check this now), we set the coeff to ZERO
        if (ind(igk) .gt. 0) then
          wfnkq%zkq(naddv+igk,ispin) = zin(ind(igk),ispin)*ph(igk)
        else
          wfnkq%zkq(naddv+igk,ispin) = ZERO
        endif
      enddo
      ! ZL: original lines
!      wfnkq%zkq(1+naddv:wfnkq%nkpt+naddv,ispin) = &
!        zin(ind(1:wfnkq%nkpt),ispin)*ph(1:wfnkq%nkpt)
    enddo

! In spinor case, we must rotate spinors according to spinor rotation matrix umtrx

#ifdef CPLX
      if (kp%nspinor.eq.2) then
        SAFE_ALLOCATE(zkqtemp, (wfnkq%nkpt,kp%nspin*kp%nspinor))
        call susymmetries(crys%bvec(:,:),syms%mtrx(:,:,it),umtrx,it)
        zkqtemp = matmul(wfnkq%zkq(1+naddv:wfnkq%nkpt+naddv,:),umtrx)
        wfnkq%zkq(1+naddv:wfnkq%nkpt+naddv,:) = zkqtemp
        SAFE_DEALLOCATE(zkqtemp)
      endif
#endif

! BAB: We check if the norm differs appreciably from unity.
! There is no longer a need to further normalize the vector.

    if (peinf%check_norms) then
      do ispin=1,sig%nspin
        call compute_norm(xnorm(iin,ispin),ispin,wfnkq%nkpt,kp%nspinor,wfnkq%zkq(:,:))
        if(abs(xnorm(iin,ispin) - 1d0) > TOL_Small) then
          write(0,*) 'Bad norm: irk=',irk,'ispin=',ispin,'norm=',xnorm(iin,ispin)
          call die('genwf_mpi: Bad norm')
        endif
      enddo
    endif
  enddo
  SAFE_DEALLOCATE(ekin)
  SAFE_DEALLOCATE(isorti)
  SAFE_DEALLOCATE(ind)
  SAFE_DEALLOCATE(ph)
  
!---------------------------------------
! end of loop over wavefunctions

  POP_SUB(genwf_mpi)
  
  return
end subroutine genwf_mpi


subroutine time_reversal(wfnkq, gvec)
  type (wfnkqstates), intent(inout) :: wfnkq
  type (gspace), intent(in) :: gvec

  integer :: mgv(3)
  integer :: ikg
  integer :: ispin

  ! ZL: no PUSH/POP because called too frequently

  ! ZL: TR only supports nspin*nspinor=1
  ispin = 1

  do ikg=1,wfnkq%nkpt
    mgv(:) = -gvec%components(:,wfnkq%isrtkq(ikg))
    call findvector(wfnkq%isrtkq(ikg), mgv, gvec)
  enddo
  if (any(wfnkq%isrtkq(1:wfnkq%nkpt)==0)) then
    call die("Master G-list is not large enough for time-reversal.", only_root_writes=.false.)
  endif
  ! ZL: for all bands distributed on each process
  wfnkq%zkq(:,ispin) = MYCONJG(wfnkq%zkq(:,ispin))

  return
end subroutine time_reversal


subroutine shift_phonq_g0(g0, wfnk, gvec)
  integer, intent(in) :: g0(3)
  type (wfnkstates), intent(inout) :: wfnk
  type (gspace), intent(in) :: gvec

  integer :: gmg0(3) ! G minus G0, the transform is G=>G-G0
  integer :: ikg

  ! ZL: no PUSH/POP for simplicity

  ! wfn(k_out) = wfn(k_in)
  ! k_out = k_in + g0
  ! G_out = G_in - g0

  do ikg=1,wfnk%nkpt
    gmg0(:) = gvec%components(:,wfnk%isrtk(ikg)) - g0(:)  ! ZL: MUST be consistent with input.f90, in which g0 is defined
    call findvector(wfnk%isrtk(ikg), gmg0, gvec)
  enddo
  if (any(wfnk%isrtk(1:wfnk%nkpt)==0)) then
    call die("Master G-list is not large enough for shift_phonq_g0.", only_root_writes=.false.)
  endif

  return
end subroutine shift_phonq_g0

end module genwf_mpi_m
