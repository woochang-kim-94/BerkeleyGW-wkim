!-------------------------------------------------------------------------------
!
!   Routines for kgrid.x
!   Generates a list of k-points with a small q-shift for espresso input file
!   inspired by espresso-3.2.3/pwtools/kpoints.f
!   written by G. Samsonidze (March 2008)
!
!-------------------------------------------------------------------------------

#include "f_defs.h"

module kgrid_routines_m

  use global_m
  use misc_m

  implicit none
  
  private

  public :: fold, reduce, sort, dqunfold, dqsubgrp, dqfold, dqsort, write_sym, &
    write_espresso, write_octopus

contains

!-------------------------------------------------------------------------------

  subroutine fold(nkf,kpoint,kweight,kfold,ksymmetry,nr,r,trs,eps)
    integer, intent(in) :: nkf
    real(DP), intent(in) :: kpoint(:,:) !< (3,nkf)
    integer, intent(inout) :: kweight(:) !< (nkf)
    integer, intent(inout) :: kfold(:) !< (nkf)
    integer, intent(inout) :: ksymmetry(:) !< (nkf)
    integer, intent(in) :: nr
    integer, intent(in) :: r(:,:,:) !< (3,3,48)
    logical, intent(in) :: trs !< time-reversal symmetry
    real(DP), intent(in) :: eps
    
    integer :: ik1,ik2,ikf,ir,gpt(3)
    real(DP) :: k(3)
    
    PUSH_SUB(fold)
    
    do ik1 = 1, nkf
! DAS: a cleaner way to write the loop, removing warning about infinite loop
      ir_loop: do ir = 1, nr
        k(1:3) = matmul(dble(r(1:3, 1:3, ir)), kpoint(1:3, ik1))
        call k_range(k,gpt,eps)
! DAS: a cleaner way to write the loop, removing warning about infinite loop
        do ik2 = 1, ik1 - 1
          if(all(abs(k(1:3)-kpoint(1:3,ik2)-dble(nint(k(1:3)-kpoint(1:3,ik2)))).lt.eps) .or. &
            (trs .and. all(abs(k(1:3)+kpoint(1:3,ik2)-dble(nint(k(1:3)+kpoint(1:3,ik2)))).lt.eps))) then
            kweight(ik1)=0
            ikf=ik2
            do while(kweight(ikf).eq.0)
              ikf=kfold(ikf)
            enddo
            kweight(ikf)=kweight(ikf)+1
            kfold(ik1)=ikf
            ksymmetry(ik1)=ir
            exit ir_loop
          endif
        enddo
      enddo ir_loop
    enddo
    
    POP_SUB(fold)
    return
    
  end subroutine fold
  
!-------------------------------------------------------------------------------

  subroutine reduce(nkf,kweight,kfold)
    integer, intent(in) :: nkf
    integer, intent(inout) :: kweight(:),kfold(:) !< (nkf)
    
    integer :: ik1,ik2
    integer, allocatable :: ws(:),fs(:)
    
    PUSH_SUB(reduce)
    
    SAFE_ALLOCATE(ws, (nkf))
    SAFE_ALLOCATE(fs, (nkf))
    do ik1=1,nkf
      if (kweight(ik1).gt.0) then
        do ik2=1,nkf
          if (kfold(ik2).eq.ik1.or.ik2.eq.ik1) then
            ws(ik2)=0
            fs(ik2)=ik1
          endif
        enddo
        ws(ik1)=kweight(ik1)
        fs(ik1)=0
      endif
    enddo
    kweight(1:nkf)=ws(1:nkf)
    kfold(1:nkf)=fs(1:nkf)
    SAFE_DEALLOCATE(ws)
    SAFE_DEALLOCATE(fs)
    
    POP_SUB(reduce)
    return
    
  end subroutine reduce

!-------------------------------------------------------------------------------

  subroutine sort(nkf,nkr,kpoint,kweight,kindex,eps)
    integer, intent(in) :: nkf,nkr
    real(DP), intent(in) :: kpoint(:,:) !< (3,nkf)
    integer, intent(inout) :: kweight(:) !< (nkf)
    integer, intent(out) :: kindex(:) !< (nkf)
    real(DP), intent(in) :: eps
    
    integer :: ik1,ik2,ki1,ki2
    
    PUSH_SUB(sort)
    
    ik1=0
    do ik2=1,nkf
      if (kweight(ik2).gt.0) then
        ik1=ik1+1
        kindex(ik1)=ik2
      endif
    enddo
    do ik1=1,nkr
      do ik2=1,nkr-1
        if ((kpoint(1,kindex(ik2+1)).lt. &
          kpoint(1,kindex(ik2))-eps).or. &
          (abs(kpoint(1,kindex(ik2+1))- &
          kpoint(1,kindex(ik2))).lt.eps.and. &
          kpoint(2,kindex(ik2+1)).lt. &
          kpoint(2,kindex(ik2))-eps).or. &
          (abs(kpoint(1,kindex(ik2+1))- &
          kpoint(1,kindex(ik2))).lt.eps.and. &
          abs(kpoint(2,kindex(ik2+1))- &
          kpoint(2,kindex(ik2))).lt.eps.and. &
          kpoint(3,kindex(ik2+1)).lt. &
          kpoint(3,kindex(ik2))-eps)) then
          ki1=kindex(ik2)
          ki2=kindex(ik2+1)
          kindex(ik2)=ki2
          kindex(ik2+1)=ki1
        endif
      enddo
    enddo
    
    POP_SUB(sort)
    return
    
  end subroutine sort
  
!-------------------------------------------------------------------------------

  subroutine dqunfold(nkr,kr,nkf,kf,kfw,nr,r,trs,eps)
    integer, intent(in) :: nkr !< number of kpoints before unfolding
    real(DP), intent(in) :: kr(:,:) !< (3,nkr) coords of kpoints before unfolding
    integer, intent(out) :: nkf
    real(DP), intent(out) :: kf(:,:) !< (3,nkf)
    real(DP), intent(out) :: kfw(:) !< (nkf)
    integer, intent(in) :: nr !< number of symmetries with which to unfold
    integer, intent(in) :: r(:,:,:) !< (3,3,48) symmetry matrices
    logical, intent(in) :: trs !< time-reversal symmetry
    real(DP), intent(in) :: eps
    
    integer :: ikr,ikf,ir,gpt(3)
    real(DP) :: k(3),minusk(3)
    logical :: add_kpt, add_kpt_trs
    
    PUSH_SUB(dqunfold)
    
    nkf=0
    do ikr=1,nkr
      ! see whether it is unique, with and without trs
      do ir=1,nr
        k(1:3) = matmul(dble(r(1:3, 1:3, ir)), kr(1:3, ikr))
        minusk(1:3) = -k(1:3)
        call k_range(k,gpt,eps)
        call k_range(minusk,gpt,eps)

        add_kpt = .true.
        add_kpt_trs = trs .and. .not. all(abs(k(1:3)-minusk(1:3)).lt.eps)
        if (nkf.gt.0) then
          do ikf=1,nkf
            if (all(abs(kf(1:3,ikf)-k(1:3)).lt.eps)) add_kpt = .false.
            if (all(abs(kf(1:3,ikf)-minusk(1:3)).lt.eps)) add_kpt_trs=.false.
          enddo
        endif
        if(add_kpt) then
          nkf=nkf+1
          kf(1:3,nkf)=k(1:3)
          kfw(nkf)=1.0d0
        endif
        if(add_kpt_trs) then
          nkf=nkf+1
          kf(1:3,nkf)=minusk(1:3)
          kfw(nkf)=1.0d0
        endif
      enddo
    enddo

    POP_SUB(dqunfold)
    return
    
  end subroutine dqunfold

!-------------------------------------------------------------------------------
  
  subroutine dqsubgrp(dq,nr,r,nrq,rq,syms,trs,eps)
    real(DP), intent(in) :: dq(:) !< (3)
    integer, intent(in) :: nr
    integer, intent(in) :: r(:,:,:) !< (3,3,48)
    integer, intent(out) :: nrq
    integer, intent(out) :: rq(:,:,:) !< (3,3,48)
    logical, intent(out) :: syms(:) !< (48)
    logical, intent(out) :: trs !< time-reversal symmetry
    real(DP), intent(in) :: eps
    
    integer :: ir,irq
    real(DP) :: k(3)
    integer :: rq_inv(3,3)
    
    PUSH_SUB(dqsubgrp)
    
    irq=0
    do ir=1,nr
      syms(ir)=.false.
      k(1:3) = matmul(dble(r(1:3, 1:3, ir)), dq(1:3))
      if (all(abs(k(1:3)-dq(1:3)).lt.eps)) then
        irq=irq+1
        rq(1:3,1:3,irq)=r(1:3,1:3,ir)
        syms(ir)=.true.
      endif
    enddo
    nrq=irq
    do irq=1,nrq
      call invert_matrix_int(rq(1:3, 1:3, irq), rq_inv)
      rq(1:3, 1:3, irq) = rq_inv(1:3, 1:3)
    enddo
    
    ! is dq invariant under time-reversal symmetry?
    trs = all(abs(dq(:)) < eps .or. abs(abs(dq(:)) - 0.5d0) < eps)

    POP_SUB(dqsubgrp)
    return
  
  end subroutine dqsubgrp

!-------------------------------------------------------------------------------

  subroutine dqfold(nkf,kf,kfw,nkq,kq,kqw,nrq,rq,trs,eps)
    integer, intent(in) :: nkf
    real(DP), intent(in) :: kf(:,:) !< (3,nkf)
    real(DP), intent(in) :: kfw(:) !< (nkf)
    integer, intent(inout) :: nkq
    real(DP), intent(out) :: kq(:,:) !< (3,nkq)
    real(DP), intent(out) :: kqw(:) !< (nkq)
    integer, intent(in) :: nrq
    integer, intent(in) :: rq(:,:,:) !< (3,3,48)
    logical, intent(in) :: trs !< time-reversal symmetry
    real(DP), intent(in) :: eps
    
    integer :: ikf,ikq,irq,gpt(3)
    real(DP) :: k(3), minusk(3)
  
    PUSH_SUB(dqfold)
    
    nkq=0
    ikf_loop: do ikf=1,nkf
      if (ikf.gt.1) then
        do irq=1,nrq
          k(1:3) = matmul(rq(1:3, 1:3, irq), kf(1:3, ikf))
          minusk(1:3) = -k(1:3)
          call k_range(k,gpt,eps)
          call k_range(minusk,gpt,eps)
          do ikq=1,nkq
            if (all(abs(kq(1:3,ikq)-k(1:3)).lt.eps) .or. &
              (trs .and. all(abs(kq(1:3,ikq)-minusk(1:3)).lt.eps))) then
              kqw(ikq)=kqw(ikq)+kfw(ikf)
              cycle ikf_loop
            endif
          enddo
        enddo
      endif
      nkq=nkq+1
      kq(1:3, nkq) = kf(1:3, ikf)
      kqw(nkq)=kfw(ikf)
    enddo ikf_loop
  
    POP_SUB(dqfold)
    return
  
  end subroutine dqfold

!-------------------------------------------------------------------------------

  subroutine dqsort(nkq,kq,kqw,eps)
    integer, intent(in) :: nkq
    real(DP), intent(inout) :: kq(:,:) !< (3,nkq)
    real(DP), intent(inout) :: kqw(:) !< (nkq)
    real(DP), intent(in) :: eps
    
    integer :: ik1,ik2,imin
    real(DP) :: kmin(3),kwmin
    
    PUSH_SUB(dqsort)
    
    do ik1=1,nkq-1
      imin=0
      kmin(1:3)=INF
      do ik2=ik1,nkq
        if ((kq(1,ik2).lt.kmin(1)-eps).or. &
          (abs(kq(1,ik2)-kmin(1)).lt.eps.and. &
          kq(2,ik2).lt.kmin(2)-eps).or. &
          (abs(kq(1,ik2)-kmin(1)).lt.eps.and. &
          abs(kq(2,ik2)-kmin(2)).lt.eps.and. &
          kq(3,ik2).lt.kmin(3)-eps)) then
          imin=ik2
          kmin(1:3)=kq(1:3,ik2)
          kwmin=kqw(ik2)
        endif
      enddo
      if (imin.ne.ik1) then
        kq(1:3,imin)=kq(1:3,ik1)
        kqw(imin)=kqw(ik1)
        kq(1:3,ik1)=kmin(1:3)
        kqw(ik1)=kwmin
      endif
    enddo
    
    POP_SUB(dqsort)
    return
    
  end subroutine dqsort

!-------------------------------------------------------------------------------

  subroutine write_sym(name, symbol, nr, sym_ops, frac_trans, spacegroup)
    character(len=*),  intent(in) :: name
    character(len=*),  intent(in) :: symbol
    integer,           intent(in) :: nr !< number of symmetry operations
    integer,           intent(in) :: sym_ops(:,:,:) !< (3,3,nr)
    real(DP),          intent(in) :: frac_trans(:,:) !< (3,nr)
    integer, optional, intent(in) :: spacegroup

    integer :: ir, ii, jj

    PUSH_SUB(write_sym)

    write(23,'(a)')
    write(23,'(2a)') "symmetries ", TRUNC(name)
    write(23,'(i3)') nr
    write(23,'(a,i3,a,a)') 'Space group ', spacegroup, ', symbol ', symbol
    do ir=1,nr
      write(23,'("r",i2.2," =",3(2x,3i4),2x,3f12.6)') &
        ir,((sym_ops(jj,ii,ir),jj=1,3),ii=1,3),(frac_trans(ii,ir),ii=1,3)
    enddo
    if (nr.eq.0) write(23,'("none")')

    POP_SUB(write_sym)
    return

  end subroutine write_sym

  subroutine write_espresso(iunit, nk, coords, weights)
    integer,  intent(in) :: iunit
    integer,  intent(in) :: nk
    real(DP), intent(in) :: coords(:,:) !< (3, nk)
    real(DP), intent(in) :: weights(:) !< (nk)

    integer :: ik

    PUSH_SUB(write_espresso)

    write(iunit,'(a)') "K_POINTS crystal"
    write(iunit,'(i5)') nk

    do ik = 1, nk
      write(iunit,'(3f13.9,f6.1)') coords(1:3, ik), weights(ik)
    enddo

    POP_SUB(write_espresso)
  end subroutine write_espresso

  subroutine write_octopus(iunit, nk, coords, weights, kgrid, shift)
    integer,  intent(in) :: iunit
    integer,  intent(in) :: nk
    real(DP), intent(in) :: coords(:,:) !< (3, nk)
    real(DP), intent(in) :: weights(:) !< (nk)
    integer,  intent(in) :: kgrid(:) !< (3)
    real(DP), intent(in) :: shift(:) !< (3); shift = dk + dq*nk.

    integer :: ik

    PUSH_SUB(write_octopus)

    write(iunit,'(a)') "%KGrid"
    write(iunit,'(i10,a,i10,a,i10)')       kgrid(1), ' | ', kgrid(2), ' | ', kgrid(3)
    write(iunit,'(f10.8,a,f10.8,a,f10.8)') shift(1), ' | ', shift(2), ' | ', shift(3)
    write(iunit,'(a)') "%"

    write(iunit,'(a)') "%KPointsReduced"

    do ik = 1, nk
      write(iunit,'(f6.1,3(" | ",f13.9))') weights(ik), coords(1:3, ik)
    enddo

    write(iunit,'(a)') "%"

    POP_SUB(write_octopus)
  end subroutine write_octopus

end module kgrid_routines_m
