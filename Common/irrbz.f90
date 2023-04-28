#include "f_defs.h"

module irrbz_m

  use global_m
  use misc_m

  implicit none

  private

  public :: irrbz

contains

  subroutine irrbz(syms,nfk,fk,nrq,neq,indrq,rq,nq,qq,itnrq,kg0,nfix)
    type (symmetry), intent(in) :: syms
    integer, intent(in) :: nfk
    real(DP), intent(in) :: fk(:,:) !< (3,nfk)
    integer, intent(out) :: nrq
    integer, optional, intent(out) :: neq(:) !< nrq
    integer, optional, intent(out) :: indrq(:) !< nrq
    real(DP), optional, intent(out) :: rq(:,:) !< (3, nrq)
    integer, optional, intent(in) :: nq
    real(DP), optional, intent(in) :: qq(:,:)  !< (3, nq)
    integer, optional, intent(out) :: itnrq(:) !< nrq
    integer, optional, intent(out) :: kg0(:,:) !< (3, nrq)
    !> don`t try to fold back the first nfix points (defaults to 1)
    integer, optional, intent(in) :: nfix
  
    integer :: i,j,irq,ik,it,iflag,iq,kg(3),nq_read,line_skip,iqsave,iostat_c,nfix_
    integer, allocatable :: indrk(:)
    real(DP) :: qk(3)
    real(DP), allocatable :: rq_read(:,:)
    logical :: use_ntran ! true for Sigma, false for Epsilon
    logical :: irrbz_read
  
!
!  i   nq,qq              q-points for which eps**(-1)(q) is known
!  i   nfk,fk             k-points in full BZ
!  o   nrq,rq             q-points in irr BZ with respect to subgroup
!  o   neq                number of points equivalent to rq
!  o   indrq,itnrq,kg0    index, transformation, Umklapp required for
!                           rq = r(q) + kg0
!
    PUSH_SUB(irrbz)

    nfix_ = 1
    if (present(nfix)) nfix_ = nfix
  
    ! all or nothing
    if(present(rq) .and. present(nq) .and. present(qq) .and. present(itnrq) .and. present(kg0) .neqv. &
       present(rq) .or.  present(nq) .or.  present(qq) .or.  present(itnrq) .or.  present(kg0)) then
       call die("irrbz internal error: bad parameters")
    endif
    use_ntran = present(qq)
    if(use_ntran .and. .not. (present(neq) .and. present(indrq))) then
      call die("irrbz internal error: bad parameters 2")
    endif
  
    if(size(fk, 1) /= 3) then
      write(0,*) 'size(fk, 1) = ', size(fk, 1)
      call die("irrbz internal error: fk must have first dimension = 3")
    endif
    if(size(fk, 2) /= nfk) then
      write(0,*) 'nfk = ', nfk, 'size(fk, 2) = ', size(fk, 2)
      call die("irrbz internal error: fk must have second dimension = nfk")
    endif

    irrbz_read = .false.
    if (.not. use_ntran) then
      call open_file(unit=21,file='irrbz.dat',status='old',iostat=iostat_c)
      if (iostat_c==0) irrbz_read=.true.
    endif
    if (irrbz_read) then
      write(6,'(3x,a)') 'Reading the Irreducible Brillouin Zone from irrbz.dat'
      read(21,*) nq_read
      SAFE_ALLOCATE(rq_read,(3,nq_read))
      read(21,*) (rq_read(1:3,iq), iq = 1,nq_read)
  
      !Compare and find which rq_read corresponds to the syms%rq 
      do iq=1,nq_read
        if (all(abs(rq_read(1:3,iq)-syms%rq(1:3)) .lt. TOL_Small)) then
           iqsave=iq
           exit
        endif
        write(0,'(3x,a)') 'rq not found in irrbz.dat',syms%rq(1:3)
        call die("rq not found in irrbz.dat")
      enddo
      SAFE_DEALLOCATE(rq_read)
     
      !Now just go to the line number from where we can read all the information
      line_skip = 3*(iqsave-1)
      do i=1,line_skip
        read(21,*)
      enddo
      read(21,*) nrq
      if (present(neq)) then
        read(21,*) (neq(i), i=1,nrq)
      else
        read(21,*) 
      endif
      if (present(indrq)) then
        read(21,*) (indrq(i), i=1,nrq)
      else
        read(21,*) 
      endif
      call close_file(21)
    
    else
  
      SAFE_ALLOCATE(indrk, (nfk))
  
!   initialize number of points in irr. BZ
!
      nrq = 0
!
!   loop over k-points in full zone
!
      ik_loop: do ik=1,nfk
        if (ik > nfix_) then
!
!       loop over transformation matrices
!
          do it=1,syms%ntranq
            qk(:) = matmul(dble(syms%mtrx(:,:,syms%indsub(it))),fk(:,ik))
            call k_range(qk,kg,TOL_Small)
  
            if(.not. present(neq)) cycle
!
!         compare to other k-points in the irr. BZ with respect to qvec
!
            do irq=1,nrq
              if (all(abs(fk(1:3,indrk(irq))-qk(1:3)) .lt. TOL_Small)) then
                neq(irq) = neq(irq) + 1
                cycle ik_loop
              endif
            enddo
          enddo ! loop over ntranq tranformation
  
        endif
  
        nrq = nrq + 1
        indrk(nrq) = ik
        if(use_ntran) rq(1:3,nrq) = fk(1:3,ik)
        if(present(neq)) neq(nrq) = 1
  
        if(use_ntran) then
!
!     find qq to which rq is equivalent
!
          iflag = 0
          iq_loop: do iq=1,nq
            do it=1,syms%ntran
              qk(:) = matmul(dble(syms%mtrx(:,:,it)),qq(:,iq))
              call k_range(qk,kg,TOL_Small)
            
              if (all(abs(rq(1:3,nrq)-qk(1:3)) .lt. TOL_Small)) then
                iflag = 1
                indrq(nrq) = iq
                itnrq(nrq) = it
                kg0(1:3,nrq) = kg(1:3)
! there may be more than one matching sym op, but better be only one q-vector!
!              exit iq_loop
              endif
            enddo
          enddo iq_loop
  !
  !     end of loop over q
  
          if (iflag .eq. 0) then
            write(0,'(a,3f12.6)') 'rq = ', rq(1:3,nrq)
            call die('irrbz: q/rq mismatch')
          endif
  
        else
          if(present(indrq)) indrq(nrq) = ik
! we do not need these, but this is what they would be
!        kg0(1:3,nrq) = kg(1:3)
!        itnrq(nrq) = it
        endif ! use_ntran
  
      enddo ik_loop !end loop over full BZ
    endif !whether read or generated
  
    if (peinf%inode .eq. 0 .and. use_ntran) then
      ! FHJ: We only print neq, everything else depends on which is the first fq
      ! related to rq.
      write(6,'(1x,a,i0)') 'Number of q-points in the irreducible BZ(k) (nrq): ', nrq
      if (peinf%verb_medium) then
        write(6,'(/6x,a,5x,a)') 'q-point rq (irr. BZ)', '#eq/fBZ'
        write(6,'(1x,29("-"),1x,7("-"))')
        write(6,'(3(1x,f9.6),1x,i7)') (rq(1:3,j), neq(j), j=1,nrq)
      endif
    endif
  
    SAFE_DEALLOCATE(indrk)
    POP_SUB(irrbz)
  
    return
  end subroutine irrbz

end module irrbz_m
