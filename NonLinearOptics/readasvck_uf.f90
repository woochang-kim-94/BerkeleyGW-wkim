!=============================================================================
!
! Routines:
!
! (1) readasvck_uf()    Orginally By JRD        Last Modified 6/6/2008 (JRD)
!
!     This routine reads in eigenvectors and calculates (by calling oscstrength)
!     of the ultrafast oscillator strength of a system.
!
!=============================================================================

#include "f_defs.h"

module readasvck_uf_m

  use global_m
  use oscstrength_m
  use evecs_m
  implicit none

  private
  public :: readasvck_uf

contains

subroutine readasvck_uf(win,nmat,neig,s1,xct,oscarray,enarray,nk,ikmax,nkpt,k)
  type (windowinfo), intent(in) :: win
  integer, intent(in) :: nmat
  integer, intent(inout) :: neig
  SCALAR, intent(in) :: s1(:) !< (nmat)
  type (xctinfo), intent(inout) :: xct
  real(DP), intent(out) :: oscarray(:,:) !< (nmat,win%nstates)
  real(DP), intent(out) :: enarray(:,:) !< (nmat,win%nstates)
  integer, intent(out) :: nk
  integer, intent(out) :: ikmax(:) !< (nmat)
  integer, intent(in) :: nkpt
  real(DP), intent(out) :: k(:,:) !< (3,nkpt)

  integer :: ns, nc, nv, i, m, unit, c, v, s, j
  integer :: ncounter, nblock, nmatold
  integer :: iflagfound,jj
  integer :: ik, cmax(nmat), vmax(nmat)
  real(DP) :: e, w, e0(win%nstates), wmax,osc
  real(DP), allocatable :: kosc(:)
  SCALAR, allocatable :: A(:,:,:,:)
  SCALAR :: osct
  type(evecs_t) :: evecs

  PUSH_SUB(readasvck_uf)

  unit = 10
  iflagfound=0

  call evecs%read_file(master_only=.true.)
  
  ns = evecs%ns
  nv = evecs%nv
  nc = evecs%nc
  nk = evecs%nk
  nmatold = evecs%nmat
  nblock = nv*nc*ns
  k(:,:) = evecs%kpts(:,:)

!------------- I Processor is 0 --------------------------------------------

  if (peinf%inode.eq.0) then

    call open_file(11,file='osc.dat',form='formatted',status='replace')
    call open_file(12,file='kosc.dat',form='formatted',status='replace')
    call open_file(13,file='kosc3D.dat',form='formatted',status='replace')
    call open_file(14,file='maxes.dat',form='formatted',status='replace')


! JRD: Debug

    write(6,*) ' '
    write(6,*) "------------------------------------------"
    write(6,*) "Reading eigenvectors"
    write(6,*) " "
    write(6,*) 'ns,nv,nc,nk = ',ns,nv,nc,nk
    write(6,*) xct%nspin,xct%nvb_fi,xct%ncb_fi,xct%nkpt_fi
    write(6,*) 'nmat, nmatold:',nmat, nmatold
    write(6,*) 'S0 Energy = ', win%evalue
    write(6,*) ' '

    SAFE_ALLOCATE(kosc, (nk))
    SAFE_ALLOCATE(A, (ns,nv,nc,nk))
    ncounter=0

    do i=1,nmatold
      
      do jj=1,win%nstates
        e = evecs%evals(jj)
        if (abs(e-win%estates(jj)) < TOL_Small) then
          e0(jj)=e
          if (win%istates(jj) .ne. i) then
            write(0,*) i,win%istates(jj)
            call die('Mismatch of States')
          endif
          write(6,*) 'State Found',jj,e0(jj)
        endif
      enddo
      
    enddo
    
    ncounter=1
    
    do i=1,nmatold
      
      e = evecs%evals(i)
      
      if (e < win%emax)  then
        m = 0
        ncounter=ncounter+1
        write(6,*) ' '
        write(6,*) ' In energy window: ', ncounter
        A(:,:,:,:) = evecs%Avc(:,:,:,:,i)
        
        wmax = 0.0d0
        
        do c=1,nc
          do v=1,nv
            w = 0.0d0
            do ik=1,nk
              do s=1,ns
                w = w + abs(A(s,v,c,ik))**2
                if (abs(A(s,v,c,ik))**2 > wmax) then
                  wmax = abs(A(s,v,c,ik))**2
                  ikmax(ncounter-1) = ik
                  cmax(ncounter-1) = c
                  vmax(ncounter-1) = v
                endif
              enddo
            enddo
          enddo
        enddo
        write(14,*) ncounter-1,e-e0,ikmax(ncounter-1), &
          vmax(ncounter-1),cmax(ncounter-1)
        
        write(6,*) ' '
        write(6,*) ' Eigenvalue #', ncounter
        
        do jj = 1, win%nstates
          if (win%istates(jj) .ne. i) then
            call oscstrength(nmat,nk,ns,xct,evecs%Avc(:,:,:,:,jj),A,s1,osct)
            osc = abs( osct )**2
            write(6,*) 'jj, osc, e-e0(jj) = ', jj, osc, e-e0(jj)
            oscarray(ncounter-1,jj)=osc
            enarray(ncounter-1,jj)=e-e0(jj)
          else
            enarray(ncounter-1,jj)=1d0
            oscarray(ncounter-1,jj)=0d0
          endif
        enddo
        write(11,'(2i8,3e16.8)') i,ncounter-1,k(3,ikmax(ncounter-1)), &
          e-e0(1), oscarray(ncounter-1,1)
      endif
    enddo
    
    do i = 1, nk
      kosc(i)=0D0
    end do
    
    do i = 1, ncounter-1
      kosc(ikmax(i))=kosc(ikmax(i))+oscarray(i,1)
      write(13,*) vmax(i),cmax(i),k(3,ikmax(i)),enarray(i,1), &
        oscarray(i,1)
    end do
    
    do i = 1, nk
      write(12,*) (k(j,i),j=1,3), kosc(i)
    end do
    
    write(6,*) 'Number of EVals in Window', ncounter-1
    neig=ncounter-1
    
    !call close_file(unit)
    call close_file(11)
    call close_file(12)
    call close_file(13)
    call close_file(14)
  end if
  
#ifdef MPI
  call MPI_BCAST(neig,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  
  call evecs%free()
  if (peinf%inode.eq.0) then
    SAFE_DEALLOCATE(A)
  endif
  
  POP_SUB(readasvck_uf)
  
  return
end subroutine readasvck_uf

end module readasvck_uf_m
