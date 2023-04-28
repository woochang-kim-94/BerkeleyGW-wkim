!==============================================================================================
!
! Utilities:
!
! (1) ploteps           Originally By JRD               Last Modified 6/30/2008 (JRD)
!
!==============================================================================================

#include "f_defs.h"

program ploteps

  use global_m
  implicit none

  integer, allocatable :: isrtq(:),oldx(:),oldy(:),oldz(:)
  integer, allocatable :: isrtold(:)
  real(DP) :: qvec(3),qk(3),gmax_in,garbaged
  integer :: ii, jj, nold, nge, nmtx, nrq0

  SCALAR, allocatable :: eps(:,:)

  character :: ajname*6, adate*11

  write(*,*) 'Welcome to Epsilon Plotting'

  call open_file(unit=10,file='eps0mat',form='unformatted',status='old')
!      call open_file(unit=11,file='epsmat',form='unformatted',status='old')
  call open_file(unit=42,file='eps_surface.out',status='replace')
  call open_file(unit=43,file='eps_diagonal.out',status='replace')

!---------------------------
! Now read eps0mat

  read(10)
  read(10) ii
  if (ii.ne.0) then
    call die('Full frequency dependence not supported')
  endif
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10) nold
  read(10) nge
  rewind(10)
!      deallocate(oldx)
  SAFE_ALLOCATE(oldx, (nold))
!      deallocate(oldy)
  SAFE_ALLOCATE(oldy, (nold))
!      deallocate(oldz)
  SAFE_ALLOCATE(oldz, (nold))
!      deallocate(ekold)
!      allocate(ekold(nge))
!      deallocate(isrtold)
  SAFE_ALLOCATE(isrtold, (nge))
!      deallocate(isrtq)
!      allocate(isrtq(nge))
!      deallocate(isrtq)
  SAFE_ALLOCATE(isrtq, (nge))

  read(10) ajname,adate
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10) gmax_in
  read(10) nrq0,(qvec(ii),ii=1,3)
  read(10) nold,(oldx(ii),oldy(ii),oldz(ii),ii=1,nold)
  
! JRD: This should be changed later

  if(nrq0.gt.1) then
    call die("There is more than one q-point in epsmat", only_root_writes = .true.)
  endif

! Read q->0 dielectric matrix

  read(10) nge,nmtx,(isrtold(ii),jj,ii=1,nge)
  write(*,*) 'nge, nmtx =', nge,nmtx
  SAFE_ALLOCATE(eps, (nmtx,nmtx))
  write(*,*) ' '
  write(*,*) 'Allocated eps!', nmtx
  write(*,*) ' '
  read(10) (garbaged,ii=1,nge)
  read(10) (qk(ii),ii=1,3)
  do jj = 1, nmtx
    read(10) (eps(ii,jj),ii=1,nmtx)
  enddo
!      read(10) (totalreal,ii=1,nge)
  write(*,*) 'read in eps'
  isrtq=0

  call close_file(10)
!      call close_file(11)

  do ii = 1, nmtx
    do jj = 1, nmtx
      write(42,*) ii, jj, dble(eps(ii,jj))
      if (ii .eq. jj) then
        write(43,*) ii,oldx(isrtold(ii)),oldy(isrtold(ii)), &
          oldz(isrtold(ii)),dble(eps(ii,jj))
      endif
    end do
    write(42,*) ' '
  end do

  write(*,*) 'Finalizing'
  
end program ploteps
