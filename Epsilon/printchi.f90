!==============================================================================================
!
! Utilities:
!
! (1) printchi           Originally By JRD               Last Modified 6/30/2008 (JRD)
!
!==============================================================================================

#include "f_defs.h"

program printchi

  use global_m
  implicit none

  integer, allocatable :: irow(:), isrtx(:)
  real(DP) :: tnp(200, 3)
  real(DP), allocatable :: ekinx(:)
  integer :: ii, jj, nmtx, ntranq, mtrx(200, 3, 3), kg0(200, 3), ttt

  SCALAR, allocatable :: eps(:,:)

  character :: ajname*6, adate*11
  character :: aheadinput*6
  integer :: ps, ps2, np, i, k, j, n, nband
  integer :: freq_dep, nfreq, kgrid(3), nrk, ng, nFFTgridpts, FFTgrid(3), nq, indexq0
  real(DP) :: initfreq, deltafreq, brdning, ecuts, bdot(3, 3)
  integer, allocatable :: kgg(:, :), index_vec(:)
  real(DP), allocatable :: qpt(:, :)



  write(6,*) 'Welcome to Chi Printing'
  call open_file(unit=10,file='chi0mat',form='unformatted',status='old')
  ttt=1
  !       call open_file(unit=10,file='chi0mat',form='unformatted',status='old')
  !  ttt=nq-1! nq is number of q points in epsilon.inp      
  call open_file(unit=11,file='chi.gconv',status='replace')

  !       call open_file(unit=6, file='temp.dat')

  !       write(6, *) 1111

  !--------------------------
  ! Now read chi0mat

  !       read(10)
  !       read(10) ii
  !       if (ii.ne.0) then
  !         call die('Full frequency dependence not supported')
  !       endif
  !       read(10)
  !       read(10)
  !       read(10)
  !       read(10)
  !       read(10)
  !       read(10)
  !       read(10) nold
  !       read(10)
  !       read(10)
  !       read(10) nge
  !       rewind(10)
  ! !      deallocate(oldx)
  !       SAFE_ALLOCATE(oldx, (nold))
  ! !      deallocate(oldy)
  !       SAFE_ALLOCATE(oldy, (nold))
  ! !      deallocate(oldz)
  !       SAFE_ALLOCATE(oldz, (nold))
  ! !      deallocate(ekold)
  ! !      allocate(ekold(nge))
  ! !      deallocate(isrtold)
  !       SAFE_ALLOCATE(isrtold, (nge))
  ! !      deallocate(isrtq)
  ! !      allocate(isrtq(nge))
  ! !      deallocate(isrtq)
  !       SAFE_ALLOCATE(isrtq, (nge))

  read(10) aheadinput,ajname,adate
  write(11, *) 'line-1'
  write(11, *) aheadinput,ajname,adate
  read(10) freq_dep, nFreq, InitFreq, &
       DeltaFreq,Brdning
  write(11, *) 'line-2'
  write(11, *) freq_dep, nFreq, InitFreq, &
       DeltaFreq,Brdning
  read(10) (kgrid(i),i=1,3)
  write(11, *) 'line-3'
  write(11, *) (kgrid(i),i=1,3)
  read(10)
  write(11, *) 'line-4'
  write(11, *)
  read(10)
  write(11, *) 'line-5'      
  write(11, *)
  read(10)
  write(11, *) 'line-6'      
  write(11, *)
  read(10) ecuts, nband
  write(11, *) 'line-7'      
  write(11, *) ecuts, nband
  read(10) nrk ! also invflag is here, but not needed
  write(11, *) 'line-8'      
  write(11, *) nrk, 1 ! invflag
  read(10) ng, nFFTgridpts
  SAFE_ALLOCATE(kgg, (3, ng))
  SAFE_ALLOCATE(index_vec, (nFFTgridpts))
  backspace(10)
  read(10) ng,nFFTgridpts,(FFTgrid(i),i=1,3), &
       ((kgg(j,i),j=1,3),i=1,ng), &
       ((bdot(i,j),j=1,3),i=1,3),(index_vec(i),i=1,nFFTgridpts)
  write(11, *) 'line-9'      
  write(11, *) ng,nFFTgridpts,(FFTgrid(i),i=1,3), &
       ((kgg(j,i),j=1,3),i=1,ng), &
       ((bdot(i,j),j=1,3),i=1,3),(index_vec(i),i=1,nFFTgridpts)
  read(10) nq
  SAFE_ALLOCATE(qpt, (3, nq))
  backspace(10)
  read(10) nq,indexq0,((qpt(j,i),j=1,3),i=1,nq)
  write(11, *) 'line-10'      
  write(11, *) nq,indexq0,((qpt(j,i),j=1,3),i=1,nq)
  SAFE_DEALLOCATE(kgg)
  SAFE_DEALLOCATE(index_vec)
  SAFE_DEALLOCATE(qpt)

  !       write(6, *) 'temp1'

  ! Read q->0 dielectric matrix
  do ps2=1, ttt
     mtrx=0
     tnp=0.0d0
     kg0=0
     read(10) ntranq
     backspace(10)

     read(10) ntranq, (((mtrx(n,i,j),i=1,3),j=1,3), &
          (tnp(n,k), &
          kg0(n,k),k=1,3),n=1,ntranq)
     write(11, *) 'line-11', ps2
     write(11, *) ntranq, (((mtrx(n,i,j),i=1,3),j=1,3), &
          (tnp(n,k), &
          kg0(n,k),k=1,3),n=1,ntranq)
     read(10) nmtx
     backspace(10)
     SAFE_ALLOCATE(isrtx, (ng))
     SAFE_ALLOCATE(ekinx, (ng))
     SAFE_ALLOCATE(irow, (nmtx))
     read(10) nmtx, np,(isrtx(i),ekinx(i),i=1 &
          ,ng),(irow(i),i=1,nmtx)
     write(11, *) 'line-12', ps2
     write(11, *) nmtx, np,(isrtx(i),ekinx(i),i=1 &
          ,ng),(irow(i),i=1,nmtx)
     SAFE_ALLOCATE(eps, (nmtx,nmtx))
     write(11, *) 'line-13', ps2
     do jj = 1, nmtx
        read(10) (eps(ii,jj),ii=1,nmtx)
        do ps=1, nmtx
           write(11,*) jj, eps(ps,jj)
        end do
     enddo
     FLUSH(11)
     SAFE_DEALLOCATE(eps)
     SAFE_DEALLOCATE(isrtx)
     SAFE_DEALLOCATE(ekinx)
     SAFE_DEALLOCATE(irow)
  end do
  !      read(10) (totalreal,ii=1,nge)
  call close_file(10)
  call close_file(11)

end program printchi
