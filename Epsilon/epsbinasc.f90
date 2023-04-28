!=========================================================================
!
! Utilities:
!
! (1) epsbinasc()  Originally by JLL       Last Modified 5/5/2008 (JRD)
!
!     This utility converts a binary epsmat file to an ascii epsmat file.
!     It uses the input file epsconv.inp.
!
!=========================================================================

#include "f_defs.h"

program epsbinasc

  use global_m
  implicit none

  character :: ajname*6,adate*11,outfile*20
  real(DP) :: ecuts1,ecuts2
  real(DP), allocatable :: dFreqGrid(:)
  complex(DPC), allocatable :: dFreqBrd(:)
  integer :: iunit,ig,im,iq,istart,j,jm, &
    nfiles,ng1,ng2,ngq1,ngq2,nmtx1,nmtx2,nq,nqtot, &
    ii,qgrid(3),freq_dep,nFreq,ijk
  
  character*120, allocatable :: filename(:)
  real(DP), allocatable :: ekin(:),q2(:,:)
  SCALAR, allocatable :: eps(:)
  integer, allocatable :: isort1(:),isort2(:),kx(:),ky(:),kz(:)

  logical :: file_exists

  write(6,*) 'This routine should only be used on non HDF5-based epsmat files.'
  write(6,*) 'If you built with HDF5 support (which is ideal), or otherwise'
  write(6,*) 'have HDF5-based epsmat files, you can use the h5dump command'
  write(6,*) 'to see your data in ascii format.'

  call open_file(55,file='epsconv.inp',form='formatted',status='old')
  read(55,*) nqtot
  SAFE_ALLOCATE(q2, (3,nqtot))
  read(55,'(a20)') outfile
  write(6,*) 'Output -> ',outfile
  write(6,*)

!--------------------------
! Find maximal values, and check consistency between
! the input file and the epsmat files...

  istart=1
  ng1=0
  ngq1=0
  nmtx1=0
  read(55,*) nfiles
  SAFE_ALLOCATE(filename, (nfiles))
  do iunit=1,nfiles
    read(55,'(a20)') filename(iunit)
#ifdef HDF5
    ! stop without error to allow the testsuite to run when compiled with -DHDF5
    INQUIRE(FILE=filename(iunit), EXIST=file_exists)
    if (.not. file_exists) then
      write(0,'(a)') 'WARNING: File "' // trim(filename(iunit)) // '" does not exist. Stopping.'
      stop
    endif
#endif
    write(6,*) 'Checking file ',TRUNC(filename(iunit))
    call open_file(unit=11,file=filename(iunit),form='unformatted',status='old')
    
    read(11) ajname,adate
    if(ajname /= 'chiGG0') then
      call die("Incorrect header '" // ajname // "' (must be 'chiGG0') in file '" // TRUNC(filename(iunit)) // "'")
    endif
    read(11) freq_dep,nFreq
    if (freq_dep.ne.0) then
      call die('epsbinasc: freq_dep')
    endif
    read(11) (qgrid(ii),ii=1,3)
    if (freq_dep .eq. 2) then
      SAFE_ALLOCATE(dFreqGrid,(nFreq))
      SAFE_ALLOCATE(dFreqBrd,(nFreq))
      read(11) (dFreqGrid(ijk),ijk=1,nFreq),(dFreqBrd(ijk),ijk=1,nFreq)
    else
      read(11)
    endif
    read(11)
    read(11)
    read(11) ecuts2
    
    if(iunit == 1) then
      ecuts1 = ecuts2
    else
      if(ecuts2.ne.ecuts1) then
        write(0,*) 'The cut-off in previous file (',ecuts1,') does not match ', &
          'the one in file ',TRUNC(filename(iunit)),' (',ecuts2,').'
        call die('epsbinasc cutoff mismatch')
      endif
    endif

    read(11) nq,((q2(j,iq),j=1,3),iq=istart,istart+nq-1)

    read(11) ng2
    if(ng1.eq.0) ng1=ng2
    if(ng2.ne.ng1) then
      call die('The number of G-vectors differs in epsmat files.')
    endif

    do iq=istart,istart+nq-1
      read(11) ngq2,nmtx2
      if(ngq1.lt.ngq2) ngq1=ngq2
      if(nmtx1.lt.nmtx2) nmtx1=nmtx2
      read(11)
      read(11) (q2(j,iq),j=1,3)
      do jm =1, nmtx2
        read(11)
      enddo
    enddo

    istart=istart+nq
    call close_file(11)
  enddo
  if(istart-1.ne.nqtot) then
    write(0,*) 'found = ', istart - 1, ' expected ', nqtot
    call die('Number of q-vectors found differs from number in input file.')
  endif

  SAFE_ALLOCATE(kx, (ng1))
  SAFE_ALLOCATE(ky, (ng1))
  SAFE_ALLOCATE(kz, (ng1))
  call open_file(unit=11,file=filename(1),form='unformatted',status='old')
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11) ng1,(kx(ig),ky(ig),kz(ig),ig=1,ng1)
  call close_file(11)
  SAFE_ALLOCATE(isort1, (ng1))
  SAFE_ALLOCATE(isort2, (ng1))
  SAFE_ALLOCATE(ekin, (ngq1))
  SAFE_ALLOCATE(eps, (nmtx1))
  
  call open_file(unit=12,file=outfile,form='formatted',status='replace')
  write(12,'(1x,a6,1x,a11)') ajname,adate
  write(12,*) freq_dep,nFreq
  write(12,*) (qgrid(ii),ii=1,3)
  if (freq_dep .eq. 2) then
    write(12,*) (dFreqGrid(ii),ii=1,nFreq),(dFreqBrd(ii),ii=1,nFreq)
  else
    write(12,*)
  endif
  write(12,*)
  write(12,*)
  write(12,*) ecuts1
  write(12,*) nqtot,((q2(j,iq),j=1,3),iq=1,nqtot)
  write(12,*) ng1,(kx(ig),ky(ig),kz(ig),ig=1,ng1)

  write(6,*)
  istart=1
  do iunit=1,nfiles
    write(6,*) 'Dealing with file ',TRUNC(filename(iunit))
    call open_file(unit=11,file=filename(iunit),form='unformatted',status='old')
    
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11) nq
    read(11)
    do iq=istart,istart+nq-1
      write(6,'(a,f9.6,3x,f9.6,3x,f9.6)') ' -> q=',(q2(j,iq),j=1,3)
      read(11) ngq2,nmtx2,(isort1(ig),isort2(ig),ig=1,ngq2)
      write(12,*) ngq2,nmtx2,(isort1(ig),isort2(ig),ig=1,ngq2)
      read(11) (ekin(ig),ig=1,ngq2)
      write(12,*) (ekin(ig),ig=1,ngq2)
      read(11)
      write(12,*) (q2(j,iq),j=1,3)
      do jm = 1, nmtx2
        read(11) (eps(im),im=1,nmtx2)
        write(12,*) (eps(im),im=1,nmtx2)
      enddo
    enddo

    istart=istart+nq
    call close_file(11)
  enddo
  call close_file(12)

end program epsbinasc
