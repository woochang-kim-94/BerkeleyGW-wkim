!=========================================================================
!
! Utilities:
!
! (1) eps0sym     Originally by MJ      Last Modified 1/12/2010 (MJ)
!
!     This utility symmetrizes eps0mat file.
!
!=========================================================================

#include "f_defs.h"

program eps0sym

  use global_m
  use epsread_hdf5_m
  use epswrite_hdf5_m
#ifdef HDF5
  use hdf5
#endif
  use write_matrix_m

  implicit none
  
  character :: ajname*6,adate*11,outfile*80,infile*80
  real(DP)  :: ecuts
  real(DP), allocatable :: dFreqGrid(:), epsdiag(:,:,:)
  complex(DPC), allocatable :: dFreqBrd(:)
  integer   :: i,ig,j, &
    ng,ngq,nmtx, &
    ii,qgrid(3),freq_dep,nFreq,nfreq_imag,&
    jj,kk,ll,nargs,imap,iout, &
    nq, gx, gy, gz, gmax, &
    ig0, igx, igy, igz, jgx, jgy, jgz, jout  
  real(DP), allocatable :: ekin(:)
  real(DP) :: q(3,1), qk(3,1), errorMax
  SCALAR, allocatable :: eps(:,:), tempeps(:,:)
  SCALAR :: errorMaxTemp
  integer, allocatable :: isort(:),isorti(:),kx(:),ky(:),kz(:)
  integer, allocatable :: minusgidx(:), map(:), old(:,:)
  integer :: nmtx0_of_q(1), isize, error
  logical :: use_hdf5
#ifdef HDF5
  integer(HID_T) :: file_id
#endif
  
  write(6,*) MYFLAVOR // ' version is used to symmetrize file'
  
!------------------
! Get file names from command-line arguments

  nargs = command_argument_count()
  
  if (nargs .ne. 2) then
    call die('Usage: eps0sym eps0mat_in eps0mat_out')
  endif
  
  call get_command_argument(1,infile)
  call get_command_argument(2,outfile)
  
  use_hdf5 = .false.
#ifdef HDF5

  call h5open_f(error)

  call h5fopen_f(trim(infile), H5F_ACC_RDONLY_F, file_id, error)
  if(error == 0) use_hdf5 = .true.

  if(use_hdf5) then
    call h5fclose_f(file_id, error)
    call system("cp "//trim(infile)//" "//trim(outfile))
    call h5fopen_f(trim(outfile), H5F_ACC_RDWR_F, file_id, error)
    call read_eps_grid_sizes_hdf5(ngq, nq, ecuts, nfreq, nfreq_imag, nmtx, qgrid, freq_dep, infile)
    if (freq_dep.ne.0) then
      call die('eps0sym: Full frequency not supported')
    endif
    ng = ngq

    call read_eps_qgrid_hdf5(nq,q,nmtx0_of_q,infile)
    qk = q
  else

#endif
  
    call open_file(unit=11,file=TRUNC(infile),form='unformatted',status='old')

!--------------------------
!  Read in the eps0mat file
!

    read(11) ajname,adate
    read(11) freq_dep,nFreq
    if (freq_dep.ne.0) then
      call die('eps0sym: Full frequency not supported')
    endif
    read(11) (qgrid(ii),ii=1,3)
    if (freq_dep .eq. 2) then
      SAFE_ALLOCATE(dFreqGrid,(nFreq))
      SAFE_ALLOCATE(dFreqBrd,(nFreq))
      read(11) (dFreqGrid(i),i=1,nFreq),(dFreqBrd(i),i=1,nFreq)
    else
      read(11)
    endif
    read(11)
    read(11)
    read(11) ecuts
    read(11) nq,(q(j,1),j=1,3)
    !if (nq .ne. 1) then
    !  call die('This only works for the q->0 point.')
    !endif
    read(11) ng
    read(11) ngq,nmtx
    
    call close_file(11)

#ifdef HDF5
  endif
#endif

  SAFE_ALLOCATE(kx, (ng))
  SAFE_ALLOCATE(ky, (ng))
  SAFE_ALLOCATE(kz, (ng))
  SAFE_ALLOCATE(isort, (ngq))
  SAFE_ALLOCATE(isorti, (ngq))
  SAFE_ALLOCATE(ekin, (ngq))
  SAFE_ALLOCATE(eps, (nmtx,nmtx))
 
#ifdef HDF5
  if(use_hdf5) then

    SAFE_ALLOCATE(old,(3,ngq))
    call read_eps_old_gvecs_hdf5(ngq,old,infile)
    kx(:)=old(1,:)
    ky(:)=old(2,:)
    kz(:)=old(3,:)

    call read_eps_gvecsofq_hdf5(ngq,isort,isorti,ekin,1,infile)
    
    call read_eps_matrix_ser_hdf5(eps, nmtx, 1, 1, infile)

  else
#endif
 
    call open_file(unit=11,file=TRUNC(infile),form='unformatted',status='old')
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11) ng,(kx(i),ky(i),kz(i),i=1,ng)
    read(11) ngq,nmtx,(isort(ig),isorti(ig),ig=1,ngq)
    read(11) (ekin(ig),ig=1,ngq)
    read(11) (qk(j,1),j=1,3)
    do jj = 1, nmtx
      read(11) (eps(ii,jj),ii=1,nmtx)
    enddo
    call close_file(11)

#ifdef HDF5
  endif
#endif
  
  ! Since we want the q=0 dielectric function but we have the
  ! q0<>0 but small dielectric, we can average over q0 and -q0
  ! to get a better dielectric (linear terms in q0 will be canceled)

! Calculate the maximum
  gmax = 0
  do ii = 1,ng
    if (abs(kx(ii)) > gmax) gmax = abs(kx(ii))
    if (abs(ky(ii)) > gmax) gmax = abs(ky(ii))
    if (abs(kz(ii)) > gmax) gmax = abs(kz(ii))
  enddo
! Create a map
!  write(6,*) gmax
  SAFE_ALLOCATE(map, ((2*gmax+1)*(2*gmax+1)*(2*gmax+1)))
  map = 0
  do ii = 1, ngq
    iout = isort(ii)
    if (iout.eq.0) cycle
    gx = kx(iout)
    gy = ky(iout)
    gz = kz(iout)
    imap = ((gx+gmax)*(2*gmax+1)+gy+gmax)*(2*gmax+1)+gz+gmax+1
    map(imap) = ii
  enddo

! For each g, find -g in the gvector list
! and also mark which ii corresponds to G=0

  SAFE_ALLOCATE(minusgidx, (nmtx))
  minusgidx = 0
  do ii=1,nmtx
    iout = isort(ii)
    gx = kx(iout)
    gy = ky(iout)
    gz = kz(iout)
    imap = ((-gx+gmax)*(2*gmax+1)-gy+gmax)*(2*gmax+1)-gz+gmax+1
    minusgidx(ii) = map(imap)
    if (gx .eq. 0 .and. gy .eq. 0 .and. gz .eq. 0) ig0 = ii
  enddo
  
!  do ii = 1, min(100,nmtx)
!    iout = isort(ii)
!    gx = kx(iout)
!    gy = ky(iout)
!    gz = kz(iout)
!    mgx = kx(isort(minusgidx(ii)))
!    mgy = ky(isort(minusgidx(ii)))
!    mgz = kz(isort(minusgidx(ii)))
!    write(6,'(7i4)') ii, gx, gy , gz, mgx, mgy, mgz
!  enddo

! Set the wings to zero
! This is as per Baldereschi and Tosatti, PRB 17, 4710 (1978)
! This is perhaps not the correct thing to do. One should 
! still symmetrize epsilon - but the wings should come out
! whatever they need to be automatically. What is given in that
! paper by Baldereschi and Tosatti is at q=0 and the dielectric
! function is weird at that point. What is needed in the GW
! code is the average in the minibz...
  !do ii=1,nmtx
  !  do jj=1,nmtx
  !    if (ii .eq. ig0 .and. jj .eq. ig0) cycle
  !    if (ii .eq. ig0) eps(ii,jj) = 0.0d0
  !    if (jj .eq. ig0) eps(ii,jj) = 0.0d0
  !  enddo
  !enddo
  
! Copy eps(q0) into a temporary

  SAFE_ALLOCATE(tempeps, (nmtx,nmtx))
  tempeps = eps

  errorMax=0D0
  errorMaxTemp=0D0
  
  write(6,'(5a4,3x,3a4,1x)',advance='no') "ig", "ig'", "Gx", "Gy", "Gz", "G'x", "G'y", "G'z"
  ! handle the fact that the two components of the complex numbers will be written in the complex case
#ifdef CPLX  
  write(6,'(6a16)') "Re eps(G,G')", "Im eps(G,G')", "Re eps(-G,-G')", "-Im eps(-G,-G')", "Re diff", "Im diff"
#else
  write(6,'(3a16)') "eps(G,G')", "eps(-G,-G')*", "difference"
#endif

!  do ii = 1, min(100,nmtx)
!    do jj = 1, min(100,nmtx)
  do ii = 1, nmtx
    do jj = 1, nmtx
      iout = isort(ii)
      igx = kx(iout)
      igy = ky(iout)
      igz = kz(iout)
      jout = isort(jj)
      jgx = kx(jout)
      jgy = ky(jout)
      jgz = kz(jout)
      kk = minusgidx(ii)
      ll = minusgidx(jj)
      if ((kk .le. nmtx) .and. (ll .le. nmtx)) then
        errorMaxTemp=eps(ii,jj)-MYCONJG(eps(kk,ll))
        if (abs(errorMaxTemp) .gt. errorMax) errorMax = abs(errorMaxTemp)
        if (ii .lt. 20 .and. jj .lt. 20) then
          write(6,'(5i4,3x,3i4,1x,6E16.5)') ii,jj,igx,igy,igz,jgx,jgy,jgz,eps(ii,jj),MYCONJG(eps(kk,ll)),errorMaxTemp
        endif
      endif
    enddo
  enddo
 
  write(6,*) "Error = eps(G,G')-eps*(-G,-G')"
  write(6,'("The max error in your matrix is",E12.5)') errorMax
  write(6,*) "Symmetrizing the matrix"

! Now add in contribution from -q0 which means conjg(eps(-g,-gp))

  do ii=1,nmtx
    do jj=1,nmtx

      kk = minusgidx(ii)
      ll = minusgidx(jj)
      if ((kk .le. nmtx) .and. (ll .le. nmtx)) then
        tempeps(ii,jj) = tempeps(ii,jj) + MYCONJG(eps(kk,ll))
      endif
    enddo
  enddo
! Average over q0 and -q0 and put back into eps

  eps = 0.5d0*tempeps
  SAFE_DEALLOCATE(tempeps)
  SAFE_DEALLOCATE(minusgidx)
  
#ifdef HDF5
  if(use_hdf5) then
    SAFE_DEALLOCATE(old)
    call write_gvec_indices_hdf(ng,isort,isorti,ekin,1,outfile)

    isize=SCALARSIZE

    call write_matrix_ser_hdf(eps, nmtx, 1, 1, outfile)

    SAFE_ALLOCATE(epsdiag,(isize,nmtx,1))

    do ii = 1, nmtx
      epsdiag(1,ii,1) = dble(eps(ii,ii))
#ifdef CPLX
      epsdiag(2,ii,1) = IMAG(eps(ii,ii))
#endif
    enddo
    
    call write_matrix_diagonal_hdf(epsdiag, nmtx, 1, isize, outfile)

    SAFE_DEALLOCATE(epsdiag)

    call h5close_f(error)
  else
#endif

    call open_file(unit=20,file=TRUNC(outfile),form='unformatted',status='replace')

    ajname='chiGG0'
    write(20) ajname,adate
    write(20) freq_dep,nFreq
    write(20) (qgrid(ii),ii=1,3)
    if (freq_dep .eq. 2) then
      write(20) (dFreqGrid(i),i=1,nFreq),(dFreqBrd(i),i=1,nFreq)
    else
      write(20)
    endif
    write(20)
    write(20)
    write(20) ecuts
    write(20) nq,(q(j,1),j=1,3)
    write(20) ng,(kx(ig),ky(ig),kz(ig),ig=1,ng)
    write(20) ngq,nmtx,(isort(ig),isorti(ig),ig=1,ngq)
    write(20) (ekin(ig),ig=1,ngq)
    write(20) (qk(j,1),j=1,3)
    do jj = 1, nmtx
      write(20) (eps(ii,jj),ii=1,nmtx)
    enddo
  
    call close_file(20)

#ifdef HDF5
  endif
#endif

end program eps0sym
