!=========================================================================
!
! Routines:
!
! (1) offdiag           Originally By gsm        Last Modified 5/8/2009 (gsm)
!
! This routine reads in sigma_hp.log file, builds the Sigma matrix for each
! k-point and for each finite difference point (Ecor - dE, Ecor, Ecor + dE)
! according to Eq. (6) of Rohlfing & Louie PRB 62 4927, diagonalizes it with
! LAPACK and writes out the eigenvalues. This is a serial program, no MPI.
! If toff = -1/+1 the Hermitian matrix is constructed from the lower/upper
! triangle. No input file or command-line arguments are needed.
!
!=========================================================================

#include "f_defs.h"

program offdiag

  use global_m
  use lapack_m
  implicit none

  integer :: iunit,ierr,spin,ii,jj,ll,ik
  integer :: freq_dep,bmin,bmax,loff,toff,fdf
  integer :: info,lda,ldvl,ldvr,lwork,nband,iw,nstart,nend
  real(DP) :: elda,ecor,exch,sx,ch,sig,vxc,eqp0,eqp1,z1,z2
  SCALAR :: sig3,vxc3
  real(DP) :: max6(6)
  integer :: idx6(6,2)
  character*256 :: fl,str,tag
  integer, allocatable :: isort(:)
  SCALAR, allocatable :: ham(:,:),alda(:,:),vl(:,:),vr(:,:),work(:)
#ifdef CPLX
  real(DP) :: x2,sx2,ch2,sig2,vxc2,i2,j2,l2
  real(DP), allocatable :: rwork(:)
  complex(DPC), allocatable :: ww(:)
#else
  real(DP), allocatable :: wi(:),wr(:)
#endif  

  iunit=21
  fl = "sigma_hp.log"
  max6(:) = 0.0d0
  idx6(:,:)=0

  freq_dep=-2
  bmin=0
  bmax=0
  loff=-3
  toff=-2
  fdf=-3
  call open_file(unit=iunit,file=fl,status='old',form='formatted')
  ierr=0
  do while (ierr.eq.0)
    read(iunit,'(a)',iostat=ierr) str
    if (str(2:21).eq."frequency_dependence") then
      read(str(22:),*) freq_dep
      write(6,'(a,i2)') 'frequency_dependence = ', freq_dep
    endif
    if (str(2:11).eq."band_index") then
      read(str(12:),*) bmin, bmax
      write(6,'(a,i6,a,i6)') 'band min = ', bmin, ', band max = ', bmax
    endif
    if (str(2:13).eq."sigma_matrix") then
      read(str(14:),*) loff, toff
      write(6,'(a,i6,a,i6)') 'loff = ', loff, ', toff = ', toff
    endif
    if (str(2:23).eq."finite_difference_form") then
      read(str(24:),*) fdf
      write(6,*) 'finite difference form = ', fdf
    endif
  enddo
  call close_file(unit=iunit)
  if(freq_dep.lt.-1.or.freq_dep.gt.2) call die("unknown frequency dependence")
  if(bmin.lt.1.or.bmax.lt.bmin) call die("bmin out of range")
  if(loff.lt.-2.or.(loff.gt.0.and.loff.lt.bmin).or.loff.gt.bmax) call die("loff out of range")
  if(loff < 0) write(0,'(a)') 'WARNING: Sigma is not Hermitian unless all matrix elements are evaluated at the same energy.'
  if(toff.lt.-1.or.toff.gt.1) call die("toff out of range")
  if(fdf.lt.-2.or.fdf.gt.2) call die("fdf out of range")
  if (freq_dep .eq. 2) call die("Full frequency is not supported")
  
  if(loff == 0) call die("You need to use sigma_matrix in the sigma run to be able to use offdiag.")

  if (fdf.eq.-1) then
    nstart = 1
    nend = 2
  elseif (fdf.eq.0) then
    nstart = 1
    nend = 3
  elseif (fdf.eq.1) then
    nstart = 2
    nend = 3
  else
    nstart = 2
    nend = 2
  endif
  
  nband=bmax-bmin+1
  lda=nband
  ldvl=1
  ldvr=1
  lwork=2*nband
  
  SAFE_ALLOCATE(isort, (nband))
  SAFE_ALLOCATE(ham, (lda,nband))
  SAFE_ALLOCATE(alda, (lda,nband))
  SAFE_ALLOCATE(vl, (ldvl,nband))
  SAFE_ALLOCATE(vr, (ldvr,nband))
  SAFE_ALLOCATE(work, (lwork))
  lwork=-1
#ifdef CPLX
  SAFE_ALLOCATE(ww, (nband))
  SAFE_ALLOCATE(rwork, (2*nband))
  call zgeev('N','N',nband,ham,lda,ww,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
  if (info.eq.0) lwork=int(dble(work(1)))
#else
  SAFE_ALLOCATE(wr, (nband))
  SAFE_ALLOCATE(wi, (nband))
  call dgeev('N','N',nband,ham,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
  if (info.eq.0) lwork=int(work(1))
#endif
  if (lwork.lt.1) then
    lwork=2*nband
  else
    SAFE_DEALLOCATE(work)
    SAFE_ALLOCATE(work, (lwork))
  endif
  
  alda(:,:)=0.0d0
  
  call open_file(unit=iunit,file=fl,status='old',form='formatted')
  write(6,*) "Reading matrix elements from " // trim(fl)
  ierr=0
  do while (ierr.eq.0)
    read(iunit,'(a)',iostat=ierr) str
    if (str(42:45).eq."ik =") then
      read(str(47:49),*) ik
      read(str(57:),*) spin
      write(6,205) ik, spin
      read(iunit,*)
      read(iunit,*)
      do
        read(iunit,'(a)') str
        if (len(trim(str)).eq.0) exit
        read(str,*) ii,elda,ecor,exch,sx,ch,sig,vxc,eqp0,eqp1
        alda(ii-bmin+1,ii-bmin+1)=elda
      enddo
      do iw = nstart, nend
        ham=alda
        if (iw.eq.1) write(6,501)
        if (iw.eq.2) write(6,502)
        if (iw.eq.3) write(6,503)
        read(iunit,*)
        read(iunit,*)
        do
          read(iunit,'(a)') str
          if (len(trim(str)).eq.0) exit
          read(str,*) ii,jj,ll,tag,exch,sx,ch,sig,vxc
          if(tag(1:4) /= 'real') call die("Incorrect tag " // TRUNC(tag) // " found in place of 'real'.")
#ifdef CPLX
          read(iunit,'(a)') str
          read(str,*) i2,j2,l2,tag,x2,sx2,ch2,sig2,vxc2
          if(tag(1:4) /= 'imag') call die("Incorrect tag " // TRUNC(tag) // " found in place of 'imag'.")
          sig3=CMPLX(sig,sig2)
          vxc3=CMPLX(vxc,vxc2)
#else
          sig3=sig
          vxc3=vxc
#endif
          ham(ii-bmin+1,jj-bmin+1)=ham(ii-bmin+1,jj-bmin+1)+sig3-vxc3
          if (ii .eq. jj) then
            if (max6(1) .lt. abs(sig3)) then
              max6(1) = abs(sig3)
              idx6(1,1) = ii
            endif
            if (max6(2) .lt. abs(vxc3)) then
              max6(2) = abs(vxc3)
              idx6(2,1) = ii
            endif
            if (max6(3) .lt. abs(sig3 - vxc3)) then
              max6(3) = abs(sig3 - vxc3)
              idx6(3,1) = ii
            endif
          else
            if (max6(4) .lt. abs(sig3)) then
              max6(4) = abs(sig3)
              idx6(4,1) = ii
              idx6(4,2) = jj
            endif
            if (max6(5) .lt. abs(vxc3)) then
              max6(5) = abs(vxc3)
              idx6(5,1) = ii
              idx6(5,2) = jj
            endif
            if (max6(6) .lt. abs(sig3 - vxc3)) then
              max6(6) = abs(sig3 - vxc3)
              idx6(6,1) = ii
              idx6(6,2) = jj
            endif
          endif
        enddo
        !
        ! construct the Hermitian matrix from the lower triangle
        !
        if (toff.eq.-1) then
          do ii=1,nband
            do jj=ii+1,nband
              ham(ii,jj)=MYCONJG(ham(jj,ii))
            enddo
          enddo
        endif
        !
        ! construct the Hermitian matrix from the upper triangle
        !
        if (toff.eq.1) then
          do ii=1,nband
            do jj=1,ii-1
              ham(ii,jj)=MYCONJG(ham(jj,ii))
            enddo
          enddo
        endif
        !
        ! diagonalize with LAPACK
        !
#ifdef CPLX
        call zgeev('N','N',nband,ham,lda,ww,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
#else
        call dgeev('N','N',nband,ham,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
#endif
        !
        ! sort and output eigenvalues
        !
        if (info.eq.0) then
          do ii=1,nband
            isort(ii)=ii
          enddo
          do ii=1,nband-1
            ll=0
            z1=1.0d6
            do jj=ii,nband
#ifdef CPLX
              z2=dble(ww(isort(jj)))
#else
              z2=wr(isort(jj))
#endif
              if (z2.lt.z1) then
                ll=jj
                z1=z2
              endif
            enddo
            if (ll.gt.0) then
              jj=isort(ii)
              isort(ii)=isort(ll)
              isort(ll)=jj
            endif
          enddo
          do ii=1,nband
#ifdef CPLX
            write(6,209) ii,ww(isort(ii))
#else
            write(6,209) ii,wr(isort(ii)),wi(isort(ii))
#endif
          enddo
          write(6,*)
          write(6,701)max6(1),idx6(1,1)
          write(6,702)max6(2),idx6(2,1)
          write(6,703)max6(3),idx6(3,1)
          write(6,704)max6(4),idx6(4,1),idx6(4,2)
          write(6,705)max6(5),idx6(5,1),idx6(5,2)
          write(6,706)max6(6),idx6(6,1),idx6(6,2)
        else
          write(0,'(a,i6)') 'LAPACK error code = ', info
          call die("Failed to diagonalize Sigma matrix.")
        endif
      enddo
    endif
  enddo
  call close_file(unit=iunit)
  
  SAFE_DEALLOCATE(isort)
  SAFE_DEALLOCATE(ham)
  SAFE_DEALLOCATE(alda)
  SAFE_DEALLOCATE(vl)
  SAFE_DEALLOCATE(vr)
  SAFE_DEALLOCATE(work)
#ifdef CPLX
  SAFE_DEALLOCATE(ww)
  SAFE_DEALLOCATE(rwork)
#else
  SAFE_DEALLOCATE(wr)
  SAFE_DEALLOCATE(wi)
#endif

205 format(/,1x,"ik =",i3,1x,"s =",i2)
209 format(1x,i4,2f12.6)
501 format(/,1x,"Sig(Eo - dE)",/)
502 format(/,1x,"Sig(Eo)",/)
503 format(/,1x,"Sig(Eo + dE)",/)
701 format(1x,"max_n  |<n|Sig|n>|       =",f10.6," eV at n =",i4)
702 format(1x,"max_n  |<n|Vxc|n>|       =",f10.6," eV at n =",i4)
703 format(1x,"max_n  |<n|Sig - Vxc|n>| =",f10.6," eV at n =",i4)
704 format(1x,"max_nm |<n|Sig|m>|       =",f10.6," eV at n =",i4," m =",i4)
705 format(1x,"max_nm |<n|Vxc|m>|       =",f10.6," eV at n =",i4," m =",i4)
706 format(1x,"max_nm |<n|Sig - Vxc|m>| =",f10.6," eV at n =",i4," m =",i4)

end program offdiag
