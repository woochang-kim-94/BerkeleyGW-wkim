!===============================================================================
!
! Utilities:
!
! (1) epsinvomega      Originally By gsm      Last Modified 11/17/2010 (gsm)
!
! Serial code that plots frequency dependence of Plasmon-Pole, Retarded,
! Advanced epsilon inverse for a given q, G, G` vectors. Input parameters
! are read from file epsinvomega.inp in the working directory.
!
! Note that epsInvPP constructed with full-frequency epsmat will be
! different from epsInvPP constructed with static epsmat because of
! finite broadening used in full-frequency epsmat.
!
! FIXME: use WFN_RHO_VXC for reading RHO
!
! FIXME: [2013-01-15 gsm] compact form of GPP model doesn`t work when
! wtilde2 < 0 which is often the case for G .ne. G`; replace with the
! explicit expression from Hybertsen & Louie (see commented out code
! in the GPP section of epsinvomega.f90)
!
!===============================================================================

#include "f_defs.h"

program epsinvomega

  use global_m
  use epsread_hdf5_m
#ifdef HDF5
  use hdf5
#endif

  implicit none

  integer :: i,j,k,iq,mq,ig,igp,iunit
  complex(DPC) :: omega
  
  character*256, parameter :: fninp = "epsinvomega.inp"
  character*256 :: fneps,fnrho,fngpp,fnffr,fnffa
  character :: infile*80
  real(DP) :: q(3)
  integer :: g(3),gp(3),gmgp(3)
  
  integer :: freq_dep,nFreq_tmp,nFreq,nfreq_imag,ii,nq,ng,nmtx,kgrid(3)
  real(DP) :: dDeltaFreq,dBrdning,ecuts,delta
  real(DP), allocatable :: qpt(:,:)
  real(DP), allocatable :: ekin(:)
  real(DP), allocatable :: dFreqGrid(:)
  integer, allocatable :: isrtx(:)
  integer, allocatable :: isorti(:)
  integer, allocatable :: nmtx_of_q(:)
  SCALAR, allocatable :: eps(:,:)
  complex(DPC), allocatable :: dFreqBrd(:)
  complex(DPC), allocatable :: epsPP(:)
  complex(DPC), allocatable :: epsR(:)
  complex(DPC), allocatable :: epsA(:)
  complex(DPC), allocatable :: epsRTemp(:,:)
  complex(DPC), allocatable :: epsATemp(:,:)
  
  real(DP) :: bvec(3,3),bdot(3,3),blat,celvol,reccelvol,ebind
  integer :: nvecs,nspin,qgrid(3)
  
  integer nproc_para,num_gvec, error
  complex(DPC) :: epsStatic
  real(DP) :: rho0
  SCALAR :: rhogmgp
  integer, allocatable :: gvec(:,:)
  SCALAR, allocatable :: xcdum(:,:)
  
  real(DP) :: wp2,qg(3),qgp(3),qgqg,qgqgp,lambda,phi
  SCALAR :: Omega2,wtilde2,epsggp,I_epsggp,eps_static,eps_dynamic
  complex(DPC) :: wtilde2_temp
!  real(DP) :: epsPPRe,epsPPIm
!  complex(DPC) :: wtilde,ampl

!-----------------------------
! read input file

  write(6,'(/,1x,"reading",1x,a,1x,"file",/)')trim(fninp)
  
  call open_file(55,file=trim(fninp),form='formatted',status='old')
  
  read(55,'(a)') fneps
  read(55,'(a)') fnrho
  read(55,*) (q(i),i=1,3)
  read(55,*) (g(i),i=1,3)
  read(55,*) (gp(i),i=1,3)
  read(55,*) nFreq
  read(55,*) dDeltaFreq
  read(55,*) dBrdning
  read(55,*) ebind
  read(55,'(a)') fngpp
  read(55,'(a)') fnffr
  read(55,'(a)') fnffa
  
  call close_file(55)
  
  write(6,'(2a)')         "     eps  file  = ", trim(fneps)
  write(6,'(2a)')         "     rho  file  = ", trim(fnrho)
  write(6,'(a, 3f15.12)') "     q  vector  = ", q(1:3)
  write(6,'(a, 3i4)')     "     G  vector  = ", g(1:3)
  write(6,'(a, 3i4)')     "     G' vector  = ", gp(1:3)
  write(6,'(a, i6)')      "     nFreq      = ", nFreq
  write(6,'(a, f7.3)')    "     dDeltaFreq = ", dDeltaFreq
  write(6,'(a, f7.3)')    "     dBrdning   = ", dBrdning
  write(6,'(a, f7.3)')    "     ebind      = ", ebind
  write(6,'(2a)')         "     GPP file   = ", trim(fngpp)
  write(6,'(2a)')         "     FFR file   = ", trim(fnffr)
  write(6,'(2a,/)')       "     FFA file   = ", trim(fnffa)
  
  gmgp(:)=g(:)-gp(:)
  
!-----------------------------
! read eps file

  write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fneps)

#ifdef HDF5

  call h5open_f(error)

  infile = trim(fneps)

  call read_eps_grid_sizes_hdf5(ng, nq, ecuts, nfreq_tmp, nfreq_imag, nmtx, qgrid, freq_dep, infile)
  if (freq_dep==2 .and. nFreq_tmp/=nFreq) nFreq = nFreq_tmp

  SAFE_ALLOCATE(dFreqGrid,(nFreq))
  SAFE_ALLOCATE(dFreqBrd,(nFreq))
  SAFE_ALLOCATE(qpt, (3,nq))
  SAFE_ALLOCATE(nmtx_of_q, (nq))
  SAFE_ALLOCATE(gvec, (3,ng))

  call read_eps_qgrid_hdf5(nq,qpt,nmtx_of_q,infile)
  if (freq_dep .eq. 2) then
    call read_eps_freqgrid_hdf5(nFreq,dFreqGrid,dFreqBrd,infile)
  else
    do i=1,nFreq
      dFreqGrid(i)=dDeltaFreq*dble(i-1)
      dFreqBrd(i)=dBrdning
    enddo
  endif

  call read_eps_old_gvecs_hdf5(ng,gvec,infile)

#else

  iunit=12
  call open_file(unit=iunit,file=trim(fneps),form='unformatted',status='old')

  read(iunit)
  read(iunit) freq_dep,ii
  if (freq_dep.eq.2) then
    nFreq_tmp=ii
  endif
  if (freq_dep==2 .and. nFreq_tmp/=nFreq) nFreq = nFreq_tmp
  read(iunit) (kgrid(i),i=1,3)
  SAFE_ALLOCATE(dFreqGrid,(nFreq))
  SAFE_ALLOCATE(dFreqBrd,(nFreq))
  if (freq_dep.eq.2) then
    read(iunit) (dFreqGrid(i),i=1,nFreq),(dFreqBrd(i),i=1,nFreq)
    if (nFreq.gt.1) dDeltaFreq=dFreqGrid(2)-dFreqGrid(1)
    dBrdning=IMAG(dFreqBrd(1))
  else
    do i=1,nFreq
      dFreqGrid(i)=dDeltaFreq*dble(i-1)
      dFreqBrd(i)=dBrdning
    enddo
    read(iunit)
  endif
  read(iunit)
  read(iunit)
  read(iunit) ecuts
  read(iunit) nq
  read(iunit) ng
  
  rewind(iunit)
  
  SAFE_ALLOCATE(qpt, (3,nq))
  SAFE_ALLOCATE(gvec, (3,ng))
  
  read(iunit)
  read(iunit)
  read(iunit)
  read(iunit)
  read(iunit)
  read(iunit)
  read(iunit)
  read(iunit) nq,((qpt(j,i),j=1,3),i=1,nq)
  read(iunit) ng,((gvec(j,i),j=1,3),i=1,ng)

#endif
  
  ig=0
  igp=0
  do i=1,ng
    if (gvec(1,i).eq.g(1).and.gvec(2,i).eq.g(2).and. &
      gvec(3,i).eq.g(3)) ig=i
    if (gvec(1,i).eq.gp(1).and.gvec(2,i).eq.gp(2).and. &
      gvec(3,i).eq.gp(3)) igp=i
  enddo
  if (ig.eq.0) then
    call die("cannot find G vector in file " // trim(fneps))
  endif
  if (igp.eq.0) then
    call die("cannot find G' vector in file " // trim(fneps))
  endif
  
  mq=-1
  do iq=1,nq
    if (abs(q(1)-qpt(1,iq)).lt.TOL_Zero.and. &
      abs(q(2)-qpt(2,iq)).lt.TOL_Zero.and. &
      abs(q(3)-qpt(3,iq)).lt.TOL_Zero) then
      mq=iq-1
      exit
    endif
  enddo
  if (mq.eq.-1) then
    call die("cannot find q vector in file " // trim(fneps))
  endif

#ifdef HDF5

  SAFE_ALLOCATE(isrtx, (ng))
  SAFE_ALLOCATE(isorti, (ng))
  SAFE_ALLOCATE(ekin, (ng))
  call read_eps_gvecsofq_hdf5(ng,isrtx,isorti,ekin,mq+1,infile)
  nmtx = nmtx_of_q(mq+1)
  SAFE_DEALLOCATE(ekin)
  SAFE_DEALLOCATE(nmtx_of_q)

#else
  
  do iq=1,mq
    read(iunit) ng,nmtx
    read(iunit)
    read(iunit)
    if (freq_dep.eq.0) then
      do j=1,nmtx
        read(iunit)
      enddo
    endif
    if (freq_dep.eq.2) then
      do j=1,nmtx
        do i=1,nmtx
          read(iunit)
#ifdef CPLX
          read(iunit)
#endif
        enddo
      enddo
    endif
  enddo
  
  SAFE_ALLOCATE(isrtx, (ng))
  SAFE_ALLOCATE(isorti, (ng))
  
  read(iunit) ng,nmtx,(isrtx(i),isorti(i),i=1,ng)
  read(iunit)
  read(iunit)
 
#endif

  ig=isorti(ig)
  igp=isorti(igp)
  if (ig.eq.0) then
    call die("cannot find G vector in file " // trim(fneps))
  endif
  if (igp.eq.0) then
    call die("cannot find G' vector in file " // trim(fneps))
  endif
  
  SAFE_ALLOCATE(epsPP, (nFreq))
  SAFE_ALLOCATE(epsR, (nFreq))
  SAFE_ALLOCATE(epsA, (nFreq))
  
#ifdef HDF5

  if (freq_dep.eq.0) then
    SAFE_ALLOCATE(eps, (nmtx,1))
    call read_eps_matrix_col_hdf5(eps,igp,1,nmtx,mq+1,1,infile)
    epsStatic=eps(ig,1)
    SAFE_DEALLOCATE(eps)
  endif
  if (freq_dep.eq.2) then
    SAFE_ALLOCATE(epsRTemp, (nFreq,nmtx))
#ifdef CPLX
    SAFE_ALLOCATE(epsATemp, (nFreq,nmtx))
    call read_eps_matrix_col_f_hdf5(epsRTemp,nFreq,igp,nmtx,mq+1,1,infile,advanced=epsATemp)
    epsA(:)=epsATemp(:,ig)
    SAFE_DEALLOCATE(epsATemp)
#else
    call read_eps_matrix_col_f_hdf5(epsRTemp,nFreq,igp,nmtx,mq+1,1,infile)
#endif
    epsR(:)=epsRTemp(:,ig)
    SAFE_DEALLOCATE(epsRTemp)
    epsStatic=epsR(1)
  endif

#else

  if (freq_dep.eq.0) then
    SAFE_ALLOCATE(eps, (nmtx,1))
    do j=1,nmtx
      if (j.eq.igp) then
        read(iunit) (eps(i,1),i=1,nmtx)
        epsStatic=eps(ig,1)
      else
        read(iunit)
      endif
    enddo
    SAFE_DEALLOCATE(eps)
  endif
  if (freq_dep.eq.2) then
    do j=1,nmtx
      do i=1,nmtx
        if (j.eq.igp.and.i.eq.ig) then
          read(iunit) (epsR(k),k=1,nFreq)
          epsStatic=epsR(1)
        else
          read(iunit) 
        endif
      enddo
#ifdef CPLX
      do i=1,nmtx
        if (j.eq.igp.and.i.eq.ig) then
          read(iunit) (epsA(k),k=1,nFreq)
        else
          read(iunit)
        endif
      enddo
#endif
    enddo
#ifndef CPLX
    do i=1,nFreq
      epsA(i)=CONJG(epsR(i))
    enddo
#endif
  endif

  call close_file(iunit)

#endif
  
  SAFE_DEALLOCATE(isrtx)
  SAFE_DEALLOCATE(isorti)
  
  SAFE_DEALLOCATE(qpt)
  SAFE_DEALLOCATE(gvec)
  
  write(6,'(a,i6)')     "     omega num  = ", nFreq
  write(6,'(a,f7.3)')   "     omega step = ", dDeltaFreq
  write(6,'(a,f7.3,/)') "     omega brd  = ", dBrdning

!-----------------------------
! read rho file

  write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fnrho)
  
  iunit=95
  call open_file(unit=iunit,file=trim(fnrho),form='unformatted',status='old')
  
  rho0=0.0d0
  rhogmgp=ZERO

  read(iunit)
  read(iunit) nspin,nvecs
  read(iunit)
  read(iunit) celvol
  read(iunit) reccelvol, blat, ((bvec(i,j),i=1,3),j=1,3), ((bdot(i,j),i=1,3),j=1,3)
  read(iunit)
  read(iunit)
  read(iunit)
  SAFE_ALLOCATE(gvec, (3, nvecs))
  SAFE_ALLOCATE(xcdum, (nvecs, nspin))
  read(iunit) nproc_para
  ig = 1
  do i=1,nproc_para
    read(iunit) num_gvec
    read(iunit) ((gvec(k, j), k = 1, 3), j = ig, ig + num_gvec - 1)
    ig = ig + num_gvec
  enddo
  read(iunit) nproc_para
  ig = 1
  do i=1,nproc_para
    read(iunit) num_gvec
    read(iunit) ((xcdum(j, k), j = ig, ig + num_gvec - 1), k = 1, nspin)
    ig = ig + num_gvec
  enddo
  do j=1,nvecs
    if (gvec(1,j).eq.0.and.gvec(2,j).eq.0.and.gvec(3,j).eq.0) then
      do k=1,nspin
        rho0=rho0+dble(xcdum(j,k))
      enddo
    endif
    if (gvec(1,j).eq.gmgp(1).and.gvec(2,j).eq.gmgp(2).and.gvec(3,j).eq.gmgp(3)) then
      do k=1,nspin
        rhogmgp=rhogmgp+xcdum(j,k)
      enddo
    endif
  enddo
  SAFE_DEALLOCATE(gvec)
  SAFE_DEALLOCATE(xcdum)
    
  if (abs(rho0).le.TOL_Zero) then
    call die("cannot find rho(0) in file " // trim(fnrho))
  endif
  if (abs(rhogmgp).le.TOL_Zero) then
    call die("cannot find rho(G-G') in file " // trim(fnrho))
  endif
  
  call close_file(iunit)
  
  write(6,'(5x,"cel vol =",e20.12)') celvol
  write(6,'(5x,"rho(0) =",f10.3,/)') rho0

!-----------------------------
! construct generalized plasmon pole model

  write(6,'(1x,"constructing GPP model",/)')
  
  wp2=ryd*ryd*16.0d0*PI_D*rho0/celvol
  qg(:)=q(:)+dble(g(:))
  qgp(:)=q(:)+dble(gp(:))
  qgqg=dot_product(qg,matmul(bdot,qg))
  qgqgp=dot_product(qg,matmul(bdot,qgp))
  if (abs(qgqg) .lt. TOL_Zero) call die("GPP model diverges")
  Omega2=wp2*qgqgp/qgqg*rhogmgp/rho0
#ifdef CPLX
  epsggp = epsStatic
#else
  epsggp = dble(epsStatic)
#endif
  if (all(g(1:3) .eq. gp(1:3))) then
    delta = 1.0d0
  else
    delta = 0.0d0
  endif
  I_epsggp = delta - epsggp
  if (abs(I_epsggp) .lt. TOL_Small) call die("GPP model diverges")
#ifdef CPLX
! Complex GPP [PRB 40, 3162 (1989)]
  wtilde2_temp = Omega2 / I_epsggp
  lambda = abs(wtilde2_temp)
  if (lambda .lt. TOL_Small) call die("GPP model diverges")
  phi = atan2(IMAG(wtilde2_temp), dble(wtilde2_temp))
  if (abs(cos(phi)) .lt. TOL_Small) call die("GPP model diverges")
  wtilde2 = lambda / cos(phi)
  Omega2 = Omega2 * CMPLX(1.0d0, -tan(phi))
#else
! Real GPP [PRB 34, 5390 (1986)]
  wtilde2 = Omega2 / I_epsggp
  if (abs(wtilde2) .lt. TOL_Small) call die("GPP model diverges")
#endif
!  wtilde=dble(sqrt(COMPLEXIFY(wtilde2)))
!  ampl=-0.5d0*PI_D*sqrt(COMPLEXIFY(Omega2/wtilde2))
  do i=1,nFreq
    omega=dFreqGrid(i)
!    epsPPRe=delta+dble(Omega2/(omega**2-wtilde2))
!    epsPPIm=dble(ampl)/sqrt(PI_D)/dBrdning* &
!     (exp(-(omega-wtilde)**2/dBrdning**2) &
!     -exp(-(omega+wtilde)**2/dBrdning**2))
!    epsPP(i)=CMPLX(epsPPRe,epsPPIm)
! Instead of using the above, we write real and imaginary in a compact 
! form by adding dbrdning in denominator. The imaginary part should be 
! a sharp (delta function) at the plasmon frequency.
    epsPP(i)=delta+Omega2/(omega**2-wtilde2-CMPLX(0D0,dBrdning))
  enddo

  write(6,'(5x,"plasma frequency =",f10.3," eV",/)') sqrt(wp2)

!-----------------------------
! write generalized plasmon pole file

  write(6,'(1x,"writing",1x,a,1x,"file",/)')trim(fngpp)
    
  call open_file(unit=7,file=trim(fngpp),form='formatted',status='replace')
    
  do i=1,nFreq
    omega=dFreqGrid(i)
    write(7,100)omega,epsPP(i)
  enddo
    
  call close_file(7)

!-----------------------------
! write full frequency files

  if (freq_dep.eq.2) then
    write(6,'(1x,"writing",1x,a,1x,"and",1x,a,1x,"files",/)') &
      trim(fnffr),trim(fnffa)
    
    call open_file(unit=8,file=trim(fnffr),form='formatted',status='replace')
    call open_file(unit=9,file=trim(fnffa),form='formatted',status='replace')
    
    eps_static=epsR(1)
    eps_dynamic=eps_static
    
    do i=1,nFreq
      omega=dFreqGrid(i)+dFreqBrd(i)
      write(8,100)omega,epsR(i)
      write(9,100)omega,epsA(i)
      if (i .gt. 1) then
        eps_dynamic=eps_dynamic-(2.0d0/PI_D)*(dFreqGrid(i)-dFreqGrid(i-1))* &
          IMAG(epsR(i))*((1.0d0/omega)-1.0d0/(omega+ebind))
      endif
    enddo
    write(6,*) 'Static Head:', eps_static, 'Dynamic Head:', eps_dynamic
    write(6,*)
      
    call close_file(8)
    call close_file(9)
  endif

!-----------------------------
! deallocate and finish
  
  SAFE_DEALLOCATE(epsPP)
  SAFE_DEALLOCATE(epsR)
  SAFE_DEALLOCATE(epsA)
  SAFE_DEALLOCATE(dFreqGrid)
  SAFE_DEALLOCATE(dFreqBrd)
  
100 format(4f25.15)
  
end program epsinvomega
