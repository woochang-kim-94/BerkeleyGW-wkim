!===============================================================================
!
! Utilities:
!
! (1) epsomega      Originally By gsm      Last Modified 11/18/2010 (gsm)
!
! Parallel code that plots frequency dependence of Plasmon-Pole, Retarded,
! Advanced epsilon for a given q, G, G` vectors. Input parameters are read
! from file epsomega.inp in the working directory.
!
! Note that epsPP constructed with full-frequency epsmat will be
! different from epsPP constructed with static epsmat because of
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

program epsomega

  use global_m
  use input_utils_m
  use inversion_m
  use read_matrix_m
  use scalapack_m
  use epsread_hdf5_m
#ifdef HDF5
  use hdf5
#endif

  implicit none
  
  logical :: sflag
  integer :: i,j,k,iq,mq,ig,igp,jg,jgp,iunit,ierr
  integer :: icurr,icol,irow,icolm,irowm,icols,irows
  integer :: iout,nFFTgridpts,FFTgrid(3)
  real(DP) :: omega
  
  character*256, parameter :: fninp = "epsomega.inp"
  character*256 :: fneps,fnrho,fngpp,fnffr,fnffa
  real(DP) :: q(3)
  integer :: g(3),gp(3)

  character :: infile*80
  
  type(gspace) :: gvec
  integer :: freq_dep,nFreq,nfreq_imag,ii,nq,ng,nmtx,kgrid(3)
  real(DP) :: dDeltaFreq,dBrdning,ecuts,delta,qvec(3)
  real(DP), allocatable :: qpt(:,:)
  real(DP), allocatable :: ekin(:)
  real(DP), allocatable :: dFreqGrid(:)
  integer, allocatable :: isrtx(:)
  integer, allocatable :: isorti(:)
  integer, allocatable :: nmtx_of_q(:)
  SCALAR, allocatable :: eps(:,:)
  complex(DPC), allocatable :: dFreqBrd(:)
  complex(DPC), allocatable :: epsPP(:,:,:)
  complex(DPC), allocatable :: epsR(:,:,:)
  complex(DPC), allocatable :: epsA(:,:,:)
  complex(DPC), allocatable :: epsAux(:,:)
  SCALAR, allocatable :: rhog(:)  
  real(DP) :: bvec(3,3),bdot(3,3),blat,celvol,reccelvol
  integer :: nvecs,nspin
  
  integer :: nproc_para,num_gvec,gx,gy,gz,qgrid(3),error
  real(DP) :: rho0
  integer, allocatable :: gvec_rho(:,:)
  SCALAR, allocatable :: xcdum(:,:)
  
  real(DP) :: wp2,qg(3),qgp(3),qgqg,qgqgp,lambda,phi
  SCALAR :: Omega2,wtilde2,epsggp,I_epsggp
  complex(DPC) :: wtilde2_temp

  type (scalapack) :: scal

  call peinfo_init()

!-----------------------------
! read input file

  if (peinf%inode.eq.0) then
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
    write(6,'(2a)')         "     GPP file   = ", trim(fngpp)
    write(6,'(2a)')         "     FFR file   = ", trim(fnffr)
    write(6,'(2a,/)')       "     FFA file   = ", trim(fnffa)
  endif
  
#ifdef MPI
  call MPI_Bcast(fneps,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(fnrho,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(q,3,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(g,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(gp,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(fngpp,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(fnffr,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(fnffa,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
#endif

!-----------------------------
! read eps file

  if (peinf%inode.eq.0) then
    write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fneps)

#ifdef HDF5

    ! JRD: There is no equivalent on read_matrix_d
    ! yet. So, the beginings of HDF5 support are here.
    ! But, it is not complete
    call die("No Support for HDF5 Yet. epsinvomega routine does support HDF5.")

    call h5open_f(error)

    infile = trim(fneps)

    call read_eps_grid_sizes_hdf5(ng, nq, ecuts, nfreq, nfreq_imag, nmtx, qgrid, freq_dep, infile)
    SAFE_ALLOCATE(dFreqGrid,(nFreq))
    SAFE_ALLOCATE(dFreqBrd,(nFreq))
    SAFE_ALLOCATE(qpt, (3,nq))
    SAFE_ALLOCATE(nmtx_of_q, (nq))
    SAFE_ALLOCATE(gvec%components, (3,ng))

    call read_eps_qgrid_hdf5(nq,qpt,nmtx_of_q,infile)
    if (freq_dep .ne. 0) then
      call read_eps_freqgrid_hdf5(nFreq,dFreqGrid,dFreqBrd,infile)
    endif

#else
    
    iunit=12
    call open_file(unit=iunit,file=trim(fneps),form='unformatted',status='old')
    
    read(iunit)
    read(iunit) freq_dep,ii
    if (freq_dep.ne.0) then
      nFreq=ii
    endif
    read(iunit) (kgrid(i),i=1,3)
    SAFE_ALLOCATE(dFreqGrid,(nFreq))
    SAFE_ALLOCATE(dFreqBrd,(nFreq))
    if (freq_dep.ne.0) then
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
    SAFE_ALLOCATE(gvec%components, (3,ng))
    
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit) nq,((qpt(j,i),j=1,3),i=1,nq)
    read(iunit) ng,((gvec%components(j,i),j=1,3),i=1,ng)

#endif
    
    ig=0
    igp=0
    do i=1,ng
      if (gvec%components(1,i).eq.g(1).and.gvec%components(2,i).eq.g(2).and. &
        gvec%components(3,i).eq.g(3)) ig=i
      if (gvec%components(1,i).eq.gp(1).and.gvec%components(2,i).eq.gp(2).and. &
        gvec%components(3,i).eq.gp(3)) igp=i
    enddo
    
    mq=-1
    do iq=1,nq
      if (abs(q(1)-qpt(1,iq)).lt.TOL_Zero.and. &
        abs(q(2)-qpt(2,iq)).lt.TOL_Zero.and. &
        abs(q(3)-qpt(3,iq)).lt.TOL_Zero) then
        mq=iq-1
        exit
      endif
    enddo
    
    SAFE_DEALLOCATE(qpt)
    
    write(6,'(a,i6)')     "     omega num  = ", nFreq
  endif
  
#ifdef MPI
  call MPI_Bcast(freq_dep,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(nFreq,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(dDeltaFreq,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(dBrdning,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(mq,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(ig,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(igp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(ng,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  if (peinf%inode.ne.0) then
    SAFE_ALLOCATE(dFreqGrid,(nFreq))
    SAFE_ALLOCATE(dFreqBrd,(nFreq))
    SAFE_ALLOCATE(gvec%components, (3,ng))
  endif
  call MPI_Bcast(dFreqGrid,nFreq,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(dFreqBrd,nFreq,MPI_COMPLEX_DPC,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(gvec%components,3*ng,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  
  if (mq.eq.-1) then
    call die("cannot find q vector in file " // trim(fneps))
  endif
  if (ig.eq.0) then
    call die("cannot find G vector in file " // trim(fneps))
  endif
  if (igp.eq.0) then
    call die("cannot find G' vector in file " // trim(fneps))
  endif
  
  if (peinf%inode.eq.0) then

#ifdef HDF5

    SAFE_ALLOCATE(isrtx, (ng))
    SAFE_ALLOCATE(isorti, (ng))
    SAFE_ALLOCATE(ekin, (ng))
    call read_eps_gvecsofq_hdf5(ng,isrtx,isorti,ekin,mq+1,infile)
    nmtx = nmtx_of_q(mq+1)
    SAFE_DEALLOCATE(ekin)
    SAFE_DEALLOCATE(isorti)
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
      else
        do j=1,nmtx
          do i=1,nmtx
            read(iunit)
#ifdef CPLX
            if(freq_dep.ne.3) then
              read(iunit)
            endif
#endif
          enddo
        enddo
      endif
    enddo
    
    SAFE_ALLOCATE(isrtx, (ng))
    SAFE_ALLOCATE(isorti, (ng))
    isrtx(:)=0
    isorti(:)=0
    
    read(iunit) ng,nmtx,(isrtx(i),isorti(i),i=1,ng)
    read(iunit)
    read(iunit)
    
    ig=isorti(ig)
    igp=isorti(igp)
    
    SAFE_DEALLOCATE(isorti)

#endif

  endif
  
#ifdef MPI
  call MPI_Bcast(nmtx,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(ig,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(igp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  if (peinf%inode.ne.0) then
    SAFE_ALLOCATE(isrtx, (ng))
  endif
  call MPI_Bcast(isrtx,ng,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  
  if (ig.eq.0) then
    call die("cannot find G vector in file " // trim(fneps))
  endif
  if (igp.eq.0) then
    call die("cannot find G' vector in file " // trim(fneps))
  endif
  
  call blacs_setup(scal, nmtx, .true.)

!-----------------------------
! read distributed matrices

  SAFE_ALLOCATE(eps, (scal%npr,scal%npc))
  SAFE_ALLOCATE(epsPP, (nFreq,scal%npr,scal%npc))
  if (freq_dep.ne.0) then
    SAFE_ALLOCATE(epsR, (nFreq,scal%npr,scal%npc))
    SAFE_ALLOCATE(epsA, (nFreq,scal%npr,scal%npc))
  endif
  SAFE_ALLOCATE(epsAux, (scal%npr,scal%npc))
  
  if (freq_dep.eq.0) then
    call read_matrix_d(scal,eps,nmtx,iunit)
  else
    peinf%rank_mtxel=0
    call read_matrix_f(scal, nFreq, nFreq, epsR, nmtx, 1, iunit, advanced=epsA)
#ifndef CPLX
    do icolm=1,scal%npc
      do irowm=1,scal%npr
        do i=1,nFreq
          epsA(i,irowm,icolm)=CONJG(epsR(i,irowm,icolm))
        enddo
      enddo
    enddo
#else
    if (freq_dep.eq.3) then
      do icolm=1,scal%npc
        do irowm=1,scal%npr
          do i=1,nFreq
            epsA(i,irowm,icolm)=CONJG(epsR(i,irowm,icolm))
          enddo
        enddo
      enddo
    endif
#endif
    do icolm=1,scal%npc
      do irowm=1,scal%npr
#ifdef CPLX
        eps(irowm,icolm)=epsR(1,irowm,icolm)
#else
        eps(irowm,icolm)=dble(epsR(1,irowm,icolm))
#endif
      enddo
    enddo
  endif
  
  if (peinf%inode.eq.0) then
    call close_file(iunit)
  endif
  
  if (peinf%inode.eq.0) then
    write(6,'(a,i6)')     "     omega num  = ", nFreq
    write(6,'(a,f7.3)')   "     omega step = ", dDeltaFreq
    write(6,'(a,f7.3,/)') "     omega brd  = ", dBrdning
  endif

!-----------------------------
! read rho file

  SAFE_ALLOCATE(rhog, (ng))
  rhog(:)=ZERO
  
  if (peinf%inode.eq.0) then
    write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fnrho)
    
    iunit=95
    call open_file(unit=iunit,file=trim(fnrho),form='unformatted',status='old')
    
    rho0=0.0d0

    read(iunit)
    read(iunit) nspin,nvecs
    read(iunit) (FFTgrid(i),i=1,3)
  endif
! (gsm) this is ugly but the only way I see to quickly fix the mess
! with gvec structure and gvec_index subroutine
#ifdef MPI
  call MPI_Bcast(FFTgrid,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  gvec%ng=ng
  gvec%FFTgrid=FFTgrid
  call gvec_index(gvec)
  nFFTgridpts=gvec%nFFTgridpts
  if (peinf%inode.eq.0) then
    read(iunit) celvol
    read(iunit) reccelvol, blat, ((bvec(i,j),i=1,3),j=1,3), ((bdot(i,j),i=1,3),j=1,3)
    read(iunit)
    read(iunit)
    read(iunit)
    SAFE_ALLOCATE(gvec_rho, (3, nvecs))
    SAFE_ALLOCATE(xcdum, (nvecs, nspin))
    read(iunit) nproc_para
    jg = 1
    do i=1,nproc_para
      read(iunit) num_gvec
      read(iunit) ((gvec_rho(k, j), k = 1, 3), j = jg, jg + num_gvec - 1)
      jg = jg + num_gvec
    enddo
    read(iunit) nproc_para
    jg = 1
    do i=1,nproc_para
      read(iunit) num_gvec
      read(iunit) ((xcdum(j, k), j = jg, jg + num_gvec - 1), k = 1, nspin)
      jg = jg + num_gvec
    enddo
    ierr=0
    do j=1,nvecs
      if (gvec_rho(1,j).eq.0.and.gvec_rho(2,j).eq.0.and.gvec_rho(3,j).eq.0) then
        do k=1,nspin
          rho0=rho0+dble(xcdum(j,k))
        enddo
      endif
      iout=((gvec_rho(1,j)+FFTgrid(1)/2)*FFTgrid(2)+gvec_rho(2,j)+FFTgrid(2)/2)* &
        FFTgrid(3)+gvec_rho(3,j)+FFTgrid(3)/2+1
      if (iout.ge.1.and.iout.le.nFFTgridpts) then
        iout=gvec%index_vec(iout)
        if (iout.ge.1.and.iout.le.ng) then
          rhog(iout)=ZERO
          do k=1,nspin
            rhog(iout)=rhog(iout) + xcdum(j,k)
          enddo
        else
          ierr=1
        endif
      else
        ierr=1
      endif
    enddo
    SAFE_DEALLOCATE(gvec_rho)
    SAFE_DEALLOCATE(xcdum)

    call close_file(iunit)
    
    write(6,'(5x,"cel vol =",e20.12)') celvol
    write(6,'(5x,"rho(0) =",f10.3,/)') rho0
  endif
  
#ifdef MPI
  call MPI_Bcast(bdot,9,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(celvol,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(FFTgrid,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(rho0,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(rhog,ng,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#endif
  
  if (ierr.ne.0) then
    call die("unknown G vector in file " // trim(fnrho))
  endif
  if (abs(rho0).le.TOL_Zero) then
    call die("cannot find rho(0) in file " // trim(fnrho))
  endif

!-----------------------------
! construct generalized plasmon pole model

  if (peinf%inode.eq.0) then
    write(6,'(1x,"constructing GPP model",/)')
  endif
  
  epsPP(:,:,:)=(0.0d0,0.0d0)
  wp2=ryd*ryd*16.0d0*PI_D*rho0/celvol
  sflag=.false.
  irows=0
  icols=0
  icurr=0
  do jgp=1,nmtx
    icol=mod(int(((jgp-1)/scal%nbl)+TOL_Small),scal%npcol)
    do jg=1,nmtx
      irow=mod(int(((jg-1)/scal%nbl)+TOL_Small),scal%nprow)
      if (irow.eq.scal%myprow.and.icol.eq.scal%mypcol) then
        icurr=icurr+1
        icolm=int((icurr-1)/scal%npr+TOL_Small)+1
        irowm=mod((icurr-1),scal%npr)+1
        if (jg.eq.ig.and.jgp.eq.igp) then
          sflag=.true.
          irows=irowm
          icols=icolm
        endif
        qg(:)=q(:)+dble(gvec%components(:,isrtx(jg)))
        qgp(:)=q(:)+dble(gvec%components(:,isrtx(jgp)))
        qgqg=dot_product(qg,matmul(bdot,qg))
        qgqgp=dot_product(qg,matmul(bdot,qgp))
        if (abs(qgqg).lt.TOL_Zero) cycle
        gx=gvec%components(1,isrtx(jg))-gvec%components(1,isrtx(jgp))
        gy=gvec%components(2,isrtx(jg))-gvec%components(2,isrtx(jgp))
        gz=gvec%components(3,isrtx(jg))-gvec%components(3,isrtx(jgp))
        iout=((gx+FFTgrid(1)/2)*FFTgrid(2)+gy+FFTgrid(2)/2)* &
          FFTgrid(3)+gz+FFTgrid(3)/2+1
        if (iout.lt.1.or.iout.gt.nFFTgridpts) cycle
        iout=gvec%index_vec(iout)
        if (iout.lt.1.or.iout.gt.ng) cycle
        Omega2=wp2*qgqgp/qgqg*rhog(iout)/rho0
        epsggp=eps(irowm,icolm)
        if (jg.eq.jgp) then
          delta=1.0d0
        else
          delta=0.0d0
        endif
        I_epsggp=delta-epsggp
        if (abs(I_epsggp).lt.TOL_Small) cycle
#ifdef CPLX
! Complex GPP [PRB 40, 3162 (1989)]
        wtilde2_temp = Omega2 / I_epsggp
        lambda = abs(wtilde2_temp)
        if (lambda .lt. TOL_Small) cycle
        phi = atan2(IMAG(wtilde2_temp), dble(wtilde2_temp))
        if (abs(cos(phi)) .lt. TOL_Small) cycle
        wtilde2 = lambda / cos(phi)
        Omega2 = Omega2 * CMPLX(1.0d0, -tan(phi))
#else
! Real GPP [PRB 34, 5390 (1986)]
        wtilde2 = Omega2 / I_epsggp
        if (abs(wtilde2) .lt. TOL_Small) cycle
#endif
        do i=1,nFreq
          omega=dFreqGrid(i)
          epsPP(i,irowm,icolm)=delta+Omega2/(omega**2-wtilde2-CMPLX(0D0,dBrdning))
        enddo
      endif
    enddo
  enddo
  
  if (peinf%inode.eq.0) then
    write(6,'(5x,"plasma frequency =",f10.3," eV",/)') sqrt(wp2)
  endif

!-----------------------------
! invert matrices
  
  if (peinf%inode.eq.0) then
    write(6,'(1x,"inverting matrices",/)')
  endif

#if defined USESCALAPACK
  do i=1,nFreq
    epsAux(:,:) = epsPP(i,:,:)
    call zinvert_with_scalapack(nmtx, scal, epsAux)
    epsPP(i,:,:) = epsAux(:,:)
  enddo
  if(freq_dep .ne. 0) then
    do i=1,nFreq
      epsAux(:,:) = epsR(i,:,:)
      call zinvert_with_scalapack(nmtx, scal, epsAux)
      epsR(i,:,:) = epsAux(:,:)
      epsAux(:,:) = epsA(i,:,:)
      call zinvert_with_scalapack(nmtx, scal, epsAux)
      epsA(i,:,:) = epsAux(:,:)
    enddo
  endif
#else
  do i=1,nFreq
    epsAux(:,:) = epsPP(i,:,:)
    call zinvert_serial(nmtx, epsAux)
    epsPP(i,:,:) = epsAux(:,:)
  enddo
  if(freq_dep .ne. 0) then
    do i=1,nFreq
      epsAux(:,:) = epsR(i,:,:)
      call zinvert_serial(nmtx, epsAux)
      epsR(i,:,:) = epsAux(:,:)
      epsAux(:,:) = epsA(i,:,:)
      call zinvert_serial(nmtx, epsAux)
      epsA(i,:,:) = epsAux(:,:)
    enddo
  endif
#endif

!-----------------------------
! write generalized plasmon pole file

  if (sflag) then
    write(6,'(1x,"writing",1x,a,1x,"file",/)')trim(fngpp)
    
    call open_file(unit=7,file=trim(fngpp),form='formatted',status='replace')
    
    do i=1,nFreq
      omega=dFreqGrid(i)
      write(7,100)omega,epsPP(i,irows,icols)
    enddo
    
    call close_file(7)
  endif

!-----------------------------
! write full frequency files

  if (sflag) then
    if (freq_dep.ne.0) then
      write(6,'(1x,"writing",1x,a,1x,"and",1x,a,1x,"files",/)') &
        trim(fnffr),trim(fnffa)
      
      call open_file(unit=8,file=trim(fnffr),form='formatted',status='replace')
      call open_file(unit=9,file=trim(fnffa),form='formatted',status='replace')
      
      do i=1,nFreq
        omega=dFreqGrid(i)
        write(8,100)omega,epsR(i,irows,icols)
        write(9,100)omega,epsA(i,irows,icols)
      enddo
      
      call close_file(8)
      call close_file(9)
    endif
  endif

!-----------------------------
! deallocate and finish
  
  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE(isrtx)
  SAFE_DEALLOCATE_P(gvec%index_vec)
  SAFE_DEALLOCATE(eps)
  SAFE_DEALLOCATE(epsPP)
  SAFE_DEALLOCATE(dFreqGrid)
  SAFE_DEALLOCATE(dFreqBrd)
  if (freq_dep.ne.0) then
    SAFE_DEALLOCATE(epsR)
    SAFE_DEALLOCATE(epsA)
  endif
  SAFE_DEALLOCATE(epsAux)
  SAFE_DEALLOCATE(rhog)
  
#ifdef MPI
  call MPI_Finalize(mpierr)
#endif
  
100 format(3f25.15)
  
end program epsomega
