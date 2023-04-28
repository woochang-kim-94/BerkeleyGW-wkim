!================================================================================
!
! Routines:
!
! (1) epscopy()  Originally By MLT               Last Modified: 5/5/2008 (JRD)
!
!     This routine reads in epsmat/eps0mat.
!
!     Input: crys,gvec,syms types
!            xct%ecute
!            xct%ecutg
!            is_background
!
!     Output: qg type
!             INT_EPS_* files
!
!================================================================================

#include "f_defs.h"

module epscopy_m

  use checkbz_m
  use epsread_hdf5_m
  use fullbz_m
  use global_m
  use misc_m
  implicit none

  private

  public :: &
    epscopy, &
    tddft_bz_gen

contains

subroutine epscopy(crys,gvec,syms,qg,xct,is_background)
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (symmetry), intent(in) :: syms
  type (grid), intent(out) :: qg
  type (xctinfo), intent(inout) :: xct
  logical, intent(in) :: is_background ! if true, read background dielectric matrix instead of epsmat/eps0mat

!---------------------------
! From units 10 and 11

  integer :: igamma,nrq0,nrq1,nmtx,ng,ngmax,iowner
  real(DP) :: q0(3),q0t(3,1),qk(3)
  real(DP), allocatable :: q1(:,:),eknq(:)
  SCALAR, allocatable :: epscol(:),tempepsdiag(:,:),eps(:,:)

!---------------------------
! Local stuff

  character :: ajname*6,adate*10
  character :: ajname2*6,adate2*10
  character :: tmpfn*16
  integer :: ii,jj,kk,ll,nold,gg(3),ngt,nmtxmax0,nmtxmax1,nmtxt,j,i, igp, igp_loc
  integer :: iout, irq, qgrid(3), idummy, idummy2, ig, freq_dep
  real(DP) :: gmax_in,emax,qshift(3)

  integer, allocatable :: nmtx_of_q(:),nmtx0_of_q(:)
  integer, allocatable :: isrtold(:),isrtinvdummy(:),oldg(:,:)
  real(DP), allocatable :: ekold(:)
  logical :: skip_checkbz

  character :: filenameh5*80
  character :: filenameh50*80
  character :: filename*80
  character :: filename0*80

  logical :: file_exists

  PUSH_SUB(epscopy)

  qgrid(:)=0
  SAFE_ALLOCATE(eknq, (gvec%ng))

!----------------- Read information for inverse dielectric matrix for q->0 unit10 --

  filenameh5 = 'epsmat.h5'
  filenameh50 = 'eps0mat.h5'
  filename = 'epsmat'
  filename0 = 'eps0mat'
!#BEGIN_INTERNAL_ONLY
  if(is_background) then
    filenameh5 = 'epsmat_bg.h5'
    filenameh50 = 'eps0mat_bg.h5'
    filename = 'epsmat_bg'
    filename0 = 'eps0mat_bg'
  endif
!#END_INTERNAL_ONLY

  if(peinf%inode.eq.0) then

#ifdef HDF5
    if(xct%use_hdf5) then

      call read_eps_grid_sizes_hdf5(nold,nrq0,gmax_in,idummy,idummy2,&
        nmtxmax0,qgrid,freq_dep,filenameh50)
      ng=nold
      if(is_background) then
        xct%nmtxmax_bg = nmtxmax0
        xct%ecute_bg = gmax_in
      else
        xct%nmtxmax = nmtxmax0
      endif

! XXX Check sanity below. Make sure nold doesn`t change
! XXX NEED igamma

      INQUIRE(FILE=filenameh5, EXIST=file_exists)
      if ( file_exists ) then
        igamma = 0 
      else
        igamma = 1
      endif
      
      if(igamma.eq.0) then
        call read_eps_grid_sizes_hdf5(nold,nrq1,gmax_in,idummy,idummy2,&
          nmtxmax1,qgrid,freq_dep,filenameh5)
        if(is_background) then
          if (nmtxmax1 .gt. xct%nmtxmax_bg) xct%nmtxmax_bg = nmtxmax1
        else
          if (nmtxmax1 .gt. xct%nmtxmax) xct%nmtxmax = nmtxmax1
        endif
      else
        nrq1=0
      endif
      
      ngmax=nold

    else
#endif

      call open_file(unit=10,file=filename0,form='unformatted',status='old')
      call open_file(unit=11,file=filename,form='unformatted',status='old',iostat=igamma)
    
      read(10) 
      read(10) ii
      if (ii.ne.0) call die('epscopy: freq_dep')
      read(10)
      read(10)
      read(10)
      read(10)
      read(10)
      read(10)
      read(10) nold
      read(10) ng, nmtx
      call close_file(10)
      
      if(is_background) then
        xct%nmtxmax_bg = nmtx
      else
        xct%nmtxmax = nmtx
      endif
      
      if(igamma.eq.0) then
        read(11) 
        read(11) ii
        if (ii.ne.0) call die('epscopy: freq_dep')
        read(11)
        read(11)
        read(11)
        read(11)
        read(11)
        read(11) nrq1
        read(11) nold
        ngmax= 0
        do ii=1,nrq1
          read(11) ngt, nmtxt
          read(11)
          read(11)
          do jj = 1, nmtxt
            read(11)
          enddo
          if (ngt.gt.ngmax) ngmax= ngt
          if(is_background) then
            if (nmtxt.gt.xct%nmtxmax_bg) xct%nmtxmax_bg= nmtxt
          else
            if (nmtxt.gt.xct%nmtxmax) xct%nmtxmax= nmtxt
          endif
        enddo
        call close_file(11)
      else
        nrq1=0
      endif

#ifdef HDF5
    endif
#endif

    if(ng > gvec%ng) then
      write(0,*) 'Read from epsmat ng = ', ng, ' > gvec%ng = ', gvec%ng
      call die("epscopy: read illegal ng from epsmat")
    endif

    SAFE_ALLOCATE(oldg, (3, nold))
    SAFE_ALLOCATE(isrtold, (ng))
    SAFE_ALLOCATE(ekold, (ng))
    
#ifdef HDF5
    if(xct%use_hdf5) then

      ajname='chiGG0'
      adate = 'nodate'

      SAFE_ALLOCATE(nmtx0_of_q,(nrq0))
      call read_eps_qgrid_hdf5(nrq0,q0t,nmtx0_of_q,filenameh50)
      q0(:)=q0t(:,1)

      call read_eps_old_gvecs_hdf5(nold,oldg,filenameh50)

    else
#endif

      call open_file(unit=10,file=filename0,form='unformatted',status='old')    
      
      read(10) ajname,adate
      read(10)
      read(10) (qgrid(ii),ii=1,3)
      read(10)
      read(10)
      read(10)
      read(10) gmax_in
      read(10) nrq0,(q0(ii),ii=1,3)
      read(10) nold,(oldg(1:3,ii),ii=1,nold)

      if (is_background) xct%ecute_bg=gmax_in
#ifdef HDF5
    endif
#endif

    if(nrq0.gt.1) then
      call die("There is more than one q-point in eps0mat.", only_root_writes = .true.)
    endif
    if(xct%use_hdf5) then
      write(6,'(1x,a)') 'Epsilon matrix for q->0 read from eps0mat.h5:'
    else
      write(6,'(1x,a)') 'Epsilon matrix for q->0 read from eps0mat:'
    endif
    write(6,'(1x,a,i0)') '- Number of q-points: ', nrq0
    write(6,'(1x,a,f0.2)') '- Dielectric cutoff (Ry): ', gmax_in
    if (peinf%verb_high) then
      write(6,'(1x,a,3(1x,f10.6))') '- q0 vector:', q0(1:3)
    endif
    write(6,'()')
  endif

#ifdef MPI
  call MPI_BCAST(nrq0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_BCAST(nrq1,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  if(is_background) then
    call MPI_BCAST(xct%nmtxmax_bg,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(xct%ecute_bg, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
  else
    call MPI_BCAST(xct%nmtxmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  endif
  call MPI_BCAST(igamma,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_BCAST(q0,3,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif

  ! FHJ: Block size for column distribution
  xct%nb = 1
  xct%ngpown = NUMROC(xct%nmtxmax, xct%nb, peinf%inode, 0, peinf%npes)
  xct%ngpown_max = NUMROC(xct%nmtxmax, xct%nb, 0, 0, peinf%npes)
  xct%ngpown_bg = NUMROC(xct%nmtxmax_bg, xct%nb, peinf%inode, 0, peinf%npes)
  xct%ngpown_max_bg = NUMROC(xct%nmtxmax_bg, xct%nb, 0, 0, peinf%npes)

  if (is_background) then
    SAFE_ALLOCATE(xct%nmtxa_bg, (nrq1+1))
  else
    SAFE_ALLOCATE(xct%nmtxa, (nrq1+1))
  endif
  if (is_background) then
    SAFE_ALLOCATE(xct%isrtqi_bg, (gvec%ng,nrq1+1))
  elseif (peinf%inode .eq. 0 .or. xct%bLowComm) then
    SAFE_ALLOCATE(xct%isrtqi, (gvec%ng,nrq1+1))
  endif
  if (is_background) then
    SAFE_ALLOCATE(xct%epscol_bg, (xct%nmtxmax_bg,xct%nmtxmax_bg,nrq1+1))
  elseif (xct%bLowComm) then
    SAFE_ALLOCATE(xct%epscol, (xct%nmtxmax,xct%nmtxmax,nrq1+1))
  else 
    SAFE_ALLOCATE(xct%epscol, (xct%nmtxmax,xct%ngpown_max,nrq1+1))
  endif

  if (is_background) then
    SAFE_ALLOCATE(xct%epsdiag_bg, (xct%nmtxmax_bg,nrq1+1))
  else
    SAFE_ALLOCATE(xct%epsdiag, (xct%nmtxmax,nrq1+1))
  endif
! Read q->0 dielectric matrix

  if(peinf%inode.eq.0) then

#ifdef HDF5
    if(xct%use_hdf5) then

      SAFE_ALLOCATE(isrtinvdummy, (ng))
      call read_eps_gvecsofq_hdf5(ng,isrtold,isrtinvdummy,ekold,1,filenameh50)
      SAFE_DEALLOCATE(isrtinvdummy)
      nmtx = nmtx0_of_q(1)
      qk = q0

    else
#endif

      read(10) ng,nmtx,(isrtold(ii),jj,ii=1,ng)
      read(10) (ekold(ii),ii=1,ng)
      read(10) (qk(ii),ii=1,3)

#ifdef HDF5
    endif
#endif

    if (is_background) then
      xct%isrtqi_bg(:,1)= 0
    else
      xct%isrtqi(:,1)= 0
    endif

!---------------------
! Sort the eps. matrix elements according to gvec%.
! Emax is some large energy, bigger than xct%ecute (but it does not
! need to be as large as the ekmax used to write epsmat/eps0mat).
! Check if the value of emax is OK.

    emax=xct%ecutg

    do ii=1,ng
      if (ekold(isrtold(ii)).lt.emax) then
        gg(1:3)=oldg(1:3, isrtold(ii))
        call findvector(iout,gg,gvec)
        if (iout.gt.gvec%ng) call die('epscopy: iout > ng')
        if (iout.le.0) call die('epscopy: iout <= 0')

!--------------------------------
! isrtqi has the sorting of G-vectors from eps to gvec%components

        if (peinf%inode .eq. 0) then
          if (is_background) then
            xct%isrtqi_bg(iout,1)=ii
          else
            xct%isrtqi(iout,1)=ii
          endif
        endif
      endif
    enddo
    SAFE_DEALLOCATE(isrtold)
    SAFE_DEALLOCATE(ekold)
    SAFE_DEALLOCATE(oldg)
  endif

  xct%q0vec=q0
  q0 = 0.0d0

! JRD: Write header of INT_EPS_*

  irq=0
  
#ifdef MPI
  call MPI_BCAST(nmtx,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif

  if (is_background) then
    xct%nmtxa_bg(1)=nmtx
  else
    xct%nmtxa(1)=nmtx
  endif
!------------------------------
! JRD: Finally Read In Eps 

#ifdef HDF5
  if(xct%use_hdf5) then

    if (xct%bLowComm .or. is_background) then
! JRD: XXX This is slower than need be. Should use parallel read + MPI_ALL_ALL
      SAFE_ALLOCATE(eps,(nmtx,nmtx))
      if (peinf%inode .eq. 0) then
        call read_eps_matrix_ser_hdf5(eps,nmtx,1,1,filenameh50)
      endif
#ifdef MPI
      call MPI_BCAST(eps(1,1),nmtx*nmtx,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#endif
      if (is_background) then
        xct%epscol_bg(1:nmtx,1:nmtx,1) = eps(1:nmtx,1:nmtx)
      else
        xct%epscol(1:nmtx,1:nmtx,1) = eps(1:nmtx,1:nmtx)
      endif
      SAFE_DEALLOCATE(eps)
    else
      call read_eps_matrix_par_hdf5(xct%epscol(:,:,1), &
        xct%nb, peinf%inode, peinf%npes, nmtx, 1, 1, filenameh50)
    endif

! XXX  Print Constant

  else
#endif

    SAFE_ALLOCATE(epscol, (nmtx))
    do igp =1,nmtx
      if (peinf%inode .eq. 0) then
        read(10) (epscol(ii),ii=1,nmtx)

        ! Report on dielectric constant
        if (igp .eq. 1) then
          write(6,'()')
          write(6,*) 'Head of epsilon inverse : ', epscol(1)
          write(6,'()')
          
          if(dble(epscol(1))<1d-3 .and. xct%iscreen==SCREEN_SEMICOND .and. peinf%inode==0) then
            write(0,'(a)') 'WARNING: You are using semiconductor screening, but the'
            write(0,'(a)') 'head of epsilon inverse is very small and seems metallic.'
          endif
        endif
      endif

      ! Write dielectric matrix column for q->0 to unit17 unformatted
#ifdef MPI
      call MPI_BCAST(epscol,nmtx,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#endif
      igp_loc = INDXG2L(igp, xct%nb, peinf%inode, 0, peinf%npes)
      iowner = INDXG2P(igp, xct%nb, peinf%inode, 0, peinf%npes)
      if (is_background) then
        xct%epscol_bg(1:nmtx,igp,1) = epscol(:)
      elseif (xct%bLowComm ) then
        xct%epscol(1:nmtx,igp,1) = epscol(:)
      else
        if (iowner==peinf%inode) then
          xct%epscol(1:nmtx,igp_loc,1) = epscol(:)
        endif
      endif

      if (peinf%inode .eq. 0) then
        if (is_background) then
          xct%epsdiag_bg(igp,1) = epscol(igp)  
        else
          xct%epsdiag(igp,1) = epscol(igp)  
        endif
      endif
    enddo
    SAFE_DEALLOCATE(epscol)
    
    if (peinf%inode .eq. 0) then
      call close_file(10)
    endif
    
#ifdef HDF5
  endif
#endif

!----------------- Read dielectric matrices from unit11 for q<>0 --------------------

  if(igamma.ne.0) then
    nrq1=0
  else

! Have to allocate oldg again...

    if(peinf%inode.eq.0) then
      SAFE_ALLOCATE(oldg, (3, nold))
      SAFE_ALLOCATE(isrtold, (ngmax))
      SAFE_ALLOCATE(ekold, (ngmax))
      SAFE_ALLOCATE(q1, (3,nrq1))
#ifdef HDF5
      if(xct%use_hdf5) then
        SAFE_ALLOCATE(nmtx_of_q,(nrq1))
        call read_eps_qgrid_hdf5(nrq1,q1,nmtx_of_q,filenameh5)
        
        call read_eps_old_gvecs_hdf5(nold,oldg,filenameh5)
      else
#endif
        call open_file(unit=11,file=filename,form='unformatted',status='old')
        read(11) ajname2,adate2
        read(11)
        read(11) (qgrid(ii),ii=1,3)
        read(11)
        read(11)
        read(11)
        read(11) gmax_in
        read(11) nrq1,((q1(ii,jj),ii=1,3),jj=1,nrq1)
        read(11) nold,(oldg(1:3,ii),ii=1,nold)
#ifdef HDF5
      endif
#endif

    endif

#ifdef MPI
    if(peinf%inode.ne.0) then
      SAFE_ALLOCATE(q1, (3,nrq1))
    endif
    call MPI_BCAST(q1,3*nrq1,  MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif

    if(peinf%inode.eq.0) then
      if(xct%use_hdf5) then
        write(6,'(1x,a)') 'Epsilon matrix for q/=0 read from epsmat.h5:'
      else
        write(6,'(1x,a)') 'Epsilon matrix for q/=0 read from epsmat:'
      endif
      write(6,'(1x,a,i0)') '- Number of q-points: ', nrq1
      write(6,'(1x,a,f0.2)') '- Dielectric cutoff (Ry): ', gmax_in
      if (peinf%verb_high) then
        write(6,'(1x,a)') '- Q-points:'
        write(6,'(1(2x,3(1x,f10.6)))') q1(1:3,1:nrq1)
      endif
    endif
    
  endif
  if (is_background) then
    if (nrq1+1 .ne. qg%nr) call die('Background dielectric matrix must have the same number of qpoints as epsmat')
  else
    qg%nr=nrq1+1
  endif
  if (.not.is_background) then
    SAFE_ALLOCATE(qg%r, (3,qg%nr))
    qg%r(1:3,1)=q0
    if(nrq1.ne.0) then
      qg%r(1:3,2:qg%nr)=q1(1:3,1:nrq1)
      SAFE_DEALLOCATE(q1)
    endif
  else
    if (.not. all(qg%r(1:3,2:qg%nr).eq.q1(1:3,1:nrq1))) then
      call die('Background dielectric matrix must have the same qpoints as epsmat')
    endif
  endif


! Read inverse dielectric matrices from unit11 for q<>0

  if(igamma == 0) then

    do irq=1,nrq1
      if(peinf%inode.eq.0) then
        isrtold=0
        ekold=0.d0
        if (is_background) then
          xct%isrtqi_bg(:,irq+1)=0
        else
          xct%isrtqi(:,irq+1)=0
        endif
#ifdef HDF5
        if(xct%use_hdf5) then
          SAFE_ALLOCATE(isrtinvdummy, (nold))
          call read_eps_gvecsofq_hdf5(nold,isrtold,isrtinvdummy,ekold,irq,filenameh5)
          SAFE_DEALLOCATE(isrtinvdummy)
          nmtx = nmtx_of_q(irq)
          qk(:) = qg%r(:,irq+1)
        else
#endif
          read(11) ng,nmtx,(isrtold(ii),jj,ii=1,ng)
          read(11) (ekold(ii),ii=1,ng)
          read(11) (qk(ii),ii=1,3)
#ifdef HDF5
        endif
#endif
        
        do ii=1,ng
          if (ekold(isrtold(ii)).lt.emax) then
            gg(1:3)=oldg(1:3, isrtold(ii))
            call findvector(iout,gg,gvec)
            if (iout.gt.gvec%ng) call die('epscopy: iout > ng')
            if (iout.le.0) call die('epscopy: iout <= 0')
            if (is_background) then
              xct%isrtqi_bg(iout,irq+1)=ii
            else
              xct%isrtqi(iout,irq+1)=ii
            endif
            eknq(ii)=ekold(isrtold(ii))
          endif
        enddo
      endif
#ifdef MPI
      call MPI_BCAST(nmtx,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
      
      if (is_background) then
        xct%nmtxa_bg(irq+1)=nmtx
      else
        xct%nmtxa(irq+1)=nmtx
      endif
#ifdef HDF5
      if(xct%use_hdf5) then

        
        if (xct%bLowComm .or. is_background) then
! JRD: XXX This is slower than need be. Should use parallel read + MPI_ALL_ALL
          SAFE_ALLOCATE(eps,(nmtx,nmtx))
          if (peinf%inode .eq. 0) then
            call read_eps_matrix_ser_hdf5(eps,nmtx,irq,1,filenameh5)
          endif
#ifdef MPI
          call MPI_BCAST(eps(1,1),nmtx*nmtx,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#endif
          if (is_background) then
            xct%epscol_bg(1:nmtx,1:nmtx,irq+1) = eps(1:nmtx,1:nmtx)
          else
            xct%epscol(1:nmtx,1:nmtx,irq+1) = eps(1:nmtx,1:nmtx)
          endif
          SAFE_DEALLOCATE(eps)
        else
          call read_eps_matrix_par_hdf5(xct%epscol(:,:,irq+1), &
            xct%nb, peinf%inode, peinf%npes, nmtx, irq, 1, filenameh5)
        endif

      else
#endif
        SAFE_ALLOCATE(epscol, (nmtx))

        do igp=1,nmtx
          if (peinf%inode .eq. 0) then
            read(11) (epscol(ii),ii=1,nmtx)
            if (is_background) then
              xct%epsdiag_bg(igp,irq+1)=epscol(igp)
            else
              xct%epsdiag(igp,irq+1)=epscol(igp)
            endif
          endif
#ifdef MPI
          call MPI_BCAST(epscol,nmtx,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#endif
          igp_loc = INDXG2L(igp, xct%nb, peinf%inode, 0, peinf%npes)
          iowner = INDXG2P(igp, xct%nb, peinf%inode, 0, peinf%npes)
          if (is_background) then
            xct%epscol_bg(1:nmtx,igp,irq+1) = epscol(:)
          elseif (xct%bLowComm) then
            xct%epscol(1:nmtx,igp,irq+1) = epscol(:)
          else
            if (iowner==peinf%inode) then
              xct%epscol(1:nmtx,igp_loc,irq+1) = epscol(:)
            endif
          endif

        enddo

        SAFE_DEALLOCATE(epscol)
#ifdef HDF5
      endif
#endif

    enddo
    
    if (peinf%inode .eq. 0) then
      SAFE_DEALLOCATE(isrtold)
      SAFE_DEALLOCATE(ekold)
      SAFE_DEALLOCATE(oldg)
      if(.not. xct%use_hdf5) call close_file(11)
    endif
    
  endif
  
  SAFE_DEALLOCATE(eknq)

!-------------------------
! Generate full brillouin zone from irreducible wedge
! rq -> fq

  if (.not. is_background) then
    call timacc(7,1)
    call fullbz(crys,syms,qg,syms%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
    qshift(:)=0.0d0
    if (igamma.ne.0) then
      tmpfn='eps0mat'
    else
      tmpfn='epsmat'
    endif
    if (.not. skip_checkbz) then
      call checkbz(qg%nf,qg%f,qgrid,qshift,crys%bdot,tmpfn,'q',.true.,xct%freplacebz,xct%fwritebz)
    endif
    call timacc(7,2)
    
    if (peinf%verb_high .and. peinf%inode==0) then
      write(6,'(1x,a)') '- Unfolded q-points:'
      do ii=1,qg%nf
        write(6,'(2x,3(1x,f10.6))') qg%f(:,ii)
      enddo
    endif
  endif
  
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD,mpierr)
  ! Doing the communication once here. Saves us a communication each time mtxel_kernel
  ! is called
  if (is_background) then
    call MPI_BCAST(xct%isrtqi_bg,gvec%ng*(nrq1+1),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  elseif (xct%bLowComm) then
    call MPI_BCAST(xct%isrtqi,gvec%ng*(nrq1+1),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  endif

  if (xct%use_hdf5 .and. .not. xct%bLowComm .and. .not. is_background) then
! XXX This should be ALL to ALL not ALL REDUCE
! XXX IT SHOULD ACTUALLY USE THE NEW matrix_diagonal ROUTINE INSTEAD
    xct%epsdiag=0D0
    SAFE_ALLOCATE(tempepsdiag, (xct%nmtxmax,nrq1+1))
    do igp_loc = 1, xct%ngpown
      igp = INDXL2G(igp_loc, xct%nb, peinf%inode, 0, peinf%npes)
      do irq = 1, nrq1+1
        tempepsdiag(igp, irq) = xct%epscol(igp, igp_loc, irq)
      enddo
    enddo
    call MPI_ALLREDUCE(tempepsdiag(1,1),xct%epsdiag(1,1),xct%nmtxmax*(nrq1+1),MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
    SAFE_DEALLOCATE(tempepsdiag)
  elseif (is_background) then
    call MPI_BCAST(xct%epsdiag_bg,xct%nmtxmax*(nrq1+1),MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
  else
    call MPI_BCAST(xct%epsdiag,xct%nmtxmax*(nrq1+1),MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
  endif
#endif
  
  if (xct%use_hdf5 .and. peinf%inode .eq. 0) then
    SAFE_DEALLOCATE(nmtx0_of_q)
    if (igamma.eq.0) then
      SAFE_DEALLOCATE(nmtx_of_q)
    endif
  endif
  
  POP_SUB(epscopy)
  
  return
end subroutine epscopy

subroutine tddft_bz_gen(crys,syms,qg,xct)
  type (crystal), intent(in) :: crys
  type (symmetry), intent(in) :: syms
  type (grid), intent(out) :: qg
  type (xctinfo), intent(inout) :: xct

  real(DP) :: qshift(3)
  character :: tmpfn*16
  logical :: skip_checkbz
  integer :: ii
!-------------------------
! Generate full brillouin zone from irreducible wedge
! rq -> fq

  PUSH_SUB(tddft_bz_gen)
  call timacc(7,1)
  call fullbz(crys,syms,qg,syms%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  qshift(:)=0.0d0
  tmpfn='tddft_q'

  if (.not. skip_checkbz) then
    call checkbz(qg%nf,qg%f,xct%qgrid,qshift,crys%bdot,tmpfn,'q',.true.,xct%freplacebz,xct%fwritebz)
  endif
  call timacc(7,2)
  
  if (peinf%verb_high .and. peinf%inode==0) then
    write(6,'(/,1x,a)') 'Unfolded q-points:'
    do ii=1,qg%nf
      write(6,'(2x,3(1x,f10.6))') qg%f(:,ii)
    enddo
  endif
  POP_SUB(tddft_bz_gen)

  return
end subroutine tddft_bz_gen

end module epscopy_m
