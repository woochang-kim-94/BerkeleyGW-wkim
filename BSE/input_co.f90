!==========================================================================
!
! Routines:
!
! (1) input_co()        Originally By ?         Last Modified 4/19/2009 (gsm)
!
!     input: crys, gvec,  syms, xct, flagbz types
!
!     output: kg, distgwfco types
!
!     Reads in the coarse-grid wavefunctions from file WFN_co and
!     distributes them between processors. The k-point grid is stored in kg.
!
!==========================================================================

#include "f_defs.h"

module input_co_m

  use checkbz_m
  use eqpcor_m
  use fullbz_m
  use global_m
#ifdef HDF5
  use hdf5
#endif
  use input_utils_m
  use io_utils_m
  use misc_m
  use scissors_m
  use wfn_rho_vxc_io_m
  use wfn_io_hdf5_m
  implicit none

  private

  public :: &
    input_co

contains

subroutine input_co(kp,kp_co,crys,gvec,kg,kgq,syms,xct,flagbz,distgwfco,eqp)
  type (kpoints), intent(in) :: kp
  type (kpoints), intent(out) :: kp_co
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (grid), intent(out) :: kg,kgq
  type (symmetry), intent(in) :: syms
  type (xctinfo), intent(inout) :: xct
  integer, intent(in) :: flagbz
  type (tdistgwf), intent(out) :: distgwfco
  type (eqpinfo), intent(inout) :: eqp

  type (crystal) :: crys_co
  type (symmetry) :: syms_co
  type (wavefunction) :: wfnv,wfnc
  type (kpoints) :: kpq
  character :: filenamev*20,filenamec*20
  character :: tmpfn*16
  character :: fncor*32
  integer :: iunit_v,iunit_c
  integer :: irk,irks,ik,ikq,umk
  integer :: ii,jj,kk,is,isp,minband,maxband
  real(DP) :: kt(3),div,tol,qq_temp(3),delta
  integer, allocatable :: indxk(:)
  SCALAR, allocatable :: cg(:,:), cgarray(:)

  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvec_co, gvec_kpt
  logical :: skip_checkbz, broken_degeneracy, dont_read
  logical :: dont_die_consistency
  type(progress_info) :: prog_info

  ! WFN HDF5 stuff
#ifdef HDF5
  integer(HID_T) :: file_id
  integer(HID_T) :: plist_id
#endif
  integer, allocatable :: gvec_kpt_components_all(:,:)
  integer, allocatable :: ib_size_array(:)
  integer :: istart, ngktot
  integer :: iband_min, iband_max, tot_nbands, band_block, min_band_block
  integer :: ib_first, ib_last, ib_size, ib_size_max
  integer :: ipe, npes_hdf5, ii_start, ii_end, ii_loc
  SCALAR, allocatable :: my_cg_all(:,:,:), ipe_cg_all(:,:,:)
  integer :: error

  PUSH_SUB(input_co)

!-------------------------
! Print to stdout

  sheader = 'WFN'
  iflavor = 0
  if ( xct%use_wfn_hdf5 ) then
    if (peinf%inode==0) write(6,'(a)') ' Reading header of WFN_co.h5'
#ifdef HDF5
    call read_hdf5_header_type('WFN_co.h5', sheader, iflavor, kp_co, gvec_co, syms_co, crys_co)
#endif
  else
    if (peinf%inode.eq.0) call open_file(unit=25,file='WFN_co',form='unformatted',status='old')
    call read_binary_header_type(25, sheader, iflavor, kp_co, gvec_co, syms_co, crys_co, warn = .false., &
      dont_warn_kgrid=.true.)
  end if
  call check_trunc_kpts(xct%icutv, kp_co)

  call check_header('WFN_fi', kp, gvec, syms, crys, 'WFN_co', kp_co, gvec_co, syms_co, crys_co, is_wfn = .true.)

  if(xct%skipinterp) then
    ! we check this first to be sure the second comparison will not segfault
    if(kp%nrk /= kp_co%nrk .or. any(kp%kgrid(1:3) /= kp_co%kgrid(1:3))) then
      call die("Cannot skip interpolation if coarse and fine grids differ.")
    endif
    if(any(abs(kp%rk(1:3, 1:kp%nrk) - kp_co%rk(1:3, 1:kp%nrk)) > TOL_Zero)) then
      call die("Cannot skip interpolation if coarse and fine k-points differ.")
    endif
  endif

  SAFE_ALLOCATE(gvec_co%components, (3, gvec%ng))
  if( xct%use_wfn_hdf5 ) then
#ifdef HDF5
    call read_hdf5_gvectors('WFN_co.h5', gvec%ng, gvec_co%components)
#endif
  else
    call read_binary_gvectors(25, gvec%ng, gvec%ng, gvec_co%components)
  end if
  SAFE_DEALLOCATE_P(gvec_co%components)

!-----------------------------------------------------------------------
! Read eqp_co.dat for possible interpolation

  SAFE_ALLOCATE(eqp%evshift_co, (xct%nvb_co,kp_co%nrk,kp_co%nspin))
  eqp%evshift_co=0D0
  SAFE_ALLOCATE(eqp%ecshift_co, (xct%ncb_co,kp_co%nrk,kp_co%nspin))
  eqp%ecshift_co=0D0
  
  fncor = ''

  if(xct%eqp_corrections .and. xct%skipinterp) fncor = 'eqp.dat'
  ! not interpolating, coarse and fine grids are identical.
  ! must correct this grid or wreck occupations. we need kp_co%el for efermi below.
  ! DAS: this is a hack to make things work. Better is never to read anything of the coarse grid if not interpolating!

  ! Note: eqp_co_corrections and skipinterp are incompatible and blocked in inread
  if(xct%eqp_co_corrections) then
    fncor = 'eqp_co.dat'
    xct%inteqp = .true.
  endif

  SAFE_ALLOCATE(kp_co%elda, (kp_co%mnband, kp_co%nrk, kp_co%nspin))
  kp_co%el(:,:,:) = kp_co%el(:,:,:) - xct%avgpot / ryd
  kp_co%elda(:,:,:) = kp_co%el(:,:,:)

  minband = 1
  maxband = kp_co%mnband
  if(trim(fncor) /= '') then
    ! FIXME: for metals this is asking for a few more bands than actually needed on some k-points
    minband = minval(kp_co%ifmax(:,:)-xct%nvb_co+1)
    maxband = maxval(kp_co%ifmax(:,:)+xct%ncb_co)
    call eqpcor(fncor,peinf%inode,peinf%npes,kp_co, &
      minband,maxband,0,0,kp_co%el,eqp%evshift_co,eqp%ecshift_co,1,0)
  endif

  if(xct%eqp_co_corrections .and. xct%eqp_corrections) xct%inteqp = .false.
  ! if we have the fine-grid QP energies, we do not need to interpolate from the fine grid

  ! scissor shift is only needed for consistency with Fermi level shift here
  ! since 'interpolating' to fine grid is the same as just applying scissor shift directly
  call scissors_shift(kp_co, eqp%scis, eqp%spl_tck)

  ! FHJ: We don`t allow inconsistent ifmin/ifmax fields. The only exception is
  ! if we are running inteqp with unrestricted interpolation and we are not 
  ! changing the FE. In this case, we freeze ifmin/ifmax from LDA, which is
  ! fine b/c unrestricted_transf doesn`t distinguish val/cond states.
  dont_die_consistency = xct%inteqp.and.xct%unrestricted_transf.and.&
    xct%rfermi.and.(dabs(xct%efermi_input)<TOL_SMALL)
  ! JRD: If we are using eqp_co.dat, then this is our best estimate of the true fermi qp fermi
  ! energy, so we update it. If we don`t update it, we may also mess up the occupations of the
  ! coarse grid because qp efermi of coarse grid may not be lda efermi of fine grids.
  if(xct%eqp_co_corrections) then
    call find_efermi(xct%rfermi, xct%efermi, xct%efermi_input, kp_co, maxband, minband, &
      "coarse grid", should_search = .true., should_update = .true., write7 = .false., &
      dont_die_consistency=dont_die_consistency)
  else
    call find_efermi(xct%rfermi, xct%efermi, xct%efermi_input, kp_co, maxband, minband, &
      "coarse grid", should_search = .true., should_update = .false., write7 = .false., &
      dont_die_consistency=dont_die_consistency)
  endif

  ! now we call again to initialize the eqp arrays
  if(xct%eqp_co_corrections) then
    call eqpcor(fncor,peinf%inode,peinf%npes,kp_co,0,0, &
      xct%nvb_co,xct%ncb_co,kp_co%el,eqp%evshift_co,eqp%ecshift_co,1,2,dont_write=.true.)
  endif

  if(any(kp_co%ifmax(:,:) == 0)) & 
    call die("BSE codes cannot handle a system where some k-points have no occupied bands.", only_root_writes = .true.) 

  kp_co%nvband=minval(kp_co%ifmax(:,:)-kp_co%ifmin(:,:))+1
  kp_co%ncband=kp_co%mnband-maxval(kp_co%ifmax(:,:))

!----------------------------------------------------------------
! (gsm) check whether the requested number of bands
!       is available in the wavefunction file

  if(xct%nvb_co .gt. kp_co%nvband) then
    call die("The requested number of valence bands is not available in WFN_co.", only_root_writes = .true.)
  endif
  if(xct%ncb_co .gt. kp_co%ncband) then
    call die("The requested number of conduction bands is not available in WFN_co.", only_root_writes = .true.)
  endif

! DAS: degenerate subspace check

  if (peinf%inode.eq.0) then
    ! MHN
    !if(xct%ncb_co .eq. kp_co%ncband) then
    !  call die("You must provide one more conduction band in WFN_co in order to assess degeneracy.")
    !endif
    broken_degeneracy = .false.
    !do jj = 1, kp_co%nspin
    !  do ii = 1, kp_co%nrk
    !    if(kp_co%ifmax(ii, jj) - xct%nvb_co > 0) then
    !      ! no need to compare against band 0 if all valence are included
    !      if(abs(kp_co%elda(kp_co%ifmax(ii, jj) - xct%nvb_co + 1, ii, jj) &
    !        - kp_co%elda(kp_co%ifmax(ii, jj) - xct%nvb_co, ii, jj)) .lt. TOL_Degeneracy) then
    !        broken_degeneracy = .true.
    !      endif
    !    endif
    !  enddo
    !enddo

    if(broken_degeneracy) then
      if(xct%degeneracy_check_override) then
        write(0,'(a)') &
          "WARNING: Selected number of valence bands breaks degenerate subspace in WFN_co. " // &
          "Run degeneracy_check.x for allowable numbers."
        write(0,*)
      else
        write(0,'(a)') &
          "Run degeneracy_check.x for allowable numbers, or use keyword " // &
          "degeneracy_check_override to run anyway (at your peril!)."
        call die("Selected number of valence bands breaks degenerate subspace in WFN_co.")
      endif
    endif
    
    broken_degeneracy = .false.
    do jj = 1, kp_co%nspin
      do ii = 1, kp_co%nrk
        if(abs(kp_co%elda(kp_co%ifmax(ii, jj) + xct%ncb_co, ii, jj) &
          - kp_co%elda(kp_co%ifmax(ii, jj) + xct%ncb_co + 1, ii, jj)) .lt. TOL_Degeneracy) then
          broken_degeneracy = .true.
        endif
      enddo
    enddo

    if(broken_degeneracy) then
      if(xct%degeneracy_check_override) then
        write(0,'(a)') &
          "WARNING: Selected number of conduction bands breaks degenerate subspace in WFN_co. " // &
          "Run degeneracy_check.x for allowable numbers."
        write(0,*)
      else
        write(0,'(a)') &
          "Run degeneracy_check.x for allowable numbers, or use keyword " // &
          "degeneracy_check_override to run anyway (at your peril!)."
        call die("Selected number of conduction bands breaks degenerate subspace in WFN_co.")
      endif
    endif
  endif

  SAFE_DEALLOCATE(kp_co%elda)

!-----------------------------------------------------------------------
!     Read k-points from file kpoints_co (if it exists) or from WFN_co
!     Array indxk has the same meaning as in input

  if (xct%read_kpoints) then
    if (peinf%inode.eq.0) then
      call open_file(9,file='kpoints_co',form='formatted',status='old')
      read(9,*) kg%nr
      SAFE_ALLOCATE(kg%r, (3,kg%nr))
      do ii=1,kg%nr
        read(9,*) (kg%r(jj,ii),jj=1,3),div
        kg%r(:,ii) = kg%r(:,ii)/div
      enddo
      call close_file(9)
    endif ! node 0
#ifdef MPI
    call MPI_BCAST(kg%nr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if (peinf%inode.ne.0) then
      SAFE_ALLOCATE(kg%r, (3,kg%nr))
    endif
    call MPI_BCAST(kg%r,3*kg%nr,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif
    tol = TOL_Small
    SAFE_ALLOCATE(indxk, (kg%nr))
    indxk=0
    do jj=1,kg%nr
      do ii=1,kp_co%nrk
        kt(:) = kg%r(:,jj) - kp_co%rk(:,ii)
        if (all(abs(kt(1:3)).lt.tol)) then
          if (indxk(jj).ne.0) then
            if (peinf%inode.eq.0) write(6,996) jj,indxk(jj),kg%r(:,jj)
          endif
          indxk(jj)=ii
        endif
      enddo
      if (indxk(jj).eq.0) then
        if (peinf%inode.eq.0) write(0,995) kg%r(:,jj)
      endif
    enddo
  else
    kg%nr=kp_co%nrk
    SAFE_ALLOCATE(kg%r, (3,kg%nr))
    kg%r(1:3,1:kg%nr)=kp_co%rk(1:3,1:kp_co%nrk)
    SAFE_ALLOCATE(indxk, (kg%nr))
    do ii=1,kg%nr
      indxk(ii)=ii
    enddo
  endif
#ifdef MPI
  call MPI_BCAST(indxk,kg%nr,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
996 format(1x,'WARNING: Multiple definition of k-point',2i4,3f10.6)
995 format(1x,'WARNING: Could not find k-point',3f10.6,1x,'in WFN_co')
  

!-----------------------------------------------------------------------
!     Generate full Brillouin zone from irreducible wedge, rk -> fk

  if (flagbz.eq.1) then
    call fullbz(crys,syms,kg,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  else
    call fullbz(crys,syms,kg,syms%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  endif
  tmpfn='WFN_co'
  if (.not. skip_checkbz .and. .not.xct%patched_sampling) then
    call checkbz(kg%nf,kg%f,kp_co%kgrid,kp_co%shift,crys%bdot, &
      tmpfn,'k',.true.,xct%freplacebz,xct%fwritebz)
  endif
  
  if (flagbz.eq.0.and.peinf%inode.eq.0) write(6,801)
  if (flagbz.eq.1.and.peinf%inode.eq.0) write(6,802)
801 format(1x,'Using symmetries to expand the coarse-grid sampling')
802 format(1x,'No symmetries used in the coarse-grid sampling')

  if (xct%nkpt_co.ne.kg%nf) then
   if(peinf%inode == 0) write(0,994) xct%nkpt_co,kg%nf
994 format('The given number of points in the coarse grid (',i4, &
      ') does not match the number of points in file WFN_co after unfolding (',i4,').')
    call die('If you are sure WFN_co is correct, please change the .inp file and try again.', only_root_writes = .true.)
  endif

!------------------------------------------------------------------------
! If there is a finite center-of-mass momentum, Q, find mapping between k
! and k+Q

  SAFE_ALLOCATE(xct%indexq,(kg%nf))
  if (xct%qflag.eq.1) then
    do ik=1,kg%nf
      xct%indexq(ik) = ik
    enddo
  endif

  SAFE_ALLOCATE(kgq%f,(3,kg%nf))
  SAFE_ALLOCATE(kgq%kg0,(3,kg%nf))
  SAFE_ALLOCATE(kgq%indr,(kg%nf))
! DYQ: Use a center-of-mass momentum commensurate with the kgrid of WFN_co, when qflag=2
! Find the index of the point in kg corresponding to kg%f+xct%finiteq 
! and store the index mapping in indexq. kgq is the grid kg shifted by Q.
  if (xct%qflag .eq. 2) then
    do ik=1,kg%nf
      kgq%f(:,ik) = kg%f(:,ik) + xct%finiteq(:) !k+Q 
      ikq = 0
      delta=0.1d0
      do while ((delta .gt. TOL_Small) .and. (ikq.lt.kg%nf))
        ikq = ikq+1
        qq_temp(:) = kgq%f(:,ik) - kg%f(:,ikq)
        do jj=1,3
          qq_temp(jj) = qq_temp(jj) - anint( qq_temp(jj) )
        enddo
        delta=sqrt((qq_temp(1))**2+(qq_temp(2))**2+(qq_temp(3))**2)
      enddo
      ! With patched sampling, k+Q might fall outside the patch 
      if (delta.gt.TOL_Small .and. xct%patched_sampling_co) then
        if(peinf%inode.eq.0) then
          write(6,*) '  Could not find point equivalent to ', (kg%f(ii,ik),ii=1,3)
          write(6,*) '  Skipping point.'
        endif
        xct%indexq(ik)= 0
        kgq%indr(ik) = ik
      elseif (delta.gt.TOL_Small) then
        if(peinf%inode.eq.0) then
          write(0,*) '  Could not find point equivalent to ', (kg%f(ii,ik),ii=1,3)
        endif
        call die('Finite momentum not commensurate with kgrid of WFN_co',only_root_writes = .true.)
      else
        xct%indexq(ik)=ikq
! kg%f and kgq%f may differ by a lattice vector near the zone edge
        do jj=1,3
          umk = nint( kgq%f(jj,ik)-kg%f(jj,ikq))
          kgq%kg0(jj,ik) = kg%kg0(jj,ikq) + umk
        enddo
        kgq%indr(ik)=kg%indr(ikq)
      endif
    enddo
  endif !qflag=2

!-----------------------------------------------------------------------
! Initialization of distributed wavefunctions

  distgwfco%nk=kg%nr
  distgwfco%ngm=kp_co%ngkmax
  distgwfco%ns=kp_co%nspin
  distgwfco%nspinor=kp_co%nspinor
  distgwfco%nv=xct%nvb_co
  distgwfco%nc=xct%ncb_co
    
  ! FHJ: Use standard BLACS distribution for G-vectors
  distgwfco%block_sz = DIVUP(distgwfco%ngm, peinf%npes)
  ! ngl = local number of G-vectors that I own.
  distgwfco%ngl = NUMROC(distgwfco%ngm, distgwfco%block_sz, peinf%inode, 0, peinf%npes)
  ! Local to global index translation: ig_g = ig_l + tgl
  distgwfco%tgl = distgwfco%block_sz * peinf%inode
    
  SAFE_ALLOCATE(distgwfco%ng, (distgwfco%nk))
  SAFE_ALLOCATE(distgwfco%isort, (distgwfco%ngl,distgwfco%nk))
  if (xct%qflag.ne.0) then
    SAFE_ALLOCATE(distgwfco%zv, (distgwfco%ngl,distgwfco%nv,distgwfco%ns*distgwfco%nspinor,distgwfco%nk))
  endif
  SAFE_ALLOCATE(distgwfco%zc, (distgwfco%ngl,distgwfco%nc,distgwfco%ns*distgwfco%nspinor,distgwfco%nk))
   
  distgwfco%ng(:)=0
  distgwfco%isort(:,:)=0
  if (xct%qflag.ne.0) then
    distgwfco%zv(:,:,:,:)=ZERO
  endif
  distgwfco%zc(:,:,:,:)=ZERO

!-----------------------------------------------------------------------
! Read the wavefunctions and distribute 

  SAFE_ALLOCATE(wfnv%isort, (gvec%ng))
  wfnv%nband=xct%nvb_co
  wfnv%nspin=kp_co%nspin
  wfnv%nspinor=kp_co%nspinor
  wfnc%nband=xct%ncb_co
  wfnc%nspin=kp_co%nspin
  wfnc%nspinor=kp_co%nspinor
  
  if ( xct%use_wfn_hdf5  ) then
    if (peinf%inode==0) write(6,*)
    if (peinf%inode==0) write(6,'(a)') ' Reading HDF5 wavefuntion (WFN_co.h5)'
    ngktot = SUM(kp_co%ngk)
    SAFE_ALLOCATE( gvec_kpt_components_all, (3,ngktot) )
#ifdef HDF5
    call read_hdf5_wfn_gvectors('WFN_co.h5', gvec_kpt_components_all, ngktot)
#endif

    iband_min      = MINVAL(kp_co%ifmax(:,:)) - xct%nvb_co + 1
    iband_min      = MAX(iband_min, 1)
    iband_max      = MAXVAL(kp_co%ifmax(:,:)) + xct%ncb_co
    iband_max      = MIN(iband_max, kp_co%mnband)
    tot_nbands     = iband_max - iband_min + 1
    band_block     = (tot_nbands + peinf%npes - 1) / peinf%npes
    ! read at least 128 Mb per MPI task
    min_band_block = 128.0D+00 / ( dble(ngktot) * dble(kp_co%nspin*kp_co%nspinor) * 16.0D+00 /1024.0D+00/1024D+00 )
    min_band_block = MAX(min_band_block, 1)
    if ( xct%wfn_hdf5_min_band_block > 0 )  min_band_block = xct%wfn_hdf5_min_band_block
    if ( min_band_block > band_block ) then
      band_block = min_band_block
    end if

    ib_first = iband_min + band_block * peinf%inode
    ib_last  = min(ib_first + band_block - 1, iband_max)
    ib_size  = ib_last - ib_first + 1
    if ( ib_size < 1 ) then
      ! don`t read
      ib_first = -1
      ib_last  = -1
      ib_size  =  0
    end if
    SAFE_ALLOCATE( ib_size_array, (peinf%npes) )
    ib_size_array = 0
    ib_size_array(peinf%inode+1) = ib_size

    !XXX Debug
    !XXX write(*,*) peinf%inode, iband_min, iband_max, tot_nbands, band_block, min_band_block, ib_first, ib_last, ib_size

#ifdef HDF5
#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, mpierr)
    call h5fopen_f('WFN_co.h5', H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id,error)
    call MPI_AllReduce(MPI_IN_PLACE, ib_size_array, peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
    call h5fopen_f('WFN_co.h5', H5F_ACC_RDONLY_F, file_id, error)
#endif
#endif

    ib_size_max = MAXVAL( ib_size_array )
    SAFE_ALLOCATE( my_cg_all, (ngktot, kp_co%nspin*kp_co%nspinor, MAX(ib_size,1)) )
    SAFE_ALLOCATE( ipe_cg_all, (ngktot, kp_co%nspin*kp_co%nspinor, ib_size_max) )
    my_cg_all = ZERO

#ifdef HDF5
    ! dummy allocation
    SAFE_ALLOCATE( peinf%does_it_ownv, (1, 1) )
    peinf%does_it_ownv = .true.
    call read_hdf5_bands_block(file_id, kp_co, ib_size_max, ib_size, &
                               peinf%does_it_ownv, ib_first, my_cg_all, is_bse = .true.)

    SAFE_DEALLOCATE_P(peinf%does_it_ownv)

    call h5fclose_f(file_id, error)
#endif

    ! allocate only once
    if (xct%qflag.ne.0) then
      SAFE_ALLOCATE(wfnv%cg, (MAXVAL(kp_co%ngk),wfnv%nband,wfnv%nspin*wfnv%nspinor))
    endif
    SAFE_ALLOCATE(wfnc%cg, (MAXVAL(kp_co%ngk),wfnc%nband,wfnc%nspin*wfnc%nspinor))
    SAFE_ALLOCATE(cgarray, (MAXVAL(kp_co%ngk)) )
  end if

  if ( xct%use_wfn_hdf5 ) then
    call progress_init(prog_info, 'reading wavefunctions (WFN_co.h5)', 'MPI task', peinf%npes)
  else
    call progress_init(prog_info, 'reading wavefunctions (WFN_co)', 'k-point', kp_co%nrk)
  end if

  npes_hdf5 = 0
  if ( xct%use_wfn_hdf5 ) npes_hdf5 = peinf%npes-1

do ipe = 0, npes_hdf5
  if ( xct%use_wfn_hdf5 ) then
    call progress_step(prog_info, ipe+1)
    if ( ib_size_array(ipe+1) > 0) then
     ipe_cg_all = ZERO
     if ( ipe == peinf%inode ) then
       ipe_cg_all(:, :, 1:ib_size) = my_cg_all(:, :, 1:ib_size)
     end if
#ifdef MPI
     call MPI_Bcast(ipe_cg_all, ngktot * kp_co%nspin*kp_co%nspinor * ib_size_max, &
                    MPI_SCALAR, ipe, MPI_COMM_WORLD, mpierr)
#endif
    else
      ! no elements on this MPI task, cycle
      cycle
    end if
  end if

  istart = 1
  do irk=1,kp_co%nrk
    if (.not. xct%use_wfn_hdf5) call progress_step(prog_info, irk)
    irks = 0
    do ii=1,kg%nr
      if (indxk(ii) == irk) then
        irks=ii
        exit
      endif
    enddo

    SAFE_ALLOCATE(gvec_kpt%components, (3, kp_co%ngk(irk)))
    if ( xct%use_wfn_hdf5 ) then
      gvec_kpt%components(:,:) = gvec_kpt_components_all(1:3, istart:istart+kp_co%ngk(irk)-1)
    else
      call read_binary_gvectors(25, kp_co%ngk(irk), kp_co%ngk(irk), gvec_kpt%components)
    end if

    SAFE_ALLOCATE(cg, (kp_co%ngk(irk),kp_co%nspin*kp_co%nspinor))
    if(irks > 0) then
      do ii = 1, kp_co%ngk(irk)
        call findvector(wfnv%isort(ii), gvec_kpt%components(:, ii), gvec)
        if (wfnv%isort(ii) == 0) call die('Could not find g-vector.')
      enddo
      
      wfnv%ng=kp_co%ngk(irk)
      wfnc%ng=kp_co%ngk(irk)
      if ( xct%use_wfn_hdf5 ) then
        ! for HDF5 we allocate on all MPI tasks only once 
        cgarray = ZERO
      else
        if(peinf%inode == 0) then
          if (xct%qflag.ne.0) then
            SAFE_ALLOCATE(wfnv%cg, (wfnv%ng,wfnv%nband,wfnv%nspin*wfnv%nspinor))
          endif
          SAFE_ALLOCATE(wfnc%cg, (wfnc%ng,wfnc%nband,wfnc%nspin*wfnc%nspinor))
          SAFE_ALLOCATE(cgarray, (kp_co%ngk(irk)))
        endif
      end if
    endif
    SAFE_DEALLOCATE_P(gvec_kpt%components)

! Loop over the bands
    !XXX  do ii=1,kp_co%mnband
    ii_start = 1
    ii_end   = kp_co%mnband
    if ( xct%use_wfn_hdf5 ) then
      ii_start = iband_min + band_block * ipe
      ii_end   = min(ii_start+ band_block - 1, iband_max)
    end if
    ii_loc = 0
    do ii = ii_start, ii_end
      ii_loc = ii_loc + 1

      ! FHJ: Don`t bother reading the WFNs if this is not a band we want.
      dont_read = (ii<=minval(kp_co%ifmax(irk,:))-xct%nvb_co) .or. &
                  (ii>maxval(kp_co%ifmax(irk,:)+xct%ncb_co))

      if ( xct%use_wfn_hdf5 ) then
        cg = ZERO
        if ( .not. dont_read ) then
          cg (:,:) = ipe_cg_all(istart:istart+kp_co%ngk(irk)-1, 1:kp_co%nspin*kp_co%nspinor, ii_loc)
        end if
      else
        call read_binary_data(25, kp_co%ngk(irk), kp_co%ngk(irk), &
          kp_co%nspin*kp_co%nspinor, cg, dont_read=dont_read, bcast=.false.)
      end if

      if(irks == 0) cycle
        
      if(peinf%inode == 0 .or. xct%use_wfn_hdf5) then
        do is=1, kp_co%nspin
          if (ii .gt. (kp_co%ifmax(irk,is)-xct%nvb_co) .and. ii .le. (kp_co%ifmax(irk,is)+xct%ncb_co)) then

            do isp=1, kp_co%nspinor
              do kk=1, kp_co%ngk(irk)
                cgarray(kk)=cg(kk, is*isp)
              end do
              if (peinf%verb_debug) then
                write(6,'(a, 3i7, 2(f18.13))') 'input_co', irks, ii, is*isp, cgarray(1)
              endif
              if (xct%qflag.ne.0) then
                if ((ii.le.kp_co%ifmax(irk,is)).and. &
                  (ii.gt.kp_co%ifmax(irk,is)-xct%nvb_co)) &
                  wfnv%cg(1:wfnv%ng,kp_co%ifmax(irk,is)-ii+1,is*isp)=cgarray(1:wfnv%ng)
              endif
              if ((ii.gt.kp_co%ifmax(irk,is)).and. &
                (ii.le.kp_co%ifmax(irk,is)+xct%ncb_co)) &
                wfnc%cg(1:wfnc%ng,ii-kp_co%ifmax(irk,is),is*isp)=cgarray(1:wfnc%ng)
            enddo
            call checknorm('WFN_co',ii,irks,kp_co%ngk(irk),is,kp_co%nspinor,cg(:,:))
          end if
        enddo ! loop over spins
      endif
      if (ii>maxval(kp_co%ifmax)+xct%ncb_co .and. irk==kp_co%nrk) then
        exit
      endif

    enddo ! ii (loop over bands)
    
    SAFE_DEALLOCATE(cg)
    if(peinf%inode == 0 .and. (.not.(xct%use_wfn_hdf5))) then
      SAFE_DEALLOCATE(cgarray)
    endif

    if ( .not. xct%use_wfn_hdf5 ) then
#ifdef MPI
      if (peinf%inode.ne.0) then
        if (xct%qflag.ne.0) then
          SAFE_ALLOCATE(wfnv%cg, (wfnv%ng,wfnv%nband,wfnv%nspin*wfnv%nspinor))
        endif
        SAFE_ALLOCATE(wfnc%cg, (wfnc%ng,wfnc%nband,wfnc%nspin*wfnc%nspinor))
      endif
      if (xct%qflag.ne.0) then
        call MPI_BCAST(wfnv%cg(1,1,1),wfnv%ng*wfnv%nband*wfnv%nspin*wfnv%nspinor, &
          MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
      endif
      call MPI_BCAST(wfnc%cg(1,1,1),wfnc%ng*wfnc%nband*wfnc%nspin*wfnc%nspinor, &
        MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#endif
    end if ! .not. HDF5

    distgwfco%ng(irks)=wfnv%ng
    do ii=1,distgwfco%ngl
      if (ii+distgwfco%tgl.le.wfnv%ng) &
        distgwfco%isort(ii,irks)=wfnv%isort(ii+distgwfco%tgl)
    enddo
    if (xct%qflag.ne.0) then
      if ( xct%use_wfn_hdf5 ) then
        !
        do is=1, kp_co%nspin
          do isp=1, kp_co%nspinor
            kk = is * isp
            do jj = ii_start, MIN(ii_end, kp_co%ifmax(irk,is))
              if ((jj.le.kp_co%ifmax(irk,is)).and. &
                  (jj.gt.kp_co%ifmax(irk,is)-xct%nvb_co)) then
                do ii=1,distgwfco%ngl
                  if (ii+distgwfco%tgl.le.wfnv%ng) then
                    distgwfco%zv(ii,kp_co%ifmax(irk,is)-jj+1,kk,irks)=wfnv%cg(ii+distgwfco%tgl,kp_co%ifmax(irk,is)-jj+1,kk)
                  endif
                enddo
              end if
            end do
          end do
        end do
        !
      else
        do kk=1,distgwfco%ns*distgwfco%nspinor
          do jj=1,distgwfco%nv
            do ii=1,distgwfco%ngl
              if (ii+distgwfco%tgl.le.wfnv%ng) then
                distgwfco%zv(ii,jj,kk,irks)=wfnv%cg(ii+distgwfco%tgl,jj,kk)
              endif
            enddo
          enddo
        enddo
      end if
    endif
    if ( xct%use_wfn_hdf5 ) then
      !
      do is=1, kp_co%nspin
        do isp=1, kp_co%nspinor
          kk = is * isp
          do jj = MAX(ii_start,kp_co%ifmax(irk,is)+1), ii_end
            if ((jj.gt.kp_co%ifmax(irk,is)).and. &
                (jj.le.kp_co%ifmax(irk,is)+xct%ncb_co)) then
              do ii=1,distgwfco%ngl
                if (ii+distgwfco%tgl.le.wfnv%ng) then
                  distgwfco%zc(ii,jj-kp_co%ifmax(irk,is),kk,irks)=wfnc%cg(ii+distgwfco%tgl,jj-kp_co%ifmax(irk,is),kk)
                endif
              enddo
            end if
          end do
        end do
      end do
      !
    else
      do kk=1,distgwfco%ns*distgwfco%nspinor
        do jj=1,distgwfco%nc
          do ii=1,distgwfco%ngl
            if (ii+distgwfco%tgl.le.wfnv%ng) then
              distgwfco%zc(ii,jj,kk,irks)=wfnc%cg(ii+distgwfco%tgl,jj,kk)
            endif
          enddo
        enddo
      enddo
    end if

    if ( .not.(xct%use_wfn_hdf5) ) then
      if (xct%qflag.ne.0)  then
        SAFE_DEALLOCATE_P(wfnv%cg)
      endif
      SAFE_DEALLOCATE_P(wfnc%cg)
    end if

    istart = istart + kp_co%ngk(irk)
  enddo ! irk (loop over k-points)
end do ! ipe loop for HDF5 case
  call progress_free(prog_info)

  if ( xct%use_wfn_hdf5  ) then
    SAFE_DEALLOCATE( gvec_kpt_components_all )
    SAFE_DEALLOCATE( ib_size_array )
    SAFE_DEALLOCATE( my_cg_all )
    SAFE_DEALLOCATE( ipe_cg_all )
    SAFE_DEALLOCATE( cgarray )
    if (xct%qflag.ne.0)  then
      SAFE_DEALLOCATE_P(wfnv%cg)
    endif
    SAFE_DEALLOCATE_P(wfnc%cg)
  end if

  SAFE_DEALLOCATE(indxk)
  SAFE_DEALLOCATE_P(wfnv%isort)  
  
  if (peinf%inode.eq.0) then
    write(6,'(/,1x,a)') 'Coarse-grid wavefunctions read from file WFN_co:'
    write(6,'(1x,a,i0)') '- Number of k-points in irreducible BZ: ', kg%nr
    write(6,'(1x,a,i0)') '- Number of k-points in full BZ: ', kg%nf
    if (peinf%verb_high) then
      write(6,'(1x,a)') '- Listing all k-points:'
      write(6,'(1(2x,3(1x,f10.6)))') (kg%r(:,jj), jj=1,kg%nr)
    endif
    if ( .not. xct%use_wfn_hdf5 ) call close_file(25)
  endif ! node 0

  if (xct%qflag.ne.0) then
    SAFE_DEALLOCATE(kp_co%rk)
    SAFE_DEALLOCATE(kp_co%ifmin)
    SAFE_DEALLOCATE(kp_co%ifmax)
    SAFE_DEALLOCATE(kp_co%el)
  endif
  
  POP_SUB(input_co)

  return
end subroutine input_co

end module input_co_m
