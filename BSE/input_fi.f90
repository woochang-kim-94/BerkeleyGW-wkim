#include "f_defs.h"

module input_fi_m

contains

!> Read and initialize data in preparation for BSE calculation. That is:
!! - k-points grid
!! - quasiparticle corrections
!! - k-points distribution among PEs
!! - G-vectors
!! - wavefunctions on the fine grid
!!
!!    output: crys, gvec, kg, syms, eqp, xct, int_wavefunction types
!!            INT_CWFN_* files (if flag%vm.ne.1.or. .not. flag%read_dtmat)
!!
subroutine input_fi(crys,gvec,kg,kp,syms,eqp,xct,flag, &
  omega_plasma,is_diag,intwfnc,read_all_bands)

  use checkbz_m
  use distrib_m
  use eqpcor_m
  use fullbz_m
  use global_m
#ifdef HDF5
  use hdf5
#endif
  use input_utils_m
  use intpts_m
  use io_utils_m
  use misc_m
  use scissors_m
  use sort_m
  use wfn_rho_vxc_io_m
  use wfn_io_hdf5_m
  implicit none

  type (crystal), intent(out) :: crys
  type (gspace), intent(out) :: gvec
  type (grid), intent(out) :: kg
  type (kpoints), intent(out) :: kp
  type (symmetry), intent(out) :: syms
  type (eqpinfo), intent(inout) :: eqp
  type (xctinfo), intent(inout) :: xct
  type (flags), intent(inout) :: flag
  real(DP), intent(out) :: omega_plasma
  logical, intent(in) :: is_diag
  type (int_wavefunction), intent(out) :: intwfnc
  logical, intent(in),optional :: read_all_bands  !< Read all bands instead of
                                                  !! only conduction bands

  type (wavefunction) :: wfnc
  character :: filenamec*20
  character :: tmpfn*16
  character :: fncor*32
  integer :: iunit_c,iwrite

  integer :: ii,jj,kk,nn,nmat,is,isp,ik,irk,irk_kp,ib,ib_kp
  integer :: dest,tag,ic,iv,iwritetotal,ijk
  integer :: irks
  integer :: ncb, nvb
  real(DP) :: diffvol,vcell,kt(3),div
  real(DP) :: test,tol,qtot
  real(DP), allocatable :: eltr(:),kltr(:,:),ek_tmp(:)
  integer, allocatable :: kcvltr(:,:),indxk(:),k_tmp(:,:)
  integer, allocatable :: index(:),isend(:),iwriteik(:)
  SCALAR, allocatable :: cg(:,:), cgarray(:)

  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvec_kpt
  integer :: last_ng, last_ng_match, last_ikt
  logical :: skip_checkbz, broken_degeneracy, dont_read
  !> Number of bands that were corrected (via scissors or eqp). We can`t find 
  !! the FE using more bands than this number.
  integer :: bands_cor, minband
  type(progress_info) :: prog_info
  logical :: read_all_bands_

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

  if (present(read_all_bands)) then
    read_all_bands_ = read_all_bands
  else
    read_all_bands_ = .false.
  end if

!-----------------------------------------------------------------------
! Read info for crystal from WFN_fi

  PUSH_SUB(input)

  sheader = 'WFN'
  iflavor = 0
  if ( xct%use_wfn_hdf5 ) then
    if (peinf%inode==0) write(6,'(a)') ' Reading header of WFN_fi.h5'
#ifdef HDF5
    call read_hdf5_header_type('WFN_fi.h5', sheader, iflavor, kp, gvec, syms, crys)
#endif
  else
    if(peinf%inode == 0) call open_file(25,file='WFN_fi',form='unformatted',status='old')
    call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys, &
      dont_warn_kgrid=.true.)
  end if
  call check_trunc_kpts(xct%icutv, kp)

!-----------------------------------------------------------------------
!     Read info for g-vectors from WFN_fi
!

  call logit('input:  reading gvec info')

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  if( xct%use_wfn_hdf5 ) then
#ifdef HDF5
    call read_hdf5_gvectors('WFN_fi.h5', gvec%ng, gvec%components)
#endif
  else
    call read_binary_gvectors(25, gvec%ng, gvec%ng, gvec%components)
  end if

  call get_volume(vcell,crys%bdot)
  diffvol=abs(crys%celvol-vcell)
  if (diffvol.gt.TOL_Small) then
    call die('volume mismatch', only_root_writes = .true.)
  endif

  xct%nspin = kp%nspin
  xct%nspinor = kp%nspinor
  if (xct%nspinor==2) then
    ! For a local fields calculations with spinors, overwrite the krnl flag and set direct_fact to
    ! zero out direct term
    if (xct%krnl==2) then
      xct%direct_fact = 0d0
    endif
    flag%krnl = 3
    xct%krnl = 3
    if (xct%theory.eq.1) call die('TDDFT not compatible with spinors')
  endif

  if(any(kp%ifmax(:,:) == 0)) &
    call die("BSE codes cannot handle a system where some k-points have no occupied bands.", only_root_writes = .true.)

  kp%nvband=minval(kp%ifmax(:,:)-kp%ifmin(:,:))+1
  kp%ncband=kp%mnband-maxval(kp%ifmax(:,:))

!-----------------------------------------------------------
! Manual override of band numbers

  if (xct%vmax.ne.0) then
    kp%nvband = xct%vmax-xct%vmin+1
    kp%ncband = kp%mnband-xct%vmax
    if (peinf%inode.eq.0) then
      write(6,*)
      write(6,*) '*** Overwrite min/max occupied state for fine-grid wavefunctions'
      write(6,*) '*** kp%nvband =',kp%nvband,' kp%ncband =',kp%ncband
      write(6,*)
    endif
  endif

  ! Define the actual number of conduction and valence bands used.
  ! These include the whole set of bands if it was requested.
  if (read_all_bands_) then
    nvb = kp%mnband
    ncb = kp%mnband
  else
    nvb = xct%nvb_fi
    ncb = xct%ncb_fi
  end if

!-----------------------------------------------------------
! Manual override of regular grid
! it matters only if you want to perform minicell averages

  if (xct%rgrid(1).ne.0) kp%kgrid = xct%rgrid

!----------------------------------------------------------------
! (gsm) check whether the requested number of bands
!       is available in the wavefunction file

! we only use the conduction bands from WFN_fi
! if flag%vm.eq.1.and.flag%read_dtmat we don`t need them at all

  if (flag%vm .ne. 1 .or. .not. flag%read_dtmat) then

    if (xct%ncb_fi.gt.kp%ncband .and. .not. read_all_bands_) then
      call die("The requested number of conduction bands is not available in WFN_fi.")
    endif

  endif

!----------------------------------------------------------------------
!     Read the k-point sampling from kpoints (if it exists) or from
!     WFN_fi. In either case, this sampling will define the irreducible
!     Brillouin zone.
!     Also, set up the correspondance between the irreducible k-points
!     in the grid object and the kpoints object. 

  if (xct%read_kpoints) then
    if (peinf%inode.eq.0) then
      call open_file(9,file='kpoints',form='formatted',status='old')
      read(9,*) kg%nr
      SAFE_ALLOCATE(kg%r, (3,kg%nr))
      do ii=1,kg%nr
        read(9,*) (kg%r(jj,ii),jj=1,3),div
        kg%r(:,ii) = kg%r(:,ii)/div
      enddo
      call close_file(9)
    endif
#ifdef MPI
    call MPI_BCAST(kg%nr,   1,     MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if(peinf%inode.ne.0) then
      SAFE_ALLOCATE(kg%r, (3,kg%nr))
    endif
    call MPI_BCAST(kg%r,    3*kg%nr,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif

!     indxk : stores the correspondence between k-points kg%r and kp%rk
!     (it is used to select the set of wavefunctions to be stored)
!     tol : tolerance in the coordinates of k-points

    tol = TOL_Small
    SAFE_ALLOCATE(indxk, (kg%nr))
    indxk=0
    do jj=1,kg%nr
      do ii=1,kp%nrk
        kt(:) = kg%r(:,jj) - kp%rk(:,ii)
        if (all(abs(kt(1:3)).lt.tol)) then
          if (indxk(jj).ne.0) write(0,*) 'WARNING: multiple definition of k-point',jj,indxk(jj),kg%r(:,jj)
          indxk(jj)=ii
        endif
      enddo
      
!     If some k-point listed in kg%r is not found in WFN_fi, indxk
!     will store zero. Later, the job will stop in genwf.

      if (indxk(jj).eq.0) write(0,*) 'WARNING: could not find vector ',kg%r(:,jj),' in WFN_fi'
    enddo
  else
    kg%nr=kp%nrk
    SAFE_ALLOCATE(kg%r, (3,kg%nr))
    kg%r(1:3,1:kg%nr)=kp%rk(1:3,1:kp%nrk)
    SAFE_ALLOCATE(indxk, (kg%nr))
    do ii=1,kg%nr
      indxk(ii) = ii
    enddo
  endif

!-----------------------------------------------------------------------
!     Generate full brillouin zone from irreducible wedge, rk -> fk
!
!     If flag%bz0.eq.1, only Identity will be used as
!     symmetry operation. In this case, kg%r (irreducible BZ) and kg%f
!     (full BZ) will be identical.
!
  if (flag%bz0.eq.1 .and. .not. xct%is_absorption) then
    ! in the inteqp code, we want to leave the k-points alone, if symmetries are not being used.
    call fullbz(crys,syms,kg,1,skip_checkbz,wigner_seitz=.false.,paranoid=.false.,do_nothing=.true.)
  else if (flag%bz0.eq.1) then
    call fullbz(crys,syms,kg,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  else
    call fullbz(crys,syms,kg,syms%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  endif
  tmpfn='WFN_fi'
  if (.not. skip_checkbz .and. .not.xct%patched_sampling) then
    call checkbz(kg%nf,kg%f,kp%kgrid,kp%shift,crys%bdot, &
      tmpfn,'k',.true.,xct%freplacebz,xct%fwritebz)
  endif
!
  if (flag%bz0.eq.0.and.peinf%inode.eq.0) write(6,801)
  if (flag%bz0.eq.1.and.peinf%inode.eq.0) write(6,802)
801 format(1x,'Using symmetries to expand the fine-grid sampling')
802 format(1x,'No symmetries used in the fine grid-sampling')
!
  xct%nkpt_fi = kg%nf
  xct%nktotal = kg%nf
  if (xct%patched_sampling) then
    ! FHJ: the main assumption of the patched sampling is that the volume of
    ! each mini-BZ is given by the number of k-points before decimating the WFN.
    ! Still, shouldn`t this line be used always, even for regular sampling?
    xct%nktotal = product(kp%kgrid)
  endif
  
  if (peinf%inode .eq. 0) then
    write(6,"(1x,a,3(1x,i0))") 'K-grid:', kp%kgrid
  endif
  
!-----------------------------------------------------------------------
! If quasi-particle correction requested, read the corrected
! qp energies from file (in eV)

  SAFE_ALLOCATE(kp%elda, (kp%mnband, kp%nrk, kp%nspin))
  kp%el(:,:,:) = kp%el(:,:,:) - xct%avgpot / ryd
  kp%elda(:,:,:) = kp%el(:,:,:)
  call scissors_shift(kp, eqp%scis, eqp%spl_tck)

  bands_cor = kp%mnband
  minband = 1
  if(xct%eqp_corrections) then
    fncor='eqp.dat'
    bands_cor = maxval(kp%ifmax(:,:)) + xct%ncb_fi
    minband = minval(kp%ifmax(:,:)-xct%nvb_fi+1)
    ! FIXME: for metals this is asking for a few more bands than actually needed on some k-points
    call eqpcor(fncor,peinf%inode,peinf%npes,kp,&
      minval(kp%ifmax(:,:)-xct%nvb_fi+1),bands_cor,0,0,kp%el,kp%el,kp%el,1,0)
  endif

  if(xct%unrestricted_transf) then
    call find_efermi(xct%rfermi, xct%efermi, xct%efermi_input, kp, bands_cor, minband, &
      "fine grid", should_search = .true., should_update = .true., write7 = .false.,  dont_die_consistency = .true.)
  else
    call find_efermi(xct%rfermi, xct%efermi, xct%efermi_input, kp, bands_cor, minband, &
      "fine grid", should_search = .true., should_update = .true., write7 = .false.,  dont_die_consistency = .false.)
  endif

  call assess_degeneracies(kp, kp%el(kp%mnband, :, :), kp%mnband - 1, xct%efermi, TOL_Degeneracy)

  call calc_qtot(kp, crys%celvol, xct%efermi, qtot, omega_plasma, write7 = .false.)

! DAS: degenerate subspace check

  ! degeneracy does not matter for WFN_fi in inteqp
  if((flag%vm .ne. 1 .or. .not. flag%read_dtmat) .and. peinf%inode == 0 .and. xct%is_absorption) then
    if(xct%ncb_fi.eq.kp%ncband) then
      call die("You must provide one more conduction band in WFN_fi in order to assess degeneracy.")
    endif
    broken_degeneracy = .false.
    do jj = 1, kp%nspin
      do ii = 1, kp%nrk
        if(abs(kp%elda(kp%ifmax(ii, jj) + xct%ncb_fi, ii, jj) &
          - kp%elda(kp%ifmax(ii, jj) + xct%ncb_fi + 1, ii, jj)) .lt. TOL_Degeneracy) then
          broken_degeneracy = .true.
        endif
      enddo
    enddo

    if(broken_degeneracy) then
      if(xct%degeneracy_check_override) then
        write(0,'(a)') &
          "WARNING: Selected number of conduction bands breaks degenerate subspace in WFN_fi. " // &
          "Run degeneracy_check.x for allowable numbers."
        write(0,*)
      else
        write(0,'(a)') &
          "Run degeneracy_check.x for allowable numbers, or use keyword " // &
          "degeneracy_check_override to run anyway (at your peril!)."
        call die("Selected number of conduction bands breaks degenerate subspace in WFN_fi.")
      endif
    endif
  endif

!-----------------------------------------------------------------------
!     Work with energy arrays: ec,ev,eclda,evlda (all in Ryd!)
!     xct%evqp          QP valence energies
!     xct%ecqp          QP conduction energies
!     xct%evlda         LDA valence energies (used only
!        with momentum operator)
!     xct%eclda         LDA conduction energies (idem)
!      The label of bands in xct%evqp and xct%evlda is again reversed:
!      from higher to lower energies.

  if (read_all_bands_) then
    eqp%band_ordering = 1
  else
    eqp%band_ordering = 0
  end if

  SAFE_ALLOCATE(eqp%evqp, (nvb,xct%nkpt_fi,xct%nspin))
  SAFE_ALLOCATE(eqp%ecqp, (ncb,xct%nkpt_fi,xct%nspin))
  SAFE_ALLOCATE(eqp%evlda, (nvb,xct%nkpt_fi,xct%nspin))
  SAFE_ALLOCATE(eqp%eclda, (ncb,xct%nkpt_fi,xct%nspin))

! FHJ: loop thru indices on `eqp` grid, and then find the
!      corresponding labels on the `kp` grid
  do is=1,xct%nspin
    do irk=1,xct%nkpt_fi
      !irkq = index of q-pt on `kp` grid
      irk_kp = indxk(kg%indr(irk))

      ! If we read all bands, then both vb and cb are counted from bottom up 
      if (read_all_bands_) then

        do ib_kp=1, kp%mnband
          eqp%evqp (ib_kp, irk, is) = kp%el  (ib_kp, irk_kp, is)
          eqp%evlda(ib_kp, irk, is) = kp%elda(ib_kp, irk_kp, is)
          eqp%ecqp (ib_kp, irk, is) = kp%el  (ib_kp, irk_kp, is)
          eqp%eclda(ib_kp, irk, is) = kp%elda(ib_kp, irk_kp, is)
        enddo

      else

        !loop thru bands of `kp` grid
        do ib_kp=1, kp%mnband

          !ib = band index on `eqp` grid
          ib = kp%ifmax(irk_kp,is) - ib_kp + 1

          if (ib > 0) then
            !if ib is one of the selected valence bands
            if (ib <= xct%nvb_fi) then
              eqp%evqp (ib, irk, is) = kp%el  (ib_kp, irk_kp, is)
              eqp%evlda(ib, irk, is) = kp%elda(ib_kp, irk_kp, is)
            endif
          else
            ib = 1 - ib !transform negative val. index into positive cond. index
            !if ib is one of the selected conduction bands
            if (ib <= xct%ncb_fi) then
              eqp%ecqp (ib, irk, is) = kp%el  (ib_kp, irk_kp, is)
              eqp%eclda(ib, irk, is) = kp%elda(ib_kp, irk_kp, is)
            endif
          endif

        enddo
      endif
    enddo
  enddo
      
  if (peinf%inode.eq.0) then
    call scissors_write(6, eqp%scis)
  endif

!-----------------------------------------------------------------------
!       Distribute kpoints among the PEs
!
  call logit('input:  calling distrib')
  call distrib(xct,gvec%FFTgrid,kp%kgrid,is_diag)

!-----------------------------------------------------------------------
! gsm: store xct%ifmax for generating eqp/eqp_q.dat in intwfn

  SAFE_ALLOCATE(xct%ifmax, (xct%nkpt_fi,xct%nspin))
  xct%ifmax(:,:) = kp%ifmax(kg%indr(:),:)

!-----------------------------------------------------------------------
!     Order g-vectors with respect to their kinetic energy

  call logit('input:  reordering gvecs')
  SAFE_ALLOCATE(index, (gvec%ng))
  SAFE_ALLOCATE(gvec%ekin, (gvec%ng))
  call kinetic_energies(gvec, crys%bdot, gvec%ekin)
  call sortrx(gvec%ng, gvec%ekin, index, gvec = gvec%components)

  SAFE_ALLOCATE(ek_tmp, (gvec%ng))
  ek_tmp = gvec%ekin
  SAFE_ALLOCATE(k_tmp, (3,gvec%ng))
  k_tmp = gvec%components
  do ii=1,gvec%ng
    gvec%ekin(ii) = ek_tmp(index(ii))
    gvec%components(:,ii) = k_tmp(:,index(ii))
  enddo
  SAFE_DEALLOCATE(ek_tmp)
  SAFE_DEALLOCATE(k_tmp)
  SAFE_DEALLOCATE(index)

  call gvec_index(gvec)
  
!-----------------------------------------------------------------------
!     Read the wavefunctions and create INT_CWFN_*
!
  if (flag%read_dtmat .and. flag%vm /= 0) then
    if (peinf%inode.eq.0) write(6,*) ' Bypassing INT_CWFN_*'
  else
    
    call logit('input:  reading WFN_fi')
    wfnc%nband=ncb
    wfnc%nspin=kp%nspin
    wfnc%nspinor=kp%nspinor

  if ( xct%use_wfn_hdf5  ) then
      if (peinf%inode==0) write(6,*)
      if (peinf%inode==0) write(6,'(a)') ' Reading HDF5 wavefuntion (WFN_fi.h5)'
      ngktot = SUM(kp%ngk)
      SAFE_ALLOCATE( gvec_kpt_components_all, (3,ngktot) )
#ifdef HDF5
      call read_hdf5_wfn_gvectors('WFN_fi.h5', gvec_kpt_components_all, ngktot)
#endif

      if ( read_all_bands_ ) then
        iband_min = 1
        iband_max = kp%mnband
      else
        iband_min      = MINVAL(kp%ifmax(:,:)) + 1 !XXXX test this
        iband_min      = MAX(iband_min, 1)
        iband_max      = MAXVAL(kp%ifmax(:,:)) + xct%ncb_fi
        iband_max      = MIN(iband_max, kp%mnband)
      end if
      tot_nbands     = iband_max - iband_min + 1
      band_block     = (tot_nbands + peinf%npes - 1) / peinf%npes
      ! read at least 128 Mb per MPI task
      min_band_block = 128.0D+00 / ( dble(ngktot) * dble(kp%nspin*kp%nspinor) * 16.0D+00 /1024.0D+00/1024D+00 )
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
      !XXX write(2000+peinf%inode,*) peinf%inode, iband_min, iband_max, tot_nbands, band_block, min_band_block, ib_first, ib_last, ib_size
      !XXX flush(2000+peinf%inode)
#ifdef HDF5
#ifdef MPI
       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
       call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, mpierr)
       call h5fopen_f('WFN_fi.h5', H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
       call h5pclose_f(plist_id,error)
       call MPI_AllReduce(MPI_IN_PLACE, ib_size_array, peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
       call h5fopen_f('WFN_fi.h5', H5F_ACC_RDONLY_F, file_id, error)
#endif
#endif

       ib_size_max = MAXVAL( ib_size_array )
       SAFE_ALLOCATE( my_cg_all, (ngktot, kp%nspin*kp%nspinor, MAX(ib_size,1)) )
       SAFE_ALLOCATE( ipe_cg_all, (ngktot, kp%nspin*kp%nspinor, ib_size_max) )
       my_cg_all = ZERO

#ifdef HDF5
       ! dummy allocation
       SAFE_ALLOCATE( peinf%does_it_ownv, (1, 1) )
       peinf%does_it_ownv = .true.
       call read_hdf5_bands_block(file_id, kp, ib_size_max, ib_size, &
                                  peinf%does_it_ownv, ib_first, my_cg_all, is_bse = .true.)

       SAFE_DEALLOCATE_P(peinf%does_it_ownv)
       call h5fclose_f(file_id, error)
#endif

    end if

    SAFE_ALLOCATE(intwfnc%ng, (peinf%ikt(peinf%inode+1)))
    SAFE_ALLOCATE(intwfnc%isort, (gvec%ng,peinf%ikt(peinf%inode+1)))
    SAFE_ALLOCATE(intwfnc%cgk, (kp%ngkmax,wfnc%nband,kp%nspin*kp%nspinor,peinf%ikt(peinf%inode+1)))
    intwfnc%nspinor=wfnc%nspinor

    SAFE_ALLOCATE(wfnc%isort, (gvec%ng))
    SAFE_ALLOCATE(isend, (peinf%npes))
  if ( xct%use_wfn_hdf5 ) then
    call progress_init(prog_info, 'reading wavefunctions (WFN_fi.h5)', 'MPI task', peinf%npes)
  else
    call progress_init(prog_info, 'reading wavefunctions (WFN_fi)', 'k-point', kp%nrk)
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
       call MPI_Bcast(ipe_cg_all, ngktot * kp%nspin*kp%nspinor * ib_size_max, &
                      MPI_SCALAR, ipe, MPI_COMM_WORLD, mpierr)
#endif
      else
        ! no elements on this MPI task, cycle
        cycle
      end if
    end if

    wfnc%isort=0
    last_ikt=-1
    last_ng=-1
    last_ng_match=-1
    istart = 1
    do irk=1,kp%nrk
      if (.not. xct%use_wfn_hdf5) call progress_step(prog_info, irk)
      irks = 0
      do ii=1,kg%nr
        if (irk.eq.indxk(ii)) then
          irks=ii
          exit
        endif
      enddo

      wfnc%ng = kp%ngk(irk)

!   FHJ: Realloc arrays. Note that we can`t do something like
!        "if (wfnc%ng>last_ng)"  b/c fortran complains at read_binary_gvectors 
!        if the vectors are not exactly wfnc%ng big.
      if(wfnc%ng/=last_ng) then
        if(last_ng/=-1) then
          SAFE_DEALLOCATE_P(gvec_kpt%components)
          SAFE_DEALLOCATE(cg)
        endif
        SAFE_ALLOCATE(gvec_kpt%components, (3, wfnc%ng))
        SAFE_ALLOCATE(cg, (wfnc%ng, kp%nspin*kp%nspinor))
        last_ng = wfnc%ng
      endif

      if ( xct%use_wfn_hdf5 ) then
        gvec_kpt%components(:,:) = gvec_kpt_components_all(1:3, istart:istart+kp%ngk(irk)-1)
      else
        call read_binary_gvectors(25, wfnc%ng, wfnc%ng, gvec_kpt%components)
      end if

      if(irks > 0) then
        do ii = 1, kp%ngk(irk)
          call findvector(wfnc%isort(ii), gvec_kpt%components(:, ii), gvec)
          if(wfnc%isort(ii) == 0) call die('input: could not find gvec')
        enddo

!   FHJ: Realloc arrays.
        if(wfnc%ng/=last_ng_match) then
          if(last_ng_match/=-1) then
            SAFE_DEALLOCATE_P(wfnc%cg)
            SAFE_DEALLOCATE(cgarray)
          endif
          SAFE_ALLOCATE(wfnc%cg, (wfnc%ng, wfnc%nband, wfnc%nspin*wfnc%nspinor))
          SAFE_ALLOCATE(cgarray, (wfnc%ng))
          last_ng_match = wfnc%ng
        endif
        if(peinf%ikt(peinf%inode+1)/=last_ikt) then
          if(last_ikt/=-1) then
            SAFE_DEALLOCATE(iwriteik)
          endif
          SAFE_ALLOCATE(iwriteik,(peinf%ikt(peinf%inode+1)))
          last_ikt = peinf%ikt(peinf%inode+1)
        endif

!         Determine which PEs will write the wavefunctions for this k-point
        iwrite=0
        iwritetotal=0
        iwriteik=0
        do ii=1, peinf%ikt(peinf%inode+1)
          if(kg%indr(peinf%ik(peinf%inode+1,ii)).eq.irks) then
            iwritetotal=iwritetotal+1
            iwriteik(iwritetotal)=ii
            iwrite=1
          endif
        enddo

        if ( .not.  xct%use_wfn_hdf5 ) then
!         Determine to which PEs the wavefunctions for this k-point
!         need to be sent...
          isend=0
          if(peinf%inode.eq.0) then
            do jj=2,peinf%npes
              do ii=1, peinf%ikt(jj)
                if(kg%indr(peinf%ik(jj,ii)).eq.irks) then
                  isend(jj)=1
                  exit
                endif
              enddo
            enddo
          endif
!
        end if ! hdf5
      endif
!
!         Loop over the bands
!
      ii_start = 1
      ii_end   = kp%mnband
      if ( xct%use_wfn_hdf5 ) then
        ii_start = iband_min + band_block * ipe
        ii_end   = min(ii_start+ band_block - 1, iband_max)
      end if
      ii_loc = 0
      do ii = ii_start, ii_end
        ii_loc = ii_loc + 1

        ! Skip reading the WFNs if this is not a band we want.
        dont_read = (((ii<=minval(kp%ifmax(irks,:))) .or. &
                      (ii>maxval(kp%ifmax(irks,:))+xct%ncb_fi)) .and. &
                     .not. read_all_bands_)

        if ( xct%use_wfn_hdf5 ) then
          cg = ZERO
          if ( .not. dont_read ) then
            cg (:,:) = ipe_cg_all(istart:istart+kp%ngk(irk)-1, 1:kp%nspin*kp%nspinor, ii_loc)
          end if
        else
          call read_binary_data(25, kp%ngk(irk), kp%ngk(irk), kp%nspin*kp%nspinor, &
            cg, dont_read=dont_read, bcast=.false.)
        end if

        if(irks == 0) cycle

        do is=1, kp%nspin

          ! Skip bands out of range for this spin polarization
          if( .not. (ii.gt.kp%ifmax(irks,is).and. &
             ii.le.kp%ifmax(irks,is)+xct%ncb_fi) .and. &
             .not. read_all_bands_) then
            cycle
          end if

          if (read_all_bands_) then
            ib = ii
          else
            ib = ii - kp%ifmax(irks,is)
          end if

          do isp=1, kp%nspinor
            if ( xct%use_wfn_hdf5 ) then
              cgarray(1:kp%ngk(irk)) = cg(1:kp%ngk(irk), is*isp)
            else
              if (peinf%inode.eq.0) then
                cgarray(1:kp%ngk(irk))=cg(1:kp%ngk(irk), is*isp)
                if (peinf%verb_debug) then
                  write(6,'(a, 3(1x,i0), 2(1x,f18.13))') 'input', irks, ii, is*isp, cgarray(1)
                endif
              endif
#ifdef MPI
              if(peinf%inode.eq.0) then
                do jj=2,peinf%npes
                  if(isend(jj).eq.1) then
                    dest=jj-1
                    tag=1000+dest*2-(is-1)
                    call MPI_SEND(cgarray,kp%ngk(irk),MPI_SCALAR, &
                      dest,tag,MPI_COMM_WORLD,mpierr)
                  endif
                enddo
              else
                if(iwrite.eq.1) then
                  tag=1000+peinf%inode*2-(is-1)
                  call MPI_RECV(cgarray,kp%ngk(irk),MPI_SCALAR,0,tag, &
                    MPI_COMM_WORLD,mpistatus,mpierr)
                endif
              endif
#endif
            end if ! hdf5
            if(iwrite.eq.1) then
              wfnc%cg(1:wfnc%ng,ib,is*isp) = cgarray
            endif
          enddo
          if(iwrite.eq.1) then
            !FHJ: We only Send/Recv and check one spin at a time
            call checknorm('WFN_fi',ii,irks,kp%ngk(irk),is,kp%nspinor,&
                         wfnc%cg(:,ib,:))
          endif
        enddo
      enddo
      !XX Increment here, a lot of "if(irks == 0) cycle" statement
      istart = istart + kp%ngk(irk)
      if(irks == 0) cycle
      
      if(iwrite.eq.1) then
        
        do ijk = 1, iwritetotal
          intwfnc%ng(iwriteik(ijk))=wfnc%ng
          intwfnc%isort(:,iwriteik(ijk))=wfnc%isort(:)

          if ( xct%use_wfn_hdf5 ) then
            do ii = ii_start, ii_end
              do is=1, kp%nspin
                ! Skip bands out of range for this spin polarization
                if( .not. (ii.gt.kp%ifmax(irks,is).and. &
                   ii.le.kp%ifmax(irks,is)+xct%ncb_fi) .and. &
                   .not. read_all_bands_) then
                  cycle
                end if
                if (read_all_bands_) then
                  ib = ii
                else
                  ib = ii - kp%ifmax(irks,is)
                end if
                do isp=1, kp%nspinor
                  intwfnc%cgk(1:wfnc%ng,ib,is*isp,iwriteik(ijk))=wfnc%cg(1:wfnc%ng,ib,is*isp)
                end do
              end do
            end do
          else
            intwfnc%cgk(1:wfnc%ng,:,:,iwriteik(ijk))=wfnc%cg(1:wfnc%ng,:,:)
          end if

        enddo
      endif

    enddo !end loop over k-points

    if(last_ikt/=-1) then
      SAFE_DEALLOCATE(iwriteik)
    endif
    if(last_ng_match/=-1) then
      SAFE_DEALLOCATE_P(wfnc%cg)
      SAFE_DEALLOCATE(cgarray)
    endif
    if(last_ng/=-1) then
      SAFE_DEALLOCATE_P(gvec_kpt%components)
      SAFE_DEALLOCATE(cg)
    endif

  end do ! ipe loop for HDF5 case
    call progress_free(prog_info)

    SAFE_DEALLOCATE(isend)
    SAFE_DEALLOCATE_P(wfnc%isort)

    if ( xct%use_wfn_hdf5  ) then
      SAFE_DEALLOCATE( gvec_kpt_components_all )
      SAFE_DEALLOCATE( ib_size_array )
      SAFE_DEALLOCATE( my_cg_all )
      SAFE_DEALLOCATE( ipe_cg_all )
    end if

    !XX if(last_ikt/=-1) then
    !XX   SAFE_DEALLOCATE(iwriteik)
    !XX endif
    !XX if(last_ng_match/=-1) then
    !XX   SAFE_DEALLOCATE_P(wfnc%cg)
    !XX   SAFE_DEALLOCATE(cgarray)
    !XX endif
    !XX if(last_ng/=-1) then
    !XX   SAFE_DEALLOCATE_P(gvec_kpt%components)
    !XX   SAFE_DEALLOCATE(cg)
    !XX endif

  endif
  
  if ( .not. xct%use_wfn_hdf5 ) then
    if(peinf%inode.eq.0) call close_file(25)
  end if

  SAFE_DEALLOCATE(indxk)

!-----------------------------------------------------------------------
! Print out some information

  call write_transitions()

  POP_SUB(input)  
  return

contains

  subroutine write_transitions()

    PUSH_SUB(input_fi.write_transitions)
    
    if(peinf%inode.eq.0) then
      write(6,'(/,1x,a)') 'Conduction wavefunctions read from file WFN_fi:'
      write(6,'(1x,a,i0)') '- Number of k-points in irreducible BZ: ', kg%nr
      write(6,'(1x,a,i0)') '- Number of k-points in full BZ: ', kg%nf
      if (peinf%verb_high) then
        write(6,'(1x,a)') '- Listing all k-points:'
        write(6,'(1(2x,3(1x,f10.6)))') (kg%r(:,jj), jj=1,kg%nr)
      endif
    endif ! node 0

    !! hang on to k-points to check in input_co that grids are the same
    !if(.not. xct%skipinterp) then
    !  SAFE_DEALLOCATE(kp%rk)
    !endif
    !SAFE_DEALLOCATE(kp%ifmin)
    !SAFE_DEALLOCATE(kp%ifmax)
    !SAFE_DEALLOCATE(kp%el)

!-----------------------------------------------------------------------
!     Print out information about the nn lowest energy transitions for
!     spin up states. The array eltr has the lowest nmat energy
!     transitions, with corresponding bands/k-points in kltr, kcvltr.
!     This can also be used to locate the location in BZ
!     of all transitions within a predefined energy range.
!
    if(peinf%inode.eq.0) then
      nmat=100
      SAFE_ALLOCATE(eltr, (nmat))
      SAFE_ALLOCATE(kltr, (3,nmat))
      SAFE_ALLOCATE(kcvltr, (3,nmat))
      kltr=0.d0
      kcvltr=0
      eltr(:)=1.d8
      eltr(1)=eqp%ecqp(1,1,1) - eqp%evqp(1,1,1)
      kltr(:,1)=kg%f(:,1)
      kcvltr(:,1)=1
      
      do ik=1,xct%nkpt_fi
        do iv=1,xct%nvb_fi
          ic_loop: do ic=1,xct%ncb_fi
            if (ik.eq.1 .and. ic.eq.1 .and. iv.eq.1) cycle ic_loop
            test=eqp%ecqp(ic,ik,1) - eqp%evqp(iv,ik,1)
            do jj=1,nmat
              if (test.lt.eltr(jj)) then
                do kk=nmat-1,jj,-1
                  eltr(kk+1) = eltr(kk)
                  kltr(:,kk+1)=kltr(:,kk)
                  kcvltr(:,kk+1)=kcvltr(:,kk)
                enddo
                eltr(jj)=test
                kltr(:,jj)=kg%f(:,ik)
                kcvltr(1,jj)=ik
                kcvltr(2,jj)=ic
                kcvltr(3,jj)=iv
                cycle ic_loop
              endif
            enddo
          enddo ic_loop
        enddo
      enddo
        
      eltr(:)=eltr(:)*ryd
      nn = 5
      write(6,'(/1x,a/)') 'Lowest-energy independent-particle transitions:'
      write(6,'(1x,15x,a7,12x,3(1x,a6),1x,a15)') 'k-point', 'ik', 'ic', 'iv', 'energy (eV)'
      write(6,'(3x,32("-"),3(1x,6("-")),1x,15("-"))')
      do jj=1,nn
        if (kcvltr(1,jj)>0) then
          write(6,'(2x,3(1x,f10.6),3(1x,i6),1x,f15.6)') &
            kltr(1:3,jj), kcvltr(1:3,jj), eltr(jj)
        endif
      enddo
      SAFE_DEALLOCATE(eltr)
      SAFE_DEALLOCATE(kltr)
      SAFE_DEALLOCATE(kcvltr)
    endif

    POP_SUB(input_fi.write_transitions)
    return
      
  end subroutine write_transitions

end subroutine input_fi

end module input_fi_m
