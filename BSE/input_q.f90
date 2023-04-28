#include "f_defs.h"

module input_q_m

contains

!>    Read data from file WFN_fi or WFNq_fi and initialize variables
!!    (only called if flag%vm.ne.1.or.flag%dtm.ne.1)
!!
!!    input: crys, gvec, kg, kp, syms, xct types
!!
!!    output: kgq type
!!            indexq(1:xct%nkpt_fi) - mapping of points in shifted grid
!!            INT_VWFNQ_* files
subroutine input_q(kp,crys,gvec,kg,kgq,kpq,syms,xct, &
  indexq,eqp,flag,intwfnv,read_all_bands)
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

  type (crystal), intent(in) :: crys
  type (kpoints), intent(in) :: kp
  type (gspace), intent(in) :: gvec
  type (grid), intent(in) :: kg
  type (grid), intent(out) :: kgq
  type (kpoints), intent(out) :: kpq
  type (symmetry), intent(in) :: syms
  type (xctinfo), intent(inout) :: xct
  type (flags), intent(in) :: flag
  integer, intent(out) :: indexq(xct%nkpt_fi)
  type (eqpinfo), intent(inout) :: eqp
  type (int_wavefunction), intent(out) :: intwfnv      
  logical, intent(in),optional :: read_all_bands  !< Read all bands instead of
                                                  !! only valence bands

  type (symmetry) :: symsq
  type (crystal) :: crysq
  type (wavefunction) :: wfnv
  character :: wfnq0*16
  character :: filenamev*20
  character :: tmpfn*16
  character :: fncor*32
  character(len=64) :: tmpstr
  integer :: iunit_v,iwrite
  integer :: ncb, nvb
  integer :: ii,jj,kk,ik,ikq,irkq,is,ispinor,ib,ibq, irk
  integer :: dest,tag,dummygvec(1,1)
  integer :: iwritetotal,ijk,minband
  real(DP) :: delta,qq(3),tol
  integer, allocatable :: isend(:),iwriteik(:)
  SCALAR, allocatable :: cg(:,:), cgarray(:)

  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvecq, gvec_kpt
  logical :: irkq_match, dont_read
  integer :: last_ng, last_ng_match, last_ikt
  logical :: skip_checkbz, broken_degeneracy
  !> Number of bands that were corrected (via scissors or eqp). We can`t find 
  !! the FE using more bands than this number.
  integer :: bands_cor
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

  PUSH_SUB(input_q)

  if (present(read_all_bands)) then
    read_all_bands_ = read_all_bands
  else
    read_all_bands_ = .false.
  end if

  if (flag%opr .eq. 0) then
    wfnq0 = 'WFNq_fi'
    fncor='eqp_q.dat'
  else
    wfnq0 = 'WFN_fi'
    fncor='eqp.dat'
  endif

  sheader = 'WFN'
  iflavor = 0
  if ( xct%use_wfn_hdf5 ) then
    wfnq0 = TRIM(wfnq0)//'.h5'
    if (peinf%inode==0) write(6,'(a)') ' Reading header of '//TRIM(wfnq0)
#ifdef HDF5
    call read_hdf5_header_type(wfnq0, sheader, iflavor, kpq, gvecq, symsq, crysq)
#endif
  else
    if(peinf%inode == 0) call open_file(26,file=wfnq0,form='unformatted',status='old')
    call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, symsq, crysq, &
      warn = .false., dont_warn_kgrid = .true.)
  end if

! When using the velocity operator and q shift in a truncated direction,
! WFNq_fi will have k-points with nonzero coordinates in that direction, so
! we should not issue a warning in that case. --DAS
!  call check_trunc_kpts(xct%icutv, kpq)
    
  call check_header('WFN_fi', kp, gvec, syms, crys, wfnq0, kpq, gvecq, symsq, crysq, is_wfn = .true.)

  if(any(kp%kgrid(1:3) /= kpq%kgrid(1:3))) then
    if(peinf%inode == 0) then
      write(0,*) 'WFN_fi  kgrid = ', kp%kgrid(1:3)
      write(0,*) 'WFNq_fi kgrid = ', kpq%kgrid(1:3)
    endif
    call die('kgrids for WFN_fi and WFNq_fi must be the same', only_root_writes = .true.)
  endif

  SAFE_ALLOCATE(kpq%elda, (kpq%mnband, kpq%nrk, kpq%nspin))
  kpq%el(:,:,:) = kpq%el(:,:,:) - xct%avgpot / ryd
  kpq%elda(1:kpq%mnband, 1:kpq%nrk, 1:kpq%nspin) = kpq%el(1:kpq%mnband, 1:kpq%nrk, 1:kpq%nspin)
  call scissors_shift(kpq, eqp%scis, eqp%spl_tck)

!-----------------------------------------------------------------------
! If quasi-particle correction requested, read the corrected
! qp energies from file (in eV)

  minband = 1
  bands_cor = kpq%mnband
  if(xct%eqp_corrections) then
    if (flag%opr .eq. 0) then !eqp_q.dat
      bands_cor = maxval(kpq%ifmax(:,:))
    else !eqp.dat
      bands_cor = maxval(kpq%ifmax(:,:)) + xct%ncb_fi
    endif
    bands_cor = min(bands_cor, kpq%mnband)
    minband = minval(kpq%ifmax(:,:)-xct%nvb_fi+1)
    ! FIXME: for metals this is asking for a few more bands than actually needed on some k-points
    call eqpcor(fncor,peinf%inode,peinf%npes,kpq, &
      minband,bands_cor,0,0,kpq%el,kpq%el,kpq%el,1,0)
  endif

  if(xct%unrestricted_transf) then
    call find_efermi(xct%rfermi, xct%efermi, xct%efermi_input, kpq, bands_cor, minband, &
      "fine grid", should_search = .true., should_update = .true., write7 = .false.,  dont_die_consistency = .true.)
  else
    call find_efermi(xct%rfermi, xct%efermi, xct%efermi_input, kpq, bands_cor, minband, &
      "fine grid", should_search = .true., should_update = .true., write7 = .false.,  dont_die_consistency = .false.)
  endif

  if( xct%use_wfn_hdf5 ) then
    ! I assume nothing to be done here
#ifdef HDF5
    ! call read_hdf5_gvectors()
#endif
  else
    call read_binary_gvectors(26, gvec%ng, gvec%ng, dummygvec, dont_read = .true.)
  end if

  if(any(kpq%ifmax(:,:) == 0)) & 
    call die("BSE codes cannot handle a system where some k-points have no occupied bands.", only_root_writes = .true.) 
 
  kpq%nvband=minval(kpq%ifmax(:,:)-kpq%ifmin(:,:))+1
  kpq%ncband=kpq%mnband-maxval(kpq%ifmax(:,:))

!----------------------------------
! Manual override of band numbers
!
  if (xct%vmax.ne.0) then
    kpq%nvband = xct%vmax-xct%vmin+1
    kpq%ncband = kpq%mnband-xct%vmax
    if (peinf%inode.eq.0) then
      write(6,*)
      write(6,*) '*** Overwrite min/max occupied state for fine grid wavefunctions'
      write(6,*) '*** kpq%nvband =',kpq%nvband,' kpq%ncband =',kpq%ncband
      write(6,*)
    endif
  endif

  ! Define the actual number of conduction and valence bands used.
  ! These include the whole set of bands if it was requested.
  if (read_all_bands_) then
    nvb = kpq%mnband
    ncb = kpq%mnband
  else
    nvb = xct%nvb_fi
    ncb = xct%ncb_fi
  end if

!----------------------------------------------------------------
! (gsm) check whether the requested number of bands
!       is available in the wavefunction file

! we only use the valence bands from WFNq_fi
! if flag%vm.eq.1.and.flag%dtm.eq.1 we don`t need them at all
! subroutine input_q is only called if flag%vm.ne.1.or.flag%dtm.ne.1

  if (xct%nvb_fi .gt. kpq%nvband .and. .not. read_all_bands_) then
    call die("The requested number of valence bands is not available in " // TRUNC(wfnq0) // ".")
  endif

! DAS: degenerate subspace check

  ! degeneracy does not matter for WFNq_fi in inteqp
  if (peinf%inode.eq.0 .and. xct%is_absorption) then
    broken_degeneracy = .false.
    do jj = 1, kpq%nspin
      do ii = 1, kpq%nrk
        if(kpq%ifmax(ii, jj) - xct%nvb_fi > 0) then
          ! no need to compare against band 0 if all valence are included
          if(abs(kpq%elda(kpq%ifmax(ii, jj) - xct%nvb_fi + 1, ii, jj) &
            - kpq%elda(kpq%ifmax(ii, jj) - xct%nvb_fi, ii, jj)) .lt. TOL_Degeneracy) then
            broken_degeneracy = .true.
          endif
        endif
      enddo
    enddo

    if(broken_degeneracy) then
      if(xct%degeneracy_check_override) then
        write(0,'(a)') &
          "WARNING: Selected number of valence bands breaks degenerate subspace in " // TRUNC(wfnq0) // &
          ". Run degeneracy_check.x for allowable numbers."
        write(0,*)
      else
        write(0,'(a)') &
          "Run degeneracy_check.x for allowable numbers, or use keyword " // &
          "degeneracy_check_override to run anyway (at your peril!)."
        call die("Selected number of valence bands breaks degenerate subspace in " // TRUNC(wfnq0) // ".")
      endif
    endif
  endif

!-----------------------------------------------------------------------
!     Define shifted irreducible BZ, kgq%r, and define the shift vector
!     (if not done before). Make sure it is right!
!
  kgq%nr=kpq%nrk
  SAFE_ALLOCATE(kgq%r, (3,kgq%nr))
  kgq%r(1:3,1:kgq%nr)=kpq%rk(1:3,1:kpq%nrk)
  xct%qshift= sqrt( DOT_PRODUCT(xct%shift(:) + xct%finiteq(:), &
    MATMUL(crys%bdot,xct%shift(:) +xct%finiteq(:))) )
  if (.not. xct%read_kpoints .and. abs(xct%qshift).lt.TOL_Zero) then
    xct%shift(:)=kgq%r(:,1)-kg%r(:,1)
    xct%qshift= sqrt( DOT_PRODUCT(xct%shift(:) +xct%finiteq(:), &
      MATMUL(crys%bdot,xct%shift(:) +xct%finiteq(:))) )
  endif
  ! We only need q_shift for the velocity operator
  if(flag%opr == 0) then
    if(peinf%inode.eq.0) write(6,90) xct%shift(:),xct%qshift
90 format(1x,'Shift vector : ',3f9.5,2x,'Length =',f8.5,/)
    if (xct%qflag.ne.1) then
      if(peinf%inode.eq.0) write(6,91) xct%finiteq(:),xct%qshift
91 format(1x,'Exciton momentum : ',3f9.5,2x,'Length =',f8.5,/)
    endif

    if(xct%qshift < TOL_Small) call die("q-shift vector may not be zero.", only_root_writes = .true.)
  endif

!-----------------------------------------------------------------------
!     Define the polarization vector (used only with momentum operator and absorption)
!
  xct%lpol=sqrt(DOT_PRODUCT(xct%pol,MATMUL(crys%bdot,xct%pol)))
  if (abs(xct%lpol).lt.TOL_Zero) then
    xct%pol = xct%shift(:)
    xct%lpol=sqrt(DOT_PRODUCT(xct%pol,MATMUL(crys%bdot,xct%pol)))
  endif
  if(flag%opr==1 .and. xct%is_absorption .and. xct%npol==1) then
    if(peinf%inode.eq.0) write(6,92) xct%pol(:),xct%lpol
92 format(1x,'Polarization vector : ',3f9.5,2x,'Length =',f8.5,/)
    if(xct%lpol < TOL_Small) call die("Polarization vector may not be zero.", only_root_writes = .true.)
  endif
  
!-----------------------------------------------------------------------
!     Generate full Brillouin zone from irreducible wedge, rk -> fk
!
  if (flag%bzq.eq.1 .and. .not. xct%is_absorption) then
    ! in the inteqp code, we want to leave the k-points alone, if symmetries are not being used.
    call fullbz(crys,symsq,kgq,1,skip_checkbz,wigner_seitz=.false.,paranoid=.false.,do_nothing=.true.)
  else if (flag%bzq.eq.1) then
    call fullbz(crys,symsq,kgq,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  else
    call fullbz(crys,symsq,kgq,symsq%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  endif
  xct%nkptq_fi = kgq%nf
  tmpfn=wfnq0
  if (.not. skip_checkbz .and. .not.xct%patched_sampling) then
    call checkbz(kgq%nf,kgq%f,kpq%kgrid,kpq%shift,crys%bdot, &
      tmpfn,'k',.true.,xct%freplacebz,xct%fwritebz)
  endif

  if (flag%bzq.eq.0.and.peinf%inode.eq.0) write(6,801)
  if (flag%bzq.eq.1.and.peinf%inode.eq.0) write(6,802)
801 format(1x,'Using symmetries to expand the shifted-grid sampling')
802 format(1x,'No symmetries used in the shifted-grid sampling')

  SAFE_ALLOCATE(xct%ifmaxq, (kgq%nf,xct%nspin))
  xct%ifmaxq(:,:)=kpq%ifmax(kgq%indr(:),:)

!
!-----------------------------------------------------------------------
!     Find correspondence with fk from WFN_fi
!
!     indexq : correspondence between a k-point in the full BZ, kg%f, and
!       its shifted vector, in kgq%f. xct%finiteq is subtracted from 
!       kgq%f to account for center-of-mass momentum
!     tol : tolerance
!
  tol = TOL_Small

  do ik=1,kg%nf
    ikq=0
    delta=0.1d0
    do while((delta.gt.tol).and.(ikq.lt.kgq%nf))
      ikq=ikq+1      
      qq(:) = kg%f(:,ik)-(kgq%f(:,ikq)-xct%shift(:)-xct%finiteq(:))
      do kk=1,3
        qq(kk) = qq(kk) - anint( qq(kk) )
      enddo
      delta=sqrt(sum(qq(1:3)**2))
    enddo
    ! With patched sampling, k+Q might fall outside the patch 
    if(delta.gt.tol .and. xct%patched_sampling) then
      if(peinf%inode.eq.0) then
        write(6,*) '  Could not find point equivalent to ', (kg%f(ii,ik),ii=1,3)
        write(6,*) '  Skipping point.'
      endif
      xct%indexq_fi(ik)= 0
    elseif (delta.gt.tol) then
      if(peinf%inode.eq.0) then
        write(0,*) 'Could not find point equivalent to ', (kg%f(ii,ik),ii=1,3)
      endif
      call die('k-point mismatch between WFN_fi and ' // TRUNC(wfnq0), only_root_writes = .true.)
    else
      !
      !     make sure that kgq%f(:,ikq)-kg%f(:,ik) = shift vector
      !     Near the zone edge, they may differ by a lattice vector
      !
      do jj=1,3
        ii = nint( kgq%f(jj,ikq)-kg%f(jj,ik) )
        kgq%f(jj,ikq) = kgq%f(jj,ikq) - dble(ii)
        kgq%kg0(jj,ikq) = kgq%kg0(jj,ikq) - ii
      enddo
      xct%indexq_fi(ik)=ikq
      indexq(ik) = ikq
    endif
  enddo

  if(kgq%nf /= kg%nf) then
    if(peinf%inode == 0) write(0,'(a,i7,a,i7,a)') trim(wfnq0)//' unfolds to ', kgq%nf, &
      ' k-points; WFN_fi unfolds to ', kg%nf, ' k-points.'
    call die("WFNq_fi and WFN_fi must have the same number of k-points in the unfolded BZ.", only_root_writes = .true.)
  endif
  ! Each WFN_fi point only needs one WFNq_fi point. In principle it would be ok to have extra unused WFNq_fi points, but in fact
  ! the subsequent code assumes that there are the same number of points, so we will check that here rather than allow mysterious
  ! segfaults or incorrect results later. --DAS

!-----------------------------------------------------------------------
!     Read the wavefunctions and create INT_VWFNQ_*
!

! JRD: Debugging
!      if (peinf%inode .eq. 0) then
!         write(6,*) 'Creating INT_VWFNQ_* files'
!      end if

  wfnv%nband=nvb
  wfnv%nspin=kpq%nspin
  wfnv%nspinor=kpq%nspinor
  
  SAFE_ALLOCATE(intwfnv%ng, (peinf%ikt(peinf%inode+1)))
  SAFE_ALLOCATE(intwfnv%isort, (gvec%ng,peinf%ikt(peinf%inode+1)))
  SAFE_ALLOCATE(intwfnv%cgk, (kpq%ngkmax,wfnv%nband,kpq%nspin*kpq%nspinor,peinf%ikt(peinf%inode+1)))
  intwfnv%nspinor=wfnv%nspinor

  if ( xct%use_wfn_hdf5  ) then
    if (peinf%inode==0) write(6,*)
    if (peinf%inode==0) write(6,'(a)') ' Reading HDF5 wavefuntion '//TRIM(wfnq0)
    ngktot = SUM(kpq%ngk)
    SAFE_ALLOCATE( gvec_kpt_components_all, (3,ngktot) )
#ifdef HDF5
    call read_hdf5_wfn_gvectors(wfnq0, gvec_kpt_components_all, ngktot)
#endif

    if ( read_all_bands_ ) then
      iband_min = 1
      iband_max = kpq%mnband
    else
      iband_min = MINVAL(kpq%ifmax(:,:)) - xct%nvb_fi + 1
      iband_min = MAX(iband_min, 1)
      iband_max = MAXVAL(kpq%ifmax(:,:))
      iband_max = MIN(iband_max, kpq%mnband)
    end if
    tot_nbands     = iband_max - iband_min + 1
    band_block     = (tot_nbands + peinf%npes - 1) / peinf%npes
    ! read at least 128 Mb per MPI task
    min_band_block = 128.0D+00 / ( dble(ngktot) * dble(kpq%nspin*kpq%nspinor) * 16.0D+00 /1024.0D+00/1024D+00 )
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
    !XXX write(1000+peinf%inode,*) peinf%inode, iband_min, iband_max, tot_nbands, band_block, min_band_block, ib_first, ib_last, ib_size
    !XXX flush(1000)
#ifdef HDF5
#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, mpierr)
    call h5fopen_f(wfnq0, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id,error)
    call MPI_AllReduce(MPI_IN_PLACE, ib_size_array, peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
    call h5fopen_f(wfnq0, H5F_ACC_RDONLY_F, file_id, error)
#endif
#endif

    ib_size_max = MAXVAL( ib_size_array )
    SAFE_ALLOCATE( my_cg_all, (ngktot, kpq%nspin*kpq%nspinor, MAX(ib_size,1)) )
    SAFE_ALLOCATE( ipe_cg_all, (ngktot, kpq%nspin*kpq%nspinor, ib_size_max) )
    my_cg_all = ZERO

#ifdef HDF5
    ! dummy allocation
    SAFE_ALLOCATE( peinf%does_it_ownv, (1, 1) )
    peinf%does_it_ownv = .true.
    call read_hdf5_bands_block(file_id, kpq, ib_size_max, ib_size, &
                               peinf%does_it_ownv, ib_first, my_cg_all, is_bse = .true.)

    SAFE_DEALLOCATE_P(peinf%does_it_ownv)

    call h5fclose_f(file_id, error)
#endif

  end if

  SAFE_ALLOCATE(wfnv%isort, (gvec%ng))
  SAFE_ALLOCATE(isend, (peinf%npes))

  write(tmpstr,'(3a)') 'reading wavefunctions (', trim(wfnq0), ')'
  if ( xct%use_wfn_hdf5 ) then
    call progress_init(prog_info, trim(tmpstr), 'MPI task', peinf%npes)
  else
    call progress_init(prog_info, trim(tmpstr), 'k-point', kpq%nrk)
  end if

npes_hdf5 = 0
if ( xct%use_wfn_hdf5 ) npes_hdf5 = peinf%npes-1

! loop over process only in case of wfn_hdf5
do ipe = 0, npes_hdf5
  if ( xct%use_wfn_hdf5 ) then
    call progress_step(prog_info, ipe+1)
    if ( ib_size_array(ipe+1) > 0) then
     ipe_cg_all = ZERO
     if ( ipe == peinf%inode ) then
       ipe_cg_all(:, :, 1:ib_size) = my_cg_all(:, :, 1:ib_size)
     end if
#ifdef MPI
       call MPI_Bcast(ipe_cg_all, ngktot * kpq%nspin*kpq%nspinor * ib_size_max, &
                      MPI_SCALAR, ipe, MPI_COMM_WORLD, mpierr)
#endif
    else
      ! no elements on this MPI task, cycle
      cycle
    end if
  end if

  wfnv%isort=0
  last_ng=-1
  last_ng_match=-1
  last_ikt=-1
  istart = 1
  do irkq = 1, kpq%nrk
    if ( .not. xct%use_wfn_hdf5 ) call progress_step(prog_info, irkq)
    irkq_match = .false.
    ! Find if a point corresponding to the reduced point, irkq, read from WFNq_fi exists on kg
    do ii=1,kg%nf
      if (xct%qflag.eq.1) then
        if (irkq == kgq%indr(indexq(ii))) then
          irkq_match = .true.
          exit
        endif
      else
        if(xct%indexq_fi(ii) .eq.0) cycle
        if (irkq == kgq%indr(xct%indexq_fi(ii))) then
          irkq_match = .true.
          exit
        endif
      endif
    enddo

    wfnv%ng = kpq%ngk(irkq)

! FHJ: Realloc arrays. Note that we can`t do something like
!      "if (wfnv%ng>last_ng)"  b/c fortran complains at read_binary_gvectors 
!      if the vectors are not exactly wfnv%ng big.
    if(wfnv%ng/=last_ng) then
      if(last_ng/=-1) then
        SAFE_DEALLOCATE_P(gvec_kpt%components)
        SAFE_DEALLOCATE(cg)
      endif
      SAFE_ALLOCATE(gvec_kpt%components, (3, wfnv%ng))
      SAFE_ALLOCATE(cg, (wfnv%ng, kpq%nspin*kpq%nspinor))
      last_ng = wfnv%ng
    endif

    if ( xct%use_wfn_hdf5 ) then
       gvec_kpt%components(:,:) = gvec_kpt_components_all(1:3, istart:istart+kpq%ngk(irkq)-1)
    else
       call read_binary_gvectors(26, wfnv%ng, wfnv%ng, gvec_kpt%components)
    end if

!    Skip this k-point if there is no k-point in kg%f that
!    corresponds to it
    if(irkq_match) then
      do ii = 1, kpq%ngk(irkq)
        call findvector(wfnv%isort(ii), gvec_kpt%components(:, ii), gvec)
        if(wfnv%isort(ii) == 0) call die('input_q: could not find gvec')
      enddo

! FHJ: Realloc arrays.
      if(wfnv%ng/=last_ng_match) then
        if(last_ng_match/=-1) then
          SAFE_DEALLOCATE_P(wfnv%cg)
          SAFE_DEALLOCATE(cgarray)
        endif
        SAFE_ALLOCATE(wfnv%cg, (wfnv%ng, wfnv%nband, wfnv%nspin*wfnv%nspinor))
        SAFE_ALLOCATE(cgarray, (wfnv%ng))
        last_ng_match = wfnv%ng
      endif        
      if(peinf%ikt(peinf%inode+1)/=last_ikt) then
        if(last_ikt/=-1) then
          SAFE_DEALLOCATE(iwriteik)
        endif
        SAFE_ALLOCATE(iwriteik,(peinf%ikt(peinf%inode+1)))
        last_ikt = peinf%ikt(peinf%inode+1)
      endif

!       Determine which PEs will write the valence bands for this k-point
      iwrite=0
      iwritetotal=0
      iwriteik=0
      do ii=1, peinf%ikt(peinf%inode+1) ! loop over # of k-points belonging to proc.
        ! Write wfn if the point being read corresponds to a point on kg that belongs to the proc.
        if (xct%indexq_fi(peinf%ik(peinf%inode+1,ii)).eq.0) cycle
        if(kgq%indr(xct%indexq_fi(peinf%ik(peinf%inode+1,ii))) == irkq) then
          iwritetotal=iwritetotal+1
          iwriteik(iwritetotal)=ii
          iwrite=1
        endif
      enddo

      if ( .not. xct%use_wfn_hdf5 ) then 
!       Determine to which PEs the valence bands for this k-point
!       need to be sent...
        isend=0
        if(peinf%inode.eq.0) then
          do jj=2,peinf%npes
            do ii=1, peinf%ikt(jj)
              if (xct%indexq_fi(peinf%ik(jj,ii)).eq.0) cycle
              if(kgq%indr(xct%indexq_fi(peinf%ik(jj,ii))) == irkq) then
                isend(jj)=1
                exit
              endif
            enddo
          enddo
        endif
!     
      end if ! hdf5
    endif ! if(irkq_match)
!
!       Loop over the bands
!
    ii_start = 1
    ii_end   = kpq%mnband
    if ( xct%use_wfn_hdf5 ) then
      ii_start = iband_min + band_block * ipe
      ii_end   = min(ii_start + band_block - 1, iband_max)
    end if
    ii_loc = 0
    do ii = ii_start, ii_end
      ii_loc = ii_loc + 1

      ! Skip reading the WFNs if this is not a band we want.
      dont_read = (((ii<=minval(kpq%ifmax(irkq,:))-xct%nvb_fi) .or. &
                   (ii>maxval(kpq%ifmax(irkq,:)))) .and. &
                    .not. read_all_bands_)

      if ( xct%use_wfn_hdf5 ) then
        cg = ZERO
        if ( .not. dont_read ) then
          cg (:,:) = ipe_cg_all(istart:istart+kpq%ngk(irkq)-1, 1:kpq%nspin*kpq%nspinor, ii_loc)
        end if
      else
        call read_binary_data(26, kpq%ngk(irkq), kpq%ngk(irkq), &
          kpq%nspin*kpq%nspinor, cg, dont_read=dont_read, bcast=.false.)
      end if

      if(.not. irkq_match) cycle

      do is=1, kpq%nspin

        ! Skip bands out of range for this spin polarization
        if (.not. (ii.le.kpq%ifmax(irkq,is).and. &
            ii.gt.kpq%ifmax(irkq,is)-xct%nvb_fi) .and. &
            .not. read_all_bands_) then
          cycle
        end if

          if (read_all_bands_) then
            ! Counting up
            ib = ii
          else
            ! Counting down from VBM
            ib = kpq%ifmax(irkq,is) - ii + 1
          end if

        do ispinor=1, kpq%nspinor
          if ( xct%use_wfn_hdf5 ) then
            cgarray(1:kpq%ngk(irkq))=cg(1:kpq%ngk(irkq), is*ispinor)
          else
            if (peinf%inode.eq.0) then
              cgarray(1:kpq%ngk(irkq))=cg(1:kpq%ngk(irkq), is*ispinor)
              if (peinf%verb_debug) then
                write(6,'(a,3i0,2(f18.13))') 'input_q', irkq, ii, is*ispinor, cgarray(1)
              endif
            endif
#ifdef MPI
            if(peinf%inode.eq.0) then
              do jj=2,peinf%npes
                if(isend(jj).eq.1) then
                  dest=jj-1
                  tag=1000+dest*2-(is-1)
                  call MPI_SEND(cgarray,kpq%ngk(irkq),MPI_SCALAR, &
                    dest,tag,MPI_COMM_WORLD,mpierr)
                endif
              enddo
            else
              if(iwrite.eq.1) then
                tag=1000+peinf%inode*2-(is-1)
                call MPI_RECV(cgarray,kpq%ngk(irkq),MPI_SCALAR,0,tag, &
                  MPI_COMM_WORLD,mpistatus,mpierr)
              endif
            endif
#endif
          end if ! hdf5
          if(iwrite.eq.1) then
            wfnv%cg(1:wfnv%ng,ib,is*ispinor) = cgarray
          endif
        enddo
        if(iwrite.eq.1) then
          !FHJ: We only Send/Recv and check one spin at a time
          call checknorm(wfnq0,ii,irkq,kpq%ngk(irkq),is,kpq%nspinor,&
                       wfnv%cg(:,ib,:))
        endif
      end do
    enddo

    istart = istart + kpq%ngk(irkq)
    if(.not. irkq_match) cycle
    
    if(iwrite.eq.1) then
      do ijk = 1, iwritetotal
        intwfnv%ng(iwriteik(ijk))=wfnv%ng
        intwfnv%isort(:,iwriteik(ijk))=wfnv%isort(:)

        if ( xct%use_wfn_hdf5 ) then
          do ii = ii_start, ii_end
            do is=1, kpq%nspin
              ! Skip bands out of range for this spin polarization
              if (.not. (ii.le.kpq%ifmax(irkq,is).and. &
                  ii.gt.kpq%ifmax(irkq,is)-xct%nvb_fi) .and. &
                  .not. read_all_bands_) then
                cycle
              end if
              if (read_all_bands_) then
                ! Counting up
                ib = ii
              else
                ! Counting down from VBM
                ib = kpq%ifmax(irkq,is) - ii + 1
              end if
              do ispinor=1, kpq%nspinor
                intwfnv%cgk(1:wfnv%ng,ib,is*ispinor,iwriteik(ijk))=wfnv%cg(1:wfnv%ng,ib,is*ispinor)
              end do
            end do
          end do
        else
          intwfnv%cgk(1:wfnv%ng,:,:,iwriteik(ijk))=wfnv%cg(1:wfnv%ng,:,:)
        end if

      enddo
    endif
  enddo

  if(last_ikt/=-1) then
    SAFE_DEALLOCATE(iwriteik)
  endif
  if(last_ng_match/=-1) then
    SAFE_DEALLOCATE_P(wfnv%cg)
    SAFE_DEALLOCATE(cgarray)
  endif
  if(last_ng/=-1) then
    SAFE_DEALLOCATE_P(gvec_kpt%components)
    SAFE_DEALLOCATE(cg)
  endif

end do ! ipe loop for HDF5 case
  call progress_free(prog_info)

  SAFE_DEALLOCATE(isend)
  SAFE_DEALLOCATE_P(wfnv%isort)

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
  !XX   SAFE_DEALLOCATE_P(wfnv%cg)
  !XX   SAFE_DEALLOCATE(cgarray)
  !XX endif
  !XX if(last_ng/=-1) then
  !XX   SAFE_DEALLOCATE_P(gvec_kpt%components)
  !XX   SAFE_DEALLOCATE(cg)
  !XX endif
      
!
! End loop over k-points
!--------------------------------------------------------------------------------
           
! JRD: For Finite Q we need valence energies
! FHJ: loop thru indices on `eqp` grid, and then find the
!      corresponding labels on the shifted grid

  if (read_all_bands_) then
    eqp%band_ordering = 1
  else
    eqp%band_ordering = 0
  end if

  do is=1,xct%nspin
    do irk=1,xct%nkpt_fi
      !irkq = index of q-pt on shifted grid
      irkq = kgq%indr(indexq(irk))

      if (.not. read_all_bands_) then

        !loop thru bands of shifted grid
        do ibq=1, kpq%mnband

          !ib = band index on `eqp` grid
          ib = kpq%ifmax(irkq,is) - ibq + 1

          if ( (ib > 0) .and. (ib <= xct%nvb_fi) ) then
            eqp%evqp (ib, irk, is) = kpq%el(ibq, irkq, is)
            eqp%evlda(ib, irk, is) = kpq%elda(ibq, irkq, is) 
          endif

        enddo
      endif
    enddo
  enddo

!-----------------------------
! Deallocate

  SAFE_DEALLOCATE(kpq%ifmin)
  SAFE_DEALLOCATE(kpq%ifmax)
  SAFE_DEALLOCATE(kpq%rk)
  SAFE_DEALLOCATE(kpq%el)
  SAFE_DEALLOCATE(kpq%elda)


! JRD: Debugging
!      if (peinf%inode .eq. 0) then
!         write(6,*) 'Deallocated arrays'
!      end if

!-----------------------------
! Write out info about xtal

  if(peinf%inode.eq.0) then
    write(6,'(/,1x,3a)') 'Valence wavefunctions read from file ', TRUNC(wfnq0), ':'
    write(6,'(1x,a,i0)') '- Number of k-points in irreducible BZ: ', kgq%nr
    write(6,'(1x,a,i0)') '- Number of k-points in full BZ: ', kgq%nf
    if (peinf%verb_high) then
      write(6,'(1x,a)') '- Listing all k-points:'
      write(6,'(1(2x,3(1x,f10.6)))') (kgq%r(:,jj), jj=1,kgq%nr)
    endif
    if ( .not. xct%use_wfn_hdf5 ) call close_file(26)
  endif ! node 0

  POP_SUB(input_q)

  return
end subroutine input_q

end module input_q_m
