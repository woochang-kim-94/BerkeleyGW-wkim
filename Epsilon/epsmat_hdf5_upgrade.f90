!=============================================================================
!
!  Utilities:
!
!    epsmat_hdf5_upgrade  Originally by FHJ    Last Modified 08/2016 (FHJ)
!
!     This utility upgrades an epsmat.h5 file from BerkeleyGW-1.1-beta
!     (version 1) to the version used in BerkeleyGW-1.2.0 (version 3).
!
!=============================================================================

#include "f_defs.h"

program epsmat_hdf5_upgrade

  use global_m
  use hdf5
  use hdf5_io_m
  use write_matrix_m
  use epswrite_hdf5_m
  use wfn_rho_vxc_io_m
  use wfn_io_hdf5_m
  use input_utils_m

  implicit none

  character(len=128) :: fname_eps, fname_wfn, tmpstr
  integer :: nspin, ng
  integer :: version, intinfo(9)
  integer :: iw, narg, error
  integer, allocatable :: gvecs(:,:), isorts(:,:,:)
  logical, allocatable :: qpt_done(:)
  real(DP), allocatable :: freqs_tmp(:,:), freqgrid(:)
  logical :: skip_wfn, mf_header_exists, nspin_exists, exists
  integer(HID_T) :: file_id
  type(mf_header_t) :: mf
  type(polarizability) :: pol

  narg = command_argument_count()
  if (narg/=4) then
    write(0,*)
    write(0,*) 'Usage: epsmat_hdf5_upgrade.flavor.x epsmat.h5 WFN icutv [itruncval]'
    write(0,*) '==================================================================='
    write(0,*)
    write(0,*) 'This utility will upgrade an epsmat.h5 file from version 1'
!#BEGIN_INTERNAL_ONLY
#if 1
    write(0,*) '(generated prior to BerkeleyGW r6583 / BerkeleyGW 1.2.0 stable)'

#else
!#END_INTERNAL_ONLY
    write(0,*) '(generated prior to BerkeleyGW 1.2.0)'
!#BEGIN_INTERNAL_ONLY
#endif
    write(0,*) 'to the version 3, currently used in BerkeleyGW.'
!#END_INTERNAL_ONLY
    write(0,*)
    write(0,*) 'Required arguments:'
    write(0,*) '-------------------'
    write(0,*) '  - epsmat.h5: epsmat HDF5 file in previous (beta) version'
    write(0,*) '  - WFN: WFN file from which we`ll extract the header'
!#BEGIN_INTERNAL_ONLY
#if 1
    write(0,*) '      You may also specify a single dash (-) as an argument if your'
    write(0,*) '      epsmat.h5 was generated after r6358 (~Aug 2014). This is'
    write(0,*) '      the case for BerkeleyGW 1.1-beta2, but not for BerkeleyGW 1.1-beta!'
#else
!#END_INTERNAL_ONLY
    write(0,*) '      You may also specify a single dash (-) as an argument, but only'
    write(0,*) '      if you generated your epsmat.h5 with BerkeleyGW 1.1-beta2, and not'
    write(0,*) '      BerkeleyGW 1.1-beta!'
!#BEGIN_INTERNAL_ONLY
#endif
!#END_INTERNAL_ONLY
    write(0,*) '  - icutv: truncation flag used in the calculation. Options are:'
    write(0,*) '      0 -> no truncation'
    write(0,*) '      2 -> spherical_truncation'
    write(0,*) '      4 -> cell_wire_truncation'
    write(0,*) '      5 -> cell_box_truncation'
    write(0,*) '      6 -> cell_slab_truncation'
    write(0,*) '  - itruncval: truncation radius, only needed if icutv==2.'
    write(0,*)
    write(0,*) 'Warning:'
    write(0,*) '--------'
    write(0,*) '  The upgrade process happens in place and modifies the source file.'
    write(0,*) '  The process is instantaneous, but because there is always a chance that'
    write(0,*) '  something will go wrong, BACK UP YOUR EPSMAT.H5 FILE FIRST!'
    write(0,*)
    stop
  endif
  call get_command_argument(1, fname_eps)
  call get_command_argument(2, fname_wfn)
  call get_command_argument(3, tmpstr)
  read(tmpstr,*) pol%icutv
  call get_command_argument(4, tmpstr)
  read(tmpstr,*) pol%truncval(1)

  write(6,*)
  write(6,*) 'Epsmat HDF5 Upgrade Utility'
  write(6,*) '==========================='
  write(6,*)
  call h5open_f(error)

  skip_wfn = TRUNC(fname_wfn)=="-"
  ! FHJ: Read WFN file, get the header
  if (skip_wfn) then
    write(6,*) 'We will not read any WFN file.'
  else
    write(6,*) 'Reading WFN file '//TRUNC(fname_wfn)//'.'
!#BEGIN_INTERNAL_ONLY
    if (index(fname_wfn,'.h5')/=0) then
      call read_hdf5_mf_header(TRUNC(fname_wfn), mf, sheader='WFN', iflavor=SCALARSIZE, &
        warn=.false., dont_warn_kgrid=.true.)
      SAFE_ALLOCATE(mf%gvec%components, (3, mf%gvec%ng))
      call read_hdf5_gvectors(TRUNC(fname_wfn), mf%gvec%ng, mf%gvec%components)
    else
!#END_INTERNAL_ONLY
      call open_file(unit=25,file=TRUNC(fname_wfn),form='unformatted',status='old')
      call read_mf_header(25, mf, sheader='WFN', iflavor=SCALARSIZE, &
        warn=.false., dont_warn_kgrid=.true.)
      SAFE_ALLOCATE(mf%gvec%components, (3, mf%gvec%ng))
      call read_binary_gvectors(25, mf%gvec%ng, mf%gvec%ng, mf%gvec%components)
      call close_file(25)
!#BEGIN_INTERNAL_ONLY
    endif
!#END_INTERNAL_ONLY
  endif

  write(6,*) 'Reading epsmat file '//TRUNC(fname_eps)//'.'
  call open_file(99, trim(fname_eps), status='old')
  call close_file(99)
  call h5fopen_f(trim(fname_eps), H5F_ACC_RDONLY_F, file_id, error)
  write(6,*) 'Performing consistency checks.'
  version = -1
  call hdf5_read_int(file_id, 'versionnumber', version, error)
  if (version/=1.or.error/=0) then
    write(0,*)
    write(0,*) 'ERROR: source file ',trim(fname_eps),' has an incorrect version tag.'
    write(0,*) '  expecting: ', 1
    write(0,*) '  read: ', version
    write(0,*) '  error flag: ', error
    write(0,*)
    call die('Invalid file format in '//trim(fname_eps)//'.')
  endif

  !FHJ: mf_header is a recent addition, so it might not be present in older HDF5 file
  call h5lexists_f(file_id, 'mf_header', mf_header_exists, error)
  if (error/=0) call die('Could not search for mf_header.')
  if (mf_header_exists) then
    write(6,*) 'File '//TRUNC(fname_eps)//' contains mean-field information.'
  else
    write(6,*) 'File '//TRUNC(fname_eps)//' does not contain mean-field information.'
    if (skip_wfn) &
      call die('Could not find mf_header in '//trim(fname_eps)//'. You`ll need to specify the WFN argument.')
  endif
  !FHJ: nspin is a recent addition, so it might not be present in older HDF5 file
  call h5lexists_f(file_id, 'nspin', nspin_exists, error)
  if (error/=0) call die('Could not search for nspin.')
  if (nspin_exists) then
    call hdf5_read_int(file_id, 'nspin', nspin)
  else
    if (skip_wfn) &
      call die('Could not find nspin in '//trim(fname_eps)//'. You`ll need to specify the WFN argument.')
  endif

  !FHJ: from the old epsmat.h5.spec file:
  !Dataset: intinfo
  !Rank: 1
  !Dims(1): 9
  !Value(1): 0 for REAL and 1 for CPLX
  !Value(2): # of qpoints
  !Value(3): freq_dep
  !Value(4): # of frequencies
  !Value(5): # of gvectors
  !Value(6): Max of nmtx(q) over q
  !Value(7): kgrid(1)
  !Value(8): kgrid(2)
  !Value(9): kgrid(3)
  call hdf5_read_double(file_id, 'dblinfo', pol%ecuts)
  call hdf5_read_int_array(file_id, 'intinfo', (/9/), intinfo)
  pol%nq = intinfo(2)
  pol%freq_dep = intinfo(3)
  pol%nfreq = intinfo(4)
  ng = intinfo(5)
  if (pol%nfreq<1) pol%nfreq = 1
  if (intinfo(1)+1/=SCALARSIZE) call die('Wrong flavor.')
  if (.not.skip_wfn) then
    ! FHJ: Only perform this consistency check if we are using an external WFN file
    if (nspin_exists) then
      if (nspin/=mf%kp%nspin) &
        call die('nspin mismatch between '//trim(fname_wfn)//' and '//trim(fname_eps)//'.')
    else
      nspin = mf%kp%nspin
    endif
    write(6,*) 'Checking consistency between G-spaces from '//&
      TRUNC(fname_eps)//' and '//TRUNC(fname_wfn)//'.'
    if (ng/=mf%gvec%ng) &
      call die('ng mismatch between '//trim(fname_wfn)//' and '//trim(fname_eps)//'.')
    SAFE_ALLOCATE(gvecs, (3,ng))
    call hdf5_read_int_array(file_id, 'gvecs', (/3,ng/), gvecs)
    if (any(mf%gvec%components/=gvecs)) &
      call die('gvectors mismatch between '//trim(fname_wfn)//' and '//trim(fname_eps)//'.')
    SAFE_DEALLOCATE(gvecs)
  endif
  call h5fclose_f(file_id, error)

  if (.not.mf_header_exists) then
    ! FHJ: Only write mf_header if it wasn`t there in the first place.
    write(6,*) 'Writing mean-field header.'
    call setup_hdf5_mf_file(trim(fname_eps), create_file=.false.)
    !mf%sheader = 'WFN'
    !mf%iflavor = SCALARSIZE
    call write_hdf5_mf_header(trim(fname_eps), mf)
    call write_hdf5_gvectors(trim(fname_eps), ng, mf%gvec%components)
  endif

  call h5fopen_f(trim(fname_eps), H5F_ACC_RDWR_F, file_id, error)
  if (mf_header_exists) then
    call h5lexists_f(file_id, 'mf_header/flavor', exists, error)
    ! FHJ: Some old mf_headers didn`t have the flavor/versionnumber.
    if (.not.exists.or.error/=0) call hdf5_write_int(file_id, 'mf_header/flavor', SCALARSIZE)
    ! FHJ: version hard coded to 1, we don`t want to auto-update it for now.
    call h5lexists_f(file_id, 'mf_header/versionnumber', exists, error)
    if (.not.exists.or.error/=0) call hdf5_write_int(file_id, 'mf_header/versionnumber', VER_WFN_HDF5)
  else
    ! FHJ: Re-read mf-header
    call read_hdf5_mf_header(trim(fname_eps), mf, sheader='WFN', iflavor=SCALARSIZE, &
      warn=.false., dont_warn_kgrid=.true.)
  endif



  write(6,*) 'Setting up new file layout.'
  call hdf5_create_group(file_id, 'eps_header')
  call hdf5_create_group(file_id, 'eps_header/params')
  call hdf5_create_group(file_id, 'eps_header/qpoints')
  call hdf5_create_group(file_id, 'eps_header/freqs')
  call hdf5_create_group(file_id, 'eps_header/gspace')
  call hdf5_create_group(file_id, 'mats')
  call hdf5_write_int(file_id, 'eps_header/versionnumber', VER_EPS_HDF5)
  call hdf5_write_int(file_id, 'eps_header/flavor', SCALARSIZE)

  ! FHJ: some defaults..
  if (index(fname_eps,'epsmat')/=0 .or. index(fname_eps,'eps0mat')/=0) then
    pol%matrix_type = 0
  elseif (index(fname_eps,'chimat')/=0 .or. index(fname_eps,'chi0mat')/=0) then
    pol%matrix_type = 2
  else
    write(0,*)
    write(0,*) 'WARNING: could not determine if we have a chimat or epsmat file.'
    write(0,*) 'Assuming we are dealing with an epsmat file.'
    write(0,*)
    pol%matrix_type = 0
  endif
  call eps_setup_sizes(pol, SCALARSIZE, nspin)
  pol%efermi = 0d0
  pol%intraband_flag = 0
  pol%intraband_overlap_min = 0.9d0
  pol%subsample = .false.

  ! FHJ: General datasets
  write(6,*) 'Rewriting general datasets.'
  call hdf5_write_int(file_id, 'eps_header/params/matrix_type', pol%matrix_type)
  call hdf5_write_logical(file_id, 'eps_header/params/has_advanced', pol%has_advanced)
  call hdf5_write_int(file_id, 'eps_header/params/nmatrix', pol%nmatrix)
  call hdf5_write_int(file_id, 'eps_header/params/matrix_flavor', pol%matrix_flavor)
  call hdf5_write_int(file_id, 'eps_header/params/icutv', pol%icutv)
  call hdf5_write_double(file_id, 'eps_header/params/ecuts', pol%ecuts)
  call hdf5_write_int(file_id, 'eps_header/params/nband', -1)
  call hdf5_write_double(file_id, 'eps_header/params/efermi', pol%efermi/ryd)
  call hdf5_write_int(file_id, 'eps_header/params/intraband_flag', pol%intraband_flag)
  call hdf5_write_double(file_id, 'eps_header/params/intraband_overlap_min', pol%intraband_overlap_min)
  call hdf5_write_logical(file_id, 'eps_header/params/subsample', pol%subsample)

  ! FHJ: Q-points-related datasets
  write(6,*) 'Rewriting q-points-related datasets.'
  call hdf5_write_int(file_id, 'eps_header/qpoints/nq', pol%nq)
  call move_dset('qpoints', 'eps_header/qpoints/qpts')
  call hdf5_write_int_array(file_id, 'eps_header/qpoints/qgrid', (/3/), intinfo(7:9))
  call h5lexists_f(file_id, 'qpt_done', exists, error)
  if (exists) then
    call move_dset('qpt_done', 'eps_header/qpoints/qpt_done')
  else
    SAFE_ALLOCATE(qpt_done, (pol%nq))
    qpt_done(:) = .true.
    call hdf5_write_logical_array(file_id, 'eps_header/qpoints/qpt_done', (/pol%nq/), qpt_done)
    SAFE_DEALLOCATE(qpt_done)
  endif

  ! FHJ: Frequency-related datasets
  write(6,*) 'Rewriting frequency-related datasets.'
  call hdf5_write_int(file_id, 'eps_header/freqs/freq_dep', pol%freq_dep)
  call hdf5_write_int(file_id, 'eps_header/freqs/nfreq', pol%nfreq)
  if (pol%freq_dep==0) then
    call hdf5_write_double_array(file_id, 'eps_header/freqs/freqs', (/2,pol%nfreq/), (/0d0,0d0/))
  else
    SAFE_ALLOCATE(freqgrid, (pol%nfreq))
    SAFE_ALLOCATE(freqs_tmp, (2,pol%nfreq))
    call hdf5_read_double_array(file_id, 'freqs', (/pol%nfreq/), freqgrid)
    call hdf5_read_double_array(file_id, 'freqbrds', (/2,pol%nfreq/), freqs_tmp)
    do iw=1,pol%nfreq
      freqs_tmp(1,iw) = freqs_tmp(1,iw) + freqgrid(iw)
    enddo
    call hdf5_write_double_array(file_id, 'eps_header/freqs/freqs', (/2, pol%nfreq/), freqs_tmp)
    SAFE_DEALLOCATE(freqgrid)
    SAFE_DEALLOCATE(freqs_tmp)
  endif

  ! FHJ: G-vectors-related datasets
  write(6,*) 'Rewriting G-vectors-related datasets.'
  call move_dset('nmtx-of-q', 'eps_header/gspace/nmtx')
  call hdf5_write_int(file_id, 'eps_header/gspace/nmtx_max',  intinfo(6))
  call move_dset('q-gvec-ekin', 'eps_header/gspace/ekin')
  SAFE_ALLOCATE(isorts, (ng,2,pol%nq))
  call hdf5_read_int_array(file_id, 'q-gvec-index', shape(isorts), isorts)
  call hdf5_write_int_array(file_id, 'eps_header/gspace/gind_eps2rho', &
    (/ng,pol%nq/), isorts(:,1,:), error)
  call hdf5_write_int_array(file_id, 'eps_header/gspace/gind_rho2eps', &
    (/ng,pol%nq/), isorts(:,2,:), error)
  SAFE_DEALLOCATE(isorts)
  ! FHJ: Some old mf_headers don`t contain the gvec%components
  call h5lexists_f(file_id, 'mf_header/gspace/components', exists, error)
  if (.not.exists.or.error/=0) call move_dset('gvecs', 'mf_header/gspace/components')

  ! FHJ: Move large matrix datasets
  write(6,*) 'Moving matrices.'
  call move_dset('matrix', 'mats/matrix')
  call move_dset('matrix-diagonal', 'mats/matrix-diagonal')

  write(6,*) 'Cleaning up old structures.'
  call h5ldelete_f(file_id, 'q-gvec-index', error)
  call h5ldelete_f(file_id, 'dblinfo', error)
  call h5ldelete_f(file_id, 'intinfo', error)
  call delete_if_exists('freqs') ! FHJ: Not present in freq_dep==0
  call delete_if_exists('freqbrds') ! Same
  call delete_if_exists('gvecs')
  call delete_if_exists('info')
  call delete_if_exists('nspin')
  call delete_if_exists('versionnumber')

  call h5fclose_f(file_id, error)
  write(6,*) 'All done!'
  write(6,*)
  call h5close_f(error)

contains

  subroutine move_dset(src, dest)
    character(len=*), intent(in) :: src, dest

    PUSH_SUB(move_dset)

    call h5lcreate_hard_f(file_id, src, file_id, dest, error)
    call h5ldelete_f(file_id, src, error)

    POP_SUB(move_dset)

  end subroutine move_dset

  subroutine delete_if_exists(dset)
    character(len=*), intent(in) :: dset

    PUSH_SUB(delete_if_exists)

    exists = .false.
    call h5lexists_f(file_id, dset, exists, error)
    if (exists) call h5ldelete_f(file_id, dset, error)

    POP_SUB(delete_if_exists)

  end subroutine delete_if_exists


end program epsmat_hdf5_upgrade
