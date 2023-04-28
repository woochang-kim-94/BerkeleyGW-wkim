!=========================================================================
!
! Included from file kernel_io.F90.
! You are expected to understand this. --FHJ
!
!=========================================================================

#ifdef READ
  #define READWRITE(x) read ## x
  #define HDF5_READWRITE_INT() hdf5_read_int
  #define HDF5_READWRITE_DOUBLE() hdf5_read_double
  #define HDF5_READWRITE_INT_ARRAY() hdf5_read_int_array
  #define HDF5_READWRITE_DOUBLE_ARRAY() hdf5_read_double_array
  #define INTENT out
  #define OPEN_MODE 'r'
#else
  #define READWRITE(x) write ## x
  #define HDF5_READWRITE_INT() hdf5_write_int
  #define HDF5_READWRITE_DOUBLE() hdf5_write_double
  #define HDF5_READWRITE_INT_ARRAY() hdf5_write_int_array
  #define HDF5_READWRITE_DOUBLE_ARRAY() hdf5_write_double_array
  #define INTENT in
  #define OPEN_MODE 'rw'
#endif

#define READWRITE_KERNELHEADER_HDF5 READWRITE(_kernel_header_hdf5)

!> Defines a subroutine with the template {read,write}_kernel_header_hdf5
!! Note: this routine doesn`t create any groups. So, call
!! setup_kernel_hdf5 before write_kernel_hdf5.
!! Unlike WFN files, you can`t specify the flavor manually when check it,
!! it is always compared against SCALARSIZE.
subroutine READWRITE_KERNELHEADER_HDF5(fname, kern, check_version, check_flavor, rw_mf_header)
  character(len=*), intent(in) :: fname
  type(kernel_header_t), intent(INTENT) :: kern !< kernel_header_t type
  !> If .true. and reading the file, we`ll make sure the file version is correct.
  !! Default is .true.
  logical, intent(in), optional :: check_version
  !> If .true. and reading the file, we`ll make sure the flavor is correct.
  !! Default is .true.
  logical, intent(in), optional :: check_flavor
  !> Should we read or write the mf header stuff? Default is .true.
  logical, intent(in), optional :: rw_mf_header

  integer(HID_T) :: file_id
  logical :: check_version_, check_flavor_, rw_mf_header_, mf_exists
  integer :: error, dims(2)

  PUSH_SUB(READWRITE_KERNELHEADER_HDF5)

  check_version_ = .true.
  if (present(check_version)) check_version_ = check_version
  check_flavor_ = .true.
  if (present(check_flavor)) check_flavor_ = check_flavor
  rw_mf_header_ = .true.
  if (present(rw_mf_header)) rw_mf_header_ = rw_mf_header

#ifdef READ
  call hdf5_open_file(trim(fname), OPEN_MODE, file_id)
  mf_exists = hdf5_exists(file_id, '/mf_header')
  call hdf5_close_file(file_id)
  rw_mf_header_ = rw_mf_header_ .and. mf_exists
  if (rw_mf_header_) then
    call read_hdf5_mf_header(fname, kern%mf)
    SAFE_ALLOCATE(kern%mf%gvec%components, (3,kern%mf%gvec%ng))
    call read_hdf5_gvectors(fname, kern%mf%gvec%ng, kern%mf%gvec%components)
  else
    kern%mf%gvec%components => Null()
  endif
#else
  rw_mf_header_ = rw_mf_header_ .and. associated(kern%mf%gvec%components)
  if (rw_mf_header_) then
    call setup_hdf5_mf_file(fname, create_file=.false.)
    call write_hdf5_mf_header(fname, kern%mf)
    call write_hdf5_gvectors(fname, kern%mf%gvec%ng, kern%mf%gvec%components)
  endif
#endif

  call hdf5_open_file(trim(fname), OPEN_MODE, file_id)

  ! FHJ: Read header stuff. Note that kern%version and kern%iflavor must be set correctly,
  ! for example, by calling xctinfo_to_kernel_header.
  call HDF5_READWRITE_INT()(file_id, 'bse_header/versionnumber', kern%version)
#ifdef READ
  if (check_version_) then
    call hdf5_require_version(file_id, 'bse_header/versionnumber', VER_BSE_HDF5, fname)
  endif
#endif
  call HDF5_READWRITE_INT()(file_id, 'bse_header/flavor', kern%iflavor)
#ifdef READ
  if (check_flavor_) then
    call hdf5_require_flavor(file_id, 'bse_header/flavor', SCALARSIZE, fname)
  endif
#endif

  ! FHJ: GENERAL PARAMS
  call HDF5_READWRITE_INT()(file_id, 'bse_header/params/screening', kern%iscreen)
  call HDF5_READWRITE_INT()(file_id, 'bse_header/params/icutv', kern%icutv)
  call HDF5_READWRITE_DOUBLE()(file_id, 'bse_header/params/ecuts', kern%ecuts)
  call HDF5_READWRITE_DOUBLE()(file_id, 'bse_header/params/ecutg', kern%ecutg)
  call HDF5_READWRITE_DOUBLE()(file_id, 'bse_header/params/efermi', kern%efermi)
  call HDF5_READWRITE_INT()(file_id, 'bse_header/params/theory', kern%theory)
  call HDF5_READWRITE_INT()(file_id, 'bse_header/params/nblocks', kern%nblocks)
  call HDF5_READWRITE_INT()(file_id, 'bse_header/params/storage', kern%storage)
  call HDF5_READWRITE_INT()(file_id, 'bse_header/params/nmat', kern%nmat)
  ! FHJ: BANDS
  call HDF5_READWRITE_INT()(file_id, 'bse_header/bands/nvb', kern%nvb)
  call HDF5_READWRITE_INT()(file_id, 'bse_header/bands/ncb', kern%ncb)
  call HDF5_READWRITE_INT()(file_id, 'bse_header/bands/n1b', kern%n1b)
  call HDF5_READWRITE_INT()(file_id, 'bse_header/bands/n2b', kern%n2b)
  ! FHJ: These are really duplicate from mf_header...
  call HDF5_READWRITE_INT()(file_id, 'bse_header/bands/ns', kern%ns)
  call HDF5_READWRITE_INT()(file_id, 'bse_header/bands/nspinor', kern%nspinor)
  ! FHJ: K-POINTS
  call HDF5_READWRITE_INT()(file_id, 'bse_header/kpoints/nk', kern%nk)
#ifdef READ
  SAFE_ALLOCATE(kern%kpts, (3,kern%nk))
#endif
  call HDF5_READWRITE_DOUBLE_ARRAY()(file_id, 'bse_header/kpoints/kpts', (/3,kern%nk/), kern%kpts)
  call HDF5_READWRITE_INT_ARRAY()(file_id, 'bse_header/kpoints/kgrid', (/3/), kern%kgrid)
  call HDF5_READWRITE_INT()(file_id, 'bse_header/kpoints/qflag', kern%qflag)
  call HDF5_READWRITE_DOUBLE_ARRAY()(file_id, 'bse_header/kpoints/exciton_Q_shift', (/3/), &
                                     kern%exciton_Q_shift)

  call hdf5_close_file(file_id)

  POP_SUB(READWRITE_KERNELHEADER_HDF5)

end subroutine READWRITE_KERNELHEADER_HDF5

! these undefs prevent lots of cpp warnings
#undef READWRITE_KERNELHEADER_HDF5
#undef READWRITE
#undef HDF5_READWRITE_INT
#undef HDF5_READWRITE_DOUBLE
#undef HDF5_READWRITE_INT_ARRAY
#undef HDF5_READWRITE_DOUBLE_ARRAY
#undef INTENT
#undef OPEN_MODE
