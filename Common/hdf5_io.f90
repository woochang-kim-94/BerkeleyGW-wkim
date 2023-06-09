!>=========================================================================
!!
!! Module:
!!
!! hdf5_io_m        Originally by FHJ     Last Modified 08/2019 (FHJ)
!!
!!     High-level routines to read and write data in HDF5 format.
!!     This module imports all functions from hdf5_io_safe_m and
!!     hdf5_io_data_m, so you only need to use this module.
!!
!!     Note that high-level HDF5 routines typically take and optional
!!     `error` argument that should be omitted unless you want to explicitly
!!     check for non-zero exist status. If the `error` optional arguments
!!     are not passed and the underlying HDF5 routine returns a nonzero
!!     error, the code will die instead of continuing (which is the default
!!     behavior of the HDF5 library).
!!
!!=========================================================================

#include "f_defs.h"

module hdf5_io_m
  use, intrinsic :: iso_c_binding
  use global_m
  use hdf5_io_safe_m
  use hdf5_io_data_m
#ifdef HDF5
  use hdf5

  implicit none

  ! FHJ: Note that all subroutines are public in this file, but we only
  ! explicitly export a few of them in mdoule hdf5_io_m

contains

!> Make sure that an hdf5 file has the correct version number.
subroutine hdf5_require_version(loc_id, dset_name, req_version, fname, allow_greater)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in) :: req_version !< version to require
  character(LEN=*), intent(in) :: fname !< file name, for debugging purposes
  !> allow the file to have a version greater than req_version? Defaults to .false.
  logical, intent(in), optional :: allow_greater

  integer :: file_version, errcode
  logical :: allow_greater_

  PUSH_SUB(hdf5_require_version)

  file_version = -1
  call hdf5_read_int(loc_id, dset_name, file_version, errcode)
  allow_greater_ = .false.
  if (present(allow_greater)) allow_greater_ = allow_greater
  if (file_version<req_version .or. &
    (file_version/=req_version.and..not.allow_greater_) .or. errcode/=0) then
    if (peinf%inode==0) then
      write(0,*)
      write(0,*) 'ERROR: Incorrect version in file ', trim(fname),' while reading ',trim(dset_name)
      write(0,*) '       Expecting: ', req_version
      write(0,*) '       Got: ', file_version
      write(0,*) '       Errcode: ', errcode
      write(0,*) 'Your file was probably generated by an older version of BerkeleyGW and'
      write(0,*) 'is now obsolete. Consult the documentation and use the appropriate converter.'
      write(0,*)
    endif
    call die("Wrong version for file '"+trim(fname)+"'.", only_root_writes=.true.)
  endif

  POP_SUB(hdf5_require_version)

end subroutine hdf5_require_version


!> Make sure that an hdf5 file has the correct flavor.
subroutine hdf5_require_flavor(loc_id, dset_name, req_flavor, fname)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in) :: req_flavor !< flavor to require
  character(LEN=*), intent(in) :: fname !< file name, for debugging purposes

  integer :: file_flavor, errcode

  PUSH_SUB(hdf5_require_flavor)

  file_flavor = -1
  call hdf5_read_int(loc_id, dset_name, file_flavor, errcode)
  if (file_flavor/=req_flavor.or.errcode/=0) then
    if (peinf%inode==0) then
      write(0,*)
      write(0,*) 'ERROR: Incorrect flavor in file ', trim(fname), ' while reading ',trim(dset_name)
      write(0,*) '       Expecting: ', req_flavor
      write(0,*) '       Got: ', file_flavor
      write(0,*) '       Errcode: ', errcode
      write(0,*) 'You are probably linking the wrong file or running the BerkeleyGW binary with'
      write(0,*) 'the wrong flavor.'
      write(0,*)
    endif
    call die("Wrong flavor in file "+trim(fname)+"'.", only_root_writes=.true.)
  endif

  POP_SUB(hdf5_require_flavor)

end subroutine hdf5_require_flavor


!> Create an HDF5 file, optionally capturing the error message.
!! Pass paralell_io=.true. if you want to open the file with parallel IO support,
!! and optionally set the communicator `comm` (defaults to MPI_COMM_WORLD)
subroutine hdf5_create_file(filename, file_id, parallel_io, comm, error)
  character(len=*), intent(in) :: filename
  integer(hid_t), intent(out) :: file_id
  logical, optional, intent(in) :: parallel_io
  integer, optional, intent(in) :: comm
  integer, optional, intent(out) :: error

  integer :: errcode, comm_
#ifdef MPI
  integer(HID_T) :: plist_id
#endif
  logical :: parallel_io_

  PUSH_SUB(hdf5_create_file)

  parallel_io_ = .false.
  if (present(parallel_io)) parallel_io_ = parallel_io
#ifdef MPI
  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm
#else
  comm_ = 0
#endif

#ifdef MPI
  if (parallel_io_) then
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, errcode)
    call hdf5_check_error('h5pcreate_f', errcode, error)
    if (present(error)) then; if (error/=0) return; endif

    call h5pset_fapl_mpio_f(plist_id, comm_, MPI_INFO_NULL, errcode)
    call hdf5_check_error('h5pset_fapl_mpio_f', errcode, error)
    if (present(error)) then; if (error/=0) return; endif

    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, errcode, access_prp=plist_id)
    call hdf5_check_error('h5fcreate_f', errcode, error)
    if (present(error)) then; if (error/=0) return; endif

    call h5pclose_f(plist_id, errcode)
    call hdf5_check_error('h5pclose_f', errcode, error)
    if (present(error)) then; if (error/=0) return; endif
  else
#endif
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, errcode)
    call hdf5_check_error('h5fcreate_f', errcode, error)
    if (present(error)) then; if (error/=0) return; endif
#ifdef MPI
  endif
#endif
  POP_SUB(hdf5_create_file)

end subroutine hdf5_create_file


!> Open an HDF5 file, optionally capturing the error message.
!! Valid modes are: `r` to read, `w` to `write`, and `rw` or `r+` to read and write.
!! Note that passing `w` is the same as calling `hdf5_create_file`.
!! Pass paralell_io=.true. if you want to open the file with parallel IO support,
!! and optionally set the communicator `comm` (defaults to MPI_COMM_WORLD)
subroutine hdf5_open_file(filename, mode, file_id, parallel_io, comm, error)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: mode !< r, w, rw = r+
  integer(hid_t), intent(out) :: file_id
  logical, optional, intent(in) :: parallel_io
  integer, optional, intent(in) :: comm
  integer, optional, intent(out) :: error

  integer :: errcode, access_flag, comm_
#ifdef MPI
  integer(HID_T) :: plist_id
#endif
  logical :: parallel_io_, exists

  PUSH_SUB(hdf5_open_file)

  select case(mode)
    case ('w', 'W')
      call hdf5_create_file(filename, file_id, parallel_io, comm, error)
      POP_SUB(hdf5_open_file)
      return
   case ('r', 'R')
      access_flag = H5F_ACC_RDONLY_F
   case ('r+', 'R+', 'rw', 'RW')
      access_flag = H5F_ACC_RDWR_F
   case default
      call die('Unknown HDF5 open mode: '//mode, only_root_writes=.true.)
  end select
  parallel_io_ = .false.
  if (present(parallel_io)) parallel_io_ = parallel_io
#ifdef MPI
  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm
#else
  comm_ = 0
#endif

  inquire(file=filename, exist=exists)
  if (.not.exists) then
    if (present(error)) then
      error = -1
      POP_SUB(hdf5_open_file)
    endif
    call die('Input HDF5 file "'//filename//'" does not exist.', only_root_writes=.true.)
  endif

#ifdef MPI
  if (parallel_io_) then
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, errcode)
    call hdf5_check_error('h5pcreate_f', errcode, error)
    if (present(error)) then; if (error/=0) return; endif

    call h5pset_fapl_mpio_f(plist_id, comm_, MPI_INFO_NULL, errcode)
    call hdf5_check_error('h5pset_fapl_mpio_f', errcode, error)
    if (present(error)) then; if (error/=0) return; endif

    call h5fopen_f(filename, access_flag, file_id, errcode, access_prp=plist_id)
    call hdf5_check_error('h5fopen_f', errcode, error)
    if (present(error)) then; if (error/=0) return; endif

    call h5pclose_f(plist_id, errcode)
    call hdf5_check_error('h5pclose_f', errcode, error)
    if (present(error)) then; if (error/=0) return; endif
  else
#endif
    call h5fopen_f(filename, access_flag, file_id, errcode)
    call hdf5_check_error('h5fopen_f', errcode, error)
    if (present(error)) then; if (error/=0) return; endif
#ifdef MPI
  endif
#endif
  POP_SUB(hdf5_open_file)

end subroutine hdf5_open_file


!> Close an file, optionally capturing the error message.
subroutine hdf5_close_file(file_id, error)
  integer(hid_t), intent(inout) :: file_id
  integer, optional, intent(out) :: error

  integer :: errcode

  call h5fclose_f(file_id, errcode)
  call hdf5_check_error('h5fclose_f', errcode, error)
  file_id = -1

end subroutine hdf5_close_file


!> Create an empty dataset, optionally capturing the error message.
subroutine hdf5_create_dset(loc_id, dset_name, dtype, dims, error)
  integer(HID_T), intent(in) :: loc_id !< hdf5 file id
  character(LEN=*), intent(in) :: dset_name !< hdf5 dataset name
  integer(HID_T), intent(in) :: dtype !< hdf5 data type
  integer, intent(in) :: dims(:) !< dimensions of the array
  integer, intent(out), optional :: error !< error code

  integer(HSIZE_T) :: hdims(size(dims))
  integer(HID_T) :: dset_id
  integer(HID_T) :: dspace

  PUSH_SUB(hdf5_create_dset)

  hdims(:) = dims(:)
  call safe_h5screate_simple(size(dims), hdims, dspace, error)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5dcreate(loc_id, dset_name, dtype, dspace, dset_id, error)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5dclose(dset_id, error)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5sclose(dspace, error)
  if (present(error)) then; if (error/=0) return; endif

  POP_SUB(hdf5_create_dset)

end subroutine hdf5_create_dset


!> Create an empty group, optionally capturing the error message.
subroutine hdf5_create_group(loc_id, group_name, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: group_name !< HDF5 group name
  integer, optional, intent(out) :: error

  integer(HID_T) :: group_id
  integer :: errcode

  PUSH_SUB(hdf5_create_group)

  call h5gcreate_f(loc_id, group_name, group_id, errcode)
  call hdf5_check_error('h5gcreate_f', errcode, error)
  if (present(error)) then; if (error/=0) return; endif
  call h5gclose_f(group_id, errcode)
  call hdf5_check_error('h5gclose_f', errcode, error)
  if (present(error)) then; if (error/=0) return; endif

  POP_SUB(hdf5_create_group)

end subroutine hdf5_create_group


!> Check whether the path to `link_name` exists.
logical function hdf5_exists(loc_id, link_name)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: link_name !< HDF5 group/dataset name

  integer :: error

  call h5lexists_f(loc_id, link_name, hdf5_exists, error)
  if (error/=0) hdf5_exists = .false.

end function hdf5_exists
#endif

end module hdf5_io_m
