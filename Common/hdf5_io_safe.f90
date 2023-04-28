!>=========================================================================
!!
!! Module:
!!
!! hdf5_io_safe_m        Originally by FHJ     Last Modified 08/2019 (FHJ)
!!
!!     Wrapper around low-level HDF5 routines with safe error checking.
!!
!!=========================================================================

#include "f_defs.h"

module hdf5_io_safe_m
  use, intrinsic :: iso_c_binding
  use global_m
#ifdef HDF5
  use hdf5

  implicit none

  ! FHJ: Note that all subroutines are public in this file, but we only
  ! explicitly export a few of them in mdoule hdf5_io_m

contains


!===============================================================================
! Low-level HDF5 routines that capture the error and die, if necessary.
! We don`t export these subrotuines and users should not call these directly.
!===============================================================================

subroutine hdf5_check_error(label, errcode, error_out)
  character(len=*), intent(in) :: label
  integer, intent(in) :: errcode
  integer, intent(out), optional :: error_out

  PUSH_SUB(hdf5_check_error)

  if (errcode/=0) then
    if (present(error_out)) then
      ! No need to die if the error is being captured.
      if (peinf%verb_debug) then
        write(0,'(3a,i0,".")') 'WARNING: Call to HDF5`s function ', label, ' returned ', errcode
      endif
    else
      ! Error code is not being captured, so die.
      write(0,'(3a,i0,".")') 'ERROR: Call to HDF5`s function ', label, ' returned ', errcode
      call die('Got nonzero error code from HDF5`s function '//trim(label), &
        only_root_writes=.true.)
    endif
  endif
  if (present(error_out)) error_out = errcode

  POP_SUB(hdf5_check_error)

end subroutine hdf5_check_error


subroutine safe_h5screate_simple(rank, dims, space_id, hdferr, maxdims)
  integer, intent(in) :: rank             ! number of dataspace dimensions
  integer(hsize_t), intent(in) :: dims(*) ! array with current dimension sizes
  integer(hid_t), intent(out) :: space_id ! dataspace identifier
  integer, optional, intent(out) :: hdferr! error code
                                          ! 0 on success and -1 on failure
  integer(hsize_t), optional, intent(in) :: maxdims(*)
                                          ! array with the maximum
                                          ! dimension sizes

  integer :: errcode


  PUSH_SUB(safe_h5screate_simple)
  call h5screate_simple_f(rank, dims, space_id, errcode, maxdims)
  call hdf5_check_error('h5screate_simple_f', errcode, hdferr)
  POP_SUB(safe_h5screate_simple)

end subroutine safe_h5screate_simple


subroutine safe_h5dcreate(loc_id, name, type_id, space_id, dset_id, &
     hdferr, dcpl_id, lcpl_id, dapl_id)
  integer(hid_t), intent(in) :: loc_id   ! file or group identifier
  character(len=*), intent(in) :: name   ! name of the dataset
  integer(hid_t), intent(in) :: type_id  ! datatype identifier
  integer(hid_t), intent(in) :: space_id ! dataspace identifier
  integer(hid_t), intent(out) :: dset_id ! dataset identifier
  integer, optional, intent(out) :: hdferr! error code
                                         ! 0 on success and -1 on failure
  integer(hid_t), optional, intent(in) :: dcpl_id
                                         ! dataset creation property list
  integer(hid_t), optional, intent(in) :: lcpl_id
                                         ! link creation property list
  integer(hid_t), optional, intent(in) :: dapl_id
                                         ! dataset access property list

  integer :: errcode

  PUSH_SUB(safe_h5dcreate)
  call h5dcreate_f(loc_id, name, type_id, space_id, dset_id, &
     errcode, dcpl_id, lcpl_id, dapl_id)
  call hdf5_check_error('h5dcreate_f', errcode, hdferr)
  POP_SUB(safe_h5dcreate)

end subroutine safe_h5dcreate


subroutine safe_h5dopen(loc_id, name, dset_id, hdferr, dapl_id)
  integer(hid_t), intent(in) :: loc_id   ! file or group identifier
  character(len=*), intent(in) :: name   ! name of the dataset
  integer(hid_t), intent(out) :: dset_id ! dataset identifier
  integer, optional, intent(out) :: hdferr         ! error code:
                                         ! 0 on success and -1 on failure
  integer(hid_t), optional, intent(in) :: dapl_id
                                         ! dataset access property list

  integer :: errcode

  PUSH_SUB(safe_h5dopen)
  call h5dopen_f(loc_id, name, dset_id, errcode, dapl_id)
  call hdf5_check_error('h5dcreate_f', errcode, hdferr)
  POP_SUB(safe_h5dopen)

end subroutine safe_h5dopen


subroutine safe_h5dclose(dset_id, hdferr)
  integer(hid_t), intent(in) :: dset_id
  integer, optional, intent(out) :: hdferr

  integer :: errcode

  PUSH_SUB(safe_h5dclose)
  call h5dclose_f(dset_id, errcode)
  call hdf5_check_error('h5dclose_f', errcode, hdferr)
  POP_SUB(safe_h5dclose)

end subroutine safe_h5dclose


subroutine safe_h5dget_space(dataset_id, dataspace_id, hdferr)
  integer(hid_t), intent(in) :: dataset_id      ! dataset identifier
  integer(hid_t), intent(out) :: dataspace_id   ! dataspace identifier
  integer, optional, intent(out) :: hdferr       ! error code
                                                ! 0 on success and -1 on failure

  integer :: errcode

  PUSH_SUB(safe_h5dget_space)
  call h5dget_space_f(dataset_id, dataspace_id, errcode)
  call hdf5_check_error('h5dget_space_f', errcode, hdferr)
  POP_SUB(safe_h5dget_space)

end subroutine safe_h5dget_space


subroutine safe_h5sselect_hyperslab(space_id, operator, start, count, &
  hdferr, stride, block)
  integer(hid_t), intent(in) :: space_id  ! dataspace identifier
  integer, intent(in) :: operator         ! flag, valid values are:
                                          !    h5s_select_set_f
                                          !    h5s_select_or_f
  integer(hsize_t), dimension(*), intent(in) :: start
                                          ! offset of start of hyperslab
  integer(hsize_t), dimension(*), intent(in) :: count
                                          ! number of blocks to select
                                          ! from dataspace
  integer, optional, intent(out) :: hdferr! error code
                                          ! 0 on success and -1 on failure
  integer(hsize_t), dimension(:), optional, intent(in) :: stride
                                          ! array of how many elements to
                                          ! move in each direction
  integer(hsize_t), dimension(:), optional, intent(in) :: block
                                          ! size of the element block

  integer :: errcode

  PUSH_SUB(safe_h5sselect_hyperslab)
  call h5sselect_hyperslab_f(space_id, operator, start, count, errcode, &
                             stride, block)
  call hdf5_check_error('h5sselect_hyperslab_f', errcode, hdferr)
  POP_SUB(safe_h5sselect_hyperslab)

end subroutine safe_h5sselect_hyperslab


subroutine safe_h5sclose(space_id, hdferr)
  integer(hid_t), intent(in) :: space_id
  integer, optional, intent(out) :: hdferr

  integer :: errcode

  PUSH_SUB(safe_h5sclose)
  call h5sclose_f(space_id, errcode)
  call hdf5_check_error('h5sclose_f', errcode, hdferr)
  POP_SUB(safe_h5sclose)

end subroutine safe_h5sclose


subroutine safe_h5pcreate(classtype, prp_id, hdferr)
  integer(hid_t), intent(in) :: classtype
  integer(hid_t), intent(out) :: prp_id
  integer, optional, intent(out) :: hdferr

  integer :: errcode

  PUSH_SUB(safe_h5pcreate)
  call h5pcreate_f(classtype, prp_id, errcode)
  call hdf5_check_error('h5pcreate_f', errcode, hdferr)
  POP_SUB(safe_h5pcreate)

end subroutine safe_h5pcreate


subroutine safe_h5pset_dxpl_mpio(prp_id, data_xfer_mode, hdferr)
  integer(hid_t), intent(in) :: prp_id
  integer, intent(in) :: data_xfer_mode
  integer, optional, intent(out) :: hdferr

  integer :: errcode

  PUSH_SUB(safe_h5pset_dxpl_mpio)
#ifdef MPI
  call h5pset_dxpl_mpio_f(prp_id, data_xfer_mode, errcode)
  call hdf5_check_error('h5pset_dxpl_mpio_f', errcode, hdferr)
#else
  if (present(hdferr)) hdferr = 0
#endif
  POP_SUB(safe_h5pset_dxpl_mpio)

end subroutine safe_h5pset_dxpl_mpio


subroutine safe_h5pclose(prp_id, hdferr)
  integer(hid_t), intent(in) :: prp_id
  integer, optional, intent(out) :: hdferr

  integer :: errcode

  PUSH_SUB(safe_h5pclose)
  call h5pclose_f(prp_id, errcode)
  call hdf5_check_error('h5pclose_f', errcode, hdferr)
  POP_SUB(safe_h5pclose)

end subroutine safe_h5pclose

#endif
end module hdf5_io_safe_m
