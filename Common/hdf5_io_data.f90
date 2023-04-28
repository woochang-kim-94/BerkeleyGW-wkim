! File hdf5_io_data.f90 automatically created from hdf5_io_data.f90p by mako_preprocess.py.
! Do not edit the resulting file (hdf5_io_data.f90) directly!

!>=========================================================================
!!
!! Module:
!!
!! hdf5_io_data_m        Originally by FHJ     Last Modified 08/2019 (FHJ)
!!
!!     High-level routines to read and write data in HDF5 format.
!!     This module uses python Mako template to generate Fortran code.
!!     Make sure to run the `mako_preprocess.py` if you edit the source
!!     .f90p file, and don`t edit the generated .f90 directly.
!!
!!     In the call signatures below:
!!     - `loc_id` is the HDF5 file handle,
!!     - `dset_name` is the path to the dataset
!!     - `error` is an optional argument with the error code, which should
!!       typically *not* be passed unless you want to do the error checking
!!       yourself and don`t want the code to die in case the HDF5 library
!!       fails.
!!
!!     Routines:
!!     - hdf5_{read,write}_{int,double,complex,logical,scalar}
!!       (loc_id, dset_name, buf, error):
!!       integer(HID_T), intent(in) :: loc_id
!!       character(len=*), intent(in) :: dset_name
!!       {data_type}, intent({intent}), dimension(*) :: buf
!!       integer, intent(out), optional :: error
!!
!!       Reads/writes a single element store in buffer `buf`.
!!       When writing, the dataset in `dset_name` may or may not exist, but if
!!       it already exist, it must have the correct type and shape.
!!
!!
!!     - hdf5_{read,write}_{int,double,complex,logical,scalar}_array
!!       (loc_id, dset_name, dims, buf, error):
!!       integer(HID_T), intent(in) :: loc_id
!!       character(len=*), intent(in) :: dset_name
!!       integer, intent(in) :: dims(:)
!!       {data_type}, intent({intent}), dimension(*) :: buf
!!       integer, intent(out), optional :: error
!!
!!       Reads/writes an array of size `dims` from/to buffer `buf`.
!!       This routine writes the full array, and not part thereof.
!!       When writing, the dataset in `dset_name` may or may not exist, but if
!!       it already exist, it must have the correct type and shape.
!!
!!
!!     - hdf5_{read,write}_{int,double,complex,logical,scalar}_hyperslab
!!       (loc_id, dset_name, rw_count, offset, buf, error)
!!       integer(HID_T), intent(in) :: loc_id
!!       character(len=*), intent(in) :: dset_name
!!       integer, intent(in) :: rw_count(:)
!!       integer, intent(in) :: offset(:)
!!       {data_type}, intent({intent}), dimension(*) :: buf
!!       integer, intent(out), optional :: error
!!
!!       Reads/writes part of an array stored in `buf`, where;
!!       * `rw_count` is the number of elements read/written from/to
!!         the daset for each dimension. This does not need to be the same
!!         as the size of the dataset, but it is the overall size of the
!!         output/input array `buf`
!!       * `offset` is the offset for each dimension, relative to the HDF5
!!         dataset when reading/writing from/to the HDF5 file.
!!
!!       Unlike the previous subroutines, the dataset must already exist
!!       before calling it, since the size of the dataset in file may be
!!       different than `rw_count`.
!!
!!     Notes:
!!     - Logical datasets are stored as integers.
!!       * When writing, we map .false.=>0 and .true.=>1
!!       * When reading, we map 0=>.false. and nonzero=>.true.
!!
!!     - Because there is no native HDF5 support for complex data, we
!!       read/write complex dataset as double dataset by appending one
!!       extra fast dimension with size=2.
!!
!!     - There are interfaces for scalar dataset, which either read/write
!!       complex or double data depending on the type of `buf`. However,
!!       to make the number of dimensions the same for double or complex
!!       scalar data, the scalar routines will **also append an extra
!!       dimension with size=1** when reading/writing double data.
!!       The code will not add an extra dimension if you call the
!!       non-scalar routine to write a double buffer.
!!
!!       For example, suppose you want to write `buf(1:N)`:
!!       * If `buf` is `real(DP)`:
!!         - `write_double_array` writes a real dataset of dimension [N].
!!         - `write_scalar_array` writes a real dataset of dimension [1,N].
!!       * If `buf` is `complex(DPC)`:
!!         - Both `write_complex_array` and `write_scalar_array` write
!!           real datasets of dimension [2,N].
!!
!!=========================================================================

#include "f_defs.h"

module hdf5_io_data_m
  use, intrinsic :: iso_c_binding
  use global_m
  use hdf5_io_safe_m
#ifdef HDF5
  use hdf5

  implicit none

  private

  public :: &
    hdf5_read_int, &
    hdf5_read_int_array, &
    hdf5_read_int_hyperslab, &
    hdf5_write_int, &
    hdf5_write_int_array, &
    hdf5_write_int_hyperslab, &
    hdf5_read_double, &
    hdf5_read_double_array, &
    hdf5_read_double_hyperslab, &
    hdf5_write_double, &
    hdf5_write_double_array, &
    hdf5_write_double_hyperslab, &
    hdf5_read_complex, &
    hdf5_read_complex_array, &
    hdf5_read_complex_hyperslab, &
    hdf5_write_complex, &
    hdf5_write_complex_array, &
    hdf5_write_complex_hyperslab, &
    hdf5_read_scalar, &
    hdf5_read_scalar_array, &
    hdf5_read_scalar_hyperslab, &
    hdf5_write_scalar, &
    hdf5_write_scalar_array, &
    hdf5_write_scalar_hyperslab, &
    hdf5_read_logical, &
    hdf5_read_logical_array, &
    hdf5_read_logical_hyperslab, &
    hdf5_write_logical, &
    hdf5_write_logical_array, &
    hdf5_write_logical_hyperslab


  interface hdf5_read_scalar
    module procedure hdf5_read_scalar_double, hdf5_read_scalar_complex
  end interface

  interface hdf5_read_scalar_array
    module procedure &
hdf5_read_scalar_double_array_1, hdf5_read_scalar_complex_array_1, &
hdf5_read_scalar_double_array_2, hdf5_read_scalar_complex_array_2, &
hdf5_read_scalar_double_array_3, hdf5_read_scalar_complex_array_3, &
hdf5_read_scalar_double_array_4, hdf5_read_scalar_complex_array_4, &
hdf5_read_scalar_double_array_5, hdf5_read_scalar_complex_array_5, &
hdf5_read_scalar_double_array_6, hdf5_read_scalar_complex_array_6, &
hdf5_read_scalar_double_array_7, hdf5_read_scalar_complex_array_7
  end interface

  interface hdf5_read_scalar_hyperslab
    module procedure &
hdf5_read_scalar_double_hyperslab_1, hdf5_read_scalar_complex_hyperslab_1, &
hdf5_read_scalar_double_hyperslab_2, hdf5_read_scalar_complex_hyperslab_2, &
hdf5_read_scalar_double_hyperslab_3, hdf5_read_scalar_complex_hyperslab_3, &
hdf5_read_scalar_double_hyperslab_4, hdf5_read_scalar_complex_hyperslab_4, &
hdf5_read_scalar_double_hyperslab_5, hdf5_read_scalar_complex_hyperslab_5, &
hdf5_read_scalar_double_hyperslab_6, hdf5_read_scalar_complex_hyperslab_6, &
hdf5_read_scalar_double_hyperslab_7, hdf5_read_scalar_complex_hyperslab_7
  end interface

  interface hdf5_write_scalar
    module procedure hdf5_write_scalar_double, hdf5_write_scalar_complex
  end interface

  interface hdf5_write_scalar_array
    module procedure &
hdf5_write_scalar_double_array_1, hdf5_write_scalar_complex_array_1, &
hdf5_write_scalar_double_array_2, hdf5_write_scalar_complex_array_2, &
hdf5_write_scalar_double_array_3, hdf5_write_scalar_complex_array_3, &
hdf5_write_scalar_double_array_4, hdf5_write_scalar_complex_array_4, &
hdf5_write_scalar_double_array_5, hdf5_write_scalar_complex_array_5, &
hdf5_write_scalar_double_array_6, hdf5_write_scalar_complex_array_6, &
hdf5_write_scalar_double_array_7, hdf5_write_scalar_complex_array_7
  end interface

  interface hdf5_write_scalar_hyperslab
    module procedure &
hdf5_write_scalar_double_hyperslab_1, hdf5_write_scalar_complex_hyperslab_1, &
hdf5_write_scalar_double_hyperslab_2, hdf5_write_scalar_complex_hyperslab_2, &
hdf5_write_scalar_double_hyperslab_3, hdf5_write_scalar_complex_hyperslab_3, &
hdf5_write_scalar_double_hyperslab_4, hdf5_write_scalar_complex_hyperslab_4, &
hdf5_write_scalar_double_hyperslab_5, hdf5_write_scalar_complex_hyperslab_5, &
hdf5_write_scalar_double_hyperslab_6, hdf5_write_scalar_complex_hyperslab_6, &
hdf5_write_scalar_double_hyperslab_7, hdf5_write_scalar_complex_hyperslab_7
  end interface

contains


!################################################
!# High-level interfaces for general data types #
!################################################


!-------------------------------
! Routines for int data type
!-------------------------------

!> read int value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_read_int(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(inout), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  integer, pointer :: buf_1d(:)

  PUSH_SUB(hdf5_read_int)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_read_int_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_read_int)

end subroutine hdf5_read_int

!> read int array to/from an HDF5 file.
subroutine hdf5_read_int_array(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  integer, intent(inout), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_int_array)
  call hdf5_read_int_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_int_array)

end subroutine hdf5_read_int_array

!> read section (hyperslab) of int array to/from an HDF5 file.
subroutine hdf5_read_int_hyperslab(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  integer, intent(inout), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_int_hyperslab)
  call hdf5_read_int_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_int_hyperslab)

end subroutine hdf5_read_int_hyperslab

!> write int value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_write_int(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  integer, pointer :: buf_1d(:)

  PUSH_SUB(hdf5_write_int)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_write_int_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_write_int)

end subroutine hdf5_write_int

!> write int array to/from an HDF5 file.
subroutine hdf5_write_int_array(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  integer, intent(in), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_int_array)
  call hdf5_write_int_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_int_array)

end subroutine hdf5_write_int_array

!> write section (hyperslab) of int array to/from an HDF5 file.
subroutine hdf5_write_int_hyperslab(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  integer, intent(in), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_int_hyperslab)
  call hdf5_write_int_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_int_hyperslab)

end subroutine hdf5_write_int_hyperslab

!-------------------------------
! Routines for double data type
!-------------------------------

!> read double value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_read_double(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  real(DP), intent(inout), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  real(DP), pointer :: buf_1d(:)

  PUSH_SUB(hdf5_read_double)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_read_double_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_read_double)

end subroutine hdf5_read_double

!> read double array to/from an HDF5 file.
subroutine hdf5_read_double_array(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(inout), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_double_array)
  call hdf5_read_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_double_array)

end subroutine hdf5_read_double_array

!> read section (hyperslab) of double array to/from an HDF5 file.
subroutine hdf5_read_double_hyperslab(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(inout), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_double_hyperslab)
  call hdf5_read_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_double_hyperslab)

end subroutine hdf5_read_double_hyperslab

!> write double value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_write_double(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  real(DP), intent(in), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  real(DP), pointer :: buf_1d(:)

  PUSH_SUB(hdf5_write_double)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_write_double_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_write_double)

end subroutine hdf5_write_double

!> write double array to/from an HDF5 file.
subroutine hdf5_write_double_array(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(in), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_double_array)
  call hdf5_write_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_double_array)

end subroutine hdf5_write_double_array

!> write section (hyperslab) of double array to/from an HDF5 file.
subroutine hdf5_write_double_hyperslab(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(in), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_double_hyperslab)
  call hdf5_write_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_double_hyperslab)

end subroutine hdf5_write_double_hyperslab

!-------------------------------
! Routines for complex data type
!-------------------------------

!> read complex value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_read_complex(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  complex(DPC), intent(inout), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  complex(DPC), pointer :: buf_1d(:)

  PUSH_SUB(hdf5_read_complex)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_read_complex_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_read_complex)

end subroutine hdf5_read_complex

!> read complex array to/from an HDF5 file.
subroutine hdf5_read_complex_array(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(inout), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_complex_array)
  call hdf5_read_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_complex_array)

end subroutine hdf5_read_complex_array

!> read section (hyperslab) of complex array to/from an HDF5 file.
subroutine hdf5_read_complex_hyperslab(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(inout), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_complex_hyperslab)
  call hdf5_read_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_complex_hyperslab)

end subroutine hdf5_read_complex_hyperslab

!> write complex value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_write_complex(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  complex(DPC), intent(in), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  complex(DPC), pointer :: buf_1d(:)

  PUSH_SUB(hdf5_write_complex)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_write_complex_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_write_complex)

end subroutine hdf5_write_complex

!> write complex array to/from an HDF5 file.
subroutine hdf5_write_complex_array(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(in), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_complex_array)
  call hdf5_write_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_complex_array)

end subroutine hdf5_write_complex_array

!> write section (hyperslab) of complex array to/from an HDF5 file.
subroutine hdf5_write_complex_hyperslab(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(in), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_complex_hyperslab)
  call hdf5_write_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_complex_hyperslab)

end subroutine hdf5_write_complex_hyperslab

!-------------------------------
! Routines for logical data type
!-------------------------------

!> read logical value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_read_logical(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  logical, intent(inout), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  logical, pointer :: buf_1d(:)

  PUSH_SUB(hdf5_read_logical)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_read_logical_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_read_logical)

end subroutine hdf5_read_logical

!> read logical array to/from an HDF5 file.
subroutine hdf5_read_logical_array(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  logical, intent(inout), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_logical_array)
  call hdf5_read_logical_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_logical_array)

end subroutine hdf5_read_logical_array

!> read section (hyperslab) of logical array to/from an HDF5 file.
subroutine hdf5_read_logical_hyperslab(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  logical, intent(inout), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_logical_hyperslab)
  call hdf5_read_logical_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_logical_hyperslab)

end subroutine hdf5_read_logical_hyperslab

!> write logical value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_write_logical(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  logical, intent(in), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  logical, pointer :: buf_1d(:)

  PUSH_SUB(hdf5_write_logical)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_write_logical_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_write_logical)

end subroutine hdf5_write_logical

!> write logical array to/from an HDF5 file.
subroutine hdf5_write_logical_array(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  logical, intent(in), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_logical_array)
  call hdf5_write_logical_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_logical_array)

end subroutine hdf5_write_logical_array

!> write section (hyperslab) of logical array to/from an HDF5 file.
subroutine hdf5_write_logical_hyperslab(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  logical, intent(in), dimension(*) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_logical_hyperslab)
  call hdf5_write_logical_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_logical_hyperslab)

end subroutine hdf5_write_logical_hyperslab

!#######################################################
!# End of high-level interfaces for general data types #
!#######################################################


!###############################################
!# High-level interfaces for scalar data types #
!###############################################


!-------------------------------
! Routines for ('scalar', 'double') data type
!-------------------------------

!> read ('scalar', 'double') value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_read_scalar_double(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  real(DP), intent(inout), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  real(DP), pointer :: buf_1d(:)

  PUSH_SUB(hdf5_read_scalar_double)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_read_scalar_double)

end subroutine hdf5_read_scalar_double

!> read ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_array_1(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(inout), dimension(:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_array_1)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_double_array_1)

end subroutine hdf5_read_scalar_double_array_1

!> read section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_hyperslab_1(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(inout), dimension(:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_hyperslab_1)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_double_hyperslab_1)

end subroutine hdf5_read_scalar_double_hyperslab_1
!> read ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_array_2(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(inout), dimension(:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_array_2)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_double_array_2)

end subroutine hdf5_read_scalar_double_array_2

!> read section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_hyperslab_2(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(inout), dimension(:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_hyperslab_2)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_double_hyperslab_2)

end subroutine hdf5_read_scalar_double_hyperslab_2
!> read ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_array_3(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(inout), dimension(:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_array_3)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_double_array_3)

end subroutine hdf5_read_scalar_double_array_3

!> read section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_hyperslab_3(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(inout), dimension(:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_hyperslab_3)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_double_hyperslab_3)

end subroutine hdf5_read_scalar_double_hyperslab_3
!> read ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_array_4(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(inout), dimension(:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_array_4)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_double_array_4)

end subroutine hdf5_read_scalar_double_array_4

!> read section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_hyperslab_4(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(inout), dimension(:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_hyperslab_4)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_double_hyperslab_4)

end subroutine hdf5_read_scalar_double_hyperslab_4
!> read ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_array_5(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(inout), dimension(:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_array_5)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_double_array_5)

end subroutine hdf5_read_scalar_double_array_5

!> read section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_hyperslab_5(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(inout), dimension(:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_hyperslab_5)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_double_hyperslab_5)

end subroutine hdf5_read_scalar_double_hyperslab_5
!> read ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_array_6(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(inout), dimension(:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_array_6)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_double_array_6)

end subroutine hdf5_read_scalar_double_array_6

!> read section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_hyperslab_6(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(inout), dimension(:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_hyperslab_6)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_double_hyperslab_6)

end subroutine hdf5_read_scalar_double_hyperslab_6
!> read ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_array_7(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(inout), dimension(:,:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_array_7)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_double_array_7)

end subroutine hdf5_read_scalar_double_array_7

!> read section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_read_scalar_double_hyperslab_7(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(inout), dimension(:,:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_double_hyperslab_7)
  call hdf5_read_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_double_hyperslab_7)

end subroutine hdf5_read_scalar_double_hyperslab_7

!> write ('scalar', 'double') value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_write_scalar_double(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  real(DP), intent(in), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  real(DP), pointer :: buf_1d(:)

  PUSH_SUB(hdf5_write_scalar_double)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_write_scalar_double)

end subroutine hdf5_write_scalar_double

!> write ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_array_1(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(in), dimension(:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_array_1)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_double_array_1)

end subroutine hdf5_write_scalar_double_array_1

!> write section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_hyperslab_1(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(in), dimension(:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_hyperslab_1)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_double_hyperslab_1)

end subroutine hdf5_write_scalar_double_hyperslab_1
!> write ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_array_2(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(in), dimension(:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_array_2)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_double_array_2)

end subroutine hdf5_write_scalar_double_array_2

!> write section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_hyperslab_2(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(in), dimension(:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_hyperslab_2)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_double_hyperslab_2)

end subroutine hdf5_write_scalar_double_hyperslab_2
!> write ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_array_3(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(in), dimension(:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_array_3)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_double_array_3)

end subroutine hdf5_write_scalar_double_array_3

!> write section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_hyperslab_3(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(in), dimension(:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_hyperslab_3)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_double_hyperslab_3)

end subroutine hdf5_write_scalar_double_hyperslab_3
!> write ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_array_4(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(in), dimension(:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_array_4)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_double_array_4)

end subroutine hdf5_write_scalar_double_array_4

!> write section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_hyperslab_4(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(in), dimension(:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_hyperslab_4)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_double_hyperslab_4)

end subroutine hdf5_write_scalar_double_hyperslab_4
!> write ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_array_5(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(in), dimension(:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_array_5)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_double_array_5)

end subroutine hdf5_write_scalar_double_array_5

!> write section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_hyperslab_5(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(in), dimension(:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_hyperslab_5)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_double_hyperslab_5)

end subroutine hdf5_write_scalar_double_hyperslab_5
!> write ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_array_6(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(in), dimension(:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_array_6)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_double_array_6)

end subroutine hdf5_write_scalar_double_array_6

!> write section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_hyperslab_6(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(in), dimension(:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_hyperslab_6)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_double_hyperslab_6)

end subroutine hdf5_write_scalar_double_hyperslab_6
!> write ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_array_7(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  real(DP), intent(in), dimension(:,:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_array_7)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_double_array_7)

end subroutine hdf5_write_scalar_double_array_7

!> write section (hyperslab) of ('scalar', 'double') array to/from an HDF5 file.
subroutine hdf5_write_scalar_double_hyperslab_7(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  real(DP), intent(in), dimension(:,:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_double_hyperslab_7)
  call hdf5_write_scalar_double_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_double_hyperslab_7)

end subroutine hdf5_write_scalar_double_hyperslab_7

!-------------------------------
! Routines for ('scalar', 'complex') data type
!-------------------------------

!> read ('scalar', 'complex') value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_read_scalar_complex(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  complex(DPC), intent(inout), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  complex(DPC), pointer :: buf_1d(:)

  PUSH_SUB(hdf5_read_scalar_complex)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_read_scalar_complex)

end subroutine hdf5_read_scalar_complex

!> read ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_array_1(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(inout), dimension(:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_array_1)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_complex_array_1)

end subroutine hdf5_read_scalar_complex_array_1

!> read section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_hyperslab_1(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(inout), dimension(:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_hyperslab_1)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_complex_hyperslab_1)

end subroutine hdf5_read_scalar_complex_hyperslab_1
!> read ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_array_2(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(inout), dimension(:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_array_2)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_complex_array_2)

end subroutine hdf5_read_scalar_complex_array_2

!> read section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_hyperslab_2(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(inout), dimension(:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_hyperslab_2)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_complex_hyperslab_2)

end subroutine hdf5_read_scalar_complex_hyperslab_2
!> read ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_array_3(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(inout), dimension(:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_array_3)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_complex_array_3)

end subroutine hdf5_read_scalar_complex_array_3

!> read section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_hyperslab_3(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(inout), dimension(:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_hyperslab_3)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_complex_hyperslab_3)

end subroutine hdf5_read_scalar_complex_hyperslab_3
!> read ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_array_4(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(inout), dimension(:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_array_4)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_complex_array_4)

end subroutine hdf5_read_scalar_complex_array_4

!> read section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_hyperslab_4(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(inout), dimension(:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_hyperslab_4)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_complex_hyperslab_4)

end subroutine hdf5_read_scalar_complex_hyperslab_4
!> read ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_array_5(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(inout), dimension(:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_array_5)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_complex_array_5)

end subroutine hdf5_read_scalar_complex_array_5

!> read section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_hyperslab_5(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(inout), dimension(:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_hyperslab_5)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_complex_hyperslab_5)

end subroutine hdf5_read_scalar_complex_hyperslab_5
!> read ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_array_6(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(inout), dimension(:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_array_6)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_complex_array_6)

end subroutine hdf5_read_scalar_complex_array_6

!> read section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_hyperslab_6(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(inout), dimension(:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_hyperslab_6)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_complex_hyperslab_6)

end subroutine hdf5_read_scalar_complex_hyperslab_6
!> read ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_array_7(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(inout), dimension(:,:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_array_7)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_read_scalar_complex_array_7)

end subroutine hdf5_read_scalar_complex_array_7

!> read section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_read_scalar_complex_hyperslab_7(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(inout), dimension(:,:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_read_scalar_complex_hyperslab_7)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_read_scalar_complex_hyperslab_7)

end subroutine hdf5_read_scalar_complex_hyperslab_7

!> write ('scalar', 'complex') value (rank-0 array) to/from an HDF5 file.
subroutine hdf5_write_scalar_complex(loc_id, dset_name, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  complex(DPC), intent(in), target :: buf !< data buffer
  integer, intent(out), optional :: error !< HDF5 error code

  complex(DPC), pointer :: buf_1d(:)

  PUSH_SUB(hdf5_write_scalar_complex)
  call c_f_pointer(c_loc(buf), buf_1d, [1])
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, [-1_HSIZE_T], buf_1d, error=error)
  POP_SUB(hdf5_write_scalar_complex)

end subroutine hdf5_write_scalar_complex

!> write ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_array_1(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(in), dimension(:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_array_1)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_complex_array_1)

end subroutine hdf5_write_scalar_complex_array_1

!> write section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_hyperslab_1(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(in), dimension(:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_hyperslab_1)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_complex_hyperslab_1)

end subroutine hdf5_write_scalar_complex_hyperslab_1
!> write ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_array_2(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(in), dimension(:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_array_2)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_complex_array_2)

end subroutine hdf5_write_scalar_complex_array_2

!> write section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_hyperslab_2(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(in), dimension(:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_hyperslab_2)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_complex_hyperslab_2)

end subroutine hdf5_write_scalar_complex_hyperslab_2
!> write ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_array_3(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(in), dimension(:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_array_3)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_complex_array_3)

end subroutine hdf5_write_scalar_complex_array_3

!> write section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_hyperslab_3(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(in), dimension(:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_hyperslab_3)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_complex_hyperslab_3)

end subroutine hdf5_write_scalar_complex_hyperslab_3
!> write ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_array_4(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(in), dimension(:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_array_4)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_complex_array_4)

end subroutine hdf5_write_scalar_complex_array_4

!> write section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_hyperslab_4(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(in), dimension(:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_hyperslab_4)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_complex_hyperslab_4)

end subroutine hdf5_write_scalar_complex_hyperslab_4
!> write ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_array_5(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(in), dimension(:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_array_5)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_complex_array_5)

end subroutine hdf5_write_scalar_complex_array_5

!> write section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_hyperslab_5(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(in), dimension(:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_hyperslab_5)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_complex_hyperslab_5)

end subroutine hdf5_write_scalar_complex_hyperslab_5
!> write ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_array_6(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(in), dimension(:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_array_6)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_complex_array_6)

end subroutine hdf5_write_scalar_complex_array_6

!> write section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_hyperslab_6(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(in), dimension(:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_hyperslab_6)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_complex_hyperslab_6)

end subroutine hdf5_write_scalar_complex_hyperslab_6
!> write ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_array_7(loc_id, dset_name, dims, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  complex(DPC), intent(in), dimension(:,:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_array_7)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(dims,kind=HSIZE_T), buf, error=error)
  POP_SUB(hdf5_write_scalar_complex_array_7)

end subroutine hdf5_write_scalar_complex_array_7

!> write section (hyperslab) of ('scalar', 'complex') array to/from an HDF5 file.
subroutine hdf5_write_scalar_complex_hyperslab_7(loc_id, dset_name, rw_count, offset, buf, error)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(len=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from or write to the dataset for each dimension.
  !> This does not need to be the same as the size of the dataset.
  integer, intent(in) :: rw_count(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offset(:)
  complex(DPC), intent(in), dimension(:,:,:,:,:,:,:) :: buf !< data buffer
  integer, intent(out), optional :: error !< error code

  PUSH_SUB(hdf5_write_scalar_complex_hyperslab_7)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, int(rw_count,kind=HSIZE_T), buf, error, offsetf=offset)
  POP_SUB(hdf5_write_scalar_complex_hyperslab_7)

end subroutine hdf5_write_scalar_complex_hyperslab_7

!######################################################
!# End of high-level interfaces for scalar data types #
!######################################################


!###############################################
!# LOW LEVEL routines -- INT + DOUBLE versions #
!###############################################


!> read int to/from an HDF5 file.
subroutine hdf5_read_int_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  integer, intent(inout), dimension(*) :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)


  integer(HSIZE_T) :: hcountf(size(countf)) !< Count for file dataspace
  integer(HSIZE_T) :: hcountm(1) !< Count for memory dataspace
  integer(HSIZE_T) :: hoffsetf(size(countf)) !< Offset for file filespace_id
  integer(HID_T) :: dset_id, filespace_id, memspace_id
  integer :: errcode, rank
  logical :: exists

  PUSH_SUB(hdf5_read_int_lowlevel)

  if (present(error)) error = 0
  call h5lexists_f(loc_id, dset_name, exists, errcode)
  call hdf5_check_error('h5lexists_f', merge(0,1,exists.and.errcode==0), error)
  if (present(error)) then; if (error/=0) return; endif

  ! FHJ: Create or get file dataspace
  if (present(offsetf)) then
    hoffsetf(:) = offsetf(:)
  else
    hoffsetf(:) = 0
  endif
  hcountf(:) = countf(:)
  rank = size(hcountf)
  if (any(countf<1)) then
    rank = 0
    ! Note: hcountf and hoffsetf are not referenced if rank==0
  endif

  if (.not.exists) then
    ! FHJ: Create file dataspace
    call safe_h5screate_simple(rank, hcountf, filespace_id, error)
    if (present(error)) then; if (error/=0) return; endif
    ! FHJ: Create dataset
    call safe_h5dcreate(loc_id, dset_name, H5T_NATIVE_INTEGER, filespace_id, dset_id, error)
    if (present(error)) then; if (error/=0) return; endif
  else
    ! FHJ: Open dataset
    call safe_h5dopen(loc_id, dset_name, dset_id, error)
    if (present(error)) then; if (error/=0) return; endif
    ! FHJ: Get file dataspace
    call safe_h5dget_space(dset_id, filespace_id, error)
    if (present(error)) then; if (error/=0) return; endif
  endif

  ! FHJ: Select hyperslab from file
  if (rank>0) then
    call safe_h5sselect_hyperslab(filespace_id, H5S_SELECT_SET_F, hoffsetf, hcountf, error)
    if (present(error)) then; if (error/=0) return; endif
  endif

  ! FHJ: Create flat memory filespace_id
  hcountm(1) = max(1, product(countf))
  call safe_h5screate_simple(1, hcountm, memspace_id, error)
  if (present(error)) then; if (error/=0) return; endif

  ! FHJ: read filespace_id
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, buf, hcountm, errcode, &
    memspace_id, filespace_id)
  call hdf5_check_error('h5dread_f', errcode, error)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5sclose(memspace_id, errcode)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5sclose(filespace_id, errcode)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5dclose(dset_id, errcode)
  if (present(error)) then; if (error/=0) return; endif

  POP_SUB(hdf5_read_int_lowlevel)

end subroutine hdf5_read_int_lowlevel


!> read double to/from an HDF5 file.
subroutine hdf5_read_double_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  real(DP), intent(inout), dimension(*) :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)


  integer(HSIZE_T) :: hcountf(size(countf)) !< Count for file dataspace
  integer(HSIZE_T) :: hcountm(1) !< Count for memory dataspace
  integer(HSIZE_T) :: hoffsetf(size(countf)) !< Offset for file filespace_id
  integer(HID_T) :: dset_id, filespace_id, memspace_id
  integer :: errcode, rank
  logical :: exists

  PUSH_SUB(hdf5_read_double_lowlevel)

  if (present(error)) error = 0
  call h5lexists_f(loc_id, dset_name, exists, errcode)
  call hdf5_check_error('h5lexists_f', merge(0,1,exists.and.errcode==0), error)
  if (present(error)) then; if (error/=0) return; endif

  ! FHJ: Create or get file dataspace
  if (present(offsetf)) then
    hoffsetf(:) = offsetf(:)
  else
    hoffsetf(:) = 0
  endif
  hcountf(:) = countf(:)
  rank = size(hcountf)
  if (any(countf<1)) then
    rank = 0
    ! Note: hcountf and hoffsetf are not referenced if rank==0
  endif

  if (.not.exists) then
    ! FHJ: Create file dataspace
    call safe_h5screate_simple(rank, hcountf, filespace_id, error)
    if (present(error)) then; if (error/=0) return; endif
    ! FHJ: Create dataset
    call safe_h5dcreate(loc_id, dset_name, H5T_NATIVE_DOUBLE, filespace_id, dset_id, error)
    if (present(error)) then; if (error/=0) return; endif
  else
    ! FHJ: Open dataset
    call safe_h5dopen(loc_id, dset_name, dset_id, error)
    if (present(error)) then; if (error/=0) return; endif
    ! FHJ: Get file dataspace
    call safe_h5dget_space(dset_id, filespace_id, error)
    if (present(error)) then; if (error/=0) return; endif
  endif

  ! FHJ: Select hyperslab from file
  if (rank>0) then
    call safe_h5sselect_hyperslab(filespace_id, H5S_SELECT_SET_F, hoffsetf, hcountf, error)
    if (present(error)) then; if (error/=0) return; endif
  endif

  ! FHJ: Create flat memory filespace_id
  hcountm(1) = max(1, product(countf))
  call safe_h5screate_simple(1, hcountm, memspace_id, error)
  if (present(error)) then; if (error/=0) return; endif

  ! FHJ: read filespace_id
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, hcountm, errcode, &
    memspace_id, filespace_id)
  call hdf5_check_error('h5dread_f', errcode, error)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5sclose(memspace_id, errcode)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5sclose(filespace_id, errcode)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5dclose(dset_id, errcode)
  if (present(error)) then; if (error/=0) return; endif

  POP_SUB(hdf5_read_double_lowlevel)

end subroutine hdf5_read_double_lowlevel


!> write int to/from an HDF5 file.
subroutine hdf5_write_int_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  integer, intent(in), dimension(*) :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)


  integer(HSIZE_T) :: hcountf(size(countf)) !< Count for file dataspace
  integer(HSIZE_T) :: hcountm(1) !< Count for memory dataspace
  integer(HSIZE_T) :: hoffsetf(size(countf)) !< Offset for file filespace_id
  integer(HID_T) :: dset_id, filespace_id, memspace_id
  integer :: errcode, rank
  logical :: exists

  PUSH_SUB(hdf5_write_int_lowlevel)

  if (present(error)) error = 0
  call h5lexists_f(loc_id, dset_name, exists, errcode)
  if (.not.exists .and. present(offsetf)) then
    ! FHJ: The developer should manually create the dataset, because countf
    ! is *NOT* the total size of the dataset, so we don`t have enough into to
    ! create the dataset.
    call die('Internal error: cannot automatically create dataset with "offsetf" argument.', &
      only_root_writes=.true.)
  endif

  ! FHJ: Create or get file dataspace
  if (present(offsetf)) then
    hoffsetf(:) = offsetf(:)
  else
    hoffsetf(:) = 0
  endif
  hcountf(:) = countf(:)
  rank = size(hcountf)
  if (any(countf<1)) then
    rank = 0
    ! Note: hcountf and hoffsetf are not referenced if rank==0
  endif

  if (.not.exists) then
    ! FHJ: Create file dataspace
    call safe_h5screate_simple(rank, hcountf, filespace_id, error)
    if (present(error)) then; if (error/=0) return; endif
    ! FHJ: Create dataset
    call safe_h5dcreate(loc_id, dset_name, H5T_NATIVE_INTEGER, filespace_id, dset_id, error)
    if (present(error)) then; if (error/=0) return; endif
  else
    ! FHJ: Open dataset
    call safe_h5dopen(loc_id, dset_name, dset_id, error)
    if (present(error)) then; if (error/=0) return; endif
    ! FHJ: Get file dataspace
    call safe_h5dget_space(dset_id, filespace_id, error)
    if (present(error)) then; if (error/=0) return; endif
  endif

  ! FHJ: Select hyperslab from file
  if (rank>0) then
    call safe_h5sselect_hyperslab(filespace_id, H5S_SELECT_SET_F, hoffsetf, hcountf, error)
    if (present(error)) then; if (error/=0) return; endif
  endif

  ! FHJ: Create flat memory filespace_id
  hcountm(1) = max(1, product(countf))
  call safe_h5screate_simple(1, hcountm, memspace_id, error)
  if (present(error)) then; if (error/=0) return; endif

  ! FHJ: write filespace_id
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, buf, hcountm, errcode, &
    memspace_id, filespace_id)
  call hdf5_check_error('h5dwrite_f', errcode, error)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5sclose(memspace_id, errcode)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5sclose(filespace_id, errcode)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5dclose(dset_id, errcode)
  if (present(error)) then; if (error/=0) return; endif

  POP_SUB(hdf5_write_int_lowlevel)

end subroutine hdf5_write_int_lowlevel


!> write double to/from an HDF5 file.
subroutine hdf5_write_double_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  real(DP), intent(in), dimension(*) :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)


  integer(HSIZE_T) :: hcountf(size(countf)) !< Count for file dataspace
  integer(HSIZE_T) :: hcountm(1) !< Count for memory dataspace
  integer(HSIZE_T) :: hoffsetf(size(countf)) !< Offset for file filespace_id
  integer(HID_T) :: dset_id, filespace_id, memspace_id
  integer :: errcode, rank
  logical :: exists

  PUSH_SUB(hdf5_write_double_lowlevel)

  if (present(error)) error = 0
  call h5lexists_f(loc_id, dset_name, exists, errcode)
  if (.not.exists .and. present(offsetf)) then
    ! FHJ: The developer should manually create the dataset, because countf
    ! is *NOT* the total size of the dataset, so we don`t have enough into to
    ! create the dataset.
    call die('Internal error: cannot automatically create dataset with "offsetf" argument.', &
      only_root_writes=.true.)
  endif

  ! FHJ: Create or get file dataspace
  if (present(offsetf)) then
    hoffsetf(:) = offsetf(:)
  else
    hoffsetf(:) = 0
  endif
  hcountf(:) = countf(:)
  rank = size(hcountf)
  if (any(countf<1)) then
    rank = 0
    ! Note: hcountf and hoffsetf are not referenced if rank==0
  endif

  if (.not.exists) then
    ! FHJ: Create file dataspace
    call safe_h5screate_simple(rank, hcountf, filespace_id, error)
    if (present(error)) then; if (error/=0) return; endif
    ! FHJ: Create dataset
    call safe_h5dcreate(loc_id, dset_name, H5T_NATIVE_DOUBLE, filespace_id, dset_id, error)
    if (present(error)) then; if (error/=0) return; endif
  else
    ! FHJ: Open dataset
    call safe_h5dopen(loc_id, dset_name, dset_id, error)
    if (present(error)) then; if (error/=0) return; endif
    ! FHJ: Get file dataspace
    call safe_h5dget_space(dset_id, filespace_id, error)
    if (present(error)) then; if (error/=0) return; endif
  endif

  ! FHJ: Select hyperslab from file
  if (rank>0) then
    call safe_h5sselect_hyperslab(filespace_id, H5S_SELECT_SET_F, hoffsetf, hcountf, error)
    if (present(error)) then; if (error/=0) return; endif
  endif

  ! FHJ: Create flat memory filespace_id
  hcountm(1) = max(1, product(countf))
  call safe_h5screate_simple(1, hcountm, memspace_id, error)
  if (present(error)) then; if (error/=0) return; endif

  ! FHJ: write filespace_id
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buf, hcountm, errcode, &
    memspace_id, filespace_id)
  call hdf5_check_error('h5dwrite_f', errcode, error)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5sclose(memspace_id, errcode)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5sclose(filespace_id, errcode)
  if (present(error)) then; if (error/=0) return; endif
  call safe_h5dclose(dset_id, errcode)
  if (present(error)) then; if (error/=0) return; endif

  POP_SUB(hdf5_write_double_lowlevel)

end subroutine hdf5_write_double_lowlevel



!########################################
!# LOW LEVEL routine -- LOGICAL version #
!########################################


!> read logical to/from an HDF5 file.
!! Note that this just maps the logical data to integers, and calls hdf5_read_int_lowlevel.
subroutine hdf5_read_logical_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  logical, intent(inout), dimension(*) :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)

  integer :: buf_int(max(1_HSIZE_T, product(countf)))
  integer(HSIZE_T) :: siz

  PUSH_SUB(hdf5_read_logical_lowlevel)

  siz = max(1_HSIZE_T, product(countf))
  call hdf5_read_int_lowlevel(loc_id, dset_name, countf, buf_int, error=error, offsetf=offsetf)
  buf(1:siz) = buf_int(:) /= 0

  POP_SUB(hdf5_read_logical_lowlevel)

end subroutine hdf5_read_logical_lowlevel

!> write logical to/from an HDF5 file.
!! Note that this just maps the logical data to integers, and calls hdf5_write_int_lowlevel.
subroutine hdf5_write_logical_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  logical, intent(in), dimension(*) :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)

  integer :: buf_int(max(1_HSIZE_T, product(countf)))
  integer(HSIZE_T) :: siz

  PUSH_SUB(hdf5_write_logical_lowlevel)

  siz = max(1_HSIZE_T, product(countf))
  where(buf(1:siz))
    buf_int = 1
  elsewhere
    buf_int = 0
  endwhere
  call hdf5_write_int_lowlevel(loc_id, dset_name, countf, buf_int, error=error, offsetf=offsetf)

  POP_SUB(hdf5_write_logical_lowlevel)

end subroutine hdf5_write_logical_lowlevel




!> read ('scalar', 'double') to/from an HDF5 file.
!! We map a scalar array of type double to another double array by adding
!! an extra dimension of size=1. Adding this extra dimension makes the shape
!! of the array consistent between scalar datasets of type complex and real.
subroutine hdf5_read_scalar_double_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  real(DP), intent(inout), dimension(*), target :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)

  integer :: rank_double
  integer(HSIZE_T) :: size_double
  integer(HSIZE_T) :: countf_double(size(countf)+1)
  integer :: offsetf_double(size(countf)+1)

  PUSH_SUB(hdf5_read_scalar_double_lowlevel)

  countf_double(1) = 1
  countf_double(2:) = countf
  offsetf_double(:) = 0
  if (any(countf<1)) then
    ! Buf is a 0-d array, i.e., a number, stored in the HDF5 file as a rank-1 array.
    rank_double = 1
    size_double = 1
  else
    ! Buf is a n-d array, stored in the HDF5 file as a rank-(n+1) array.
    rank_double = 1 + size(countf)
    size_double = 1_HSIZE_T * product(countf)
  endif

  ! Can pass double array `buf` directly when calling hdf5_read_double_lowlevel
  if (present(offsetf)) then
    offsetf_double(2:) = offsetf
    call hdf5_read_double_lowlevel(loc_id, dset_name, countf_double(1:rank_double), &
      buf, error=error, offsetf=offsetf_double(1:rank_double))
  else
    call hdf5_read_double_lowlevel(loc_id, dset_name, countf_double(1:rank_double), &
      buf, error=error)
  endif

  POP_SUB(hdf5_read_scalar_double_lowlevel)

end subroutine hdf5_read_scalar_double_lowlevel
!> read ('scalar', 'complex') to/from an HDF5 file.
!! We map a scalar array of type complex to another double array by adding
!! an extra dimension of size=2. Adding this extra dimension makes the shape
!! of the array consistent between scalar datasets of type complex and real.
subroutine hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  complex(DPC), intent(inout), dimension(*), target :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)

  integer :: rank_double
  integer(HSIZE_T) :: size_double
  integer(HSIZE_T) :: countf_double(size(countf)+1)
  integer :: offsetf_double(size(countf)+1)
  real(DP), pointer :: buf_double(:)

  PUSH_SUB(hdf5_read_scalar_complex_lowlevel)

  countf_double(1) = 2
  countf_double(2:) = countf
  offsetf_double(:) = 0
  if (any(countf<1)) then
    ! Buf is a 0-d array, i.e., a number, stored in the HDF5 file as a rank-1 array.
    rank_double = 1
    size_double = 2
  else
    ! Buf is a n-d array, stored in the HDF5 file as a rank-(n+1) array.
    rank_double = 1 + size(countf)
    size_double = 2_HSIZE_T * product(countf)
  endif

  ! Need to cast complex array `buf` into a double array before calling hdf5_read_double_lowlevel
  call c_f_pointer(c_loc(buf), buf_double, [size_double])
  if (present(offsetf)) then
    offsetf_double(2:) = offsetf
    call hdf5_read_double_lowlevel(loc_id, dset_name, countf_double(1:rank_double), &
      buf_double, error=error, offsetf=offsetf_double(1:rank_double))
  else
    call hdf5_read_double_lowlevel(loc_id, dset_name, countf_double(1:rank_double), &
      buf_double, error=error)
  endif

  POP_SUB(hdf5_read_scalar_complex_lowlevel)

end subroutine hdf5_read_scalar_complex_lowlevel

!> read complex to/from an HDF5 file.
!! This subroutine just calls hdf5_read_scalar_complex_lowlevel directly, which appends an extra dimension.
!! The behavior of this subroutine is slightly different from that of hdf5_read_double_lowlevel,
!! which does not add an extra dimension. That is because it is possible to
!! natively store real numbers in HDF5 files, but not complex numbers.
!! So, while both hdf5_read_scalar_complex_lowlevel and
!! hdf5_read_complex_lowlevel add an extra dimension to
!! the array, only hdf5_read_scalar_double_lowlevel adds the
!! extra dimension (to make the number of dimensions the same as that in
!! get_subroutine_label(read_or_write, ('scalar','double'), 'lowlevel'). However,
!! hdf5_read_double_lowlevel will not add
!! an extra dimension.
subroutine hdf5_read_complex_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  complex(DPC), intent(inout), dimension(*), target :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)

  PUSH_SUB(hdf5_read_complex_lowlevel)
  call hdf5_read_scalar_complex_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  POP_SUB(hdf5_read_complex_lowlevel)

end subroutine hdf5_read_complex_lowlevel

!> write ('scalar', 'double') to/from an HDF5 file.
!! We map a scalar array of type double to another double array by adding
!! an extra dimension of size=1. Adding this extra dimension makes the shape
!! of the array consistent between scalar datasets of type complex and real.
subroutine hdf5_write_scalar_double_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  real(DP), intent(in), dimension(*), target :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)

  integer :: rank_double
  integer(HSIZE_T) :: size_double
  integer(HSIZE_T) :: countf_double(size(countf)+1)
  integer :: offsetf_double(size(countf)+1)

  PUSH_SUB(hdf5_write_scalar_double_lowlevel)

  countf_double(1) = 1
  countf_double(2:) = countf
  offsetf_double(:) = 0
  if (any(countf<1)) then
    ! Buf is a 0-d array, i.e., a number, stored in the HDF5 file as a rank-1 array.
    rank_double = 1
    size_double = 1
  else
    ! Buf is a n-d array, stored in the HDF5 file as a rank-(n+1) array.
    rank_double = 1 + size(countf)
    size_double = 1_HSIZE_T * product(countf)
  endif

  ! Can pass double array `buf` directly when calling hdf5_write_double_lowlevel
  if (present(offsetf)) then
    offsetf_double(2:) = offsetf
    call hdf5_write_double_lowlevel(loc_id, dset_name, countf_double(1:rank_double), &
      buf, error=error, offsetf=offsetf_double(1:rank_double))
  else
    call hdf5_write_double_lowlevel(loc_id, dset_name, countf_double(1:rank_double), &
      buf, error=error)
  endif

  POP_SUB(hdf5_write_scalar_double_lowlevel)

end subroutine hdf5_write_scalar_double_lowlevel
!> write ('scalar', 'complex') to/from an HDF5 file.
!! We map a scalar array of type complex to another double array by adding
!! an extra dimension of size=2. Adding this extra dimension makes the shape
!! of the array consistent between scalar datasets of type complex and real.
subroutine hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  complex(DPC), intent(in), dimension(*), target :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)

  integer :: rank_double
  integer(HSIZE_T) :: size_double
  integer(HSIZE_T) :: countf_double(size(countf)+1)
  integer :: offsetf_double(size(countf)+1)
  real(DP), pointer :: buf_double(:)

  PUSH_SUB(hdf5_write_scalar_complex_lowlevel)

  countf_double(1) = 2
  countf_double(2:) = countf
  offsetf_double(:) = 0
  if (any(countf<1)) then
    ! Buf is a 0-d array, i.e., a number, stored in the HDF5 file as a rank-1 array.
    rank_double = 1
    size_double = 2
  else
    ! Buf is a n-d array, stored in the HDF5 file as a rank-(n+1) array.
    rank_double = 1 + size(countf)
    size_double = 2_HSIZE_T * product(countf)
  endif

  ! Need to cast complex array `buf` into a double array before calling hdf5_write_double_lowlevel
  call c_f_pointer(c_loc(buf), buf_double, [size_double])
  if (present(offsetf)) then
    offsetf_double(2:) = offsetf
    call hdf5_write_double_lowlevel(loc_id, dset_name, countf_double(1:rank_double), &
      buf_double, error=error, offsetf=offsetf_double(1:rank_double))
  else
    call hdf5_write_double_lowlevel(loc_id, dset_name, countf_double(1:rank_double), &
      buf_double, error=error)
  endif

  POP_SUB(hdf5_write_scalar_complex_lowlevel)

end subroutine hdf5_write_scalar_complex_lowlevel

!> write complex to/from an HDF5 file.
!! This subroutine just calls hdf5_write_scalar_complex_lowlevel directly, which appends an extra dimension.
!! The behavior of this subroutine is slightly different from that of hdf5_write_double_lowlevel,
!! which does not add an extra dimension. That is because it is possible to
!! natively store real numbers in HDF5 files, but not complex numbers.
!! So, while both hdf5_write_scalar_complex_lowlevel and
!! hdf5_write_complex_lowlevel add an extra dimension to
!! the array, only hdf5_write_scalar_double_lowlevel adds the
!! extra dimension (to make the number of dimensions the same as that in
!! get_subroutine_label(read_or_write, ('scalar','double'), 'lowlevel'). However,
!! hdf5_write_double_lowlevel will not add
!! an extra dimension.
subroutine hdf5_write_complex_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimension.
  !! Pass (/-1/) to read/write a scalar rank-0 array.
  integer(HSIZE_T), intent(in) :: countf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  complex(DPC), intent(in), dimension(*), target :: buf
  integer, intent(out), optional :: error !< error code
  !> Offset when reading dataset from file.
  integer, intent(in), optional :: offsetf(:)

  PUSH_SUB(hdf5_write_complex_lowlevel)
  call hdf5_write_scalar_complex_lowlevel(loc_id, dset_name, countf, buf, error, offsetf)
  POP_SUB(hdf5_write_complex_lowlevel)

end subroutine hdf5_write_complex_lowlevel

#endif
end module hdf5_io_data_m
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
