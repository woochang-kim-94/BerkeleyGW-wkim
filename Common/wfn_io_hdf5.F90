!>=========================================================================
!!
!! Module:
!!
!! (1) wfn_io_hdf5_m     Originally by JIM     Last Modified 4/25/2012 (JIM)
!!
!!     Routines to read and write wavefunctions in HDF5 format.
!!     The code is generated through repeated inclusion of a file with
!!     different preprocessor definitions each time. You are not expected to
!!     understand this. Consult the resulting .p.f file for clarity.
!!
!!=========================================================================

#include "f_defs.h"

module wfn_io_hdf5_m
  use global_m
#ifdef HDF5
  use hdf5_io_m
  use hdf5

  implicit none

  private

  public :: &
    read_hdf5_header_type   , &
    write_hdf5_header_type  , &
    read_hdf5_gvectors      , &
    write_hdf5_gvectors     , &
    read_hdf5_wfn_gvectors  , &
    write_hdf5_wfn_gvectors , &
    read_hdf5_band_real     , &
    write_hdf5_band_real    , &
    read_hdf5_band_complex  , &
    write_hdf5_band_complex , &
    read_hdf5_band          , &
    write_hdf5_band         , &
    read_hdf5_bands_block   , &
    setup_hdf5_mf_file      , &
    setup_hdf5_wfn_file     , &
    read_hdf5_mf_header     , &
    write_hdf5_mf_header

  interface read_hdf5_band
    module procedure read_hdf5_band_real, read_hdf5_band_complex
  end interface

  interface write_hdf5_band
    module procedure write_hdf5_band_real, write_hdf5_band_complex
  end interface

contains

!===============================================================================
#define READ
!
#define TEMP_HEADER
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_HEADER
!
#define TEMP_WFN_GVEC
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_WFN_GVEC
!
#define TEMP_WFN_DATA
#define TEMP_COMPLEX
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_COMPLEX
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_WFN_DATA
!===============================================================================
#undef READ
!
#define TEMP_HEADER
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_HEADER
!
#define TEMP_WFN_GVEC
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_WFN_GVEC
!
#define TEMP_WFN_DATA
#define TEMP_COMPLEX
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_COMPLEX
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_WFN_DATA
!
#define TEMP_OTHER
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_OTHER
!===============================================================================


  !> Create the appropriate structures that hold info about the mean-field
  !! calculation in an HDF5 file. No data is actually written.
  subroutine setup_hdf5_mf_file(fname, create_file)
    character(len=*), intent(in) :: fname
    logical, intent(in), optional :: create_file

    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer :: error
    logical :: create_file_

    PUSH_SUB(setup_hdf5_mf_file)

    create_file_ = .true.
    if (present(create_file)) create_file_ = create_file
    if (create_file_) then
      call hdf5_create_file(fname, file_id)
    else
      call hdf5_open_file(fname, 'rw', file_id)
    endif

    call h5gcreate_f(file_id, '/mf_header', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/gspace', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/symmetry', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/crystal', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/kpoints', group_id, error)
    call h5gclose_f(group_id, error)

    call hdf5_close_file(file_id)

    POP_SUB(setup_hdf5_mf_file)

  end subroutine setup_hdf5_mf_file


  !> Create the appropriate structures that hold info about the WFN
  !! coefficients in an HDF5 file. No data is actually written.
  subroutine setup_hdf5_wfn_file(fname, iflavor, kp)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: iflavor
    type(kpoints), intent(in) :: kp

    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer(HID_T) :: dataspace_id
    integer(HID_T) :: dataset_id
    integer(HSIZE_T) :: a3(4)
    integer :: error

    PUSH_SUB(setup_hdf5_wfn_file)

    call setup_hdf5_mf_file(fname)

    call hdf5_open_file(fname, 'rw', file_id)
    call h5gcreate_f(file_id, '/wfns', group_id, error)
    call h5gclose_f(group_id, error)
    if (iflavor/=1 .and. iflavor/=2) then
      write(0,*) 'ERROR: got iflavor=', iflavor
      call die('Internal error: invalid flavor in setup_hdf5_wfn_file.', only_root_writes=.true.)
    endif
    a3(1) = iflavor
    a3(2) = sum(kp%ngk)
    a3(3) = kp%nspin*kp%nspinor
    a3(4) = kp%mnband
    CALL h5screate_simple_f(4, a3, dataspace_id, error)
    call h5dcreate_f(file_id, '/wfns/coeffs', H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
    call h5dclose_f(dataset_id, error)
    CALL h5sclose_f(dataspace_id, error)
    call hdf5_close_file(file_id)

    POP_SUB(setup_hdf5_wfn_file)

  end subroutine setup_hdf5_wfn_file


  !> A high-level wrapper for write_*_header* functions
  subroutine write_hdf5_mf_header(fname, mf)
    character(len=*), intent(in) :: fname
    type(mf_header_t), intent(in) :: mf

    character(len=3) :: sheader

    PUSH_SUB(write_hdf5_mf_header)

    sheader = mf%sheader
    call write_hdf5_header_type(fname, sheader, mf%iflavor, &
                                mf%kp, mf%gvec, mf%syms, mf%crys)

    POP_SUB(write_hdf5_mf_header)

  end subroutine write_hdf5_mf_header


  !> A high-level wrapper for read_*_header* functions
  !! Note: the optional fields are ignored for now.
  subroutine read_hdf5_mf_header(fname, mf, iflavor, sheader, warn, dont_warn_kgrid)
    character(len=*), intent(in) :: fname
    type(mf_header_t), intent(out) :: mf
    integer, intent(in), optional :: iflavor
    character(len=3), intent(in), optional :: sheader
    logical, intent(in), optional :: warn
    logical, intent(in), optional :: dont_warn_kgrid

    PUSH_SUB(read_hdf5_mf_header)

    if (present(sheader)) then
      mf%sheader = sheader
    else
      mf%sheader = 'GET'
    endif
    if (present(iflavor)) then
      mf%iflavor = iflavor
    else
      mf%iflavor = -1
    endif
    call read_hdf5_header_type(fname, mf%sheader, mf%iflavor, &
      !mf%kp, mf%gvec, mf%syms, mf%crys, version=mf%version, sdate=mf%sdate, stime=mf%stime)
      mf%kp, mf%gvec, mf%syms, mf%crys)

    !FHJ: FIXME - Implement in WFN file
    mf%sheader = 'WFN'
    mf%stime = 'N/A'
    mf%sdate = 'N/A'

    POP_SUB(read_hdf5_mf_header)

  end subroutine read_hdf5_mf_header

#endif

end module wfn_io_hdf5_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
