!>=========================================================================
!!
!!  Module:
!!
!!  epsread_hdf5_m     Originally by JRD     Last Modified 12/2014 (FHJ)
!!
!!    Routines to read header info and matrices from epsmat files in
!!    HDF5 format.
!!
!!=========================================================================

#include "f_defs.h"

module epsread_hdf5_m
#ifdef HDF5

  use global_m
  use hdf5
  use hdf5_io_m
  use scalapack_m
  !use io_utils_m
  use timing_m, only: timing => common_timing
  implicit none

  private

  public :: &
    read_eps_params_hdf5, &
    read_eps_grid_sizes_hdf5, &
    read_eps_subspace_info, &
    read_eps_matrix_diagonal_hdf5, &
    read_eps_qgrid_hdf5, &
    read_eps_sub_neig_of_q, &
    read_eps_freqgrid_hdf5, &
    read_eps_old_gvecs_hdf5, &
    read_eps_gvecsofq_hdf5, &
    read_vcoul_hdf5, &
    read_eps_matrix_col_freq_f_hdf5, &
    read_eps_matrix_col_f_hdf5, &
    read_eps_matrix_col_hdf5, &
    read_eps_matrix_ser_hdf5, &
    read_eps_matrix_ser_f_hdf5, &
    read_eps_matrix_par_hdf5, &
    read_eps_matrix_par_f_hdf5
#if defined MPI && defined USESCALAPACK
  public :: read_eps_matrix_par_hdf5_general, &
            read_eps_matrix_par_hdf5_general_f
#endif

contains


subroutine read_eps_params_hdf5(pol, name, nband)
  type(polarizability), intent(inout) :: pol
  character(len=*), intent(in) :: name
  integer, intent(out), optional :: nband

  integer(HID_T) :: file_id       ! File identifier
  integer :: error

  PUSH_SUB(read_eps_params_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  call hdf5_read_int(file_id, 'eps_header/params/matrix_type', pol%matrix_type)
  call hdf5_read_logical(file_id, 'eps_header/params/has_advanced', pol%has_advanced)
  call hdf5_read_int(file_id, 'eps_header/params/nmatrix', pol%nmatrix)
  call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', pol%matrix_flavor)
  call hdf5_read_int(file_id, 'eps_header/params/icutv', pol%icutv)
  call hdf5_read_double(file_id, 'eps_header/params/ecuts', pol%ecuts)
  if (present(nband)) &
    call hdf5_read_int(file_id, 'eps_header/params/nband', nband)
  call hdf5_read_double(file_id, 'eps_header/params/efermi', pol%efermi)
  pol%efermi = pol%efermi*ryd
  call hdf5_read_int(file_id, 'eps_header/params/intraband_flag', pol%intraband_flag)
  call hdf5_read_double(file_id, 'eps_header/params/intraband_overlap_min', pol%intraband_overlap_min)
  call hdf5_read_logical(file_id, 'eps_header/params/subsample', pol%subsample)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_params_hdf5)

end subroutine read_eps_params_hdf5

!===================================================================================

subroutine read_eps_grid_sizes_hdf5(ng, nq, ecuts, nfreq, nfreq_imag, nmtxmax, qgrid, freq_dep, name)
  integer, intent(out) :: ng
  integer, intent(out) :: nq
  real(DP), intent(out) :: ecuts
  integer, intent(out) :: nfreq
  integer, intent(out) :: nfreq_imag
  integer, intent(out) :: nmtxmax
  integer, intent(out) :: qgrid(3)
  integer, intent(out) :: freq_dep
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier
  integer :: version, flavor

  PUSH_SUB(read_eps_grid_sizes_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  call hdf5_require_version(file_id, 'eps_header/versionnumber', VER_EPS_HDF5, trim(name))
  call hdf5_require_version(file_id, 'mf_header/versionnumber', VER_WFN_HDF5, trim(name))
  call hdf5_require_flavor(file_id, 'eps_header/flavor', SCALARSIZE, trim(name))
  call hdf5_require_flavor(file_id, 'mf_header/flavor', SCALARSIZE, trim(name))

  call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq)
  call hdf5_read_int(file_id, 'mf_header/gspace/ng', ng)
  call hdf5_read_int(file_id, 'eps_header/freqs/nfreq', nfreq)
  if (hdf5_exists(file_id, 'eps_header/freqs/nfreq_imag')) then
    call hdf5_read_int(file_id, 'eps_header/freqs/nfreq_imag', nfreq_imag)
  else
    nfreq_imag = 0
  endif
  call hdf5_read_int(file_id, 'eps_header/freqs/freq_dep', freq_dep)
  call hdf5_read_int(file_id, 'eps_header/gspace/nmtx_max', nmtxmax)
  call hdf5_read_int_array(file_id, 'eps_header/qpoints/qgrid', (/3/), qgrid)
  call hdf5_read_double(file_id, 'eps_header/params/ecuts', ecuts)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_grid_sizes_hdf5)

end subroutine read_eps_grid_sizes_hdf5

!====================================================================================

subroutine read_eps_subspace_info(is_subspace, matrix_in_subspace_basis, keep_full_eps_static, &
                                  neig_sub, chi_eigenvalue_cutoff, name)
  logical, intent(out) :: is_subspace, matrix_in_subspace_basis, keep_full_eps_static
  integer, intent(out) :: neig_sub
  real(dp) :: chi_eigenvalue_cutoff
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier

  PUSH_SUB(read_eps_subspace_info)

  is_subspace = .false.
  matrix_in_subspace_basis = .false.
  keep_full_eps_static = .false.
  neig_sub = 0
  chi_eigenvalue_cutoff = 0.0d0

  call hdf5_open_file(trim(name), 'r', file_id)
  ! check if the subspace flag exists
  if (hdf5_exists(file_id, 'eps_header/params/subspace')) then
    call hdf5_read_logical(file_id, 'eps_header/params/subspace', is_subspace)
    if (is_subspace) then
      call hdf5_read_logical(file_id, 'eps_header/subspace/matrix_in_subspace_basis', matrix_in_subspace_basis)
      call hdf5_read_logical(file_id, 'eps_header/subspace/keep_full_eps_static', keep_full_eps_static)
      call hdf5_read_double(file_id, 'eps_header/subspace/eps_eigenvalue_cutoff', chi_eigenvalue_cutoff)
      call hdf5_read_int(file_id, 'eps_header/subspace/neig_max',  neig_sub)
    end if
  end if
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_subspace_info)

end subroutine

!====================================================================================

subroutine read_eps_matrix_diagonal_hdf5(nmtx, iq, epsdiag, name)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  real(DP), intent(out) :: epsdiag(:,:) !< (matrix_flavor,nmtx)
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier
  integer :: matrix_flavor

  PUSH_SUB(read_eps_matrix_diagonal_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', matrix_flavor)
  if (matrix_flavor/=size(epsdiag,1)) then
    write(0,*) 'ERROR: Got size(epsdiag,1)=',size(epsdiag,1),', but matrix_flavor=',matrix_flavor
    call die('Internal error in read_eps_matrix_diagonal_hdf5', &
      only_root_writes=.true.)
  endif
  call hdf5_read_double_hyperslab(file_id, 'mats/matrix-diagonal', &
    (/matrix_flavor,nmtx,1/), (/0,0,iq-1/), epsdiag)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_matrix_diagonal_hdf5)

end subroutine read_eps_matrix_diagonal_hdf5

!====================================================================================

subroutine read_eps_qgrid_hdf5(nq, qpts, nmtx, name)
  integer, intent(in) :: nq
  real(DP), intent(inout) :: qpts(:,:) !< (3,nq)
  integer, intent(out) :: nmtx(:) !< (nq)
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier

  PUSH_SUB(read_eps_qgrid_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  call hdf5_read_double_array(file_id, 'eps_header/qpoints/qpts', (/3,nq/), qpts)
  call hdf5_read_int_array(file_id, 'eps_header/gspace/nmtx', (/nq/), nmtx)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_qgrid_hdf5)

end subroutine read_eps_qgrid_hdf5

!===================================================================================

subroutine read_eps_sub_neig_of_q(nq, neig, name)

  integer, intent(in)  :: nq
  integer, intent(out) :: neig(:) !< (nq)
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier

  PUSH_SUB(read_eps_sub_neig_of_q)

  call hdf5_open_file(trim(name), 'r', file_id)
  call hdf5_read_int_array(file_id, 'eps_header/subspace/neig', (/nq/), neig)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_sub_neig_of_q)

end subroutine read_eps_sub_neig_of_q

!===================================================================================

subroutine read_eps_freqgrid_hdf5(nfreq, dFreqGrid, dFreqBrd, name)
  integer, intent(in) :: nfreq
  real(DP), intent(out) :: dFreqGrid(:) !< (nfreq)
  complex(DPC), intent(out) :: dFreqBrd(:) !< (nfreq)
  character(len=*), intent(in) :: name

  real(DP) :: freqs_tmp(2,nfreq)
  integer :: iw
  integer(HID_T) :: file_id       ! File identifier

  PUSH_SUB(read_eps_freqgrid_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  call hdf5_read_double_array(file_id, 'eps_header/freqs/freqs', (/2,nfreq/), freqs_tmp)
  do iw=1,nfreq
    dFreqGrid(iw) = freqs_tmp(1,iw)
    dFreqBrd(iw) = CMPLX(0,freqs_tmp(2,iw))
  enddo
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_freqgrid_hdf5)

end subroutine read_eps_freqgrid_hdf5

!=================================================================================

subroutine read_eps_old_gvecs_hdf5(ng, gvecs, name)
  integer, intent(in) :: ng
  integer, intent(out) :: gvecs(:,:) !< (3,ng)
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier

  PUSH_SUB(read_eps_old_gvecs_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  call hdf5_read_int_array(file_id, 'mf_header/gspace/components', (/3,ng/), gvecs)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_old_gvecs_hdf5)

end subroutine read_eps_old_gvecs_hdf5

!====================================================================================

subroutine read_eps_gvecsofq_hdf5(ng, gind_eps2rho, gind_rho2eps, ekin, iq, name)
  integer, intent(in) :: ng !< Number of G-vectors
  integer, intent(out) :: gind_eps2rho(:) !< (ng)
  integer, intent(out) :: gind_rho2eps(:) !< (ng)
  real(DP), intent(out) :: ekin(:) !< (ng)
  integer, intent(in) :: iq
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id
  integer :: countf(2), offsetf(2), ng_real

  PUSH_SUB(read_eps_gvecsofq_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  ng_real = min(ng, size(gind_eps2rho))
  countf(:) = [ng_real, 1]
  offsetf(:) = [0, iq-1]
  call hdf5_read_int_hyperslab(file_id, 'eps_header/gspace/gind_eps2rho', &
    countf, offsetf, gind_eps2rho)
  call hdf5_read_int_hyperslab(file_id, 'eps_header/gspace/gind_rho2eps', &
    countf, offsetf, gind_rho2eps)
  call hdf5_read_double_hyperslab(file_id, 'eps_header/gspace/ekin', &
    countf, offsetf, ekin)

  call hdf5_close_file(file_id)

  POP_SUB(read_eps_gvecsofq_hdf5)

end subroutine read_eps_gvecsofq_hdf5

!====================================================================================

subroutine read_vcoul_hdf5(vcoul, iq, name)
  real(DP), intent(out) :: vcoul(:) !< (nmtx(iq))
  integer, intent(in) :: iq
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id
  integer :: nmtx(1)
  integer :: countf(2), offsetf(2)

  PUSH_SUB(read_vcoul_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  call hdf5_read_int_hyperslab(file_id, 'eps_header/gspace/nmtx', (/1/), (/iq-1/), nmtx)
  if (size(vcoul)<nmtx(1)) &
    call die('Internal error: array vcoul too small in read_vcoul_hdf.', only_root_writes=.true.)
  countf(:) = (/nmtx(1), 1/)
  offsetf(:) = (/0, iq-1/)
  call hdf5_read_double_hyperslab(file_id, 'eps_header/gspace/vcoul', &
    countf, offsetf, vcoul)
  call hdf5_close_file(file_id)

  POP_SUB(read_vcoul_hdf5)

end subroutine read_vcoul_hdf5

!===========================================================================================

!> Reads a specific column (Gp) and frequency of the FF dielectric matrix.
!! FIXME: Ideally, we should have a version that reads all Gs and Gps for a
!! given frequency.
subroutine read_eps_matrix_col_freq_f_hdf5(retarded, nfreq, ifreq, igp, nmtx, iq, is, name, advanced)
  integer, intent(in) :: iq
  integer, intent(in) :: is
  integer, intent(in) :: nfreq
  integer, intent(in) :: ifreq
  integer, intent(in) :: nmtx
  integer, intent(in) :: igp
  complex(DPC), intent(out) :: retarded(:) !< (nmtx)
  character(len=*), intent(in) :: name
  complex(DPC), optional, intent(out) :: advanced(:) !< (nmtx)

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: data_id       ! Property list identifier
  integer(HID_T) :: dataspace        ! Property list identifier
  integer(HID_T) :: memspace        ! Property list identifier
  integer :: error
  integer(HSIZE_T) :: count(6), offset(6)

  real(DP), allocatable :: data(:,:,:,:,:,:)
  integer :: nmatrix_per_spin, nspin, buf_sz, version
  logical :: has_advanced

  PUSH_SUB(read_eps_matrix_col_freq_f_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)

  ! FHJ: the default is never to read the advanced matrix, unless the file
  ! version is <3 (on which case we didn`t store the Coulomb interaction)
  call hdf5_read_int(file_id, 'eps_header/versionnumber', version)
  call hdf5_read_int(file_id, 'eps_header/params/nmatrix', nmatrix_per_spin)
  call hdf5_read_int(file_id, 'mf_header/kpoints/nspin', nspin)
  nmatrix_per_spin = nmatrix_per_spin / nspin
  has_advanced = .false.
  buf_sz = 1
  if (version<3) then
    call hdf5_read_logical(file_id, 'eps_header/params/has_advanced', has_advanced)
    if (present(advanced) .and. .not.has_advanced) &
      call die('Inconsistent epsmat file: version<3, but no advanced matrix', only_root_writes=.true.)
  endif
  if (has_advanced) buf_sz = 2

  call h5dopen_f(file_id, 'mats/matrix', data_id, error)
  call h5dget_space_f(data_id,dataspace,error)

  count(:) = (/2, nmtx, 1, 1, buf_sz, 1/)
  call h5screate_simple_f(6, count, memspace, error)
  offset(:) = (/0, 0, igp-1, ifreq-1, nmatrix_per_spin*(is-1), iq-1/)
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

! XXX - this is bad because we double memory usage. Oh well.
  SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
  call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
  retarded(1:nmtx) = CMPLX(data(1,1:nmtx,1,1,1,1),data(2,1:nmtx,1,1,1,1))
  if (has_advanced .and. present(advanced)) then
    advanced(1:nmtx) = CMPLX(data(1,1:nmtx,1,1,2,1),data(2,1:nmtx,1,1,2,1))
  endif
  SAFE_DEALLOCATE(data)
  call h5sclose_f(memspace, error)

  if (.not.has_advanced .and. present(advanced)) call get_advanced_from_retarded()

  call h5sclose_f(dataspace, error)
  call h5dclose_f(data_id, error)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_matrix_col_freq_f_hdf5)

contains

  !> FHJ: This routine does the magic of generating advanced epsilon/chi out of
  !! the retarded epsilon/chi.
  !! We use the fact that W_A = (W_R)^H to write
  !! epsinv_A = (epsinv_R * v)^H * (1/v)
  subroutine get_advanced_from_retarded()
    real(DP) :: vcoul(nmtx)
    integer :: matrix_type, ig

    PUSH_SUB(read_eps_matrix_col_f_hdf5.get_advanced_from_retarded)

    call hdf5_read_int(file_id, 'eps_header/params/matrix_type', matrix_type)
    vcoul(:) = 1d0
    ! FHJ: There is no Coulomb interaction in a chimat file.
    if (matrix_type<2) then
      call hdf5_read_double_hyperslab(file_id, 'eps_header/gspace/vcoul', &
        (/nmtx,1/), (/0,iq-1/), vcoul)
    endif

    count(:) = (/2, 1, nmtx, 1, buf_sz, 1/)
    call h5screate_simple_f(6, count, memspace, error)
    offset(:) = (/0, igp-1, 0, ifreq-1, nmatrix_per_spin*(is-1), iq-1/)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

    ! XXX - this is bad because we double memory usage. Oh well.
    SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))

    ! JRD XXX These reads have horrendous locality
    ! FHJ: you are right, they have terrible I/O locality. But this is only
    ! for the serial code.

    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    ! FHJ: Note the complex conjugation here
    advanced(1:nmtx) = CMPLX(data(1,1,1:nmtx,1,1,1),-data(2,1,1:nmtx,1,1,1))
    ! FHJ: Correct asymmetry from Coulomb interaction
    do ig = 1, nmtx
      advanced(ig) = advanced(ig) * (vcoul(ig)/vcoul(igp))
    enddo
    SAFE_DEALLOCATE(data)
    call h5sclose_f(memspace, error)

    POP_SUB(read_eps_matrix_col_f_hdf5.get_advanced_from_retarded)

  end subroutine get_advanced_from_retarded


end subroutine read_eps_matrix_col_freq_f_hdf5

!===========================================================================================

!> Reads a specific column (Gp) and all frequencies of the FF dielectric matrix.
!! FIXME: Because of the new way the file is stored, this function is very
!! ineficient. Use `read_eps_matrix_col_freq_f_hdf5` if possible, for now.
subroutine read_eps_matrix_col_f_hdf5(retarded, nFreq, igp, nmtx, iq, is, name, advanced)
  integer, intent(in) :: iq
  integer, intent(in) :: is
  integer, intent(in) :: nFreq
  integer, intent(in) :: nmtx
  integer, intent(in) :: igp
  complex(DPC), intent(out) :: retarded(:,:) !< (nmtx,nFreq)
  character(len=*), intent(in) :: name
  complex(DPC), optional, intent(out) :: advanced(:,:) !< (nmtx,nFreq)

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: data_id       ! Property list identifier
  integer(HID_T) :: dataspace        ! Property list identifier
  integer(HID_T) :: memspace        ! Property list identifier
  integer :: error
  integer(HSIZE_T) :: count(6), offset(6)

  real(DP), allocatable :: data(:,:,:,:,:,:)
  integer :: nmatrix_per_spin, nspin, buf_sz, version
  logical :: has_advanced

  PUSH_SUB(read_eps_matrix_col_f_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)

  ! FHJ: the default is never to read the advanced matrix, unless the file
  ! version is <3 (on which case we didn`t store the Coulomb interaction)
  call hdf5_read_int(file_id, 'eps_header/versionnumber', version)
  call hdf5_read_int(file_id, 'eps_header/params/nmatrix', nmatrix_per_spin)
  call hdf5_read_int(file_id, 'mf_header/kpoints/nspin', nspin)
  nmatrix_per_spin = nmatrix_per_spin / nspin
  has_advanced = .false.
  buf_sz = 1
  if (version<3) then
    call hdf5_read_logical(file_id, 'eps_header/params/has_advanced', has_advanced)
    if (present(advanced) .and. .not.has_advanced) &
      call die('Inconsistent epsmat file: version<3, but no advanced matrix', only_root_writes=.true.)
  endif
  if (has_advanced) buf_sz = 2

  call h5dopen_f(file_id, 'mats/matrix', data_id, error)
  call h5dget_space_f(data_id,dataspace,error)

  count(:) = (/2, nmtx, 1, nFreq, buf_sz, 1/)
  call h5screate_simple_f(6, count, memspace, error)
  offset(:) = (/0, 0, igp-1, 0, nmatrix_per_spin*(is-1), iq-1/)
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

! XXX - this is bad because we double memory usage. Oh well.
  SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
  call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
  retarded(1:nmtx, 1:nFreq) = CMPLX(data(1,1:nmtx,1,1:nFreq,1,1),data(2,1:nmtx,1,1:nFreq,1,1))
  if (has_advanced .and. present(advanced)) then
    advanced(1:nmtx,1:nFreq) = CMPLX(data(1,1:nmtx,1,1:nFreq,2,1),data(2,1:nmtx,1,1:nFreq,2,1))
  endif
  SAFE_DEALLOCATE(data)
  call h5sclose_f(memspace, error)

  if (.not.has_advanced .and. present(advanced)) call get_advanced_from_retarded()

  call h5sclose_f(dataspace, error)
  call h5dclose_f(data_id, error)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_matrix_col_f_hdf5)

contains

  !> FHJ: This routine does the magic of generating advanced epsilon/chi out of
  !! the retarded epsilon/chi.
  !! We use the fact that W_A = (W_R)^H to write
  !! epsinv_A = (epsinv_R * v)^H * (1/v)
  subroutine get_advanced_from_retarded()
    real(DP) :: vcoul(nmtx)
    integer :: matrix_type, ig

    PUSH_SUB(read_eps_matrix_col_f_hdf5.get_advanced_from_retarded)

    call hdf5_read_int(file_id, 'eps_header/params/matrix_type', matrix_type)
    vcoul(:) = 1d0
    ! FHJ: There is no Coulomb interaction in a chimat file.
    if (matrix_type<2) then
      call hdf5_read_double_hyperslab(file_id, 'eps_header/gspace/vcoul', &
        (/nmtx,1/), (/0,iq-1/), vcoul)
    endif

    count(:) = (/2, 1, nmtx, nFreq, buf_sz, 1/)
    call h5screate_simple_f(6, count, memspace, error)
    offset(:) = (/0, igp-1, 0, 0, nmatrix_per_spin*(is-1), iq-1/)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

    ! XXX - this is bad because we double memory usage. Oh well.
    SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))

    ! JRD XXX These reads have horrendous locality
    ! FHJ: you are right, they have terrible I/O locality. But this is only
    ! for the serial code.

    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    ! FHJ: Note the complex conjugation here
    advanced(1:nmtx,1:nFreq) = CMPLX(data(1,1,1:nmtx,1:nFreq,1,1),-data(2,1,1:nmtx,1:nFreq,1,1))
    ! FHJ: Correct asymmetry from Coulomb interaction
    do ig = 1, nmtx
      advanced(ig,1:nFreq) = advanced(ig,1:nFreq) * (vcoul(ig)/vcoul(igp))
    enddo
    SAFE_DEALLOCATE(data)
    call h5sclose_f(memspace, error)

    POP_SUB(read_eps_matrix_col_f_hdf5.get_advanced_from_retarded)

  end subroutine get_advanced_from_retarded


end subroutine read_eps_matrix_col_f_hdf5

!===========================================================================================

subroutine read_eps_matrix_col_hdf5(eps,j,ncol,nmtx,iq,is,name)
  integer, intent(in) :: iq
  integer, intent(in) :: is
  integer, intent(in) :: nmtx
  integer, intent(in) :: j
  integer, intent(in) :: ncol     ! Reading a block of column from j to j+ncol-1
  SCALAR, intent(out) :: eps(:,:) !< (nmtx,ncol)
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: data_id       ! Property list identifier
  integer(HID_T) :: dataspace        ! Property list identifier
  integer(HID_T) :: memspace        ! Property list identifier
  integer :: error, rank
  integer(HSIZE_T) :: count(6), offset(6)

  real(DP), allocatable :: data(:,:,:,:,:,:)

  PUSH_SUB(read_eps_matrix_col_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  call h5dopen_f(file_id, 'mats/matrix', data_id, error)
  call h5dget_space_f(data_id, dataspace, error)

  rank = 6
  count(1) = SCALARSIZE
  count(2) = nmtx
  !XXX count(3) = 1
  count(3) = ncol
  count(4) = 1
  count(5) = 1 !mat
  count(6) = 1 !iq
  call h5screate_simple_f(rank, count, memspace, error)

  offset(:) = 0
  offset(3) = j - 1
  offset(5) = is - 1
  offset(6) = iq - 1
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

! XXX - this is bad because we double memory usage. Oh well.
  SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
  call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
!  write(6,*) 'About to die?', nmtx, data(1,1,:,1,1,1)
  eps(1:nmtx,1:ncol) = SCALARIFY2(data(1,1:nmtx,1:ncol,1,1,1),data(2,1:nmtx,1:ncol,1,1,1))
  SAFE_DEALLOCATE(data)

  call h5sclose_f(memspace, error)
  call h5sclose_f(dataspace, error)
  call h5dclose_f(data_id, error)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_matrix_col_hdf5)

end subroutine read_eps_matrix_col_hdf5

!===========================================================================================

subroutine read_eps_matrix_ser_hdf5(eps,nmtx,iq,is,name)
  integer, intent(in) :: iq
  integer, intent(in) :: is
  integer, intent(in) :: nmtx
  SCALAR, intent(out) :: eps(:,:) !< (nmtx,nmtx)
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: data_id       ! Property list identifier
  integer(HID_T) :: dataspace        ! Property list identifier
  integer(HID_T) :: memspace        ! Property list identifier
  integer :: error, rank
  integer(HSIZE_T) :: count(6), offset(6)

  real(DP), allocatable :: data(:,:,:,:,:,:)

  PUSH_SUB(read_eps_matrix_ser_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  call h5dopen_f(file_id, 'mats/matrix', data_id, error)
  call h5dget_space_f(data_id, dataspace, error)

  rank = 6
  count(1) = SCALARSIZE
  count(2) = nmtx
  count(3) = nmtx
  count(4) = 1
  count(5) = 1 !mat
  count(6) = 1 !iq
  call h5screate_simple_f(rank, count, memspace, error)

  offset(:) = 0
  offset(5) = is - 1
  offset(6) = iq - 1
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

! XXX - this is bad because we double memory usage. Oh well.
  SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
  call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
  eps(1:nmtx,1:nmtx) = SCALARIFY2(data(1,1:nmtx,1:nmtx,1,1,1),data(2,1:nmtx,1:nmtx,1,1,1))
  SAFE_DEALLOCATE(data)

  call h5sclose_f(memspace, error)
  call h5sclose_f(dataspace, error)
  call h5dclose_f(data_id, error)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_matrix_ser_hdf5)

end subroutine read_eps_matrix_ser_hdf5

!===========================================================================================

!> Reads the whole dielectric matrix in serial.
subroutine read_eps_matrix_ser_f_hdf5(eps,nmtx,nfreq,iq,is,name)
  complex(DPC), intent(out) :: eps(:,:,:) !< (nmtx,nmtx,nfreq)
  integer, intent(in) :: nmtx
  integer, intent(in) :: nfreq
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: data_id       ! Property list identifier
  integer(HID_T) :: dataspace        ! Property list identifier
  integer(HID_T) :: memspace        ! Property list identifier
  integer :: error, rank
  integer(HSIZE_T) :: count(6), offset(6)

  real(DP), allocatable :: data(:,:,:,:,:,:)

  PUSH_SUB(read_eps_matrix_ser_hdf5)

  call hdf5_open_file(trim(name), 'r', file_id)
  call h5dopen_f(file_id, 'mats/matrix', data_id, error)
  call h5dget_space_f(data_id, dataspace, error)

  rank = 6
  count(1) = SCALARSIZE
  count(2) = nmtx
  count(3) = nmtx
  count(4) = nfreq
  count(5) = 1 !mat
  count(6) = 1 !iq
  call h5screate_simple_f(rank, count, memspace, error)

  offset(:) = 0
  offset(5) = is - 1
  offset(6) = iq - 1
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

! XXX - this is bad because we double memory usage. Oh well.
  SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
  call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
  eps(1:nmtx,1:nmtx,1:nfreq) = CMPLX(data(1,1:nmtx,1:nmtx,1:nfreq,1,1),data(2,1:nmtx,1:nmtx,1:nfreq,1,1))
  SAFE_DEALLOCATE(data)

  call h5sclose_f(memspace, error)
  call h5sclose_f(dataspace, error)
  call h5dclose_f(data_id, error)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_matrix_ser_f_hdf5)

end subroutine read_eps_matrix_ser_f_hdf5

!===========================================================================================

subroutine read_eps_matrix_par_hdf5(eps, nb, rank, npes, nmtx, iq, is, name)
  SCALAR, intent(inout) :: eps(:,:) !< (neps,ngpown)
  integer, intent(in) :: nb !< block size
  integer, intent(in) :: rank !< processor rank for column distribution
  integer, intent(in) :: npes !< number of processors over which we distribute
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: data_id       ! Property list identifier
  integer(HID_T) :: plist_id       ! Property list identifier
  integer(HID_T) :: dataspace        ! Property list identifier
  integer(HID_T) :: memspace        ! Property list identifier
  integer :: error
  integer :: comm, info
  integer(HSIZE_T) :: count(6), offset(6), countm(6)
  real(DP), allocatable :: data(:,:,:,:,:,:)
  integer :: ngpown_max, igp, igp_loc

  PUSH_SUB(read_eps_matrix_par_hdf5)

  ngpown_max = NUMROC(nmtx, nb, 0, 0, npes)
  SAFE_ALLOCATE(data,(SCALARSIZE,nmtx,1,1,1,1))

  call hdf5_open_file(trim(name), 'r', file_id, parallel_io=.true.)

  do igp_loc = 1, ngpown_max

    igp = INDXL2G(igp_loc, nb, rank, 0, npes)
    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    call h5dget_space_f(data_id, dataspace, error)

! JRD: The commented code in this routine represents efforts to use a single HDF5 read call
! with an appropriate block and stride for each proc. In general, it currently appears to
! perform worse (though this may be related to size of matrix. So, it is commented until
! further investigation.

    countm(1) = SCALARSIZE
    countm(2) = nmtx
    countm(3) = 1
    countm(4) = 1
    countm(5) = 1
    countm(6) = 1
    call h5screate_simple_f(6, countm, memspace, error)

    count(1) = SCALARSIZE
    count(2) = nmtx
    count(3) = 1
    count(4) = 1
    count(5) = 1
    count(6) = 1
    ! Construct data and offset

    if (igp <= nmtx) then
    !if (ngpown .ne. 0 .and. my_igp .le. nmtx) then
      offset(1)=0
      offset(2)=0
      offset(3)=igp-1
      offset(4)=0
      offset(5)=is-1
      offset(6)=iq-1
      ! Select hyperslab
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

    else

      call h5sselect_none_f(memspace, error)
      call h5sselect_none_f(dataspace, error)

    endif

    ! Create property list for collective dataset read
#ifdef MPI
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    ! Collectively read the file
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, dataspace, &
                      xfer_prp = plist_id)
#else
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, dataspace)
#endif
    if (igp <= nmtx) then
      !    if (ngpown .ne. 0 .and. my_igp .le. nmtx) then
      !XXX PROBABLY NEED THREADED LOOP HERE
      eps(1:nmtx,igp_loc) = SCALARIFY2(data(1,1:nmtx,1,1,1,1),data(2,1:nmtx,1,1,1,1))
    endif
    call h5sclose_f(memspace, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(data_id, error)
  enddo

  call hdf5_close_file(file_id)
  SAFE_DEALLOCATE(data)

  POP_SUB(read_eps_matrix_par_hdf5)

end subroutine read_eps_matrix_par_hdf5

!==========================================================================================

subroutine read_eps_matrix_par_f_hdf5(retarded, nb, pool_comm, my_pool, npools, &
  nmtx, nFreq, iq, is, name, advanced)
  complex(DPC), intent(inout) :: retarded(:,:,:) !< (neps,ngpown,nFreq)
  integer, intent(in) :: nb !< block size
  integer, intent(in) :: pool_comm !< MPI comm for each pool
  integer, intent(in) :: my_pool !< my pool, starting from 0 (0=no pools)
  integer, intent(in) :: npools !< number of pools (1=no pools).
  integer, intent(in) :: nmtx
  integer, intent(in) :: nFreq
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name
  complex(DPC), optional, intent(inout) :: advanced(:,:,:) !< (nFreq,neps,ngpown)

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: data_id       ! Property list identifier
  integer(HID_T) :: plist_id       ! Property list identifier
  integer(HID_T) :: dataspace        ! Property list identifier
  integer(HID_T) :: memspace        ! Property list identifier
  integer :: error
  integer(HSIZE_T) :: count(6), offset(6)
  real(DP), allocatable :: data(:,:,:,:,:,:)

  integer :: pool_rank !< processor rank for column distribution
  integer :: npes_pool !< number of processors over which we distribute
  integer :: ngpown_max, igp, igp_loc, nmatrix_per_spin, nspin, buf_sz, version
  logical :: want_advanced, read_advanced
  !type(progress_info) :: prog_info !< a user-friendly progress report

  integer :: iw

  PUSH_SUB(read_eps_matrix_par_f_hdf5)

#ifdef MPI
  call MPI_Comm_rank(pool_comm, pool_rank, mpierr)
  call MPI_Comm_size(pool_comm, npes_pool, mpierr)
  ! FHJ: We need the following BCast for processors left over out of the pool
  call MPI_Bcast(npes_pool, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#else
  pool_rank = 0
  npes_pool = 1
#endif
  call hdf5_open_file(trim(name), 'r', file_id, parallel_io=.true.)

  want_advanced = present(advanced)
  ! FHJ: the default is never to read the advanced matrix, unless the file
  ! version is <3 (on which case we didn`t store the Coulomb interaction)
  call hdf5_read_int(file_id, 'eps_header/versionnumber', version)
  call hdf5_read_int(file_id, 'eps_header/params/nmatrix', nmatrix_per_spin)
  call hdf5_read_int(file_id, 'mf_header/kpoints/nspin', nspin)
  nmatrix_per_spin = nmatrix_per_spin / nspin
  read_advanced = .false.
  buf_sz = 1
  if (version<3) then
    call hdf5_read_logical(file_id, 'eps_header/params/has_advanced', read_advanced)
    if (SCALARSIZE==2 .and. .not.read_advanced) &
      call die('Inconsistent epsmat file: version<3, but no advanced matrix', only_root_writes=.true.)
  endif
  if (read_advanced) buf_sz = 2

  ngpown_max = NUMROC(nmtx, nb, 0, 0, npes_pool)
  SAFE_ALLOCATE(data,(2,nmtx,1,nFreq,buf_sz,1))

  call logit('Reading HDF5 file')
  !call progress_init(prog_info, 'reading '//trim(name), &
  !  'distributed column', ngpown_max)
! JRD XXX, we should do this as a single read instead of loop!!!
  do igp_loc = 1, ngpown_max
    !call progress_step(prog_info)
    if (my_pool>=0) then
      igp = INDXL2G(igp_loc, nb, pool_rank, 0, npes_pool)
    else
      igp = nmtx + 1
    endif
    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    call h5dget_space_f(data_id,dataspace,error)

    count(:) = (/2, nmtx, 1, nFreq, buf_sz, 1/)
    call h5screate_simple_f(6, count, memspace, error)

    ! Construct data and offset
    if (igp <= nmtx) then
      offset = (/0, 0, igp-1, 0, nmatrix_per_spin*(is-1), iq-1/)
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
    else
      call H5sselect_none_f(memspace,error);
      call H5sselect_none_f(dataspace,error);
    endif

! Create property list for collective dataset read
! Read is serial for now

#ifdef MPI
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, &
                      xfer_prp = plist_id)
#else
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
#endif

    if (igp <= nmtx) then
      retarded(1:nmtx,igp_loc,1:nFreq) = CMPLX(data(1,1:nmtx,1,1:nFreq,1,1),data(2,1:nmtx,1,1:nFreq,1,1))
      if (want_advanced .and. read_advanced) then
        advanced(1:nmtx,igp_loc,1:nFreq) = CMPLX(data(1,1:nmtx,1,1:nFreq,2,1),data(2,1:nmtx,1,1:nFreq,2,1))
      endif
    endif

    call h5sclose_f(memspace, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(data_id, error)
  enddo

  !call progress_free(prog_info)
  SAFE_DEALLOCATE(data)

  if (want_advanced .and. .not.read_advanced) call get_advanced_from_retarded()

  call hdf5_close_file(file_id)

  POP_SUB(read_eps_matrix_par_f_hdf5)

contains

  !> FHJ: This routine does the magic of generating advanced epsilon/chi out of
  !! the retarded epsilon/chi.
  !! We use the fact that W_A = (W_R)^H to write
  !! epsinv_A = (epsinv_R * v)^H * (1/v)
  subroutine get_advanced_from_retarded()
    real(DP) :: vcoul(nmtx)
    integer :: ig, igp, igp_loc, matrix_type
#ifdef MPI
    complex(DPC) :: send_buf(ngpown_max,ngpown_max,nFreq,npes_pool), recv_buf(ngpown_max,ngpown_max,nFreq,npes_pool)
    integer :: req_send(npes_pool), req_recv(npes_pool)
    integer :: ig_loc, ngown_ipe, ngpown, ipe
#endif

    PUSH_SUB(read_eps_matrix_par_f_hdf5.get_advanced_from_retarded)

    call logit('Getting advanced matrix from retarded')
    call hdf5_read_int(file_id, 'eps_header/params/matrix_type', matrix_type)
    vcoul(:) = 1d0
    ! FHJ: There is no Coulomb interaction in a chimat file.
    if (matrix_type<2) then
      call hdf5_read_double_hyperslab(file_id, 'eps_header/gspace/vcoul', &
        (/nmtx,1/), (/0,iq-1/), vcoul)
    endif

    if (my_pool<0) then
      POP_SUB(read_eps_matrix_par_f_hdf5.get_advanced_from_retarded)
      return
    endif

#ifdef MPI
    ngpown = NUMROC(nmtx, nb, pool_rank, 0, npes_pool)

    call logit('Preparing send buffers')
    ! FHJ: Prepare the send buffer. The send buffer will retarded^H, but with the
    ! columns of retarded^H (rows of retarded) reordered from ig=1..nmtx to irow=(ig_loc,ipe).
    ! We do this to distribute all the matrix elements with a single MPI_Scatter.
    send_buf = (0d0, 0d0)
    do ipe = 0, npes_pool-1
      ngown_ipe = NUMROC(nmtx, nb, ipe, 0, npes_pool)
      ! JRD XXX Loops not ordered well. FHJ this well be killer on systems with lots of G vecs
      ! FHJ: We are doing two matrix transpositions. We will *never* get good locality.
      ! The indices are optimized here for MPI communication. I actually timed this
      ! and saw that, with the opposite index order, MPI communication kills.
      ! For large ng, we are still paying the price of a large matrix transposition,
      ! which is still cheaper than lots of unordered MPI_BCasts`s or MPI_Gather`s.
      do ig_loc = 1, ngown_ipe
        ig = INDXL2G(ig_loc, nb, ipe, 0, npes_pool)
        send_buf(1:ngpown,ig_loc,1:nFreq,ipe+1) = CONJG(retarded(ig,1:ngpown,1:nFreq))
      enddo
    enddo

    call logit('Setting up non-blocking communication')
    call timing%start(timing%epscopy_comm)
    ! FHJ: Each processor scatters its retarded^H matrix, stored in send_buf, to all
    ! all other processors. Each PE gets nFreq*ngpown_max**2 complex numbers,
    ! stored as a (ngpown_max,ngpown_max) matrix of vectors of length nFreq.
    ! FHJ: Open non-blocking receiving buffers
    do ipe = 0, npes_pool-1
      call MPI_Irecv(recv_buf(1,1,1,ipe+1), nFreq*ngpown_max**2, MPI_COMPLEX_DPC, ipe, &
        pool_rank*npes_pool + ipe, pool_comm, req_recv(ipe+1), mpierr)
    enddo
    ! FHJ: Send data in non-blocking way
    do ipe = 0, npes_pool-1
      call MPI_Isend(send_buf(1,1,1,ipe+1), nFreq*ngpown_max**2, MPI_COMPLEX_DPC, ipe, &
        ipe*npes_pool + pool_rank, pool_comm, req_send(ipe+1), mpierr)
    enddo

    call logit('Receiving data')
    do while (.true.)
      ! FHJ: Work on whichever buffer arrived first.
      call MPI_Waitany(npes_pool, req_recv, ipe, MPI_STATUS_IGNORE, mpierr)
      if (ipe==MPI_UNDEFINED) exit
      ipe = ipe - 1
      ! FHJ: Processor "ipe" was kind enough to reorder its rows of retarded such in
      ! a way that it corresponds to all columns of advanced that we own, and in the
      ! correct order. However, each "ipe" originally owned a subset of all
      ! columns of retarded, so here we received a subset of the rows of advanced.
      ngown_ipe = NUMROC(nmtx, nb, ipe, 0, npes_pool)
      ! JRD: XXX this loop is horrifying, FHJ - what were you thinking??
      ! FHJ: Again, we are doing a **two matrix transpositions in parallel**.
      ! We will *never* get good locality.
      ! Seriously, before complaining so much about this line, did you see
      ! how efficient the MPI communication is -- which has a much smaller bandwitdh
      ! than this in-memory transposition?! And how we send a small number of
      ! MPI messages?! And how much all this is *much faster* than reading
      ! data from disk twice, which has an even smaller bandwidth?!
      ! There is no free lunch, you either pay the price of MPI communication
      ! + disk I/O or memory transposition.
      ! Memory transposition is still much cheaper, and I timed it.
      ! And BTW: this code is faster than BLACS for the matrix transposition.
      do iw = 1, nFreq
        do igp_loc = 1, ngpown
          do ig_loc = 1, ngown_ipe
            ig = INDXL2G(ig_loc, nb, ipe, 0, npes_pool)
            advanced(ig,igp_loc,iw) = recv_buf(ig_loc,igp_loc,iw,ipe+1)
          enddo
        enddo
      enddo
    enddo

    call logit('Releasing send buffer')
    call MPI_Waitall(npes_pool, req_send, MPI_STATUSES_IGNORE, mpierr)
    call timing%stop(timing%epscopy_comm)

    call logit('Fixing asymmetry in the Coulomb interaction')
    ! FHJ: Finally, fix the asymmetry in the Coulomb interaction.
    ! JRD XXX Should OMP here
    ! FHJ: ok, but note that nFreq might be small.
    do iw = 1, nFreq
      do igp_loc = 1, ngpown
        igp = INDXL2G(igp_loc, nb, pool_rank, 0, npes_pool)
        do ig = 1, nmtx
          advanced(ig,igp_loc,iw) = advanced(ig,igp_loc,iw) * (vcoul(ig)/vcoul(igp))
        enddo
      enddo
    enddo
#else
    ! FHJ: serial is just slightly simpler
    do igp = 1, nmtx
      do ig = 1, nmtx
        advanced(ig,igp,1:nFreq) = CONJG(retarded(igp,ig,1:nFreq)) * (vcoul(ig)/vcoul(igp))
      enddo
    enddo
#endif

    POP_SUB(read_eps_matrix_par_f_hdf5.get_advanced_from_retarded)

  end subroutine get_advanced_from_retarded

end subroutine read_eps_matrix_par_f_hdf5

!==========================================================================================
! Routine used in SUBSPACE method, can be used for other application also
#if defined MPI && defined USESCALAPACK

subroutine read_eps_matrix_par_hdf5_general(scal, matrix, nmtx, nmtx_col, iq, is, name, set_name, &
                                            comm)
  type(scalapack), intent(in) :: scal
  complex(DPC), intent(inout) :: matrix(:,:) !< (local_row,local_col)
  integer, intent(in) :: nmtx, nmtx_col  ! global row and col of the matrix to read in
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: set_name ! where to read from
  integer, intent(in) :: comm ! MPI communicator

  integer :: ii, jj, error, rank
  real(DP), allocatable :: data(:,:,:,:,:,:)

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
  integer(HID_T) :: plist_id      ! Property list identifier for parallel IO

  integer(HSIZE_T) :: count(6), countm(6), offset(6), offsetm(6), stride(6), block(6)
  integer(HSIZE_T) :: countr(6), offsetr(6), strider(6), blockr(6)

  integer :: info, rowremainder, colremainder

  PUSH_SUB(read_eps_matrix_par_hdf5_general)

  ! synchronize subgroup
  call MPI_barrier(comm, mpierr)

  rank=6
  countm(1) = 2
  countm(2) = scal%npr
  countm(3) = scal%npc
  countm(4) = 1
  countm(5) = 1
  countm(6) = 1

  offsetm(:) = 0

  count(1) = 1
  count(2) = scal%npr/scal%nbl
  count(3) = scal%npc/scal%nbl
  count(4) = 1
  count(5) = 1
  count(6) = 1

  block(1) = 2
  block(2) = scal%nbl
  block(3) = scal%nbl
  block(4) = 1
  block(5) = 1
  block(6) = 1

  offset(1) = 0
  offset(2) = scal%myprow*scal%nbl
  offset(3) = scal%mypcol*scal%nbl
  offset(4) = 0
  offset(5) = is-1
  offset(6) = iq-1

  stride(1) = 1
  stride(2) = scal%nprow*scal%nbl
  stride(3) = scal%npcol*scal%nbl
  stride(4) = 1
  stride(5) = 1
  stride(6) = 1

  SAFE_ALLOCATE(data,(countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))
  data = 0.0d0

  call hdf5_open_file(trim(name), 'r', file_id, parallel_io=.true., comm=comm)

  call h5dopen_f(file_id, trim(set_name), dset_id, error)
  call h5dget_space_f(dset_id, filespace, error)

  call h5screate_simple_f(rank, countm, memspace, error)
  !XXX call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)

  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)

  ! Bottom Rows
  rowremainder = mod(scal%npr,scal%nbl)
  if (rowremainder .ne. 0) then
    offsetr=offset
    countr=count
    blockr=block
    strider=stride
    offsetr(2)=nmtx-rowremainder
    countr(2)=rowremainder
    blockr(2)=1
    strider(2)=1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
    !write(6,*) peinf%inode, "I have the bottom row", rowremainder, scal%npc
  endif

  ! Right Columns
  colremainder = mod(scal%npc,scal%nbl)
  if (colremainder .ne. 0) then
    offsetr=offset
    countr=count
    blockr=block
    strider=stride
    offsetr(3)=nmtx_col - colremainder
    countr(3)=colremainder
    blockr(3)=1
    strider(3)=1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
    !write(6,*) peinf%inode, "I have the right column", colremainder, scal%npr
    ! Bottom Corner of Matrix
    if (rowremainder .ne. 0) then
      offsetr=offset
      countr=count
      blockr=block
      strider=stride
      offsetr(2)=nmtx-rowremainder
      countr(2)=rowremainder
      blockr(2)=1
      strider(2)=1
      offsetr(3)=nmtx_col - colremainder
      countr(3)=colremainder
      blockr(3)=1
      strider(3)=1
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
      !write(6,*) peinf%inode, "I have bottom both"
    endif
  endif

  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, &
                  xfer_prp = plist_id)

  do jj = 1, scal%npc
    do ii = 1, scal%npr
      matrix(ii,jj) = CMPLX(data(1,ii,jj,1,1,1),data(2,ii,jj,1,1,1))
    enddo
  enddo

  SAFE_DEALLOCATE(data)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(memspace, error)
  call h5sclose_f(filespace, error)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_matrix_par_hdf5_general)

end subroutine read_eps_matrix_par_hdf5_general

subroutine read_eps_matrix_par_hdf5_general_f(scal, matrix, nmtx, nmtx_col, Nfreq, iq, is, name, set_name, &
                                              comm)
  type(scalapack), intent(in) :: scal
  complex(DPC), intent(inout) :: matrix(:,:,:) !< (local_row,local_col,Nfreq)
  integer, intent(in) :: nmtx, nmtx_col  ! global row and col of the matrix to read in
  integer, intent(in) :: nfreq
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: set_name ! where to read from
  integer, intent(in) :: comm ! MPI communicator

  integer :: ii, jj, kk, error, rank
  real(DP), allocatable :: data(:,:,:,:,:,:)

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
  integer(HID_T) :: plist_id      ! Property list identifier for parallel IO

  integer(HSIZE_T) :: count(6), countm(6), offset(6), offsetm(6), stride(6), block(6)
  integer(HSIZE_T) :: countr(6), offsetr(6), strider(6), blockr(6)

  integer :: info, rowremainder, colremainder

  PUSH_SUB(read_eps_matrix_par_hdf5_general_f)

  ! even if possible lets keep it simple for now and assume nmtx = nmtx_col
  ! meaning that we consider only square subspace matrices

  ! MPI_INFO_NULL can be used for info in place of of any actual information
  ! (also when using subgroups, hopefully)
  info = MPI_INFO_NULL

  ! synchronize subgroup
  call MPI_barrier(comm,mpierr)

  rank=6
  countm(1) = 2
  countm(2) = scal%npr
  countm(3) = scal%npc
  countm(4) = Nfreq
  countm(5) = 1
  countm(6) = 1

  offsetm(:) = 0

  count(1) = 1
  count(2) = scal%npr/scal%nbl
  count(3) = scal%npc/scal%nbl
  count(4) = 1
  count(5) = 1
  count(6) = 1

  block(1) = 2
  block(2) = scal%nbl
  block(3) = scal%nbl
  block(4) = Nfreq
  block(5) = 1
  block(6) = 1

  offset(1) = 0
  offset(2) = scal%myprow*scal%nbl
  offset(3) = scal%mypcol*scal%nbl
  offset(4) = 0
  offset(5) = is-1
  offset(6) = iq-1

  stride(1) = 1
  stride(2) = scal%nprow*scal%nbl
  stride(3) = scal%npcol*scal%nbl
  stride(4) = 1
  stride(5) = 1
  stride(6) = 1

  SAFE_ALLOCATE(data,(countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))
  data = 0.0d0

  call hdf5_open_file(trim(name), 'r', file_id, parallel_io=.true., comm=comm)
  call h5dopen_f(file_id, trim(set_name), dset_id, error)
  call h5dget_space_f(dset_id,filespace,error)

  call h5screate_simple_f(rank, countm, memspace, error)
  !XXX call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)

  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)

  ! Bottom Rows
  rowremainder = mod(scal%npr,scal%nbl)
  if (rowremainder .ne. 0) then
    offsetr=offset
    countr=count
    blockr=block
    strider=stride
    offsetr(2)=nmtx-rowremainder
    countr(2)=rowremainder
    blockr(2)=1
    strider(2)=1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
    !write(6,*) peinf%inode, "I have the bottom row", rowremainder, scal%npc
  endif

  ! Right Columns
  colremainder = mod(scal%npc,scal%nbl)
  if (colremainder .ne. 0) then
    offsetr=offset
    countr=count
    blockr=block
    strider=stride
    offsetr(3)=nmtx-colremainder
    countr(3)=colremainder
    blockr(3)=1
    strider(3)=1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
    !write(6,*) peinf%inode, "I have the right column", colremainder, scal%npr
    ! Bottom Corner of Matrix
    if (rowremainder .ne. 0) then
      offsetr=offset
      countr=count
      blockr=block
      strider=stride
      offsetr(2)=nmtx-rowremainder
      countr(2)=rowremainder
      blockr(2)=1
      strider(2)=1
      offsetr(3)=nmtx-colremainder
      countr(3)=colremainder
      blockr(3)=1
      strider(3)=1
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
      !write(6,*) peinf%inode, "I have bottom both"
    endif
  endif

  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, &
                  xfer_prp = plist_id)

  do kk =1, Nfreq
    do jj = 1, scal%npc
      do ii = 1, scal%npr
        matrix(ii,jj,kk) = CMPLX(data(1,ii,jj,kk,1,1),data(2,ii,jj,kk,1,1))
      enddo
    enddo
  end do

  SAFE_DEALLOCATE(data)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(memspace, error)
  call h5sclose_f(filespace, error)
  call hdf5_close_file(file_id)

  POP_SUB(read_eps_matrix_par_hdf5_general_f)

end subroutine read_eps_matrix_par_hdf5_general_f

#endif
! END SUBSPACE specific routine
!==========================================================================================

#endif
end module epsread_hdf5_m
