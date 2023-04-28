!>=========================================================================
!!
!!  Module:
!!
!!  epswrite_hdf5_m     Originally by JRD     Last Modified 12/2014 (FHJ)
!!
!!    Routines to write header info for epsmat files in HDF5 format.
!!
!!=========================================================================

#include "f_defs.h"

module epswrite_hdf5_m
#ifdef HDF5
  use hdf5
  use global_m
  use hdf5_io_m
  use wfn_io_hdf5_m
  implicit none

  private

  public :: &
    eps_hdf5_setup, &
    set_subspace_neigen_q, &
    set_qpt_done, &
    is_qpt_done

contains


subroutine eps_hdf5_setup(kp, gvec, syms, crys, pol, qgrid, nq, qpts, nmtx, &
  nmtx_max, nband, name, restart)
  type(kpoints), intent(in) :: kp
  type(gspace), intent(in) :: gvec
  type(symmetry), intent(in) :: syms
  type(crystal), intent(in) :: crys
  type(polarizability), intent(in) :: pol
  integer, intent(in) :: qgrid(3)
  integer, intent(in) :: nq
  real(DP), intent(in) :: qpts(:,:) !< (3,nq)
  integer, intent(in) :: nmtx(:) !< (nq)
  integer, intent(in) :: nmtx_max
  integer, intent(in) :: nband
  character(len=*), intent(in) :: name
  logical, intent(inout), optional :: restart

  integer(HID_T) :: file_id
  integer(HID_T) :: dspace
  integer(HID_T) :: dset_id
  integer(HSIZE_T) :: dims(6)

  integer :: ii, intdata(1)
  integer :: neigen_of_q(nq)
  logical :: qpt_done(nq)
  real(DP), allocatable :: tmp(:,:)
  real(DP) :: freqs_tmp(2,pol%nfreq)
  real(DP) :: realdata(1)
  logical :: restart_, file_exists, file_ok
  character(len=3) :: sheader='WFN'

  PUSH_SUB(eps_hdf5_setup)

  restart_=.false.
  if (present(restart)) restart_ = restart

  ! FHJ: try to restart the calculation, if possible and desired.
  ! We ignore the restart flags if the file doesn`t exist. However, the code
  ! dies if the file exists and is incompatible or looks corrupted.
  if (restart_) then
    call try_restart()
    if (file_ok) then
      POP_SUB(eps_hdf5_setup)
      return
    end if
  endif

  ! FHJ: Set up file: write MF header and create groups
  write(6,'(1x,2a)') "Initializing ", trim(name)
  call setup_hdf5_mf_file(trim(name))
  call write_hdf5_header_type(trim(name), sheader, SCALARSIZE, kp, gvec, syms, crys)
  call write_hdf5_gvectors(trim(name), gvec%ng, gvec%components)
  call hdf5_open_file(trim(name), 'rw', file_id)
  call hdf5_create_group(file_id, 'eps_header')
  call hdf5_create_group(file_id, 'eps_header/params')
  call hdf5_create_group(file_id, 'eps_header/qpoints')
  call hdf5_create_group(file_id, 'eps_header/freqs')
  call hdf5_create_group(file_id, 'eps_header/gspace')
  call hdf5_create_group(file_id, 'mats')

  if( pol%subspace ) then
    call hdf5_create_group(file_id, 'eps_header/subspace')
  endif

  call hdf5_write_int(file_id, 'eps_header/versionnumber', VER_EPS_HDF5)
  call hdf5_write_int(file_id, 'eps_header/flavor', SCALARSIZE)

  ! FHJ: General datasets
  call hdf5_write_int(file_id, 'eps_header/params/matrix_type', pol%matrix_type)
  call hdf5_write_logical(file_id, 'eps_header/params/has_advanced', pol%has_advanced)
  call hdf5_write_int(file_id, 'eps_header/params/nmatrix', pol%nmatrix)
  call hdf5_write_int(file_id, 'eps_header/params/matrix_flavor', pol%matrix_flavor)
  call hdf5_write_int(file_id, 'eps_header/params/icutv', pol%icutv)
  call hdf5_write_double(file_id, 'eps_header/params/ecuts', pol%ecuts)
  call hdf5_write_int(file_id, 'eps_header/params/nband', nband)
  call hdf5_write_double(file_id, 'eps_header/params/efermi', pol%efermi/ryd)
  call hdf5_write_int(file_id, 'eps_header/params/intraband_flag', pol%intraband_flag)
  call hdf5_write_double(file_id, 'eps_header/params/intraband_overlap_min', pol%intraband_overlap_min)
  call hdf5_write_logical(file_id, 'eps_header/params/subsample', pol%subsample)
  call hdf5_write_logical(file_id, 'eps_header/params/subspace', pol%subspace)

  ! FHJ: Q-points-related datasets
  qpt_done(:) = .false.
  call hdf5_write_int(file_id, 'eps_header/qpoints/nq', nq)
  call hdf5_write_double_array(file_id, 'eps_header/qpoints/qpts', (/3,nq/), qpts)
  call hdf5_write_int_array(file_id, 'eps_header/qpoints/qgrid', (/3/), qgrid)
  call hdf5_write_logical_array(file_id, 'eps_header/qpoints/qpt_done', (/nq/), qpt_done)

  ! FHJ: Frequency-related datasets
  call hdf5_write_int(file_id, 'eps_header/freqs/freq_dep', pol%freq_dep)
  call hdf5_write_int(file_id, 'eps_header/freqs/nfreq', pol%nfreq)
  call hdf5_write_int(file_id, 'eps_header/freqs/nfreq_imag', pol%nfreq_imag)
  do ii=1,pol%nfreq !FHJ: TODO - have an unique complex freq grid in the future!
    freqs_tmp(1,ii) = pol%dFreqGrid(ii) + dble(pol%dFreqBrd(ii))
    freqs_tmp(2,ii) = IMAG(pol%dFreqBrd(ii))
  enddo
  call hdf5_write_double_array(file_id, 'eps_header/freqs/freqs', (/2, pol%nfreq/), freqs_tmp)

  ! FHJ: G-vectors-related datasets
  call hdf5_write_int_array(file_id, 'eps_header/gspace/nmtx', (/nq/), nmtx)
  call hdf5_write_int(file_id, 'eps_header/gspace/nmtx_max',  nmtx_max)
  call hdf5_create_dset(file_id, 'eps_header/gspace/ekin', H5T_NATIVE_DOUBLE, (/gvec%ng,nq/))
  ! FHJ: Epsmat version 3+ includes the bare Coulomb interaction as well.
  ! We use this to get the advanced epsilon from the retarded matrix.
  call hdf5_create_dset(file_id, 'eps_header/gspace/vcoul', H5T_NATIVE_DOUBLE, (/nmtx_max,nq/))
  call hdf5_create_dset(file_id, 'eps_header/gspace/gind_eps2rho', H5T_NATIVE_INTEGER, (/gvec%ng,nq/))
  call hdf5_create_dset(file_id, 'eps_header/gspace/gind_rho2eps', H5T_NATIVE_INTEGER, (/gvec%ng,nq/))

  !DVF: subspace approximation-related datasets
  if( pol%subspace ) then
    call hdf5_write_logical(file_id, 'eps_header/subspace/keep_full_eps_static', pol%keep_full_eps_static)
    call hdf5_write_logical(file_id, 'eps_header/subspace/matrix_in_subspace_basis', pol%matrix_in_subspace_basis)
    call hdf5_write_double(file_id, 'eps_header/subspace/eps_eigenvalue_cutoff', pol%chi_eigenvalue_cutoff)
    call hdf5_write_int(file_id, 'eps_header/subspace/neig_max',  pol%neig_sub_input)
    ! use same strategy as qpt_done, since the number of eigenvector will never exceed neig_sub_input, we
    ! initially set all of them to neig_sub_input and we update the actual value when calculated
    neigen_of_q(:) = pol%neig_sub_input
    call hdf5_write_int_array(file_id, 'eps_header/subspace/neig', (/nq/), neigen_of_q)
  endif

  ! FHJ: Matrix-elements-related datasets
  call hdf5_create_dset(file_id, 'mats/matrix', H5T_NATIVE_DOUBLE, &
    (/pol%matrix_flavor,nmtx_max,nmtx_max,pol%nfreq,pol%nmatrix,nq/))
  call hdf5_create_dset(file_id, 'mats/matrix-diagonal', H5T_NATIVE_DOUBLE, &
    (/pol%matrix_flavor,nmtx_max,nq/))

  if( pol%subspace .and. pol%matrix_in_subspace_basis ) then
    !
    ! here pol%nmatrix should always be ONE, this may change if we want to write
    ! the subspace chi, which can be spin polarized, which not sure if make any sense.
    !
    call hdf5_create_dset(file_id, 'mats/matrix_fulleps0', H5T_NATIVE_DOUBLE, &
       (/pol%matrix_flavor,nmtx_max,nmtx_max,1,pol%nmatrix,nq/))
    !XXX
    !XXX actual (max) size
    call hdf5_create_dset(file_id, 'mats/matrix_eigenvec', H5T_NATIVE_DOUBLE, &
        (/pol%matrix_flavor,nmtx_max,pol%neig_sub_input,1,pol%nmatrix,nq/))
    !XXX
    !XXX full matrix
    !XXX call hdf5_create_dset(file_id, 'mats/matrix_eigenvec', H5T_NATIVE_DOUBLE, &
    !XXX     (/pol%matrix_flavor,nmtx_max,nmtx_max,1,pol%nmatrix,nq/))
    !XXX
    call hdf5_create_dset(file_id, 'mats/matrix_subspace', H5T_NATIVE_DOUBLE, &
       (/pol%matrix_flavor,pol%neig_sub_input,pol%neig_sub_input,pol%nfreq,pol%nmatrix,nq/))
  end if

  call hdf5_close_file(file_id)

  POP_SUB(eps_hdf5_setup)

contains

  subroutine try_restart()
    integer :: nspin_old, nq_old
    real(DP) :: qpts_old(3,nq)
    integer :: error, freq_dep_old, nfreq_old
    integer :: ng_old, gvecs_old(3,gvec%ng), nmtx_old(nq), matrix_flavor_old, nmatrix_old
    integer :: neig_max_old
    logical :: subspace_old, matrix_in_subspace_basis_old, keep_full_eps_static_old
    logical :: subspace_error

    PUSH_SUB(eps_hdf5_setup.try_restart)

    subspace_error = .false.
    write(6,'(1x,2a)') "Trying to restart file ", trim(name)
    inquire(file=trim(name), exist=file_exists)
    if (.not.file_exists) then
      file_ok = .false.
      write(0,'(3a)') "WARNING: file ", trim(name), " doesn`t exist. We`ll start the calculation from scratch."
      restart_ = .false.
      if (present(restart)) restart = .false.

      POP_SUB(eps_hdf5_setup.try_restart)
      return
    endif

    call hdf5_open_file(trim(name), 'r', file_id, error=error)
    if (error/=0) then
      write(0,'(3a,i0,a)') "ERROR: ", trim(name), " is not a valid HDF5 file (error code: ", error," )."
      write(0,'(a)') 'Make sure the file is not corrupted'
      call die("file "//trim(name)//" looks corrupted", only_root_writes=.true.)
    endif

    ! FHJ: Consistency check.
    file_ok = hdf5_exists(file_id, 'eps_header/qpoints/qpt_done')
    ! We only execute the loop below once, but the loop structure is
    ! convenient because it allows us to use the `exit` statement.
    do while(file_ok)
      call hdf5_read_int(file_id, 'mf_header/kpoints/nspin', nspin_old, error)
      if (error/=0) exit
      call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq_old, error)
      if (error/=0) exit
      call hdf5_read_double_array(file_id, 'eps_header/qpoints/qpts', (/3,nq/), qpts_old, error)
      if (error/=0) exit
      call hdf5_read_int(file_id, 'eps_header/freqs/freq_dep', freq_dep_old, error)
      if (error/=0) exit
      call hdf5_read_int(file_id, 'eps_header/freqs/nfreq', nfreq_old, error)
      if (error/=0) exit
      call hdf5_read_int(file_id, 'mf_header/gspace/ng', ng_old, error)
      if (error/=0) exit
      call hdf5_read_int_array(file_id, 'mf_header/gspace/components', (/3,gvec%ng/), gvecs_old, error)
      if (error/=0) exit
      call hdf5_read_int_array(file_id, 'eps_header/gspace/nmtx', (/nq/), nmtx_old, error)
      if (error/=0) exit
      call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', matrix_flavor_old, error)
      if (error/=0) exit
      call hdf5_read_int(file_id, 'eps_header/params/nmatrix', nmatrix_old, error)
      if (error/=0) exit
      exit
    enddo
    file_ok = (error==0) .and. (nspin_old==kp%nspin .and. nq_old==nq .and. &
      all(dabs(qpts_old-qpts)<TOL_ZERO) .and. freq_dep_old==pol%freq_dep .and. &
      nfreq_old==pol%nfreq .and. ng_old==gvec%ng .and. all(gvecs_old==gvec%components) .and. &
      all(nmtx_old==nmtx) .and. matrix_flavor_old==pol%matrix_flavor .and. &
      nmatrix_old==pol%nmatrix)

    if (.not.file_ok) then
      write(0,*)
      write(0,'(3a)') "ERROR: restart file ", trim(name), " is incompatible with the current calculation."
      write(0,*) 'Values from file vs. calculation:'
      write(0,*) 'nspin:', nspin_old, kp%nspin
      write(0,*) 'nq:', nq_old, nq
      write(0,*) 'qpts (same?):', all(dabs(qpts_old-qpts)<TOL_ZERO)
      write(0,*) 'freq_dep:', freq_dep_old, pol%freq_dep
      write(0,*) 'nfreq:', nfreq_old, pol%nfreq
      write(0,*) 'ng:', ng_old, gvec%ng
      write(0,*) 'gvecs (same?):', all(gvecs_old==gvec%components)
      write(0,*) 'nmtx (same?):', all(nmtx_old==nmtx)
      write(0,*) 'matrix_flavor:', matrix_flavor_old, pol%matrix_flavor
      write(0,*) 'nmatrix:', nmatrix_old, pol%nmatrix
      write(0,*)
      write(0,*) 'NOTE: you should only trust the first pair of values that disagree.'
      write(0,*)
      call die("restart file "//trim(name)//" is incompatible with the current calculation.", &
        only_root_writes=.true.)
    endif

    ! Additional check in case of subspace
    ! Check the existence of the parameter (this ensure that the routine will
    ! work also for old .h5 files)
    subspace_old = .false.
    if (hdf5_exists(file_id, 'eps_header/params/subspace')) then
      call hdf5_read_logical(file_id, 'eps_header/params/subspace', subspace_old, error)
    endif

    ! eps_header/params/subspace exists, check new and old flags
    file_ok = (error==0) .and. (pol%subspace .eqv. subspace_old)
    if (.not.file_ok) then
      write(0,'(3a)') "ERROR: restart file ", trim(name), " is incompatible with the current calculation."
      if (error/=0) then
        write(6,'(1x,a)') "Dataset eps_header/params/subspace is corrupted."
      elseif (pol%subspace) then
        write(6,'(1x,a)') "Trying to restart a subspace calc from a non-subspace calc"
      else
        write(6,'(1x,a)') "Trying to restart a non-subspace calc from a subspace calc"
      endif
      write(0,*)
      call die("restart file "//trim(name)//" is incompatible with the current calculation.", &
        only_root_writes=.true.)
    endif

    if (subspace_old) then
      ! if the flags are equivalent and we are doing a subspace calculation
      ! check the remaining parameters
      do while (.true.)
        call hdf5_read_int(file_id, 'eps_header/subspace/neig_max', &
          neig_max_old, error)
         if (error/=0) exit
        call hdf5_read_logical(file_id, 'eps_header/subspace/matrix_in_subspace_basis', &
          matrix_in_subspace_basis_old, error)
         if (error/=0) exit
        call hdf5_read_logical(file_id, 'eps_header/subspace/keep_full_eps_static', &
          keep_full_eps_static_old, error)
        if (error/=0) exit
        exit
      enddo
      file_ok = (error==0) &
                .and. (neig_max_old == pol%neig_sub_input) &
                .and. (matrix_in_subspace_basis_old .eqv. pol%matrix_in_subspace_basis) &
                .and. (keep_full_eps_static_old .eqv. pol%keep_full_eps_static)
      if (.not.file_ok) then
        write(0,'(3a)') "ERROR: restart file ", trim(name), " is incompatible with the current calculation."
        write(0,*) "Problem is most likely due to incompatible .h5 files within subspace approximation"
        write(0,*) 'Values from file vs. calculation:'
        write(0,*) 'neig_max:', neig_max_old, pol%neig_sub_input
        write(0,*) 'matrix_in_subspace_basis:', matrix_in_subspace_basis_old, pol%matrix_in_subspace_basis
        write(0,*) 'keep_full_eps_static:', keep_full_eps_static_old, pol%keep_full_eps_static
        write(0,*)
        write(0,*) 'NOTE: you should only trust the first pair of values that disagree.'
        write(0,*)
        call die("restart file "//trim(name)//" is incompatible with the current calculation.", &
          only_root_writes=.true.)
      endif
    endif
    call hdf5_close_file(file_id)
    ! FHJ: File *seems* alright, we don`t have to initialize it
    write(6,'(1x,2a)') "Everything looks ok: restarting file ", trim(name)

    POP_SUB(eps_hdf5_setup.try_restart)
  end subroutine try_restart

end subroutine eps_hdf5_setup


subroutine setup_subspace_mats_hdf5(fname, pol,iq)
  character(len=*), intent(in) :: fname
  integer, intent(in) :: iq
  type(polarizability), intent(in) :: pol

  integer(HID_T) :: file_id
  character(len=20) :: subspace_name
! DVF: we create different subgroups in the mats group for the
! subspace matrices and basis vectors at each q-point. This is
! necessary because we don`t know the size of the subspace basis
! at each q-point before we actually do the calculation like we do
! for the reciprocal space basis. So, we need to allocate the matrices
! right after we have run scalapack/ELPA, and then create the needed
! subgroups and dsets. That is the purpose of this routine.

! The above description assumes we specify the eigenvector tolerance and
! not the number of eigenvectors. Specifying the tolerance is best practice
! since the eigenvector tolerance is directly related to the convergence of
! one`s calculation, and specifying the tolerance requires no prior knowledge
! of the size of the g-space, and hence is more automatic for users.

  PUSH_SUB(setup_subspace_mats_hdf5)

  call hdf5_open_file(trim(fname), 'rw', file_id)
  ! DVF: subspace matrix-elements-related datasets
  write(subspace_name,'(a14,i4)') "mats/subspace_",iq
  call hdf5_create_group(file_id, trim(subspace_name))
  call hdf5_create_dset(file_id, trim(subspace_name) // '/matrix_sub', H5T_NATIVE_DOUBLE, &
    (/pol%matrix_flavor,pol%neig_sub,pol%neig_sub,pol%nfreq,pol%nmatrix,1/))
  call hdf5_create_dset(file_id, trim(subspace_name) // '/basis_sub', H5T_NATIVE_DOUBLE, &
    (/pol%matrix_flavor,pol%nmtx,pol%neig_sub,1,pol%nmatrix,1/))
  call hdf5_close_file(file_id)

  POP_SUB(setup_subspace_mats_hdf5)

end subroutine setup_subspace_mats_hdf5

subroutine set_subspace_neigen_q(fname, iq, neigen)
  character(len=*), intent(in) :: fname
  integer, intent(in) :: iq, neigen

  integer(HID_T) :: file_id
  integer :: nq
  integer, allocatable :: neigen_of_q(:)

  PUSH_SUB(set_subspace_neigen_q)

  call hdf5_open_file(trim(fname), 'rw', file_id)
  call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq)
  SAFE_ALLOCATE(neigen_of_q, (nq))

  call hdf5_read_int_array(file_id, 'eps_header/subspace/neig', (/nq/), neigen_of_q)
  neigen_of_q(iq) = neigen
  call hdf5_write_int_array(file_id, 'eps_header/subspace/neig', (/nq/), neigen_of_q)
  call hdf5_close_file(file_id)
  SAFE_DEALLOCATE(neigen_of_q)

  POP_SUB(set_subspace_neigen_q)

end subroutine set_subspace_neigen_q

subroutine set_qpt_done(fname, iq)
  character(len=*), intent(in) :: fname
  integer, intent(in) :: iq

  integer(HID_T) :: file_id
  integer :: nq
  logical, allocatable :: qpt_done(:)

  PUSH_SUB(set_qpt_done)

  call hdf5_open_file(trim(fname), 'rw', file_id)
  call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq)
  SAFE_ALLOCATE(qpt_done, (nq))
  call hdf5_read_logical_array(file_id, 'eps_header/qpoints/qpt_done', (/nq/), qpt_done)
  qpt_done(iq) = .true.
  call hdf5_write_logical_array(file_id, 'eps_header/qpoints/qpt_done', (/nq/), qpt_done)
  call hdf5_close_file(file_id)
  SAFE_DEALLOCATE(qpt_done)

  POP_SUB(set_qpt_done)

end subroutine set_qpt_done


logical function is_qpt_done(fname, iq)
  character(len=*), intent(in) :: fname
  integer, intent(in) :: iq

  integer(HID_T) :: file_id
  integer :: nq
  logical, allocatable :: qpt_done(:)

  PUSH_SUB(is_qpt_done)

  call hdf5_open_file(trim(fname), 'r', file_id)
  call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq)
  SAFE_ALLOCATE(qpt_done, (nq))
  call hdf5_read_logical_array(file_id, 'eps_header/qpoints/qpt_done', (/nq/), qpt_done)
  is_qpt_done = qpt_done(iq)
  call hdf5_close_file(file_id)
  SAFE_DEALLOCATE(qpt_done)

  POP_SUB(is_qpt_done)

end function is_qpt_done

#endif
end module epswrite_hdf5_m
