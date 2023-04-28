!==============================================================================
!
! Module gkq_m
!
! Originally by GKA (2018)
!
! Objects holding electron-phonon coupling matrix elements.
!
!==============================================================================

#include "f_defs.h"

module gkq_m

#ifdef HDF5
use hdf5
use hdf5_io_m
#endif

use global_m

implicit none

public

! =========================================================================== !

!> electron-phonon coupling matrix elements
!! (electron state basis)
type gkq_t

  ! Dimensions
  integer :: ns       !< Number of spins
  integer :: nk       !< Number of k-points
  integer :: nband    !< Number of bands at k
  integer :: mband    !< Number of bands at k+q
  integer :: nq       !< Number of q-points
  integer :: nat      !< number of atoms
  integer :: ndir=3   !< number of spacial (reduced) directions
  integer :: nmode    !< nat * ndir
  integer :: nq_me    !< Number of q-points own by the local worker

  ! Parallel info
  logical :: is_master = .true.
  logical :: is_active = .true.
  integer :: iq_start_me = 1 !< index of the first q-point
  integer :: qpt_mpi_comm    !< qpt communicator
  integer :: qpt_mpi_info    !< qpt info

  !> masses, in atomic mass units
  !! amu(nat)
  real(DP), allocatable :: amu(:)

  !> k-points, in reduced coordinates
  !! kpts(3,nk)
  real(DP), allocatable :: kpts(:,:)

  !> q-points, in reduced coordinates
  !! qpts(3,nq_me)
  real(DP), allocatable :: qpts(:,:)

  !> Electron-phonon coupling matrix elemnts, in atom/direction basis [Hartree]
  !! gkq(mband,nband,nk,ns,ndir,nat,nq)
  SCALAR, allocatable :: g_at(:,:,:,:,:,:,:)

  !> Electron-phonon coupling matrix elemnts, in mode basis [Hartree]
  !! gkq(mband,nband,nk,ns,nmode,nq)
  SCALAR, allocatable :: g_nu(:,:,:,:,:,:)

  contains

  ! Core procedures
  procedure, pass :: init => init_gkq
  procedure, pass :: alloc => alloc_gkq
  procedure, pass :: free => free_gkq
  procedure, pass :: copy => copy_gkq
  procedure, pass :: compare => compare_gkq
  procedure, pass :: setup_paral => setup_paral_gkq

  ! Setting procedures
  procedure, pass :: set_amu
  procedure, pass :: set_qpts_local => set_qpts_local_gkq
  procedure, pass :: set_qpts_global => set_qpts_global_gkq
  procedure, pass :: set_kpts => set_kpts_gkq

#ifdef HDF5
  ! Writing procedures
  procedure, pass :: write_hdf5 => write_gkq_hdf5
  procedure, pass :: create_and_write_gkq_header_hdf5
  procedure, pass :: write_gkq_matrix_par_hdf5
  procedure, pass :: write_gkq_matrix_ser_hdf5
  procedure, pass :: write_gkq_global_arrays_hdf5
  procedure, pass :: write_gkq_local_arrays_par_hdf5
  procedure, pass :: write_gkq_local_arrays_ser_hdf5

  ! Reading procedures
  procedure, pass :: read_hdf5 => read_gkq_hdf5
  procedure, pass :: read_and_broadcast_gkq_header_hdf5
  procedure, pass :: read_gkq_header_hdf5
  procedure, pass :: broadcast_gkq_dimensions
  procedure, pass :: read_and_broadcast_gkq_data_hdf5
  procedure, pass :: read_and_broadcast_gkq_global_arrays_hdf5
  procedure, pass :: read_gkq_local_arrays_par_hdf5
  procedure, pass :: read_gkq_local_arrays_ser_hdf5
  procedure, pass :: read_gkq_matrix_par_hdf5
  procedure, pass :: read_gkq_matrix_ser_hdf5
  procedure, pass :: read_and_broadcast_amu
  procedure, pass :: read_and_broadcast_kpts
#endif

end type gkq_t

! =========================================================================== !
contains
! =========================================================================== !

!> Initialize dimensions of the arrays.
subroutine init_gkq(this, ns, nk, nband, mband, nq, nat, nq_me)
  class(gkq_t), intent(inout) :: this
  integer, intent(in) :: ns       !< Number of spins
  integer, intent(in) :: nk       !< Number of k-points
  integer, intent(in) :: nband    !< Number of bands at k
  integer, intent(in) :: mband    !< Number of bands at k+q
  integer, intent(in) :: nq       !< Number of q-points
  integer, intent(in) :: nat      !< number of atoms
  integer, intent(in),optional :: nq_me    !< Number of q-points

  PUSH_SUB(init_gkq)

  this%ns = ns
  this%nk = nk
  this%nband = nband
  this%mband = mband
  this%nq = nq
  this%nat = nat
  this%nmode = this%nat * this%ndir

  if (present(nq_me)) then
    this%nq_me = nq_me
  else
    this%nq_me = this%nq
  end if

  POP_SUB(init_gkq)

end subroutine init_gkq

! =========================================================================== !

!> Setup the parallel layout of the data.
!! At the moment, only parallelism over q-points is supported.
subroutine setup_paral_gkq(this, peinf, np_qpt)
  class(gkq_t), intent(inout) :: this
  type(peinfo), intent(in) :: peinf     !< info on parallelization
  integer, intent(in) :: np_qpt         !< number of q-point groups

  integer :: inode
  integer :: nq_me_max

  PUSH_SUB(setup_paral_gkq)

  inode = peinf%inode
  this%is_master = (inode .eq. 0)

  ! Maximum number of q-point per processor
  nq_me_max = this%nq / np_qpt + min(1, mod(this%nq, np_qpt))

  this%iq_start_me = 1 + inode * nq_me_max

  ! Compute local number of q-point
  if (this%iq_start_me + nq_me_max - 1 .le. this%nq) then
    this%nq_me = nq_me_max
  else if (this%iq_start_me .le. this%nq) then
    this%nq_me = this%nq - this%iq_start_me + 1
  else
    this%nq_me = 0
  end if

  this%is_active = (this%nq_me .ne. 0)

#ifdef MPI
  this%qpt_mpi_comm = MPI_COMM_WORLD
  this%qpt_mpi_info = MPI_INFO_NULL
#endif

  POP_SUB(setup_paral_gkq)

end subroutine setup_paral_gkq

! =========================================================================== !

!> Allocate memory.
!! It is assumed that setup_paral has been called prior to this routine.
subroutine alloc_gkq(this)
  class(gkq_t), intent(inout) :: this

  PUSH_SUB(alloc_gkq)

  if (.not. allocated(this%amu)) then
    SAFE_ALLOCATE(this%amu, (this%nat))
  end if

  if (.not. allocated(this%kpts)) then
    SAFE_ALLOCATE(this%kpts, (3, this%nk))
  end if

  if (.not. allocated(this%qpts)) then
    SAFE_ALLOCATE(this%qpts, (3, this%nq_me))
  end if

  !if (.not. allocated(this%g_at)) then
  !  SAFE_ALLOCATE(this%g_at, (this%mband, this%nband, this%nk, this%ns, this%ndir, this%nat, this%nq_me))
  !end if

  if (.not. allocated(this%g_nu)) then
    SAFE_ALLOCATE(this%g_nu, (this%mband, this%nband, this%nk, this%ns, this%nmode, this%nq_me))
  end if

  POP_SUB(alloc_gkq)

end subroutine alloc_gkq

! =========================================================================== !

!> Free memory.
subroutine free_gkq(this)

  class(gkq_t), intent(inout) :: this

  PUSH_SUB(free_gkq)

  SAFE_DEALLOCATE(this%amu)
  SAFE_DEALLOCATE(this%kpts)
  SAFE_DEALLOCATE(this%qpts)
  !SAFE_DEALLOCATE(this%g_at)
  SAFE_DEALLOCATE(this%g_nu)

  POP_SUB(free_gkq)

end subroutine free_gkq

! =========================================================================== !

!> Copy object onto another.
subroutine copy_gkq(this, other)
  class(gkq_t), intent(inout) :: this
  class(gkq_t), intent(out) :: other

  PUSH_SUB(copy_gkq)

  other%ns = this%ns
  other%nk = this%nk
  other%nband = this%nband
  other%mband = this%mband
  other%nq = this%nq
  other%nat = this%nat
  other%ndir = this%ndir
  other%nmode = this%nmode

  other%is_master = this%is_master
  other%is_active = this%is_active
  other%nq_me = this%nq_me
  other%iq_start_me  = this%iq_start_me
  other%qpt_mpi_comm  = this%qpt_mpi_comm
  other%qpt_mpi_info  = this%qpt_mpi_info

  call other%alloc()

  other%amu = this%amu
  other%kpts = this%kpts
  other%qpts = this%qpts
  !other%g_at = this%g_at
  other%g_nu = this%g_nu

  POP_SUB(copy_gkq)

end subroutine copy_gkq

! =========================================================================== !

!> Compare a gkq object against another, and report any discrepancy.
!! Mostly for testing purposes.
subroutine compare_gkq(this, other, ndiff, verbose)
  class(gkq_t), intent(in) :: this
  class(gkq_t), intent(in) :: other
  integer, intent(out) :: ndiff
  logical, intent(in), optional :: verbose

  logical :: verbose_ = .True.
  real(DP), parameter :: tol=1.0d-5
  real(DP) :: checksum
  integer :: unt=6

  if (present(verbose)) verbose_ = verbose

  PUSH_SUB(compare_gkq)

  ndiff = 0
  checksum = abs(sum(this%amu) - sum(other%amu))
  if (checksum .gt. tol) then
    if (verbose_) write(unt,*) 'Error in checksum of: amu. ','delta=',checksum
    ndiff = ndiff + 1
  end if
  checksum = abs(sum(this%kpts) - sum(other%kpts))
  if (checksum .gt. tol) then
    if (verbose_) write(unt,*) 'Error in checksum of: kpts. ','delta=',checksum
    ndiff = ndiff + 1
  end if
  checksum = abs(sum(this%qpts) - sum(other%qpts))
  if (checksum .gt. tol) then
    if (verbose_) write(unt,*) 'Error in checksum of: qpts. ','delta=',checksum
    ndiff = ndiff + 1
  end if
  !checksum = abs(sum(this%g_at) - sum(other%g_at))
  checksum = abs(sum(this%g_nu) - sum(other%g_nu))
  if (checksum .gt. tol) then
    if (verbose_) write(unt,*) 'Error in checksum of: g_nu. ','delta=',checksum
    ndiff = ndiff + 1
  end if

  POP_SUB(compare_gkq)

end subroutine compare_gkq

! =========================================================================== !

!> Set value for the masses.
subroutine set_amu(this, amu)
  class(gkq_t), intent(inout) :: this
  real(DP), intent(in) :: amu(this%nat)

  PUSH_SUB(set_amu)

  this%amu = amu

  POP_SUB(set_amu)

end subroutine set_amu

! =========================================================================== !

!> Set value for the k-points.
subroutine set_kpts_gkq(this, kpts)
  class(gkq_t), intent(inout) :: this
  real(DP), intent(in) :: kpts(3*this%nk)

  PUSH_SUB(set_kpts_gkq)

  this%kpts = reshape(kpts, (/3, this%nk/))

  POP_SUB(set_kpts_gkq)

end subroutine set_kpts_gkq

! =========================================================================== !

!> Set values for the local q-points.
!! It is assumed that setup_paral has been called prior to this routine.
subroutine set_qpts_local_gkq(this, qpts)
  class(gkq_t), intent(inout) :: this
  real(DP), intent(in) :: qpts(3*this%nq_me)

  PUSH_SUB(set_qpts_local_gkq)

  this%qpts = reshape(qpts, (/3, this%nq_me/))

  POP_SUB(set_qpts_local_gkq)

end subroutine set_qpts_local_gkq

! =========================================================================== !

!> Set values for all q-points.
subroutine set_qpts_global_gkq(this, qpts)
  class(gkq_t), intent(inout) :: this
  real(DP), intent(in) :: qpts(3, this%nq)

  integer :: iq, jq

  PUSH_SUB(set_qpts_global_gkq)

  if (this%is_active) then
    do iq=1,this%nq_me
      jq = this%iq_start_me + iq - 1
      this%qpts(:,iq) = qpts(:,jq)
    end do
  end if

  POP_SUB(set_qpts_global_gkq)

end subroutine set_qpts_global_gkq

! =========================================================================== !
! Writing routines
! =========================================================================== !

#ifdef HDF5

!> Write the full gkq object into a file (main writing routine).
subroutine write_gkq_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id       ! File identifier
  integer :: error

  PUSH_SUB(write_gkq_hdf5)

  if (this%is_master) then
    ! Create the file and write dimensions
    call this%create_and_write_gkq_header_hdf5(fname)

    ! Write the arrays that are not distributed
    call this%write_gkq_global_arrays_hdf5(fname)
  end if

#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  ! Verify that the file exists (make sure that an error is raised otherwise)
  call open_file(99, trim(fname), status='old')
  call close_file(99)


  ! Write the main gkq matrix and the distributed arrays
#ifdef MPI
  call this%write_gkq_local_arrays_par_hdf5(fname)
  call this%write_gkq_matrix_par_hdf5(fname)
#else
  call this%write_gkq_local_arrays_ser_hdf5(fname)
  call this%write_gkq_matrix_ser_hdf5(fname)
#endif

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif

  POP_SUB(write_gkq_hdf5)

end subroutine write_gkq_hdf5

! =========================================================================== !

!> Create the file and write dimensions.
subroutine create_and_write_gkq_header_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id
  integer :: error

  PUSH_SUB(create_and_write_gkq_header_hdf5)

  ! Create a new file using default properties.
  call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, error)

  ! Create the groups
  call hdf5_create_group(file_id, 'gkq_header', error)
  call hdf5_create_group(file_id, 'gkq_data', error)

  ! Write dimensions
  call hdf5_write_int(file_id, 'gkq_header/ns',  this%ns, error)
  call hdf5_write_int(file_id, 'gkq_header/nk',  this%nk, error)
  call hdf5_write_int(file_id, 'gkq_header/nband',  this%nband, error)
  call hdf5_write_int(file_id, 'gkq_header/mband',  this%mband, error)
  call hdf5_write_int(file_id, 'gkq_header/nq',  this%nq, error)
  call hdf5_write_int(file_id, 'gkq_header/nat',  this%nat, error)
  call hdf5_write_int(file_id, 'gkq_header/ndir',  this%ndir, error)
  call hdf5_write_int(file_id, 'gkq_header/nmode',  this%nmode, error)

  ! Create the datasets for various arrays
  call hdf5_create_dset(file_id, 'gkq_data/kpts', H5T_NATIVE_DOUBLE, &
    (/3,this%nk/), error)

  call hdf5_create_dset(file_id, 'gkq_data/qpts', H5T_NATIVE_DOUBLE, &
    (/3,this%nq/), error)

  ! Create the dataset for the main matrix
  !call hdf5_create_dset(file_id, 'gkq_data/g_at', H5T_NATIVE_DOUBLE, &
  !  (/SCALARSIZE,this%mband,this%nband,this%nk,this%ns,this%nmode,this%nq/),error)

  call hdf5_create_dset(file_id, 'gkq_data/g_nu', H5T_NATIVE_DOUBLE, &
    (/SCALARSIZE,this%mband,this%nband,this%nk,this%ns,this%nmode,this%nq/),error)

  call h5fclose_f(file_id, error)

  POP_SUB(create_and_write_gkq_header_hdf5)

end subroutine create_and_write_gkq_header_hdf5

! =========================================================================== !

!> Write the global arrays of qkq.
subroutine write_gkq_global_arrays_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id       ! File identifier
  integer :: error

  PUSH_SUB(write_gkq_global_arrays_hdf5)

  if (.not. this%is_master) return

  call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error)

  ! -------------------------------
  ! Write arrays
  call hdf5_write_double_array(file_id, 'gkq_data/amu', &
                               (/this%nat/), this%amu, error)

  call hdf5_write_double_array(file_id, 'gkq_data/kpts', &
                               (/3,this%nk/), this%kpts, error)

  ! -------------------------------

  call h5fclose_f(file_id, error)

  POP_SUB(write_gkq_global_arrays_hdf5)

end subroutine write_gkq_global_arrays_hdf5

! =========================================================================== !

!> Write local (distributed) arrays using parallel hdf5.
subroutine write_gkq_local_arrays_par_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id       ! File identifier
  integer :: error
  integer(HID_T) :: dset_id     ! Dataset identifier
  integer(HID_T) :: filespace   ! Dataspace identifier in file
  integer(HID_T) :: memspace    ! Dataspace identifier in mem
  integer(HID_T) :: dspace_id   ! Dataspace identifier in mem
  integer(HID_T) :: plist_id    ! Property list identifyer
  integer(HSIZE_T) :: dims(2), fulldims(2), offset(2), offsetm(2)
  integer :: rank
  integer :: iband,jband,ik,is,idir,iat,iq,imode
  integer :: comm, info

  real(DP), allocatable :: data(:,:)

  PUSH_SUB(write_gkq_local_arrays_par_hdf5)

  comm = this%qpt_mpi_comm
  info = this%qpt_mpi_info

  ! =============================== !
  ! Set up the data to be written
  ! =============================== !

  rank = 2

  ! Dimensions of the full array
  fulldims(1) = 3
  fulldims(2) = this%nq

  ! Dimensions of the local portion of the array
  dims(1) = 3
  dims(2) = this%nq_me

  ! Offset for writing the portion of the array
  offset(:) = 0
  offset(2) = this%iq_start_me - 1

  ! Offset for accessing the data buffer
  !offsetm(:) = 0

  if (.not. this%is_active) then
    dims(:) = 1
    offset(:) = 0
  end if

  ! Copy the matrix into the data buffer
  SAFE_ALLOCATE(data,(dims(1), dims(2)))

  if (this%is_active) then

    do iq = 1, this%nq_me
      do idir = 1, 3
        data(idir,iq) = this%qpts(idir,iq)
      end do
    end do
  else
    data = 1.0_dp
  end if

  ! =============================== !
  ! Open file for parrallel IO
  ! =============================== !

  ! Create the property list to specify that we open the file for parallel IO
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
#ifdef MPI
  call h5pset_fapl_mpio_f(plist_id, comm, info, error)
#endif
  ! Open the file
  call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error, access_prp=plist_id)

  ! Close the property list. It is no longer needed.
  call h5pclose_f(plist_id, error)

  ! =============================== !
  ! Open dataset and setup memory space and file space
  ! =============================== !


  ! Create the memspace array
  call h5screate_simple_f(rank, dims, memspace, error)
  if (this%is_active) then
    ! Specify that we are going to use the full array
    call h5sselect_all_f(memspace, error)
  else
    ! Specify that we are not using any part of the array
    call h5sselect_none_f(memspace, error)
  end if

  ! Open the dataspace
  call h5dopen_f(file_id, 'gkq_data/qpts', dset_id, error)

  ! Create the filespace
  call h5dget_space_f(dset_id, filespace, error)


  ! Specify which part of the file we are writing
  if (this%is_active) then
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
  else
    call h5sselect_none_f(filespace, error)
  end if

  ! =============================== !
  ! Write the main data
  ! =============================== !

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
#ifdef MPI
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#endif
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                  memspace, filespace, xfer_prp=plist_id)

  ! =============================== !
  ! Close everything
  ! =============================== !

  call h5dclose_f(dset_id, error)  ! Dataset
  call h5sclose_f(memspace, error)  ! Memory space
  call h5sclose_f(filespace, error) ! File space
  call h5fclose_f(file_id, error)  ! File

  ! Free memory
  SAFE_DEALLOCATE(data)

  POP_SUB(write_gkq_local_arrays_par_hdf5)

end subroutine write_gkq_local_arrays_par_hdf5


!> Write local (distributed) arrays using serial hdf5.
subroutine write_gkq_local_arrays_ser_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id       ! File identifier
  integer :: error
  integer(HID_T) :: dset_id     ! Dataset identifier
  integer(HID_T) :: filespace   ! Dataspace identifier in file
  integer(HID_T) :: memspace    ! Dataspace identifier in mem
  integer(HID_T) :: dspace_id   ! Dataspace identifier in mem
  integer(HID_T) :: plist_id    ! Property list identifyer
  integer(HSIZE_T) :: dims(2), fulldims(2), offset(2), offsetm(2)
  integer :: rank
  integer :: iband,jband,ik,is,idir,iat,iq,imode

  real(DP), allocatable :: data(:,:)

  PUSH_SUB(write_gkq_local_arrays_ser_hdf5)

  ! =============================== !
  ! Set up the data to be written
  ! =============================== !

  rank = 2

  ! Dimensions of the full array
  fulldims(1) = 3
  fulldims(2) = this%nq

  ! Dimensions of the local portion of the array
  dims(1) = 3
  dims(2) = this%nq_me

  ! Offset for writing the portion of the array
  offset(:) = 0
  offset(2) = this%iq_start_me - 1

  ! Offset for accessing the data buffer
  !offsetm(:) = 0

  if (.not. this%is_active) then
    dims(:) = 1
    offset(:) = 0
  end if

  ! Copy the matrix into the data buffer
  SAFE_ALLOCATE(data,(dims(1), dims(2)))

  if (this%is_active) then

    do iq = 1, this%nq
      do idir = 1, 3
        data(idir,iq) = this%qpts(idir,iq)
      end do
    end do
  else
    data = 1.0_dp
  end if

  ! =============================== !
  ! Done setting up the data
  ! =============================== !

  ! Open the file
  call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error)

  ! Open the dataspace
  call h5dopen_f(file_id, 'gkq_data/qpts', dset_id, error)

  ! Write data
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error)

  ! Close
  call h5dclose_f(dset_id, error)  ! Dataset
  call h5fclose_f(file_id, error)  ! File

  ! Free memory
  SAFE_DEALLOCATE(data)

  POP_SUB(write_gkq_local_arrays_ser_hdf5)

end subroutine write_gkq_local_arrays_ser_hdf5

! =========================================================================== !

!> Write the gkq matrix in serial, using hdf5.
subroutine write_gkq_matrix_ser_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer :: iband,jband,ik,is,idir,iat,iq,imode
  integer :: rank
  integer :: error
  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
  integer(HID_T) :: dspace_id     ! Dataspace identifier in mem
  integer(HSIZE_T) :: dims(7), offset(7), offsetm(7)
  real(DP), allocatable :: data(:,:,:,:,:,:,:)

  PUSH_SUB(write_gkq_matrix_ser_hdf5)

  ! Note: Fortran2003 does not allow for arrays with rank larger than 7
  !       Hence we combine the dimensions iat and idir into imode
  rank = 7
  dims(1) = SCALARSIZE
  dims(2) = this%mband
  dims(3) = this%nband
  dims(4) = this%nk
  dims(5) = this%ns
  dims(6) = this%nmode
  dims(7) = this%nq

  offset(:) = 0
  offsetm(:) = 0

  SAFE_ALLOCATE(data,(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))

  do iq = 1, this%nq
    do iat = 1, this%nat
      do idir = 1, this%ndir
        imode = (iat-1) * this%ndir + idir
        do is = 1, this%ns
          do ik = 1, this%nk
            do iband = 1, this%nband
              do jband = 1, this%mband
                data(1,jband,iband,ik,is,imode,iq) = &
!&                 dble(this%g_at(jband,iband,ik,is,idir,iat,iq))
&                 dble(this%g_nu(jband,iband,ik,is,imode,iq))
#ifdef CPLX
                data(2,jband,iband,ik,is,imode,iq) = &
!&                 IMAG(this%g_at(jband,iband,ik,is,idir,iat,iq))
&                 IMAG(this%g_nu(jband,iband,ik,is,imode,iq))
#endif
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  ! Open the file
  call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error)

  ! Open the dataspace
  !call h5dopen_f(file_id, 'gkq_data/g_at', dset_id, error)
  call h5dopen_f(file_id, 'gkq_data/g_nu', dset_id, error)

  ! Write data
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error)

  call h5dclose_f(dset_id, error)  ! Dataset
  call h5fclose_f(file_id, error)  ! File

  ! Free memory
  SAFE_DEALLOCATE(data)

  POP_SUB(write_gkq_matrix_ser_hdf5)

end subroutine write_gkq_matrix_ser_hdf5

! =========================================================================== !

!> Write the gkq matrix in parallel (qpoint), using hdf5.
subroutine write_gkq_matrix_par_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id     ! File identifier
  integer :: error
  integer(HID_T) :: dset_id     ! Dataset identifier
  integer(HID_T) :: filespace   ! Dataspace identifier in file
  integer(HID_T) :: memspace    ! Dataspace identifier in mem
  integer(HID_T) :: plist_id    ! Property list identifyer
  integer(HSIZE_T) :: dims(7), fulldims(7), offset(7), offsetm(7)
  integer :: rank
  integer :: iband,jband,ik,is,idir,iat,iq,imode
  integer :: comm, info
  real(DP), allocatable :: data(:,:,:,:,:,:,:)

  PUSH_SUB(write_gkq_matrix_par_hdf5)

  comm = this%qpt_mpi_comm
  info = this%qpt_mpi_info

  ! =============================== !
  ! Set up the data to be written
  ! =============================== !

  ! Note: Fortran2003 does not allow for arrays with rank larger than 7
  !       Hence we combine the dimensions iat and idir into imode
  rank = 7

  ! Dimensions of the full matrix
  fulldims(1) = SCALARSIZE
  fulldims(2) = this%mband
  fulldims(3) = this%nband
  fulldims(4) = this%nk
  fulldims(5) = this%ns
  fulldims(6) = this%nmode
  fulldims(7) = this%nq

  ! Dimensions of the local portion of the matrix
  dims(1) = SCALARSIZE
  dims(2) = this%mband
  dims(3) = this%nband
  dims(4) = this%nk
  dims(5) = this%ns
  dims(6) = this%nmode
  dims(7) = this%nq_me

  ! Offset for writing the portion of the matrix
  offset(:) = 0
  offset(7) = this%iq_start_me - 1

  ! Offset for accessing the data buffer
  !offsetm(:) = 0

  if (.not. this%is_active) then
    dims(:) = 1
    offset(:) = 0
  end if

  ! Copy the matrix into the data buffer
  SAFE_ALLOCATE(data,(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))

  if (this%is_active) then

    do iq = 1, this%nq_me
      do iat = 1, this%nat
        do idir = 1, this%ndir
          imode = (iat-1) * this%ndir + idir
          do is = 1, this%ns
            do ik = 1, this%nk
              do iband = 1, this%nband
                do jband = 1, this%mband
                  data(1,jband,iband,ik,is,imode,iq) = &
                    dble(this%g_nu(jband,iband,ik,is,imode,iq))
#ifdef CPLX
                  data(2,jband,iband,ik,is,imode,iq) = &
                    IMAG(this%g_nu(jband,iband,ik,is,imode,iq))
#endif
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  else
    data = 0.0_dp
  end if

  ! =============================== !
  ! done setting up the data
  ! =============================== !

  ! Create the property list to specify that we open the file for parallel IO
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
#ifdef MPI
  call h5pset_fapl_mpio_f(plist_id, comm, info, error)
#endif
  ! Open the file, using the property list to specify the parallel IO
  call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error, access_prp=plist_id)

  ! Close the property list. It is no longer needed.
  call h5pclose_f(plist_id, error)


  ! Create the memspace array
  call h5screate_simple_f(rank, dims, memspace, error)
  if (this%is_active) then
    ! Specify that we are going to use the full array
    call h5sselect_all_f(memspace, error)
  else
    ! Specify that we are not using any part of the array
    call h5sselect_none_f(memspace, error)
  end if


  ! Open the dataspace
  !call h5dopen_f(file_id, 'gkq_data/g_at', dset_id, error)
  call h5dopen_f(file_id, 'gkq_data/g_nu', dset_id, error)

  ! Create the filespace
  call h5dget_space_f(dset_id, filespace, error)


  ! Specify which part of the file we are writing
  if (this%is_active) then
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
  else
    call h5sselect_none_f(filespace, error)
  end if


  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
#ifdef MPI
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#endif
  ! =============================== !
  ! Write the main data
  ! =============================== !

  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                  memspace, filespace, xfer_prp=plist_id)

  ! =============================== !
  ! Done writing the main data
  ! =============================== !

  ! Close
  call h5dclose_f(dset_id, error)  ! Dataset
  call h5sclose_f(memspace, error)  ! Memory space
  call h5sclose_f(filespace, error) ! File space
  call h5fclose_f(file_id, error)  ! File

  ! Free memory
  SAFE_DEALLOCATE(data)

  POP_SUB(write_gkq_matrix_par_hdf5)

end subroutine write_gkq_matrix_par_hdf5

! =========================================================================== !
! Reading routines
! =========================================================================== !


!> Read the gkq file in parallel
!! It is assumed that setup_paral has been called prior to this routine.
!! This mean the dimensions must be known prior to reading the file.
!! If this is not possible, then setup_paral must be called
!! after reading the dimensions.
subroutine read_gkq_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id       ! File identifier
  integer :: error

  PUSH_SUB(read_gkq_hdf5)

  ! Read and broadcast dimensions
  call this%read_and_broadcast_gkq_header_hdf5(fname)

  ! Allocate arrays
  call this%free()
  call this%alloc()

  ! Read and broadcast data
  call this%read_and_broadcast_gkq_data_hdf5(fname)

  POP_SUB(read_gkq_hdf5)

end subroutine read_gkq_hdf5

! =========================================================================== !

!> Read the dimensions then broadcast the result to everyone.
subroutine read_and_broadcast_gkq_header_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id
  integer :: error

  PUSH_SUB(read_and_broadcast_gkq_header_hdf5)

  if (this%is_master) then
    call this%read_gkq_header_hdf5(fname)
  end if
  call this%broadcast_gkq_dimensions()

  POP_SUB(read_and_broadcast_gkq_header_hdf5)

end subroutine read_and_broadcast_gkq_header_hdf5

! =========================================================================== !

!> Read the dimensions
subroutine read_gkq_header_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id
  integer :: error

  PUSH_SUB(read_gkq_header_hdf5)

  call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)

  ! Read dimensions
  call hdf5_read_int(file_id, 'gkq_header/ns',  this%ns, error)
  call hdf5_read_int(file_id, 'gkq_header/nk',  this%nk, error)
  call hdf5_read_int(file_id, 'gkq_header/nband',  this%nband, error)
  call hdf5_read_int(file_id, 'gkq_header/mband',  this%mband, error)
  call hdf5_read_int(file_id, 'gkq_header/nq',  this%nq, error)
  call hdf5_read_int(file_id, 'gkq_header/nat',  this%nat, error)
  call hdf5_read_int(file_id, 'gkq_header/ndir',  this%ndir, error)
  call hdf5_read_int(file_id, 'gkq_header/nmode',  this%nmode, error)

  call h5fclose_f(file_id, error)

  POP_SUB(read_gkq_header_hdf5)

end subroutine read_gkq_header_hdf5

! =========================================================================== !

!> Broadcast the dimensions to all workers in the group.
subroutine broadcast_gkq_dimensions(this, comm)

  class(gkq_t), intent(inout) :: this
  integer, intent(in), optional :: comm

  integer :: ndim=7
  integer :: dims(7)
  integer :: comm_

  PUSH_SUB(broadcast_gkq_dimensions)

#ifdef MPI
  if (present(comm)) then
    comm_ = comm
  else
    comm_ = MPI_COMM_WORLD
  end if
#else
  comm_ = 0
#endif

  if (this%is_master) then
     dims = (/this%mband, this%nband, this%nk, this%ns, this%ndir, this%nat, this%nq/)
  end if

#ifdef MPI
  call MPI_BCAST(dims, ndim, MPI_INTEGER, 0, comm_, mpierr)
  call MPI_BARRIER(comm_, mpierr)
#endif

  if (.not. this%is_master) then
     this%mband = dims(1)
     this%nband = dims(2)
     this%nk = dims(3)
     this%ns = dims(4)
     this%ndir = dims(5)
     this%nat = dims(6)
     this%nq = dims(7)
     this%nmode = this%nat * 3
  end if

  POP_SUB(broadcast_gkq_dimensions)

end subroutine broadcast_gkq_dimensions

! =========================================================================== !

!> Read the arrays in serial or parallel, and broadcast those read in serial.
subroutine read_and_broadcast_gkq_data_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id       ! File identifier
  integer :: error

  PUSH_SUB(read_and_broadcast_gkq_data_hdf5)

  ! Read and broadcast arrays that are not distributed
  call this%read_and_broadcast_gkq_global_arrays_hdf5(fname)

  ! Read the main gkq matrix and the local arrays
#ifdef MPI
  call this%read_gkq_local_arrays_par_hdf5(fname)
  call this%read_gkq_matrix_par_hdf5(fname)
#else
  call this%read_gkq_local_arrays_ser_hdf5(fname)
  call this%read_gkq_matrix_ser_hdf5(fname)
#endif

  POP_SUB(read_and_broadcast_gkq_data_hdf5)

end subroutine read_and_broadcast_gkq_data_hdf5

! =========================================================================== !

!> Master reads the global (non-distributed) arrays and broadcast it.
subroutine read_and_broadcast_gkq_global_arrays_hdf5(this, fname, comm)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname
  integer, intent(in), optional :: comm

  integer :: comm_
  integer(HID_T) :: file_id
  integer :: error

  PUSH_SUB(read_and_broadcast_gkq_global_arrays_hdf5)

#ifdef MPI
  if (present(comm)) then
    comm_ = comm
  else
    comm_ = MPI_COMM_WORLD
  end if
#else
  comm_ = 0
#endif

  if (this%is_master) then
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
  end if

  call this%read_and_broadcast_amu(file_id, comm=comm_)
  call this%read_and_broadcast_kpts(file_id, comm=comm_)

  if (this%is_master) then
    call h5fclose_f(file_id, error)
  end if

  POP_SUB(read_and_broadcast_gkq_global_arrays_hdf5)

end subroutine read_and_broadcast_gkq_global_arrays_hdf5

! =========================================================================== !

!> Master reads amu and broadcasts the array.
subroutine read_and_broadcast_amu(this, file_id, comm)
  class(gkq_t), intent(inout) :: this
  integer(HID_T), intent(in) :: file_id
  integer, intent(in), optional :: comm

  integer :: comm_
  integer :: error

  PUSH_SUB(read_and_broadcast_amu)

#ifdef MPI
  if (present(comm)) then
    comm_ = comm
  else
    comm_ = MPI_COMM_WORLD
  end if
#else
  comm_ = 0
#endif

  if (this%is_master) then
    call hdf5_read_double_array(file_id, 'gkq_data/amu',(/this%nat/), this%amu, error)
  end if

#ifdef MPI
  call MPI_BCAST(this%amu, this%nat, MPI_DOUBLE, 0, comm_, mpierr)
  call MPI_BARRIER(comm_,mpierr)
#endif

  POP_SUB(read_and_broadcast_amu)

end subroutine read_and_broadcast_amu

! =========================================================================== !

!> Master reads amu.
subroutine read_and_broadcast_kpts(this, file_id, comm)
  class(gkq_t), intent(inout) :: this
  integer(HID_T), intent(in) :: file_id       ! File identifier
  integer, intent(in), optional :: comm

  integer :: comm_
  integer :: rank=1
  integer :: error

  PUSH_SUB(read_and_broadcast_kpts)

#ifdef MPI
  if (present(comm)) then
    comm_ = comm
  else
    comm_ = MPI_COMM_WORLD
  end if
#else
  comm_ = 0
#endif

  ! Master reads
  if (this%is_master) then

    ! Read arrays
    call hdf5_read_double_array(file_id, 'gkq_data/kpts', &
    &                            (/3,this%nk/), this%kpts, error)

  end if

  ! Broadcast
#ifdef MPI
  call MPI_BCAST(this%kpts, 3*this%nk, MPI_DOUBLE, 0, comm_, mpierr)
  call MPI_BARRIER(comm_,mpierr)
#endif

  POP_SUB(read_and_broadcast_kpts)

end subroutine read_and_broadcast_kpts

! =========================================================================== !

!> Read the gkq matrix in serial, using hdf5.
subroutine read_gkq_matrix_ser_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer :: iband,jband,ik,is,idir,iat,iq,imode
  integer :: rank
  integer :: error
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
  integer(HSIZE_T) :: dims(7), offset(7), offsetm(7)
  real(DP), allocatable :: data(:,:,:,:,:,:,:)

  PUSH_SUB(read_gkq_matrix_ser_hdf5)

  ! Note: Fortran2003 does not allow for arrays with rank larger than 7
  !       Hence we combine the dimensions iat and idir into imode
  rank = 7
  dims(1) = SCALARSIZE
  dims(2) = this%mband
  dims(3) = this%nband
  dims(4) = this%nk
  dims(5) = this%ns
  dims(6) = this%nmode
  dims(7) = this%nq

  offset(:) = 0

  SAFE_ALLOCATE(data,(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))

  ! Open the file
  call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error)

  !call h5dopen_f(file_id, 'gkq_data/g_at', dset_id, error)
  call h5dopen_f(file_id, 'gkq_data/g_nu', dset_id, error)
  call h5dget_space_f(dset_id,filespace,error)
  call h5screate_simple_f(rank, dims, memspace, error)
  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, memspace, filespace)

  do iq = 1, this%nq
    do iat = 1, this%nat
      do idir = 1, this%ndir
        imode = (iat-1) * this%ndir + idir
        do is = 1, this%ns
          do ik = 1, this%nk
            do iband = 1, this%nband
              do jband = 1, this%mband
                !this%g_at(jband,iband,ik,is,idir,iat,iq) = &
                this%g_nu(jband,iband,ik,is,imode,iq) = &
                  SCALARIFY2(data(1,jband,iband,ik,is,imode,iq), data(2,jband,iband,ik,is,imode,iq))
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5dclose_f(dset_id,error)
  call h5fclose_f(file_id, error)

  SAFE_DEALLOCATE(data)

  POP_SUB(read_gkq_matrix_ser_hdf5)

end subroutine read_gkq_matrix_ser_hdf5

! =========================================================================== !

!> Read the gkq matrix in parallel, using hdf5.
subroutine read_gkq_matrix_par_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id
  integer :: error
  integer :: iband,jband,ik,is,idir,iat,iq,imode
  integer :: rank
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: plist_id      ! Property list identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
  integer(HSIZE_T) :: dims(7), fulldims(7), offset(7), offsetm(7)
  integer :: info, comm
  real(DP), allocatable :: data(:,:,:,:,:,:,:)

  PUSH_SUB(read_gkq_matrix_par_hdf5)

  comm = this%qpt_mpi_comm
  info = this%qpt_mpi_info

  ! Note: Fortran2003 does not allow for arrays with rank larger than 7
  !       Hence we combine the dimensions iat and idir into imode
  rank = 7

  ! Dimensions of the full matrix
  fulldims(1) = SCALARSIZE
  fulldims(2) = this%mband
  fulldims(3) = this%nband
  fulldims(4) = this%nk
  fulldims(5) = this%ns
  fulldims(6) = this%nmode
  fulldims(7) = this%nq

  ! Dimensions of the local portion of the matrix
  dims(1) = SCALARSIZE
  dims(2) = this%mband
  dims(3) = this%nband
  dims(4) = this%nk
  dims(5) = this%ns
  dims(6) = this%nmode
  dims(7) = this%nq_me

  ! Offset for reading the portion of the matrix
  offset(:) = 0
  offset(7) = this%iq_start_me - 1

  if (.not. this%is_active) then
    dims(:) = 1
    offset(:) = 0
  end if

  SAFE_ALLOCATE(data,(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))

  ! Create the property list to specify that we open the file for parallel IO
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
#ifdef MPI
  call h5pset_fapl_mpio_f(plist_id, comm, info, error)
#endif
  ! Open the file, using the property list to specify the parallel IO
  call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, error, access_prp=plist_id)

  ! Close the property list. It is no longer needed.
  call h5pclose_f(plist_id, error)

  ! Create the memspace array
  call h5screate_simple_f(rank, dims, memspace, error)
  if (this%is_active) then
    ! Specify that we are going to use the full array
    call h5sselect_all_f(memspace, error)
  else
    ! Specify that we are not using any part of the array
    call h5sselect_none_f(memspace, error)
  end if

  ! Open the dataspace
  !call h5dopen_f(file_id, 'gkq_data/g_at', dset_id, error)
  call h5dopen_f(file_id, 'gkq_data/g_nu', dset_id, error)

  ! Create the filespace
  call h5dget_space_f(dset_id, filespace, error)

  ! Specify which part of the file we are reading
  if (this%is_active) then
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
  else
    call h5sselect_none_f(filespace, error)
  end if

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
#ifdef MPI
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#endif
  ! =============================== !
  ! Read the main data
  ! =============================== !

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                  memspace, filespace, xfer_prp=plist_id)

  ! =============================== !
  ! Done reading the main data
  ! =============================== !

  ! Close
  call h5dclose_f(dset_id, error)  ! Dataset
  call h5sclose_f(memspace, error)  ! Memory space
  call h5sclose_f(filespace, error) ! File space
  call h5fclose_f(file_id, error)  ! File

  ! Store the data
  if (this%is_active) then
    do iq = 1, this%nq_me
      do iat = 1, this%nat
        do idir = 1, this%ndir
          imode = (iat-1) * this%ndir + idir
          do is = 1, this%ns
            do ik = 1, this%nk
              do iband = 1, this%nband
                do jband = 1, this%mband
                  !this%g_at(jband,iband,ik,is,idir,iat,iq) = &
                  this%g_nu(jband,iband,ik,is,imode,iq) = &
                    SCALARIFY2(data(1,jband,iband,ik,is,imode,iq), data(2,jband,iband,ik,is,imode,iq))
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end if

  SAFE_DEALLOCATE(data)

  POP_SUB(read_gkq_matrix_par_hdf5)

end subroutine read_gkq_matrix_par_hdf5

! =========================================================================== !

!> Read the local (distributed) arrays using parallel hdf5.
subroutine read_gkq_local_arrays_par_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id       ! File identifier
  integer :: error
  integer(HID_T) :: dset_id     ! Dataset identifier
  integer(HID_T) :: filespace   ! Dataspace identifier in file
  integer(HID_T) :: memspace    ! Dataspace identifier in mem
  integer(HID_T) :: dspace_id   ! Dataspace identifier in mem
  integer(HID_T) :: plist_id    ! Property list identifyer
  integer(HSIZE_T) :: dims(2), fulldims(2), offset(2), offsetm(2)
  integer :: rank
  integer :: iband,jband,ik,is,idir,iat,iq,imode
  integer :: comm, info
  real(DP), allocatable :: data(:,:)

  PUSH_SUB(read_gkq_local_arrays_par_hdf5)

  comm = this%qpt_mpi_comm
  info = this%qpt_mpi_info

  rank = 2

  ! Dimensions of the full array
  fulldims(1) = 3
  fulldims(2) = this%nq

  ! Dimensions of the local portion of the array
  dims(1) = 3
  dims(2) = this%nq_me

  ! Offset for writing the portion of the array
  offset(:) = 0
  offset(2) = this%iq_start_me - 1

  ! Offset for accessing the data buffer
  !offsetm(:) = 0

  if (.not. this%is_active) then
    dims(:) = 1
    offset(:) = 0
  end if

  ! Allocate memory
  SAFE_ALLOCATE(data,(dims(1), dims(2)))

  ! Create the property list to specify that we open the file for parallel IO
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
#ifdef MPI
  call h5pset_fapl_mpio_f(plist_id, comm, info, error)
#endif
  ! Open the file, using the property list to specify the parallel IO
  call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error, access_prp=plist_id)

  ! Close the property list. It is no longer needed.
  call h5pclose_f(plist_id, error)


  ! Create the memspace array
  call h5screate_simple_f(rank, dims, memspace, error)
  if (this%is_active) then
    ! Specify that we are going to use the full array
    call h5sselect_all_f(memspace, error)
  else
    ! Specify that we are not using any part of the array
    call h5sselect_none_f(memspace, error)
  end if


  ! Open the dataspace
  call h5dopen_f(file_id, 'gkq_data/qpts', dset_id, error)

  ! Create the filespace
  call h5dget_space_f(dset_id, filespace, error)


  ! Specify which part of the file we are writing
  if (this%is_active) then
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
  else
    call h5sselect_none_f(filespace, error)
  end if


  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
#ifdef MPI
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#endif
  ! =============================== !
  ! Read the main data
  ! =============================== !

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                  memspace, filespace, xfer_prp=plist_id)

  if (this%is_active) then
    this%qpts = data
  end if

  ! =============================== !
  ! Done reading the main data
  ! =============================== !

  ! Close
  call h5dclose_f(dset_id, error)  ! Dataset
  call h5sclose_f(memspace, error)  ! Memory space
  call h5sclose_f(filespace, error) ! File space
  call h5fclose_f(file_id, error)  ! File

  ! Free memory
  SAFE_DEALLOCATE(data)

  POP_SUB(read_gkq_local_arrays_par_hdf5)

end subroutine read_gkq_local_arrays_par_hdf5


!> Read the local (distributed) arrays using serial hdf5.
subroutine read_gkq_local_arrays_ser_hdf5(this, fname)
  class(gkq_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id       ! File identifier
  integer :: error
  integer(HID_T) :: dset_id     ! Dataset identifier
  integer(HID_T) :: filespace   ! Dataspace identifier in file
  integer(HID_T) :: memspace    ! Dataspace identifier in mem
  integer(HID_T) :: dspace_id   ! Dataspace identifier in mem
  integer(HID_T) :: plist_id    ! Property list identifyer
  integer(HSIZE_T) :: dims(2), fulldims(2), offset(2), offsetm(2)
  integer :: rank
  integer :: iband,jband,ik,is,idir,iat,iq,imode
  real(DP), allocatable :: data(:,:)

  PUSH_SUB(read_gkq_local_arrays_ser_hdf5)

  rank = 2

  ! Dimensions of the full array
  fulldims(1) = 3
  fulldims(2) = this%nq

  ! Dimensions of the local portion of the array
  dims(1) = 3
  dims(2) = this%nq_me

  ! Offset for writing the portion of the array
  offset(:) = 0
  offset(2) = this%iq_start_me - 1

  ! Offset for accessing the data buffer
  !offsetm(:) = 0

  if (.not. this%is_active) then
    dims(:) = 1
    offset(:) = 0
  end if

  ! Allocate memory
  SAFE_ALLOCATE(data,(dims(1), dims(2)))

  ! Open the file
  call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error)

  ! Open the dataspace
  call h5dopen_f(file_id, 'gkq_data/qpts', dset_id, error)

  ! Create the filespace
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_all_f(filespace, error)

  ! Create the memspace array
  call h5screate_simple_f(rank, dims, memspace, error)

  ! =============================== !
  ! Read the main data
  ! =============================== !

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                  memspace, filespace)

  this%qpts = data

  ! =============================== !
  ! Done reading the main data
  ! =============================== !

  ! Close
  call h5dclose_f(dset_id, error)  ! Dataset
  call h5sclose_f(memspace, error)  ! Memory space
  call h5sclose_f(filespace, error) ! File space
  call h5fclose_f(file_id, error)  ! File

  ! Free memory
  SAFE_DEALLOCATE(data)

  POP_SUB(read_gkq_local_arrays_ser_hdf5)

end subroutine read_gkq_local_arrays_ser_hdf5

#endif

end module gkq_m
