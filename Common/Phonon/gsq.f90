!==============================================================================
!
! Module gsq_m
!
! Originally by GKA (2018)
!
!==============================================================================

#include "f_defs.h"

module gsq_m

#ifdef HDF5
use hdf5
use hdf5_io_m
#endif

use global_m

implicit none

public

! =========================================================================== !

!> exciton-phonon coupling matrix elements
type gsq_t

  ! Dimensions
  integer :: ns       !< Number of spins
  integer :: nqxct    !< Number of Q-points (exciton momentum)
  integer :: nband    !< Number of bands at Q
  integer :: mband    !< Number of bands at Q+q
  integer :: nq       !< Number of q-points (phonons)
  integer :: nat      !< number of atoms
  integer :: ndir=3   !< number of spacial (reduced) directions
  integer :: nmode    !< nat * ndir
  integer :: nq_me    !< Number of q-points own by the local worker

  ! Parallel info
  logical :: is_master
  logical :: is_active
  integer :: iq_start_me    !< index of the first q-point
  integer :: qpt_mpi_comm   !< qpt communicator
  integer :: qpt_mpi_info   !< qpt info

  !> masses, in atomic mass units
  !! amu(nat)
  real(DP), allocatable :: amu(:)

  !> Exciton Q-points, in reduced coordinates
  !! qpts_xct(3,nk)
  real(DP), allocatable :: qpts_xct(:,:)

  !> q-points, in reduced coordinates
  !! qpts(3,nq_me)
  real(DP), allocatable :: qpts(:,:)

  !> Electron-phonon coupling matrix elemnts, in atom/direction basis [Hartree]
  !! gsq(mband,nband,nqxct,ns,ndir,nat,nq)
  SCALAR, allocatable :: gsq_at(:,:,:,:,:,:,:)

  !> Electron-phonon coupling matrix elemnts, in mode basis [Hartree]
  !! gsq(mband,nband,nqxct,ns,nmode,nq)
  SCALAR, allocatable :: gsq_nu(:,:,:,:,:,:)

  contains

  ! Core procedures
  procedure, pass :: init => init_gsq
  procedure, pass :: alloc => alloc_gsq
  procedure, pass :: free => free_gsq
  procedure, pass :: copy => copy_gsq
  procedure, pass :: compare => compare_gsq
  procedure, pass :: setup_paral => setup_paral_gsq

  ! Setting procedures
  procedure, pass :: set_amu
  procedure, pass :: set_qpts_local => set_qpts_local_gsq
  procedure, pass :: set_qpts_global => set_qpts_global_gsq
  procedure, pass :: set_kpts => set_kpts_gsq

end type gsq_t

! =========================================================================== !
contains
! =========================================================================== !

!> Initialize dimensions of the arrays.
subroutine init_gsq(this, ns, nk, nband, mband, nq, nat, nq_me)
  class(gsq_t), intent(inout) :: this
  integer, intent(in) :: ns       !< Number of spins
  integer, intent(in) :: nk       !< Number of k-points
  integer, intent(in) :: nband    !< Number of bands at k
  integer, intent(in) :: mband    !< Number of bands at k+q
  integer, intent(in) :: nq       !< Number of q-points
  integer, intent(in) :: nat      !< number of atoms
  integer, intent(in),optional :: nq_me    !< Number of q-points

  PUSH_SUB(init_gsq)

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

  POP_SUB(init_gsq)

end subroutine init_gsq

! =========================================================================== !

!> Setup the parallel layout of the data.
!! At the moment, only parallelism over q-points is supported.
subroutine setup_paral_gsq(this, peinf, np_qpt)
  class(gsq_t), intent(inout) :: this
  type(peinfo), intent(in) :: peinf     !< info on parallelization
  integer, intent(in) :: np_qpt         !< number of q-point groups

  integer :: inode
  integer :: nq_me_max

  PUSH_SUB(setup_paral_gsq)

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

  this%qpt_mpi_comm = MPI_COMM_WORLD
  this%qpt_mpi_info = MPI_INFO_NULL

  POP_SUB(setup_paral_gsq)

end subroutine setup_paral_gsq

! =========================================================================== !

!> Allocate memory.
!! It is assumed that setup_paral has been called prior to this routine.
subroutine alloc_gsq(this)
  class(gsq_t), intent(inout) :: this

  PUSH_SUB(alloc_gsq)

  if (.not. allocated(this%amu)) then
    SAFE_ALLOCATE(this%amu, (this%nat))
  end if

  if (.not. allocated(this%kpts)) then
    SAFE_ALLOCATE(this%kpts, (3, this%nk))
  end if

  if (.not. allocated(this%qpts)) then
    SAFE_ALLOCATE(this%qpts, (3, this%nq_me))
  end if

  if (.not. allocated(this%gsq_at)) then
    SAFE_ALLOCATE(this%gsq_at, (this%mband, this%nband, this%nk, this%ns, this%ndir, this%nat, this%nq_me))
  end if

  if (.not. allocated(this%gsq_nu)) then
    SAFE_ALLOCATE(this%gsq_nu, (this%mband, this%nband, this%nk, this%ns, this%nmode, this%nq_me))
  end if

  POP_SUB(alloc_gsq)

end subroutine alloc_gsq

! =========================================================================== !

!> Free memory.
subroutine free_gsq(this)

  class(gsq_t), intent(inout) :: this

  PUSH_SUB(free_gsq)

  SAFE_DEALLOCATE(this%amu)
  SAFE_DEALLOCATE(this%kpts)
  SAFE_DEALLOCATE(this%qpts)
  SAFE_DEALLOCATE(this%gsq_at)
  SAFE_DEALLOCATE(this%gsq_nu)

  POP_SUB(free_gsq)

end subroutine free_gsq

! =========================================================================== !

!> Copy object onto another.
subroutine copy_gsq(this, other)
  class(gsq_t), intent(inout) :: this
  class(gsq_t), intent(out) :: other

  PUSH_SUB(copy_gsq)

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
  other%gsq_at = this%gsq_at
  other%gsq_nu = this%gsq_nu

  POP_SUB(copy_gsq)

end subroutine copy_gsq

! =========================================================================== !

!> Compare a gsq object against another, and report any discrepancy.
!! Mostly for testing purposes.
subroutine compare_gsq(this, other, ndiff, verbose)
  class(gsq_t), intent(in) :: this
  class(gsq_t), intent(in) :: other
  integer, intent(out) :: ndiff
  logical, intent(in), optional :: verbose

  logical :: verbose_ = .True.
  real(DP), parameter :: tol=1.0d-5
  real(DP) :: checksum
  integer :: unt=6

  if (present(verbose)) verbose_ = verbose

  PUSH_SUB(compare_gsq)

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
  checksum = abs(sum(this%gsq_at) - sum(other%gsq_at))
  if (checksum .gt. tol) then
    if (verbose_) write(unt,*) 'Error in checksum of: gsq_at. ','delta=',checksum
    ndiff = ndiff + 1
  end if

  POP_SUB(compare_gsq)

end subroutine compare_gsq

! =========================================================================== !

!> Set value for the masses.
subroutine set_amu(this, amu)
  class(gsq_t), intent(inout) :: this
  real(DP), intent(in) :: amu(this%nat)

  PUSH_SUB(set_amu)

  this%amu = amu

  POP_SUB(set_amu)

end subroutine set_amu

! =========================================================================== !

!> Set value for the k-points.
subroutine set_kpts_gsq(this, kpts)
  class(gsq_t), intent(inout) :: this
  real(DP), intent(in) :: kpts(3*this%nk)

  PUSH_SUB(set_kpts_gsq)

  this%kpts = reshape(kpts, (/3, this%nk/))

  POP_SUB(set_kpts_gsq)

end subroutine set_kpts_gsq

! =========================================================================== !

!> Set values for the local q-points.
!! It is assumed that setup_paral has been called prior to this routine.
subroutine set_qpts_local_gsq(this, qpts)
  class(gsq_t), intent(inout) :: this
  real(DP), intent(in) :: qpts(3*this%nq_me)

  PUSH_SUB(set_qpts_local_gsq)

  this%qpts = reshape(qpts, (/3, this%nq_me/))

  POP_SUB(set_qpts_local_gsq)

end subroutine set_qpts_local_gsq

! =========================================================================== !

!> Set values for all q-points.
subroutine set_qpts_global_gsq(this, qpts)
  class(gsq_t), intent(inout) :: this
  real(DP), intent(in) :: qpts(3, this%nq)

  integer :: iq, jq

  PUSH_SUB(set_qpts_global_gsq)

  if (this%is_active) then
    do iq=1,this%nq_me
      jq = this%iq_start_me + iq - 1
      this%qpts(:,iq) = qpts(:,jq)
    end do
  end if

  POP_SUB(set_qpts_global_gsq)

end subroutine set_qpts_global_gsq

end module gsq_m
