!==============================================================================
!
! Module evecs_m
!
! Created by GKA (2018).
!
! Objects holding the exciton wavefunctions, with corresponding IO routines.
!
!==============================================================================


#include "f_defs.h"

module evecs_m

#ifdef HDF5
use hdf5
use hdf5_io_m
use wfn_io_hdf5_m
#endif
use global_m
use message_m

implicit none

public

! =========================================================================== !

!> BSE eigenvectors and eigenvalues
type evecs_t

  ! Dimensions
  integer :: ns       !< Number of spins
  integer :: nv       !< Number of valence bands
  integer :: nc       !< Number of conduction bands
  integer :: nk       !< Number of k-points
  integer :: nmat     !< Size of the BSE hamiltonian (maximum number of eigenvalues)
                      !! We expect nmat=ns*nk*nc*nv 
  integer :: usize    !< Size of the eigenvectors (nmat or 2*nmat)
  integer :: neig     !< Total number of eigenvectors
  integer :: meig=1   !< Number of eigenvectors held locally (array dimension)
  integer :: nq=1     !< Number of q-points (hard-coded)
  integer :: krnl=1   !< Kernel type
                      !! krnl = 0 --> spin triplet kernel, direct kernel only (only allow for nspin = 1)
                      !!        1 --> spin singlet kernel (default)
                      !!        2 --> local-fields + RPA, exchange kernel only
                      !!        3 --> spinor kernel

  logical :: tda = .true.            !< Whether Tamm-Dancoff approximation is used
  logical :: has_Avc = .false.       !< Whether the object holds eigenvectors as Avc 
  logical :: has_flat_Avc = .false.  !< Whether the object holds the flat eigenvectors
                                     !! Note: eigenvectors are named u when they are stored
                                     !! as flat arrays, and Avc when they are stored as (s,v,c,k).

  type(mf_header_t) :: mf_header     !< The mean field header
  logical :: has_mf_header = .false. !< Whether the mean field header is initialized

  !> q-points, or center-of-mass momentum of exciton (qpt=0 for optical excitons).
  !! qpts(3,nq)
  real(DP), allocatable :: qpts(:,:)

  !> k-points
  !! kpts(3,nk)
  real(DP), allocatable :: kpts(:,:)

  !> BSE eigenvalues. This array is not distributed. Each node has all eigenvalues
  !! evals(neig)
  real(DP), allocatable :: evals(:)

  !> Right eigenvectors (flat)
  !! u_r(usize,meig)
  SCALAR, allocatable :: u_r(:,:)

  !> Left eigenvectors (flat)
  !! u_l(usize,meig)
  SCALAR, allocatable :: u_l(:,:)

  ! GKA: At the moment, the number of q-points (nq) is hard-coded to 1, and
  !      for simplicity, the e-h coefficient arrays do not have an nq dimension.
  !      In the hdf5 file, however, they do have this extra dimension.

  !> Right e-h coefficients
  !! Avc(ns,nv,nc,nk,meig)
  SCALAR, allocatable :: Avc(:,:,:,:,:)

  !> Left e-h coefficients
  !! Avc_l(ns,nv,nc,nk,meig)
  SCALAR, allocatable :: Avc_l(:,:,:,:,:)

  !> Right e-h coefficients, deexcitation component
  !! Bvc_r(ns,nv,nc,nk,meig)
  SCALAR, allocatable :: Bvc_r(:,:,:,:,:)

  !> Left e-h coefficients, deexcitation component
  !! Bvc_l(ns,nv,nc,nk,meig)
  SCALAR, allocatable :: Bvc_l(:,:,:,:,:)

  ! Parallel info
  logical :: is_master = .true.           !< Whether we are on the master node
  logical :: is_distributed = .false.     !< Whether each node owns a different set of eigenvalues
  integer :: npes=1                       !< Number of processors
  integer :: inode=0                      !< Rank of this node
#ifdef MPI
  integer :: comm=MPI_COMM_WORLD          !< MPI communicator
  integer :: mpi_info=MPI_INFO_NULL       !< MPI info
#else
  integer :: comm=0                       !< MPI communicator
  integer :: mpi_info=0                   !< MPI info
#endif
  integer :: global_meig                  !< The maximum value of meig accross processors
  integer, allocatable :: global_ieigs(:) !< The global eigenvalue indices of those owned locally
  integer, allocatable :: who_owns(:)     !< The index of the node holding all eigenvector
  integer, allocatable :: local_ieigs(:)  !< The local index on the node of all eigenvector
                                          !! This is a global array.

  ! IO parameters
  integer :: file_version=1     !< HDF5 file version
  logical :: is_open = .false.  !< Whether the file is open for sequential reading
  integer :: nwrite=0           !< Number of eigenvectors to be written
  integer :: ieig_read=0        !< Reading counter
  integer :: uur = 11           !< Unit number for right eigenvectors
  integer :: uul = 12           !< Unit number for left eigenvectors

  character(30) :: default_fname_bin='eigenvectors'
  character(30) :: default_fname_hdf5='eigenvectors.h5'

#ifdef HDF5
  integer(HID_T) :: file_id
  logical :: use_hdf5 = .true.
#else
  logical :: use_hdf5 = .false.
#endif


  ! ------------------------------------------------------------------------- !
  contains
  ! ------------------------------------------------------------------------- !

  ! Initialization
  procedure :: init_from_xctinfo => evecs_init_from_xctinfo
  procedure :: init => evecs_init

  ! Allocation / deallocation
  procedure :: alloc => evecs_alloc
  procedure :: alloc_global_arrays => evecs_alloc_global_arrays
  procedure :: alloc_local_arrays => evecs_alloc_local_arrays
  procedure :: free => evecs_free

  ! I/O interface
  procedure :: write_file => evecs_write
  procedure :: read_file => evecs_read
  procedure :: open_read => evecs_open_read
  procedure :: close_file => evecs_close
  procedure :: read_header => evecs_read_header
  procedure :: read_next_eigenvector => evecs_read_next_eigenvector

  ! Output
  procedure :: print_out_header_info => evecs_print_out_header_info
  procedure :: print_out_dimensions => evecs_print_out_dimensions

  ! Arrays manipulation
  procedure :: reshape_Avc => evecs_reshape_Avc

  ! Parallel setup
  procedure :: set_mpi_comm => evecs_set_mpi_comm
  procedure :: setup_paral_default => evecs_setup_paral_default
  procedure :: setup_paral_distribution => evecs_setup_paral_distribution
  procedure :: get_BSE_paral => evecs_get_BSE_paral

  ! MPI communication
  procedure :: broadcast_header => evecs_broadcast_header
  procedure :: broadcast_global_arrays => evecs_broadcast_global_arrays
  procedure :: broadcast_local_arrays => evecs_broadcast_local_arrays

  ! Testing
  procedure :: copy => evecs_copy
  procedure :: compare => evecs_compare

  ! Binary I/O
  procedure :: evecs_write_bin
  procedure :: evecs_read_bin
  procedure :: evecs_read_and_broadcast_bin
  procedure :: evecs_open_read_bin
  procedure :: evecs_close_bin
  procedure :: evecs_read_header_bin
  procedure :: evecs_read_dimensions_bin
  procedure :: evecs_count_neig_and_rewind_bin
  procedure :: evecs_read_single_eigenvector_bin
  procedure :: evecs_skip_eigenvectors_bin
  procedure :: evecs_read_and_distribute_bin

  ! HDF5 I/O
#ifdef HDF5
  procedure :: evecs_write_hdf5
  procedure :: evecs_create_and_write_header_hdf5
  procedure :: evecs_write_global_arrays_hdf5
  procedure :: evecs_write_local_arrays_ser_hdf5
  procedure :: evecs_write_local_arrays_par_hdf5
  procedure :: evecs_write_single_eigenvector_par_hdf5

  procedure :: evecs_read_hdf5
  procedure :: evecs_open_read_hdf5
  procedure :: evecs_close_hdf5
  procedure :: evecs_read_header_hdf5
  procedure :: evecs_read_global_arrays_hdf5
  procedure :: evecs_read_local_arrays_ser_hdf5
  procedure :: evecs_read_local_arrays_par_hdf5
  procedure :: evecs_read_single_eigenvector_par_hdf5
#endif

end type evecs_t

! =========================================================================== !
contains
! =========================================================================== !

!> Initialize manually, mostly for testing purposes.
subroutine evecs_init(this, ns, nk, nv, nc, neig, meig, tda, mf_header)
  class(evecs_t), intent(inout) :: this
  integer, intent(in) :: ns   !< Number of spin
  integer, intent(in) :: nk   !< Number of k-points
  integer, intent(in) :: nv   !< Number of valence bands
  integer, intent(in) :: nc   !< Number of conduction bands
  integer, intent(in), optional :: neig !< Number of eigenvalues
  integer, intent(in), optional :: meig !< Number of eigenvectors for array dimension
                                        !! Cannot be larger than neig
  logical, intent(in), optional :: tda  !< Use Tamm-Dancoff approximation
  type (mf_header_t), intent(in), optional :: mf_header !< The mean field header

  PUSH_SUB(evecs_init)

  this%ns = ns
  this%nk = nk
  this%nv = nv
  this%nc = nc

  if (present(neig)) then
    this%neig = neig
  else
    this%neig = this%ns * this%nk * this%nv * this%nc
  end if

  if (present(meig)) then
    this%meig = meig
  else
    this%meig = this%neig
  end if

  this%tda = .true.
  if (present(tda)) this%tda = tda

  this%nmat = this%ns * this%nk * this%nv * this%nc
  if (this%tda) then
    this%usize = this%nmat
  else
    this%usize = 2 * this%nmat
  end if

  SAFE_ALLOCATE(this%kpts, (3, this%nk))
  this%kpts = 0.0D0

  SAFE_ALLOCATE(this%qpts, (3, this%nq))
  this%qpts = 0.0D0

  this%ieig_read = 0
  this%nwrite = 0

  if (present(mf_header)) then
    this%mf_header = mf_header
    this%has_mf_header = .true.
  end if

  ! Get parallel info
  call this%setup_paral_default()

  POP_SUB(evecs_init)

end subroutine evecs_init

! =========================================================================== !

!> Initialize the eigenvectors using the xctinfo object typically found
!! in a BSE calculation.
subroutine evecs_init_from_xctinfo(this, xct, kg, nmat, neig, meig, mf_header)
  class(evecs_t), intent(inout) :: this
  type (xctinfo), intent(in) :: xct   !< Information on the BSE calculation
  type (grid), intent(in) :: kg       !< Information on the k-point grid
  integer, intent(in) :: nmat
  integer, intent(in) :: neig
  integer, intent(in) :: meig
  type (mf_header_t), intent(in), optional :: mf_header !< The mean field header

  integer :: ii, ik

  PUSH_SUB(evecs_init_from_xctinfo)

  this%ieig_read = 0

  this%nmat = nmat
  this%neig = neig
  this%meig = meig
  this%nwrite = neig

  this%ns = xct%nspin
  this%nv = xct%nvb_fi
  this%nc = xct%ncb_fi
  this%nk = xct%nkpt_fi
  this%tda = xct%tda
  this%krnl = xct%krnl

  this%use_hdf5 = (xct%use_hdf5 .and. xct%use_hdf5_output)

  SAFE_ALLOCATE(this%kpts, (3, this%nk))
  this%kpts = kg%f

  SAFE_ALLOCATE(this%qpts, (3, this%nq))
  this%qpts(:,:) = 0.0D0
  if (xct%qflag == 1) then
    continue
  else if (xct%qflag == 2 .or. xct%qflag == 0) then
    this%qpts(:,1) = xct%finiteq
  end if

  if (present(mf_header)) then
    this%mf_header = mf_header
    this%has_mf_header = .true.
  end if

  ! Get parallel info
  call this%get_BSE_paral()

  POP_SUB(evecs_init_from_xctinfo)

end subroutine evecs_init_from_xctinfo

! =========================================================================== !

!> Copy object into a new instance.
!! For testing purposes.
subroutine evecs_copy(this, other)
  class(evecs_t), intent(inout) :: this
  class(evecs_t), intent(out) :: other
  PUSH_SUB(evecs_copy)

  other%ns = this%ns
  other%nk = this%nk
  other%nc = this%nc
  other%nv = this%nv
  other%neig = this%neig
  other%meig = this%meig

  other%nmat = this%nmat
  other%usize = this%usize

  other%tda = this%tda
  other%has_Avc = this%has_Avc
  other%has_flat_Avc = this%has_flat_Avc

  other%mf_header = this%mf_header
  other%has_mf_header = this%has_mf_header

  other%is_master = this%is_master
  other%is_distributed = this%is_distributed
  other%global_meig = this%global_meig

  other%is_open = this%is_open
  other%nwrite = this%nwrite
  other%ieig_read = this%ieig_read
  other%uur = this%uur
  other%uul = this%uul
  other%default_fname_bin = this%default_fname_bin
  other%default_fname_hdf5 = this%default_fname_hdf5
  other%use_hdf5 = this%use_hdf5

  call other%alloc(with_Avc=.true.)

  other%qpts = this%qpts
  if (allocated(this%kpts) .and. allocated(other%kpts)) other%kpts = this%kpts
  if (allocated(this%evals) .and. allocated(other%evals)) other%evals = this%evals
  if (allocated(this%u_r) .and. allocated(other%u_r)) other%u_r = this%u_r
  if (allocated(this%u_l) .and. allocated(other%u_l)) other%u_l = this%u_l
  if (allocated(this%Avc) .and. allocated(other%Avc)) other%Avc = this%Avc
  if (allocated(this%Avc_l) .and. allocated(other%Avc_l)) other%Avc_l = this%Avc_l
  if (allocated(this%Bvc_r) .and. allocated(other%Bvc_r)) other%Bvc_r = this%Bvc_r
  if (allocated(this%Bvc_l) .and. allocated(other%Bvc_l)) other%Bvc_l = this%Bvc_l

  if (allocated(this%global_ieigs) .and. allocated(other%global_ieigs)) other%global_ieigs = this%global_ieigs
  if (allocated(this%who_owns) .and. allocated(other%who_owns)) other%who_owns = this%who_owns
  if (allocated(this%local_ieigs) .and. allocated(other%local_ieigs)) other%local_ieigs = this%local_ieigs

  POP_SUB(evecs_copy)
end subroutine evecs_copy

! =========================================================================== !

!> Compare a evecs object against another, and report any discrepancy.
!! For testing purposes.
subroutine evecs_compare(this, other, ndiff, verbose)
  class(evecs_t), intent(in) :: this
  class(evecs_t), intent(in) :: other
  integer, intent(out) :: ndiff
  logical, intent(in), optional :: verbose

  logical :: verbose_ = .True.
  real(DP), parameter :: tol=1.0d-5
  real(DP) :: checksum
  integer :: unt=6

  if (present(verbose)) verbose_ = verbose

  PUSH_SUB(evecs_compare)

  ndiff = 0

  checksum = abs(sum(this%qpts) - sum(other%qpts))
  if (checksum .gt. tol) then
    if (verbose_) write(unt,*) 'Error in checksum of: qpt. ','delta=',checksum
    ndiff = ndiff + 1
  end if

  if (allocated(this%kpts) .and. allocated(other%kpts)) then
    checksum = abs(sum(this%kpts) - sum(other%kpts))
    if (checksum .gt. tol) then
      if (verbose_) write(unt,*) 'Error in checksum of: kpts. ','delta=',checksum
      ndiff = ndiff + 1
    end if
  end if

  if (allocated(this%Avc) .and. allocated(other%Avc)) then
    checksum = abs(sum(this%Avc) - sum(other%Avc))
    if (checksum .gt. tol) then
      if (verbose_) write(unt,*) 'Error in checksum of: Avc. ','delta=',checksum
      ndiff = ndiff + 1
    end if
  end if

  if (allocated(this%Avc_l) .and. allocated(other%Avc_l)) then
    checksum = abs(sum(this%Avc_l) - sum(other%Avc_l))
    if (checksum .gt. tol) then
      if (verbose_) write(unt,*) 'Error in checksum of: Avc_l. ','delta=',checksum
      ndiff = ndiff + 1
    end if
  end if

  if (allocated(this%Bvc_r) .and. allocated(other%Bvc_r)) then
    checksum = abs(sum(this%Bvc_r) - sum(other%Bvc_r))
    if (checksum .gt. tol) then
      if (verbose_) write(unt,*) 'Error in checksum of: Bvc_r. ','delta=',checksum
      ndiff = ndiff + 1
    end if
  end if

  if (allocated(this%Bvc_l) .and. allocated(other%Bvc_l)) then
    checksum = abs(sum(this%Bvc_l) - sum(other%Bvc_l))
    if (checksum .gt. tol) then
      if (verbose_) write(unt,*) 'Error in checksum of: Bvc_l. ','delta=',checksum
      ndiff = ndiff + 1
    end if
  end if

  POP_SUB(evecs_compare)

end subroutine evecs_compare

! =========================================================================== !

!> Pass the MPI communicator
subroutine evecs_set_mpi_comm(this, comm)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: comm  !< MPI communicator

  PUSH_SUB(evecs_set_mpi_comm)

#ifdef MPI
  this%comm = MPI_COMM_WORLD
  this%mpi_info = MPI_INFO_NULL
  if (present(comm)) then
    this%comm = comm
  end if
  call MPI_Comm_rank(this%comm, this%inode, mpierr)
  call MPI_Comm_size(this%comm, this%npes, mpierr)
#endif

  this%is_master = (this%inode .eq. 0)

  POP_SUB(evecs_set_mpi_comm)

end subroutine evecs_set_mpi_comm

! =========================================================================== !

!> Setup the default parallel scheme, where each node owns a copy
!! of all eigenvectors. No information on the dimensions is needed.
subroutine evecs_setup_paral_default(this, comm, neig)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: comm  !< MPI communicator
  integer, intent(in), optional :: neig  !< Number of eigenvalues

  integer :: ieig

  PUSH_SUB(evecs_setup_paral_default)

  this%is_distributed = .false.

  call this%set_mpi_comm(comm)

  this%global_meig = this%meig
  if (present(neig)) this%neig = neig
  ! No distribution: Master owns everything (but evecs might be broadcast)
  if (this%neig .gt. 0) then
    if (allocated(this%global_ieigs)) then
      SAFE_DEALLOCATE(this%global_ieigs)
    end if

    if (allocated(this%local_ieigs)) then
      SAFE_DEALLOCATE(this%local_ieigs)
    end if

    if (allocated(this%who_owns)) then
      SAFE_DEALLOCATE(this%who_owns)
    end if

    SAFE_ALLOCATE(this%global_ieigs, (this%neig))
    SAFE_ALLOCATE(this%local_ieigs, (this%neig))
    SAFE_ALLOCATE(this%who_owns, (this%neig))

    this%who_owns = 0  !< The index of the node holding each eigenvector
    this%local_ieigs = 0  !< The local index on the node of all eigenvector
    this%global_ieigs = 0  !< The global index of each local eigenvector

    do ieig=1,this%neig
      this%local_ieigs(ieig) = ieig
      this%global_ieigs(ieig) = ieig
    enddo
    
  end if

  POP_SUB(evecs_setup_paral_default)

end subroutine evecs_setup_paral_default

! =========================================================================== !

!> Initialize parallel info by relying on peinfo, which has been set up for BSE
!! Assume that neig is known
subroutine evecs_get_BSE_paral(this)
  class(evecs_t), intent(inout) :: this

  integer :: ieig_local, ieig_global, peadd

  PUSH_SUB(evecs_get_BSE_paral)

  call this%setup_paral_default()

  if (allocated(peinf%neig)) then
    this%meig = peinf%neig(peinf%inode+1)
  end if
  this%global_meig = this%meig
#ifdef MPI
  call MPI_REDUCE(this%meig, this%global_meig, 1, MPI_INT, MPI_MAX, 0, this%comm, mpierr)
  call MPI_BCAST(this%global_meig, 1, MPI_INTEGER, 0, this%comm, mpierr)
#endif

  if (.not. allocated(this%global_ieigs)) then
    SAFE_ALLOCATE(this%global_ieigs, (this%global_meig))
  end if

  if (.not. allocated(this%local_ieigs)) then
    SAFE_ALLOCATE(this%local_ieigs, (this%neig))
  end if

  if (.not. allocated(this%who_owns)) then
    SAFE_ALLOCATE(this%who_owns, (this%neig))
  end if

  this%who_owns = 0  !< The index of the node holding each eigenvector
  this%local_ieigs = 0  !< The local index on the node of all eigenvector
  this%global_ieigs = 0  !< The global index of each local eigenvector

  if (allocated(peinf%peig)) then

    this%is_distributed = .true.

    do ieig_local=1,this%meig
      this%global_ieigs(ieig_local) = peinf%peig(peinf%inode+1, ieig_local)
    end do

    do peadd=1,peinf%npes
      do ieig_local=1,peinf%neig(peadd)
        ieig_global = peinf%peig(peadd, ieig_local)
        this%who_owns(ieig_global) = peadd - 1
        this%local_ieigs(ieig_global) = ieig_local
      end do
    end do

  else

    this%is_distributed = .false.

    do ieig_local=1,this%meig
      this%global_ieigs(ieig_local) = ieig_local
      this%local_ieigs(ieig_local) = ieig_local
    end do

  end if

  POP_SUB(evecs_get_BSE_paral)

end subroutine evecs_get_BSE_paral

! =========================================================================== !

!> Setup the distribution of eigenvectors among all nodes.
!! The number of eigenvalues must be known prior to calling this function
!! or be given explicitly.
!! If there are N nodes in the group numbered from 0 to N-1
!! and there are M eigenvalues numbered from 1 to M,
!! then eigenvalue i will be owned by node mod(i-1, N)
subroutine evecs_setup_paral_distribution(this, neig, comm)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: neig  !< Number of eigenvectors
  integer, intent(in), optional :: comm  !< MPI communicator

  integer :: ieig, ieig_local
  integer, allocatable :: tmp_local_ieigs(:)

  PUSH_SUB(evecs_setup_paral_distribution)

  call this%set_mpi_comm(comm=comm)

  if (present(neig)) then
    this%neig = neig
  else if (this%neig .eq. 0) then
    call die('evecs_setup_paral_distribution was called without specifying neig, but neig is unknown.')
  end if

  if (allocated(this%global_ieigs)) then
    SAFE_DEALLOCATE(this%global_ieigs)
  end if

  if (allocated(this%local_ieigs)) then
    SAFE_DEALLOCATE(this%local_ieigs)
  end if

  if (allocated(this%who_owns)) then
    SAFE_DEALLOCATE(this%who_owns)
  end if

  ! This array only needs to be of length meig
  ! but we are using neig as an upper bound.
  SAFE_ALLOCATE(this%global_ieigs, (this%neig))
  SAFE_ALLOCATE(this%local_ieigs, (this%neig))
  SAFE_ALLOCATE(this%who_owns, (this%neig))

  this%who_owns = 0  !< The index of the node holding each eigenvector
  this%local_ieigs = 0  !< The local index on the node of all eigenvector
  this%global_ieigs = 0  !< The global index of each local eigenvector

  ! Manual distribution of the eigenvectors
  ieig_local = 0
  do ieig=1,this%neig
    this%who_owns(ieig) = mod(ieig-1, this%npes)
    if (this%who_owns(ieig) .eq. this%inode) then
      ieig_local = ieig_local + 1
      this%local_ieigs(ieig) = ieig_local
      this%global_ieigs(ieig_local) = ieig
    endif
  enddo

  this%meig = ieig_local

  ! Figure out maximum number of eigenvectors locally stored
  this%global_meig = this%meig
#ifdef MPI
  call MPI_REDUCE(this%meig, this%global_meig, 1, MPI_INT, MPI_MAX, 0, this%comm, mpierr)
  call MPI_BCAST(this%global_meig, 1, MPI_INTEGER, 0, this%comm, mpierr)
#endif

  ! Communicate all local indices
#ifdef MPI
  SAFE_ALLOCATE(tmp_local_ieigs, (this%neig))
  call MPI_ALLREDUCE(this%local_ieigs, tmp_local_ieigs, this%neig, &
                     MPI_INTEGER, MPI_SUM, this%comm, mpierr)
  this%local_ieigs = tmp_local_ieigs
  SAFE_DEALLOCATE(tmp_local_ieigs)
#endif

  this%is_distributed = .true.

  POP_SUB(evecs_setup_paral_distribution)

end subroutine evecs_setup_paral_distribution

! =========================================================================== !

!> Allocate memory
subroutine evecs_alloc(this, with_Avc, with_flat_Avc)
  class(evecs_t), intent(inout) :: this
  logical, intent(in), optional :: with_Avc  !< Allocate e-h coefficients
  logical, intent(in), optional :: with_flat_Avc  !< Allocate u arrays

  PUSH_SUB(evecs_alloc)

  call this%alloc_global_arrays()
  call this%alloc_local_arrays(with_Avc=with_Avc, with_flat_Avc=with_flat_Avc)

  POP_SUB(evecs_alloc)

end subroutine evecs_alloc


!> Allocate memory
subroutine evecs_alloc_global_arrays(this)
  class(evecs_t), intent(inout) :: this

  PUSH_SUB(evecs_alloc_global_arrays)


  ! Global arrays
  ! -------------
  ! k-points
  if (.not. allocated(this%kpts)) then
    SAFE_ALLOCATE(this%kpts, (3, this%nk))
  end if

  ! Eigenvalues
  if (.not. allocated(this%evals)) then
    SAFE_ALLOCATE(this%evals, (this%neig))
  end if

  if (.not. allocated(this%global_ieigs)) then
    SAFE_ALLOCATE(this%global_ieigs, (this%meig))
  end if

  ! q-points
  if (.not. allocated(this%qpts)) then
    SAFE_ALLOCATE(this%qpts, (3, this%nq))
    this%qpts = 0.0D0
  end if

  POP_SUB(evecs_alloc_global_arrays)

end subroutine evecs_alloc_global_arrays

! =========================================================================== !

!> Allocate memory
subroutine evecs_alloc_local_arrays(this, with_Avc, with_flat_Avc)
  class(evecs_t), intent(inout) :: this
  logical, intent(in), optional :: with_Avc  !< Whether Avc arrays should be allocated
  logical, intent(in), optional :: with_flat_Avc  !< Whether u arrays should be allocated

  logical :: with_Avc_, with_flat_Avc_  !< Whether u arrays should be allocated

  PUSH_SUB(evecs_alloc_local_arrays)

  with_Avc_ = .false.
  with_flat_Avc_ = .true.
  if (present(with_Avc)) then
    with_Avc_ = with_Avc
  end if
  if (present(with_flat_Avc)) then
    with_flat_Avc_ = with_flat_Avc
  end if

  ! Flat Avc
  if ( this%tda ) then
    this%usize = this%nmat
    if (with_flat_Avc_) then
      SAFE_ALLOCATE(this%u_r, (this%usize,this%meig))
    end if
  else
    this%usize = 2*this%nmat
    if (with_flat_Avc_) then
      SAFE_ALLOCATE(this%u_r, (this%usize, this%meig))
      SAFE_ALLOCATE(this%u_l, (this%usize, this%meig))
    end if
  end if

  ! Avc
  if (with_Avc_) then
    if ( this%tda ) then
      SAFE_ALLOCATE(this%Avc, (this%ns,this%nv,this%nc,this%nk,this%meig))
    else
      SAFE_ALLOCATE(this%Avc, (this%ns,this%nv,this%nc,this%nk, this%meig))
      SAFE_ALLOCATE(this%Avc_l, (this%ns,this%nv,this%nc,this%nk, this%meig))
      SAFE_ALLOCATE(this%Bvc_r, (this%ns,this%nv,this%nc,this%nk, this%meig))
      SAFE_ALLOCATE(this%Bvc_l, (this%ns,this%nv,this%nc,this%nk, this%meig))
    end if
  end if

  POP_SUB(evecs_alloc_local_arrays)

end subroutine evecs_alloc_local_arrays

! =========================================================================== !

!> Free memory
subroutine evecs_free(this)
  class(evecs_t), intent(inout) :: this

  PUSH_SUB(evecs_free)

  SAFE_DEALLOCATE(this%kpts)
  SAFE_DEALLOCATE(this%qpts)
  SAFE_DEALLOCATE(this%evals)
  SAFE_DEALLOCATE(this%u_r)
  SAFE_DEALLOCATE(this%u_l)
  SAFE_DEALLOCATE(this%Avc)
  SAFE_DEALLOCATE(this%Avc_l)
  SAFE_DEALLOCATE(this%Bvc_r)
  SAFE_DEALLOCATE(this%Bvc_l)
  SAFE_DEALLOCATE(this%global_ieigs)
  SAFE_DEALLOCATE(this%local_ieigs)
  SAFE_DEALLOCATE(this%who_owns)

  ! Reinitialize some variables
  this%has_Avc = .false.
  this%has_flat_Avc = .false.
  this%is_distributed = .false.
  this%ieig_read = 0
  this%nwrite = 0

  POP_SUB(evecs_free)

end subroutine evecs_free

! =========================================================================== !

!> Transform the flat eigenvectors into the Avc components
subroutine evecs_reshape_Avc(this, force)
  class(evecs_t), intent(inout) :: this
  logical, intent(in), optional :: force  !< Do it even if done before

  integer :: ieig,iflat,ik,ic,iv,is
  logical :: force_  !< Do it even if done before

  PUSH_SUB(evecs_reshape_Avc)

  if (this%has_Avc) then
    force_ = .false.
    if (present(force)) then
      force_ = force
    end if
    if (.not. force_) return
  end if

  ! local worker may not even have the data
  if (.not. allocated(this%u_r)) then
    this%has_Avc = .true.
    return
  end if

  if (.not. allocated(this%Avc)) then
    call this%alloc_local_arrays(with_Avc=.true., with_flat_Avc=.false.)
  end if

  if ( this%tda ) then
    do ieig=1,this%meig
      iflat = 0
      do ik=1,this%nk
        do ic=1,this%nc
          do iv=1,this%nv
            do is=1,this%ns
              iflat = iflat + 1
              this%Avc(is,iv,ic,ik,ieig) = this%u_r(iflat,ieig)
            enddo
          enddo
        enddo
      enddo
    enddo
  else
    do ieig=1,this%meig
      iflat = 0
      do ik=1,this%nk
        do ic=1,this%nc
          do iv=1,this%nv
            do is=1,this%ns
              iflat = iflat + 1
              this%Avc(is,iv,ic,ik,ieig) = this%u_r(iflat,ieig)
              this%Bvc_r(is,iv,ic,ik,ieig) = this%u_r(this%nmat+iflat,ieig)
              this%Avc_l(is,iv,ic,ik,ieig) = this%u_l(iflat,ieig)
              this%Bvc_l(is,iv,ic,ik,ieig) = this%u_l(this%nmat+iflat,ieig)
            enddo
          enddo
        enddo
      enddo
    enddo
  end if !tda

  this%has_Avc = .true.

  POP_SUB(evecs_reshape_Avc)

end subroutine evecs_reshape_Avc

! =========================================================================== !

!> Broadcast dimensions and k-points
subroutine evecs_broadcast_header(this, comm)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: comm  !< MPI communicator

  integer :: dims(11)
  integer :: tda

  PUSH_SUB(evecs_broadcast_header)

  call this%set_mpi_comm(comm)

  if (this%is_master) then
    if (this%tda) then
      tda = 1
    else
      tda = 0
    end if
    dims(:) = [this%nmat, this%usize, this%neig, this%ns, this%nk, this%nc, this%nv, this%meig, this%nq, this%krnl, tda]
  end if

  ! Broadcast
#ifdef MPI
  call MPI_BCAST(dims, 11, MPI_INTEGER, 0, this%comm, mpierr)
  call MPI_BARRIER(this%comm, mpierr)
#endif

  if (.not. this%is_master) then
    this%nmat = dims(1)
    this%usize = dims(2)
    this%neig = dims(3)
    this%ns = dims(4)
    this%nk = dims(5)
    this%nc = dims(6)
    this%nv = dims(7)
    this%meig = dims(8)
    this%nq = dims(9)
    this%krnl = dims(10)
    tda = dims(11)
    if (tda.eq.0) then
      this%tda = .false.
    else
      this%tda = .true.
    end if
  end if

  POP_SUB(evecs_broadcast_header)

end subroutine evecs_broadcast_header

! =========================================================================== !

!> Broadcast the arrays that should be owned by all processors
!! Assume that the arrays are already allocated
subroutine evecs_broadcast_global_arrays(this)
  class(evecs_t), intent(inout) :: this

  PUSH_SUB(evecs_broadcast_global_arrays)

#ifdef MPI
  call MPI_BARRIER(this%comm,mpierr)
  call MPI_BCAST(this%qpts, 3*this%nq, MPI_REAL_DP, 0, this%comm, mpierr)
  call MPI_BCAST(this%kpts, 3*this%nk, MPI_REAL_DP, 0, this%comm, mpierr)
  call MPI_BCAST(this%evals, this%neig, MPI_REAL_DP, 0, this%comm, mpierr)
  call MPI_BARRIER(this%comm,mpierr)
#endif

  POP_SUB(evecs_broadcast_global_arrays)

end subroutine evecs_broadcast_global_arrays

! =========================================================================== !

!> Broadcast the eigenvectors
subroutine evecs_broadcast_local_arrays(this)
  class(evecs_t), intent(inout) :: this

  integer :: dims

  PUSH_SUB(evecs_broadcast_local_arrays)

#ifdef MPI
  if (this%has_Avc) then
    dims = this%ns * this%nv * this%nc * this%nk * this%meig * this%nq
    call MPI_BCAST(this%Avc, dims , MPI_SCALAR, 0, this%comm, mpierr)
    if (.not. this%tda) then
      call MPI_BCAST(this%Avc_l, dims , MPI_SCALAR, 0, this%comm, mpierr)
      call MPI_BCAST(this%Bvc_r, dims , MPI_SCALAR, 0, this%comm, mpierr)
      call MPI_BCAST(this%Bvc_l, dims , MPI_SCALAR, 0, this%comm, mpierr)
    end if
  end if

  if (this%has_flat_Avc) then
    dims = this%usize * this%meig
    call MPI_BCAST(this%u_r, dims , MPI_SCALAR, 0, this%comm, mpierr)
    if (.not. this%tda) then
      call MPI_BCAST(this%u_l, dims , MPI_SCALAR, 0, this%comm, mpierr)
    end if
  end if

  call MPI_BARRIER(this%comm, mpierr)
#endif

  this%is_distributed = .false.

  POP_SUB(evecs_broadcast_local_arrays)

end subroutine evecs_broadcast_local_arrays

! =========================================================================== !
! Interface routines for binary / HDF5 selection
! =========================================================================== !

!> Main reading routine
subroutine evecs_read(this, fname, neig_read, ieig_offset, distribute, master_only, comm)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname
  integer, intent(in), optional :: neig_read    !< Number of eigenvectors to be read
  integer, intent(in), optional :: ieig_offset  !< Number of eigenvectors to skip befor reading
                                                !< The first eigenvector index is thus ieig_offset+1
  logical, intent(in), optional :: distribute   !< Distribute eigenvectors among all processors
  logical, intent(in), optional :: master_only  !< Only master reads eivenvectors
  integer, intent(in), optional :: comm         !< MPI communicator

  PUSH_SUB(evecs_read)

  if (this%use_hdf5) then
#ifdef HDF5
    call this%evecs_read_hdf5(fname, neig_read=neig_read, ieig_offset=ieig_offset, &
      distribute=distribute, master_only=master_only, comm=comm)
#endif
  else
    call this%evecs_read_bin(fname=fname, neig_read=neig_read, ieig_offset=ieig_offset, &
      with_Avc=.true., distribute=distribute, master_only=master_only, comm=comm)
  end if

  POP_SUB(evecs_read)

end subroutine evecs_read

! =========================================================================== !

!> Open file for sequential reading
subroutine evecs_open_read(this, fname, file_id)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname  !< File name
#ifdef HDF5
  integer(HID_T), intent(out), optional :: file_id  !< File ID for sequential read
#else
  integer, intent(out), optional :: file_id  !< File ID for sequential read
#endif
  PUSH_SUB(evecs_open_read)
  if (this%use_hdf5) then
#ifdef HDF5
    call this%evecs_open_read_hdf5(fname, file_id)
#endif
  else
    call this%evecs_open_read_bin(fname)
  end if
  POP_SUB(evecs_open_read)
end subroutine evecs_open_read

! =========================================================================== !

!> Interface function for reading header
subroutine evecs_read_header(this, fname, file_id, broadcast)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname  !< File name
#ifdef HDF5
  integer(HID_T), intent(inout), optional :: file_id  !< File ID for HDF5 sequential read
#else
  integer, intent(inout), optional :: file_id  !< File ID for HDF5 sequential read
#endif
  logical, intent(in), optional :: broadcast  !< Broadcast after

  PUSH_SUB(evecs_read_header)

  if (this%use_hdf5) then
#ifdef HDF5
    call this%evecs_read_header_hdf5(fname, file_id)
    call this%evecs_read_global_arrays_hdf5(fname)
#endif
  else
    call this%evecs_read_header_bin()
  end if

  if (present(broadcast)) then
    if (broadcast) then
      call this%broadcast_header()
      call this%broadcast_global_arrays()
    end if
  end if

  POP_SUB(evecs_read_header)

end subroutine evecs_read_header

! =========================================================================== !

!> Close file after sequential reading
subroutine evecs_close(this, file_id)
  class(evecs_t), intent(inout) :: this
#ifdef HDF5
  integer(HID_T), intent(inout), optional :: file_id
#else
  integer, intent(in), optional :: file_id
#endif
  PUSH_SUB(evecs_close)
  if (this%use_hdf5) then
#ifdef HDF5
      call this%evecs_close_hdf5(file_id)
#endif
  else
    call this%evecs_close_bin()
  end if
  POP_SUB(evecs_close)
end subroutine evecs_close

! =========================================================================== !

!> Read a single eigenvalue and eigenvector in serial
!! and store it in position 1 of the evals and evecs arrays
!!
!! If one wants needs this function to be compatible with both binary and hdf5,
!! then one should make sure the file is open beforehand and file_id is passed.
!! Else, for a call that only needs to be compatible with hdf5
!! the file_id can be omitted and the file will be closed after
!!
!! This is the sequential reading mode, which minimizes memory requirements.
!! It is used in summariza_eigenvectors.
subroutine evecs_read_next_eigenvector(this, fname, file_id)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname  !< File name
#ifdef HDF5
  integer(HID_T), intent(inout), optional :: file_id  !< File ID for HDF5 sequential read
#else
  integer, intent(inout), optional :: file_id  !< File ID for HDF5 sequential read
#endif

  PUSH_SUB(evecs_read_next_eigenvector)
  if (this%use_hdf5) then
#ifdef HDF5
    call this%evecs_read_local_arrays_ser_hdf5(fname=fname, file_id=file_id, neig_read=1, ieig_offset=this%ieig_read)
    this%ieig_read = this%ieig_read + 1
    this%evals(1) = this%evals(this%ieig_read)
#endif
  else
    call this%evecs_read_single_eigenvector_bin()
  end if

  POP_SUB(evecs_read_next_eigenvector)

end subroutine evecs_read_next_eigenvector


! =========================================================================== !
! Binary reading routines
! =========================================================================== !

!> Read eigenvectors in binary format
subroutine evecs_read_bin(this, fname, neig_read, ieig_offset, &
                          with_Avc, distribute, master_only, comm)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname
  integer, intent(in), optional :: neig_read    !< Number of eigenvectors to be read
                                                !! It is optional only if HDF5 is used
  integer, intent(in), optional :: ieig_offset  !< number of eigenvectors to skip befor reading
                                                !! the first eigenvector index is thus ieig_offset+1
                                                !! will not be used if eigenvectors are distributed
  logical, intent(in), optional :: with_Avc     !< Reshape into Avc
  logical, intent(in), optional :: distribute   !< spread the eigenvectors among processors
  logical, intent(in), optional :: master_only  !< Only master reads and holds eigenvectors
  integer, intent(in), optional :: comm         !< MPI communicator

  character(len=512) :: fname_

  PUSH_SUB(evecs_read_bin)

  ! Setup the distribution if requested
  if (present(distribute)) then
    if (distribute) then
      this%is_distributed = distribute
    endif
  endif

  if (this%is_distributed) then
    if (present(ieig_offset)) then
      call die('Cannot read and distribute binary file with ieig_offset')
    end if
    call this%evecs_read_and_distribute_bin(fname, neig_read=neig_read, &
                                            with_Avc=.true., comm=comm)
  else
    call this%evecs_read_and_broadcast_bin(fname, neig_read=neig_read, ieig_offset=ieig_offset, &
                    with_Avc=.true., master_only=master_only, comm=comm)
  end if

  POP_SUB(evecs_read_bin)
end subroutine evecs_read_bin

! =========================================================================== !

!> Open the binary eigenvectors file(s) for reading
subroutine evecs_open_read_bin(this, rootname)
  class(evecs_t), intent(inout) :: this
  character(*), intent(in), optional :: rootname  !< File name, to be appended by _r or _l
                                                  !! for non-tda.

  character(60) :: fname

  PUSH_SUB(evecs_open_read_bin)

  if (present(rootname)) then
    fname = trim(rootname)
  else
    fname = this%default_fname_bin
  end if

  if (this%tda) then 
    call open_file(unit=this%uur,file=trim(fname),form='unformatted',status='old')
  else
    call open_file(unit=this%uul,file=trim(trim(fname)//'_l'),form='unformatted',status='old')
    call open_file(unit=this%uur,file=trim(trim(fname)//'_r'),form='unformatted',status='old')
  end if

  this%is_open = .true.

  POP_SUB(evecs_open_read_bin)

end subroutine evecs_open_read_bin

! =========================================================================== !

!> Read the dimensions from binary file(s) then count number of eigenvectors
subroutine evecs_read_header_bin(this, neig)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: neig  !< Provide the number of eigenvectors
                                         !! so it needs not be counted

  PUSH_SUB(evecs_read_header_bin)

  call this%evecs_read_dimensions_bin()
  if (present(neig)) then
    this%neig = neig
    this%meig = neig
  else
    call this%evecs_count_neig_and_rewind_bin()
  end if

  if (allocated(this%qpts)) then
    SAFE_DEALLOCATE(this%qpts)
  end if
  SAFE_ALLOCATE(this%qpts, (3,this%nq))
  this%qpts = 0.0D0

  POP_SUB(evecs_read_header_bin)
end subroutine evecs_read_header_bin

! =========================================================================== !

!> Read the dimensions from binary file(s)
subroutine evecs_read_dimensions_bin(this)
  class(evecs_t), intent(inout) :: this

  real(DP), allocatable :: kpts_l(:,:)
  integer :: ns_, nv_, nc_, nk_, nmat_
  integer :: ik

  PUSH_SUB(evecs_read_dimensions_bin)

  if (.not. this%is_open) then
    call this%evecs_open_read_bin()
  end if

  read(this%uur) this%ns
  read(this%uur) this%nv
  read(this%uur) this%nc
  read(this%uur) this%nk
  this%nmat = this%ns * this%nv * this%nc * this%nk
  this%neig = this%nmat  !< Temporary value
  this%meig = this%neig
  if ( this%tda ) then
    this%usize = this%nmat
  else
    this%usize = 2*this%nmat
  end if
  if (.not. this%tda ) then
    read(this%uul) ns_
    read(this%uul) nv_
    read(this%uul) nc_
    read(this%uul) nk_
    nmat_ = ns_*nv_*nc_*nk_
    if( this%nmat .ne. nmat_) then
      call die('Error found: nmat should be equal for left and right eigenvalue matrices')
    end if
  end if

  if (allocated(this%kpts)) then
    SAFE_DEALLOCATE(this%kpts)
  end if
  SAFE_ALLOCATE(this%kpts, (3,this%nk))
  read(this%uur) this%kpts(:,:)
  if( .not. this%tda ) then
    SAFE_ALLOCATE(kpts_l, (3,this%nk))
    read(this%uul) kpts_l(:,:)
!   Check:
    do ik=1,this%nk
      if( abs(this%kpts(1,ik)-kpts_l(1,ik)) > 1.d06) then
        if( abs(this%kpts(2,ik)-kpts_l(2,ik)) > 1.d06) then
          if( abs(this%kpts(3,ik)-kpts_l(3,ik)) > 1.d06) then
            call die('Inconsistency in k-points found in left and right eigenvalues')
          end if
        end if
      end if 
    end do
    SAFE_DEALLOCATE(kpts_l)
  end if

  POP_SUB(evecs_read_dimensions_bin)

end subroutine evecs_read_dimensions_bin

! =========================================================================== !

!> Read a single eigenvalue and eigenvector(s)
subroutine evecs_read_single_eigenvector_bin(this, ieig_evec, ieig_eval)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: ieig_evec  !< Index of the eigenvector for storage
  integer, intent(in), optional :: ieig_eval  !< Index of the eigenvalue for storage

  integer :: ieig_evec_, ieig_eval_
  real(DP) :: eval_l

  PUSH_SUB(evecs_read_single_eigenvector_bin)

  ieig_eval_ = 1
  ieig_evec_ = 1
  if (present(ieig_eval)) ieig_eval_ = ieig_eval
  if (present(ieig_evec)) ieig_evec_ = ieig_evec

  if ( this%tda ) then
     read(this%uur) this%evals(ieig_eval_)
     read(this%uur) this%u_r(:,ieig_evec_)
  else
     read(this%uur) this%evals(ieig_eval_)
     read(this%uur) this%u_r(:,ieig_evec_)
     read(this%uul) eval_l
     read(this%uul) this%u_l(:,ieig_evec_)

     if( abs(eval_l - this%evals(ieig_eval_))>1.d-12) then
        call die('Inconsistency in energies in left and right eigenfuncion files')
     end if
  end if

  this%has_Avc = .false.

  POP_SUB(evecs_read_single_eigenvector_bin)

end subroutine evecs_read_single_eigenvector_bin

! =========================================================================== !

!> Skip eigenvectors in the file, but read and store the eigenvalues.
subroutine evecs_skip_eigenvectors_bin(this, nskip)
  class(evecs_t), intent(inout) :: this
  integer, intent(in) :: nskip  !< Number of eigenvector to skip

  integer :: ii
  real(DP) :: eval_l

  PUSH_SUB(evecs_skip_eigenvectors_bin)

  if ( this%tda ) then
    do ii=1,nskip
      read(this%uur)
      read(this%uur)
    end do
  else
    do ii=1,nskip
      read(this%uur)
      read(this%uur)
      read(this%uul)
      read(this%uul)
    end do
  end if

  POP_SUB(evecs_skip_eigenvectors_bin)

end subroutine evecs_skip_eigenvectors_bin

! =========================================================================== !

!> Count number of eigenvectors in binary file then rewind
!! at the position after header.
subroutine evecs_count_neig_and_rewind_bin(this)
  class(evecs_t), intent(inout) :: this

  integer :: ii, neig
  integer :: iostatus
  real(dp) :: eval

  PUSH_SUB(evecs_count_neig_and_rewind_bin)

  neig = 0
  if ( this%tda ) then
    do ii=1,this%nmat
      read(this%uur, iostat=iostatus) eval
      if (iostatus .eq. 0) then
        neig = neig + 1
      else if (iostatus .lt. 0) then
        exit 
      else
        call die('Error counting eigenvectors')
      end if
      read(this%uur)
    end do
  else
    do ii=1,this%nmat
      read(this%uur, iostat=iostatus) eval
      if (iostatus .eq. 0) then
        neig = neig + 1
      else if (iostatus .lt. 0) then
        exit 
      else
        call die('Error counting eigenvectors')
      end if
      read(this%uur)
      read(this%uul)
      read(this%uul)
    end do
  end if

  ! Now rewind the file
  rewind(this%uur)
  if (.not. this%tda) then
    rewind(this%uul)
  end if

  call this%evecs_read_dimensions_bin()
  this%neig = neig
  this%meig = neig

  POP_SUB(evecs_count_neig_and_rewind_bin)

end subroutine evecs_count_neig_and_rewind_bin

! =========================================================================== !

!> Binary reading of all eigenvectors in serial and broadcast all eigenvectors
!! to every node.
subroutine evecs_read_and_broadcast_bin(this, rootname, neig_read, ieig_offset, with_Avc, master_only, comm)
  class(evecs_t), intent(inout) :: this
  character(*), intent(in), optional :: rootname
  integer, intent(in), optional :: neig_read    !< Number of eigenvectors to be read
  integer, intent(in), optional :: ieig_offset  !< Number of eigenvectors to skip befor reading
                                                !< The first eigenvector index is thus ieig_offset+1
  logical, intent(in), optional :: with_Avc
  logical, intent(in), optional :: master_only  !< Only master reads eivenvectors
  integer, intent(in), optional :: comm         !< MPI communicator

  logical :: master_only_
  logical :: with_Avc_
  integer :: ieig
  integer :: nread, nskip

  PUSH_SUB(evecs_read_and_broadcast_bin)

  call this%set_mpi_comm(comm=comm)

  with_Avc_ = .true.
  if (present(with_Avc)) with_Avc_ = with_Avc

  if (present(master_only)) then
    master_only_ = master_only
  else
    master_only_ = .false.
  end if
  if (master_only_) this%is_distributed = .false.

  if (this%is_master) then

    call this%evecs_open_read_bin(rootname)
    call this%evecs_read_header_bin()

    nskip = 0
    if (present(ieig_offset)) nskip = ieig_offset

    if (present(neig_read)) then
      nread = neig_read
      this%neig = nread
      this%meig = nread
    else if (nskip > 0) then
      nread = this%neig - nskip
      this%neig = nread
      this%meig = nread
    else
      nread = this%neig
    end if

  end if

  call this%broadcast_header(comm=comm)

  call this%setup_paral_default(comm=comm)

  ! Allocate memory
  if (master_only_) then
    call this%alloc_global_arrays()
    if (this%is_master) then
      call this%alloc_local_arrays(with_Avc=with_Avc_)
    end if
  else
    call this%alloc(with_Avc=with_Avc_)
  end if

  if (this%is_master) then

    if (nskip > 0) then
      call this%evecs_skip_eigenvectors_bin(nskip)
    end if

    do ieig=1,nread
      call this%evecs_read_single_eigenvector_bin(ieig_evec=ieig, ieig_eval=ieig)
    end do

    call this%evecs_close_bin()

  end if

  this%has_flat_Avc = .true.

#ifdef MPI
  call this%broadcast_global_arrays()
  if (.not. master_only_) then
    call this%broadcast_local_arrays()
  end if
#endif

  if (with_Avc_) then
    call this%reshape_Avc()
  end if

  POP_SUB(evecs_read_and_broadcast_bin)

end subroutine evecs_read_and_broadcast_bin

! =========================================================================== !

!> Binary reading of all eigenvectors in serial then distribute eigenvectors
!! among the nodes.
subroutine evecs_read_and_distribute_bin(this, rootname, neig_read, with_Avc, comm)
  class(evecs_t), intent(inout) :: this
  character(*), intent(in), optional :: rootname
  integer, intent(in), optional :: neig_read    !< Number of eigenvectors to be read
  logical, intent(in), optional :: with_Avc     !< Reshape into Avc
  integer, intent(in), optional :: comm         !< MPI communicator

  logical :: with_Avc_
  integer :: ieig, local_ieig
  integer :: dest, tag
  SCALAR, allocatable :: tmp_u_r(:,:), tmp_u_l(:,:)

  PUSH_SUB(evecs_read_and_distribute_bin)

  call this%set_mpi_comm(comm=comm)

  if (this%is_master) then
    call this%evecs_open_read_bin(rootname)
    call this%evecs_read_header_bin()
  end if

  call this%broadcast_header(comm=comm)

  call this%setup_paral_distribution(neig=neig_read, comm=comm)

  with_Avc_ = .true.
  if (present(with_Avc)) with_Avc_ = with_Avc

  call this%alloc(with_Avc=with_Avc_)

  ! Allocate temporary array for storage
  if (this%is_master) then
    SAFE_ALLOCATE(tmp_u_r, (this%usize, this%meig))
    if (.not. this%tda) then
      SAFE_ALLOCATE(tmp_u_l, (this%usize, this%meig))
    end if
  end if

  ! Read eigenvectors
  do ieig=1,this%neig

    if (this%is_master) then
      call this%evecs_read_single_eigenvector_bin(ieig_evec=1, ieig_eval=ieig)
    end if

    ! Communicate the eigenvector to the appropriate node
    dest = this%who_owns(ieig)
    local_ieig = this%local_ieigs(ieig)
    tag = 2000 + dest

    ! If destination is master, store in a temporary array
    if ((dest .eq. 0) .and. this%is_master) then
      tmp_u_r(:,local_ieig) = this%u_r(:,1)

    ! Otherwise, master sends the eigenvector
    else if (this%is_master) then
#ifdef MPI
      call MPI_SEND(this%u_r(:,1),this%usize, MPI_SCALAR, &                                                                            
                    dest, tag, this%comm, mpierr)
#endif

    ! node receives the eigenvector
    else if (this%inode .eq. dest) then

#ifdef MPI
      call MPI_RECV(this%u_r(:,local_ieig),this%usize, MPI_SCALAR, &
                    0, tag, this%comm ,mpistatus, mpierr)
#endif

    end if

    if (.not. this%tda) then
      tag = 3000 + dest
      ! If destination is master, store in a temporary array
      if (dest .eq. 0) then
        tmp_u_l(:,local_ieig) = this%u_l(:,1)

      ! Otherwise, master sends the eigenvector
      else if (this%is_master) then
#ifdef MPI
        call MPI_SEND(this%u_l(:,1),this%usize, MPI_SCALAR, &                                                                            
                      dest, tag, this%comm, mpierr)
#endif

      ! node receives the eigenvector
      else if (this%inode .eq. dest) then

#ifdef MPI
        call MPI_RECV(this%u_l(:,local_ieig),this%usize, MPI_SCALAR, &
                      0, tag, this%comm ,mpistatus, mpierr)
#endif

      end if
    end if

  end do

  if (this%is_master) then
    call this%evecs_close_bin()

    ! Copy eigenvector from temporary array
    do ieig=1,this%meig
      this%u_r(:,ieig) = tmp_u_r(:,ieig)
    end do
    SAFE_DEALLOCATE(tmp_u_r)
    if (.not. this%tda) then
      do ieig=1,this%meig
        this%u_l(:,ieig) = tmp_u_l(:,ieig)
      end do
      SAFE_DEALLOCATE(tmp_u_l)
    end if

  end if

  if (with_Avc_) then
    call this%reshape_Avc()
  end if

  call this%broadcast_global_arrays()

  this%has_flat_Avc = .true.

  POP_SUB(evecs_read_and_distribute_bin)

end subroutine evecs_read_and_distribute_bin

! =========================================================================== !

!> Close the binary eigenvectors file(s)
subroutine evecs_close_bin(this)
  class(evecs_t), intent(inout) :: this

  PUSH_SUB(evecs_close_bin)

  call close_file(this%uur)
  if ( .not. this%tda) call close_file(this%uul)
  this%is_open = .false.

  POP_SUB(evecs_close_bin)

end subroutine evecs_close_bin

! =========================================================================== !
! Writing routines
! =========================================================================== !

!> Generic routine to write BSE eigenvectors, whether binary or HDF5 
subroutine evecs_write(this, nwrite, fname)
  class(evecs_t), intent(inout) :: this
  integer, intent(inout), optional :: nwrite
  character(len=*), intent(in), optional :: fname

  PUSH_SUB(evecs_write)

  if (this%use_hdf5) then
#ifdef HDF5
    call this%evecs_write_hdf5(nwrite, fname=fname) 
#endif
  else
    call this%evecs_write_bin(nwrite, fname=fname) 
  end if

  POP_SUB(evecs_write)

end subroutine evecs_write

! =========================================================================== !

!> Write the binary eigenvectors file(s)
subroutine evecs_write_bin(this, nwrite, fname)
  class(evecs_t), intent(inout) :: this
  integer, intent(inout), optional :: nwrite  !< Number of eigenvectors to write
  character(len=*), intent(in), optional :: fname  !< File name

  integer :: nwrite_
  integer :: ii,jj,ik
  integer :: ieig, ieig_local
  integer :: who_owns
  integer :: rank_r, rank_l
  integer :: tag
  logical :: io_r, io_l
  character(len=512) :: fname_
  character(len=512) :: fname_r, fname_l
  SCALAR, allocatable :: u_r(:), u_l(:)

  PUSH_SUB(evecs_write_bin)

  if (present(fname)) then
    fname_ = fname
  else
    fname_ = this%default_fname_bin
  end if

! Who can do io?
! FHJ: we do i/o for the right and left eigenvectors using different MPI ranks.
! The hope is to get better load balance on a lustre FS. The optimal solution,
! however, is to use HDF5.
  rank_r = 0
  io_r = (this%inode==rank_r)
  rank_l = this%npes-1
  io_l = ((this%inode==rank_l) .and. (.not. this%tda))

  nwrite_ = this%neig
  if (present(nwrite)) then
    if (nwrite.lt.0) then
      nwrite = this%neig
    else
      nwrite_ = nwrite
    end if
  end if
  this%nwrite = nwrite_

  ! Temporary arrays for communication
  SAFE_ALLOCATE(u_r, (this%usize))
  if (.not. this%tda) then
    SAFE_ALLOCATE(u_l, (this%usize))
  end if

  ! Open the file we will write to and write header information
  if (io_r) then
    write(6,*)
    if (.not. this%tda) then
      !fname_r = "eigenvectors_r"
      fname_r = trim(fname_)//"_r"
      write(6,'(1x,a,i0,a)') 'Writing ',nwrite_, ' right eigenvectors to file "'//trim(fname_r)//'"'
    else
      !fname_r = "eigenvectors"
      fname_r = trim(fname_)
      write(6,'(1x,a,i0,a)') 'Writing ',nwrite_, ' eigenvectors to file "'//trim(fname_r)//'"'
    endif
    write(6,'(1x,a,i0)') 'Length of each vector: ', this%usize
    write(6,*)
    call open_file(unit=this%uur,file=trim(fname_r),form='unformatted',status='replace')
    write(this%uur) this%ns
    write(this%uur) this%nv
    write(this%uur) this%nc
    write(this%uur) this%nk
    write(this%uur) ((this%kpts(jj,ik),jj=1,3),ik=1,this%nk)
  endif
  if (io_l) then
    !fname_l = "eigenvectors_l"
    fname_r = trim(fname_)//"_l"
    write(6,*)
    write(6,'(1x,a,i0,a)') 'Writing ',nwrite_, ' left eigenvectors to file "eigenvectors_l"'
    write(6,'(1x,a,i0)') 'Length of each vector: ', this%usize
    write(6,*)
    call open_file(unit=this%uul,file=trim(fname_l),form='unformatted',status='replace')
    write(this%uul) this%ns
    write(this%uul) this%nv
    write(this%uul) this%nc
    write(this%uul) this%nk
    write(this%uul) ((this%kpts(jj,ik),jj=1,3),ik=1,this%nk)
  endif

  ! Loop over states to be written
  do ieig=1,nwrite_

    ! Figure out which processor and column state ieig belongs to

    if (this%is_distributed) then
      who_owns = this%who_owns(ieig)
      ieig_local = this%local_ieigs(ieig)
    else
      who_owns = 0
      ieig_local = ieig
    end if

    ! Get the coeffs for state ieig into A (on all processors)
    if (this%inode==who_owns) then
      u_r(:) = this%u_r(:,ieig_local)
#ifdef MPI
      if (who_owns/=rank_r) then
        tag = ieig
        call MPI_Send(u_r(1), this%usize, MPI_SCALAR, rank_r, tag, &
          this%comm, mpierr)
      endif
#endif
      if (.not. this%tda) then
        u_l(:) = this%u_l(:,ieig_local)
#ifdef MPI
        if (who_owns/=rank_l) then
          tag = ieig + nwrite_
          call MPI_Send(u_l(1), this%usize, MPI_SCALAR, rank_l, tag, &
            this%comm, mpierr)
        endif
#endif
      endif
    endif
#ifdef MPI
    if (io_r .and. who_owns/=rank_r) then
      tag = ieig
      call MPI_Recv(u_r(1), this%usize, MPI_SCALAR, who_owns, tag, &
        this%comm, MPI_STATUS_IGNORE, mpierr)
    endif
    if (io_l .and. who_owns/=rank_l) then
      tag = ieig + nwrite_
      call MPI_Recv(u_l(1), this%usize, MPI_SCALAR, who_owns, tag, &
        this%comm, MPI_STATUS_IGNORE, mpierr)
    endif
#endif

    ! Write to file
    if (io_r) then
      write(this%uur) this%evals(ieig)
      write(this%uur) (u_r(ii),ii=1,this%usize)
    endif
    if (io_l) then
      write(this%uul) this%evals(ieig)
      write(this%uul) (u_l(ii),ii=1,this%usize)
    endif

  end do
  
  if (io_r) call close_file(this%uur)
  if (io_l) call close_file(this%uul)

  SAFE_DEALLOCATE(u_r)
  SAFE_DEALLOCATE(u_l)
  
  POP_SUB(evecs_write_bin)
  
end subroutine evecs_write_bin

! =========================================================================== !
! HDF5 writing routines
! =========================================================================== !

#ifdef HDF5

!> Write the hdf5 eigenvectors file(s)
subroutine evecs_write_hdf5(this, nwrite, fname)
  class(evecs_t), intent(inout) :: this
  integer, intent(inout), optional :: nwrite
  character(len=*), intent(in), optional :: fname

  integer :: nwrite_
  character(len=512) :: fname_

  PUSH_SUB(evecs_write_hdf5)

  if (present(fname)) then
    fname_ = fname
  else
    fname_ = this%default_fname_hdf5
  end if

  nwrite_ = this%neig
  if (present(nwrite)) then
    if (nwrite .lt. 0) then
      nwrite = this%neig
    else
      nwrite_ = nwrite
    end if
  end if
  this%nwrite = nwrite_

  if (this%is_master) then

    ! Create the file and write dimensions
    call this%evecs_create_and_write_header_hdf5(fname_)

    ! Write the arrays that are not distributed
    call this%evecs_write_global_arrays_hdf5(fname_)

  end if

  ! Write the distributed eigenvectors
#ifdef MPI
  if (this%is_distributed) then
    call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    call this%evecs_write_local_arrays_par_hdf5(fname_)
    call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
  else
    if (this%is_master) then
      call this%evecs_write_local_arrays_ser_hdf5(fname_)
    end if
  end if
#else
  call this%evecs_write_local_arrays_ser_hdf5(fname_)
#endif
  
  POP_SUB(evecs_write_hdf5)

end subroutine evecs_write_hdf5

! =========================================================================== !

!> Create the file and write header dimensions.
subroutine evecs_create_and_write_header_hdf5(this, fname)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in) :: fname  !< File name

  integer(HID_T) :: file_id
  integer :: error
  integer, parameter :: rank=7
  integer :: dims(rank)

  PUSH_SUB(evecs_create_and_write_header_hdf5)

  ! Create a new file using default properties.
  call hdf5_open_file(trim(fname), 'w', file_id)

  ! Create the groups
  call hdf5_create_group(file_id, 'exciton_header')
  call hdf5_create_group(file_id, 'exciton_data')
  call hdf5_create_group(file_id, 'exciton_header/params')
  call hdf5_create_group(file_id, 'exciton_header/kpoints')

  ! Write dimensions
  call hdf5_write_int(file_id, 'exciton_header/versionnumber',  this%file_version)
  call hdf5_write_int(file_id, 'exciton_header/flavor',  SCALARSIZE)
  call hdf5_write_int(file_id, 'exciton_header/params/bse_hamiltonian_size',  this%nmat)
  call hdf5_write_int(file_id, 'exciton_header/params/evec_sz',  this%usize)
  call hdf5_write_int(file_id, 'exciton_header/params/nevecs',  this%nwrite)
  call hdf5_write_int(file_id, 'exciton_header/params/ns',  this%ns)
  call hdf5_write_int(file_id, 'exciton_header/params/nc',  this%nc)
  call hdf5_write_int(file_id, 'exciton_header/params/nv',  this%nv)
  call hdf5_write_int(file_id, 'exciton_header/params/spin_kernel',  this%krnl)
  call hdf5_write_int(file_id, 'exciton_header/kpoints/nk',  this%nk)
  call hdf5_write_int(file_id, 'exciton_header/kpoints/nQ',  this%nq)

  call hdf5_write_logical(file_id, 'exciton_header/params/use_tda',  this%tda)
  
  ! Create the datasets for various arrays
  call hdf5_create_dset(file_id, 'exciton_header/kpoints/kpts', H5T_NATIVE_DOUBLE, [3, this%nk])
  call hdf5_create_dset(file_id, 'exciton_header/kpoints/exciton_Q_shifts', H5T_NATIVE_DOUBLE, [3, this%nq])

  call hdf5_create_dset(file_id, 'exciton_data/eigenvalues', H5T_NATIVE_DOUBLE, [this%nwrite])

  dims = [SCALARSIZE, this%ns, this%nv, this%nc, this%nk, this%nwrite, this%nq]

  call hdf5_create_dset(file_id, 'exciton_data/eigenvectors', &
                        h5t_native_double, dims)

  call hdf5_create_dset(file_id, 'exciton_data/eigenvectors_left', &
                        h5t_native_double, dims)

  call hdf5_create_dset(file_id, 'exciton_data/eigenvectors_deexcitation', &
                        h5t_native_double, dims)

  call hdf5_create_dset(file_id, 'exciton_data/eigenvectors_deexcitation_left', &
                        h5t_native_double, dims)

  call hdf5_close_file(file_id)

  ! Write mean field header
  if (this%has_mf_header) then
    call setup_hdf5_mf_file(fname, create_file=.false.)
    call write_hdf5_mf_header(fname, this%mf_header)
  end if

  POP_SUB(evecs_create_and_write_header_hdf5)

end subroutine evecs_create_and_write_header_hdf5

! =========================================================================== !

!> Serial writing of data owned by all processors
subroutine evecs_write_global_arrays_hdf5(this, fname)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in) :: fname  !< File name

  integer(HID_T) :: file_id

  PUSH_SUB(evecs_write_global_arrays_hdf5)

  if (.not. this%is_master) return

  call hdf5_open_file(trim(fname), 'rw', file_id)

  ! Write arrays
  call hdf5_write_double_array(file_id, 'exciton_header/kpoints/kpts', &
                               [3, this%nk], this%kpts)
  call hdf5_write_double_array(file_id, 'exciton_header/kpoints/exciton_Q_shifts', &
                               [3, this%nq], this%qpts)
  call hdf5_write_double_array(file_id, 'exciton_data/eigenvalues', &
                               [this%nwrite], this%evals(1:this%nwrite))

  call hdf5_close_file(file_id)

  POP_SUB(evecs_write_global_arrays_hdf5)
  
end subroutine evecs_write_global_arrays_hdf5

! =========================================================================== !

!> Serial writing of the data
subroutine evecs_write_local_arrays_ser_hdf5(this, fname)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in) :: fname  !< File name

  integer, parameter :: rank=6
  integer :: dims(rank)
  integer :: is, iv, ic, ik, ieig
  integer :: error
  integer(HID_T) :: file_id       ! File identifier

  PUSH_SUB(evecs_write_local_arrays_ser_hdf5)

  if (.not. this%is_master) return

  if (.not. this%has_Avc) then
    call this%reshape_Avc()
  end if

  dims = [this%ns, this%nv, this%nc, this%nk, this%nwrite, this%nq]

  call hdf5_open_file(trim(fname), 'rw', file_id)

  call hdf5_write_scalar_array(file_id, 'exciton_data/eigenvectors', &
                               dims, this%Avc(:,:,:,:,1:this%nwrite))

  if (.not. this%tda) then

    call hdf5_write_scalar_array(file_id, 'exciton_data/eigenvectors_left', &
                                 dims, this%Avc_l(:,:,:,:,1:this%nwrite))

    call hdf5_write_scalar_array(file_id, 'exciton_data/eigenvectors_deexcitation', &
                                 dims, this%Bvc_r(:,:,:,:,1:this%nwrite))

    call hdf5_write_scalar_array(file_id, 'exciton_data/eigenvectors_deexcitation_left', &
                                 dims, this%Bvc_l(:,:,:,:,1:this%nwrite))

  end if

  call hdf5_close_file(file_id)

  POP_SUB(evecs_write_local_arrays_ser_hdf5)
  
end subroutine evecs_write_local_arrays_ser_hdf5

! =========================================================================== !

!> Parallel writing of the data
subroutine evecs_write_local_arrays_par_hdf5(this, fname)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in) :: fname  !< File name

  integer :: meig_owned, ieig_local, ieig_global
  integer :: error
  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: plist_id    ! Property list identifyer

  PUSH_SUB(evecs_write_local_arrays_par_hdf5)

  if (.not. this%has_Avc) then
    call this%reshape_Avc()
  end if

  ! Open file in parallel io mode
  call hdf5_open_file(trim(fname), 'rw', file_id, parallel_io=.true., comm=this%comm)

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
#ifdef MPI
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#endif

  ! Write data
  if (this%tda) then
    do ieig_local=1,this%global_meig
      ieig_global = this%global_ieigs(ieig_local)
      call this%evecs_write_single_eigenvector_par_hdf5(ieig_local, ieig_global, file_id, plist_id, 1)
    end do
  else
    do ieig_local=1,this%global_meig
      ieig_global = this%global_ieigs(ieig_local)
      call this%evecs_write_single_eigenvector_par_hdf5(ieig_local, ieig_global, file_id, plist_id, 1)
      call this%evecs_write_single_eigenvector_par_hdf5(ieig_local, ieig_global, file_id, plist_id, 2)
      call this%evecs_write_single_eigenvector_par_hdf5(ieig_local, ieig_global, file_id, plist_id, 3)
      call this%evecs_write_single_eigenvector_par_hdf5(ieig_local, ieig_global, file_id, plist_id, 4)
    end do
  end if

  ! Close things
  call h5pclose_f(plist_id, error)  ! Property list
  call hdf5_close_file(file_id)  ! Close file

  POP_SUB(evecs_write_local_arrays_par_hdf5)
  
end subroutine evecs_write_local_arrays_par_hdf5

! =========================================================================== !

!> Write a single eigenvector
subroutine evecs_write_single_eigenvector_par_hdf5(this, ieig_local, ieig_global, file_id, plist_id, which)
  class(evecs_t), intent(inout) :: this
  integer, intent(in) :: ieig_local  !< Index of the eigenvector in the local array
  integer, intent(in) :: ieig_global !< Index of the eigenvector in the global array
                                     !! (i.e. actual eigenvalue index)
  integer(HID_T), intent(inout) :: file_id     !< File identifier
  integer(HID_T), intent(in) :: plist_id    !< Property list identifyer
  integer, intent(in) :: which  !< Which kind of eigenvector to write
                                !! 1 = Avc
                                !! 2 = Avc_l
                                !! 3 = Bvc_r
                                !! 4 = Bvc_l

  integer :: is, ik, ic, iv, iq
  integer :: error
  integer, parameter :: rank=7
  real(DP), allocatable :: data(:,:,:,:,:,:,:)
  integer(HSIZE_T) :: dims(rank), offset(rank)
  character(len=100) :: dataspace_name
  logical :: dummywrite
  integer(HID_T) :: dset_id     ! Dataset identifier
  integer(HID_T) :: filespace   ! Dataspace identifier in file
  integer(HID_T) :: memspace    ! Dataspace identifier in mem
  integer(HID_T) :: dspace_id   ! Dataspace identifier in mem

  PUSH_SUB(evecs_write_single_eigenvector_par_hdf5)

  ! All nodes must participate in collective write action
  ! even if they have nothing to write.
  ! This is the condition for not having to actually write data
  dummywrite = (ieig_local > this%meig .or. ieig_global > this%nwrite .or. ieig_global < 1)

  ! =============================== !
  ! Set up the data to be written
  ! =============================== !

  if (dummywrite) then
    dims = 1
  else
    dims = [SCALARSIZE, this%ns, this%nv, this%nc, this%nk, 1, this%nq]
  end if
  offset(:) = 0
  if (.not. dummywrite) then
    offset(rank-1) = ieig_global - 1
  end if


  if (which == 1) then
    dataspace_name = 'exciton_data/eigenvectors'
  else if (which == 2) then
    dataspace_name = 'exciton_data/eigenvectors_left'
  else if (which == 3) then
    dataspace_name = 'exciton_data/eigenvectors_deexcitation'
  else if (which == 4) then
    dataspace_name = 'exciton_data/eigenvectors_deexcitation_left'
  end if

  ! Copy the matrix into the data buffer
  SAFE_ALLOCATE(data, (dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))

  if (dummywrite) then
    data = 0.0_dp
  else
    ! A more elegant way would be to use pointer to select the array, but hey
    if (which == 1) then
      iq = 1
      do ik = 1, this%nk
        do ic = 1, this%nc
          do iv = 1, this%nv
            do is = 1, this%ns
#ifdef CPLX
              data(1,is,iv,ic,ik,1,iq) = &
                dble(this%Avc(is,iv,ic,ik,ieig_local))
              data(2,is,iv,ic,ik,1,iq) = &
                IMAG(this%Avc(is,iv,ic,ik,ieig_local))
#else
              data(1,is,iv,ic,ik,1,iq) = &
                dble(this%Avc(is,iv,ic,ik,ieig_local))
#endif
            end do
          end do
        end do
      end do
    else if (which == 2) then
      iq = 1
      do ik = 1, this%nk
        do ic = 1, this%nc
          do iv = 1, this%nv
            do is = 1, this%ns
#ifdef CPLX
              data(1,is,iv,ic,ik,1,iq) = &
                dble(this%Avc_l(is,iv,ic,ik,ieig_local))
              data(2,is,iv,ic,ik,1,iq) = &
                IMAG(this%Avc_l(is,iv,ic,ik,ieig_local))
#else
              data(1,is,iv,ic,ik,1,iq) = &
                dble(this%Avc_l(is,iv,ic,ik,ieig_local))
#endif
            end do
          end do
        end do
      end do
    else if (which == 3) then
      iq = 1
      do ik = 1, this%nk
        do ic = 1, this%nc
          do iv = 1, this%nv
            do is = 1, this%ns
#ifdef CPLX
              data(1,is,iv,ic,ik,1,iq) = &
                dble(this%Bvc_r(is,iv,ic,ik,ieig_local))
              data(2,is,iv,ic,ik,1,iq) = &
                IMAG(this%Bvc_r(is,iv,ic,ik,ieig_local))
#else
              data(1,is,iv,ic,ik,1,iq) = &
                dble(this%Bvc_r(is,iv,ic,ik,ieig_local))
#endif
            end do
          end do
        end do
      end do
    else if (which == 4) then
      iq = 1
      do ik = 1, this%nk
        do ic = 1, this%nc
          do iv = 1, this%nv
            do is = 1, this%ns
#ifdef CPLX
              data(1,is,iv,ic,ik,1,iq) = &
                dble(this%Bvc_l(is,iv,ic,ik,ieig_local))
              data(2,is,iv,ic,ik,1,iq) = &
                IMAG(this%Bvc_l(is,iv,ic,ik,ieig_local))
#else
              data(1,is,iv,ic,ik,1,iq) = &
                dble(this%Bvc_l(is,iv,ic,ik,ieig_local))
#endif
            end do
          end do
        end do
      end do
    end if
  end if


  ! =============================== !
  ! Done setting up the data
  ! =============================== !

  ! Create the memspace array
  call h5screate_simple_f(rank, dims, memspace, error)
  if (dummywrite) then
    ! Specify that we are not using any part of the array
    call h5sselect_none_f(memspace, error)
  else
    ! Specify that we are going to use the full array
    call h5sselect_all_f(memspace, error)
  end if

  ! Open the dataspace
  call h5dopen_f(file_id, trim(dataspace_name), dset_id, error)

  ! Create the filespace
  call h5dget_space_f(dset_id, filespace, error)

  ! Specify which part of the file we are writing
  if (dummywrite) then
    call h5sselect_none_f(filespace, error)
  else
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
  end if

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

  ! Free memory
  SAFE_DEALLOCATE(data)

  POP_SUB(evecs_write_single_eigenvector_par_hdf5)
  
end subroutine evecs_write_single_eigenvector_par_hdf5

! =========================================================================== !
! HDF5 reading routines
! =========================================================================== !

!> Read the hdf5 eigenvectors file(s)
subroutine evecs_read_hdf5(this, fname, file_id, neig_read, ieig_offset, &
                           distribute, master_only, comm)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname  !< File name
  integer(HID_T), intent(inout), optional :: file_id !< File identifier
  integer, intent(in), optional :: neig_read    !< Number of eigenvectors to be read
  integer, intent(in), optional :: ieig_offset  !< number of eigenvectors to skip befor reading
                                                !! the first eigenvector index is thus ieig_offset+1
                                                !! will not be used if eigenvectors are distributed
  logical, intent(in), optional :: distribute   !< spread the eigenvectors among processors
  logical, intent(in), optional :: master_only  !< Only master reads eivenvectors
  integer, intent(in), optional :: comm         !< mpi communicator

  character(len=512) :: fname_

  logical :: master_only_

  PUSH_SUB(evecs_read_hdf5)

  ! Do we split the eigenvectors among processors?
  ! this%is_distributed = .false.
#ifdef MPI
  if (present(distribute)) then
    this%is_distributed = distribute
  end if
#endif

  if (present(master_only)) then
    master_only_ = master_only
  else
    master_only_ = .false.
  end if
  if (master_only_) this%is_distributed = .false.

  ! Read header
  if (this%is_master) then
    call this%evecs_read_header_hdf5(fname)
  end if

#ifdef MPI
  call this%broadcast_header(comm=comm)
#endif

  ! Setup the eigenvectors distribution
  if (this%is_distributed) then
    call this%setup_paral_distribution(neig=neig_read, comm=comm)
  else
    call this%setup_paral_default(neig=neig_read)
    this%meig = this%neig
    if (present(neig_read)) then
      this%meig = neig_read
    end if
  end if


  ! Allocate memory
  if (master_only_) then
    call this%alloc_global_arrays()
    if (this%is_master) then
      call this%alloc_local_arrays(with_Avc=.true., with_flat_Avc=.false.)
    end if
  else
    call this%alloc(with_Avc=.true., with_flat_Avc=.false.)
  end if

  this%has_Avc = .true.

  ! Read and broadcast global arrays
  if (this%is_master) then
    call this%evecs_read_global_arrays_hdf5(fname, ieig_offset=ieig_offset)
  end if

#ifdef MPI
  call this%broadcast_global_arrays()
#endif

  if (this%is_distributed) then
    ! Parallel reading of the eigenvectors
    call this%evecs_read_local_arrays_par_hdf5(fname)
  else
    ! Master reads the arrays
    if (this%is_master) then
      call this%evecs_read_local_arrays_ser_hdf5(fname, file_id, neig_read=neig_read, ieig_offset=ieig_offset)
    end if
#ifdef MPI
    if (.not. master_only_) then
      call this%broadcast_local_arrays()
    end if
#endif
  end if
  
  POP_SUB(evecs_read_hdf5)

end subroutine evecs_read_hdf5

! =========================================================================== !

!> Open eigenvector file for sequential reading in HDF5
subroutine evecs_open_read_hdf5(this, fname, file_id)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname   !< File name
  integer(HID_T), intent(out), optional :: file_id  !< File identifier

  integer(HID_T) :: file_id_
  integer :: error
  character(len=512) :: fname_

  PUSH_SUB(evecs_open_read_hdf5)

  if (present(fname)) then
    fname_ = trim(fname)
  else
    fname_ = trim(this%default_fname_hdf5)
  end if

  call hdf5_open_file(trim(fname_), 'r', file_id_)
  this%file_id = file_id_

  if (present(file_id)) then
    file_id = file_id_
  end if

  this%is_open = .true.

  POP_SUB(evecs_open_read_hdf5)

end subroutine evecs_open_read_hdf5

! =========================================================================== !

!> Close eigenvector file after sequential reading
subroutine evecs_close_hdf5(this, file_id)
  class(evecs_t), intent(inout) :: this
  integer(HID_T), intent(inout), optional :: file_id  !< File identifier

  integer :: error

  PUSH_SUB(evecs_close_hdf5)

  if (present(file_id)) then
    call hdf5_close_file(file_id)
  else if (this%is_open) then
    call hdf5_close_file(this%file_id)
  end if
  this%is_open = .false.

  POP_SUB(evecs_close_hdf5)

end subroutine evecs_close_hdf5

! =========================================================================== !

!> Read the dimensions
subroutine evecs_read_header_hdf5(this, fname, file_id)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname  !< File name
  integer(HID_T), intent(inout), optional :: file_id  !< File identifier

  character(len=512) :: fname_
  integer(HID_T) :: file_id_
  integer :: scalarsize_read
  logical :: tda_read, file_has_mf_header

  PUSH_SUB(evecs_read_header_hdf5)

  if (.not. this%is_master) return

  if (present(fname)) then
    fname_ = fname
  else
    fname_ = this%default_fname_hdf5
  end if

  if (present(file_id)) then
    file_id_ = file_id
  else
    file_id_ = this%file_id
  end if

  if (.not. this%is_open) then
    call hdf5_open_file(trim(fname_), 'r', file_id_)
  end if

  ! Read dimensions
  call hdf5_read_int(file_id_, 'exciton_header/flavor',  scalarsize_read)
  call hdf5_read_int(file_id_, 'exciton_header/params/bse_hamiltonian_size',  this%nmat)
  call hdf5_read_int(file_id_, 'exciton_header/params/evec_sz',  this%usize)
  call hdf5_read_int(file_id_, 'exciton_header/params/nevecs',  this%neig)
  call hdf5_read_int(file_id_, 'exciton_header/params/ns',  this%ns)
  call hdf5_read_int(file_id_, 'exciton_header/params/nc',  this%nc)
  call hdf5_read_int(file_id_, 'exciton_header/params/nv',  this%nv)
  call hdf5_read_int(file_id_, 'exciton_header/params/spin_kernel',  this%krnl)
  call hdf5_read_int(file_id_, 'exciton_header/kpoints/nk',  this%nk)
  call hdf5_read_int(file_id_, 'exciton_header/kpoints/nQ',  this%nq)

  call hdf5_read_logical(file_id_, 'exciton_header/params/use_tda',  tda_read)

  this%meig = this%neig
  this%tda = tda_read

  file_has_mf_header = hdf5_exists(file_id_, "mf_header")

  ! Close the file only if it wasnt open before.
  if (.not. this%is_open) then
    call hdf5_close_file(file_id_)
  end if

  if (file_has_mf_header) then
    call read_hdf5_mf_header(fname_, this%mf_header)
    this%has_mf_header = .true.
  end if

  POP_SUB(evecs_read_header_hdf5)

end subroutine evecs_read_header_hdf5


! =========================================================================== !

!> Read global arrays (owned by all processors)
subroutine evecs_read_global_arrays_hdf5(this, fname, ieig_offset)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname  !< File name
  integer, intent(in), optional :: ieig_offset  !< number of eigenvectors to skip befor reading
                                                !! the first eigenvector index is thus ieig_offset+1

  integer :: neig_read
  integer :: ieig_offset_
  integer(HID_T) :: file_id
  character(len=512) :: fname_
  real(DP), allocatable :: evals_read(:)

  PUSH_SUB(evecs_read_global_arrays_hdf5)

  if (.not. this%is_master) return

  call this%alloc_global_arrays()

  if (present(ieig_offset)) then
    ieig_offset_ = ieig_offset
  else
    ieig_offset_ = 0
  end if
  neig_read = this%neig + ieig_offset_
  SAFE_ALLOCATE(evals_read, (neig_read))

  if (present(fname)) then
    fname_ = fname
  else
    fname_ = this%default_fname_hdf5
  end if

  if (this%is_open) then
    file_id = this%file_id
  else
    call hdf5_open_file(trim(fname_), 'r', file_id)
  end if

  call hdf5_read_double_array(file_id, 'exciton_header/kpoints/kpts', &
                              [3, this%nk], this%kpts)
  call hdf5_read_double_array(file_id, 'exciton_header/kpoints/exciton_Q_shifts', &
                              [3, this%nq], this%qpts)
  call hdf5_read_double_array(file_id, 'exciton_data/eigenvalues', &
                              [neig_read], evals_read)

  this%evals(:) = evals_read(ieig_offset_+1:ieig_offset_+this%neig)

  SAFE_DEALLOCATE(evals_read)

  if (.not. this%is_open) then
    call hdf5_close_file(file_id)
  end if

  POP_SUB(evecs_read_global_arrays_hdf5)

end subroutine evecs_read_global_arrays_hdf5


! =========================================================================== !

!> Read the exciton eigenvectors in serial mode.
subroutine evecs_read_local_arrays_ser_hdf5(this, fname, file_id, neig_read, ieig_offset)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname  !< File name
  integer(HID_T), intent(inout), optional :: file_id  !< File identifier
  integer, intent(in), optional :: neig_read    !< Number of eigenvectors to be read
  integer, intent(in), optional :: ieig_offset  !< Number of eigenvectors to skip befor reading
                                                !< The first eigenvector index is thus ieig_offset+1

  integer, parameter :: rank = 6
  integer :: read_count(rank)
  integer :: offset(rank)
  integer(HID_T) :: file_id_
  integer :: error
  integer :: neig
  character(len=512) :: fname_

  PUSH_SUB(evecs_read_local_arrays_ser_hdf5)

  if (.not. this%is_master) return

  if (present(fname)) then
    fname_ = fname
  else
    fname_ = this%default_fname_hdf5
  end if

  if (present(file_id)) then
    file_id_ = file_id
  else
    file_id_ = this%file_id
  end if

  if (.not. this%is_open) then
    call hdf5_open_file(trim(fname_), 'r', file_id_)
  end if

  if (present(neig_read)) then
    neig = neig_read
  else
    neig = this%meig
  end if

  read_count = [this%ns, this%nv, this%nc, this%nk, neig, this%nq]
  offset = 0
  if (present(ieig_offset)) offset(rank-1) = ieig_offset

  call hdf5_read_scalar_hyperslab(file_id_, 'exciton_data/eigenvectors', &
                read_count, offset, this%Avc(:,:,:,:,1:neig), error)

  if (.not. this%tda) then

    call hdf5_read_scalar_hyperslab(file_id_, 'exciton_data/eigenvectors_left', &
                  read_count, offset, this%Avc_l(:,:,:,:,1:neig), error)

    call hdf5_read_scalar_hyperslab(file_id_, 'exciton_data/eigenvectors_deexcitation', &
                  read_count, offset, this%Bvc_r(:,:,:,:,1:neig), error)

    call hdf5_read_scalar_hyperslab(file_id_, 'exciton_data/eigenvectors_deexcitation_left', &
                  read_count, offset, this%Bvc_l(:,:,:,:,1:neig), error)

  end if

  ! Close the file if it wasnt open at the start.
  if (.not. this%is_open) then
    call hdf5_close_file(file_id_)
  end if

  this%has_Avc = .true.

  POP_SUB(evecs_read_local_arrays_ser_hdf5)

end subroutine evecs_read_local_arrays_ser_hdf5

! =========================================================================== !

!> Parallel reading of the eigenvectors
subroutine evecs_read_local_arrays_par_hdf5(this, fname)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname  !< File name

  integer :: ieig_local, ieig_global
  integer :: error
  character(len=512) :: fname_
  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: plist_id    ! Property list identifyer

  PUSH_SUB(evecs_read_local_arrays_par_hdf5)

  if (present(fname)) then
    fname_ = fname
  else
    fname_ = this%default_fname_hdf5
  end if


  ! Open file in parallel io mode
  call hdf5_open_file(trim(fname_), 'r', file_id, parallel_io=.true., comm=this%comm)

  ! Create property list for collective dataset read
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
#ifdef MPI
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#endif

  if (this%tda) then
    do ieig_local=1,this%global_meig
      ieig_global = this%global_ieigs(ieig_local)
      call this%evecs_read_single_eigenvector_par_hdf5(ieig_local, ieig_global, file_id, plist_id, 1)
    end do
  else
    do ieig_local=1,this%global_meig
      ieig_global = this%global_ieigs(ieig_local)
      call this%evecs_read_single_eigenvector_par_hdf5(ieig_local, ieig_global, file_id, plist_id, 1)
      call this%evecs_read_single_eigenvector_par_hdf5(ieig_local, ieig_global, file_id, plist_id, 2)
      call this%evecs_read_single_eigenvector_par_hdf5(ieig_local, ieig_global, file_id, plist_id, 3)
      call this%evecs_read_single_eigenvector_par_hdf5(ieig_local, ieig_global, file_id, plist_id, 4)
    end do
  end if

  ! Close things
  call h5pclose_f(plist_id, error)  ! Property list
  call hdf5_close_file(file_id)  ! Close file

  POP_SUB(evecs_read_local_arrays_par_hdf5)

end subroutine evecs_read_local_arrays_par_hdf5

! =========================================================================== !

!> Read a single eigenvector
subroutine evecs_read_single_eigenvector_par_hdf5(this, ieig_local, ieig_global, file_id, plist_id, which)
  class(evecs_t), intent(inout) :: this
  integer, intent(in) :: ieig_local  !< Index of the eigenvector in the local array
  integer, intent(in) :: ieig_global !< Index of the eigenvector in the global array
                                     !! (i.e. actual eigenvalue index)
  integer(HID_T), intent(inout) :: file_id   !< File identifier
  integer(HID_T), intent(in) :: plist_id  !< Property list identifyer
  integer, intent(in) :: which  !< Which kind of eigenvector to read
                                !! 1 = Avc
                                !! 2 = Avc_l
                                !! 3 = Bvc_r
                                !! 4 = Bvc_l

  integer :: iq, is, ik, ic, iv
  integer :: error
  integer, parameter :: rank=7
  real(DP), allocatable :: data(:,:,:,:,:,:,:)
  integer(HSIZE_T) :: dims(rank), offset(rank)
  character(len=100) :: dataspace_name
  logical :: dummyread
  integer(HID_T) :: dset_id     ! Dataset identifier
  integer(HID_T) :: filespace   ! Dataspace identifier in file
  integer(HID_T) :: memspace    ! Dataspace identifier in mem
  integer(HID_T) :: dspace_id   ! Dataspace identifier in mem

  PUSH_SUB(evecs_read_single_eigenvector_par_hdf5)

  ! All nodes must participate in collective read action
  ! even if they have nothing to read.
  ! This is the condition for not having to actually read data
  dummyread = (ieig_local > this%meig .or. ieig_global > this%neig .or. ieig_global < 1)

  ! =============================== !
  ! Set up data buffer
  ! =============================== !

  if (dummyread) then
    dims = 1
  else
    dims = [SCALARSIZE, this%ns, this%nv, this%nc, this%nk, 1, this%nq]
  end if
  offset(:) = 0
  if (.not. dummyread) then
    offset(rank-1) = ieig_global - 1
  end if

  SAFE_ALLOCATE(data, (dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))

  if (which == 1) then
    dataspace_name = 'exciton_data/eigenvectors'
  else if (which == 2) then
    dataspace_name = 'exciton_data/eigenvectors_left'
  else if (which == 3) then
    dataspace_name = 'exciton_data/eigenvectors_deexcitation'
  else if (which == 4) then
    dataspace_name = 'exciton_data/eigenvectors_deexcitation_left'
  end if

  ! =============================== !
  ! Done setting up data buffer
  ! =============================== !

  ! Create the memspace array
  call h5screate_simple_f(rank, dims, memspace, error)
  if (dummyread) then
    ! Specify that we are not using any part of the array
    call h5sselect_none_f(memspace, error)
  else
    ! Specify that we are going to use the full array
    call h5sselect_all_f(memspace, error)
  end if

  ! Open the dataspace
  call h5dopen_f(file_id, trim(dataspace_name), dset_id, error)

  ! Create the filespace
  call h5dget_space_f(dset_id, filespace, error)

  ! Specify which part of the file we are reading
  if (dummyread) then
    call h5sselect_none_f(filespace, error)
  else
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
  end if

  ! =============================== !
  ! Read the main data
  ! =============================== !

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error, &
                  memspace, filespace, xfer_prp=plist_id)

  if (error .ne. 0) then
    call die('Problem reading eigenvector')
  end if

  ! =============================== !
  ! Done reading the main data
  ! =============================== !

  ! Close
  call h5dclose_f(dset_id, error)  ! Dataset
  call h5sclose_f(memspace, error)  ! Memory space
  call h5sclose_f(filespace, error) ! File space


  ! Copy the buffer into Avc
  if (.not. dummyread) then
    ! A more elegant way would be to use pointer to select the array, but hey
    if (which == 1) then
      iq = 1
      do ik = 1, this%nk
        do ic = 1, this%nc
          do iv = 1, this%nv
            do is = 1, this%ns
#ifdef CPLX
              this%Avc(is,iv,ic,ik,ieig_local) = &
                CMPLX(data(1,is,iv,ic,ik,1,iq), data(2,is,iv,ic,ik,1,iq))
#else
              this%Avc(is,iv,ic,ik,ieig_local) = data(1,is,iv,ic,ik,1,iq)
#endif
            end do
          end do
        end do
      end do
    else if (which == 2) then
      iq = 1
      do ik = 1, this%nk
        do ic = 1, this%nc
          do iv = 1, this%nv
            do is = 1, this%ns
#ifdef CPLX
              this%Avc_l(is,iv,ic,ik,ieig_local) = &
                CMPLX(data(1,is,iv,ic,ik,1,iq), data(2,is,iv,ic,ik,1,iq))
#else
              this%Avc_l(is,iv,ic,ik,ieig_local) = data(1,is,iv,ic,ik,1,iq)
#endif
            end do
          end do
        end do
      end do
    else if (which == 3) then
      iq = 1
      do ik = 1, this%nk
        do ic = 1, this%nc
          do iv = 1, this%nv
            do is = 1, this%ns
#ifdef CPLX
              this%Bvc_r(is,iv,ic,ik,ieig_local) = &
                CMPLX(data(1,is,iv,ic,ik,1,iq), data(2,is,iv,ic,ik,1,iq))
#else
              this%Bvc_r(is,iv,ic,ik,ieig_local) = data(1,is,iv,ic,ik,1,iq)
#endif
            end do
          end do
        end do
      end do
    else if (which == 4) then
      iq = 1
      do ik = 1, this%nk
        do ic = 1, this%nc
          do iv = 1, this%nv
            do is = 1, this%ns
#ifdef CPLX
              this%Bvc_l(is,iv,ic,ik,ieig_local) = &
                CMPLX(data(1,is,iv,ic,ik,1,iq), data(2,is,iv,ic,ik,1,iq))
#else
              this%Bvc_l(is,iv,ic,ik,ieig_local) = data(1,is,iv,ic,ik,1,iq)
#endif
            end do
          end do
        end do
      end do
    end if
  end if

  ! Free memory
  SAFE_DEALLOCATE(data)

  POP_SUB(evecs_read_single_eigenvector_par_hdf5)
  
end subroutine evecs_read_single_eigenvector_par_hdf5

#endif

! =========================================================================== !
! Output routines
! =========================================================================== !

!> Print out the information read in the header
subroutine evecs_print_out_header_info(this, unt)
  class(evecs_t), intent(inout) :: this
  integer, intent(in) :: unt  !< The unit for writing

  integer :: ik

  PUSH_SUB(evecs_print_out_header_info)

  if (this%is_master) then
    ! Print out information found in file
    write(unt,'(a)')
    write(unt,'(a,4i5)') ' ns, nv, nc, nk = ',this%ns,this%nv,this%nc,this%nk
    write(unt,'(a,i8)') ' nmat = ',this%nmat
    write(unt,'(a)')
    write(unt,'(a)') 'kpoints follow:'
    do ik=1, this%nk
      write(unt,'(i5,3f10.5)') ik, this%kpts(:,ik)
    enddo
  end if

  POP_SUB(evecs_print_out_header_info)

end subroutine evecs_print_out_header_info

! =========================================================================== !

!> Print out current dimensions of the object
subroutine evecs_print_out_dimensions(this, unt)
  class(evecs_t), intent(inout) :: this
  integer, intent(in) :: unt  !< The unit for writing

  PUSH_SUB(evecs_print_out_dimensions)

  if (this%is_master) then
    write(unt,'(a)') ' ---------------------------------------'
    write(unt,'(a)') ' BSE eigenvectors dimensions'
    write(unt,'(a,4i5)') ' ns, nv, nc, nk = ',this%ns,this%nv,this%nc,this%nk
    write(unt,'(a,i8)') ' nmat = ',this%nmat
    write(unt,'(a,i8)') ' neig = ',this%neig
    write(unt,'(a,i8)') ' meig = ',this%meig
    write(unt,'(a)') ' ---------------------------------------'
    write(unt,'(a)')
  end if

  POP_SUB(evecs_print_out_dimensions)

end subroutine evecs_print_out_dimensions

! =========================================================================== !

end module evecs_m
