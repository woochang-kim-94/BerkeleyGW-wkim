#include "f_defs.h"

module kpoint_pool_m

  use global_m
  implicit none

  private

    type kpoint_pool_t
      !> BLACS context wherein all processors in this pool are distributed in
      !! a single column. There are no leftover processors.
      integer :: cntxt_c = -1
      !> BLACS context wherein all processors in this pool are distributed in
      !! a single row. There are no leftover processors.
      integer :: cntxt_r = -1
      !> BLACS context wherein processors in this pool are distributed in a 2D
      !! grid. There might be leftover processors that are not part of the context.
      integer :: cntxt_2d = -1
      integer :: npools = 1
      integer :: ipool = 0
      integer :: inode = 0
      integer :: npes = 0
      integer :: npes_diag = 0
      integer :: group = 0
      integer :: comm = 0
      integer :: iunit = 6
    contains
      procedure :: setup => kpoint_pool_setup
      procedure :: cleanup => kpoint_pool_cleanup
    endtype kpoint_pool_t

  public :: kpoint_pool_t


contains


!> Note: this routine must be called by ALL MPI PROCESSORS!
subroutine kpoint_pool_setup(this, nk)
  class(kpoint_pool_t), intent(inout) :: this
  integer, intent(in) :: nk

#ifdef MPI
  integer, allocatable :: ranks(:,:)
  integer :: world_group, ipool, ii, tmp_cntxt
  integer :: nprow, npcol

  PUSH_SUB(kpoint_pool_setup)

  if (this%npools<=0) this%npools = nk
  this%npes = peinf%npes / this%npools
  if (this%npes<1) then
    this%npools = peinf%npes
    this%npes = 1
    if (peinf%inode==0) &
      write(0,'(1x,a,i0)') 'WARNING: Overriding number_pools to ', this%npools
  endif

  ! Create one MPI group per k-point pool
  this%ipool = peinf%inode / this%npes
  SAFE_ALLOCATE(ranks, (this%npes,1))
  ranks(:,1) = (/ (ii, ii = 0, this%npes-1) /) + this%ipool*this%npes
  if (this%ipool < this%npools) then
    call MPI_Comm_Group(MPI_COMM_WORLD, world_group, mpierr)
    call MPI_Group_Incl(world_group, this%npes, ranks(:,1), this%group, mpierr)
    call MPI_Comm_Create(MPI_COMM_WORLD, this%group, this%comm, mpierr)
    call MPI_Comm_rank(this%comm, this%inode, mpierr)
  else
    call MPI_Comm_Create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, this%comm, mpierr)
    this%inode = -1
    this%ipool = -1
  endif
  SAFE_DEALLOCATE(ranks)

  ! Create 1D BLACS contexts. They have to be created globally!
  SAFE_ALLOCATE(ranks, (1,this%npes))
  do ipool = 0, this%npools-1
    ranks(1,:) = (/ (ii, ii=0, this%npes-1) /)
    ranks = ranks + ipool*this%npes
    call blacs_get(-1, 0, tmp_cntxt)
    call blacs_gridmap(tmp_cntxt, ranks, 1, 1, this%npes) ! This will modify cntxt
    if (ipool==this%ipool) this%cntxt_c = tmp_cntxt
  enddo
  SAFE_DEALLOCATE(ranks)
  SAFE_ALLOCATE(ranks, (this%npes,1))
  do ipool = 0, this%npools-1
    ranks(:,1) = (/ (ii, ii=0, this%npes-1) /)
    ranks = ranks + ipool*this%npes
    call blacs_get(-1, 0, tmp_cntxt)
    call blacs_gridmap(tmp_cntxt, ranks, this%npes, this%npes, 1) ! This will modify cntxt
    if (ipool==this%ipool) this%cntxt_r = tmp_cntxt
  enddo
  SAFE_DEALLOCATE(ranks)

  ! Create 2D BLACS contexts. They have to be created globally!
  nprow = int(sqrt(dble(this%npes)))
  npcol = this%npes/nprow
  this%npes_diag = nprow*npcol
  SAFE_ALLOCATE(ranks, (nprow, npcol))
  do ipool = 0, this%npools-1
    do ii = 0, nprow*npcol-1
      ranks(mod(ii,nprow)+1, ii/nprow+1) = ii
    enddo
    ranks = ranks + ipool*this%npes
    call blacs_get(-1, 0, tmp_cntxt)
    call blacs_gridmap(tmp_cntxt, ranks(1,1), nprow, nprow, npcol) ! This will modify cntxt
    if (ipool==this%ipool) this%cntxt_2d = tmp_cntxt
  enddo
  SAFE_DEALLOCATE(ranks)

  POP_SUB(kpoint_pool_setup)
#endif

end subroutine kpoint_pool_setup


subroutine kpoint_pool_cleanup(this)
  class(kpoint_pool_t), intent(inout) :: this

#ifdef MPI
  PUSH_SUB(kpoint_pool_cleanup)

  if (this%cntxt_c>=0) call blacs_gridexit(this%cntxt_c)
  if (this%cntxt_r>=0) call blacs_gridexit(this%cntxt_r)
  if (this%cntxt_2d>=0) call blacs_gridexit(this%cntxt_2d)
  this%cntxt_c = -1
  this%cntxt_r = -1
  this%cntxt_2d = -1

  POP_SUB(kpoint_pool_cleanup)
#endif

end subroutine kpoint_pool_cleanup

end module kpoint_pool_m
