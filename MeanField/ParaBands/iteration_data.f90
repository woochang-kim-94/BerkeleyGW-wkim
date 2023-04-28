#include "f_defs.h"
#if defined MPI && !defined USESCALAPACK
  #error ScaLAPACK is required for MPI builds.
#endif

module iteration_data_m

  use global_m

  implicit none

  private

    type iteration_data_t
      logical :: is_working = .true. !< Whether I`m working
      logical :: is_working_leader = .true. ! Whether I`m working and the inode==0 of my pool
      integer :: working_group = -1
      integer :: working_comm = -1
      integer :: working_inode = -1
      integer :: leader_group = -1
      integer :: leader_comm = -1
      integer :: leader_inode = -1
    contains
      procedure :: setup_comm => iteration_data_setup_comm
      procedure :: free_comm => iteration_data_free_comm
    endtype iteration_data_t

  public :: iteration_data_t

contains

subroutine iteration_data_setup_comm(this, should_work, kp_inode)
  class(iteration_data_t), intent(inout) :: this
  logical, intent(in) :: should_work
  integer, intent(in) :: kp_inode

#ifdef MPI
  logical :: in_group(peinf%npes)
  integer :: ipe, ninclude, group_ranks(peinf%npes), world_group

  this%is_working = should_work
  group_ranks = 0
  ! Create MPI group for everyone that is_working
  call MPI_Allgather(this%is_working, 1, MPI_LOGICAL, in_group, 1, MPI_LOGICAL, &
    MPI_COMM_WORLD, mpierr)
  ninclude = 0
  do ipe = 0, peinf%npes-1
    if (in_group(ipe+1)) then
      ninclude = ninclude + 1
      group_ranks(ninclude) = ipe
    endif
  enddo
  if (this%is_working) then
    call MPI_Comm_group(MPI_COMM_WORLD, world_group, mpierr)
    call MPI_Group_incl(world_group, ninclude, group_ranks, this%working_group, mpierr)
    call MPI_Comm_create(MPI_COMM_WORLD, this%working_group, this%working_comm, mpierr)
    call MPI_Comm_rank(this%working_comm, this%working_inode, mpierr)
  else
    call MPI_Comm_Create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, this%working_comm, mpierr)
  endif

  this%is_working_leader = should_work .and. kp_inode==0
  group_ranks = 0
  ! Create MPI group for everyone that is_working and kp_inode==0
  call MPI_Allgather(this%is_working_leader, 1, MPI_LOGICAL, in_group, 1, MPI_LOGICAL, &
    MPI_COMM_WORLD, mpierr)
  ninclude = 0
  do ipe = 0, peinf%npes-1
    if (in_group(ipe+1)) then
      ninclude = ninclude + 1
      group_ranks(ninclude) = ipe
    endif
  enddo
  if (this%is_working_leader) then
    call MPI_Comm_group(MPI_COMM_WORLD, world_group, mpierr)
    call MPI_Group_incl(world_group, ninclude, group_ranks, this%leader_group, mpierr)
    call MPI_Comm_create(MPI_COMM_WORLD, this%leader_group, this%leader_comm, mpierr)
    call MPI_Comm_rank(this%leader_comm, this%leader_inode, mpierr)
  else
    call MPI_Comm_Create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, this%leader_comm, mpierr)
  endif
#endif

end subroutine iteration_data_setup_comm


subroutine iteration_data_free_comm(this)
  class(iteration_data_t), intent(inout) :: this

#ifdef MPI
  if (this%is_working) then
    call MPI_Group_free(this%working_group, mpierr)
    call MPI_Comm_free(this%working_comm, mpierr)
  endif
  this%working_group = -1
  this%working_comm = -1
  this%working_inode = -1
  if (this%is_working_leader) then
    call MPI_Group_free(this%leader_group, mpierr)
    call MPI_Comm_free(this%leader_comm, mpierr)
  endif
  this%leader_group = -1
  this%leader_comm = -1
  this%leader_inode = -1
#endif

end subroutine iteration_data_free_comm

end module iteration_data_m
