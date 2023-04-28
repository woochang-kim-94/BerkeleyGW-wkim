! File bgw_mpi.f90 automatically created from bgw_mpi.f90p by mako_preprocess.py.
! Do not edit the resulting file (bgw_mpi.f90) directly!

!>=========================================================================
!!
!! Module:
!!
!! bgw_mpi_m        Originally by FHJ     Last Modified 03/2018 (FHJ)
!!
!!     Wrappers around MPI functions. All functions work both with
!!     MPI and serial runs. Right now, only MPI_Bcast is supported.
!!     TODO: implement wrappers for MPI_Send, MPI_Recv, gather, etc.
!!
!!     This module uses python Mako template to generate Fortran code.
!!     Make sure to run the `mako_preprocess.py` if you edit the source
!!     .f90p file, and don`t edit the generated .f90 directly.
!!
!!=========================================================================

#include "f_defs.h"


module bgw_mpi_m
  use, intrinsic :: iso_c_binding
  use global_m

  implicit none

  interface bgw_bcast
     module procedure &
       bgw_bcast_str, &
       bgw_bcast_int_0, &
       bgw_bcast_int_1, &
       bgw_bcast_int_2, &
       bgw_bcast_int_3, &
       bgw_bcast_int_4, &
       bgw_bcast_int_5, &
       bgw_bcast_int_6, &
       bgw_bcast_int_7, &
       bgw_bcast_real_0, &
       bgw_bcast_real_1, &
       bgw_bcast_real_2, &
       bgw_bcast_real_3, &
       bgw_bcast_real_4, &
       bgw_bcast_real_5, &
       bgw_bcast_real_6, &
       bgw_bcast_real_7, &
       bgw_bcast_complex_0, &
       bgw_bcast_complex_1, &
       bgw_bcast_complex_2, &
       bgw_bcast_complex_3, &
       bgw_bcast_complex_4, &
       bgw_bcast_complex_5, &
       bgw_bcast_complex_6, &
       bgw_bcast_complex_7, &
       bgw_bcast_logical_0, &
       bgw_bcast_logical_1, &
       bgw_bcast_logical_2, &
       bgw_bcast_logical_3, &
       bgw_bcast_logical_4, &
       bgw_bcast_logical_5, &
       bgw_bcast_logical_6, &
       bgw_bcast_logical_7
  end interface
  interface bgw_gather
     module procedure &
       bgw_gather_int_0, &
       bgw_gather_int_1, &
       bgw_gather_int_2, &
       bgw_gather_int_3, &
       bgw_gather_int_4, &
       bgw_gather_int_5, &
       bgw_gather_int_6, &
       bgw_gather_int_7, &
       bgw_gather_real_0, &
       bgw_gather_real_1, &
       bgw_gather_real_2, &
       bgw_gather_real_3, &
       bgw_gather_real_4, &
       bgw_gather_real_5, &
       bgw_gather_real_6, &
       bgw_gather_real_7, &
       bgw_gather_complex_0, &
       bgw_gather_complex_1, &
       bgw_gather_complex_2, &
       bgw_gather_complex_3, &
       bgw_gather_complex_4, &
       bgw_gather_complex_5, &
       bgw_gather_complex_6, &
       bgw_gather_complex_7, &
       bgw_gather_logical_0, &
       bgw_gather_logical_1, &
       bgw_gather_logical_2, &
       bgw_gather_logical_3, &
       bgw_gather_logical_4, &
       bgw_gather_logical_5, &
       bgw_gather_logical_6, &
       bgw_gather_logical_7
  end interface
  interface bgw_allgather
     module procedure &
       bgw_allgather_int_0, &
       bgw_allgather_int_1, &
       bgw_allgather_int_2, &
       bgw_allgather_int_3, &
       bgw_allgather_int_4, &
       bgw_allgather_int_5, &
       bgw_allgather_int_6, &
       bgw_allgather_int_7, &
       bgw_allgather_real_0, &
       bgw_allgather_real_1, &
       bgw_allgather_real_2, &
       bgw_allgather_real_3, &
       bgw_allgather_real_4, &
       bgw_allgather_real_5, &
       bgw_allgather_real_6, &
       bgw_allgather_real_7, &
       bgw_allgather_complex_0, &
       bgw_allgather_complex_1, &
       bgw_allgather_complex_2, &
       bgw_allgather_complex_3, &
       bgw_allgather_complex_4, &
       bgw_allgather_complex_5, &
       bgw_allgather_complex_6, &
       bgw_allgather_complex_7, &
       bgw_allgather_logical_0, &
       bgw_allgather_logical_1, &
       bgw_allgather_logical_2, &
       bgw_allgather_logical_3, &
       bgw_allgather_logical_4, &
       bgw_allgather_logical_5, &
       bgw_allgather_logical_6, &
       bgw_allgather_logical_7
  end interface
  interface bgw_allgatherv
     module procedure &
       bgw_allgatherv_int_0, &
       bgw_allgatherv_int_1, &
       bgw_allgatherv_int_2, &
       bgw_allgatherv_int_3, &
       bgw_allgatherv_int_4, &
       bgw_allgatherv_int_5, &
       bgw_allgatherv_int_6, &
       bgw_allgatherv_int_7, &
       bgw_allgatherv_real_0, &
       bgw_allgatherv_real_1, &
       bgw_allgatherv_real_2, &
       bgw_allgatherv_real_3, &
       bgw_allgatherv_real_4, &
       bgw_allgatherv_real_5, &
       bgw_allgatherv_real_6, &
       bgw_allgatherv_real_7, &
       bgw_allgatherv_complex_0, &
       bgw_allgatherv_complex_1, &
       bgw_allgatherv_complex_2, &
       bgw_allgatherv_complex_3, &
       bgw_allgatherv_complex_4, &
       bgw_allgatherv_complex_5, &
       bgw_allgatherv_complex_6, &
       bgw_allgatherv_complex_7, &
       bgw_allgatherv_logical_0, &
       bgw_allgatherv_logical_1, &
       bgw_allgatherv_logical_2, &
       bgw_allgatherv_logical_3, &
       bgw_allgatherv_logical_4, &
       bgw_allgatherv_logical_5, &
       bgw_allgatherv_logical_6, &
       bgw_allgatherv_logical_7
  end interface

  public :: bgw_bcast, bgw_gather, bgw_allgather, bgw_allgatherv, &
    bgw_barrier, bgw_comm_size

contains


!==============================================================================
! Auxiliary routines
!==============================================================================


subroutine bgw_mpi_error(bgw_routine, mpi_routine, ierr, comm, root)
  character(len=*), intent(in) :: bgw_routine
  character(len=*), intent(in) :: mpi_routine
  integer, intent(in) :: ierr
  integer, intent(in), optional :: comm
  integer, intent(in), optional :: root

  character(len=256) :: error_str
  integer :: comm_, root_

  PUSH_SUB(bgw_mpi_error)

#ifdef MPI
  root_ = get_optional_arg(0, root)
  comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
  if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
    comm_/=MPI_COMM_WORLD) then
    write(0,'(/3a)') 'Error in ', trim(bgw_routine), ':'
    write(0,'(3a,i0/)') 'In call to ', trim(mpi_routine), ', got error=', ierr
  endif
  write(error_str,'(a,i0,4a)') 'Got error=', ierr, ' in ', &
    trim(bgw_routine), ' when calling ', trim(mpi_routine)
  call die(trim(error_str), only_root_writes=comm_==MPI_COMM_WORLD)
#endif

  POP_SUB(bgw_mpi_error)

end subroutine bgw_mpi_error


!==============================================================================
! High-level routines
!==============================================================================

!------------------------------------------------------------------------------
! Comm_size
!------------------------------------------------------------------------------

integer function bgw_comm_size(comm)
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_comm_size)

  bgw_comm_size = 1
#ifdef MPI
  comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
  call MPI_Comm_size(comm_, bgw_comm_size, ierr)
  if (ierr/=0) call bgw_mpi_error('bgw_comm_size', 'MPI_Comm_size', ierr, comm)
#endif

  POP_SUB(bgw_comm_size)

end function bgw_comm_size

!------------------------------------------------------------------------------
! Comm_rank
!------------------------------------------------------------------------------

integer function bgw_comm_rank(comm)
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_comm_rank)

  bgw_comm_rank = 0
#ifdef MPI
  comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
  call MPI_Comm_rank(comm_, bgw_comm_rank, ierr)
  if (ierr/=0) call bgw_mpi_error('bgw_comm_rank', 'MPI_Comm_rank', ierr, comm)
#endif

  POP_SUB(bgw_comm_rank)

end function bgw_comm_rank

!------------------------------------------------------------------------------
! Barrier
!------------------------------------------------------------------------------

subroutine bgw_barrier(comm)
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_barrier)

#ifdef MPI
  if (peinf%npes>1) then
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Barrier(comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_barrier', 'MPI_Barrier', ierr, comm)
  endif
#endif

  POP_SUB(bgw_barrier)

end subroutine bgw_barrier

!------------------------------------------------------------------------------
! Bcast
!------------------------------------------------------------------------------

subroutine bgw_bcast_int_0(x, root, comm)
  integer, intent(inout), target :: x
  integer, intent(in), optional :: root, comm

  integer, pointer :: x_1d(:)

  PUSH_SUB(bgw_bcast_int_0)
  call c_f_pointer(c_loc(x), x_1d, [1])
  call bgw_bcast_int_lowlevel(x_1d, 1, root=root, comm=comm)
  POP_SUB(bgw_bcast_int_0)

end subroutine bgw_bcast_int_0

subroutine bgw_bcast_int_1(x, root, comm)
  integer, intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_int_1)
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_int_1)

end subroutine bgw_bcast_int_1

subroutine bgw_bcast_int_2(x, root, comm)
  integer, intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_int_2)
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_int_2)

end subroutine bgw_bcast_int_2

subroutine bgw_bcast_int_3(x, root, comm)
  integer, intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_int_3)
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_int_3)

end subroutine bgw_bcast_int_3

subroutine bgw_bcast_int_4(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_int_4)
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_int_4)

end subroutine bgw_bcast_int_4

subroutine bgw_bcast_int_5(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_int_5)
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_int_5)

end subroutine bgw_bcast_int_5

subroutine bgw_bcast_int_6(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_int_6)
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_int_6)

end subroutine bgw_bcast_int_6

subroutine bgw_bcast_int_7(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_int_7)
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_int_7)

end subroutine bgw_bcast_int_7

subroutine bgw_bcast_real_0(x, root, comm)
  real(DP), intent(inout), target :: x
  integer, intent(in), optional :: root, comm

  real(DP), pointer :: x_1d(:)

  PUSH_SUB(bgw_bcast_real_0)
  call c_f_pointer(c_loc(x), x_1d, [1])
  call bgw_bcast_real_lowlevel(x_1d, 1, root=root, comm=comm)
  POP_SUB(bgw_bcast_real_0)

end subroutine bgw_bcast_real_0

subroutine bgw_bcast_real_1(x, root, comm)
  real(DP), intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_real_1)
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_real_1)

end subroutine bgw_bcast_real_1

subroutine bgw_bcast_real_2(x, root, comm)
  real(DP), intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_real_2)
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_real_2)

end subroutine bgw_bcast_real_2

subroutine bgw_bcast_real_3(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_real_3)
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_real_3)

end subroutine bgw_bcast_real_3

subroutine bgw_bcast_real_4(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_real_4)
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_real_4)

end subroutine bgw_bcast_real_4

subroutine bgw_bcast_real_5(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_real_5)
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_real_5)

end subroutine bgw_bcast_real_5

subroutine bgw_bcast_real_6(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_real_6)
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_real_6)

end subroutine bgw_bcast_real_6

subroutine bgw_bcast_real_7(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_real_7)
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_real_7)

end subroutine bgw_bcast_real_7

subroutine bgw_bcast_complex_0(x, root, comm)
  complex(DPC), intent(inout), target :: x
  integer, intent(in), optional :: root, comm

  complex(DPC), pointer :: x_1d(:)

  PUSH_SUB(bgw_bcast_complex_0)
  call c_f_pointer(c_loc(x), x_1d, [1])
  call bgw_bcast_complex_lowlevel(x_1d, 1, root=root, comm=comm)
  POP_SUB(bgw_bcast_complex_0)

end subroutine bgw_bcast_complex_0

subroutine bgw_bcast_complex_1(x, root, comm)
  complex(DPC), intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_complex_1)
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_complex_1)

end subroutine bgw_bcast_complex_1

subroutine bgw_bcast_complex_2(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_complex_2)
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_complex_2)

end subroutine bgw_bcast_complex_2

subroutine bgw_bcast_complex_3(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_complex_3)
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_complex_3)

end subroutine bgw_bcast_complex_3

subroutine bgw_bcast_complex_4(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_complex_4)
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_complex_4)

end subroutine bgw_bcast_complex_4

subroutine bgw_bcast_complex_5(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_complex_5)
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_complex_5)

end subroutine bgw_bcast_complex_5

subroutine bgw_bcast_complex_6(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_complex_6)
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_complex_6)

end subroutine bgw_bcast_complex_6

subroutine bgw_bcast_complex_7(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_complex_7)
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_complex_7)

end subroutine bgw_bcast_complex_7

subroutine bgw_bcast_logical_0(x, root, comm)
  logical, intent(inout), target :: x
  integer, intent(in), optional :: root, comm

  logical, pointer :: x_1d(:)

  PUSH_SUB(bgw_bcast_logical_0)
  call c_f_pointer(c_loc(x), x_1d, [1])
  call bgw_bcast_logical_lowlevel(x_1d, 1, root=root, comm=comm)
  POP_SUB(bgw_bcast_logical_0)

end subroutine bgw_bcast_logical_0

subroutine bgw_bcast_logical_1(x, root, comm)
  logical, intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_logical_1)
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_logical_1)

end subroutine bgw_bcast_logical_1

subroutine bgw_bcast_logical_2(x, root, comm)
  logical, intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_logical_2)
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_logical_2)

end subroutine bgw_bcast_logical_2

subroutine bgw_bcast_logical_3(x, root, comm)
  logical, intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_logical_3)
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_logical_3)

end subroutine bgw_bcast_logical_3

subroutine bgw_bcast_logical_4(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_logical_4)
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_logical_4)

end subroutine bgw_bcast_logical_4

subroutine bgw_bcast_logical_5(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_logical_5)
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_logical_5)

end subroutine bgw_bcast_logical_5

subroutine bgw_bcast_logical_6(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_logical_6)
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_logical_6)

end subroutine bgw_bcast_logical_6

subroutine bgw_bcast_logical_7(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  PUSH_SUB(bgw_bcast_logical_7)
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
  POP_SUB(bgw_bcast_logical_7)

end subroutine bgw_bcast_logical_7


!------------------------------------------------------------------------------
! Gather
!------------------------------------------------------------------------------

subroutine bgw_gather_int_0(x_in, x_out, root, comm)
  integer, intent(in), target :: x_in
  integer, intent(out), target :: x_out
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  integer, pointer :: x_in_1d(:)
  integer, pointer :: x_out_1d(:)

  PUSH_SUB(bgw_gather_int_0)
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_gather_int_lowlevel(x_in_1d, x_out_1d, 1, root, comm=comm)
  POP_SUB(bgw_gather_int_0)

end subroutine bgw_gather_int_0

subroutine bgw_gather_int_1(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:)
  integer, intent(out) :: x_out(:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_int_1)
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_int_1)

end subroutine bgw_gather_int_1

subroutine bgw_gather_int_2(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:)
  integer, intent(out) :: x_out(:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_int_2)
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_int_2)

end subroutine bgw_gather_int_2

subroutine bgw_gather_int_3(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:,:)
  integer, intent(out) :: x_out(:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_int_3)
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_int_3)

end subroutine bgw_gather_int_3

subroutine bgw_gather_int_4(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_int_4)
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_int_4)

end subroutine bgw_gather_int_4

subroutine bgw_gather_int_5(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_int_5)
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_int_5)

end subroutine bgw_gather_int_5

subroutine bgw_gather_int_6(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_int_6)
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_int_6)

end subroutine bgw_gather_int_6

subroutine bgw_gather_int_7(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_int_7)
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_int_7)

end subroutine bgw_gather_int_7

subroutine bgw_gather_real_0(x_in, x_out, root, comm)
  real(DP), intent(in), target :: x_in
  real(DP), intent(out), target :: x_out
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  real(DP), pointer :: x_in_1d(:)
  real(DP), pointer :: x_out_1d(:)

  PUSH_SUB(bgw_gather_real_0)
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_gather_real_lowlevel(x_in_1d, x_out_1d, 1, root, comm=comm)
  POP_SUB(bgw_gather_real_0)

end subroutine bgw_gather_real_0

subroutine bgw_gather_real_1(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:)
  real(DP), intent(out) :: x_out(:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_real_1)
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_real_1)

end subroutine bgw_gather_real_1

subroutine bgw_gather_real_2(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:)
  real(DP), intent(out) :: x_out(:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_real_2)
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_real_2)

end subroutine bgw_gather_real_2

subroutine bgw_gather_real_3(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:,:)
  real(DP), intent(out) :: x_out(:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_real_3)
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_real_3)

end subroutine bgw_gather_real_3

subroutine bgw_gather_real_4(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_real_4)
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_real_4)

end subroutine bgw_gather_real_4

subroutine bgw_gather_real_5(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_real_5)
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_real_5)

end subroutine bgw_gather_real_5

subroutine bgw_gather_real_6(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_real_6)
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_real_6)

end subroutine bgw_gather_real_6

subroutine bgw_gather_real_7(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_real_7)
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_real_7)

end subroutine bgw_gather_real_7

subroutine bgw_gather_complex_0(x_in, x_out, root, comm)
  complex(DPC), intent(in), target :: x_in
  complex(DPC), intent(out), target :: x_out
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  complex(DPC), pointer :: x_in_1d(:)
  complex(DPC), pointer :: x_out_1d(:)

  PUSH_SUB(bgw_gather_complex_0)
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_gather_complex_lowlevel(x_in_1d, x_out_1d, 1, root, comm=comm)
  POP_SUB(bgw_gather_complex_0)

end subroutine bgw_gather_complex_0

subroutine bgw_gather_complex_1(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:)
  complex(DPC), intent(out) :: x_out(:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_complex_1)
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_complex_1)

end subroutine bgw_gather_complex_1

subroutine bgw_gather_complex_2(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:)
  complex(DPC), intent(out) :: x_out(:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_complex_2)
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_complex_2)

end subroutine bgw_gather_complex_2

subroutine bgw_gather_complex_3(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_complex_3)
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_complex_3)

end subroutine bgw_gather_complex_3

subroutine bgw_gather_complex_4(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_complex_4)
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_complex_4)

end subroutine bgw_gather_complex_4

subroutine bgw_gather_complex_5(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_complex_5)
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_complex_5)

end subroutine bgw_gather_complex_5

subroutine bgw_gather_complex_6(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_complex_6)
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_complex_6)

end subroutine bgw_gather_complex_6

subroutine bgw_gather_complex_7(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_complex_7)
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_complex_7)

end subroutine bgw_gather_complex_7

subroutine bgw_gather_logical_0(x_in, x_out, root, comm)
  logical, intent(in), target :: x_in
  logical, intent(out), target :: x_out
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  logical, pointer :: x_in_1d(:)
  logical, pointer :: x_out_1d(:)

  PUSH_SUB(bgw_gather_logical_0)
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_gather_logical_lowlevel(x_in_1d, x_out_1d, 1, root, comm=comm)
  POP_SUB(bgw_gather_logical_0)

end subroutine bgw_gather_logical_0

subroutine bgw_gather_logical_1(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:)
  logical, intent(out) :: x_out(:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_logical_1)
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_logical_1)

end subroutine bgw_gather_logical_1

subroutine bgw_gather_logical_2(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:)
  logical, intent(out) :: x_out(:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_logical_2)
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_logical_2)

end subroutine bgw_gather_logical_2

subroutine bgw_gather_logical_3(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:,:)
  logical, intent(out) :: x_out(:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_logical_3)
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_logical_3)

end subroutine bgw_gather_logical_3

subroutine bgw_gather_logical_4(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_logical_4)
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_logical_4)

end subroutine bgw_gather_logical_4

subroutine bgw_gather_logical_5(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_logical_5)
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_logical_5)

end subroutine bgw_gather_logical_5

subroutine bgw_gather_logical_6(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_logical_6)
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_logical_6)

end subroutine bgw_gather_logical_6

subroutine bgw_gather_logical_7(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_gather_logical_7)
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
  POP_SUB(bgw_gather_logical_7)

end subroutine bgw_gather_logical_7


!------------------------------------------------------------------------------
! Allgather
!------------------------------------------------------------------------------

subroutine bgw_allgather_int_0(x_in, x_out, comm)
  integer, intent(in), target :: x_in
  integer, intent(out), target :: x_out
  integer, intent(in), optional :: comm

  integer, pointer :: x_in_1d(:)
  integer, pointer :: x_out_1d(:)

  PUSH_SUB(bgw_allgather_int_0)
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgather_int_lowlevel(x_in_1d, x_out_1d, 1, comm=comm)
  POP_SUB(bgw_allgather_int_0)

end subroutine bgw_allgather_int_0

subroutine bgw_allgather_int_1(x_in, x_out, comm)
  integer, intent(in) :: x_in(:)
  integer, intent(out) :: x_out(:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_int_1)
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_int_1)

end subroutine bgw_allgather_int_1

subroutine bgw_allgather_int_2(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:)
  integer, intent(out) :: x_out(:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_int_2)
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_int_2)

end subroutine bgw_allgather_int_2

subroutine bgw_allgather_int_3(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:,:)
  integer, intent(out) :: x_out(:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_int_3)
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_int_3)

end subroutine bgw_allgather_int_3

subroutine bgw_allgather_int_4(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_int_4)
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_int_4)

end subroutine bgw_allgather_int_4

subroutine bgw_allgather_int_5(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_int_5)
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_int_5)

end subroutine bgw_allgather_int_5

subroutine bgw_allgather_int_6(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_int_6)
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_int_6)

end subroutine bgw_allgather_int_6

subroutine bgw_allgather_int_7(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_int_7)
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_int_7)

end subroutine bgw_allgather_int_7

subroutine bgw_allgather_real_0(x_in, x_out, comm)
  real(DP), intent(in), target :: x_in
  real(DP), intent(out), target :: x_out
  integer, intent(in), optional :: comm

  real(DP), pointer :: x_in_1d(:)
  real(DP), pointer :: x_out_1d(:)

  PUSH_SUB(bgw_allgather_real_0)
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgather_real_lowlevel(x_in_1d, x_out_1d, 1, comm=comm)
  POP_SUB(bgw_allgather_real_0)

end subroutine bgw_allgather_real_0

subroutine bgw_allgather_real_1(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:)
  real(DP), intent(out) :: x_out(:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_real_1)
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_real_1)

end subroutine bgw_allgather_real_1

subroutine bgw_allgather_real_2(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:)
  real(DP), intent(out) :: x_out(:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_real_2)
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_real_2)

end subroutine bgw_allgather_real_2

subroutine bgw_allgather_real_3(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:,:)
  real(DP), intent(out) :: x_out(:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_real_3)
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_real_3)

end subroutine bgw_allgather_real_3

subroutine bgw_allgather_real_4(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_real_4)
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_real_4)

end subroutine bgw_allgather_real_4

subroutine bgw_allgather_real_5(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_real_5)
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_real_5)

end subroutine bgw_allgather_real_5

subroutine bgw_allgather_real_6(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_real_6)
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_real_6)

end subroutine bgw_allgather_real_6

subroutine bgw_allgather_real_7(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_real_7)
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_real_7)

end subroutine bgw_allgather_real_7

subroutine bgw_allgather_complex_0(x_in, x_out, comm)
  complex(DPC), intent(in), target :: x_in
  complex(DPC), intent(out), target :: x_out
  integer, intent(in), optional :: comm

  complex(DPC), pointer :: x_in_1d(:)
  complex(DPC), pointer :: x_out_1d(:)

  PUSH_SUB(bgw_allgather_complex_0)
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgather_complex_lowlevel(x_in_1d, x_out_1d, 1, comm=comm)
  POP_SUB(bgw_allgather_complex_0)

end subroutine bgw_allgather_complex_0

subroutine bgw_allgather_complex_1(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:)
  complex(DPC), intent(out) :: x_out(:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_complex_1)
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_complex_1)

end subroutine bgw_allgather_complex_1

subroutine bgw_allgather_complex_2(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:)
  complex(DPC), intent(out) :: x_out(:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_complex_2)
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_complex_2)

end subroutine bgw_allgather_complex_2

subroutine bgw_allgather_complex_3(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_complex_3)
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_complex_3)

end subroutine bgw_allgather_complex_3

subroutine bgw_allgather_complex_4(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_complex_4)
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_complex_4)

end subroutine bgw_allgather_complex_4

subroutine bgw_allgather_complex_5(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_complex_5)
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_complex_5)

end subroutine bgw_allgather_complex_5

subroutine bgw_allgather_complex_6(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_complex_6)
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_complex_6)

end subroutine bgw_allgather_complex_6

subroutine bgw_allgather_complex_7(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_complex_7)
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_complex_7)

end subroutine bgw_allgather_complex_7

subroutine bgw_allgather_logical_0(x_in, x_out, comm)
  logical, intent(in), target :: x_in
  logical, intent(out), target :: x_out
  integer, intent(in), optional :: comm

  logical, pointer :: x_in_1d(:)
  logical, pointer :: x_out_1d(:)

  PUSH_SUB(bgw_allgather_logical_0)
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgather_logical_lowlevel(x_in_1d, x_out_1d, 1, comm=comm)
  POP_SUB(bgw_allgather_logical_0)

end subroutine bgw_allgather_logical_0

subroutine bgw_allgather_logical_1(x_in, x_out, comm)
  logical, intent(in) :: x_in(:)
  logical, intent(out) :: x_out(:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_logical_1)
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_logical_1)

end subroutine bgw_allgather_logical_1

subroutine bgw_allgather_logical_2(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:)
  logical, intent(out) :: x_out(:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_logical_2)
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_logical_2)

end subroutine bgw_allgather_logical_2

subroutine bgw_allgather_logical_3(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:,:)
  logical, intent(out) :: x_out(:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_logical_3)
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_logical_3)

end subroutine bgw_allgather_logical_3

subroutine bgw_allgather_logical_4(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_logical_4)
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_logical_4)

end subroutine bgw_allgather_logical_4

subroutine bgw_allgather_logical_5(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_logical_5)
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_logical_5)

end subroutine bgw_allgather_logical_5

subroutine bgw_allgather_logical_6(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_logical_6)
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_logical_6)

end subroutine bgw_allgather_logical_6

subroutine bgw_allgather_logical_7(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in), optional :: comm

  PUSH_SUB(bgw_allgather_logical_7)
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
  POP_SUB(bgw_allgather_logical_7)

end subroutine bgw_allgather_logical_7


!------------------------------------------------------------------------------
! Allgatherv
!------------------------------------------------------------------------------

subroutine bgw_allgatherv_int_0(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in), target :: x_in
  integer, intent(out), target :: x_out
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: x_in_1d(:)
  integer, pointer :: x_out_1d(:)
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_int_0)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    siz = 1
  else
    siz = 1
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgatherv_int_lowlevel(x_in_1d, x_out_1d, 1, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_int_0)

end subroutine bgw_allgatherv_int_0

subroutine bgw_allgatherv_int_1(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:)
  integer, intent(out) :: x_out(:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_int_1)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_int_1)

end subroutine bgw_allgatherv_int_1

subroutine bgw_allgatherv_int_2(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:)
  integer, intent(out) :: x_out(:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_int_2)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_int_2)

end subroutine bgw_allgatherv_int_2

subroutine bgw_allgatherv_int_3(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:,:)
  integer, intent(out) :: x_out(:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_int_3)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_int_3)

end subroutine bgw_allgatherv_int_3

subroutine bgw_allgatherv_int_4(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_int_4)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_int_4)

end subroutine bgw_allgatherv_int_4

subroutine bgw_allgatherv_int_5(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_int_5)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_int_5)

end subroutine bgw_allgatherv_int_5

subroutine bgw_allgatherv_int_6(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_int_6)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_int_6)

end subroutine bgw_allgatherv_int_6

subroutine bgw_allgatherv_int_7(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_int_7)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_int_7)

end subroutine bgw_allgatherv_int_7

subroutine bgw_allgatherv_real_0(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in), target :: x_in
  real(DP), intent(out), target :: x_out
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  real(DP), pointer :: x_in_1d(:)
  real(DP), pointer :: x_out_1d(:)
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_real_0)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    siz = 1
  else
    siz = 1
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgatherv_real_lowlevel(x_in_1d, x_out_1d, 1, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_real_0)

end subroutine bgw_allgatherv_real_0

subroutine bgw_allgatherv_real_1(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:)
  real(DP), intent(out) :: x_out(:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_real_1)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_real_1)

end subroutine bgw_allgatherv_real_1

subroutine bgw_allgatherv_real_2(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:)
  real(DP), intent(out) :: x_out(:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_real_2)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_real_2)

end subroutine bgw_allgatherv_real_2

subroutine bgw_allgatherv_real_3(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:,:)
  real(DP), intent(out) :: x_out(:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_real_3)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_real_3)

end subroutine bgw_allgatherv_real_3

subroutine bgw_allgatherv_real_4(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_real_4)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_real_4)

end subroutine bgw_allgatherv_real_4

subroutine bgw_allgatherv_real_5(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_real_5)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_real_5)

end subroutine bgw_allgatherv_real_5

subroutine bgw_allgatherv_real_6(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_real_6)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_real_6)

end subroutine bgw_allgatherv_real_6

subroutine bgw_allgatherv_real_7(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_real_7)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_real_7)

end subroutine bgw_allgatherv_real_7

subroutine bgw_allgatherv_complex_0(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in), target :: x_in
  complex(DPC), intent(out), target :: x_out
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  complex(DPC), pointer :: x_in_1d(:)
  complex(DPC), pointer :: x_out_1d(:)
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_complex_0)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    siz = 1
  else
    siz = 1
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgatherv_complex_lowlevel(x_in_1d, x_out_1d, 1, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_complex_0)

end subroutine bgw_allgatherv_complex_0

subroutine bgw_allgatherv_complex_1(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:)
  complex(DPC), intent(out) :: x_out(:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_complex_1)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_complex_1)

end subroutine bgw_allgatherv_complex_1

subroutine bgw_allgatherv_complex_2(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:)
  complex(DPC), intent(out) :: x_out(:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_complex_2)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_complex_2)

end subroutine bgw_allgatherv_complex_2

subroutine bgw_allgatherv_complex_3(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_complex_3)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_complex_3)

end subroutine bgw_allgatherv_complex_3

subroutine bgw_allgatherv_complex_4(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_complex_4)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_complex_4)

end subroutine bgw_allgatherv_complex_4

subroutine bgw_allgatherv_complex_5(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_complex_5)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_complex_5)

end subroutine bgw_allgatherv_complex_5

subroutine bgw_allgatherv_complex_6(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_complex_6)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_complex_6)

end subroutine bgw_allgatherv_complex_6

subroutine bgw_allgatherv_complex_7(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_complex_7)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_complex_7)

end subroutine bgw_allgatherv_complex_7

subroutine bgw_allgatherv_logical_0(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in), target :: x_in
  logical, intent(out), target :: x_out
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  logical, pointer :: x_in_1d(:)
  logical, pointer :: x_out_1d(:)
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_logical_0)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    siz = 1
  else
    siz = 1
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgatherv_logical_lowlevel(x_in_1d, x_out_1d, 1, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_logical_0)

end subroutine bgw_allgatherv_logical_0

subroutine bgw_allgatherv_logical_1(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:)
  logical, intent(out) :: x_out(:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_logical_1)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_logical_1)

end subroutine bgw_allgatherv_logical_1

subroutine bgw_allgatherv_logical_2(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:)
  logical, intent(out) :: x_out(:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_logical_2)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_logical_2)

end subroutine bgw_allgatherv_logical_2

subroutine bgw_allgatherv_logical_3(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:,:)
  logical, intent(out) :: x_out(:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_logical_3)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_logical_3)

end subroutine bgw_allgatherv_logical_3

subroutine bgw_allgatherv_logical_4(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_logical_4)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_logical_4)

end subroutine bgw_allgatherv_logical_4

subroutine bgw_allgatherv_logical_5(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_logical_5)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_logical_5)

end subroutine bgw_allgatherv_logical_5

subroutine bgw_allgatherv_logical_6(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_logical_6)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_logical_6)

end subroutine bgw_allgatherv_logical_6

subroutine bgw_allgatherv_logical_7(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm

  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe

  PUSH_SUB(bgw_allgatherv_logical_7)

  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif

  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    SAFE_ALLOCATE(recvcounts_, (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif

  if (present(displs)) then
    displs_ => displs
  else
    SAFE_ALLOCATE(displs_, (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif

  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)

  if (.not.present(recvcounts)) then
    SAFE_DEALLOCATE_P(recvcounts_)
  endif
  if (.not.present(displs)) then
    SAFE_DEALLOCATE_P(displs_)
  endif

  POP_SUB(bgw_allgatherv_logical_7)

end subroutine bgw_allgatherv_logical_7



!==============================================================================
! Low-level implementations
!==============================================================================

!------------------------------------------------------------------------------
! Bcast: str
!------------------------------------------------------------------------------


subroutine bgw_bcast_str(x, root, comm)
  character(len=*), intent(inout) :: x
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, ierr

  PUSH_SUB(bgw_bcast_str)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = get_optional_arg(0, root)
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Bcast(x, len(x), MPI_CHARACTER, root_, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_bcast_str', 'MPI_Bcast', ierr, comm, root)
  endif
#endif

  POP_SUB(bgw_bcast_str)

end subroutine bgw_bcast_str


!------------------------------------------------------------------------------
! Bcast: int, real, complex, logical
!------------------------------------------------------------------------------


subroutine bgw_bcast_int_lowlevel(x, siz, root, comm)
  integer, intent(inout) :: x(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, ierr

  PUSH_SUB(bgw_bcast_int_lowlevel)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = get_optional_arg(0, root)
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Bcast(x, siz, MPI_INTEGER, root_, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_bcast_int_lowlevel', 'MPI_Bcast', ierr, comm, root)
  endif
#endif

  POP_SUB(bgw_bcast_int_lowlevel)

end subroutine bgw_bcast_int_lowlevel

subroutine bgw_bcast_real_lowlevel(x, siz, root, comm)
  real(DP), intent(inout) :: x(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, ierr

  PUSH_SUB(bgw_bcast_real_lowlevel)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = get_optional_arg(0, root)
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Bcast(x, siz, MPI_REAL_DP, root_, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_bcast_real_lowlevel', 'MPI_Bcast', ierr, comm, root)
  endif
#endif

  POP_SUB(bgw_bcast_real_lowlevel)

end subroutine bgw_bcast_real_lowlevel

subroutine bgw_bcast_complex_lowlevel(x, siz, root, comm)
  complex(DPC), intent(inout) :: x(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, ierr

  PUSH_SUB(bgw_bcast_complex_lowlevel)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = get_optional_arg(0, root)
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Bcast(x, siz, MPI_COMPLEX_DPC, root_, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_bcast_complex_lowlevel', 'MPI_Bcast', ierr, comm, root)
  endif
#endif

  POP_SUB(bgw_bcast_complex_lowlevel)

end subroutine bgw_bcast_complex_lowlevel

subroutine bgw_bcast_logical_lowlevel(x, siz, root, comm)
  logical, intent(inout) :: x(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, ierr

  PUSH_SUB(bgw_bcast_logical_lowlevel)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = get_optional_arg(0, root)
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Bcast(x, siz, MPI_LOGICAL, root_, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_bcast_logical_lowlevel', 'MPI_Bcast', ierr, comm, root)
  endif
#endif

  POP_SUB(bgw_bcast_logical_lowlevel)

end subroutine bgw_bcast_logical_lowlevel


!------------------------------------------------------------------------------
! Gather: int, real, complex, logical
!------------------------------------------------------------------------------

subroutine bgw_gather_int_lowlevel(x_in, x_out, siz, root, comm)
  integer, intent(in) :: x_in(*)
  integer, intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_gather_int_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Gather(x_in, siz, MPI_INTEGER, &
      x_out, siz, MPI_INTEGER, comm_, root, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_gather_int_lowlevel', 'MPI_Allgather', ierr, comm, root)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_gather_int_lowlevel)

end subroutine bgw_gather_int_lowlevel

subroutine bgw_gather_real_lowlevel(x_in, x_out, siz, root, comm)
  real(DP), intent(in) :: x_in(*)
  real(DP), intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_gather_real_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Gather(x_in, siz, MPI_REAL_DP, &
      x_out, siz, MPI_REAL_DP, comm_, root, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_gather_real_lowlevel', 'MPI_Allgather', ierr, comm, root)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_gather_real_lowlevel)

end subroutine bgw_gather_real_lowlevel

subroutine bgw_gather_complex_lowlevel(x_in, x_out, siz, root, comm)
  complex(DPC), intent(in) :: x_in(*)
  complex(DPC), intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_gather_complex_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Gather(x_in, siz, MPI_COMPLEX_DPC, &
      x_out, siz, MPI_COMPLEX_DPC, comm_, root, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_gather_complex_lowlevel', 'MPI_Allgather', ierr, comm, root)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_gather_complex_lowlevel)

end subroutine bgw_gather_complex_lowlevel

subroutine bgw_gather_logical_lowlevel(x_in, x_out, siz, root, comm)
  logical, intent(in) :: x_in(*)
  logical, intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in) :: root
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_gather_logical_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Gather(x_in, siz, MPI_LOGICAL, &
      x_out, siz, MPI_LOGICAL, comm_, root, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_gather_logical_lowlevel', 'MPI_Allgather', ierr, comm, root)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_gather_logical_lowlevel)

end subroutine bgw_gather_logical_lowlevel


!------------------------------------------------------------------------------
! Allgather: int, real, complex, logical
!------------------------------------------------------------------------------

subroutine bgw_allgather_int_lowlevel(x_in, x_out, siz, comm)
  integer, intent(in) :: x_in(*)
  integer, intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_allgather_int_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Allgather(x_in, siz, MPI_INTEGER, &
      x_out, siz, MPI_INTEGER, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_allgather_int_lowlevel', 'MPI_Allgather', ierr, comm)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_allgather_int_lowlevel)

end subroutine bgw_allgather_int_lowlevel

subroutine bgw_allgather_real_lowlevel(x_in, x_out, siz, comm)
  real(DP), intent(in) :: x_in(*)
  real(DP), intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_allgather_real_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Allgather(x_in, siz, MPI_REAL_DP, &
      x_out, siz, MPI_REAL_DP, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_allgather_real_lowlevel', 'MPI_Allgather', ierr, comm)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_allgather_real_lowlevel)

end subroutine bgw_allgather_real_lowlevel

subroutine bgw_allgather_complex_lowlevel(x_in, x_out, siz, comm)
  complex(DPC), intent(in) :: x_in(*)
  complex(DPC), intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_allgather_complex_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Allgather(x_in, siz, MPI_COMPLEX_DPC, &
      x_out, siz, MPI_COMPLEX_DPC, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_allgather_complex_lowlevel', 'MPI_Allgather', ierr, comm)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_allgather_complex_lowlevel)

end subroutine bgw_allgather_complex_lowlevel

subroutine bgw_allgather_logical_lowlevel(x_in, x_out, siz, comm)
  logical, intent(in) :: x_in(*)
  logical, intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_allgather_logical_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Allgather(x_in, siz, MPI_LOGICAL, &
      x_out, siz, MPI_LOGICAL, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_allgather_logical_lowlevel', 'MPI_Allgather', ierr, comm)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_allgather_logical_lowlevel)

end subroutine bgw_allgather_logical_lowlevel


!------------------------------------------------------------------------------
! Allgatherv: int, real, complex, logical
!------------------------------------------------------------------------------

subroutine bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts, displs, comm)
  integer, intent(in) :: x_in(*)
  integer, intent(out) :: x_out(*)
  integer, intent(in) :: siz, recvcounts(*), displs(*)
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_allgatherv_int_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Allgatherv(x_in, siz, MPI_INTEGER, &
      x_out, recvcounts, displs, MPI_INTEGER, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_allgatherv_int_lowlevel', 'MPI_Allgatherv', ierr, comm)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_allgatherv_int_lowlevel)

end subroutine bgw_allgatherv_int_lowlevel

subroutine bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(*)
  real(DP), intent(out) :: x_out(*)
  integer, intent(in) :: siz, recvcounts(*), displs(*)
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_allgatherv_real_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Allgatherv(x_in, siz, MPI_REAL_DP, &
      x_out, recvcounts, displs, MPI_REAL_DP, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_allgatherv_real_lowlevel', 'MPI_Allgatherv', ierr, comm)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_allgatherv_real_lowlevel)

end subroutine bgw_allgatherv_real_lowlevel

subroutine bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(*)
  complex(DPC), intent(out) :: x_out(*)
  integer, intent(in) :: siz, recvcounts(*), displs(*)
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_allgatherv_complex_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Allgatherv(x_in, siz, MPI_COMPLEX_DPC, &
      x_out, recvcounts, displs, MPI_COMPLEX_DPC, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_allgatherv_complex_lowlevel', 'MPI_Allgatherv', ierr, comm)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_allgatherv_complex_lowlevel)

end subroutine bgw_allgatherv_complex_lowlevel

subroutine bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts, displs, comm)
  logical, intent(in) :: x_in(*)
  logical, intent(out) :: x_out(*)
  integer, intent(in) :: siz, recvcounts(*), displs(*)
  integer, intent(in), optional :: comm

  integer :: comm_, ierr

  PUSH_SUB(bgw_allgatherv_logical_lowlevel)

  if (peinf%npes>1) then
#ifdef MPI
    comm_ = get_optional_arg(MPI_COMM_WORLD, comm)
    call MPI_Allgatherv(x_in, siz, MPI_LOGICAL, &
      x_out, recvcounts, displs, MPI_LOGICAL, comm_, ierr)
    if (ierr/=0) call bgw_mpi_error('bgw_allgatherv_logical_lowlevel', 'MPI_Allgatherv', ierr, comm)
#endif
  else
    x_out(1:siz) = x_in(1:siz)
  endif

  POP_SUB(bgw_allgatherv_logical_lowlevel)

end subroutine bgw_allgatherv_logical_lowlevel


end module bgw_mpi_m
