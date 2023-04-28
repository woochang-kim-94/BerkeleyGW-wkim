! File bgw_mpi.f90 automatically created from bgw_mpi.f90p by mako_preprocess.py.
! Do not edit the resulting file (bgw_mpi.f90) directly!

#include "f_defs.h"


module bgw_mpi_m

  use global_m

  implicit none

  interface bgw_bcast
    module procedure &
      bgw_bcast_int0, &
      bgw_bcast_int1, &
      bgw_bcast_int2, &
      bgw_bcast_int3, &
      bgw_bcast_int4, &
      bgw_bcast_int5, &
      bgw_bcast_int6, &
      bgw_bcast_int7, &
      bgw_bcast_real0, &
      bgw_bcast_real1, &
      bgw_bcast_real2, &
      bgw_bcast_real3, &
      bgw_bcast_real4, &
      bgw_bcast_real5, &
      bgw_bcast_real6, &
      bgw_bcast_real7, &
      bgw_bcast_complex0, &
      bgw_bcast_complex1, &
      bgw_bcast_complex2, &
      bgw_bcast_complex3, &
      bgw_bcast_complex4, &
      bgw_bcast_complex5, &
      bgw_bcast_complex6, &
      bgw_bcast_complex7, &
      bgw_bcast_logical0, &
      bgw_bcast_logical1, &
      bgw_bcast_logical2, &
      bgw_bcast_logical3, &
      bgw_bcast_logical4, &
      bgw_bcast_logical5, &
      bgw_bcast_logical6, &
      bgw_bcast_logical7
  end interface

  public :: bgw_bcast

contains


subroutine bgw_bcast_int0(x, root, comm)
  integer, intent(inout) :: x
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_int0)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, 1, MPI_INTEGER, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_int0:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_int0', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_int0)

end subroutine bgw_bcast_int0


subroutine bgw_bcast_int1(x, root, comm)
  integer, intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_int1)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_INTEGER, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_int1:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_int1', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_int1)

end subroutine bgw_bcast_int1


subroutine bgw_bcast_int2(x, root, comm)
  integer, intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_int2)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_INTEGER, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_int2:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_int2', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_int2)

end subroutine bgw_bcast_int2


subroutine bgw_bcast_int3(x, root, comm)
  integer, intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_int3)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_INTEGER, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_int3:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_int3', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_int3)

end subroutine bgw_bcast_int3


subroutine bgw_bcast_int4(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_int4)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_INTEGER, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_int4:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_int4', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_int4)

end subroutine bgw_bcast_int4


subroutine bgw_bcast_int5(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_int5)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_INTEGER, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_int5:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_int5', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_int5)

end subroutine bgw_bcast_int5


subroutine bgw_bcast_int6(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_int6)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_INTEGER, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_int6:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_int6', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_int6)

end subroutine bgw_bcast_int6


subroutine bgw_bcast_int7(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_int7)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_INTEGER, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_int7:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_int7', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_int7)

end subroutine bgw_bcast_int7


subroutine bgw_bcast_real0(x, root, comm)
  real(DP), intent(inout) :: x
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_real0)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, 1, MPI_REAL_DP, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_real0:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_real0', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_real0)

end subroutine bgw_bcast_real0


subroutine bgw_bcast_real1(x, root, comm)
  real(DP), intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_real1)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_REAL_DP, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_real1:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_real1', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_real1)

end subroutine bgw_bcast_real1


subroutine bgw_bcast_real2(x, root, comm)
  real(DP), intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_real2)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_REAL_DP, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_real2:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_real2', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_real2)

end subroutine bgw_bcast_real2


subroutine bgw_bcast_real3(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_real3)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_REAL_DP, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_real3:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_real3', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_real3)

end subroutine bgw_bcast_real3


subroutine bgw_bcast_real4(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_real4)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_REAL_DP, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_real4:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_real4', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_real4)

end subroutine bgw_bcast_real4


subroutine bgw_bcast_real5(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_real5)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_REAL_DP, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_real5:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_real5', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_real5)

end subroutine bgw_bcast_real5


subroutine bgw_bcast_real6(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_real6)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_REAL_DP, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_real6:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_real6', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_real6)

end subroutine bgw_bcast_real6


subroutine bgw_bcast_real7(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_real7)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_REAL_DP, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_real7:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_real7', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_real7)

end subroutine bgw_bcast_real7


subroutine bgw_bcast_complex0(x, root, comm)
  complex(DPC), intent(inout) :: x
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_complex0)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, 1, MPI_COMPLEX_DPC, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_complex0:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_complex0', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_complex0)

end subroutine bgw_bcast_complex0


subroutine bgw_bcast_complex1(x, root, comm)
  complex(DPC), intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_complex1)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_COMPLEX_DPC, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_complex1:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_complex1', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_complex1)

end subroutine bgw_bcast_complex1


subroutine bgw_bcast_complex2(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_complex2)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_COMPLEX_DPC, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_complex2:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_complex2', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_complex2)

end subroutine bgw_bcast_complex2


subroutine bgw_bcast_complex3(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_complex3)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_COMPLEX_DPC, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_complex3:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_complex3', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_complex3)

end subroutine bgw_bcast_complex3


subroutine bgw_bcast_complex4(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_complex4)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_COMPLEX_DPC, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_complex4:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_complex4', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_complex4)

end subroutine bgw_bcast_complex4


subroutine bgw_bcast_complex5(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_complex5)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_COMPLEX_DPC, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_complex5:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_complex5', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_complex5)

end subroutine bgw_bcast_complex5


subroutine bgw_bcast_complex6(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_complex6)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_COMPLEX_DPC, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_complex6:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_complex6', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_complex6)

end subroutine bgw_bcast_complex6


subroutine bgw_bcast_complex7(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_complex7)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_COMPLEX_DPC, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_complex7:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_complex7', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_complex7)

end subroutine bgw_bcast_complex7


subroutine bgw_bcast_logical0(x, root, comm)
  logical, intent(inout) :: x
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_logical0)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, 1, MPI_LOGICAL, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_logical0:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_logical0', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_logical0)

end subroutine bgw_bcast_logical0


subroutine bgw_bcast_logical1(x, root, comm)
  logical, intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_logical1)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_LOGICAL, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_logical1:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_logical1', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_logical1)

end subroutine bgw_bcast_logical1


subroutine bgw_bcast_logical2(x, root, comm)
  logical, intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_logical2)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_LOGICAL, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_logical2:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_logical2', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_logical2)

end subroutine bgw_bcast_logical2


subroutine bgw_bcast_logical3(x, root, comm)
  logical, intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_logical3)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_LOGICAL, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_logical3:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_logical3', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_logical3)

end subroutine bgw_bcast_logical3


subroutine bgw_bcast_logical4(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_logical4)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_LOGICAL, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_logical4:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_logical4', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_logical4)

end subroutine bgw_bcast_logical4


subroutine bgw_bcast_logical5(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_logical5)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_LOGICAL, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_logical5:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_logical5', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_logical5)

end subroutine bgw_bcast_logical5


subroutine bgw_bcast_logical6(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_logical6)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_LOGICAL, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_logical6:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_logical6', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_logical6)

end subroutine bgw_bcast_logical6


subroutine bgw_bcast_logical7(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm

  integer :: root_, comm_, error

  PUSH_SUB(bgw_bcast_logical7)

#ifdef MPI
  if (peinf%npes>1) then
    root_ = 0
    comm_ = MPI_COMM_WORLD
    if (present(root)) root_ = root
    if (present(comm)) comm_ = comm
    call MPI_Bcast(x, size(x), MPI_LOGICAL, root_, comm_, error)
    if (error/=0) then
      if ((peinf%inode==root_.and.comm_==MPI_COMM_WORLD) .or. &
        comm_/=MPI_COMM_WORLD) then
        write(0,'(/a)') 'Error in bgw_bcast_logical7:'
        write(0,'(a,i0/)') 'In call to MPI_Bcast, got error=',error
      endif
      call die('Got error/=0 for MPI_Bcast in bgw_bcast_logical7', &
        only_root_writes=comm_==MPI_COMM_WORLD)
    endif
  endif
#endif

  POP_SUB(bgw_bcast_logical7)

end subroutine bgw_bcast_logical7


end module bgw_mpi_m
