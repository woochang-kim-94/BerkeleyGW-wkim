#include "f_defs.h"
#if defined MPI && !defined USESCALAPACK
  #error ScaLAPACK is required for MPI builds.
#endif

module distribution_m

  use global_m
  use blas_m
  use lapack_m
  use scalapack_m
  use kpoint_pool_m, only: kpoint_pool_t

  implicit none

  private

    type distrib_mat_t
      integer :: cntxt
      integer :: desc(9)
      !> Number of rows in distributed matrix
      integer :: M
      !> Number of cols in distributed matrix
      integer :: N
      !> Number of rows in local matrix buffer
      integer :: Ml
      !> Number of cols in local matrix buffer
      integer :: Nl
      !> Number of rows of the distributed matrix that I own (<=Ml)
      integer :: Mown
      !> Number of cols of the distributed matrix that I own (<=Nl)
      integer :: Nown
      !> Block size
      integer :: Mb, Nb
      !> Number of processors in each row of the *processor grid*
      integer :: nprow
      !> Number of processors in each col of the *processor grid*
      integer :: npcol
      !> My row number in the blacs *processor grid*
      integer :: myprow
      !> My col number in the blacs *processor grid*
      integer :: mypcol
    contains
      procedure, private :: setup_serial_mat => distrib_setup_serial_mat
      procedure :: setup_mat_1d => distrib_setup_mat_1d
      procedure :: setup_mat_2d => distrib_setup_mat_2d
    endtype distrib_mat_t

  public :: distrib_mat_t

contains

subroutine distrib_setup_serial_mat(dm, M, N)
  class(distrib_mat_t), intent(inout) :: dm
  integer, intent(in) :: M
  integer, intent(in) :: N

  PUSH_SUB(distrib_setup_serial_mat)

  dm%nprow = 1
  dm%npcol = 1
  dm%myprow = 0
  dm%mypcol = 0
  dm%M = M
  dm%N = N
  dm%Mb = M
  dm%Nb = N
  dm%Mown = M
  dm%Nown = N
  dm%Ml = M
  dm%Nl = N

  PUSH_SUB(distrib_setup_serial_mat)

end subroutine distrib_setup_serial_mat


subroutine distrib_setup_mat_1d(dm, M, N, kpp, row_or_col)
  class(distrib_mat_t), intent(inout) :: dm
  integer, intent(in) :: M
  integer, intent(in) :: N
  type(kpoint_pool_t), intent(in) :: kpp
  character, optional, intent(in) :: row_or_col

  integer :: info
  character :: row_or_col_

  PUSH_SUB(distrib_setup_mat_1d)

  row_or_col_ = 'c'
  if (present(row_or_col)) row_or_col_ = row_or_col

#ifdef MPI
  ! Global matrix size
  dm%N = N
  dm%M = M
  if (row_or_col_=='c'.or.row_or_col_=='C') then
    dm%cntxt = kpp%cntxt_c
    ! Blacs processor grid
    dm%nprow = 1
    dm%npcol = kpp%npes
    dm%myprow = 0
    dm%mypcol = kpp%inode
    ! Block sizes
    dm%Nb = DIVUP(N, kpp%npes)
    dm%Mb = dm%Nb
    ! Local matrix buffer size
    dm%Nl = dm%Nb
    dm%Ml = M
    ! Number of rows/cols I own
    dm%Nown = numroc(N, dm%Nb, kpp%inode, 0, kpp%npes)
    dm%Mown = M
  elseif (row_or_col_=='r'.or.row_or_col_=='R') then
    dm%cntxt = kpp%cntxt_r
    ! Blacs processor grid
    dm%npcol = 1
    dm%nprow = kpp%npes
    dm%mypcol = 0
    dm%myprow = kpp%inode
    ! Block sizes
    dm%Mb = DIVUP(M, kpp%npes)
    dm%Nb = dm%Mb
    ! Local matrix buffer size
    dm%Ml = dm%Mb
    dm%Nl = N
    ! Number of rows/cols I own
    dm%Mown = numroc(M, dm%Mb, kpp%inode, 0, kpp%npes)
    dm%Nown = N
  else
    call die('Unknown value for row_of_col: '//row_or_col_, only_root_writes=.true.)
  endif

  call descinit(dm%desc, dm%M, dm%N, dm%Mb, dm%Nb, 0, 0, dm%cntxt, dm%Ml, info)
  if (info/=0) call die('got info/=0 in descinit')
#else
  call dm%setup_serial_mat(M, N)
#endif

  POP_SUB(distrib_setup_mat_1d)

end subroutine distrib_setup_mat_1d


subroutine distrib_setup_mat_2d(dm, M, N, block_size, kpp)
  class(distrib_mat_t), intent(inout) :: dm
  integer, intent(in) :: M
  integer, intent(in) :: N
  integer, intent(in) :: block_size
  type(kpoint_pool_t), intent(in) :: kpp

  PUSH_SUB(distrib_setup_mat_2d)

#ifdef MPI
  dm%cntxt = kpp%cntxt_2d
  ! Block sizes
  dm%Nb = block_size
  dm%Mb = block_size
  ! Blacs processor grid
  dm%nprow = int(sqrt(dble(kpp%npes)))
  dm%npcol = kpp%npes/dm%nprow
  call blacs_gridinfo(dm%cntxt, dm%nprow, dm%npcol, dm%myprow, dm%mypcol)
  if (kpp%inode==0) then
    write(kpp%iunit,'(/,1x,4(a,i0))') 'BLACS processor grid: ', &
      dm%nprow, ' x ', dm%npcol, '; BLOCKSIZE = ', dm%Mb,' x ', dm%Nb
    write(kpp%iunit,'(1x,a,i0,/)') 'Number of idle processors: ', kpp%npes - dm%nprow*dm%npcol
    FLUSH(kpp%iunit)
  endif
  ! Global matrix size
  dm%N = N
  dm%M = M
  ! Number of rows/cols I own
  if (dm%cntxt>=0) then
    dm%Mown = numroc(M, dm%Mb, dm%myprow, 0, dm%nprow)
    dm%Nown = numroc(N, dm%Nb, dm%mypcol, 0, dm%npcol)
  else
    dm%Mown = 0
    dm%Nown = 0
  endif
  ! Local matrix buffer size
  dm%Ml = max(1, dm%Mown)
  dm%Nl = max(1, dm%Nown)
  call descset(dm%desc, dm%M, dm%N, dm%Mb, dm%Nb, 0, 0, dm%cntxt, dm%Ml)
#else
  call dm%setup_serial_mat(M, N)
#endif

  POP_SUB(distrib_setup_mat_2d)

end subroutine distrib_setup_mat_2d


end module distribution_m
