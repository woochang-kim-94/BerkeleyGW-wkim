!================================================================================
!
! Modules:
!
! (1) scalapack_hl_m      Originally by FHJ 06/2020
!
!     High-level routines for ScaLAPACK (e.g., diagonaliation, etc.)
!
!================================================================================

#include "f_defs.h"

module scalapack_hl_m

#ifdef USESCALAPACK
  use global_m
  use scalapack_m
  implicit none

  private

  interface realloc_query_workspace
    module procedure realloc_query_workspace_int, realloc_query_workspace_real, realloc_query_workspace_cplx
  end interface realloc_query_workspace

  public :: &
    diagonalize_scalapack

  interface diagonalize_scalapack
    module procedure &
        diagonalize_scalapack_real, diagonalize_scalapack_cplx
  end interface

  interface diagonalize_scalapack_x
    module procedure &
        diagonalize_scalapack_x_real, diagonalize_scalapack_x_cplx
  end interface

#ifdef USEMR3
  interface diagonalize_scalapack_r
    module procedure &
        diagonalize_scalapack_r_real, diagonalize_scalapack_r_cplx
  end interface
#endif

  interface diagonalize_scalapack_d
    module procedure &
        diagonalize_scalapack_d_real, diagonalize_scalapack_d_cplx
  end interface

contains


!-------------------------------------------------------------------------------
! Internal routines
!-------------------------------------------------------------------------------

subroutine check_lapack_ev_result(routine_str, info, nwant, nfound, nzfound, is_query)
  character(len=*), intent(in) :: routine_str
  integer, optional, intent(in) :: info
  integer, optional, intent(in) :: nwant
  integer, optional, intent(in) :: nfound
  integer, optional, intent(in) :: nzfound
  logical, optional, intent(in) :: is_query

  character(len=64) :: routname
  character(len=128) :: errmsg

  PUSH_SUB(check_lapack_ev_result)

  if (present(is_query)) then
    if (is_query) write(routname,'(a)') trim(routine_str)//' (query mode)'
  else
    routname = trim(routine_str)
  endif
  if (info/=0) then
    if ((trim(routname)=='PZHEEVX') .and. info==2) then
      ! Reorthogonalization problem, but not critical.
      if (peinf%inode==0) then
        write(0,'(/,a)') 'WARNING: got info=2 from call to PZHEEVX'
        write(0,'(a)') '(failed to reorthogonalize eigenvectors).'
        write(0,'(a)') 'You may want to rerun BerkeleyGW with a different ScaLAPACK solver'
        write(0,'(a/)') 'using the environment variable BGW_SCALAPACK_ALGO'
      endif
    else
      ! Fatal error
      if (peinf%inode==0) then
        write(0,'(/,a)') 'ERROR: Diagonalization with ', trim(routname),' failed:'
        if (present(nfound)) then
          if (nfound/=nwant) write(0,'(3(a,i0))') 'Only ', nfound, ' out of ', nwant, ' eigvals found.'
        endif
        if (present(nzfound)) then
          if (nzfound/=nwant) write(0,'(3(a,i0))') 'Only ', nfound, ' out of ', nwant, ' eigvecs found.'
        endif
        write(0,'(a,i0,/)') 'Got info=', info
      endif
      write(errmsg,'(3a)') trim(routname), ' failed with info=', info
      call die(errmsg, only_root_writes=.true.)
    endif
  endif

  POP_SUB(check_lapack_ev_result)

end subroutine check_lapack_ev_result


subroutine safe2int(real_in, int_out, rescaling)
  real(DP), intent(in) :: real_in
  integer, intent(out) :: int_out
  real(DP), intent(in), optional :: rescaling

  character(len=128) :: scalapack_mem_rescale_str
  real(DP) :: x, scalapack_mem_rescale
  integer :: err

  PUSH_SUB(safe2int)

  x = real_in
  if (present(rescaling)) x = x * rescaling
  call get_environment_variable("BGW_SCALAPACK_MEM_RESCALE", scalapack_mem_rescale_str)
  read(scalapack_mem_rescale_str,*,iostat=err) scalapack_mem_rescale
  if (err==0) x = x * scalapack_mem_rescale

  if (peinf%verb_high .and. peinf%inode==0) then
    write(6,'(/1x,a,f0.3/)') 'Workspace memory rescaling factor: ', rescaling
  endif
  x = x + 16d0

  int_out = huge(int_out)
  if (x > 0d0) then
    if (x < dble(int_out)) then
      int_out = int(x)
    endif
  endif

  POP_SUB(safe2int)

end subroutine safe2int


subroutine realloc_query_workspace_int(buf, length, rescaling)
  integer, allocatable, intent(inout) :: buf(:)
  integer, intent(out) :: length
  real(DP), intent(in), optional :: rescaling

  PUSH_SUB(realloc_query_workspace_int)

  if (.not.(allocated(buf))) then
    call die('INTERNAL ERROR: buffer `buf` not allocated in call to realloc_query_workspace_int!')
  endif

  call safe2int(dble(buf(1)), length, rescaling)
  SAFE_DEALLOCATE(buf)
  SAFE_ALLOCATE(buf, (length))

  POP_SUB(realloc_query_workspace_int)

end subroutine realloc_query_workspace_int


subroutine realloc_query_workspace_real(buf, length, rescaling)
  real(DP), allocatable, intent(inout) :: buf(:)
  integer, intent(out) :: length
  real(DP), intent(in), optional :: rescaling

  PUSH_SUB(realloc_query_workspace_real)

  if (.not.(allocated(buf))) then
    call die('INTERNAL ERROR: buffer `buf` not allocated in call to realloc_query_workspace_real!')
  endif

  call safe2int(dble(buf(1)), length, rescaling)
  SAFE_DEALLOCATE(buf)
  SAFE_ALLOCATE(buf, (length))

  POP_SUB(realloc_query_workspace_real)

end subroutine realloc_query_workspace_real


subroutine realloc_query_workspace_cplx(buf, length, rescaling)
  complex(DPC), allocatable, intent(inout) :: buf(:)
  integer, intent(out) :: length
  real(DP), intent(in), optional :: rescaling

  PUSH_SUB(realloc_query_workspace_cplx)

  if (.not.(allocated(buf))) then
    call die('INTERNAL ERROR: buffer `buf` not allocated in call to realloc_query_workspace_real!')
  endif

  call safe2int(dble(buf(1)), length, rescaling)
  SAFE_DEALLOCATE(buf)
  SAFE_ALLOCATE(buf, (length))

  POP_SUB(realloc_query_workspace_cplx)

end subroutine realloc_query_workspace_cplx



!-------------------------------------------------------------------------------
! bissection algorithm
!-------------------------------------------------------------------------------



!> Diagonalize a matrix with ScaLAPACK using the bissection algorithm.
!! Implementation for real type.
subroutine diagonalize_scalapack_x_real(desc, Neig, matrix, eigenval, eigen_vect, &
                                        iunit)
  integer :: desc(9) !< ScaLAPACK/BLACS descriptor. See documentation in DESCINIT
  integer, intent(in) :: Neig !< Number of eigenvalues/eigenvectors to compute
  real(DP), intent(inout) :: matrix(*) !< Local buffer for the 2D distributed matrix
  real(DP), intent(out) :: eigenval(*) !< Eigenvalues. Only the first Neig will be written.
  real(DP), intent(inout) :: eigen_vect(*) !< Local buffer for the 2D distributed eigenvectors
  integer, intent(in), optional :: iunit !< Unit to write output. Defaults to 6.

  ! BLACS stuff
  integer :: ctxt, N, Mb, Nb
  integer :: nprow, npcol, myprow, mypcol, Ml, Nl

  real(DP), parameter :: scalapack_mem_rescale = 1d0
  integer :: info, inode, npes, iunit_

  real(DP), allocatable :: work(:)
  integer :: lwork, liwork
  real(DP) :: en_t(desc(3))
  integer, allocatable :: iwork(:)
  integer :: nfound, nzfound

  character :: range_char
  real(DP), parameter :: orfac = 5d-7
  real(DP), allocatable :: gap(:)
  integer, allocatable :: ifail(:), iclustr(:)

  PUSH_SUB(diagonalize_scalapack_x_real)

  iunit_ = 6
  if (present(iunit)) iunit_ = iunit
  ctxt = desc(2)
  if (desc(3)/=desc(4)) call die('Matrix is not square')
  N = desc(4)
  Mb = desc(5)
  Nb = desc(6)
  if (Neig == N) then
    range_char = 'A'
  else
    range_char = 'I'
  endif

  if (ctxt>=0) then

    call blacs_gridinfo(ctxt, nprow, npcol, myprow, mypcol)
    Ml = NUMROC(N, Mb, myprow, 0, nprow) ! # of local rows I own
    Nl = NUMROC(N, Nb, mypcol, 0, npcol) ! # of local columns I own

    nfound = Neig
    nzfound = Neig
    call blacs_pinfo(inode, npes)
    if (inode==0) then
      write(iunit_,'(1x,2a)') 'Using diagonalization routine: PDSYEVX'
      FLUSH(iunit_)
    endif

    ! FHJ: Call in query mode
    SAFE_ALLOCATE(iclustr, (max(1,2*nprow*npcol)))
    SAFE_ALLOCATE(gap, (max(1,nprow*npcol)))
    SAFE_ALLOCATE(ifail, (N))

    SAFE_ALLOCATE(work, (10))
    SAFE_ALLOCATE(iwork, (10))

    lwork = -1
    liwork = -1

    call scalapack_call()

    call check_lapack_ev_result('PDSYEVX', info, Neig, is_query=.true.)
    call realloc_query_workspace(work, lwork, scalapack_mem_rescale)
    call realloc_query_workspace(iwork, liwork, scalapack_mem_rescale)

    ! FHJ: Call scalapack for real
    if (inode==0) then
      write(iunit_,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', N
      FLUSH(iunit_)
    endif

    call scalapack_call()

    if (inode==0) then
      write(iunit_,'(1x,a)') 'Done ScaLAPACK diagonalization'
      FLUSH(iunit_)
    endif
    call check_lapack_ev_result('PDSYEVX', info, Neig, nfound, nzfound)

    SAFE_DEALLOCATE(work)
    SAFE_DEALLOCATE(iwork)

    ! FHJ: Copy eigenvalues/eigenvectors
    eigenval(1:Neig) = en_t(1:Neig)
    SAFE_DEALLOCATE(iclustr)
    SAFE_DEALLOCATE(gap)
    SAFE_DEALLOCATE(ifail)
  endif

  POP_SUB(diagonalize_scalapack_x_real)

contains

  subroutine scalapack_call

    PUSH_SUB(diagonalize_scalapack_x_real.scalapack_call)

    call PDSYEVX &
      ('V', range_char, 'U', N, matrix, 1, 1, desc, &
      0d0, 0d0, 1, Neig, 0d0, nfound, &
      nzfound, en_t, orfac, eigen_vect, 1, 1, desc, work, lwork, &
      iwork, liwork, ifail, iclustr, gap, info)

    POP_SUB(diagonalize_scalapack_x_real.scalapack_call)

  end subroutine scalapack_call

end subroutine diagonalize_scalapack_x_real



!> Diagonalize a matrix with ScaLAPACK using the bissection algorithm.
!! Implementation for cplx type.
subroutine diagonalize_scalapack_x_cplx(desc, Neig, matrix, eigenval, eigen_vect, &
                                        iunit)
  integer :: desc(9) !< ScaLAPACK/BLACS descriptor. See documentation in DESCINIT
  integer, intent(in) :: Neig !< Number of eigenvalues/eigenvectors to compute
  complex(DPC), intent(inout) :: matrix(*) !< Local buffer for the 2D distributed matrix
  real(DP), intent(out) :: eigenval(*) !< Eigenvalues. Only the first Neig will be written.
  complex(DPC), intent(inout) :: eigen_vect(*) !< Local buffer for the 2D distributed eigenvectors
  integer, intent(in), optional :: iunit !< Unit to write output. Defaults to 6.

  ! BLACS stuff
  integer :: ctxt, N, Mb, Nb
  integer :: nprow, npcol, myprow, mypcol, Ml, Nl

  real(DP), parameter :: scalapack_mem_rescale = 1d0
  integer :: info, inode, npes, iunit_

  complex(DPC), allocatable :: work(:)
  integer :: lwork, liwork
  real(DP), allocatable :: rwork(:)
  integer :: lrwork
  real(DP) :: en_t(desc(3))
  integer, allocatable :: iwork(:)
  integer :: nfound, nzfound

  character :: range_char
  real(DP), parameter :: orfac = 5d-7
  real(DP), allocatable :: gap(:)
  integer, allocatable :: ifail(:), iclustr(:)

  PUSH_SUB(diagonalize_scalapack_x_cplx)

  iunit_ = 6
  if (present(iunit)) iunit_ = iunit
  ctxt = desc(2)
  if (desc(3)/=desc(4)) call die('Matrix is not square')
  N = desc(4)
  Mb = desc(5)
  Nb = desc(6)
  if (Neig == N) then
    range_char = 'A'
  else
    range_char = 'I'
  endif

  if (ctxt>=0) then

    call blacs_gridinfo(ctxt, nprow, npcol, myprow, mypcol)
    Ml = NUMROC(N, Mb, myprow, 0, nprow) ! # of local rows I own
    Nl = NUMROC(N, Nb, mypcol, 0, npcol) ! # of local columns I own

    nfound = Neig
    nzfound = Neig
    call blacs_pinfo(inode, npes)
    if (inode==0) then
      write(iunit_,'(1x,2a)') 'Using diagonalization routine: PZHEEVX'
      FLUSH(iunit_)
    endif

    ! FHJ: Call in query mode
    SAFE_ALLOCATE(iclustr, (max(1,2*nprow*npcol)))
    SAFE_ALLOCATE(gap, (max(1,nprow*npcol)))
    SAFE_ALLOCATE(ifail, (N))

    SAFE_ALLOCATE(work, (10))
    SAFE_ALLOCATE(rwork, (10))
    lrwork = -1
    SAFE_ALLOCATE(iwork, (10))

    lwork = -1
    liwork = -1

    call scalapack_call()

    call check_lapack_ev_result('PZHEEVX', info, Neig, is_query=.true.)
    call realloc_query_workspace(work, lwork, scalapack_mem_rescale)
    ! Need to add (CLUSTERSIZE-1)*N to rwork to make eigenvectors orthogonal
    ! Assume that CLUSTERSIZE ~ 10
    rwork(1) = rwork(1) + 10d0*N
    call realloc_query_workspace(rwork, lrwork, scalapack_mem_rescale)
    call realloc_query_workspace(iwork, liwork, scalapack_mem_rescale)

    ! FHJ: Call scalapack for real
    if (inode==0) then
      write(iunit_,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', N
      FLUSH(iunit_)
    endif

    call scalapack_call()

    if (inode==0) then
      write(iunit_,'(1x,a)') 'Done ScaLAPACK diagonalization'
      FLUSH(iunit_)
    endif
    call check_lapack_ev_result('PZHEEVX', info, Neig, nfound, nzfound)

    SAFE_DEALLOCATE(work)
    SAFE_DEALLOCATE(rwork)
    SAFE_DEALLOCATE(iwork)

    ! FHJ: Copy eigenvalues/eigenvectors
    eigenval(1:Neig) = en_t(1:Neig)
    SAFE_DEALLOCATE(iclustr)
    SAFE_DEALLOCATE(gap)
    SAFE_DEALLOCATE(ifail)
  endif

  POP_SUB(diagonalize_scalapack_x_cplx)

contains

  subroutine scalapack_call

    PUSH_SUB(diagonalize_scalapack_x_cplx.scalapack_call)

    call PZHEEVX &
      ('V', range_char, 'U', N, matrix, 1, 1, desc, &
      0d0, 0d0, 1, Neig, 0d0, nfound, &
      nzfound, en_t, orfac, eigen_vect, 1, 1, desc, work, lwork, &
      rwork, lrwork, &
      iwork, liwork, ifail, iclustr, gap, info)

    POP_SUB(diagonalize_scalapack_x_cplx.scalapack_call)

  end subroutine scalapack_call

end subroutine diagonalize_scalapack_x_cplx




!-------------------------------------------------------------------------------
! MR3 algorithm
!-------------------------------------------------------------------------------



#ifdef USEMR3
!> Diagonalize a matrix with ScaLAPACK using the MR3 algorithm.
!! Implementation for real type.
subroutine diagonalize_scalapack_r_real(desc, Neig, matrix, eigenval, eigen_vect, &
                                        iunit)
  integer :: desc(9) !< ScaLAPACK/BLACS descriptor. See documentation in DESCINIT
  integer, intent(in) :: Neig !< Number of eigenvalues/eigenvectors to compute
  real(DP), intent(inout) :: matrix(*) !< Local buffer for the 2D distributed matrix
  real(DP), intent(out) :: eigenval(*) !< Eigenvalues. Only the first Neig will be written.
  real(DP), intent(inout) :: eigen_vect(*) !< Local buffer for the 2D distributed eigenvectors
  integer, intent(in), optional :: iunit !< Unit to write output. Defaults to 6.

  ! BLACS stuff
  integer :: ctxt, N, Mb, Nb
  integer :: nprow, npcol, myprow, mypcol, Ml, Nl

  real(DP), parameter :: scalapack_mem_rescale = 1d0
  integer :: info, inode, npes, iunit_

  real(DP), allocatable :: work(:)
  integer :: lwork, liwork
  real(DP) :: en_t(desc(3))
  integer, allocatable :: iwork(:)
  integer :: nfound, nzfound

  character :: range_char

  PUSH_SUB(diagonalize_scalapack_r_real)

  iunit_ = 6
  if (present(iunit)) iunit_ = iunit
  ctxt = desc(2)
  if (desc(3)/=desc(4)) call die('Matrix is not square')
  N = desc(4)
  Mb = desc(5)
  Nb = desc(6)
  if (Neig == N) then
    range_char = 'A'
  else
    range_char = 'I'
  endif

  if (ctxt>=0) then

    call blacs_gridinfo(ctxt, nprow, npcol, myprow, mypcol)
    Ml = NUMROC(N, Mb, myprow, 0, nprow) ! # of local rows I own
    Nl = NUMROC(N, Nb, mypcol, 0, npcol) ! # of local columns I own

    nfound = Neig
    nzfound = Neig
    call blacs_pinfo(inode, npes)
    if (inode==0) then
      write(iunit_,'(1x,2a)') 'Using diagonalization routine: PDSYEVR'
      FLUSH(iunit_)
    endif

    ! FHJ: Call in query mode

    SAFE_ALLOCATE(work, (10))
    SAFE_ALLOCATE(iwork, (10))

    lwork = -1
    liwork = -1

    call scalapack_call()

    call check_lapack_ev_result('PDSYEVR', info, Neig, is_query=.true.)
    call realloc_query_workspace(work, lwork, scalapack_mem_rescale)
    call realloc_query_workspace(iwork, liwork, scalapack_mem_rescale)

    ! FHJ: Call scalapack for real
    if (inode==0) then
      write(iunit_,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', N
      FLUSH(iunit_)
    endif

    call scalapack_call()

    if (inode==0) then
      write(iunit_,'(1x,a)') 'Done ScaLAPACK diagonalization'
      FLUSH(iunit_)
    endif
    call check_lapack_ev_result('PDSYEVR', info, Neig, nfound, nzfound)

    SAFE_DEALLOCATE(work)
    SAFE_DEALLOCATE(iwork)

    ! FHJ: Copy eigenvalues/eigenvectors
    eigenval(1:Neig) = en_t(1:Neig)
  endif

  POP_SUB(diagonalize_scalapack_r_real)

contains

  subroutine scalapack_call

    PUSH_SUB(diagonalize_scalapack_r_real.scalapack_call)

    call PDSYEVR &
      ('V', range_char, 'U', N, matrix, 1, 1, desc, &
      0d0, 0d0, 1, Neig, 0d0, nfound, &
      nzfound, en_t, eigen_vect, 1, 1, desc, work, lwork, &
      iwork, liwork, info)

    POP_SUB(diagonalize_scalapack_r_real.scalapack_call)

  end subroutine scalapack_call

end subroutine diagonalize_scalapack_r_real

#endif


#ifdef USEMR3
!> Diagonalize a matrix with ScaLAPACK using the MR3 algorithm.
!! Implementation for cplx type.
subroutine diagonalize_scalapack_r_cplx(desc, Neig, matrix, eigenval, eigen_vect, &
                                        iunit)
  integer :: desc(9) !< ScaLAPACK/BLACS descriptor. See documentation in DESCINIT
  integer, intent(in) :: Neig !< Number of eigenvalues/eigenvectors to compute
  complex(DPC), intent(inout) :: matrix(*) !< Local buffer for the 2D distributed matrix
  real(DP), intent(out) :: eigenval(*) !< Eigenvalues. Only the first Neig will be written.
  complex(DPC), intent(inout) :: eigen_vect(*) !< Local buffer for the 2D distributed eigenvectors
  integer, intent(in), optional :: iunit !< Unit to write output. Defaults to 6.

  ! BLACS stuff
  integer :: ctxt, N, Mb, Nb
  integer :: nprow, npcol, myprow, mypcol, Ml, Nl

  real(DP), parameter :: scalapack_mem_rescale = 1d0
  integer :: info, inode, npes, iunit_

  complex(DPC), allocatable :: work(:)
  integer :: lwork, liwork
  real(DP), allocatable :: rwork(:)
  integer :: lrwork
  real(DP) :: en_t(desc(3))
  integer, allocatable :: iwork(:)
  integer :: nfound, nzfound

  character :: range_char

  PUSH_SUB(diagonalize_scalapack_r_cplx)

  iunit_ = 6
  if (present(iunit)) iunit_ = iunit
  ctxt = desc(2)
  if (desc(3)/=desc(4)) call die('Matrix is not square')
  N = desc(4)
  Mb = desc(5)
  Nb = desc(6)
  if (Neig == N) then
    range_char = 'A'
  else
    range_char = 'I'
  endif

  if (ctxt>=0) then

    call blacs_gridinfo(ctxt, nprow, npcol, myprow, mypcol)
    Ml = NUMROC(N, Mb, myprow, 0, nprow) ! # of local rows I own
    Nl = NUMROC(N, Nb, mypcol, 0, npcol) ! # of local columns I own

    nfound = Neig
    nzfound = Neig
    call blacs_pinfo(inode, npes)
    if (inode==0) then
      write(iunit_,'(1x,2a)') 'Using diagonalization routine: PZHEEVR'
      FLUSH(iunit_)
    endif

    ! FHJ: Call in query mode

    SAFE_ALLOCATE(work, (10))
    SAFE_ALLOCATE(rwork, (10))
    lrwork = -1
    SAFE_ALLOCATE(iwork, (10))

    lwork = -1
    liwork = -1

    call scalapack_call()

    call check_lapack_ev_result('PZHEEVR', info, Neig, is_query=.true.)
    call realloc_query_workspace(work, lwork, scalapack_mem_rescale)
    call realloc_query_workspace(rwork, lrwork, scalapack_mem_rescale)
    call realloc_query_workspace(iwork, liwork, scalapack_mem_rescale)

    ! FHJ: Call scalapack for real
    if (inode==0) then
      write(iunit_,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', N
      FLUSH(iunit_)
    endif

    call scalapack_call()

    if (inode==0) then
      write(iunit_,'(1x,a)') 'Done ScaLAPACK diagonalization'
      FLUSH(iunit_)
    endif
    call check_lapack_ev_result('PZHEEVR', info, Neig, nfound, nzfound)

    SAFE_DEALLOCATE(work)
    SAFE_DEALLOCATE(rwork)
    SAFE_DEALLOCATE(iwork)

    ! FHJ: Copy eigenvalues/eigenvectors
    eigenval(1:Neig) = en_t(1:Neig)
  endif

  POP_SUB(diagonalize_scalapack_r_cplx)

contains

  subroutine scalapack_call

    PUSH_SUB(diagonalize_scalapack_r_cplx.scalapack_call)

    call PZHEEVR &
      ('V', range_char, 'U', N, matrix, 1, 1, desc, &
      0d0, 0d0, 1, Neig, 0d0, nfound, &
      nzfound, en_t, eigen_vect, 1, 1, desc, work, lwork, &
      rwork, lrwork, &
      iwork, liwork, info)

    POP_SUB(diagonalize_scalapack_r_cplx.scalapack_call)

  end subroutine scalapack_call

end subroutine diagonalize_scalapack_r_cplx

#endif



!-------------------------------------------------------------------------------
! divide & conquer algorithm
!-------------------------------------------------------------------------------



!> Diagonalize a matrix with ScaLAPACK using the divide & conquer algorithm.
!! Implementation for real type.
subroutine diagonalize_scalapack_d_real(desc, Neig, matrix, eigenval, eigen_vect, &
                                        iunit)
  integer :: desc(9) !< ScaLAPACK/BLACS descriptor. See documentation in DESCINIT
  integer, intent(in) :: Neig !< Number of eigenvalues/eigenvectors to compute
  real(DP), intent(inout) :: matrix(*) !< Local buffer for the 2D distributed matrix
  real(DP), intent(out) :: eigenval(*) !< Eigenvalues. Only the first Neig will be written.
  real(DP), intent(inout) :: eigen_vect(*) !< Local buffer for the 2D distributed eigenvectors
  integer, intent(in), optional :: iunit !< Unit to write output. Defaults to 6.

  ! BLACS stuff
  integer :: ctxt, N, Mb, Nb
  integer :: nprow, npcol, myprow, mypcol, Ml, Nl

  real(DP), parameter :: scalapack_mem_rescale = 1d0
  integer :: info, inode, npes, iunit_

  real(DP), allocatable :: work(:)
  integer :: lwork, liwork
  real(DP) :: en_t(desc(3))
  integer, allocatable :: iwork(:)
  integer :: nfound, nzfound

  real(DP), allocatable :: eigen_vect_t(:)

  PUSH_SUB(diagonalize_scalapack_d_real)

  iunit_ = 6
  if (present(iunit)) iunit_ = iunit
  ctxt = desc(2)
  if (desc(3)/=desc(4)) call die('Matrix is not square')
  N = desc(4)
  Mb = desc(5)
  Nb = desc(6)

  if (ctxt>=0) then

    call blacs_gridinfo(ctxt, nprow, npcol, myprow, mypcol)
    Ml = NUMROC(N, Mb, myprow, 0, nprow) ! # of local rows I own
    Nl = NUMROC(N, Nb, mypcol, 0, npcol) ! # of local columns I own

    nfound = Neig
    nzfound = Neig
    call blacs_pinfo(inode, npes)
    if (inode==0) then
      write(iunit_,'(1x,2a)') 'Using diagonalization routine: PDSYEVD'
      FLUSH(iunit_)
    endif

    ! FHJ: Call in query mode
    SAFE_ALLOCATE(eigen_vect_t, (max(1,Ml*Nl)))

    SAFE_ALLOCATE(work, (10))
    SAFE_ALLOCATE(iwork, (10))

    lwork = -1
    liwork = -1

    call scalapack_call()

    call check_lapack_ev_result('PDSYEVD', info, Neig, is_query=.true.)
    call realloc_query_workspace(work, lwork, scalapack_mem_rescale)
    call realloc_query_workspace(iwork, liwork, scalapack_mem_rescale)

    ! FHJ: Call scalapack for real
    if (inode==0) then
      write(iunit_,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', N
      FLUSH(iunit_)
    endif

    call scalapack_call()

    if (inode==0) then
      write(iunit_,'(1x,a)') 'Done ScaLAPACK diagonalization'
      FLUSH(iunit_)
    endif
    call check_lapack_ev_result('PDSYEVD', info, Neig, nfound, nzfound)

    SAFE_DEALLOCATE(work)
    SAFE_DEALLOCATE(iwork)

    ! FHJ: Copy eigenvalues/eigenvectors
    eigenval(1:Neig) = en_t(1:Neig)
    ! FHJ: Copy eigenvectors
    call pdgemr2d(N, Neig, eigen_vect_t, 1, 1, desc, eigen_vect, 1, 1, desc, ctxt)
    SAFE_DEALLOCATE(eigen_vect_t)
  endif

  POP_SUB(diagonalize_scalapack_d_real)

contains

  subroutine scalapack_call

    PUSH_SUB(diagonalize_scalapack_d_real.scalapack_call)

    call PDSYEVD &
      ('V', 'U', N, matrix, 1, 1, desc, &
      en_t, eigen_vect_t, 1, 1, desc, work, lwork, &
      iwork, liwork, info)

    POP_SUB(diagonalize_scalapack_d_real.scalapack_call)

  end subroutine scalapack_call

end subroutine diagonalize_scalapack_d_real



!> Diagonalize a matrix with ScaLAPACK using the divide & conquer algorithm.
!! Implementation for cplx type.
subroutine diagonalize_scalapack_d_cplx(desc, Neig, matrix, eigenval, eigen_vect, &
                                        iunit)
  integer :: desc(9) !< ScaLAPACK/BLACS descriptor. See documentation in DESCINIT
  integer, intent(in) :: Neig !< Number of eigenvalues/eigenvectors to compute
  complex(DPC), intent(inout) :: matrix(*) !< Local buffer for the 2D distributed matrix
  real(DP), intent(out) :: eigenval(*) !< Eigenvalues. Only the first Neig will be written.
  complex(DPC), intent(inout) :: eigen_vect(*) !< Local buffer for the 2D distributed eigenvectors
  integer, intent(in), optional :: iunit !< Unit to write output. Defaults to 6.

  ! BLACS stuff
  integer :: ctxt, N, Mb, Nb
  integer :: nprow, npcol, myprow, mypcol, Ml, Nl

  real(DP), parameter :: scalapack_mem_rescale = 1d0
  integer :: info, inode, npes, iunit_

  complex(DPC), allocatable :: work(:)
  integer :: lwork, liwork
  real(DP), allocatable :: rwork(:)
  integer :: lrwork
  real(DP) :: en_t(desc(3))
  integer, allocatable :: iwork(:)
  integer :: nfound, nzfound

  complex(DPC), allocatable :: eigen_vect_t(:)

  PUSH_SUB(diagonalize_scalapack_d_cplx)

  iunit_ = 6
  if (present(iunit)) iunit_ = iunit
  ctxt = desc(2)
  if (desc(3)/=desc(4)) call die('Matrix is not square')
  N = desc(4)
  Mb = desc(5)
  Nb = desc(6)

  if (ctxt>=0) then

    call blacs_gridinfo(ctxt, nprow, npcol, myprow, mypcol)
    Ml = NUMROC(N, Mb, myprow, 0, nprow) ! # of local rows I own
    Nl = NUMROC(N, Nb, mypcol, 0, npcol) ! # of local columns I own

    nfound = Neig
    nzfound = Neig
    call blacs_pinfo(inode, npes)
    if (inode==0) then
      write(iunit_,'(1x,2a)') 'Using diagonalization routine: PZHEEVD'
      FLUSH(iunit_)
    endif

    ! FHJ: Call in query mode
    SAFE_ALLOCATE(eigen_vect_t, (max(1,Ml*Nl)))

    SAFE_ALLOCATE(work, (10))
    SAFE_ALLOCATE(rwork, (10))
    lrwork = -1
    SAFE_ALLOCATE(iwork, (10))

    lwork = -1
    liwork = -1

    call scalapack_call()

    call check_lapack_ev_result('PZHEEVD', info, Neig, is_query=.true.)
    call realloc_query_workspace(work, lwork, scalapack_mem_rescale)
    call realloc_query_workspace(rwork, lrwork, scalapack_mem_rescale)
    call realloc_query_workspace(iwork, liwork, scalapack_mem_rescale)

    ! FHJ: Call scalapack for real
    if (inode==0) then
      write(iunit_,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', N
      FLUSH(iunit_)
    endif

    call scalapack_call()

    if (inode==0) then
      write(iunit_,'(1x,a)') 'Done ScaLAPACK diagonalization'
      FLUSH(iunit_)
    endif
    call check_lapack_ev_result('PZHEEVD', info, Neig, nfound, nzfound)

    SAFE_DEALLOCATE(work)
    SAFE_DEALLOCATE(rwork)
    SAFE_DEALLOCATE(iwork)

    ! FHJ: Copy eigenvalues/eigenvectors
    eigenval(1:Neig) = en_t(1:Neig)
    ! FHJ: Copy eigenvectors
    call pzgemr2d(N, Neig, eigen_vect_t, 1, 1, desc, eigen_vect, 1, 1, desc, ctxt)
    SAFE_DEALLOCATE(eigen_vect_t)
  endif

  POP_SUB(diagonalize_scalapack_d_cplx)

contains

  subroutine scalapack_call

    PUSH_SUB(diagonalize_scalapack_d_cplx.scalapack_call)

    call PZHEEVD &
      ('V', 'U', N, matrix, 1, 1, desc, &
      en_t, eigen_vect_t, 1, 1, desc, work, lwork, &
      rwork, lrwork, &
      iwork, liwork, info)

    POP_SUB(diagonalize_scalapack_d_cplx.scalapack_call)

  end subroutine scalapack_call

end subroutine diagonalize_scalapack_d_cplx





!> High-level wrapper to diagonalize a matrix with ScaLAPACK.
!! This subroutine accepts matricies in arbitrary BLACS distribution,
!! as long as matrix and eigen_vect are distributed in the same way.
!! A 2D block-cyclic layout is recommended. All the parallelization
!! information is stored in the BLACS descriptor `desc`. Optionally, one
!! can specify the algorithm (0=>bissection, 1=>MR3, or 2=>D&C) via `algo`
subroutine diagonalize_scalapack_real(desc, Neig, matrix, eigenval, eigen_vect, &
                                      iunit, algo)
  integer :: desc(9) !< ScaLAPACK/BLACS descriptor. See documentation in DESCINIT
  integer, intent(in) :: Neig !< Number of eigenvalues/eigenvectors to compute
  real(DP), intent(inout) :: matrix(*) !< Local buffer for the 2D distributed matrix
  real(DP), intent(out) :: eigenval(*) !< Eigenvalues. Only the first Neig will be written.
  real(DP), intent(inout) :: eigen_vect(*) !< Local buffer for the 2D distributed eigenvectors
  integer, intent(in), optional :: iunit !< Unit to write output. Defaults to 6.
  integer, intent(in), optional :: algo !< algorithm to use. 0=bissection, 1=MR3, 2=D&C

  character(len=64) :: scalapack_algo_str
  integer :: err, algo_

  PUSH_SUB(diagonalize_scalapack_real)

  algo_ = -1
  if (present(algo)) algo_ = algo
  if (algo_==-1) then
    call get_environment_variable("BGW_SCALAPACK_ALGO", scalapack_algo_str)
    read(scalapack_algo_str,*,iostat=err) algo_
    if (err/=0) algo_ = -1
  endif

#ifndef USEMR3
  ! If we did not compile with MR3 support
  if (algo_==1) algo_=0
#endif

  select case(algo_)
    case (-1,0)
      call diagonalize_scalapack_x_real(desc, Neig, matrix, eigenval, eigen_vect, iunit)
#ifdef USEMR3
    case (1)
      call diagonalize_scalapack_r_real(desc, Neig, matrix, eigenval, eigen_vect, iunit)
#endif
    case (2)
      call diagonalize_scalapack_d_real(desc, Neig, matrix, eigenval, eigen_vect, iunit)
    case default
      call die('Invalid algorithm in diagonalize_scalapack_real')
  end select

  POP_SUB(diagonalize_scalapack_real)

end subroutine diagonalize_scalapack_real



!> High-level wrapper to diagonalize a matrix with ScaLAPACK.
!! This subroutine accepts matricies in arbitrary BLACS distribution,
!! as long as matrix and eigen_vect are distributed in the same way.
!! A 2D block-cyclic layout is recommended. All the parallelization
!! information is stored in the BLACS descriptor `desc`. Optionally, one
!! can specify the algorithm (0=>bissection, 1=>MR3, or 2=>D&C) via `algo`
subroutine diagonalize_scalapack_cplx(desc, Neig, matrix, eigenval, eigen_vect, &
                                      iunit, algo)
  integer :: desc(9) !< ScaLAPACK/BLACS descriptor. See documentation in DESCINIT
  integer, intent(in) :: Neig !< Number of eigenvalues/eigenvectors to compute
  complex(DPC), intent(inout) :: matrix(*) !< Local buffer for the 2D distributed matrix
  real(DP), intent(out) :: eigenval(*) !< Eigenvalues. Only the first Neig will be written.
  complex(DPC), intent(inout) :: eigen_vect(*) !< Local buffer for the 2D distributed eigenvectors
  integer, intent(in), optional :: iunit !< Unit to write output. Defaults to 6.
  integer, intent(in), optional :: algo !< algorithm to use. 0=bissection, 1=MR3, 2=D&C

  character(len=64) :: scalapack_algo_str
  integer :: err, algo_

  PUSH_SUB(diagonalize_scalapack_cplx)

  algo_ = -1
  if (present(algo)) algo_ = algo
  if (algo_==-1) then
    call get_environment_variable("BGW_SCALAPACK_ALGO", scalapack_algo_str)
    read(scalapack_algo_str,*,iostat=err) algo_
    if (err/=0) algo_ = -1
  endif

#ifndef USEMR3
  ! If we did not compile with MR3 support
  if (algo_==1) algo_=0
#endif

  select case(algo_)
    case (-1,0)
      call diagonalize_scalapack_x_cplx(desc, Neig, matrix, eigenval, eigen_vect, iunit)
#ifdef USEMR3
    case (1)
      call diagonalize_scalapack_r_cplx(desc, Neig, matrix, eigenval, eigen_vect, iunit)
#endif
    case (2)
      call diagonalize_scalapack_d_cplx(desc, Neig, matrix, eigenval, eigen_vect, iunit)
    case default
      call die('Invalid algorithm in diagonalize_scalapack_cplx')
  end select

  POP_SUB(diagonalize_scalapack_cplx)

end subroutine diagonalize_scalapack_cplx


#endif

end module scalapack_hl_m
