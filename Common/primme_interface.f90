#include "f_defs.h"
#if defined MPI && !defined USESCALAPACK
  #error ScaLAPACK is required for MPI builds.
#endif

module primme_interface_m

#ifdef USEPRIMME
#ifdef MPI
  use, intrinsic :: iso_c_binding
  use global_m
  use blas_m
  use primme_bindings_m
  use bgw_mpi_m
#ifdef ACC
  use cublas, only: cublasDgemv, cublasZgemv
#endif
  use nvtx_m

  implicit none

  private

    abstract interface
      subroutine matvec_iface(N, Nown, nk, ldx, ldy, x_loc, y_loc, recvcounts, displs)
        import
        integer, intent(in) :: N, Nown, nk, ldx, ldy
        SCALAR, intent(in) :: x_loc(ldx,nk)
        SCALAR, intent(out) :: y_loc(ldy,nk)
        integer, intent(in) :: recvcounts(*), displs(*)
      end subroutine
    end interface

    abstract interface
      subroutine init_evecs_iface(nev, Nown_loc, offset_loc, &
        evecs_loc, num_restart_read)
        import
        integer, intent(in) :: nev, Nown_loc, offset_loc
        SCALAR, intent(inout) :: evecs_loc(:,:)
        integer, intent(out) :: num_restart_read
      end subroutine
    end interface

    type primme_vars_t
      logical :: use_gpu
      logical :: verbose
      integer :: N, Nown
      integer, allocatable :: recvcounts(:), displs(:)
      real(DP), allocatable :: precond_priv(:)
      real(DP) :: ritz_shift
      real(DP) :: ritz_shift_min
      SCALAR, pointer :: matrix_loc(:,:) => NULL()
      procedure(matvec_iface), nopass, pointer :: matvec => NULL()
      procedure(init_evecs_iface), nopass, pointer :: init_evecs => NULL()
    endtype primme_vars_t

    type(primme_vars_t) :: primme_vars
    SCALAR, allocatable :: matrix_loc_gpu(:,:)
    !$acc declare create(matrix_loc_gpu)

  include 'primme_f90.inc'

  public :: &
    solve_primme_implicit, solve_primme_explicit, set_init_evecs_func


contains


subroutine set_init_evecs_func(init_evecs_func)
  procedure(init_evecs_iface) :: init_evecs_func

  PUSH_SUB(set_init_evecs_func)
  primme_vars%init_evecs => init_evecs_func
  POP_SUB(set_init_evecs_func)

end subroutine set_init_evecs_func


subroutine primme_init_params(primme, N, Nown, nev, apply_precond, &
    tol, primme_method, max_basis_size, max_block_size)!, comm, npes, inode)
  class(primme_t), intent(out) :: primme
  integer, intent(in) :: N, Nown, nev
  logical, optional, intent(in) :: apply_precond
  real(DP), optional, intent(in) :: tol
  integer, optional, intent(in) :: primme_method, max_basis_size, max_block_size
  !integer, optional, intent(in) :: comm, npes, inode

  logical :: apply_precond_
  integer :: maxMatvecs=100000, ld
  integer :: printLevel
  integer(primme_int) :: tmp

  PUSH_SUB(primme_init_params)

  call primme%init()

  call primme%set(PRIMME_n, N)
  call primme%set(PRIMME_nLocal, Nown)
  call primme%set(PRIMME_ldevecs, max(1,Nown))

  call primme%set(PRIMME_numEvals, nev)
  call primme%set(PRIMME_eps, get_optional_arg(1d-6, tol))
  call primme%set(PRIMME_maxBasisSize, max(8, get_optional_arg(256, max_basis_size)))
  call primme%set(PRIMME_maxBlockSize, max(4, get_optional_arg(4, max_block_size)))
  apply_precond_ = get_optional_arg(.false., apply_precond)
  if (apply_precond_) then
    call primme%set(PRIMME_correctionParams_precondition, 1)
    call primme%set(PRIMME_applyPreconditioner, precond)
    primme_vars%ritz_shift = 0.1d0 / ryd
    primme_vars%ritz_shift_min = 1d0 / TOL_ZERO / ryd
    !primme_vars%ritz_shift = ip%get_double('ritz_shift', 0.1d0) / ryd
    !primme_vars%ritz_shift_min = ip%get_double('ritz_shift_min', 1d0/TOL_ZERO) / ryd
  else
    call primme%set(PRIMME_correctionParams_precondition, 0)
  endif

  call primme%set(PRIMME_target, PRIMME_smallest)
  printLevel = 3
  if (peinf%verb_high) printLevel = 4
  call primme%set(PRIMME_printLevel, printLevel)
  !call primme%set(PRIMME_maxMatvecs, maxMatvecs)
  !call primme%set(PRIMME_restartingParams_scheme, PRIMME_thick)

  call primme%set(PRIMME_matrixMatvec, matvec_para_wrapper)
  !call primme%set(PRIMME_commInfo, get_default_arg(MPI_COMM_WORLD, comm))
  !call primme%set(PRIMME_numProcs, get_default_arg(peinf%npes, npes))
  !call primme%set(PRIMME_procID, get_default_arg(peinf%inode, inode))
  call primme%set(PRIMME_commInfo, MPI_COMM_WORLD)
  call primme%set(PRIMME_numProcs, peinf%npes)
  call primme%set(PRIMME_procID, peinf%inode)
  call primme%set(PRIMME_globalSumReal, sum_para)
  call primme%set(PRIMME_broadcastReal, bcast_para)
  call primme%set_method(PRIMME_DEFAULT_METHOD)
  !call primme%set_method(PRIMME_DYNAMIC)
  !call primme%set_method(PRIMME_DEFAULT_MIN_MATVECS)

  POP_SUB(primme_init_params)

end subroutine primme_init_params


!==============================================================================
! Parallel routines
!==============================================================================


!> Implement matvec for a dense distributed matrix
subroutine matvec_dense(N, Nown, nk, ldx, ldy, x_loc, y_loc, recvcounts, displs)
  integer, intent(in) :: N, Nown, nk, ldx, ldy
  SCALAR, intent(in) :: x_loc(ldx,nk)
  SCALAR, intent(out) :: y_loc(ldy,nk)
  integer, intent(in) :: recvcounts(*), displs(*)

  integer :: ik, ipe
  SCALAR :: x(N), y(N)

  !print *, 'inode=',peinf%inode, 'N=',N, 'Nown=', Nown, 'ldx,ldy=',ldx,ldy

#ifdef DEBUG
  ! FIXME:
  FLUSH(0)
  FLUSH(6)
  call bgw_barrier()
  do ipe = 0, peinf%npes-1
    if (peinf%inode==ipe) then
      print *, 'matvec_dense: inode=', peinf%inode, 'ldx=', ldx, 'Nown=', Nown
      FLUSH(0)
      FLUSH(6)
    endif
    call bgw_barrier()
  enddo
#endif

  call nvtxStartRange('matvec_dense', 1)

  !call timacc(52,1)
  do ik = 1, nk
    !call timacc(57,1)
    call MPI_Allgatherv(x_loc(:,ik), Nown, MPI_SCALAR, x, recvcounts, displs, &
      MPI_SCALAR, MPI_COMM_WORLD, mpierr)
    !call timacc(57,2)
    if (Nown>0) then
#ifdef ACC
      if (primme_vars%use_gpu) then
        !$acc data copyin(x) copyout(y_loc(1:N,ik))
        !$acc host_data use_device(matrix_loc_gpu,x,y_loc)
#ifdef CPLX
        call cublasZgemv(&
#else
        call cublasDgemv(&
#endif
        'C', N, Nown, ONE, matrix_loc_gpu, &
          N, x, 1, ZERO, y_loc(:,ik), 1)
        !$acc end host_data
        !$acc end data
      else
#endif
        call X(gemv)('C', N, Nown, ONE, primme_vars%matrix_loc, &
          N, x, 1, ZERO, y_loc(:,ik), 1)
#ifdef ACC
      endif
#endif
    else
      y_loc(:,ik) = ZERO
    endif
  enddo
  !call timacc(52,2)
  call nvtxEndRange()

#ifdef DEBUG
  print *, 'matvec_dense.end: inode=', peinf%inode
  FLUSH(0)
  FLUSH(6)
  call bgw_barrier()
#endif

end subroutine matvec_dense


!> Generic matvec function called by PRIMME. Relies on a user-defined matvec
!! stored in p%matvec, which can point to the default matvec_dense.
subroutine matvec_para_wrapper(x, ldx, y, ldy, k, primme_ptr, ierr)
  integer(primme_int), intent(in) :: ldx
  integer(primme_int), intent(in) :: ldy
  SCALAR, intent(in) :: x(ldx,k)
  SCALAR, intent(out) :: y(ldy,k)
  integer, intent(in) :: k
  type(c_ptr), intent(in) :: primme_ptr
  integer, intent(out) :: ierr

  integer :: ik

  ierr = 0
  associate(p=>primme_vars)
  call p%matvec(p%N, p%Nown, k, int(ldx), int(ldy), &
    x(:,:), y(:,:), p%recvcounts, p%displs)
  endassociate

end subroutine matvec_para_wrapper


subroutine precond(x, ldx, y, ldy, k, primme_ptr, ierr)
  integer(primme_int), intent(in) :: ldx, ldy
  SCALAR, intent(in) :: x(ldx,k)
  SCALAR, intent(out) :: y(ldy,k)
  integer, intent(in) :: k
  type(c_ptr), intent(in) :: primme_ptr
  integer, intent(out) :: ierr

  real(DP) :: shift
  integer :: ik

  ! no push/pop, called too often
  ierr = 0
  associate(pv=>primme_vars)
  do ik = 1, k
    call primme_get_prec_shift_f77(primme_ptr, ik, shift)
    ! FHJ: We add an extra shift because Ritz value it an upper bound to the
    ! eigenvalue.
    shift = min(shift - pv%ritz_shift, pv%ritz_shift_min)
    if (peinf%inode==0 .and. peinf%verb_medium) write(6,'(1x,a,i0,a,f0.6,a)') 'ik = ', ik, ' shift = ', shift, ' Ry'
    y(:,ik) = x(:,ik) / (pv%precond_priv(:) - shift)
  enddo
  endassociate
  ! no push/pop, called too often

end subroutine precond


subroutine sum_para(x, y, k, primme_ptr, ierr)
  integer, intent(in) :: k
  real(DP), target :: x(k)
  real(DP), target :: y(k)
  type(c_ptr), intent(in) :: primme_ptr
  integer, intent(out) :: ierr

  real(kind=8), allocatable :: y_loc(:)
  integer :: k_max, k_min

  k_max = k
  call MPI_Allreduce(MPI_IN_PLACE, k_max, 1, MPI_INTEGER, MPI_MAX, &
                     MPI_COMM_WORLD, ierr)
  k_min = k
  call MPI_Allreduce(MPI_IN_PLACE, k_min, 1, MPI_INTEGER, MPI_MIN, &
                     MPI_COMM_WORLD, ierr)

  if ( k_max /= k_min ) then
    SAFE_ALLOCATE ( y_loc, ( k_max ) )
    y_loc = 0.0D+00
    if (c_associated(c_loc(x),c_loc(y))) then
      y_loc(1:k) = y(1:k)
    else
      y_loc(1:k) = x(1:k)
    end if
    call MPI_Allreduce(MPI_IN_PLACE, y_loc, k_max, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    y(1:k) = y_loc(1:k)
    SAFE_DEALLOCATE(y_loc)
  else

    if (c_associated(c_loc(x),c_loc(y))) then
      call MPI_Allreduce(MPI_IN_PLACE, y, k, MPI_DOUBLE_PRECISION, MPI_SUM, &
        MPI_COMM_WORLD, ierr)
    else
      call MPI_Allreduce(x, y, k, MPI_DOUBLE_PRECISION, MPI_SUM, &
        MPI_COMM_WORLD, ierr)
    endif

  end if

  !XX MDB Original buggy implementation: if k has not the same value across all processes
  !XX ! no push/pop, called too often
  !XX !call timacc(56,1)
  !XX if (c_associated(c_loc(x),c_loc(y))) then
  !XX   call MPI_Allreduce(MPI_IN_PLACE, y, k, MPI_DOUBLE_PRECISION, MPI_SUM, &
  !XX     MPI_COMM_WORLD, ierr)
  !XX else
  !XX   call MPI_Allreduce(x, y, k, MPI_DOUBLE_PRECISION, MPI_SUM, &
  !XX     MPI_COMM_WORLD, ierr)
  !XX endif
  !XX !call timacc(56,2)
  !XX ! no push/pop, called too often

end subroutine sum_para


subroutine bcast_para(buf, count, primme_ptr, ierr)
  real(DP), intent(inout) :: buf(*)
  integer, intent(in) :: count
  type(c_ptr), intent(in) :: primme_ptr
  integer, intent(out) :: ierr

  ! no push/pop, called too often
  !call MPI_Bcast(buf(1), count, MPI_DOUBLE_PRECISION, 0, comm, mpierr)
  call MPI_Bcast(buf(1), count, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
  ierr = mpierr
  ! no push/pop, called too often

end subroutine bcast_para



!> Use PRIMME to find eigenpairs given a matvec.
!!
!! The input local matrix must be distributed in a block-column BLACS
!! layout. There is also a constraint on the block size.
!!
!! \param N [in] Rank of the matrix
!! \param Nown [in] Number of local columns of the matrix I own.
!! \param nev [in] Number of eigenvectors to compute
!! \param evals [out] (nev) Eigenvalues
!! \param evecs_loc [out] (N,nev) Eigenvectors
!! \param matvec [in] The matrix-vector function. See interface for information.
!! \param mat_diag [in] (N) Diagonal entries of the matrix.
!! \param tol [in] PRIMME tolerance for finding eigenvalue, optional.
!! \param primme_method [in] PRIMME method to use, optional.
!! \param max_basis_size [in] PRIMME maximum basis size, optional.
!! \param max_block_size [in] PRIMME maximum block size for a single matvec, optional.
!! \param ritz_shift [in] Amount by which we shift the ritz values when
!!        applying the preconditioner.
subroutine solve_primme_implicit(N, Nown, nev, evals, evecs_loc, matvec, &
    matrix_diag, tol, primme_method, max_basis_size, max_block_size, ritz_shift, use_gpu, verbose)!, comm, npes, inode)
  integer, intent(in) :: N
  integer, intent(in) :: Nown
  integer, intent(in) :: nev
  real(DP), intent(out) :: evals(nev)
  SCALAR, intent(out) :: evecs_loc(max(1,Nown),nev)
  procedure(matvec_iface) :: matvec
  real(DP), optional, intent(in) :: matrix_diag(max(1,Nown))
  real(DP), optional, intent(in) :: tol
  integer, optional, intent(in) :: primme_method, max_basis_size, max_block_size
  real(DP), optional, intent(in) :: ritz_shift
  logical, optional, intent(in) :: use_gpu
  logical, optional, intent(in) :: verbose
  !integer, optional, intent(in) :: comm, npes, inode

  type(primme_t) :: primme
  integer :: ipe, ieval, num_restart
  real(DP) :: rnorms(nev), t0, t1, w0, w1

  PUSH_SUB(solve_primme_implicit)

  if (peinf%inode==0) write(*,'(/1x,a)') 'Preparing diagonalization with PRIMME'

  primme_vars%verbose = get_optional_arg(.true., verbose)
  primme_vars%N = N
  primme_vars%Nown = Nown
  primme_vars%matvec => matvec
  evecs_loc = ZERO
  evals(:) = 0d0

  associate(p=>primme_vars)
  SAFE_ALLOCATE(p%recvcounts, (peinf%npes))
  SAFE_ALLOCATE(p%displs, (peinf%npes))
  SAFE_ALLOCATE(p%precond_priv, (max(1,Nown)))
  call bgw_allgather([Nown], p%recvcounts)
  p%displs(1) = 0
  do ipe = 1, peinf%npes-1
    p%displs(ipe+1) = p%displs(ipe) + p%recvcounts(ipe)
  enddo
  !print *, 'solve_primme_implicit: inode=', peinf%inode, 'displs=', p%displs, 'recvcounts=', p%recvcounts

  p%precond_priv = 0d0
  if (present(matrix_diag)) then
    if (Nown>0) p%precond_priv(1:Nown) = matrix_diag(1:Nown)
  endif

  num_restart = 0
!subroutine primme_init_params(primme, N, Nown, nev, apply_precond, &
!    tol, primme_method, max_basis_size, max_block_size)!, comm, npes, inode)
  call primme_init_params(primme, N, Nown, nev, &
    present(matrix_diag), tol, primme_method, max_basis_size, max_block_size)!, comm, npes, inode)
  if (associated(p%init_evecs)) then
    call p%init_evecs(nev, Nown, p%displs(peinf%inode+1), &
      evecs_loc, num_restart)
    if (num_restart > 0) then
      if (peinf%inode==0) write(6,*) 'Starting vectors: from file'
      call primme%set(PRIMME_initSize, num_restart)
    endif
  endif
  if (num_restart==0) then
    if (peinf%inode==0) write(6,*) 'Starting vectors: random'
  endif

#if 0
  ! DEBUGGING!
  block
    SCALAR :: x_loc(Nown), y_loc(Nown)
    if (peinf%npes==1) then
      x_loc = [1d0, 2d0, 3d0]
    else
      x_loc = dble(peinf%inode + 1)
    endif
    !x_loc = x_loc / sqrt(1d0 + 2d0**2 + 3d0**2)
    !call matvec_dense(N, Nown, 1, Nown, Nown, x_loc, y_loc, p%recvcounts, p%displs)
    !x_loc = y_loc
    call matvec_dense(N, Nown, 1, Nown, Nown, x_loc, y_loc, p%recvcounts, p%displs)
    print *, peinf%inode, y_loc
  end block
  call bgw_barrier()
  block
    integer, parameter :: k=2
    integer :: i
    real(DP) :: x(k)
    x = [(10*i + peinf%inode, i=1,k)]
    print *, 'A', peinf%inode, x
    call sum_para(x, x, k, c_null_ptr, mpierr)
    print *, 'B', peinf%inode, x
  end block
#endif

  if (peinf%inode==0 .and. primme_vars%verbose) then
    print '(/1x,a/)', 'PRIMME interal parameters:'
    call primme_display_params_f77(primme%ptr)
    print '(/1x,a/)', 'Calling PRIMME solver'
    FLUSH(6)
  endif
  call bgw_barrier()

  call timget(t0, w0)
  call primme%solve_eigs(evals, evecs_loc, rnorms)
  call timget(t1, w1)
  if (peinf%inode==0) print '(/1x,a,f0.3,a/)', &
    'Successfully called PRIMME solver after ', w1-w0, ' s.'
  call bgw_barrier()

  endassociate

  ! ----------------------------------------------------------------
  ! Reporting results
  !call primme_display_stats_f77(primme)
  !call primmetop_get_member_f77(primme, PRIMME_aNorm, aNorm)
  !print*, 'Tolerance used: ',epsil, ' Esimated norm(A):', aNorm

  call primme%free()

  !if (params%kp_inode==0) call close_primme_file(primme%fp, params%cur_ik)

  SAFE_DEALLOCATE(primme_vars%recvcounts)
  SAFE_DEALLOCATE(primme_vars%displs)
  SAFE_DEALLOCATE(primme_vars%precond_priv)

  POP_SUB(solve_primme_implicit)

end subroutine solve_primme_implicit


!> Use PRIMME to find eigenpairs given a dense, distributed matrix.
!!
!! The input local matrix must be distributed in a block-column BLACS
!! layout. There is also a constraint on the block size.
!!
!! \param N [in] Rank of the matrix
!! \param Nown [in] Number of local columns of the matrix I own.
!! \param nev [in] Number of eigenvectors to compute
!! \param evals [out] (nev) Eigenvalues
!! \param evecs_loc [out] (N,nev) Eigenvectors
!! \param matrix_loc [in] (N,Nown) BLACS block-column distributed matrix.
!! \param tol [in] PRIMME tolerance for finding eigenvalue, optional.
!! \param primme_method [in] PRIMME method to use, optional.
!! \param max_basis_size [in] PRIMME maximum basis size, optional.
!! \param max_block_size [in] PRIMME maximum block size for a single matvec, optional.
!! \param ritz_shift [in] Amount by which we shift the ritz values when
!!        applying the preconditioner.
subroutine solve_primme_explicit(N, Nown, nev, evals, evecs_loc, matrix_loc, &
    tol, primme_method, max_basis_size, max_block_size, ritz_shift, use_gpu, verbose)!, comm, npes, inode)
  integer, intent(in) :: N
  integer, intent(in) :: Nown
  integer, intent(in) :: nev
  real(DP), intent(out) :: evals(nev)
  SCALAR, intent(out) :: evecs_loc(max(1,Nown),nev)
  SCALAR, target, intent(in) :: matrix_loc(N,max(1,Nown))
  real(DP), optional, intent(in) :: tol
  integer, optional, intent(in) :: primme_method, max_basis_size, max_block_size
  real(DP), optional, intent(in) :: ritz_shift
  logical, optional, intent(in) :: use_gpu
  logical, optional, intent(in) :: verbose
  !integer, optional, intent(in) :: comm, npes, inode

  integer :: ii, recvcounts(peinf%npes), offset
  real(DP) :: mat_diag(max(1,Nown)), matrix_diag(max(1,Nown))

  PUSH_SUB(solve_primme_explicit)

  call bgw_allgather([Nown], recvcounts)
  if (peinf%inode==0) print '(1x,a,i0,", ",i0)', &
      'PRIMME interface: max./min. number of rows stored by (first, last) PE = ', &
      recvcounts(1), recvcounts(peinf%npes)
  offset = 0
  if (peinf%inode>0) offset = sum(recvcounts(1:peinf%inode))

  matrix_diag(:) = 0d0
  matrix_diag(1:Nown) = [(dble(matrix_loc(ii+offset,ii)), ii=1,Nown)]

  primme_vars%use_gpu = get_optional_arg(peinf%ngpus>0, use_gpu)
  primme_vars%use_gpu = .false.
  if (primme_vars%use_gpu) then
    SAFE_ALLOCATE(matrix_loc_gpu, (size(matrix_loc,1),size(matrix_loc,2)))
    matrix_loc_gpu(:,:) = matrix_loc(:,:)
    !$acc update device(matrix_loc_gpu)
  else
    primme_vars%matrix_loc => matrix_loc
  endif

  call solve_primme_implicit(N, Nown, nev, evals, evecs_loc, matvec_dense, &
    matrix_diag, tol, primme_method, max_basis_size, max_block_size, ritz_shift, use_gpu, verbose)!, comm, npes, inode)

  if (primme_vars%use_gpu) then
    SAFE_DEALLOCATE(matrix_loc_gpu)
  else
    primme_vars%matrix_loc => NULL()
  endif

  POP_SUB(solve_primme_explicit)

end subroutine solve_primme_explicit

#endif
#endif

end module primme_interface_m
