#include "f_defs.h"
#if defined MPI && !defined USESCALAPACK
  #error ScaLAPACK is required for MPI builds.
#endif

module diag_primme_m

#ifdef USEPRIMME

  use iso_c_binding
  use global_m
  use blas_m
  use lapack_m
  use scalapack_m
  use inread_m,       only: pb_params_t
  use primme_m
  use diag_scalapack_m
  use distribution_m, only: distrib_mat_t
#ifdef USEELPA
  use diag_elpa_m
#endif

  implicit none

  private

  SCALAR, allocatable :: ham_d_priv_t(:,:)
  real(DP), allocatable :: precond_priv(:)
  integer :: ngk, ngk_own, offset_g, comm, npes, inode
  integer, allocatable :: recvcounts(:)
  real(DP) :: ritz_shift, ritz_shift_min
  include 'primme_eigs_f90.inc'
  !include 'primme_f77.h'

  public :: &
#ifdef MPI
    diag_primme_para, matvec_para, &
#endif
    diag_primme_serial, matvec_serial, precond

contains

subroutine precond(x, ldx, y, ldy, k, primme_ptr, ierr)
  integer(primme_int), intent(in) :: ldx, ldy
  SCALAR, intent(in) :: x(ldx,*)
  SCALAR, intent(out) :: y(ldy,*)
  integer, intent(in) :: k
  type(c_ptr), intent(in) :: primme_ptr
  integer, intent(out) :: ierr

  real(DP) :: shift
  integer :: ik

  ! no push/pop, called too often
  ierr = 0
  y(:,1:k) = ZERO
  do ik = 1, k
    call primme_get_prec_shift_f77(primme_ptr, ik, shift, ierr)
    ! FHJ: We add an extra shift because Ritz value it an upper bound to the
    ! eigenvalue.
    shift = min(shift - ritz_shift, ritz_shift_min)
    y(1:ngk_own,ik) = x(1:ngk_own,ik) / (precond_priv(1:ngk_own) - shift)
  enddo
  ! no push/pop, called too often

end subroutine precond

subroutine primme_init_params(primme, params, dm_ham, dm_wfn)
  class(primme_t), intent(out) :: primme
  type(pb_params_t), intent(in) :: params
  class(distrib_mat_t), intent(in) :: dm_ham
  class(distrib_mat_t), intent(in) :: dm_wfn

  integer :: printLevel=3, maxMatvecs=100000
  integer :: ierr

  PUSH_SUB(primme_init_params)

  call primme%init()
  ritz_shift = params%primme_ritz_shift
  ritz_shift_min = params%primme_ritz_shift_min

  if (inode==0) then
    call open_primme_file(primme%fp, params%cur_ik, params%kpp%ipool)
    call primme%set(PRIMME_outputFile, primme%fp, ierr)
  endif

  call primme%set(PRIMME_n, dm_ham%M, ierr)
  call primme%set(PRIMME_nLocal, dm_ham%Nown, ierr)
  call primme%set(PRIMME_numEvals, dm_wfn%N, ierr)
  call primme%set(PRIMME_maxBasisSize, params%primme_max_basis_sz, ierr)
  call primme%set(PRIMME_maxBlockSize, params%primme_max_block_sz, ierr)
  call primme%set(PRIMME_correctionParams_precondition, 1, ierr)
  call primme%set(PRIMME_applyPreconditioner, precond, ierr)

  call primme%set(PRIMME_eps, params%primme_tol, ierr)
  call primme%set(PRIMME_target, PRIMME_smallest, ierr)
  call primme%set(PRIMME_printLevel, printLevel, ierr)
  call primme%set(PRIMME_maxMatvecs, maxMatvecs, ierr)
  !call primme%set(PRIMME_restartingParams_scheme, PRIMME_thick, ierr)
  !call primme%set(PRIMME_aNorm, fnorm)

#ifdef MPI
  if (npes>1) then
    call primme%set(PRIMME_matrixMatvec, matvec_para, ierr)
    call primme%set(PRIMME_commInfo, comm, ierr)
    call primme%set(PRIMME_numProcs, npes, ierr)
    call primme%set(PRIMME_procID, inode, ierr)
    call primme%set(PRIMME_globalSumReal, sum_para, ierr)
    call primme%set(PRIMME_broadcastReal, bcast_para, ierr)
  else
#endif
    call primme%set(PRIMME_matrixMatvec, matvec_serial, ierr)
#ifdef MPI
  endif
#endif
  call primme%set_method(params%solver_alg-100, ierr)

  if (ierr/=0) then
    write(0,'(/a,i0,a/)') 'ERROR: Got error message ', ierr, ' when setting PRIMME method.'
    call die('Could not set PRIMME method')
  endif

  POP_SUB(primme_init_params)

end subroutine primme_init_params

!==============================================================================
! Serial routines
!==============================================================================


subroutine matvec_serial(x, ldx, y, ldy, k, primme_ptr, ierr)
  integer(primme_int), intent(in) :: ldx
  integer(primme_int), intent(in) :: ldy
  SCALAR, intent(in) :: x(ldx,*)
  SCALAR, intent(out) :: y(ldy,*)
  integer, intent(in) :: k
  type(c_ptr), intent(in) :: primme_ptr
  integer, intent(out) :: ierr

  integer :: ldx_, ldy_

  ! no push/pop, called too often
  ierr = 0
  ldx_ = int(ldx)
  ldy_ = int(ldy)
  call X(GEMM)('T', 'N', ldx_, k, ldx_, ONE, ham_d_priv_t, ldx_, x, ldx_, ZERO, y, ldy_)
  ! no push/pop, called too often

end subroutine matvec_serial


subroutine diag_primme_serial(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  SCALAR, intent(in) :: ham_d(dm_ham%Ml,dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  SCALAR, intent(out) :: wfn_d(dm_wfn%Ml,dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  type(primme_t) :: primme
  integer :: ierr, ig
  real(DP) :: rnorms(dm_wfn%N)

  PUSH_SUB(diag_primme_serial)

  wfn_d(:,:) = ZERO
  en(:) = 0d0

  SAFE_ALLOCATE(ham_d_priv_t, (dm_ham%Nl,dm_ham%Ml))
  ham_d_priv_t = transpose(ham_d)
  ngk = dm_ham%Ml
  ngk_own = ngk
  SAFE_ALLOCATE(precond_priv, (dm_ham%Ml))
  do ig = 1, ngk_own
    precond_priv(ig) = dble(ham_d_priv_t(ig,ig))
  enddo
  !SAFE_DEALLOCATE(ham_d)
  !SAFE_ALLOCATE(ham_d, (1,1))

  inode = 0
  npes = 1
  call primme_init_params(primme, params, dm_ham, dm_wfn)
  !call primme_display_params_f77(primme)
  call X(primme_f77)(en, wfn_d, rnorms, primme%ptr, ierr)

  ! ----------------------------------------------------------------
  ! Reporting results
  if (inode==0) then
    if (ierr==0) then
      write(params%kpp%iunit,'(1x,a)') 'PRIMME has returned successfully'
      FLUSH(params%kpp%iunit)
    else
       write(0,'(1x,a,i0)') 'ERROR: PRIMME returned with error: ', ierr
       call die('PRIMME returned with error')
    endif
  endif
  !call primme_display_stats_f77(primme)
  !call primmetop_get_member_f77(primme, PRIMME_aNorm, aNorm)
  !print*, 'Tolerance used: ',epsil, ' Esimated norm(A):', aNorm

  call primme%free()
  if (inode==0) call close_primme_file(primme%fp, params%cur_ik)

  SAFE_DEALLOCATE(ham_d_priv_t)
  SAFE_DEALLOCATE(precond_priv)

  POP_SUB(diag_primme_serial)

end subroutine diag_primme_serial


!==============================================================================
! Parallel routines
!==============================================================================

#ifdef MPI

subroutine matvec_para(x, ldx, y, ldy, k, primme_ptr, ierr)
  integer(primme_int), intent(in) :: ldx
  integer(primme_int), intent(in) :: ldy
  SCALAR, intent(in) :: x(ldx,k)
  SCALAR, intent(out) :: y(ldy,k)
  integer, intent(in) :: k
  type(c_ptr), intent(in) :: primme_ptr
  integer, intent(out) :: ierr

  SCALAR :: y2(k,ngk), y3(k,ngk_own), yl(ngk,k)
  integer :: recvs(npes), ldy_min

  ! no push/pop, called too often
  ierr = 0

#if 0
  y = ZERO
  y2 = ZERO
  yl = ZERO
  if (ngk_own>0) then
    !call X(GEMM)('T', 'N', k, ngk, ngk_own, ONE, x, ldx, ham_d_priv_t, ngk_own, ZERO, y2, k)
    call X(GEMM)('T', 'N', ngk, k, ngk_own, ONE, ham_d_priv_t, ngk_own, x, ldx, ZERO, yl, ngk)
  endif
  call MPI_Allreduce(MPI_IN_PLACE, yl, ngk*k, MPI_SCALAR, MPI_SUM, comm, mpierr)
  if (ngk_own>0) then
    !y3(1:k,1:ngk_own) = y2(1:k, offset_g+1:offset_g+ngk_own)
    !y(1:ngk_own,1:k) = transpose(y3)
    y(1:ngk_own, 1:k) = yl(offset_g+1:offset_g+ngk_own, 1:k)
  endif
#else
  ! Want: global y = H x
  ! H_priv => Local cols that I own: Hl = H(1:N,l)
  ! N => Ngk
  ! O => ngk_Own
  ! l: 1..O
  ! [ y_1 ] = [H_1 ... H_O] [ x_1 ]
  ! [ ... ]                 [ ... ]
  ! [ y_N ]                 [ x_O ]
  ! H_priv_t => transpose of H_priv = H(l,1:N)
  ! [yt_1 ... yt_N] = [xt_1 ... xt_O] [ Ht_1 ]
  !                                   [ ...  ]
  !                                   [ Ht_O ]
  ! Apply xl^T * Hl^T to get yl^T, then reduce yl^T to y^T
  ! Work with the transpose form so that we can lump everything into a
  ! contiguous chuck of memory and use a single call to MPI_Reduce_scatter.
  y2(:,:) = ZERO
  y3(:,:) = ZERO
  y(:,1:k) = ZERO
  call X(GEMM)('T', 'N', k, ngk, ngk_own, ONE, x, ldx, ham_d_priv_t, ngk_own, ZERO, y2, k)
  recvs(1:npes) = recvcounts(1:npes)*k
  call MPI_Reduce_scatter(y2, y3, recvs, MPI_SCALAR, MPI_SUM, comm, mpierr)
  y(1:ngk_own,1:k) = transpose(y3)
  ! no push/pop, called too often
#endif

end subroutine matvec_para


subroutine sum_para(sendbuf, recvbuf, count, primme_ptr, ierr)
  real(DP), intent(in), target :: sendbuf(*)
  real(DP), intent(out), target :: recvbuf(*)
  integer, intent(in) :: count
  type(c_ptr), intent(in) :: primme_ptr
  integer, intent(out) :: ierr

  ! no push/pop, called too often
  if (c_associated(c_loc(sendbuf), c_loc(recvbuf))) then
    call MPI_Allreduce(MPI_IN_PLACE, recvBuf(1), count, MPI_DOUBLE_PRECISION, MPI_SUM, comm, mpierr)
  else
    call MPI_Allreduce(sendBuf(1), recvBuf(1), count, MPI_DOUBLE_PRECISION, MPI_SUM, comm, mpierr)
  endif
  ierr = mpierr
  ! no push/pop, called too often

end subroutine sum_para


subroutine bcast_para(buf, count, primme_ptr, ierr)
  real(DP), intent(inout) :: buf(*)
  integer, intent(in) :: count
  type(c_ptr), intent(in) :: primme_ptr
  integer, intent(out) :: ierr

  ! no push/pop, called too often
  call MPI_Bcast(buf(1), count, MPI_DOUBLE_PRECISION, 0, comm, mpierr)
  ierr = mpierr
  ! no push/pop, called too often

end subroutine bcast_para


subroutine diag_primme_para(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  !SCALAR, intent(inout)  :: ham_d(dm_ham%Ml,dm_ham%Nl)
  SCALAR, intent(inout), allocatable :: ham_d(:,:)
  class(distrib_mat_t), intent(in) :: dm_wfn
  SCALAR, intent(out) :: wfn_d(dm_wfn%Ml,dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  type(distrib_mat_t) :: dm_evecs
  SCALAR :: sendbuf(dm_ham%Nb), recvbuf(dm_ham%Nb*params%kpp%npes)
  SCALAR :: evecs(dm_ham%Nl,dm_wfn%N)
  SCALAR :: evecs_guess(dm_ham%Nl,dm_wfn%N)
  !SCALAR, allocatable :: evecs_guess(:,:)
  type(primme_t) :: primme
  integer :: ierr, ig, ig_l, ig_g, ib, ib_l, ipe, init_size
  integer :: color
  real(DP) :: rnorms(dm_wfn%N)

  PUSH_SUB(diag_primme_para)

  ngk = dm_ham%M
  ngk_own = dm_ham%Nown

  inode = params%kpp%inode
  ! PRIMME doesn`t work when a processor doesn`t onn any data. So,
  ! let`s split the MPI communicator and only include the processors with data.
  ! If inode==-1, the processor is not in the group.
  color = 1
  if (ngk_own==0) color = 2
  inode = params%kpp%inode
  call MPI_Comm_split(params%kpp%comm, color, inode, comm, mpierr)
  if (params%kpp%inode==0) call MPI_Comm_size(comm, npes, mpierr)
  call MPI_Bcast(npes, 1, MPI_INTEGER, 0, params%kpp%comm, mpierr)
  if (color==2) inode = -1

  wfn_d(:,:) = ZERO
  en(:) = 0d0

  if (inode>-1) then
    !ASSERT(dm_ham%M==dm_ham%Ml)
    SAFE_ALLOCATE(ham_d_priv_t, (ngk_own,ngk))
    ham_d_priv_t = ZERO
    if (ngk_own > 0) then
      ham_d_priv_t = transpose(ham_d(1:ngk, 1:ngk_own))
    endif

    SAFE_ALLOCATE(precond_priv, (ngk_own))
    precond_priv(:) = 0d0
    offset_g = inode*dm_ham%Nb
    do ig_l = 1, ngk_own
      precond_priv(ig_l) = dble(ham_d_priv_t(ig_l, offset_g+ig_l))
    enddo

    SAFE_ALLOCATE(recvcounts, (params%kpp%npes))
    do ipe = 0, params%kpp%npes-1
      recvcounts(ipe+1) = NUMROC(dm_ham%N, dm_ham%Nb, ipe, 0, params%kpp%npes)
    enddo

    call primme_init_params(primme, params, dm_ham, dm_wfn)
    !call primme_display_params_f77(primme)
  endif

  init_size = 0
  evecs = ZERO
  if (params%primme_exact_diag_size>0) then
    !params%primme_exact_diag_size = min(params%primme_exact_diag_size, dm_ham%M)
    init_size = min(params%primme_exact_diag_size, dm_wfn%N)
    !SAFE_ALLOCATE(evecs_guess, (dm_ham%Nl,init_size))
    call dm_evecs%setup_mat_1d(dm_ham%M, init_size, params%kpp, row_or_col='r')
    call get_guess_evecs(params%primme_exact_diag_size, &
      dm_ham, ham_d, dm_evecs, evecs_guess, en, params)
    evecs(1:dm_ham%Nl,1:init_size) = evecs_guess(1:dm_ham%Nl,1:init_size)
    !evecs = evecs_guess(1:dm_ham%M,1:init_size)
    if (inode>-1) call primme%set(PRIMME_initSize, init_size, ierr)
  endif

   SAFE_DEALLOCATE(ham_d)
   SAFE_ALLOCATE(ham_d, (1,1))

  ! Call PRIMME
  if (inode>-1) call X(primme_f77)(en, evecs, rnorms, primme%ptr, ierr)
  call MPI_Barrier(params%kpp%comm, mpierr)

  sendbuf = ZERO
  do ib = 1, dm_wfn%N
    ipe = INDXG2P(ib, dm_wfn%Nb, 0, 0, params%kpp%npes)
    ib_l = INDXG2L(ib, dm_wfn%Nb, 0, 0, params%kpp%npes)
    sendbuf(1:ngk_own) = evecs(:,ib)
    call MPI_Gather(sendbuf, dm_ham%Nb, MPI_SCALAR, recvbuf, dm_ham%Nb, MPI_SCALAR, ipe, params%kpp%comm, ierr)
    if (ipe==params%kpp%inode) then
      wfn_d(:,ib_l) = recvbuf(1:ngk)
    endif
  enddo

  ! ----------------------------------------------------------------
  ! Reporting results
  if (params%kpp%inode==0) then
    if (ierr==0) then
      write(params%kpp%iunit,'(1x,a)') 'PRIMME has returned successfully'
      FLUSH(params%kpp%iunit)
    else
       write(0,'(1x,a,i0)') 'ERROR: PRIMME returned with error: ', ierr
       call die('PRIMME returned with error.')
    endif
  endif
  if (params%kpp%inode==0) then
    write(params%kpp%iunit,'(1x,a)') 'First eigenvector written to evec.dat:'
    call open_file(file='evec.dat', unit=666, form='formatted', status='replace')
#ifdef CPLX
    write(666,'(1x,f12.6,1x,f12.6)') wfn_d(:,1)
#else
    write(666,'(1x,f12.6)') wfn_d(:,1)
#endif
    call close_file(666)
  endif
  !call primme_display_stats_f77(primme)
  !call primmetop_get_member_f77(primme, PRIMME_aNorm, aNorm)
  !print*, 'Tolerance used: ',epsil, ' Esimated norm(A):', aNorm

  call primme%free()
  if (params%kpp%inode==0) call close_primme_file(primme%fp, params%cur_ik)

  if (inode>-1) then
    SAFE_DEALLOCATE(ham_d_priv_t)
    SAFE_DEALLOCATE(precond_priv)
    SAFE_DEALLOCATE(recvcounts)
  endif

  POP_SUB(diag_primme_para)

end subroutine diag_primme_para


!> Perform diagonalization of the dense Hamiltonian with a smaller cutoff,
!! add a random perturbation to the eigenvectors up to the original cutoff,
!! and perform a QR orthogonalization on the perturbed eigenvectors.
subroutine get_guess_evecs(M_sub, dm_ham, ham_d, dm_evecs, evecs, en, params)
  integer, intent(in) :: M_sub !< size of subspace Hamiltonian
  class(distrib_mat_t), intent(in) :: dm_ham
  SCALAR, intent(in) :: ham_d(dm_ham%Ml, dm_ham%Nl) !< original Hamiltonian, any distrib.
  !> Need dm_evecs%M=dm_ham%M.
  class(distrib_mat_t), intent(in) :: dm_evecs
  !> evecs up to full Ham. size, and including only number of columns that we
  !! care about for the guess size. Typically, dm_evecs%N << M_sub << dm_ham%M.
  SCALAR, intent(out) :: evecs(dm_evecs%Ml,dm_evecs%Nl)
  real(DP), intent(out) :: en(dm_evecs%N)
  type(pb_params_t), intent(in) :: params

  type(distrib_mat_t) :: dm_ham_sub, dm_wfn_sub
  SCALAR, allocatable :: ham_d_sub(:,:), wfn_d_sub(:,:)
  SCALAR  :: evecs_rnd(dm_evecs%Ml,dm_evecs%Nl)
  integer :: iseed(4)

  if (params%kpp%inode==0) then
    write(params%kpp%iunit,'(1x,a,i0)') 'Performing diagonalizing on smaller Hamiltonian with N=', M_sub
    write(params%kpp%iunit,'(1x,a,i0,a)') 'Keeping only lowest ', dm_evecs%N, ' states.'
    FLUSH(params%kpp%iunit)
  endif

  call dm_ham_sub%setup_mat_1d(M_sub, M_sub, params%kpp)
  SAFE_ALLOCATE(ham_d_sub, (dm_ham_sub%Ml, dm_ham_sub%Nl))
  call pX(gemr2d)(dm_ham_sub%M, dm_ham_sub%N, &
    ham_d, 1, 1, dm_ham%desc, &
    ham_d_sub, 1, 1, dm_ham_sub%desc, &
    dm_ham%cntxt)

  call dm_wfn_sub%setup_mat_1d(M_sub, dm_evecs%N, params%kpp)
  SAFE_ALLOCATE(wfn_d_sub, (dm_wfn_sub%Ml, dm_wfn_sub%Nl))

  en(:) = 0d0
#ifdef USEELPA
  call diag_elpa_para(dm_ham_sub, ham_d_sub, dm_wfn_sub, wfn_d_sub, &
    en(1:dm_wfn_sub%N), params)
#else
  call diag_scalapack_para_x(dm_ham_sub, ham_d_sub, dm_wfn_sub, wfn_d_sub, &
    en(1:dm_wfn_sub%N), params)
#endif

  if (params%kpp%inode==0) then
    write(params%kpp%iunit,'(1x,a,f9.3,a,f9.3)') &
      'Eigenvalues of the small Hamiltonian (Ry):', &
      en(1), ' ... ', en(dm_wfn_sub%N)
    FLUSH(params%kpp%iunit)
  endif

  ! Copy eigenvectors into evecs
  evecs(:,:) = ZERO
  call pX(gemr2d)(dm_wfn_sub%M, dm_wfn_sub%N, &
    wfn_d_sub, 1, 1, dm_wfn_sub%desc, &
    evecs, 1, 1, dm_evecs%desc, &
    dm_ham%cntxt)

  ! Add random perturbation
  iseed(:) = [1,2,3,4]
  iseed(1) = peinf%inode
  call zlarnv(4, iseed, dm_evecs%Ml*dm_evecs%Nl, evecs_rnd(1,1))
  evecs(:,:) = evecs(:,:) + evecs_rnd(:,:) / dm_evecs%M * 1d-4

  ! Re-orthogonalize evecs
  call do_qr(dm_evecs, evecs, params)

  SAFE_DEALLOCATE(wfn_d_sub)
  SAFE_DEALLOCATE(ham_d_sub)

end subroutine get_guess_evecs


!> Performs a QR orthogonalization of the matrix `mat` without pivotting.
subroutine do_qr(dm_mat, mat, params)
  class(distrib_mat_t), intent(in) :: dm_mat
  SCALAR, intent(inout) :: mat(dm_mat%Ml, dm_mat%Nl)
  type(pb_params_t), intent(in) :: params

  SCALAR :: tau(min(dm_mat%M,dm_mat%N)), work_tmp(4)
  SCALAR, allocatable :: work(:)
  integer :: lwork, info

  lwork = -1
  call pX(geqrf)(dm_mat%M, dm_mat%N, mat, 1, 1, dm_mat%desc, tau, work_tmp, lwork, info)
  if (info/=0) call die('Got info/=0 for pXgeqrf in query mode', only_root_writes=.true.)
  lwork = int(work_tmp(1))
  SAFE_ALLOCATE(work, (lwork))
  call pX(geqrf)(dm_mat%M, dm_mat%N, mat, 1, 1, dm_mat%desc, tau, work, lwork, info)
  if (info/=0) call die('Got info/=0 for pXgeqrf', only_root_writes=.true.)
  SAFE_DEALLOCATE(work)

#ifdef CPLX
#define PXGQR pzungqr
#else
#define PXGQR pdorgqr
#endif
  lwork = -1
  call PXGQR(dm_mat%M, dm_mat%N, dm_mat%N, mat, 1, 1, dm_mat%desc, tau, work_tmp, lwork, info)
  if (info/=0) call die('Got info/=0 for pzungqr in query mode', only_root_writes=.true.)
  lwork = int(work_tmp(1))
  SAFE_ALLOCATE(work, (lwork))
  call PXGQR(dm_mat%M, dm_mat%N, dm_mat%N, mat, 1, 1, dm_mat%desc, tau, work, lwork, info)
  if (info/=0) call die('Got info/=0 for pzungqr', only_root_writes=.true.)
  SAFE_DEALLOCATE(work)

end subroutine do_qr

#endif

#endif

end module diag_primme_m
