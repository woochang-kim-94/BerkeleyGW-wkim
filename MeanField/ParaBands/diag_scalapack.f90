#include "f_defs.h"
#if defined MPI && !defined USESCALAPACK
  #error ScaLAPACK is required for MPI builds.
#endif
#ifdef CPLX
#define XHE(x) ZHE ## x
#define pXHE(x) PZHE ## x
#else
#define XHE(x) DSY ## x
#define pXHE(x) PDSY ## x
#endif

module diag_scalapack_m

  use global_m
  use blas_m
  use lapack_m
  use scalapack_m
  use inread_m,       only: pb_params_t
  use distribution_m, only: distrib_mat_t

  implicit none

  private

  interface realloc_query_workspace
    module procedure realloc_query_workspace_int, realloc_query_workspace_real, realloc_query_workspace_cplx
  end interface realloc_query_workspace

  public :: &
#ifdef MPI
  diag_scalapack_para_x, diag_scalapack_para_d, diag_scalapack_para_r, &
#endif
  diag_scalapack_serial_x, diag_scalapack_serial_d, diag_scalapack_serial_r

contains


subroutine get_lapack_ev_routname(suffix, routname, params)
  character(len=*), intent(in) :: suffix
  character(len=64), intent(out) :: routname
  type(pb_params_t), intent(in) :: params

  character(len=1)  :: tmpstr

  PUSH_SUB(get_lapack_ev_routname)

  tmpstr = ''
  if (params%kpp%npes>1) write(tmpstr,'(a)') 'P'
#ifdef CPLX
  write(routname,'(a)') trim(tmpstr)//'ZHEEV'//trim(suffix)
#else
  write(routname,'(a)') trim(tmpstr)//'DSYEV'//trim(suffix)
#endif

  POP_SUB(get_lapack_ev_routname)

end subroutine get_lapack_ev_routname


subroutine print_lapack_ev_routine(suffix, params)
  character(len=*), intent(in) :: suffix
  type(pb_params_t), intent(in) :: params

  character(len=64) :: routname

  PUSH_SUB(print_lapack_ev_routine)

  call get_lapack_ev_routname(suffix, routname, params)
  if (params%kpp%inode==0) then
    write(params%kpp%iunit,'(1x,2a)') 'Using diagonalization routine: ', trim(routname)
    FLUSH(params%kpp%iunit)
  endif

  POP_SUB(print_lapack_ev_routine)

end subroutine print_lapack_ev_routine


subroutine check_lapack_ev_result(params, suffix, info, nwant, nfound, nzfound, is_query)
  type(pb_params_t), intent(in) :: params
  character(len=*), intent(in) :: suffix
  integer, optional, intent(in) :: info
  integer, optional, intent(in) :: nwant
  integer, optional, intent(in) :: nfound
  integer, optional, intent(in) :: nzfound
  logical, optional, intent(in) :: is_query

  character(len=64) :: tmpname, routname
  character(len=128) :: errmsg

  PUSH_SUB(check_lapack_ev_result)

  call get_lapack_ev_routname(suffix, tmpname, params)
  if (present(is_query)) then
    if (is_query) write(routname,'(a)') trim(tmpname)//' (query mode)'
  else
    routname = tmpname
  endif
  if (info==2 .and. trim(routname)=='PZHEEVX') then
    if (params%kpp%inode==0) then
      write(0,'(/,a)') 'WARNING: got info=2 from call to PZHEEVX'
      write(0,'(a)') '(failed to reorthogonalize eigenvectors).'
      write(0,'(a)') 'You may want to rerun BerkeleyGW with a different ScaLAPACK solver'
      write(0,'(a/)') 'using the environment variable BGW_SCALAPACK_ALGO'
    endif
  elseif (info/=0) then
    if (params%kpp%inode==0) then
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

  POP_SUB(check_lapack_ev_result)

end subroutine check_lapack_ev_result


subroutine safe2int(real_in, int_out, rescaling)
  real(DP), intent(in) :: real_in
  integer, intent(out) :: int_out
  real(DP), intent(in), optional :: rescaling

  real(DP) :: x

  PUSH_SUB(safe2int)

  !call get_environment_variable("BGW_SCALAPACK_MEM_RESCALE", scalapack_mem_rescale_str)
  !read(scalapack_mem_rescale_str,*,iostat=err) scalapack_mem_rescale
  !if (err/=0) then
  !  scalapack_mem_rescale = 1d0
  !endif
  x = real_in
  if (present(rescaling)) x = x * rescaling
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




!==============================================================================
! Serial routines
!==============================================================================


subroutine diag_scalapack_serial_x(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  SCALAR, intent(in) :: ham_d(dm_ham%Ml,dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  SCALAR, intent(out) :: wfn_d(dm_wfn%Ml,dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  SCALAR, allocatable :: work(:)
  integer :: iwork(5*dm_ham%M)
  integer :: lwork, ifail(dm_ham%M), nfound, info
#ifdef CPLX
  real(DP) :: rwork(7*dm_ham%M)
#endif
  real(DP) :: en_t(dm_ham%M)

  PUSH_SUB(diag_scalapack_serial_x)

  call print_lapack_ev_routine('X', params)
  wfn_d(:,:) = ZERO
  en(:) = 0d0

  ! Query mode.
  SAFE_ALLOCATE(work, (10))
  lwork = -1
  call XHE(evx)('V', 'I', 'L', dm_ham%M, ham_d, dm_ham%M, 0d0, 0d0, 1, dm_wfn%N, 0d0, &
    nfound, en_t, wfn_d, dm_wfn%M, work, lwork, &
#ifdef CPLX
    rwork, &
#endif
   iwork, ifail, info)
  call check_lapack_ev_result(params, 'X', info, dm_wfn%N, is_query=.true.)

  call realloc_query_workspace(work, lwork, params%scalapack_mem_rescale)
  ! Call for realz
  call XHE(evx)('V', 'I', 'L', dm_ham%M, ham_d, dm_ham%M, 0d0, 0d0, 1, dm_wfn%N, 0d0, &
    nfound, en_t, wfn_d, dm_wfn%M, work, lwork, &
#ifdef CPLX
    rwork, &
#endif
   iwork, ifail, info)
  call check_lapack_ev_result(params, 'X', info, dm_wfn%N, nfound)

  SAFE_DEALLOCATE(work)

  en(1:dm_wfn%N) = en_t(1:dm_wfn%N)

  POP_SUB(diag_scalapack_serial_x)

end subroutine diag_scalapack_serial_x


subroutine diag_scalapack_serial_d(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  SCALAR, intent(in) :: ham_d(dm_ham%Ml,dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  SCALAR, intent(out) :: wfn_d(dm_wfn%Ml,dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  SCALAR, allocatable :: work(:)
  integer, allocatable :: iwork(:)
  integer :: lwork, liwork
#ifdef CPLX
  real(DP), allocatable :: rwork(:)
  integer :: lrwork
#endif
  real(DP) :: en_t(dm_ham%M)
  integer :: info

  PUSH_SUB(diag_scalapack_serial_d)

  call print_lapack_ev_routine('D', params)
  wfn_d(:,:) = ZERO
  en(:) = 0d0

  ! Query mode.
  SAFE_ALLOCATE(work, (10))
  SAFE_ALLOCATE(iwork, (10))
#ifdef CPLX
  SAFE_ALLOCATE(rwork, (10))
  lrwork = -1
#endif
  lwork = -1
  liwork = -1
  call XHE(evd)('V', 'L', dm_ham%M, ham_d, dm_ham%M, en_t, work, lwork, &
#ifdef CPLX
    rwork, lrwork, &
#endif
   iwork, liwork, info)
  call check_lapack_ev_result(params, 'D', info, dm_wfn%N, is_query=.true.)

  call realloc_query_workspace(work, lwork, params%scalapack_mem_rescale)
  call realloc_query_workspace(iwork, liwork, params%scalapack_mem_rescale)
#ifdef CPLX
  call realloc_query_workspace(rwork, lrwork, params%scalapack_mem_rescale)
#endif
  ! Call for realz
  call XHE(evd)('V', 'L', dm_ham%M, ham_d, dm_ham%M, en_t, work, lwork, &
#ifdef CPLX
    rwork, lrwork, &
#endif
   iwork, liwork, info)
  call check_lapack_ev_result(params, 'D', info, dm_wfn%N)

  SAFE_DEALLOCATE(work)
  SAFE_DEALLOCATE(iwork)
#ifdef CPLX
  SAFE_DEALLOCATE(rwork)
#endif

  ! Note: ham might be shorter than wfn (M_psi<=M_wfn), but the wfn might
  ! be thinner than ham (N_psi>=N_wfn)
  wfn_d(1:dm_ham%Ml, 1:dm_wfn%Nl) = ham_d(1:dm_ham%Ml, 1:dm_wfn%Nl)
  en(1:dm_wfn%N) = en_t(1:dm_wfn%N)

  POP_SUB(diag_scalapack_serial_d)

end subroutine diag_scalapack_serial_d


subroutine diag_scalapack_serial_r(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  SCALAR, intent(in) :: ham_d(dm_ham%Ml,dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  SCALAR, intent(out) :: wfn_d(dm_wfn%Ml,dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  SCALAR, allocatable :: work(:)
  integer, allocatable :: iwork(:)
  integer :: lwork, liwork
#ifdef CPLX
  real(DP), allocatable :: rwork(:)
  integer :: lrwork
#endif
  real(DP) :: en_t(dm_ham%M)
  integer :: isuppz(2*max(1,dm_wfn%N)), info, nfound

  PUSH_SUB(diag_scalapack_serial_r)

  call print_lapack_ev_routine('R', params)
  wfn_d(:,:) = ZERO
  en(:) = 0d0

  ! Query mode.
  SAFE_ALLOCATE(work, (10))
  SAFE_ALLOCATE(iwork, (10))
#ifdef CPLX
  SAFE_ALLOCATE(rwork, (10))
  lrwork = -1
#endif
  lwork = -1
  liwork = -1
  call XHE(evr)('V', 'I', 'L', dm_ham%M, ham_d, dm_ham%M, 0d0, 0d0, 1, dm_wfn%N, 0d0, &
    nfound, en_t, wfn_d, dm_wfn%M, isuppz, work, lwork, &
#ifdef CPLX
    rwork, lrwork, &
#endif
   iwork, liwork, info)
  call check_lapack_ev_result(params, 'R', info, dm_wfn%N, is_query=.true.)

  call realloc_query_workspace(work, lwork, params%scalapack_mem_rescale)
  call realloc_query_workspace(iwork, liwork, params%scalapack_mem_rescale)
#ifdef CPLX
  call realloc_query_workspace(rwork, lrwork, params%scalapack_mem_rescale)
#endif
  ! Call for realz
  call XHE(evr)('V', 'I', 'L', dm_ham%M, ham_d, dm_ham%M, 0d0, 0d0, 1, dm_wfn%N, 0d0, &
    nfound, en_t, wfn_d, dm_wfn%M, isuppz, work, lwork, &
#ifdef CPLX
    rwork, lrwork, &
#endif
   iwork, liwork, info)
  call check_lapack_ev_result(params, 'R', info, dm_wfn%N, nfound)

  SAFE_DEALLOCATE(work)
  SAFE_DEALLOCATE(iwork)
#ifdef CPLX
  SAFE_DEALLOCATE(rwork)
#endif

  en(1:dm_wfn%N) = en_t(1:dm_wfn%N)

  POP_SUB(diag_scalapack_serial_r)

end subroutine diag_scalapack_serial_r


!==============================================================================
! Parallel routines
!==============================================================================


#ifdef MPI
subroutine diag_scalapack_para_x(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  SCALAR, intent(inout) :: ham_d(dm_ham%Ml, dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  SCALAR, intent(out) :: wfn_d(dm_wfn%Ml, dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  type(distrib_mat_t) :: dm_2d
  SCALAR, allocatable :: work(:)
  integer :: lwork, liwork
#ifdef CPLX
  real(DP), allocatable :: rwork(:)
  integer :: lrwork
#endif
  SCALAR, allocatable :: ham_2d(:), wfn_2d(:)
  integer :: info, nfound, nzfound
  real(DP) :: en_t(dm_ham%M), orfac
  real(DP), allocatable :: gap(:)
  integer, allocatable :: ifail(:), iclustr(:), iwork(:)

  PUSH_SUB(diag_scalapack_para_x)

  call print_lapack_ev_routine('X', params)
  wfn_d(:,:) = ZERO
  en(:) = 0d0

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices ham_2d, wfn_2d.
  ! Some processors might get excluded for the sake of load balancing.
  call dm_2d%setup_mat_2d(dm_ham%M, dm_ham%N, params%block_sz, params%kpp)

  SAFE_ALLOCATE(ham_2d, (dm_2d%Ml * dm_2d%Nl))
  ham_2d(:) = ZERO
  SAFE_ALLOCATE(wfn_2d, (dm_2d%Ml * dm_2d%Nl))
  wfn_2d(:) = ZERO

  ! FHJ: Copy matrix ham_d from 1d cyclic layout into ham_2d (2d cyclic layout)
  call pX(gemr2d)(dm_ham%M, dm_ham%N, ham_d, 1, 1, dm_ham%desc, &
    ham_2d, 1, 1, dm_2d%desc, dm_ham%cntxt)

  if (dm_2d%cntxt>=0) then
    ! FHJ: Call p*evx in query mode
    orfac = 5d-7
    SAFE_ALLOCATE(iclustr, (2*dm_2d%nprow * dm_2d%npcol))
    SAFE_ALLOCATE(gap, (dm_2d%nprow * dm_2d%npcol))
    SAFE_ALLOCATE(ifail, (dm_2d%N))
    SAFE_ALLOCATE(work, (10))
#ifdef CPLX
    SAFE_ALLOCATE(rwork, (10))
    lrwork = -1
#endif
    SAFE_ALLOCATE(iwork, (10))

    lwork = -1
    liwork = -1
    call pXHE(evx) &
      ('V', 'I', 'L', dm_2d%N, ham_2d, 1, 1, dm_2d%desc, &
      0d0, 0d0, 1, dm_wfn%N, 0d0, nfound, &
      nzfound, en_t, orfac, wfn_2d, 1, 1, dm_2d%desc, work, lwork, &
#ifdef CPLX
      rwork, lrwork, &
#endif
      iwork, liwork, ifail, iclustr, gap, info)
    call check_lapack_ev_result(params, 'X', info, dm_wfn%N, is_query=.true.)

    call realloc_query_workspace(work, lwork, params%scalapack_mem_rescale)
#ifdef CPLX
    ! Need to add (CLUSTERSIZE-1)*N to rwork to make eigenvectors orthogonal
    ! Assume that CLUSTERSIZE ~ 10
    rwork(1) = rwork(1) + 10d0*dm_2d%N
    call realloc_query_workspace(rwork, lrwork, params%scalapack_mem_rescale)
#endif
    call realloc_query_workspace(iwork, liwork, params%scalapack_mem_rescale)

    ! FHJ: Call p*evx for realz
    if (params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', dm_2d%N
      FLUSH(params%kpp%iunit)
    endif
    call pXHE(evx) &
      ('V', 'I', 'L', dm_2d%N, ham_2d, 1, 1, dm_2d%desc, &
      0d0, 0d0, 1, dm_wfn%N, 0d0, nfound, &
      nzfound, en_t, orfac, wfn_2d, 1, 1, dm_2d%desc, work, lwork, &
#ifdef CPLX
      rwork, lrwork, &
#endif
      iwork, liwork, ifail, iclustr, gap, info)
    if (params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a)') 'Done ScaLAPACK diagonalization'
      FLUSH(params%kpp%iunit)
    endif
    call check_lapack_ev_result(params, 'X', info, dm_wfn%N, nfound, nzfound)

    SAFE_DEALLOCATE(work)
#ifdef CPLX
    SAFE_DEALLOCATE(rwork)
#endif
    SAFE_DEALLOCATE(iwork)

    ! FHJ: Copy eigenvalues/eigenvectors
    en(1:dm_wfn%N) = en_t(1:dm_wfn%N)
    SAFE_DEALLOCATE(iclustr)
    SAFE_DEALLOCATE(gap)
    SAFE_DEALLOCATE(ifail)
  endif !dm_2d%cntxt>=0

  call MPI_Barrier(params%kpp%comm, mpierr)
  ! FHJ: PEs outside the grid don`t have the evals!
  call MPI_Bcast(en, size(en), MPI_DOUBLE_PRECISION, 0, params%kpp%comm, mpierr)
  ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
  call logit('Redistributing wavefunctions', params%kpp%inode==0, params%kpp%iunit)
  call pX(gemr2d)(dm_wfn%M, dm_wfn%N, wfn_2d, 1, 1, dm_2d%desc, wfn_d, 1, 1, dm_wfn%desc, dm_wfn%cntxt)

  ! FHJ: Cleanup
  SAFE_DEALLOCATE(ham_2d)
  SAFE_DEALLOCATE(wfn_2d)
  call MPI_Barrier(params%kpp%comm, mpierr)

  POP_SUB(diag_scalapack_para_x)

end subroutine diag_scalapack_para_x


subroutine diag_scalapack_para_d(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  SCALAR, intent(inout) :: ham_d(dm_ham%Ml, dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  SCALAR, intent(out) :: wfn_d(dm_wfn%Ml, dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  type(distrib_mat_t) :: dm_2d
  SCALAR, allocatable :: work(:)
  integer :: lwork, liwork
#ifdef CPLX
  real(DP), allocatable :: rwork(:)
  integer :: lrwork
#endif
  SCALAR, allocatable :: ham_2d(:), wfn_2d(:)
  integer :: info
  real(DP) :: en_t(dm_ham%M)
  integer, allocatable ::  iwork(:)

  PUSH_SUB(diag_scalapack_para_d)

  call print_lapack_ev_routine('D', params)
  wfn_d(:,:) = ZERO
  en(:) = 0d0

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices ham_2d, wfn_2d.
  ! Some processors might get excluded for the sake of load balancing.
  call dm_2d%setup_mat_2d(dm_ham%M, dm_ham%N, params%block_sz, params%kpp)

  SAFE_ALLOCATE(ham_2d, (dm_2d%Ml * dm_2d%Nl))
  ham_2d(:) = ZERO
  SAFE_ALLOCATE(wfn_2d, (dm_2d%Ml * dm_2d%Nl))
  wfn_2d(:) = ZERO

  ! FHJ: Copy matrix ham_d from 1d cyclic layout into ham_2d (2d cyclic layout)
  call pX(gemr2d)(dm_ham%M, dm_ham%N, ham_d, 1, 1, dm_ham%desc, &
    ham_2d, 1, 1, dm_2d%desc, dm_ham%cntxt)

  if (dm_2d%cntxt>=0) then
    ! FHJ: Call p*evx in query mode
    SAFE_ALLOCATE(work, (10))
#ifdef CPLX
    SAFE_ALLOCATE(rwork, (10))
    lrwork = -1
#endif
    SAFE_ALLOCATE(iwork, (10))

    lwork = -1
    liwork = -1
    call pXHE(evd) &
      ('V', 'L', dm_2d%N, ham_2d, 1, 1, dm_2d%desc, &
      en_t, wfn_2d, 1, 1, dm_2d%desc, work, lwork, &
#ifdef CPLX
      rwork, lrwork, &
#endif
      iwork, liwork, info)
    call check_lapack_ev_result(params, 'D', info, dm_wfn%N, is_query=.true.)

    call realloc_query_workspace(work, lwork, params%scalapack_mem_rescale)
#ifdef CPLX
    call realloc_query_workspace(rwork, lrwork, params%scalapack_mem_rescale)
#endif
    call realloc_query_workspace(iwork, liwork, params%scalapack_mem_rescale)

    ! FHJ: Call p*evx for realz
    if (params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', dm_2d%N
      FLUSH(params%kpp%iunit)
    endif
    call pXHE(evd) &
      ('V', 'L', dm_2d%N, ham_2d, 1, 1, dm_2d%desc, &
      en_t, wfn_2d, 1, 1, dm_2d%desc, work, lwork, &
#ifdef CPLX
      rwork, lrwork, &
#endif
      iwork, liwork, info)
    if (params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a)') 'Done ScaLAPACK diagonalization'
      FLUSH(params%kpp%iunit)
    endif
    call check_lapack_ev_result(params, 'D', info, dm_wfn%N)

    SAFE_DEALLOCATE(work)
#ifdef CPLX
    SAFE_DEALLOCATE(rwork)
#endif
    SAFE_DEALLOCATE(iwork)

    ! FHJ: Copy eigenvalues/eigenvectors
    en(1:dm_wfn%N) = en_t(1:dm_wfn%N)
  endif !dm_2d%cntxt>=0

  call MPI_Barrier(params%kpp%comm, mpierr)
  ! FHJ: PEs outside the grid don`t have the evals!
  call MPI_Bcast(en, size(en), MPI_DOUBLE_PRECISION, 0, params%kpp%comm, mpierr)
  ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
  call logit('Redistributing wavefunctions', params%kpp%inode==0, params%kpp%iunit)
  call pX(gemr2d)(dm_wfn%M, dm_wfn%N, wfn_2d, 1, 1, dm_2d%desc, wfn_d, 1, 1, dm_wfn%desc, dm_wfn%cntxt)

  ! FHJ: Cleanup
  SAFE_DEALLOCATE(ham_2d)
  SAFE_DEALLOCATE(wfn_2d)
  call MPI_Barrier(params%kpp%comm, mpierr)

  POP_SUB(diag_scalapack_para_d)

end subroutine diag_scalapack_para_d


subroutine diag_scalapack_para_r(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  SCALAR, intent(inout) :: ham_d(dm_ham%Ml, dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  SCALAR, intent(out) :: wfn_d(dm_wfn%Ml, dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  type(distrib_mat_t) :: dm_2d
  SCALAR, allocatable :: work(:)
  integer :: lwork, liwork
#ifdef CPLX
  real(DP), allocatable :: rwork(:)
  integer :: lrwork
#endif
  SCALAR, allocatable :: ham_2d(:), wfn_2d(:)
  integer :: info, nfound, nzfound
  real(DP) :: en_t(dm_ham%M)
  integer, allocatable ::  iwork(:)

  PUSH_SUB(diag_scalapack_para_r)

#ifdef USEMR3
  call print_lapack_ev_routine('R', params)
  wfn_d(:,:) = ZERO
  en(:) = 0d0

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices ham_2d, wfn_2d.
  ! Some processors might get excluded for the sake of load balancing.
  call dm_2d%setup_mat_2d(dm_ham%M, dm_ham%N, params%block_sz, params%kpp)

  SAFE_ALLOCATE(ham_2d, (dm_2d%Ml * dm_2d%Nl))
  ham_2d(:) = ZERO
  SAFE_ALLOCATE(wfn_2d, (dm_2d%Ml * dm_2d%Nl))
  wfn_2d(:) = ZERO

  ! FHJ: Copy matrix ham_d from 1d cyclic layout into ham_2d (2d cyclic layout)
  call pX(gemr2d)(dm_ham%M, dm_ham%N, ham_d, 1, 1, dm_ham%desc, &
    ham_2d, 1, 1, dm_2d%desc, dm_ham%cntxt)

  if (dm_2d%cntxt>=0) then
    ! FHJ: Call p*evx in query mode
    SAFE_ALLOCATE(work, (10))
#ifdef CPLX
    SAFE_ALLOCATE(rwork, (10))
    lrwork = -1
#endif
    SAFE_ALLOCATE(iwork, (10))

    lwork = -1
    liwork = -1
    call pXHE(evr) &
      ('V', 'I', 'L', dm_2d%N, ham_2d, 1, 1, dm_2d%desc, &
      0d0, 0d0, 1, dm_wfn%N, nfound, &
      nzfound, en_t, wfn_2d, 1, 1, dm_2d%desc, work, lwork, &
#ifdef CPLX
      rwork, lrwork, &
#endif
      iwork, liwork, info)
    call check_lapack_ev_result(params, 'R', info, dm_wfn%N, is_query=.true.)

    call realloc_query_workspace(work, lwork, params%scalapack_mem_rescale)
#ifdef CPLX
    call realloc_query_workspace(rwork, lrwork, params%scalapack_mem_rescale)
#endif
    call realloc_query_workspace(iwork, liwork, params%scalapack_mem_rescale)

    ! FHJ: Call p*evx for realz
    if (params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', dm_2d%N
      FLUSH(params%kpp%iunit)
    endif
    call pXHE(evr) &
      ('V', 'I', 'L', dm_2d%N, ham_2d, 1, 1, dm_2d%desc, &
      0d0, 0d0, 1, dm_wfn%N, nfound, &
      nzfound, en_t, wfn_2d, 1, 1, dm_2d%desc, work, lwork, &
#ifdef CPLX
      rwork, lrwork, &
#endif
      iwork, liwork, info)
    if (params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a)') 'Done ScaLAPACK diagonalization'
      FLUSH(params%kpp%iunit)
    endif
    call check_lapack_ev_result(params, 'R', info, dm_wfn%N, nfound)

    SAFE_DEALLOCATE(work)
#ifdef CPLX
    SAFE_DEALLOCATE(rwork)
#endif
    SAFE_DEALLOCATE(iwork)

    ! FHJ: Copy eigenvalues/eigenvectors
    en(1:dm_wfn%N) = en_t(1:dm_wfn%N)
  endif !dm_2d%cntxt>=0

  call MPI_Barrier(params%kpp%comm, mpierr)
  ! FHJ: PEs outside the grid don`t have the evals!
  call MPI_Bcast(en, size(en), MPI_DOUBLE_PRECISION, 0, params%kpp%comm, mpierr)
  ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
  call logit('Redistributing wavefunctions', params%kpp%inode==0, params%kpp%iunit)
  call pX(gemr2d)(dm_wfn%M, dm_wfn%N, wfn_2d, 1, 1, dm_2d%desc, wfn_d, 1, 1, dm_wfn%desc, dm_wfn%cntxt)

  ! FHJ: Cleanup
  SAFE_DEALLOCATE(ham_2d)
  SAFE_DEALLOCATE(wfn_2d)
  call MPI_Barrier(params%kpp%comm, mpierr)
#else
  call diag_scalapack_para_x(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
#endif

  POP_SUB(diag_scalapack_para_r)

end subroutine diag_scalapack_para_r
#endif

end module diag_scalapack_m
