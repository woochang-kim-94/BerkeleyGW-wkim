#include "f_defs.h"
#if defined MPI && !defined USESCALAPACK
  #error ScaLAPACK is required for MPI builds.
#endif
#if defined USEELPA && !defined MPI
  #error MPI is required for ELPA builds.
#endif


module diag_elpa_m
#ifdef USEELPA

  use, intrinsic :: iso_c_binding
  use global_m
  use scalapack_m
  use inread_m,       only: pb_params_t
  use distribution_m, only: distrib_mat_t
  use elpa,           only: elpa_t, elpa_init, elpa_ok, elpa_allocate, &
                            elpa_deallocate, elpa_uninit, elpa_initialized, &
                            ELPA_SOLVER_2STAGE
  implicit none

  private

  interface diag_elpa_para
    module procedure diag_elpa_para_s, diag_elpa_para_d, diag_elpa_para_c, diag_elpa_para_z
  end interface diag_elpa_para

  public :: &
    diag_elpa_para

contains


subroutine assert_elpa_ok(elpa_err)
  integer, intent(in) :: elpa_err
  character(len=128) :: tmpstr

  PUSH_SUB(assert_elpa_ok)

  if (elpa_err/=elpa_ok) then
    write(tmpstr,'(a,i0)') 'Elpa assertion error: ', elpa_err
    call die(tmpstr)
  endif

  POP_SUB(assert_elpa_ok)

end subroutine assert_elpa_ok



subroutine diag_elpa_para_s(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  real, intent(inout) :: ham_d(dm_ham%Ml, dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  real, intent(out) :: wfn_d(dm_wfn%Ml, dm_wfn%Nl)
  real, intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  type(distrib_mat_t) :: dm_2d
  real, allocatable :: ham_2d(:,:), wfn_2d(:,:)
  real :: en_t(dm_ham%M)
  integer :: comm_rows, comm_cols, ranks(dm_ham%nprow * dm_ham%npcol)
  integer :: ipe, group_size, orig_group, active_group, active_comm
  logical :: should_include(params%kpp%npes)
  class(elpa_t), pointer :: elpa
  integer :: elpa_err

  PUSH_SUB(diag_elpa_para_s)

  if (params%kpp%inode==0) then
    write(params%kpp%iunit,'(1x,a)') 'Solving with ELPA'
    FLUSH(params%kpp%iunit)
  endif

  !if (.not.elpa_initialized()) then
  !  if (elpa_init(20170403)/=elpa_ok) call die("ELPA API version not supported")
  !endif
  if (elpa_init(20170403)/=elpa_ok) call die("ELPA API version not supported")

  wfn_d(:,:) = ZERO
  en(:) = 0d0

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices ham_2d, wfn_2d.
  ! Some processors might get excluded for the sake of load balancing.
  call dm_2d%setup_mat_2d(dm_ham%M, dm_ham%N, params%block_sz, params%kpp)

  ! FHJ: Note that the distributed matrix wfn_2d is dm_ham%M x dm%ham%N = Hamiltonian size.
  ! ELPA requires this. However, we only copy the first dm_wfn%N columns from wfn_2d to
  ! wfn_d at the end.
  SAFE_ALLOCATE(ham_2d, (dm_2d%Ml, dm_2d%Nl))
  ham_2d = ZERO
  SAFE_ALLOCATE(wfn_2d, (dm_2d%Ml, dm_2d%Nl))
  wfn_2d = ZERO
  ! FHJ: Copy matrix ham_d from 1d cyclic layout into ham_2d (2d cyclic layout)
  call psgemr2d(dm_ham%M, dm_ham%N, ham_d, 1, 1, dm_ham%desc, &
    ham_2d, 1, 1, dm_2d%desc, dm_ham%cntxt)

  ! FHJ: Create an MPI Group and Communicator for processors inside the Blacs grid
  call logit('Determining active processors', params%kpp%inode==0, params%kpp%iunit)
  should_include(:) = .false.
  should_include(params%kpp%inode+1) = dm_2d%cntxt>=0
  call MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, should_include, &
    1, MPI_LOGICAL, params%kpp%comm, mpierr)

  group_size = 0
  do ipe = 0, params%kpp%npes-1
    if (should_include(ipe+1)) then
      group_size = group_size + 1
      ranks(group_size) = ipe
    endif
  enddo
  if (dm_2d%cntxt>=0 .and. (group_size /= dm_2d%nprow*dm_2d%npcol)) then
    write(0,'(/a)') 'ERROR: Inconsistent processors in group'
    write(0,'(2(a,i0)/)') 'group_size=', group_size, &
      ' nprow*npcol=', dm_2d%nprow*dm_2d%npcol
    call die('Inconsistent processors in group')
  endif

  call logit('Creating communicators', params%kpp%inode==0, params%kpp%iunit)
  call MPI_Comm_group(params%kpp%comm, orig_group, mpierr)
  call MPI_Group_incl(orig_group, group_size, ranks, active_group, mpierr)
  call MPI_Comm_create(params%kpp%comm, active_group, active_comm, mpierr)
  if ((dm_2d%cntxt>=0) .neqv. (active_comm/=MPI_COMM_NULL)) then
    write(0,'(/a)') 'ERROR: Inconsistent group created:'
    write(0,'(4(a,i0)/)') &
      'kp%inode=', params%kpp%inode, ' group_size=', group_size, &
      ' dm_2d%cntxt=', dm_2d%cntxt, ' active_comm=', active_comm
    call die('Inconsistent group')
  endif

  if (dm_2d%cntxt>=0) then
    ! FHJ: Call ELPA
    if (params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a,i0)') 'Beginning ELPA diagonalization. Size: ', dm_2d%N
      FLUSH(params%kpp%iunit)
    endif

    call logit('Calling elpa_allocate')
    elpa => elpa_allocate()
    call logit('ok')

    call elpa%set("na", dm_2d%N, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("nev", dm_wfn%N, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("local_nrows", dm_2d%Ml, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("local_ncols", dm_2d%Nl, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("nblk", dm_2d%Nb, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("mpi_comm_parent", active_comm, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("process_row", dm_2d%myprow, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("process_col", dm_2d%mypcol, elpa_err)
    call assert_elpa_ok(elpa_err)

    call assert_elpa_ok(elpa%setup())
    call elpa%set("solver", ELPA_SOLVER_2STAGE, elpa_err)
    call assert_elpa_ok(elpa_err)

    call logit('Calling elpa solver')
    call elpa%eigenvectors(ham_2d, en_t, wfn_2d, elpa_err)
    call assert_elpa_ok(elpa_err)
    call logit('ok')
    call MPI_Barrier(active_comm, mpierr)

    call logit('Calling elpa_deallocate')
    call elpa_deallocate(elpa)
    call logit('ok')
    call MPI_Barrier(active_comm, mpierr)

    ! FHJ: Copy eigenvalues/eigenvectors
    en(1:dm_wfn%N) = en_t(1:dm_wfn%N)

    call MPI_Comm_free(active_comm, mpierr)
    call MPI_Group_free(active_group, mpierr)
  endif !dm_2d%cntxt>=0

  call MPI_Barrier(params%kpp%comm, mpierr)
  ! FHJ: PEs outside the grid don`t have the evals!
  call MPI_Bcast(en, size(en), MPI_DOUBLE_PRECISION, 0, params%kpp%comm, mpierr)
  ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
  call logit('Redistributing wavefunctions', params%kpp%inode==0, params%kpp%iunit)
  call psgemr2d(dm_wfn%M, dm_wfn%N, wfn_2d, 1, 1, dm_2d%desc, wfn_d, 1, 1, dm_wfn%desc, dm_wfn%cntxt)

  ! FHJ: Cleanup
  SAFE_DEALLOCATE(ham_2d)
  SAFE_DEALLOCATE(wfn_2d)
  call MPI_Group_free(orig_group, mpierr)
  call elpa_uninit()
  call MPI_Barrier(params%kpp%comm, mpierr)

  POP_SUB(diag_elpa_para_s)

end subroutine diag_elpa_para_s

subroutine diag_elpa_para_d(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  real(DP), intent(inout) :: ham_d(dm_ham%Ml, dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  real(DP), intent(out) :: wfn_d(dm_wfn%Ml, dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  type(distrib_mat_t) :: dm_2d
  real(DP), allocatable :: ham_2d(:,:), wfn_2d(:,:)
  real(DP) :: en_t(dm_ham%M)
  integer :: comm_rows, comm_cols, ranks(dm_ham%nprow * dm_ham%npcol)
  integer :: ipe, group_size, orig_group, active_group, active_comm
  logical :: should_include(params%kpp%npes)
  class(elpa_t), pointer :: elpa
  integer :: elpa_err

  PUSH_SUB(diag_elpa_para_d)

  if (params%kpp%inode==0) then
    write(params%kpp%iunit,'(1x,a)') 'Solving with ELPA'
    FLUSH(params%kpp%iunit)
  endif

  !if (.not.elpa_initialized()) then
  !  if (elpa_init(20170403)/=elpa_ok) call die("ELPA API version not supported")
  !endif
  if (elpa_init(20170403)/=elpa_ok) call die("ELPA API version not supported")

  wfn_d(:,:) = ZERO
  en(:) = 0d0

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices ham_2d, wfn_2d.
  ! Some processors might get excluded for the sake of load balancing.
  call dm_2d%setup_mat_2d(dm_ham%M, dm_ham%N, params%block_sz, params%kpp)

  ! FHJ: Note that the distributed matrix wfn_2d is dm_ham%M x dm%ham%N = Hamiltonian size.
  ! ELPA requires this. However, we only copy the first dm_wfn%N columns from wfn_2d to
  ! wfn_d at the end.
  SAFE_ALLOCATE(ham_2d, (dm_2d%Ml, dm_2d%Nl))
  ham_2d = ZERO
  SAFE_ALLOCATE(wfn_2d, (dm_2d%Ml, dm_2d%Nl))
  wfn_2d = ZERO
  ! FHJ: Copy matrix ham_d from 1d cyclic layout into ham_2d (2d cyclic layout)
  call pdgemr2d(dm_ham%M, dm_ham%N, ham_d, 1, 1, dm_ham%desc, &
    ham_2d, 1, 1, dm_2d%desc, dm_ham%cntxt)

  ! FHJ: Create an MPI Group and Communicator for processors inside the Blacs grid
  call logit('Determining active processors', params%kpp%inode==0, params%kpp%iunit)
  should_include(:) = .false.
  should_include(params%kpp%inode+1) = dm_2d%cntxt>=0
  call MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, should_include, &
    1, MPI_LOGICAL, params%kpp%comm, mpierr)

  group_size = 0
  do ipe = 0, params%kpp%npes-1
    if (should_include(ipe+1)) then
      group_size = group_size + 1
      ranks(group_size) = ipe
    endif
  enddo
  if (dm_2d%cntxt>=0 .and. (group_size /= dm_2d%nprow*dm_2d%npcol)) then
    write(0,'(/a)') 'ERROR: Inconsistent processors in group'
    write(0,'(2(a,i0)/)') 'group_size=', group_size, &
      ' nprow*npcol=', dm_2d%nprow*dm_2d%npcol
    call die('Inconsistent processors in group')
  endif

  call logit('Creating communicators', params%kpp%inode==0, params%kpp%iunit)
  call MPI_Comm_group(params%kpp%comm, orig_group, mpierr)
  call MPI_Group_incl(orig_group, group_size, ranks, active_group, mpierr)
  call MPI_Comm_create(params%kpp%comm, active_group, active_comm, mpierr)
  if ((dm_2d%cntxt>=0) .neqv. (active_comm/=MPI_COMM_NULL)) then
    write(0,'(/a)') 'ERROR: Inconsistent group created:'
    write(0,'(4(a,i0)/)') &
      'kp%inode=', params%kpp%inode, ' group_size=', group_size, &
      ' dm_2d%cntxt=', dm_2d%cntxt, ' active_comm=', active_comm
    call die('Inconsistent group')
  endif

  if (dm_2d%cntxt>=0) then
    ! FHJ: Call ELPA
    if (params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a,i0)') 'Beginning ELPA diagonalization. Size: ', dm_2d%N
      FLUSH(params%kpp%iunit)
    endif

    call logit('Calling elpa_allocate')
    elpa => elpa_allocate()
    call logit('ok')

    call elpa%set("na", dm_2d%N, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("nev", dm_wfn%N, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("local_nrows", dm_2d%Ml, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("local_ncols", dm_2d%Nl, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("nblk", dm_2d%Nb, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("mpi_comm_parent", active_comm, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("process_row", dm_2d%myprow, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("process_col", dm_2d%mypcol, elpa_err)
    call assert_elpa_ok(elpa_err)

    call assert_elpa_ok(elpa%setup())
    call elpa%set("solver", ELPA_SOLVER_2STAGE, elpa_err)
    call assert_elpa_ok(elpa_err)

    call logit('Calling elpa solver')
    call elpa%eigenvectors(ham_2d, en_t, wfn_2d, elpa_err)
    call assert_elpa_ok(elpa_err)
    call logit('ok')
    call MPI_Barrier(active_comm, mpierr)

    call logit('Calling elpa_deallocate')
    call elpa_deallocate(elpa)
    call logit('ok')
    call MPI_Barrier(active_comm, mpierr)

    ! FHJ: Copy eigenvalues/eigenvectors
    en(1:dm_wfn%N) = en_t(1:dm_wfn%N)

    call MPI_Comm_free(active_comm, mpierr)
    call MPI_Group_free(active_group, mpierr)
  endif !dm_2d%cntxt>=0

  call MPI_Barrier(params%kpp%comm, mpierr)
  ! FHJ: PEs outside the grid don`t have the evals!
  call MPI_Bcast(en, size(en), MPI_DOUBLE_PRECISION, 0, params%kpp%comm, mpierr)
  ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
  call logit('Redistributing wavefunctions', params%kpp%inode==0, params%kpp%iunit)
  call pdgemr2d(dm_wfn%M, dm_wfn%N, wfn_2d, 1, 1, dm_2d%desc, wfn_d, 1, 1, dm_wfn%desc, dm_wfn%cntxt)

  ! FHJ: Cleanup
  SAFE_DEALLOCATE(ham_2d)
  SAFE_DEALLOCATE(wfn_2d)
  call MPI_Group_free(orig_group, mpierr)
  call elpa_uninit()
  call MPI_Barrier(params%kpp%comm, mpierr)

  POP_SUB(diag_elpa_para_d)

end subroutine diag_elpa_para_d

subroutine diag_elpa_para_c(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  complex, intent(inout) :: ham_d(dm_ham%Ml, dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  complex, intent(out) :: wfn_d(dm_wfn%Ml, dm_wfn%Nl)
  real, intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  type(distrib_mat_t) :: dm_2d
  complex, allocatable :: ham_2d(:,:), wfn_2d(:,:)
  real :: en_t(dm_ham%M)
  integer :: comm_rows, comm_cols, ranks(dm_ham%nprow * dm_ham%npcol)
  integer :: ipe, group_size, orig_group, active_group, active_comm
  logical :: should_include(params%kpp%npes)
  class(elpa_t), pointer :: elpa
  integer :: elpa_err

  PUSH_SUB(diag_elpa_para_c)

  if (params%kpp%inode==0) then
    write(params%kpp%iunit,'(1x,a)') 'Solving with ELPA'
    FLUSH(params%kpp%iunit)
  endif

  !if (.not.elpa_initialized()) then
  !  if (elpa_init(20170403)/=elpa_ok) call die("ELPA API version not supported")
  !endif
  if (elpa_init(20170403)/=elpa_ok) call die("ELPA API version not supported")

  wfn_d(:,:) = ZERO
  en(:) = 0d0

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices ham_2d, wfn_2d.
  ! Some processors might get excluded for the sake of load balancing.
  call dm_2d%setup_mat_2d(dm_ham%M, dm_ham%N, params%block_sz, params%kpp)

  ! FHJ: Note that the distributed matrix wfn_2d is dm_ham%M x dm%ham%N = Hamiltonian size.
  ! ELPA requires this. However, we only copy the first dm_wfn%N columns from wfn_2d to
  ! wfn_d at the end.
  SAFE_ALLOCATE(ham_2d, (dm_2d%Ml, dm_2d%Nl))
  ham_2d = ZERO
  SAFE_ALLOCATE(wfn_2d, (dm_2d%Ml, dm_2d%Nl))
  wfn_2d = ZERO
  ! FHJ: Copy matrix ham_d from 1d cyclic layout into ham_2d (2d cyclic layout)
  call pcgemr2d(dm_ham%M, dm_ham%N, ham_d, 1, 1, dm_ham%desc, &
    ham_2d, 1, 1, dm_2d%desc, dm_ham%cntxt)

  ! FHJ: Create an MPI Group and Communicator for processors inside the Blacs grid
  call logit('Determining active processors', params%kpp%inode==0, params%kpp%iunit)
  should_include(:) = .false.
  should_include(params%kpp%inode+1) = dm_2d%cntxt>=0
  call MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, should_include, &
    1, MPI_LOGICAL, params%kpp%comm, mpierr)

  group_size = 0
  do ipe = 0, params%kpp%npes-1
    if (should_include(ipe+1)) then
      group_size = group_size + 1
      ranks(group_size) = ipe
    endif
  enddo
  if (dm_2d%cntxt>=0 .and. (group_size /= dm_2d%nprow*dm_2d%npcol)) then
    write(0,'(/a)') 'ERROR: Inconsistent processors in group'
    write(0,'(2(a,i0)/)') 'group_size=', group_size, &
      ' nprow*npcol=', dm_2d%nprow*dm_2d%npcol
    call die('Inconsistent processors in group')
  endif

  call logit('Creating communicators', params%kpp%inode==0, params%kpp%iunit)
  call MPI_Comm_group(params%kpp%comm, orig_group, mpierr)
  call MPI_Group_incl(orig_group, group_size, ranks, active_group, mpierr)
  call MPI_Comm_create(params%kpp%comm, active_group, active_comm, mpierr)
  if ((dm_2d%cntxt>=0) .neqv. (active_comm/=MPI_COMM_NULL)) then
    write(0,'(/a)') 'ERROR: Inconsistent group created:'
    write(0,'(4(a,i0)/)') &
      'kp%inode=', params%kpp%inode, ' group_size=', group_size, &
      ' dm_2d%cntxt=', dm_2d%cntxt, ' active_comm=', active_comm
    call die('Inconsistent group')
  endif

  if (dm_2d%cntxt>=0) then
    ! FHJ: Call ELPA
    if (params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a,i0)') 'Beginning ELPA diagonalization. Size: ', dm_2d%N
      FLUSH(params%kpp%iunit)
    endif

    call logit('Calling elpa_allocate')
    elpa => elpa_allocate()
    call logit('ok')

    call elpa%set("na", dm_2d%N, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("nev", dm_wfn%N, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("local_nrows", dm_2d%Ml, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("local_ncols", dm_2d%Nl, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("nblk", dm_2d%Nb, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("mpi_comm_parent", active_comm, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("process_row", dm_2d%myprow, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("process_col", dm_2d%mypcol, elpa_err)
    call assert_elpa_ok(elpa_err)

    call assert_elpa_ok(elpa%setup())
    call elpa%set("solver", ELPA_SOLVER_2STAGE, elpa_err)
    call assert_elpa_ok(elpa_err)

    call logit('Calling elpa solver')
    call elpa%eigenvectors(ham_2d, en_t, wfn_2d, elpa_err)
    call assert_elpa_ok(elpa_err)
    call logit('ok')
    call MPI_Barrier(active_comm, mpierr)

    call logit('Calling elpa_deallocate')
    call elpa_deallocate(elpa)
    call logit('ok')
    call MPI_Barrier(active_comm, mpierr)

    ! FHJ: Copy eigenvalues/eigenvectors
    en(1:dm_wfn%N) = en_t(1:dm_wfn%N)

    call MPI_Comm_free(active_comm, mpierr)
    call MPI_Group_free(active_group, mpierr)
  endif !dm_2d%cntxt>=0

  call MPI_Barrier(params%kpp%comm, mpierr)
  ! FHJ: PEs outside the grid don`t have the evals!
  call MPI_Bcast(en, size(en), MPI_DOUBLE_PRECISION, 0, params%kpp%comm, mpierr)
  ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
  call logit('Redistributing wavefunctions', params%kpp%inode==0, params%kpp%iunit)
  call pcgemr2d(dm_wfn%M, dm_wfn%N, wfn_2d, 1, 1, dm_2d%desc, wfn_d, 1, 1, dm_wfn%desc, dm_wfn%cntxt)

  ! FHJ: Cleanup
  SAFE_DEALLOCATE(ham_2d)
  SAFE_DEALLOCATE(wfn_2d)
  call MPI_Group_free(orig_group, mpierr)
  call elpa_uninit()
  call MPI_Barrier(params%kpp%comm, mpierr)

  POP_SUB(diag_elpa_para_c)

end subroutine diag_elpa_para_c

subroutine diag_elpa_para_z(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  complex(DPC), intent(inout) :: ham_d(dm_ham%Ml, dm_ham%Nl)
  class(distrib_mat_t), intent(in) :: dm_wfn
  complex(DPC), intent(out) :: wfn_d(dm_wfn%Ml, dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  type(distrib_mat_t) :: dm_2d
  complex(DPC), allocatable :: ham_2d(:,:), wfn_2d(:,:)
  real(DP) :: en_t(dm_ham%M)
  integer :: comm_rows, comm_cols, ranks(dm_ham%nprow * dm_ham%npcol)
  integer :: ipe, group_size, orig_group, active_group, active_comm
  logical :: should_include(params%kpp%npes)
  class(elpa_t), pointer :: elpa
  integer :: elpa_err

  PUSH_SUB(diag_elpa_para_z)

  if (params%kpp%inode==0) then
    write(params%kpp%iunit,'(1x,a)') 'Solving with ELPA'
    FLUSH(params%kpp%iunit)
  endif

  !if (.not.elpa_initialized()) then
  !  if (elpa_init(20170403)/=elpa_ok) call die("ELPA API version not supported")
  !endif
  if (elpa_init(20170403)/=elpa_ok) call die("ELPA API version not supported")

  wfn_d(:,:) = ZERO
  en(:) = 0d0

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices ham_2d, wfn_2d.
  ! Some processors might get excluded for the sake of load balancing.
  call dm_2d%setup_mat_2d(dm_ham%M, dm_ham%N, params%block_sz, params%kpp)

  ! FHJ: Note that the distributed matrix wfn_2d is dm_ham%M x dm%ham%N = Hamiltonian size.
  ! ELPA requires this. However, we only copy the first dm_wfn%N columns from wfn_2d to
  ! wfn_d at the end.
  SAFE_ALLOCATE(ham_2d, (dm_2d%Ml, dm_2d%Nl))
  ham_2d = ZERO
  SAFE_ALLOCATE(wfn_2d, (dm_2d%Ml, dm_2d%Nl))
  wfn_2d = ZERO
  ! FHJ: Copy matrix ham_d from 1d cyclic layout into ham_2d (2d cyclic layout)
  call pzgemr2d(dm_ham%M, dm_ham%N, ham_d, 1, 1, dm_ham%desc, &
    ham_2d, 1, 1, dm_2d%desc, dm_ham%cntxt)

  ! FHJ: Create an MPI Group and Communicator for processors inside the Blacs grid
  call logit('Determining active processors', params%kpp%inode==0, params%kpp%iunit)
  should_include(:) = .false.
  should_include(params%kpp%inode+1) = dm_2d%cntxt>=0
  call MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, should_include, &
    1, MPI_LOGICAL, params%kpp%comm, mpierr)

  group_size = 0
  do ipe = 0, params%kpp%npes-1
    if (should_include(ipe+1)) then
      group_size = group_size + 1
      ranks(group_size) = ipe
    endif
  enddo
  if (dm_2d%cntxt>=0 .and. (group_size /= dm_2d%nprow*dm_2d%npcol)) then
    write(0,'(/a)') 'ERROR: Inconsistent processors in group'
    write(0,'(2(a,i0)/)') 'group_size=', group_size, &
      ' nprow*npcol=', dm_2d%nprow*dm_2d%npcol
    call die('Inconsistent processors in group')
  endif

  call logit('Creating communicators', params%kpp%inode==0, params%kpp%iunit)
  call MPI_Comm_group(params%kpp%comm, orig_group, mpierr)
  call MPI_Group_incl(orig_group, group_size, ranks, active_group, mpierr)
  call MPI_Comm_create(params%kpp%comm, active_group, active_comm, mpierr)
  if ((dm_2d%cntxt>=0) .neqv. (active_comm/=MPI_COMM_NULL)) then
    write(0,'(/a)') 'ERROR: Inconsistent group created:'
    write(0,'(4(a,i0)/)') &
      'kp%inode=', params%kpp%inode, ' group_size=', group_size, &
      ' dm_2d%cntxt=', dm_2d%cntxt, ' active_comm=', active_comm
    call die('Inconsistent group')
  endif

  if (dm_2d%cntxt>=0) then
    ! FHJ: Call ELPA
    if (params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a,i0)') 'Beginning ELPA diagonalization. Size: ', dm_2d%N
      FLUSH(params%kpp%iunit)
    endif

    call logit('Calling elpa_allocate')
    elpa => elpa_allocate()
    call logit('ok')

    call elpa%set("na", dm_2d%N, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("nev", dm_wfn%N, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("local_nrows", dm_2d%Ml, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("local_ncols", dm_2d%Nl, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("nblk", dm_2d%Nb, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("mpi_comm_parent", active_comm, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("process_row", dm_2d%myprow, elpa_err)
    call assert_elpa_ok(elpa_err)
    call elpa%set("process_col", dm_2d%mypcol, elpa_err)
    call assert_elpa_ok(elpa_err)

    call assert_elpa_ok(elpa%setup())
    call elpa%set("solver", ELPA_SOLVER_2STAGE, elpa_err)
    call assert_elpa_ok(elpa_err)

    call logit('Calling elpa solver')
    call elpa%eigenvectors(ham_2d, en_t, wfn_2d, elpa_err)
    call assert_elpa_ok(elpa_err)
    call logit('ok')
    call MPI_Barrier(active_comm, mpierr)

    call logit('Calling elpa_deallocate')
    call elpa_deallocate(elpa)
    call logit('ok')
    call MPI_Barrier(active_comm, mpierr)

    ! FHJ: Copy eigenvalues/eigenvectors
    en(1:dm_wfn%N) = en_t(1:dm_wfn%N)

    call MPI_Comm_free(active_comm, mpierr)
    call MPI_Group_free(active_group, mpierr)
  endif !dm_2d%cntxt>=0

  call MPI_Barrier(params%kpp%comm, mpierr)
  ! FHJ: PEs outside the grid don`t have the evals!
  call MPI_Bcast(en, size(en), MPI_DOUBLE_PRECISION, 0, params%kpp%comm, mpierr)
  ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
  call logit('Redistributing wavefunctions', params%kpp%inode==0, params%kpp%iunit)
  call pzgemr2d(dm_wfn%M, dm_wfn%N, wfn_2d, 1, 1, dm_2d%desc, wfn_d, 1, 1, dm_wfn%desc, dm_wfn%cntxt)

  ! FHJ: Cleanup
  SAFE_DEALLOCATE(ham_2d)
  SAFE_DEALLOCATE(wfn_2d)
  call MPI_Group_free(orig_group, mpierr)
  call elpa_uninit()
  call MPI_Barrier(params%kpp%comm, mpierr)

  POP_SUB(diag_elpa_para_z)

end subroutine diag_elpa_para_z
#endif

end module diag_elpa_m
