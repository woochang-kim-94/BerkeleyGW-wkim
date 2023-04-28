#include "f_defs.h"
#if defined MPI && !defined USESCALAPACK
  #error ScaLAPACK is required for MPI builds.
#endif
#if defined USEELPA && !defined MPI
  #error MPI is required for ELPA builds.
#endif

module diag_driver_m

  use global_m
  use inread_m,       only: pb_params_t
  use distribution_m, only: distrib_mat_t
  use diag_scalapack_m
#ifdef USEELPA
  use diag_elpa_m
#endif
#ifdef USEPRIMME
  use diag_primme_m
#endif

  implicit none

  private

  public :: check_valid_solver, diag_driver

contains


subroutine check_valid_solver(params, solver_name)
  type(pb_params_t), intent(in) :: params
  character(len=*), intent(out) :: solver_name

  character(len=4) :: pxhe
  character(len=3) :: xhe

  PUSH_SUB(check_valid_solver)

  solver_name = 'UNKNOWN'
#ifdef CPLX
  pxhe = 'PZHE'
  xhe = 'ZHE'
#else
  pxhe = 'PDSY'
  xhe = 'DSY'
#endif

#ifdef MPI
  if (params%kpp%npes>1) then
    select case (params%solver_alg)
      case (-2)
        solver_name = 'Dummy'
      case (-1)
#ifdef USEELPA
        solver_name = 'ELPA'
#else
        solver_name = pxhe//'EVX'
#endif
      case (0)
        solver_name = pxhe//'EVX'
      case (1)
        solver_name = pxhe//'EVD'
      case (2)
#ifdef USEMR3
        solver_name = pxhe//'EVR'
#else
        call die('Not compiled with support for MR3', only_root_writes=.true.)
#endif
      case (10)
#ifdef USEELPA
        solver_name = 'ELPA'
#else
        call die('Not compiled with support for ELPA', only_root_writes=.true.)
#endif
      case (100:115)
#ifdef USEPRIMME
        solver_name = 'PRIMME'
#else
        call die('Not compiled with support for PRIMME', only_root_writes=.true.)
#endif
      case default
        if (params%kpp%inode==0) &
          write(0,'(a,i0)') 'ERROR: Got invalid solver_algorithm: ', params%solver_alg
        call die('Invalid option for solver_algorithm', only_root_writes=.true.)
    endselect
  else
#endif
    select case (params%solver_alg)
      case (-2)
        solver_name = 'Dummy'
      case (-1)
        solver_name = xhe//'EVX'
      case (0)
        solver_name = xhe//'EVX'
      case (1)
        solver_name = xhe//'EVD'
      case (2)
        solver_name = xhe//'EVR'
      case (10)
#ifdef USEELPA
        solver_name = 'ELPA'
#else
        call die('Not compiled with support for ELPA', only_root_writes=.true.)
#endif
      case (100:115)
#ifdef USEPRIMME
        solver_name = 'PRIMME'
#else
        call die('Not compiled with support for PRIMME', only_root_writes=.true.)
#endif
      case default
        if (params%kpp%inode==0) &
          write(0,'(a,i0)') 'ERROR: Got invalid solver_algorithm: ', params%solver_alg
        call die('Invalid option for solver_algorithm', only_root_writes=.true.)
    endselect
#ifdef MPI
  endif
#endif

  POP_SUB(check_valid_solver)

end subroutine check_valid_solver


subroutine diag_driver(dm_ham, ham_d, dm_wfn, wfn_d, en, params)
  class(distrib_mat_t), intent(in) :: dm_ham
  !SCALAR, intent(inout) :: ham_d(dm_ham%Ml, dm_ham%Nl)
  SCALAR, intent(inout), allocatable :: ham_d(:,:)
  class(distrib_mat_t), intent(in) :: dm_wfn
  SCALAR, intent(out) :: wfn_d(dm_wfn%Ml, dm_wfn%Nl)
  real(DP), intent(out) :: en(dm_wfn%N)
  type(pb_params_t), intent(in) :: params

  integer :: ib
#ifdef USEELPA
  integer :: nprow_max
#endif
  logical, save :: warned_elpa_block_size=.false.

  PUSH_SUB(diag_driver)

#ifdef MPI
  if (params%kpp%npes>1) then
    select case (params%solver_alg)
      case (-2)
        ! DOS ~ E^(1/2) => E_i = Emax * (i/N)^(2/3)
        en(:) = [ (60d0 * (dble(ib)/dm_wfn%N)**(2d0/3d0), ib=1,dm_wfn%N) ]
        wfn_d(:,:) = ONE
      case (-1)
#ifdef USEELPA
        nprow_max = int(ceiling(sqrt(dble(params%kpp%npes))))
        if (params%block_sz*(nprow_max - 1) >= dm_ham%N) then
          if (params%kpp%inode==0 .and. .not.warned_elpa_block_size) then
            write(0,'(/1x,a,i0/)') &
              'WARNING: switching from ELPA to ScaLAPACK solver due to large block size'
            warned_elpa_block_size = .true.
          endif
          call diag_scalapack_para_x(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
        else
          call diag_elpa_para(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
        endif
#else
        call diag_scalapack_para_x(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
#endif
      case (0)
        call diag_scalapack_para_x(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
      case (1)
        call diag_scalapack_para_d(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
      case (2)
#ifdef USEMR3
        call diag_scalapack_para_r(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
#else
        call die('Not compiled with support for MR3', only_root_writes=.true.)
#endif
      case (10)
#ifdef USEELPA
        nprow_max = int(ceiling(sqrt(dble(params%kpp%npes))))
        if (params%block_sz*(nprow_max - 1) >= dm_ham%N) then
          if (params%kpp%inode==0 .and. .not.warned_elpa_block_size) then
            write(0,'(/1x,a,i0/)') &
              'WARNING: switching from ELPA to ScaLAPACK solver due to large block size'
            warned_elpa_block_size = .true.
          endif
          call diag_scalapack_para_x(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
        else
          call diag_elpa_para(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
        endif
#else
        call die('Not compiled with support for ELPA', only_root_writes=.true.)
#endif
      case (100:115)
#ifdef USEPRIMME
        call diag_primme_para(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
#else
        call die('Not compiled with support for PRIMME', only_root_writes=.true.)
#endif
      case default
        if (params%kpp%inode==0) &
          write(0,'(a,i0)') 'ERROR: Got solver_alg: ', params%solver_alg
        call die('Invalid option for the solver_alg', only_root_writes=.true.)
    endselect
  else
#endif
    select case (params%solver_alg)
      case (-2)
        ! DOS ~ E^(1/2) => E_i = Emax * (i/N)^(2/3)
        en(:) = [ (60d0 * (dble(ib)/dm_wfn%N)**(2d0/3d0), ib=1,dm_wfn%N) ]
        wfn_d(:,:) = ONE
      case (-1)
        call diag_scalapack_serial_x(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
      case (0)
        call diag_scalapack_serial_x(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
      case (1)
        call diag_scalapack_serial_d(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
      case (2)
        call diag_scalapack_serial_r(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
      case (10)
#ifdef USEELPA
        call diag_elpa_para(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
#else
        call die('Not compiled with support for ELPA', only_root_writes=.true.)
#endif
      case (100:115)
#ifdef USEPRIMME
        call diag_primme_serial(dm_ham, ham_d, dm_wfn, wfn_d(:,:), en(:), params)
#else
        call die('Not compiled with support for PRIMME', only_root_writes=.true.)
#endif
      case default
        if (params%kpp%inode==0) &
          write(0,'(a,i0)') 'ERROR: Got solver_alg: ', params%solver_alg
        call die('Invalid option for the solver_alg', only_root_writes=.true.)
    endselect
#ifdef MPI
  endif
#endif

  POP_SUB(diag_driver)

end subroutine diag_driver


end module diag_driver_m
