#include "f_defs.h"

! WPH (12 July 2021):  This module is currently a nightmare
! from a code smell standpoint; everything is a copy-paste of a
! small number of subroutine templates which, in principle,
! should be united.
! The problem is the way that OpenACC handles data clauses in
! Fortran.  OpenACC requires that it knows the bounds of arrays
! via interfaces in order to offload them via data regions in the
! Fortran style, so using assumed size arrays doesn't work (I
! tried.)
! One way to get around this would be to adopt a C-style
! approach where we always explictly pass the bounds of arrays
! and are explicit in data clauses which data ranges are being
! offloaded/updated, in which case we should be able to use
! base subroutines with assumed size arrays wrapped by wrapper
! subroutines with assumed shape arrays to play nicely with
! module procedure interfaces.
! The C-style approach would require a lot of modification in
! the main code, as we would need to modify every data movement
! call, and there would still be some copy-paste to support
! int, real, and complex data types (until we can use assumed
! type arguments from Fortran 2018, that is).

module accel_memory_m

  use algos_common_m
  use typedefs_m

#if defined (OPENACC) || defined (OMP_TARGET)
  use openacc
#endif

  implicit none

  private

  public :: accel_enter_data_map_alloc
  public :: accel_enter_data_map_to
  public :: accel_exit_data_map_from
  public :: accel_exit_data_map_delete
  public :: accel_update_to
  public :: accel_update_from

  interface accel_enter_data_map_alloc
    module procedure accel_enter_data_map_alloc_int_1d
    module procedure accel_enter_data_map_alloc_real_1d
    module procedure accel_enter_data_map_alloc_cplx_1d
    module procedure accel_enter_data_map_alloc_int_2d
    module procedure accel_enter_data_map_alloc_real_2d
    module procedure accel_enter_data_map_alloc_cplx_2d
    module procedure accel_enter_data_map_alloc_cplx_3d
  end interface

  interface accel_enter_data_map_to
    module procedure accel_enter_data_map_to_int_1d
    module procedure accel_enter_data_map_to_int_2d
    module procedure accel_enter_data_map_to_real_1d
    module procedure accel_enter_data_map_to_cplx_1d
    module procedure accel_enter_data_map_to_real_2d
    module procedure accel_enter_data_map_to_cplx_2d
    module procedure accel_enter_data_map_to_real_3d
    module procedure accel_enter_data_map_to_cplx_3d
    module procedure accel_enter_data_map_to_real_4d
    module procedure accel_enter_data_map_to_cplx_4d
    module procedure accel_enter_data_map_to_cplx_5d
  end interface

  interface accel_exit_data_map_from
    module procedure accel_exit_data_map_from_real_2d
    module procedure accel_exit_data_map_from_cplx_2d
    module procedure accel_exit_data_map_from_cplx_3d
    module procedure accel_exit_data_map_from_cplx_5d
  end interface

  interface accel_exit_data_map_delete
    module procedure accel_exit_data_map_delete_int_1d
    module procedure accel_exit_data_map_delete_real_1d
    module procedure accel_exit_data_map_delete_cplx_1d
    module procedure accel_exit_data_map_delete_int_2d
    module procedure accel_exit_data_map_delete_real_2d
    module procedure accel_exit_data_map_delete_cplx_2d
    module procedure accel_exit_data_map_delete_real_3d
    module procedure accel_exit_data_map_delete_cplx_3d
    module procedure accel_exit_data_map_delete_real_4d
    module procedure accel_exit_data_map_delete_cplx_4d
    module procedure accel_exit_data_map_delete_cplx_5d
  end interface

  interface accel_update_to
    module procedure accel_update_to_int_1d
    module procedure accel_update_to_int_2d
    module procedure accel_update_to_real_1d
    module procedure accel_update_to_cplx_1d
    module procedure accel_update_to_real_1d_subset
    module procedure accel_update_to_cplx_1d_subset
    module procedure accel_update_to_real_2d
    module procedure accel_update_to_cplx_2d
    module procedure accel_update_to_real_3d
    module procedure accel_update_to_cplx_3d
    module procedure accel_update_to_real_4d
    module procedure accel_update_to_cplx_4d
  end interface

  interface accel_update_from
    module procedure accel_update_from_real_1d
    module procedure accel_update_from_cplx_1d
    module procedure accel_update_from_real_1d_subset
    module procedure accel_update_from_cplx_1d_subset
    module procedure accel_update_from_real_2d
    module procedure accel_update_from_cplx_2d
    module procedure accel_update_from_real_3d
    module procedure accel_update_from_cplx_3d
    module procedure accel_update_from_real_4d
    module procedure accel_update_from_cplx_4d
  end interface

contains

  subroutine accel_enter_data_map_alloc_int_1d(array, algo)
    implicit none

    integer, intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data create(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(alloc: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_alloc", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_alloc_real_1d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data create(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(alloc: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_alloc", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_alloc_cplx_1d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data create(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(alloc: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_alloc", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_alloc_int_2d(array, algo)
    implicit none

    integer, intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data create(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(alloc: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_alloc", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_alloc_real_2d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data create(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(alloc: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_alloc", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_alloc_cplx_2d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data create(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(alloc: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_alloc", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_alloc_cplx_3d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data create(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(alloc: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_alloc", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_int_1d(array, algo)
    implicit none

    integer, intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_int_2d(array, algo)
    implicit none

    integer, intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_real_1d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_cplx_1d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_real_2d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_cplx_2d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_real_3d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_cplx_3d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_real_4d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
     !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_cplx_4d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_enter_data_map_to_cplx_5d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc enter data copyin(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target enter data map(to: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_enter_data_map_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_from_real_2d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data copyout(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(from: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_from_cplx_2d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data copyout(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(from: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_from_cplx_3d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data copyout(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(from: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_from_cplx_5d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data copyout(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(from: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_int_1d(array, algo)
    implicit none

    integer, intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_int_2d(array, algo)
    implicit none

    integer, intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_real_1d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_cplx_1d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_real_2d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_cplx_2d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_real_3d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_cplx_3d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_real_4d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_cplx_4d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_exit_data_map_delete_cplx_5d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc exit data delete(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target exit data map(delete: array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_exit_data_map_delete", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_int_1d(array, algo)
    implicit none

    integer, intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update device(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_int_2d(array, algo)
    implicit none

    integer, intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update device(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_real_1d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update device(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_cplx_1d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update device(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_real_1d_subset(array, start, length, algo)
    implicit none

    real(DP), intent(inout) :: array(:)
    integer, intent(in) :: start
    integer, intent(in) :: length
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array(start:length))
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_cplx_1d_subset(array, start, length, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:)
    integer, intent(in) :: start
    integer, intent(in) :: length
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array(start:length))
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_from_real_1d_subset(array, start, length, algo)
    implicit none

    real(DP), intent(inout) :: array(:)
    integer, intent(in) :: start
    integer, intent(in) :: length
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update from(array(start:length))
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_from_cplx_1d_subset(array, start, length, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:)
    integer, intent(in) :: start
    integer, intent(in) :: length
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update from(array(start:length))
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_real_2d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update device(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_cplx_2d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update device(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_real_3d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update device(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_cplx_3d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update device(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_real_4d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update device(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_to_cplx_4d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update device(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update to(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_to", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_from_real_1d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update self(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update from(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_from_cplx_1d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update self(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update from(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_from_real_2d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update self(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update from(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_from_cplx_2d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update self(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update from(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_from_real_3d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update self(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update from(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_from_cplx_3d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update self(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update from(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_from_real_4d(array, algo)
    implicit none

    real(DP), intent(inout) :: array(:, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update self(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update from(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_from", &
               only_root_writes = .true.)
    end select
  end subroutine

  subroutine accel_update_from_cplx_4d(array, algo)
    implicit none

    complex(DPC), intent(inout) :: array(:, :, :, :)
    integer(kind(CPU_ALGO)), intent(in) :: algo

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      !$acc update self(array)
#else
      call die_algos("OpenACC")
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      !$omp target update from(array)
#else
      call die_algos("OpenMP Target")
#endif
    case (CPU_ALGO)
    case default
      call die("Invald algorithm for accel_update_from", &
               only_root_writes = .true.)
    end select
  end subroutine

end module accel_memory_m
