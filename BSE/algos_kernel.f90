#include "f_defs.h"

! Module to control selection of various accelerated algos in Kernel
! Any build-specific preprocessor directives should be as high-level
! as possible, i.e. no explicit framework-based calls
module algos_kernel_m

  use ISO_C_BINDING

  use accel_memory_m
  use algos_common_m
  use message_m
  use nrtype_m, only: DP, DPC
  use peinfo_m, only: peinf
  use push_pop_m

  implicit none

  private

  ! Publicly exposed subroutines
  public :: CPU_ALGO, OPENACC_ALGO, OMP_TARGET_ALGO
  public :: output_algos
  public :: set_algos_to_cpu
  public :: set_algos_to_best_available_gpu
  public :: verify_gpu_settings
  public :: initialize_gpu
  public :: algos_inread
  public :: allocate_accel_kernel, deallocate_accel_kernel

  ! General control variables
  integer(kind(CPU_ALGO)), public :: mtxel_algo
  integer(kind(CPU_ALGO)), public :: w_sum_algo
  integer(kind(CPU_ALGO)), public :: g_sum_algo

  ! Variables for OpenACC acceleration
  ! WPH:  We put the OpenACC variables here, rather than in their own separate
  !       OpenACC module, because OpenACC is valid CPU code.
  SCALAR, target, allocatable, public :: tempb_accel(:,:,:,:), temph_accel(:,:,:)
  SCALAR, target, allocatable, public :: tempw_accel(:,:,:,:)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     __  ____  _ ___ __  _           !
  !    / / / / /_(_) (_) /_(_)__  _____ !
  !   / / / / __/ / / / __/ / _ \/ ___/ !
  !  / /_/ / /_/ / / / /_/ /  __(__  )  !
  !  \____/\__/_/_/_/\__/_/\___/____/   !
  !                                     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_algos()
    PUSH_SUB(output_algos)

    if (peinf%inode.eq.0) then
      call output_algos_common()
      write(6,*)
      write(6,'(1x,a)')   'Algorithms used:'
      write(6,'(1x,a,a)')  '- mtxel : ', trim(get_algo_str(mtxel_algo))
      write(6,'(1x,a,a)')  '- w_sum : ', trim(get_algo_str(w_sum_algo))
      write(6,'(1x,a,a)')  '- g_sum : ', trim(get_algo_str(g_sum_algo))
    endif

    POP_SUB(output_algos)
  end subroutine

  subroutine set_algos_to_cpu()
    PUSH_SUB(set_algos_to_cpu)

    call set_algos_to_cpu_common()
    mtxel_algo = CPU_ALGO
    g_sum_algo = CPU_ALGO
    w_sum_algo = CPU_ALGO

    POP_SUB(set_algos_to_cpu)
  end subroutine

  subroutine set_algos_to_best_available_gpu()
    PUSH_SUB(set_algos_to_best_available_gpu)

    call set_algos_to_best_available_gpu_common()
#ifdef OMP_TARGET
    mtxel_algo = OMP_TARGET_ALGO
    g_sum_algo = OMP_TARGET_ALGO
    w_sum_algo = OMP_TARGET_ALGO
#endif
#ifdef OPENACC
    mtxel_algo = OPENACC_ALGO
    g_sum_algo = OPENACC_ALGO
    w_sum_algo = OPENACC_ALGO
#endif

    POP_SUB(set_algos_to_best_available_gpu)
  end subroutine

  subroutine verify_gpu_settings()
    PUSH_SUB(verify_gpu_settings)

    call verify_gpu_settings_common()

    ! First check if the desired algorithm is compiled in
#ifndef OMP_TARGET
    if (mtxel_algo == OMP_TARGET_ALGO .or. &
        w_sum_algo == OMP_TARGET_ALGO .or. &
        g_sum_algo == OMP_TARGET_ALGO) then
      call die_algos("OpenMP Target")
    end if
#endif
#ifndef OPENACC
    if (mtxel_algo == OPENACC_ALGO .or. &
        w_sum_algo == OPENACC_ALGO .or. &
        g_sum_algo == OPENACC_ALGO) then
      call die_algos("OpenACC")
    end if
#endif

    ! Then check if we`re using the GPU
    if (.not. use_gpu .and. (mtxel_algo /= CPU_ALGO .or. &
                             w_sum_algo /= CPU_ALGO .or. &
                             g_sum_algo /= CPU_ALGO)) then
      call die("You have specified one or more GPU-accelerated algorithms, but&
               & the GPU is not enabled.", only_root_writes=.true.)
    end if

    ! Finally, check for algorithms that don`t play nicely with one another
    if (w_sum_algo /= g_sum_algo) then
      call die("You have specified different algorithms for w_sum() and g_sum().&
               & They must use the same algorithm.", only_root_writes=.true.)
    end if

    POP_SUB(verify_gpu_settings)
  end subroutine

  subroutine initialize_gpu(algo_type)
    integer(kind(CPU_ALGO)), intent(in) :: algo_type
    PUSH_SUB(initialize_gpu)
    call initialize_gpu_common(algo_type)
    POP_SUB(initialize_gpu)
  end subroutine initialize_gpu

  subroutine algos_inread(keyword, line, found)
    character(len=*), intent(in) :: keyword, line
    logical, intent(out) :: found

    character*256 :: errmsg

    PUSH_SUB(algos_inread)

    call algos_inread_common(keyword, line, found)

    if (.not. found) then
      found = .true.
      if (trim(keyword) .eq. 'use_gpu') then
        if (trim(line) .eq. '.true.') then
#if defined(OPENACC) || defined(OMP_TARGET)
          call set_algos_to_best_available_gpu()
#else
          write (errmsg,'(a)') 'GPU acceleration not compiled into this executable, use_gpu keyword is invalid.'
          call die(errmsg, only_root_writes = .true.)
#endif
        else if(trim(line) .eq. '.false.') then
          call set_algos_to_cpu()
        else
          write(errmsg,'(3a)') 'Unexpected parameter "', trim(line), '" for keyword "use_gpu" was found in kernel.inp.'
          call die(errmsg, only_root_writes = .true.)
        end if
      else if(trim(keyword) .eq. 'mtxel_algo') then
        mtxel_algo = get_algo(trim(line))
      else if(trim(keyword) .eq. 'g_sum_algo') then
        g_sum_algo = get_algo(trim(line))
      else if(trim(keyword) .eq. 'w_sum_algo') then
        w_sum_algo = get_algo(trim(line))
      else
        found = .false.
      end if
    end if

    POP_SUB(algos_inread)
  end subroutine algos_inread

  subroutine allocate_accel_kernel(bsedbody, bsedhead, bsedwing, n_ele, n_spin)

    SCALAR, allocatable, intent(out) :: bsedbody(:, :, :)
    SCALAR, allocatable, intent(out) :: bsedhead(:, :, :)
    SCALAR, allocatable, intent(out) :: bsedwing(:, :, :)

    integer, intent(in) :: n_ele
    integer, intent(in) :: n_spin

    PUSH_SUB(allocate_accel_kernel)

    allocate(bsedbody(n_ele, n_spin, n_spin))
    allocate(bsedhead(n_ele, n_spin, n_spin))
    allocate(bsedwing(n_ele, n_spin, n_spin))
    bsedbody(:,:,:) = ZERO
    bsedhead(:,:,:) = ZERO
    bsedwing(:,:,:) = ZERO

    call accel_enter_data_map_to(bsedbody, g_sum_algo)
    call accel_enter_data_map_to(bsedhead, g_sum_algo)
    call accel_enter_data_map_to(bsedwing, g_sum_algo)

    POP_SUB(allocate_accel_kernel)
  end subroutine allocate_accel_kernel

  subroutine deallocate_accel_kernel(bsedbody, bsedhead, bsedwing)

    SCALAR, allocatable, intent(out) :: bsedbody(:, :, :)
    SCALAR, allocatable, intent(out) :: bsedhead(:, :, :)
    SCALAR, allocatable, intent(out) :: bsedwing(:, :, :)

    PUSH_SUB(deallocate_accel_kernel)

    call accel_exit_data_map_delete(bsedbody, g_sum_algo)
    call accel_exit_data_map_delete(bsedhead, g_sum_algo)
    call accel_exit_data_map_delete(bsedwing, g_sum_algo)

    SAFE_DEALLOCATE(bsedhead)
    SAFE_DEALLOCATE(bsedbody)
    SAFE_DEALLOCATE(bsedwing)

    POP_SUB(deallocate_accel_kernel)
  end subroutine deallocate_accel_kernel

end module algos_kernel_m
