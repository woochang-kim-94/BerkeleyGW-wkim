#include "f_defs.h"

! Module to control selection of various accelerated algos in Kernel
! Any build-specific preprocessor directives should be as high-level
! as possible, i.e. no explicit framework-based calls
module algos_sigma_m

  use ISO_C_BINDING

  use algos_common_m
  use message_m
  use peinfo_m, only: peinf
  use push_pop_m

  implicit none

  private

  ! Publicly exposed members
  public :: CPU_ALGO, OPENACC_ALGO, OMP_TARGET_ALGO
  public :: output_algos
  public :: set_algos_to_cpu
  public :: set_algos_to_best_available_gpu
  public :: verify_gpu_settings
  public :: initialize_gpu
  public :: algos_inread
  public :: use_gpu

  ! General control variables
  integer(kind(CPU_ALGO)), public :: mtxel_algo
  integer(kind(CPU_ALGO)), public :: sigma_gpp_algo

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
      write(6,'(1x,a,a)')  '- mtxel     : ', trim(get_algo_str(mtxel_algo))
      write(6,'(1x,a,a)')  '- sigma_gpp : ', trim(get_algo_str(sigma_gpp_algo))
    endif
    POP_SUB(output_algos)
  end subroutine

  subroutine set_algos_to_cpu()
    PUSH_SUB(set_algos_to_cpu)
    call set_algos_to_cpu_common()
    mtxel_algo = CPU_ALGO
    sigma_gpp_algo = CPU_ALGO
    POP_SUB(set_algos_to_cpu)
  end subroutine

  subroutine set_algos_to_best_available_gpu()
    PUSH_SUB(set_algos_to_best_available_gpu)
    call set_algos_to_best_available_gpu_common()
#ifdef OMP_TARGET
    mtxel_algo = OMP_TARGET_ALGO
    sigma_gpp_algo = OMP_TARGET_ALGO
#endif
#ifdef OPENACC
    mtxel_algo = OPENACC_ALGO
    sigma_gpp_algo = OPENACC_ALGO
#endif
    POP_SUB(set_algos_to_best_available_gpu)
  end subroutine

  subroutine verify_gpu_settings()
    PUSH_SUB(verify_gpu_settings)
    ! the main body of the code.
    call verify_gpu_settings_common()

    ! First check if the desired algorithm is compiled in
#ifndef OPENACC
    if (mtxel_algo == OPENACC_ALGO .or. &
        sigma_gpp_algo == OPENACC_ALGO) then
      call die_algos("OpenACC")
    end if
#endif

#ifndef OMP_TARGET
    if (mtxel_algo == OMP_TARGET_ALGO .or. &
        sigma_gpp_algo == OMP_TARGET_ALGO) then
      call die_algos("OpenMP Target")
    end if
#endif

    ! Then check if we're using the GPU
    if (.not. use_gpu .and. (mtxel_algo /= CPU_ALGO .or. &
                             sigma_gpp_algo /= CPU_ALGO)) then
      call die("You have specified one or more GPU-accelerated algorithms, but&
               & the GPU is not enabled.", only_root_writes=.true.)
    end if

    ! Finally, check for algorithms that don't play nicely with one another
    ! None so far...
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
        elseif(trim(line) .eq. '.false.') then
          call set_algos_to_cpu()
        else
          write(errmsg, '(3a)') 'Unexpected parameter "', trim(line), '" for keyword "use_gpu" was found in kernel.inp.'
          call die(errmsg, only_root_writes = .true.)
        end if
      else if (trim(keyword) .eq. 'mtxel_algo') then
        mtxel_algo = get_algo(trim(line))
      else if (trim(keyword) .eq. 'sigma_gpp_algo') then
        sigma_gpp_algo = get_algo(trim(line))
      else
        found = .false.
      end if
    end if

    POP_SUB(algos_inread)
  end subroutine algos_inread

end module algos_sigma_m
