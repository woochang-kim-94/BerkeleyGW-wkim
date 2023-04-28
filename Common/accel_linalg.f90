#include "f_defs.h"

module accel_linalg_m

  use ISO_C_BINDING

  use algos_common_m
  use blas_m
  use message_m,  only: die
  use nrtype_m

#if defined (OPENACC) || defined (OMP_TARGET)
  use cublas
  use cudafor
  use openacc
#endif

  implicit none

  private

  logical :: common_blas_handle_exist = .false.
#if defined (OPENACC) || defined (OMP_TARGET)
  type(cublasHandle) :: cublas_handle_common
#endif

  public :: accel_dgemm
  public :: accel_zgemm
  public :: accel_xgemm

  public :: accel_dgemv
  public :: accel_zgemv
  public :: accel_xgemv

  public :: create_common_blas_handle, destroy_common_blas_handle
  public :: set_queue_to_common_blas_handle

  ! WPH: Implementation of X(gemm) as a generic subroutine
  interface accel_xgemm
    module procedure accel_dgemm_2d
    module procedure accel_dgemm_3d_3d_1d
    module procedure accel_dgemm_2d_5d_5d
    module procedure accel_zgemm_2d
    module procedure accel_zgemm_3d_3d_1d
    module procedure accel_zgemm_2d_5d_5d
  end interface

  interface accel_dgemm
    module procedure accel_dgemm_2d
    module procedure accel_dgemm_3d_3d_1d
    module procedure accel_dgemm_2d_5d_5d
  end interface

  interface accel_zgemm
    module procedure accel_zgemm_2d
    module procedure accel_zgemm_3d_3d_1d
    module procedure accel_zgemm_2d_5d_5d
  end interface

  ! WPH: Implementation of X(gemv) as a generic subroutine
  interface accel_xgemv
    module procedure accel_zgemv_3d_3d_1d
    module procedure accel_dgemv_3d_3d_1d
  end interface

  interface accel_dgemv
    module procedure accel_dgemv_3d_3d_1d
  end interface

  interface accel_zgemv
    module procedure accel_zgemv_3d_3d_1d
  end interface

contains

  subroutine create_common_blas_handle(algo)
    integer(kind(CPU_ALGO)) :: algo
    integer :: cuErr
   
    common_blas_handle_exist = .false.
    select case (algo)
    case (OPENACC_ALGO, OMP_TARGET_ALGO)
#if defined (OPENACC) || defined (OMP_TARGET)
      cuErr = cublasCreate(cublas_handle_common)
      common_blas_handle_exist = .true.
#else
      call die("OpenACC/OpenMP-target not "&
            &"compiled into this executable", only_root_writes = .true.)
#endif
    case default
      ! do nothing
    end select

  end subroutine

  subroutine destroy_common_blas_handle(algo)
    integer(kind(CPU_ALGO)) :: algo
    integer :: cuErr

    if ( common_blas_handle_exist ) then
      select case (algo)
      case (OPENACC_ALGO, OMP_TARGET_ALGO)
#if defined (OPENACC) || defined (OMP_TARGET)
        cuErr = cublasDestroy(cublas_handle_common)
#else
        call die("OpenACC/OpenMP-target not "&
              &"compiled into this executable", only_root_writes = .true.)
#endif
      case default
        ! do nothing
      end select
      ! set it back to false
      common_blas_handle_exist = .false.
    end if

  end subroutine

  subroutine set_queue_to_common_blas_handle(algo, queue)
    integer(kind(CPU_ALGO)) :: algo
    integer :: queue
    integer :: cuErr

    if ( common_blas_handle_exist ) then
      select case (algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        cuErr = cublasSetStream(cublas_handle_common, acc_get_cuda_stream(queue))
#else
        call die("Can not set stream/queue, OpenACC "&
              &"compiled into this executable", only_root_writes = .true.)
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        ! no queue for now with OpenMP target
#else
#endif
      case default
        ! do nothing
      end select
    end if

  end subroutine

  subroutine accel_dgemm_core(transa, transb, &
                              m, n, k, &
                              alpha, &
                              a, lda, &
                              b, ldb, &
                              beta, &
                              c, ldc, &
                              algo)
    implicit none

    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double precision :: alpha, beta
    double precision :: a(*), b(*), c(*)
    integer(kind(CPU_ALGO)) :: algo

#if defined (OPENACC) || defined (OMP_TARGET)
    integer :: cuErr
    type(cublasHandle) :: cublas_handle
#endif

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      cuErr = cublasCreate(cublas_handle)
      !$acc host_data use_device(a, b, c)
      call cublasDgemm(transa, transb, &
                       m, n, k, &
                       alpha, a, lda, b, ldb, &
                       beta, c, ldc)
      !$acc end host_data
      cuErr = cublasDestroy(cublas_handle)
#else
      call die("OpenACC version of accel_dgemm requested, but OpenACC not "&
            &"compiled into this executable", only_root_writes = .true.)
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      cuErr = cublasCreate(cublas_handle)
      !$omp target data use_device_ptr(a, b, c)
      call cublasDgemm(transa, transb, &
                       m, n, k, &
                       alpha, a, lda, b, ldb, &
                       beta, c, ldc)
      !$omp end target data
      cuErr = cublasDestroy(cublas_handle)
#else
      call die("OpenMP Target version of accel_dgemm requested, but OpenMP Target not "&
            &"compiled into this executable", only_root_writes = .true.)
#endif
    case (CPU_ALGO)
      call dgemm(transa, transb, &
                 m, n, k, &
                 alpha, a, lda, b, ldb, &
                 beta, c, ldc)
    case default
      call die("Invald algorithm for accel_dgemm", only_root_writes = .true.)
    end select

  end subroutine accel_dgemm_core

  subroutine accel_zgemm_core(transa, transb, &
                              m, n, k, &
                              alpha, &
                              a, lda, &
                              b, ldb, &
                              beta, &
                              c, ldc, &
                              algo)
    implicit none

    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double complex :: alpha, beta
    double complex :: a(*), b(*), c(*)
    integer(kind(CPU_ALGO)) :: algo

#if defined (OPENACC) || defined (OMP_TARGET)
    integer :: cuErr
    type(cublasHandle) :: cublas_handle
#endif

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      if ( common_blas_handle_exist ) then
        !$acc host_data use_device(a, b, c)
        cuErr = cublasZgemm_v2(cublas_handle_common, &
                               cublasConvertCharToV2(transa), cublasConvertCharToV2(transb), &
                               m, n, k, &
                               alpha, a, lda, b, ldb, &
                               beta, c, ldc)
        !$acc end host_data
      else
        cuErr = cublasCreate(cublas_handle)
        !$acc host_data use_device(a, b, c)
        call cublasZgemm(transa, transb, &
                         m, n, k, &
                         alpha, a, lda, b, ldb, &
                         beta, c, ldc)
        !$acc end host_data
        cuErr = cublasDestroy(cublas_handle)
      end if
#else
      call die("OpenACC version of accel_zgemm requested, but OpenACC not "&
            &"compiled into this executable", only_root_writes = .true.)
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      if ( common_blas_handle_exist ) then
        !$omp target data use_device_ptr(a, b, c)
        cuErr = cublasZgemm_v2(cublas_handle_common, &
                               cublasConvertCharToV2(transa), cublasConvertCharToV2(transb), &
                               m, n, k, &
                               alpha, a, lda, b, ldb, &
                               beta, c, ldc)
        !$omp end target data
      else
        cuErr = cublasCreate(cublas_handle)
        !$omp target data use_device_ptr(a, b, c)
        call cublasZgemm(transa, transb, &
                         m, n, k, &
                         alpha, a, lda, b, ldb, &
                         beta, c, ldc)
        !$omp end target data
        cuErr = cublasDestroy(cublas_handle)
      end if
#else
      call die("OpenMP Target version of accel_zgemm requested, but OpenMP Target not "&
            &"compiled into this executable", only_root_writes = .true.)
#endif
    case (CPU_ALGO)
      call zgemm(transa, transb, &
                 m, n, k, &
                 alpha, a, lda, b, ldb, &
                 beta, c, ldc)
    case default
      call die("Invald algorithm for accel_zgemm", only_root_writes = .true.)
    end select

  end subroutine accel_zgemm_core

  subroutine accel_dgemv_core(trans, &
                              m, n, &
                              alpha, &
                              a, lda, &
                              x, incx, &
                              beta, &
                              y, incy, &
                              algo)
    implicit none

    double precision :: alpha, beta
    character :: trans
    integer :: incx, incy, lda, m, n
    double precision :: a(*), x(*), y(*)
    integer(kind(CPU_ALGO)) :: algo

#if defined (OPENACC) || defined (OMP_TARGET)
    integer :: cuErr
    type(cublasHandle) :: cublas_handle
#endif

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      cuErr = cublasCreate(cublas_handle)
      !$acc host_data use_device(a, x, y)
      call cublasDgemv(trans, &
                       m, n, &
                       alpha, a, lda, x, incx, &
                       beta, y, incy)
      !$acc end host_data
      cuErr = cublasDestroy(cublas_handle)
#else
      call die("OpenACC version of accel_zgemv requested, but OpenACC not "&
            &"compiled into this executable", only_root_writes = .true.)
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      cuErr = cublasCreate(cublas_handle)
      !$omp target data use_device_ptr(a, x, y)
      call cublasDgemv(trans, &
                       m, n, &
                       alpha, a, lda, x, incx, &
                       beta, y, incy)
      !$omp end target data
      cuErr = cublasDestroy(cublas_handle)
#else
      call die("OpenMP Target version of accel_zgemv requested, but OpenMP Target not "&
            &"compiled into this executable", only_root_writes = .true.)
#endif
    case (CPU_ALGO)
      call dgemv(trans, &
                 m, n, &
                 alpha, a, lda, x, incx, &
                 beta, y, incy)
    case default
      call die("Invald algorithm for accel_zgemv", only_root_writes = .true.)
    end select

  end subroutine accel_dgemv_core

  subroutine accel_zgemv_core(trans, &
                              m, n, &
                              alpha, &
                              a, lda, &
                              x, incx, &
                              beta, &
                              y, incy, &
                              algo)
    implicit none

    double complex :: alpha, beta
    character :: trans
    integer :: incx, incy, lda, m, n
    double complex :: a(*), x(*), y(*)
    integer(kind(CPU_ALGO)) :: algo

#if defined (OPENACC) || defined (OMP_TARGET)
    integer :: cuErr
    type(cublasHandle) :: cublas_handle
#endif

    select case (algo)
    case (OPENACC_ALGO)
#ifdef OPENACC
      cuErr = cublasCreate(cublas_handle)
      !$acc host_data use_device(a, x, y)
      call cublasZgemv(trans, &
                       m, n, &
                       alpha, a, lda, x, incx, &
                       beta, y, incy)
      !$acc end host_data
      cuErr = cublasDestroy(cublas_handle)
#else
      call die("OpenACC version of accel_zgemv requested, but OpenACC not "&
            &"compiled into this executable", only_root_writes = .true.)
#endif
    case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
      cuErr = cublasCreate(cublas_handle)
      !$omp target data use_device_ptr(a, x, y)
      call cublasZgemv(trans, &
                       m, n, &
                       alpha, a, lda, x, incx, &
                       beta, y, incy)
      !$omp end target data
      cuErr = cublasDestroy(cublas_handle)
#else
      call die("OpenMP Target version of accel_zgemv requested, but OpenMP Target not "&
            &"compiled into this executable", only_root_writes = .true.)
#endif
    case (CPU_ALGO)
      call zgemv(trans, &
                 m, n, &
                 alpha, a, lda, x, incx, &
                 beta, y, incy)
    case default
      call die("Invald algorithm for accel_zgemv", only_root_writes = .true.)
    end select

  end subroutine accel_zgemv_core

  subroutine accel_zgemm_2d(transa, transb, &
                            m, n, k, &
                            alpha, &
                            a, lda, &
                            b, ldb, &
                            beta, &
                            c, ldc, &
                            algo)
    implicit none

    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double complex :: alpha, beta
    double complex :: a(:,:), b(:,:), c(:,:)
    integer(kind(CPU_ALGO)) :: algo

    call accel_zgemm_core(transa, transb, &
                          m, n, k, &
                          alpha, &
                          a, lda, &
                          b, ldb, &
                          beta, &
                          c, ldc, &
                          algo)

  end subroutine accel_zgemm_2d

  subroutine accel_dgemm_2d(transa, transb, &
                            m, n, k, &
                            alpha, &
                            a, lda, &
                            b, ldb, &
                            beta, &
                            c, ldc, &
                            algo)
    implicit none

    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double precision :: alpha, beta
    double precision :: a(:,:), b(:,:), c(:,:)
    integer(kind(CPU_ALGO)) :: algo

    call accel_dgemm_core(transa, transb, &
                          m, n, k, &
                          alpha, &
                          a, lda, &
                          b, ldb, &
                          beta, &
                          c, ldc, &
                          algo)

  end subroutine accel_dgemm_2d

  subroutine accel_zgemm_3d_3d_1d(transa, transb, &
                                  m, n, k, &
                                  alpha, &
                                  a, lda, &
                                  b, ldb, &
                                  beta, &
                                  c, ldc, &
                                  algo)
    implicit none

    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double complex :: alpha, beta
    double complex :: a(:,:,:), b(:,:,:), c(:)
    integer(kind(CPU_ALGO)) :: algo

    call accel_zgemm_core(transa, transb, &
                          m, n, k, &
                          alpha, &
                          a, lda, &
                          b, ldb, &
                          beta, &
                          c, ldc, &
                          algo)

  end subroutine accel_zgemm_3d_3d_1d

  subroutine accel_dgemm_3d_3d_1d(transa, transb, &
                                  m, n, k, &
                                  alpha, &
                                  a, lda, &
                                  b, ldb, &
                                  beta, &
                                  c, ldc, &
                                  algo)
    implicit none

    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double precision :: alpha, beta
    double precision :: a(:,:,:), b(:,:,:), c(:)
    integer(kind(CPU_ALGO)) :: algo

    call accel_dgemm_core(transa, transb, &
                          m, n, k, &
                          alpha, &
                          a, lda, &
                          b, ldb, &
                          beta, &
                          c, ldc, &
                          algo)

  end subroutine accel_dgemm_3d_3d_1d

  subroutine accel_zgemm_2d_5d_5d(transa, transb, &
                                  m, n, k, &
                                  alpha, &
                                  a, lda, &
                                  b, ldb, &
                                  beta, &
                                  c, ldc, &
                                  algo)
    implicit none

    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double complex :: alpha, beta
    double complex :: a(:,:), b(:,:,:,:,:), c(:,:,:,:,:)
    integer(kind(CPU_ALGO)) :: algo

    call accel_zgemm_core(transa, transb, &
                          m, n, k, &
                          alpha, &
                          a, lda, &
                          b, ldb, &
                          beta, &
                          c, ldc, &
                          algo)

  end subroutine accel_zgemm_2d_5d_5d

  subroutine accel_dgemm_2d_5d_5d(transa, transb, &
                                  m, n, k, &
                                  alpha, &
                                  a, lda, &
                                  b, ldb, &
                                  beta, &
                                  c, ldc, &
                                  algo)
    implicit none

    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double precision :: alpha, beta
    double precision :: a(:,:), b(:,:,:,:,:), c(:,:,:,:,:)
    integer(kind(CPU_ALGO)) :: algo

    call accel_dgemm_core(transa, transb, &
                          m, n, k, &
                          alpha, &
                          a, lda, &
                          b, ldb, &
                          beta, &
                          c, ldc, &
                          algo)

  end subroutine accel_dgemm_2d_5d_5d

  subroutine accel_dgemv_3d_3d_1d(trans, &
                                  m, n, &
                                  alpha, &
                                  a, lda, &
                                  x, incx, &
                                  beta, &
                                  y, incy, &
                                  algo)
    implicit none

    double precision :: alpha, beta
    character :: trans
    integer :: incx, incy, lda, m, n
    double precision :: a(:,:,:), x(:,:,:), y(:)
    integer(kind(CPU_ALGO)) :: algo

    call accel_dgemv_core(trans, &
                          m, n, &
                          alpha, &
                          a, lda, &
                          x, incx, &
                          beta, &
                          y, incy, &
                          algo)

  end subroutine accel_dgemv_3d_3d_1d

  subroutine accel_zgemv_3d_3d_1d(trans, &
                                  m, n, &
                                  alpha, &
                                  a, lda, &
                                  x, incx, &
                                  beta, &
                                  y, incy, &
                                  algo)
    implicit none

    double complex :: alpha, beta
    character :: trans
    integer :: incx, incy, lda, m, n
    double complex :: a(:,:,:), x(:,:,:), y(:)
    integer(kind(CPU_ALGO)) :: algo

    call accel_zgemv_core(trans, &
                          m, n, &
                          alpha, &
                          a, lda, &
                          x, incx, &
                          beta, &
                          y, incy, &
                          algo)

  end subroutine accel_zgemv_3d_3d_1d

end module accel_linalg_m
