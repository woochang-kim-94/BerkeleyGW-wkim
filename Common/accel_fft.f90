#include "f_defs.h"

module accel_fft_m

  use ISO_C_BINDING

  use accel_memory_m
  use algos_common_m
  use message_m
  use nrtype_m
  use peinfo_m,   only: peinf
  use push_pop_m

#if defined (OPENACC) || defined (OMP_TARGET)
  use openacc
  use cufft
#endif

  implicit none

  private

  public :: allocate_accel_mtxel_sig, deallocate_accel_mtxel_sig
  public :: allocate_accel_mtxel_eps, deallocate_accel_mtxel_eps
  public :: allocate_accel_mtxel_ker, deallocate_accel_mtxel_ker
  public :: accel_zero_box, accel_put_into_fftbox, accel_box_multiply, accel_get_from_fftbox
  public :: accel_run_fft
  public :: accel_zero_box_batched, accel_box_multiply_batched, accel_run_fft_batched, &
            accel_get_from_fftbox_batched, accel_put_into_fftbox_batched

  interface accel_put_into_fftbox_batched
    module procedure accel_put_into_fftbox_real_batched_ker
    module procedure accel_put_into_fftbox_cplx_batched_ker
    module procedure accel_put_into_fftbox_real_batched_eps
    module procedure accel_put_into_fftbox_cplx_batched_eps
  end interface

  interface accel_put_into_fftbox
    module procedure accel_put_into_fftbox_real
    module procedure accel_put_into_fftbox_cplx
  end interface

  interface accel_get_from_fftbox_batched
    module procedure accel_get_from_fftbox_batched_real
    module procedure accel_get_from_fftbox_batched_cplx
  end interface

  interface accel_get_from_fftbox
    module procedure accel_get_from_fftbox_real
    module procedure accel_get_from_fftbox_cplx
  end interface

  ! Type definition
  type :: fft_box_type
    complex(DPC), allocatable :: box(:,:,:)
  end type

  type :: fft_vec_type
    complex(DPC), allocatable :: vec(:)
  end type

  type :: fft_vec_type_scalar
    SCALAR, allocatable :: vec(:)
  end type

  type :: accel_mtxel_sig_type
    ! gvec info
    integer, allocatable :: gvec_comp(:,:)
    integer, allocatable :: gvec_isrtx_k(:), gvec_isrtx_q(:), gvec_isrtx_mtxel(:)
    ! eqp aux function (outer WFN)
    complex(DPC), allocatable :: eqp_box(:,:,:), eqp_vec(:)
    ! bands aux function (inner WFN)
    type(fft_vec_type), allocatable :: bands_vec(:)
    ! output (MTXEL)
    type(fft_box_type), allocatable :: mtxel_box(:)
    type(fft_vec_type), allocatable :: mtxel_vec(:)
    ! fft plan
    integer              :: eqp_plan
    integer, allocatable :: mtxel_plan(:)   ! do we need many plans??
    ! host stuff
    complex(DPC), allocatable :: fftbox_aux(:,:,:)
    SCALAR, allocatable       :: vec_aux(:)
    ! control size
    integer :: Nfft(3)
    integer :: mtxel_band_block_size
  end type accel_mtxel_sig_type

  type :: accel_mtxel_eps_type
    ! control size
    integer :: Nfft(3)
    ! ordering information
    integer, allocatable :: gvec_comp(:,:)
    integer, allocatable :: inner_sort(:)
    integer, allocatable :: outer_sort(:)
    integer, allocatable :: result_sort(:)
    ! FFT plan
    integer :: outer_plan
    integer :: inner_plan_batched
    integer, allocatable :: inner_plan(:)
    ! FFT boxes
    complex(DPC), dimension(:,:,:), allocatable :: outer_fftbox
    complex(DPC), dimension(:,:,:), allocatable :: inner_fftbox
    complex(DPC), dimension(:,:,:,:), allocatable :: inner_fftbox_batched
    type(fft_box_type), allocatable :: inner_box(:)
    ! Bands
    SCALAR, allocatable :: outer_bands(:)
    SCALAR, allocatable :: outer_bands_2d(:, :)
    SCALAR, allocatable :: inner_bands(:)
    SCALAR, allocatable :: inner_bands_2d(:, :)
    ! Result
    SCALAR, dimension(:), allocatable :: result_array
    SCALAR, dimension(:,:), allocatable :: result_array_batched
    SCALAR, allocatable :: result_array_2d(:, :)
    type(fft_vec_type_scalar), allocatable :: result_a(:)
    ! block size if using streams
    integer :: mtxel_band_block_size
    integer :: n_ffts_per_batch  !< Number of FFTs to run in each batched FFT invocation
  end type accel_mtxel_eps_type

  type(accel_mtxel_sig_type), public :: accel_mtxel_sig
  type(accel_mtxel_eps_type), public :: accel_mtxel_eps

contains

  ! WPH: Adding the box as future-proofing for MKL offload
  subroutine accel_create_fft_plan(plan, Nfft, in_box, out_box)
    integer, intent(out) :: plan
    integer, intent(in) :: Nfft(3)
    complex(DPC), intent(in) :: in_box(:, :, :)
    complex(DPC), intent(in) :: out_box(:, :, :)

    integer :: ierr

#if defined (OPENACC) || defined (OMP_TARGET)
    ierr = cufftPlan3D(plan, Nfft(3), Nfft(2), Nfft(1), CUFFT_Z2Z)
#else
    call die_algos("OpenACC")
#endif
  end subroutine accel_create_fft_plan


  subroutine accel_create_fft_plan_batched(plan, Nfft, in_box, out_box, size_batch)
    integer, intent(out) :: plan
    integer, intent(in) :: Nfft(3)
    complex(DPC), intent(in) :: in_box(:, :, :, :)
    complex(DPC), intent(in) :: out_box(:, :, :, :)
    integer, intent(in) :: size_batch

    integer, dimension(3) :: Nfft_plan
    integer :: ierr, iodist

    iodist = Nfft(1) * Nfft(2) * Nfft(3)
    ! cuFFT plan creation reverses the order of the FFT grid
    Nfft_plan(1) = Nfft(3)
    Nfft_plan(2) = Nfft(2)
    Nfft_plan(3) = Nfft(1)

#if defined (OPENACC) || defined (OMP_TARGET)
    ierr = cufftPlanMany(plan, 3, Nfft_plan, &
                         in_box, 1, iodist, &
                         out_box, 1, iodist, &
                         CUFFT_Z2Z, size_batch)
#else
    call die_algos("OpenACC")
#endif
  end subroutine accel_create_fft_plan_batched

  subroutine accel_destroy_fft_plan(plan)
    integer, intent(out) :: plan

    integer :: ierr

#if defined (OPENACC) || defined (OMP_TARGET)
    ierr = cufftDestroy(plan)
#else
    call die_algos("OpenACC")
#endif
  end subroutine

  subroutine accel_run_fft_batched(box, plan, fsign, algo)
    complex(DPC) :: box(:,:,:,:)
    integer :: plan, fsign
    integer(kind(CPU_ALGO)) :: algo

    integer :: ierr

    PUSH_SUB(accel_run_fft_batched)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if ( fsign > 0 ) then
          !$acc host_data use_device(box)
          iErr = cufftExecZ2Z(plan, box, box, CUFFT_INVERSE)
          !$acc end host_data
          call cudaDeviceSynchronize()
        else
          call die("fsign < 0 not implemented", only_root_writes = .true.)
        end if
#else
        call die_algos("OpenACC")
#endif
      case(OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        if ( fsign > 0 ) then
          !$omp target data use_device_ptr(box)
          iErr = cufftExecZ2Z(plan, box, box, CUFFT_INVERSE)
          !$omp end target data
          call cudaDeviceSynchronize()
        else
          call die("fsign < 0 not implemented", only_root_writes = .true.)
        end if
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_run_fft", only_root_writes = .true.)
    end select

    POP_SUB(accel_run_fft_batched)
  end subroutine accel_run_fft_batched

  subroutine accel_run_fft(box, plan, fsign, algo)
    complex(DPC) :: box(:,:,:)
    integer :: plan, fsign
    integer(kind(CPU_ALGO)) :: algo

    integer :: ierr

    PUSH_SUB(accel_run_fft)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if ( fsign > 0 ) then
          !$acc host_data use_device(box)
          iErr = cufftExecZ2Z(plan, box, box, CUFFT_INVERSE)
          !$acc end host_data
        else
          !$acc host_data use_device(box)
          iErr = cufftExecZ2Z(plan, box, box, CUFFT_FORWARD)
          !$acc end host_data
        end if
#else
        call die_algos("OpenACC")
#endif
      case(OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        if ( fsign > 0 ) then
          !$omp target data use_device_ptr(box)
          iErr = cufftExecZ2Z(plan, box, box, CUFFT_INVERSE)
          !$omp end target data
        else
          !$omp target data use_device_ptr(box)
          iErr = cufftExecZ2Z(plan, box, box, CUFFT_FORWARD)
          !$omp end target data
        end if
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_run_fft", only_root_writes = .true.)
    end select

    POP_SUB(accel_run_fft)
  end subroutine accel_run_fft

  subroutine accel_zero_box_batched(box, Nfft, how_many, algo, queue)
    complex(DPC) :: box(:,:,:,:)
    integer :: Nfft(3), how_many
    integer(kind(CPU_ALGO)) :: algo
    integer, optional :: queue

    integer :: i_batch, ix, iy, iz

    PUSH_SUB(accel_zero_box_batched)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          call die("queue not implemented", only_root_writes = .true.)
        else
          !$acc parallel present(box) vector_length(512)
          !$acc loop gang vector collapse(4)
          do i_batch = 1, how_many
            do iz = 1, Nfft(3)
              do iy = 1, Nfft(2)
                do ix = 1, Nfft(1)
                  box(ix,iy,iz,i_batch) = CMPLX(0.0d0,0.0d0)
                enddo
              enddo
            enddo
          enddo
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(4)
        do i_batch = 1, how_many
          do iz = 1, Nfft(3)
            do iy = 1, Nfft(2)
              do ix = 1, Nfft(1)
                box(ix,iy,iz,i_batch) = CMPLX(0.0d0,0.0d0)
              enddo
            enddo
          enddo
        enddo
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_zero_box", only_root_writes = .true.)
    end select

    POP_SUB(accel_zero_box_batched)
  end subroutine accel_zero_box_batched

  subroutine accel_zero_box(box, Nfft, algo, queue)
    complex(DPC) :: box(:,:,:)
    integer :: Nfft(3)
    integer(kind(CPU_ALGO)) :: algo
    integer, optional :: queue

    integer :: ix, iy, iz

    PUSH_SUB(accel_zero_box)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          !$acc parallel present(box,Nfft) async(queue) vector_length(512)
          !$acc loop gang vector collapse(3) 
          do iz = 1, Nfft(3)
            do iy = 1, Nfft(2)
              do ix = 1, Nfft(1)
                box(ix,iy,iz) = CMPLX(0.0d0,0.0d0)
              enddo
            enddo
          enddo
          !$acc end parallel
        else
          !$acc parallel present(box) vector_length(512)
          !$acc loop gang vector collapse(3) 
          do iz = 1, Nfft(3)
            do iy = 1, Nfft(2)
              do ix = 1, Nfft(1)
                box(ix,iy,iz) = CMPLX(0.0d0,0.0d0)
              enddo
            enddo
          enddo
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(3)
        do iz = 1, Nfft(3)
          do iy = 1, Nfft(2)
            do ix = 1, Nfft(1)
              box(ix,iy,iz) = CMPLX(0.0d0,0.0d0)
            enddo
          enddo
        enddo
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_zero_box", only_root_writes = .true.)
    end select

    POP_SUB(accel_zero_box)
  end subroutine accel_zero_box

  subroutine accel_box_multiply_batched(box_out, box_in, Nfft, conj, how_many, algo, queue)
    complex(DPC) :: box_out(:,:,:,:), box_in(:,:,:)
    integer :: Nfft(3), conj, how_many
    integer(kind(CPU_ALGO)) :: algo
    integer, optional :: queue

    integer :: i_batch, ix, iy, iz

    PUSH_SUB(accel_box_multiply_batched)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if ( conj > 0 ) then
          if (present(queue)) then
            call die("queue not implemented", only_root_writes = .true.)
          else
            !$acc parallel present(box_out, box_in, Nfft) vector_length(512)
            !$acc loop gang vector collapse(4)
            do i_batch = 1, how_many
              do iz = 1, Nfft(3)
                do iy = 1, Nfft(2)
                  do ix = 1, Nfft(1)
                    box_out(ix,iy,iz,i_batch) = CONJG( box_in(ix,iy,iz) ) * box_out(ix,iy,iz,i_batch)
                  enddo
                enddo
              enddo
            enddo
            !$acc end parallel
          end if
        else
          call die("conj < 0 not implemented", only_root_writes = .true.)
        end if
#else
        call die_algos("OpenACC")
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        if ( conj > 0 ) then
          !$omp target teams distribute parallel do collapse(4)
          do i_batch = 1, how_many
            do iz = 1, Nfft(3)
              do iy = 1, Nfft(2)
                do ix = 1, Nfft(1)
                  box_out(ix,iy,iz,i_batch) = CONJG( box_in(ix,iy,iz) ) * box_out(ix,iy,iz,i_batch)
                enddo
              enddo
            enddo
          enddo
          !$omp end target teams distribute parallel do
        else
          call die("conj < 0 not implemented", only_root_writes = .true.)
        end if
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_box_multiply", only_root_writes = .true.)
    end select

    POP_SUB(accel_box_multiply_batched)
  end subroutine accel_box_multiply_batched

  subroutine accel_box_multiply(box_out, box_in, Nfft, conj, algo, queue)
    complex(DPC) :: box_out(:,:,:), box_in(:,:,:)
    integer :: Nfft(3), conj
    integer(kind(CPU_ALGO)) :: algo
    integer, optional :: queue

    integer :: ix, iy, iz

    PUSH_SUB(accel_box_multiply)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if ( conj > 0 ) then
          if (present(queue)) then
            !$acc parallel present(box_out, box_in, Nfft) async(queue) vector_length(512)
            !$acc loop gang vector collapse(3) 
            do iz = 1, Nfft(3)
              do iy = 1, Nfft(2)
                do ix = 1, Nfft(1)
                  box_out(ix,iy,iz) = CONJG( box_in(ix,iy,iz) ) * box_out(ix,iy,iz)
                enddo
              enddo
            enddo
            !$acc end parallel
          else
            !$acc parallel present(box_out, box_in, Nfft) vector_length(512)
            !$acc loop gang vector collapse(3)
            do iz = 1, Nfft(3)
              do iy = 1, Nfft(2)
                do ix = 1, Nfft(1)
                  box_out(ix,iy,iz) = CONJG( box_in(ix,iy,iz) ) * box_out(ix,iy,iz)
                enddo
              enddo
            enddo
            !$acc end parallel
          end if
        else
          if (present(queue)) then
            !$acc parallel present(box_out, box_in, Nfft) async(queue) vector_length(512)
            !$acc loop gang vector collapse(3)
            do iz = 1, Nfft(3)
              do iy = 1, Nfft(2)
                do ix = 1, Nfft(1)
                  box_out(ix,iy,iz) = box_in(ix,iy,iz) * box_out(ix,iy,iz)
                enddo
              enddo
            enddo
            !$acc end parallel
          else
            !$acc parallel present(box_out, box_in, Nfft) vector_length(512)
            !$acc loop gang vector collapse(3)
            do iz = 1, Nfft(3)
              do iy = 1, Nfft(2)
                do ix = 1, Nfft(1)
                  box_out(ix,iy,iz) = box_in(ix,iy,iz) * box_out(ix,iy,iz)
                enddo
              enddo
            enddo
            !$acc end parallel
          end if
        end if
#else
        call die_algos("OpenACC")
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        if ( conj > 0 ) then
          !$omp target teams distribute parallel do collapse(3)
          do iz = 1, Nfft(3)
            do iy = 1, Nfft(2)
              do ix = 1, Nfft(1)
                box_out(ix,iy,iz) = CONJG( box_in(ix,iy,iz) ) * box_out(ix,iy,iz)
              enddo
            enddo
          enddo
          !$omp end target teams distribute parallel do
        else
          !$omp target teams distribute parallel do collapse(3)
          do iz = 1, Nfft(3)
            do iy = 1, Nfft(2)
              do ix = 1, Nfft(1)
                box_out(ix,iy,iz) = box_in(ix,iy,iz) * box_out(ix,iy,iz)
              enddo
            enddo
          enddo
          !$omp end target teams distribute parallel do
        end if
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_box_multiply", only_root_writes = .true.)
    end select

    POP_SUB(accel_box_multiply)
  end subroutine accel_box_multiply

  ! Used to populate initial FFT boxes for Kernel, where entire wavefunctions are populated at once
  ! and "batches" refer to batches over wavefunctions
  subroutine accel_put_into_fftbox_real_batched_ker(Ng, vec, g_comp, g_index, box, Nfft, alpha, offset, how_many, algo, queue)
    integer :: Ng, Nfft(3)
    real(DP) :: vec(:,:)
    integer :: g_comp(:,:), g_index(:)
    complex(DPC) :: box(:,:,:,:)
    real(DP) :: alpha
    integer :: offset
    integer :: how_many
    integer(kind(CPU_ALGO)) :: algo
    integer, optional :: queue

    integer :: iii, i_batch, bidx(3)

    PUSH_SUB(accel_put_into_fftbox_batched_ker)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          call die("queue not implemented", only_root_writes = .true.)
        else
          !$acc parallel private(bidx) present(g_comp,g_index,vec,box,Nfft) vector_length(512)
          !$acc loop gang vector collapse(2)
          do i_batch = 1, how_many
            do iii=1, Ng

              bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

              if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
              if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
              if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

              box(bidx(1), bidx(2), bidx(3), i_batch) = vec(iii, offset + i_batch) * alpha
            end do
          end do
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case(OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(bidx) collapse(2)
        do i_batch = 1, how_many
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            box(bidx(1), bidx(2), bidx(3), i_batch) = vec(iii, offset + i_batch) * alpha
          end do
        end do
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_put_into_fftbox", only_root_writes = .true.)
    end select

    POP_SUB(accel_put_into_fftbox_batched_ker)
  end subroutine accel_put_into_fftbox_real_batched_ker

  ! Used to populate initial FFT boxes for Epsilon, where pieces of wavefunctions (i.e. slices
  ! of the wavefunction array) are populated and "batches" refer to batches over these slices
  subroutine accel_put_into_fftbox_real_batched_eps(Ng, vec, g_comp, g_index, box, Nfft, alpha, offset, how_many, algo, queue)
    integer :: Ng, Nfft(3)
    real(DP) :: vec(:)
    integer :: g_comp(:,:), g_index(:)
    complex(DPC) :: box(:,:,:,:)
    real(DP) :: alpha
    integer :: offset(:)
    integer :: how_many
    integer(kind(CPU_ALGO)) :: algo
    integer, optional  :: queue

    integer :: iii, i_batch, bidx(3)

    PUSH_SUB(accel_put_into_fftbox_batched_eps)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          call die("queue not implemented", only_root_writes = .true.)
        else
          !$acc parallel loop private(bidx) present(g_comp,g_index,vec,box,Nfft,offset) collapse(2)
          do i_batch = 1, how_many
            do iii=1, Ng

              bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

              if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
              if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
              if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

              box(bidx(1), bidx(2), bidx(3), i_batch) = vec(offset(i_batch) + iii) * alpha
            end do
          end do
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case(OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(bidx) collapse(2)
        do i_batch = 1, how_many
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            box(bidx(1), bidx(2), bidx(3), i_batch) = vec(offset(i_batch) + iii) * alpha
          end do
        end do
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_put_into_fftbox", only_root_writes = .true.)
    end select

    POP_SUB(accel_put_into_fftbox_batched_eps)
  end subroutine accel_put_into_fftbox_real_batched_eps

  subroutine accel_put_into_fftbox_real(Ng, vec, g_comp, g_index, box, Nfft, alpha, algo, queue)
    integer :: Ng, Nfft(3)
    real(DP) :: vec(:)
    integer :: g_comp(:,:), g_index(:)
    complex(DPC) :: box(:,:,:)
    real(DP) :: alpha
    integer(kind(CPU_ALGO)) :: algo
    integer, optional  :: queue

    integer :: iii, bidx(3)

    PUSH_SUB(accel_put_into_fftbox)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          !$acc parallel loop private(bidx) present(g_comp,g_index,vec,box,Nfft) async(queue)
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            box(bidx(1),bidx(2),bidx(3)) = vec(iii) * alpha
          end do
          !$acc end parallel
        else
          !$acc parallel loop private(bidx) present(g_comp,g_index,vec,box,Nfft)
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            box(bidx(1),bidx(2),bidx(3)) = vec(iii) * alpha
          end do
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case(OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(bidx)
        do iii=1, Ng

          bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

          if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
          if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
          if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

          box(bidx(1),bidx(2),bidx(3)) = vec(iii) * alpha
        end do
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_put_into_fftbox", only_root_writes = .true.)
    end select

    POP_SUB(accel_put_into_fftbox)
  end subroutine accel_put_into_fftbox_real

  ! Used to populate initial FFT boxes for Kernel, where entire wavefunctions are populated at once
  ! and "batches" refer to batches over wavefunctions
  subroutine accel_put_into_fftbox_cplx_batched_ker(Ng, vec, g_comp, g_index, box, Nfft, alpha, offset, how_many, algo, queue)
    integer :: Ng, Nfft(3)
    complex(DPC) :: vec(:,:)
    integer :: g_comp(:,:), g_index(:)
    complex(DPC) :: box(:,:,:,:)
    real(DP) :: alpha
    integer :: offset
    integer :: how_many
    integer(kind(CPU_ALGO)) :: algo
    integer, optional :: queue

    integer :: iii, i_batch, bidx(3)

    PUSH_SUB(accel_put_into_fftbox_batched_ker)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          call die("queue not implemented", only_root_writes = .true.)
        else
          !$acc parallel private(bidx) present(g_comp,g_index,vec,box,Nfft) vector_length(512)
          !$acc loop gang vector collapse(2)
          do i_batch = 1, how_many
            do iii=1, Ng

              bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

              if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
              if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
              if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

              box(bidx(1), bidx(2), bidx(3), i_batch) = vec(iii, offset + i_batch) * alpha
            end do
          end do
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case(OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(bidx) collapse(2)
        do i_batch = 1, how_many
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            box(bidx(1), bidx(2), bidx(3), i_batch) = vec(iii, offset + i_batch) * alpha
          end do
        end do
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_put_into_fftbox", only_root_writes = .true.)
    end select

    POP_SUB(accel_put_into_fftbox_batched_ker)
  end subroutine accel_put_into_fftbox_cplx_batched_ker

  ! Used to populate initial FFT boxes for Epsilon, where pieces of wavefunctions (i.e. slices
  ! of the wavefunction array) are populated and "batches" refer to batches over these slices
  subroutine accel_put_into_fftbox_cplx_batched_eps(Ng, vec, g_comp, g_index, box, Nfft, alpha, offset, how_many, algo, queue)
    integer :: Ng, Nfft(3)
    complex(DPC) :: vec(:)
    integer :: g_comp(:,:), g_index(:)
    complex(DPC) :: box(:,:,:,:)
    real(DP) :: alpha
    integer :: offset(:)
    integer :: how_many
    integer(kind(CPU_ALGO)) :: algo
    integer, optional :: queue

    integer :: iii, i_batch, bidx(3)

    PUSH_SUB(accel_put_into_fftbox_batched)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          call die("queue not implemented", only_root_writes = .true.)
        else
          !$acc parallel private(bidx) present(g_comp,g_index,vec,box,Nfft,offset) vector_length(512)
          !$acc loop gang vector collapse(2)
          do i_batch = 1, how_many
            do iii=1, Ng

              bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

              if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
              if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
              if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

              box(bidx(1), bidx(2), bidx(3), i_batch) = vec(offset(i_batch) + iii) * alpha
            end do
          end do
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case(OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(bidx) collapse(2)
        do i_batch = 1, how_many
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            box(bidx(1), bidx(2), bidx(3), i_batch) = vec(offset(i_batch) + iii) * alpha
          end do
        end do
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_put_into_fftbox", only_root_writes = .true.)
    end select

    POP_SUB(accel_put_into_fftbox_batched)
  end subroutine accel_put_into_fftbox_cplx_batched_eps

  subroutine accel_put_into_fftbox_cplx(Ng, vec, g_comp, g_index, box, Nfft, alpha, algo, queue)
    integer :: Ng, Nfft(3)
    complex(DPC) :: vec(:)
    integer :: g_comp(:,:), g_index(:)
    complex(DPC) :: box(:,:,:)
    real(DP) :: alpha
    integer(kind(CPU_ALGO)) :: algo
    integer, optional :: queue

    integer :: iii, bidx(3)

    PUSH_SUB(accel_put_into_fftbox)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          !$acc parallel private(bidx) present(g_comp,g_index,vec,box,Nfft) async(queue) vector_length(512)
          !$acc loop gang vector 
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            box(bidx(1),bidx(2),bidx(3)) = vec(iii) * alpha
          end do
          !$acc end parallel
        else
          !$acc parallel private(bidx) present(g_comp,g_index,vec,box,Nfft) vector_length(512)
          !$acc loop gang vector
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            box(bidx(1),bidx(2),bidx(3)) = vec(iii) * alpha
          end do
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case(OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(bidx)
        do iii=1, Ng

          bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

          if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
          if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
          if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

          box(bidx(1),bidx(2),bidx(3)) = vec(iii) * alpha
        end do
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_put_into_fftbox", only_root_writes = .true.)
    end select

    POP_SUB(accel_put_into_fftbox)
  end subroutine accel_put_into_fftbox_cplx

  subroutine accel_get_from_fftbox_batched_real(Ng, vec, g_comp, g_index, box, Nfft, alpha, how_many, algo, queue)
    integer :: Ng, Nfft(3)
    real(DP) :: vec(:,:)
    integer :: g_comp(:,:), g_index(:)
    complex(DPC) :: box(:,:,:,:)
    real(DP) :: alpha
    integer(kind(CPU_ALGO)) :: algo
    integer :: how_many
    integer, optional  :: queue

    integer :: iii, i_batch, bidx(3)

    PUSH_SUB(accel_get_from_fftbox_batched)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          call die("queue not implemented", only_root_writes = .true.)
        else
          !$acc parallel loop private(bidx) present(g_comp,g_index,vec,box,Nfft) collapse(2)
          do i_batch = 1, how_many
            do iii=1, Ng

              bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

              if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
              if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
              if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

              vec(iii, i_batch) = box(bidx(1), bidx(2), bidx(3), i_batch) * alpha
            end do
          end do
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(bidx) collapse(2)
        do i_batch = 1, how_many
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            vec(iii, i_batch) = box(bidx(1), bidx(2), bidx(3), i_batch) * alpha
          end do
        end do
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_get_from_fftbox", only_root_writes = .true.)
    end select

    POP_SUB(accel_get_from_fftbox_batched)
  end subroutine accel_get_from_fftbox_batched_real

  subroutine accel_get_from_fftbox_real(Ng, vec, g_comp, g_index, box, Nfft, alpha, algo, queue)
    integer :: Ng, Nfft(3)
    real(DP) :: vec(:)
    integer :: g_comp(:,:), g_index(:)
    complex(DPC) :: box(:,:,:)
    real(DP) :: alpha
    integer(kind(CPU_ALGO)) :: algo
    integer, optional  :: queue

    integer :: iii, bidx(3)

    PUSH_SUB(accel_get_from_fftbox)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          !$acc parallel loop private(bidx) present(g_comp,g_index,vec,box,Nfft) async(queue)
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            vec(iii) = box(bidx(1),bidx(2),bidx(3)) * alpha
          end do
          !$acc end parallel
        else
          !$acc parallel loop private(bidx) present(g_comp,g_index,vec,box,Nfft)
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            vec(iii) = box(bidx(1),bidx(2),bidx(3)) * alpha
          end do
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(bidx)
        do iii=1, Ng

          bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

          if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
          if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
          if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

          vec(iii) = box(bidx(1),bidx(2),bidx(3)) * alpha
        end do
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_get_from_fftbox", only_root_writes = .true.)
    end select

    POP_SUB(accel_get_from_fftbox)
  end subroutine accel_get_from_fftbox_real

  subroutine accel_get_from_fftbox_batched_cplx(Ng, vec, g_comp, g_index, box, Nfft, alpha, how_many, algo, queue)
    integer :: Ng, Nfft(3)
    complex(DPC) :: vec(:,:)
    integer :: g_comp(:,:), g_index(:)
    complex(DPC) :: box(:,:,:,:)
    real(DP) :: alpha
    integer :: how_many
    integer(kind(CPU_ALGO)) :: algo
    integer, optional :: queue

    integer :: iii, i_batch, bidx(3)

    PUSH_SUB(accel_get_from_fftbox)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          call die("queue not implemented", only_root_writes = .true.)
        else
          !$acc parallel private(bidx) present(g_comp,g_index,vec,box,Nfft) vector_length(512)
          !$acc loop gang vector collapse(2)
          do i_batch = 1, how_many
            do iii=1, Ng

              bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

              if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
              if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
              if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

              vec(iii, i_batch) = box(bidx(1), bidx(2), bidx(3), i_batch) * alpha
            end do
          end do
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(bidx) collapse(2)
        do i_batch = 1, how_many
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            vec(iii, i_batch) = box(bidx(1), bidx(2), bidx(3), i_batch) * alpha
          end do
        end do
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_get_from_fftbox", only_root_writes = .true.)
    end select

    POP_SUB(accel_get_from_fftbox_batched)
  end subroutine accel_get_from_fftbox_batched_cplx

  subroutine accel_get_from_fftbox_cplx(Ng, vec, g_comp, g_index, box, Nfft, alpha, algo, queue)
    integer :: Ng, Nfft(3)
    complex(DPC) :: vec(:)
    integer :: g_comp(:,:), g_index(:)
    complex(DPC) :: box(:,:,:)
    real(DP) :: alpha
    integer(kind(CPU_ALGO)) :: algo
    integer, optional :: queue

    integer :: iii, bidx(3)

    PUSH_SUB(accel_get_from_fftbox)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        if (present(queue)) then
          !$acc parallel private(bidx) present(g_comp,g_index,vec,box,Nfft) async(queue) vector_length(512)
          !$acc loop gang vector
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            vec(iii) = box(bidx(1),bidx(2),bidx(3)) * alpha
          end do
          !$acc end parallel
        else
          !$acc parallel private(bidx) present(g_comp,g_index,vec,box,Nfft) vector_length(512)
          !$acc loop gang vector
          do iii=1, Ng

            bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

            if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
            if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
            if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

            vec(iii) = box(bidx(1),bidx(2),bidx(3)) * alpha
          end do
          !$acc end parallel
        end if
#else
        call die_algos("OpenACC")
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(bidx)
        do iii=1, Ng

          bidx(1:3) = g_comp(1:3, g_index(iii)) + 1

          if (g_comp(1,g_index(iii)) < 0) bidx(1) = Nfft(1) + bidx(1)
          if (g_comp(2,g_index(iii)) < 0) bidx(2) = Nfft(2) + bidx(2)
          if (g_comp(3,g_index(iii)) < 0) bidx(3) = Nfft(3) + bidx(3)

          vec(iii) = box(bidx(1),bidx(2),bidx(3)) * alpha
        end do
        !$omp end target teams distribute parallel do
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for accel_get_from_fftbox", only_root_writes = .true.)
    end select

    POP_SUB(accel_get_from_fftbox)
  end subroutine accel_get_from_fftbox_cplx

  subroutine allocate_accel_mtxel_eps_common(Nfft, size_inner, size_sort, &
                                             size_result, size_batch, &
                                             size_inner_bands, algo)
    implicit none

    integer, intent(in) :: Nfft(3)
    integer, intent(in) :: size_inner
    integer, intent(in) :: size_sort
    integer, intent(in) :: size_result
    integer, intent(in) :: size_inner_bands
    integer, intent(in) :: size_batch
    integer, intent(in) :: algo

    integer :: n1, n2, n3
    complex(DPC), dimension(:,:,:,:), allocatable :: fftbox2

    PUSH_SUB(allocate_accel_mtxel_eps_ommon)

    if (algo == CPU_ALGO) return

    accel_mtxel_eps%n_ffts_per_batch = size_batch

    n1 = Nfft(1)
    n2 = Nfft(2)
    n3 = Nfft(3)

    SAFE_ALLOCATE (accel_mtxel_eps%gvec_comp, (3, size_sort))
    accel_mtxel_eps%gvec_comp = 0

    SAFE_ALLOCATE (accel_mtxel_eps%result_sort, (size_sort))
    accel_mtxel_eps%result_sort = 0

    SAFE_ALLOCATE (accel_mtxel_eps%inner_sort, (size_sort))
    accel_mtxel_eps%inner_sort = 0

    SAFE_ALLOCATE (accel_mtxel_eps%outer_sort, (size_sort))
    accel_mtxel_eps%outer_sort = 0

    SAFE_ALLOCATE(accel_mtxel_eps%inner_fftbox, (n1, n2, n3))
    accel_mtxel_eps%inner_fftbox = cmplx(0.0D+00, 0.0D+00)

    SAFE_ALLOCATE(accel_mtxel_eps%outer_fftbox, (n1, n2, n3))
    accel_mtxel_eps%outer_fftbox = cmplx(0.0D+00, 0.0D+00)

    SAFE_ALLOCATE (accel_mtxel_eps%inner_bands_2d, (size_inner, size_inner_bands))
    accel_mtxel_eps%inner_bands_2d = ZERO

    accel_mtxel_eps%Nfft = Nfft

    SAFE_ALLOCATE(accel_mtxel_eps%inner_fftbox_batched, (n1, n2, n3, size_batch))
    accel_mtxel_eps%inner_fftbox_batched = cmplx(0.0D+00, 0.0D+00)

    SAFE_ALLOCATE (accel_mtxel_eps%result_array_batched, (size_result, size_batch))
    accel_mtxel_eps%result_array_batched = ZERO

    call accel_enter_data_map_to(accel_mtxel_eps%gvec_comp, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%result_sort, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%inner_sort, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%outer_sort, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%result_array_batched, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%inner_fftbox, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%inner_fftbox_batched, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%outer_fftbox, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%inner_bands_2d, algo)
    call accel_enter_data_map_alloc(accel_mtxel_eps%Nfft, algo)

    call accel_create_fft_plan(accel_mtxel_eps%outer_plan, Nfft, &
                               accel_mtxel_eps%outer_fftbox, accel_mtxel_eps%outer_fftbox)
    ! TODO:  Figure out why dummy array works but the real one doesn't!
    call accel_create_fft_plan_batched(accel_mtxel_eps%inner_plan_batched, Nfft, &
                                       fftbox2, fftbox2, size_batch)
  end subroutine allocate_accel_mtxel_eps_common

  subroutine deallocate_accel_mtxel_eps_common(algo)
      implicit none

      integer, intent(in) :: algo

      if (algo == CPU_ALGO) return

      accel_mtxel_eps%n_ffts_per_batch = 0

      call accel_exit_data_map_delete(accel_mtxel_eps%gvec_comp, algo)
      call accel_exit_data_map_delete(accel_mtxel_eps%result_sort, algo)
      call accel_exit_data_map_delete(accel_mtxel_eps%inner_sort, algo)
      call accel_exit_data_map_delete(accel_mtxel_eps%outer_sort, algo)
      call accel_exit_data_map_delete(accel_mtxel_eps%result_array_batched, algo)
      call accel_exit_data_map_delete(accel_mtxel_eps%inner_fftbox, algo)
      call accel_exit_data_map_delete(accel_mtxel_eps%inner_fftbox_batched, algo)
      call accel_exit_data_map_delete(accel_mtxel_eps%outer_fftbox, algo)
      call accel_exit_data_map_delete(accel_mtxel_eps%inner_bands_2d, algo)
      call accel_exit_data_map_delete(accel_mtxel_eps%Nfft, algo)

      call accel_destroy_fft_plan(accel_mtxel_eps%outer_plan)
      call accel_destroy_fft_plan(accel_mtxel_eps%inner_plan_batched)

      deallocate(accel_mtxel_eps%gvec_comp)
      deallocate(accel_mtxel_eps%result_sort)
      deallocate(accel_mtxel_eps%inner_sort)
      deallocate(accel_mtxel_eps%outer_sort)
      deallocate(accel_mtxel_eps%result_array_batched)
      deallocate(accel_mtxel_eps%inner_fftbox)
      deallocate(accel_mtxel_eps%inner_fftbox_batched)
      deallocate(accel_mtxel_eps%outer_fftbox)
      deallocate(accel_mtxel_eps%inner_bands_2d)
  end subroutine deallocate_accel_mtxel_eps_common

  subroutine allocate_accel_mtxel_ker(Nfft, gvec_components, size_outer, &
                                      size_inner, size_sort, &
                                      size_result, size_outer_bands, &
                                      size_inner_bands, size_batch, algo)
    implicit none

    integer, intent(in) :: Nfft(3)
    integer, intent(in) :: gvec_components(:, :)
    integer, intent(in) :: size_outer
    integer, intent(in) :: size_inner
    integer, intent(in) :: size_sort
    integer, intent(in) :: size_result
    integer, intent(in) :: size_outer_bands
    integer, intent(in) :: size_inner_bands
    integer, intent(in) :: size_batch
    integer, intent(in) :: algo

    if (algo == CPU_ALGO) return

    call allocate_accel_mtxel_eps_common(Nfft, size_inner, size_sort, size_result, &
                                         size_batch, size_inner_bands, algo)

    SAFE_ALLOCATE(accel_mtxel_eps%result_array_2d, (size_result, size_inner_bands))
    accel_mtxel_eps%result_array_2d = 0
    SAFE_ALLOCATE(accel_mtxel_eps%outer_bands_2d, (size_outer, size_outer_bands))
    accel_mtxel_eps%outer_bands_2d = 0

    call accel_enter_data_map_to(accel_mtxel_eps%result_array_2d, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%outer_bands_2d, algo)

    accel_mtxel_eps%gvec_comp = 0
    accel_mtxel_eps%gvec_comp(1:3, 1:size_sort) = gvec_components(1:3, 1:size_sort)
    call accel_update_to(accel_mtxel_eps%gvec_comp, algo)
  end subroutine allocate_accel_mtxel_ker

  subroutine deallocate_accel_mtxel_ker(algo)
      implicit none

      integer, intent(in) :: algo

      if (algo == CPU_ALGO) return

      call deallocate_accel_mtxel_eps_common(algo)

      call accel_exit_data_map_delete(accel_mtxel_eps%result_array_2d, algo)
      call accel_exit_data_map_delete(accel_mtxel_eps%outer_bands_2d, algo)

      deallocate(accel_mtxel_eps%result_array_2d)
      deallocate(accel_mtxel_eps%outer_bands_2d)
  end subroutine deallocate_accel_mtxel_ker

  subroutine allocate_accel_mtxel_eps(Nfft, size_outer, size_inner, size_sort, &
                                    size_result, size_spin, size_batch, algo)
    implicit none

    integer, intent(in) :: Nfft(3)
    integer, intent(in) :: size_outer
    integer, intent(in) :: size_inner
    integer, intent(in) :: size_sort
    integer, intent(in) :: size_result
    integer, intent(in) :: size_spin
    integer, intent(in) :: size_batch
    integer, intent(in) :: algo

    integer :: n1, n2, n3
    integer :: band_block_size, n1_loc
    integer :: ierr

    PUSH_SUB(allocate_accel_mtxel_eps)

    call allocate_accel_mtxel_eps_common(Nfft, size_inner, size_sort, size_result, &
                                         size_batch, size_spin, algo)

    n1 = Nfft(1)
    n2 = Nfft(2)
    n3 = Nfft(3)

    band_block_size = accel_mtxel_eps%mtxel_band_block_size

    SAFE_ALLOCATE (accel_mtxel_eps%result_array, (size_result))
    accel_mtxel_eps%result_array = ZERO

    SAFE_ALLOCATE (accel_mtxel_eps%outer_bands, (size_outer))
    accel_mtxel_eps%outer_bands = ZERO

    SAFE_ALLOCATE (accel_mtxel_eps%inner_bands, (size_inner))
    accel_mtxel_eps%inner_bands = ZERO

    call accel_enter_data_map_to(accel_mtxel_eps%outer_bands, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%inner_bands, algo)
    call accel_enter_data_map_to(accel_mtxel_eps%result_array, algo)

    ALLOCATE(accel_mtxel_eps%result_a   (band_block_size) )
    ALLOCATE(accel_mtxel_eps%inner_box  (band_block_size) )
    ALLOCATE(accel_mtxel_eps%inner_plan (band_block_size) )
    do n1_loc = 1, band_block_size
      ALLOCATE(accel_mtxel_eps%result_a(n1_loc)%vec (size_result) )
      accel_mtxel_eps%result_a(n1_loc)%vec  = cmplx(0.0D+00, 0.0D+00)

      ALLOCATE(accel_mtxel_eps%inner_box(n1_loc)%box (n1, n2, n3) )
      accel_mtxel_eps%inner_box(n1_loc)%box = cmplx(0.0D+00, 0.0D+00)
    end do
    accel_mtxel_eps%inner_plan = 0

    if (algo == OPENACC_ALGO) then
#ifdef OPENACC
      !$acc enter data copyin(accel_mtxel_eps%result_a, accel_mtxel_eps%inner_box)
      do n1_loc = 1, band_block_size
        !$acc enter data copyin( accel_mtxel_eps%result_a(n1_loc)%vec )
        !$acc enter data copyin( accel_mtxel_eps%inner_box(n1_loc)%box )
      end do
#else
        call die_algos("OpenACC")
#endif
    end if

#if defined (OPENACC) || defined (OMP_TARGET)
    if ( algo == OPENACC_ALGO ) then
      ierr = cufftSetStream(accel_mtxel_eps%outer_plan, acc_get_cuda_stream(1))
    end if
    do n1_loc = 1, band_block_size
      call accel_create_fft_plan(accel_mtxel_eps%inner_plan(n1_loc), Nfft, &
                                 accel_mtxel_eps%inner_fftbox, accel_mtxel_eps%inner_fftbox)
      if ( algo == OPENACC_ALGO ) then
        ierr = cufftSetStream(accel_mtxel_eps%inner_plan(n1_loc), acc_get_cuda_stream(n1_loc + 1))
      end if
    end do
#endif

    POP_SUB(allocate_accel_mtxel_eps)
  end subroutine allocate_accel_mtxel_eps

  subroutine deallocate_accel_mtxel_eps(algo)
    implicit none

    integer, intent(in) :: algo

    integer :: n1_loc
    integer :: ierr

    PUSH_SUB(deallocate_accel_mtxel_eps)

    call deallocate_accel_mtxel_eps_common(algo)

    call accel_exit_data_map_delete(accel_mtxel_eps%outer_bands, algo)
    call accel_exit_data_map_delete(accel_mtxel_eps%inner_bands, algo)
    call accel_exit_data_map_delete(accel_mtxel_eps%result_array, algo)

    SAFE_DEALLOCATE(accel_mtxel_eps%outer_bands)
    SAFE_DEALLOCATE(accel_mtxel_eps%inner_bands)
    SAFE_DEALLOCATE(accel_mtxel_eps%result_array)

    if (algo == OPENACC_ALGO) then
#ifdef OPENACC
        do n1_loc = 1, accel_mtxel_eps%mtxel_band_block_size
          !$acc exit data delete( accel_mtxel_eps%result_a(n1_loc)%vec )
          !$acc exit data delete( accel_mtxel_eps%inner_box(n1_loc)%box )
        end do
        !$acc exit data delete(accel_mtxel_eps%result_a)
        !$acc exit data delete(accel_mtxel_eps%inner_box)
        !$acc exit data delete(accel_mtxel_eps)
#else
        call die_algos("OpenACC")
#endif
    end if

    do n1_loc = 1, accel_mtxel_eps%mtxel_band_block_size
      DEALLOCATE(accel_mtxel_eps%result_a(n1_loc)%vec )
      DEALLOCATE(accel_mtxel_eps%inner_box(n1_loc)%box )
    end do
    DEALLOCATE(accel_mtxel_eps%result_a )
    DEALLOCATE(accel_mtxel_eps%inner_box )
    DEALLOCATE(accel_mtxel_eps%inner_plan)

    POP_SUB(deallocate_accel_mtxel_eps)
  end subroutine deallocate_accel_mtxel_eps

  subroutine allocate_accel_mtxel_sig ( Nfft, NG_q, NG_k, ng, ncoul, algo )
    integer, intent(in)      :: Nfft(3), NG_q, NG_k, ng, ncoul
    integer, intent(in)      :: algo

    integer :: nbands, n1, n2, n3
    integer :: nloc, ierr

    PUSH_SUB(allocate_accel_mtxel_sig)

    nbands = accel_mtxel_sig%mtxel_band_block_size
    n1 = Nfft(1)
    n2 = Nfft(2)
    n3 = Nfft(3)

    accel_mtxel_sig%Nfft = Nfft
    SAFE_ALLOCATE ( accel_mtxel_sig%gvec_comp, (3,ng) )
    SAFE_ALLOCATE ( accel_mtxel_sig%gvec_isrtx_k, (ng) )
    SAFE_ALLOCATE ( accel_mtxel_sig%gvec_isrtx_q, (ng) )
    SAFE_ALLOCATE ( accel_mtxel_sig%gvec_isrtx_mtxel, (ng) )
    accel_mtxel_sig%gvec_comp        = 0
    accel_mtxel_sig%gvec_isrtx_k     = 0
    accel_mtxel_sig%gvec_isrtx_q     = 0
    accel_mtxel_sig%gvec_isrtx_mtxel = 0

    SAFE_ALLOCATE ( accel_mtxel_sig%bands_vec, (nbands) )
    SAFE_ALLOCATE ( accel_mtxel_sig%mtxel_vec, (nbands) )
    SAFE_ALLOCATE ( accel_mtxel_sig%mtxel_box, (nbands) )
    SAFE_ALLOCATE ( accel_mtxel_sig%mtxel_plan, (nbands) )

    SAFE_ALLOCATE ( accel_mtxel_sig%eqp_vec, (NG_q) )
    accel_mtxel_sig%eqp_vec = CMPLX(0.0D+00,0.0D+00)
    SAFE_ALLOCATE ( accel_mtxel_sig%eqp_box, (n1,n2,n3) )
    accel_mtxel_sig%eqp_box = CMPLX(0.0D+00,0.0D+00)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        !$acc enter data copyin( accel_mtxel_sig )
        !$acc enter data copyin( accel_mtxel_sig%gvec_comp )
        !$acc enter data copyin( accel_mtxel_sig%gvec_isrtx_k, accel_mtxel_sig%gvec_isrtx_q, accel_mtxel_sig%gvec_isrtx_mtxel )
        !$acc enter data copyin( accel_mtxel_sig%bands_vec, accel_mtxel_sig%mtxel_box, accel_mtxel_sig%mtxel_vec )
        !$acc enter data copyin( accel_mtxel_sig%eqp_vec, accel_mtxel_sig%eqp_box )
#else
        call die_algos("OpenACC")
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target enter data map (to: accel_mtxel_sig )
        !$omp target enter data map (to: accel_mtxel_sig%gvec_comp)
        !$omp target enter data map (to: accel_mtxel_sig%gvec_isrtx_k, accel_mtxel_sig%gvec_isrtx_q, accel_mtxel_sig%gvec_isrtx_mtxel )
        !$omp target enter data map (to: accel_mtxel_sig%bands_vec, accel_mtxel_sig%mtxel_box, accel_mtxel_sig%mtxel_vec )
        !$omp target enter data map (to: accel_mtxel_sig%eqp_vec, accel_mtxel_sig%eqp_box )
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for allocate_accel_mtxel_sig", only_root_writes = .true.)
    end select

    ! is this reverse order?
#if defined (OPENACC) || defined (OMP_TARGET)
    call accel_create_fft_plan(accel_mtxel_sig%eqp_plan, Nfft, &
                               accel_mtxel_sig%eqp_box, accel_mtxel_sig%eqp_box)
    if (algo == OPENACC_ALGO) then
      ierr = cufftSetStream(accel_mtxel_sig%eqp_plan, acc_get_cuda_stream(1))
    end if
#endif

    do nloc = 1, nbands
      SAFE_ALLOCATE ( accel_mtxel_sig%bands_vec(nloc)%vec, (NG_k) )
      SAFE_ALLOCATE ( accel_mtxel_sig%mtxel_vec(nloc)%vec, (ncoul) )
      SAFE_ALLOCATE ( accel_mtxel_sig%mtxel_box(nloc)%box, (n1,n2,n3) )
      accel_mtxel_sig%bands_vec(nloc)%vec = CMPLX(0.0D+00,0.0D+00)
      accel_mtxel_sig%mtxel_vec(nloc)%vec = CMPLX(0.0D+00,0.0D+00)
      accel_mtxel_sig%mtxel_box(nloc)%box = CMPLX(0.0D+00,0.0D+00)
      select case(algo)
        case (OPENACC_ALGO)
#ifdef OPENACC
          !$acc enter data copyin( accel_mtxel_sig%bands_vec(nloc)%vec, accel_mtxel_sig%mtxel_vec(nloc)%vec )
          !$acc enter data copyin( accel_mtxel_sig%mtxel_box(nloc)%box )
#else
        call die_algos("OpenACC")
#endif
        case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
          !$omp target enter data map (to: accel_mtxel_sig%bands_vec(nloc)%vec, accel_mtxel_sig%mtxel_vec(nloc)%vec )
          !$omp target enter data map (to: accel_mtxel_sig%mtxel_box(nloc)%box )
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for allocate_accel_mtxel_sig", only_root_writes = .true.)
      end select

#if defined (OPENACC) || defined (OMP_TARGET)
      call accel_create_fft_plan(accel_mtxel_sig%mtxel_plan(nloc), Nfft, &
                                 accel_mtxel_sig%mtxel_box(nloc)%box, &
                                 accel_mtxel_sig%mtxel_box(nloc)%box)
      if (algo == OPENACC_ALGO) then
        ierr = cufftSetStream(accel_mtxel_sig%mtxel_plan(nloc), acc_get_cuda_stream(nloc+1))
      end if
#endif
    end do

    POP_SUB(allocate_accel_mtxel_sig)
  end subroutine allocate_accel_mtxel_sig

  subroutine deallocate_accel_mtxel_sig(algo)
    implicit none

    integer, intent(in) :: algo

    integer :: nloc, nbands, ierr

    PUSH_SUB(deallocate_accel_mtxel_sig)

    nbands = accel_mtxel_sig%mtxel_band_block_size

    do nloc = 1, nbands
      select case(algo)
        case (OPENACC_ALGO)
#ifdef OPENACC
          !$acc exit data delete ( accel_mtxel_sig%bands_vec(nloc)%vec, accel_mtxel_sig%mtxel_vec(nloc)%vec )
          !$acc exit data delete ( accel_mtxel_sig%mtxel_box(nloc)%box )
#else
        call die_algos("OpenACC")
#endif
        case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
          !$omp target exit data map (delete: accel_mtxel_sig%bands_vec(nloc)%vec, accel_mtxel_sig%mtxel_vec(nloc)%vec )
          !$omp target exit data map (delete: accel_mtxel_sig%mtxel_box(nloc)%box )
#else
        call die_algos("OpenMP Target")
#endif
        case default
          call die("Invald algorithm for deallocate_accel_mtxel_sig", only_root_writes = .true.)
      end select
      SAFE_DEALLOCATE ( accel_mtxel_sig%bands_vec(nloc)%vec )
      SAFE_DEALLOCATE ( accel_mtxel_sig%mtxel_vec(nloc)%vec )
      SAFE_DEALLOCATE ( accel_mtxel_sig%mtxel_box(nloc)%box )
      call accel_destroy_fft_plan(accel_mtxel_sig%mtxel_plan(nloc))
    end do

    call accel_destroy_fft_plan(accel_mtxel_sig%eqp_plan)

    select case(algo)
      case (OPENACC_ALGO)
#ifdef OPENACC
        !$acc exit data delete( accel_mtxel_sig%gvec_comp )
        !$acc exit data delete( accel_mtxel_sig%gvec_isrtx_k, accel_mtxel_sig%gvec_isrtx_q, accel_mtxel_sig%gvec_isrtx_mtxel )
        !$acc exit data delete( accel_mtxel_sig%bands_vec, accel_mtxel_sig%mtxel_box, accel_mtxel_sig%mtxel_vec )
        !$acc exit data delete( accel_mtxel_sig%eqp_vec, accel_mtxel_sig%eqp_box )
        !$acc exit data delete( accel_mtxel_sig )
#else
        call die_algos("OpenACC")
#endif
      case (OMP_TARGET_ALGO)
#ifdef OMP_TARGET
        !$omp target exit data map (delete: accel_mtxel_sig%gvec_comp )
        !$omp target exit data map (delete: accel_mtxel_sig%gvec_isrtx_k, accel_mtxel_sig%gvec_isrtx_q, accel_mtxel_sig%gvec_isrtx_mtxel )
        !$omp target exit data map (delete: accel_mtxel_sig%bands_vec, accel_mtxel_sig%mtxel_box, accel_mtxel_sig%mtxel_vec )
        !$omp target exit data map (delete: accel_mtxel_sig%eqp_vec, accel_mtxel_sig%eqp_box )
        !$omp target exit data map (delete: accel_mtxel_sig )
#else
        call die_algos("OpenMP Target")
#endif
      case default
        call die("Invald algorithm for deallocate_accel_mtxel_sig", only_root_writes = .true.)
    end select

    SAFE_DEALLOCATE ( accel_mtxel_sig%gvec_comp )
    SAFE_DEALLOCATE ( accel_mtxel_sig%gvec_isrtx_k )
    SAFE_DEALLOCATE ( accel_mtxel_sig%gvec_isrtx_q )
    SAFE_DEALLOCATE ( accel_mtxel_sig%gvec_isrtx_mtxel )

    SAFE_DEALLOCATE ( accel_mtxel_sig%bands_vec )
    SAFE_DEALLOCATE ( accel_mtxel_sig%mtxel_vec )
    SAFE_DEALLOCATE ( accel_mtxel_sig%mtxel_box )
    SAFE_DEALLOCATE ( accel_mtxel_sig%mtxel_plan )

    SAFE_DEALLOCATE ( accel_mtxel_sig%eqp_vec )
    SAFE_DEALLOCATE ( accel_mtxel_sig%eqp_box )

    POP_SUB(deallocate_accel_mtxel_sig)
  end subroutine deallocate_accel_mtxel_sig

end module accel_fft_m
