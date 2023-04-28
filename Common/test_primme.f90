#include "f_defs.h"

program test_primme

  use global_m
  use bgw_mpi_m
  use primme_interface_m
  use write_program_header_m
  implicit none

  logical, parameter :: &
    dump_matrix = .false., &
    test_lapack = .false., &
    test_gpu = .true., &
    test_nogpu = .true.

  integer :: N, Nev, Nb, Nown, ii, jj, jj_l, ipe
  integer :: lwork, lrwork, info
  character(len=16) :: Nstr

  SCALAR, allocatable :: work(:), evecs_lapack(:,:)
  real(DP), allocatable :: rwork(:), evals_lapack(:)

  real(DP), allocatable :: evals_primme(:)
  SCALAR, allocatable :: evecs_primme(:,:)
  SCALAR, allocatable :: matrix_loc(:,:)

  call peinfo_init()
  if (peinf%inode==0) call write_program_header('test_primme', .false.)
  flush(6)

  N = 100
  if (command_argument_count()>0) then
    call get_command_argument(1, Nstr)
    read(Nstr,*) N
  endif

  if (peinf%inode==0) print *, 'N =', N
  FLUSH(6)

  Nb = DIVUP(N, peinf%npes)
  Nown = NUMROC(N, Nb, peinf%inode, 0, peinf%npes)
  Nev = 10

  SAFE_ALLOCATE(matrix_loc, (N,max(1,Nown)))

  !omega = exp(-2d0 * PI_D * (0d0,1d0) / N)
  matrix_loc(:,:) = ZERO
  call bgw_barrier()
  do jj_l = 1, Nown
    jj = indxl2g(jj_l, Nb, peinf%inode, 0, peinf%npes)
    !do ii = 1, N
    !  matrix_loc(ii,jj_l) = 1d0 / (abs(ii - jj) + 1d0)
    !  !matrix_loc(ii,jj_l) = omega**((ii-1)*(jj-1))
    !  !if (ii==jj) matrix_loc(ii,jj_l) = ii
    !enddo
    matrix_loc(jj,jj_l) = 100d0/(jj+1)
  enddo
  !matrix_loc = matrix_loc / sqrt(dble(N))
  call bgw_barrier()

  if (peinf%inode==0) then
    write(6,*) 'N =', N
    write(6,*) 'Nb =', Nb
    write(6,*) 'Nown =', Nown
  endif
  FLUSH(6)

  if (dump_matrix) then
    if (peinf%inode==0) then
      write(6,*) 'matrix_loc.T:'
      FLUSH(6)
    endif
    do ipe = 1, peinf%npes
      call bgw_barrier()
      if (peinf%inode==ipe-1) then
        FLUSH(6)
        print *
        print *, 'inode =', peinf%inode
        do jj_l = 1, Nown
#ifdef CPLX
          write(6,'(*(f5.2,",",f5.2,:,1x))') matrix_loc(:,jj_l)
#else
          write(6,'(*(f5.2,:,1x))') matrix_loc(:,jj_l)
#endif
        enddo
        FLUSH(6)
      endif
      call bgw_barrier()
    enddo
  endif

  if (test_lapack) then
    if (peinf%inode==0) then
      write(*,*)
      write(*,*) '*****************************************************************'
      write(*,*) 'LAPACK SOLVER'
      write(*,*) '*****************************************************************'
      FLUSH(6)
    endif

    if (peinf%inode==0 .and. peinf%npes==1) then
      lwork = 3*N - 1
      lrwork = 3*N - 2

      SAFE_ALLOCATE(work, (lwork))
      SAFE_ALLOCATE(rwork, (lrwork))
      SAFE_ALLOCATE(evecs_lapack, (N,N))
      SAFE_ALLOCATE(evals_lapack, (N))

      write(6,*) 'Eigenvalues (LAPACK):'
      FLUSH(6)
      evecs_lapack = matrix_loc
  #ifdef CPLX
      call zheev('N', 'U', N, evecs_lapack, N, evals_lapack, &
        work, lwork, rwork, info)
  #else
      call dsyev('N', 'U', N, evecs_lapack, N, evals_lapack, &
        work, lwork, info)
  #endif
      write(6,*) evals_lapack(1:nev)
      FLUSH(6)

      SAFE_DEALLOCATE(evals_lapack)
      SAFE_DEALLOCATE(evecs_lapack)
      SAFE_DEALLOCATE(rwork)
      SAFE_DEALLOCATE(work)
    endif
  endif

  do ii = 1, 2
    if (ii==1 .and. .not. test_nogpu) cycle
    if (ii==2 .and. .not. test_gpu) cycle

    if (peinf%inode==0) then
      write(*,*)
      write(*,*) '*****************************************************************'
      if (ii==1) then
        write(*,*) 'PRIMME SOLVER (no GPU)'
      else
        write(*,*) 'PRIMME SOLVER (GPU)'
      endif
      write(*,*) '*****************************************************************'
      FLUSH(6)
    endif

    SAFE_ALLOCATE(evals_primme, (Nev))
    SAFE_ALLOCATE(evecs_primme, (max(1,Nown),Nev))

    call solve_primme_explicit(N, Nown, nev, evals_primme, evecs_primme, &
      matrix_loc, use_gpu=(ii==2), verbose=.true.)
    FLUSH(6)

    if (peinf%inode==0) then
      write(6,*) 'Eigenvalues (PRIMME):'
      write(6,*) evals_primme
      FLUSH(6)
    endif

    SAFE_DEALLOCATE(evals_primme)
    SAFE_DEALLOCATE(evecs_primme)
  enddo

  SAFE_DEALLOCATE(matrix_loc)

  call MPI_Finalize(mpierr)

end program test_primme
