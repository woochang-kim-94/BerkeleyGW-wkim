!===============================================================================
!
! Routines:
!   Contains utilities subroutines/functions
!
! (1) kindex()
!     find a given kpoint in the list, return index
!
!===============================================================================

#include "f_defs.h"

module utility_m

  use global_m
  use timing_m, only: timing => tdgw_timing

  implicit none

  private

  public :: &
    kindex, &
    commutation, &
    generalized_commutation, &
    gaussian, &
    mpi_dispatcher, &
    mpi_dispatcher_one_dim

contains

!=========================================================================
!> given a kpoint coordinate and a list of k, return its index in the list
!> if two kpoints are differed by an umklapp G vector, they are considered the same kpoint
!> fractional coordinates
!  taken from ../Sigma/sympert_utils.f90
integer function kindex(kpt, klist, nk, skip_check)
  real(DP), intent(in) :: kpt(3)
  real(DP), allocatable, intent(in) :: klist(:,:)
  integer, intent(in) :: nk
  logical, intent(in), optional :: skip_check
  !
  integer :: ik
  real(DP) :: kdiff(3), knorm
  logical :: check

  ! no PUSH/POP because called too frequently

  check = .true.
  if (present(skip_check)) then
    if (skip_check) then
      check = .false.
    endif
  endif

  kindex = 0  ! initialize with an invalid value

  do ik = 1, nk
    kdiff(:) = kpt(:) - klist(:, ik)
    kdiff(:) = kdiff(:) - nint(kdiff(:))
    knorm = sqrt(dot_product(kdiff, kdiff))

    if (knorm .lt. TOL_Small) then
      kindex = ik
      exit
    endif
  enddo

  if (check) then
    if (kindex .eq. 0) then
      call die('kindex: kpt not found in provided klist.')
    endif
  endif

  return

end function kindex

!=========================================================================
!> get commutator of matrix_A and matrix_B, commutator = [matA, matB] = matA.matB - matB.matA
subroutine commutation(dimens, matA, matB, commutator)

  integer, intent(in) :: dimens
  complex(DPC), intent(in) :: matA(:,:)
  complex(DPC), intent(in) :: matB(:,:)
  complex(DPC), intent(inout) :: commutator(:,:)
  !
  integer :: i, j, k

  PUSH_SUB(commutation)

  call timing%start(timing%commutation)

  commutator(:,:) = (0.0d0,0.0d0)

  do i = 1, dimens
    do j = 1, dimens
      do k = 1, dimens
        commutator(i, j) = commutator(i, j) & 
                          + matA(i, k) * matB(k, j) - matB(i, k) * matA(k, j)
      enddo ! k
    enddo ! j
  enddo ! i

  call timing%stop(timing%commutation)

  POP_SUB(commutation)

end subroutine

!=========================================================================
!> define generalized commutator of matrix_A, matrix_B, matrix_C and matrix_D
!    commutator = matA.matB - matD.matC
subroutine generalized_commutation(dimens, matA, matB, matC, matD, commutator)

  integer, intent(in) :: dimens
  complex(DPC), intent(in) :: matA(:,:)
  complex(DPC), intent(in) :: matB(:,:)
  complex(DPC), intent(in) :: matC(:,:)
  complex(DPC), intent(in) :: matD(:,:)
  complex(DPC), intent(inout) :: commutator(:,:)
  !
  integer :: i, j, k

  PUSH_SUB(commutation)

  call timing%start(timing%commutation)

  commutator(:,:) = (0.0d0,0.0d0)

  do i = 1, dimens
    do j = 1, dimens
      do k = 1, dimens
        commutator(i, j) = commutator(i, j) & 
                          + matA(i, k) * matB(k, j) - matD(i, k) * matC(k, j)
      enddo ! k
    enddo ! j
  enddo ! i

  call timing%stop(timing%commutation)

  POP_SUB(commutation)

end subroutine


!=========================================================================
!> gaussian function normalized to one
real(DP) function gaussian(x, mu, sigma)
  real(DP), intent(in) :: x
  real(DP), intent(in) :: mu
  real(DP), intent(in) :: sigma

  ! no PUSH/POP because called too frequently

  gaussian = exp(-0.5d0 * (x - mu) * (x - mu) / sigma / sigma) / sigma / sqrt(2.0d0 * PI_D)

  return

end function gaussian

!=========================================================================
!> mpi dispatcher for (k, k`) or (k, q`) pairs
!  the outer most q is not supposed to be parallelized
!  p below labels prime points
subroutine mpi_dispatcher(nk, np, kstart, kend, pstart, pend)
  integer, intent(in)  :: nk
  integer, intent(in)  :: np
  integer, intent(out) :: kstart
  integer, intent(out) :: kend
  integer, intent(out) :: pstart
  integer, intent(out) :: pend
  !
  integer :: std_k_npes     ! number of processors per k for two-level parallelization
  integer :: my_k_npes      ! number of processors per k being used
  integer :: inode_my_k     ! local node index
  integer :: std_np, my_np  ! distributed np of each processor

  PUSH_SUB(mpi_dispatcher)

  ! initialize output with invalid indices since some cores will be idle
  kstart = 0
  kend = -1
  pstart = 0
  pend = -1

  ! choose one-level vs. two-level parallelization according to total number of cores
  if (peinf%npes <= nk) then
    ! one-level parallelization over nk only
    if (peinf%inode == 0) write(*,*) 'One-level parallelization.'

    ! distribute nk over global npes
    call mpi_dispatcher_one_dim(peinf%inode, peinf%npes, nk, kstart, kend)

    ! no distribution for np points
    pstart = 1
    pend = np

  else
    ! two-level parallelization over (nk, np)
    if (peinf%inode == 0) write(*,*) 'Two-level parallelization.'

    ! here there is no need to deal with remaining cores
    ! since all k have the same np points and the computing time is determined by the slowest ik
    ! for non-perfect division, uneven distribution of ncores over nk will not speed up the calculation

    if (mod(peinf%npes, nk) .ne. 0) then
      if(peinf%inode == 0) then
        write(*,'(1x,a)')     'WARNING: Non-ideal distribuiton of all computing cores in two-level parallelization.'
        write(*,'(1x,a,i10)') '         Number of idle cores:', mod(peinf%npes, nk)
        write(*,'(1x,a)')     'Use integer multiples of numbers of kpoints as number of cores for optimal distribution of 1st level.'
      endif
    endif

    std_k_npes = peinf%npes / nk   ! integer division

    ! only work within the subset of cores within nk usage
    if (peinf%inode < std_k_npes * nk) then
      ! determine for each core, which ik it belongs to
      ! note: kstart = ik, kend = ik, for the loop
      kstart = peinf%inode / std_k_npes + 1
      kend = kstart
  
      ! now prepare to distribute np points for each ik
      if (std_k_npes > np) then
        ! if for each ik, we have more cores than np
        ! we only use the first np cores
        my_k_npes = np
  
        if (peinf%inode == 0) write(*,*) 'WARNING: 2nd-level in parallelization is not using all allocated cores.'
      else
        my_k_npes = std_k_npes
      endif

      ! only work within the subset of cores within my_k_npes
      inode_my_k = mod(peinf%inode, std_k_npes)  ! inode counts from zero

      if (inode_my_k < my_k_npes) then
        call mpi_dispatcher_one_dim(inode_my_k, my_k_npes, np, pstart, pend)
      endif
    endif

  endif ! choice between one-level vs. two-level

  POP_SUB(mpi_dispatcher)

end subroutine mpi_dispatcher


!=========================================================================
!> mpi dispatcher for one dimension distribution
subroutine mpi_dispatcher_one_dim(inode, npes, nk, kstart, kend)
  integer, intent(in)  :: inode
  integer, intent(in)  :: npes
  integer, intent(in)  :: nk
  integer, intent(out) :: kstart
  integer, intent(out) :: kend
  !
  integer :: std_nk
  integer, allocatable :: count_my_nk(:)  ! distributed nk of each processor
  integer :: ipe, ik

  PUSH_SUB(mpi_dispatcher_one_dim)

  if (npes > nk) then
    call die('Too many cores to be distributed.')  ! this should never happen... controled from outside call
  endif

  SAFE_ALLOCATE(count_my_nk, (npes))

  count_my_nk(:) = 0

  ! this is just used to count, not the actual kpoints assignment to nodes
  do ik = 1, nk
    ipe = mod(ik-1, npes) + 1
    count_my_nk(ipe) = count_my_nk(ipe) + 1
  enddo

  kstart = sum(count_my_nk(:inode)) + 1  ! Fortran counting convention in do loop
  kend = sum(count_my_nk(:inode+1))

  SAFE_DEALLOCATE(count_my_nk)

  POP_SUB(mpi_dispatcher_one_dim)

end subroutine mpi_dispatcher_one_dim

end module utility_m
