!===============================================================================
! 
! Types:
! (1) measurement
!
! Routines:
! (1) init_measurement()
!     initialize measurement object
! 
! (2) take_measurement()
!     take measurement of a given td_equation object
!
!===============================================================================

#include "f_defs.h"

module measurement_m

  use global_m
  use timing_m, only: timing => tdgw_timing
  use equation_of_motion_m

  implicit none
  
  private

  public :: &
    measurement, &
    init_measurement, &
    take_measurement, &
    write_measurement

  type measurement
    integer :: iq                                     ! i-th q point label
    integer :: nt                                     ! number of time steps
    real(DP), allocatable :: t_series(:)              ! time series (nt)
    real(DP), allocatable :: efield(:)                ! efield as function of time (nt)
    complex(DPC), allocatable :: occupation(:,:,:)    ! occupation (nb, nk, nt)
    complex(DPC), allocatable :: polarization(:)      ! polarization (nt)

  end type measurement

contains

!=========================================================================
!> initialize measurement object
subroutine init_measurement(dyn, tdevol, measure)

  type(dynaminfo), intent(in) :: dyn
  type(td_equation), intent(in) :: tdevol
  type(measurement), intent(inout) :: measure

  PUSH_SUB(init_measurement)

  ! allocate arrays
  SAFE_ALLOCATE(measure%t_series, (dyn%nt))
  SAFE_ALLOCATE(measure%efield, (dyn%nt))
  SAFE_ALLOCATE(measure%occupation, (tdevol%nb, tdevol%nk, dyn%nt))
  SAFE_ALLOCATE(measure%polarization, (dyn%nt))

  measure%iq = tdevol%iq   ! initialize iq label

  ! take measurement for equilibrium tdevol
  if (tdevol%itlabel .ne. 1) then
    call die('Measurement must be initialized with equilibrium tdevol (it = 0).')
  endif

  call take_measurement(dyn, tdevol, measure)

  POP_SUB(init_measurement)

end subroutine init_measurement

!=========================================================================
!> take measurement of given td_equation object
subroutine take_measurement(dyn, tdevol, measure)

  type(dynaminfo), intent(in) :: dyn
  type(td_equation), intent(in) :: tdevol
  type(measurement), intent(inout) :: measure

  integer :: ib, jb, ik, it
  complex :: pol
  complex :: imag_one

  PUSH_SUB(take_measurement)

  call timing%start(timing%measurement)

  imag_one = (0.0d0, 1.0d0)

  ! check iq consistency
  if (measure%iq .ne. tdevol%iq) then
    call die('Mismatch iq in measure and tdevol.')
  endif

  ! convert itlabel to integer
  it = int(tdevol%itlabel)

  if (abs(tdevol%itlabel - it) .gt. TOL_Small) then
    call die('Measurement cannot be taken for non-integer itlabel value')
  endif

  ! record t_series
  measure%t_series(it) = (it - 1) * dyn%dt  ! time axis origin is zero

  ! measure efield
  measure%efield(it) = tdevol%efield

  ! measure occupation
  do ik = 1, tdevol%nk
    do ib = 1, tdevol%nb
      measure%occupation(ib, ik, it) = -imag_one * tdevol%gf(ib, ib, ik)
    enddo ! ib
  enddo ! ik

  ! measure polarization
  ! P = -i \hbar \sum_{nm,k} r_{nm,k} G_{mn,k}
  pol = (0.0d0, 0.0d0)

  do ik = 1, tdevol%nk
    do ib = 1, tdevol%nb
      do jb = 1, tdevol%nb
        pol = pol - imag_one * tdevol%dipole(ib, jb, ik) * tdevol%gf(jb, ib, ik)
      enddo ! jb
    enddo ! ib
  enddo ! ik

  measure%polarization(it) = pol

  call timing%stop(timing%measurement)

  POP_SUB(take_measurement)

end subroutine take_measurement

!=========================================================================
!> write measurement of given td_equation object
subroutine write_measurement(dyn, measure)

  type(dynaminfo), intent(in) :: dyn
  type(measurement), intent(in) :: measure

  character*256 :: pol_fname
  character*256 :: iq_str
  real(DP) :: time
  integer :: it

  PUSH_SUB(write_measurement)

  call timing%start(timing%output)

  write(iq_str,*) measure%iq
  pol_fname = 'bse_'//trim(adjustl(iq_str))//'/polarization'//'.dat'

  call open_file(21, file=pol_fname, form='formatted', status='replace')

  do it = 1, dyn%nt
    time = it * dyn%dt  ! consistent with efield definition in proparation_q.f90
                        ! guaranteed integer as protected by take_measurement converting itlabel to it

    write(21, '(f20.3, 2f20.15)') time, real(measure%polarization(it)), aimag(measure%polarization(it))
  enddo

  call close_file(21)

  call timing%stop(timing%output)

  POP_SUB(write_measurement)

end subroutine write_measurement

end module measurement_m
