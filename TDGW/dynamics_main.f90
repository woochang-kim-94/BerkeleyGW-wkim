!===============================================================================
!
! Routines:
!
! (1) dynamics            Originally: Yang-Hao Chan, td_cohsex code (python/fortran)
!                                     (2016-2020)
!                         BerkeleyGW adaption: Zhenglu Li (2021)
!
! This is the main routine for the time-dependent GW (TDGW) code
!
! Original version by Yang-Hao Chan (YHC)
! -> python td_cohsex code for TD-aGW (Q=0) and a fortran version
! -> gauge rotation for intraband dipole in nonlinear optics
! -> ref: Chan, Qiu, da Jornada, Louie, PNAS 22, e1906938118 (2021)
!
! Restructured and extended by Zhenglu Li (ZL):
! -> adapted into BerkeleyGW package coding scheme
! -> reorganized file and data structures
! -> generalized version for finite Q
!
! TODO list:
! -> include coherent phonon dynamics (ZL)
! -> include electron-phonon scattering (ZL)
! -> extend insulating systems to include metals: best ways is probably to 
!       get rid of the prep_kernel.py step and directly read bsemat.h5,
!       and rearrange the head, wing, body, for semiconductor, metal, and graphene
!
!===============================================================================

#include "f_defs.h"

program dynamics

  use global_m
  use write_program_header_m
  use hdf5
  use io_utils_m
  use timing_m, only: common_timing, timing => tdgw_timing
  use input_dynamics_m
  use equation_of_motion_m
  use propagation_m
  use propagation_q_m
  use measurement_m

  implicit none

!----------------------
! declare variables

  type(progress_info) :: prog_info  ! progress report
  type(dynaminfo) :: dyn            ! control parameters
  type(td_equation), dimension(:), allocatable :: tdevols   ! array for q of tdevol
  type(measurement), dimension(:), allocatable :: measures  ! array for q of measure

  integer :: iq                 ! qpoint label
  integer :: it                 ! time step label

!---------------- Begin Program ----------------

  call timing%init()
  call common_timing%init()
  call timing%start(timing%total)

  call peinfo_init()

  call write_program_header('TDGW-Dynamics', .false.)

  call read_input(dyn)

  ! TODO: distribute tdevols and measures
  SAFE_ALLOCATE(tdevols, (dyn%n_finite_q))
  SAFE_ALLOCATE(measures, (dyn%n_finite_q))

!=========================================================================
!> initiliazation

  ! TODO: to be parallelized over q
  do iq = 1, dyn%n_finite_q
    ! read interaction matrix elements
    call load_bse_kernel(dyn, iq, tdevols(iq))

    ! read quasiparticle energies and optical dipole matrxi elements
    call load_eqp_vmtxel(dyn, tdevols(iq))

    ! initialize green`s function
    call init_green_function(tdevols(iq))

    ! initialize hamiltonian
    call init_hamiltonian(tdevols(iq))

    ! solve bse hamiltonian
    call solve_bse_hamiltonian(tdevols(iq))

    ! initialize measurement, which will also take measurement at it = 1
    call init_measurement(dyn, tdevols(iq), measures(iq))
  enddo

!=========================================================================
!> main loop over time steps

  call progress_init(prog_info, 'time propagation', 'time step', dyn%nt)
  call progress_step(prog_info)  ! this one step is for the above initialization process

  ! it = 1 records equilibrium and represents time = -inf, and has been done in init_measurement
  ! the real time dynamics starts, i.e. efield starts to apply, when it = 2
  do it = 2, dyn%nt
    call progress_step(prog_info)

    ! propagate to next time step
    call propagator_q(dyn, tdevols)

    do iq = 1, dyn%n_finite_q
      call take_measurement(dyn, tdevols(iq), measures(iq))
    enddo
  enddo ! it

  ! all nodes save all q
  if (peinf%inode == 0) then
    do iq = 1, dyn%n_finite_q
      call write_measurement(dyn, measures(iq))
    enddo
  endif

  call progress_free(prog_info)

  call timing%stop(timing%total)
  call timing%print(common_timing)

end program dynamics
