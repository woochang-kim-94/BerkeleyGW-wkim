!===============================================================================
! 
! Routines:
!     Time propagators to give next iteration
!
! (1) propagator()
!     propagator interface
!
! (2) euler_propagator()
!     simplest first-order procedure
!
! (3) runge_kutta_propagaor()
!     4th order runge-kutta algorithm
!
!===============================================================================

#include "f_defs.h"

module propagation_m

  use global_m
  use timing_m, only: timing => tdgw_timing
  use equation_of_motion_m
  use utility_m
  use measurement_m

  implicit none

  private

  public :: &
    propagator, &
    euler_propagator

  contains

!=========================================================================
!> propagator interface
subroutine propagator(dyn, tdevol)

  type(dynaminfo), intent(in) :: dyn
  type(td_equation), intent(inout) :: tdevol

  PUSH_SUB(propagator)

  call timing%start(timing%propagator)

  if (dyn%prop_algo == 1) then
    call euler_propagator(dyn, tdevol)
  elseif (dyn%prop_algo == 2) then
    call runge_kutta_propagator(dyn, tdevol)
  else
    call die('Invalid prop_algo input flag.')
  endif

  call timing%stop(timing%propagator)

  POP_SUB(propagator)

end subroutine propagator

!=========================================================================
!> euler propagator (equivalent to 1st order runge-kutta)
subroutine euler_propagator(dyn, tdevol, commute_save)

  type(dynaminfo), intent(in) :: dyn
  type(td_equation), intent(inout) :: tdevol
  complex(DPC), intent(inout), optional :: commute_save(:,:,:)  ! save commutators
  !
  complex(DPC), allocatable :: h_int(:,:)       ! interacting hamiltonian including equilibrium h, field, and interaction
  complex(DPC), allocatable :: commutator(:,:)  ! commutator = i*dG/dt for one kpoint
  integer :: ik
  complex(DPC) :: imag_one
  real(DP) :: efield

  PUSH_SUB(euler_propagator)

  imag_one = (0.0d0,1.0d0)

  ! H(t) = H(U(t),G(t))
  ! time dependence in H is in external varying field and green`s function
  efield = dyn%eamp * gaussian(tdevol%itlabel * dyn%dt, dyn%tpulse, dyn%ewidth)
  tdevol%efield = efield

  SAFE_ALLOCATE(h_int, (tdevol%nb, tdevol%nb))
  SAFE_ALLOCATE(commutator, (tdevol%nb, tdevol%nb))

  do ik = 1, tdevol%nk
    ! prepare h_int
    h_int(:,:) = tdevol%hamiltonian(:,:,ik) - efield * tdevol%dipole(:,:,ik) + tdevol%interaction(:,:,ik)

    ! compute commutator
    call commutation(tdevol%nb, h_int, tdevol%gf(:,:,ik), commutator)

    ! update green`s function of ik 
    tdevol%gf(:,:,ik) = tdevol%gf(:,:,ik) + commutator(:,:) * dyn%dt / imag_one

    ! if request saving, collect commutators here
    if (present(commute_save)) then
      commute_save(:,:,ik) = commutator(:,:)
    endif
  enddo

  SAFE_DEALLOCATE(h_int)
  SAFE_DEALLOCATE(commutator)

  ! update interaction for all kpts
  call get_interaction(tdevol)

  ! update time step label
  tdevol%itlabel = tdevol%itlabel + 1.0d0

  !write(*,*) 'finishing one iteration of euler propagation'

  POP_SUB(euler_propagator)

end subroutine euler_propagator

!=========================================================================
!> runge-kutta propagator (4th order implementation), default algorithm
subroutine runge_kutta_propagator(dyn, tdevol)

  type(dynaminfo), intent(in) :: dyn
  type(td_equation), intent(inout) :: tdevol
  !
  type(td_equation) :: tdev
  complex(DPC), allocatable :: commu1(:,:,:), commu2(:,:,:), commu3(:,:,:), commu4(:,:,:)
  complex(DPC) :: imag_one
  real(DP) :: efield
  integer :: it

  PUSH_SUB(runge_kutta_propagator)

  imag_one = (0.0d0,1.0d0)

  ! Runge-Kutta algorithm
  ! in our data structure and euler_propagator implementaiton, we finish updating interaction and gf on exit
  ! here we map one-to-one correspondence with Runge-Kutta (4th order) method
  !   ref: https://en.wikipedia.org/wiki/Runge-Kutta_methods
  !
  !           |       Runge-Kutta notation         |         TD-aGW notation
  !--------------------------------------------------------------------------------
  !   e.o.m.  |        dy/dt = f(t, y(t))          |    dG/dt = -i [H(G(t)), G(t)]
  ! variables |             y(t)                   |              G(t)
  !           |           f(t, y(t))               |      -i [H(G(t)), G(t)]
  !           |              h                     |               dt

  ! prepare four commutators corresponding to i*k1, i*k2, i*k3, i*k4
  SAFE_ALLOCATE(commu1, (tdevol%nb, tdevol%nb, tdevol%nk))
  SAFE_ALLOCATE(commu2, (tdevol%nb, tdevol%nb, tdevol%nk))
  SAFE_ALLOCATE(commu3, (tdevol%nb, tdevol%nb, tdevol%nk))
  SAFE_ALLOCATE(commu4, (tdevol%nb, tdevol%nb, tdevol%nk))

  commu1(:,:,:) = (0.0d0,0.0d0)
  commu2(:,:,:) = (0.0d0,0.0d0)
  commu3(:,:,:) = (0.0d0,0.0d0)
  commu4(:,:,:) = (0.0d0,0.0d0)

  ! ==================== step 1 ===================
  !> k1 = f(t, y(t))        <=>    k1 = dG/dt|_t = -i [H(G(t)), G(t)]
  !
  call init_td_equation_from_copy(tdevol, tdev)

  ! compute commu1
  call euler_propagator(dyn, tdev, commute_save=commu1)

  ! save this efield value
  efield = tdev%efield

  call free_td_equation(tdev)

  ! ==================== step 2 ====================
  ! k2 = f(t+0.5, y1(t+0.5))  <=>  k2 = dG`/dt|_{t+0.5} = -i [H`(G`(t+0.5)), G`(t+0.5)]
  ! y1 = G`(t+0.5) = y(t) + h * k1 / 2.0
  ! here G` denotes that this is a euler extrapolated quantity with specific input
  !
  call init_td_equation_from_copy(tdevol, tdev)

  ! U(t+0.5)
  tdev%itlabel = tdev%itlabel + 0.5d0
  ! G`(t+0.5) using commu1
  tdev%gf(:,:,:) = tdev%gf(:,:,:) + commu1(:,:,:) * dyn%dt / imag_one / 2.0d0
  ! H`(G`(t+0.5))
  call get_interaction(tdev)

  ! compute commu2
  call euler_propagator(dyn, tdev, commute_save=commu2)

  call free_td_equation(tdev)

  ! ==================== step 3 ====================
  ! k3 = f(t+0.5, y2(t+0.5))  <=>  k3 = dG``/dt|_{t+0.5} = -i [H``(G``(t+0.5)), G``(t+0.5)]
  ! y2 = G``(t+0.5) = y(t) + h * k2 / 2.0
  ! note G`` is differently constructed as G`, even though they are both at the same time snapshot
  !     
  call init_td_equation_from_copy(tdevol, tdev)

  ! U(t+0.5)
  tdev%itlabel = tdev%itlabel + 0.5d0
  ! G``(t+0.5) using commu2
  tdev%gf(:,:,:) = tdev%gf(:,:,:) + commu2(:,:,:) * dyn%dt / imag_one / 2.0d0
  ! H``(G``(t+0.5))
  call get_interaction(tdev)

  ! compute commu3
  call euler_propagator(dyn, tdev, commute_save=commu3)

  call free_td_equation(tdev)

  ! ==================== step 4 ====================
  ! k4 = f(t+1, y3(t+1))  <=>  k4 = dG/dt|_{t+1} = -i [H(G(t+1)), G(t+1)]
  ! y3 = G(t+1) = y(t) + h * k3
  !
  call init_td_equation_from_copy(tdevol, tdev)

  ! U(t+1)
  tdev%itlabel = tdev%itlabel + 1.0d0
  ! G(t+1) using commu3
  tdev%gf(:,:,:) = tdev%gf(:,:,:) + commu3(:,:,:) * dyn%dt / imag_one
  ! H(G(t+1))
  call get_interaction(tdev)

  ! commute commu4
  call euler_propagator(dyn, tdev, commute_save=commu4)

  call free_td_equation(tdev)

  ! ==================== update ====================
  ! Runge-Kutta (4th-order) update
  !
  ! y(t+1) = y(t) + (1/6) * h * (k1 + 2*k2 + 2*k3 + k4)
  !
  tdevol%efield = efield

  ! update green`s function
  tdevol%gf(:,:,:) = tdevol%gf(:,:,:) + (1.0d0 / 6.0d0) * (dyn%dt / imag_one) & 
                   & * (commu1(:,:,:) + 2.0d0 * commu2(:,:,:) + 2.0d0 * commu3(:,:,:) + commu4(:,:,:))

  ! get interaction
  call get_interaction(tdevol)

  ! update time step label 
  tdevol%itlabel = tdevol%itlabel + 1.0d0

  ! deallocate
  SAFE_DEALLOCATE(commu1)
  SAFE_DEALLOCATE(commu2)
  SAFE_DEALLOCATE(commu3)
  SAFE_DEALLOCATE(commu4)

  POP_SUB(runge_kutta_propagator)

end subroutine runge_kutta_propagator

end module propagation_m
