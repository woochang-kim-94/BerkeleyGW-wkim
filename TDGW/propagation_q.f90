!===============================================================================
! 
! Routines:
!     Time propagators to give next iteration
!     All q`s co-evolve
!
! (1) propagator_q()
!     propagator interface
!
! (2) euler_propagator_q()
!     simplest first-order procedure
!
! (3) runge_kutta_propagaor_q()
!     4th order runge-kutta algorithm
!
!===============================================================================

#include "f_defs.h"

module propagation_q_m

  use global_m
  use timing_m, only: timing => tdgw_timing
  use equation_of_motion_m
  use utility_m
  use measurement_m

  implicit none

  private

  public :: &
    propagator_q, &
    euler_propagator_q

  contains

!=========================================================================
!> propagator interface
subroutine propagator_q(dyn, tdevols)

  type(dynaminfo), intent(in) :: dyn
  type(td_equation), intent(inout) :: tdevols(:)

  PUSH_SUB(propagator)

  call timing%start(timing%propagator)

  if (dyn%prop_algo == 1) then
    call euler_propagator_q(dyn, tdevols)
  elseif (dyn%prop_algo == 2) then
    call runge_kutta_propagator_q(dyn, tdevols)
  else
    call die('Invalid prop_algo input flag.')
  endif

  call timing%stop(timing%propagator)

  POP_SUB(propagator)

end subroutine propagator_q

!=========================================================================
!> euler propagator (equivalent to 1st order runge-kutta)
!  parallelization is done withing this subroutine
subroutine euler_propagator_q(dyn, tdevols, commute_save)

  type(dynaminfo), intent(in) :: dyn
  type(td_equation), intent(inout) :: tdevols(:)
  complex(DPC), intent(inout), optional :: commute_save(:,:,:,:)  ! save commutators (nb, nb, nk, nq)
  !
  complex(DPC), allocatable :: h_int_A(:,:), gf_B(:,:), h_int_C(:,:), gf_D(:,:) 
                             ! interacting hamiltonian including equilibrium h, field, and interaction
                             ! matA, matB, matC, matD for generalized_commutation 
  complex(DPC), allocatable :: commutator(:,:), commutators(:,:,:,:)  ! commutator = i*dG/dt for all kpoints and qpoints
  integer :: ncount
  integer :: ik, ikpqp, iq, iqp, iqmqp
  real(DP) :: qmqp(3)  ! Q-Q` point
  real(DP) :: kpqp(3)  ! k+Q` point
  complex(DPC) :: imag_one
  real(DP) :: efield

  PUSH_SUB(euler_propagator_q)

  imag_one = (0.0d0,1.0d0)

  ! H(t) = H(U(t),G(t))
  ! time dependence in H is in external varying field and green`s function
  ! only for iq_efield tdevol
  do iq = 1, dyn%n_finite_q
    ! hermitian light field made of two conjugate components
    if ((iq == dyn%iq_efield(1)) .or. (iq == dyn%iq_efield(2))) then
      efield = dyn%eamp * gaussian(tdevols(iq)%itlabel * dyn%dt, dyn%tpulse, dyn%ewidth)
      tdevols(iq)%efield = efield  ! other iq (other than efield_iq) have zero efield all the time
    else
      tdevols(iq)%efield = 0.0d0
    endif
  enddo

  SAFE_ALLOCATE(h_int_A, (tdevols(1)%nb, tdevols(1)%nb))
  SAFE_ALLOCATE(gf_B,    (tdevols(1)%nb, tdevols(1)%nb))
  SAFE_ALLOCATE(h_int_C, (tdevols(1)%nb, tdevols(1)%nb))
  SAFE_ALLOCATE(gf_D,    (tdevols(1)%nb, tdevols(1)%nb))

  SAFE_ALLOCATE(commutator, (tdevols(1)%nb, tdevols(1)%nb))
  SAFE_ALLOCATE(commutators, (tdevols(1)%nb, tdevols(1)%nb, tdevols(1)%nk, dyn%n_finite_q))

  commutators(:,:,:,:) = (0.0d0,0.0d0)

  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  ! first we work with terms involving sum q`
  do iq = 1, dyn%n_finite_q
    ! here we have (k, q`) loops, we can impose parallelization directly and gather at the end without double counting terms
    ! we first work with sum over q` terms since we can directly allgather as they start as zeros
    ! here we parallelize nk and nq`
    do ik = tdevols(iq)%mpi_k_start, tdevols(iq)%mpi_k_end  ! serial version: do ik = 1, tdevols(iq)%nk
      ! then each (ik, iq) component sums over all q` points
      do iqp = tdevols(iq)%mpi_qp_start, tdevols(iq)%mpi_qp_end  ! serial version: do iqp = 1, dyn%n_finite_q
        !> == first part: terms without involving sum q`
        !     we put is under iqp loop to take advantage of the parallelization
        !     we take the treatment with delta_{q`,0} i.e. iqp == 1
        !!
        !     The reason why we do not put it outside the iqp loop is because we have mpi_allreduce later
        !     this contribution can only be counted once <=> controlled by delta_{q`,0}, i.e. iqp == 1
        !     our current implementation works for both one-level and two-level parallelization allreduce
        if (iqp == 1) then
          ! first get commutators of hamiltonian part as it does not involve sum over q`
          h_int_A(:,:) = tdevols(iq)%hamiltonian(:,:,ik)
          gf_B(:,:) = tdevols(iq)%gf(:,:,ik)
    
          ! standard commutation
          commutator(:,:) = (0.0d0,0.0d0)
          call commutation(tdevols(iq)%nb, h_int_A, gf_B, commutator)
          ! accumulate to commutators
          commutators(:,:,ik,iq) = commutators(:,:,ik,iq) + commutator(:,:)
        endif ! iqp == 1
        
        !> == second part: terms involving sum q`
        ! Q-Q` point
        qmqp(:) = tdevols(iq)%qpt(:) - tdevols(iqp)%qpt(:)
        ! find iqmqp index for Q-Q` in qpts
        iqmqp = kindex(qmqp, dyn%qpts, dyn%n_finite_q, skip_check=.true.)

        ! skip non-existing points
        if (iqmqp == 0) then
          !if (ik == 1) write(*,'(1x,a,3f10.6,a)') 'Warning: cannot find (', qmqp, ') q`-point in propagator summation.'
          cycle
        endif

        ! k+Q` point
        kpqp(:) = tdevols(iq)%kpts(:,ik) + tdevols(iqp)%qpt(:)
        ! find ikpqp index for k+Q` in kpts (of iqmqp), guaranteed to find
        ikpqp = kindex(kpqp, tdevols(iqmqp)%kpts, tdevols(iqmqp)%nk)

        !===> arrange four matrices matA, matB, matC, matD for generalized_commutation
        !=> 1) matA: h_int_A of (k, Q`)
        h_int_A(:,:) = - tdevols(iqp)%efield * tdevols(iqp)%dipole(:,:,ik) + tdevols(iqp)%interaction(:,:,ik)
 
        !=> 2) matB: gf_B of (k+Q`, Q-Q`)
        gf_B(:,:) = tdevols(iqmqp)%gf(:,:,ikpqp)

        !=> 3) matC: h_int_C of (k+Q`, Q-Q`)
        h_int_C(:,:) = - tdevols(iqmqp)%efield * tdevols(iqmqp)%dipole(:,:,ikpqp) + tdevols(iqmqp)%interaction(:,:,ikpqp)

        !=> 4) matD: gf_D of (k, Q`)
        gf_D(:,:) = tdevols(iqp)%gf(:,:,ik)

        ! compute commutator
        commutator(:,:) = (0.0d0,0.0d0)
        call generalized_commutation(tdevols(iqp)%nb, h_int_A, gf_B, h_int_C, gf_D, commutator)
        ! commutator of (ik, iq) accumates over iqp
        commutators(:,:,ik,iq) = commutators(:,:,ik,iq) + commutator(:,:)

      enddo ! iqp
    enddo ! ik
  enddo ! iq

  ! collect data
  ncount = tdevols(1)%nb * tdevols(1)%nb * tdevols(1)%nk * dyn%n_finite_q
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  call timing%start(timing%communication)
  ! MPI_IN_PLACE requires intent of recvbuf to be inout: commutators is declared local as variable
  call MPI_Allreduce(MPI_IN_PLACE, commutators, ncount, MPI_COMPLEX_DPC, MPI_SUM, MPI_COMM_WORLD, mpierr)
  call timing%stop(timing%communication)

  ! update all green`s functions after commutators are calculated
  do iq = 1, dyn%n_finite_q
    ! update green`s functions
    tdevols(iq)%gf(:,:,:) = tdevols(iq)%gf(:,:,:) + commutators(:,:,:,iq) * dyn%dt / imag_one

    ! update interaction
    call get_interaction(tdevols(iq))

    ! update time step label
    tdevols(iq)%itlabel = tdevols(iq)%itlabel + 1.0d0
  enddo ! iq

  ! if request saving, collect commutators here
  if (present(commute_save)) then
    commute_save(:,:,:,:) = commutators(:,:,:,:)
  endif

  SAFE_DEALLOCATE(h_int_A)
  SAFE_DEALLOCATE(gf_B)
  SAFE_DEALLOCATE(h_int_C)
  SAFE_DEALLOCATE(gf_D)

  SAFE_DEALLOCATE(commutator)
  SAFE_DEALLOCATE(commutators)

  POP_SUB(euler_propagator_q)

end subroutine euler_propagator_q


!=========================================================================
!> runge-kutta propagator (4th order implementation), default algorithm
subroutine runge_kutta_propagator_q(dyn, tdevols)

  type(dynaminfo), intent(in) :: dyn
  type(td_equation), intent(inout) :: tdevols(:)
  !
  type(td_equation), allocatable :: tdevs(:)
  complex(DPC), allocatable :: commus1(:,:,:,:), commus2(:,:,:,:), commus3(:,:,:,:), commus4(:,:,:,:)
  complex(DPC) :: imag_one
  real(DP), allocatable :: efields(:)
  integer :: it, iq

  PUSH_SUB(runge_kutta_propagator_q)

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

  SAFE_ALLOCATE(tdevs, (dyn%n_finite_q))
  SAFE_ALLOCATE(efields, (dyn%n_finite_q))

  ! prepare four commutators corresponding to i*k1, i*k2, i*k3, i*k4
  SAFE_ALLOCATE(commus1, (tdevols(1)%nb, tdevols(1)%nb, tdevols(1)%nk, dyn%n_finite_q))
  SAFE_ALLOCATE(commus2, (tdevols(1)%nb, tdevols(1)%nb, tdevols(1)%nk, dyn%n_finite_q))
  SAFE_ALLOCATE(commus3, (tdevols(1)%nb, tdevols(1)%nb, tdevols(1)%nk, dyn%n_finite_q))
  SAFE_ALLOCATE(commus4, (tdevols(1)%nb, tdevols(1)%nb, tdevols(1)%nk, dyn%n_finite_q))

  commus1(:,:,:,:) = (0.0d0,0.0d0)
  commus2(:,:,:,:) = (0.0d0,0.0d0)
  commus3(:,:,:,:) = (0.0d0,0.0d0)
  commus4(:,:,:,:) = (0.0d0,0.0d0)

  ! ==================== step 1 ===================
  !> k1 = f(t, y(t))        <=>    k1 = dG/dt|_t = -i [H(G(t)), G(t)]
  !
  do iq = 1, dyn%n_finite_q
    call init_td_equation_from_copy(tdevols(iq), tdevs(iq))
  enddo

  ! compute commus1
  call euler_propagator_q(dyn, tdevs, commute_save=commus1)

  do iq = 1, dyn%n_finite_q
    ! save this efield value
    efields(iq) = tdevs(iq)%efield
  
    call free_td_equation(tdevs(iq))
  enddo

  ! ==================== step 2 ====================
  ! k2 = f(t+0.5, y1(t+0.5))  <=>  k2 = dG`/dt|_{t+0.5} = -i [H`(G`(t+0.5)), G`(t+0.5)]
  ! y1 = G`(t+0.5) = y(t) + h * k1 / 2.0
  ! here G` denotes that this is a euler extrapolated quantity with specific input
  !
  do iq = 1, dyn%n_finite_q
    call init_td_equation_from_copy(tdevols(iq), tdevs(iq))
  
    ! U(t+0.5)
    tdevs(iq)%itlabel = tdevs(iq)%itlabel + 0.5d0
    ! G`(t+0.5) using commus1
    tdevs(iq)%gf(:,:,:) = tdevs(iq)%gf(:,:,:) + commus1(:,:,:,iq) * dyn%dt / imag_one / 2.0d0
    ! H`(G`(t+0.5))
    call get_interaction(tdevs(iq))
  enddo

  ! compute commus2
  call euler_propagator_q(dyn, tdevs, commute_save=commus2)

  do iq = 1, dyn%n_finite_q
    call free_td_equation(tdevs(iq))
  enddo

  ! ==================== step 3 ====================
  ! k3 = f(t+0.5, y2(t+0.5))  <=>  k3 = dG``/dt|_{t+0.5} = -i [H``(G``(t+0.5)), G``(t+0.5)]
  ! y2 = G``(t+0.5) = y(t) + h * k2 / 2.0
  ! note G`` is differently constructed as G`, even though they are both at the same time snapshot
  !     
  do iq = 1, dyn%n_finite_q
    call init_td_equation_from_copy(tdevols(iq), tdevs(iq))
  
    ! U(t+0.5)
    tdevs(iq)%itlabel = tdevs(iq)%itlabel + 0.5d0
    ! G``(t+0.5) using commus2
    tdevs(iq)%gf(:,:,:) = tdevs(iq)%gf(:,:,:) + commus2(:,:,:,iq) * dyn%dt / imag_one / 2.0d0
    ! H``(G``(t+0.5))
    call get_interaction(tdevs(iq))
  enddo

  ! compute commus3
  call euler_propagator_q(dyn, tdevs, commute_save=commus3)

  do iq = 1, dyn%n_finite_q
    call free_td_equation(tdevs(iq))
  enddo

  ! ==================== step 4 ====================
  ! k4 = f(t+1, y3(t+1))  <=>  k4 = dG/dt|_{t+1} = -i [H(G(t+1)), G(t+1)]
  ! y3 = G(t+1) = y(t) + h * k3
  !
  do iq = 1, dyn%n_finite_q
    call init_td_equation_from_copy(tdevols(iq), tdevs(iq))
  
    ! U(t+1)
    tdevs(iq)%itlabel = tdevs(iq)%itlabel + 1.0d0
    ! G(t+1) using commus3
    tdevs(iq)%gf(:,:,:) = tdevs(iq)%gf(:,:,:) + commus3(:,:,:,iq) * dyn%dt / imag_one
    ! H(G(t+1))
    call get_interaction(tdevs(iq))
  enddo

  ! commute commus4
  call euler_propagator_q(dyn, tdevs, commute_save=commus4)

  do iq = 1, dyn%n_finite_q
    call free_td_equation(tdevs(iq))
  enddo

  ! ==================== update ====================
  ! Runge-Kutta (4th-order) update
  !
  ! y(t+1) = y(t) + (1/6) * h * (k1 + 2*k2 + 2*k3 + k4)
  !
  do iq = 1, dyn%n_finite_q
    tdevols(iq)%efield = efields(iq)

    ! update green`s function
    tdevols(iq)%gf(:,:,:) = tdevols(iq)%gf(:,:,:) + (1.0d0 / 6.0d0) * (dyn%dt / imag_one) & 
                        & * (commus1(:,:,:,iq) + 2.0d0 * commus2(:,:,:,iq) + 2.0d0 * commus3(:,:,:,iq) + commus4(:,:,:,iq))

    ! get interaction
    call get_interaction(tdevols(iq))

    ! update time step label 
    tdevols(iq)%itlabel = tdevols(iq)%itlabel + 1.0d0
  enddo

  ! deallocate
  SAFE_DEALLOCATE(tdevs)
  SAFE_DEALLOCATE(efields)
  SAFE_DEALLOCATE(commus1)
  SAFE_DEALLOCATE(commus2)
  SAFE_DEALLOCATE(commus3)
  SAFE_DEALLOCATE(commus4)

  POP_SUB(runge_kutta_propagator_q)

end subroutine runge_kutta_propagator_q

end module propagation_q_m
