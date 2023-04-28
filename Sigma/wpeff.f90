!=============================================================================
!
! Routines:
!
! (1) wpeff()    Originally By SIB              Last Modified 5/1/2008 (JRD)
!
!     This routine computes effective plasma frequencies (squared).
!
!     Specifically, for each G` (igp) it computes for all G (ig),
!     the quantity in formula (31) on page 5396 of
!     Hybertsen and Louie, PRB vol 34, #8 (1986) given by
!
!     Omega(G,G`)^2 = wp^2 * [rho(G-G`)/rho(0)] * (q+G).(q+G`)*vc(q+G)/(8pi)
!
!     (vc(q+G) is the Coulomb interaction), and places the numbers
!     into wpmtx(ig).  Units are eV^2.
!
!     Computes effective plasma freq matrix-
!     wp(g,g`) = wp**2 (rho(g-g`)/rho(0)) (q+g).(q+g`) vc(q+g)/8pi
!
!     isrtrq             index array, (rq+g)**2
!     nelec,wpsq         rho(g=0), feg plasma frequency^2
!     rho                rho(g)
!
!=============================================================================

#include "f_defs.h"

module wpeff_m

  use global_m
  use misc_m
  implicit none

  private

  public :: &
    wpeff

contains

subroutine wpeff(crys, gvec, wpg, sig, neps, isrtrq, igp, ncouls, wpmtx, nspin, qk,vcoul,coulfact)
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (wpgen), intent(in) :: wpg
  type (siginfo), intent(in) :: sig
  integer, intent(in) :: neps
  integer, intent(in) :: isrtrq(neps), igp, ncouls
  SCALAR, intent(out) :: wpmtx(neps)
  integer, intent(in) :: nspin
  real(DP), intent(in) :: qk(3), vcoul(ncouls), coulfact

  SCALAR :: rho_g_minus_gp
  integer :: ig, igadd, igpadd, kadd, gg(3)
  real(DP) :: qg(3), qgp(3)
  logical, save :: warned=.false.

!-----------------------------
! SIB: stuff that is remembered between calls in order to speed
! things up
! 'save' is deprecated, replace with an object... DAS

  !logical, save :: first_call = .true.
  !logical, save :: q_is_not_zero
  !real(DP), save :: qk_old(3), fact, qnorm2
  !real(DP), allocatable, save :: precalc(:,:)
  real(DP) :: fact, qnorm2
  logical :: q_is_not_zero
  real(DP), allocatable :: precalc(:,:)

  ! no push_sub, called too frequently

! WARNING: if Fermi level is adjusted, the charge density here does not change,
! which will make the plasma frequency inconsistent. -- DAS

  if((abs(sig%efermi_input) > TOL_Zero .or. .not. sig%rfermi) &
    .and. (peinf%inode == 0) .and..not.warned) then
    write(0,'(a)') "WARNING: GPP plasma frequency does not respond to manual adjustment of the Fermi level."
    warned = .true.
  endif

!--------------- Begin Calculation ------------------------------------------------------


!------------------------
! SIB: Calculate stuff that is the same between calls for a given qk.
! The precalculated table contains:
!
! precalc(:,ig) = wp**2/rho(0)*(vc(q+G)/8pi)*(q+G)

  !if (sum(abs(qk-qk_old)) > 1.0d-12) then
    SAFE_ALLOCATE(precalc, (3,ncouls))
    precalc = 0.0d0

! Common factor

    fact=(wpg%wpsq(1)+(nspin-1)*wpg%wpsq(2))/(wpg%nelec(1)+(nspin-1)*wpg%nelec(2))

! Square length of qk

    qnorm2=dot_product(qk,matmul(crys%bdot,qk))

! Set a flag if q is not "zero"

    if (qnorm2 > 1.0d-12) then
      q_is_not_zero = .true.
    else
      q_is_not_zero = .false.
    endif


!---------- Loop over g and tabulate -------------------------------

    do ig=1,ncouls

! Get q+g and |q+g|^2

      igadd=isrtrq(ig)

! If g<>0 or q<>0, we just calculate the formula

      if (igadd.ne.1 .or. q_is_not_zero) then
        qg(:)=gvec%components(:,igadd)+qk(:)
        precalc(:,ig) = fact*vcoul(ig)*qg(:)/coulfact
      else

! If g=q=0, we have to avoid dividing by zero;
! we handle this special case separately below

        precalc(:,ig) = 0.0d0
      endif

    enddo ! g loop
  !endif


!----------------- Here starts the main calculation! ----------------------------


! Get q+gp

  igpadd=isrtrq(igp)
  qgp(:)=gvec%components(:,igpadd)+qk(:)


!!----------- Loop over g and calculate Omega^2(g,gp) -------------------------


  do ig=1,ncouls
    igadd=isrtrq(ig)
    
    wpmtx(ig) = ZERO

! Compute address of g-gp, and if it is a vector for which
! we have the density, get rho(g-gp); if out of bounds,
! skip this g
    
    gg(:) = gvec%components(:, igadd) - gvec%components(:, igpadd)
    call findvector(kadd, gg, gvec)
    if(kadd.eq.0) cycle
    rho_g_minus_gp = sum(wpg%rho(kadd,1:nspin))

! Using precalculated stuff, assemble together wpmtx
! if g<>0 or q<>0, we can just do what the formula says

    if (igadd.ne.1 .or. q_is_not_zero) then
      wpmtx(ig) = dot_product(qgp(:),matmul(crys%bdot,precalc(:,ig)))*rho_g_minus_gp

! The special case q=g=0

    else

! If gp=0, then wpmtx=wp**2 times a possible
! (1-cos(0))=0 if truncating

      if (igpadd.eq.1) then
        wpmtx(ig) = fact*rho_g_minus_gp
!            if (sig%icutv/=TRUNC_NONE) wpmtx(ig) = 0.0d0

! When q=g=0, gp<>0, the result is just set to zero
! (it is infinite for true coulomb interaction and zero
! if truncated).

      else
        wpmtx(ig) = 0.0d0
      endif
    endif

  enddo !! g loop

  SAFE_DEALLOCATE(precalc)

  return
end subroutine wpeff

end module wpeff_m
