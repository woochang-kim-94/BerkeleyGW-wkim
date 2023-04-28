!===========================================================================
!
! Routines:
!
! (1) rqstar()          Originally By SIB       Last Modified 9/1/2010 (DAS)
!
!     Takes vector rq (in irr. zone) and applies all the symmetries
!     of current q-vector to it, and compiles a list of the unique vectors
!     thus generated.  On exit,
!       nstar = number of such vectors
!       insdt = indices of symmetries which generate them, where indices
!               are in the list of those of q (subgroup).
!
!===========================================================================

#include "f_defs.h"

module rqstar_m

  use global_m
  use misc_m
  implicit none

  private

  public :: &
    rqstar

contains

subroutine rqstar(syms,nstar,indst,rq)
  type (symmetry), intent(in) :: syms
  integer, intent(out) :: nstar,indst(*)
  real(DP), intent(in) :: rq(3)

  integer :: it,istar,gpt(3)
  real(DP) :: qk(3),rqs(3,48)
  logical :: found

!-------------- loop over elements of subgroup ---------------------------

  PUSH_SUB(rqstar)

  nstar = 0
  do it = 1, syms%ntranq

! Rotate rq, so qk = syms%mtrix(indsub(it))*rq,
! and ensure qk(i) is between 0 and 1.

    qk(1:3) = matmul(syms%mtrx(1:3, 1:3, syms%indsub(it)), rq(1:3))
    call k_range(qk(1:3), gpt(1:3), TOL_Small)

! Compare to other elements of star to see if it is already present
    found = .false.
    do istar = 1, nstar
      if (all(abs(qk(1:3) - rqs(1:3, istar)) .lt. TOL_Small)) then
        found = .true.
        exit
      endif
    enddo
    
    if(.not. found) then
! Store new element of star and rotation which gives it
      nstar = nstar + 1
      rqs(1:3, nstar) = qk(1:3)
      indst(nstar) = it
    endif
  enddo

!-------------- end loop over elements of subgroup -----------------------

  POP_SUB(rqstar)
  
  return
end subroutine rqstar

end module rqstar_m
