!==========================================================================
!
! Routines:
!
! (1) find_kpt_match()    Originally by SIB
!
!   Look for rkq in the list of kpoints rotated by the symmetries.
!   If found, it puts its index in ikrkq, itqq is the index
!   of the symmetry that worked, and kgqq is an "umklapp" vector (i.e.
!   integer 3-vector) so that rkq = symm*kvec + kgqq.
!
!==========================================================================

#include "f_defs.h"

module find_kpt_match_m

  use global_m
  implicit none

  private
  public :: find_kpt_match

contains

  subroutine find_kpt_match(kp_point, syms, rkq, ikrkq, itqq, kgqq)
    type(kpoints), intent(in) :: kp_point
    type(symmetry), intent(in) :: syms
    real(DP), intent(in) :: rkq(3)
    integer, intent(out) :: ikrkq
    integer, intent(out) :: itqq
    integer, intent(out) :: kgqq(3)

    integer :: ik, itq, ii
    real(DP) :: qk(3), del(3)

    PUSH_SUB(find_kpt_match)

    ikrkq = 0
    ik_loop: do ik = 1, kp_point%nrk
      do itq = 1, syms%ntran
        qk(1:3) = matmul(syms%mtrx(1:3, 1:3, itq), kp_point%rk(1:3, ik))
        do ii = 1, 3
          del(ii) = rkq(ii) - qk(ii)
          if(del(ii).ge.0.0d0) kgqq(ii) = del(ii) + TOL_Small
          if(del(ii).lt.0.0d0) kgqq(ii) = del(ii) - TOL_Small
        enddo
        if(all(abs(del(1:3)-kgqq(1:3)) .lt. TOL_Small)) then
          ikrkq=ik
          itqq = itq
          exit ik_loop
        endif
      enddo
    enddo ik_loop
    
    POP_SUB(find_kpt_match)
    return

  end subroutine find_kpt_match

end module find_kpt_match_m
