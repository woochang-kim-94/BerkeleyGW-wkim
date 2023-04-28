!==========================================================================
!
! Routines:
!
! 1. subgrp()   Originally By ?                 Last Modified 6/12/2008 (JRD)
!
!    Determines a subgroup of the symmetry group that preserves a q-vector.
!
!===========================================================================

#include "f_defs.h"

module subgrp_m

  use global_m
  use misc_m
  implicit none

  private

  public :: &
    subgrp

contains

subroutine subgrp(qq, syms)
  type (symmetry), intent(inout) :: syms
  real(DP), intent(in) :: qq(3)

  integer ::  it, kg(3)
  real(DP) :: qk(3), dqk(3)

  PUSH_SUB(subgrp)
  
!---------------------
! Loop over transformations testing for r(q) = q + kg0

  syms%rq = qq
  syms%ntranq = 0
  do it = 1, syms%ntran
    qk(1:3) = matmul(syms%mtrx(1:3, 1:3, it), qq(1:3))
    dqk(1:3) = qk(1:3) - qq(1:3)
    call k_range(dqk(1:3), kg(1:3), TOL_Small)
    
    if (all(abs(dqk(1:3)) .lt. TOL_Small)) then
      !--------------------
      ! Store index of element of subgroup
      syms%ntranq = syms%ntranq + 1
      syms%indsub(syms%ntranq) = it
      syms%kgzero(1:3, syms%ntranq) = -kg(1:3)
    endif
  enddo
  
  POP_SUB(subgrp)
  
  return
end subroutine subgrp

end module subgrp_m
