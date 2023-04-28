!===============================================================================
!
! Routines:
!
! (1) mtxel_jdos()         Originally By MJ               Last Modified 22/7/2010 MJ
!
!     input: s1 and nmat
!
!     output: s1  matrix element of the JDOS operator
!
!     Calculates the JDOS operator 
!     exp(phi(ic,iv))  where phi(ic,iv) = random number between 0 and 2pi
!
!===============================================================================

#include "f_defs.h"

module mtxel_jdos_m

  use global_m
  use random_m
  implicit none

  private

  public :: &
    mtxel_jdos

contains

subroutine mtxel_jdos(s1,nmat)
  integer, intent(in) :: nmat
  SCALAR, intent(inout) :: s1(nmat)
  
  integer :: i
  real(DP) :: rnd,sq2
  complex(DPC) :: zi
  SCALAR :: random

  PUSH_SUB(mtxel_jdos)

  s1 = ZERO
  zi = (0.0d0,1.0d0)
  sq2 = sqrt(2.0d0)

!----------------------------------
! Calculate s1(iv) = exp(phi)

  do i = 1, nmat
    call genrand_real4(rnd)
#ifdef CPLX
    random = exp(zi*2.0d0*PI_D*rnd)
#else
    random = sq2*cos(2.0d0*PI_D*rnd)
#endif
    s1(i) = random
  enddo  

  POP_SUB(mtxel_jdos)

  return
end subroutine mtxel_jdos

end module mtxel_jdos_m
