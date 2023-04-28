!===============================================================================
!
! Routines:
!
! 1. groupk()        Originally By gsm       Last Modified 8/31/2010 (gsm)
!
!    Generates symmetries of the k-point.
!
!===============================================================================

#include "f_defs.h"

module groupk_m

  use global_m
  implicit none

  private

  public :: groupk

contains

subroutine groupk(kvec,nsyml,rotl,nsymk,rotk)
  real(DP), intent(in) :: kvec(:) !< (3)
  integer, intent(in) :: nsyml
  integer, intent(in) :: rotl(:,:,:) !< (3,3,48)
  integer, intent(out) :: nsymk
  integer, intent(out) :: rotk(:,:,:) !< (3,3,48)
  
  integer :: isyml
  real(DP) :: krot(3)
  
  PUSH_SUB(groupk)

  nsymk=0
  do isyml=1,nsyml
    krot(1:3) = matmul(dble(rotl(1:3, 1:3, isyml)), kvec(1:3))
    if (all(abs(krot(1:3)-kvec(1:3)).lt.TOL_Small)) then
      nsymk=nsymk+1
      rotk(1:3, 1:3, nsymk) = rotl(1:3, 1:3, isyml)
    endif
  enddo
  
  POP_SUB(groupk)
  return  
end subroutine groupk

end module groupk_m
