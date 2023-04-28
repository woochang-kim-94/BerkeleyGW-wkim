!=============================================================================
!
! Routines:
!
! (1) distrib()         Originally By MLT       Last Modified 6/2008 FJR
!
!     Distributes kpoints among the processors
!
!=============================================================================

#include "f_defs.h"

module distrib_m

  use global_m
  implicit none

  private

  public :: &
    distrib

contains

subroutine distrib(xct)
  type (xctinfo), intent(in) :: xct
  
  integer :: ik,ipe
  
  PUSH_SUB(distrib)

  peinf%nkpe=xct%nn/peinf%npes
  if(xct%nn-peinf%npes*peinf%nkpe.gt.0) peinf%nkpe=peinf%nkpe+1
  
  SAFE_ALLOCATE(peinf%ik, (peinf%npes,peinf%nkpe))
  SAFE_ALLOCATE(peinf%ikt, (peinf%npes))
  
  peinf%ik=0
  peinf%ikt=0
  
  ipe=0
  do ik=1,xct%nn
    ipe=ipe+1
    if(ipe.eq.peinf%npes+1) ipe=1
    peinf%ikt(ipe)=peinf%ikt(ipe)+1
    peinf%ik(ipe,peinf%ikt(ipe))=ik
  enddo
  
  if(peinf%inode.eq.0) then
    write(6,70)
70  format(/,1x,'Number of kpoints treated by each PE:')
    write(6,80) peinf%ikt
80  format(3x,8i4)
  endif
  
  POP_SUB(distrib)
  
  return
end subroutine distrib

end module distrib_m
