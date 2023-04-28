!============================================================================
!
! Routines:
!
! (1) distrib()         Originally By JRD       Last Modified (JRD) 5/29/2008
!
!     Distributes k-points among processors.  It is best if number of kpoints
!     divides number of processor evenly, otherwise who knows...
!
!     nk = # of kpts
!     np = # of processors
!
!     peinf%nkpe = floor(nk/np) = # of kpts per proc
!     peinf%ikt(p) = # of kpoints belonging to proc p
!     peinf%ik(p,i) = kpt label (1..xct%nkpt_fi) of ith kpts on proc p
!
!     input:  xct%nkpt_fi
!             xct%ncb_fi
!             xct%nvb_fi
!             xct%neps
!             xct%nspin
!             peinf type
!
!      output: peinf type
!
!============================================================================

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

!      write(6,*) 'In distrib xct%nkpt_fi = ',xct%nkpt_fi
!      write(6,*) 'peinf%npes =',peinf%npes


!------------- Distribute work ----------------------------------------------


  PUSH_SUB(distrib)
  
  peinf%nkpe=xct%nkpt_fi/peinf%npes
  if(xct%nkpt_fi-peinf%npes*peinf%nkpe.gt.0) peinf%nkpe=peinf%nkpe+1
  
  SAFE_ALLOCATE(peinf%ik, (peinf%npes,peinf%nkpe))
  SAFE_ALLOCATE(peinf%ikt, (peinf%npes))
  
  peinf%ik=0
  peinf%ikt=0
  
  ipe=0
  do ik=1,xct%nkpt_fi
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
    write(6,*) ' '
    if (peinf%nkpe*peinf%npes.ne.xct%nkpt_fi) then
      write(0,*) '...k-points not evenly distributed among PEs.'
      write(0,*) 'This job will not run with optimum load balance.'
    endif
  endif
  
  POP_SUB(distrib)
  
  return
end subroutine distrib

end module distrib_m
