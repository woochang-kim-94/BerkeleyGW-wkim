!============================================================================
!
! Routines:
!
! (1) oscstrength()     Originally By JRD       Last Modified 6/6/2008 (JRD)
!
!     Calculates and exciton-exciton transition oscilator strength between S`
!     (Af) and S (A).
!
!============================================================================

#include "f_defs.h"

module oscstrength_m

  use global_m
  implicit none
  
  private
  public :: oscstrength

contains

subroutine oscstrength(nmat,nk,ns,xct,Af,A,s1,osc)
  integer, intent(in) :: nmat,nk,ns  
  type (xctinfo), intent(inout) :: xct
  SCALAR, intent(in) :: Af(:,:,:,:) !< (ns,xct%nvb_fi,xct%ncb_fi,nk)
  SCALAR, intent(in) :: A(:,:,:,:) !< (ns,xct%nvb_fi,xct%ncb_fi,nk)
  SCALAR, intent(in) :: s1(:) !< (nmat)
  SCALAR, intent(out) :: osc

  integer :: iv1,iv2,ic1,ic2,ik,is,ikcvs
  
  SCALAR :: tmp,s1tmp

  PUSH_SUB(oscstrength)

  osc=0d0
  tmp=0d0

  do is=1,xct%nspin
    do ik=1,xct%nkpt_fi
      do iv2=1,xct%nvb_fi
        do iv1=1,xct%nvb_fi
          
          ikcvs=is + (iv1 - 1 + (iv2 - 1 + (ik - 1)*(xct%ncb_fi+ &
            xct%nvb_fi))*(xct%nvb_fi+xct%ncb_fi))*xct%nspin
          
          s1tmp = s1(ikcvs)
          
          do ic1=1,xct%ncb_fi
            
            tmp = tmp + MYCONJG(Af(is,xct%nvb_fi-iv1+1,ic1,ik))* &
              A(is,xct%nvb_fi-iv2+1,ic1,ik)*s1tmp
          enddo
        enddo
      enddo
    enddo
  enddo
  
  do is=1,xct%nspin
    do ik=1,xct%nkpt_fi
      do ic2=1,xct%ncb_fi
        do ic1=1,xct%ncb_fi
          
          ikcvs=is + ((xct%nvb_fi+ic1) - 1 + ((xct%nvb_fi+ic2-1) + &
            (ik - 1)*(xct%ncb_fi+xct%nvb_fi))* &
            (xct%nvb_fi+xct%ncb_fi))*xct%nspin
          
          s1tmp = s1(ikcvs) 
          
          do iv1=1,xct%nvb_fi
            
            tmp = tmp + MYCONJG(Af(is,iv1,ic1,ik))* &
              A(is,iv1,ic2,ik)*s1tmp
          enddo
        enddo
      enddo
    enddo
  enddo
  
  osc= tmp
  
  POP_SUB(oscstrength)
  
  return
end subroutine oscstrength

end module oscstrength_m
