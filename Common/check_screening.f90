!============================================================================
!
! Routines:
!
! (1) check_screening_trunc   Originally by JRD     Last Modified: 2/09/2009 (JRD)
!
!     Die if screening, truncation, and q0vec are not set consistently.
!
!============================================================================

#include "f_defs.h"

module check_screening_m

  use global_m
  implicit none

  private

  public :: check_screening_trunc

contains

subroutine check_screening_trunc(itruncflag,iscreen,q0vec,bdot)
  integer, intent(in) :: itruncflag
  integer, intent(in) :: iscreen
  real(DP), intent(in) :: q0vec(3)
  real(DP), intent(in) :: bdot(3,3)
  
  real(DP) :: q0len

  PUSH_SUB(check_screening_trunc)

  q0len = sqrt(DOT_PRODUCT(q0vec,MATMUL(bdot,q0vec)))

  if (iscreen==SCREEN_METAL .and. q0len<TOL_SMALL) then
    if (peinf%inode == 0) then
      write(0,'()')
      write(0,'(a)') 'ERROR: cannot use metallic screening with q0 = 0.'
      write(0,'(a)') 'You should either specify a nonzero q0->0 vector or use another screening flag.'
      write(0,'()')
    endif
    call die('Inconsistent screening, trunction, or q0 vector')
  endif

  if ((itruncflag==TRUNC_NONE .or. itruncflag==TRUNC_WIRE .or. &
    itruncflag==TRUNC_SLAB) .and. q0len<TOL_SMALL) then
    if (peinf%inode == 0) then
      write(0,'()')
      write(0,'(a)') 'ERROR: the input truncation flag indicates that the Coloumb interaction v(q0)'
      write(0,'(a)') 'diverges for q0->0. However, you have q0 exactly zero.'
      write(0,'(a)') 'You should always specify a *nonzero* q0->0 vector unless you have 0D'
      write(0,'(a)') 'truncation, i.e., spherical or box truncation.'
      write(0,'()')
    endif
    call die('Inconsistent screening, trunction, or q0 vector')
  endif
  
  if ((itruncflag/=TRUNC_NONE .and. itruncflag/=TRUNC_WIRE .and. &
    itruncflag/=TRUNC_SLAB) .and. iscreen/=SCREEN_METAL .and. q0len>=TOL_SMALL) then
    if (peinf%inode == 0) then
      write(0,'()')
      write(0,'(a)') 'ERROR: the input truncation '
      write(0,*) 'You want semiconductor or graphene screening with truncation'
      write(0,*) 'but specified nonzero q0vec!!'
    endif
    call die('Inconsistent Screening', only_root_writes = .true.)
  endif
  
  POP_SUB(check_screening_trunc)
  
end subroutine check_screening_trunc

end module check_screening_m
