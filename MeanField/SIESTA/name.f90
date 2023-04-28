
!===============================================================================
!
! Routines:
!
! 1. xvfilename()        Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Returns SIESTA file name SystemLabel.XV
!
! 2. kpfilename()        Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Returns SIESTA file name SystemLabel.KP
!
! 3. eigfilename()       Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Returns SIESTA file name SystemLabel.EIG
!
! 4. wfnrfilename()      Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Returns SIESTA file name SystemLabel{.K$k}.WF$n{.UP|.DOWN}{.REAL|.IMAG}.cube
!
! 5. rhorfilename()      Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Returns SIESTA file name SystemLabel.RHO{.UP|.DN}.cube
!
! 6. vtrfilename()       Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Returns SIESTA file name SystemLabel.VT{.UP|.DN}.cube
!
! 7. vhrfilename()       Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Returns SIESTA file name SystemLabel.VH.cube
!
!===============================================================================

#include "f_defs.h"

module name_m

  use global_m

  implicit none

  private

  public :: xvfilename, kpfilename, eigfilename, &
    wfnrfilename, rhorfilename, vtrfilename, vhrfilename

contains

!-------------------------------------------------------------------------------

subroutine xvfilename(systemlabel,filename)

  character(len=256), intent(in) :: systemlabel
  character(len=256), intent(out) :: filename
  
  PUSH_SUB(xvfilename)
  
  filename = TRUNC(systemlabel) // '.XV'
  
  POP_SUB(xvfilename)

  return

end subroutine xvfilename

!-------------------------------------------------------------------------------

subroutine kpfilename(systemlabel,filename)
  
  character(len=256), intent(in) :: systemlabel
  character(len=256), intent(out) :: filename
  
  PUSH_SUB(kpfilename)
  
  filename = TRUNC(systemlabel) // '.KP'
  
  POP_SUB(kpfilename)
  
  return
  
end subroutine kpfilename

!-------------------------------------------------------------------------------

subroutine eigfilename(systemlabel,filename)
  
  character(len=256), intent(in) :: systemlabel
  character(len=256), intent(out) :: filename
  
  PUSH_SUB(eigfilename)
  
  filename = TRUNC(systemlabel) // '.EIG'
  
  POP_SUB(eigfilename)
  
  return
  
end subroutine eigfilename

!-------------------------------------------------------------------------------

subroutine wfnrfilename(ik,is,ib,ip,ns,np,systemlabel,filename)
  
  integer, intent(in) :: ik,is,ib,ip,ns,np
  character(len=256), intent(in) :: systemlabel
  character(len=256), intent(out) :: filename
  
  character(len=16) :: c1,c2,s1,s2,s3,s4
  
  PUSH_SUB(wfnrfilename)
  
  write(c1,301)ik
  write(c2,301)ib
  
  if (np .eq. 1) then
    s1 = ''
  else
    s1 = '.K' // TRUNC(c1)
  endif
  
  s2 = '.WF' // TRUNC(c2)
  
  if (ns .eq. 1) then
    s3 = ''
  else
    if (is .eq. 1) then
      s3 = '.UP'
    elseif (is .eq. 2) then
      s3 = '.DOWN'
    else
      s3 = ''
    endif
  endif
  
  if (np .eq. 1) then
    s4 = ''
  else
    if (ip .eq. 1) then
      s4 = '.REAL'
    elseif (ip .eq. 2) then
      s4 = '.IMAG'
    else
      s4 = ''
    endif
  endif
  
  filename = TRUNC(systemlabel) // TRUNC(s1) // TRUNC(s2) // &
    TRUNC(s3) // TRUNC(s4) // '.cube'
  
  POP_SUB(wfnrfilename)
  
  return
  
301 format(i16)
  
end subroutine wfnrfilename

!-------------------------------------------------------------------------------

subroutine rhorfilename(is,ns,systemlabel,filename)
  
  integer, intent(in) :: is,ns
  character(len=256), intent(in) :: systemlabel
  character(len=256), intent(out) :: filename
  
  character(len=16) :: s1
  
  PUSH_SUB(rhorfilename)
  
  if (ns .eq. 1) then
    s1 = ''
  else
    if (is .eq. 1) then
      s1 = '.UP'
    elseif (is .eq. 2) then
      s1 = '.DN'
    else
      s1 = ''
    endif
  endif
  
  filename = TRUNC(systemlabel) // '.RHO' // TRUNC(s1) // '.cube'
  
  POP_SUB(rhorfilename)
  
  return
  
end subroutine rhorfilename

!-------------------------------------------------------------------------------

subroutine vtrfilename(is,ns,systemlabel,filename)
  
  integer, intent(in) :: is,ns
  character(len=256), intent(in) :: systemlabel
  character(len=256), intent(out) :: filename
  
  character(len=16) :: s1
  
  PUSH_SUB(vtrfilename)
  
  if (ns .eq. 1) then
    s1 = ''
  else
    if (is .eq. 1) then
      s1 = '.UP'
    elseif (is .eq. 2) then
      s1 = '.DN'
    else
      s1 = ''
    endif
  endif
  
  filename = TRUNC(systemlabel) // '.VT' // TRUNC(s1) // '.cube'
  
  POP_SUB(vtrfilename)
  
  return
  
end subroutine vtrfilename

!-------------------------------------------------------------------------------

subroutine vhrfilename(systemlabel,filename)
  
  character(len=256), intent(in) :: systemlabel
  character(len=256), intent(out) :: filename
  
  PUSH_SUB(vhrfilename)
  
  filename = TRUNC(systemlabel) // '.VH' // '.cube'
  
  POP_SUB(vhrfilename)
  
  return
  
end subroutine vhrfilename

!-------------------------------------------------------------------------------

end module name_m

