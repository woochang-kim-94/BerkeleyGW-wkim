#include "f_defs.h"

module print_error_m

  use global_m
  implicit none

  private

  public :: &
      print_error


contains

!This routine is very simple now, probably this has to be substituted by a BGW routine.
   subroutine print_error(a,b,n)
    IMPLICIT NONE

    CHARACTER(LEN=*)    :: a, b
    INTEGER, INTENT(IN) :: n
    INTEGER :: ip, nproc, mpime
! no push/pop since called too frequently.

    if(n==1) then    
      WRITE (0,'("ERROR in ",a,":")') trim(a)
    else 
      WRITE (0,'("WARNING in ",a,":")') trim(a)
    end if
    WRITE (0,'(a)') trim(b)
    if(n==1) then
! no push/pop since called too frequently.
      call die("Stopping")
    end if

! no push/pop since called too frequently.
  end subroutine print_error


end module print_error_m
