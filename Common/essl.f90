!================================================================================
!
! Modules:
!
! (1) essl_m      Originally By DAS      Last Modified 2/2/2011 (das)
!
!     Interfaces for ESSL functions.
!     Every ESSL function used in the code should be listed here, and this
!     module should be used in every routine containing ESSL calls to ensure
!     the argument types are correct.
!
!     Note that routines can take extra arguments like "*100" which directs to
!     go to line 100 if error 100 occurs. This is an "alternate return specifier"
!     considered an obsolescent feature in Fortran 95 and Fortran 90.
!     It does not appear to be possible to declare in an interface in a way
!     that will not produce compilation warnings.
!
! http://publib.boulder.ibm.com/infocenter/clresctr/vxrx/index.jsp?topic=/com.ibm.cluster.essl43.guideref.doc/am501mst_xtoc.html
!
! In the idiotic tradition of IBM programming, ESSL now bombs
! if the matrix is singular regardless even when calling
! supposedly non-bombing routines like dgeicd... so instead
! we have to catch the error #2103 which means the matrix is singular.
! Beautiful, no!?!?
!
!================================================================================

#include "f_defs.h"

module essl_m

  public ! only interfaces in this module

#ifdef USEESSL
  interface
    subroutine einfo(icode, inf1, inf2)
      implicit none
      integer, intent(in) :: icode
      integer, optional, intent(out) :: inf1
      integer, optional, intent(out) :: inf2
    end subroutine einfo
  end interface

  interface
    subroutine errset(ierno, inoal, inomes, itrace, iusadr, irange)
      implicit none
      integer, intent(in) :: ierno
      integer, intent(in) :: inoal
      integer, intent(in) :: inomes
      integer, intent(in) :: itrace
      integer, intent(in) :: iusadr
      integer, intent(in) :: irange
    end subroutine errset
  end interface

  interface
    subroutine dgef(a, lda, n, ipvt, *)
      implicit none
      real*8, intent(inout) :: a !< (lda,n)
      integer, intent(in) :: lda
      integer, intent(in) :: n
      integer, intent(out) :: ipvt !< (n)
    end subroutine dgef
  end interface

  interface
    subroutine zgef(a, lda, n, ipvt)
      implicit none
      complex*16, intent(inout) :: a !< (lda,n)
      integer, intent(in) :: lda
      integer, intent(in) :: n
      integer, intent(out) :: ipvt !< (n)
    end subroutine zgef
  end interface

  interface
    subroutine dgesm(trans, a, lda, n, ipvt, bx, ldb, nrhs)
      implicit none
      character, intent(in) :: trans
      real*8, intent(in) :: a !< (lda, n)
      integer, intent(in) :: lda
      integer, intent(in) :: n
      integer, intent(in) :: ipvt !< (n)
      real*8, intent(inout) :: bx !< (lda, nrhs)
      integer, intent(in) :: ldb
      integer, intent(in) :: nrhs
    end subroutine dgesm
  end interface

  interface
    subroutine zgesm(trans, a, lda, n, ipvt, bx, ldb, nrhs)
      implicit none
      character, intent(in) :: trans
      complex*16, intent(in) :: a !< (lda, n)
      integer, intent(in) :: lda
      integer, intent(in) :: n
      integer, intent(in) :: ipvt !< (n)
      complex*16, intent(inout) :: bx !< (lda, nrhs)
      integer, intent(in) :: ldb
      integer, intent(in) :: nrhs
    end subroutine zgesm
  end interface

#endif

end module essl_m
