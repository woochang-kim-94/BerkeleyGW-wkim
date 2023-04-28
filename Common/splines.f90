!===============================================================================
!
! MODULE:
!
!   splines_m     Originally by FHJ     Last Modified 10/24/2011 (FHJ)
!
!>    Routines to evaluate spline curves.
!
! DESCRIPTION:
!
!     Main functions were adapted from P. Dierckx`s FITPACK (BSD License):
!       http://www.netlib.org/dierckx/
!       http://nalag.cs.kuleuven.be/research/topics/fitpack.shtml
!       http://www.scipy.org/doc/api_docs/SciPy.interpolate.fitpack.html
!
!===============================================================================

#include "f_defs.h"

module splines_m

  use global_m
  implicit none

  private

  public :: splev, splev_shift

contains

  !The following functions were downloaded from:
  ! http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=dierckx%2Fsplev.f
  !The code was converted to F90 using Alan Miller`s 'to_f90' utility:
  ! http://jblevins.org/mirror/amiller/to_f90.f90


  !> Subroutine fpbspl evaluates the (k+1) non-zero b-splines of degree k at
  !! t(l) <= x < t(l+1) using the stable recurrence relation of de Boor and Cox.
  subroutine fpbspl(t,n,k,x,l,h)
    real(DP), intent(in)   :: t(:) !< (n)
    integer, intent(in)    :: n
    integer, intent(in)    :: k
    real(DP), intent(in)   :: x
    integer, intent(in)    :: l
    real(DP), intent(out)  :: h(6)

    !  ..local scalars..
    real(DP) :: f,one
    integer :: i,j,li,lj
    !  ..local arrays..
    real(DP) :: hh(5)
    !  ..

    PUSH_SUB(fpbspl)

    one = 0.1d+1
    h(1) = one
    do j=1,k
      do i=1,j
        hh(i) = h(i)
      end do
      h(1) = 0.0d0
      do i=1,j
        li = l+i
        lj = li-j
        f = hh(i)/(t(li)-t(lj))
        h(i) = h(i)+f*(t(li)-x)
        h(i+1) = f*(x-t(lj))
      end do
    end do

    POP_SUB(fpbspl)
  end subroutine fpbspl

  !> Subroutine splev evaluates in a number of points x(i),i=1,2,...,m
  !!   a spline s(x) of degree k, given in its b-spline representation.
  !!
  !! Restrictions:
  !!   \li m >= 1
  !!   \li t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
  !!
  !!   \param t  array, length n, which contains the position of the knots.
  !!   \param n  integer, giving the total number of knots of s(x).
  !!   \param c  array, length n, which contains the b-spline coefficients.
  !!   \param k  integer, giving the degree of s(x).
  !!   \param x  array, length m, which contains the points where s(x) must
  !!     be evaluated.
  !!   \param m  integer, giving the number of points where s(x) must be
  !!     evaluated.
  !!   \param y  array,length m, giving the value of s(x) at the different
  !!     points.
  !!   \param ier error flag:
  !!     \li ier = 0:  normal return
  !!     \li ier = 10: invalid input data (see restrictions)
  subroutine splev(t,n,c,k,x,y,m,ier)
    !  Other subroutines required: fpbspl.

    !  References :
    !    de boor c  : on calculating with b-splines, j. approximation theory
    !                 6 (1972) 50-62.
    !    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
    !                 applics 10 (1972) 134-149.
    !    dierckx p. : curve and surface fitting with splines, monographs on
    !                 numerical analysis, oxford university press, 1993.

    !  Author :
    !    p.dierckx
    !    dept. computer science, k.u.leuven
    !    celestijnenlaan 200a, b-3001 heverlee, belgium.
    !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be

    !  Latest Update : march 1987
    !  Removed go to usage, DAS Jan 2012

    !  ..scalar arguments..

    real(DP), intent(in)      :: t(:) !< (n)
    integer, intent(in)       :: n
    real(DP), intent(in)      :: c(:) !< (n)
    integer, intent(in)       :: k
    real(DP), intent(in)      :: x(:) !< (m)
    real(DP), intent(out)     :: y(:) !< (m)
    integer, intent(in)       :: m
    integer, intent(out)      :: ier

    !  ..array arguments..

    !  ..local scalars..
    integer :: i,j,k1,l,ll,l1,nk1
    real(DP) :: arg,sp,tb,te
    !  ..local array..
    real(DP) :: h(6)

    PUSH_SUB(splev)

    !  ..
    !  before starting computations a data check is made. if the input data
    !  are invalid control is immediately repassed to the calling program.
    ier = 10
    if(m-1 < 0) then
      POP_SUB(splev)
      return
    endif
    do i=2,m
      if(x(i) < x(i-1)) then
        POP_SUB(splev)
        return
      endif
    end do
    ier = 0
    !  fetch tb and te, the boundaries of the approximation interval.
    k1 = k+1
    nk1 = n-k1
    tb = t(k1)
    te = t(nk1+1)
    l = k1
    l1 = l+1
    !  main loop for the different points.
    do  i=1,m
    !  fetch a new x-value arg.
      arg = x(i)
      if(arg < tb) arg = tb
      if(arg > te) arg = te
    !  search for knot interval t(l) <= arg < t(l+1)
      40    if(arg < t(l1) .or. l == nk1) go to 50
      l = l1
      l1 = l+1
      go to 40
    !  evaluate the non-zero b-splines at arg.
      50    call fpbspl(t,n,k,arg,l,h)
    !  find the value of s(x) at x=arg.
      sp = 0.
      ll = l-k1
      do  j=1,k1
        ll = ll+1
        sp = sp+c(ll)*h(j)
      end do
      y(i) = sp
    end do

    POP_SUB(splev)
    return
  end subroutine splev
  
  !> Similar to splines_m::splev, but modifies the energy by using the scissors operators
  !!   described by tck.
  !!   \param tck b-spline coefficient for scissors operators, in eV.
  !!   \param E LDA energy (will be overwritten by the corrected one).
  !!   \param ryd optional parameter that tells whether E is in Ryd.
  !!     Default is <tt>.true.</tt>. Note: tck is always assumed to be in eV.
  !!
  !!   \sa typedefs_m::spline_tck
  subroutine splev_shift(tck, E, in_ryd)
    type(spline_tck), intent(in)  :: tck
    real(DP), intent(inout)       :: E
    logical, optional, intent(in) :: in_ryd
  
    real(DP) :: x_arr(1)
    real(DP) :: y_arr(1)
    integer :: m, ier
    logical :: conv_units

    PUSH_SUB(splev_shift)

    conv_units = .true.
    if (present(in_ryd)) then
      if (.not. in_ryd) conv_units = .false.
    endif

    m=1 
    if (conv_units) then
      x_arr(1) = E*ryd
    else
      x_arr(1) = E
    endif
    call splev(tck%t,tck%n,tck%c,tck%k,x_arr,y_arr,m,ier)
    if (ier==0) then
      if (conv_units) then
        E = E + y_arr(1)/ryd
      else
        E = E + y_arr(1)
      endif
    else
      write(0,'(a)') 'WARNING: Error doing spline interpolation. Ignoring correction!'
    endif

    POP_SUB(splev_shift)
  end subroutine splev_shift

end module splines_m
