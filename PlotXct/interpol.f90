!==============================================================================
!
! Routines:
!
! (1) interpol()        Originally By MLT       Last Modified 6/2008 FJR
!
!     Computes the density at any point r by linear interpolation
!     inside the eight vertices of the surrounding cube
!
!==============================================================================

#include "f_defs.h"

module interpol_m

  use global_m
  implicit none

  private

  public :: &
    interpol

contains

subroutine interpol(rhole,nfft,wfnvint,fftboxv)
  real(DP), intent(in) :: rhole(3)
  integer, intent(in) :: nfft(3)
  complex(DPC), intent(out) :: wfnvint
  complex(DPC), intent(in) :: fftboxv(nfft(1),nfft(2),nfft(3))
  
  integer :: ir(3),pr(3)
  real(DP) :: rr(3),x(3),shift

!--------------------------
! Shift to [0,1) cell

  PUSH_SUB(interpol)
  
  shift = anint( 10.d0 + maxval(rhole) )
  rr(1:3) = rhole(1:3) + shift - aint( rhole(1:3) + shift )
  
! Define the surrouding cube
!   ir .le. rr*nfft .lt. pr ; 0 .lt. ir,pr .le. nfft

  ir(1:3) = int(rr(1:3) * nfft(1:3))
  
  pr(1:3) = mod(ir(1:3) + 1, nfft(1:3))
  
  x(1:3) = rr(1:3) * nfft(1:3) - ir(1:3)
  
! Periodic boundary conditions

  if(ir(1) == 0) ir(1) = nfft(1)
  if(ir(2) == 0) ir(2) = nfft(2)
  if(ir(3) == 0) ir(3) = nfft(3)
  
  if(pr(1) == 0) pr(1) = nfft(1)
  if(pr(2) == 0) pr(2) = nfft(2)
  if(pr(3) == 0) pr(3) = nfft(3)
  
! Interpolate

  wfnvint=0.0d0
  wfnvint=wfnvint+fftboxv(ir(1),ir(2),ir(3))*(1.0-x(1))*(1.0-x(2))*(1.0-x(3))
  wfnvint=wfnvint+fftboxv(pr(1),ir(2),ir(3))*x(1)*(1.0-x(2))*(1.0-x(3))
  wfnvint=wfnvint+fftboxv(ir(1),pr(2),ir(3))*(1.0-x(1))*x(2)*(1.0-x(3))
  wfnvint=wfnvint+fftboxv(ir(1),ir(2),pr(3))*(1.0-x(1))*(1.0-x(2))*x(3)
  wfnvint=wfnvint+fftboxv(pr(1),pr(2),ir(3))*x(1)*x(2)*(1.0-x(3))
  wfnvint=wfnvint+fftboxv(ir(1),pr(2),pr(3))*(1.0-x(1))*x(2)*x(3)
  wfnvint=wfnvint+fftboxv(pr(1),ir(2),pr(3))*x(1)*(1.0-x(2))*x(3)
  wfnvint=wfnvint+fftboxv(pr(1),pr(2),pr(3))*x(1)*x(2)*x(3)
  
  POP_SUB(interpol)
  
  return
end subroutine interpol

end module interpol_m
