!==============================================================================
!
! Routines:
!
! (1) absh()    Originally By MLT               Last Modified 7/1/2008 (JRD)
!
!     Calculates absorption spectra in haydock code.
!
!     input: mmts type (has the Haydock coefficients, a_n b_n)
!            nspin    number of spins
!            eta      energy resolution
!
!     output: file "absorption_eh.dat"
!
!     Calculate the dielectric function (real and imaginary parts),
!     from coefficients mmts%an, mmts%bn
!
!     omega = frequency, given in eV
!     eta = energy broadening, given in eV
!
!==============================================================================

#include "f_defs.h"

module absh_m

  use global_m
  implicit none

  private

  public :: &
    absh

contains

subroutine absh(mmts,nspin,nspinor,eta)
  type (mmtsinfo), intent(in) :: mmts
  integer, intent(in) :: nspin, nspinor
  real(DP), intent(in) :: eta

  integer :: ih,iw,nwstep
  real(DP) :: w0,wfin,omega,pref,eps1,eps2
  complex(DPC) :: sum,zo

  pref = 16.d0 * PI_D**2 /(mmts%vol * dble(nspin) * dble(nspinor))

  PUSH_SUB(absh)

  w0 =  0.d0
  wfin = 20.d0
  nwstep = 2000
  
  call open_file(10,file='absorption_eh.dat',form='formatted',status='replace')

  write(10,*) "# Column 1: omega"
  write(10,*) "# Column 2: eps2(omega)"
  write(10,*) "# Column 3: eps1(omega)"

  do iw=0,nwstep
    eps2 = 0.d0
    eps1 = 0.d0
    omega = w0 + (wfin - w0) * dble(iw) / dble(nwstep)

!------------------------------------
! Absorption contribution

    zo = CMPLX(omega,eta)/RYD

! Assume that all coefficients above mmts%nmax are constant...
!
!        sum = sqrt( ((zo - mmts%an(mmts%nmax) )/2.d0)**2 - mmts%bn(mmts%nmax) )
!
! or zero... (it doesn`t matter if mmts%nmax is sufficiently high)

    sum = CMPLX(0.d0,0.d0)
    do ih= mmts%nmax, 2, -1
      sum = mmts%bn(ih-1)/( zo - mmts%an(ih) - sum)
    enddo
    sum = CMPLX(mmts%norm / PI_D, 0.d0)/( zo - mmts%an(1) - sum)
    eps2 = -pref*IMAG(sum)
    eps1 = -pref*dble(sum) + 1.d0
    
!-------------------------------------
! Emission contribution

    zo = CMPLX(-omega,eta)/RYD

! Assume that all coefficients above mmts%nmax are constant...
!
!         sum = sqrt( ((zo - mmts%an(mmts%nmax) )/2.d0)**2 -
!     >                       mmts%bn(mmts%nmax) )
!
! or zero...

    sum = CMPLX(0.d0,0.d0)
    do ih= mmts%nmax, 2, -1
      sum = mmts%bn(ih-1)/( zo - mmts%an(ih) - sum)
    enddo
    sum = CMPLX(mmts%norm / PI_D, 0.d0)/( zo - mmts%an(1) - sum)
    eps2 = pref*IMAG(sum) + eps2
    eps1 = -pref*dble(sum) + eps1
    write(10,100) omega,eps2,eps1
  enddo
  
100 format(3f16.9)
  
  call close_file(10)
  
  POP_SUB(absh)
  
  return
end subroutine absh

end module absh_m
