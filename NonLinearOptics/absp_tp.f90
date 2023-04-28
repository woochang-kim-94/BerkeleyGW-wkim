!===================================================================================
!
! Routines:
!
! (1) absp_tp() Originally By JRD               Last Modified 6/6/2008 (JRD)
!
!     input: nspin     number of spins
!          : nspinor   number of spinors
!            eta       energy resolution
!            neig      number of excitonic states
!            cs        transition matrix element of each state
!            en        energy eigenvalue of each state
!            vol       crystal volume (cell volume times # of k-points)
!
!     output: file "absorption_nl.dat"
!
!     Calculate the nonlinear absorption coefficient as a function of photon 
!     energy in the presence of N photons per unit volume.  
!     See Bassani "Electronic states and optical transitions in solids" pg. 165
!
!     Each delta function is replaced by a Gaussian or Lorentzian peak,
!     delta(x) -> exp(-x^2/(2*D^2))/(sqrt(2*pi)*D), D = eta
!     delta(x) -> (D/pi)/(x^2+D^2), D = eta
!
!     The prefactor is 64 pi^3 / c * omega^2 (where a second factor of 2 comes from the 
!     order at which the photons are absorbed - and another factor of 2 comes 
!     from real to complex representation). We divide by nspin so that spin
!     polarized and un-polarized calculations give the same answer.
!
!     omega = frequency, given in eV
!     eta = energy broadening, given in eV
!
!===================================================================================

#include "f_defs.h"

module absp_tp_m

  use global_m
  implicit none

  private

  public :: &
    absp_tp

contains

subroutine absp_tp(nspin,nspinor,eta,neig,cs,en,vol, &
  nmat,flag)
  type (flags) :: flag
  integer :: nspin,nspinor,neig,nmat
  real(DP) :: eta,vol,en(nmat)
  SCALAR :: tmp1,tmp2
  real(DP) :: cs(neig)

!----------------------------------
! Local variables

  integer :: ii,iemax,iw,nwstep
  real(DP) :: emin,emax,eps2
  real(DP) :: omega,fac,fac2,sum2,pref
  
  PUSH_SUB(absp_tp)
  
  pref = 64.d0 * PI_D**3 / ( vol * dble(nspin * nspinor) * LIGHTSPEED) 
  
  emin= minval(en)
!      emax = maxval(en)
  emax = en(neig)
!      emin = max(emin - 10.d0 * eta, 0.d0)
!      emax = emax + 10.d0 * eta
!      nwstep = 1000
  emin = 0.d0
  iemax = int(emax + 10.d0 * eta) + 1
  emax = dble(iemax)
  nwstep = 100 * iemax

  call open_file(10,file='absorption_nl.dat',form='formatted',status='replace')

  write(10,*) "# Column 1: omega"
  write(10,*) "# Column 2: alpha(omega)"
  write(10,*) " "
  write(10,*) " "

  do iw=0,nwstep
    
    eps2 = 0.d0
    
    omega = emin + (emax - emin) * dble(iw) / dble(nwstep)
    sum2 = 0.d0

!------------------------------
! Loop over Final states

    do ii=1,neig
      tmp1 = 0.d0
      tmp2 = 0.d0

!------------------------------
! Absorption contribution

      fac = omega - en(ii)
      if (flag%lor .eq. 0) then
        fac2 = exp(-fac**2 / (2.d0 * eta**2))
        fac2 = fac2 / (sqrt(2.d0 * PI_D) * eta)
      else
        fac2 = (eta / PI_D) / (fac**2 + eta**2)
      endif

!----------------------------
! Emission contribution
! JRD: This should be included to give 0 at zero frequency....
!
!          fac = -omega - en(ii)
!          if (flag%lor .eq. 0) then
!            fac2 = -exp(-fac**2 / (2.d0 * eta**2))
!            fac2 = fac2 / (sqrt(2.d0 * PI_D) * eta)
!          else
!            fac2 = -(eta / PI_D) / (fac**2 + eta**2)
!          endif

! JRD: omega still power of 2.  See Bassani 5-34

      sum2 = sum2 + cs(ii)*fac2*ryd / &
        ((en(ii)/ryd)**2)
    enddo

! This factor in front may not be correct

    eps2 = pref * sum2
    write(10,100) omega,eps2
  enddo
  
  call close_file(10)
  
  POP_SUB(absp_tp)

  return

100 format(3f22.9)

end subroutine absp_tp

end module absp_tp_m
