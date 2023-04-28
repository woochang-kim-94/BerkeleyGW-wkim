!===================================================================================
!
! Routines:
!
! (1) absp_uf() Originally By JRD               Last Modified 6/6/2008 (JRD)
!
!     input: eta       energy resolution
!            spin      number of spins
!            spinor    number of spinors
!            neig      number of excitonic states
!            cs        transition matrix element of each state
!            en        energy eigenvalue of each state
!            vol       crystal volume (cell volume times # of k-points)
!
!     output: file "absorption_nl.dat"
!
!     Calculate the nonlinear absorption as a function of photon energy.
!
!     Each delta function is replaced by a Gaussian or Lorentzian peak,
!     delta(x) -> exp(-x^2/(2*D^2))/(sqrt(2*pi)*D), D = eta
!     delta(x) -> (D/pi)/(x^2+D^2), D = eta
!
!     See the comment about the prefactor in file BSE/absp0.f90
!
!     omega = frequency, given in eV
!     eta = energy broadening, given in eV
!
!===================================================================================

#include "f_defs.h"

module absp_uf_m

  use global_m
  implicit none

  private

  public :: &
    absp_uf

contains

subroutine absp_uf(win,eta,spin,spinor,neig,cs,en,vol,nmat,flag)
  type (flags) :: flag
  type (windowinfo) :: win 
  integer :: spin,spinor,neig,nmat
  real(DP) :: eta,vol,en(nmat,win%nstates), &
    cs(nmat,win%nstates)

!----------------------------------
! Local variables

  integer :: ii,iemax,iw,nwstep,jj
  real(DP) :: emin,emax,eps2,dos
  real(DP) :: omega,fac,fac2,sum,sum1,sum2,pref
  
  pref = 16.d0 * PI_D**2 / (vol * dble(spin*spinor))
  
  PUSH_SUB(absp_uf)
  
  emin= minval(en)
  emax = en(neig,1)
  emin = 0.d0
  iemax = int(emax + 10.d0 * eta) + 1
  emax = dble(iemax)
  nwstep = 100 * iemax
  
  call open_file(10,file='absorption_nl.dat',form='formatted',status='replace')
  
  write(10,*) "# Column 1: omega"
  write(10,*) "# Column 2: alpha(omega)"
  write(10,*) "# Column 3: JDOS "
  write(10,*) " "
  
  do iw=0,nwstep
    eps2 = 0.d0

!------------------------------
! Absorption contribution

    omega = emin + (emax - emin) * dble(iw) / dble(nwstep)
    sum = 0.d0
    sum1 = 0.d0
    sum2 = 0.d0
    do ii=1,neig
      do jj = 1, win%nstates
        fac = omega - en(ii,jj)
        if ( flag%lor .eq. 0) then
          fac2 = exp(-fac**2 / (2.d0 * eta**2))
          fac2 = fac2 / (sqrt(2.d0 * PI_D) * eta)
        else
          fac2 = (eta / PI_D) / (fac**2 + eta**2)
        endif
        sum2 = sum2 + cs(ii,jj)*fac2*ryd &
          / ((en(ii,jj)/ryd)**2)
        sum = sum + fac2 * ryd * win%cstates(jj)
      enddo
    enddo
    eps2 = pref * sum2
    dos = sum / (PI_D * dble(neig))

!----------------------------
! Emission contribution

    sum1 = 0.d0
    sum2 = 0.d0
    do ii=1,neig
      do jj = 1, win%nstates
        fac = -omega - en(ii,jj)
        if (flag%lor .eq. 0) then
          fac2 = -exp(-fac**2 / (2.d0 * eta**2))
          fac2 = fac2 / (sqrt(2.d0 * PI_D) * eta)
        else
          fac2 = -(eta / PI_D) / (fac**2 + eta**2)
        endif
        sum2 = sum2 + cs(ii,jj)*fac2*ryd*win%cstates(jj) &
          / ((en(ii,jj)/ryd)**2)
      enddo
    enddo
    eps2 = eps2 + pref * sum2
    write(10,100) omega,eps2,dos
  enddo
  
  call close_file(10)
  
  POP_SUB(absp_uf)
  
  return

100 format(3f16.9)

end subroutine absp_uf

end module absp_uf_m
