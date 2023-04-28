!===============================================================================
!
! Routines:
!
! (1) absp()            Originally By: MLT      Last Modified   6/5/2008 (JRD)
!
!     input: eta,sigma,gamma energy resolution, given in eV
!            neig            number of excitonic states
!            cs              transition matrix elements, given in au
!            en              energy eigenvalue of each state, given in eV
!            vol             crystal volume (cell volume times # of k-points)
!            omega_plasma    plasma frequency, given in Ryd
!            flag            flags type
!
!     output: file "absorption_eh.dat"
!
!     Calculate the absorption with electron-hole interaction,
!     using eq. (26) of Rohlfing & Louie PRB 62 4927 (2000),
!     and the density of exciton states, given in Ryd^-1.
!
!     epsilon_1 always uses Lorentzian broadening.
!     epsilon_2 uses Lorentzian, Gaussian or Voigt, depending on flag%lor
!     delta(x) -> (eta/pi)/(x^2+eta^2)
!     delta(x) -> exp(-x^2/(2*eta^2))/(sqrt(2*pi)*eta)
!     delta(x) -> Voigt(x,sigma,gamma)
!
!     See the comment about the prefactor in file absp0.f90
!
!     omega = frequency, given in eV
!     eta,sigma,gamma = energy broadening, given in eV
!
!===============================================================================

#include "f_defs.h"

module absp_m

  use global_m
  use misc_m
  implicit none

  private

  public :: &
    absp

contains

subroutine absp(xct,neig,cs,en,vol,omega_plasma,flag,ipol)
  type(xctinfo), intent(in) :: xct
  integer, intent(in) :: neig
  real(DP), intent(in) :: cs(neig), en(neig), vol, omega_plasma
  type (flags), intent(in) :: flag
  integer, intent(in) :: ipol

  integer :: ii,iemax,iw,nwstep
  real(DP) :: emin,emax,eps1,eps2,dos
  real(DP) :: omega,fac,fac1,fac2,sum0,sum1,sum2,pref,domega
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
  
  PUSH_SUB(absp)
  
  pref = 16.d0 * PI_D**2 / (vol * dble(xct%nspinor))
  ! the factor 1/sqrt(nspin) in absp0 is handled by the normalization of the exciton here
  
!----------------------------
! Check the f-sum rule using BSE eigenvalues.
! Exact value of the sum rule:
! sum1 = (pi / 2.d0) * (plasma frequency)^2

  if (xct%tda) then
    sum1 = sum(cs(1:neig)*en(1:neig))
  else
    sum1 = sum(cs*en,en>0)
  endif
  sum1 = sum1 * pref * ryd
  if (ipol==1) then
    write(6,'()')
    write(6,*) 'Plasma Frequency : ',omega_plasma*ryd,' eV'
  endif
  write(6,'()')
  if (xct%npol==1) then
    if (flag%opr==0) then
      write(6,'(1x,a,3(f10.6,1x))') 'Polarization: ', xct%shift
    else
      write(6,'(1x,a,3(f10.6,1x))') 'Polarization: ', xct%pol
    endif
  else
    write(6,'(1x,a)') 'Polarization: '//suffix(ipol)
  endif
  if (omega_plasma.lt.TOL_Small) then
    write(6,*) 'Sum rule (BSE) : ',sum1,' eV^2'
  else
    sum1 = sum1 / (0.5d0 * PI_D * omega_plasma**2 * ryd**2)
    write(6,*) 'Sum rule (BSE) : ',sum1
  endif
  !emin = minval(en)
  emin = 0.d0
  emax = maxval(en)
  iemax = int(emax + 10.d0 * xct%eta) + 1
  emax = dble(iemax)
  domega = xct%delta_frequency
  nwstep = int(iemax / domega)


  if (xct%npol==1) then
    fname = 'absorption_eh.dat'
  else
    fname = 'absorption_'//suffix(ipol)//'_eh.dat'
  endif
  call open_file(10,file=trim(fname),form='formatted',status='replace')
  
  write(10,*) "# Column 1: omega"
  write(10,*) "# Column 2: eps2(omega)"
  write(10,*) "# Column 3: eps1(omega)"
  write(10,*) "# Column 4: JDOS(omega)"

  do iw=0,nwstep
    eps2 = 0.d0
    eps1 = 0.d0

!-----------------------------
! Absorption contribution

    omega = emin + (emax - emin) * dble(iw) / dble(nwstep)
    sum0 = 0.d0
    sum1 = 0.d0
    sum2 = 0.d0
    do ii=1,neig
      if (.not.xct%tda.and.en(ii)<0) cycle
      fac = omega - en(ii)
      fac1 = (-fac / PI_D)/(fac**2 + xct%eta**2)
      sum1 = sum1 + cs(ii) * fac1 * ryd
      
      if (flag%lor .eq. 0) then
        fac2 = (xct%eta / PI_D) / (fac**2 + xct%eta**2)
      else if (flag%lor .eq. 1) then
        fac2 = exp(-fac**2 / (2.d0 * xct%eta**2)) / (sqrt(2.d0 * PI_D) * xct%eta)
      else
        fac2 = voigt(fac, xct%sigma, xct%gamma)
      endif
      
      sum2 = sum2 + cs(ii) * fac2 * ryd
      sum0 = sum0 + fac2
    enddo
    eps1 = 1.d0 + pref*sum1
    eps2 = pref*sum2
    dos = sum0*xct%nspin/dble(neig)

!---------------------------
! Emission contribution
! These correspond to peaks in the spectra at negative omega used so
! eps2 gets negative contribution, so it vanishes at zero frequency.

    sum1 = 0.d0
    sum2 = 0.d0
    do ii=1,neig
      if (.not.xct%tda.and.en(ii)<0) cycle
      fac = -omega - en(ii)
      fac1 = (-fac / PI_D) / (fac**2 + xct%eta**2)
      sum1 = sum1 + cs(ii) * fac1 * ryd
      
      if (flag%lor .eq. 0) then
        fac2 = -(xct%eta / PI_D) / (fac**2 + xct%eta**2)
      else if (flag%lor .eq. 1) then
        fac2 = -exp(-fac**2 / (2.d0 * xct%eta**2)) / (sqrt(2.d0 * PI_D) * xct%eta)
      else
        fac2 = -voigt(fac, xct%sigma, xct%gamma)
      endif
      
      sum2 = sum2 + cs(ii) * fac2 * ryd
    enddo
    eps1 = eps1 + pref * sum1
    eps2 = eps2 + pref * sum2
    write(10,100)omega,eps2,eps1,dos
  enddo
  
  call close_file(10)
  
  POP_SUB(absp)
  
  return
  
100 format(4f16.9)
  
end subroutine absp

end module absp_m
