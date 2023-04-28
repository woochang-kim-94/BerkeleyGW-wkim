!===============================================================================
!
! Routines:
!
! (1) absp0()           Originally By MLT       Last Modified 6/5/2008 (JRD)
!
!     input: eqp             eqpinfo type
!            xct             xctinfo type
!            s0              velocity or momentum matrix elements
!            vol             crystal volume (cell volume times # of k-points)
!            omega_plasma    plasma frequency, given in Ryd
!            flag            flags type
!
!     output: file "absorption_noeh.dat"
!
!     Calculate the absorption without electron-hole interaction,
!     using eq. (27) of Rohlfing & Louie PRB 62 4927 (2000).
!     Absorption and emission contributions are calculated separately.
!     The causal, retarded dielectric function is assumed.
!
!     For the time-ordered dielectric function, both real and imaginary
!     parts are even functions of frequency. For the causal, retarded
!     dielectric function, the real part is an even function and the
!     imaginary part is an odd function of frequency.
!
!     NOTE: I am not sure if the numerical factor in eqs. (26) and (27)
!     is right. There should be a 4*Pi^2 instead of 16*Pi. (Murilo)
!
!     REPLY: There was a typo in eqs. (26) and (27), the correct factor
!     is 16*Pi^2. The code, however, was correct, an additional factor
!     of Pi was hidden in the delta function (gaussian and lorentzian).
!     I put it in the prefactor. There is a factor of 4*Pi^2 in eq.
!     (6.48) of Yu & Cardona "Fundamentals of Semiconductors" because
!     they use the real representation for electromagnetic waves.
!     Using the complex representation and including spin,
!     their factor 4*Pi^2 becomes our 16*Pi^2. (gsm)
!
!     RE-REPLY: There is another typo in eqs. (26) and (27). Including spin,
!     the prefactor should be 8*Pi^2*e^2, where e^2 = 2 because we use 
!     Rydberg units. If we follow RMP 74, 2 (2002), eqn (4.12):
!       eps(q, w) = 1 - v(q) \sum |<n1|e^{-iq.r}|n2>*A|^2 / (w - E + i*eta)
!     But v(q) = 4*Pi*e^2 / q^2 and Im[1/(x + i*eta)] ~ pi * delta(x).
!     If we multiply by the spin degeneracy, we get the prefactor 8*Pi^2*e^2.
!     Another nice reference is PRB 73, 045112 (2006), eqn (13).
!     The bottom line is that the code was right, but the extra factor of two
!     comes from the e^2 and not from the real/complex representation. (FHJ)
!
!     Each delta function is replaced by a Lorentzian, Gaussian or Voigt peak:
!     delta(x) -> (xct%eta/pi)/(x^2+xct%eta^2)
!     delta(x) -> exp(-x^2/(2*xct%eta^2))/(sqrt(2*pi)*xct%eta)
!     delta(x) -> Voigt(x,xct%sigma,xct%gamma)
!
!     omega = frequency, given in eV
!     xct%eta,xct%sigma,xct%gamma = energy broadening, given in eV
!
!===============================================================================

#include "f_defs.h"

module absp0_m

  use absp_io_m
  use global_m
  use misc_m
  implicit none

  private

  public :: &
    absp0

contains

subroutine absp0(eqp,xct,s0,vol,omega_plasma,flag,ipol)
  type (eqpinfo), intent(in) :: eqp
  type (xctinfo), intent(in) :: xct
  SCALAR, intent(in) :: s0(xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin)
  real(DP), intent(in) :: vol, omega_plasma
  type (flags), intent(in) :: flag
  integer, intent(in) :: ipol

  integer :: ic,iv,ik,ikcvs,is,iemax,iw,nwstep
  real(DP) :: emin,emax,eps1,eps2,dos
  real(DP) :: omega,omegalda,fac1,fac2,sum1,sum2,pref,fac,domega
  real(DP), allocatable :: eps1_w(:), eps2_w(:), dos_w(:)
  integer :: neig
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
  
  PUSH_SUB(absp0)
  
  pref = 16.d0 * PI_D**2 / (vol * dble(xct%nspin*xct%nspinor))
  neig = xct%nkpt_fi * xct%ncb_fi * xct%nvb_fi * xct%nspin

!----------------------------
! Check the f-sum rule using DFT and DFT+GW eigenvalues.
! The sum rule differs because of nonlocal contributions
! in the GW corrections.
! Exact value of the sum rule (nonlocal contributions neglected):
! sum1 = (pi / 2.d0) * (plasma frequency)^2

  emin = 1.d10
  emax = 0.d0
  sum1 = 0.d0
  sum2 = 0.d0
  do ik=1,xct%nkpt_fi
    do ic=1,xct%ncb_fi
      do iv=1,xct%nvb_fi
        do is=1,xct%nspin
          ikcvs = bse_index(ik, ic, iv, is, xct)
          if (xct%qflag.ne.2) then
            omega = eqp%ecqp(ic,ik,is) - eqp%evqp(iv,ik,is)
            omegalda = eqp%eclda(ic,ik,is) - eqp%evlda(iv,ik,is)
          else
            if (xct%indexq_fi(ik).eq.0 .and. xct%patched_sampling) cycle
            omega = eqp%ecqp(ic,ik,is) - eqp%evqp(iv,xct%indexq_fi(ik),is)
            omegalda = eqp%eclda(ic,ik,is) - eqp%evlda(iv,xct%indexq_fi(ik),is)
          endif
          if (omega.lt.emin) emin = omega
          if (omega.gt.emax) emax = omega

          sum1 = sum1 + omega * abs(s0(ikcvs))**2
          sum2 = sum2 + omegalda * abs(s0(ikcvs))**2

!            ediff = eqp%eclda(ic,ik,is) - eqp%evlda(iv,ik,is)
!            sum1 = sum1 + omega * abs(s0(ikcvs))**2 / ediff**2
!            ediff = eqp%eclda(ic,ik,is) - eqp%evlda(iv,ik,is)
!            sum2 = sum2 + omegalda * abs(s0(ikcvs))**2 / ediff**2
        enddo
      enddo
    enddo
  enddo
  sum1 = sum1 * pref * ryd**2
  sum2 = sum2 * pref * ryd**2
  if (ipol==1) then
    write(6,*)
    write(6,*) 'Plasma Frequency : ',omega_plasma*ryd,' eV'
    write(6,*)
  endif
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
    write(6,*) 'Sum rule (DFT) : ',sum2,' eV^2'
    write(6,*) 'Sum rule (GW)  : ',sum1,' eV^2'
    write(6,*)
  else
    sum1 = sum1 / (0.5d0 * PI_D * omega_plasma**2 * ryd**2)
    sum2 = sum2 / (0.5d0 * PI_D * omega_plasma**2 * ryd**2)
    write(6,*) 'Sum rule (DFT) : ',sum2
    write(6,*) 'Sum rule (GW)  : ',sum1
    write(6,*)
  endif
  !emin = emin * ryd
  emin = 0.d0
  emax = emax * ryd
  iemax = int(emax + 10.d0 * xct%eta) + 1
  emax = dble(iemax)
  domega = xct%delta_frequency
  nwstep = int(iemax / domega)
 
  call write_eigenvalues_noeh(xct,neig,vol,eqp,s0,ipol)

  if (xct%npol==1) then
    fname = 'absorption_noeh.dat'
  else
    fname = 'absorption_'//suffix(ipol)//'_noeh.dat'
  endif
  call open_file(10,file=trim(fname),form='formatted',status='replace')

  write(10,*) "# Column 1: omega"
  write(10,*) "# Column 2: eps2(omega)"
  write(10,*) "# Column 3: eps1(omega)"
  write(10,*) "# Column 4: JDOS(omega)"

  SAFE_ALLOCATE(eps1_w, (nwstep+1))
  SAFE_ALLOCATE(eps2_w, (nwstep+1))
  SAFE_ALLOCATE(dos_w, (nwstep+1))

  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP PRIVATE(iw,eps1,eps2,dos,omega,sum1,sum2,ik,ic,iv,is,ikcvs,fac,fac1,fac2)
  do iw = 0, nwstep
    eps1 = 0.d0
    eps2 = 0.d0
    dos = 0.d0
    omega = emin + (emax - emin) * dble(iw) / dble(nwstep)

!----------------------------
! Absorption contribution

    sum1 = 0.d0
    sum2 = 0.d0
    ikcvs = 0
    do ik=1,xct%nkpt_fi
      do ic=1,xct%ncb_fi
        do iv=1,xct%nvb_fi
          do is=1,xct%nspin
            ikcvs = ikcvs + 1
            if (xct%qflag.ne.2) then
              fac = omega / ryd - eqp%ecqp(ic,ik,is) + eqp%evqp(iv,ik,is)
            else
              if (xct%indexq_fi(ik).eq.0 .and. xct%patched_sampling) cycle
              fac = omega / ryd - eqp%ecqp(ic,ik,is) + eqp%evqp(iv,xct%indexq_fi(ik),is)
            endif
!          ediff = eqp%eclda(ic,ik,is) - eqp%evlda(iv,ik,is)
            fac1 = (-fac / PI_D) / (fac**2 + (xct%eta / ryd)**2)
            sum1 = sum1 + abs(s0(ikcvs))**2 * fac1

!            sum1 = sum1 + abs(s0(ikcvs))**2 * fac1 / ediff**2
            
            if (flag%lor.eq.0) then
              fac2 = (xct%eta / ryd / PI_D) / (fac**2 + (xct%eta / ryd)**2)
            else if (flag%lor.eq.1) then
              fac2 = exp(-fac**2 / (2.d0 * (xct%eta / ryd)**2)) / (sqrt(2.d0 * PI_D) * xct%eta / ryd)
            else
              fac2 = voigt(fac, xct%sigma, xct%gamma)
            endif
            
            sum2 = sum2 + abs(s0(ikcvs))**2 * fac2
!            sum2 = sum2 + abs(s0(ikcvs))**2 * fac2 / ediff**2
            dos = dos + fac2
          enddo
        enddo
      enddo
    enddo
    eps1 = 1.d0 + pref * sum1
    eps2 = pref * sum2
    dos = dos / (ryd * dble(neig))
    
!--------------------------
! Emission contribution
! eps2 gets negative contribution, so it vanishes at zero frequency

    sum1 = 0.d0
    sum2 = 0.d0
    ikcvs = 0
    do ik=1,xct%nkpt_fi
      do ic=1,xct%ncb_fi
        do iv=1,xct%nvb_fi
          do is=1,xct%nspin
            ikcvs = ikcvs + 1
            if (xct%qflag.ne.2) then
              fac = -omega / ryd - eqp%ecqp(ic,ik,is) + eqp%evqp(iv,ik,is)
            else
              if (xct%patched_sampling .and. xct%indexq_fi(ik).eq.0) cycle
              fac = -omega / ryd - eqp%ecqp(ic,ik,is) + eqp%evqp(iv,xct%indexq_fi(ik),is)
            endif
!          ediff = eqp%eclda(ic,ik,is) - eqp%evlda(iv,ik,is)
            fac1 = (-fac / PI_D) / (fac**2 + (xct%eta/ryd)**2)
            sum1 = sum1 + abs(s0(ikcvs))**2 * fac1
!            sum1 = sum1 + abs(s0(ikcvs))**2 * fac1 / ediff**2
            
            if (flag%lor .eq. 0) then
              fac2 = -(xct%eta / ryd / PI_D) / (fac**2 + (xct%eta / ryd)**2)
            else if (flag%lor.eq.1) then
              fac2 = -exp(-fac**2 / (2.d0 * (xct%eta / ryd)**2)) / (sqrt(2.d0 * PI_D) * xct%eta / ryd)
            else
              fac2 = -voigt(fac, xct%sigma, xct%gamma)
            endif
            
            sum2 = sum2 + abs(s0(ikcvs))**2 * fac2
!            sum2 = sum2 + abs(s0(ikcvs))**2 * fac2 / ediff**2
          enddo
        enddo
      enddo
    enddo
!        if (flag%opr.eq.1) then
!          sum1 = sum1 / (omega / ryd)**2
!          sum2 = sum2 / (omega / ryd)**2
!        endif
    eps1 = eps1 + pref * sum1
    eps2 = eps2 + pref * sum2
    eps1_w(iw+1) = eps1
    eps2_w(iw+1) = eps2
    dos_w(iw+1) = dos
  enddo
  !$OMP END PARALLEL DO

  do iw = 0, nwstep
    omega = emin + (emax - emin) * dble(iw) / dble(nwstep)
    write(10,'(4f16.9)') omega, eps2_w(iw+1), eps1_w(iw+1), dos_w(iw+1)
  enddo
  SAFE_DEALLOCATE(eps1_w)
  SAFE_DEALLOCATE(eps2_w)
  SAFE_DEALLOCATE(dos_w)
  
  call close_file(10)
  
  POP_SUB(absp0)
  
  return
  
end subroutine absp0

end module absp0_m
