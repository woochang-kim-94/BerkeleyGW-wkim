!==================================================================================
!
! Routines:
!
! (1) absp3d()          Originally By JRD               Last Modified 6/6/2008 (JRD)
!
!     Calculate absorption vs energy vs kpoint for ultrafast.
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
!==================================================================================

#include "f_defs.h"

module absp3d_m

  use global_m
  implicit none

  private

  public :: &
    absp3d

contains

subroutine absp3d(nspin,nspinor,eta,neig,cs,en,vol,nmat, &
  nk,ikmax,keta,knx,kny,knz,k,flag)
  type (flags) :: flag
  
  integer :: nspin,nspinor,neig,nmat
  integer :: ikmax(nmat)
  real(DP) :: eta,vol,en(nmat),cs(nmat)
  integer :: ii,iemax,iw,nwstep,jj,nk,knx,kny,knz,kk,ll
  real(DP) :: emin,emax,eps1,eps2,dos,fac3,fac4
  real(DP) :: omega,fac,fac2,sum,sum1,sum2,pref
  real(DP) :: k(3,nk),keta,kr(3)
  
  pref = 16.d0 * PI_D**2 / (vol * dble(nspin*nspinor))
  
  PUSH_SUB(absp3d)
  
  emin= minval(en)
  emax = maxval(en)
!      emin = max(emin - 10.d0 * eta, 0.d0)
!      emax = emax + 10.d0 * eta
!      nwstep = 1000
  emin = 0.d0
  iemax = int(emax + 10.d0 * eta) + 1
  emax = dble(iemax)
  nwstep = 10 * iemax

  call open_file(11,file='absorp3D',form='formatted',status='replace')
!      call open_file(12,file='eps3D',form='unformatted')
!      call open_file(13,file='dos3D',form='unformatted')
!      call open_file(12,file='eps3D')
!      call open_file(13,file='dos3D')

  write(6,*) 'In absp3d!'
  write(6,*) ''

!      do ii = 1, neig
!        write(6,*) 'eig,ikmax,en',ii,ikmax(ii),en(ii)
!      end do

  write(6,*) ''
  
  do iw=0,nwstep
    do jj=1,knx
      do kk=1,kny
        do ll=1,knz
          
          eps2 = 0.d0
          eps1 = 0.d0

!----------------------------
! Absorption contribution

          omega = emin + (emax - emin) *dble(iw) / dble(nwstep)
          if (knx .eq. 1) then
            kr(1) = k(1,1)
          else
            kr(1) = dble(jj) / dble(knx) - 0.5d0
          end if
          if (kny .eq. 1) then
            kr(2) = k(2,1)
          else
            kr(2) = dble(kk) / dble(kny) - 0.5d0
          end if
          if (knz .eq. 1) then
            kr(3) = k(3,1)
          else
            kr(3) = dble(ll) / dble(knz) - 0.5d0
          end if
          
          sum = 0.d0
          sum1 = 0.d0
          sum2 = 0.d0
          do ii=1,neig
            fac = omega - en(ii)
            fac3 = (k(1,ikmax(ii)) - kr(1))**2 + &
              (k(2,ikmax(ii)) - kr(2))**2 + &
              (k(3,ikmax(ii)) - kr(3))**2
            if (flag%lor .eq. 0) then
              fac2 = exp(-fac**2 / (2.d0 * eta**2))
              fac2 = fac2 / (sqrt(2.d0 * PI_D) * eta)
              fac4 = exp(-fac3**2 / (2.d0 * keta**2))
              fac4 = fac4 / (sqrt(2.d0 * PI_D) * keta)
            else
              fac2 = (eta / PI_D) / (fac**2 + eta**2)
              fac4 = (keta / PI_D) / (fac3**2 + keta**2)
            endif
            sum2 = sum2 + cs(ii) * fac2 * fac4 * ryd
            sum = sum + fac2 * fac4 * ryd
          enddo

          eps2 = pref * sum2

          dos = sum / (PI_D * dble(neig))

!------------------------------
! Emission contribution

          sum1 = 0.d0
          sum2 = 0.d0
          do ii=1,neig
            fac = -omega - en(ii)
            fac3 = (k(1,ikmax(ii)) - kr(1))**2 + &
              (k(2,ikmax(ii)) - kr(2))**2 + &
              (k(3,ikmax(ii)) - kr(3))**2
            if (flag%lor .eq. 0) then
              fac2 = -exp(-fac**2 / (2.d0 * eta**2))
              fac2 = fac2 / (sqrt(2.d0 * PI_D) * eta)
              fac4 = exp(-fac3**2 /(2.d0 * keta**2))
              fac4 = fac4 / (sqrt(2.d0 * PI_D) * keta)
            else
              fac2 = (eta / PI_D) / (fac**2 + eta**2)
              fac4 = (keta / PI_D) / (fac3**2 + keta**2)
            endif
            
            sum2 = sum2 + cs(ii) * fac2 * fac4 * ryd
          enddo
          eps2 = eps2 + pref * sum2
          
          write(11,100) kr(1),kr(2),kr(3),omega,eps2,dos

!        eps23D(jj,kk,ll,iw)=eps2
!        dos3D(jj,kk,ll,iw)=dos

        enddo
      enddo
    enddo
    
  enddo
  
  call close_file(11)
  call close_file(12)
  call close_file(13)
  
  POP_SUB(absp3d)
  
  return
  
100 format(6f16.9)

end subroutine absp3d

end module absp3d_m
