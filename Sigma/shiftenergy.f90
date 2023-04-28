!===========================================================================
!
! Routines()
!
! (1) shiftenergy()     Originally by ?         Last Edited: Apr/2016 (FHJ)
!
!     Computes and symmetrizes the quasiparticle spectrum, and
!     perform linear interpolation/extrapolation to solve Dyson`s equation.
!     Note that, since we may use more than two frequency points to evaluate
!     Sigma(omega), our "eqp1" in general goes beyond the linear correction
!     scheme in Eq. (37) of Hybertsen & Louie PRB.
!
!===========================================================================

#include "f_defs.h"

module shiftenergy_m

  use global_m
  implicit none

  private

  public :: &
    shiftenergy

contains

subroutine shiftenergy(sig, wfnk, alda, asx, ach, achcor, ach_n1, achcor_n1, &
  ax, efsto, asig, enew, zrenorm, nfreqgpp, ncore_excl, akih) ! ZL: akih attached as optional
  type (siginfo), intent(in) :: sig
  type (wfnkstates), intent(in) :: wfnk
  integer, intent(in) :: nfreqgpp
  SCALAR, intent(inout) :: & 
    alda(sig%ndiag+sig%noffdiag,sig%nspin), & !< vxc
    asx(nfreqgpp,sig%ndiag+sig%noffdiag,sig%nspin), & !< sx
    ach(nfreqgpp,sig%ndiag+sig%noffdiag,sig%nspin)    !< ch
  complex(DPC), intent(inout) :: achcor(sig%ndiag+sig%noffdiag,sig%nspin) !< static remainder
  SCALAR, intent(inout) :: ach_n1(sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin)
  SCALAR, intent(inout) :: achcor_n1(sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin)
  SCALAR, intent(inout) :: ax(sig%ndiag+sig%noffdiag,sig%nspin) !< x
  real(DP), intent(out) :: efsto(sig%ndiag,sig%nspin) !< eqp0
  SCALAR, intent(inout) :: asig(sig%ndiag+sig%noffdiag,sig%nspin) !< sig
  real(DP), intent(out) :: enew(sig%ndiag,sig%nspin)  !< eqp1
  real(DP), intent(out) :: zrenorm(sig%ndiag,sig%nspin) !< Znk
  integer, intent(in) :: ncore_excl !< number of core states excluded  
  SCALAR, intent(inout), optional :: akih(sig%ndiag+sig%noffdiag,sig%nspin) !< kih

  integer :: ii,jj,istart,istop,nsubs,ispin,iwlda0,iwlda1,iwlda2,iw
  integer :: ndeg(sig%ndiag)
  real(DP) :: fact,dek,dele,diffmin,delta_omega,diff,e_lk,freq0
  SCALAR :: aldai,axi, akihi ! ZL: add akihi
  SCALAR, dimension(nfreqgpp) :: asigi, asxi, achi
  SCALAR, dimension(sig%ntband) :: ach_n1i, achcor_n1i
  complex(DPC) :: achcori
  
  PUSH_SUB(shiftenergy)

  ! ZL: check optional akih
  if(sig%use_kihdat) then
    if (.not.present(akih)) call die('In shiftenergy, to use KIH, pass akih into this subroutine.')
  endif
  
  do ispin=1,sig%nspin
    
    nsubs = 1
    ndeg(nsubs) = 1
    ! This loop is used to find degeneracy
    do ii=2,sig%ndiag
! DVF : ncore_excl has to be substracted here because wfnk%elda is defined in the 
! read_wavefunction subroutine in input.f90 to be referenced to the case with 
! no core states. Same with wfnk%el/elda below. 
      dek = wfnk%elda(sig%diag(ii)-ncore_excl,ispin) - wfnk%elda(sig%diag(ii-1)-ncore_excl,ispin)
      if (abs(dek)<sig%tol .and. sig%symmetrize) then
        ! Band ii is degenerate to ii-1
        ndeg(nsubs) = ndeg(nsubs) + 1
      else
        ! Start a new degenerate subspace
        nsubs = nsubs + 1
        ndeg(nsubs) = 1
      endif
    enddo
    
    istop = 0
    ! loop degeneracy
    do ii=1,nsubs
      istart = istop + 1
      istop = istart + ndeg(ii) - 1
      aldai = ZERO
      akihi = ZERO ! ZL: for KIH
      axi = ZERO
      asxi = ZERO
      achi = ZERO
      achcori = (0.0d0, 0.0d0)
      ach_n1i(:) = ZERO
      achcor_n1i(:) = ZERO
      ! start real loop over bands
      do jj=istart,istop
        ! average VXC in the same degen group: combine first
        aldai = aldai + alda(jj,ispin)
        ! "i" is intermediate steps in Ry
        ! ZL: do the same for KIH as VXC
        akihi = akihi + akih(jj,ispin)
        ! Note: while using Vxc or KIH, the other array would just be zeros
        axi = axi + ax(jj,ispin)
        asxi(:) = asxi(:) + asx(:,jj,ispin)
        achi(:) = achi(:) + ach(:,jj,ispin)
        achcori = achcori + achcor(jj,ispin)
        ach_n1i(:) = ach_n1i(:) + ach_n1(:,jj,ispin)
        achcor_n1i(:) = achcor_n1i(:) + achcor_n1(:,jj,ispin)
      enddo

      fact = ryd / dble(ndeg(ii))
      do jj=istart,istop
        ! then AVERAGE
        ! common quantities in eV
        alda(jj,ispin) = aldai * fact
        ! ZL: do the same for KIH
        akih(jj,ispin) = akihi * fact
        ax(jj,ispin) = axi * fact
        asx(:,jj,ispin) = asxi(:) * fact
        ach(:,jj,ispin) = achi(:) * fact
        achcor(jj,ispin) = achcori * fact
        ach_n1(:,jj,ispin) = ach_n1i(:) * fact
        achcor_n1(:,jj,ispin) = achcor_n1i(:) * fact
        ! intermediate the quantity, sigma in Ry, unnecessary i
        asigi(:) = ax(jj,ispin) + asx(:,jj,ispin) + ach(:,jj,ispin)


        ! FHJ: Determine the expansion points for the linear interpolation: iwlda0
        ! is where we evaluate eqp0, and iwlda1, iwlda2 are the expansion points
        ! we used in the finite difference linearization (ie, used to compute the
        ! derivatives)
        if (sig%fdf==-3) then
          ! FHJ: Find iw closest to e_lk
          iwlda1 = 1
          iwlda2 = 2
          diffmin = INF
          e_lk = wfnk%ek(sig%diag(jj)-ncore_excl,ispin)
          ! FHJ: Figure out starting frequency for freq. grid
          if (sig%freq_grid_shift<2) then
            freq0 = sig%freqevalmin
          else
            freq0 = e_lk - sig%freqevalstep*(sig%nfreqeval-1)/2
          endif
          do iw=1,sig%nfreqeval
            diff = abs(freq0 + (iw-1)*sig%freqevalstep - e_lk)
            if (diff .lt. diffmin) then
              diffmin=diff
              iwlda2=iwlda1
              iwlda1=iw
            endif
          enddo
          iwlda0 = iwlda1
          delta_omega = (iwlda2 - iwlda1)*sig%freqevalstep
        else
          !iwlda0 is always 2 here
          iwlda0 = 2
          iwlda1 = 2
          iwlda2 = 2
          select case (sig%fdf)
            case (-1)
              iwlda1 = 1
              iwlda2 = 2
            case (0)
              iwlda1 = 1
              iwlda2 = 3
            case (1,2)
              iwlda1 = 2
              iwlda2 = 3
          end select
          delta_omega = (iwlda2 - iwlda1)*sig%dw
        endif
        asig(jj,ispin) = asigi(iwlda0)

!FHJ: Perform linear interpolation/extrapolation for all finite difference
!     schemes at once. enew is the interpolated off-shell answer = "eqp1",
!     and efsto is the on-shell answer = "eqp0".

        ! ZL: Vxc vs. KIH
        if(.not. sig%use_kihdat) then
          efsto(jj,ispin) = wfnk%elda(sig%diag(jj)-ncore_excl,ispin) - &
            alda(jj,ispin) + asig(jj,ispin)
        else
          efsto(jj,ispin) = akih(jj,ispin) + asig(jj,ispin)
        endif
        dele = efsto(jj,ispin) - wfnk%ek(sig%diag(jj)-ncore_excl,ispin)

        if (iwlda1/=iwlda2) then
          enew(jj,ispin) = efsto(jj,ispin) + &
            (asigi(iwlda2)-asigi(iwlda1)) / (delta_omega - asigi(iwlda2) + asigi(iwlda1)) * dele
          zrenorm(jj, ispin) = 1d0 / (1d0 - (asigi(iwlda2) - asigi(iwlda1))/delta_omega)
        else
          enew(jj,ispin) = efsto(jj,ispin)
          zrenorm(jj, ispin) = 1d0
        endif

      enddo ! jj
    enddo ! ii
  enddo ! ispin

  POP_SUB(shiftenergy)
  
  return
end subroutine shiftenergy

end module shiftenergy_m
