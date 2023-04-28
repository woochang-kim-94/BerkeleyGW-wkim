!===========================================================================
!
! Routines()
!
! (1) shiftenergy_dyn() Originally by ?         Last Edited: 5/12/2008 (JRD)
!
!     Computes and symmetrizes the quasiparticle spectrum
!
!===========================================================================

#include "f_defs.h"

module shiftenergy_dyn_m

  use global_m
  implicit none

  private

  public :: &
    shiftenergy_dyn

contains

! FHJ: For some reason, we dropped the "D"s for the variable names in these
! subroutine. This is confusing.
subroutine shiftenergy_dyn(sig, wfnk, alda, asx, ach, ach_cor, ach_corb, &
  ach2, achcor, ach_n1, achcor_n1, ax, efsto, asig, enew, enew_nosr, &
  neqp1, neqp1_nosr, ikn, kp, ncore_excl, akih) ! ZL: akih attached as optional
  type (siginfo), intent(in) :: sig
  type (wfnkstates), intent(in) :: wfnk
  SCALAR, intent(inout) :: alda(sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(inout) :: &
    asx(sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin), &
    ach(sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin), &
    ach_cor(sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin), &
    ach_corb(sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin), &
    ach2(sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin), &
    achcor(sig%ndiag+sig%noffdiag,sig%nspin), &
    ach_n1(sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin)
  SCALAR, intent(inout) :: achcor_n1(sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin)
  SCALAR, intent(inout) :: ax(sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(out) :: efsto(sig%ndiag,sig%nspin)
  complex(DPC), intent(inout) :: asig(sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(out) :: enew(sig%ndiag,sig%nspin)
  complex(DPC), intent(out) :: enew_nosr(sig%ndiag,sig%nspin)
  integer, intent(out) :: neqp1(sig%ndiag,sig%nspin)
  integer, intent(out) :: neqp1_nosr(sig%ndiag,sig%nspin)
  integer, intent(in) :: ikn
  type (kpoints), intent(in) :: kp
  integer, intent(in) :: ncore_excl !< number of core states excluded  
  SCALAR, intent(inout), optional :: akih(sig%ndiag+sig%noffdiag,sig%nspin) !< kih

  logical, save :: warned=.false.
  integer :: iwlda,iw
  integer :: ii,jj,istart,istop,nl,iflag,ispin
  integer :: ndeg(sig%ndiag)
  real(DP) :: fact,dek,eval,diff,diffmin,e_lk,freq0
  real(DP) :: specsum(sig%nfreqeval)
  SCALAR :: aldai,axi, akihi  ! ZL: add akihi
  SCALAR :: achcor_n1i(sig%ntband)
  complex(DPC) :: ach_n1i(sig%ntband)
  complex(DPC) :: asxi(sig%nfreqeval), achi(sig%nfreqeval), asigt(sig%nfreqeval)
  complex(DPC) :: ach_cori(sig%nfreqeval)
  complex(DPC) :: ach_corbi(sig%nfreqeval)
  complex(DPC) :: asigt_cor(sig%nfreqeval)
  complex(DPC) :: asigt_corb(sig%nfreqeval)
  complex(DPC) :: ach2i(sig%nfreqeval), asigt2(sig%nfreqeval), achcori
  complex(DPC) :: fmin(sig%nfreqeval), fmin_nosr(sig%nfreqeval)
  real(DP) :: freqs(sig%nfreqeval)
  
! JRD: CHANGE THIS ROUTINE A----LOT!!! :) 

  PUSH_SUB(shiftenergy_dyn)

  ! ZL: check optional akih
  if(sig%use_kihdat) then
    if (.not.present(akih)) call die('In shiftenergy_dyn, to use KIH, pass akih into this subroutine.')
  endif

  do ispin=1,sig%nspin
    
    nl=1
    ndeg(nl)=1
    do ii=2,sig%ndiag
      iflag=0
! DVF : ncore_excl has to be substracted here because wfnk%ek is defined in the 
! read_wavefunction subroutine in input.f90 to be referenced to the case with 
! no core states. Same with wfnk%elda below. 
      dek = wfnk%elda(sig%diag(ii)-ncore_excl,ispin) - wfnk%elda(sig%diag(ii-1)-ncore_excl,ispin)
      if(abs(dek) .lt. sig%tol .and. sig%symmetrize) iflag=1
      if (iflag.eq.0) nl=nl+1
      if (iflag.eq.0) ndeg(nl)=1
      if (iflag.eq.1) ndeg(nl)=ndeg(nl)+1
    enddo
    
    specsum = 0D0
    
    istop = 0
    do ii=1,nl
      istart = istop + 1
      istop = istart + ndeg(ii) - 1
      aldai =0.0d0
      akihi =0.0d0  ! ZL: kih
      axi =0.0d0
      asxi(:) = 0.0d0
      achi(:) = 0.0d0
      ach_cori(:) = 0.0d0
      ach_corbi(:) = 0.0d0
      ach2i(:) = 0.0d0
      achcori = 0.0d0
      ach_n1i(:) = (0d0, 0d0)
      achcor_n1i(:) = ZERO
      do jj=istart,istop
        aldai = aldai + alda(jj,ispin)
        akihi = akihi + akih(jj,ispin) ! kih
        axi = axi + ax(jj,ispin)
        asxi(:) = asxi(:) + asx(:,jj,ispin)
        achi(:) = achi(:) + ach(:,jj,ispin)
        ach_cori(:) = ach_cori(:) + ach_cor(:,jj,ispin)
        ach_corbi(:) = ach_corbi(:) + ach_corb(:,jj,ispin)
        ach2i(:) = ach2i(:) + ach2(:,jj,ispin)
        achcori = achcori + achcor(jj,ispin)
        ach_n1i(:) = ach_n1i(:) + ach_n1(:,jj,ispin)
        achcor_n1i(:) = achcor_n1i(:) + achcor_n1(:,jj,ispin)
      enddo
      
      fact = ryd / dble(ndeg(ii))
      do jj=istart,istop
        alda(jj,ispin) = aldai * fact
        akih(jj,ispin) = akihi * fact  ! kih
        ax(jj,ispin) = axi * fact
        asx(:,jj,ispin) = asxi(:) * fact
        ach(:,jj,ispin) = achi(:) * fact
        ach_cor(:,jj,ispin) = ach_cori(:) * fact
        ach_corb(:,jj,ispin) = ach_corbi(:) * fact
        ach2(:,jj,ispin) = ach2i(:) * fact
        achcor(jj,ispin) = achcori * fact
        ach_n1(:,jj,ispin) = ach_n1i(:) * fact
        achcor_n1(:,jj,ispin) = achcor_n1i(:) * fact
        if (sig%freq_dep_method .eq. 2) then
          asigt(:) = ax(jj,ispin) + asx(:,jj,ispin) + ach(:,jj,ispin) + achcor(jj,ispin)
          asigt_cor(:) = ach(:,jj,ispin) + achcor(jj,ispin)
          asigt_corb(:) = asx(:,jj,ispin) + ach(:,jj,ispin) + achcor(jj,ispin)
          asigt2(:) = asx(:,jj,ispin)
        else
          asigt(:) = ax(jj,ispin) + asx(:,jj,ispin) + ach(:,jj,ispin) + achcor(jj,ispin)
          asigt_cor(:) = ax(jj,ispin) + ach_cor(:,jj,ispin) + achcor(jj,ispin)
          asigt_corb(:) = ax(jj,ispin) + ach_corb(:,jj,ispin) + achcor(jj,ispin)
          asigt2(:) = ax(jj,ispin) + asx(:,jj,ispin) + ach2(:,jj,ispin) + achcor(jj,ispin)
        endif

! JRD: Find iw closest to e_lk
        
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
            iwlda=iw
          endif
        enddo
        
!            write(6,*) 'ShiftEnergy - iwlda', j, e_lk, iwlda 

! gsm: asig shall not contain static remainder, it is added in write_result_dyn/write_result_dyn_hp
!            asig(jj,ispin) = asigt(iwlda) - achcor(jj,ispin)
        asig(jj,ispin) = ax(jj,ispin) + asx(iwlda,jj,ispin) + ach(iwlda,jj,ispin)

! JRD: Write out Sigma(omega)

        if (peinf%inode .eq. 0) then
          if (sig%freq_dep_method .eq. 2) then
            write(8000,2002) 
          else
            write(8000,2001) 
          endif
          do iw=1,sig%nfreqeval
            eval = freq0 + (iw-1)*sig%freqevalstep
            freqs(iw) = eval

            if (IMAG(asigt(iw))>TOL_ZERO .and. sig%freq_dep_method/=2) then
              if (.not.warned) then
                write(0,'(a,3f12.6)') 'WARNING: You have a positive imaginary Sigma at k = ', kp%rk(:,ikn)
                write(0,'(a,i6,a,f12.6,a,2f12.6)') &
                  'band calc ', jj, ', evaluation energy ', eval, ' : Im Sigma = ', IMAG(asigt(iw))
                warned = .true.
              endif
            endif
! ZL: the line below was commented out in the trunk, so do not use.  But if it is used in the future,
!     note that one should add the if-statement to differ Vxc (alda) or KIH (akih)
!
!                eval2 = wfnk%elda(sig%diag(jj),ispin) + dble(asigt(iw)) - alda(jj,ispin)
!
! CHP: spectral function is not very meaningful unless one knows the final Fermi energy
!      (i.e., after the GW correction) conserving the Fermi sphere volume. In practice,
!      one has to (1) read the calculated real part of the self energy, (2) find out
!      what is the new Fermi energy (manually), (3) subtract this value from the
!      real part of the self energies, (4) and calculate the spectral function.
!
!      Also, EQP(k,w) is not physical.  EQP is a function only of k which can
!      be obtained by solving EQP(k)=E0(k)+Re(Sig(k,EQP(k))).  This process is
!      straightforward once one knows Re(Sig(k,w)).
!
!      Printed energy is now with respect to the Fermi energy (sig%efermi).
!
!                spectral = (1D0/PI_D) * abs(IMAG(asigt(iw))) / ( (eval - eval2)**2 + IMAG(asigt(iw))**2)
!                specsum(iw) = specsum(iw) + spectral
!                write(8000,2000) kp%rk(:,ikn), ispin, sig%diag(jj), iw, eval, eval2, &
!                 dble(asigt(iw)),IMAG(asigt(iw)),spectral
            if (sig%freq_dep_method .eq. 2) then
              write(8000,3000) kp%rk(:,ikn), ispin, sig%diag(jj), iw, eval, &
                dble(asigt(iw)),IMAG(asigt(iw)),dble(asigt2(iw)),IMAG(asigt2(iw)), &
                dble(asigt_corb(iw)),IMAG(asigt_corb(iw)),dble(asigt_cor(iw)),IMAG(asigt_cor(iw))
            else
              write(8000,2000) kp%rk(:,ikn), ispin, sig%diag(jj), iw, eval, &
                dble(asigt(iw)),IMAG(asigt(iw)),dble(asigt_cor(iw)),IMAG(asigt_cor(iw))
            endif

! CHP: IM(asigt2) was obtained by using a zero energy broadening in the energy
!      denominator for the self energy evaluation.  Thus, IM(asigt2) vanishes
!      at the Fermi level, which is physically correct: the scattering rate of
!      a quasiparticle at the Fermi surface is zero.  However, (1) this routine
!      currently is meaningful only for systems having inversion symmetry
!      and (2) IM(asigt2) does not satisfy the Kramers-Kronig relation with
!      RE(asigt) for obvious reasons.

          enddo
          write(8000,*) ''
              
          ! This is Eqp0 w/o SR
          ! ZL: Vxc vs. KIH
          if(.not.sig%use_kihdat) then
            efsto(jj,ispin) = wfnk%elda(sig%diag(jj)-ncore_excl,ispin) - &
              alda(jj,ispin) + asig(jj,ispin)
          else
            efsto(jj,ispin) = akih(jj,ispin) + asig(jj,ispin)
          endif
          ! FHJ: Finding Eqp1 is equiv to finding the roots of:
          ! f(w) := e_mf + sigma(w) - vxc - w
          if (sig%freq_dep_method==2) then
            ! Using: Sigma = (X+COR) + SR = Res + Int + SR
            ! ZL: Vxc vs. KIH
            if(.not.sig%use_kihdat) then
              fmin(:) = wfnk%elda(sig%diag(jj)-ncore_excl,ispin) + asigt(:) - alda(jj,ispin) - freqs(:)
            else
              fmin(:) = asigt(:) + akih(jj,ispin) - freqs(:)
            endif
          else
            ! Using: Sigma = (X+COR) + SR. NOTE: we compute Eqp0 with SX+CH!!
            ! ZL: Vxc vs. KIH
            if(.not.sig%use_kihdat) then
              fmin(:) = wfnk%elda(sig%diag(jj)-ncore_excl,ispin) + asigt_cor(:) - alda(jj,ispin) - freqs(:)
            else
              fmin(:) = asigt_cor(:) + akih(jj,ispin) - freqs(:)
            endif
          endif
          fmin_nosr(:) = fmin(:) - achcor(jj,ispin)

          call get_eqp1(sig%nfreqeval, freqs, fmin, efsto(jj,ispin)+achcor(jj,ispin), &
            enew(jj,ispin), neqp1(jj,ispin))
          call get_eqp1(sig%nfreqeval, freqs, fmin_nosr, efsto(jj,ispin), &
            enew_nosr(jj,ispin), neqp1_nosr(jj,ispin))
        endif

      enddo ! jj
    enddo ! ii
  enddo ! ispin

!      if (peinf%inode .eq. 0) then
!        do iw = 1, sig%nfreqeval
!          eval = freq0 + (iw-1)*sig%freqevalstep
!          write(8001,2002) kp%rk(:,ikn), eval, specsum(iw)
!        enddo
!          write(8001,*) ''
!      endif

! CHP: not to printout the spectral function

! 2000 format(3F12.5,2x,3i4,2x,7F12.5)
! 2001 format("#",6x,"kx",10x,"ky",10x,"kz",7x,"spn",1x,"bnd",2x,"iw",8x,"Ew",10x,"EQP", &
!      5x,"RE(SIGMA)",3x,"IM(SIGMA)",6x,"SPEC")
! 2002 format(3F12.5,2x,2F12.5)
 2000 format(3F12.5,2x,3i4,2x,5F12.5) !RA
 3000 format(3F12.5,2x,3i4,2x,9F12.5) !CD
 2001 format("#",6x,"kx",10x,"ky",10x,"kz",7x,"spn",1x,"bnd",2x,"iw",8x,"Ew",3x, &
      "   Re(SX+CH)", "   Im(SX+CH)", "   Re(X+Cor)", "   Im(X+Cor)")
 2002 format("#",6x,"kx",10x,"ky",10x,"kz",7x,"spn",1x,"bnd",2x,"iw",8x,"Ew",3x, &
      "   Re(X+Cor)", "   Im(X+Cor)", "     Re(Res)", "     Im(Res)", &
      "     Re(Cor)", "     Im(Cor)", "     Re(Int)", "     Im(Int)")

  POP_SUB(shiftenergy_dyn)

  return

contains

  !> FHJ: Find the roots of fmin(w) to get Eqp1. The strategy is to find all
  !! segments where fmin goes from a positive value to a negative value. We
  !! don`t consider solutions where dSigma/dE>0 because they typically have
  !! a large imag. part associated. If we find multiple solutions, we keep
  !! the one closest to Eqp0. If we don`t find any solution, we try to extrap.
  !! fmin. If there is still no solution, we use Eqp0.
  subroutine get_eqp1(nf, freqs, fmin, eqp0, eqp1, neqp1)
    integer, intent(in) :: nf
    real(DP), intent(in) :: freqs(nf)
    complex(DPC), intent(in) :: fmin(nf)
    complex(DPC), intent(in) :: eqp0
    complex(DPC), intent(out) :: eqp1
    integer, intent(out) :: neqp1

    integer :: iw, nsols, isol
    complex(DPC) :: solns(nf)
    real(DP) :: rfmin(nf), ifmin(nf), dw, rsoln

    PUSH_SUB(shiftenergy_dyn.get_eqp1)

    if (nf<2) then
      eqp1 = eqp0
      neqp1 = 0
      POP_SUB(shiftenergy_dyn.get_eqp1)
      return
    endif

    nsols = 0
    isol = 1 ! Default is to use first solution
    solns(:) = (0d0, 0d0)
    rfmin(:) = dble(fmin(:))
    ifmin(:) = AIMAG(fmin(:))
    dw = freqs(2) - freqs(1)
    do iw=2, nf
      if (rfmin(iw-1)>0 .and. rfmin(iw)<0) then
        nsols = nsols + 1
        ! FHJ: This works for the real and imag. parts
    !    solns(nsols) = freqs(iw-1) + (fmin(iw-1)*dw)/(rfmin(iw-1)-rfmin(iw))
        ! JWJ: Linear interpolation of real part and imaginary part
        rsoln = freqs(iw-1) + (rfmin(iw-1)*dw)/(rfmin(iw-1)-rfmin(iw))
        solns(nsols) = CMPLX(rsoln,((freqs(iw)-rsoln)*ifmin(iw-1)+(rsoln-freqs(iw-1))*ifmin(iw))/dw)
        ! FHJ: Is this solution better than the other ones?
        if (abs(dble(solns(nsols)-eqp0)) < abs(dble(solns(isol)-eqp0)) .or. nsols==1) then
          isol = nsols
        endif
      endif
    enddo
    neqp1 = nsols
    if (nsols==0) then
      if (rfmin(1)<0 .and. rfmin(1)>rfmin(2)) then
        ! FHJ: Extrap. left
        !  solns(1) = (fmin(2)*freqs(1) - fmin(1)*freqs(2))/(rfmin(2) - rfmin(1))
        ! JWJ: Fixed imag part Extrap. left
          rsoln = (rfmin(2)*freqs(1) - rfmin(1)*freqs(2))/(rfmin(2) - rfmin(1))
          solns(1) = CMPLX(rsoln,((freqs(2)-rsoln)*ifmin(1)+(rsoln-freqs(1))*ifmin(2))/dw)
          neqp1 = -1
      elseif (rfmin(nf)>0 .and. rfmin(nf-1)>rfmin(nf)) then
        ! FHJ: Extrap. right
        !  solns(1) = (fmin(nf-1)*freqs(nf) - fmin(nf)*freqs(nf-1))/(rfmin(nf-1) - rfmin(nf))
        ! JWJ: Fixed imag part Extrap. right
          rsoln = (rfmin(nf-1)*freqs(nf) - rfmin(nf)*freqs(nf-1))/(rfmin(nf-1) - rfmin(nf))
          solns(1) = CMPLX(rsoln,((freqs(nf)-rsoln)*ifmin(nf-1)+(rsoln-freqs(nf-1))*ifmin(nf))/dw)
          neqp1 = -2
      else
        ! FHJ: I give up, just use Eqp0
        solns(nsols) = eqp0
        neqp1 = 0
      endif
    endif
    eqp1 = solns(isol)

    POP_SUB(shiftenergy_dyn.get_eqp1)

  end subroutine get_eqp1

end subroutine shiftenergy_dyn

end module shiftenergy_dyn_m
