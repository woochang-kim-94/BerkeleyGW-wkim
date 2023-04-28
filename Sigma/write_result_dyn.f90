!===========================================================================
!
! Routines:
!
! (1) write_result_dyn()        Originally By ?         Last Modified 7/3/2008 (JRD)
!
!     Writes the quasiparticle spectrum to the output
!
!===========================================================================

#include "f_defs.h"

module write_result_dyn_m

  use global_m
  implicit none

  private
  public :: write_result_dyn

contains

subroutine write_result_dyn(kp, wfnk, sig, ax, asx, ach, ach_cor, &
  achcor, asig, alda, efsto, enew, enew_nosr, neqp1, neqp1_nosr, ikn,ncore_excl, akih)
  type (kpoints), intent(in) :: kp
  type (wfnkstates), intent(in) :: wfnk
  type (siginfo), intent(in) :: sig
  SCALAR, intent(in) :: ax(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(in) :: asx(:,:,:) !< (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(in) :: ach(:,:,:) !< (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(in) :: ach_cor(:,:,:) !< (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(in) :: achcor(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(in) :: asig(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin)
  SCALAR, intent(in) :: alda(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(in) :: efsto(:,:) !< (sig%ndiag,sig%nspin)
  complex(DPC), intent(in) :: enew(:,:) !< (sig%ndiag,sig%nspin)
  complex(DPC), intent(in) :: enew_nosr(:,:) !< (sig%ndiag,sig%nspin)
  integer, intent(in) :: neqp1(:,:) !< (sig%ndiag,sig%nspin)
  integer, intent(in) :: neqp1_nosr(:,:) !< (sig%ndiag,sig%nspin)
  integer, intent(in) :: ikn
  integer, intent(in) :: ncore_excl !< number of core states excluded  
  SCALAR, intent(in), optional :: akih(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin), ZL

  integer :: iw, i, j, ispin
  integer, allocatable :: iwlda(:)
  real(DP) :: diff, diffmin, e_lk, freq0

  PUSH_SUB(write_result_dyn)

  ! ZL: check optional akih
  if(sig%use_kihdat) then
    if (.not.present(akih)) call die('In write_result_dyn, to use KIH, pass akih into this subroutine.')
  endif

  SAFE_ALLOCATE(iwlda, (sig%ndiag))
  do ispin=1,sig%nspin

! Sigma Diagonal
! JRD: Find iw closest to e_lk
! DVF : ncore_excl has to be substracted here because wfnk%ek is defined in the 
! read_wavefunction subroutine in input.f90 to be referenced to the case with 
! no core states. 
    do i = 1, sig%ndiag
      diffmin = INF
      e_lk = wfnk%ek(sig%diag(i)-ncore_excl,ispin)
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
          iwlda(i)=iw
        endif
      enddo
    enddo
    call write_qp(6)
  enddo ! ispin
  SAFE_DEALLOCATE(iwlda)
  
  POP_SUB(write_result_dyn)
  
  return

contains

  subroutine write_qp(iunit)
    integer, intent(in) :: iunit

    character(len=9) :: soln_str

    write(iunit,979) (kp%rk(j,ikn),j=1,3),ikn,sig%spin_index(ispin)
979 format(/,7x,"k =",3f10.6,1x,"ik =",i4,1x,"spin =",i2)

    if (sig%freq_dep_method==2) then
      write(iunit,*)
      if(.not.sig%use_kihdat) then
        write(iunit,900) &
          "n", "Emf", "Eo", "Vxc", "X", "Re Cor", "Re Eqp0", "Re Eqp1", "Soln"
        write(iunit,900) &
          '', '', '', '', '', "Im Cor", "Im Eqp0", "Im Eqp1"
      else
        write(iunit,900) &
          "n", "Emf", "Eo", "KIH", "X", "Re Cor", "Re Eqp0", "Re Eqp1", "Soln"
        write(iunit,900) &
          '', '', '', '', '', "Im Cor", "Im Eqp0", "Im Eqp1"
      endif
    else
      write(iunit,*)
      if(.not.sig%use_kihdat) then
        write(iunit,900) &
          "n", "Emf", "Eo", "Vxc", "X", "Re Cor", "Re Eqp0", "Re Eqp1", "Soln"
        write(iunit,900) &
          '', '', '', '', '', "Im Cor", "Im Eqp0", "Im Eqp1"
      else
        write(iunit,900) &
          "n", "Emf", "Eo", "KIH", "X", "Re Cor", "Re Eqp0", "Re Eqp1", "Soln"
        write(iunit,900) &
          '', '', '', '', '', "Im Cor", "Im Eqp0", "Im Eqp1"
      endif
    endif
    do i = 1, sig%ndiag
      select case (neqp1(i,ispin))
        case (-2)
          soln_str = 'extrap+'
        case (-1)
          soln_str = 'extrap-'
        case (0)
          soln_str = 'NO_SOLN!'
        case (1)
          soln_str = 'unique'
        case default
          write(soln_str,'("MULT:",i0)') neqp1(i,ispin)
      endselect
      if(.not.sig%use_kihdat) then
        write(iunit,901) sig%diag(i), wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
          wfnk%ek(sig%diag(i)-ncore_excl,ispin), dble(alda(i,ispin)), &
          dble(ax(i,ispin)), dble(ach_cor(iwlda(i),i,ispin)+achcor(i,ispin)), &
          dble(efsto(i,ispin)+achcor(i,ispin)), dble(enew(i,ispin)), trim(soln_str)
        write(iunit,902) IMAG(ach_cor(iwlda(i),i,ispin)), &
          IMAG(efsto(i,ispin)), IMAG(enew(i,ispin))
      else
        write(iunit,901) sig%diag(i), wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
          wfnk%ek(sig%diag(i)-ncore_excl,ispin), dble(akih(i,ispin)), &
          dble(ax(i,ispin)), dble(ach_cor(iwlda(i),i,ispin)+achcor(i,ispin)), &
          dble(efsto(i,ispin)+achcor(i,ispin)), dble(enew(i,ispin)), trim(soln_str)
        write(iunit,902) IMAG(ach_cor(iwlda(i),i,ispin)), &
          IMAG(efsto(i,ispin)), IMAG(enew(i,ispin))
      endif
    enddo
900 format(a6,8(a9))
901 format(i6,7(f9.3),a9)
902 format(6x,4(9x),3(f9.3))

! Sigma Off-Diagonal

    if (sig%noffdiag.gt.0) then
      if (sig%freq_dep_method .eq. 2) then
        if(.not.sig%use_kihdat) then
          write(iunit,959)
        else
          write(iunit,9591)
        endif
      else
        if(.not.sig%use_kihdat) then
          write(iunit,969)
        else
          write(iunit,9691)
        endif
      endif
969 format(/,3x,"n",3x,"m",3x,"l",10x,6x,"Vxc",8x,"X",5x,"SX-X",7x,"CH",6x,"Sig")
9691 format(/,3x,"n",3x,"m",3x,"l",10x,6x,"KIH",8x,"X",5x,"SX-X",7x,"CH",6x,"Sig")
959 format(/,3x,"n",3x,"m",3x,"l",10x,6x,"Vxc",8x,"X",6x,"Res",6x,"Int",6x,"Sig")
9591 format(/,3x,"n",3x,"m",3x,"l",10x,6x,"KIH",8x,"X",6x,"Res",6x,"Int",6x,"Sig")

      do i=sig%ndiag+1,sig%ndiag+sig%noffdiag
        iw=iwlda(sig%offmap(i-sig%ndiag, 3))

        if(.not.sig%use_kihdat) then
          write(iunit,968) sig%off1(i-sig%ndiag), sig%off2(i-sig%ndiag), &
            sig%off3(i-sig%ndiag), dble(alda(i,ispin)), dble(ax(i,ispin)), &
            dble(asx(iw,i,ispin)), dble(ach(iw,i,ispin)+achcor(i,ispin)), &
            dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin))
        else
          write(iunit,968) sig%off1(i-sig%ndiag), sig%off2(i-sig%ndiag), &
            sig%off3(i-sig%ndiag), dble(akih(i,ispin)), dble(ax(i,ispin)), &
            dble(asx(iw,i,ispin)), dble(ach(iw,i,ispin)+achcor(i,ispin)), &
            dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin))
        endif
968 format(3i4,3x,"real",3x,5f9.3)

        ! FHJ: why do we add IMAG(achcor) only for the real case?!!
#ifdef CPLX
        if(.not.sig%use_kihdat) then
          write(iunit,967) sig%off1(i-sig%ndiag), sig%off2(i-sig%ndiag), &
            sig%off3(i-sig%ndiag), IMAG(alda(i,ispin)), IMAG(ax(i,ispin)), &
            IMAG(asx(iw,i,ispin)), IMAG(ach(iw,i,ispin)), &
            IMAG(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin))
        else
          write(iunit,967) sig%off1(i-sig%ndiag), sig%off2(i-sig%ndiag), &
            sig%off3(i-sig%ndiag), IMAG(akih(i,ispin)), IMAG(ax(i,ispin)), &
            IMAG(asx(iw,i,ispin)), IMAG(ach(iw,i,ispin)), &
            IMAG(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin))
        endif
        ! ZL: note that KIH only supports complex flavor now
        ! so the below lines remain unchanged.
        ! In fact, alda (Vxc) has been set to zeros by hand below
#else
        write(iunit,967) sig%off1(i-sig%ndiag), sig%off2(i-sig%ndiag), &
          sig%off3(i-sig%ndiag), 0d0, 0d0, IMAG(asx(iw,i,ispin)), &
          IMAG(ach(iw,i,ispin)+achcor(i,ispin)), &
          IMAG(asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin))
#endif
967 format(3i4,3x,"imag",3x,5f9.3)
      enddo
    endif

  end subroutine write_qp

end subroutine write_result_dyn

end module write_result_dyn_m
