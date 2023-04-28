!===========================================================================
!
! Routines:
!
! (1) write_result_dyn_hp()        Originally By ?         Last Modified 7/3/2008 (JRD)
!
!     Writes the quasiparticle spectrum to the output
!
!===========================================================================

#include "f_defs.h"

module write_result_dyn_hp_m

  use global_m
  implicit none

  private
  public :: write_result_dyn_hp

contains

subroutine write_result_dyn_hp(kp, wfnk, sig, ax, asx, ach, ach_cor, achcor, &
  asig, alda, efsto, enew, enew_nosr, neqp1, neqp1_nosr, ikn, ncore_excl, akih)
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
  character(len=12) :: soln_str, soln_str_nosr
  
  PUSH_SUB(write_result_dyn_hp)

  ! ZL: check optional akih
  if(sig%use_kihdat) then
    if (.not.present(akih)) call die('In write_result_dyn_hp, to use KIH, pass akih into this subroutine.')
  endif

  do ispin=1,sig%nspin

! Sigma Diagonal

    SAFE_ALLOCATE(iwlda, (sig%ndiag))

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

    write(8,979) (kp%rk(j,ikn),j=1,3),ikn,sig%spin_index(ispin)
979 format(7x,"k =",3f10.6,1x,"ik =",i4,1x,"spin =",i2)
    
    if(sig%exact_ch==0) then
      if (sig%freq_dep_method==2) then
        write(8,*)
        if(.not.sig%use_kihdat) then
          write(8,900) &
            "n", "Emf", "Eo", "X", "Re Res", "Re Int", "Re Sig", "Re Cor", "Vxc", "Re Eqp0", "Re Eqp1", "Soln"
          write(8,900) &
            '', '', '', '', "Im Res", "Im Int", "Im Sig", "Im Cor", '', "Im Eqp0", "Im Eqp1"
        else
          write(8,900) &
            "n", "Emf", "Eo", "X", "Re Res", "Re Int", "Re Sig", "Re Cor", "KIH", "Re Eqp0", "Re Eqp1", "Soln"
          write(8,900) &
            '', '', '', '', "Im Res", "Im Int", "Im Sig", "Im Cor", '', "Im Eqp0", "Im Eqp1"
        endif
      else
        write(8,*)
        if(.not.sig%use_kihdat) then
          write(8,900) &
            "n", "Emf", "Eo", "X", "Re SX-X", "Re ch", "Re Sig", "Re Cor", "Vxc", "Re Eqp0", "Re Eqp1", "Soln"
          write(8,900) &
            '', '', '', '', "Im SX-X", "Im ch", "Im Sig", "Im Cor", '', "Im Eqp0", "Im Eqp1"
        else
          write(8,900) &
            "n", "Emf", "Eo", "X", "Re SX-X", "Re ch", "Re Sig", "Re Cor", "KIH", "Re Eqp0", "Re Eqp1", "Soln"
          write(8,900) &
            '', '', '', '', "Im SX-X", "Im ch", "Im Sig", "Im Cor", '', "Im Eqp0", "Im Eqp1"
        endif
      endif
      do i = 1, sig%ndiag
        call get_soln_str(neqp1(i,ispin), soln_str)
        if(.not.sig%use_kihdat) then
          write(8,901) sig%diag(i), wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
            wfnk%ek(sig%diag(i)-ncore_excl,ispin), dble(ax(i,ispin)), dble(asx(iwlda(i),i,ispin)), &
            dble(ach(iwlda(i),i,ispin)+achcor(i,ispin)), dble(asig(i,ispin)+achcor(i,ispin)), &
            dble(ach_cor(iwlda(i),i,ispin)+achcor(i,ispin)), dble(alda(i,ispin)), &
            dble(efsto(i,ispin)+achcor(i,ispin)), dble(enew(i,ispin)), trim(soln_str)
          write(8,902) IMAG(asx(iwlda(i),i,ispin)), &
            IMAG(ach(iwlda(i),i,ispin)), IMAG(asig(i,ispin)), &
            IMAG(ach_cor(iwlda(i),i,ispin)), IMAG(efsto(i,ispin)), IMAG(enew(i,ispin))
        else
          write(8,901) sig%diag(i), wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
            wfnk%ek(sig%diag(i)-ncore_excl,ispin), dble(ax(i,ispin)), dble(asx(iwlda(i),i,ispin)), &
            dble(ach(iwlda(i),i,ispin)+achcor(i,ispin)), dble(asig(i,ispin)+achcor(i,ispin)), &
            dble(ach_cor(iwlda(i),i,ispin)+achcor(i,ispin)), dble(akih(i,ispin)), &
            dble(efsto(i,ispin)+achcor(i,ispin)), dble(enew(i,ispin)), trim(soln_str)
          write(8,902) IMAG(asx(iwlda(i),i,ispin)), &
            IMAG(ach(iwlda(i),i,ispin)), IMAG(asig(i,ispin)), &
            IMAG(ach_cor(iwlda(i),i,ispin)), IMAG(efsto(i,ispin)), IMAG(enew(i,ispin))
        endif
      enddo
    else         
      if (sig%freq_dep_method==2) then
        write(8,*)
        if(.not.sig%use_kihdat) then
          write(8,900) &
            "n", "Emf", "Eo", "X", "Re Res", "Re Int", "Re Sig", "Re Cor", "Vxc", &
            "Re Eqp0", "Re Eqp1", "Soln", "Re Int`", "Re Sig`", "Re Cor`", "Re Eqp0`", "Re Eqp1`", "Soln`"
          write(8,900) &
            '', '', '', '', "Im Res", "Im Int", "Im Sig", "Im Cor", '', "Im Eqp0", "Im Eqp1", &
            '', '', '', '', '', "Im Eqp1`"
        else
          write(8,900) &
            "n", "Emf", "Eo", "X", "Re Res", "Re Int", "Re Sig", "Re Cor", "KIH", &
            "Re Eqp0", "Re Eqp1", "Soln", "Re Int`", "Re Sig`", "Re Cor`", "Re Eqp0`", "Re Eqp1`", "Soln`"
          write(8,900) &
            '', '', '', '', "Im Res", "Im Int", "Im Sig", "Im Cor", '', "Im Eqp0", "Im Eqp1", &
            '', '', '', '', '', "Im Eqp1`"
        endif ! not kih
      else
        write(8,*)
        if(.not.sig%use_kihdat) then
          write(8,900) &
            "n", "Emf", "Eo", "X", "Re SX-X", "Re CH", "Re Sig", "Re Cor", "Vxc", &
            "Re Eqp0", "Re Eqp1", "Soln", "Re CH`", "Re Sig`", "Re Cor`", "Re Eqp0`", "Re Eqp1`", "Soln`"
          write(8,900) &
            '', '', '', '', "Im SX-X", "Im CH", "Im Sig", "Im Cor", '', "Im Eqp0", "Im Eqp1", &
            '', '', '', '', '', "Im Eqp1`"
        else
          write(8,900) &
            "n", "Emf", "Eo", "X", "Re SX-X", "Re CH", "Re Sig", "Re Cor", "KIH", &
            "Re Eqp0", "Re Eqp1", "Soln", "Re CH`", "Re Sig`", "Re Cor`", "Re Eqp0`", "Re Eqp1`", "Soln`"
          write(8,900) &
            '', '', '', '', "Im SX-X", "Im CH", "Im Sig", "Im Cor", '', "Im Eqp0", "Im Eqp1", &
            '', '', '', '', '', "Im Eqp1`"
        endif ! not kih
      endif
      do i = 1, sig%ndiag
        call get_soln_str(neqp1(i,ispin), soln_str)
        call get_soln_str(neqp1_nosr(i,ispin), soln_str_nosr)
        if(.not.sig%use_kihdat) then
          write(8,901) sig%diag(i), wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
            wfnk%ek(sig%diag(i)-ncore_excl,ispin), dble(ax(i,ispin)), dble(asx(iwlda(i),i,ispin)), &
            dble(ach(iwlda(i),i,ispin)+achcor(i,ispin)), dble(asig(i,ispin)+achcor(i,ispin)), &
            dble(ach_cor(iwlda(i),i,ispin)+achcor(i,ispin)), dble(alda(i,ispin)), &
            dble(efsto(i,ispin)+achcor(i,ispin)), dble(enew(i,ispin)), trim(soln_str), &
            dble(ach(iwlda(i),i,ispin)), dble(asig(i,ispin)), dble(ach_cor(iwlda(i),i,ispin)), &
            dble(efsto(i,ispin)), dble(enew_nosr(i,ispin)), trim(soln_str_nosr)
          write(8,902) IMAG(asx(iwlda(i),i,ispin)), &
            IMAG(ach(iwlda(i),i,ispin)), IMAG(asig(i,ispin)), &
            IMAG(ach_cor(iwlda(i),i,ispin)), IMAG(efsto(i,ispin)), IMAG(enew(i,ispin)), &
            IMAG(enew_nosr(i,ispin))
        else
          write(8,901) sig%diag(i), wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
            wfnk%ek(sig%diag(i)-ncore_excl,ispin), dble(ax(i,ispin)), dble(asx(iwlda(i),i,ispin)), &
            dble(ach(iwlda(i),i,ispin)+achcor(i,ispin)), dble(asig(i,ispin)+achcor(i,ispin)), &
            dble(ach_cor(iwlda(i),i,ispin)+achcor(i,ispin)), dble(akih(i,ispin)), &
            dble(efsto(i,ispin)+achcor(i,ispin)), dble(enew(i,ispin)), trim(soln_str), &
            dble(ach(iwlda(i),i,ispin)), dble(asig(i,ispin)), dble(ach_cor(iwlda(i),i,ispin)), &
            dble(efsto(i,ispin)), dble(enew_nosr(i,ispin)), trim(soln_str_nosr)
          write(8,902) IMAG(asx(iwlda(i),i,ispin)), &
            IMAG(ach(iwlda(i),i,ispin)), IMAG(asig(i,ispin)), &
            IMAG(ach_cor(iwlda(i),i,ispin)), IMAG(efsto(i,ispin)), IMAG(enew(i,ispin)), &
            IMAG(enew_nosr(i,ispin))
        endif ! not kih
      enddo
    endif
900 format(a6,17(a12))
901 format(i6,10(f12.6),a12,5(f12.5),a12)
902 format(6x,3(12x),4(f12.6),12x,2(f12.6),5(12x),f12.6)

! Sigma Off-Diagonal

    if (sig%noffdiag.gt.0) then
      if (sig%freq_dep_method .eq. 2) then
        if(sig%exact_ch == 0) then
          if(.not.sig%use_kihdat) then
            write(8,959)
          else
            write(8,9591)
          endif ! not kih
        else
          if(.not.sig%use_kihdat) then
            write(8,960)
          else
            write(8,9601)
          endif ! not kih
        endif
      else
        if(sig%exact_ch == 0) then
          if(.not.sig%use_kihdat) then
            write(8,969)
          else
            write(8,9691)
          endif ! not kih
        else
          if(.not.sig%use_kihdat) then
            write(8,970)
          else
            write(8,9701)
          endif ! not kih
        endif
      endif
      do i=sig%ndiag+1,sig%ndiag+sig%noffdiag
        iw=iwlda(sig%offmap(i-sig%ndiag, 3))
        if (sig%exact_ch .eq. 0) then
          if(.not.sig%use_kihdat) then
            write(8,968) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
              sig%off3(i-sig%ndiag),dble(ax(i,ispin)),dble(asx(iw,i,ispin)), &
              dble(ach(iw,i,ispin)+achcor(i,ispin)),dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
              dble(alda(i,ispin))
          else
            write(8,968) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
              sig%off3(i-sig%ndiag),dble(ax(i,ispin)),dble(asx(iw,i,ispin)), &
              dble(ach(iw,i,ispin)+achcor(i,ispin)),dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
              dble(akih(i,ispin))
          endif ! not kih
        else
          if(.not.sig%use_kihdat) then
            write(8,971) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
              sig%off3(i-sig%ndiag),dble(ax(i,ispin)),dble(asx(iw,i,ispin)), &
              dble(ach(iw,i,ispin)+achcor(i,ispin)),dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
              dble(alda(i,ispin)), dble(ach(iw,i,ispin)), dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin))
          else
            write(8,971) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
              sig%off3(i-sig%ndiag),dble(ax(i,ispin)),dble(asx(iw,i,ispin)), &
              dble(ach(iw,i,ispin)+achcor(i,ispin)),dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
              dble(akih(i,ispin)), dble(ach(iw,i,ispin)), dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin))
          endif ! not kih
        endif
#ifdef CPLX
        if(.not.sig%use_kihdat) then
          write(8,967) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
            sig%off3(i-sig%ndiag),IMAG(ax(i,ispin)),IMAG(asx(iw,i,ispin)), &
            IMAG(ach(iw,i,ispin)+achcor(i,ispin)),IMAG(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
            IMAG(alda(i,ispin))
        else
          write(8,967) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
            sig%off3(i-sig%ndiag),IMAG(ax(i,ispin)),IMAG(asx(iw,i,ispin)), &
            IMAG(ach(iw,i,ispin)+achcor(i,ispin)),IMAG(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
            IMAG(akih(i,ispin))
        endif ! not kih
#else
972 format(3i4,3x,"imag",3x,7f12.6)
        write(8,972) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
          sig%off3(i-sig%ndiag),0.0d0,IMAG(asx(iw,i,ispin)), &
          IMAG(ach(iw,i,ispin)+achcor(i,ispin)),IMAG(asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
          0.0d0
#endif
      enddo
    endif
969 format(/,3x,"n",3x,"m",3x,"l",21x,"X",8x,"SX-X",10x,"CH",9x, &
      "Sig",9x,"Vxc")
9691 format(/,3x,"n",3x,"m",3x,"l",21x,"X",8x,"SX-X",10x,"CH",9x, &
      "Sig",9x,"KIH")
959 format(/,3x,"n",3x,"m",3x,"l",21x,"X",9x,"Res",9x,"Int",9x, &
      "Sig",9x,"Vxc")
9591 format(/,3x,"n",3x,"m",3x,"l",21x,"X",9x,"Res",9x,"Int",9x, &
      "Sig",9x,"KIH")
968 format(3i4,3x,"real",3x,5f12.6)
967 format(3i4,3x,"imag",3x,5f12.6)
971 format(3i4,3x,"real",3x,7f12.6)
970 format(/,3x,"n",3x,"m",3x,"l",21x,"X",8x,"SX-X",10x,"CH",9x, &
      "Sig",9x,"Vxc",9x,"CH`",8x,"Sig`")
9701 format(/,3x,"n",3x,"m",3x,"l",21x,"X",8x,"SX-X",10x,"CH",9x, &
      "Sig",9x,"KIH",9x,"CH`",8x,"Sig`")
960 format(/,3x,"n",3x,"m",3x,"l",21x,"X",9x,"Res",9x,"Int",9x, &
      "Sig",9x,"Vxc",8x,"Int`",8x,"Sig`")
9601 format(/,3x,"n",3x,"m",3x,"l",21x,"X",9x,"Res",9x,"Int",9x, &
      "Sig",9x,"KIH",8x,"Int`",8x,"Sig`")
   
    SAFE_DEALLOCATE(iwlda)
    
    write(8,*)
    
  enddo ! ispin
  
  POP_SUB(write_result_dyn_hp)  
  return

contains

  subroutine get_soln_str(neqp, str)
    integer, intent(in) :: neqp
    character(len=*), intent(out) :: str

    PUSH_SUB(write_result_dyn_hp.get_soln_str)

    select case (neqp)
      case (-2)
        str = 'extrap+'
      case (-1)
        str = 'extrap-'
      case (0)
        str = 'NO_SOLN!'
      case (1)
        str = 'unique'
      case default
        write(soln_str,'("MULT:",i0)') neqp
    endselect

    POP_SUB(write_result_dyn_hp.get_soln_str)

  end subroutine get_soln_str

end subroutine write_result_dyn_hp

end module write_result_dyn_hp_m
