!===========================================================================
!
! Routines:
!
! (1) ch_converge()    Originally By ?         Last Modified May/2016 (FHJ)
!
!     Write the convergence of the Coulomb-hole term as a function of bands.
!
!===========================================================================

#include "f_defs.h"

module ch_converge_m

  use global_m
  implicit none

  private
  public :: ch_converge

contains

! FHJ: Note that the input arrays ach_n1 and achcor_n1 are all complex, even
! though they may be computed as SCALAR. We do this to simply the call to this
! subroutine by FF and non-FF calculations.
subroutine ch_converge(kp, sig, ach_n1, achcor_n1, ikn)
  type (kpoints), intent(in) :: kp
  type (siginfo), intent(in) :: sig
  complex(DPC), intent(in) :: ach_n1(:,:,:) !< (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(in) :: achcor_n1(:,:,:)  !< (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin)
  integer, intent(in) :: ikn

  integer :: i, j, ivbm, icbm, ispin
  !> FHJ: CH without static remainder
  complex(DPC), dimension(sig%ntband) :: chp_v, chp_c, chp_d
  !> FHJ: Static remainder
  complex(DPC), dimension(sig%ntband) :: sr_v, sr_c
  !> FHJ: CH with static remainder
  complex(DPC), dimension(sig%ntband) :: ch_v, ch_c, ch_d
  complex(DPC) :: extrap_v, extrap_c, extrap_d, extrapp_v, extrapp_c, extrapp_d
  character(len=3) :: str
  character(len=9) :: str_v, str_c, strp_v, strp_c
  
  PUSH_SUB(ch_converge)

  if (sig%freq_dep==-1 .or. (sig%freq_dep==0 .and. sig%exact_ch==1)) then
    PUSH_SUB(ch_converge)
    return
  endif

  ! FHJ: figure out which are the VBM/CBM. This is only used for sig%fullConvLog==0.
  ivbm = 0
  icbm = 0
! DVF : sig%ncore_excl has to be substracted here because matrix elements are referenced to
! case with core states, while nvband has been redefined in input.f90 to exclude core states. 
  do i = 1, sig%ndiag
    if (sig%diag(i)-sig%ncore_excl .eq.sig%nvband) ivbm=i
    if (sig%diag(i)-sig%ncore_excl .eq.sig%nvband+1) icbm=i
  enddo

  ! FHJ: What we really compute is the CH term for GPP and RA calculations, and
  ! the integral contribution to COR for CD calculations.
  str = 'CH'
  if (sig%freq_dep==2 .and. sig%freq_dep_method==2) then
    str = 'Int'
  endif

  do ispin=1,sig%nspin
    write(127,'(a,3(1x,f10.6),a,i1)') '# k =', (kp%rk(j,ikn),j=1,3), ' spin = ',sig%spin_index(ispin)
    write(127,'(3a)') '# <ib|Sigma_', trim(str), &
      '[1..Nb]|ib> (eV), partial and symmetrized sum with up to Nb "inner" bands'

    if (sig%fullConvLog==0) then
      write(127,'(a)') "# Reporting only VBM and CBM, use 'full_ch_conv_log 1' keyword for full report."

      ! FHJ: Calculate CH(band), where "chp" is CH w/o SR.
      chp_v(:) = 0d0
      sr_v(:) = 0d0
      ch_v(:) = 0d0
      if (ivbm.ne.0) then
        call cumsum(ach_n1(:,ivbm,ispin), chp_v)
        call cumsum(achcor_n1(:,ivbm,ispin), sr_v)
        ch_v = chp_v + sr_v
      endif
      chp_c(:) = 0d0
      sr_c(:) = 0d0
      ch_c(:) = 0d0
      if (icbm.ne.0) then
        call cumsum(ach_n1(:,icbm,ispin), chp_c)
        call cumsum(achcor_n1(:,icbm,ispin), sr_c)
        ch_c = chp_c + sr_c
      endif
      chp_d = chp_c - chp_v
      ch_d = ch_c - ch_v

      ! FHJ: 1/N extrapolation and error estimate
      extrap_v = extrap(ch_v)
      extrap_c = extrap(ch_c)
      extrap_d = extrap(ch_d)
      extrapp_v = extrap(chp_v)
      extrapp_c = extrap(chp_c)
      extrapp_d = extrap(chp_d)
      write(str_v,'(2a)') trim(str), '(vbm)'
      write(str_c,'(2a)') trim(str), '(cbm)'
      write(strp_v,'(2a)') trim(str), '`(vbm)'
      write(strp_c,'(2a)') trim(str), '`(cbm)'
      write(127,'(a,3(1x,a),3(1x,es16.8))') '# 1/N extrap for', &
        trim(str_v), trim(str_c), 'diff =', &
        dble((/extrap_v, extrap_c, extrap_d/))
      write(127,'(a,3(1x,a),3(1x,es16.8))') '# Error est. for', &
        trim(str_v), trim(str_c), 'diff =', &
        dble((/ch_v(sig%ntband)-extrap_v, ch_c(sig%ntband)-extrap_c, ch_d(sig%ntband)-extrap_d/))

      if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
        ! FHJ: if we computed the SR correction, we also report the "primed" values w/o SR
        write(127,'(a,3(1x,a),3(1x,es16.8))') '# 1/N extrap for', &
          trim(strp_v), trim(strp_c), 'diff =', &
          dble((/extrapp_v, extrapp_c, extrapp_d/))
        write(127,'(a,3(1x,a),3(1x,es16.8))') '# Error est. for', &
          trim(strp_v), trim(strp_c), 'diff =', &
          dble((/chp_v(sig%ntband)-extrapp_v, chp_c(sig%ntband)-extrapp_c, chp_d(sig%ntband)-extrapp_d/))
        write(127,'(a1,a6,6(1x,a16))') '#', 'Nb' , &
          trim(str_v), trim(str_c), 'diff', &
          trim(strp_v), trim(strp_c), 'diff`'
      else
        write(127,'(a1,a6,3(1x,a16))') '#', 'Nb' , &
          trim(str_v), trim(str_c), 'diff'
      endif

      ! FHJ: Write CH(band)
      if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
        do j = 1, sig%ntband
          write(127,'(1x,i6,6(1x,es16.8))') j, dble((/ch_v(j), ch_c(j), ch_d(j), chp_v(j), chp_c(j), chp_d(j)/))
        enddo
      else
        do j = 1, sig%ntband
          write(127,'(1x,i6,3(1x,es16.8))') j, dble((/ch_v(j), ch_c(j), ch_d(j)/))
        enddo
      endif

    else
      ! FHJ: Full report for all bands, Re and Im parts

      if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
        ! FHJ: Please, keep the double slashes here, because Open64 segfaults
        ! if we use the plus sign for string concatenation.
        write(127,'(a1,a5,1x,a6,4(1x,a16))') '#', 'ib', 'Nb' , &
          'Re['//trim(str)//']', 'Im['//trim(str)//']', 'Re['//trim(str)//'`]', 'Im['//trim(str)//'`]'
      else
        ! FHJ: Please, keep the double slashes here, because Open64 segfaults
        ! if we use the plus sign for string concatenation.
        write(127,'(a1,a5,1x,a6,2(1x,a16))') '#', 'ib', 'Nb' , &
          'Re['//trim(str)//']', 'Im['//trim(str)//']'
      endif

      do i = 1, sig%ndiag
        ! FHJ: Calculate CH(band), where "chp" is CH w/o SR.
        call cumsum(ach_n1(:,i,ispin), chp_v)
        call cumsum(achcor_n1(:,i,ispin), sr_v)
        ch_v = chp_v + sr_v

        ! FHJ: 1/N extrapolation and error estimate
        extrap_v = extrap(ch_v)
        extrapp_v = extrap(chp_v)
        if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
          write(127,'(a1,i5,1x,a6,4(1x,es16.8))') '#', i, 'extrap', extrap_v, extrapp_v
          write(127,'(a1,i5,1x,a6,4(1x,es16.8))') '#', i, 'error', &
            ch_v(sig%ntband)-extrap_v, chp_v(sig%ntband)-extrapp_v
        else
          write(127,'(a1,i5,1x,a6,2(1x,es16.8))') '#', i, 'extrap', extrap_v
          write(127,'(a1,i5,1x,a6,2(1x,es16.8))') '#', i, 'error', ch_v(sig%ntband)-extrap_v
        endif

        ! FHJ: Write CH(band)
        if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
          do j = 1, sig%ntband
            write(127,'(i6,1x,i6,4(1x,es16.8))') i, j, ch_v(j), chp_v(j)
          enddo
        else
          do j = 1, sig%ntband
            write(127,'(i6,1x,i6,2(1x,es16.8))') i, j, ch_v(j)
          enddo
        endif
      enddo
    endif ! sig%fullConvLog

  enddo ! ispin
  FLUSH(127)

  POP_SUB(ch_converge)

  return

contains

  subroutine cumsum(x, y)
    complex(DPC), intent(in) :: x(:)
    complex(DPC), intent(out) :: y(:)

    PUSH_SUB(ch_converge.cumsum)

    y(1) = x(1)
    do j = 2, sig%ntband
      y(j) = x(j) + y(j-1)
    enddo

    POP_SUB(ch_converge.cumsum)

  end subroutine cumsum

  complex(DPC) function extrap(ch)
    complex(DPC), intent(in) :: ch(:)
    integer :: idis

    PUSH_SUB(ch_converge.extrap)

    extrap = 0d0
    idis = sig%ntband / 10
    if (idis>0) then
      extrap = (dble(sig%ntband) * ch(sig%ntband) - &
        dble(sig%ntband - idis) * ch(sig%ntband - idis)) / dble(idis)
    endif

    POP_SUB(ch_converge.extrap)

  end function extrap

end subroutine ch_converge

end module ch_converge_m
