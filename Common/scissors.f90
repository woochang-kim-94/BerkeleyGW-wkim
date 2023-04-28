!===============================================================================
!
! Modules:
!
! scissors_m   Originally By DAS
!
!   Routines for scissors corrections to mean-field eigenvalues.
!
!===============================================================================

#include "f_defs.h"

module scissors_m

  use global_m
  use splines_m

  implicit none

  private

  public ::               &
    scissors_shift,        &
    scissors_write,       &
    scissors_zero,        &
    scissors_function,    &
    scissors_inread

contains

  !---------------------------------------------------------------------------------------------------
  !> Scissors operator to get quasi-particle spectrum from LDA eigenvalues
  !! evs, ecs, ev0, ec0 are supposed to be in eV
  !! kp%el is taken to be in Ry however
  !! spl_tck (spline coefficients) in Ry
  subroutine scissors_shift(kp, scis, spl_tck)
    type(kpoints), intent(inout) :: kp
    type(scissors_t), intent(in) :: scis
    type(spline_tck), optional, intent(in) :: spl_tck

    integer :: is, irk, ib
    logical :: has_spline, use_spline
    real(DP) :: spline_emin, spline_emax

    PUSH_SUB(scissors_shift)

    has_spline=.false.
    if (present(spl_tck)) then
      if (spl_tck%n>0) then
        has_spline = .true.
        spline_emin = spl_tck%t(1)/ryd
        spline_emax = spl_tck%t(spl_tck%n)/ryd
      endif
    endif

    do is=1,kp%nspin
      do irk=1,kp%nrk
        do ib=1,kp%mnband

          use_spline=.false.
          if (has_spline) then
            if ( (kp%el(ib,irk,is) > spline_emin) .and. &
                 (kp%el(ib,irk,is) < spline_emax) ) then
                 use_spline=.true.
            endif
          endif

          if (use_spline) then
            ! use spline interpolation for scissors ops
            call splev_shift( spl_tck, kp%el(ib,irk,is) )
          else
            ! use regular (linear) scissors ops
            if(ib <= kp%ifmax(irk,is)) then
              kp%el(ib,irk,is) = scissors_function(kp%el(ib, irk, is), scis%val)
            else
              kp%el(ib,irk,is) = scissors_function(kp%el(ib, irk, is), scis%cond)
            endif
          endif

        enddo ! ib (band)
      enddo ! irk (kpoint)
    enddo ! is (spin)
    
    POP_SUB(scissors_shift)
    return
  end subroutine scissors_shift

  !---------------------------------------------------------------------------------------------------
  subroutine scissors_write(iunit, scis, suffix)
    integer, intent(in) :: iunit !< unit to which to write
    type(scissors_t), intent(in) :: scis
    character(len=*), optional, intent(in) :: suffix !< identifier to be printed

    PUSH_SUB(scissors_write)

    write(iunit,'(/1x,a)',advance='no') "Scissors parameters"
    if(present(suffix)) then
      write(iunit,'(3a)') ' (', TRUNC(suffix), '):'
    else
      write(iunit,'(a)') ":"
    endif
    write(iunit,'(1x,3(a,f8.4))') &
      '- Valence:    es = ', scis%val%es,  ' eV, e0 = ', scis%val%e0,  ' eV, edel = ', scis%val%edel
    write(iunit,'(1x,3(a,f8.4))') &
      '- Conduction: es = ', scis%cond%es, ' eV, e0 = ', scis%cond%e0, ' eV, edel = ', scis%cond%edel

    POP_SUB(scissors_write)
    return
  end subroutine scissors_write

  !---------------------------------------------------------------------------------------------------
  subroutine scissors_zero(scis)
    type(scissors_t), intent(out) :: scis

    PUSH_SUB(scissors_zero)

    scis%val%es    = 0d0
    scis%val%e0    = 0d0
    scis%val%edel  = 0d0
    scis%cond%es   = 0d0
    scis%cond%e0   = 0d0
    scis%cond%edel = 0d0

    POP_SUB(scissors_zero)
    return
  end subroutine scissors_zero

  !---------------------------------------------------------------------------------------------------
  real(DP) function scissors_function(energy, scis)
    real(DP), intent(in) :: energy !< should be in Ry
    type(sub_scissors_t), intent(in) :: scis
    
    PUSH_SUB(scissors_function)
    scissors_function = energy + scis%es / ryd + scis%edel * (energy - scis%e0 / ryd)

    POP_SUB(scissors_function)
    return
  end function scissors_function

  !---------------------------------------------------------------------------------------------------
  subroutine scissors_inread(keyword, line, scis, found, suffix)
    character(len=*), intent(in) :: keyword, line
    type(scissors_t), intent(inout) :: scis
    logical, intent(out) :: found
    character(len=*), optional, intent(in) :: suffix

    integer :: iostat, trunc_length
    character*120 :: trunc_keyword

    PUSH_SUB(scissors_inread)

    if(present(suffix)) then
      if(len(suffix) > 0) then
        trunc_length = len(trim(keyword))
        if(keyword(trunc_length-len(suffix)+1:trunc_length) == suffix) then
          trunc_keyword = keyword(1:trunc_length-len(suffix))
        else
          found = .false.
          POP_SUB(scissors_inread)
          return
        endif
      endif
    else
      trunc_keyword = keyword
    endif

    found = .true.
    if(trim(trunc_keyword).eq.'cvfit') then
      read(line,*,iostat = iostat) scis%val%es,scis%val%e0,scis%val%edel, &
        scis%cond%es,scis%cond%e0,scis%cond%edel
    elseif(trim(trunc_keyword).eq.'evs') then
      read(line,*,iostat = iostat) scis%val%es
    elseif(trim(trunc_keyword).eq.'ev0') then
      read(line,*,iostat = iostat) scis%val%e0
    elseif(trim(trunc_keyword).eq.'evdel') then
      read(line,*,iostat = iostat) scis%val%edel
    elseif(trim(trunc_keyword).eq.'ecs') then
      read(line,*,iostat = iostat) scis%cond%es
    elseif(trim(trunc_keyword).eq.'ec0') then
      read(line,*,iostat = iostat) scis%cond%e0
    elseif(trim(trunc_keyword).eq.'ecdel') then
      read(line,*,iostat = iostat) scis%cond%edel
    else
      found = .false.
    endif

    if(found .and. iostat /= 0) call die( &
      'Unexpected characters were found while reading the value for the keyword ' &
      // trim(keyword) // '. ', only_root_writes = .true.)

    POP_SUB(scissors_inread)
    return
  end subroutine scissors_inread

end module scissors_m
