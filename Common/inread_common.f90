!===============================================================================
!
! Modules:
!
! inread_common_m   Originally By FHJ
!
!   A first attempt to unify the inread routines. Right now, this module
!   only implements consistency checks and warning/error messages.
!
!===============================================================================

#include "f_defs.h"

module inread_common_m

  use global_m
  implicit none

  private

  public ::                   &
    check_bounds_nkq,         &
    check_bounds_nbands,      &
    check_consistency_nbands, &
    try_inread_truncation,    &
    try_inread_screening,     &
    print_truncation_summary

contains


  !> FHJ: Makes sure the number of {q,k}-points is less than MAX_KPTS.
  !! Call this subroutine just after you read the keyword `number_{k,q}-points`.
  subroutine check_bounds_nkq(nkq, k_or_q, keyword)
    !> Number of {k,q}-points expected to be read (eg, pol%nq)
    integer, intent(in) :: nkq
    character(len=1), intent(in) :: k_or_q !< Either "k" or "q"
    character(len=*), intent(in) :: keyword !< a keyword, such as `number_qpoints`

    PUSH_SUB(check_bounds_nkq)

    if (nkq>MAX_KPTS) then
      write(0,*)
      write(0,'(/,a)') 'ERROR: The number of '//&
        k_or_q//'-points specified in the keyword `'//keyword//'` is '
      write(0,'(a,i0)') ' larger than the maximum: MAX_KPTS=',MAX_KPTS
      write(0,'(a,/)') ' Either use less '//&
        k_or_q//'-points or increase MAX_KPTS in Common/nrtype.f90'
      write(0,*)
      call die('Too many '//k_or_q//'-points. Increase MAX_KPTS in Common/nrtypes.f90.')
    endif

    POP_SUB(check_bounds_nkq)

  end subroutine check_bounds_nkq

  !> FHJ: Makes sure nb<MAX_BANDS
  !! Call this subroutine just after you read the keyword `number_bands`.
  subroutine check_bounds_nbands(nb, keyword)
    integer, intent(in) :: nb
    character(len=*), intent(in) :: keyword !< a keyword, such as `number_qpoints`

    PUSH_SUB(check_bounds_nbands)

    if(nb>MAX_BANDS) then
      write(0,*)
      write(0,'(a)') 'ERROR: The number of bands specified in the keyword `'//keyword//'` is larger '
      write(0,'(a,i0)') ' than the maximum: MAX_BANDS=',MAX_BANDS
      write(0,'(a)') ' Either use less bands or increase MAX_BANDS in Common/nrtype.f90'
      write(0,*)
      call die("Too many bands. Increase MAX_BANDS in Common/nrtypes.f90.")
    endif

    POP_SUB(check_bounds_nbands)

  end subroutine check_bounds_nbands


  !> FHJ: Makes sure nb<MAX_BANDS and that nb!=0 if the keyword is required.
  !! Call this after you parse the whole input file.
  subroutine check_consistency_nbands(nb, is_required)
    integer, intent(in) :: nb !< number of bands
    logical, intent(in) :: is_required !< dies if nb==0

    PUSH_SUB(check_consistency_nbands)
    
    if (peinf%inode>0) then
      POP_SUB(check_consistency_nbands)
      return
    endif

    if(is_required .and. nb<1) then
      call die("The keyword `number_bands` could not be found.", only_root_writes = .true.)
    endif

    call check_bounds_nbands(nb, 'number_bands')

    POP_SUB(check_consistency_nbands)

  end subroutine check_consistency_nbands


  !> Try to read input options related to the truncation of the Coulomb
  !! potential
  logical function try_inread_truncation(keyword, line, icutv, truncval)
    character(len=*), intent(in) :: keyword, line
    integer, intent(inout) :: icutv
    real(DP), intent(inout) :: truncval

    integer :: ierr
    logical :: found

    PUSH_SUB(try_inread_truncation)

    ierr = 0
    found = .true.
    if(trim(keyword)=='spherical_truncation') then
      icutv = TRUNC_SPHERICAL
    elseif(trim(keyword)=='cell_wire_truncation') then
      icutv = TRUNC_WIRE
    elseif(trim(keyword)=='cell_box_truncation') then
      icutv = TRUNC_BOX
    elseif(trim(keyword)=='cell_slab_truncation') then
      icutv = TRUNC_SLAB
!#BEGIN_INTERNAL_ONLY
    elseif(trim(keyword)=='supercell_box_truncation') then
      icutv = TRUNC_SUPERCELL
!#END_INTERNAL_ONLY
    elseif(trim(keyword).eq.'coulomb_truncation_x') then
      read(line,*,iostat=ierr) truncval
    elseif(trim(keyword).eq.'coulomb_truncation_radius') then
      read(line,*,iostat=ierr) truncval
    else
      found = .false.
    endif

    if (found .and. ierr/=0) call die( &
      'Unexpected characters were found while reading the value for the keyword ' &
      // trim(keyword) // '. ', only_root_writes = .true.)

    try_inread_truncation = found

    POP_SUB(try_inread_truncation)

  end function try_inread_truncation

  !> Try to read input options related to the analytical behavior for the q->0
  !! limit of W
  logical function try_inread_screening(keyword, line, iscreen)
    character(len=*), intent(in) :: keyword, line
    integer, intent(inout) :: iscreen

    integer :: ierr
    logical :: found

    PUSH_SUB(try_inread_screening)

    ierr = 0
    found = .true.
    if(trim(keyword)=='screening_semiconductor') then
      iscreen = SCREEN_SEMICOND
    elseif(trim(keyword)=='screening_graphene') then
      iscreen = SCREEN_GRAPHENE
    elseif(trim(keyword)=='screening_metal') then
      iscreen = SCREEN_METAL
    else
      found = .false.
    endif

    if (found .and. ierr/=0) call die( &
      'Unexpected characters were found while reading the value for the keyword ' &
      // trim(keyword) // '. ', only_root_writes = .true.)

    try_inread_screening = found

    POP_SUB(try_inread_screening)

  end function try_inread_screening


  subroutine print_truncation_summary(icutv, truncval, iunit)
    integer, intent(in) :: icutv
    real(DP), intent(in) :: truncval
    integer, intent(in), optional :: iunit

    integer :: iunit_

    PUSH_SUB(print_truncation_summary)

    iunit_ = 6
    if (present(iunit)) iunit_ = iunit
    if (peinf%inode/=0) then
      POP_SUB(print_truncation_summary)
      return
    endif

    select case (icutv)
      case (TRUNC_NONE)
        write(iunit_,'(1x,a)') 'We are using no truncation'
      case (TRUNC_SPHERICAL)
        write(iunit_,'(1x,a)') 'We are using a truncated Coulomb interaction: spherical'
        write(iunit_,'(1x,"r_cut = ",f10.6," Bohr")') truncval
      case (TRUNC_WIRE)
        write(iunit_,'(1x,a)') 'We are using a truncated Coulomb interaction: cell wire'
      case (TRUNC_BOX)
        write(iunit_,'(1x,a)') 'We are using a truncated Coulomb interaction: cell box'
      case (TRUNC_SLAB)
        write(iunit_,'(1x,a)') 'We are using a truncated Coulomb interaction: cell slab'
      case default
        call die('Unknown truncation flag', only_root_writes=.true.)
      write(iunit_,'()')
    endselect

    POP_SUB(print_truncation_summary)

  end subroutine print_truncation_summary


end module inread_common_m
