#include "f_defs.h"

!*************************************************************************
Module TagHandlerMod
!*************************************************************************
! The module contains the public functions:
!    TagHandlerInit(tagHandler, fileno, filename)
!                              Opens the file filename with unit number
!                              fileno, then sets up the case insensitivity.
!                              Copies the filename.
!    TagHandlerDestroy         Closes the unit fileno, deallocates
!                              the copy of filename.
!    FindTag(tagname, fterr)   Searches the file filename from the 
!                              begining to find the tag tagname.  It
!                              returns fterr = 0 if the tag was found.
!                              The unit fileno will be set at
!                              immediately after tagname.  Search is
!                              case insensitive.
!   StringCompare(str1, str2)  Logical function that returns true if
!                              str1 has a subset equal to str2.  The
!                              comparison is case insensitive.
!*************************************************************************
! The module contains the public variables:
!    fileno                    The unit number of the file that taghandler
!                              will search.
!*************************************************************************
! The module contains the private functions:
!    CopyString(str, pstr)     Allocates space in pstr equal to the size of
!                              str then copies the contents of str into pstr.
!    CUpperCase(c)             Converts the character c to uppercase.
!    CopyToUppercase(str,pstr) Like CopyString except that it converts the
!                              string to uppercase letters.
!*************************************************************************
! The module contains the private variables:
!    ucase                     The integer index of A in the ascii sequence.
!    lcase                     The integer index of a in the ascii sequence.
!    filename                  A copy of the parameter to TagHandlerInit.
!*************************************************************************
! Rules:                       Tags must be contain only alphabetic characters.
!                              Tags are case insensitive.
!*************************************************************************

  use push_pop_m
  implicit none

  integer              :: ucase, lcase
 
  type tagHandlerT
    integer            :: fileno
    character, pointer :: filename(:)
  end type tagHandlerT

  private ucase
  private lcase
  private CopyString
  private CUpperCase
  private CopyToUpperCase

  contains

!*************************************************************************
  Subroutine CopyString(str, pstr)
!*************************************************************************

    implicit none

    character*(*), intent(in) :: str
    character, pointer        :: pstr(:)
    integer                   :: error, i, clen

    PUSH_SUB(CopyString)

    clen = len(str)
    allocate(pstr(clen), stat = error)
    if(error /= 0) then
      stop 'Error allocating space in CopyString.'
    end if

    do i = 1, clen
      pstr(i) = str(i:i)
    end do

    POP_SUB(CopyString)
  end Subroutine CopyString

!*************************************************************************
  Function CUpperCase(c)
!*************************************************************************

    implicit none

    character, intent(in) :: c
    character             :: CUpperCase
    integer               :: i

    ! no push_sub, called too frequently

    i = ichar(c)
    CUpperCase = char(mod(i,lcase) + int(i/lcase) * ucase)

  end function CUpperCase

!*************************************************************************
  Subroutine CopyToUpperCase(str, pstr)
!*************************************************************************

    implicit none

    character*(*), intent(in) :: str
    character, pointer        :: pstr(:)
    integer                   :: clen, i

    PUSH_SUB(CopyToUpperCase)

    call CopyString(str, pstr)
! Changed below since we know from routine above that lbound = 1.
! And if it weren`t, the loop below would fail despite this 'precaution.' --DAS
!    clen = ubound(pstr,1) - lbound(pstr,1) + 1
    clen = ubound(pstr,1)

    do i = 1, clen
      pstr(i) = CUpperCase(pstr(i))
    end do

    POP_SUB(CopyToUpperCase)
  end Subroutine CopyToUpperCase

!*************************************************************************
  Function StringCompare(str1, str2)
!*************************************************************************

    implicit none

    character*(*), intent(in) :: str1, str2
    character, pointer        :: ustr1(:), ustr2(:)
    integer                   :: clen1, clen2, i, match, error
    logical                   :: StringCompare

    PUSH_SUB(StringCompare)

    StringCompare = .False.
    call CopyToUpperCase(str1, ustr1)
    call CopyToUpperCase(str2, ustr2)
    clen1 = len(str1)
    clen2 = len(str2)
    if(clen1 .lt. 1 .or. clen2 .lt. 1) then
      stop 'Error in StringCompare, string has zero length.'
    end if
    match = 1
    i = 1

    do while(match .le. clen2 .and. i .le. clen1)
      if(ustr1(i) .eq. ustr2(match)) then
        match = match + 1
      else
        match = 1
      end if
      i = i + 1
    end do

    if(match .gt. clen2) StringCompare = .True.

    deallocate(ustr1, ustr2, stat = error)
    if(error /= 0) stop 'Error deallocating ustr1 in StingCompare.'

    POP_SUB(StringCompare)
  end Function StringCompare

!*************************************************************************
  Subroutine TagHandlerInit(tagHandler, no, fname)
!*************************************************************************

    use message_m
    implicit none

    type( tagHandlerT ), pointer :: tagHandler
    integer, intent(in)          :: no
    character*(*), intent(in)    :: fname
    integer                      :: error

    PUSH_SUB(TagHandlerInit)

    if(associated( tagHandler )) stop 'Error, tagHandler is already allocated.'

    allocate( tagHandler , stat = error )
    if( error /= 0 ) stop 'Error, tagHandler allocation failed.'

    call open_file(unit = no, file = fname, status = 'OLD', FORM = 'formatted')

    taghandler%fileno = no
    call CopyString(fname, taghandler%filename)
    ucase = ichar('A')
    lcase = ichar('a')

    if(lcase .lt. ucase) then

      error = lcase
      lcase = ucase
      ucase = error

    end if

    POP_SUB(TagHandlerInit)
  end Subroutine TagHandlerInit

!*************************************************************************
  Subroutine TagHandlerDestroy( tagHandler )
!*************************************************************************

    use message_m
    implicit none

    type( tagHandlerT ), pointer :: tagHandler
    integer :: error

    PUSH_SUB(TagHandlerDestroy)

    if( .not. associated( tagHandler )) stop 'Error; tagHandler not allocated.'

    call close_file(unit = taghandler%fileno)

!   close_file handles this issue now. --DAS
!    if(error /= 0) then
!      stop 'Error closing file in TagHandler.'
!    end if

!BUS ERROR begin
    deallocate(tagHandler%filename, stat = error)

    if(error /= 0) then
      stop 'Error deallocating filename in TagHandler.'
    end if

    deallocate( tagHandler, stat = error )
    
    if(error /= 0) then
      stop 'Error deallocating TagHandler.'
    end if

    nullify( tagHandler )
!BUS ERROR end

    POP_SUB(TagHandlerDestroy)
  end Subroutine TagHandlerDestroy

!*************************************************************************
  Subroutine FindTag(tagHandler, tag, fterr)
!*************************************************************************

    implicit none

    type( tagHandlerT ), pointer :: tagHandler
    character*(*), intent(in)    :: tag
    integer, intent(out)         :: fterr
    character, pointer           :: uppercasetag(:)
    character                    :: c
    integer                      :: error, clen, match

    PUSH_SUB(FindTag)

    if( .not. associated( tagHandler )) stop 'Error, tagHandler not allocated in FindTag.'
    fterr = 0
    clen = len(tag)
    rewind(unit = tagHandler%fileno, iostat = error)

    if(error /= 0) then
      stop 'Error rewinding file in FindTag'
      fterr = 2
    end if

    match = 1
    
    call CopyToUpperCase(tag, uppercasetag)

    do while(match .le. clen)
      read(unit = tagHandler%fileno, FMT = '(A1)', advance = 'NO', &
          & iostat = error, end = 100) c
      if(CUpperCase(c) .eq. uppercasetag(match)) then
        match = match + 1
      else
        match = 1
      end if
    end do

100 if(match .le. clen) then
      print *, 'Warning: tag ', tag, ' could not be found in the input file.'
      fterr = 1
    end if

    deallocate(uppercasetag, stat = error)
    if(error /= 0) then
      stop 'Error, in TagHandler, deallocating uppercasetag.'
    end if

    POP_SUB(FindTag)
  end Subroutine FindTag

end Module TagHandlerMod

