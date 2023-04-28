
!-------------------------------------------------------------------------------
!
!     The primitive iolib modules
!
!     The following integer markers have been used:
!
!     Int      Object:
!
!       1      gspace (not the parallel one)
!       2      isrt arrays
!       3      wave functions
!       4      kpoint
!       5      level
!       6      bands
!       7      crystal structure
!       8      symmetry info
!
!     If you add a new one, make sure you update iolib_skip()
!     in iolib_m
!
!-------------------------------------------------------------------------------

#include "f_defs.h"

module iolib_m

  use global_m
  implicit none

contains  

  !
  !     =============== generic routine to find certain marker ===========
  !
  subroutine iolib_find_marker(imark, iunit, info)
    !
    !     INPUT:
    !
    integer, intent(in) :: imark, &    ! integer marker
         iunit                         ! the unit to be browsed
    !
    !     OUTPUT:
    !
    integer, intent(out) :: info       ! 1 if marker found, 0 otherwise
    !
    !     ------------- local variables
    !
    integer :: marker, ierr, iostat

    PUSH_SUB(iolib_find_marker)

    rewind(iunit)  
    !      write(0,*) 'searching for marker:', imark
    do while (.true.)  
       read(iunit, *, iostat=iostat) marker  
       if(iostat < 0) exit
       !         write(0,*) 'found marker:', marker
       ! found the right object
       if (marker == imark) then  
          info = 1  
          exit
       else  
          ! found other object. skip
          call iolib_skip(iunit, marker, ierr)  
          if (ierr /= 0) then
            info = 0
            exit
          endif
       end if
    end do

    POP_SUB(iolib_find_marker)

    return  

  end subroutine iolib_find_marker
  !
  !     ============= routine to skip an object. has to be updated =======
  !
  !
  !
  !
  subroutine iolib_skip(iunit, imark, ierr)
    !
    !     INPUT:
    !
    integer, intent(in) :: iunit, &    ! the io unit to be used
         imark                   ! the marker code of the object to be skipped
    !
    !     OUTPUT:
    !
    integer, intent(out) :: ierr  ! error code. if error, then 1. Zero otherwise
    !
    !     ------- local variables
    !
    integer :: i, iostat
    !
    !     take whatever action is necessary to skip the object
    !

    PUSH_SUB(iolib_skip)

    select case(imark)  
    case default  
       write(0, *) 'iolib: illegal marker: ', imark  
       call die('the linked iolib probably not up to date')
    case(1)  
       ! ------ how to skip object 1
       read(iunit, *, iostat=iostat)
       read(iunit, *, iostat=iostat)
       read(iunit, *, iostat=iostat)  
       read(iunit, *, iostat=iostat)  
       read(iunit, *, iostat=iostat)  
       read(iunit, *, iostat=iostat)  
       read(iunit, *, iostat=iostat)  
       read(iunit, *, iostat=iostat)  
    case(2)  
       ! ------ how to skip object 2
       do i = 1, 4  
          read(iunit, *, iostat=iostat)  
       end do
    case(3)  
       ! ------ how to skip object 3
       do i = 1, 7  
          read(iunit, *, iostat=iostat)  
       end do
    case(4)  
       ! ------ how to skip object 4
       do i = 1, 10  
          read(iunit, *, iostat=iostat)  
       end do
    case(5)  
       ! ------ how to skip object 5
       do i = 1, 10  
          read(iunit, *, iostat=iostat)  
       end do
    case(6)  
       ! ------ how to skip object 6
       do i = 1, 6  
          read(iunit, *, iostat=iostat)  
       end do
    case(7)  
       ! ------ how to skip object 7
       do i = 1, 11  
          read(iunit, *, iostat=iostat)  
       end do
    case(8)  
       ! ------ how to skip object 8
       do i = 1, 5  
          read(iunit, *, iostat=iostat)  
       end do
    end select

    ierr = iostat

    POP_SUB(iolib_skip)
    return  

  end subroutine iolib_skip

end module iolib_m
