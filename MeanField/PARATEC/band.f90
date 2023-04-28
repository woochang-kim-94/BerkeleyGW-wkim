
#include "f_defs.h"

module band_m

  use global_m
  use iolib_m
  implicit none

  type band 
    integer :: nrk, min(2), max, nspin
    integer, pointer :: nband(:,:), ifmax(:,:)
    real(DP), pointer :: energy(:,:,:), ekn(:,:,:), occup(:,:,:)
  end type band

contains

  subroutine read_band(fname,b)
    character(len=*), intent(in) :: fname      ! the filename
    type(band), intent(inout) :: b

    integer :: j, k, l, info
    integer, parameter :: mymarker = 6 
    integer :: nrk, max, nspin

    PUSH_SUB(read_band)

    call open_file(unit = 21, file = fname, status = 'unknown', &
      form = 'formatted')
    call iolib_find_marker(mymarker, 21, info)
    if (info == 0) then
      write(0,*) 'There is a problem'
      POP_SUB(read_band)
      return  
    end if
    read(21, *, end = 102) b%nrk, b%min, b%max, b%nspin
    nrk = b%nrk
    nspin = b%nspin
    max = b%max
    SAFE_ALLOCATE(b%nband, (nrk, nspin)); b%nband = 0
    SAFE_ALLOCATE(b%occup, (max, nrk, nspin)); b%occup = 0d0
    SAFE_ALLOCATE(b%energy, (max, nrk, nspin)); b%energy = 0d0
    SAFE_ALLOCATE(b%ekn, (max, nrk, nspin)); b%ekn = 0d0
    SAFE_ALLOCATE(b%ifmax, (nrk, nspin)); b%ifmax = 0


    read(21, *, end = 102) ((b%nband(k, j), k = 1, nrk), j = 1, nspin)
    read(21, *, end = 102) ((b%ifmax(k, j), k = 1, nrk), j = 1, nspin)
    read(21, *, end = 102) (((b%energy(l, k, j), l = 1, max), &
      k = 1, nrk), j = 1, nspin)
    read(21, *, end = 102) (((b%ekn(l, k, j), l = 1, max), &
      k = 1, nrk), j = 1, nspin)
    read(21, *, end = 102) (((b%occup(l, k, j), l = 1, max), &
      k = 1, nrk), j = 1, nspin)
    call close_file(21)

    POP_SUB(read_band)
    return

102 continue
    write(0, *) 'iolib: reached eof when reading bands'

    POP_SUB(read_band)

    return
  end subroutine read_band

  subroutine write_band(fname, struc)

    implicit none      ! implicit? no!

    character(len=*), intent(in) :: fname      ! the filename
    type(band), intent(in) :: struc

    !     ----------- local variables
    integer :: j, k, l
    integer, parameter :: mymarker = 6

    PUSH_SUB(write_band)

    call open_file(unit = 21, file = fname, position = 'rewind', status = 'unknown', &
      form = 'formatted')

    write(21, *) mymarker, ' band info'
    write(21, *) struc%nrk, struc%min, struc%max, struc%nspin
    write(21, *) ((struc%nband(k, j), k = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) ((struc%ifmax(k, j), k = 1, struc%nrk), j = 1, struc%nspin)

    write(21, *) (((struc%energy(l, k, j), l = 1, struc%max), &
      k = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) (((struc%ekn(l, k, j), l = 1, struc%max), &
      k = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) (((struc%occup(l, k, j), l = 1, struc%max), &
      k = 1, struc%nrk), j = 1, struc%nspin)
    call close_file(21)

    POP_SUB(write_band)
    return  

  end subroutine write_band

end module band_m
