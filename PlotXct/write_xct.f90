!==============================================================================
!
! Routines()
!
! (1) write_xct()       Originally By MLT       Last Modified 6/2008 FJR
!
!     Writes out the wavefunction.
!
!=============================================================================

#include "f_defs.h"

module write_xct_m

  use global_m
  use plotxct_common_m

  implicit none

  private

  public :: write_xct

contains

subroutine write_xct(pxct, crys, nfft, scfft)
  type(plotxct_t), intent(in) :: pxct
  type(crystal), intent(in) :: crys
  integer, intent(in) :: nfft(3)
  complex(DPC), intent(inout) :: scfft(:,:,:) !< (nw_eff(1),nw_eff(2),nw_eff(3))
  
  integer :: i1,i2,i3
  integer :: nw(3), nw_int(3), nw_eff(3)
  character :: fn*30
  real(DP) :: rh_a0(3), int_val

  PUSH_SUB(write_xct)
  
  nw(:) = nfft(:) * pxct%nsuper(:)
  nw_eff = (nw + pxct%downsample - 1) / pxct%downsample
  !FHJ: nw_int is the size nw_eff but also taking into consideration the 
  !     integration over the appropriate axes.
  where (pxct%int_dim)
    nw_int = 1
  elsewhere
    nw_int = nw_eff
  endwhere

  !FHJ: note that the scfft array is product(nw_eff(1:3)) big, and we should
  !     never use nw directly!
  rh_a0 = matmul(crys%avec, pxct%rhole) * crys%alat
  
!-----------------------------
! Print output data

  write(fn,'(a,i6.6,a,i1.1,a)') 'xct.', pxct%pstate, '_s', pxct%ispin, '.a3Dr'
  
  write(6,'(1x,3a,3(1x,i0))') 'Writing data set to "', trim(fn), '" with nw =', nw_int(:)
  
  call open_file(30,file=fn,form='formatted',status='replace')
  
  write(30, '(a,i6)') '# ie = ', pxct%pstate
  write(30, '(a,f12.4,a)') '# e  = ', pxct%e_S, ' eV'
  write(30, '(a,3f14.6,a)') '# rh = ', rh_a0, ' Bohr'
  write(30, '(a,i1)') '# spin = ', pxct%ispin
  write(30, '(a)') '#'
  write(30, '(a)') '# unit cell'
  write(30, '(a,3f14.6,a)') '# a1 = ', pxct%nsuper(1) * crys%avec(1:3,1) * crys%alat, ' Bohr'
  write(30, '(a,3f14.6,a)') '# a2 = ', pxct%nsuper(2) * crys%avec(1:3,2) * crys%alat, ' Bohr'
  write(30, '(a,3f14.6,a)') '# a3 = ', pxct%nsuper(3) * crys%avec(1:3,3) * crys%alat, ' Bohr'
  write(30, '(a)') '#'
  ! only write 1/8 = (1/2)^3 of the full data to save disk space
  write(30, '(a,3i8)') '# ni = ', nw_int(:)
  write(30, '(a)') '#'
  if (pxct%only_psi2) then
    write(30, '(a)') '#    abs^2'
  else
    write(30, '(a)') '#    real         imag        abs^2'
  endif

  if(.not.any(pxct%int_dim)) then
    !FHJ : Not integrating dimensions, so nw_int = nw_eff.

    if (pxct%only_psi2) then
      do i3 = 1, nw_eff(3)
        do i2 = 1, nw_eff(2)
          do i1 = 1, nw_eff(1)
            write(30, '(3es13.5)') ABS2(scfft(i1, i2, i3))
          enddo
          write(30,*)
        enddo
        write(30,*)
      enddo
    else
      do i3 = 1, nw_eff(3)
        do i2 = 1, nw_eff(2)
          do i1 = 1, nw_eff(1)
            write(30, '(3es13.5)') scfft(i1, i2, i3),ABS2(scfft(i1, i2, i3))
          enddo
          write(30,*)
        enddo
        write(30,*)
      enddo
    endif

  else
    ! FHJ : Integrating on one/two dimensions
    ! The innermost index must be the one(s) we are integrating!

    if (count(pxct%int_dim)==1) then !integrate 1 dim

      if (pxct%int_dim(1)) then
        do i3 = 1, nw_eff(3)
          do i2 = 1, nw_eff(2)
            int_val = 0.d0
            do i1 = 1, nw_eff(1)
              int_val = int_val + ABS2(scfft(i1,i2,i3))
            enddo
            scfft(1,i2,i3) = int_val
          enddo
        enddo
      else if (pxct%int_dim(2)) then
        do i3 = 1, nw_eff(3)
          do i1 = 1, nw_eff(1)
            int_val = 0.d0
            do i2 = 1, nw_eff(2)
              int_val = int_val + ABS2(scfft(i1,i2,i3))
            enddo
            scfft(i1,1,i3) = int_val
          enddo
        enddo
      else
        do i2 = 1, nw_eff(2)
          do i1 = 1, nw_eff(1)
            int_val = 0.d0
            do i3 = 1, nw_eff(3)
              int_val = int_val + ABS2(scfft(i1,i2,i3))
            enddo
            scfft(i1,i2,1) = int_val
          enddo
        enddo
      endif

    else !integrate 2 dims

      if (.not.pxct%int_dim(3)) then
        do i3 = 1, nw_eff(3)
          int_val = 0.d0
          do i2 = 1, nw_eff(2)
            do i1 = 1, nw_eff(1)
              int_val = int_val + ABS2(scfft(i1,i2,i3))
            enddo
          enddo
          scfft(1,1,i3) = int_val
        enddo
      else if (.not.pxct%int_dim(2)) then
        do i2 = 1, nw_eff(2)
          int_val = 0.d0
          do i3 = 1, nw_eff(3)
            do i1 = 1, nw_eff(1)
              int_val = int_val + ABS2(scfft(i1,i2,i3))
            enddo
          enddo
          scfft(1,i2,1) = int_val
        enddo
      else
        do i1 = 1, nw_eff(1)
          int_val = 0.d0
          do i3 = 1, nw_eff(3)
            do i2 = 1, nw_eff(2)
              int_val = int_val + ABS2(scfft(i1,i2,i3))
            enddo
          enddo
          scfft(i1,1,1) = int_val
        enddo
      endif

    endif !count(pxct%int_dim)==1

    !FHJ: Yes, we are using an array of complex to store doubles, but this is
    ! better than allocating yet another array!
    do i3 = 1, nw_int(3)
      do i2 = 1, nw_int(2)
        do i1 = 1, nw_int(1)
          write(30, '(1es13.5)') dble(scfft(i1, i2, i3))
        enddo
        write(30,*)
      enddo
      write(30,*)
    enddo

  endif
  
  call close_file(30)
  
  POP_SUB(write_xct)
  
  return
end subroutine write_xct

end module write_xct_m
