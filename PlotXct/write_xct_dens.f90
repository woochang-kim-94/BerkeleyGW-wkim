!==============================================================================
!
! Routines()
!
! (1) write_xct_dens()       By PWD         Last Modified PWD 4/19/2014
!        from write_xct() Originally By MLT       
!
!     Writes excited states density.
!
!=============================================================================

#include "f_defs.h"

module write_xct_dens_m

  use global_m
  use plotxct_common_m

  implicit none

  private

  public :: write_xct_dens

contains

subroutine write_xct_dens(pxct, crys, nfft, typeS, scfft)
  type(plotxct_t), intent(in) :: pxct
  type(crystal), intent(in) :: crys
  integer, intent(in) :: nfft(3)
  character (len=*), intent(in) :: typeS
  complex(DPC), intent(inout) :: scfft(:,:,:,:) !< (xctstate,nw_eff(1),nw_eff(2) ! nw_eff(3))
  
  integer :: ixct,i1,i2,i3
  integer :: nw(3) 
  character :: fn*30
  real(DP) :: rh_a0(3), int_val

  PUSH_SUB(write_xct_dens)
  
  nw(:) = nfft(:) * pxct%nsuper(:)
  !nw_eff = (nw + pxct%downsample - 1) / pxct%downsample

  
!-----------------------------
! Print output data

  do ixct = pxct%pstateLow-pxct%offset,pxct%pstateHigh-pxct%offset

  write(fn,'(a,i6.6,a,i1.1,a,a,a)') 'xct.', ixct+pxct%offset, '_s', pxct%ispin, '_t', typeS, '.a3Dr'
  
  write(6,'(a,a,a,3i8)') 'Writing data set to "', trim(fn), '" with nw = ',nw(:)
  
  call open_file(30,file=fn,form='formatted',status='replace')
  
  write(30, '(a,i6)') '# ie = ', ixct+pxct%offset
  write(30, '(a,f12.4,a)') '# e  = ', pxct%e_Ses(ixct), ' eV'
  write(30, '(a,3f14.6,a)') '# rh = ', 0,0,0, ' Bohr'
  write(30, '(a,i1)') '# spin = ', pxct%ispin
  write(30, '(a)') '#'
  write(30, '(a)') '# unit cell'
  write(30, '(a,3f14.6,a)') '# a1 = ', pxct%nsuper(1) * crys%avec(1:3,1) * crys%alat, ' Bohr'
  write(30, '(a,3f14.6,a)') '# a2 = ', pxct%nsuper(2) * crys%avec(1:3,2) * crys%alat, ' Bohr'
  write(30, '(a,3f14.6,a)') '# a3 = ', pxct%nsuper(3) * crys%avec(1:3,3) * crys%alat, ' Bohr'
  write(30, '(a)') '#'
  ! only write 1/8 = (1/2)^3 of the full data to save disk space
  write(30, '(a,3i8)') '# ni = ', nw(:)
  write(30, '(a)') '#'
  write(30, '(a)') '#    real         imag        abs^2'

  if(.not.any(pxct%int_dim)) then
    !FHJ : Not integrating dimensions, so nw_int = nw_eff.

      do i3 = 1, nw(3)
        do i2 = 1, nw(2)
          do i1 = 1, nw(1)
            write(30, '(3es13.5)') scfft(ixct,i1, i2, i3),ABS2(scfft(ixct,i1, i2, i3))
          enddo
          write(30,*)
        enddo
        write(30,*)
      enddo

  endif
  
  call close_file(30)
  end do

  POP_SUB(write_xct_dens)
  
  return
end subroutine write_xct_dens

end module write_xct_dens_m
