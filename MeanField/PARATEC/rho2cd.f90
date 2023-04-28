!-------------------------------------------------------------------------------
!
!   rho2cd.f90
!   converts BerkeleyGW file RHO to PARATEC file CD
!   written by M. Jain (May 2010)
!
!-------------------------------------------------------------------------------

#include "f_defs.h"

program rho2cd

  use global_m
  use wfn_rho_vxc_io_m
  implicit none

  character(len=24) :: bdate
  character(len=8) :: btime
  character(len=11) :: cdate
  character(len=14) :: stime

  type(kpoints) :: kp
  type(gspace) :: gvec
  type(symmetry) :: syms
  type(crystal) :: crys

  integer :: ig, is
  SCALAR, allocatable :: rho(:,:)

  integer :: iflavor
  character(len=3) :: sheader
  
  call open_file(95,file='RHO',form='unformatted',status='old')
  call open_file(21,file='CD',form='unformatted',status='replace')
  
  write(6,*) 'Opening RHO'

  sheader = 'RHO'
  iflavor = 0
  call read_binary_header_type(95, sheader, iflavor, kp, gvec, syms, crys, warn = .false.)

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  call read_binary_gvectors(95, gvec%ng, gvec%ng, gvec%components)

  SAFE_ALLOCATE(rho, (gvec%ng, kp%nspin))
  call read_binary_data(95, gvec%ng, gvec%ng, kp%nspin, rho)

  call date_time(cdate, stime)
  write(bdate, '(a11,13x)') cdate
  write(btime, '(a8)') stime(1:8)
  write(21) bdate, btime
  
  write(21) 1, kp%nspin
  
  do ig = 1, gvec%ng
    write(21) gvec%components(1:3, ig)
    write(21) (COMPLEXIFY(rho(ig, is)), is = 1, kp%nspin)
  enddo
  write(21) -1234567, 0, 0
  ! the purpose of this last line is a profound mystery
  
  call dealloc_header_type(sheader, crys, kp)
  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE(rho)

  call close_file(95)
  call close_file(21)
  
end program rho2cd
