!=============================================================================
!
! Utilities:
!
! (1) wfnmix_spinor_rho_vxc    Originally By BAB       Last Modified 1/23/2015 (BAB)
!
!     Based on wfnmix_spinor_rho, for extremely narrow usage with spinor testsuite
!     Simply changes the spinor designation of RHO or VXC file.
!
!     The purpose is to debug the spinor functionality of Sigma GPP, or for
!     similar test calculations at the BSE level.
!
!     FHJ: I think this file is wrong. We are not doubling kp%el, kp%occ, etc.
!
!==============================================================================

#include "f_defs.h"

program wfnmix_spinor_rho_vxc

  use blas_m
  use global_m
  use misc_m
  use wfn_rho_vxc_io_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp,kps
  type(wpgen) :: wpg
  type(siginfo) :: sig
  type(gspace) :: gvec

  character(len=80)  :: infile, outfile
  integer :: nargs
  
  character(len=3) :: sheader
  integer :: iflavor
  SCALAR, allocatable :: vxctemp(:,:)

#ifndef CPLX
  call die('Can not use wfnmix_spinor in real case')
#endif

!---------------------------------
! Get file names from command-line arguments

  nargs = command_argument_count()
  
  if (nargs .ne. 3) then
    call die('Usage: wfnmix_spinor_rho_vxc RHO/VXC RHO1/VXC1 RHO2/VXC2')
  endif
  
  call get_command_argument(1,sheader)
  call get_command_argument(2,infile)
  call get_command_argument(3,outfile)

  if ( (sheader .ne. 'RHO') .and. (sheader .ne. 'VXC')) then
    call die('sheader must be RHO or VXC')
  endif

! Open units

  call open_file(unit=7,file=TRUNC(infile),form='unformatted',status='old')
  call open_file(unit=8,file=TRUNC(outfile),form='unformatted',status='replace')

  if (sheader .eq. 'RHO') then
    write(6,*) 'Converting file ',TRUNC(infile),' into spinor charge density file'
  else
    write(6,*) 'Converting file ',TRUNC(infile),' into spinor VXC file'
  endif

! Read/write data

  iflavor = 0
  call read_binary_header_type(7, sheader, iflavor, kp, gvec, syms, crys, warn = .false.)
  if (iflavor.eq.1) then
    call die('RHO/VXC must be type CPLX')
  endif

! Allocating arrays for kps type

  if (kp%nspinor.eq.2) then 
    call die('Input RHO/VXC file already has two spinor components')
  endif

  kps%nspinor=2
  kps%nspin=kp%nspin
  if (kps%nspin.ne.1) then
    call die('nspin not equal to unity')
  endif
  kps%mnband=2*kp%mnband
  kps%nvband=2*kp%nvband
  kps%ncband=2*kp%ncband
  kps%nrk=kp%nrk
  kps%kgrid(:)=kp%kgrid(:)
  kps%shift(:)=kp%shift(:)
  kps%ecutwfc=kp%ecutwfc

  ! Now write new header info
  call write_binary_header_type(8, sheader, iflavor, kps, gvec, syms, crys)

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))

  call read_binary_gvectors(7, gvec%ng, gvec%ng, gvec%components)
  call write_binary_gvectors(8, gvec%ng, gvec%ng, gvec%components)

! Output info

  if (sheader.eq.'RHO') then

    SAFE_ALLOCATE(wpg%rho, (gvec%ng,kp%nspin))

    call read_binary_data(7, gvec%ng, gvec%ng, kp%nspin, wpg%rho)
    call write_binary_data(8, gvec%ng, gvec%ng, kps%nspin, wpg%rho)

    SAFE_DEALLOCATE_P(wpg%rho)

  elseif (sheader.eq.'VXC') then

    SAFE_ALLOCATE(sig%vxc, (gvec%ng,kp%nspin))
    call read_binary_data(7, gvec%ng, gvec%ng, kp%nspin, sig%vxc)

    SAFE_ALLOCATE(vxctemp, (gvec%ng,kps%nspin))
    vxctemp(:,:) = ZERO

    ! from previous call die statements, sig%vxc is guaranteed to only have one spin component
    vxctemp(:,1) = sig%vxc(:,1)

    SAFE_DEALLOCATE_P(sig%vxc)

    call write_binary_data(8, gvec%ng, gvec%ng, kps%nspin, vxctemp)
    SAFE_DEALLOCATE(vxctemp)

  endif

  call close_file(7)
  call close_file(8)

  SAFE_DEALLOCATE_P(gvec%components)
  ! FHJ: The deallotion bellow is wrong and triggers a runtime error with NAG.
  !      Someone should fix this in the future.
  !SAFE_DEALLOCATE_P(kps%ngk)
  !SAFE_DEALLOCATE_P(kps%w)
  !SAFE_DEALLOCATE_P(kps%rk)
  !SAFE_DEALLOCATE_P(kps%ifmin)
  !SAFE_DEALLOCATE_P(kps%ifmax)
  !SAFE_DEALLOCATE_P(kps%el)
  !SAFE_DEALLOCATE_P(kps%occ)

  call dealloc_header_type(sheader, crys, kp)
  !call dealloc_header_type(sheader, crys, kps)
  
  write(6,*) 'Done '
  
end program wfnmix_spinor_rho_vxc
