!=============================================================================
!
! Utilities:
!
! (1) real2complex         Originally By BAB      Last Modified 05/13/2015 (BAB)
!
!     Converts WFN/RHO/VXC files, Real to Complex flavor.
!
!     This program is useful for testing spinor functionality in the testsuite,
!     converting pre-generated real RHO and VXC files to complex flavor.
!
!==============================================================================

#include "f_defs.h"

program real2complex

  use global_m
  use wfn_rho_vxc_io_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp
  type(gspace) :: gvec
  character*3 :: sheader, sformat
  character*7 :: sflavor
  character*11 :: inform, outform
  character*256 :: infile, outfile, usage
  integer :: ik, ib, iflavor, iflavor_new
  integer :: nargs
  real(DP), pointer :: dwfn(:,:)
  complex(DPC), pointer :: zwfn(:,:)
  logical :: informat, outformat

  usage = 'Usage: real2complex.x infile outfile'

! Get file names from command-line arguments

  nargs = command_argument_count()

  if (nargs < 2) then
    call die(usage)
  endif

! mpiexec from mpich1 may add 4 extra arguments to the list
  if (nargs > 2) then
    write(0,'(a,i3,a)') 'WARNING: ', nargs, ' arguments found, only first 2 being used.'
  endif

  call get_command_argument(1, infile)
  call get_command_argument(2, outfile)

! Open units

  informat = .false.
  inform = 'unformatted'
  outformat = .false.
  outform = 'unformatted'

  call open_file(unit=7,file=TRUNC(infile),form=inform,status='old')
  call open_file(unit=8,file=TRUNC(outfile),form=outform,status='replace')

  write(6,'(/,3x,a,/)') 'Converting file ' // TRUNC(infile) // ' from ' // &
    TRUNC(inform) // ' to ' // TRUNC(outform)

  sheader = 'GET'
  iflavor = -1
  call read_header_type(7, informat, sheader, iflavor, kp, gvec, syms, crys, dont_warn_kgrid = .true.)

  if (iflavor .eq. 1) then
    sflavor = RFLAVOR
    iflavor_new = 2
  else
    call die('Input file neads to be REAL flavor')
  endif

! Output info

  write(6,'(3x,"File header:",17x,a,/)') sheader // '-' // TRUNC(sflavor)
  write(6,'(3x,"Crystal volume:",f32.14)') crys%celvol
  write(6,'(3x,"Number of G-vectors:",i12)') gvec%ng
  write(6,'(3x,"Number of spins:",i16)') kp%nspin
  if (sheader .eq. 'WFN') then
    write(6,'(3x,"Number of bands:",i16)') kp%mnband
    write(6,'(3x,"Number of k-points:",i13)') kp%nrk
  endif
  write(6,*)

  call write_header_type(8, outformat, sheader, iflavor_new, kp, gvec, syms, crys)

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))

  call read_gvectors(7, informat, gvec%ng, gvec%ng, gvec%components)
  call write_gvectors(8, outformat, gvec%ng, gvec%ng, gvec%components)

  if (sheader .eq. 'WFN') then
    SAFE_ALLOCATE(dwfn, (kp%ngkmax, kp%nspin))
    SAFE_ALLOCATE(zwfn, (kp%ngkmax, kp%nspin))
  else
    SAFE_ALLOCATE(dwfn, (gvec%ng, kp%nspin))
    SAFE_ALLOCATE(zwfn, (gvec%ng, kp%nspin))
  endif

  if (sheader .eq. 'WFN') then
    do ik = 1, kp%nrk
      call read_gvectors(7, informat, kp%ngk(ik), kp%ngkmax, gvec%components)
      call write_gvectors(8, outformat, kp%ngk(ik), kp%ngkmax, gvec%components)

      do ib = 1, kp%mnband
        call read_real_data(7, informat, kp%ngk(ik), kp%ngkmax, kp%nspin, dwfn)
        zwfn=COMPLEXIFY(dwfn)
        call write_complex_data(8, outformat, kp%ngk(ik), kp%ngkmax, kp%nspin, zwfn)
      enddo
    enddo
  else
    call read_real_data(7, informat, gvec%ng, gvec%ng, kp%nspin, dwfn)
    zwfn=COMPLEXIFY(dwfn)
    call write_complex_data(8, outformat, gvec%ng, gvec%ng, kp%nspin, zwfn)
  endif

  SAFE_DEALLOCATE_P(dwfn)
  SAFE_DEALLOCATE_P(zwfn)

  SAFE_DEALLOCATE_P(gvec%components)

  call dealloc_header_type(sheader, crys, kp)

  call close_file(7)
  call close_file(8)

  write(6,'(3x,"Done",/)')

end program real2complex
