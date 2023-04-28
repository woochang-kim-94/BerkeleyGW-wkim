!=============================================================================
!
! Utilities:
!
! (1) mf_convert         Originally By DAS      Last Modified 10/17/2010 (DAS)
!
!     Converts WFN/RHO/VXC files, Real/Complex flavor, binary <--> ascii
!
!==============================================================================

#include "f_defs.h"

program mf_convert

  use global_m
  use wfn_rho_vxc_io_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp
  type(gspace) :: gvec
  character*3 :: sheader, sformat
  character*7 :: sflavor
  character*11 :: inform, outform, inaccess, outaccess
  character*256 :: infile, outfile, usage
  integer :: ik, ib, iflavor
  integer :: nargs
  real(DP), pointer :: dwfn(:,:)
  complex(DPC), pointer :: zwfn(:,:)
  logical :: informat, outformat

  usage = 'Usage: mf_convert.x A2B|B2A infile outfile'

! Get file names from command-line arguments

  nargs = command_argument_count()

  if (nargs < 3) then
    call die(usage)
  endif

! mpiexec from mpich1 may add 4 extra arguments to the list
  if (nargs > 3) then
    write(0,'(a,i3,a)') 'WARNING: ', nargs, ' arguments found, only first 3 being used.'
  endif

  call get_command_argument(1, sformat)
  call get_command_argument(2, infile)
  call get_command_argument(3, outfile)

! Open units

  if (sformat .eq. 'A2B') then
    informat = .true.
    inform = 'formatted'
    inaccess = 'stream'
    outformat = .false.
    outform = 'unformatted'
    outaccess = 'sequential'
  elseif (sformat .eq. 'B2A') then
    informat = .false.
    inform = 'unformatted'
    inaccess = 'sequential'
    outformat = .true.
    outform = 'formatted'
    outaccess = 'stream'
  else  
    call die(usage)
  endif

  call open_file(unit=7,file=TRUNC(infile),form=inform,status='old',access=inaccess)
  call open_file(unit=8,file=TRUNC(outfile),form=outform,status='replace',access=outaccess)

  write(6,'(/,3x,a,/)') 'Converting file ' // TRUNC(infile) // ' from ' // &
    TRUNC(inform) // ' to ' // TRUNC(outform)

  sheader = 'GET'
  iflavor = -1
  call read_header_type(7, informat, sheader, iflavor, kp, gvec, syms, crys, dont_warn_kgrid = .true.)

  if (iflavor .eq. 1) then
    sflavor = RFLAVOR
  else
    sflavor = CFLAVOR
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

  call write_header_type(8, outformat, sheader, iflavor, kp, gvec, syms, crys)

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))

  call read_gvectors(7, informat, gvec%ng, gvec%ng, gvec%components)
  call write_gvectors(8, outformat, gvec%ng, gvec%ng, gvec%components)

  if(iflavor == 1) then
    if (sheader .eq. 'WFN') then
      SAFE_ALLOCATE(dwfn, (kp%ngkmax, kp%nspin))
    else
      SAFE_ALLOCATE(dwfn, (gvec%ng, kp%nspin))
    endif
  else
    if (sheader .eq. 'WFN') then
      SAFE_ALLOCATE(zwfn, (kp%ngkmax, kp%nspin*kp%nspinor))
    else
      SAFE_ALLOCATE(zwfn, (gvec%ng, kp%nspin))
    endif
  endif

  if (sheader .eq. 'WFN') then
    do ik = 1, kp%nrk
      call read_gvectors(7, informat, kp%ngk(ik), kp%ngkmax, gvec%components)
      call write_gvectors(8, outformat, kp%ngk(ik), kp%ngkmax, gvec%components)

      do ib = 1, kp%mnband
        if(iflavor == 1) then
          call read_real_data(7, informat, kp%ngk(ik), kp%ngkmax, kp%nspin, dwfn)
          call write_real_data(8, outformat, kp%ngk(ik), kp%ngkmax, kp%nspin, dwfn)
        else
          call read_complex_data(7, informat, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, zwfn)
          call write_complex_data(8, outformat, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, zwfn)
        endif
      enddo
    enddo
  else
    if(iflavor == 1) then
      call read_real_data(7, informat, gvec%ng, gvec%ng, kp%nspin, dwfn)
      call write_real_data(8, outformat, gvec%ng, gvec%ng, kp%nspin, dwfn)
    else
      call read_complex_data(7, informat, gvec%ng, gvec%ng, kp%nspin, zwfn)
      call write_complex_data(8, outformat, gvec%ng, gvec%ng, kp%nspin, zwfn)
    endif
  endif

  if(iflavor == 1) then
    SAFE_DEALLOCATE_P(dwfn)
  else
    SAFE_DEALLOCATE_P(zwfn)
  endif

  SAFE_DEALLOCATE_P(gvec%components)

  call dealloc_header_type(sheader, crys, kp)

  call close_file(7)
  call close_file(8)

  write(6,'(3x,"Done",/)')

end program mf_convert
