!===============================================================================
!
! Utilities:
!
! (1) hdf2wfn         Originally By JIM      Last Modified 4/23/2012 (JIM)
!
!     Converts HDF5 WFN to ascii
!
!===============================================================================

#include "f_defs.h"

program hdf2wfn
  use global_m
  use hdf5
  use io_utils_m
  use wfn_io_hdf5_m
  use wfn_rho_vxc_io_m

  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp
  type(gspace) :: gvec
  character(len=3) :: sheader
  character(len=5) :: file_fmt
  character(len=256) :: infile
  character(len=256) :: outfile
  integer :: ik, ib, iflavor
  real(DP), allocatable :: dwfn(:,:)
  complex(DPC), allocatable :: zwfn(:,:)
  integer, pointer :: gvec_all(:,:)
  integer :: ioffsetk
  integer :: error, narg
  logical :: is_fmt
  type(progress_info) :: prog_info !< a user-friendly progress report

!-------------------------------
! JIM: Open HDF interface
#ifdef HDF5
  call h5open_f(error)
#endif

  narg = command_argument_count()
  if (narg/=3) then
    call die('Usage: hdf2wfn.x ASCII|BIN wfn_in.h5 wfn_out')
  endif
  call get_command_argument(1, file_fmt)
  call get_command_argument(2, infile)
  call get_command_argument(3, outfile)

  if (file_fmt=='ASCII') then
    call open_file(unit=8, file=trim(outfile), form='formatted', status='replace')
    is_fmt = .true.
  elseif (file_fmt=='BIN  ') then
    call open_file(unit=8, file=trim(outfile), form='unformatted', status='replace')
    is_fmt = .false.
  else
    call die('Invalid format "'//trim(file_fmt)//'". Options are ASCII or BIN.')
  endif

  write(*,*) 'Reading/writing header and G vectors from ', trim(infile),' to ', trim(outfile)
  iflavor = -1
  call read_hdf5_header_type(trim(infile), sheader, iflavor, kp, gvec, syms, crys)
  sheader = 'WFN'
  call write_header_type(8, is_fmt, sheader, iflavor, kp, gvec, syms, crys)

  SAFE_ALLOCATE(gvec%components, (3,gvec%ng))
  SAFE_ALLOCATE(gvec_all, (3,sum(kp%ngk)))
  if (iflavor == 1) then
    SAFE_ALLOCATE(dwfn, (kp%ngkmax,kp%nspin*kp%nspinor))
  else
    SAFE_ALLOCATE(zwfn, (kp%ngkmax,kp%nspin*kp%nspinor))
  endif

  call read_hdf5_gvectors(trim(infile), gvec%ng, gvec%components)
  call write_gvectors(8, is_fmt, gvec%ng, gvec%ng, gvec%components)

  write(*,*) 'Reading/writing wavefunctions G vectors from ', trim(infile),' to ', trim(outfile)

  call read_hdf5_wfn_gvectors(trim(infile), gvec_all, sum(kp%ngk))

  write(*,*) 'Reading/writing wavefunctions from ', trim(infile), ' to ', trim(outfile)

  ioffsetk = 0
  call progress_init(prog_info, 'conversion of wavefunctions', 'k-points', kp%nrk)
  do ik = 1, kp%nrk
    call progress_step(prog_info)
    gvec%components(:,:kp%ngk(ik)) = gvec_all(:,ioffsetk+1:ioffsetk+kp%ngk(ik))
    call write_gvectors(8, is_fmt, kp%ngk(ik), kp%ngkmax, gvec%components)
    do ib = 1, kp%mnband
      ! write(*,'(8X,A,I6)') 'band', ib
      if (iflavor == 1) then
        call read_hdf5_band_real(trim(infile), dwfn(1:kp%ngk(ik),:), kp%ngk(ik), &
          kp%nspin*kp%nspinor, ioffsetk, ib-1)
        call write_real_data(8, is_fmt, kp%ngk(ik), kp%ngkmax, kp%nspin, dwfn)
      else
        call read_hdf5_band_complex(trim(infile), zwfn(1:kp%ngk(ik),:), kp%ngk(ik), &
          kp%nspin*kp%nspinor, ioffsetk, ib-1)
        call write_complex_data(8, is_fmt, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, zwfn)
      endif
    enddo
    ioffsetk = ioffsetk + kp%ngk(ik)
  enddo
  call progress_free(prog_info)

  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE_P(gvec_all)
  call dealloc_header_type(sheader, crys, kp)

  if(iflavor == 1) then
    SAFE_DEALLOCATE(dwfn)
  else
    SAFE_DEALLOCATE(zwfn)
  endif

!-------------------------------
! JIM: Close HDF interface
#ifdef HDF5
  call h5close_f(error)
#endif

  call close_file(8)

  write(*,*) 'Conversion from HDF5 to binary file complete!'

end program hdf2wfn
