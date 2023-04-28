!===============================================================================
!
! Utilities:
!
! (1) wfn2hdf         Originally By JIM      Last Modified 4/23/2012 (JIM)
!
!     Converts binary WFN to HDF5
!
!===============================================================================

#include "f_defs.h"

program wfn2hdf
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
  integer :: ioffsetk
  integer, allocatable :: wfn_gvec_all(:,:)
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
    call die('Usage: wfn2hdf.x ASCII|BIN wfn_in wfn_out.h5')
  endif
  call get_command_argument(1, file_fmt)
  call get_command_argument(2, infile)
  call get_command_argument(3, outfile)

  if (file_fmt=='ASCII' .or. file_fmt=='ascii') then
    call open_file(unit=7, file=trim(infile), form='formatted', status='old')
    is_fmt = .true.
  elseif (file_fmt=='BIN  ' .or. file_fmt=='bin ') then
    call open_file(unit=7, file=trim(infile), form='unformatted', status='old')
    is_fmt = .false.
  else
    call die('Invalid format "'//trim(file_fmt)//'". Options are ASCII or BIN.')
  endif

  write(*,*) 'Reading header and gvectors from ', trim(infile)

  sheader = 'GET'
  iflavor = -1
  call read_header_type(7, is_fmt, sheader, iflavor, kp, gvec, syms, crys, dont_warn_kgrid=.true.)
  if (sheader .ne. 'WFN') call die("wfn2hdf only works with WFN files")

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  SAFE_ALLOCATE(wfn_gvec_all, (3,sum(kp%ngk)))

  call read_gvectors(7, is_fmt, gvec%ng, gvec%ng, gvec%components)

  if (iflavor == 1) then
    SAFE_ALLOCATE(dwfn, (kp%ngkmax,kp%nspin*kp%nspinor))
  else
    SAFE_ALLOCATE(zwfn, (kp%ngkmax,kp%nspin*kp%nspinor))
  endif

  call setup_hdf5_wfn_file(trim(outfile), iflavor, kp)

  write(*,*) 'Writing header and gvectors to ', trim(outfile)

  call write_hdf5_header_type(trim(outfile), sheader, iflavor, kp, gvec, syms, crys)
  call write_hdf5_gvectors(trim(outfile), gvec%ng, gvec%components)

  write(*,*) 'Reading/writing wavefunctions from ', trim(infile), ' to ', trim(outfile)

  ioffsetk = 0
  call progress_init(prog_info, 'conversion of wavefunctions', 'k-points', kp%nrk)
  do ik = 1, kp%nrk
    call progress_step(prog_info)
    call read_gvectors(7, is_fmt, kp%ngk(ik), kp%ngkmax, gvec%components)
    wfn_gvec_all(:, ioffsetk+1:ioffsetk+kp%ngk(ik))=gvec%components(:,1:kp%ngk(ik))
    do ib = 1, kp%mnband
      ! write(*,'(8X,A,I6)') 'band', ib
      if (iflavor == 1) then
        call read_real_data(7, is_fmt, kp%ngk(ik), kp%ngkmax, kp%nspin, dwfn)
        call write_hdf5_band_real(trim(outfile), dwfn(1:kp%ngk(ik),:), kp%ngk(ik), kp%nspin,ioffsetk, ib-1)
      else
        call read_complex_data(7, is_fmt, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, zwfn)
        call write_hdf5_band_complex(trim(outfile), zwfn(1:kp%ngk(ik),:), kp%ngk(ik), kp%nspin*kp%nspinor,ioffsetk,ib-1)
      endif
    enddo
    ioffsetk = ioffsetk + kp%ngk(ik)
  enddo
  call progress_free(prog_info)

  call write_hdf5_wfn_gvectors(trim(outfile), wfn_gvec_all, sum(kp%ngk))

  call dealloc_header_type(sheader, crys, kp)
  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE(wfn_gvec_all)

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

  call close_file(7)
  write(*,*) 'Conversion from binary to HDF5 file complete!'

end program wfn2hdf
