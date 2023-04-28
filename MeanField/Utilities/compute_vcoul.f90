!=============================================================================
!
! Utilities:
!
! compute_vocul          Originally By FHJ      Last Modified 07/2019 (FHJ)
!
!     Compute the mini-BZ averaged Coulomb potential.
!
!==============================================================================

#include "f_defs.h"

program compute_vcoul

  use hdf5
  use hdf5_io_m
  use global_m
  use vcoul_generator_m
  use wfn_io_hdf5_m
  use input_utils_m
  use sort_m

  implicit none

  !type(input_parser_t) :: input_parser
  integer(HID_T) :: file_id

  character(len=256) :: fname_in
  type(mf_header_t) :: mf
  integer :: icutv, iscreen
  real(DP) :: truncval(3), q0vec(3), qq(3), avgcut, oneoverq
  integer :: ig
  integer :: iparallel, iwritecoul
  real(DP), allocatable :: vcoul(:)
  type (twork_scell) :: work_scell
  logical :: has_epsmat

  integer :: error
  SCALAR :: epshead, wcoul
  real(DP) :: epshead_tmp(1)
  real(DP), allocatable :: epsheads(:)

  avgcut = TOL_ZERO
  iparallel = 1
  iwritecoul = 0
  truncval(:) = 0d0

  call peinfo_init()
  peinf%jobtypeeval = 1
  call h5open_f(error)
  if (error/=0) stop

  call open_file(file='compute_vcoul.inp', unit=666, form='formatted', status='old')
  read(666,*) fname_in
  read(666,*) icutv
  read(666,*) iscreen
  call close_file(666)

  ! Read mf header from eps0mat
  call read_hdf5_mf_header(trim(fname_in), mf)
  SAFE_ALLOCATE(mf%gvec%components, (3,mf%gvec%ng))
  call read_hdf5_gvectors(trim(fname_in), mf%gvec%ng, mf%gvec%components)
  call gvec_index(mf%gvec)

  if (peinf%inode==0) then
    write(*,*) 'Number of G-vectors:', mf%gvec%ng
    write(*,*) 'Cell volume:', mf%crys%celvol
    write(*,*) 'icutv:', icutv
    write(*,*) 'iscreen:', iscreen
  endif

  if (peinf%inode==0) then
    ! Setup output file
    write(*,*) 'Writing header'
    mf%sheader='MBZ'
    call setup_hdf5_mf_file('vcoul.h5')
    call write_hdf5_mf_header('vcoul.h5', mf)
    call write_hdf5_gvectors('vcoul.h5', mf%gvec%ng, mf%gvec%components)
    call hdf5_open_file('vcoul.h5', 'rw', file_id)
    write(*,*) 'Creating group'
    call hdf5_create_group(file_id, '/vcoul')
  endif

  ! Compute Coulomb potential and write to file
  SAFE_ALLOCATE(vcoul, (mf%gvec%ng))

  if (peinf%inode==0) then
    write(*,*) 'Computing Coulomb potential...'
    write(*,'()')
  endif

  epshead = ONE
  wcoul = ZERO
  call vcoul_generator(icutv, truncval, mf%gvec, mf%crys%bdot, mf%crys%celvol, &
    1, mf%gvec%ng, [(ig, ig=1,mf%gvec%ng)], iscreen, [0d0,0d0,0d0], &
    [0d0,0d0,0d0], vcoul, iwritecoul, iparallel, avgcut, oneoverq, &
    mf%kp%kgrid, epshead, work_scell, .true., wcoul)

  if (peinf%inode==0) then
    write(*,'(1x,a,f0.6)') '- Vcoul(1): ', dble(vcoul(1))
    call hdf5_write_double_array(file_id, '/vcoul/vcoul', [mf%gvec%ng], vcoul)
  endif

  if (peinf%inode==0) then
    write(*,*) 'File "vcoul.h5" written!'
    call hdf5_close_file(file_id)
  endif
  call h5close_f(error)
  if (error/=0) stop

end program compute_vcoul
