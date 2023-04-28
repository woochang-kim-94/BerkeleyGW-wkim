!===============================================================================
!
! extract_vcoul:
!
! Computes the v(q+G) using the same mBZ averaging as BerkeleyGW.  (FHJ/2020)
!
! No input parameters are required. However:
! 1) You need at least an eps0mat.h5 file in your directory, from which
!    we read the q points and the G sphere. An epsmat.h5 may also be present.
! 2) We hardcode the screening model to that of a semiconductor.
!
!===============================================================================

#include "f_defs.h"

program extract_vcoul

  use hdf5
  use hdf5_io_m
  use global_m
  use vcoul_generator_m
  use wfn_io_hdf5_m
  use input_utils_m
  use sort_m

  implicit none

  !type(input_parser_t) :: input_parser
  integer(HID_T) :: file_id, group_id, dataspace_id, dataset_id
  integer(HSIZE_T) :: dset_sz(2)

  type(mf_header_t) :: mf
  integer :: icutv, iscreen, ng
  real(DP) :: truncval(3), q0vec(3), qq(3), avgcut, oneoverq
  integer :: nf, iq, nq, qgrid(3)
  integer :: iparallel, iwritecoul
  real(DP), allocatable :: vcoul(:), qpts(:,:)
  real(DP) :: g2_max, cutoff
  integer, allocatable :: isort(:)
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

  ! Read mf header from eps0mat
  call read_hdf5_mf_header('eps0mat.h5', mf)
  SAFE_ALLOCATE(mf%gvec%components, (3,mf%gvec%ng))
  call read_hdf5_gvectors('eps0mat.h5', mf%gvec%ng, mf%gvec%components)
  call gvec_index(mf%gvec)

  ! Figure out the G-sphere. We use the wavefunction cutoff
  ! plus 2 G-vectors for the cutoff.
  associate(bdot=>mf%crys%bdot, gvec=>mf%gvec)
    ! Maximum lenghts of the reciprocal lattice vectors
    g2_max = max(max(bdot(1,1), bdot(2,2)), bdot(3,3))
    ! Factor of 2d0 is for the extra padding
    cutoff = (sqrt(gvec%ecutrho/4d0) + sqrt(g2_max)*2d0)**2
    ! And now find the cutoff
    SAFE_ALLOCATE(isort, (gvec%ng))
    SAFE_ALLOCATE(gvec%ekin, (gvec%ng))
    call kinetic_energies(gvec, bdot, gvec%ekin)
    call sortrx(gvec%ng, gvec%ekin, isort, gvec=gvec%components)
    ng = gcutoff(gvec%ng, gvec%ekin, isort, cutoff)
  endassociate
  if (peinf%inode==0) then
    write(*,'(1x,a,f0.2,a)') 'Charge-density cutoff: ', mf%gvec%ecutrho, ' Ry'
    write(*,'(1x,a,f0.2,a)') 'WFN cutoff: ', mf%gvec%ecutrho/4d0, ' Ry'
    write(*,'(1x,a,f0.2,a)') 'v(q+G) cutoff: ', cutoff, ' Ry'
  endif

  ! Read q-points from eps0mat
  call h5fopen_f('eps0mat.h5', H5F_ACC_RDONLY_F, file_id, error)
  if (error/=0) stop
  call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq, error)
  if (error/=0) stop
  if (nq/=1) call die('Only support eps0mat.h5 with a single q-point.')
  call hdf5_read_double_array(file_id, 'eps_header/qpoints/qpts', [3,1], q0vec, error)
  if (error/=0) stop
  ! Note: we shouldn`t trust the qgrid from eps0mat.h5!
  call hdf5_read_int_array(file_id, 'eps_header/qpoints/qgrid', [3], qgrid, error)
  if (error/=0) stop
  call hdf5_read_int(file_id, 'eps_header/params/icutv', icutv, error)
  if (error/=0) stop
  call hdf5_read_double_hyperslab(file_id, 'mats/matrix-diagonal', [1,1,1], &
    [0,0,0], epshead_tmp, error)
  if (error/=0) stop
  ! FIXME
  iscreen = SCREEN_SEMICOND
  nf = product(qgrid)
  call h5fclose_f(file_id, error)
  if (error/=0) stop

  inquire(file='epsmat.h5', exist=has_epsmat)
  ! Read q-points from epsmat
  if (has_epsmat) then
    call h5fopen_f('epsmat.h5', H5F_ACC_RDONLY_F, file_id, error)
    if (error/=0) stop
    call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq, error)
    if (error/=0) stop
    nq = nq + 1
    SAFE_ALLOCATE(qpts, (3,nq))
    SAFE_ALLOCATE(epsheads, (nq))
    epsheads(1) = epshead_tmp(1)
    do iq = 2, nq
      call hdf5_read_double_hyperslab(file_id, 'mats/matrix-diagonal', [1,1,1], &
        [0,0,iq-2], epshead_tmp, error)
      epsheads(iq) = epshead_tmp(1)
    enddo
    call hdf5_read_double_array(file_id, 'eps_header/qpoints/qpts', [3,nq-1], qpts(:,2:nq), error)
    if (error/=0) stop
    call hdf5_read_int_array(file_id, 'eps_header/qpoints/qgrid', [3], qgrid, error)
    if (error/=0) stop
    call h5fclose_f(file_id, error)
    if (error/=0) stop
  else
    nq = 1
    SAFE_ALLOCATE(qpts, (3,nq))
    SAFE_ALLOCATE(epsheads, (nq))
    epsheads(1) = epshead_tmp(1)
  endif
  nf = product(qgrid)
  qpts(:,1) = 0d0
  if (peinf%inode==0) then
    write(*,*) 'Number of G-vectors:', ng
    write(*,*) 'Number of q-points in the full BZ:', nf
    write(*,*) 'Number of q-points:', nq
    write(*,*) 'Q grid:', qgrid
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
    call h5fopen_f('vcoul.h5', H5F_ACC_RDWR_F, file_id, error)
    if (error/=0) stop

    write(*,*) 'Creating group'
    call h5gcreate_f(file_id, '/vcoul', group_id, error)
    if (error/=0) stop
    call h5gclose_f(group_id, error)
    if (error/=0) stop
    write(*,*) 'Writing isort'
    call hdf5_write_int_array(file_id, '/vcoul/isort', [ng], isort(1:ng), error)
    if (error/=0) stop
    write(*,*) 'Writing nq'
    call hdf5_write_int(file_id, '/vcoul/nq', nq, error)
    if (error/=0) stop
    write(*,*) 'Writing q0vec'
    call hdf5_write_double_array(file_id, '/vcoul/q0vec', [3], q0vec, error)
    if (error/=0) stop
    write(*,*) 'Writing qpts'
    call hdf5_write_double_array(file_id, '/vcoul/qpts', [3,nq], qpts, error)
    if (error/=0) stop

    write(*,*) 'Preparing dataset vcoul'
    dset_sz(1) = ng
    dset_sz(2) = nq
    call h5screate_simple_f(2, dset_sz(:2), dataspace_id, error)
    if (error/=0) stop
    call h5dcreate_f(file_id, '/vcoul/vcoul', H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
    if (error/=0) stop
    call h5dclose_f(dataset_id, error)
    if (error/=0) stop
    call h5sclose_f(dataspace_id, error)
    if (error/=0) stop

    write(*,*) 'Preparing dataset wcoul'
    dset_sz(1) = nq
    call h5screate_simple_f(1, dset_sz(:1), dataspace_id, error)
    if (error/=0) stop
    call h5dcreate_f(file_id, '/vcoul/wcoul', H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
    if (error/=0) stop
    call h5dclose_f(dataset_id, error)
    if (error/=0) stop
    call h5sclose_f(dataspace_id, error)
    if (error/=0) stop

    write(*,*) 'Preparing dataset oneoverq'
    dset_sz(1) = nq
    call h5screate_simple_f(1, dset_sz(:1), dataspace_id, error)
    if (error/=0) stop
    call h5dcreate_f(file_id, '/vcoul/oneoverq', H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
    if (error/=0) stop
    call h5dclose_f(dataset_id, error)
    if (error/=0) stop
    call h5sclose_f(dataspace_id, error)
    if (error/=0) stop
  endif

  ! Compute Coulomb potential and write to file
  SAFE_ALLOCATE(vcoul, (ng))

  if (peinf%inode==0) then
    write(*,*) 'Starting the loop...'
    write(*,'()')
  endif

  do iq = 1, nq
    qq = qpts(:,iq)
    epshead = epsheads(iq) * ONE
    if (peinf%inode==0) then
      write(*,'(1x,a,i0,a,i0,a,3(1x,f9.6))') 'Dealing with iq=', iq, '/', nq, ' ; q=', qq
      write(*,'(1x,a,f0.6)') '- Head of epsilon: ', epsheads(iq)
    endif

    wcoul = ZERO
    call vcoul_generator(icutv, truncval, mf%gvec, &
      mf%crys%bdot, mf%crys%celvol, nf, ng, isort, iscreen, qq, q0vec, &
      vcoul, iwritecoul, iparallel, avgcut, oneoverq, &
      qgrid, epshead, work_scell, .true., wcoul)

    if (peinf%inode==0) then
      write(*,'(1x,a,f0.6)') '- Vcoul(1): ', dble(vcoul(1))
      write(*,'(1x,a,f0.6)') '- Wcoul: ', dble(wcoul)
      write(*,'()')
      call hdf5_write_double_hyperslab(file_id, '/vcoul/vcoul', &
        [ng,1], [0,iq-1], vcoul, error)
      if (error/=0) stop
      call hdf5_write_double_hyperslab(file_id, '/vcoul/wcoul', &
        [1], [iq-1], [dble(wcoul)], error)
      if (error/=0) stop
      call hdf5_write_double_hyperslab(file_id, '/vcoul/oneoverq', &
        [1], [iq-1], [oneoverq], error)
      if (error/=0) stop
    endif
  enddo

  if (peinf%inode==0) then
    write(*,*) 'File "vcoul.h5" written!'
    call h5fclose_f(file_id, error)
    if (error/=0) stop
  endif
  call h5close_f(error)
  if (error/=0) stop

end program extract_vcoul
