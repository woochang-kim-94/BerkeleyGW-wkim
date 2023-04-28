#include "f_defs.h"

module wfn_io_m

  use global_m
  use misc_m,            only: findvector
  use check_inversion_m, only: check_inversion
  use input_utils_m,     only: gvec_index
  use wfn_rho_vxc_io_m,  only: read_mf_header, read_binary_gvectors, read_binary_data
#ifdef HDF5
  use hdf5
  use hdf5_io_m
  use wfn_io_hdf5_m,     only: setup_hdf5_wfn_file, write_hdf5_mf_header, &
                               write_hdf5_gvectors, write_hdf5_wfn_gvectors
#endif
  use inread_m,          only: pb_params_t
  use distribution_m,    only: distrib_mat_t
  use iteration_data_m,  only: iteration_data_t

  implicit none

  private

  public :: &
#ifdef HDF5
    wfn_write_header, &
    wfn_write_kpt, &
    is_kpt_done, &
    set_kpt_done, &
#endif
    wfn_read


contains


#ifdef HDF5
logical function is_kpt_done(fname, ik)
  character(len=*), intent(in) :: fname
  integer, intent(in) :: ik

  integer(HID_T) :: file_id
  integer :: nrk
  logical, allocatable :: kpt_done(:)

  PUSH_SUB(is_kpt_done)

  call hdf5_open_file(trim(fname), 'r', file_id)
  call hdf5_read_int(file_id, 'mf_header/kpoints/nrk', nrk)
  SAFE_ALLOCATE(kpt_done, (nrk))
  call hdf5_read_logical_array(file_id, '/parabands/kpt_done', (/nrk/), kpt_done)
  is_kpt_done = kpt_done(ik)
  call hdf5_close_file(file_id)
  SAFE_DEALLOCATE(kpt_done)

  POP_SUB(is_kpt_done)

end function is_kpt_done


subroutine set_kpt_done(fname, ik)
  character(len=*), intent(in) :: fname
  integer, intent(in) :: ik

  integer(HID_T) :: file_id
  integer :: nrk
  logical, allocatable :: kpt_done(:)

  PUSH_SUB(set_qpt_done)

  call hdf5_open_file(trim(fname), 'rw', file_id)
  call hdf5_read_int(file_id, 'mf_header/kpoints/nrk', nrk)
  SAFE_ALLOCATE(kpt_done, (nrk))
  call hdf5_read_logical_array(file_id, '/parabands/kpt_done', (/nrk/), kpt_done)
  kpt_done(ik) = .true.
  call hdf5_write_logical_array(file_id, '/parabands/kpt_done', (/nrk/), kpt_done)
  call hdf5_close_file(file_id)
  SAFE_DEALLOCATE(kpt_done)


  POP_SUB(set_kpt_done)

end subroutine set_kpt_done
#endif


subroutine wfn_read(mf, params, isort)
  type(mf_header_t), intent(out) :: mf
  type(pb_params_t), intent(inout) :: params
  integer, allocatable, intent(inout) :: isort(:,:)

  integer, allocatable :: gvectmp(:,:)
  SCALAR, allocatable :: wfn_buf(:,:)
  integer :: ig, ib, ik, iout

  PUSH_SUB(wfn_read)

!------------------------
! Read header of input wavefunction file
  if (peinf%inode.eq.0) then
    call open_file(7, file=params%fname_wfn_in, status='old', form='unformatted')
  endif ! peinf%inode.eq.0
  call read_mf_header(7, mf, iflavor=SCALARSIZE, sheader='WFN', warn=.false., dont_warn_kgrid=.true.)
  call check_inversion(SCALARSIZE, mf%syms%ntran, mf%syms%mtrx, &
    mf%kp%nspin, .true., .true., tnp=mf%syms%tnp)

!------------------------
! Read master list of G-vectors
  if (peinf%inode==0) then
    write(6,'(1x,a)') 'Reading wavefunctions from file '//TRUNC(params%fname_wfn_in)
  endif

  SAFE_ALLOCATE(mf%gvec%components, (3, mf%gvec%ng))
  call read_binary_gvectors(7, mf%gvec%ng, mf%gvec%ng, mf%gvec%components)

!------------------------
! Allocate and compute index_vec indices
  call gvec_index(mf%gvec)

!------------------------
! Loop over k-points and read G-vectors/WFNs
  SAFE_ALLOCATE(gvectmp, (3,mf%kp%ngkmax))
  SAFE_ALLOCATE(wfn_buf, (mf%kp%ngkmax,mf%kp%nspin*mf%kp%nspinor))
  SAFE_ALLOCATE(isort, (mf%kp%ngkmax,mf%kp%nrk))
  isort(:,:) = 0
  do ik=1,mf%kp%nrk
    ! Read g-vectors for current k-point, determine their indices in list of
    ! G-vectors, and distribute index arrays over processors
    call read_binary_gvectors(7, mf%kp%ngk(ik), mf%kp%ngkmax, gvectmp)
    do ig = 1, mf%kp%ngk(ik)
      call findvector(iout, gvectmp(:,ig), mf%gvec)
      isort(ig,ik) = iout
    enddo
    if (any(isort(1:mf%kp%ngk(ik),ik)==0)) then
      if (peinf%inode==0) write(0,*) ' ik=',ik
      call die('failed to find G-vector for k-point', only_root_writes=.true.)
    endif
    if (ik<mf%kp%nrk) then
      do ib = 1, mf%kp%mnband
        call read_binary_data(7, mf%kp%ngk(ik), mf%kp%ngkmax, mf%kp%nspin*mf%kp%nspinor, &
          wfn_buf, dont_read=.true.)
      enddo
    endif
  enddo ! ik
  if (peinf%inode==0) then
    call close_file(7)
  endif
  SAFE_DEALLOCATE(gvectmp)
  SAFE_DEALLOCATE(wfn_buf)

  POP_SUB(wfn_read)

end subroutine wfn_read


#ifdef HDF5
! This routine is actually executed in serial
! If isort is not allocated, we will generate it on-the-fly from
! the kpoint_grid derived type.
subroutine wfn_write_header(mf, params, isort)
  type(mf_header_t), intent(in) :: mf
  type(pb_params_t), intent(in) :: params
  integer, intent(in), allocatable :: isort(:,:)

  integer :: gvecs_all(3,sum(mf%kp%ngk))
  integer :: offset_gk, ngk, ik
  integer :: isortk(mf%gvec%ng)
  integer(HID_T) :: file_id
  logical :: kpt_done(mf%kp%nrk)

  if (peinf%inode==0) then
    write(6,'(/1x,3a/)') 'Setting up wavefunctions file ', &
      TRUNC(params%fname_wfn_out), ' in HDF5 format.'

    call setup_hdf5_wfn_file(TRUNC(params%fname_wfn_out), SCALARSIZE, mf%kp)
    call write_hdf5_mf_header(TRUNC(params%fname_wfn_out), mf)
    call write_hdf5_gvectors(TRUNC(params%fname_wfn_out), mf%gvec%ng, mf%gvec%components)
    offset_gk = 0
    gvecs_all = 0
    do ik = 1, mf%kp%nrk
      ngk = mf%kp%ngk(ik)
      if (allocated(isort)) then
        isortk(1:ngk) = isort(1:ngk,ik)
      else
        call die('isort not allocated, but not using kpg')
      endif
      gvecs_all(1:3, offset_gk+1:offset_gk+ngk) = mf%gvec%components(1:3, isortk(1:ngk))
      offset_gk = offset_gk + ngk
    enddo
    call write_hdf5_wfn_gvectors(TRUNC(params%fname_wfn_out), gvecs_all, sum(mf%kp%ngk))

    call hdf5_open_file(trim(params%fname_wfn_out), 'rw', file_id)
    kpt_done(:) = .false.
    call hdf5_create_group(file_id, '/parabands')
    call hdf5_write_logical_array(file_id, '/parabands/kpt_done', &
      [mf%kp%nrk], kpt_done)
    call hdf5_close_file(file_id)
  endif

end subroutine wfn_write_header


subroutine wfn_write_kpt(mf, params, iter_data, dm_wfn, ik, wfn_d, en)
  type(mf_header_t), intent(in) :: mf
  type(pb_params_t), intent(in) :: params
  type(iteration_data_t), intent(in) :: iter_data
  class(distrib_mat_t), intent(in) :: dm_wfn
  integer, intent(in) :: ik
  SCALAR, intent(in) :: wfn_d(:,:,:)
  real(DP), intent(in) :: en(:,:)

  ! Need this barrier otherwise MPI ranks may get here before rank==0 sets up
  ! the output WFN file.
#ifdef MPI
  call MPI_Barrier(iter_data%working_comm, mpierr)
#endif
  call wfn_write_kpt_direct(mf, params, iter_data, dm_wfn, ik, wfn_d, en)

end subroutine wfn_write_kpt


subroutine wfn_write_kpt_direct(mf, params, iter_data, dm_wfn, ik, wfn_d, en)
  type(mf_header_t), intent(in) :: mf
  type(pb_params_t), intent(in) :: params
  type(iteration_data_t), intent(in) :: iter_data
  class(distrib_mat_t), intent(in) :: dm_wfn
  integer, intent(in) :: ik
  SCALAR, intent(in) :: wfn_d(:,:,:)
  real(DP), intent(in) :: en(:,:)

  integer(HID_T) :: file_id
#ifdef MPI
  integer(HID_T) :: plist_id
#endif
  integer(HID_T) :: dset_id
  integer(HID_T) :: filespace
  integer(HID_T) :: memspace
  integer(HSIZE_T) :: count(4), offsetf(4)
  integer :: count_(4)
  integer :: error
  real(DP), pointer, contiguous :: wfn_d2(:,:,:,:)
  integer :: offsetk
  integer*8 :: clock_count
  real(DP) :: clock_inv_rate, t0, t1
  integer :: is

  if (iter_data%is_working_leader) then
    write(params%kpp%iunit,'(1x,3a)') &
      'Writing wavefunctions to file ', trim(params%fname_wfn_out), '.'
    FLUSH(params%kpp%iunit)
  endif
  call system_clock(count_rate=clock_count)
  clock_inv_rate = 1d0/clock_count
  call system_clock(count=clock_count)
  t0 = clock_count * clock_inv_rate


  ! Write wavefunctions
  !--------------------
  if (iter_data%is_working) then
    call hdf5_open_file(trim(params%fname_wfn_out), 'rw', file_id, &
      parallel_io=.true., comm=iter_data%working_comm)

    call safe_h5dopen(file_id, 'wfns/coeffs', dset_id)
    call safe_h5dget_space(dset_id, filespace)
    count(1) = SCALARSIZE
    count(2) = dm_wfn%Mown/mf%kp%nspinor
    count(3) = mf%kp%nspin*mf%kp%nspinor
    count(4) = dm_wfn%Nown
    call safe_h5screate_simple(4, count, memspace)

    ! FHJ: offset the file data depending on the k-point
    offsetk = 0
    if (ik>1) offsetk = sum(mf%kp%ngk(1:ik-1))
    offsetf(:) = 0
    offsetf(2) = offsetk
    offsetf(4) = params%kpp%inode * dm_wfn%Nb
    call safe_h5sselect_hyperslab(filespace, H5S_SELECT_SET_F, offsetf, count)

    if (mf%kp%nspin==2) then
      ! FHJ: if nspin==1 (including nspinor==2), we simply point wfn_d2 to wfn_d
      SAFE_ALLOCATE(wfn_d2, (SCALARSIZE, dm_wfn%Ml, mf%kp%nspin*mf%kp%nspinor, dm_wfn%Nl))
      do is = 1, mf%kp%nspin
        wfn_d2(1,:,is,:) = dble(wfn_d(:,:,is))
#ifdef CPLX
        wfn_d2(2,:,is,:) = IMAG(wfn_d(:,:,is))
#endif
      enddo
    else
      count_ = int(count, kind=kind(ik))
      call wfn2real_ptr(wfn_d, wfn_d2, count_)
    endif

#ifdef MPI
    call safe_h5pcreate(H5P_DATASET_XFER_F, plist_id)
    if (params%wfn_io_mpiio_mode==0) then
      call safe_h5pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE_F)
    elseif (params%wfn_io_mpiio_mode==1) then
      call safe_h5pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT_F)
    endif
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wfn_d2(:,:,:,:), &
      count, error, memspace, filespace, xfer_prp=plist_id)
    call safe_h5pclose(plist_id)
#else
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wfn_d2(:,:,:,:), &
      count, error, memspace, filespace)
#endif
    if (error/=0) call die('Error writing wavefunctions to HDF5 file')
    call safe_h5sclose(memspace)
    call safe_h5sclose(filespace)
    call safe_h5dclose(dset_id)
    call hdf5_close_file(file_id)
    if (params%kpp%inode==0) then
      call system_clock(count=clock_count)
      t1 = clock_count * clock_inv_rate
      write(params%kpp%iunit,'(1x,a)') 'Done writing wavefunctions.'
      write(params%kpp%iunit,'(1x,a,f0.3,a/)') 'Time elapsed: ', t1-t0, ' s.'
      FLUSH(params%kpp%iunit)
    endif
  endif !should_work

  if (iter_data%is_working_leader) then
    write(params%kpp%iunit,'(1x,3a)') &
      'Writing energies to file ', trim(params%fname_wfn_out), '.'
    FLUSH(params%kpp%iunit)
    t0 = t1
  endif

  ! Now, update energies and occupations. Only the kpp%inode==0 of each pool participates.
  if (iter_data%is_working_leader) then
    call hdf5_open_file(trim(params%fname_wfn_out), 'rw', file_id, &
      parallel_io=.true., comm=iter_data%leader_comm)


    ! Write energies (el)
    !--------------------
    call safe_h5dopen(file_id, 'mf_header/kpoints/el', dset_id)
    call safe_h5dget_space(dset_id, filespace)
    count(1) = dm_wfn%N
    count(2) = 1
    count(3) = mf%kp%nspin
    call safe_h5screate_simple(3, count(1:3), memspace)
    offsetf(1) = 0
    offsetf(2) = ik-1
    offsetf(3) = 0
    call safe_h5sselect_hyperslab(filespace, H5S_SELECT_SET_F, offsetf, count(1:3))

#ifdef MPI
    call safe_h5pcreate(H5P_DATASET_XFER_F, plist_id)
    if (params%wfn_io_mpiio_mode==0) then
      call safe_h5pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE_F)
    elseif (params%wfn_io_mpiio_mode==1) then
      call safe_h5pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT_F)
    endif
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, en(:,:), &
      count(1:3), error, memspace, filespace, xfer_prp=plist_id)
    call safe_h5pclose(plist_id)
#else
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, en(:,:), &
      count(1:3), error, memspace, filespace)
#endif
    if (error/=0) call die('Error writing energies to HDF5 file')
    call safe_h5sclose(memspace)
    call safe_h5sclose(filespace)
    call safe_h5dclose(dset_id)


    ! Write occupations (occ)
    !------------------------
    call safe_h5dopen(file_id, 'mf_header/kpoints/occ', dset_id)
    call safe_h5dget_space(dset_id, filespace)
    count(1) = dm_wfn%N
    count(2) = 1
    count(3) = mf%kp%nspin
    call safe_h5screate_simple(3, count(1:3), memspace)
    offsetf(1) = 0
    offsetf(2) = ik-1
    offsetf(3) = 0
    call safe_h5sselect_hyperslab(filespace, H5S_SELECT_SET_F, offsetf, count(1:3))

#ifdef MPI
    call safe_h5pcreate(H5P_DATASET_XFER_F, plist_id)
    if (params%wfn_io_mpiio_mode==0) then
      call safe_h5pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE_F)
    elseif (params%wfn_io_mpiio_mode==1) then
      call safe_h5pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT_F)
    endif
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mf%kp%occ(:,ik:ik,:), &
      count(1:3), error, memspace, filespace, xfer_prp=plist_id)
    call safe_h5pclose(plist_id)
#else
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mf%kp%occ(:,ik:ik,:), &
      count(1:3), error, memspace, filespace)
#endif
    if (error/=0) call die('Error writing occupations to HDF5 file')
    call safe_h5sclose(memspace)
    call safe_h5sclose(filespace)
    call safe_h5dclose(dset_id)


    ! Write occupations (ifmax)
    !--------------------------
    call safe_h5dopen(file_id, 'mf_header/kpoints/ifmax', dset_id)
    call safe_h5dget_space(dset_id, filespace)
    count(1) = 1
    count(2) = mf%kp%nspin
    call safe_h5screate_simple(2, count(1:2), memspace)
    offsetf(1) = ik-1
    offsetf(2) = 0
    call safe_h5sselect_hyperslab(filespace, H5S_SELECT_SET_F, offsetf(1:2), count(1:2))

#ifdef MPI
    call safe_h5pcreate(H5P_DATASET_XFER_F, plist_id)
    if (params%wfn_io_mpiio_mode==0) then
      call safe_h5pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE_F)
    elseif (params%wfn_io_mpiio_mode==1) then
      call safe_h5pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT_F)
    endif
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, mf%kp%ifmax(ik:ik,:), &
      count(1:2), error, memspace, filespace, xfer_prp=plist_id)
    call safe_h5pclose(plist_id)
#else
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, mf%kp%ifmax(ik:ik,:), &
      count(1:2), error, memspace, filespace)
#endif
    if (error/=0) call die('Error writing max. occupied states to HDF5 file')
    call safe_h5sclose(memspace)
    call safe_h5sclose(filespace)
    call safe_h5dclose(dset_id)

    call hdf5_close_file(file_id)

    call system_clock(count=clock_count)
    t1 = clock_count * clock_inv_rate
    write(params%kpp%iunit,'(1x,a)') 'Done writing energies and occupations.'
    write(params%kpp%iunit,'(1x,a,f0.3,a/)') 'Time elapsed: ', t1-t0, ' s.'
    FLUSH(params%kpp%iunit)
  endif

  if (mf%kp%nspin==2) then
    SAFE_DEALLOCATE_P(wfn_d2)
  endif

contains

  subroutine wfn2real_ptr(wfn_d, wfn_d2, count)
    SCALAR, target, intent(in) :: wfn_d(:,:,:)
    real(DP), pointer, intent(out) :: wfn_d2(:,:,:,:)
    integer, intent(in) :: count(:)

    type(c_ptr) :: ptr

    ptr = c_loc(wfn_d)
    call c_f_pointer(ptr, wfn_d2, count)

  end subroutine wfn2real_ptr

end subroutine wfn_write_kpt_direct
#endif

end module wfn_io_m
