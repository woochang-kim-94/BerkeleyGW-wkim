!=================================================================================
!
! Module write_matrix_m
!
! (1) write_matrix_d()          Originally by JRD       Last Modified 5/1/2008 (JRD)
!
! This program writes a distributed matrix like chimat or epsmat to file.
!
! (2) write_matrix_f()          Originally by JRD       Last Modified 2/5/2009 (CHP)
!
! Modification of write_matrix_d for full-frequency.
!
!=================================================================================

#include "f_defs.h"

module write_matrix_m

  use, intrinsic :: iso_c_binding
  use global_m
#ifdef HDF5
  use hdf5
#endif
  use hdf5_io_m
  use scalapack_m
  use io_utils_m
  use timing_m, only: timing => common_timing

  implicit none

  private

  public :: &
    write_matrix_d, &
    write_matrix_d_sub, &
    write_matrix_f
#ifdef HDF5
  public :: &
    write_matrix_ser_hdf, &
    write_matrix_f_ser_hdf, &
    write_matrix_d_hdf, &
    write_matrix_f_hdf, &
    write_matrix_diagonal_hdf, &
    write_gvec_indices_hdf, &
    write_vcoul_hdf
#ifdef USESCALAPACK
  public :: &
    write_matrix_d_par_hdf, &
    write_matrix_d_par_hdf_sub, &
    write_matrix_f_par_hdf
#endif
#endif

contains

!===================================================================================

subroutine write_matrix_d(scal,matrix,nmtx,iunit)
  type(scalapack), intent(in) :: scal
  SCALAR, intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iunit

  integer :: ii, jj
#ifdef USESCALAPACK
  SCALAR, allocatable :: tempcol(:),tempcol2(:)
  integer :: irow, icol, irowm, icolm
  integer :: icurr
#endif
  type(progress_info) :: prog_info !< a user-friendly progress report

  PUSH_SUB(write_matrix_d)

  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nmtx, iunit
    write(6,*)
  endif

#ifdef USESCALAPACK
  SAFE_ALLOCATE(tempcol, (nmtx))
  SAFE_ALLOCATE(tempcol2, (nmtx))

  icurr=0

  call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  do jj = 1, nmtx
    call progress_step(prog_info, jj)
!        if (peinf%inode .eq. 0) then
!          write(6,*) ' In loop: ', ii
!        endif
    icol=MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
    tempcol=0d0
    if (icol .eq. scal%mypcol) then
      do ii = 1, nmtx
        irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
        if (irow .eq. scal%myprow) then
          icurr=icurr+1
          icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
          irowm=MOD((icurr-1),scal%npr)+1
          tempcol(ii)=matrix(irowm,icolm)

!                if (icolm .gt. scal%npc .or. irowm.gt.scal%npr) then
!                  write(6,*) 'Error: ', scal%npc,scal%npr,icolm,irowm
!                endif

        endif
      enddo
    endif
    if (peinf%inode .eq. 0) then
      tempcol2=0d0
    endif
    call MPI_REDUCE(tempcol,tempcol2,nmtx,MPI_SCALAR,MPI_SUM,0, &
      MPI_COMM_WORLD,mpierr)
    if (peinf%inode .eq. 0) then
      write(iunit) (tempcol2(ii),ii=1,nmtx)
    endif

    call MPI_barrier(MPI_COMM_WORLD,mpierr)

  enddo
  call progress_free(prog_info)

  SAFE_DEALLOCATE(tempcol)
  SAFE_DEALLOCATE(tempcol2)

!      if(peinf%inode .eq. 0) then
!        write(6,*) ' Done Writing chimat: '
!      endif

#else

  if (peinf%inode .eq. 0) then
    call progress_init(prog_info, 'writing matrix', 'column', nmtx)
    do jj = 1, nmtx
      call progress_step(prog_info, jj)
      write(iunit) (matrix(ii, jj), ii = 1, nmtx)
    enddo
    call progress_free(prog_info)
  endif

#endif

  POP_SUB(write_matrix_d)

  return
end subroutine write_matrix_d

subroutine write_matrix_d_sub(scal,matrix,nmtx,iunit,neig_sub)
  type(scalapack), intent(in) :: scal
  complex(DPC), intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iunit
  integer, intent(in), optional :: neig_sub

  integer :: ii, jj, nmtx_col
#ifdef USESCALAPACK
  complex(DPC), allocatable :: tempcol(:),tempcol2(:)
  integer :: irow, icol, irowm, icolm
  integer :: icurr
#endif
  type(progress_info) :: prog_info !< a user-friendly progress report

  PUSH_SUB(write_matrix_d_sub)

  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nmtx, iunit
    write(6,*)
  endif

  ! neig_sub allows to write only neig_sub columns of the actual matrix
  nmtx_col = nmtx
  IF(PRESENT(neig_sub)) nmtx_col = neig_sub

#ifdef USESCALAPACK
  SAFE_ALLOCATE(tempcol, (nmtx))
  SAFE_ALLOCATE(tempcol2, (nmtx))

  icurr=0

  call progress_init(prog_info, 'writing matrix', 'column', nmtx_col)
  do jj = 1, nmtx_col
    call progress_step(prog_info, jj)
!        if (peinf%inode .eq. 0) then
!          write(6,*) ' In loop: ', ii
!        endif
    icol=MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
    tempcol=0d0
    if (icol .eq. scal%mypcol) then
      do ii = 1, nmtx
        irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
        if (irow .eq. scal%myprow) then
          icurr=icurr+1
          icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
          irowm=MOD((icurr-1),scal%npr)+1
          tempcol(ii)=matrix(irowm,icolm)

!                if (icolm .gt. scal%npc .or. irowm.gt.scal%npr) then
!                  write(6,*) 'Error: ', scal%npc,scal%npr,icolm,irowm
!                endif

        endif
      enddo
    endif
    if (peinf%inode .eq. 0) then
      tempcol2=0d0
    endif
    call MPI_REDUCE(tempcol,tempcol2,nmtx,MPI_COMPLEX_DPC,MPI_SUM,0, &
      MPI_COMM_WORLD,mpierr)
    if (peinf%inode .eq. 0) then
      write(iunit) (tempcol2(ii),ii=1,nmtx)
    endif

    call MPI_barrier(MPI_COMM_WORLD,mpierr)

  enddo
  call progress_free(prog_info)

  if (peinf%inode .eq. 0) then
    ! write empty rows for the missin column so sigma will work also with previously
    ! generated epsmat (this can go in the future)
    do jj = nmtx_col + 1, nmtx
      write(iunit)
    end do
  end if

  SAFE_DEALLOCATE(tempcol)
  SAFE_DEALLOCATE(tempcol2)

!      if(peinf%inode .eq. 0) then
!        write(6,*) ' Done Writing chimat: '
!      endif

#else

  if (peinf%inode .eq. 0) then
    call progress_init(prog_info, 'writing matrix', 'column', nmtx_col)
    do jj = 1, nmtx_col
      call progress_step(prog_info, jj)
      write(iunit) (matrix(ii, jj), ii = 1, nmtx)
    enddo
    call progress_free(prog_info)
    !XXXX
    do jj = nmtx_col + 1, nmtx
      write(iunit)
    end do
    !XXXX
  endif

#endif

  POP_SUB(write_matrix_d_sub)

  return
end subroutine write_matrix_d_sub

!=================================================================================

subroutine write_matrix_f(scal,nfreq,retarded,nmtx,iunit,nfreq_group,advanced)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq
  complex(DPC), intent(in) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iunit
  integer, intent(in) :: nfreq_group
  complex(DPC), optional, intent(in) :: advanced(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)

  integer :: ii, jj, ifreq
#ifdef USESCALAPACK
  complex(DPC), allocatable :: tempcolR(:,:),tempcolR2(:,:)
  complex(DPC), allocatable :: tempcolA(:,:),tempcolA2(:,:)
  integer :: irow, icol, irowm, icolm, freq_grp_ind, ifreq_para
#endif
  type(progress_info) :: prog_info !< a user-friendly progress report
  logical :: has_advanced

  PUSH_SUB(write_matrix_f)

  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nfreq, nmtx, iunit
    write(6,*)
  endif

  has_advanced = present(advanced)
#ifdef USESCALAPACK
  SAFE_ALLOCATE(tempcolR, (nfreq,nmtx))
  SAFE_ALLOCATE(tempcolR2, (nfreq,nmtx))
#ifdef CPLX
  SAFE_ALLOCATE(tempcolA, (nfreq,nmtx))
  SAFE_ALLOCATE(tempcolA2, (nfreq,nmtx))
#endif

  call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  do jj = 1, nmtx
    call progress_step(prog_info, jj)
    icol = indxg2p(jj, scal%nbl, scal%mypcol, 0, scal%npcol)
    tempcolR(:,:) = (0d0, 0d0)
#ifdef CPLX
    tempcolA(:,:) = (0d0, 0d0)
#endif
    if (icol==scal%mypcol) then
      icolm = indxg2l(jj, scal%nbl, scal%mypcol, 0, scal%npcol)
      do irowm = 1, scal%npr
        ii = indxl2g(irowm, scal%nbl, scal%myprow, 0, scal%nprow)
        do ifreq = 1, nfreq
          freq_grp_ind = mod(ifreq-1,nfreq_group)
          ifreq_para = (ifreq+nfreq_group-1)/nfreq_group
          if (freq_grp_ind .eq. peinf%rank_mtxel) then
            tempcolR(ifreq,ii) = retarded(irowm,icolm,ifreq_para)
#ifdef CPLX
            if (has_advanced) then
              tempcolA(ifreq,ii) = advanced(irowm,icolm,ifreq_para)
            endif
#endif
          endif
        enddo
      enddo
    endif
    if (peinf%inode .eq. 0) then
      tempcolR2(:,:) = (0d0, 0d0)
#ifdef CPLX
      if (has_advanced) then
        tempcolA2(:,:) = (0d0, 0d0)
      endif
#endif
    endif
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(tempcolR(1,1),tempcolR2(1,1),nfreq*nmtx, &
      MPI_COMPLEX_DPC,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
#ifdef CPLX
    if (has_advanced) then
      call MPI_Reduce(tempcolA(1,1),tempcolA2(1,1),nfreq*nmtx, &
        MPI_COMPLEX_DPC,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
    endif
#endif
    if (peinf%inode .eq. 0) then
      do ii = 1, nmtx
        write(iunit) (tempcolR2(ifreq,ii),ifreq=1,nfreq)
      enddo
#ifdef CPLX
      if (has_advanced) then
        do ii = 1, nmtx
          write(iunit) (tempcolA2(ifreq,ii),ifreq=1,nfreq)
        enddo
      else
        do ii = 1, nmtx
          write(iunit)
        enddo
      endif
#endif
    endif
    call MPI_Barrier(MPI_COMM_WORLD,mpierr)
  enddo
  call progress_free(prog_info)

  SAFE_DEALLOCATE(tempcolR)
  SAFE_DEALLOCATE(tempcolR2)
#ifdef CPLX
  SAFE_DEALLOCATE(tempcolA)
  SAFE_DEALLOCATE(tempcolA2)
#endif

#else

  if(peinf%inode .eq. 0) then
    call progress_init(prog_info, 'writing matrix', 'column', nmtx)
    do jj = 1, nmtx
      call progress_step(prog_info, jj)
      do ii = 1, nmtx
        write(iunit) (retarded(ii, jj, ifreq), ifreq= 1, nfreq)
      enddo
#ifdef CPLX
      if (has_advanced) then
        do ii = 1, nmtx
          write(iunit) (advanced(ii, jj, ifreq),ifreq = 1, nfreq)
        enddo
      else
        do ii = 1, nmtx
          write(iunit)
        enddo
      endif
#endif
    enddo
    call progress_free(prog_info)
  endif

#endif

  POP_SUB(write_matrix_f)

  return
end subroutine write_matrix_f

#ifdef HDF5

!========================================================================
! JRD: The HDF5 Equivalents of the above routines.
!========================================================================

!==========================================================================================

subroutine write_matrix_diagonal_hdf(epsdiag,nmtx,iq,isize,name)
  real(DP), intent(in) :: epsdiag(:,:,:) !< (isize,nmtx,1)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  integer, intent(in) :: isize
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem

  integer(HSIZE_T) :: dims(3), offset(3), offsetm(3)

  integer :: error, rank

  PUSH_SUB(write_matrix_diagonal_hdf)

  call hdf5_open_file(trim(name), 'rw', file_id)

! Write Array

  rank = 3
  dims(1) = isize
  dims(2) = nmtx
  dims(3) = 1
  offset(1) = 0
  offset(2) = 0
  offset(3) = iq - 1
  offsetm(:) = 0

  call h5dopen_f(file_id, 'mats/matrix-diagonal', dset_id, error)
  call h5screate_simple_f(rank, dims, memspace, error)
  call h5dget_space_f(dset_id,filespace,error)
  call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, dims, error)
  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, epsdiag, dims, error, memspace, filespace)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(memspace, error)
  call h5sclose_f(filespace, error)
  call hdf5_close_file(file_id)

  POP_SUB(write_matrix_diagonal_hdf)

end subroutine write_matrix_diagonal_hdf

!===================================================================================

subroutine write_matrix_ser_hdf(matrix,nmtx,iq,is,name)
  SCALAR, intent(in) :: matrix(:,:) !< (nmtx,nmtx)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name

  integer :: error, rank, ii, jj

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem

  integer(HSIZE_T) :: count(6), offset(6), offsetm(6)

  real(DP), allocatable :: data(:,:,:,:,:,:)

  PUSH_SUB(write_matrix_ser_hdf)

  ! FHJ - FIXME: complex/scalar interface is bad
  call hdf5_open_file(trim(name), 'rw', file_id)

  rank=6
  count(1) = SCALARSIZE
  count(2) = nmtx
  count(3) = nmtx
  count(4) = 1
  count(5) = 1
  count(6) = 1

  offset(:) = 0
  offset(5) = is - 1
  offset(6) = iq - 1

  offsetm(:) = 0

  SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))

  do jj = 1, nmtx
  do ii = 1, nmtx
    data(1,ii,jj,1,1,1) = dble(matrix(ii,jj))
#ifdef CPLX
    data(2,ii,jj,1,1,1) = IMAG(matrix(ii,jj))
#endif
  enddo
  enddo

  call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
  call h5screate_simple_f(rank, count, memspace, error)
  call h5dget_space_f(dset_id,filespace,error)
  call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, count, error)
  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, filespace)

  SAFE_DEALLOCATE(data)

  call h5dclose_f(dset_id, error)
  call h5sclose_f(memspace, error)
  call h5sclose_f(filespace, error)
  call hdf5_close_file(file_id)

  POP_SUB(write_matrix_ser_hdf)

end subroutine write_matrix_ser_hdf

!========================================================================

subroutine write_matrix_f_ser_hdf(nfreq, retarded, nmtx, iq, is, name)
  integer, intent(in) :: nfreq
  complex(DPC), intent(in) :: retarded(:,:,:) !< (nmtx,nmtx,nfreq)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name

  integer :: ii, jj, error, rank
  real(DP), allocatable :: data(:,:,:,:,:,:)
  type(progress_info) :: prog_info !< a user-friendly progress report

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem

  integer(HSIZE_T) :: count(6), offset(6), offsetm(6)

  PUSH_SUB(write_matrix_f_ser_hdf)

! DVF: this routine was built off of write_matrix_f_hdf to do the serial
! writing of an hdf format matrix. This is needed for epsmat_old2hdf5.f90

  rank=6
  count(1) = 2
  count(2) = nmtx
  count(3) = 1
  count(4) = nfreq
  count(5) = 1
  count(6) = 1

  SAFE_ALLOCATE(data, (count(1),count(2),count(3),count(4),count(5),count(6)))

  call hdf5_open_file(trim(name), 'rw', file_id)

  call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
  call h5screate_simple_f(rank, count, memspace, error)
  call h5dget_space_f(dset_id,filespace,error)

  call progress_init(prog_info, 'writing matrix', 'column', nmtx)
! JRD XXX the fact that jj is not outer loop presents a bit of challenge
! but this serial routine is just for legacy support
  do jj = 1, nmtx
    call progress_step(prog_info, jj)
      data(1,:,1,:,1,1)=dble(retarded(:,jj,:))
      data(2,:,1,:,1,1)=IMAG(retarded(:,jj,:))

    offsetm(:) = 0
    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, count, error)

    offset(1)=0
    offset(2)=0
    offset(3)=jj-1
    offset(4)=0
    offset(5)=is-1
    offset(6)=iq-1

    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, filespace)
  enddo
  call progress_free(prog_info)

  SAFE_DEALLOCATE(data)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(memspace, error)
  call h5sclose_f(filespace, error)
  call hdf5_close_file(file_id)

  POP_SUB(write_matrix_f_ser_hdf)

  return

end subroutine write_matrix_f_ser_hdf

!========================================================================

subroutine write_matrix_d_hdf(scal,matrix,nmtx,iq,is,name)
  type(scalapack), intent(in) :: scal
  SCALAR, intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name

  integer :: ii, jj, error, size, rank
#ifdef USESCALAPACK
  real(DP), allocatable :: datatmp(:,:,:,:,:,:)
  integer :: irow, icol, irowm, icolm
  integer :: icurr
#endif
  real(DP), allocatable :: data(:,:,:,:,:,:)
  type(progress_info) :: prog_info !< a user-friendly progress report

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
#ifdef USESCALAPACK
!  integer(HID_T) :: plist_id      ! Property list identifier for parallel IO
!                                 ! Not used yet...
#endif

  integer(HSIZE_T) :: count(6), offset(6), offsetm(6)

  PUSH_SUB(write_matrix_d_hdf)

  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nmtx
    write(6,*)
  endif

! XXX: For now, we will still have only proc 0 write...
! We should changes this to parallel writes. But doing
! this effectively from the scalapack, block cyclic layout
! seems a bit tricky. So, ignoring for now...

  rank=6
  count(1) = SCALARSIZE
  count(2) = nmtx
  count(3) = 1
  count(4) = 1
  count(5) = 1
  count(6) = 1

  if (peinf%inode .eq. 0) then
    SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
    call hdf5_open_file(trim(name), 'rw', file_id)
    call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
    call h5screate_simple_f(rank, count, memspace, error)
    call h5dget_space_f(dset_id,filespace,error)
  endif

#ifdef USESCALAPACK
  SAFE_ALLOCATE(datatmp, (count(1),count(2),count(3),count(4),count(5),count(6)))
  icurr=0
#endif

  call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  do jj = 1, nmtx

    call progress_step(prog_info, jj)
#ifdef USESCALAPACK

    call timing%start(timing%eps_i_o_comm)

    icol=MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
    datatmp=0d0
    if (icol .eq. scal%mypcol) then
      do ii = 1, nmtx
        irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
        if (irow .eq. scal%myprow) then
          icurr=icurr+1
          icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
          irowm=MOD((icurr-1),scal%npr)+1
          datatmp(1,ii,1,1,1,1)=dble(matrix(irowm,icolm))
#ifdef CPLX
          datatmp(2,ii,1,1,1,1)=IMAG(matrix(irowm,icolm))
#endif
        endif
      enddo
    endif
    if (peinf%inode .eq. 0) then
      data=0d0
    endif

! XXX This is a big waste of communication. Should be fixed when do
! parallel IO.

    size = nmtx * SCALARSIZE

    call MPI_REDUCE(datatmp,data,size,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      MPI_COMM_WORLD,mpierr)

    call timing%stop(timing%eps_i_o_comm)

#else

    if (peinf%inode .eq. 0) then
      do ii = 1, nmtx
        data(1,ii,1,1,1,1) = dble(matrix(ii,jj))
#ifdef CPLX
        data(2,ii,1,1,1,1) = IMAG(matrix(ii,jj))
#endif
      enddo
    endif

#endif

    call timing%start(timing%eps_i_o_io)

    if (peinf%inode .eq. 0) then

      offsetm(:) = 0
      call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, count, error)

      offset(1)=0
      offset(2)=0
      offset(3)=jj-1
      offset(4)=0
      offset(5)=is-1
      offset(6)=iq-1

      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, filespace)

    endif

    call timing%stop(timing%eps_i_o_io)

#ifdef USESCALAPACK
    call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif

  enddo
  call progress_free(prog_info)

#ifdef USESCALAPACK
  SAFE_DEALLOCATE(datatmp)
#endif

  if (peinf%inode .eq. 0) then
    SAFE_DEALLOCATE(data)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(memspace, error)
    call h5sclose_f(filespace, error)
    call hdf5_close_file(file_id)
  endif

  POP_SUB(write_matrix_d_hdf)

  return
end subroutine write_matrix_d_hdf

!========================================================================

#ifdef USESCALAPACK

subroutine write_matrix_d_par_hdf(scal,matrix,nmtx,iq,is,name)
  type(scalapack), intent(in) :: scal
  SCALAR, intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name

  integer :: ii, jj, error, rank
  real(DP), allocatable :: data(:,:,:,:,:,:)

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
  integer(HID_T) :: plist_id      ! Property list identifier for parallel IO

  integer(HSIZE_T) :: count(6), countm(6), offset(6), offsetm(6), stride(6), block(6)
  integer(HSIZE_T) :: countr(6), offsetr(6), strider(6), blockr(6)

  integer :: info, rowremainder, colremainder

  PUSH_SUB(write_matrix_d_par_hdf)

  call timing%start(timing%eps_i_o_comm)

! JRD: We need a barrier here or else parallel file opening gets mixed up with
! peinf%inode 0 opening the file to write the diagonal (which is called first).
  call MPI_barrier(MPI_COMM_WORLD,mpierr)

  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nmtx
    write(6,*)
  endif

! JRD Should be ok with npr and npc = 0

  !if (scal%npr .eq. 0 .or. scal%npc .eq. 0) then
  !  write(6,*) peinf%inode,"Zero npr or npc!!", scal%npr, scal%npc
  !endif

  rank=6
  countm(1) = SCALARSIZE
  countm(2) = scal%npr
  countm(3) = scal%npc
  countm(4) = 1
  countm(5) = 1
  countm(6) = 1

  offsetm(:) = 0

  count(1) = 1
  count(2) = scal%npr/scal%nbl
  count(3) = scal%npc/scal%nbl
  count(4) = 1
  count(5) = 1
  count(6) = 1

  block(1) = SCALARSIZE
  block(2) = scal%nbl
  block(3) = scal%nbl
  block(4) = 1
  block(5) = 1
  block(6) = 1

  offset(1) = 0
  offset(2) = scal%myprow*scal%nbl
  offset(3) = scal%mypcol*scal%nbl
  offset(4) = 0
  offset(5) = is-1
  offset(6) = iq-1

  stride(1) = 1
  stride(2) = scal%nprow*scal%nbl
  stride(3) = scal%npcol*scal%nbl
  stride(4) = 1
  stride(5) = 1
  stride(6) = 1

  call hdf5_open_file(trim(name), 'rw', file_id, parallel_io=.true.)

  SAFE_ALLOCATE(data,(countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))
!XXX create data can we avoid duplication?
  do jj = 1, scal%npc
  !do jj = 1, scal%npc - mod(scal%npc,scal%nbl)
    do ii = 1, scal%npr
    !do ii = 1, scal%npr - mod(scal%npr,scal%nbl)
        data(1,ii,jj,1,1,1) = dble(matrix(ii,jj))
#ifdef CPLX
        data(2,ii,jj,1,1,1) = IMAG(matrix(ii,jj))
#endif
    enddo
  enddo

  call h5screate_simple_f(rank, countm, memspace, error)
  call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)

  call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
  call h5dget_space_f(dset_id,filespace,error)

  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)

! Add in remainders, in case scal%nbl doesnt perfectly divide nmtx

! Bottom Rows
  rowremainder = mod(scal%npr,scal%nbl)
  if (rowremainder .ne. 0) then
    offsetr=offset
    countr=count
    blockr=block
    strider=stride
    offsetr(2)=nmtx-rowremainder
    countr(2)=rowremainder
    blockr(2)=1
    strider(2)=1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
    !write(6,*) peinf%inode, "I have the bottom row", rowremainder, scal%npc
  endif

! Right Columns
  colremainder = mod(scal%npc,scal%nbl)
  if (colremainder .ne. 0) then
    offsetr=offset
    countr=count
    blockr=block
    strider=stride
    offsetr(3)=nmtx-colremainder
    countr(3)=colremainder
    blockr(3)=1
    strider(3)=1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
    !write(6,*) peinf%inode, "I have the right column", colremainder, scal%npr
! Bottom Corner of Matrix
    if (rowremainder .ne. 0) then
      offsetr=offset
      countr=count
      blockr=block
      strider=stride
      offsetr(2)=nmtx-rowremainder
      countr(2)=rowremainder
      blockr(2)=1
      strider(2)=1
      offsetr(3)=nmtx-colremainder
      countr(3)=colremainder
      blockr(3)=1
      strider(3)=1
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
      !write(6,*) peinf%inode, "I have bottom both"
    endif
  endif

  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  call timing%stop(timing%eps_i_o_comm)
  call timing%start(timing%eps_i_o_io)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, &
                      xfer_prp = plist_id)
  call timing%stop(timing%eps_i_o_io)
  call h5pclose_f(plist_id, error)

  SAFE_DEALLOCATE(data)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(memspace, error)
  call h5sclose_f(filespace, error)
  call hdf5_close_file(file_id)

  POP_SUB(write_matrix_d_par_hdf)

  return
end subroutine write_matrix_d_par_hdf

!========================================================================================
! subspace routine for eigenvectors and full epsilon zero
subroutine write_matrix_d_par_hdf_sub(scal, matrix, nmtx, nmtx_col, iq, is, name, set_name)
  type(scalapack), intent(in) :: scal
  complex(DPC), intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx, nmtx_col
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: set_name

  integer :: ii, jj, error, rank
  real(DP), allocatable :: data(:,:,:,:,:,:)

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
  integer(HID_T) :: plist_id      ! Property list identifier for parallel IO

  integer(HSIZE_T) :: count(6), countm(6), offset(6), offsetm(6), stride(6), block(6)
  integer(HSIZE_T) :: countr(6), offsetr(6), strider(6), blockr(6)

  integer :: rowremainder, colremainder

  PUSH_SUB(write_matrix_d_par_hdf_sub)

  call timing%start(timing%eps_i_o_comm)

! JRD: We need a barrier here or else parallel file opening gets mixed up with
! peinf%inode 0 opening the file to write the diagonal (which is called first).
  call MPI_barrier(MPI_COMM_WORLD,mpierr)

  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nmtx, 'x', nmtx_col
    write(6,*)
  endif

! JRD Should be ok with npr and npc = 0

  !if (scal%npr .eq. 0 .or. scal%npc .eq. 0) then
  !  write(6,*) peinf%inode,"Zero npr or npc!!", scal%npr, scal%npc
  !endif

  rank=6
  countm(1) = 2
  countm(2) = scal%npr
  countm(3) = scal%npc
  countm(4) = 1
  countm(5) = 1
  countm(6) = 1

  offsetm(:) = 0

  count(1) = 1
  count(2) = scal%npr/scal%nbl
  count(3) = scal%npc/scal%nbl
  count(4) = 1
  count(5) = 1
  count(6) = 1

  block(1) = 2
  block(2) = scal%nbl
  block(3) = scal%nbl
  block(4) = 1
  block(5) = 1
  block(6) = 1

  offset(1) = 0
  offset(2) = scal%myprow*scal%nbl
  offset(3) = scal%mypcol*scal%nbl
  offset(4) = 0
  offset(5) = is-1
  offset(6) = iq-1

  stride(1) = 1
  stride(2) = scal%nprow*scal%nbl
  stride(3) = scal%npcol*scal%nbl
  stride(4) = 1
  stride(5) = 1
  stride(6) = 1

  call hdf5_open_file(trim(name), 'rw', file_id, parallel_io=.true.)

  SAFE_ALLOCATE(data,(countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))
!XXX create data can we avoid duplication?
  do jj = 1, scal%npc
  !do jj = 1, scal%npc - mod(scal%npc,scal%nbl)
    do ii = 1, scal%npr
    !do ii = 1, scal%npr - mod(scal%npr,scal%nbl)
      data(1,ii,jj,1,1,1) = dble(matrix(ii,jj))
      data(2,ii,jj,1,1,1) = IMAG(matrix(ii,jj))
    enddo
  enddo

  call h5screate_simple_f(rank, countm, memspace, error)
  call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)

  call h5dopen_f(file_id, trim(set_name), dset_id, error)
  call h5dget_space_f(dset_id,filespace,error)

  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)

! Add in remainders, in case scal%nbl doesnt perfectly divide nmtx

! Bottom Rows
  rowremainder = mod(scal%npr,scal%nbl)
  if (rowremainder .ne. 0) then
    offsetr=offset
    countr=count
    blockr=block
    strider=stride
    offsetr(2)=nmtx-rowremainder
    countr(2)=rowremainder
    blockr(2)=1
    strider(2)=1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
    !write(6,*) peinf%inode, "I have the bottom row", rowremainder, scal%npc
  endif

! Right Columns
  colremainder = mod(scal%npc,scal%nbl)
  if (colremainder .ne. 0) then
    offsetr=offset
    countr=count
    blockr=block
    strider=stride
    offsetr(3)=nmtx_col - colremainder
    countr(3)=colremainder
    blockr(3)=1
    strider(3)=1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
    !write(6,*) peinf%inode, "I have the right column", colremainder, scal%npr
! Bottom Corner of Matrix
    if (rowremainder .ne. 0) then
      offsetr=offset
      countr=count
      blockr=block
      strider=stride
      offsetr(2)=nmtx-rowremainder
      countr(2)=rowremainder
      blockr(2)=1
      strider(2)=1
      offsetr(3)=nmtx_col-colremainder
      countr(3)=colremainder
      blockr(3)=1
      strider(3)=1
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
      !write(6,*) peinf%inode, "I have bottom both"
    endif
  endif

  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  call timing%stop(timing%eps_i_o_comm)
  call timing%start(timing%eps_i_o_io)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, &
                      xfer_prp = plist_id)
  call timing%stop(timing%eps_i_o_io)
  call h5pclose_f(plist_id, error)

  SAFE_DEALLOCATE(data)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(memspace, error)
  call h5sclose_f(filespace, error)
  call hdf5_close_file(file_id)

  POP_SUB(write_matrix_d_par_hdf_sub)

  return
end subroutine write_matrix_d_par_hdf_sub

!========================================================================================

!> Write a matrix `retarded` to file `name`: BLACS istributed version.
!! Internally, this routin may use two drivers:
!! - write_matrix_f_par_hdf_direct: uses HDF5 directly.
!! - write_matrix_f_par_hdf_redist: first redistrubutes the matrix so that
!!   each rank performs a single contiguous write to the file.
subroutine write_matrix_f_par_hdf(scal, nfreq_in_group, retarded, nmtx, &
                                  iq, is, name, nfreq_group, dset_name)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq_in_group
  complex(DPC), intent(in) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name
  integer, intent(in) :: nfreq_group
  character(len=*), intent(in), optional :: dset_name

  character(len=128) :: dset_name_

  PUSH_SUB(write_matrix_f_par_hdf)

  dset_name_ = 'mats/matrix'
  if (present(dset_name)) dset_name_ = trim(dset_name)

  if (peinf%hdf5_write_redist) then
    if (peinf%inode==0) write(6,'(1x,a)') &
      'Note: using new code to first redistribute matrix into contiguous chunks.'
    call write_matrix_f_par_hdf_redist(scal, nfreq_in_group, retarded, nmtx, &
      iq, is, name, nfreq_group, trim(dset_name_))
  else
    if (peinf%inode==0) write(6,'(1x,a)') &
      'Note: using old code to directly write the matrix to file.'
    call write_matrix_f_par_hdf_direct(scal, nfreq_in_group, retarded, nmtx, &
      iq, is, name, nfreq_group, trim(dset_name_))
  endif


  POP_SUB(write_matrix_f_par_hdf)

end subroutine write_matrix_f_par_hdf

!----------------------------------------------------------------------------------------

!> Writes a matrix `retarded` to file `name` by first redistributing the data.
!! This makes sure that each rank performs a single contiguous write to the file.
subroutine write_matrix_f_par_hdf_redist(scal, nfreq_in_group, retarded, nmtx, &
                                         iq, is, name, nfreq_group, dset_name)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq_in_group
  complex(DPC), intent(in) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name
  integer, intent(in) :: nfreq_group
  character(len=*), intent(in) :: dset_name

  type(progress_info) :: prog_info !< a user-friendly progress report
  integer :: nmtx_max, ncols_glob, ncols_loc, nb, iw, ipe, icol0_g
  integer :: icntxt_2d, icntxt_1d, desc_2d(9), desc_1d(9)
  integer :: usermap(1,peinf%npes)
  complex(DPC), allocatable, target :: mat_loc(:,:) !< (nmtx_max, ncols_loc)
  real(DPC), pointer :: mat_loc_dp(:,:,:)

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
  integer(HID_T) :: plist_id      ! Property list identifier for parallel IO
  integer :: error, info
  integer(HSIZE_T) :: count_m(6), offset_f(6)
  logical :: do_subspace, is_subspace_basis  ! check if matrices are in subspace basis

  ! FHJ: old idea:
  ! FIXME - we assume that nfreq_in_group==nfreq, ie., nfreq_group==1
  ! The global distributed matrix has dimensions (nmtx,nmtx,nfreq), and
  ! each processor owns a submatrix of size (scal%npr,scal%npc,nfreq_in_group)
  ! We will redistribute the global matrix so that each processor owns a
  ! contiguous chunk of the global matrix:
  !  matrix_loc => (nmtx_max,ncols_loc)
  ! where ncols_loc = nmtx * nfreq / npes
  !
  ! For now:
  ! Loop over frequencies, for each frequency, redistribute matrix `retarded`
  ! to `mat_loc`.

  PUSH_SUB(write_matrix_f_par_hdf_redist)

  ! Read nmtx_max from file. Sure, we could have asked this as a parameter to
  ! this subroutine..
  do_subspace = .false.
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  if (peinf%inode==0) then
    call hdf5_open_file(trim(name), 'r', file_id)
    call hdf5_read_int(file_id, 'eps_header/gspace/nmtx_max', nmtx_max)
    ! check if matrix is in subspace basis 
    call hdf5_read_logical(file_id, 'eps_header/params/subspace', do_subspace)
    if (do_subspace) then
      is_subspace_basis = .false.
      call hdf5_read_logical(file_id, 'eps_header/subspace/matrix_in_subspace_basis', is_subspace_basis)
      if (is_subspace_basis) then
        call hdf5_read_int(file_id, 'eps_header/subspace/neig_max', nmtx_max)
      end if
    end if
    ! 
    call hdf5_close_file(file_id)
  endif
  call MPI_Bcast(nmtx_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

  ! Set up size of the 1D contiguously distributed matrix
  ncols_glob = nmtx !total number of columns in distrib. matrix
  nb = DIVUP(ncols_glob, peinf%npes) !block size
  ! Number of columns I own:
  ncols_loc = NUMROC(ncols_glob, nb, peinf%inode, 0, peinf%npes)
  SAFE_ALLOCATE(mat_loc, (nmtx_max,ncols_loc)) !local matrix
  mat_loc(:,:) = ZERO
  ! first global column I own:
  icol0_g = indxl2g(1, nb, peinf%inode, 0, peinf%npes)
  ! Create a double precision pointer to the complex data, without
  ! copying any data. This is necessary because the input data is
  ! complex, but we write to a double precision dataset.
  call c_f_pointer(c_loc(mat_loc), mat_loc_dp, [2, nmtx_max, ncols_loc])

  ! FHJ: Create BLACS context for 1d cyclic matrices. Processors: all aboard!
  call blacs_get(-1, 0, icntxt_1d)
  do ipe = 1, peinf%npes
    usermap(1,ipe) = ipe - 1
  enddo
  call blacs_gridmap(icntxt_1d, usermap, 1, 1, peinf%npes) ! This will modify cntxt_1d
  call descinit(desc_1d, nmtx, nmtx, nmtx, nb, 0, 0, &
    icntxt_1d, nmtx_max, info)
  if (info/=0) call die('Error in descinit to create 1D descriptor')

  call descinit(desc_2d, nmtx, nmtx, scal%nbl, scal%nbl, 0, 0, &
    scal%icntxt, max(1,scal%npr), info)
  if (info/=0) call die('Error in descinit to create 2D descriptor')

  call hdf5_open_file(trim(name), 'rw', file_id, parallel_io=.true.)

  count_m = [2, nmtx_max, ncols_loc, 1, 1, 1]
  call h5screate_simple_f(6, count_m, memspace, error)
  call h5sselect_all_f(memspace, error)

  call h5dopen_f(file_id, trim(dset_name), dset_id, error)
  call h5dget_space_f(dset_id, filespace, error)

  ! Offset: (flavor, G, Gp, freq, matrix_index, q-point)
  offset_f(:) = [0, 0, icol0_g-1, 0, is-1, iq-1]
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
  !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  call progress_init(prog_info, 'writing matrix', 'frequency', nfreq_in_group)
  do iw = 1, nfreq_in_group
    call progress_step(prog_info)

    ! FHJ: Copy matrix `retarded` from 2d cyclic layout to `mat_loc` (1d cyclic layout)
    call timing%start(timing%eps_i_o_comm)
    call pzgemr2d(nmtx, nmtx, retarded(:,:,iw), 1, 1, desc_2d, mat_loc(:,:), &
                  1, 1, desc_1d, icntxt_1d)
    call timing%stop(timing%eps_i_o_comm)

    ! Update frequency offset and write to file
    offset_f(4) = iw - 1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_f, count_m, error)
    call timing%start(timing%eps_i_o_io)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mat_loc_dp, count_m, error, &
                    memspace, filespace, xfer_prp=plist_id)
    call timing%stop(timing%eps_i_o_io)
  enddo
  call progress_free(prog_info)

  SAFE_DEALLOCATE(mat_loc)
  call h5pclose_f(plist_id, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(memspace, error)
  call h5sclose_f(filespace, error)
  call hdf5_close_file(file_id)

  POP_SUB(write_matrix_f_par_hdf_redist)

end subroutine write_matrix_f_par_hdf_redist

!----------------------------------------------------------------------------------------

!> Directly writes a matrix `retarded` to file `name`
subroutine write_matrix_f_par_hdf_direct(scal, nfreq_in_group, retarded, nmtx, &
                                         iq, is, name, nfreq_group, dset_name)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq_in_group
  complex(DPC), intent(in) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name
  integer, intent(in) :: nfreq_group
  character(len=*), intent(in) :: dset_name

  integer :: ii, jj, error, rank
  real(DP), allocatable :: data(:,:,:,:,:,:)

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
  integer(HID_T) :: plist_id      ! Property list identifier for parallel IO

  integer(HSIZE_T) :: count(6), countm(6), offset(6), offsetm(6), stride(6), block(6)
  integer(HSIZE_T) :: countr(6), offsetr(6), strider(6), blockr(6)

  integer :: rowremainder, colremainder

  PUSH_SUB(write_matrix_f_par_hdf_direct)

  call timing%start(timing%eps_i_o_comm)

! JRD: We need a barrier here or else parallel file opening gets mixed up with
! peinf%inode 0 opening the file to write the diagonal (which is called first).
  call MPI_barrier(MPI_COMM_WORLD,mpierr)

  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nmtx
    write(6,*)
  endif

! JRD Should be ok with npr and npc = 0

  rank=6
  countm(1) = 2
  countm(2)=max(1,scal%npr)
  countm(3)=max(1,scal%npc)
  countm(4)=max(1,nfreq_in_group)
  countm(5) = 1
  countm(6) = 1

  offsetm(:) = 0

  count(1) = 1
  count(2) = scal%npr/scal%nbl
  count(3) = scal%npc/scal%nbl
  if(nfreq_group .gt. 1) then
    count(4) = nfreq_in_group
  else
    count(4) = 1
  endif
  count(5) = 1
  count(6) = 1

  block(1) = 2
  block(2) = scal%nbl
  block(3) = scal%nbl
  if(nfreq_group .gt. 1) then
    block(4) = 1
  else
    block(4) = nfreq_in_group
  endif
  block(5) = 1
  block(6) = 1

  offset(1) = 0
  offset(2) = scal%myprow*scal%nbl
  offset(3) = scal%mypcol*scal%nbl
  offset(4) = peinf%igroup_f
  offset(5) = is-1
  offset(6) = iq-1

  stride(1) = 1
  stride(2) = scal%nprow*scal%nbl
  stride(3) = scal%npcol*scal%nbl
  stride(4) = nfreq_group
  stride(5) = 1
  stride(6) = 1

  call hdf5_open_file(trim(name), 'rw', file_id, parallel_io=.true.)

  SAFE_ALLOCATE(data,(countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))
!XXX create data can we avoid duplication?
!XXX THREAD? Yes we should to get better bandwidth
  if (scal%npr/=0 .and. scal%npc/=0 .and. nfreq_in_group/=0) then
    data(1,:,:,:,1,1) = dble(retarded(:,:,:))
    data(2,:,:,:,1,1) = IMAG(retarded(:,:,:))
  endif

  call h5screate_simple_f(rank, countm, memspace, error)

 if(scal%npr*scal%npc .ne. 0) then
    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)
  else
    call H5sselect_none_f(memspace,error)
  endif

  call h5dopen_f(file_id, trim(dset_name), dset_id, error)
  call h5dget_space_f(dset_id,filespace,error)

  if(scal%npr*scal%npc .ne. 0) then
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
  else
    call H5sselect_none_f(filespace,error)
  endif

! Add in remainders, in case scal%nbl doesnt perfectly divide nmtx

! The condition that nfreq_in_group .ne. 0 is here because this
! condition will only be met by the processors that are excluded
! from the calculation when using parallel frequencies with a
! number of processors not divisible by the number of frequency
! groups (for the paranoid: this condition that some processors don`t own any frequencies
! could also be met if the number of frequency groups
! requested by the user is greater than the total number of frequencies calculated,
! but we don`t allow that, i.e. the code dies if the user makes such a request).
! They will have a row/colremainder = 0 because they own no part of the dielectric
! matrix, but we don`t want them to be involved with the hyperslab operations because
! they have no data and are spectators.

! Bottom Rows
  rowremainder = mod(scal%npr,scal%nbl)
  if (rowremainder .ne. 0 .and. nfreq_in_group .ne. 0) then
    offsetr=offset
    countr=count
    blockr=block
    strider=stride
    offsetr(2)=nmtx-rowremainder
    countr(2)=rowremainder
    blockr(2)=1
    strider(2)=1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
    !write(6,*) peinf%inode, "I have the bottom row", rowremainder, scal%npc
  endif

! Right Columns
  colremainder = mod(scal%npc,scal%nbl)
  if (colremainder .ne. 0 .and. nfreq_in_group .ne. 0) then
    offsetr=offset
    countr=count
    blockr=block
    strider=stride
    offsetr(3)=nmtx-colremainder
    countr(3)=colremainder
    blockr(3)=1
    strider(3)=1
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
    !write(6,*) peinf%inode, "I have the right column", colremainder, scal%npr
! Bottom Corner of Matrix
    if (rowremainder .ne. 0) then
      offsetr=offset
      countr=count
      blockr=block
      strider=stride
      offsetr(2)=nmtx-rowremainder
      countr(2)=rowremainder
      blockr(2)=1
      strider(2)=1
      offsetr(3)=nmtx-colremainder
      countr(3)=colremainder
      blockr(3)=1
      strider(3)=1
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
      !write(6,*) peinf%inode, "I have bottom both"
    endif
  endif

  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
  !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  call timing%stop(timing%eps_i_o_comm)
  call timing%start(timing%eps_i_o_io)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, &
                      xfer_prp = plist_id)
  call timing%stop(timing%eps_i_o_io)
  call h5pclose_f(plist_id, error)

  SAFE_DEALLOCATE(data)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(memspace, error)
  call h5sclose_f(filespace, error)
  call hdf5_close_file(file_id)

  POP_SUB(write_matrix_f_par_hdf_direct)

end subroutine write_matrix_f_par_hdf_direct

#endif


!==========================================================================================

subroutine write_matrix_f_hdf(scal, nfreq, retarded, nmtx, iq, is, name)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq
  complex(DPC), intent(in) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iq
  integer, intent(in) :: is
  character(len=*), intent(in) :: name

  integer :: ii, jj, error, size, rank
#ifdef USESCALAPACK
  real(DP), allocatable :: datatmp(:,:,:,:,:,:)
  integer :: irow, icol, irowm, icolm
  integer :: icurr
#endif
  real(DP), allocatable :: data(:,:,:,:,:,:)
  type(progress_info) :: prog_info !< a user-friendly progress report

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem

  integer(HSIZE_T) :: count(6), offset(6), offsetm(6)

  PUSH_SUB(write_matrix_f_hdf)

  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nmtx, nfreq
    write(6,*)
  endif

! XXX: For now, we will still have only proc 0 write...
! We should changes this to parallel writes. But doing
! this effectively from the scalapack, block cyclic layout
! seems a bit tricky. So, ignoring for now...

  rank=6
  count(1) = 2
  count(2) = nmtx
  count(3) = 1
  count(4) = nfreq
  count(5) = 1
  count(6) = 1

  if (peinf%inode .eq. 0) then
    SAFE_ALLOCATE(data, (count(1),count(2),count(3),count(4),count(5),count(6)))
    call hdf5_open_file(trim(name), 'rw', file_id)
    call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
    call h5screate_simple_f(rank, count, memspace, error)
    call h5dget_space_f(dset_id,filespace,error)
  endif

#ifdef USESCALAPACK
  SAFE_ALLOCATE(datatmp, (count(1),count(2),count(3),count(4),count(5),count(6)))
  icurr=0
#endif

  call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  do jj = 1, nmtx
    call progress_step(prog_info, jj)

#ifdef USESCALAPACK

    icol=MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
    datatmp=0d0
! JRD XXX The below is going to be incredibly slow. Hoping all over memory.
    if (icol .eq. scal%mypcol) then
      do ii = 1, nmtx
        irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
        if (irow .eq. scal%myprow) then
          icurr=icurr+1
          icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
          irowm=MOD((icurr-1),scal%npr)+1
          datatmp(1,ii,1,:,1,1)=dble(retarded(irowm,icolm,:))
          datatmp(2,ii,1,:,1,1)=IMAG(retarded(irowm,icolm,:))
        endif
      enddo
    endif
    if (peinf%inode .eq. 0) then
      data=0d0
    endif
! XXX This is a big waste of communication. Should be fixed when do
! parallel IO.

    size = nmtx*nfreq*2

    call MPI_REDUCE(datatmp,data,size,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      MPI_COMM_WORLD,mpierr)

#else
! JRD XXX The below is horrendously slow. Jumping all over memory
    if (peinf%inode .eq. 0) then
      do ii = 1, nmtx
        data(1,ii,1,:,1,1)=dble(retarded(ii,jj,:))
        data(2,ii,1,:,1,1)=IMAG(retarded(ii,jj,:))
      enddo
    endif

#endif

    if (peinf%inode .eq. 0) then

      offsetm(:) = 0
      call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, count, error)

      offset(1)=0
      offset(2)=0
      offset(3)=jj-1
      offset(4)=0
      offset(5)=is-1
      offset(6)=iq-1

      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, filespace)

    endif

#ifdef USESCALAPACK
    call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif


  enddo
  call progress_free(prog_info)

#ifdef USESCALAPACK
  SAFE_DEALLOCATE(datatmp)
#endif

  if (peinf%inode .eq. 0) then
    SAFE_DEALLOCATE(data)
  endif

  if (peinf%inode .eq. 0) then
    call h5dclose_f(dset_id, error)
    call h5sclose_f(memspace, error)
    call h5sclose_f(filespace, error)
    call hdf5_close_file(file_id)
  endif

  POP_SUB(write_matrix_f_hdf)

  return
end subroutine write_matrix_f_hdf

!===================================================================================

subroutine write_gvec_indices_hdf(ng, gind_eps2rho, gind_rho2eps, ekin, iq, name)
  integer, intent(in) :: ng
  integer, intent(in) :: gind_eps2rho(:) !< (ng)
  integer, intent(in) :: gind_rho2eps(:) !< (ng)
  real(DP), intent(in) :: ekin(:) !< (ng)
  integer, intent(in) :: iq
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id
  integer :: error, countf(2), offsetf(2)

  PUSH_SUB(write_gvec_indices_hdf)

  call open_file(99, trim(name), status='old')
  call close_file(99)

  call hdf5_open_file(trim(name), 'rw', file_id)
  countf(:) = (/ng, 1/)
  offsetf(:) = (/0, iq-1/)
  call hdf5_write_int_hyperslab(file_id, 'eps_header/gspace/gind_eps2rho', &
    countf, offsetf, gind_eps2rho, error)
  call hdf5_write_int_hyperslab(file_id, 'eps_header/gspace/gind_rho2eps', &
    countf, offsetf, gind_rho2eps, error)
  call hdf5_write_double_hyperslab(file_id, 'eps_header/gspace/ekin', &
    countf, offsetf, ekin, error)
  call hdf5_close_file(file_id)

  POP_SUB(write_gvec_indices_hdf)

end subroutine write_gvec_indices_hdf

!===================================================================================

subroutine write_vcoul_hdf(vcoul, iq, name)
  real(DP), intent(in) :: vcoul(:) !< (nmtx(iq))
  integer, intent(in) :: iq
  character(len=*), intent(in) :: name

  integer(HID_T) :: file_id
  integer :: countf(2), offsetf(2), nmtx(1)

  PUSH_SUB(write_vcoul_hdf)

  call hdf5_open_file(trim(name), 'rw', file_id)
  call hdf5_read_int_hyperslab(file_id, 'eps_header/gspace/nmtx', (/1/), (/iq-1/), nmtx)
  if (size(vcoul)<nmtx(1)) &
    call die('Internal error: array vcoul too small in write_vcoul_hdf.', only_root_writes=.true.)
  countf(:) = (/nmtx(1), 1/)
  offsetf(:) = (/0, iq-1/)
  call hdf5_write_double_hyperslab(file_id, 'eps_header/gspace/vcoul', &
    countf, offsetf, vcoul)
  call hdf5_close_file(file_id)

  POP_SUB(write_vcoul_hdf)

end subroutine write_vcoul_hdf


#endif

end module write_matrix_m

