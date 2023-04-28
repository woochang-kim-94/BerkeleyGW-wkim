!==============================================================================
!
! Routines:
!
! (1) bsewrite()        Originally By MLT       Last Modified 7/1/2008 (JRD)
!
!     input:  xct types
!            bsedbody,bsedhead,bsedwing,bsex
!     output: binary files "bsedmat", "bsexmat"
!
!     Write out all interaction matrices elements in "bsedmat", "bsexmat"
!
!     Since all PEs write into the same units, the access must be organized:
!     one PE at a time
!
!     It appears the imatrix code is: 1 dhead, 2 dwing, 3 dbody, 4 x
!
! (2) bse_hdf5_write()  Orgiginally By JRD      Last Modified 4/27/2012 (JRD)
!
!     hdf5 writer routine for bsemat.h5 files
!
!=============================================================================

#include "f_defs.h"

module bsewrite_m

#ifdef HDF5
  use hdf5
#endif
  use hdf5_io_m
  use global_m
  use io_utils_m
  use kernel_io_m
  use wfn_rho_vxc_io_m
  implicit none

  public :: bsewrite
#ifdef HDF5
  public :: bse_hdf5_write
#endif

  private

contains


subroutine bsewrite(xct,iownsize,bsedbody,bsedhead,bsex,kg,kp,gvec,syms,crys,&
    bsedwing,bset)
  type (xctinfo), intent(in) :: xct
  type (grid), intent(in) :: kg
  integer, intent(in) :: iownsize
  SCALAR, intent(in) :: bsedbody(:,:,:), bsedhead(:,:,:), &
    bsex(:,:,:) !< (iownsize*peinf%myown,xct%nspin,xct%nspin)
  type (kpoints), intent(in)  :: kp
  type (gspace), intent(in) :: gvec
  type (symmetry), intent(in) :: syms
  type (crystal), intent(in) :: crys
  SCALAR, intent(in), optional :: &
    bsedwing(:,:,:) !< (iownsize*peinf%myown,xct%nspin,xct%nspin)
  SCALAR, intent(in), optional :: &
    bset(:,:,:) !< (iownsize*peinf%myown,xct%nspin,xct%nspin)

  SCALAR, allocatable :: bsemt(:,:,:,:,:), bsemtt(:,:,:,:,:)
  SCALAR :: bsem

  integer :: nmatrices,imatrix,iunit,error,version
  integer :: ic,icp,ik,ikp,is1,is2,iv,ivp,it
  real(DP) :: bsedh_fro2, bsedw_fro2, bsedb_fro2, bsex_fro2, bset_fro2
  real(DP) :: bsedh_fro, bsedw_fro, bsedb_fro, bsex_fro, bset_fro
  real(DP) :: approx_sz
  type(progress_info) :: prog_info
  type(kernel_header_t) :: kern

  PUSH_SUB(bsewrite)

  approx_sz = (xct%n1b_co*xct%n2b_co*xct%nkpt_co*xct%nspin + 8d0)**2 * 8d0
  approx_sz = SCALARSIZE * approx_sz / (1024d0 * 1024d0)

  if (peinf%inode .eq. 0) then
    write(6,'(1x,a)') 'Expected size of the BSE matrices:'
    if (xct%use_hdf5) then
      write(6,'(1x,a,f0.3,a)') '  bsemat.h5: ', approx_sz*4, ' MB'
    else
      write(6,'(1x,a,f0.3,a)') '  bsexmat: ', approx_sz, ' MB'
      write(6,'(1x,a,f0.3,a)') '  bsedmat: ', approx_sz*3, ' MB'
    endif
    write(6,*)
  endif

  ! FHJ: Min. number of matrices we save is 3 => X + W_body + W_head
  nmatrices = 3
  if (present(bsedwing)) nmatrices = nmatrices + 1
  if (present(bset)) nmatrices = nmatrices + 1

  ! FHJ: Compue local contribution to the Frobenius norm of the matrices
  bsedh_fro2 = sum(MYABS2(bsedhead(:,:,:)))
  if (present(bsedwing)) bsedw_fro2 = sum(MYABS2(bsedwing(:,:,:)))
  bsedb_fro2 = sum(MYABS2(bsedbody(:,:,:)))
  bsex_fro2 = sum(MYABS2(bsex(:,:,:)))
  if (present(bset)) bset_fro2 = sum(MYABS2(bset(:,:,:)))

  version = VER_BSE_FORT
#ifdef HDF5
  if(xct%use_hdf5) then
    version = VER_BSE_HDF5
  endif
#endif
  call init_mf_header_from_types(kern%mf, 'KER', SCALARSIZE, version, kp, gvec, syms, crys)
  call xctinfo_to_kernel_header(xct, kg%f, kern, 3)
  if (present(bsedwing)) kern%nmat = kern%nmat + 1
  if (present(bset)) kern%nmat = kern%nmat + 1

#ifdef HDF5
  if(xct%use_hdf5) then
    call timacc(75,1)

    if (peinf%inode==0) then
      call setup_kernel_hdf5('bsemat.h5', kern)
      call write_kernel_header_hdf5('bsemat.h5', kern)
    endif

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

    call timacc(75,2)
    call timacc(76,1)

    if (present(bsedwing)) then
      call bse_hdf5_write(xct,iownsize,1,bsedhead)
      call bse_hdf5_write(xct,iownsize,2,bsedwing)
      call bse_hdf5_write(xct,iownsize,3,bsedbody)
      call bse_hdf5_write(xct,iownsize,4,bsex)
    else
      call bse_hdf5_write(xct,iownsize,1,bsedhead)
      call bse_hdf5_write(xct,iownsize,2,bsedbody)
      call bse_hdf5_write(xct,iownsize,3,bsex)
      call bse_hdf5_write(xct,iownsize,4,bset)
    endif

    call h5close_f(error)
    call timacc(76,2)
  else
#endif

    SAFE_ALLOCATE(bsemt,(xct%nspin,xct%nspin,xct%n1b_co,xct%n2b_co,xct%nkpt_co))
    SAFE_ALLOCATE(bsemtt,(xct%nspin,xct%nspin,xct%n1b_co,xct%n2b_co,xct%nkpt_co))
    if(peinf%inode .eq. 0 ) then
      call open_file(unit=11,file='bsedmat',form='unformatted',status='replace')
      call open_file(unit=12,file='bsexmat',form='unformatted',status='replace')
      if (present(bset)) &
        call open_file(unit=13,file='bsetmat',form='unformatted',status='replace')
    endif
    if (present(bsedwing)) then
      kern%nmat = 3
    else
      kern%nmat = 2
    endif
    call write_binary_kernel_header(11, kern)
    kern%nmat = 1
    call write_binary_kernel_header(12, kern)
    if (present(bset)) then
      kern%nmat = 1
      call write_binary_kernel_header(13, kern)
    endif

    ! FHJ: this is to generate nice output / time estimate
    call progress_init(prog_info, 'writing out BSE matrices', 'records', &
      xct%nkpt_co*nmatrices*xct%n2b_co*xct%n1b_co)
    do imatrix = 1, nmatrices
      do ikp=1,xct%nkpt_co
        do icp=1,xct%n2b_co
          do ivp=1,xct%n1b_co
            call progress_step(prog_info)
            bsemt=0.d0
            do ik=1,xct%nkpt_co
              do ic=1,xct%n2b_co
                do iv=1,xct%n1b_co

                  if (xct%icpar .eq. 0) then
                    it = peinf%wown(1,1,ikp,1,1,ik)
                    if (it .ne. 0) then
                      it = peinf%wown(1,1,ikp,1,1,ik) + xct%n1b_co*xct%n1b_co*xct%n2b_co*(icp-1) &
                        + xct%n1b_co*xct%n1b_co*(ic-1) + xct%n1b_co*(ivp-1) + iv -1
                    endif
                  else if (xct%ivpar .eq. 0) then
                    it = peinf%wown(1,icp,ikp,1,ic,ik)
                    if (it .ne. 0) then
                      it = it + xct%n1b_co*(ivp-1) + iv -1
                    endif
                  else
                    it = peinf%wown(ivp,icp,ikp,iv,ic,ik)
                  endif

                  if (it .ne. 0) then
                    do is1=1,xct%nspin
                      do is2=1,xct%nspin
                        if (present(bsedwing)) then
                          if ( imatrix .eq. 1) then
                            bsem=bsedhead(it,is1,is2)
                          else if ( imatrix .eq. 2) then
                            bsem=bsedwing(it,is1,is2)
                          else if ( imatrix .eq. 3) then
                            bsem=bsedbody(it,is1,is2)
                          else if ( imatrix .eq. 4) then
                            bsem=bsex(it,is1,is2)
                          endif
                        else if (present(bset)) then
                          if ( imatrix .eq. 1) then
                            bsem=bsedhead(it,is1,is2)
                          else if ( imatrix .eq. 2) then
                            bsem=bsedbody(it,is1,is2)
                          else if ( imatrix .eq. 3) then
                            bsem=bsex(it,is1,is2)
                          else if ( imatrix .eq. 4) then
                            bsem=bset(it,is1,is2)
                          endif
                        endif
                        bsemt(is2,is1,iv,ic,ik) = bsem
                      enddo !is2
                    enddo !is1
                  endif
                enddo !iv
              enddo !ic
            enddo !ik
            bsemtt=0D0

            call timacc(75,1)
#ifdef MPI
            call MPI_REDUCE(bsemt(1,1,1,1,1),bsemtt(1,1,1,1,1),xct%nkpt_co*xct%n2b_co*xct%n1b_co*xct%nspin*xct%nspin, &
              MPI_SCALAR,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
! JRD: This barrier may greatly decrease performance. Could hit buffer limit if don`t have it...
            call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
            call timacc(75,2)
            call timacc(76,1)
            if (peinf%inode .eq. 0) then
              if (present(bsedwing)) then
                if ( imatrix .lt. 4) then
                  iunit = 11
                else if (imatrix .eq. 4) then
                  iunit = 12
                endif
              else
                if ( imatrix .lt. 3) then
                  iunit = 11
                else if (imatrix .eq. 3) then
                  iunit = 12
                else if (imatrix .eq. 4) then
                  iunit = 13
                endif
              endif
#ifdef MPI
              write(iunit) ikp,icp,ivp, bsemtt(:,:,:,:,:)
#else
              write(iunit) ikp,icp,ivp, bsemt(:,:,:,:,:)
#endif
            endif
            call timacc(76,2)
          enddo !ivp
        enddo !icp
      enddo !ikp
    enddo !imatrix
    call progress_free(prog_info)

    if(peinf%inode .eq. 0 ) then
      call close_file(11)
      call close_file(12)
      if (present(bset)) call close_file(13)
    endif

#ifdef HDF5
  endif
#endif

! FHJ: Reduce and print the Frobenius norm from all matrices. We use Frobenius
! instead of the max norm because the Frobenius norm is invariant under unitary transf.

#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD,mpierr)
  call MPI_Reduce(bsedh_fro2, bsedh_fro, 1, MPI_REAL_DP, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
  if (present(bsedwing)) &
    call MPI_Reduce(bsedw_fro2, bsedw_fro, 1, MPI_REAL_DP, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Reduce(bsedb_fro2, bsedb_fro, 1, MPI_REAL_DP, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Reduce(bsex_fro2, bsex_fro, 1, MPI_REAL_DP, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
  if (present(bset)) &
    call MPI_Reduce(bset_fro2, bset_fro, 1, MPI_REAL_DP, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
#else
  bsedh_fro = bsedh_fro2
  if (present(bsedwing)) &
    bsedw_fro = bsedw_fro2
  bsedb_fro = bsedb_fro2
  bsex_fro = bsex_fro2
  if (present(bset)) &
    bset_fro = bset_fro2
#endif

  if (peinf%inode.eq.0) then
    write(6,*)
    write(6,'(1x,a)') 'Frobenius norm of the matrices per spin:'
    write(6,'(1x,a,es21.12e4)') '- Head : ', sqrt(bsedh_fro)/xct%nspin
    if (present(bsedwing)) &
      write(6,'(1x,a,es21.12e4)') '- Wing : ', sqrt(bsedw_fro)/xct%nspin
    write(6,'(1x,a,es21.12e4)') '- Body : ', sqrt(bsedb_fro)/xct%nspin
    write(6,'(1x,a,es21.12e4)') '- X    : ', sqrt(bsex_fro)/xct%nspin
    if (present(bset)) &
      write(6,'(1x,a,es21.12e4)') '- T    : ', sqrt(bset_fro)/xct%nspin
    write(6,*)
  endif

  if(.not. xct%use_hdf5) then
    SAFE_DEALLOCATE(bsemt)
    SAFE_DEALLOCATE(bsemtt)
  endif

  POP_SUB(bsewrite)

  return
end subroutine bsewrite

!============================================================================
!
! (2) bse_hdf5_write
!
!============================================================================

#ifdef HDF5
subroutine bse_hdf5_write(xct,iownsize,imatrix,bsemat)
  type (xctinfo), intent(in) :: xct
  integer, intent(in) :: iownsize
  integer, intent(in) :: imatrix
  SCALAR, intent(in) :: bsemat(:,:,:) !< (iownsize*peinf%myown,xct%nspin,xct%nspin)

  integer :: error, rank, ik, ikp, ic, icp, iv, ivp
  integer :: it, ii, iit, is, isp
  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
#ifdef MPI
  integer(HID_T) :: plist_id      ! Property list identifier
#endif

  integer(HSIZE_T) :: count(7), countf(7), offset(7), offsetm(7)

  real(DP), allocatable :: data(:,:,:,:,:,:,:)

  PUSH_SUB(bse_hdf5_write)

! JRD: Initialize rank
  rank = 7

! JRD: Open File and Grab dset_id for this matrix

  call hdf5_open_file('bsemat.h5', 'rw', file_id, parallel_io=.true.)
  if (xct%theory == 0) then
    if (imatrix .eq. 1) then
      call h5dopen_f(file_id, 'mats/head', dset_id, error)
    else if (imatrix .eq. 2) then
      call h5dopen_f(file_id, 'mats/wing', dset_id, error)
    else if (imatrix .eq. 3) then
      call h5dopen_f(file_id, 'mats/body', dset_id, error)
    else if (imatrix .eq. 4) then
      call h5dopen_f(file_id, 'mats/exchange', dset_id, error)
    endif
  else
    if (imatrix .eq. 1) then
      call h5dopen_f(file_id, 'mats/head', dset_id, error)
    else if (imatrix .eq. 2) then
      call h5dopen_f(file_id, 'mats/body', dset_id, error)
    else if (imatrix .eq. 3) then
      call h5dopen_f(file_id, 'mats/exchange', dset_id, error)
    else if (imatrix .eq. 4) then
      call h5dopen_f(file_id, 'mats/fxc', dset_id, error)
    endif
  endif
  ! JRD: Dumb Debugging
  !write(6,*) 'bs', peinf%inode, peinf%mypown, peinf%npown, peinf%nckpe,peinf%myown

#ifdef CPLX
  count(1) = 2
  countf(1) = 2
#else
  count(1) = 1
  countf(1) = 1
#endif

  if (xct%ivpar .eq. 0) then
    count(2) = xct%n1b_co
    count(3) = xct%n1b_co
    countf(2) = xct%n1b_co
    countf(3) = xct%n1b_co
  else
    count(3) = 1
    countf(3) = 1
    countf(2) = peinf%npown
    if (peinf%mypown .eq. 0) then
      count(2) = 1
    else
      count(2) = peinf%mypown
    endif
  endif

  if (xct%icpar .eq. 0) then
    count(4) = xct%n2b_co
    count(5) = xct%n2b_co
    countf(4) = xct%n2b_co
    countf(5) = xct%n2b_co
  else
    count(5) = 1
    countf(5) = 1
    if (xct%ivpar .eq. 1 .or. peinf%mypown .eq. 0) then
      count(4) = 1
    else
      count(4) = peinf%mypown
    endif
    if (xct%ivpar .eq. 1) then
      countf(4) = 1
    else
      countf(4) = peinf%npown
    endif
  endif

  count(7)=1
  countf(7)=1
  if (xct%icpar .eq. 1 .or. xct%ivpar .eq. 1 .or. peinf%mypown .eq. 0) then
    count(6) = 1
  else
    count(6) = peinf%mypown
  endif
  if (xct%icpar .eq. 1 .or. xct%ivpar .eq. 1) then
    countf(6) = 1
  else
    countf(6) = peinf%npown
  endif

  SAFE_ALLOCATE(data,(countf(1),countf(2),countf(3),countf(4),countf(5),countf(6),countf(7)))
  !write(6,*) 'bs2', peinf%inode, count(2),count(4),count(6)
  !write(6,*) 'bs3', peinf%inode, countf(2),countf(4),countf(6)

  do isp = 1, xct%nspin
  do is = 1, xct%nspin

! We  pool the writes here...

!  write(6,*) 'before loop', peinf%inode, peinf%npown, peinf%mypown, peinf%myown

  do ii = 1, (peinf%nckpe/peinf%npown)

! Get filespace
! Probably don't need this is if don't close above

    call h5screate_simple_f(rank, countf, memspace, error)
    call h5dget_space_f(dset_id,filespace,error)

! Construct data and offset

    iit=(ii-1)*peinf%mypown + 1

!    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
!    write(6,*) 'loop', peinf%inode, ii, iit
!    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if (peinf%mypown .ne. 0 .and. iit .le. peinf%myown) then

      it=(iit-1)*iownsize + 1

      do ikp=1,count(7)
      do ik=1,count(6)
      do icp=1,count(5)
      do ic=1,count(4)
      do ivp=1,count(3)
      do iv=1,count(2)
#ifdef CPLX
        data(1,iv,ivp,ic,icp,ik,ikp) = dble(bsemat(it,is,isp))
        data(2,iv,ivp,ic,icp,ik,ikp) = IMAG(bsemat(it,is,isp))
#else
        data(1,iv,ivp,ic,icp,ik,ikp) = bsemat(it,is,isp)
#endif
        it=it+1
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      offsetm(:) = 0
      call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, count, error)

      ik=peinf%ik(peinf%inode+1,iit)
      ikp=peinf%ikp(peinf%inode+1,iit)

      if (xct%icpar .eq. 1) then
        ic=peinf%ic(peinf%inode+1,iit)
        icp=peinf%icp(peinf%inode+1,iit)
      else
        ic=1
        icp=1
      endif

      if (xct%ivpar .eq. 1) then
        iv=peinf%iv(peinf%inode+1,iit)
        ivp=peinf%ivp(peinf%inode+1,iit)
      else
        iv=1
        ivp=1
      endif

      offset(1)=0
      offset(2)=(iv-1)
      offset(3)=(ivp-1)
      offset(4)=(ic-1)
      offset(5)=(icp-1)
      offset(6)=(ik-1)+(is-1)*xct%nkpt_co
      offset(7)=(ikp-1)+(isp-1)*xct%nkpt_co

    else

      offset(:) = 0
      call H5sselect_none_f(memspace,error);

    endif

! Select hyperslab

    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

! JRD: If this proc has no work, nullify filespace

    if (iit .gt. peinf%myown .or. peinf%mypown .eq. 0) call H5sselect_none_f(filespace,error);

!    write(6,*) 'loop3', peinf%inode, ii, iit
!    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

! Create property list for collective dataset write
! Collectively write the file

#ifdef MPI
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, countf, error, memspace, filespace, &
                      xfer_prp = plist_id)
    call h5pclose_f(plist_id, error)
#else
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, countf, error, memspace, filespace)
#endif

!    write(6,*) 'loop4', peinf%inode, ii, iit
!    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

! Close property list and filespace. Do we have to do this each time?

    call h5sclose_f(memspace, error)
    call h5sclose_f(filespace, error)

  enddo
  enddo
  enddo

! Close dataset and file

  call h5dclose_f(dset_id, error)
  call hdf5_close_file(file_id)

  SAFE_DEALLOCATE(data)

  POP_SUB(bse_hdf5_write)

  return

end subroutine bse_hdf5_write

#endif

end module bsewrite_m
