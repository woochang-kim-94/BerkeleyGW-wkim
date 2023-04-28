!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================

#ifdef READ
  #define HDF5_READ_WRITE(x) hdf5_read ## x
  #define READ_WRITE(x) read ## x
  #define INTENT out
  #define FLAVOR_INTENT inout
  #define H5D_READ_WRITE call h5dread_f
  #define H5F_FILE_ACCESS 'r'
#else
  #define HDF5_READ_WRITE(x) hdf5_write ## x
  #define READ_WRITE(x) write ## x
  #define INTENT in
  #define FLAVOR_INTENT in
  #define H5D_READ_WRITE call h5dwrite_f
  #define H5F_FILE_ACCESS 'rw'
#endif

#ifdef HDF5
  #define NAME(x) READ_WRITE(_hdf5 ## x)
#endif

#ifdef TEMP_WFN_DATA
  #ifdef TEMP_COMPLEX
    #define TEMP_SCALAR complex(DPC)
    #define MPI_TEMP_SCALAR MPI_COMPLEX_DPC
    #define LONGNAME(x) NAME(x ## _complex)
  #else
    #define TEMP_SCALAR real(DPC)
    #define MPI_TEMP_SCALAR MPI_REAL_DP
    #define LONGNAME(x) NAME(x ## _real)
  #endif
#endif

! begin read/write header
#ifdef TEMP_HEADER
subroutine NAME(_header_type)(sFileName, sheader, iflavor, kp, gvec, syms, crys)
  character(len=*), intent(in) :: sFileName
  character(len=3), intent(inout) :: sheader
  !> define type. always must be initialized. modified only if -1 on input
  !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
  integer, intent(FLAVOR_INTENT) :: iflavor
  type(kpoints), intent(INTENT) :: kp
  type(gspace), intent(INTENT) :: gvec
  type(symmetry), intent(INTENT) :: syms
  type(crystal), intent(INTENT) :: crys

  ! set values based on epsilon calculation
  logical :: is_get=.false.
  logical :: wfnflag=.true.

  PUSH_SUB(NAME(_header_type))

  if (peinf%inode == 0) then
    call READ_WRITE(_info)(TRUNC(sFileName), iflavor)
    call READ_WRITE(_kpoints)(TRUNC(sFileName), kp)
    call READ_WRITE(_gspace)(TRUNC(sFileName), gvec)
    call READ_WRITE(_symmetry)(TRUNC(sFileName), syms)
    call READ_WRITE(_crystal)(TRUNC(sFileName), crys)
  endif

#ifdef READ
#ifdef MPI
  if (peinf%npes > 1) then
    call MPI_BCAST(iflavor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(is_get, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kp%nspin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kp%nspinor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(gvec%ng, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(syms%ntran, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(syms%cell_symmetry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%nat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(gvec%ecutrho, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(gvec%FFTgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%celvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%alat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%adot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%recvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%blat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%bdot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(syms%mtrx(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(syms%tnp, 3*48, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    if (wfnflag) then
      call MPI_BCAST(kp%nrk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%mnband, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%ngkmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%ecutwfc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%kgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%shift, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    endif
  endif
#endif
  if(peinf%inode > 0 ) then
    SAFE_ALLOCATE(crys%atyp, (crys%nat))
    SAFE_ALLOCATE(crys%apos, (3, crys%nat))
    if (wfnflag) then
      SAFE_ALLOCATE(kp%ngk, (kp%nrk))
      SAFE_ALLOCATE(kp%w, (kp%nrk))
      SAFE_ALLOCATE(kp%rk, (3, kp%nrk))
      SAFE_ALLOCATE(kp%ifmin, (kp%nrk, kp%nspin))
      SAFE_ALLOCATE(kp%ifmax, (kp%nrk, kp%nspin))
      SAFE_ALLOCATE(kp%el, (kp%mnband, kp%nrk, kp%nspin))
      SAFE_ALLOCATE(kp%occ, (kp%mnband, kp%nrk, kp%nspin))
    endif
  endif
#endif

#if defined READ && defined MPI
  if (peinf%npes > 1) then
    call MPI_BCAST(crys%atyp, crys%nat, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%apos, 3*crys%nat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    if (wfnflag) then
      call MPI_BCAST(kp%ngk, kp%nrk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%w, kp%nrk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%rk(1,1), 3*kp%nrk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%ifmin(1,1), size(kp%ifmin), MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%ifmax(1,1), size(kp%ifmax), MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%el(1,1,1), size(kp%el), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%occ(1,1,1), size(kp%occ), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    endif
  endif
#endif
  POP_SUB(NAME(_header_type))
end subroutine NAME(_header_type)

subroutine READ_WRITE(_info)(sFileName, iflavor)
  character(len=*), intent(in) :: sFileName
  !> define type. always must be initialized. modified only if -1 on input
  !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
  integer, intent(FLAVOR_INTENT) :: iflavor

  integer(HID_T) :: hidFile
  integer :: iflavor2
  character(len=16) :: sflavor
  character(len=128) :: routine_name
#ifdef READ
  logical :: exists
#endif

  PUSH_SUB(READ_WRITE(_info))

  call hdf5_open_file(sFileName, H5F_FILE_ACCESS, hidFile)
  iflavor2 = iflavor
  if (iflavor==0) iflavor2 = SCALARSIZE
  routine_name = TOSTRING(READ_WRITE(_info))

#ifdef READ
  ! FHJ: Reading flavor: need some special logic b/c of iflavor==-1
  if (iflavor<-1 .or. iflavor>2) then
    write(sflavor,'(i0)') iflavor
    call die("Illegal value iflavor = " + TRUNC(sflavor) + " passed to " +&
      trim(routine_name) + ": must be -1,0,1,2.")
  endif
  if (iflavor==-1) then
    ! FHJ: read flavor directly into iflavor, make sure value is reasonable
    call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/flavor', iflavor)
    if (iflavor/=1 .and. iflavor/=2) then
      write(sflavor,'(i0)') iflavor
      call die("Illegal flavor = " + TRUNC(sflavor) + " in file " + trim(sFileName))
    endif
  else
    ! FHJ: just check flavor against iflavor2
    call hdf5_require_flavor(hidFile, 'mf_header/flavor', iflavor2, trim(sFileName))
  endif
  call hdf5_require_version(hidFile, 'mf_header/versionnumber', VER_WFN_HDF5, trim(sFileName))

#else
  ! FHJ: Writing flavor: just write adjusted value in iflavor2
  if (iflavor<0 .or. iflavor>2) then
    write(sflavor,'(i0)') iflavor
    call die("Illegal value iflavor = " + TRUNC(sflavor) + " passed to " +&
      trim(routine_name) + ": must be 0,1,2.")
  endif
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/versionnumber', VER_WFN_HDF5)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/flavor', iflavor2)
#endif

  call hdf5_close_file(hidFile)

  POP_SUB(READ_WRITE(_info))

end subroutine READ_WRITE(_info)

subroutine READ_WRITE(_gspace)(sFileName,gvec)
  character(len=*), intent(in) :: sFileName
  type(gspace), intent(INTENT) :: gvec

  integer(HID_T) :: hidFile

  PUSH_SUB(READ_WRITE(_gspace))

  call hdf5_open_file(sFileName, H5F_FILE_ACCESS, hidFile)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/gspace/ng', gvec%ng)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/gspace/ecutrho', gvec%ecutrho)
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/gspace/FFTgrid', (/3/), gvec%FFTgrid)
  call hdf5_close_file(hidFile)

  POP_SUB(READ_WRITE(_gspace))
end subroutine READ_WRITE(_gspace)

subroutine NAME(_gvectors)(sFileName, ng, gvec)
  character(len=*), intent(in) :: sFileName
  integer, intent(in) :: ng !< used size of array
  integer, intent(INTENT) :: gvec(:, :) !< (3, ng_bound)

  integer(HID_T) :: hidFile
  logical :: bcast_, dont_read_

  PUSH_SUB(NAME(_gvectors))

  dont_read_=.false.
  bcast_=.not. dont_read_

  if(peinf%inode == 0) then
    call hdf5_open_file(sFileName, H5F_FILE_ACCESS, hidFile)
    call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/gspace/components', (/3,ng/), gvec)
    call hdf5_close_file(hidFile)
  endif

#ifdef READ
#ifdef MPI
  if(peinf%npes > 1) then
    if(bcast_) then
      call MPI_BCAST(gvec(1,1), 3 * ng, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    endif
  endif
#endif
#endif

  POP_SUB(NAME(_gvectors))
end subroutine NAME(_gvectors)

subroutine READ_WRITE(_symmetry)(sFileName,syms)
  character(len=*), intent(in) :: sFileName
  type(symmetry), intent(INTENT) :: syms

  integer(HID_T) :: hidFile

  PUSH_SUB(READ_WRITE(_symmetry))

  call hdf5_open_file(sFileName, H5F_FILE_ACCESS, hidFile)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/symmetry/ntran', syms%ntran)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/symmetry/cell_symmetry', syms%cell_symmetry)
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/symmetry/mtrx', (/3, 3, 48/), syms%mtrx)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/symmetry/tnp', (/3, 48/), syms%tnp)
  call hdf5_close_file(hidFile)

  POP_SUB(READ_WRITE(_symmetry))
end subroutine READ_WRITE(_symmetry)

subroutine READ_WRITE(_crystal)(sFileName,crys)
  character(len=*), intent(in) :: sFileName
  type(crystal), intent(INTENT) :: crys

  integer(HID_T) :: hidFile

  PUSH_SUB(READ_WRITE(_crystal))

  call hdf5_open_file(sFileName, H5F_FILE_ACCESS, hidFile)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/crystal/celvol', crys%celvol)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/crystal/recvol', crys%recvol)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/crystal/alat', crys%alat)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/crystal/blat', crys%blat)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/crystal/nat', crys%nat)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/crystal/avec', (/3, 3/), crys%avec)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/crystal/bvec', (/3, 3/), crys%bvec)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/crystal/adot', (/3, 3/), crys%adot)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/crystal/bdot', (/3, 3/), crys%bdot)

#ifdef READ
  SAFE_ALLOCATE(crys%atyp, (crys%nat))
  SAFE_ALLOCATE(crys%apos, (3,crys%nat))
#endif

  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/crystal/atyp', (/crys%nat/), crys%atyp)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/crystal/apos', (/3, crys%nat/), crys%apos)
  call hdf5_close_file(hidFile)

  POP_SUB(READ_WRITE(_crystal))
end subroutine READ_WRITE(_crystal)

subroutine READ_WRITE(_kpoints)(sFileName,kp)
  character(len=*), intent(in) :: sFileName
  type(kpoints), intent(INTENT) :: kp

  integer(HID_T) :: hidFile

  PUSH_SUB(READ_WRITE(_kpoints))

  call hdf5_open_file(sFileName, H5F_FILE_ACCESS, hidFile)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/kpoints/nspin', kp%nspin)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/kpoints/nspinor', kp%nspinor)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/kpoints/nrk', kp%nrk)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/kpoints/mnband', kp%mnband)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/kpoints/ngkmax', kp%ngkmax)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/kpoints/ecutwfc', kp%ecutwfc)
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/kpoints/kgrid', (/3/), kp%kgrid)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/kpoints/shift', (/3/), kp%shift)

#ifdef READ
  SAFE_ALLOCATE(kp%ngk, (kp%nrk))
  SAFE_ALLOCATE(kp%ifmin, (kp%nrk, kp%nspin))
  SAFE_ALLOCATE(kp%ifmax, (kp%nrk, kp%nspin))
  SAFE_ALLOCATE(kp%w, (kp%nrk))
  SAFE_ALLOCATE(kp%rk, (3,kp%nrk))
  SAFE_ALLOCATE(kp%el, (kp%mnband,kp%nrk,kp%nspin))
  SAFE_ALLOCATE(kp%occ, (kp%mnband,kp%nrk,kp%nspin))
#endif

  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/kpoints/ngk', &
    (/kp%nrk/), kp%ngk)
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/kpoints/ifmin', &
    (/kp%nrk, kp%nspin/), kp%ifmin)
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/kpoints/ifmax', &
    (/kp%nrk, kp%nspin/), kp%ifmax)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/kpoints/w', &
    (/kp%nrk/), kp%w)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/kpoints/rk', &
    (/3, kp%nrk/), kp%rk)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/kpoints/el', &
    (/kp%mnband, kp%nrk, kp%nspin/), kp%el)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/kpoints/occ', &
    (/kp%mnband, kp%nrk, kp%nspin/), kp%occ)
  call hdf5_close_file(hidFile)

  POP_SUB(READ_WRITE(_kpoints))
end subroutine READ_WRITE(_kpoints)
#endif
! end read/write header

! begin read/write wfn gvectors
#ifdef TEMP_WFN_GVEC
subroutine NAME(_wfn_gvectors)(fname, gvec, ngktot)
  character(len=*) :: fname
  integer, intent(inout) :: gvec(:,:)
  integer, intent(in) :: ngktot

  integer(HID_T) :: file_id

  PUSH_SUB(NAME(_wfn_gvectors))

  if(peinf%inode == 0) then
    call hdf5_open_file(fname, H5F_FILE_ACCESS, file_id)
    call HDF5_READ_WRITE(_int_array)(file_id, 'wfns/gvecs', (/3, ngktot/), gvec)
    call hdf5_close_file(file_id)
  endif

#if defined READ && defined MPI
  if(peinf%npes > 1) then
    call MPI_Bcast(gvec(1,1), 3*ngktot, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  endif
#endif
  POP_SUB(NAME(_wfn_gvectors))
end subroutine NAME(_wfn_gvectors)
#endif
! end read/write wfn gvectors

! begin read/write wfn data
#ifdef TEMP_WFN_DATA
subroutine LONGNAME(_band)(fname, wfn, ngk, nstot, ioffsetk, ioffsetb)
  character(len=*) :: fname
  TEMP_SCALAR, intent(INTENT) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
  integer, intent(in) :: ngk
  integer, intent(in) :: nstot
  integer, intent(in) :: ioffsetk
  integer, intent(in) :: ioffsetb

#ifdef TEMP_COMPLEX
  real(DP) :: dwfn(2,ngk,nstot,1)
#else
  real(DP) :: dwfn(1,ngk,nstot,1)
#endif
  integer(HID_T) :: file_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: dataspace_id
  integer(HID_T) :: memspace_id
  integer(HSIZE_T) :: a3(4), offset(4), count(4)
  integer :: error, is

  PUSH_SUB(LONGNAME(_band))

#ifdef TEMP_COMPLEX
  a3(1) = 2
#else
  a3(1) = 1
#endif
  a3(2) = ngk
  a3(3) = nstot
  a3(4) = 1
  offset(1) = 0
  offset(2) = ioffsetk
  offset(3) = 0
  offset(4) = ioffsetb

#ifdef TEMP_COMPLEX
  count(1) = 2
#else
  count(1) = 1
#endif
  count(2) = ngk
  count(3) = nstot
  count(4) = 1

#ifndef READ
#ifdef TEMP_COMPLEX
  dwfn(1,:,:,1) = real(wfn(:,:))
  dwfn(2,:,:,1) = IMAG(wfn(:,:))
#else
  dwfn(1,:,:,1) = wfn(:,:)
#endif
#endif

  if(peinf%inode == 0) then
    call hdf5_open_file(fname, H5F_FILE_ACCESS, file_id)
    call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
    CALL h5screate_simple_f(4, count, memspace_id, error)
    call h5dget_space_f(dataset_id, dataspace_id, error)

    call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)

    H5D_READ_WRITE(dataset_id, H5T_NATIVE_DOUBLE, dwfn, a3, error, file_space_id = dataspace_id, mem_space_id = memspace_id)
    call h5dclose_f(dataset_id, error)
    call h5sclose_f(dataspace_id, error)
    call h5sclose_f(memspace_id, error)
    call hdf5_close_file(file_id)
  endif

#ifdef READ
#ifdef TEMP_COMPLEX
  wfn(:,:) = CMPLX(dwfn(1,:,:,1), dwfn(2,:,:,1))
#else
  wfn(:,:) = dwfn(1,:,:,1)
#endif
#endif

  POP_SUB(LONGNAME(_band))

end subroutine LONGNAME(_band)
#endif
! end read/write wfn data

#ifdef TEMP_OTHER
subroutine read_hdf5_bands_block(file_id, kp, nbownmax, nbownactual, does_it_ownb, ib_first, wfnsout, ioffset, is_sigma, is_bse)
  integer(HID_T), intent(in) :: file_id
  type(kpoints), intent(in) :: kp
  integer, intent(in) :: nbownmax
  integer, intent(in) :: nbownactual !< how many bands I own
  logical, intent(in) :: does_it_ownb(:,:)
  integer, intent(in) :: ib_first !< first band of the set of bands I own
  SCALAR, intent(out) :: wfnsout(:,:,:) !< (SUM(kp%ngk), kp%nspin*kp%nspinor, nbownactual)
  integer, optional, intent(in) :: ioffset
  logical, optional, intent(in) :: is_sigma, is_bse

  real(DP), allocatable :: wfndata(:,:,:,:)
  integer(HID_T) :: plist_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dataspace
  integer(HID_T) :: memspace
  integer(HSIZE_T) :: count(4), offset(4)
  integer :: error
  integer :: ipe, reader
  integer :: nread
  integer :: ngktot
  integer :: ioffset_
  integer, allocatable :: ranks(:)
  integer :: icount
  integer :: mpiworldgrp, mpigrp, bandcomm, ib, ib_, max_bands_read, bands_read, bands_read_max, is
  integer, parameter :: max_bytes_read = 536870912 ! don`t read/send more than 1/2 GB at a time
  ! FHJ: comm_style is now controled by the BGW_HDF5_BANDS_COMM_STYLE.
  ! 0=native HDF5; 1=manual group comm; 2=manual send/recvs
  integer :: comm_style
  logical :: do_read
  !real(DP) :: mem
  !integer :: nmpinode

  ! variable for HDF5 in sigma 
  logical :: my_is_sigma, my_is_bse
  integer :: inter_pool_rank, pool_rank

  PUSH_SUB(read_hdf5_bands_block)

  call logit('Reading HDF5 WFNs by blocks')

  comm_style = peinf%hdf5_bands_comm_style

  my_is_sigma = .false.
  if ( present ( is_sigma ) ) my_is_sigma = is_sigma
  my_is_bse = .false.
  if ( present ( is_bse ) )   my_is_bse = is_bse
  if ( my_is_bse .and. my_is_sigma ) then
    call die("read_hdf5_bands_block inconsistent flags: BSE and Sigma read mode both active", &
             only_root_writes=.true.)
  end if

  SAFE_ALLOCATE(ranks,(peinf%npes))
  ioffset_ = 0
  if(present(ioffset)) then
    ioffset_ = ioffset
  endif
  ngktot = SUM(kp%ngk)
  nread = peinf%npools

  call h5dopen_f(file_id, 'wfns/coeffs', dset_id, error)
  call h5dget_space_f(dset_id, dataspace, error)

  ! set comm style for sigma case
  if ( my_is_sigma )then
    ! only pool zero read
    reader = -1
    if (nbownactual>0) then
      if ( peinf%my_pool == 0) reader = peinf%inode
    end if
    ! create the interpool communicator
    inter_pool_rank = peinf%my_pool
    pool_rank = peinf%pool_rank
    if ( pool_rank < 0 .or. inter_pool_rank < 0 ) then
      ! this ensure that even if the number of pools are 
      ! not exactly dividing the total npes, the leftover 
      ! MPI tasks are placed into a pool
      inter_pool_rank = peinf%inode/peinf%npes_pool
      pool_rank = mod(peinf%inode,peinf%npes_pool)
    end if
#ifdef MPI
    call MPI_Comm_split(MPI_COMM_WORLD, pool_rank, inter_pool_rank, bandcomm, mpierr)
#endif

    ! use comm_style 1 (broadcast)
    comm_style = 1
  elseif( my_is_bse ) then 
    ! all read
    reader = -1
    if (nbownactual>0) then
      reader = peinf%inode
    end if
    comm_style = 0
  else

    ! find lowest rank PE that owns the first band of this group of bands
    reader = -1
    if (nbownactual>0) then
      do ipe = 1, peinf%npes
        if(does_it_ownb(ib_first,ipe)) then
          reader = ipe - 1
          exit
        endif
      enddo
      if (reader==-1) call die("Cannot find reader in read_hdf5_bands_block", only_root_writes=.true.)
  
#ifdef MPI
      if (comm_style==1) then
        ! We use MPI_Bcast with MPI groups
        icount = 0
        do ipe = 1, peinf%npes
          if(does_it_ownb(ib_first,ipe)) then
              icount = icount + 1
              ranks(icount) = ipe - 1
          endif
        enddo
        call MPI_Comm_Group(MPI_COMM_WORLD, mpiworldgrp, mpierr)
        call MPI_Group_Incl(mpiworldgrp, icount, ranks, mpigrp, mpierr)
        call MPI_Comm_Create(MPI_COMM_WORLD, mpigrp, bandcomm, mpierr)
      endif
#endif
    else
#ifdef MPI
      if (comm_style==1) then
        ! FHJ: Note that MPI_Comm_Create must be called by everyone in MPI_COMM_WORLD!
        call MPI_Comm_Create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, bandcomm, mpierr)
      endif
#endif
    endif

  endif ! my_is_sigma

  ! Write some infos
  if ( peinf%inode == 0 ) then
    if ( peinf%use_wfn_hdf5_independent ) then 
      write(6,*) 'Reading using H5FD_MPIO_INDEPENDENT.'
    else 
      write(6,*) 'Reading using H5FD_MPIO_COLLECTIVE.'
    end if 
    write(6,"(1x,A,i3)") 'Using communication style:', comm_style
    write(6,*) 
  end if 

  ! FHJ: read at most max_bands_read bands to avoid MPI/HDF5 buffer overflow.
  ! Here, SCALARSIZE*kp%nspin*kp%nspinor*dble(ngktot)*8d0 is the size of a
  ! single band, including G-vectors from all k-points.
  max_bands_read = min(nbownmax, &
    int(max_bytes_read/(SCALARSIZE*kp%nspin*kp%nspinor*dble(ngktot)*8d0)))
  if (max_bands_read==0) then
    max_bands_read = 1
    if (peinf%inode==0) then
      write(0,*)
      write(0,'(a)') 'WARNING: could not honor limit of '
      write(0,'(f0.3,a)') max_bytes_read/1024d0/1024d0,' MB per chunk when'
      write(0,'(a)') 'reading HDF5 WFN file. Using chunks of '
      write(0,'(f0.3,a)') (kp%nspin*kp%nspinor*SCALARSIZE*dble(ngktot)*8d0)/1024d0/1024d0,' MB.'
      write(0,*)
    endif
  endif
  !write(6,*) 'max_bands_read', max_bands_read
  SAFE_ALLOCATE(wfndata, (SCALARSIZE, ngktot, kp%nspin*kp%nspinor, max_bands_read))

  ib = 1
  do while (ib<=nbownmax)
    bands_read = max(min(nbownactual, ib-1 + max_bands_read) - ib + 1, 0)
    bands_read_max = max(min(nbownmax, ib-1 + max_bands_read) - ib + 1, 0)
    count(1) = SCALARSIZE
    count(2) = ngktot
    count(3) = kp%nspin*kp%nspinor
    count(4) = bands_read
    call h5screate_simple_f(4, count, memspace, error)
    do_read = bands_read>0.and.(peinf%inode==reader.or.comm_style==0)

    if (do_read) then
      offset(1) = 0
      offset(2) = 0
      offset(3) = 0
      offset(4) = (ib_first-1)+ioffset_+(ib-1)
      if (peinf%verb_debug .and. peinf%inode==reader) then
        write(6,'(4(a,i0,1x))') 'ib=',ib,'bands_read=',bands_read,'offset=',offset(3),'ib_first=',ib_first
      endif
    else
      offset(:) = 0
      call H5sselect_none_f(memspace,error)
    endif

    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
    if (.not.do_read) then
      call H5sselect_none_f(dataspace,error)
    endif

    if (peinf%verb_debug .and. peinf%inode==reader) then
      write(6,'(a,i0,a)') 'ib=',ib,' before read!'
      !call procmem(mem,nmpinode)
      !write(6,'(a,f0.3,a)') 'Memory available: ', mem/(1024d0**2),' MB per MPI rank.'
    endif

#ifdef MPI
    !if (peinf%inode==reader) write(6,'(a)') '[1]'
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    !if (peinf%inode==reader) write(6,'(a)') '[2]'
    if ( peinf%use_wfn_hdf5_independent ) then
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    else
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    end if
    !if (peinf%inode==reader) write(6,'(a)') '[3]'
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, wfndata(:,:,:,:), count, error, memspace, dataspace, xfer_prp=plist_id)
    !if (peinf%inode==reader) write(6,'(a)') '[4]'
    call h5pclose_f(plist_id, error)
    !if (peinf%inode==reader) write(6,'(a)') '[5]'
#else
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, wfndata(:,:,:,:), count, error, memspace, dataspace)
#endif

    if (peinf%verb_debug .and. peinf%inode==reader) then
      write(6,'(a,i0,a)') 'ib=',ib,' after read!'
    endif

    call h5sclose_f(memspace, error)
    if (do_read) then
#ifdef DEBUG
      if (peinf%verb_max .and. peinf%inode==reader) then
        do ib_ = 1, bands_read
          write(6,'(" band = ",i6,"; avg(norm) = ", f12.6)') &
            offset(3) + ib_, sum(wfndata(:,:,:,ib_)**2)/dble(kp%nrk)
        enddo
        FLUSH(6)
      endif
#endif
      do is = 1, kp%nspin*kp%nspinor
        wfnsout(:,is,ib:ib+bands_read-1) = &
          SCALARIFY2(wfndata(1,:,is,1:bands_read),wfndata(2,:,is,1:bands_read))
      enddo
    endif

#ifdef MPI
    if (bands_read>0) then
      ! FHJ: No manual distribution is necessary for comm_style==0
      if (comm_style>0) call logit('Sending bands')
      if (comm_style==2) then
        if (peinf%inode==reader) then
          do ipe = 1, peinf%npes
            if(does_it_ownb(ib_first,ipe) .and. ipe-1 .ne. peinf%inode) then
              call MPI_Send(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, &
                MPI_SCALAR, ipe-1, 0, MPI_COMM_WORLD, mpierr)
            endif
          enddo
        else
          call MPI_Recv(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, &
            MPI_SCALAR, reader, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
        endif
      elseif (comm_style==1) then
        call MPI_Bcast(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, &
          MPI_SCALAR, 0, bandcomm, mpierr)
      endif
    endif
#endif
    ib = ib + bands_read_max
  enddo

#ifdef MPI
  if (comm_style==1.and.nbownactual>0) then
    call MPI_Comm_free(bandcomm, mpierr)
    if ( .not. my_is_sigma ) call MPI_Group_free(mpigrp, mpierr)
  endif
#endif
  SAFE_DEALLOCATE(ranks)
  SAFE_DEALLOCATE(wfndata)

  call h5sclose_f(dataspace, error)
  call h5dclose_f(dset_id, error)
  POP_SUB(read_hdf5_bands_block)
end subroutine read_hdf5_bands_block
#endif
! read/write other

#undef READ_WRITE
#undef INTENT
#undef FLAVOR_INTENT
#undef NAME
#undef TEMP_SCALAR
#undef MPI_TEMP_SCALAR
#undef LONGNAME

#undef HDF5_READ_WRITE
#undef H5D_READ_WRITE
#undef H5F_FILE_ACCESS
