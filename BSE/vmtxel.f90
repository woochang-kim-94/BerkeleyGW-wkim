!==============================================================================
!
! Module vmtxel_m
!
! Originally by GKA (2018)
!
! Objects holding dipole operator.
!
!==============================================================================

#include "f_defs.h"

module vmtxel_m

#ifdef HDF5
  use hdf5
  use hdf5_io_m
#endif

  use global_m
  use mtxel_optical_m
  use misc_m, only: bse_index

  implicit none

public

! =========================================================================== !

!> dipole operator (velocity)
type vmtxel_t

  ! Objects
  type (crystal) :: crys
  type (gspace) :: gvec

  ! Dimensions
  integer :: ns       !< Number of spins
  integer :: nk       !< Number of k-points
  integer :: nband    !< Number of bands at k   (conduction bands in BSE)
  integer :: mband    !< Number of bands at k+q (valence bands in BSE)
  integer :: nmat     !< Flat dimension (ns * nk * nband * mband)
  integer :: npol     !< Number of polarization

  ! Options
  integer :: band_ordering  !< 0 = counting bands from fermi levels, with
                            !!     conduction bands going up in energy
                            !!     valence bands down in energy.
                            !! 1 = counting all bands from bottom-up,
                            !!     starting at the first valence band
                            !!     and going up in energy.
  integer :: opr  !< 0 = use velocity operator
                  !! 1 = use momentum operator

  logical :: has_velocity = .false.   !< Whether the velocity operator has
                                      !! been computed along with the
                                      !! dipole operator

  logical :: use_hdf5 = .true.  !< Use hdf5 for I/O

  ! Parallel info
  logical :: is_master

  ! Arrays

  !> The list of k-points
  !> kpt(3,nkpt)
  real(dp), allocatable :: kpt(:,:)

  !> The polarization vector.
  !> pol(3,npol)
  real(dp), allocatable :: pol(:,:)

  !> Dipole operator, for a single k-point, in band index basis.
  !! s1k(nband,mband,ns,npol)
  SCALAR, allocatable :: s1k(:,:,:,:)

  !> Dipole operator, for all k-point, in flat BSE indices basis.
  !! s1(nmat,npol)
  SCALAR, allocatable :: s1(:,:)

  !> Dipole operator, for all k-point, in band index basis.
  !! At the moment, this is treated as a temporary array for io only.
  !! rnm(nband,mband,nk,ns,npol)
  SCALAR, allocatable :: rnm(:,:,:,:,:)

  !> Velocity operator, for all k-point, in band index basis.
  !! vnm(nband,mband,nk,ns,npol)
  SCALAR, allocatable :: vnm(:,:,:,:,:)

  contains

  ! Core procedures
  procedure :: init => vmtxel_init
  procedure :: init_from_xctinfo => vmtxel_init_from_xctinfo
  procedure :: alloc => vmtxel_alloc
  procedure :: free => vmtxel_free

  ! Computation / communication
  procedure :: compute_ik_vmtxel
  procedure :: reduce => reduce_vmtxel
  procedure :: band_to_flat => vmtxel_band_to_flat
  procedure :: flat_to_band => vmtxel_flat_to_band

  ! I/O
  procedure :: write_vmtxel
  procedure :: read_vmtxel
  procedure :: write_vmtxel_bin
  procedure :: read_vmtxel_bin

#ifdef HDF5
  procedure :: write_vmtxel_hdf5
  procedure :: read_vmtxel_hdf5
  procedure :: create_and_write_vmtxel_header_hdf5
  procedure :: write_vmtxel_arrays_hdf5
  procedure :: read_and_broadcast_vmtxel_header_hdf5
  procedure :: read_vmtxel_header_hdf5
  procedure :: broadcast_vmtxel_header
  procedure :: read_and_broadcast_vmtxel_data_hdf5
#endif

end type vmtxel_t

! =========================================================================== !
contains
! =========================================================================== !

!> Initialize object, setting manually dimensions and options
subroutine vmtxel_init(this, ns, nk, nband, mband, opr, npol, &
                       band_ordering, with_velocity)
  class(vmtxel_t), intent(inout) :: this
  integer, intent(in) :: ns       !< Number of spins
  integer, intent(in) :: nk       !< Number of k-points
  integer, intent(in) :: nband    !< Number of bands at k (conduction)
  integer, intent(in) :: mband    !< Number of bands at k+q (valence)
  integer, intent(in),optional :: opr     !< 0 = use velocity operator
                                          !! 1 = use momentum operator
  integer, intent(in), optional :: npol   !< Number of polarization
  integer, intent(in), optional :: band_ordering  !< 0 = from fermi levels
                                                  !! 1 = bottom-up
  logical, intent(in), optional :: with_velocity  !< Do compute the velocity
                                                  !! operator
  !real(dp), intent(in), optional :: pol(:,:)   !< Polarization vectors

  PUSH_SUB(vmtxel_init)

  ! Get parallel info
  this%is_master = (peinf%inode.eq.0)

  this%ns = ns
  this%nk = nk
  this%nband = nband
  this%mband = mband
  this%nmat = this%ns * this%nk * this%nband * this%mband

  if (present(opr)) then
    this%opr = opr
  else
    this%opr = 1
  end if

  if (present(npol)) then
    this%npol = npol
  else
    if (this%opr .eq. 1) then
      this%npol = 3
    else
      this%npol = 1
    end if
  end if

  if (present(band_ordering)) then
    this%band_ordering = band_ordering
  else
    this%band_ordering = 0
  end if

  if (present(with_velocity)) then
    this%has_velocity = with_velocity
  else
    this%has_velocity = .false.
  end if

  POP_SUB(vmtxel_init)

end subroutine vmtxel_init

! =========================================================================== !

!> Initialize object from an xctinfo object, copying dimensions and options
!! as well as the polarization vector
subroutine vmtxel_init_from_xctinfo(this, xct, opr)
  class(vmtxel_t), intent(inout) :: this
  type(xctinfo), intent(in) :: xct
  integer, intent(in),optional :: opr     !< 0 = use velocity operator
                                          !! 1 = use momentum operator
  integer :: opr_ = 0

  PUSH_SUB(vmtxel_init_from_xctinfo)

  if (present(opr)) opr_ = opr

  call this%init(xct%nspin, xct%nkpt_fi, xct%ncb_fi, xct%nvb_fi, &
                 opr=opr_, npol=xct%npol)

  ! Set polarization vector
  if (this%npol .eq. 1) then
    SAFE_ALLOCATE(this%pol, (3, this%npol))
    this%pol(:,1) = xct%pol
  end if

  POP_SUB(vmtxel_init_from_xctinfo)

end subroutine vmtxel_init_from_xctinfo

! =========================================================================== !

!> Allocate arrays
subroutine vmtxel_alloc(this, rnm)
  class(vmtxel_t), intent(inout) :: this
  logical, intent(in), optional :: rnm

  integer :: ipol

  PUSH_SUB(vmtxel_alloc)

  SAFE_ALLOCATE(this%kpt, (3, this%nk))
  SAFE_ALLOCATE(this%s1, (this%nmat, this%npol))
  SAFE_ALLOCATE(this%s1k, (this%nband, this%mband, this%ns, this%npol))
  this%s1 = ZERO
  this%s1k = ZERO

  if (.not. allocated(this%pol)) then
    SAFE_ALLOCATE(this%pol, (3, this%npol))
    do ipol=1,this%npol
      this%pol(:,ipol) = ZERO
      this%pol(ipol,ipol) = ONE
    end do
  end if

  if (present(rnm)) then
    if (rnm) then
      SAFE_ALLOCATE(this%rnm, (this%nband, this%mband, this%nk, this%ns, this%npol))
    end if
  end if

  if (this%has_velocity) then
    SAFE_ALLOCATE(this%vnm, (this%nband, this%mband, this%nk, this%ns, this%npol))
  end if

  POP_SUB(vmtxel_alloc)

end subroutine vmtxel_alloc

! =========================================================================== !

!> Free memory
subroutine vmtxel_free(this)
  class(vmtxel_t), intent(inout) :: this

  PUSH_SUB(vmtxel_free)

  SAFE_DEALLOCATE(this%kpt)
  SAFE_DEALLOCATE(this%pol)
  SAFE_DEALLOCATE(this%s1k)
  SAFE_DEALLOCATE(this%s1)
  SAFE_DEALLOCATE(this%rnm)

  POP_SUB(vmtxel_free)

end subroutine vmtxel_free

! =========================================================================== !

!> Share matrix elements among all PEs
subroutine reduce_vmtxel(this)
  class(vmtxel_t), intent(inout) :: this

  SCALAR, allocatable :: dummy(:,:)

  PUSH_SUB(reduce_vmtxel)

  SAFE_ALLOCATE(dummy, (this%nmat, this%npol))
  dummy = this%s1
#ifdef MPI
  call MPI_ALLREDUCE(dummy(1,1), this%s1(1,1), size(dummy), &
                     MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD, mpierr)
#endif
  SAFE_DEALLOCATE(dummy)

  POP_SUB(reduce_vmtxel)

end subroutine reduce_vmtxel

! =========================================================================== !

!> Transform the dipole operator from flat BSE indices to band indices
subroutine vmtxel_flat_to_band(this)
  class(vmtxel_t), intent(inout) :: this

  type(xctinfo) :: xct_
  integer :: ipol, is, ik, ic, iv

  PUSH_SUB(vmtxel_flat_to_band)

  if (.not. allocated(this%rnm)) then
    SAFE_ALLOCATE(this%rnm, (this%nband, this%mband, this%nk, this%ns, this%npol))
  end if

  ! GKA:  FIXME Do this more elegantly
  xct_%nspin = this%ns
  xct_%ncb_fi = this%nband
  xct_%nvb_fi = this%mband

  do is=1,this%ns
    do ik=1,this%nk
      do iv=1,this%mband
        do ic=1,this%nband
          this%rnm(ic,iv,ik,is,:) = this%s1(bse_index(ik, ic, iv, is, xct_),:)
        enddo
      enddo
    enddo
  enddo

  POP_SUB(vmtxel_flat_to_band)

end subroutine vmtxel_flat_to_band


!> Transform the dipole operator from band indices to flat BSE indices
subroutine vmtxel_band_to_flat(this)
  class(vmtxel_t), intent(inout) :: this

  type(xctinfo) :: xct_
  integer :: ipol, is, ik, ic, iv

  PUSH_SUB(vmtxel_band_to_flat)

  if (.not. allocated(this%s1)) then
    SAFE_ALLOCATE(this%s1, (this%nmat, this%npol))
  end if

  ! GKA:  FIXME Do this more elegantly
  xct_%nspin = this%ns
  xct_%ncb_fi = this%nband
  xct_%nvb_fi = this%mband

  do is=1,this%ns
    do ik=1,this%nk
      do iv=1,this%mband
        do ic=1,this%nband
          this%s1(bse_index(ik, ic, iv, is, xct_),:) = this%rnm(ic,iv,ik,is,:)
        enddo
      enddo
    enddo
  enddo

  POP_SUB(vmtxel_band_to_flat)

end subroutine vmtxel_band_to_flat

! =========================================================================== !

!> Compute the dipole operator
subroutine compute_ik_vmtxel(this, ik, wfnc_fi, wfnvq_fi, gvec, qshift, &
                             crys, eqp)
  class(vmtxel_t), intent(inout) :: this
  integer, intent(in) :: ik
  type (gspace), intent(in) :: gvec
  type (wavefunction), intent(in) :: wfnc_fi
  type (wavefunction), intent(in) :: wfnvq_fi
  real(DP), intent(in), optional :: qshift
  type (crystal), intent(in), optional :: crys
  type (eqpinfo), intent(in), optional :: eqp

  type(xctinfo) :: xct_
  integer :: ipol, is, ic, iv
  real(DP) :: de

  PUSH_SUB(compute_ik_vmtxel)

  if (this%npol==1) then
    if (this%opr.eq.0) then
      call mtxel_v(wfnc_fi,wfnvq_fi,gvec,qshift,this%nband,this%mband, &
                   this%s1k(:,:,:,1))
    elseif (this%opr.eq.1) then
      call mtxel_m(crys,wfnc_fi,wfnvq_fi,gvec,eqp,this%pol(:,1), &
                   this%nband,this%mband,this%s1k(:,:,:,1),ik,.true.)
    endif
  else
    do ipol=1,3
      this%pol(:,ipol) = ZERO
      this%pol(ipol,ipol) = ONE
      if (this%opr.eq.0) then
        call mtxel_v(wfnc_fi,wfnvq_fi,gvec,qshift,this%nband,this%mband, &
                     this%s1k(:,:,:,ipol))
      else
        call mtxel_m(crys,wfnc_fi,wfnvq_fi,gvec,eqp,this%pol(:,ipol), &
                     this%nband,this%mband,this%s1k(:,:,:,ipol),ik,.true.)
      endif
    enddo
  endif

  ! GKA:  FIXME Do this more elegantly
  xct_%nspin = this%ns
  xct_%ncb_fi = this%nband
  xct_%nvb_fi = this%mband

  do is=1,this%ns
    do ic=1,this%nband
      do iv=1,this%mband
        this%s1(bse_index(ik, ic, iv, is, xct_),:) = this%s1k(ic,iv,is,:)
      enddo
    enddo
  enddo

  ! Also compute the velocity operator, along with the dipole operator
  if (this%has_velocity .and. present(eqp)) then

    do is=1,this%ns
      do ic=1,this%nband
        do iv=1,this%mband
          de = (eqp%eclda(ic,ik,is) - eqp%evlda(iv,ik,is))
          this%vnm(ic,iv,ik,is,:) = this%s1k(ic,iv,is,:) * de
        enddo
      enddo
    enddo

  end if


  POP_SUB(compute_ik_vmtxel)

end subroutine compute_ik_vmtxel

! =========================================================================== !
! Writing routines
! =========================================================================== !


!> Write the file, using the format determined with use_hdf5
subroutine write_vmtxel(this)
  class(vmtxel_t), intent(inout) :: this

  PUSH_SUB(write_vmtxel)

#ifdef HDF5
  if (this%use_hdf5) then
    call this%write_vmtxel_hdf5()
  else
    call this%write_vmtxel_bin()
  end if
#else
    call this%write_vmtxel_bin()
#endif

  POP_SUB(write_vmtxel)

end subroutine write_vmtxel

! =========================================================================== !

!> Write binary file
subroutine write_vmtxel_bin(this)
  class(vmtxel_t), intent(inout) :: this

  integer :: ipol
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)

  PUSH_SUB(write_vmtxel_bin)

  if (this%is_master) then
    write(6,'(1x,a)') 'Writing matrix elements into vmtxel'
    do ipol=1,this%npol
      if (this%npol==1) then
        fname = 'vmtxel'
      else
        fname = 'vmtxel_'//suffix(ipol)
      endif
      call open_file(16, file=trim(fname), form='unformatted', status='replace')
      write(16) this%nk, this%nband, this%mband, this%ns, this%opr
      write(16) this%s1(:,ipol)
      call close_file(16)
    enddo

  endif

  POP_SUB(write_vmtxel_bin)

end subroutine write_vmtxel_bin

! =========================================================================== !

#ifdef HDF5

!> Write hdf5 file (serial)
subroutine write_vmtxel_hdf5(this)
  class(vmtxel_t), intent(inout) :: this

  character(len=128) :: fname
  integer(HID_T) :: file_id

  PUSH_SUB(write_vmtxel_hdf5)

  fname = 'vmtxel.h5'

  if (this%is_master) then

    ! Create the file and write dimensions
    call this%create_and_write_vmtxel_header_hdf5(fname)

    ! Write the arrays that are not distributed
    call this%write_vmtxel_arrays_hdf5(fname)

  end if

  POP_SUB(write_vmtxel_hdf5)

end subroutine write_vmtxel_hdf5

! =========================================================================== !

!> Create the file and write dimensions.
subroutine create_and_write_vmtxel_header_hdf5(this, fname)
  class(vmtxel_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id

  PUSH_SUB(create_and_write_vmtxel_header_hdf5)

  !call h5open_f(error)

  ! Create a new file using default properties.
  call hdf5_create_file(trim(fname), file_id)

  ! Create the groups
  call hdf5_create_group(file_id, 'vmtxel_header')
  call hdf5_create_group(file_id, 'vmtxel_data')

  ! Write dimensions
  call hdf5_write_int(file_id, 'vmtxel_header/ns', this%ns)
  call hdf5_write_int(file_id, 'vmtxel_header/nk', this%nk)
  call hdf5_write_int(file_id, 'vmtxel_header/nband', this%nband)
  call hdf5_write_int(file_id, 'vmtxel_header/mband', this%mband)
  call hdf5_write_int(file_id, 'vmtxel_header/nmat', this%nmat)
  call hdf5_write_int(file_id, 'vmtxel_header/npol', this%npol)
  call hdf5_write_int(file_id, 'vmtxel_header/band_ordering', this%band_ordering)
  call hdf5_write_int(file_id, 'vmtxel_header/opr', this%opr)
  call hdf5_write_logical(file_id, 'vmtxel_header/has_velocity', this%has_velocity)

  call hdf5_close_file(file_id)

  !call h5close_f(error)

  POP_SUB(create_and_write_vmtxel_header_hdf5)

end subroutine create_and_write_vmtxel_header_hdf5

! =========================================================================== !

!> Write the arrays into the hdf5 file.
subroutine write_vmtxel_arrays_hdf5(this, fname)
  class(vmtxel_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id

  PUSH_SUB(write_vmtxel_arrays_hdf5)

  !call h5open_f(error)

  call hdf5_open_file(trim(fname), 'rw', file_id)

  ! -------------------------------
  ! Write arrays
  call hdf5_write_double_array(file_id, 'vmtxel_data/pol', &
                               (/3,this%npol/), this%pol)

  call hdf5_write_double_array(file_id, 'vmtxel_data/kpt', &
                               (/3,this%nk/), this%kpt)

  !call hdf5_write_complex_array(file_id, 'vmtxel_data/dipole_flat', &
  !                             (/this%nmat,this%npol/), this%s1)

  call this%flat_to_band()

  call hdf5_write_scalar_array(file_id, 'vmtxel_data/dipole', &
                (/this%nband, this%mband, this%nk, this%ns, this%npol/), &
                this%rnm)

  if (this%has_velocity) then
    call hdf5_write_scalar_array(file_id, 'vmtxel_data/velocity', &
                  (/this%nband, this%mband, this%nk, this%ns, this%npol/), &
                  this%vnm)
  end if


  ! -------------------------------

  call hdf5_close_file(file_id)

  !call h5close_f(error)

  POP_SUB(write_vmtxel_arrays_hdf5)

end subroutine write_vmtxel_arrays_hdf5


#endif

! =========================================================================== !
! Reading routines
! =========================================================================== !


!> Read the file, using the format determined with use_hdf5
subroutine read_vmtxel(this)
  class(vmtxel_t), intent(inout) :: this

  PUSH_SUB(read_vmtxel)

#ifdef HDF5
  if (this%use_hdf5) then
    call this%read_vmtxel_hdf5()
  else
    call this%read_vmtxel_bin()
  end if
#else
    call this%read_vmtxel_bin()
#endif

  POP_SUB(read_vmtxel)

end subroutine read_vmtxel

! =========================================================================== !

!> Read binary file
subroutine read_vmtxel_bin(this)
  class(vmtxel_t), intent(inout) :: this

  integer :: ii,ipol
  integer :: ic,iv,ik,is
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)

  PUSH_SUB(read_vmtxel_bin)

  if (this%is_master) then
    write(6,'(1x,a)') 'Reading matrix elements from vmtxel'
    do ipol=1,this%npol
      if (this%npol==1) then
        fname = 'vmtxel'
      else
        fname = 'vmtxel_'//suffix(ipol)
      endif
      call open_file(16, file=trim(fname), form='unformatted', status='old')
      read(16) ik,ic,iv,is,ii
      if (ik.ne.this%nk.or.ic.ne.this%nband.or.iv.ne.this%mband &
        .or.is.ne.this%ns.or.ii.ne.this%opr) then
        write(0,'(a,5i6)') 'read  : ', ik,ic,iv,is,ii
        write(0,'(a,5i6)') 'needed: ', this%nk,this%nband,this%mband,this%ns,this%opr
        call die('parameter mismatch in vmtxel')
      endif
      read(16) this%s1(:,ipol)
      call close_file(16)
    enddo
  endif

#ifdef MPI
  call MPI_BCAST(this%s1,this%npol*this%nmat,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#endif

  POP_SUB(read_vmtxel_bin)

end subroutine read_vmtxel_bin

! =========================================================================== !

#ifdef HDF5

!> Read the vmtxel file in serial then broadcast the data
subroutine read_vmtxel_hdf5(this)
  class(vmtxel_t), intent(inout) :: this

  character(len=128) :: fname
  integer(HID_T) :: file_id

  PUSH_SUB(read_vmtxel_hdf5)

  fname = 'vmtxel.h5'

  call this%free()

  ! Read and broadcast dimensions
  call this%read_and_broadcast_vmtxel_header_hdf5(trim(fname))

  ! Allocate arrays
  call this%alloc(rnm=.true.)

  ! Read and broadcast data
  call this%read_and_broadcast_vmtxel_data_hdf5(trim(fname))

  POP_SUB(read_vmtxel_hdf5)

end subroutine read_vmtxel_hdf5

! =========================================================================== !

!> Read the dimensions then broadcast the result to everyone.
subroutine read_and_broadcast_vmtxel_header_hdf5(this, fname)
  class(vmtxel_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id

  PUSH_SUB(read_and_broadcast_vmtxel_header_hdf5)

  if (this%is_master) then
    call this%read_vmtxel_header_hdf5(fname)
  end if
  call this%broadcast_vmtxel_header()

  POP_SUB(read_and_broadcast_vmtxel_header_hdf5)

end subroutine read_and_broadcast_vmtxel_header_hdf5


!> Read the dimensions
subroutine read_vmtxel_header_hdf5(this, fname)
  class(vmtxel_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id

  PUSH_SUB(read_vmtxel_header_hdf5)

  !call h5open_f(error)

  call hdf5_open_file(fname, 'r', file_id)

  ! Read dimensions
  call hdf5_read_int(file_id, 'vmtxel_header/ns', this%ns)
  call hdf5_read_int(file_id, 'vmtxel_header/nk', this%nk)
  call hdf5_read_int(file_id, 'vmtxel_header/nband', this%nband)
  call hdf5_read_int(file_id, 'vmtxel_header/mband', this%mband)
  call hdf5_read_int(file_id, 'vmtxel_header/npol', this%npol)
  call hdf5_read_int(file_id, 'vmtxel_header/band_ordering', this%band_ordering)
  call hdf5_read_int(file_id, 'vmtxel_header/opr', this%opr)
  call hdf5_read_logical(file_id, 'vmtxel_header/has_velocity', this%has_velocity)
  call hdf5_close_file(file_id)

  !call h5close_f(error)

  POP_SUB(read_vmtxel_header_hdf5)

end subroutine read_vmtxel_header_hdf5


!> Broadcast the dimensions to all workers in the group.
subroutine broadcast_vmtxel_header(this, comm)

  class(vmtxel_t), intent(inout) :: this
  integer, intent(in), optional :: comm

  integer :: ndim=8
  integer :: dims(8)
  integer :: comm_

  integer :: has_velocity

  PUSH_SUB(broadcast_vmtxel_header)

#ifdef MPI
  if (present(comm)) then
    comm_ = comm
  else
    comm_ = MPI_COMM_WORLD
  end if
#else
  comm_ = 0
#endif

  if (this%is_master) then

     if (this%has_velocity) then
       has_velocity = 1
     else
       has_velocity = 0
     end if

     dims = (/this%mband, this%nband, this%nk, this%ns, this%npol, &
              this%band_ordering, this%opr, has_velocity/)
  end if

#ifdef MPI
  call MPI_BCAST(dims, ndim, MPI_INTEGER, 0, comm_, mpierr)
  call MPI_BARRIER(comm_, mpierr)
#endif

  if (.not. this%is_master) then
     this%mband = dims(1)
     this%nband = dims(2)
     this%nk = dims(3)
     this%ns = dims(4)
     this%npol = dims(5)
     this%band_ordering = dims(6)
     this%opr = dims(7)
     if (dims(8) .eq. 1) then
       this%has_velocity = .true.
     else
       this%has_velocity = .false.
     end if
  end if

  POP_SUB(broadcast_vmtxel_header)

end subroutine broadcast_vmtxel_header

! =========================================================================== !

!> Read the arrays in serial or parallel, and broadcast those read in serial.
!! At the moment, everything is read in serial then broadcast
subroutine read_and_broadcast_vmtxel_data_hdf5(this, fname)
  class(vmtxel_t), intent(inout) :: this
  character(len=*), intent(in) :: fname

  integer(HID_T) :: file_id       ! File identifier

  integer :: comm_

  PUSH_SUB(read_and_broadcast_vmtxel_data_hdf5)

#ifdef MPI
  comm_ = MPI_COMM_WORLD
#else
  comm_ = 0
#endif

  if (this%is_master) then

    call hdf5_open_file(fname, 'r', file_id)

    ! Read arrays
    !call hdf5_read_complex_array(file_id, 'vmtxel_data/dipole_flat', &
    !                             (/this%nmat,this%npol/), this%s1)

    call hdf5_read_double_array(file_id, 'vmtxel_data/kpt', &
                (/3, this%nk/), this%kpt)

    call hdf5_read_double_array(file_id, 'vmtxel_data/pol', &
                (/3, this%npol/), this%pol)

    call hdf5_read_scalar_array(file_id, 'vmtxel_data/dipole', &
                (/this%nband, this%mband, this%nk, this%ns, this%npol/), &
                this%rnm)

    call this%band_to_flat()

    if (this%has_velocity) then
      call hdf5_read_scalar_array(file_id, 'vmtxel_data/velocity', &
                  (/this%nband, this%mband, this%nk, this%ns, this%npol/), &
                  this%vnm)
    end if


    call this%band_to_flat()

    SAFE_DEALLOCATE(this%rnm)

    call hdf5_close_file(file_id)

  end if

  ! Broadcast
#ifdef MPI
  call MPI_BCAST(this%s1, this%npol*this%nmat, MPI_SCALAR, 0, comm_, mpierr)
  call MPI_BARRIER(comm_,mpierr)
  if (this%has_velocity) then
    call MPI_BCAST(this%vnm, this%npol*this%nmat, MPI_SCALAR, 0, comm_, mpierr)
  end if
#endif

  POP_SUB(read_and_broadcast_vmtxel_data_hdf5)

end subroutine read_and_broadcast_vmtxel_data_hdf5


#endif

end module vmtxel_m
