#include "f_defs.h"

module pseudopot_vkb_m

  use global_m
  use bgw_mpi_m,        only: bgw_bcast
  use pseudopot_m,      only: pseudopot_explicit_t
  use kpoint_pool_m,    only: kpoint_pool_t
  use wfn_rho_vxc_io_m, only: read_binary_complex_data, read_binary_gvectors, check_header
  use misc_m,           only: findvector

  implicit none

  private

  type, extends(pseudopot_explicit_t) :: pseudopot_vkb_t
    integer :: nk
    integer :: ngkmax
    !> (pseudopot%nat) Maps a KB center to an atom. We use this mapping array only
    !! when we read the projectors from the file, becuase QE has a funny order.
    integer, allocatable :: center2at(:)
    !
    integer :: iunit
  contains
    ! High-level routines
    procedure :: init => pp_init
    procedure :: free => pp_free
    procedure :: prepare_kpoints => pp_prepare_kpoints
    ! Low-level routines
    procedure :: read_header => pp_read_header
    procedure :: read_beta => pp_read_beta
  endtype pseudopot_vkb_t

  public :: pseudopot_vkb_t, compare_mf_headers

contains


!--------------------------------------------------------------------------------------
! High-level routines
!--------------------------------------------------------------------------------------


subroutine pp_init(this, fname, kpp, mf)
  class(pseudopot_vkb_t), intent(out) :: this
  character(len=*), intent(in) :: fname
  type(kpoint_pool_t), intent(in) :: kpp
  type(mf_header_t), intent(in) :: mf

  PUSH_SUB(pp_init)

  call timacc(23,1)
  this%iunit = this%newunit()
  if (peinf%inode==0) then
    call open_file(this%iunit, file=fname, status='old', form='unformatted')
    write(6,'(1x,a)') 'Reading k-point-dependent Kleinman-Bylander projectors from file ' // &
      trim(fname)
  endif
  call this%read_header(this%iunit, kpp, mf)
  call timacc(23,2)

  POP_SUB(pp_init)

end subroutine pp_init


subroutine pp_prepare_kpoints(this, ik, kpp, mf, gind_k2m)
  class(pseudopot_vkb_t), intent(inout) :: this
  integer, intent(in) :: ik
  type(kpoint_pool_t), intent(in) :: kpp
  type(mf_header_t), intent(in) :: mf
  integer, intent(inout) :: gind_k2m(:)

  integer :: ipool, ik_first, ik2
  complex(DPC), allocatable :: vkb_buf(:,:,:)

  PUSH_SUB(pp_prepare_kpoints)

  ik_first = ik - kpp%ipool
  call logit('Reading VKB')
  call timacc(24,1)
  ! Get KB projectors for k-points. Everyone has to participate.
  ! For each iteration, loop over all pools and figure out the k-point for
  ! that pool. We can exit if we went past the max number of k-points
  ! (i.e., we reached an idle pool)
  SAFE_ALLOCATE(vkb_buf, (this%ngkmax,this%nkb_own,this%nspin))
  do ipool = 0, kpp%npools-1
    ik2 = ik_first + ipool
    if (ik2>mf%kp%nrk) exit
    call this%read_beta(this%iunit, kpp, mf, ik2, vkb_buf, gind_k2m)
    if (ipool == kpp%ipool) this%vkb(:,:,:) = vkb_buf(:,:,:)
  enddo
  SAFE_DEALLOCATE(vkb_buf)
  call timacc(24,2)
  call logit('Done VKB')

  POP_SUB(pp_prepare_kpoints)

end subroutine pp_prepare_kpoints


subroutine pp_free(this)
  class(pseudopot_vkb_t), intent(inout) :: this

  PUSH_SUB(pp_free)

  if (peinf%inode==0) call close_file(this%iunit)
  SAFE_DEALLOCATE(this%ityp)
  SAFE_DEALLOCATE(this%nh)
  if (this%nspinor==1) then
    SAFE_DEALLOCATE(this%deeq)
  else
    SAFE_DEALLOCATE(this%deeq_nc)
  endif
  SAFE_DEALLOCATE(this%center2at)
  SAFE_DEALLOCATE(this%vkb)

  POP_SUB(pp_free)

end subroutine pp_free


!--------------------------------------------------------------------------------------
! Low-level routines
!--------------------------------------------------------------------------------------


subroutine pp_read_header(this, iunit, kpp, mf, allocate_vkb_buf)
  class(pseudopot_vkb_t), intent(inout) :: this
  integer, intent(in) :: iunit
  type(kpoint_pool_t), intent(in) :: kpp
  type(mf_header_t), intent(in) :: mf
  logical, intent(in), optional :: allocate_vkb_buf ! Default is true

  type(mf_header_t) :: mf2
  character(len=32) :: stitle
  integer :: ii, nh, iel, iat, iat_loc, icenter
  logical :: allocate_vkb_buf_

  PUSH_SUB(pp_read_header)

  if (peinf%inode==0) then
    read(iunit) stitle, mf2%sdate, mf2%stime
    read(iunit) mf2%kp%nspin, mf2%gvec%ng, mf2%syms%ntran, mf2%syms%cell_symmetry, &
      mf2%crys%nat, mf2%gvec%ecutrho, mf2%kp%nrk, &
      this%nsp, this%nkb, this%nhm, mf2%kp%ngkmax, mf2%kp%ecutwfc
    read(iunit) mf2%gvec%FFTgrid, mf2%kp%kgrid, mf2%kp%shift
    read(iunit) mf2%crys%celvol, mf2%crys%alat, mf2%crys%avec, mf2%crys%adot
    read(iunit) mf2%crys%recvol, mf2%crys%blat, mf2%crys%bvec, mf2%crys%bdot
    read(iunit) mf2%syms%mtrx(1:3,1:3,1:mf2%syms%ntran)
    read(iunit) mf2%syms%tnp(1:3,1:mf2%syms%ntran)

    select case (mf2%kp%nspin)
      case (4)
        mf2%kp%nspin = 1
        mf2%kp%nspinor = 2
      case (1,2)
        mf2%kp%nspinor = 1
      case default
        call die('Invalid number of spins when reading '//stitle)
    endselect
    this%nspin = mf2%kp%nspin
    this%nspinor = mf2%kp%nspinor

    this%nat = mf2%crys%nat
    this%nk = mf2%kp%nrk
    this%ngkmax = mf2%kp%ngkmax

    SAFE_ALLOCATE(mf2%crys%atyp, (this%nat))
    SAFE_ALLOCATE(mf2%crys%apos, (3,this%nat))
    SAFE_ALLOCATE(mf2%kp%ngk, (this%nk))
    SAFE_ALLOCATE(mf2%kp%w, (this%nk))
    SAFE_ALLOCATE(mf2%kp%rk, (3,this%nk))

    read(iunit) (mf2%crys%apos(1:3,ii), mf2%crys%atyp(ii), ii=1,mf2%crys%nat)
    read(iunit) mf2%kp%ngk(1:this%nk)
    read(iunit) mf2%kp%w(1:this%nk)
    read(iunit) mf2%kp%rk(1:3,1:this%nk)
  endif ! peinf%inode==0

  call bgw_bcast(this%nkb)
  call bgw_bcast(this%nhm)
  call bgw_bcast(this%nsp)
  call bgw_bcast(this%nspin)
  call bgw_bcast(this%nspinor)
  call bgw_bcast(this%nat)
  call bgw_bcast(this%nk)
  call bgw_bcast(this%ngkmax)

  SAFE_ALLOCATE(this%ityp, (this%nat))
  SAFE_ALLOCATE(this%nh, (this%nsp))
  if (this%nspinor==1) then
    SAFE_ALLOCATE(this%deeq, (this%nhm,this%nhm,this%nat,this%nspin))
  else
    SAFE_ALLOCATE(this%deeq_nc, (this%nhm,this%nhm,this%nat,4))
  endif
  if (peinf%inode==0) then
    read(iunit) this%ityp(:)
    read(iunit) this%nh(:)
    if (this%nspinor==1) then
      read(iunit) this%deeq(:,:,:,:)
    else
      read(iunit) this%deeq_nc(:,:,:,:)
    endif
    SAFE_ALLOCATE(mf2%gvec%components, (3, mf2%gvec%ng))
    call read_binary_gvectors(iunit, mf2%gvec%ng, mf2%gvec%ng, mf2%gvec%components, bcast=.false.)

    if (stitle(1:11) /= 'VKB-Complex') call die('Wrong header')
    call compare_mf_headers('WFN', mf, 'VKB', mf2, .true.)

    SAFE_DEALLOCATE_P(mf2%crys%atyp)
    SAFE_DEALLOCATE_P(mf2%crys%apos)
    SAFE_DEALLOCATE(mf2%kp%ngk)
    SAFE_DEALLOCATE(mf2%kp%w)
    SAFE_DEALLOCATE(mf2%kp%rk)
    SAFE_DEALLOCATE_P(mf2%gvec%components)
  endif
  call bgw_bcast(this%ityp)
  call bgw_bcast(this%nh)
  if (this%nspinor==1) then
    call bgw_bcast(this%deeq)
  else
    call bgw_bcast(this%deeq_nc)
  endif

  ! FHJ: QE orders the KB projectors in a funny way. Instead of looping over atoms
  ! with the same order as the atomic definitions, it reorders the atoms and
  ! groups them by species. There is no computation advantage of doing so here,
  ! so we instead reoder these array from QE to the "natural" order.
  ! We use the mapping array, this%center2at, only when we read the KB projectors
  ! from file. After that, we stick with a consistent "natural" order.
  SAFE_ALLOCATE(this%center2at, (this%nat))
  icenter = 0
  do iel = 1, this%nsp
    do iat = 1, this%nat
      if (this%ityp(iat)==iel) then
        icenter = icenter + 1
        this%center2at(icenter) = iat
      endif
    enddo
  enddo

  ! FHJ: distribute atoms
  this%nat_own = NUMROC(this%nat, 1, kpp%inode, 0, kpp%npes)
  SAFE_ALLOCATE(this%atom_offset, (max(1,this%nat_own)))
  this%atom_offset(:) = 0
  this%nkb_own = 1
  do iat_loc = 1, this%nat_own
    iat = INDXL2G(iat_loc, 1, kpp%inode, 0, kpp%npes)
    iel = this%ityp(iat)
    nh = this%nh(iel)
    if (iat_loc<this%nat_own) then
      this%atom_offset(iat_loc+1) = this%atom_offset(iat_loc) + nh
    else
      this%nkb_own = this%atom_offset(iat_loc) + nh
    endif
  enddo

  allocate_vkb_buf_ = .true.
  if (present(allocate_vkb_buf)) allocate_vkb_buf_ = allocate_vkb_buf
  if (allocate_vkb_buf_) then
    SAFE_ALLOCATE(this%vkb, (this%ngkmax,this%nkb_own,this%nspin))
    this%vkb = (0d0, 0d0)
  endif

  POP_SUB(pp_read_header)

end subroutine pp_read_header


subroutine pp_read_beta(this, iunit, kpp, mf, ik, buf, gind_k2m)
  class(pseudopot_vkb_t), intent(in) :: this
  integer, intent(in) :: iunit
  type(kpoint_pool_t), intent(in) :: kpp
  type(mf_header_t), intent(in) :: mf
  integer, intent(in) :: ik
  complex(DPC), intent(out) :: buf(this%ngkmax, this%nkb_own, this%nspin)
  integer, intent(out) :: gind_k2m(:)

  complex(DPC) :: tmp_buf(this%ngkmax, this%nspin)

  integer :: gvectmp(3, this%ngkmax)
  integer :: is, ikb, iel, iat, ih, nh, ipe, iat_loc, icenter, atom_offset, ig, iout

  PUSH_SUB(pp_read_beta)

  buf = (0d0, 0d0)
  call read_binary_gvectors(iunit, mf%kp%ngk(ik), this%ngkmax, gvectmp)
  do ig = 1, mf%kp%ngk(ik)
    call findvector(iout, gvectmp(:,ig), mf%gvec)
    gind_k2m(ig) = iout
  enddo
  if (any(gind_k2m(1:mf%kp%ngk(ik))==0)) then
    if (kpp%inode==0) write(0,*) ' ik=',ik
    call die('failed to find G-vector for k-point')
  endif

  !FHJ: TODO - add consistency check between vsc/vkb gvecs and wfn isort
  do is = 1, this%nspin
    ikb = 0
    do icenter = 1, this%nat
      iat = this%center2at(icenter)
      iel = this%ityp(iat)
      nh = this%nh(iel)
      ipe = INDXG2P(iat, 1, kpp%inode, 0, kpp%npes)
      iat_loc = INDXG2L(iat, 1, kpp%inode, 0, kpp%npes)
      if (ipe==kpp%inode) atom_offset = this%atom_offset(iat_loc)
      do ih = 1, nh
        ikb = ikb + 1
        call read_binary_complex_data(iunit, mf%kp%ngk(ik), this%ngkmax, 1, &
          tmp_buf(:,is:is))
        if (ipe==kpp%inode) then
          buf(:,atom_offset+ih,is) = tmp_buf(:,is)
        endif
      enddo
    enddo ! ikb
  enddo ! is

  POP_SUB(pp_read_beta)

end subroutine pp_read_beta


!--------------------------------------------------------------------------------------
! Auxiliary routines
!--------------------------------------------------------------------------------------


subroutine compare_mf_headers(name1, mf1, name2, mf2, is_wfn)
  character(len=*) :: name1
  type(mf_header_t), intent(in) :: mf1
  character(len=*) :: name2
  type(mf_header_t), intent(in) :: mf2
  logical, intent(in) :: is_wfn

  character(len=256) :: string

  PUSH_SUB(compare_mf_headers)

  string = TRUNC(name1) // " vs. " // TRUNC(name2) // ":"
  call check_header(name1, mf1%kp, mf1%gvec, mf1%syms, mf1%crys, &
    name2, mf2%kp, mf2%gvec, mf2%syms, mf2%crys, is_wfn, .false.)

  if (is_wfn) then
    if (mf1%kp%nrk /= mf2%kp%nrk) call die(TRUNC(string)//&
      " number of k-points mismatch")
    if (any(abs(mf1%kp%rk - mf2%kp%rk) > TOL_SMALL)) call die(TRUNC(string)//&
      " k-points mismatch")
    if (any(abs(mf1%kp%w - mf2%kp%w) > TOL_SMALL)) call die(TRUNC(string)//&
      " k-weights mismatch")
    if (any(mf1%kp%kgrid /= mf2%kp%kgrid)) call die(TRUNC(string)//&
      " k-grid mismatch")
    if (any(abs(mf1%kp%shift - mf2%kp%shift) > TOL_SMALL)) call die(TRUNC(string)//&
      " k-shift mismatch")
    if (any(mf1%kp%ngk /= mf2%kp%ngk)) call die(TRUNC(string)//&
      " number of G-vectors per k-point mismatch")
    if (mf1%kp%ngkmax /= mf2%kp%ngkmax) call die(TRUNC(string)//&
      " max. number of G-vectors per k-point mismatch")
  endif
  if (any(mf1%gvec%components /= mf2%gvec%components)) call die(TRUNC(string)//&
    " G-space mismatch")

  POP_SUB(compare_mf_headers)

end subroutine compare_mf_headers


end module pseudopot_vkb_m
