#include "f_defs.h"

module pseudopot_m

  use global_m
  use kpoint_pool_m,    only: kpoint_pool_t
  use bgw_mpi_m,        only: bgw_bcast

  implicit none

  private

  type, abstract :: pseudopot_t
    integer :: nat !< Number of atoms
    integer :: nsp !< Number of species
    integer :: nkb !< Number of KB projectors
    integer :: nhm !< Max. number of KB projectors per species
    integer :: nspin !< Number of spin components
    integer :: nspinor !< Number of spinor components
    integer, allocatable :: ityp(:) !< *PP species* for a given *atom index*
    integer, allocatable :: nh(:) !< number of projectors for a given *PP species*
    !> US D matrices (ikb, jkb, iatom, ispin).
    !! ikb,jkb run from 1..nh(ityp(iatom)).
    !! The matrix elements <G s|V_NL|Gp s> become:
    !! <G,up|V_NL|G,up>     = sum_{...} vkb(...,1) deeq(...,1) vkb(...,1)^H
    !! <G,down|V_NL|G,down> = sum_{...} vkb(...,2) deeq(...,2) vkb(...,2)^H
    real(DP), allocatable :: deeq(:,:,:,:)
    !> Non-collinear US D matrices (ikb, jkb, iatom, 1..4).
    !! ikb,jkb run from 1..nh(ityp(iatom)).
    !! The matrix elements <G s|V_NL|Gp sp> become:
    !! <G,1|V_NL|G,1> = sum_{...} vkb deeq_nc(...,1) vkb^H
    !! <G,1|V_NL|G,2> = sum_{...} vkb deeq_nc(...,2) vkb^H
    !! <G,2|V_NL|G,1> = sum_{...} vkb deeq_nc(...,3) vkb^H
    !! <G,2|V_NL|G,2> = sum_{...} vkb deeq_nc(...,4) vkb^H
    complex(DPC), allocatable :: deeq_nc(:,:,:,:)
    !> Potential times KB projectors (ng, nkb, nspin)
    complex(DPC), allocatable :: vkb(:,:,:)
    real(DP) :: num_elec = -1d0
    real(DP) :: ef = 0d0
  contains
    ! High-level routines
    procedure(pp_init), deferred :: init
    procedure(pp_free), deferred :: free
    procedure(pp_prepare_kpoints), deferred :: prepare_kpoints
    procedure(pp_get_nh_for_atom), deferred :: get_nh_for_atom
    procedure(pp_get_vkb_for_atom), deferred :: get_vkb_for_atom
    ! Aux. routines
    procedure, nopass :: newunit => pp_newunit
  endtype pseudopot_t


  abstract interface
    subroutine pp_init(this, fname, kpp, mf)
      import
      class(pseudopot_t), intent(out) :: this
      character(len=*), intent(in) :: fname
      type(kpoint_pool_t), intent(in) :: kpp
      type(mf_header_t), intent(in) :: mf
    end subroutine pp_init
  end interface


  abstract interface
    subroutine pp_prepare_kpoints(this, ik, kpp, mf, gind_k2m)
      import
      class(pseudopot_t), intent(inout) :: this
      integer, intent(in) :: ik
      type(kpoint_pool_t), intent(in) :: kpp
      type(mf_header_t), intent(in) :: mf
      integer, intent(inout) :: gind_k2m(:)
    end subroutine pp_prepare_kpoints
  end interface


  abstract interface
    subroutine pp_free(this)
      import
      class(pseudopot_t), intent(inout) :: this
    end subroutine pp_free
  end interface


  abstract interface
    integer function pp_get_nh_for_atom(this, iatom)
      import
      class(pseudopot_t), intent(in) :: this
      integer, intent(in) :: iatom
    end function pp_get_nh_for_atom
  end interface


  abstract interface
    subroutine pp_get_vkb_for_atom(this, iatom, ispin, vkb_atom, kpp)
      import
      class(pseudopot_t), intent(in) :: this
      integer, intent(in) :: iatom
      integer, intent(in) :: ispin
      complex(DPC), intent(out) :: vkb_atom(:,:) !<(ngkmax,nh)
      type(kpoint_pool_t), intent(in) :: kpp
    end subroutine pp_get_vkb_for_atom
  end interface


  !> Abstract type for pseudopotentials that explicitly have the KB in memory
  !! for each k-point iteration. Right now, all PPs inherit from this type.
  type, abstract, extends(pseudopot_t) :: pseudopot_explicit_t
    !> Number of atoms that I own
    integer :: nat_own
    !> Number of KB projectors that I own
    integer :: nkb_own
    integer, allocatable :: atom_offset(:)
  contains
    procedure :: get_nh_for_atom => pp2_get_nh_for_atom
    procedure :: get_vkb_for_atom => pp2_get_vkb_for_atom
  endtype pseudopot_explicit_t


  public :: pseudopot_t, pseudopot_explicit_t


contains

!-------------------------------------------------------------------------------
! pseudopot_t

integer function pp_newunit(unit)
  integer, intent(out), optional :: unit

  integer, parameter :: LUN_MIN=201, LUN_MAX=1000
  logical :: opened
  integer :: lun

  pp_newunit = -1
  do lun = LUN_MIN, LUN_MAX
    inquire(unit=lun, opened=opened)
    if (.not.opened) then
      pp_newunit = lun
      exit
    endif
  enddo
  if (present(unit)) unit = pp_newunit

end function pp_newunit



!-------------------------------------------------------------------------------
! pseudopot_explicit_t


integer function pp2_get_nh_for_atom(this, iatom)
  class(pseudopot_explicit_t), intent(in) :: this
  integer, intent(in) :: iatom

  PUSH_SUB(pp2_get_nh_for_atom)

  pp2_get_nh_for_atom = this%nh(this%ityp(iatom))

  POP_SUB(pp2_get_nh_for_atom)

end function pp2_get_nh_for_atom


subroutine pp2_get_vkb_for_atom(this, iatom, ispin, vkb_atom, kpp)
  class(pseudopot_explicit_t), intent(in) :: this
  integer, intent(in) :: iatom
  integer, intent(in) :: ispin
  complex(DPC), intent(out) :: vkb_atom(:,:) !<(ngkmax,nh)
  type(kpoint_pool_t), intent(in) :: kpp

  integer :: nh, ipe, iatom_loc, atom_offset, ngk

  PUSH_SUB(pp2_get_vkb_for_atom)

  ! FHJ: KB projectors are distributed over the atom index round robin
  ! among ranks in the k-point group.
  nh = this%get_nh_for_atom(iatom)
  ipe = INDXG2P(iatom, 1, kpp%inode, 0, kpp%npes)
  if (ipe==kpp%inode) then
    iatom_loc = INDXG2L(iatom, 1, kpp%inode, 0, kpp%npes)
    atom_offset = this%atom_offset(iatom_loc)
    ngk = size(this%vkb, 1)
    vkb_atom(:ngk,:nh) = this%vkb(:,atom_offset+1:atom_offset+nh,ispin)
  endif
  call timacc(44,1)
  call bgw_bcast(vkb_atom, ipe, kpp%comm)
  call timacc(44,2)

  POP_SUB(pp2_get_vkb_for_atom)

end subroutine pp2_get_vkb_for_atom

end module pseudopot_m
