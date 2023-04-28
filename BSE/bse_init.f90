!==============================================================================
!
! Module bse_init_m
!
! (1) bse_init()         Originally By FHJ       Last Modified  3/24/2012 (FHJ)
!
!     Performs initializations required by absorption and inteqp calculations.
!
!==============================================================================

#include "f_defs.h"

module bse_init_m

  use global_m
#ifdef HDF5
  use hdf5
#endif
  use fullbz_m
  use intpts_m
  use wfn_rho_vxc_io_m
  use wfn_io_hdf5_m

  implicit none

  private

  public :: bse_init

contains

!> FHJ: Figures out the number of k-points in the coarse grid and the
!!      dimensionality of the system for interpolation pourposes.
!! FIXME: am I freeing all buffers correctly?
subroutine bse_init(xct,flag)
  type (xctinfo), intent(inout) :: xct
  type (flags), intent(in) :: flag

  character(len=3) :: sheader
  integer :: iflavor
  type (crystal) :: crys, crys_co
  type (gspace) :: gvec, gvec_co
  type (grid) :: kg, kg_co
  type (kpoints) :: kp, kp_co
  type (symmetry) :: syms, syms_co
  logical :: skip_checkbz

  logical :: is_periodic_old(3)
  integer :: jdim, npts_intp_kernel

  PUSH_SUB(bse_init)

  ! GKA: No additional information is needed if the eigenvalues are already computed
  if (flag%spec.eq.1) then
    return
  endif

  if (flag%read_dtmat) then
    call open_file(unit=13,file='dtmat',form='unformatted',status='old')
    read(13) xct%idimensions, xct%is_periodic(1:3), npts_intp_kernel, xct%nkpt_co
    if (peinf%inode==0) then
      if (xct%npts_intp_kernel==-1 .and. npts_intp_kernel==1) then
        write(0,*)
        write(0,'(a)') 'WARNING: dtmat was not calculated with "kernel_k_interpolation",'
        write(0,'(a)') 'so we are turning off "kernel_k_interpolation".'
        write(0,*)
      elseif (xct%npts_intp_kernel==1 .and. npts_intp_kernel>1) then
        write(0,*)
        write(0,'(a)') 'WARNING: dtmat was calculated with "kernel_k_interpolation",'
        write(0,'(a)') 'so we are turning on "kernel_k_interpolation".'
        write(0,*)
      endif
    endif
    xct%npts_intp_kernel = npts_intp_kernel
    call close_file(13)
    return
  endif

  if (peinf%inode==0) write(6,*)

  ! FHJ: Read the header of WFN_co
  sheader = 'WFN'
  iflavor = 0
  if ( xct%use_wfn_hdf5 ) then
#ifdef HDF5
   call read_hdf5_header_type('WFN_co.h5', sheader, iflavor, kp_co, gvec_co, syms_co, crys_co)
#endif
  else
    if (peinf%inode == 0) &
      call open_file(unit=25, file='WFN_co', form='unformatted', status='old')
    call read_binary_header_type(25, sheader, iflavor, kp_co, gvec_co, &
      syms_co, crys_co, warn=.false., dont_warn_kgrid=.true.)
    if (peinf%inode == 0) call close_file(25)
  end if
  kg_co%nr = kp_co%nrk
  SAFE_ALLOCATE(kg_co%r, (3,kg_co%nr))
  kg_co%r(:,:) = kp_co%rk(:,:)

  ! FHJ: We just need to unfold the k-points fast to get the number of k-points
  !      in each direction, so we won`t build the wigner seitz cell.
  if (flag%bzc==1) then
    call fullbz(crys_co, syms_co, kg_co, 1, skip_checkbz, &
                wigner_seitz=.false., paranoid=.false., do_nothing=.true.)
  else
    call fullbz(crys_co, syms_co, kg_co, syms_co%ntran, skip_checkbz, &
                wigner_seitz=.false., paranoid=.false.)
  endif

  xct%nkpt_co = kg_co%nf
  if (peinf%inode==0) then
    write(6,'(1x,a,i0)') 'Number of k-points in the coarse k-grid: ', xct%nkpt_co
  endif

  ! FHJ: find the number of periodic dimensions from the kpt sampling.
  !      Don`t trust kp%kgrid!
  call get_ndims(kg_co, xct)

  SAFE_DEALLOCATE_P(kg_co%r)
  SAFE_DEALLOCATE_P(kg_co%f)


  if (xct%skipinterp) then
    xct%npts_intp_kernel = 1
  else
    ! FHJ: Read the header of WFN, if we are interpolating
    sheader = 'WFN'
    iflavor = 0
    if ( xct%use_wfn_hdf5 ) then
#ifdef HDF5
      call read_hdf5_header_type('WFN_fi.h5', sheader, iflavor, kp, gvec, syms, crys)
#endif
    else
      if (peinf%inode == 0) &
        call open_file(25, file='WFN_fi', form='unformatted', status='old')
      call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys, &
        dont_warn_kgrid=.true.)
      if (peinf%inode == 0) call close_file(25)
    end if
    kg%nr = kp%nrk
    SAFE_ALLOCATE(kg%r, (3,kg%nr))
    kg%r(:,:) = kp%rk(:,:)

    if (flag%bz0==0.and.xct%is_absorption) then
      call fullbz(crys, syms, kg, syms%ntran, skip_checkbz, &
                  wigner_seitz=.false., paranoid=.false.)
    else
      call fullbz(crys, syms, kg, 1, skip_checkbz, &
                  wigner_seitz=.false., paranoid=.false., do_nothing=.true.)
    endif

    is_periodic_old = xct%is_periodic
    call get_ndims(kg, xct)
    xct%is_periodic = xct%is_periodic .or. is_periodic_old
    xct%idimensions = 0
    do jdim = 1, 3
      if (xct%is_periodic(jdim)) xct%idimensions = xct%idimensions + 1
    enddo

    SAFE_DEALLOCATE_P(kg%r)
    SAFE_DEALLOCATE_P(kg%f)

    if (xct%npts_intp_kernel==-1) xct%npts_intp_kernel = xct%idimensions + 1

    if (peinf%inode==0) then
      write(6,'(1x,a,i1,a)') 'A ',xct%idimensions,'-D interpolation algorithm will be employed.'
    endif
  endif ! xct%skipinterp

  POP_SUB(bse_init)

end subroutine bse_init

end module bse_init_m

