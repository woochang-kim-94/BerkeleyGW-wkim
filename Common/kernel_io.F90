!>=========================================================================
!!
!! Module:
!!
!! (1) kernel_io_m     Originally by FHJ     Last Modified 10/07/2013 (FHJ)
!!
!!     Routines to read and write kernel files (bsedmat, ...)
!!     Inspired on wfn_rho_vxc_io.F90.
!!
!!=========================================================================

#include "f_defs.h"

module kernel_io_m

#ifdef HDF5
  use hdf5
  use hdf5_io_m
  use wfn_io_hdf5_m
#endif
  use global_m
  use wfn_rho_vxc_io_m

  implicit none

  private
 !> For library usage, do not make global_m contents available
 !! to avoid namespace clashes.

  public ::                        &
    read_binary_kernel_header,     &
    write_binary_kernel_header,    &
    read_format_kernel_header,     &
    write_format_kernel_header,    &
    read_kernel_header,            &
    write_kernel_header,           &
#ifdef HDF5
    read_kernel_header_hdf5,       &
    write_kernel_header_hdf5,      &
    setup_kernel_hdf5,             &
#endif
    xctinfo_to_kernel_header,      &
    check_xctinfo_kernel_header

contains

!===============================================================================
! FHJ: Below are the routines for kernel*mat files
!===============================================================================


#define FORMATTED
#define READ
!read_formatted_kernel_header
#include "kernel_io_inc.f90"
#undef READ
!write_formatted_kernel_header
#include "kernel_io_inc.f90"
#undef FORMATTED

#define BINARY
#define READ
!read_binary_kernel_header
#include "kernel_io_inc.f90"
#undef READ
!write_binary_kernel_header
#include "kernel_io_inc.f90"
#undef BINARY

#define READ
!read_kernel_header
#include "kernel_io_inc.f90"
#undef READ
!write_kernel_header
#include "kernel_io_inc.f90"

#ifdef HDF5
#define READ
!read_kernel_header_hdf5
#include "kernel_io_hdf5_inc.f90"
#undef READ
!write_kernel_header_hdf5
#include "kernel_io_hdf5_inc.f90"

! FHJ: Prepares the groups and empty datasets of a bsemat.h5 file.
subroutine setup_kernel_hdf5(fname, kern)
  character(len=*), intent(in) :: fname
  type(kernel_header_t), intent(in) :: kern !< kernel_header_t type

  integer(HID_T) :: file_id
  integer(HID_T) :: filespace
  integer :: dims(7)

  PUSH_SUB(setup_kernel_hdf5)

  call hdf5_create_file(fname, file_id)
  ! FHJ: Create empty groups
  call hdf5_create_group(file_id, 'bse_header')
  call hdf5_create_group(file_id, 'bse_header/params')
  call hdf5_create_group(file_id, 'bse_header/bands')
  call hdf5_create_group(file_id, 'bse_header/kpoints')
  call hdf5_create_group(file_id, 'mats')
  ! FHJ: Create empty dataset for the matrix elements
  dims(1) = kern%iflavor
  dims(2:3) = kern%n1b
  dims(4:5) = kern%n2b
  dims(6:7) = kern%nk*kern%ns
  call hdf5_create_dset(file_id, 'mats/head', H5T_NATIVE_DOUBLE, dims)
  if (kern%theory==0) then
    call hdf5_create_dset(file_id, 'mats/wing', H5T_NATIVE_DOUBLE, dims)
  endif
  call hdf5_create_dset(file_id, 'mats/body', H5T_NATIVE_DOUBLE, dims)
  call hdf5_create_dset(file_id, 'mats/exchange', H5T_NATIVE_DOUBLE, dims)
  if (kern%theory==1) then
    call hdf5_create_dset(file_id, 'mats/fxc', H5T_NATIVE_DOUBLE, dims)
  endif
  call hdf5_close_file(file_id)

  POP_SUB(setup_kernel_hdf5)

end subroutine setup_kernel_hdf5

#endif


  !> Populate the non-MF part of a kernel_header_t type. We assume that the MF
  !! part, i.e., kernel%mf, is already set (although we do overwrite
  !! kernel%mf%sheader='BSE' to play safe).
  !! Unlike WFN files, you can`t specify the flavor manually, it always matches SCALARSIZE.
  subroutine xctinfo_to_kernel_header(xct, kpts, kernel, nmat)
    type(xctinfo), intent(in) :: xct
    real(DP), intent(in) :: kpts(:,:) !< (3, nk)
    type(kernel_header_t), intent(inout) :: kernel
    integer, intent(in) :: nmat !< 1 for kernelxmat, 3 for kerneldmat

    PUSH_SUB(xctinfo_to_kernel_header)

    ! Generic header
    kernel%mf%sheader = 'KER'
    kernel%version = VER_BSE_HDF5
    kernel%iflavor = SCALARSIZE

    ! General paramters
    kernel%iscreen = xct%iscreen
    kernel%icutv = xct%icutv
    kernel%ecuts = xct%ecute
    kernel%ecutg = xct%ecutg
    kernel%efermi = xct%efermi
    kernel%theory = xct%theory
    kernel%nblocks = 1
    if (xct%extended_kernel) kernel%nblocks = 4
    kernel%storage = 0 ! Hard coded for now
    kernel%nmat = nmat
    kernel%energy_loss = xct%energy_loss

    ! K-point stuff
    kernel%nk = xct%nkpt_co
    SAFE_ALLOCATE(kernel%kpts, (3,kernel%nk))
    kernel%kpts = kpts(1:3, 1:kernel%nk)
    kernel%kgrid = kernel%mf%kp%kgrid
    kernel%qflag = xct%qflag
    kernel%exciton_Q_shift = xct%finiteq
    kernel%patched_sampling = xct%patched_sampling_co

    ! Bands stuff
    kernel%nvb = xct%nvb_co
    kernel%ncb = xct%ncb_co
    kernel%n1b = xct%n1b_co
    kernel%n2b = xct%n2b_co
    kernel%ns = xct%nspin
    kernel%nspinor = xct%nspinor

    POP_SUB(xctinfo_to_kernel_header)

  end subroutine xctinfo_to_kernel_header

  subroutine check_xctinfo_kernel_header(fname, xct, kernel)
    character(len=*), intent(in) :: fname !< file name
    type(xctinfo), intent(in) :: xct
    type(kernel_header_t), intent(in) :: kernel

    integer :: nblocks

    PUSH_SUB(check_xctinfo_kernel_header)

    ! Generic variables
    ! FHJ: absorption doesn`t care about WFN cutoff (for now!)
    !call check_R('WFN cutoff', xct%ecutg, kernel%mf%kp%ecutwfc)
    call check_I('screening flag', xct%iscreen, kernel%iscreen)
    call check_I('truncation flag', xct%icutv, kernel%icutv)
    ! FHJ: absorption doesn`t care about epsilon cutoff
    !call check_R('epsilon cutoff', xct%ecute, kernel%ecuts)

    ! Specific to bsemat
    call check_I('# of k-point', xct%nkpt_co, kernel%nk)

    call check_I('# of spins', xct%nspin, kernel%ns)
    call check_I('# of spinor components', xct%nspinor, kernel%nspinor)
    call check_I('# of val. bands', xct%nvb_co, kernel%nvb)
    call check_I('# of cond. bands', xct%ncb_co, kernel%ncb)
    call check_I('# of bands 1', xct%n1b_co, kernel%n1b)
    call check_I('# of bands 2', xct%n2b_co, kernel%n2b)

    call check_I('theory level', 0, kernel%theory) ! Hard coded for now
    nblocks = 1
    if (xct%extended_kernel) nblocks = 4
    call check_I('number of transition blocks', nblocks, kernel%nblocks)
    call check_I('storage format', 0, kernel%storage) ! Hard coded for now

    POP_SUB(check_xctinfo_kernel_header)

    contains

      subroutine check_I(label, ref, got)
        character(len=*), intent(in) :: label
        integer, intent(in) :: ref
        integer, intent(in) :: got

        PUSH_SUB(check_xctinfo_kernel_header.check_I)

        if (ref/=got) then
          if (peinf%inode==0) then
            write(0,'(1x,3a)') 'ERROR: incompatible values found in file "', fname, '".'
            write(0,'(1x,2a)') 'Quantity: ', label
            write(0,'(1x,a,i0)') 'Expected: ', ref
            write(0,'(1x,a,i0)') 'Got: ', got
          endif
          call die('Incompatible values in file "'+fname+'".', &
            only_root_writes=.true.)
        endif

        POP_SUB(check_xctinfo_kernel_header.check_I)

      end subroutine check_I

      subroutine check_R(label, ref, got)
        character(len=*), intent(in) :: label
        real(DP), intent(in) :: ref
        real(DP), intent(in) :: got

        real(DP) :: tol=TOL_Small

        PUSH_SUB(check_xctinfo_kernel_header.check_R)

        if (dabs(ref-got)>tol) then
          if (peinf%inode==0) then
            write(0,'(1x,3a)') 'ERROR: incompatible values found in file "', fname, '".'
            write(0,'(1x,2a)') 'Quantity: ', label
            write(0,'(1x,a,f0.8)') 'Expected: ', ref
            write(0,'(1x,a,f0.8)') 'Got: ', got
          endif
          call die('Incompatible values in file "'+fname+'".', &
            only_root_writes=.true.)
        endif

        POP_SUB(check_xctinfo_kernel_header.check_R)

      end subroutine check_R

  end subroutine check_xctinfo_kernel_header

end module kernel_io_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
