!===============================================================================
!
! Program:
!
! (1) paraband            Originally By FHJ      Last Modified Apr/2015 (FHJ)
!
! Input is read from file split_spin.inp.
!
!===============================================================================

#include "f_defs.h"

program split_spin
  use global_m
  use inread_m,         only: pb_params_t, inread
  use wfn_rho_vxc_io_m, only: read_mf_header, read_binary_gvectors, read_binary_complex_data, &
                              write_mf_header, write_binary_gvectors, write_binary_complex_data
  use wfn_io_m,         only: wfn_read
  use io_utils_m,       only: progress_info, progress_init, progress_step, &
                              progress_free, print_dealing_with
  use pseudopot_vkb_m,  only: pseudopot_vkb_t, compare_mf_headers

  implicit none

  type(mf_header_t) :: mf_wfn, mf_vsc, mf_nospin
  complex(DPC), allocatable :: vsc(:,:), vsc_nospin(:,:), wfn_buf(:,:)
  integer, allocatable :: isort(:,:)
  integer, allocatable :: gveck_buf(:,:), gvectmp(:,:)
  complex(DPC), allocatable :: vkb_buf(:,:)
  type(pseudopot_vkb_t) :: pseudopot, pseudopot_nospin
  character(len=128) :: fname, solver_name
  character(len=32) :: stitle
  type(pb_params_t) :: params
  type(progress_info) :: prog_info
  integer :: is, ikb, ik, ii, ib


  !call peinfo_init()
  write(6,'(1x,a, a)') 'Git version: ', TOSTRING(GIT_COMMIT)
  if (peinf%npes/=1) then
    call die('You run this application with a single processor', &
    only_root_writes=.true.)
  endif

!------------------------
! Read input parameters
  write(6,'(1x,a)') 'Reading parameters from file parabands.inp'
  call inread(params)


!-----------------------
! Split WFN
  write(*,'(/1x,a)') 'Reading WFN header'
  call open_file(7, file=params%fname_wfn_in, status='old', form='unformatted')
  call read_mf_header(7, mf_wfn, iflavor=SCALARSIZE, sheader='WFN', warn=.false., dont_warn_kgrid=.true.)
  if (mf_wfn%kp%nspin/=2) then
    call die('Input WFN file is not spin polarized.')
  endif

  SAFE_ALLOCATE(mf_wfn%gvec%components, (3, mf_wfn%gvec%ng))
  call read_binary_gvectors(7, mf_wfn%gvec%ng, mf_wfn%gvec%ng, mf_wfn%gvec%components)
  SAFE_ALLOCATE(gvectmp, (3,mf_wfn%kp%ngkmax))
  SAFE_ALLOCATE(wfn_buf, (mf_wfn%kp%ngkmax,mf_wfn%kp%nspin))

  mf_nospin = mf_wfn
  associate(kp=>mf_nospin%kp)
    kp%nspin = 1
    ! By default, assignment of derived types assigns the pointers which are
    ! members of the types. But if two pointers point to the same memory
    ! and one of them gets deallocated, all pointers get deallocated.
    ! FIXME: we changed these to be allocatable now. Is this deprecated?
    SAFE_ALLOCATE(kp%ifmin, (kp%nrk,kp%nspin))
    SAFE_ALLOCATE(kp%ifmax, (kp%nrk,kp%nspin))
    SAFE_ALLOCATE(kp%el, (kp%mnband,kp%nrk,kp%nspin))
    SAFE_ALLOCATE(kp%occ, (kp%mnband,kp%nrk,kp%nspin))
  endassociate

  do is = 1, 2
    associate(kp=>mf_nospin%kp, kp_wfn=>mf_wfn%kp)
      kp%ifmin(:,1) = kp_wfn%ifmin(:,is)
      kp%ifmax(:,1) = kp_wfn%ifmax(:,is)
      kp%el(:,:,1) = kp_wfn%el(:,:,is)
      kp%occ(:,:,1) = kp_wfn%occ(:,:,is)
    endassociate
    write(fname,'(a,i0)') 'WFN_spin_', is
    call open_file(10+is, file=trim(fname), status='replace', form='unformatted')
    call write_mf_header(10+is, mf_nospin)
    call write_binary_gvectors(10+is, mf_nospin%gvec%ng, mf_nospin%gvec%ng, &
      mf_nospin%gvec%components)
  enddo

  call progress_init(prog_info, 'copying WFN', 'state', mf_wfn%kp%nrk*mf_wfn%kp%mnband*2)
  do ik = 1, mf_wfn%kp%nrk
    call read_binary_gvectors(7, mf_wfn%kp%ngk(ik), mf_wfn%kp%ngkmax, gvectmp)
    do is = 1, 2
      call write_binary_gvectors(10+is, mf_nospin%kp%ngk(ik), mf_nospin%kp%ngkmax, gvectmp)
    enddo
    do ib = 1, mf_wfn%kp%mnband
      call read_binary_complex_data(7, mf_wfn%kp%ngk(ik), mf_wfn%kp%ngkmax, mf_wfn%kp%nspin, &
        wfn_buf)
      do is = 1, 2
         call progress_step(prog_info)
        call write_binary_complex_data(10+is, mf_nospin%kp%ngk(ik), mf_nospin%kp%ngkmax, mf_nospin%kp%nspin, &
          wfn_buf(:,is:is))
      enddo
    enddo
  enddo
  call close_file(7)
  call close_file(11)
  call close_file(12)
  SAFE_DEALLOCATE(gvectmp)
  SAFE_DEALLOCATE(wfn_buf)
  call progress_free(prog_info)


!------------------------
! Split self-consistent potential file
  write(6,'(/1x,a)') 'Reading self-consistent potential from file '//TRUNC(params%fname_vsc)
  call open_file(10, file=params%fname_vsc, status='old', form='unformatted')
  call read_mf_header(10, mf_vsc, iflavor=SCALARSIZE, sheader='VSC', warn=.false., dont_warn_kgrid=.true.)
  if (mf_vsc%kp%nspin/=2) then
    call die('Input VSC file is not spin polarized.')
  endif
  SAFE_ALLOCATE(mf_vsc%gvec%components, (3, mf_vsc%gvec%ng))
  call read_binary_gvectors(10, mf_vsc%gvec%ng, mf_vsc%gvec%ng, &
    mf_vsc%gvec%components, bcast=.false.)
  call compare_mf_headers('WFN', mf_wfn, 'VSC', mf_vsc, .false.)
  SAFE_ALLOCATE(vsc, (mf_vsc%gvec%ng, mf_vsc%kp%nspin))
  call read_binary_complex_data(10, mf_vsc%gvec%ng, mf_vsc%gvec%ng, mf_vsc%kp%nspin, vsc)
  call close_file(10)

  mf_nospin = mf_vsc
  mf_nospin%kp%nspin = 1
  SAFE_ALLOCATE(vsc_nospin, (mf_nospin%gvec%ng, mf_nospin%kp%nspin))
  do is = 1, 2
    write(*,'(1x,a,i0)') 'Writing spin-unpolarized VSC for spin ', is
    vsc_nospin(:,1) = vsc(:,is)
    write(fname,'(a,i0)') 'VSC_spin_', is
    call open_file(10, file=trim(fname), status='replace', form='unformatted')
    call write_mf_header(10, mf_nospin)
    call write_binary_gvectors(10, mf_nospin%gvec%ng, mf_nospin%gvec%ng, &
      mf_nospin%gvec%components)
    call write_binary_complex_data(10, mf_nospin%gvec%ng, mf_nospin%gvec%ng, mf_nospin%kp%nspin, vsc_nospin)
    call close_file(10)
  enddo
  SAFE_DEALLOCATE(vsc_nospin)


!------------------------
! Split Kleinman-Bylander projectors file
  if (params%has_vkb) then
    call open_file(7, file=params%fname_vkb, status='old', form='unformatted')
    write(6,'(/1x,a)') 'Reading Kleinman-Bylander projectors from file ' // &
      TRUNC(params%fname_vkb)
    call pseudopot%read_header(7, params%kpp, mf_wfn, allocate_vkb_buf=.false.)

    stitle = 'VKB-Complex'
    pseudopot_nospin = pseudopot
    pseudopot_nospin%nspin = 1
    SAFE_ALLOCATE(vkb_buf, (pseudopot%ngkmax,1))
    SAFE_ALLOCATE(gveck_buf, (3,pseudopot%ngkmax))
    do is = 1, pseudopot%nspin
      ! Write from now on. We should only have pseudopot_nospin%.. below, and no pseudopot%..
      !write(*,'(1x,a,i0)') 'Writing spin-unpolarized VKB for spin ', is

      write(fname,'(a,i0)') 'VKB_spin_', is
      call open_file(10+is, file=trim(fname), status='replace', form='unformatted')
      write(10+is) stitle, mf_wfn%sdate, mf_wfn%stime
      write(10+is) 1, mf_wfn%gvec%ng, mf_wfn%syms%ntran, mf_wfn%syms%cell_symmetry, &
        mf_wfn%crys%nat, mf_wfn%gvec%ecutrho, mf_wfn%kp%nrk, &
        pseudopot_nospin%nsp, pseudopot_nospin%nkb, pseudopot_nospin%nhm, mf_wfn%kp%ngkmax, mf_wfn%kp%ecutwfc
      write(10+is) mf_wfn%gvec%FFTgrid, mf_wfn%kp%kgrid, mf_wfn%kp%shift
      write(10+is) mf_wfn%crys%celvol, mf_wfn%crys%alat, mf_wfn%crys%avec, mf_wfn%crys%adot
      write(10+is) mf_wfn%crys%recvol, mf_wfn%crys%blat, mf_wfn%crys%bvec, mf_wfn%crys%bdot
      write(10+is) mf_wfn%syms%mtrx(1:3,1:3,1:mf_wfn%syms%ntran)
      write(10+is) mf_wfn%syms%tnp(1:3,1:mf_wfn%syms%ntran)
      write(10+is) (mf_wfn%crys%apos(1:3,ii), mf_wfn%crys%atyp(ii), ii=1,mf_wfn%crys%nat)
      write(10+is) mf_wfn%kp%ngk(1:pseudopot_nospin%nk)
      write(10+is) mf_wfn%kp%w(1:pseudopot_nospin%nk)
      write(10+is) mf_wfn%kp%rk(1:3,1:pseudopot_nospin%nk)

      write(10+is) pseudopot_nospin%ityp(:)
      write(10+is) pseudopot_nospin%nh(:)
      write(10+is) pseudopot_nospin%deeq(:,:,:,is:is)
      call write_binary_gvectors(10+is, mf_wfn%gvec%ng, mf_wfn%gvec%ng, mf_wfn%gvec%components)
    enddo

    call progress_init(prog_info, 'copying VKB', 'projector', &
      pseudopot_nospin%nk*pseudopot%nspin*pseudopot%nkb)
    do ik = 1, pseudopot%nk
      call read_binary_gvectors(7, mf_wfn%kp%ngk(ik), pseudopot%ngkmax, gveck_buf(:,:))
      do is = 1, pseudopot%nspin
        call write_binary_gvectors(10+is, mf_wfn%kp%ngk(ik), pseudopot_nospin%ngkmax, gveck_buf(:,:))
        do ikb = 1, pseudopot%nkb
          call progress_step(prog_info)
          call read_binary_complex_data(7, mf_wfn%kp%ngk(ik), pseudopot%ngkmax, 1, vkb_buf(:,1:1))
          call write_binary_complex_data(10+is, mf_wfn%kp%ngk(ik), pseudopot_nospin%ngkmax, 1, vkb_buf(:,1:1))
        enddo
      enddo
    enddo
    call progress_free(prog_info)
    call close_file(7)
    call close_file(11)
    call close_file(12)

    SAFE_DEALLOCATE(gveck_buf)
    SAFE_DEALLOCATE(vkb_buf)
  endif ! params%has_vkb

  write(*,'(/1x,a)') 'All done!'

end program split_spin
