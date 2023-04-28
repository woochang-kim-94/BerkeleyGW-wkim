!=======================================================================
!
! setup_subsampling_csi  Originally By DYQ Last Modified 02/09/2015 (DYQ)
!
!
! Given a WFN_co and WFN files and a desired fine k-grid nk_fi_x nk_fi_y nk_fi_z, automatically find
! k-point and q-points required for clustered subsampling.
! Input files:
!      WFN_co: file used to generate BSE kernel on coarse grid
!      WFN: WFN file usedfor epsilon, usually the same as WFN_co
!      
! Output files:
!      kpoints_sub_*.dat: Cluster of subsampled points around each coarse points. One file
!                         is written for each coarse point on the full grid. Use to generate
!                         WFN_sub wavefunctions. Run kernel for each WFN_sub.
!      epsilon_q0s.inp: Part of epsilon.inp file used to generate epsmat compatible with subsampled
!                       k-points.
!      kpoints_wfnq.dat: Contains list of k-points needed to generate WFNq for all q points
!      kpoints_wfnq_*.dat: Separated lists of k points needed to generate WFNq for each q point in epsilon_q0s.inp
!                          kpoints_wfnq.dat should contain combined list of all kpoints in all kpoints_wfnq_*.dat files
!      subsample.inp: Part of subsample.inp
!
!
! TODO:
! - Support HDF5
! - Support subsampling along arbitrary directions (currently y and xy only)
!
!===============================================================================

#include "f_defs.h"

program setup_subsampling_csi

  use global_m
  use wfn_rho_vxc_io_m
  use random_m
  use irrbz_m
  use fullbz_m
  use checkbz_m
  use subgrp_m
  implicit none

  integer :: narg, nk_fi(3),nk_sub, nk_old, nk_new,nk_fold,nk_tot,use_syms,nsub_fact, dir
  character(len=256) :: fname_wfn_co,fname_wfn, fname_kpts,tmpstr
  character(len=5) file_fmt

  real(DP) :: dq(3), delta_co(3),delta_fi(3),qshift(3), delta_min
  real(DP) :: sub_len
  real(DP), allocatable  :: kp_sub(:,:),kpts_new(:,:),kpts_all(:,:),qshifts(:,:)
  integer :: ii,ik,iq,iq2,nk_sub_xy,nsub_max,ik_sub
  type(mf_header_t) :: mf
  logical :: skip_checkbz
  integer, allocatable :: neq(:), indrk(:),neq_all(:)
  type(grid) :: gr,gr_co
  
  PUSH_SUB(setup_subsampling_csi)

  ! Input args
  narg = command_argument_count()
  if (narg.lt.8 .or. narg.gt.12) then
    write(0,*) 'Usage: setup_subsampling_csi.x ASCII|BIN|HDF5 WFN_co WFN'
    write(0,*) '       nk_fi_1 nk_fi_2 nk_fi_3  nsub_factor direction'
    write(0,*) '       [qshift_1 qshift_2 qshift_3  use_syms]'
    write(0,*)
    write(0,*) 'Required arguments:'
    write(0,*) '  ASCII|BIN|HDF5: format of the WFN file'
    write(0,*) '  WFN_co: WFN_co file to used in the Kernel calculation'
    write(0,*) '  WFN: WFN file used in the Epsilon/Sigma calculation'
    write(0,*) '       (typically the same as WFN_co)'
    write(0,*) '  nk_fi_1 nk_fi_2 nk_fi_3: the fine k-grid desired for absorption'
    write(0,*) '                           along each reciprocal lattice direction b_i'
    write(0,*) '  nsub_factor: number of subsampled points = nk_fi/nk_co*nsub_factor'
    write(0,*) '  direction = 1 (b1-axis), 2 (b2-axis), 3 (diagonal)'
    write(0,*)
    write(0,*) 'Optional arguments:'
    write(0,*) '  qshift_1 qshift_2 qshift_3: shift applied to subsampled kpoints'
    write(0,*) '  use_syms = 1 (default), to unfold WFN_co using symmetries'
    write(0,*) '           = 0 do not unfold WFN_co using symmetries'
    write(0,*)
    write(0,*) 'Output files:'
    write(0,*) '  subsample.inp: Used for subsequent absorption calculations. See documentation'
    write(0,*) '  kpoints_sub_????.dat: Contains the k-points necessary for each subsampled'
    write(0,*) '                        WFN file, each around a different coarse k point.'
    write(0,*) '  kpoint_wfnq.dat: Contains all k-points necessary for the calculation'
    write(0,*) '                   of the dielectric matrix for all subsampled q points.'
    write(0,*) '  kpoints_wfnq_????.dat: Same as `kpoint_wfnq.dat`, but containing only the'
    write(0,*) '                         k-points necessary for each q-point calculation,'
    write(0,*) '                         associated with the q-point ????. Use these files'
    write(0,*) '                         to split your `epsilon` calculation into individual'
    write(0,*) '                         q points if you are running out of memory.'
    write(0,*) '  epsilon_q0s.inp: Contains part of the `epsilon.inp` file needed'
    write(0,*) '                   to generate epsmat for the clustered points.'
    write(0,*)
    stop
  endif
  call get_command_argument(1, file_fmt)
  call get_command_argument(2, fname_wfn_co)
  call get_command_argument(3, fname_wfn)
  call get_command_argument(4, tmpstr)
  read(tmpstr,*) nk_fi(1)
  call get_command_argument(5, tmpstr)
  read(tmpstr,*) nk_fi(2)
  call get_command_argument(6, tmpstr)
  read(tmpstr,*) nk_fi(3)
  call get_command_argument(7, tmpstr)
  read(tmpstr,*) nsub_fact
  call get_command_argument(8, tmpstr)
  read(tmpstr,*) dir
  qshift(:) = 0.d0
  use_syms=1
  if (narg.eq.9) then
    call get_command_argument(9, tmpstr)
    read(tmpstr,*) use_syms
  elseif (narg.gt.9) then
    call get_command_argument(9, tmpstr)
    read(tmpstr,*) qshift(1)
    call get_command_argument(10, tmpstr)
    read(tmpstr,*) qshift(2)
    call get_command_argument(11, tmpstr)
    read(tmpstr,*) qshift(3)
    if (narg.eq.12) then
      call get_command_argument(12, tmpstr)
      read(tmpstr,*) use_syms
    endif
    write(6,*) "qshift = ",qshift
  endif

  ! Read header of WFN_co
  if (file_fmt=='ASCII') then
    call open_file(unit=11, file=TRIM(fname_wfn_co), form='formatted', status='old')
    call read_mf_header(11, mf)
    call close_file(11)
  elseif (file_fmt=='BIN  ') then
    call open_file(unit=11, file=TRIM(fname_wfn_co), form='unformatted', status='old')
    call read_mf_header(11, mf)
    call close_file(11)
  elseif (file_fmt=='HDF5 ') then
    call die('HDF5 not implemented yet')
  else
    call die('Unknown format "'//TRIM(file_fmt)//'". Must be either ASCII, BIN or HDF5.')
  endif
  write(6,'(1x,a,3(1x,i0))') 'Read WFN_co k-grid:', mf%kp%kgrid

  ! Determine cluster size and sampling fineness from k-grid
  do ii=1,3
    delta_co(ii) = 1d0/mf%kp%kgrid(ii)
    delta_fi(ii) = 1d0/nk_fi(ii) / nsub_fact
  enddo
  ! FHJ: Do not increase fineness along b3 direction
  delta_fi(3) = 1d0
  delta_min = minval(delta_co) - minval(delta_fi)
  nsub_max = minval(nk_fi)/delta_min * nsub_fact * 2.
  nsub_max = max(nsub_max, int(delta_min / minval(delta_fi))) + 1

  ! Unfold WFN_co
  gr_co%nr = mf%kp%nrk
  SAFE_ALLOCATE(gr_co%r, (3, gr_co%nr))
  gr_co%r = mf%kp%rk
  if (use_syms.eq.1) then
    call fullbz(mf%crys, mf%syms, gr_co, mf%syms%ntran, skip_checkbz, wigner_seitz=.false., paranoid=.true.)
  else
    call fullbz(mf%crys, mf%syms, gr_co, 1, skip_checkbz, wigner_seitz=.false., paranoid=.true.)
  endif
  if (.not.skip_checkbz) call checkbz(gr_co%nf, gr_co%f, mf%kp%kgrid, mf%kp%shift, &
    mf%crys%bdot, TRUNC(fname_wfn), 'k', .false., .false., .false.)

  ! Print out list of cluster points
  SAFE_ALLOCATE(kp_sub,(3,nsub_max))
  SAFE_ALLOCATE(qshifts,(3,nsub_max))
  do ik=1,gr_co%nf
    kp_sub = 0.d0
    qshifts = 0.d0
    sub_len = 0.d0
    nk_sub = 0
    do while (sub_len<delta_min)
      nk_sub=nk_sub+1
      kp_sub(:,nk_sub) = gr_co%f(:,ik)
      !write(6,*) '!!', kp_sub(:,nk_sub)
      if (dir.eq.3) then
        kp_sub(2,nk_sub) = kp_sub(2,nk_sub) + delta_fi(2)*nk_sub
        kp_sub(1,nk_sub) = kp_sub(1,nk_sub) + delta_fi(1)*nk_sub
        qshifts(2,nk_sub) = delta_fi(2)*nk_sub
        qshifts(1,nk_sub) = delta_fi(1)*nk_sub
      else
        kp_sub(dir,nk_sub) = kp_sub(dir,nk_sub) + delta_fi(dir)*nk_sub
        qshifts(dir,nk_sub) = delta_fi(dir)*nk_sub
      endif
      sub_len=sub_len+minval(delta_fi)
    enddo
    write(fname_kpts,'(a,i4.4,a)') 'kpoints_sub_', ik, '.dat'
    call open_file(unit=12,file=fname_kpts,form='formatted',status='replace')
    write(12,'(a)') 'K_POINTS crystal'
    write(12,'(i8)') nk_sub+1
    write(12, '(3(f13.9),f6.1)') gr_co%f(:,ik) + qshift(:), 1d0
    do ik_sub=1,nk_sub
      write(12, '(3(f13.9),f6.1)') kp_sub(:,ik_sub) + qshift(:), 1d0
    enddo
    call close_file(12)
  enddo ! loop over ik_co
  write(6,'(1x,a,i0)') "Number of subsampled points per coarse point: ", nk_sub

  call open_file(unit=15,file='subsample.inp',form='formatted',status='replace')
  write(15,'(i8)') nk_sub+1
  write(15,'(i8)') gr_co%nf
  do ik=1,gr_co%nf
    write(15, '(3(f13.9))') gr_co%f(:,ik)
  enddo
  write(15,'(a)') "! List subsampled bsemat.h5 files here"
  write(15,'(a)') "! List subsampled WFN files here"
  write(15,'(a)') "! If using finite Q, list subsampled WFNq files here. Otherwise delete this line."
  call close_file(15)

  ! Read header for WFN file
  if (file_fmt=='ASCII') then
    call open_file(unit=11, file=TRIM(fname_wfn), form='formatted', status='old')
    call read_mf_header(11, mf)
    call close_file(11)
  elseif (file_fmt=='BIN  ') then
    call open_file(unit=11, file=TRIM(fname_wfn), form='unformatted', status='old')
    call read_mf_header(11, mf)
    call close_file(11)
  elseif (file_fmt=='HDF5 ') then
    call die('HDF5 not implemented yet')
  else
    call die('Unknown format "'//TRIM(file_fmt)//'". Must be either ASCII, BIN or HDF5.')
  endif
  write(6,'(1x,a,3(1x,i0))') 'Read WFN k-grid:', mf%kp%kgrid


  ! Find q-points
  call open_file(unit=13,file='epsilon_q0s.inp',form='formatted',status='replace')
  call open_file(unit=14,file='kpoints_wfnq.dat',form='formatted',status='replace')
  write(13,*)
  write(13,'(a)') 'subsample'
  write(13,'(a,i0)') 'number_qpoints ', nk_sub*2
  write(13,'(a)') 'begin qpoints'

  ! Find k-points for WFNq
  gr%nr = mf%kp%nrk
  SAFE_ALLOCATE(gr%r, (3, gr%nr))
  gr%r = mf%kp%rk
  call fullbz(mf%crys, mf%syms, gr, mf%syms%ntran, skip_checkbz, wigner_seitz=.false., paranoid=.true.)
  if (.not.skip_checkbz) call checkbz(gr%nf, gr%f, mf%kp%kgrid, mf%kp%shift, &
    mf%crys%bdot, TRUNC(fname_wfn), 'k', .false., .false., .false.)

  nk_new = gr%nf*(nk_sub)*(nk_sub+1)
  SAFE_ALLOCATE(kpts_new, (3,gr%nf))
  SAFE_ALLOCATE(indrk, (gr%nf))
  SAFE_ALLOCATE(neq, (gr%nf))    
  SAFE_ALLOCATE(kpts_all, (3,nk_new*2))
  SAFE_ALLOCATE(neq_all, (nk_new*2))    
  nk_new = 0
  nk_old = 0
  nk_tot = 0
  dq = 0
  write(6,'(1x,a,i0)') "Number of subsampled points: ", nk_sub 
  do iq=1,nk_sub
    if (dir.eq.3) then
      dq(1) = delta_fi(1)*iq
      dq(2) = delta_fi(2)*iq
    else
      dq(dir) = delta_fi(dir)*iq
    endif
      write(13,'(3(f13.9),a)') dq(:), ' 1 1'

      ! Find subgroup that leaves dq invariant
      call subgrp(dq,mf%syms)
      forall(ik=1:gr%nf) kpts_new(1:3,ik) = gr%f(1:3, ik) + dq
      nk_old = nk_new+1
      nk_new = nk_new + gr%nf
      ! Fold shifted points to irreducible wedge using subgroup symmetry
      call irrbz(mf%syms, gr%nf, kpts_new, nk_fold, neq, indrk)
      kpts_all(1:3,nk_tot+1:nk_tot+nk_fold) = kpts_new(1:3,indrk(1:nk_fold))
      neq_all(nk_tot+1:nk_tot+nk_fold) = neq(1:nk_fold)
      nk_tot = nk_tot + nk_fold

      ! write a list of kpoints corresponding to each q-point
      write(fname_kpts,'(a,i4.4,a)') 'kpoints_wfnq_', iq*2-1, '.dat'
      call open_file(unit=17,file=fname_kpts,form='formatted',status='replace')
      write(17,'(a)') 'K_POINTS crystal'
      write(17,'(i0)') nk_fold
      do ik=1,nk_fold
        write(17, '(3(f13.9),f6.1)') kpts_new(:,indrk(ik)), dble(neq(ik))
      enddo
      call close_file(17)

      ! Time-reversed q point
      if (dir.eq.3) THEN
        dq(1) = 1.0d0-delta_fi(1)*iq
        dq(2) = 1.0d0-delta_fi(2)*iq
      else
        dq(dir) = 1.0d0-delta_fi(dir)*iq
      endif
      write(13,'(3(f13.9),a)') dq(:), ' 1 1'
      ! Find subgroup that leaves dq invariant
      call subgrp(dq,mf%syms)
      forall(ik=1:gr%nf) kpts_new(1:3,ik) = gr%f(1:3, ik) + dq
      nk_old = nk_new+1
      nk_new = nk_new + gr%nf
      ! Fold shifted points to irreducible wedge using subgroup symmetry
      call irrbz(mf%syms, gr%nf, kpts_new, nk_fold, neq, indrk)
      kpts_all(1:3,nk_tot+1:nk_tot+nk_fold) = kpts_new(1:3,indrk(1:nk_fold))
      neq_all(nk_tot+1:nk_tot+nk_fold) = neq(1:nk_fold)
      nk_tot = nk_tot + nk_fold

      write(fname_kpts,'(a,i4.4,a)') 'kpoints_wfnq_', iq*2, '.dat'
      call open_file(unit=17,file=fname_kpts,form='formatted',status='replace')
      write(17,'(a)') 'K_POINTS crystal'
      write(17,'(i0)') nk_fold
      do ik=1,nk_fold
        write(17, '(3(f13.9),f6.1)') kpts_new(:,indrk(ik)), dble(neq(ik))
      enddo
      call close_file(17)      

    enddo
  write(13,'(a)') 'end'
  write(13,*)
  call close_file(13)

  write(6,'(1x,a,i0)') "Total kpoints in WFNq: ", nk_tot
  write(14,'(a)') 'K_POINTS crystal'
  write(14,'(i0)') nk_tot
  do ik=1,nk_tot
    write(14, '(3(f13.9),f6.1)') kpts_all(:,ik), dble(neq_all(ik))
  enddo
  call close_file(14)

  SAFE_DEALLOCATE(kpts_new)
  SAFE_DEALLOCATE(kpts_all)
  SAFE_DEALLOCATE(neq)
  SAFE_DEALLOCATE(neq_all)
  SAFE_DEALLOCATE(indrk)
  SAFE_DEALLOCATE(kp_sub)
  SAFE_DEALLOCATE(qshifts)

  POP_SUB(setup_subsampling_csi)

end program setup_subsampling_csi
