!===============================================================================
!
! setup_subsampling_nns     Originally By FHJ     Last Modified 07/27/2014 (FHJ)
!
!
! Given a WFN file, automatically finds the radial subsampling q-points, their
! weights, and the k-points associated to a WFNq file. The final WFNq file,
! along with WFN, can be used in Epsilon to generate the subsampled eps0mat file.
!
!===============================================================================

#include "f_defs.h"

program setup_subsampling_nns

  use global_m
  use wfn_rho_vxc_io_m
  use random_m
  use irrbz_m
  use fullbz_m
  use checkbz_m
  use wfn_io_hdf5_m
#ifdef HDF5
  use hdf5
#endif
  use subgrp_m
  implicit none

  integer :: narg
  character(len=256) :: fname_wfn, tmpstr
  character(len=5) file_fmt

  real(DP) :: dq(3), dq_inp(3)
  real(DP) :: dq_abs, qmax, bdot(3,3), dq_inp_abs, dq_cur
  real(DP) :: qq(3), qq_min(3), qq_tmp(3), len_min, len_tmp
  real(DP), allocatable  :: qran(:,:), qlen(:), qq_abs(:)
  real(DP) :: dq_abs_min, dq_abs_max, weight, degree
  integer :: nn
  integer :: i1, i2, i3, iq, nq, ii, ni(3), qgrid(3), ndims
  type(mf_header_t) :: mf
  logical :: auto_mode, skip_checkbz, dq_rel, dq_crys, is_periodic(3)
  logical, allocatable :: cond(:)
  type(grid) :: gr
  
  PUSH_SUB(setup_subsampling_nns)

  narg = command_argument_count()
  qgrid(:) = 0
  if (narg/=2.and.narg/=3.and.narg/=4.and.narg/=7.and.narg/=10.and.narg/=11.and.narg/=12.and.narg/=13) then
    write(0,'(1x,2a)') 'Usage: setup_subsamlping_nns.x ', &
#ifdef HDF5
      'ASCII|BIN|HDF5 WFN [nq] [degree] [nq1 nq2 nq3] [dq1 dq2 dq3] [dq_rel] [dq_crys] [ndims]'
#else
      'ASCII|BIN WFN [nq] [degree] [nq1 nq2 nq3] [dq1 dq2 dq3] [dq_rel] [dq_crys] [ndims]'
#endif
    write(0,1)
    write(0,1) 'Required arguments:'
#ifdef HDF5
    write(0,1) '  ASCII|BIN|HDF5: format of the WFN file'
#else
    write(0,1) '  ASCII|BIN: format of the WFN file'
#endif
    write(0,1) '  WFN: WFN file to used in the Epsilon/Sigma calculation'
    write(0,1)
    write(0,1) 'Optional arguments:'
    write(0,1) '  nq: number of radial subsampled q-points (defaults to 10)'
    write(0,1) '  degree: the radial intervals are defined by dq(i) = dq0 * i^{degree} (defaults to 1.0)'
    write(0,1) '  nq1, nq3, nq3: qgrid (defaults to WFN kgrid)'
    write(0,1) '  dq1, dq2, dq3: direction/value of the shift in crystal coords (defaults to 1 0 0)'
    write(0,1) '  dq_rel: if .true., treat the dq as the direction of the q-shift, but find' 
    write(0,1) '       the absolute value automatically depending on the Voronoi cell.'
    write(0,1) '       If .false., the dq array will be the absolute shift. (default is .true.)'
    write(0,1) '  dq_crys: if .true., the dq array is in crystal coordinates, otherwise it is in'
    write(0,1) '       cartesian coordinated. (default is .true.)'
    write(0,1) '  ndims: dimensionality of the problem:'
    write(0,1) '       ndims=3: bulk system, system periodic along all directions'
    write(0,1) '       ndims=2: quasi-2D system, system periodic along a1 and a2'
    write(0,1) '       ndims=1: quasi-2D system, system periodic along a3'
    write(0,1) '       ndims=-1: automatic, determine from qgrid. May not work for large supercells!'
    write(0,1)
1 format(1x,a)
    stop
  endif
  call get_command_argument(1, file_fmt)
  call get_command_argument(2, fname_wfn)
  nq = 10
  if (narg>2) then
    call get_command_argument(3, tmpstr)
    read(tmpstr,*) nq
  endif
  degree = 1d0
  if (narg>3) then
    call get_command_argument(4, tmpstr)
    read(tmpstr,*) degree
  endif
  if (narg>4) then
    call get_command_argument(5, tmpstr)
    read(tmpstr,*) qgrid(1)
    call get_command_argument(6, tmpstr)
    read(tmpstr,*) qgrid(2)
    call get_command_argument(7, tmpstr)
    read(tmpstr,*) qgrid(3)
  endif
  auto_mode = .true.
  if (narg>7) then
    call get_command_argument(8, tmpstr)
    read(tmpstr,*) dq_inp(1)
    call get_command_argument(9, tmpstr)
    read(tmpstr,*) dq_inp(2)
    call get_command_argument(10, tmpstr)
    read(tmpstr,*) dq_inp(3)
    auto_mode = .false.
  endif
  dq_rel = .true.
  if (narg>10) then
    call get_command_argument(11, tmpstr)
    read(tmpstr,*) dq_rel
  endif
  dq_crys = .true.
  if (narg>11) then
    call get_command_argument(12, tmpstr)
    read(tmpstr,*) dq_crys
  endif
  ndims = -1
  if (narg>12) then
    call get_command_argument(13, tmpstr)
    read(tmpstr,*) ndims
  endif

  if (file_fmt=='ASCII') then
    call open_file(unit=11, file=TRIM(fname_wfn), form='formatted', status='old')
    call read_mf_header(11, mf, dont_warn_kgrid=.true.)
    call close_file(11)
  elseif (file_fmt=='BIN  ') then
    call open_file(unit=11, file=TRIM(fname_wfn), form='unformatted', status='old')
    call read_mf_header(11, mf, dont_warn_kgrid=.true.)
    call close_file(11)
#ifdef HDF5
  elseif (file_fmt=='HDF5 ') then
    call h5open_f(i1)
    call read_hdf5_mf_header(TRUNC(fname_wfn), mf, sheader='WFN', dont_warn_kgrid=.true.)
    call h5close_f(i1)
#endif
  else
    call die('Unknown format "'//TRIM(file_fmt)//'". Must be either ASCII, BIN or HDF5.')
  endif
  write(6,'(a,3(1x,i0))') 'WFN k-grid:', mf%kp%kgrid
  if (any(qgrid<1)) qgrid = mf%kp%kgrid
  write(6,'(a,3(1x,i0))') 'Using q-grid:', qgrid
  if (abs(degree)>=1) then
    write(6,'(a,f0.3)') 'Grid degree: dq(i) = dq0 * i^{degree}, degree = ', degree
  else
    write(6,'(a,f6.3)') 'Grid degree: dq(i) = dq0 * i^{degree}, degree =', degree
  endif

  select case (ndims)
    case (3)
      is_periodic(:) = [.true., .true., .true.]
    case (2)
      is_periodic(:) = [.true., .true., .false.]
    case (1)
      is_periodic(:) = [.false., .false., .true.]
    case default
      is_periodic(:) = qgrid(:)>1
      if (all(.not.is_periodic(:))) then
        write(0,*)
        write(0,*) 'ERROR: Cannot determine periodic direction(s) because the qgrid is all 1.'
        write(0,*) '       You will need to manually set the `periodic` flag'
        write(0,*)
        call die('Cannot determine periodic direction(s)')
      endif
  endselect

  if (auto_mode) then
    where (is_periodic)
      dq_inp = 1d0
    elsewhere
      dq_inp = 0d0
    endwhere
  endif ! auto_mode
  if (dq_rel) write(6,'(a,3(1x,f5.3))') 'Using a relative shift in the direction:', dq_inp
  write(6,*)

  ! FHJ: sample Voronoi cell
  nn = 10000000
  write(6,'(a,i0,a)') 'Sampling q=0 Voronoi cell stocastically with ', nn, ' points.'
  SAFE_ALLOCATE(qran, (3,nn))
  SAFE_ALLOCATE(qlen, (nn))
  call genrand_init(put=5000)
  do iq = 1, nn
    do ii = 1, 3
      call genrand_real4(qran(ii, iq))
    enddo
  enddo
  do ii = 1, 3
    if (.not.is_periodic(ii)) qran(ii,:) = 0d0
  enddo

  ni(:) = 1
  where (.not.is_periodic)
    ni(:) = 0
  endwhere

  do ii = 1, 3
    bdot(:,ii) = mf%crys%bdot(:,ii)/qgrid(ii)
  enddo
  do ii = 1, 3
    bdot(ii,:) = bdot(ii,:)/qgrid(ii)
  enddo
  do iq = 1, nn
    len_min = INF
    qq(:) = qran(:,iq)
    do i1 = -ni(1), ni(1)
      qq_tmp(1) = qq(1) - i1
      do i2 = -ni(2), ni(2)
        qq_tmp(2) = qq(2) - i2
        do i3 = -ni(3), ni(3)
          qq_tmp(3) = qq(3) - i3
          len_tmp = DOT_PRODUCT(qq_tmp,MATMUL(bdot,qq_tmp))
          if (len_tmp<len_min) then
            len_min = len_tmp
            qq_min(:) = qq_tmp(:)
          endif
        enddo
      enddo
    enddo
    qlen(iq) = sqrt(len_min)
    qran(1:3,iq) = qq_min(1:3)
  enddo ! iq
  qmax = maxval(qlen)
  write(6,'(2x,a,es12.6,a)') 'Maximum |q| in Voronoi cell = ', qmax,' 1/bohr'


  if (.not.dq_crys) then
    ! A^T B = 2*pi*I
    ! B x_crys = x_cart
    ! x_crys = 1/(2*pi) A^T x_cart
    dq_inp = 1d0/(2d0*PI_D) * mf%crys%alat * MATMUL(transpose(mf%crys%avec), dq_inp)
  endif
  dq_inp_abs = sqrt(DOT_PRODUCT(dq_inp,MATMUL(mf%crys%bdot,dq_inp)))
  SAFE_ALLOCATE(qq_abs, (nq))

  if (dq_rel) then
    ! We want to bin the random q-points radially and find the appropriate weights
    ! |           *           |           *        ...        |
    ! 0           dq         2dq         3dq          (2nq)dq = qmax
    ! 
    ! where: * is a point where we`ll calculated epsinv(q)
    !        | is the boundary of each interval
    ! dq_abs = qmax / (2d0 * sum_(i=1)^n i^{degree})
    dq_abs = 0d0
    do iq = 1, nq
      dq_abs = dq_abs + dble(iq)**degree
    enddo
    dq_abs = qmax / (2d0 * dq_abs)
    dq = dq_inp/dq_inp_abs * dq_abs
  else
    dq = dq_inp
    dq_abs = dq_inp_abs
  endif

  write(6,'(2x,a,3(1x,es22.15),a)') 'dq = ', dq, '  (crystal units)'
  write(6,'(2x,a,3(1x,es22.15),a)') '   = ', mf%crys%blat*MATMUL(mf%crys%bvec, dq), '  1/bohr'
  write(6,'(2x,a,1x,es22.15,1x,a)') '|dq| = ', dq_abs, '  1/bohr'
  write(6,*)

  write(6,'(a,a)') 'Determining q0 weights.'
  SAFE_ALLOCATE(cond, (nn))
  dq_abs_min = 0d0
  call open_file(unit=13, file='subweights.dat', form='formatted', status='replace')
  write(13,'(i0)') nq
  do iq=1,nq
    dq_cur = 2d0*dq_abs*(iq**degree)
    dq_abs_max = dq_abs_min + dq_cur
    cond = qlen>=dq_abs_min.and.qlen<dq_abs_max
    weight = count(cond)
    qq_abs(iq) = dq_abs_min + 0.5d0*dq_cur
    write(6,'(2x,3(a,es12.6,2x))') '|q| = ', qq_abs(iq), 'weight = ', weight
    write(13,'(2(es12.6,2x))') weight, qq_abs(iq)
    dq_abs_min = dq_abs_max
  enddo
  SAFE_DEALLOCATE(cond)
  call close_file(13)
  write(6,'(2x,2a)') 'Wrote q0 weights to subweights.dat'

  call open_file(unit=13, file='epsilon_q0s.inp', form='formatted', status='replace')
  write(13,*)
  write(13,'(a)') 'subsample'
  write(13,'(a)') 'begin qpoints'
  do iq = 1, nq
    write(13,'(3(es22.15,1x),a)') qq_abs(iq)/dq_abs*dq, '1 1'
  enddo
  write(13,'(a)') 'end'
  write(13,*)
  call close_file(13)
  write(6,'(2x,a)') 'Wrote lines for epsilon.inp in epsilon_q0s.inp'
  write(6,*)

  write (6,'(a)') 'Generating k-points for WFNq file(s)'
  ! FHJ: Use symmetries to unfold BZ from WFN
  gr%nr = mf%kp%nrk
  SAFE_ALLOCATE(gr%r, (3, gr%nr))
  gr%r = mf%kp%rk
  call fullbz(mf%crys, mf%syms, gr, mf%syms%ntran, skip_checkbz, wigner_seitz=.false., paranoid=.true.)
  skip_checkbz = .true.
  if (.not.skip_checkbz) call checkbz(gr%nf, gr%f, mf%kp%kgrid, mf%kp%shift, &
    mf%crys%bdot, TRUNC(fname_wfn), 'k', .false., .false., .false.)
  ! FHJ: Find subgroup that leaves dq invariant
  call subgrp(dq, mf%syms)
  ! FHJ: For each dq, displace k-points by dq and use symmetries in the
  ! subgroup to fold the k-points
  do ii=0,nq
    call gen_kpoints_file(ii)
  enddo
  write(6,*)

  POP_SUB(setup_subsampling_nns)

  contains

    subroutine gen_kpoints_file(iq_in)
      integer, intent(in) :: iq_in

      real(DP), allocatable :: kpts_new(:,:)
      integer :: nk_new, nk_fold, iq_min, iq_max, ik
      integer, allocatable :: neq(:), indrk(:)
      character(len=64) :: fname_kpts

      POP_SUB(subsample_plan.gen_kpoints_file)

      if (iq_in>0) then
        nk_new = gr%nf
        iq_min = iq_in
        iq_max = iq_in
        write(fname_kpts,'(a,i3.3,a)') 'kpoints_', iq_in, '.dat'
      else
        nk_new = gr%nf*nq
        iq_min = 1
        iq_max = nq
        fname_kpts = 'kpoints_all.dat'
      endif

      SAFE_ALLOCATE(kpts_new, (3,nk_new))
      SAFE_ALLOCATE(indrk, (nk_new))
      SAFE_ALLOCATE(neq, (nk_new))
      nk_new = 0
      do iq=iq_min,iq_max
        qq = qq_abs(iq)/dq_abs*dq
        forall(ik=1:gr%nf) kpts_new(1:3, nk_new+ik) = gr%f(1:3, ik) + qq(1:3)
        nk_new = nk_new + gr%nf
      enddo
      call irrbz(mf%syms, nk_new, kpts_new, nk_fold, neq, indrk)

      call open_file(unit=13, file=TRUNC(fname_kpts), form='formatted', status='replace')
      write(13,'(a)') 'K_POINTS crystal'
      write(13,'(i0)') nk_fold
      do iq=1, nk_fold
        write(13, '(3(f13.9),f6.1)') kpts_new(:,indrk(iq)), dble(neq(iq))
      enddo
      call close_file(13)
      write(6,'(2x,2a)') 'Wrote kpoints to ', TRUNC(fname_kpts)

      SAFE_DEALLOCATE(kpts_new)
      SAFE_DEALLOCATE(neq)
      SAFE_DEALLOCATE(indrk)

      PUSH_SUB(subsample_plan.gen_kpoints_file)

    end subroutine gen_kpoints_file

end program setup_subsampling_nns
