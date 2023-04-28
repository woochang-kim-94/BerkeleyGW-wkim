!===============================================================================
!
! Program:
!
! (1) paraband            Originally By FHJ      Last Modified Apr/2015 (FHJ)
!
! Input is read from file parabands.inp.
!
!===============================================================================

#include "f_defs.h"
#if defined MPI && !defined USESCALAPACK
  #error ScaLAPACK is required for MPI builds.
#endif

program parabands

#ifdef HDF5
  use hdf5
#endif
  use global_m
  use misc_m,           only: procmem
  use wfn_rho_vxc_io_m, only: read_mf_header, read_binary_gvectors, read_binary_data
  use io_utils_m,       only: progress_info, progress_init, progress_step, &
                              progress_free, print_dealing_with
  use inread_m,         only: pb_params_t, inread
  use input_utils_m,    only: gvec_index
  use distribution_m,   only: distrib_mat_t
  use iteration_data_m, only: iteration_data_t
  use diag_driver_m,    only: check_valid_solver
  use wfn_io_m,         only: wfn_read
#ifdef HDF5
  use wfn_io_m,         only: wfn_write_header, wfn_write_kpt, is_kpt_done, &
                              set_kpt_done
#endif
  use pseudopot_m,      only: pseudopot_t
  use pseudopot_vkb_m,  only: pseudopot_vkb_t, compare_mf_headers
  use hamiltonian_m,    only: ham_solve
  use write_program_header_m

  implicit none

  type(mf_header_t) :: mf, mf_vsc
  SCALAR, allocatable :: wfn_d(:,:,:), vsc(:,:)
  integer, allocatable :: isort(:,:), isortk(:)
  class(pseudopot_t), allocatable :: pseudopot
  real(DP), allocatable :: occ_old(:,:,:), el_old(:,:,:), en(:,:)
  type(distrib_mat_t) :: dm_wfn
  integer :: ii, ik, nb_min, iter, ik_first, nmpinode, ind_diff(2), n_occ, is
#ifdef HDF5
  integer :: error
#endif
  logical :: should_work
  logical :: kpt_done
  integer :: nb_orig
  integer :: restart_iter
  real(DP) :: amem, max_diff
  character(len=128) :: fname, solver_name
  logical :: wfn_out_exists
  type(pb_params_t) :: params
  type(progress_info) :: prog_info
  type(iteration_data_t) :: iter_data


!------------------------
! Initialization
  call peinfo_init()
#ifdef HDF5
  call h5open_f(error)
#endif
  call timacc(0,0)
  call timacc(1,1)
  call write_program_header('ParaBands', .false.)
!#BEGIN_INTERNAL_ONLY
  if (peinf%inode==0) then
    write(6,'(1x,a, a)') 'Git version: ', TOSTRING(GIT_COMMIT)
  endif
!#END_INTERNAL_ONLY


!------------------------
! Read input parameters
  if (peinf%inode.eq.0) write(6,'(1x,a)') 'Reading parameters from file parabands.inp'
  call inread(params)
  call timacc(2,1)


!------------------------
! Open self-consistent potential file
  call timacc(22,1)
  if (peinf%inode==0) then
    write(6,'(1x,a)') 'Reading self-consistent potential from file '//TRUNC(params%fname_vsc)
    call open_file(10, file=params%fname_vsc, status='old', form='unformatted')
  endif
  call read_mf_header(10, mf_vsc, iflavor=SCALARSIZE, sheader='GET', warn=.false., dont_warn_kgrid=.true.)
  ! Legacy header was 'VXC', but we now renamed it to 'VSC'
  if (mf_vsc%sheader(1:3)/='VSC' .and. mf_vsc%sheader(1:3)/='VXC') then
    call die('Invalid header for vsc file: "'//mf_vsc%sheader//'".')
  endif

  SAFE_ALLOCATE(mf_vsc%gvec%components, (3, mf_vsc%gvec%ng))
  call read_binary_gvectors(10, mf_vsc%gvec%ng, mf_vsc%gvec%ng, &
    mf_vsc%gvec%components)
  call gvec_index(mf_vsc%gvec)
  if (params%has_vsc) then
    !call compare_mf_headers('WFN', mf, 'VSC', mf_vsc, .false.)
  endif
  !SAFE_DEALLOCATE_P(mf_vsc%gvec%components)
  SAFE_ALLOCATE(vsc, (mf_vsc%gvec%ng, mf_vsc%kp%nspin))
  call read_binary_data(10, mf_vsc%gvec%ng, mf_vsc%gvec%ng, mf_vsc%kp%nspin, vsc)
  if (peinf%inode==0) call close_file(10)
  call timacc(22,2)


!------------------------
! Read WFN header or generate that from k-point grid
  if (params%has_wfn_in) then
    ! WFN file is available
    call logit('Reading WFN header')
    call timacc(21,1)
    call wfn_read(mf, params, isort) ! Will allocate isort
    call timacc(21,2)
    nb_orig = mf%kp%mnband
    if (peinf%inode==0) then
      call compare_mf_headers('WFN', mf, 'VSC', mf_vsc, .false.)
      !SAFE_DEALLOCATE_P(mf_vsc%gvec%components)
    endif
  endif
  if (params%nb<1) params%nb = minval(mf%kp%ngk)*mf%kp%nspinor
  if (params%nb>minval(mf%kp%ngk)*mf%kp%nspinor) params%nb = minval(mf%kp%ngk)*mf%kp%nspinor
  SAFE_ALLOCATE(isortk, (mf%gvec%ng))


!------------------------
! Setup k-point pools and BLACS stuff
  call logit('Setting up pools')
  call params%kpp%setup(mf%kp%nrk)


!------------------------
! Make sure the solver is valid
  call logit('Validating solver')
  call check_valid_solver(params, solver_name)


!------------------------
! Open Kleinman-Bylander projectors file
  if (params%has_vkb) then
    SAFE_ALLOCATE_CLASS(pseudopot_vkb_t,pseudopot)
    call pseudopot%init(trim(params%fname_vkb), params%kpp, mf)
  endif ! params%has_vkb

  call timacc(2,2)

!------------------------
! Save previous MF energies and occupations. We will reallocate mf%kp%el and
! mf%kp%occ afterwards. If the user is using pseudobands, we need to overwrite
! mf%kp%mnband, and write the correct mf header to file.
  if (params%has_wfn_in) then
    nb_min = min(params%nb, mf%kp%mnband)
    SAFE_ALLOCATE(occ_old, (mf%kp%mnband, mf%kp%nrk, mf%kp%nspin))
    SAFE_ALLOCATE(el_old, (mf%kp%mnband, mf%kp%nrk, mf%kp%nspin))
    occ_old = mf%kp%occ
    el_old = mf%kp%el
  endif
  ! Will reallocate these afterwards
  SAFE_DEALLOCATE(mf%kp%occ)
  SAFE_DEALLOCATE(mf%kp%el)
  mf%kp%mnband = params%nb


!------------------------
! Cleanup old log files
  if (params%kpp%npools>1) then
    params%kpp%iunit=600
  else
    params%kpp%iunit=6
  endif
  if (params%kpp%inode==0) then
    if (params%kpp%iunit/=6) then
      write(fname,'(a,i0,a)') 'parabands_', params%kpp%ipool, '.log'
      call open_file(unit=params%kpp%iunit, file=fname, status='replace')
    endif
    if (params%solver_alg>=100 .and. params%solver_alg<=120) then
      ! Just erase the PRIIME file, because we`ll append to it afterwards
      write(fname,'(a,i0,a)') 'primme_', params%kpp%ipool, '.log'
      call open_file(unit=666, file=fname, status='replace')
      call close_file(666)
    endif
  endif

  call procmem(amem, nmpinode)

  if (peinf%inode==0) then
    write(6,'(/1x,a,f0.3,a)') 'Memory available per MPI task: ', amem/1024.0d0**2, ' MB'
    write(6,'(/1x,a)') 'Problem description:'
    write(6,'(3x,a,3(i0,1x))') 'FFT grid: ', mf%gvec%FFTgrid
    write(6,'(3x,a,i0)') 'Size of RHO G-space (gvec%ng): ', mf%gvec%ng
    write(6,'(3x,a,i0)') 'Max. number of G-vectors per k-point: ', mf%kp%ngkmax
    if (params%has_wfn_in) then
      write(6,'(3x,a,i0)') 'Original number of bands in WFN file: ', nb_orig
    endif
    write(6,'(3x,a,i0)') 'Number of spins: ', mf%kp%nspin
    write(6,'(3x,a,i0)') 'Number of spinor components: ', mf%kp%nspinor
    write(6,'(3x,a,i0)') 'Number of k-points: ', mf%kp%nrk
    if (params%has_wfn_in) then
      write(6,'(3x,a,i0)') 'Max. number of occ states: ', maxval(mf%kp%ifmax)
    endif
    write(6,'(/1x,a)') 'Summary of calculation:'
    write(6,'(3x,a,i0)') 'Number of bands to generate: ', params%nb
    write(6,'(3x,a,i0)') 'Number of k-point pools: ', params%kpp%npools
    write(6,'(3x,a,i0)') 'Number of MPI ranks per pool: ', params%kpp%npes
    write(6,'(3x,a,i0)') 'Number of idle MPI ranks building Ham.: ', &
      peinf%npes - params%kpp%npools*params%kpp%npes
    write(6,'(3x,a,i0)') 'Number of idle MPI ranks diagonalizing Ham.: ', &
      peinf%npes - params%kpp%npools*params%kpp%npes_diag
    write(6,'(3x,a,i0)') 'Matrix block size: ', params%block_sz
    write(6,'(3x,a,i0,a)') 'Diagonalization algorithm: ', params%solver_alg, &
      ' ('//trim(solver_name)//')'
    if (params%kpp%npools>1) then
      write(6,'(/1x,a,a)') 'NOTE: We are using pools, so we`ll minimize the output to keep things neat.'
      write(6,'(7x,a)') 'Detailed information on the calculation will be written to the log files'
      write(6,'(7x,a,i0,a)') 'parabands_{p}.log, one for each pool (0 <= p <= ', params%kpp%npools-1, ').'
    endif
    write(6,*)
    FLUSH(6)
  endif
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

!------------------------
! Check if wfn_out already exists: treats case of "restart from scratch"
  if (peinf%inode==0) then
    inquire(file=params%fname_wfn_out, exist=wfn_out_exists)
  endif
#ifdef MPI
  call MPI_Bcast(wfn_out_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
#endif

!------------------------
! Main part: loop over k-points and get bands
  if (params%kpp%npools>1) then
    call progress_init(prog_info, 'generating bands', 'k-point iteration', &
      (mf%kp%nrk-1)/params%kpp%npools + 1)
  endif
  ! Iterate over blocks of k-points. On each iteration we deal with kpp%npools
  ! k-points simultaneously. The max. number of iterations, properly rounded,
  ! is (mf%kp%nrk-1)/params%kpp%npools + 1

  restart_iter = 1

  do iter = 1, (mf%kp%nrk-1)/params%kpp%npools + 1
    if (params%kpp%npools>1) call progress_step(prog_info)
    ik_first = (iter-1)*params%kpp%npools + 1

    ! For each iteration, figure out what`s my k-point, and whether I
    ! should_work. I don`t have to work if I`m pool-less or if we`re past the
    ! total number of k-points (i.e., I`m idle in this iteration)
    ik = ik_first + params%kpp%ipool
    params%cur_ik = ik
    should_work = params%kpp%ipool>=0 .and. ik<=mf%kp%nrk
    call iter_data%setup_comm(should_work, params%kpp%inode)

    if (params%has_wfn_in) then
      isortk = isort(:,min(ik,mf%kp%nrk))
    endif

    if (allocated(pseudopot)) then
      ! Note that isortk may be modified if we use a VKB file.
      call pseudopot%prepare_kpoints(ik, params%kpp, mf, isortk)
    endif

#ifdef HDF5
    ! ARA: Restart feature
    if (params%has_restart .and. wfn_out_exists) then
      if (peinf%inode==0) kpt_done = is_kpt_done(TRUNC(params%fname_wfn_out), ik)
#ifdef MPI
      call MPI_Bcast(kpt_done, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
#endif
      if (kpt_done) then
        if (peinf%inode==0) then
          write(6,'(/,1x,a,/)') 'This k-point was already calculated: skipping.'
        endif
        cycle
      endif
    endif
#endif

    if (should_work) then
      if (params%kpp%inode==0) then
        call print_dealing_with(ik, mf%kp%nrk, mf%kp%rk(:,ik), 'k', iunit=params%kpp%iunit)
        write(params%kpp%iunit,'(1x,a,i0/)') &
          'Constructing and diagonalizing Hamiltonian H of size ', &
          mf%kp%ngk(ik) * mf%kp%nspinor
        FLUSH(params%kpp%iunit)
      endif

      ! Allocate distributed arrays for wavefunctions
      ! Setup 1d cyclic matrix, where the global number of columns is nb.
      call dm_wfn%setup_mat_1d(mf%kp%ngk(ik)*mf%kp%nspinor, params%nb, params%kpp)
      SAFE_ALLOCATE(wfn_d, (dm_wfn%Ml, dm_wfn%Nl, mf%kp%nspin))
      SAFE_ALLOCATE(en, (params%nb, mf%kp%nspin))

      call logit('Calling ham_solve', params%kpp%inode==0, params%kpp%iunit)
      call timacc(4,1)
      call ham_solve(params, dm_wfn, mf%gvec, mf%kp, mf%crys, isortk, &
        vsc, pseudopot, ik, en, wfn_d)
      call timacc(4,2)
      call logit('Done with ham_solve', params%kpp%inode==0, params%kpp%iunit)

      ! Check solution against input MF energies
      if (params%kpp%inode==0) then
        if (params%solver_alg/=-2) then
          if (params%has_wfn_in) then
            max_diff = maxval(abs(el_old(1:nb_min,ik,:) - en(1:nb_min,:)))
            ind_diff = maxloc(abs(el_old(1:nb_min,ik,:) - en(1:nb_min,:)))
            write(params%kpp%iunit,'(1x,a,es9.3,a,/,1x,a,2(i0,", "),i0,")."/)') &
              'Max. difference between input and output mean-field energies: ', max_diff, &
              ' Ry', 'at (ik, band, spin) = (', ik, ind_diff
            FLUSH(params%kpp%iunit)
            if (max_diff > 1d-2) then
              FLUSH(0)
              write(0,'(/,1x,a)') 'WARNING: Difference between input and output MF energies larger than 1e-2.'
              write(0,'(10x,a)') 'Make sure your input files are correct, otherwise file a bug report.'
              if (mf%kp%nspin*mf%kp%nspinor>1) then
                write(0,'(10x,a)') 'Make sure that you are not constraining the magnetization!'
              endif
              write(0,'(10x,a)') 'SOMETHING PROBABLY WENT *WRONG*!'
              write(0,'(/)')
              FLUSH(0)
            elseif (max_diff > 1d-6) then
              FLUSH(0)
              write(0,'(/,1x,a)') 'WARNING: Difference between input and output MF energies larger than 1e-6.'
              write(0,'(10x,a)') 'Are you sure your input WFNs are well converged?'
              if (mf%kp%nspin*mf%kp%nspinor>1) then
                write(0,'(10x,a)') 'Make sure that you are not constraining the magnetization!'
              endif
              write(0,'(/)')
              FLUSH(0)
            endif
          else
            write(6,'(/1x,2a/)') 'Note: will not compare mean-field energies. ', &
              'Make sure they are reasonable.'
          endif
        else
          write(params%kpp%iunit,'(1x,a/)') 'Note: will not compute energy '//&
            'differences because we are using the dummy solver.'
        endif
      endif

      if (iter==1 .or. restart_iter==1) then ! ARA: this still needs to run in restart
        ! Setup MF energies and occupations
        SAFE_ALLOCATE(mf%kp%occ, (mf%kp%mnband, mf%kp%nrk, mf%kp%nspin))
        SAFE_ALLOCATE(mf%kp%el, (mf%kp%mnband, mf%kp%nrk, mf%kp%nspin))
        mf%kp%occ = 0d0
        mf%kp%el = 0d0
        if (params%has_wfn_in) then
          ! FIXME - count number of electrons instead?
          ! NOTE: we could also figure out the occupations instead of just copying
          mf%kp%occ(1:nb_min,:,:) = occ_old(1:nb_min,:,:)
        endif

#ifdef HDF5
        if (params%has_wfn_out .and. iter==1) then
          ! Note: we write the correct occupations here, even for the case of
          ! pseudobands, but the energies are only written afterwards. We might
          ! want to change this behavior in the future if we don`t have/trust
          ! the previous MF occupations.
          call timacc(3,1)
          call timacc(31,1)
          call logit('Calling wfn_write_header')
          call wfn_write_header(mf, params, isort)
          call logit('Done with wfn_write_header')
          call timacc(31,2)
          call timacc(3,2)
        endif
#endif
      endif

    endif ! should_work

#ifdef HDF5
    if (params%has_wfn_out) then
      call logit('Calling wfn_write_kpt', params%kpp%inode==0, params%kpp%iunit)
      call timacc(3,1)
      call timacc(32,1)
      call wfn_write_kpt(mf, params, iter_data, dm_wfn, ik, wfn_d, en)
      call timacc(32,2)
      call timacc(3,2)
      call logit('Done with wfn_write_kpt', params%kpp%inode==0, params%kpp%iunit)
      if (peinf%inode==0) call set_kpt_done(TRUNC(params%fname_wfn_out), iter)
    endif
#endif
    restart_iter = restart_iter + 1

    SAFE_DEALLOCATE(wfn_d)
    SAFE_DEALLOCATE(en)
    call iter_data%free_comm()
  enddo
  if (params%kpp%npools>1) call progress_free(prog_info)
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

!------------------------
! Clean-up
  if (params%kpp%inode==0.and.params%kpp%npools>1) call close_file(params%kpp%iunit)
#ifdef MPI
  call params%kpp%cleanup()
#endif
  if (allocated(pseudopot)) then
    call pseudopot%free()
  endif
  SAFE_DEALLOCATE(isortk)
  SAFE_DEALLOCATE(occ_old)
  SAFE_DEALLOCATE(el_old)
  call timacc(1,2)

  call report_timing()

!------------------------
! Finish
#ifdef HDF5
  call h5close_f(error)
#endif
#ifdef MPI
  call MPI_Finalize(mpierr)
#endif

contains

  subroutine report_timing()
    character(len=18) :: routnam(50)
    real(DP) :: tsec(2)
    integer :: ncount

    routnam(:) = ''
    routnam(1)='TOTAL:'
    routnam(2)='INPUT:'
    routnam(3)='OUTPUT:'
    routnam(4)='HAMILTONIAN:'
    routnam(20)='-'
    routnam(21)='INPUT WFN:'
    routnam(22)='INPUT VSC:'
    routnam(23)='INPUT VKB HEAD:'
    routnam(24)='INPUT VKB DATA:'
    routnam(30)='-'
    routnam(31)='OUTPUT HEADER:'
    routnam(32)='OUTPUT WFN:'
    routnam(40)='-'
    routnam(41)='HAM KIN:'
    routnam(42)='HAM BUILD LOC:'
    routnam(43)='HAM BUILD NL:'
    routnam(44)='HAM BUILD NL COMM:'
    routnam(45)='HAM DIAG:'

    if(peinf%inode.eq.0) then
      write(6,9000)
      do ii=2,49
        if (trim(routnam(ii))=='-') then
          write(*,'()')
        elseif (len_trim(routnam(ii))>0) then
          call timacc(ii,3,tsec,ncount)
          write(6,9001) routnam(ii), tsec(1), tsec(2), ncount
        endif
      enddo
      write(6,*)
      call timacc(1,3,tsec,ncount)
      write(6,9001) routnam(1), tsec(1), tsec(2), 1
      write(6,*)
9000 format(/,18x,6x,"CPU (s)",8x,"WALL (s)",11x,"#",/)
9001 format(a18,f13.3,3x,f13.3,3x,i9)
    endif

  end subroutine report_timing

end program parabands
