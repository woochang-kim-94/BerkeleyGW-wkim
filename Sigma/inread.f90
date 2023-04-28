!===============================================================================
!
! Routines:
!
! (1) inread()  Originally By ?         Last Modified 7/8/2008 (JRD)
!
!     Reads parameters from sigma.inp for the current job.
!
!===============================================================================

#include "f_defs.h"

module inread_m

  use algos_sigma_m
  use global_m
  use inread_common_m
  use references_m
  use scissors_m
  implicit none

  private

  public :: &
    inread

contains

subroutine inread(sig)
  type (siginfo), intent(out) :: sig

  character*256 :: blockword,keyword,line,errmsg,tmpstr
  integer :: ii,jj,kk,nkpt_read,nqpt_read,ndiag_read,noff_read
  integer :: nbnmin,nbnmax,itestq,iostat
  
  real(DP) :: div,max_freq_eval
  real(DP), allocatable :: kpt_read(:,:),qpt_read(:,:)
  integer, allocatable :: off1_read(:),off2_read(:),off3_read(:)
  integer, allocatable :: off_ep_read(:) ! ZL: add for EP bra ket
  integer, allocatable :: band_occ(:),diag(:)
  integer :: spinmin,spinmax
  logical :: occ_set, q0vec_read, VXC_exists, found, old_grid_set

  integer :: ik ! ZL: for loop generating k_phonq
  real(DP), allocatable :: k_phonq_read(:)

  PUSH_SUB(inread)

#ifdef MPI
  ! Non-root nodes should wait for root to read the whole file.
  ! That way, we can be sure root gets a chance to write errors before
  ! any die call is issued by another node. Root calls MPI_Barrier below.
  if(peinf%inode /= 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

!----------------------
! Allocate arrays

  SAFE_ALLOCATE(kpt_read, (3,MAX_KPTS))
  SAFE_ALLOCATE(qpt_read, (3,MAX_KPTS))
  SAFE_ALLOCATE(off1_read, (MAX_BANDS))
  SAFE_ALLOCATE(off2_read, (MAX_BANDS))
  SAFE_ALLOCATE(off3_read, (MAX_BANDS))
  SAFE_ALLOCATE(off_ep_read, (MAX_BANDS))
  SAFE_ALLOCATE(band_occ, (MAX_BANDS))
  SAFE_ALLOCATE(diag, (MAX_BANDS))

!----------------------
! Initialize sig%igamma

  sig%igamma = 0

!----------------------
! Initialize wcoul0
  sig%wcoul0 = ZERO


!---------------------------------
! Set default values

  occ_set = .false.
  q0vec_read = .false.
  sig%invalid_gpp_mode = -1
  nkpt_read=0
  nqpt_read=0
  noff_read=0
  band_occ=0
  nbnmin=0
  nbnmax=0
  spinmin=1
  spinmax=1
  diag=0
#ifdef HDF5
  sig%use_hdf5 = .true.
#else
  sig%use_hdf5 = .false.
#endif
  sig%wfn_hdf5 = .false.
  sig%loff=0
  sig%toff=0
  sig%freq_dep=1
  sig%freq_dep_method=2
  sig%nFreq=1
  sig%exact_ch=-1
  sig%freq_grid_shift=2
  sig%freqevalmin=0d0
  old_grid_set=.false.
  sig%freqevalstep=0.2d0
  max_freq_eval=2d0
  sig%nfreqeval=1
  sig%iuseremainder=0
  sig%nkn=0
  sig%nq=0
  sig%nq0=0
  sig%nq1=0
  sig%subsample=.false.
!#BEGIN_INTERNAL_ONLY
  ! ZL: add for EP
!#END_INTERNAL_ONLY
  sig%elph=.false.
  sig%ep_bracket=0
  sig%nphonq = 0
  sig%ndiag=0
  sig%noffdiag=0
  sig%nspin=0
  sig%nfreq_imag=0
  sig%cd_int_method=0
  call scissors_zero(sig%scis)
  call scissors_zero(sig%scis_outer)
  sig%spl_tck%n=0 
  sig%spl_tck_outer%n=0
  sig%avgpot=0.0d0
  sig%avgpot_outer=0.0d0
  sig%tol = TOL_Degeneracy
  sig%use_xdat=.false.
  sig%use_vxcdat=.TRUE.
  sig%use_vxc2dat=.TRUE.
  sig%is_EXX=.FALSE.
  sig%occ_broadening_EXX=0.1 ! OAH for smearing
  sig%do_smearing_EXX=.false.
  sig%use_kihdat=.false. ! ZL: default false. Default is still to use vxc.dat, to be changed in the future TODO
  ! ZL: note, we should always keep vxcdat and VXC options but may move them from default
  !     There are other parts of code that strongly depend on VXC
  sig%ecuts=0.d0
  sig%ecutb=0.d0
  sig%ntband=0  ! ZL: number of total bands used in sigma
  sig%ncore_excl=0
  sig%nvband=0
  sig%fdf=2
  sig%dw=1.0d0
  sig%xfrac=1.0d0
  sig%gamma=0.5d0
  sig%sexcut=4.0d0
  sig%icutv=TRUNC_NONE
  sig%iscreen=SCREEN_SEMICOND
  sig%iwritecoul=0
  sig%truncval(:)=0.0d0
  sig%ncrit=0
  sig%efermi_input=0.0d0
  sig%rfermi=.true.
  sig%fullConvLog=0
  sig%qgrid(1:3)=0
  sig%avgcut=-1 !FHJ: <0 means "auto", see code below
  sig%freplacebz=.false.
  sig%fwritebz=.false.
  sig%degeneracy_check_override=.false.
  sig%offdiagsym=.true.
  sig%qgridsym=.true.
  sig%die_outside_sphere=.true.
  sig%averagew=.true.
  peinf%npools=0
  sig%eqp_corrections=.false.
  sig%eqp_outer_corrections=.false.
  sig%coulomb_mod%short_range_frac_fock=1.0d0
  sig%coulomb_mod%long_range_frac_fock=1.0d0
  sig%coulomb_mod%screening_length=0.0d0
  sig%coul_mod_flag=.false.
  sig%sigma_correction=.false.
  sig%symmetrize=.true.
  sig%do_sigma_subspace=.false.
  sig%sub_collective_redistr=.false.
  sig%tolerant_value=.false.

! Set default values for algo acceleration
  call set_algos_to_cpu()
#if defined(OPENACC) || defined(OMP_TARGET)
  ! GPU settings
  sig%accel_mtxel_band_block_size = 20
  sig%gpp_band_block_size = 16
  sig%gpp_ig_block_size = 512
  call set_algos_to_best_available_gpu()
#endif

!---------------------------------
! Never-ending loop...

  do while(0.eq.0)

! Actually the loop ends when the end of the file is reached
    
    read(55,'(a256)',iostat=iostat) line
    if(iostat < 0) exit

! Skip comment lines

    if(len_trim(line).eq.0) cycle
    if(line(1:1).eq.'#') cycle

! Determine keyword

    keyword=line(1:scan(line," ")-1)
    line=adjustl(line(scan(line," ")+1:256))

    if (trim(keyword).eq.'begin') then
      blockword=line(1:scan(line," ")-1)
      ii=0
      do while(trim(line).ne.'end')
        read(55,'(a256)',iostat=iostat) line
        if(iostat /= 0) then
          write(errmsg,'(3a)') 'The end of the file was reached while reading elements of the ', &
            trim(blockword),' block.'
          call die(errmsg, only_root_writes = .true.)
        endif
        
        if(trim(line).ne.'end') then
          ii=ii+1
          if(trim(blockword).eq.'kpoints') then
            call check_bounds_nkq(ii, 'k', 'begin kpoints')
            read(line,*,err=112) (kpt_read(jj,ii),jj=1,3),div
            kpt_read(1:3,ii)=kpt_read(1:3,ii)/div
          elseif(trim(blockword).eq.'qpoints') then
            call check_bounds_nkq(ii, 'q', 'begin qpoints')
            read(line,*,err=112) (qpt_read(jj,ii),jj=1,3),div,itestq
            qpt_read(1:3,ii)=qpt_read(1:3,ii)/div
            ! for Hartree-Fock
            if (itestq.eq.1) then
              sig%nq0 = sig%nq0 + 1
              if (.not.q0vec_read) then
                sig%q0vec(1:3)=qpt_read(1:3,ii)
                q0vec_read = .true.
              endif
            endif
          elseif(trim(blockword).eq.'diag') then
            call check_bounds_nbands(ii, 'begin diag')
            read(line,*,err=112) diag(ii)
          elseif(trim(blockword).eq.'offdiag') then
            call check_bounds_nbands(ii, 'begin offdiag')
            read(line,*,err=112) off1_read(ii),off2_read(ii),off3_read(ii)
          else
            write(errmsg,'(3a)') 'Unexpected blockword ', trim(blockword), ' was found in sigma.inp.'
            call die(errmsg, only_root_writes = .true.)
          end if
        end if
      end do
      if(trim(blockword).eq.'kpoints') then
        nkpt_read=ii
      elseif(trim(blockword).eq.'qpoints') then
        nqpt_read=ii
      elseif(trim(blockword).eq.'diag') then
        ndiag_read=ii
      elseif(trim(blockword).eq.'offdiag') then
        noff_read=ii
      endif

! Spline scissors
    elseif(trim(keyword).eq.'spline_scissors') then
      if (peinf%inode==0) then
        write(6,*) 'Reading spline coefficients from `spline_scissors.dat`'
        write(6,*)
        ! read number of pts, knots, coefficients, degree
        call open_file(20,file='spline_scissors.dat',form='formatted',status='old')
        read(20,*) sig%spl_tck%n
        SAFE_ALLOCATE(sig%spl_tck%t, (sig%spl_tck%n))
        SAFE_ALLOCATE(sig%spl_tck%c, (sig%spl_tck%n))
        read(20,*) (sig%spl_tck%t(ii), ii=1,sig%spl_tck%n)
        read(20,*) (sig%spl_tck%c(ii), ii=1,sig%spl_tck%n)
        read(20,*) sig%spl_tck%k
        call close_file(20)
      endif

! Spline scissors, outer
    elseif(trim(keyword).eq.'spline_scissors_outer') then
      if (peinf%inode==0) then
        write(6,*) 'Reading outer spline coefficients from `spline_scissors_outer.dat`'
        write(6,*)
        ! read number of pts, knots, coefficients, degree
        call open_file(20,file='spline_scissors_outer.dat',form='formatted',status='old')
        read(20,*) sig%spl_tck_outer%n
        SAFE_ALLOCATE(sig%spl_tck_outer%t, (sig%spl_tck_outer%n))
        SAFE_ALLOCATE(sig%spl_tck_outer%c, (sig%spl_tck_outer%n))
        read(20,*) (sig%spl_tck_outer%t(ii), ii=1,sig%spl_tck_outer%n)
        read(20,*) (sig%spl_tck_outer%c(ii), ii=1,sig%spl_tck_outer%n)
        read(20,*) sig%spl_tck_outer%k
        call close_file(20)
      endif

    elseif(trim(keyword).eq.'verbosity') then
      read(line,*,err=110) peinf%verbosity
! The average potential on the faces of the unit cell in the non-periodic directions for the bands in WFN_inner
    elseif(trim(keyword).eq.'avgpot') then
      read(line,*,err=110) sig%avgpot
! The average potential on the faces of the unit cell in the non-periodic directions for the bands in WFN_outer
    elseif(trim(keyword).eq.'avgpot_outer') then
      read(line,*,err=110) sig%avgpot_outer
! Frequency dependence of the inverse dielectric matrix
    elseif(trim(keyword).eq.'frequency_dependence') then
      read(line,*,err=110) sig%freq_dep
! Is calculation only of bare exchange for EXX calculation? If so, don`t need VXC or vxc.dat
    elseif(trim(keyword).eq.'for_EXX') then
      sig%is_EXX = .TRUE.
    elseif(trim(keyword).eq.'frequency_dependence_method') then
      read(line,*,err=110) sig%freq_dep_method
    elseif(trim(keyword).eq.'occ_broadening_EXX') then
      read(line,*,err=110) sig%occ_broadening_EXX
      sig%do_smearing_EXX = .true.
    elseif(trim(keyword).eq.'cd_integration_method') then
      read(line,*,err=110) sig%cd_int_method
! Frequency grid for numerical integration of full-frequency sigma (may be different from Epsilon)
    elseif(trim(keyword).eq.'init_frequency_eval') then
      read(line,*,err=110) sig%freqevalmin
      old_grid_set = .true.
    elseif(trim(keyword).eq.'delta_frequency_eval') then
      read(line,*,err=110) sig%freqevalstep
    elseif(trim(keyword).eq.'number_frequency_eval') then
      read(line,*,err=110) sig%nfreqeval
      old_grid_set = .true.
    elseif(trim(keyword).eq.'max_frequency_eval') then
      read(line,*,err=110) max_freq_eval
    elseif(trim(keyword).eq.'frequency_grid_shift') then
      read(line,*,err=110) sig%freq_grid_shift
! Use full frequency tail for better energies
    elseif(trim(keyword).eq.'use_epsilon_remainder') then
      sig%iuseremainder = 1
! Grid of Qs in Epsilon
    elseif(trim(keyword).eq.'qgrid') then
      read(line,*,err=110) sig%qgrid(1), sig%qgrid(2), sig%qgrid(3)
! Compute the exact static CH
    elseif(trim(keyword).eq.'exact_static_ch') then
      read(line,*,err=110) sig%exact_ch
! Full CH convergence logging of all calculated bands
    elseif(trim(keyword).eq.'full_ch_conv_log') then
      read(line,*,err=110) sig%fullConvLog
! Use use_xdat to skip the computation of bare exchange
    elseif(trim(keyword).eq.'use_xdat') then
      sig%use_xdat = .true.
! Use dont_use_vxcdat to compute exchange-correlation matrix elements from VXC
    elseif(trim(keyword).eq.'dont_use_vxcdat') then
      sig%use_vxcdat = .false.
    elseif(trim(keyword).eq.'dont_use_vxc2dat') then
      sig%use_vxc2dat = .false.
    ! ZL: add tolerant setting, basically for frozen-phonon GW setup, this is a
    ! hidden flag
    elseif(trim(keyword).eq.'set_tolerant') then
      sig%tolerant_value = .true.
    ! ZL: add for arbitrary functionals
    elseif(trim(keyword).eq.'use_kihdat') then
      sig%use_kihdat = .true.
    elseif(trim(keyword).eq.'dont_use_hdf5') then
      sig%use_hdf5 = .false.
    elseif(trim(keyword).eq.'number_kpoints') then
      read(line,*,err=110) sig%nkn ! FHJ: deprecated
    elseif(trim(keyword).eq.'number_qpoints') then
      read(line,*,err=110) sig%nq ! FHJ: deprecated
    elseif(trim(keyword).eq.'dont_symmetrize') then
      sig%symmetrize = .false.
    elseif(trim(keyword).eq.'subsample') then
      sig%subsample = .true.
!#BEGIN_INTERNAL_ONLY
    ! ZL: for EP calculation, under INTERNAL now
    elseif(trim(keyword).eq.'electron_phonon') then
      sig%elph = .true.
    elseif(trim(keyword).eq.'ep_bracket') then
      read(line,*,err=110) sig%ep_bracket
    elseif(trim(keyword).eq.'phonon_q') then
      sig%nphonq = 1
      read(line,*,err=110) sig%phonq(1), sig%phonq(2), sig%phonq(3)
    elseif(trim(keyword).eq.'sigma_correction') then
      sig%sigma_correction = .true.
      sig%use_vxcdat = .false.
      sig%xfrac = 0d0
    elseif(trim(keyword).eq.'dont_check_norms') then
      peinf%check_norms=.false.
!#END_INTERNAL_ONLY
    elseif(trim(keyword).eq.'use_wfn_hdf5') then
#ifdef HDF5
      sig%wfn_hdf5 = .true.
#endif
    elseif(trim(keyword).eq.'number_sigma_pools') then
      read(line,*,err=110) peinf%npools
    elseif(trim(keyword).eq.'bare_coulomb_cutoff') then
      read(line,*,err=110) sig%ecutb
    elseif(trim(keyword).eq.'cell_average_cutoff') then
      read(line,*,err=110) sig%avgcut
    elseif(trim(keyword).eq.'screened_coulomb_cutoff') then
      read(line,*,err=110) sig%ecuts
    elseif(trim(keyword).eq.'number_bands') then
      read(line,*,err=110) sig%ntband
      call check_bounds_nbands(sig%ntband, 'number_bands')
    elseif(trim(keyword).eq.'number_core_excluded') then
      read(line,*,err=110) sig%ncore_excl
    elseif(trim(keyword).eq.'band_index_min') then
      read(line,*,err=110) nbnmin
    elseif(trim(keyword).eq.'band_index_max') then
      read(line,*,err=110) nbnmax
    elseif(trim(keyword).eq.'spin_index_min') then
      read(line,*,err=110) spinmin
    elseif(trim(keyword).eq.'spin_index_max') then
      read(line,*,err=110) spinmax
    elseif(trim(keyword).eq.'finite_difference_form') then
      read(line,*,err=110) sig%fdf
    elseif(trim(keyword).eq.'finite_difference_spacing') then
      read(line,*,err=110) sig%dw
    elseif(trim(keyword).eq.'bare_exchange_fraction') then
      read(line,*,err=110) sig%xfrac
    elseif(trim(keyword).eq.'short_range_frac_fock') then
      read(line,*,err=110) sig%coulomb_mod%short_range_frac_fock
    elseif(trim(keyword).eq.'long_range_frac_fock') then
      read(line,*,err=110) sig%coulomb_mod%long_range_frac_fock
    elseif(trim(keyword).eq.'screening_length') then
      read(line,*,err=110) sig%coulomb_mod%screening_length
    elseif(trim(keyword).eq.'write_vcoul') then
      sig%iwritecoul=1
    elseif(trim(keyword).eq.'tol_degeneracy') then
      read(line,*,err=110) sig%tol
    elseif(trim(keyword).eq.'gpp_broadening') then
      read(line,*,err=110) sig%gamma
    elseif(trim(keyword).eq.'gpp_sexcutoff') then
      read(line,*,err=110) sig%sexcut
    elseif(trim(keyword).eq.'number_diag') then
      read(line,*,err=110) sig%ndiag
    elseif(trim(keyword).eq.'number_offdiag') then
      read(line,*,err=110) sig%noffdiag
    elseif(trim(keyword).eq.'sigma_matrix') then
      read(line,*,err=110) sig%loff, sig%toff
    elseif(trim(keyword).eq.'band_occupation') then
      read(line,*,err=110) (band_occ(ii),ii=1,sig%ntband)
      occ_set = .true.
    elseif(trim(keyword).eq.'number_partial_occup') then
      read(line,*,err=110) sig%ncrit
      occ_set = .true.
    elseif(trim(keyword).eq.'fermi_level') then
      read(line,*,err=110) sig%efermi_input
    elseif(trim(keyword).eq.'fermi_level_absolute') then
      sig%rfermi=.false.
    elseif(trim(keyword).eq.'fermi_level_relative') then
      sig%rfermi=.true.
    elseif(trim(keyword).eq.'fullbz_replace') then
      sig%freplacebz=.true.
    elseif(trim(keyword).eq.'fullbz_write') then
      sig%fwritebz=.true.
    elseif(trim(keyword).eq.'degeneracy_check_override') then
      sig%degeneracy_check_override=.true.
    elseif(trim(keyword).eq.'no_symmetries_offdiagonals') then
      sig%offdiagsym=.false.
    elseif(trim(keyword).eq.'no_symmetries_q_grid') then
      sig%qgridsym=.false.
    elseif(trim(keyword).eq.'use_symmetries_q_grid') then
      sig%qgridsym=.true.
    elseif(trim(keyword).eq.'die_outside_sphere') then
      sig%die_outside_sphere=.true.
    elseif(trim(keyword).eq.'ignore_outside_sphere') then
      sig%die_outside_sphere=.false.
    elseif(trim(keyword).eq.'eqp_corrections') then
      sig%eqp_corrections=.true.
    elseif(trim(keyword).eq.'eqp_outer_corrections') then
      sig%eqp_outer_corrections=.true.
    elseif(trim(keyword).eq.'skip_averagew') then
      sig%averagew = .false.
    elseif(trim(keyword).eq.'do_sigma_subspace') then
      sig%do_sigma_subspace = .true.
    elseif(trim(keyword).eq.'sub_collective_redistr') then
      sig%sub_collective_redistr = .true.
    elseif(trim(keyword).eq.'invalid_gpp_mode') then
      read(line,*,err=110) sig%invalid_gpp_mode
    ! Algo-specific GPU input
    ! WPH: I don't like these here, because algos_inread is used to set
    !      details of the GPU acceleration, but I like moving them into
    !      algos_inread and having the algos module become dependent on
    !      siginfo's internal structure even less.
    elseif(trim(keyword).eq.'accel_mtxel_band_block_size') then
      read(line,*,err=110) sig%accel_mtxel_band_block_size
    elseif(trim(keyword).eq.'accel_gpp_band_block_size') then
      read(line,*,err=110) sig%gpp_band_block_size
    elseif(trim(keyword).eq.'accel_gpp_ig_block_size') then
      read(line,*,err=110) sig%gpp_ig_block_size
    ! end GPU input
    elseif(try_inread_truncation(trim(keyword), trim(line), sig%icutv, sig%truncval(1))) then
      ! subroutine already does the job
    elseif(try_inread_screening(trim(keyword), trim(line), sig%iscreen)) then
      ! subroutine already does the job
    else
      call scissors_inread(keyword, line, sig%scis, found)
      if(.not. found) call scissors_inread(keyword, line, sig%scis_outer, found, "_outer")
      if(.not. found) call algos_inread(keyword, line, found)
      if(.not. found) then
        write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found in sigma.inp.'
        call die(errmsg, only_root_writes = .true.)
      endif
    end if
  enddo

! for the moment subspace works only in combination with MPI and SCALAPACK
#if !defined MPI || !defined USESCALAPACK
  if(sig%do_sigma_subspace) then
    if (peinf%inode==0) then
      write(0,*)
      write(0,*) 'WARNING: Static Subspace method only works with MPI and SCALAPACK.'
      write(0,*) 'Subspace method turned off.'
      write(0,*)
    endif
    sig%do_sigma_subspace = .false.
  end if
#endif

  ! entered in Ryd, stored in eV since kp%el and kp%elda are soon converted to eV.
  sig%tol = sig%tol * RYD

  call peinfo_set_verbosity()
  if (peinf%inode==0 .and. sig%nkn>0) then
    write(0,'(/,a)') 'WARNING: the `number_kpoints` flag is deprecated. The code now'
    write(0,'(a,/)') 'automatically figures out the number of q-points from the input.'
  endif
  sig%nkn = nkpt_read
  if (peinf%inode==0 .and. sig%nq>0) then
    write(0,'(/,a)') 'WARNING: the `number_qpoints` flag is deprecated. The code now'
    write(0,'(a,/)') 'automatically figures out the number of k-points from the input.'
  endif
  sig%nq = nqpt_read
  if(sig%freq_dep==-1 .and. sig%nq<1) then
    write(0,'(/,a)') 'ERROR: Hartree-Fock calculations require a list of q-points'
    call die('The `begin qpoints` block could not be found.', only_root_writes=.true.)
  endif
  sig%nq1 = sig%nq - sig%nq0

  if(any(sig%qgrid(1:3) == 0).and.sig%freq_dep.eq.-1) then
    call die('qgrid must be specified for Hartree-Fock.', only_root_writes = .true.)
  endif
  
  if(.not. q0vec_read .and. sig%freq_dep .eq. -1) then
    call die('No q0 specified in qpoints block.', only_root_writes = .true.)
  endif
  
  SAFE_ALLOCATE(sig%kpt, (3,sig%nkn))
  sig%kpt(1:3,1:sig%nkn)=kpt_read(1:3,1:sig%nkn)
!#BEGIN_INTERNAL_ONLY
  ! ZL: EP - TODO
  if(sig%elph) then
    if (sig%nphonq .ne. 1) then
      call die('ONLY nphonq=1 implemented!', only_root_writes=.true.)
    endif
    SAFE_ALLOCATE(sig%k_phonq, (3, sig%nkn * sig%nphonq)) ! ZL: only nphonq=1 considered
    SAFE_ALLOCATE(k_phonq_read, (3))
    do ik=1,sig%nkn
      k_phonq_read(1:3) = sig%kpt(1:3, ik) + sig%phonq ! ZL: k+phonq points
!      call fold_back_zero_one(k_phonq_read(1), k_phonq_read(2), k_phonq_read(3))
      sig%k_phonq(1:3, ik) = k_phonq_read(1:3)
    enddo
    SAFE_DEALLOCATE(k_phonq_read)
  endif
!#END_INTERNAL_ONLY
 
  if (sig%nq>0) then
    if (sig%freq_dep/=-1) then
      if (peinf%inode==0) &
        write(0,'(/,a,/)') 'WARNING: ignoring the qpoints block for calculations other than HF.'
    else
      SAFE_ALLOCATE(sig%qpt, (3,sig%nq))
      sig%qpt(1:3,1:sig%nq)=qpt_read(1:3,1:sig%nq)
    endif
  endif
  
  if(sig%ndiag.ne.0) then
    SAFE_ALLOCATE(sig%diag, (sig%ndiag))
    sig%diag(1:sig%ndiag)=diag(1:sig%ndiag)
  endif
  
  if(sig%ndiag.eq.0.and.nbnmin*nbnmax.ne.0) then
    sig%ndiag= nbnmax - nbnmin + 1
    ndiag_read = sig%ndiag
    SAFE_ALLOCATE(sig%diag, (sig%ndiag))
    do ii=1,sig%ndiag
      sig%diag(ii)= ii + nbnmin - 1
    enddo
  endif
  
  if(sig%ndiag.eq.0.and.nbnmin*nbnmax.eq.0) then
    write(errmsg,'(a,3i6)') 'Incomprehensible list of energy bands: ', sig%ndiag, nbnmin, nbnmax
    call die(errmsg, only_root_writes = .true.)
  endif
  
  if(ndiag_read.lt.sig%ndiag) then
    if(peinf%inode.eq.0) then
      write(0,*) 'The number of diagonal elements found in the diag block (',ndiag_read,')'
      write(0,*) '  is smaller than the one specified by the keyword number_diag (',sig%ndiag,').'
    endif
    call die("ndiag too small", only_root_writes = .true.)
  endif
  
  ! no finite-difference evaluation unless GPP
  if(sig%freq_dep /= 1 .and. sig%freq_dep /= 3) sig%fdf=-2

  if(ndiag_read.gt.sig%ndiag) then
    if(peinf%inode.eq.0) then
      write(0,887) ndiag_read, sig%ndiag, sig%ndiag
    endif
  endif
887 format(1x,"WARNING: The number of diag elements in the diag block (",i4,") is larger",/,&
      3x,"than the one specified by the keyword number_diag (",i4,").",/,&
      3x,"The program will continue using only",i4,1x,"diagonal elements.",/)
  
  if(nbnmin*nbnmax.eq.0) then
    nbnmin=max_bands
    nbnmax=-max_bands
    do ii=1,sig%ndiag
      if (sig%diag(ii).lt.nbnmin) nbnmin=sig%diag(ii)
      if (sig%diag(ii).gt.nbnmax) nbnmax=sig%diag(ii)
    enddo
  endif
  
  sig%bmin=nbnmin
  sig%bmax=nbnmax

  if(sig%loff.lt.-2.or.(sig%loff.gt.0.and.sig%loff.lt.nbnmin).or. &
    sig%loff.gt.nbnmax.or.sig%toff.lt.-1.or.sig%toff.gt.1) then
    call die("sigma_matrix parameters out of range", only_root_writes = .true.)
  endif
  
  if(sig%loff.ne.0) then
    if(noff_read.ne.0.or.sig%noffdiag.gt.0) then
      if(peinf%inode.eq.0) then
        write(0,994)
      endif
    endif
    kk=0
    do ii=nbnmin,nbnmax
      do jj=nbnmin,nbnmax
        if (sig%toff.eq.-1.and.ii.lt.jj) cycle
        if (sig%toff.eq.1.and.ii.gt.jj) cycle
        kk=kk+1
        off1_read(kk)=ii
        off2_read(kk)=jj
        if (sig%loff.eq.-1) then
          off3_read(kk)=ii
          ! ZL: add for ep, to keep track of the unselected bra or ket
          if(sig%elph) then
            off_ep_read(kk)=jj
          endif
        else if (sig%loff.eq.-2) then
          off3_read(kk)=jj
          ! ZL: add for ep
          if(sig%elph) then
            off_ep_read(kk)=ii
          endif
        else
          off3_read(kk)=sig%loff
        endif
      enddo
    enddo
    sig%noffdiag=kk
    noff_read=sig%noffdiag
  endif
994 format(1x,"WARNING: both offdiag and sigma_matrix found in sigma.inp.",/,&
      3x,"The latter overrides the former.",/)
  
  if(noff_read.lt.sig%noffdiag) then
    if(peinf%inode.eq.0) then
      write(0,997) noff_read, sig%noffdiag
    endif
    call die("offdiag error", only_root_writes = .true.)
  endif
997 format(1x,"The number of offdiag elements in the offdiag block (",i4,") is",/,&
      3x,"smaller than the one specified by the keyword number_offdiag (",i4,").")
  
  if(noff_read.gt.sig%noffdiag) then
    if(peinf%inode.eq.0) then
      write(0,996) noff_read, sig%noffdiag, sig%noffdiag
    endif
  endif
996 format(1x,"WARNING: The number of offdiag elements in the offdiag block (",i4,") is",/,&
      3x,"larger than the one specified by the keyword number_offdiag (",i4,").",/,&
      3x,"The program will continue using only",i4,1x,"offdiag elements.",/)
  
  ! allocate even if noffdiag = 0 since sunf90 -xcheck will complain when passed unallocated to read_matrix_elements
  SAFE_ALLOCATE(sig%off1, (sig%noffdiag))
  SAFE_ALLOCATE(sig%off2, (sig%noffdiag))
  SAFE_ALLOCATE(sig%off3, (sig%noffdiag))
  ! ZL: for EP
  if(sig%elph) then
    SAFE_ALLOCATE(sig%off_ep, (sig%noffdiag))
  endif
  if(sig%noffdiag.gt.0) then
    sig%off1(1:sig%noffdiag)=off1_read(1:sig%noffdiag)
    sig%off2(1:sig%noffdiag)=off2_read(1:sig%noffdiag)
    sig%off3(1:sig%noffdiag)=off3_read(1:sig%noffdiag)
    ! ZL: for EP
    if(sig%elph) then
      sig%off_ep(1:sig%noffdiag)=off_ep_read(1:sig%noffdiag)
    endif
  endif
  
  if (peinf%inode==0 .and. (sig%ecuts>TOL_ZERO .or. sig%ecutb>TOL_ZERO)) then
    write(6,'(/1x,a/)') 'NOTE: `screened_coulomb_cutoff` and `bare_coulomb_cutoff` are now optional flags.'
  endif

  if(spinmin.eq.1.and.spinmax.eq.1) then
    sig%nspin=1
    sig%spin_index(1)=1
    sig%spin_index(2)=2 ! only used for spinor case
  elseif(spinmin.eq.2.and.spinmax.eq.2) then
    sig%nspin=1
    sig%spin_index(1)=2
  elseif(spinmin.eq.1.and.spinmax.eq.2) then
    sig%nspin=2
    sig%spin_index(1)=1
    sig%spin_index(2)=2
  else
    write(errmsg,'(a,i2,a,i2,a,i2)') 'Illegal range of spin indices from ', spinmin, ' to ', spinmax
    call die(errmsg, only_root_writes = .true.)
  endif

!------------------------------
! Build the map sig%offmap
!         sig%off1(ii) = sig%diag(sig%offmap(ii,1))
!         sig%off2(ii) = sig%diag(sig%offmap(ii,2))
!         sig%off3(ii) = sig%diag(sig%offmap(ii,3))
! This is done only if sig%noffdiag .ne. 0

  if(sig%noffdiag.gt.0) then
    SAFE_ALLOCATE(sig%offmap, (sig%noffdiag,3))
    sig%offmap(1:sig%noffdiag, 1:3) = 0
    ! ZL: build for ep
    if(sig%elph) then
      SAFE_ALLOCATE(sig%offmap_ep, (sig%noffdiag))
      sig%offmap_ep(1:sig%noffdiag) = 0
    endif
    do ii=1,sig%noffdiag
      do jj=1,sig%ndiag
        if(sig%diag(jj) .eq. sig%off1(ii)) sig%offmap(ii, 1) = jj
        if(sig%diag(jj) .eq. sig%off2(ii)) sig%offmap(ii, 2) = jj
        if(sig%diag(jj) .eq. sig%off3(ii)) sig%offmap(ii, 3) = jj
        ! ZL: for EP
        if(sig%elph) then
          if(sig%diag(jj) .eq. sig%off_ep(ii)) sig%offmap_ep(ii) = jj
        endif
      enddo
      
      if(sig%offmap(ii, 1) == 0) then
        write(errmsg,'(a,i6,a)') 'Off-diagonal matrix el. requested for band ', &
          sig%off1(ii),' but no corresponding diagonal matrix el. is requested'
        call die(errmsg, only_root_writes = .true.)
      endif
      
      if(sig%offmap(ii, 2) == 0) then
        write(errmsg,'(a,i6,a)') 'Off-diagonal matrix el. requested for band ', &
          sig%off2(ii),' but no corresponding diagonal matrix el. is requested'
        call die(errmsg, only_root_writes = .true.)
      endif
      
      if(sig%offmap(ii, 3) == 0) then
        write(errmsg,'(a,i6,a)') 'Off-diagonal matrix el. requested at energy of band ', &
          sig%off3(ii),' but no corresponding diagonal matrix el. is requested'
        call die(errmsg, only_root_writes = .true.)
      endif

      ! ZL: for EP
      if(sig%elph) then
        if(sig%offmap_ep(ii) == 0) then
          write(errmsg,'(a,i6,a)') 'Off-diagonal matrix el. requested at energy of band ', &
            sig%off_ep(ii),' but no corresponding diagonal matrix el. is requested'
          call die(errmsg, only_root_writes = .true.)
        endif
      endif
    enddo
  endif ! sig%noffdiag.gt.0
 
  if (occ_set) then
    if (peinf%inode==0) then
      write(0,'(/1x,a)') 'WARNING: keywords `number_partial_occup` and `band_occupations` are deprecated.'
      write(0,'(1x,a/)') 'BerkeleyGW now figures out these parameters automatically.'
    endif

    sig%nvband = count(band_occ==1)

    if ((sig%nvband+sig%ncrit)==0) &
      call die("There are no occupied or partially occupied bands.", only_root_writes = .true.)

    if (any(band_occ/=0 .and. band_occ/=1) .and. &
      any(band_occ(2:sig%ntband)>band_occ(1:sig%ntband-1))) then
      ! FHJ: non-equilibrium occ completely disabled. Go change your mean-field!
      call die("Non-equilibrium occupations not supported.", only_root_writes = .true.)
    endif
  endif
  
  if(sig%fullConvLog<0.or.sig%fullConvLog>1) then
    call die('Invalid full_ch_conv_log', only_root_writes = .true.)
  endif

! gsm: What frequency dependence we are using?

  if(peinf%inode.eq.0) then
    if(sig%freq_dep.eq.-1) then
      write(6,699)
    elseif(sig%freq_dep.eq.0) then
      write(6,700)
    elseif(sig%freq_dep.eq.1) then
      write(6,701)
    elseif(sig%freq_dep.eq.2) then
      write(6,702)
    elseif(sig%freq_dep.eq.3) then
      write(6,703)
    else
      call die('Need to specify frequency dependence')
    endif
  endif
699 format(1x,'Using the Hartree-Fock approximation',/)
700 format(1x,'Using the static COHSEX approximation',/)
701 format(1x,'Using the Generalized Plasmon Pole model',/)
702 format(1x,'Using the full frequency-dependent inverse dielectric matrix',/)
703 format(1x,'Using the Generalized Plasmon Pole model (GN flavor)',/)

  if(peinf%inode == 0) then
    if(peinf%npes > 1) then
      write(6,803)
    else
      write(6,805)
    endif
  endif
803 format(1x,'We are communicating via MPI',/)
805 format(1x,'We are not communicating',/)

! gsm: Do we compute the exact static CH?

  if(sig%exact_ch.eq.-1) then ! set default according to freq_dep
    if(sig%freq_dep .eq. 0) then
      sig%exact_ch = 1
    else
      sig%exact_ch = 0
    endif
  endif
  if(sig%freq_dep.ne.-1) then ! HF has no CH sum so don`t write any comments about it
    if(sig%exact_ch.eq.0) then
      if(peinf%inode.eq.0) write(6,750)
    elseif(sig%exact_ch.eq.1) then
      if(sig%freq_dep.eq.0) then
        if(peinf%inode.eq.0) write(6,751)
      else
        if(peinf%inode.eq.0) write(6,752)
      endif
    else
      call die('Illegal value for exact_static_ch')
    endif
  endif
750 format(1x,'Computing CH as a partial sum over empty bands',/)
751 format(1x,'Computing the exact static CH',/)
752 format(1x,'Computing CH as a partial sum over empty bands with the static remainder',/)

  ! JRD: What screening is present?
  if(peinf%inode.eq.0) then
    if(sig%ncrit < 0) then
      call die("number_partial_occup < 0")
    endif

    select case (sig%iscreen)
      case (SCREEN_SEMICOND)
        if(sig%ncrit > 0) then
          write(0,'(a)') "WARNING: Semiconductor screening is inappropriate for number_partial_occup > 0."
          write(0,'(a)') "If a band really crosses the Fermi level, graphene or metal screening must be used."
          ! this only makes sense if ncrit was used to handle two spins, or
          ! to tune the number of conduction and valence bands for optimal parallelization.
        endif
        write(6,'(1x,a/)') 'Running with semiconductor screening'
      case (SCREEN_GRAPHENE)
        if(sig%ncrit == 0) then
          write(0,*) "WARNING: Graphene screening usually should have number_partial_occup > 0."
        endif
        write(6,'(1x,a/)') 'Running with graphene screening'
      case (SCREEN_METAL)
        if(sig%ncrit == 0) then
          write(0,*) "WARNING: Metal screening usually should have number_partial_occup > 0."
        endif
        write(6,'(1x,a/)') 'Running with metal screening'
      case default
        call die('Unknown screening type', only_root_writes=.true.)
    endselect
  endif

  call print_truncation_summary(sig%icutv, sig%truncval(1))

  ! FHJ: set cell_average_cutoff, if not set by user
  if (sig%avgcut<0) then
    sig%avgcut = TOL_ZERO
    if (sig%icutv==TRUNC_NONE .and. sig%iscreen==SCREEN_SEMICOND) then
      sig%avgcut = 1d0/TOL_ZERO
    endif
  endif
  if (peinf%inode==0) then
    write(6,'(1x,a,es10.3e3,a/)') &
      'Cutoff for Monte-Carlo average of Coulomb potential: ', sig%avgcut, ' Ry'
  endif

  if(peinf%inode.eq.0) then
    if(sig%ncrit.gt.0) then
      write(6,951)sig%ncrit
951   format(1x,'We have partially occupied bands',/, &
        1x,'number_partial_occup (ncrit) =',i3,/)
    endif
  endif

  if(sig%qgridsym .and. .not. sig%symmetrize) then
    write(errmsg,'(a,i6,a)') 'Must use no_symmetries_q_grid flag with dont_symmetrize flag'
    call die(errmsg, only_root_writes = .true.)
  endif


  if(sig%use_xdat) then
    inquire(file='x.dat', exist=sig%use_xdat)
    if(sig%use_xdat .and. peinf%inode == 0) then
      write(6,899)
899   format(1x,"Reading bare exchange matrix elements from x.dat file."/)
    endif
  endif

  ! ZL: add compatibility check for vxcdat and kihdat
  if(sig%use_vxcdat .and. sig%use_kihdat) then
    write(6,*) "You cannot use vxc.dat and kih.dat simutaneously!"
    write(6,*) "Default is to use vxc.dat (or VXC) for standard LDA/GGA starting point."
    write(6,*) "For arbitrary functional starting point, i.e. LDA/GGA/meta-GGA/Hybrid ...,"
    write(6,*) "Use kih.dat, and set dont_use_vxcdat in sigma.inp."
    call die("vxc.dat and kih.dat cannot be used simutaneously! See output instruction.")
  endif

  if(sig%use_vxcdat) then
    inquire(file='vxc.dat', exist=sig%use_vxcdat)
    if(sig%use_vxcdat .and. peinf%inode == 0) then
      write(6,898)
898   format(1x,"Reading exchange-correlation matrix elements from vxc.dat file",/)
    endif
  endif
  if(.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX .and. .not.sig%use_kihdat) then
    ! ZL: also only if not kihdat, then search for VXC file
    inquire(file='VXC', exist=VXC_exists)
    if(.not. VXC_exists) call die("VXC file is missing.", only_root_writes = .true.)
    if(peinf%inode == 0) then
      write(6,897)
897   format(1x,"Computing exchange-correlation matrix elements from VXC file",/)
    endif
  endif
  ! This is for hybrid functional like calculation (one shot)
  if (sig%freq_dep.eq.-1 .and. ((1.0d0 - sig%coulomb_mod%long_range_frac_fock > TOL_SMALL) .or. &
    (1.0d0 - sig%coulomb_mod%short_range_frac_fock > TOL_SMALL)))  then

    sig%coul_mod_flag = .true.             
    if(sig%use_vxc2dat) then
      inquire(file='vxc2.dat', exist=sig%use_vxc2dat)
      if(sig%use_vxc2dat .and. peinf%inode == 0) then
        write(6,896)
896     format(1x,"Reading exchange-correlation matrix elements from vxc2.dat file",/)
      endif
    endif
    if(.not. sig%use_vxc2dat) then
      inquire(file='VXC2', exist=VXC_exists)
      if(.not. VXC_exists) call die("VXC2 file is missing.", only_root_writes = .true.)
      if(peinf%inode == 0) then
        write(6,895)
895     format(1x,"Computing exchange-correlation matrix elements from VXC2 file",/)
      endif
    endif
  endif
  ! ZL: add check for KIH arbitrary functionals
  if(sig%use_kihdat) then
    inquire(file='kih.dat', exist=sig%use_kihdat)
    if(sig%use_kihdat .and. peinf%inode == 0) then
      write(6,890)
890   format(1x,"Reading (Kinetic energy + Ionic potential + Hartree) matrix elements from kih.dat file",/)
    endif
  endif

  !FHJ: report on the frequency dependence method
  if (peinf%inode==0) then
    select case (sig%freq_dep)
      case (-1)
        tmpstr = 'Hartree-Fock or hybrid functional approximation'
      case (0)
        tmpstr = 'Static COHSEX approximation'
      case (1)
        tmpstr = 'Hybertsen-Louie Generalized Plasmon Pole model'
      case (3)
        tmpstr = 'Godby-Needs Generalized Plasmon Pole model'
      case (2)
        select case (sig%freq_dep_method)
          case (0)
            tmpstr = 'full frequency Real Axis Integration method'
          case (2)
            tmpstr = 'full frequency Contour Deformation method'
          case default
            write(0,'(a,i0)') 'ERROR: Got sig%freq_dep_method = ', sig%freq_dep_method
            call die('Invalid option for `frequency_dependence_method` flag.')
        endselect
      case default
        write(0,'(a,i0)') 'ERROR: Got sig%freq_dep = ', sig%freq_dep
        call die('Invalid option for `frequency_dependence` flag.')
    endselect
    write(6,'(1x,a)') 'Treating W within the '//trim(tmpstr)

    if (sig%freq_dep==2 .and. sig%freq_dep_method==2) then
      if (sig%cd_int_method/=0 .and. sig%cd_int_method/=2 .and. sig%cd_int_method/=3) then
        write(0,'(a,i0)') 'ERROR: Got sig%cd_int_method = ', sig%cd_int_method
        call die('Invalid option for `cd_integration_method` flag.')
      endif
      write(6,'(1x,a,i0,a/)') 'We`ll use an integration scheme of order ', &
        sig%cd_int_method,' in the Contour Deformation method'
    else
      write(6,*)
    endif
  endif
  sig%need_advanced = SCALARSIZE==2 .and. sig%freq_dep==2 .and. sig%freq_dep_method/=2

  !FHJ: report on the kind of grid we are using
  if (sig%freq_dep==2 .or. (sig%freq_dep==1 .and. sig%fdf==-3)) then
    if (peinf%inode==0) write(6,'(1x,a)', advance='no') 'Frequency sampling: '
    if (sig%freq_grid_shift==2) then
      sig%nfreqeval = 2*int((max_freq_eval+TOL_ZERO)/sig%freqevalstep) + 1
      if (peinf%inode==0) write(6,'(a)', advance='no') 'using a QP-centered '
    elseif (sig%freq_grid_shift==1) then
      if (peinf%inode==0) write(6,'(a)', advance='no') 'using a Fermi-energy-shifted '
    elseif (sig%freq_grid_shift==0) then
      if (peinf%inode==0) write(6,'(a)', advance='no') 'using an absolute energy '
    else
      call die('Invalid value for frequency_grid_shift.')
    endif
    if (peinf%inode==0) write(6,'(a,i0,a/)') 'grid with ', sig%nfreqeval, ' points.'
  endif
  !FHJ: Make sure the user is not using incompatible flags.
  !if ((old_grid_set .and. sig%freq_grid_shift==2).and.sig%freq_dep==2) then
  if (old_grid_set .and. sig%freq_grid_shift==2) then
    if (peinf%inode==0) then
      write(0,'(/a)') 'ERROR: flags `init_frequency_eval` and/or `number_frequency_eval` were set, but'
      write(0,'(a)') 'you didn`t change the `frequency_grid_shift`. Either:'
      write(0,'(a)') ' - set `frequency_grid_shift` to 0 or 1; or'
      write(0,'(a/)') ' - don`t set the `init_frequency_eval` and `number_frequency_eval` flags.'
    endif
    call die('Inconsistent flags with frequency_grid_shift')
  endif

  call verify_gpu_settings()

  call require_reference(REF_Deslippe2012)
  call require_reference(REF_Hybertsen1986)
  ! Truncation of the Coulomb potential: slab and write
  if (sig%icutv==TRUNC_SLAB .or. sig%icutv==TRUNC_WIRE) call require_reference(REF_IsmailBeigi2016)
  ! FF algorithm by Shishkin and Kresse, as implemented by F. Liu et al
  if (sig%freq_dep_method==1) call require_reference(REF_Liu2015)
  ! Non-uniform sampling schemes
  if (sig%subsample) call require_reference(REF_Jornada2017)
  ! Subspace for G0W0
  if (sig%do_sigma_subspace) call require_reference(REF_DelBen2019Subspace)
  ! GWPT
  if (sig%elph) call require_reference(REF_Li2019)

#ifdef MPI
  ! root lets the others go after it is done reading (see beginning of function)
  if(peinf%inode == 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  ! FHJ: Broadcast spline scissors (this must be called after the barrier)
  ! FHJ: inner splines
  call MPI_BCAST(sig%spl_tck%n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  if ( sig%spl_tck%n > 0 ) then
    if(peinf%inode/=0) then
      SAFE_ALLOCATE(sig%spl_tck%t, (sig%spl_tck%n))
      SAFE_ALLOCATE(sig%spl_tck%c, (sig%spl_tck%n))
    endif
    call MPI_BCAST(sig%spl_tck%t, sig%spl_tck%n, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(sig%spl_tck%c, sig%spl_tck%n, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(sig%spl_tck%k, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  endif
  ! FHJ: outer splines
  call MPI_BCAST(sig%spl_tck_outer%n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  if ( sig%spl_tck_outer%n > 0 ) then
    if(peinf%inode/=0) then
      SAFE_ALLOCATE(sig%spl_tck_outer%t, (sig%spl_tck_outer%n))
      SAFE_ALLOCATE(sig%spl_tck_outer%c, (sig%spl_tck_outer%n))
    endif
    call MPI_BCAST(sig%spl_tck_outer%t, sig%spl_tck_outer%n, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(sig%spl_tck_outer%c, sig%spl_tck_outer%n, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(sig%spl_tck_outer%k, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  endif
#endif

!----------------------
! Free memory

  SAFE_DEALLOCATE(kpt_read)
  SAFE_DEALLOCATE(qpt_read)
  SAFE_DEALLOCATE(off1_read)
  SAFE_DEALLOCATE(off2_read)
  SAFE_DEALLOCATE(off3_read)
  SAFE_DEALLOCATE(off_ep_read)
  SAFE_DEALLOCATE(band_occ)
  SAFE_DEALLOCATE(diag)

  POP_SUB(inread)

  return
  
110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', &
      trim(keyword), '. '
  call die(errmsg, only_root_writes = .true.)
  
112 write(errmsg,'(3a)') 'Unexpected characters were found while reading elements of the ', &
      trim(blockword),' block.'
  call die(errmsg, only_root_writes = .true.)
  
end subroutine inread


! ZL: fold back kpoint into [0, 1) range
!     kpoint MUST be in fractional coordinate
subroutine fold_back_zero_one(ka, kb, kc)
  implicit none

  real(DP), intent(inout) :: ka, kb, kc
  real(DP) :: infsmall = 1.0d-5 ! ZL: k/q/phonq-points should have higher accuracy than this
  integer :: i

  PUSH_SUB(fold_back_zero_one)
  ! ZL: TODO: get rid of this function by using macros UNIT_RANGE() and TOL_SMALL

  ! ka
  do while (ka .gt. 1.0-infsmall)
    ka = ka - 1.0
  enddo
  do while (ka .lt. 0.0-infsmall)
    ka = ka + 1.0
  enddo

  ! kb
  do while (kb .gt. 1.0-infsmall)
    kb = kb - 1.0
  enddo
  do while (kb .lt. 0.0-infsmall)
    kb = kb + 1.0
  enddo

  ! kc
  do while (kc .gt. 1.0-infsmall)
    kc = kc - 1.0
  enddo
  do while (kc .lt. 0.0-infsmall)
    kc = kc + 1.0
  enddo

  POP_SUB(fold_back_zero_one)

  return
end subroutine fold_back_zero_one

end module inread_m
