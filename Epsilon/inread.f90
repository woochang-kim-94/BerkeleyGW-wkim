!==========================================================================
!
! Routines:
!
! (1) inread()  Originally By (SIB)             Last Modified 5/9/2008 (JRD)
!
!     Reads the input file and sets various values
!
!       ecuts is the epsilon_cutoff
!       pol%qpt(1:3,1:pol%nq) is the q vector information
!
!==============================================================================

#include "f_defs.h"

module inread_m

  use algos_epsilon_m
  use global_m
  use inread_common_m
  use references_m
  use scissors_m
  implicit none

  private

  public :: &
    inread

contains

subroutine inread(pol,vwfn,cwfn)
  type (polarizability), intent(out) :: pol
  type (valence_wfns), intent(out) :: vwfn
  type (conduction_wfns), intent(out) :: cwfn
  
  character*256 :: blockword,keyword,line,errmsg
  integer :: ii,itestq,jj,nqpt_read,iostat,ifreqCounter
  
  real(DP) :: zFreq, tmpFreq, freqStep, plasmaFreq,qw_sum,div
  logical :: occ_set, found, on_steroids
  integer, allocatable :: qflags(:),band_occ(:)
  real(DP), allocatable :: omega_CC(:),qpt_read(:,:),qw_read(:)

  integer :: i
  character(len=9) :: str_intra
  logical :: add_w0_freq

  PUSH_SUB(inread)
  
#ifdef MPI
  ! Non-root nodes should wait for root to read the whole file.
  ! That way, we can be sure root gets a chance to write errors before
  ! any die call is issued by another node. Root calls MPI_Barrier below.
  if(peinf%inode /= 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

!----------------------
! Allocate arrays

  SAFE_ALLOCATE(qpt_read, (3,MAX_KPTS))
  SAFE_ALLOCATE(qw_read, (MAX_KPTS))
  SAFE_ALLOCATE(qflags, (MAX_KPTS))
  SAFE_ALLOCATE(band_occ, (MAX_BANDS))

!-----------------------------
! Set default values

  nqpt_read=0
  pol%nq=0
  pol%nq0=0
  pol%nq1=0
  pol%non_uniform=.false.
  pol%subsample=.false.
  band_occ=0
  call scissors_zero(pol%scis)
  occ_set=.false.
  vwfn%nband=0
  cwfn%nband=0
  vwfn%ncore_excl=0
  pol%freq_dep=0
  pol%freq_dep_method=2
  pol%nFreq=1
  pol%nfreq_imag=15
  pol%dInitFreq=0.0d0
  pol%dDeltaFreq=-1d0
  pol%dFreqStepIncrease=1d0
  pol%dFreqCutoff1=-1d0
  pol%dFreqCutoff2=-1d0
  pol%dBrdning=0.0d0
  pol%stop_after_qpt=-1

  pol%nSFreq=1
  pol%dInitSFreq=0.0d0
  pol%dDeltaSFreq=-1d0
  pol%dSFreqStepIncrease=0d0
  pol%dSFreqCutoff1=-1d0
  pol%dSFreqCutoff2=-1d0

  pol%fullConvLog=0
  pol%icutv=TRUNC_NONE
  pol%iwritecoul=0
  pol%truncval(:)=0.d0
  pol%ecuts=0.0d0
  pol%valueq0=0
  pol%iqexactlyzero=0
  pol%ncrit=0
#ifdef HDF5
  pol%use_hdf5 = .true.
#else
  pol%use_hdf5 = .false.
#endif
  pol%efermi_input=0.0d0
  pol%rfermi=.true.
  pol%gcomm=-1
  pol%os_opt_ffts=0
  pol%n_ffts_per_batch=20
  pol%restart=.false.
  pol%min_fftgrid=.true.
  pol%lin_denominator=0d0
  pol%nfreq_group=1
  pol%wfn_hdf5=.false.
  pol%skip_epsmat_write=.false.  !WK
  pol%skip_chimat_write=.false.  !WK
  pol%write_chi_diag=.false.     !WK
  pol%energy_window=.false.      !WK
  pol%ewin_min = 0d0  !WK
  pol%ewin_max = 0d0  !WK
  pol%skip_epsilon=.false.
  pol%skip_chi=.false.
  pol%freplacebz=.false.
  pol%fwritebz=.false.
  pol%degeneracy_check_override=.false.
  peinf%npools=0
  pol%eqp_corrections=.false.
  pol%intraband_flag=0
  pol%intraband_overlap_min=0.5d0
  pol%protection_window=0
  pol%num_cond_bands_ignore=0
  pol%patched_sampling=.false.
  pol%qgrid(:) = 0
  pol%imaginary_frequency=2.0d0*ryd
  plasmaFreq = 2.0d0*ryd
! variables for subspace truncation method in epsilon
  pol%subspace = .FALSE.
  pol%chi_eigenvalue_cutoff = 1.0d-6
  pol%neig_sub_input = -1
  pol%use_elpa = .FALSE.
  pol%keep_full_eps_static = .TRUE.
  pol%matrix_in_subspace_basis = .FALSE.
  pol%do_rpa = .false.
  pol%occ_broadening = 0.1 ! partial occupations
  pol%do_smearing = .false. ! partial occupations flag
  pol%is_insulator_RPA = .false.
! variables for nonblocking scheme
  pol%nonblocking_cyclic = .false.
  pol%dont_keep_fix_buffers = .false.
  pol%sub_collective_eigen_redistr = .false.
  add_w0_freq = .false.
  pol%tda = .false.
! GPU settings
  pol%accel_mtxel_band_block_size = 20
  pol%accel_chisum_fulloffload = .true.

! Set default values for algo acceleration
  call set_algos_to_cpu()
#if defined(OPENACC) || defined(OMP_TARGET)
  ! GPU settings
  call set_algos_to_best_available_gpu()
#endif

!----------------- Never ending loop ---------------------------------------

! Actually the loop ends when the end of the file is reached


  do while(0.eq.0)

    read(55,'(a256)',iostat=iostat) line
    if(iostat < 0) exit

! Skip comment lines

    if(len_trim(line).eq.0) cycle
    if(line(1:1).eq.'#') cycle

! Determine keyword:

    keyword=line(1:scan(line," ")-1)
    line=adjustl(line(scan(line," ")+1:256))

! SIB: If we have a 'begin', then in scan the information up to 'end'
! For now, only 'qpoints' is recognized
! DVF: I need to put the do_rpa first because it changes how the q-points 
! are listed, i.e. they need weights
    if(trim(keyword).eq.'do_rpa') then
      pol%do_rpa = .true.
    elseif(trim(keyword).eq.'begin') then
      blockword=line(1:scan(line," ")-1)
      ii=0
      do while(trim(line).ne.'end')
        read(55,'(a256)',end=105) line
        if(trim(line).ne.'end') then
          ii=ii+1
          call check_bounds_nkq(ii, 'q', 'begin qpoints')
          if(trim(blockword).eq.'qpoints') then
            if(pol%do_rpa) then
              read(line,*,iostat=iostat) (qpt_read(jj,ii),jj=1,3), div,itestq,qw_read(ii)
            else
              read(line,*,iostat=iostat) (qpt_read(jj,ii),jj=1,3), div,itestq
            endif
            if(iostat /= 0) then
              write(errmsg,'(3a)') 'Unexpected characters were found while reading elements of the ', &
                trim(blockword),' block.'
              call die(errmsg, only_root_writes = .true.)
            endif
            if (itestq>0) then
              if (itestq/=1 .and. itestq/=2) &
                call die("Illegal value for last column in qpoints. May only be -1, 0, 1, 2.", &
                  only_root_writes=.true.)
              pol%nq0 = pol%nq0 + 1
              if (pol%valueq0/=0 .and. itestq/=pol%valueq0) &
                call die("All q->0 points must have the same value for last column.", only_root_writes=.true.)
              pol%valueq0=itestq
              if (all(abs(qpt_read(:,ii))<TOL_Zero)) then
                pol%iqexactlyzero = 1
                if (peinf%inode==0) &
                  write(6,'(1x,a/)') 'Doing q exactly zero'
              endif
            endif
            qpt_read(1:3,ii)=qpt_read(1:3,ii)/div
            qflags(ii) = itestq
          else
            write(errmsg,'(3a)') 'Unexpected blockword ', trim(blockword), ' was found in epsilon.inp.'
            call die(errmsg, only_root_writes = .true.)
          end if
        end if
      end do
      if(trim(blockword).eq.'qpoints') then
        nqpt_read=ii
      endif

! SIB: Other keywords than 'begin'

    elseif(trim(keyword).eq.'verbosity') then
      read(line,*,err=110) peinf%verbosity
    elseif(trim(keyword).eq.'frequency_dependence') then
      read(line,*,err=110) pol%freq_dep
    elseif(trim(keyword).eq.'frequency_dependence_method') then
      read(line,*,err=110) pol%freq_dep_method
    elseif(trim(keyword).eq.'init_frequency') then
      read(line,*,err=110) pol%dInitFreq
    elseif(trim(keyword).eq.'delta_frequency') then
      read(line,*,err=110) pol%dDeltaFreq
    elseif(trim(keyword).eq.'delta_frequency_step') then
      read(line,*,err=110) pol%dFreqStepIncrease
    elseif(trim(keyword).eq.'frequency_low_cutoff') then
      read(line,*,err=110) pol%dFreqCutoff1
    elseif(trim(keyword).eq.'frequency_high_cutoff') then
      read(line,*,err=110) pol%dFreqCutoff2
    elseif(trim(keyword).eq.'number_imaginary_freqs') then
      read(line,*,err=110) pol%nfreq_imag
    elseif(trim(keyword).eq.'plasma_freq') then
      read(line,*,err=110) plasmaFreq
    elseif(trim(keyword).eq.'broadening') then
      read(line,*,err=110) pol%dBrdning
    elseif(trim(keyword).eq.'init_sfrequency') then
      read(line,*,err=110) pol%dInitSFreq
    elseif(trim(keyword).eq.'delta_sfrequency') then
      read(line,*,err=110) pol%dDeltaSFreq
    elseif(trim(keyword).eq.'delta_sfrequency_step') then
      read(line,*,err=110) pol%dSFreqStepIncrease
    elseif(trim(keyword).eq.'sfrequency_low_cutoff') then
      read(line,*,err=110) pol%dSFreqCutoff1
    elseif(trim(keyword).eq.'sfrequency_high_cutoff') then
      read(line,*,err=110) pol%dSFreqCutoff2
    elseif(trim(keyword).eq.'full_chi_conv_log') then
      read(line,*,err=110) pol%fullConvLog
    elseif(trim(keyword).eq.'number_qpoints') then
      read(line,*,err=110) pol%nq ! FHJ: deprecated
    elseif(trim(keyword).eq.'qgrid') then
      read(line,*,err=110) pol%qgrid(1:3)
    elseif(trim(keyword).eq.'number_valence_pools') then
      read(line,*,err=110) peinf%npools
    elseif(trim(keyword).eq.'skip_epsilon') then
      pol%skip_epsilon=.true.
    elseif(trim(keyword).eq.'skip_chi') then
      pol%skip_chi=.true.
    !WK
    elseif(trim(keyword).eq.'skip_epsmat_write') then
      pol%skip_epsmat_write=.true.
    elseif(trim(keyword).eq.'skip_chimat_write') then
      pol%skip_chimat_write=.true.
    elseif(trim(keyword).eq.'write_chi_diag') then
      pol%write_chi_diag=.true.
    elseif(trim(keyword).eq.'energy_window') then
      pol%energy_window=.true.
    elseif(trim(keyword).eq.'ewin_min') then
      read(line,*,err=110) pol%ewin_min
    elseif(trim(keyword).eq.'ewin_max') then
      read(line,*,err=110) pol%ewin_max
    !WK
    elseif(trim(keyword).eq.'gcomm_elements') then
      pol%gcomm=0
    elseif(trim(keyword).eq.'gcomm_matrix') then
      pol%gcomm=-1
    elseif(trim(keyword).eq.'dont_use_hdf5') then
      pol%use_hdf5 = .false.
    elseif(trim(keyword).eq.'no_min_fftgrid') then
      pol%min_fftgrid = .false.
    elseif(trim(keyword).eq.'imaginary_frequency') then
      read(line,*,err=110) pol%imaginary_frequency
    elseif(trim(keyword).eq.'nfreq_group') then
      read(line,*,err=110) pol%nfreq_group
    elseif(trim(keyword).eq.'subsample') then
      pol%subsample = .true.
!#BEGIN_INTERNAL_ONLY
    elseif(trim(keyword).eq.'patched_sampling') then
      pol%patched_sampling = .true.
    elseif(trim(keyword).eq.'intraband_flag') then
      read(line,*,err=110) pol%intraband_flag
    elseif(trim(keyword).eq.'intraband_overlap_min') then
      read(line,*,err=110) pol%intraband_overlap_min
    elseif(trim(keyword).eq.'protection_window') then
      read(line,*,err=110) pol%protection_window
    elseif(trim(keyword).eq.'num_cond_bands_ignore') then
      read(line,*,err=110) pol%num_cond_bands_ignore
    elseif(trim(keyword).eq.'non_uniform') then
      pol%non_uniform = .true.
    elseif(trim(keyword).eq.'os_opt_ffts') then
      read(line,*,err=110) pol%os_opt_ffts
    elseif(trim(keyword).eq.'n_ffts_per_batch') then
      read(line,*,err=110) pol%n_ffts_per_batch
!#END_INTERNAL_ONLY
    elseif(trim(keyword).eq.'os_hdf5') then !FHJ: Deprecated
      if (peinf%inode==0) write(0,'(/a/)') &
        "WARNING: the flag 'os_hdf5' is deprecated. Use 'wfn_hdf5' instead."
#ifdef HDF5
      pol%wfn_hdf5 = .true.
#endif
    elseif(trim(keyword).eq.'use_wfn_hdf5') then
#ifdef HDF5
      pol%wfn_hdf5 = .true.
#endif
    elseif(trim(keyword).eq.'restart') then
      pol%restart = .true.
    elseif(trim(keyword).eq.'stop_after_qpt') then
      read(line,*,err=110) pol%stop_after_qpt
    elseif(trim(keyword).eq.'write_vcoul') then
      pol%iwritecoul=1
    elseif(trim(keyword).eq.'epsilon_cutoff') then
      read(line,*,err=110) pol%ecuts
    elseif(trim(keyword).eq.'number_bands') then
      read(line,*,err=110) cwfn%nband
      call check_bounds_nbands(cwfn%nband, 'number_bands')
    elseif(trim(keyword).eq.'band_occupation') then
      read(line,*,err=110) (band_occ(ii),ii=1,cwfn%nband)
      occ_set = .true.
    elseif(trim(keyword).eq.'number_partial_occup') then
      read(line,*,err=110) pol%ncrit
      occ_set = .true.
    elseif(trim(keyword).eq.'number_core_excluded') then
      read(line,*,err=110) vwfn%ncore_excl
    elseif(trim(keyword).eq.'fermi_level') then
      read(line,*,err=110) pol%efermi_input
    elseif(trim(keyword).eq.'fermi_level_absolute') then
      pol%rfermi=.false.
    elseif(trim(keyword).eq.'fermi_level_relative') then
      pol%rfermi=.true.
    elseif(trim(keyword).eq.'occ_broadening') then
      read(line,*,err=110) pol%occ_broadening
      pol%do_smearing = .true.
    elseif(trim(keyword).eq.'RPA_insulator') then
      pol%is_insulator_RPA = .true.
    elseif(trim(keyword).eq.'fullbz_replace') then
      pol%freplacebz=.true.
    elseif(trim(keyword).eq.'fullbz_write') then
      pol%fwritebz=.true.
!#BEGIN_INTERNAL_ONLY
    elseif(trim(keyword).eq.'lin_denominator') then
      read(line,*,err=110) pol%lin_denominator
    elseif(trim(keyword).eq.'dont_check_norms') then
      peinf%check_norms=.false.
    elseif(trim(keyword).eq.'use_tda') then
      pol%tda=.true.
!#END_INTERNAL_ONLY
    elseif(trim(keyword).eq.'degeneracy_check_override') then
      pol%degeneracy_check_override=.true.
    elseif(trim(keyword).eq.'eqp_corrections') then
      pol%eqp_corrections=.true.
!#BEGIN_INTERNAL_ONLY
    ! Algo-specific GPU input
    elseif(trim(keyword).eq.'accel_mtxel_band_block_size') then
      read(line,*,err=110) pol%accel_mtxel_band_block_size
    elseif(trim(keyword).eq.'accel_chi_summation_no_full_offload') then
      pol%accel_chisum_fulloffload = .false.
!#END_INTERNAL_ONLY
! subspace truncation
    elseif(trim(keyword).eq.'chi_eigenvalue_cutoff' .or. trim(keyword).eq.'eps_trunc_eigen') then
      read(line,*,err=110) pol%chi_eigenvalue_cutoff
      pol%subspace =.TRUE.
    elseif(trim(keyword).eq.'nbasis_subspace') then
      read(line,*,err=110) pol%neig_sub_input
      pol%subspace =.TRUE.
    elseif(trim(keyword).eq.'subspace_dont_keep_full_eps_omega0') then
      pol%keep_full_eps_static =.FALSE.
    elseif(trim(keyword).eq.'subspace_use_elpa') then
#ifdef USEELPA
      pol%use_elpa = .TRUE.
#else
      pol%use_elpa =.FALSE.
#endif
    elseif(trim(keyword).eq.'write_subspace_epsinv') then
      pol%matrix_in_subspace_basis = .TRUE.
    elseif(trim(keyword).eq.'do_rpa') then
      pol%do_rpa = .true.
    elseif(trim(keyword).eq.'comm_nonblocking_cyclic') then
      pol%nonblocking_cyclic = .true.
    elseif(trim(keyword).eq.'add_w0_freq') then
      add_w0_freq = .true.
    elseif(trim(keyword).eq.'dont_keep_fix_buffers') then
      pol%dont_keep_fix_buffers = .true.
    elseif(trim(keyword).eq.'sub_collective_eigen_redistr') then
      pol%sub_collective_eigen_redistr = .true.
    elseif(trim(keyword).eq.'read_chi_maxmem_per_block') then
      read(line,*,err=110) peinf%maxmem_read_h5_mat
    elseif(try_inread_truncation(trim(keyword), trim(line), pol%icutv, pol%truncval(1))) then
      ! subroutine already does the job
    else
      call scissors_inread(keyword, line, pol%scis, found)
      if(.not. found) call algos_inread(keyword, line, found)
      if(.not. found) then
        write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found in epsilon.inp.'
        call die(errmsg, only_root_writes = .true.)
      endif
    end if
  enddo

! End of Big Input Loop

! for the moment subspace works only in combination with MPI and SCALAPACK
#ifndef MPI
#ifndef USESCALAPACK
  if (pol%subspace) then
    if (peinf%inode==0) then
      write(0,*)
      write(0,*) 'WARNING: Static Subspace method only works with MPI and SCALAPACK.'
      write(0,*) 'Subspace method turned off.'
      write(0,*)
    endif
    pol%subspace = .false.
    pol%keep_full_eps_static = .false.
    pol%matrix_in_subspace_basis = .false.
  endif
#endif
#endif

  call peinfo_set_verbosity()
  if ((pol%freq_dep==2).and.(pol%freq_dep_method==1)) then
    if (pol%gcomm==-1) then
      pol%gcomm=0
      if (peinf%inode==0) then
        write(0,*)
        write(0,*) 'WARNING: Spectral method for polarizability does not work with gcomm_matrix.'
        write(0,*) 'Changing communication method to gcomm_elements.'
        write(0,*)
      endif
    endif
  endif

  if (pol%freq_dep==2) then
    ! Default settings for full-frequency calculations
    if (pol%freq_dep_method==3) then
      if (pol%dBrdning<0d0) pol%dBrdning = 0.25d0
      if (pol%dFreqCutoff1<0d0) pol%dFreqCutoff1 = 10d0
    else
      if (pol%dBrdning<0d0) pol%dBrdning = 0.1d0
      if (pol%dFreqCutoff1<0d0) pol%dFreqCutoff1 = 200d0
      if (pol%dFreqCutoff2<0d0) pol%dFreqCutoff2 = 4*pol%dFreqCutoff1
    endif
    if (pol%dDeltaFreq<0d0) pol%dDeltaFreq = pol%dBrdning
    if (pol%freq_dep_method==1) then
      if (pol%dDeltaSFreq<0d0) pol%dDeltaSFreq = pol%dDeltaFreq
      if (pol%dSFreqCutoff1<0d0) pol%dSFreqCutoff1 = pol%dFreqCutoff1
      if (pol%dSFreqCutoff2<0d0) pol%dSFreqCutoff2 = pol%dSFreqCutoff1
    endif
    if (pol%freq_dep_method/=2) then
      pol%nfreq_imag = 0
    endif

    if (peinf%inode==0) then
      write(6,'(1x,a)') 'Parameters for the full-frequency calculation:'
      write(6,'(1x,a,i0)') '- frequency_dependence_method: ', pol%freq_dep_method
      write(6,'(1x,a,f0.6)') '- init_frequency: ', pol%dInitFreq
      write(6,'(1x,a,f0.6)') '- broadening: ', pol%dBrdning
      write(6,'(1x,a,f0.6)') '- delta_frequency: ', pol%dDeltaFreq
      if (pol%freq_dep_method/=2) then
        write(6,'(1x,a,f0.6)') '- delta_frequency_step: ', pol%dFreqStepIncrease
      endif
      write(6,'(1x,a,f0.6)') '- frequency_low_cutoff: ', pol%dFreqCutoff1
      if (pol%freq_dep_method/=2) then
        write(6,'(1x,a,f0.6)') '- frequency_high_cutoff: ', pol%dFreqCutoff2
      else
        write(6,'(1x,a,i0)') '- number_imaginary_freqs: ', pol%nfreq_imag
      endif
      if (pol%freq_dep_method==1) then
        write(6,'(1x,a,f0.6)') '- init_sfrequency: ', pol%dInitSFreq
        write(6,'(1x,a,f0.6)') '- delta_sfrequency: ', pol%dDeltaSFreq
        write(6,'(1x,a,f0.6)') '- delta_sfrequency_step: ', pol%dSFreqStepIncrease
        write(6,'(1x,a,f0.6)') '- sfrequency_low_cutoff: ', pol%dSFreqCutoff1
        write(6,'(1x,a,f0.6)') '- sfrequency_high_cutoff: ', pol%dSFreqCutoff2
      endif
      !XXXXX
      IF(pol%subspace) THEN
        write(6,'(1x,a,e7.1)') '- using subspace truncation with EPS eigenvalues: ', &
                               pol%chi_eigenvalue_cutoff
      END IF
      !XXXXX
      write(6,'()')
    endif
  else
    pol%nfreq_imag = 0
  endif
  if(peinf%inode==0) then
    if (pol%tda) then
      write(6,'(1x,a)') 'WARNING: Polarizability is computed to be compatible with TDA'
      write(6,'(1x,a)') '         Only use if you know what you are doing.'
    endif
  endif

!-----------------------------------------------------------------------

!------------ This is where we end up if the file was read successfully -----------

  if (.not.pol%use_hdf5) pol%restart = .false.
  if (peinf%inode==0) then
    if (pol%restart) then
      write(6,'(1x,a,/)') 'We will restart the calculation from the last finished q-point.'
    else
      write(6,'(1x,a,/)') 'We will start a new calculation from scratch.'
    endif
  endif

  if(pol%nfreq_group .gt. peinf%npes) then
    write(0,'(/1x,a)') 'WARNING: Number of frequency groups cannot exceed number of processors'
    write(0,'(/,1x,2(a,i5),/)') 'Resetting nfreq_group',pol%nfreq_group,' to number of processors', peinf%npes
    pol%nfreq_group=peinf%npes
  endif

!#BEGIN_INTERNAL_ONLY
  ! FHJ: TODO - os_opt_ffts=-1 should be implemented!
  !pol%os_opt_ffts=1

  on_steroids=pol%os_opt_ffts/=0

  if (pol%os_opt_ffts<-1.or.pol%os_opt_ffts>2) &
    call die('Invalid option for os_opt_ffts.', only_root_writes = .true.)
  if (pol%os_opt_ffts==-1) &
    call die('Selected option for os_opt_ffts is not currently supported.', only_root_writes = .true.)

  if (peinf%inode==0) then
    write(6,'(1x,a)') 'Experimental Features:'
    write(6,'(1x,a,i0)') ' - FFTs Optimization Level: ', pol%os_opt_ffts
    write(6,'(1x,a,l1)') ' - Use HDF5 Parallel IO: ', pol%wfn_hdf5
#ifndef HDF5
    if(pol%wfn_hdf5)  call die("Must build code with HDF5 to use wfn_hdf5 option.", only_root_writes = .true.)
#endif    
    if (on_steroids) then
      write(0,'(/a/)') 'WARNING: Epsilon is running with experimental features!'
    endif
    write(6,'()')
  endif
!#END_INTERNAL_ONLY

  if ((pol%freq_dep .eq. 2) .and. abs(pol%dDeltaFreq) .gt. TOL_Zero) then

! JRD: Only use low_freq_cutoff for contour deformation
    if (pol%freq_dep_method .eq. 2) pol%dFreqCutoff2 = pol%dFreqCutoff1

    tmpFreq = pol%dInitFreq
    ifreqCounter = 0
    if (add_w0_freq) iFreqCounter = iFreqCounter + 1
    freqStep = pol%dDeltaFreq
    do while (tmpFreq .le. pol%dFreqCutoff2)
      ifreqCounter = ifreqCounter+1
      if (tmpFreq .lt. pol%dFreqCutoff1) then
        tmpFreq=tmpFreq+pol%dDeltaFreq
      else 
        freqstep = freqstep + pol%dFreqStepIncrease
        tmpFreq=tmpFreq+freqStep
      endif
    enddo

    if (pol%freq_dep_method .eq. 2) iFreqCounter = iFreqCounter + pol%nfreq_imag
    
    pol%nFreq = iFreqCounter

! This condition plays nicely with the condition above when nfreq_group .gt. npes, i.e. the conditions ensure
! the number of processors will always be re-set to a legitimate/sensible number 
    if(pol%nfreq_group .gt. pol%nFreq) then
      write(0,'(/1x,a)') 'WARNING: Number of frequency groups cannot exceed number of frequencies computed'
      write(0,'(/,1x,2(a,i5),/)') 'Resetting nfreq_group',pol%nfreq_group,'to number of frequencies computed', pol%nfreq
      pol%nfreq_group=peinf%npes
    endif

    !XXXXX
    if(pol%subspace) then
      if(pol%nfreq_group .lt. peinf%npes) then
        ! this is a workaround for the subspace method in order to use all proc at freq=0
        pol%nfreq_group = 1
      end if
    end if
    !XXXXX

    SAFE_ALLOCATE(pol%dFreqGrid,(pol%nFreq))
    SAFE_ALLOCATE(pol%dFreqBrd,(pol%nFreq))
    
    ifreqCounter = 0
    if (add_w0_freq) then
      iFreqCounter = iFreqCounter + 1
      pol%dFreqGrid(iFreqCounter) = 0d0
      pol%dFreqBrd(iFreqCounter) = 1.0d-4 * (0.0,1.0)
    endif
    tmpFreq = pol%dInitFreq
    freqStep = pol%dDeltaFreq
    do while (tmpFreq .le. pol%dFreqCutoff2)
      ifreqCounter = ifreqCounter+1
      pol%dFreqGrid(iFreqCounter)=tmpFreq
      pol%dFreqBrd(iFreqCounter)=pol%dBrdning*(0.0,1.0)
      
      if (tmpFreq .lt. pol%dFreqCutoff1) then
        tmpFreq=tmpFreq+pol%dDeltaFreq
      else 
        freqstep = freqstep + pol%dFreqStepIncrease
        tmpFreq=tmpFreq+freqStep
      endif

! JRD Dumb Debugging
! JRD XXX We should however probably right out the frequency grid somewhere
    !if (peinf%inode .eq. 0) write(6,*) iFreqCounter, pol%dFreqGrid(iFreqCounter), Pol%dFreqBrd(iFreqCounter)

    enddo

    if (pol%freq_dep_method .eq. 2) then
      do i = 1, pol%nfreq_imag
        iFreqCounter = iFreqCounter + 1
        zFreq = (1D0/pol%nfreq_imag)*DBLE(i-1)
        tmpFreq = -1D0 * plasmaFreq * (zFreq / (zFreq - 1D0))
        pol%dFreqGrid(iFreqCounter)=0D0
        pol%dFreqBrd(iFreqCounter)=tmpFreq*(0.0,1.0)
! JRD XXX We should write out frequency grid somewhere
        !if (peinf%inode .eq. 0) write(6,*) iFreqCounter, pol%dFreqGrid(iFreqCounter), Pol%dFreqBrd(iFreqCounter)
      enddo

      if(pol%do_rpa) then
        ! MDB these are the parameters for the CC grid
        SAFE_ALLOCATE(pol%rpa_freq_grid, (pol%nfreq_imag)) 
        SAFE_ALLOCATE(omega_CC, (pol%nfreq_imag))
        pol%rpa_freq_grid = 0.0D+00
        omega_CC = 0.0D+00
        do i = 1, pol%nfreq_imag
          omega_CC(i) = i * PI_D * 0.5D+00 / pol%nfreq_imag
        end do
        !
        do i = 1, pol%nfreq_imag - 1
          pol%rpa_freq_grid(i) = PI_D / (pol%nfreq_imag * (SIN(omega_CC(i))**2))
        end do
        pol%rpa_freq_grid(pol%nfreq_imag) = PI_D * 0.5D+00 / (pol%nfreq_imag*(SIN(omega_CC(pol%nfreq_imag))**2))
        do i = 1, pol%nfreq_imag
          omega_CC(i) = 1.0D+00 / TAN(omega_CC(i))
        end do
       ! omega_CC = omega_CC * ryd
        omega_CC = omega_CC * ryd+pol%dBrdning 
        ! copy grid
        do i = 1, pol%nfreq_imag
          pol%dFreqBrd(i + (pol%Nfreq - pol%nfreq_imag)) = omega_CC(pol%nfreq_imag-i+1) * (0.0,1.0)
        end do
        ! CC grid
      endif
    endif
  
    if(pol%freq_dep_method .eq. 1) then

      tmpFreq = pol%dInitSFreq
      ifreqCounter = 0
      freqStep = pol%dDeltaSFreq
      do while (tmpFreq .le. pol%dSFreqCutoff2)
        ifreqCounter = ifreqCounter+1
        if (tmpFreq .lt. pol%dSFreqCutoff1) then
          tmpFreq=tmpFreq+pol%dDeltaSFreq
        else
          freqstep = freqstep + pol%dSFreqStepIncrease
          tmpFreq=tmpFreq+freqStep
        endif
      enddo

      pol%nSFreq = ifreqCounter
      if (peinf%inode.eq.0) then 
        print*, 'pol%nSFreq',pol%nSFreq
      endif

      SAFE_ALLOCATE(pol%dSFreqGrid,(pol%nSFreq))

      tmpFreq = pol%dInitSFreq
      ifreqCounter = 0
      freqStep = pol%dDeltaSFreq
      do while (tmpFreq .le. pol%dSFreqCutoff2)
        ifreqCounter = ifreqCounter+1
        pol%dSFreqGrid(ifreqCounter)=tmpFreq

        if (tmpFreq .lt. pol%dSFreqCutoff1) then
          tmpFreq=tmpFreq+pol%dDeltaSFreq
        else
          freqstep = freqstep + pol%dSFreqStepIncrease
          tmpFreq=tmpFreq+freqStep
        endif
      enddo
    endif
  else if ((pol%freq_dep .eq. 2) .and. abs(pol%dDeltaFreq) .le. TOL_Zero) then
    call die("Illegal value for Delta Frequency in full frequency calculation")
  else if (pol%freq_dep .eq. 3) then
    pol%nFreq=2
    SAFE_ALLOCATE(pol%dFreqGrid,(pol%nFreq))
    SAFE_ALLOCATE(pol%dFreqBrd,(pol%nFreq))  
    pol%dFreqGrid = 0D0
    pol%dFreqBrd(1) = (0D0,0D0)      ! In Godby-Needs, the first frequency is always zero
    pol%dFreqBrd(2) = (0D0,1D0) * pol%imaginary_frequency
  else
    pol%nFreq=1
    SAFE_ALLOCATE(pol%dFreqGrid,(pol%nFreq))
    SAFE_ALLOCATE(pol%dFreqBrd,(pol%nFreq))  
    pol%dFreqGrid = 0D0
    pol%dFreqBrd = 0D0
  endif
  if (pol%freq_dep==3) then
    pol%nfreq_imag=1
  elseif (pol%freq_dep/=2) then
    pol%nfreq_imag=0
  endif

  if (peinf%inode==0 .and. pol%nq>0) then
    write(0,'(/,a)') 'WARNING: the `number_qpoints` flag is deprecated. The code now'
    write(0,'(a,/)') 'automatically figures out the number of q-points from the input.'
  endif
  pol%nq = nqpt_read

  pol%nq1 = pol%nq - pol%nq0
  if (.not.pol%subsample.and.pol%nq0>1) then
    call die('Can only have one q->0 point', only_root_writes=.true.)
  endif
  
  ! SIB: allocate polarizability q-points array pol%qpt(3,1:pol:nq)
  ! and copy what was read from qpt_read into it
  SAFE_ALLOCATE(pol%qpt, (3,pol%nq))
  pol%qpt(1:3,1:pol%nq)=qpt_read(1:3,1:pol%nq)

  SAFE_ALLOCATE(pol%qflags, (pol%nq))
  pol%qflags(1:pol%nq)=qflags(1:pol%nq)
  
  ! MDB: weights for RPA and associated array
  if(pol%do_rpa) then
    SAFE_ALLOCATE(pol%qw_rpa, (pol%nq))
    pol%qw_rpa(1:pol%nq) = qw_read(1:pol%nq)
    qw_sum = SUM(pol%qw_rpa(1:pol%nq))
    pol%qw_rpa = pol%qw_rpa / qw_sum
    SAFE_ALLOCATE(pol%E_rpa_qp, (pol%nq))
    pol%E_rpa_qp = 0.0D+00
  endif

  ! SIB: check for bad input
  if(abs(pol%ecuts).lt.TOL_Zero) then
    call die("The epsilon_cutoff keyword could not be found.", only_root_writes = .true.)
  endif

  if (occ_set) then
    if (peinf%inode==0) then
      write(0,'(/1x,a)') 'WARNING: keywords `number_partial_occup` and `band_occupations` are deprecated.'
      write(0,'(1x,a/)') 'BerkeleyGW now figures out these parameters automatically.'
    endif

    vwfn%nband = count(band_occ==1)
    if(vwfn%nband == 0 .and. pol%ncrit == 0) &
      call die("There are no occupied or partially occupied bands.", only_root_writes = .true.)
    if(vwfn%nband == cwfn%nband) &
      call die("There are no unoccupied bands.", only_root_writes = .true.)

    if (any(band_occ/=0 .and. band_occ/=1) .and. &
      any(band_occ(2:cwfn%nband)>band_occ(1:cwfn%nband-1))) then
      ! FHJ: non-equilibrium occ completely disabled. Go change your mean-field!
      call die("Non-equilibrium occupations not supported.", only_root_writes = .true.)
    endif
  endif

  if(pol%fullConvLog.lt.-1.or.pol%fullConvLog.gt.2) then
    call die('Invalid full_chi_conv_log', only_root_writes = .true.)
  endif
  
  if(pol%skip_epsilon .and. pol%skip_chi) then
    call die('Cannot skip_epsilon and skip_chi', only_root_writes = .true.)
  endif

  ! FHJ: will we need to read at least one q-point from WFNq?
  pol%need_WFNq = (pol%nq0>0.and.pol%valueq0==1.and.pol%iqexactlyzero==0)&
    .or.pol%patched_sampling.or.any(pol%qflags==-1)

!#BEGIN_INTERNAL_ONLY
  if (pol%patched_sampling.and.pol%non_uniform) then
    call die('Cannot do patched sampling and non-uniform sampling at the same time!', &
      only_root_writes=.true.)
  endif
#ifndef USEVORO
  if (pol%non_uniform) then
    call die('non_uniform flag requires compilation with -DUSEVORO flag.', &
      only_root_writes=.true.)
  endif
#endif
!#END_INTERNAL_ONLY

  if (peinf%inode==0) then
    write(6,*)
    if (pol%subsample) write(6,'(1x,a)') 'Using the subsample method to capture several q->0 points.'
    if (pol%non_uniform) then
      write(6,'(1x,a)') 'We`ll perform a non-uniform sampling of the BZ using a Voronoi decomposition scheme.'
    elseif (pol%patched_sampling) then
      write(6,'(1x,a)') 'We`ll perform a uniform sampling only on a small patch of the BZ.'
    else
      write(6,'(1x,a)') 'We`ll perform a uniform sampling of the full BZ.'
    endif
    write(6,*)
    if (pol%non_uniform) then
      call open_file(666, file='kweights.dat', form='formatted', status='replace')
      call close_file(666)
    endif
  endif

! gsm: What frequency dependence we are using?

  if(peinf%inode.eq.0) then
    if(pol%freq_dep.eq.0) then
      write(6,700)
      write(7,700)
    elseif((pol%freq_dep.eq.2).and.(pol%freq_dep_method.eq.0)) then
      write(6,702)
      write(7,702)
    elseif((pol%freq_dep.eq.2).and.(pol%freq_dep_method.eq.1)) then
      write(6,703)
      write(7,703)
    elseif((pol%freq_dep.eq.2).and.(pol%freq_dep_method.eq.2)) then
      write(6,705)
      write(7,705)
    elseif(pol%freq_dep.eq.3) then
      write(6,704)
      write(7,704)
    else
      call die('Need to specify frequency dependence', only_root_writes = .true.)
    endif
  endif
700 format(1x,'Computing the static inverse dielectric matrix',/)
702 format(1x,'Computing the full frequency-dependent inverse dielectric matrix (Adler-Wiser formulation)',/)
703 format(1x,'Computing the full frequency-dependent inverse dielectric matrix (spectral method)',/)
705 format(1x,'Computing the full frequency-dependent inverse dielectric matrix (Contour-Deformation formulation)',/)
704 format(1x,'Computing the inverse dielectric matrix for two purely imaginary frequencies (Godby-Needs formulation)',/)

! JRD: What Communication Scheme are we using?
! We should make a better default choice based on 
! whether nv*nc*nk(reduced) > nmtx*nfreq
  
  if (peinf%inode.eq.0) then
    if (pol%gcomm .eq. -1) then
      write(6,801)
      write(7,801)
      if( pol%nonblocking_cyclic .and. (pol%subspace .or. (pol%freq_dep.eq.0)) ) then
        write(6,806)
        write(7,806)
        if ( pol%subspace ) then
          if( pol%dont_keep_fix_buffers ) then
            write(6,807)
            write(7,807)
          else
            write(6,808)
            write(7,808)
          end if
        end if
      end if
      if( pol%subspace ) then
        if( pol%sub_collective_eigen_redistr ) then
          write(6,809)
          write(7,809)
        else
          write(6,810)
          write(7,810)
        end if
      end if
    else
      write(6,802)
      write(7,802)
    endif
  endif
801 format(1x,'We are using matrix communication scheme',/)
802 format(1x,'We are using element communication scheme',/)
806 format(1x,'We are using non-blocking cyclic communication scheme')
807 format(1x,'Communication buffers will NOT be statically allocated',/)
808 format(1x,'Communication buffers will be statically allocated',/)
809 format(1x,'Subspace redistribute eigenvectors with collective communication',/)
810 format(1x,'Subspace redistribute eigenvectors with point to point communication',/)

  if(peinf%inode == 0) then
    if(peinf%npes > 1) then
      write(6,803)
      write(7,803)
    else
      write(6,805)
      write(7,805)
    endif
  endif
803 format(1x,'We are communicating via MPI',/)
805 format(1x,'We are not communicating',/)

  call print_truncation_summary(pol%icutv, pol%truncval(1), iunit=6)
  call print_truncation_summary(pol%icutv, pol%truncval(1), iunit=7)

  if (peinf%inode==0 .and. pol%intraband_flag/=0) then
    if (pol%intraband_flag<0 .or. pol%intraband_flag>5) &
      call die('Invalid value for flag `intraband_flag`', only_root_writes=.true.)
    if (pol%intraband_overlap_min<0d0 .or. pol%intraband_overlap_min>1d0) &
      call die('Invalid value for flag `intraband_overlap_min`', only_root_writes=.true.)
    if (pol%intraband_flag>=1 .or. pol%intraband_flag<=2) then
      if (pol%intraband_flag==1) then
        str_intra = 'intraband'
      else
        str_intra = 'interband'
      endif
      write(6,'(/,1x,3a,f9.6,/)') 'Calculating only ', str_intra, &
        ' transitions, intraband_overlap_min = ', pol%intraband_overlap_min
      if (.not.pol%skip_epsilon) then
        write(0,'(/1x,3a)') 'WARNING: you are only calculating ', str_intra, ' transitions, but'
        write(0,'(1x,a)') 'you are not skipping epsilon. Remember that inverse dielectric matrices cannot'
        write(0,'(1x,a,/)') 'not be added together, only the polarizability!'
      endif
    else
      if (.not.pol%skip_epsilon) then
        write(0,'(/1x,3a)') 'WARNING: intraband_flag is not zero, but'
        write(0,'(1x,a)') 'you are not skipping epsilon. Remember that inverse dielectric matrices cannot'
        write(0,'(1x,a,/)') 'not be added together, only the polarizability!'
      endif
    endif
  endif
  if (pol%skip_epsilon) then
    pol%matrix_type = 2
  else
    pol%matrix_type = 0
  endif

!#BEGIN_INTERNAL_ONLY
  if (peinf%inode==0.and.pol%num_cond_bands_ignore>0) then
    write(6,'(/,1x,a,i0,/)') 'Number of conduction bands to ignore: ', pol%num_cond_bands_ignore
  endif
  if (peinf%inode==0 .and. abs(pol%lin_denominator)>TOL_ZERO) then
    write(0,'(/1x,a/)') 'WARNING: Linearized energy denominator is under development!'
    ! FHJ: LED only works with matrix comm.
    pol%gcomm = -1
    if (pol%freq_dep_method==1) then
      call die('Cannot use linearized energy denom. with the spectral method.', &
        only_root_writes=.true.)
    endif
  endif
!#END_INTERNAL_ONLY

  call verify_gpu_settings()

  call require_reference(REF_Deslippe2012)
  call require_reference(REF_Hybertsen1986)
  ! Truncation of the Coulomb potential: slab and write
  if (pol%icutv==TRUNC_SLAB .or. pol%icutv==TRUNC_WIRE) call require_reference(REF_IsmailBeigi2016)
  ! FF algorithm by Shishkin and Kresse, as implemented by F. Liu et al
  if (pol%freq_dep_method==1) call require_reference(REF_Liu2015)
  ! Non-uniform sampling schemes
  if (pol%subsample .or. pol%non_uniform) call require_reference(REF_Jornada2017)
  ! Subspace for G0W0
  if (pol%subspace) call require_reference(REF_DelBen2019Subspace)
  ! Decomposition of polarizability into interband, intraband, etc.
  if (pol%intraband_flag/=0) call require_reference(REF_Jornada2020)


!----------------------
! Free memory

  SAFE_DEALLOCATE(qpt_read)
  SAFE_DEALLOCATE(qw_read)
  SAFE_DEALLOCATE(qflags)
  SAFE_DEALLOCATE(band_occ)

#ifdef MPI
  ! root lets the others go after it is done reading (see beginning of function)
  if(peinf%inode == 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  POP_SUB(inread)
  
  return
  
105 write(errmsg,'(3a)') 'The end of the file was reached while reading elements of the ', &
      trim(blockword),' block.'
  call die(errmsg, only_root_writes = .true.)
  
110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', &
      trim(keyword), '. '
  call die(errmsg, only_root_writes = .true.)
  
end subroutine inread

end module inread_m
