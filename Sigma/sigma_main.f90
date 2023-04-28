!===============================================================================
!
! Routines:
!
! (1) sigma              Originally By MSH        Last Modified 10/5/2009 (gsm)
!
! This is the main routine for the Sigma code.  Please see the documentation
! in the README for more information on this code.
!
!===============================================================================

#include "f_defs.h"

program sigma

  use accel_fft_m, only : accel_mtxel_sig, &
                                  allocate_accel_mtxel_sig, &
                                  deallocate_accel_mtxel_sig
  use algos_sigma_m
  use ch_converge_m
  use check_screening_m
  use checkbz_m
  use checkgriduniformity_m
  use epscopy_m
  use epsread_hdf5_m
  use fftw_m
  use fixwings_m
  use fullbz_m
  use genwf_mpi_m
  use global_m
  use gmap_m
#ifdef HDF5
  use hdf5
#endif
  use input_m
  use input_outer_m
  use input_utils_m
  use io_utils_m
  use irrbz_m
  use misc_m
  use mtxel_cor_m
  use mtxel_m
  use mtxel_occ_m
  use mtxel_vxc_m
  use references_m
  use shiftenergy_dyn_m
  use shiftenergy_m
  use sort_m
  use subgrp_m
  use timing_m, only: common_timing, timing => sigma_timing
  use vcoul_generator_m
  use wfn_rho_vxc_io_m
  use write_program_header_m
  use write_result_dyn_hp_m
  use write_result_dyn_m
  use write_result_hp_m
  use write_result_m
  implicit none

!---------------------
! Derived Types

  type (crystal) :: crys
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (kpoints) :: kp
  type (siginfo) :: sig
  type (wpgen) :: wpg
  type (wfnkstates) :: wfnk,wfnkoff
  type (wfnkqstates) :: wfnkq
  type (epsmpiinfo) :: epsmpi
  type (wfnkqmpiinfo) :: wfnkqmpi
  type (wfnkmpiinfo) :: wfnkmpi
  type (twork_scell) :: work_scell
  type(progress_info) :: prog_info !< a user-friendly progress report
!------ Electron Phonon (EP) ----------
! ZL: allocate same variables to read wave functions 
! while without overwriting the useful information
! phonq/phonv: electron phonon q/v (v: nu, phonon mode band index)
! for perturbed wave functions: dqv wfn(k)
  type (crystal) :: ep_crys
  type (symmetry) :: ep_syms
  type (gspace) :: ep_gvec
  type (kpoints) :: ep_kp
  !  type (siginfo) :: ep_sig  ! sig will be re-used because somehow inread
                                 ! does not allow to be called twice
  type (wpgen) :: ep_wpg
  type (wfnkstates) :: ep_dq_wfnk
  type (wfnkqstates) :: ep_dq_wfnkq
  type (wfnkqmpiinfo) :: ep_dq_wfnkqmpi
  type (wfnkmpiinfo) :: ep_dq_wfnkmpi
  ! extra variables
  type (wfnkstates) :: wfnk_phonq, wfnk_phonq_off ! wfn(k+phonq), regular wfn at a different point
                                                  ! wfnk_phonq plyas the role of wfnk
                                                  ! wfnk_phonq_off, for dVXC mtxel
  type (wfnkmpiinfo) :: wfnk_phonq_mpi ! distributed
  type (wfnkstates) :: ep_dq_wfnk_tmp ! temporary variable
  type (wfnkmpiinfo) :: ep_dq_wfnkmpi_tmp ! temporary variable

  ! define variables associated with the second term in formalism
  type (wfnkqstates) :: wfnkq_phonq ! inner wfn(k-q+phonq)
  type (wfnkqstates) :: ep_dmq_wfnkq_phonq ! inner d_{-phonq,v} wfn(k-q+phonq), need time-reversal
                                       ! "mq" represents "minus q": "-q"

!---------------------
! k-points for the sum over BZ (the rq-points)

  integer :: nm, nrq, iout, iparallel, igp, igp_loc, igp_owner, nq0
  integer, allocatable :: neq(:),indrq(:),itnrq(:),kg0(:,:)
  real(DP), allocatable :: rq(:,:)
  type(grid) :: gr
  
!---------------------
! (k-q) kpoints involved in sigma summations
! ===> see data type cwfnkq and vwfnkq

  real(DP) :: rkq(3)
  real(DP) :: rkq_ep(3), rk_phonq(3), r_mkq_mphonq(3)  ! used for EP
  integer :: ib_trs_ep, is_trs_ep ! looper for bands and spins, in the time-reversal case

!---------------------
! Dielectric matrices
! ZL: for EP under the static screening approximation, use the same eps

  integer :: ngq, neps, nmtx, ncoul, ncoulch, ncoulb, ncouls, ngpown_q
  integer :: ngqt,nmtxt,nfreqgpp
  integer, allocatable :: isrtq(:),isrtqi(:)
  integer, allocatable :: isrtrq(:)
  integer, pointer :: isrtrqi(:)
  SCALAR, pointer :: eps(:,:)
  real(DP), allocatable :: vcoul(:), ekin(:)
  complex(DPC), pointer :: epsR(:,:,:),epsA(:,:,:)

!---------------------
! Matrix elements for sigma

  SCALAR, allocatable :: aqs(:,:), aqsaug(:,:,:,:), alda(:,:), alda2(:,:), &
    ax(:,:), asx(:,:,:), ach(:,:,:), asig(:,:), ach_n1(:,:,:), achcor_n1(:,:,:), akih(:,:)
  ! ZL: alda is for VXC, and akih(:,:) is for KIH: Kinetic+Ionic+Hartree
  SCALAR, pointer :: aqsch(:), aqsaugchd(:,:,:), aqsaugcho(:,:,:), acht_n1(:), achtcor_n1(:)
  real(DP), allocatable :: enew(:,:),efsto(:,:),zrenorm(:,:)

  complex(DPC), allocatable :: achcor(:,:)
  complex(DPC), allocatable :: asig_imag(:,:)
  complex(DPC), allocatable :: asxDyn(:,:,:), achDyn(:,:,:), achDyn_cor(:,:,:), &
    achDyn_corb(:,:,:), ach2Dyn(:,:,:), asigDyn(:,:), achD_n1(:,:,:)
  complex (DPC), pointer :: achtD_n1(:)
  complex(DPC), allocatable :: efstoDyn(:,:), enewDyn(:,:), enewDyn_nosr(:,:)
  integer, allocatable :: neqp1(:,:), neqp1_nosr(:,:)

!#BEGIN_INTERNAL_ONLY
!---------------------
! ZL: Matrix elements for electron phonon
  ! for first term: note that aqs = aqs_oner
  ! matrix elements M
  SCALAR, allocatable :: aqs_ep_dq_onel(:,:), aqsaug_ep_dq_onel(:,:,:,:)
  ! matrix elements for X, SX, CH, Sigma
  SCALAR, allocatable :: ax_ep_one(:,:), asx_ep_one(:,:,:), ach_ep_one(:,:,:), asig_ep_one(:,:)
  SCALAR, allocatable :: ach_n1_ep_one(:,:,:), achcor_n1_ep_one(:,:,:)
  complex(DPC), allocatable :: asig_imag_ep_one(:,:)
  ! aqsch is used only for exact_ch. Place holder for EP
  SCALAR, pointer :: aqsch_ep_one(:)
  complex(DPC), allocatable :: achcor_ep_one(:,:)
  ! temporary variables
  SCALAR, pointer :: asxt_ep_one(:), acht_ep_one(:)
  SCALAR :: axt_ep_one, achtcor_ep_one, asigt_imag_ep_one
  SCALAR, pointer :: acht_n1_ep_one(:), achtcor_n1_ep_one(:)

  ! for second term
  ! matrix elements M
  SCALAR, allocatable :: aqs_phonq_twol(:,:), aqs_ep_dmq_twor(:,:)
  SCALAR, allocatable :: aqsaug_phonq_twol(:,:,:,:), aqsaug_ep_dmq_twor(:,:,:,:)
  ! matrix elements for X, SX, CH, Sigma
  SCALAR, allocatable :: ax_ep_two(:,:), asx_ep_two(:,:,:), ach_ep_two(:,:,:), asig_ep_two(:,:)
  SCALAR, allocatable :: ach_n1_ep_two(:,:,:), achcor_n1_ep_two(:,:,:)
  complex(DPC), allocatable :: asig_imag_ep_two(:,:)
  SCALAR, pointer :: aqsch_ep_two(:)
  complex(DPC), allocatable :: achcor_ep_two(:,:)
  ! temporary variables
  SCALAR, pointer :: asxt_ep_two(:), acht_ep_two(:)
  SCALAR :: axt_ep_two, achtcor_ep_two, asigt_imag_ep_two
  SCALAR, pointer :: acht_n1_ep_two(:), achtcor_n1_ep_two(:)

  ! for sum of two terms
  SCALAR, allocatable :: ax_ep(:,:), asx_ep(:,:,:), ach_ep(:,:,:), asig_ep(:,:), &
                         ach_n1_ep(:,:,:), achcor_n1_ep(:,:,:), alda_ep(:,:)
  complex(DPC), allocatable :: achcor_ep(:,:), asig_imag_ep(:,:)
!#END_INTERNAL_ONLY

!----------------------
! eps distrib variables
  SCALAR, allocatable :: epstemp(:)
  
  character :: tmpstr*120
  character :: tmpfn*16
  character*20 :: fnc,fnk,fne
  character*16 :: routnam(100)
  integer :: routsrt(59),nullvec(3)
  logical :: xflag,imagvxcflag,imagxflag,found,q0flag,bExactlyZero, imagkihflag ! ZL: last one for KIH
  logical :: eqp1_warns(4) ! (GPP extrap, FF extrap, FF multiple solns, FF no soln)
  integer :: ig,i,j,k,itran,ikn,ika,ioff,error
  integer :: in,im,iw,ib,jb,idum,kg(3),jj,ii,ispin,jsp,g1,g2
  integer :: ncount,ndum,nbandi,tag,dest,source,nfold,ifold
  integer :: iwlda,irq,irq_,irq_min,n1,ierr
  integer :: s2,iunit_c,iunit_k,iunit_eps,ndv_ikn,iunit
  integer, allocatable :: ind(:), indinv(:)
  real(DP) :: fact,coulfact,weight,tempval,occ
  real(DP) :: qshift(3),oneoverq,qlen,q0len,vq(3),qk(3)
  real(DP) :: tsec(2),diffmin,diff,e_lk,avgcut,subcut,freq0
  complex(DPC), pointer :: achtDyn(:),achtDyn_cor(:),asxtDyn(:),ach2tDyn(:),achtDyn_corb(:)
  SCALAR :: achtcor,axt,epshead,asigt_imag
  SCALAR, pointer :: asxt(:), acht(:)
  SCALAR, allocatable :: ph(:)

  logical :: skip_checkbz, is_subq

!------------- ZL: Electron phonon variables --------------------------------
  logical :: ep_read ! controls reading dWFN of EP
  logical :: ep_debug ! controls debug of EP
  integer :: ik_phonq_idx ! similar role of ikn
  real(DP) :: k_phonq_coord(3) ! similar role of qk(3)
  logical :: check_norms_save
  integer :: ioff_check
  real(DP) :: Eo_save

! sigma-subspace variables --------------------------------------------------
  integer :: ipe_wing, my_pos_loc_wing, iproc_dum 
  integer :: nfreq_fixwings, ifreq
  complex(DPC), pointer :: epsR_corrections(:,:,:)

  integer, dimension(3) :: Nfft
  real(DP) :: dummy_scale

  ep_read = .false.
  ep_debug = .false.

!--------------- Begin Program -------------------------------------------------

  call peinfo_init()

#ifdef PGI
  call setvbuf3f(6,2,0)
#endif

!----------------------
! Initialize random numbers

  peinf%jobtypeeval = 1

!------------------------
! Initialize timer
  call timing%init()
  call common_timing%init()
  call timing%start(timing%total)

!------------------------
! Initialize files

  call open_file(55,file='sigma.inp',form='formatted',status='old')
  if(peinf%inode == 0) then
    call open_file(8,file='sigma_hp.log',form='formatted',status='replace')
    call open_file(30,file='eqp0.dat',form='formatted',status='replace')
    call open_file(31,file='eqp1.dat',form='formatted',status='replace')
  endif

  call write_program_header('Sigma', .false.)

!------- Read crys data and wavefunctions from WFN_inner ----------------------------

! JRD: Included in input is the inread routine which reads the
! job parameters from sigma.inp and initializes the XC potential

  call timing%start(timing%input)

!#BEGIN_INTERNAL_ONLY
  ! ZL: EP, add to permanent arguments: wfnk_phonq, wfnk_phonq_mpi
  !     They CANNOT be made optional intrisincally!
  !     Because if allocate and read them is determined inside input(), after inread()
  !     Note: the two variables will not be allocated and evaluated if not doing EP
!#END_INTERNAL_ONLY
  call input(crys,gvec,syms,kp,wpg,sig,wfnk,iunit_c,iunit_k,fnc,fnk,wfnkqmpi,wfnkmpi,wfnk_phonq,wfnk_phonq_mpi,ep_read_in=.false.)

  call output_algos()

!------------------------
! Initialize GPU
  if (sigma_gpp_algo == OPENACC_ALGO .or. mtxel_algo == OPENACC_ALGO) then
    call initialize_gpu(OPENACC_ALGO)
  end if
  if (sigma_gpp_algo == OMP_TARGET_ALGO .or. mtxel_algo == OMP_TARGET_ALGO) then
    call initialize_gpu(OMP_TARGET_ALGO)
  end if

!#BEGIN_INTERNAL_ONLY
  ! ZL: check if EP is in complex version
#ifndef CPLX
  if (sig%elph) then
    call die("Electron-phonon code MUST BE complied with COMPLEX version!")
  endif
#endif
  ! ZL: EP output information
  !     for electron-phonon calculations
  if(sig%elph .and. peinf%inode.eq.0) then
    write(6,*)
    write(6,*) "Symmetry information"
    write(6,*) "number of symmetries", syms%ntran
    write(6,*) "symmetry operations", syms%mtrx(:,:,1:syms%ntran)
    write(6,*)
    if (syms%ntran .gt. 1) then
      call die('For elph, TURN OFF ALL SYMMETRIES in ALL DFT calculations.')
    endif
    write(6,*) "number of kpoints in outer, read from sigma.inp =", sig%nkn
    write(6,*) "kpoints in outer ="
    do k=1, sig%nkn
      write(6,*) sig%kpt(1:3, k)
    enddo
    write(6,*) "number of phonon_q points =", sig%nphonq
    write(6,*) "phonon_q (supposed to be only one for now)  =", sig%phonq
    write(6,*) "k+phonon_q points mapped from sigma.inp ="
    do k=1, sig%nkn*sig%nphonq
      write(6,*) sig%k_phonq(1:3, k)
    enddo
    write(6,*)
    write(6,*) "sig%ep_bracket =", sig%ep_bracket
    write(6,*)
    write(6,*) "Fermi level in sig", sig%efermi
    write(6,*)
  endif

  if (sig%elph .and. sig%nphonq == 0) then 
    if (peinf%inode .eq. 0) write (6,*) "Linear respoinse turned on, but NO found phonq! Exiting..."
    call die('Linear respoinse turned on, but NO found phonq!', only_root_writes=.true.)
  endif

  ! ZL: for EP, set up offdiags for dSig(E)
  if(sig%elph) then
    if(sig%ep_bracket.eq.0) then
      call die('You must specify at which state, bra(-1) or ket(1), dSig(E) is evaluated!')
    endif

    if(sig%ep_bracket.eq.1) then
      if(peinf%inode.eq.0) write(6,*) 'Electron-phonon dSig(E) is evaluated at |ket> = |k> state'
    elseif (sig%ep_bracket.eq.-1) then
      if(peinf%inode.eq.0) write(6,*) 'Electron-phonon dSig(E) is evaluated at <bra| = <k+phonq| state'
    else
      call die('ep_bracket must be specified as 1 (for ket), or -1 (for bra).')
    endif
    if(peinf%inode.eq.0) write(6,*) 

    if(sig%noffdiag .gt. 0) then
      !  ZL: check if the three indices are properly set up
      !      users should know what they are doing here, but we check anyways
      if(sig%ep_bracket.eq.1) then
        do ioff_check = 1, sig%noffdiag
          if(sig%off2(ioff_check) .ne. sig%off3(ioff_check)) then
            call die('Offdiag is inconsistent with ep_bracket!')
          endif
        enddo
      elseif(sig%ep_bracket.eq.-1) then
        do ioff_check = 1, sig%noffdiag
          if(sig%off1(ioff_check) .ne. sig%off3(ioff_check)) then
            call die('Offdiag is inconsistent with ep_bracket!')
          endif
        enddo
      endif
    endif
  else
    if(sig%ep_bracket.ne.0) then
      call die('You should not set ep_bracket if you are not doing a electron-phonon calculation.')
    endif
  endif ! set up offdiags for dSig(E)

  if (sig%elph) then ! read dWFN, read first, to protect use info (if saving memory)
      ep_read = .true.
      ! ZL: ep_dq_wfnk_tmp, ep_dq_wfnkmpi_tmp will not be even allocated, they are just here to satisfy 
      !     subroutine input() argument rule - place holder
      call input(ep_crys,ep_gvec,ep_syms,ep_kp,ep_wpg,   sig,   ep_dq_wfnk,&
             iunit_c,iunit_k,fnc,fnk,ep_dq_wfnkqmpi,ep_dq_wfnkmpi, ep_dq_wfnk_tmp, ep_dq_wfnkmpi_tmp, ep_read_in=ep_read)
           ! ZL: after read in the ep_ variables, we didn`t update sigma, and 
           !     ep_kp needs some updates (although they would
           !     not be used in the following, but just in case. Note the most
           !     imporant info is ngk, which has been correctly assigned.
           !     Actually, band energies in ep_dq_wfnkmpi are not meaningful,
           !     and indeed hard to define. We are not using them anyways
           !     because we always use energies in outer wfn: wfnk and wfnk_phonq
      ep_kp%ifmin(:,:) = kp%ifmin(:,:)
      ep_kp%ifmax(:,:) = kp%ifmax(:,:)
      ep_kp%occ(:,:,:) = kp%occ(:,:,:)
      ep_kp%degeneracy(:,:,:) = kp%degeneracy(:,:,:)
  endif ! for sig%elph
  ! ZL: after reading for sig%elph, ep_gvec, ep_kp, ep_dq_wfnk, ep_dq_wfnkqmpi,
  ! ep_dq_wfnkmpi have been initialized

  ! ZL: output and check
  !     ep_gvec and gvec should be EXACTLY the same gspace
  if (sig%elph .and. peinf%inode==0) then
    write(6,*)
    write(6,*) 'gvec info', gvec%ng
    write(6,*) 'ep_gvec info', ep_gvec%ng
    write(6,*)
    if(sig%nkn.gt.1) then
      write(6,*) "kpoints in wfnkmpi"
      do k=1, sig%nkn
        write(6,*) wfnkmpi%qk(1:3, k)
      enddo
      write(6,*) "ngk", wfnkmpi%nkptotal(:)
      write(6,*)
      write(6,*) "k+phonon_q points in wfnk_phonq_mpi"
      do k=1, sig%nkn
        write(6,*) wfnk_phonq_mpi%qk(1:3, k)
      enddo
      write(6,*) "ngk", wfnk_phonq_mpi%nkptotal(:)
      write(6,*)
      write(6,*) "Printing dWFN reading information"
      write(6,*) "ep_dq_wfnkmpi will not be used in calculations, but contains info"
      write(6,*) "kpoints in ep_dq_wfnkmpi and energy (useless)"
      do k=1, sig%nkn
        write(6,*) "rk = ", ep_dq_wfnkmpi%qk(1:3, k)
  !      write(6,*) "Band energy"
  !      write(6,*) ep_dq_wfnkmpi%el(:,:,k)
      enddo
      write(6,*) "ep_dq_wfnkmpi ngk"
      write(6,*) ep_dq_wfnkmpi%nkptotal
  !    write(6,*) "Bonus: wfnk_phonq_mpi band energy"
  !    do k=1, sig%nkn
  !      write(6,*) wfnk_phonq_mpi%el(:,:,k)
  !    enddo
    endif ! sig%nkn.gt.0
    write(6,*) "Fermi level in sig", sig%efermi
    write(6,*)
  endif
!#END_INTERNAL_ONLY

  call timing%stop(timing%input)
  call timing%start(timing%input_outer)

  ! ZL: for EP, we do not use WFN_outer. Everything should be in WFN_inner and dWFN
  if (sig%elph) then
!#BEGIN_INTERNAL_ONLY
    sig%wfn_outer_present = .false.
    if(peinf%inode.eq.0) write(6,*) "Electron phonon calculation does not read WFN_outer."
!#END_INTERNAL_ONLY
  else
    call input_outer(crys,gvec,syms,kp,sig,wfnk,iunit_k,fnk,wfnkmpi)
  endif

  call timing%stop(timing%input_outer)
  SAFE_DEALLOCATE_P(sig%kpt)
  SAFE_DEALLOCATE(kp%ifmin)
  SAFE_DEALLOCATE(kp%ifmax)

!-------------------
! Initialize Various Parameters from inread

  eqp1_warns(:)=.false.

! ZL: KIH not supported with Real flavor
#ifndef CPLX
  if (sig%use_kihdat) then
    call die("KIH usage for arbitrary functional starting point MUST BE complied with COMPLEX version!")
  endif
#endif

! imaginary parts of diagonal vxc or exchange matrix elements
  imagvxcflag = .false.
  imagxflag = .false.
  imagkihflag = .true. ! ZL: always true for kih
!#BEGIN_INTERNAL_ONLY
  if(sig%elph) then
! ZL: for quasiparticle energies, imaginary parts are not necessary
!     now I turn them on because everything is complex in EP (electron phonon)
    imagvxcflag = .true.
    imagxflag = .true.
  endif
!#END_INTERNAL_ONLY

! fraction of bare exchange

  if (abs(sig%xfrac).GT.TOL_Small) then
    xflag=.true.
  else
    xflag=.false.
  endif

!---------------------
! Open Various Files 

  if (peinf%inode .eq. 0) then
    if (.not.(sig%freq_dep .eq. 0 .and. sig%exact_ch .eq. 1) .and. .not. (sig%freq_dep == -1)) then
      call open_file(127,file='ch_converge.dat',form='formatted',status='replace')
    endif

    if (sig%iwritecoul .eq. 1) then
      call open_file(19,file='vcoul',form='formatted',status='replace')
    endif

   ! This if for the hybrid functional calculations (one shot) otherwise just open x.dat
    if (sig%coul_mod_flag .and. (.not. sig%use_vxc2dat)) then
      call open_file(121,file='vxc2.dat',form='formatted',status='replace')
    else if ((.not. sig%use_xdat) .and. xflag .and. (.not. sig%coul_mod_flag)) then
      call open_file(119,file='x.dat',form='formatted',status='replace')
    endif

    if (.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX) then
      if(.not.sig%use_kihdat) then ! ZL: only meaningful to generate vxc.dat if not running with KIH
        call open_file(120,file='vxc.dat',form='formatted',status='replace')
      endif
    endif
  endif

!---------------------
! Write header of sigma_hp.log file

  if (peinf%inode.eq.0) then
    write(8,601) sig%freq_dep
    write(8,602) sig%bmin,sig%bmax
    write(8,603) sig%loff,sig%toff
    write(8,604) sig%fdf
    if(sig%fdf /= -2) write(8,605) sig%dw
    write(8,606) syms%ntran
    do itran=1,syms%ntran
      write(8,607) itran, ((syms%mtrx(i,j,itran),i=1,3),j=1,3)
    enddo
    write(8,*)
  endif
601 format(/,1x,"frequency_dependence",i4)
602 format(/,1x,"band_index",2i6)
603 format(1x,"sigma_matrix",i6,i4)
604 format(/,1x,"finite_difference_form",i4)
605 format(1x,"finite_difference_spacing",f10.6)
606 format(/,1x,"symmetries",/,1x,"ntran  =",i3)
607 format(1x,"mtrx",i2.2,1x,"=",9i3)

!---------------------
! JRD: Initialize the Full Frequency output files

  if (peinf%inode.eq.0 .and. (sig%freq_dep.eq.2 .or. (sig%fdf.eq.-3 .and. sig%freq_dep.eq.1))) then
    call open_file(8000,file='spectrum.dat',form='formatted',status='replace')
  endif

!---------------------
! Determine nq and neps
  
! JRD: This performs significantly better with hdf5

  call timing%start(timing%read_neps)

  call epscopy_init(sig, neps)
  if (sig%freq_dep/=-1) then
    ! FHJ: sig%nq and sig%qpt already defined if this is a HF calculation
    sig%nq = sig%nq0 + sig%nq1 
    SAFE_ALLOCATE(sig%qpt, (3,sig%nq))
  endif
  if (sig%nq0==0) call die('There is no q->0 point in your calculation!', only_root_writes=.true.)

  epsmpi%nb = 1
  epsmpi%ngpown = NUMROC(neps, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
  epsmpi%ngpown_max = NUMROC(neps, epsmpi%nb, 0, 0, peinf%npes_pool)

  ! define distribution for subspace matrices (MDB)
  if(sig%do_sigma_subspace) then
    ! this has to 1 (no change)
    sig%epssub%nb_sub = 1
    !
    sig%epssub%Nbas_own_max =  NUMROC(sig%neig_sub_max, sig%epssub%nb_sub, 0, 0, peinf%npes_pool)
    sig%epssub%Nbas_own = MIN(sig%epssub%Nbas_own_max*peinf%pool_rank + 1, sig%neig_sub_max)
    sig%epssub%Nbas_own = MIN(sig%epssub%Nbas_own_max*(peinf%pool_rank+1), sig%neig_sub_max) - sig%epssub%Nbas_own + 1
    ! note the block size here is the same as the epsmpi (not a good idea, this will be removed in the future)
    sig%epssub%ngpown_sub_max = NUMROC(neps, epsmpi%nb, 0, 0, peinf%npes_pool)
    sig%epssub%ngpown_sub = MIN(sig%epssub%ngpown_sub_max*peinf%pool_rank + 1, neps)
    sig%epssub%ngpown_sub = MIN(sig%epssub%ngpown_sub_max*(peinf%pool_rank+1), neps) - sig%epssub%ngpown_sub + 1
    sig%epssub%neps = neps
    !
    !XXX write(*,*) peinf%inode, peinf%pool_rank, sig%epssub%Nbas_own, sig%epssub%Nbas_own_max,&
    !XXX            sig%epssub%ngpown_sub, sig%epssub%ngpown_sub_max
  end if

!----------------------------
! Allocate arrays
! ZL: allocate wfnkq contents

! wfn(k-q) where q is the internal regular q (or p in the equation)
  SAFE_ALLOCATE(wfnkq%isrtkq, (gvec%ng))
  SAFE_ALLOCATE(wfnkq%ekq, (sig%ntband,kp%nspin))
  SAFE_ALLOCATE(alda, (sig%ndiag+sig%noffdiag,sig%nspin))
  ! ZL: add for array KIH
  SAFE_ALLOCATE(akih, (sig%ndiag+sig%noffdiag,sig%nspin))
!#BEGIN_INTERNAL_ONLY
  ! ZL: EP counterpart of wfnkq: ep_dq_wfnkq, d_phonq_v WFN(k-q) where q is p in formalism
  !   ep_gvec should be EXACTLY the same as gvec in generating wfns
  if(sig%elph) then
    ! for first term
    SAFE_ALLOCATE(ep_dq_wfnkq%isrtkq, (ep_gvec%ng))
    SAFE_ALLOCATE(ep_dq_wfnkq%ekq, (sig%ntband, kp%nspin))
    ! for second term
    SAFE_ALLOCATE(wfnkq_phonq%isrtkq, (ep_gvec%ng))
    SAFE_ALLOCATE(wfnkq_phonq%ekq, (sig%ntband, kp%nspin))
    SAFE_ALLOCATE(ep_dmq_wfnkq_phonq%isrtkq, (ep_gvec%ng))
    SAFE_ALLOCATE(ep_dmq_wfnkq_phonq%ekq, (sig%ntband, kp%nspin))
    ! for dVXC
    SAFE_ALLOCATE(alda_ep, (sig%ndiag+sig%noffdiag,sig%nspin))
  endif ! sig%elph
!#END_INTERNAL_ONLY
  
!#BEGIN_INTERNAL_ONLY
  if (sig%coul_mod_flag) then
    SAFE_ALLOCATE(alda2, (sig%ndiag+sig%noffdiag,sig%nspin))
  endif
!#END_INTERNAL_ONLY
! ZL: 'a' means ARRAY, for bare exchange
  SAFE_ALLOCATE(ax, (sig%ndiag+sig%noffdiag,sig%nspin))

!#BEGIN_INTERNAL_ONLY
! EP: perturbed bare exchange
  if (sig%elph) then
    SAFE_ALLOCATE(ax_ep_one, (sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(ax_ep_two, (sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(ax_ep, (sig%ndiag+sig%noffdiag,sig%nspin))
  endif ! sig%elph
!#END_INTERNAL_ONLY

! achcor for static remainder
  SAFE_ALLOCATE(achcor, (sig%ndiag+sig%noffdiag,sig%nspin))
  SAFE_ALLOCATE(asig_imag, (sig%ndiag+sig%noffdiag,sig%nspin))
  SAFE_ALLOCATE(achcor_n1, (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))
!#BEGIN_INTERNAL_ONLY
  if (sig%elph) then
    SAFE_ALLOCATE(achcor_ep_one, (sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(asig_imag_ep_one, (sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(achcor_n1_ep_one, (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))

    SAFE_ALLOCATE(achcor_ep_two, (sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(asig_imag_ep_two, (sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(achcor_n1_ep_two, (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))

    SAFE_ALLOCATE(achcor_ep, (sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(asig_imag_ep, (sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(achcor_n1_ep, (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))
  endif
!#END_INTERNAL_ONLY

  if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
      ! static remainder for 1 band, allocated for all bands
    SAFE_ALLOCATE(achtcor_n1, (sig%ntband))
!#BEGIN_INTERNAL_ONLY
    if (sig%elph) then
      SAFE_ALLOCATE(achtcor_n1_ep_one, (sig%ntband))
      SAFE_ALLOCATE(achtcor_n1_ep_two, (sig%ntband))
      call die('Electron phonon not set up for exact_ch yet!')
    endif
!#END_INTERNAL_ONLY
  endif
  if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
    if (sig%fdf .eq. -3) then ! fdf: finite difference form, way to eval quasi energy, -3 many freqs
      ! ZL: we do not consider unusual cases for EP now
      SAFE_ALLOCATE(asx, (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
      SAFE_ALLOCATE(ach, (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
      if (sig%elph) call die('Electron phonon not set up for fdf = -3 yet')
    else
      ! at most use 3 freqs
      SAFE_ALLOCATE(asx, (3,sig%ndiag+sig%noffdiag,sig%nspin))
      SAFE_ALLOCATE(ach, (3,sig%ndiag+sig%noffdiag,sig%nspin))
!#BEGIN_INTERNAL_ONLY
      if (sig%elph) then
        SAFE_ALLOCATE(asx_ep_one, (3,sig%ndiag+sig%noffdiag,sig%nspin))
        SAFE_ALLOCATE(ach_ep_one, (3,sig%ndiag+sig%noffdiag,sig%nspin))
        SAFE_ALLOCATE(asx_ep_two, (3,sig%ndiag+sig%noffdiag,sig%nspin))
        SAFE_ALLOCATE(ach_ep_two, (3,sig%ndiag+sig%noffdiag,sig%nspin))
        SAFE_ALLOCATE(asx_ep, (3,sig%ndiag+sig%noffdiag,sig%nspin))
        SAFE_ALLOCATE(ach_ep, (3,sig%ndiag+sig%noffdiag,sig%nspin))
      endif
!#END_INTERNAL_ONLY
    endif
    ! total self energy for final results (maybe)
    SAFE_ALLOCATE(asig, (sig%ndiag+sig%noffdiag,sig%nspin))
!#BEGIN_INTERNAL_ONLY
    ! EP
    if(sig%elph) then
      SAFE_ALLOCATE(asig_ep_one, (sig%ndiag+sig%noffdiag,sig%nspin))
      SAFE_ALLOCATE(asig_ep_two, (sig%ndiag+sig%noffdiag,sig%nspin))
      SAFE_ALLOCATE(asig_ep, (sig%ndiag+sig%noffdiag,sig%nspin))
    endif
!#END_INTERNAL_ONLY
    ! t means temperory
    SAFE_ALLOCATE(acht_n1, (sig%ntband))
    ! COH resolved for all bands, convergence usage
    SAFE_ALLOCATE(ach_n1, (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))
!#BEGIN_INTERNAL_ONLY
    if (sig%elph) then
      SAFE_ALLOCATE(acht_n1_ep_one, (sig%ntband))
      SAFE_ALLOCATE(ach_n1_ep_one, (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))
      SAFE_ALLOCATE(acht_n1_ep_two, (sig%ntband))
      SAFE_ALLOCATE(ach_n1_ep_two, (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))
      SAFE_ALLOCATE(ach_n1_ep, (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))
    endif
!#END_INTERNAL_ONLY
    ! enew eqp1
    SAFE_ALLOCATE(enew, (sig%ndiag,sig%nspin))
    ! efsto eqp0
    SAFE_ALLOCATE(efsto, (sig%ndiag,sig%nspin))
    ! Z renorm factor
    SAFE_ALLOCATE(zrenorm, (sig%ndiag,sig%nspin))
    if (sig%fdf.eq.-3) then
      nfreqgpp=sig%nfreqeval
      SAFE_ALLOCATE(asxt, (sig%nfreqeval))
      SAFE_ALLOCATE(acht, (sig%nfreqeval))
      ! ZL: if this is the case for EP, it should have happened already
      if(sig%elph) call die('Electron phonon not set up for fdf = -3 yet')
    else
      nfreqgpp=3
      SAFE_ALLOCATE(asxt, (3))
      SAFE_ALLOCATE(acht, (3))
!#BEGIN_INTERNAL_ONLY
      if (sig%elph) then
        SAFE_ALLOCATE(asxt_ep_one, (3))
        SAFE_ALLOCATE(acht_ep_one, (3))
        SAFE_ALLOCATE(asxt_ep_two, (3))
        SAFE_ALLOCATE(acht_ep_two, (3))
      endif
!#END_INTERNAL_ONLY
    endif
  endif
  ! full frequency calculation, labeled by Dyn
  if (sig%freq_dep.eq.2) then
    SAFE_ALLOCATE(asxDyn, (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(achDyn, (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(achDyn_cor, (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(achDyn_corb, (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(ach2Dyn, (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(asigDyn, (sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(achtD_n1, (sig%ntband))
    SAFE_ALLOCATE(achD_n1, (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))
    SAFE_ALLOCATE(efstoDyn, (sig%ndiag,sig%nspin))
    SAFE_ALLOCATE(enewDyn, (sig%ndiag,sig%nspin))
    SAFE_ALLOCATE(enewDyn_nosr, (sig%ndiag,sig%nspin))
    SAFE_ALLOCATE(neqp1, (sig%ndiag,sig%nspin))
    SAFE_ALLOCATE(neqp1_nosr, (sig%ndiag,sig%nspin))
    SAFE_ALLOCATE(asxtDyn, (sig%nfreqeval))
    SAFE_ALLOCATE(achtDyn, (sig%nfreqeval))
    SAFE_ALLOCATE(achtDyn_cor, (sig%nfreqeval))
    SAFE_ALLOCATE(achtDyn_corb, (sig%nfreqeval))
    SAFE_ALLOCATE(ach2tDyn, (sig%nfreqeval))
    if(sig%elph) call die('Electron phonon not set up for full frequency calculation')
  endif
  ! ZL: isrtrq is of the length of all gvec pool
  SAFE_ALLOCATE(isrtrq, (gvec%ng))
  ! GPP is freq_dep = 1, the only case we consider for EP for now
  if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then
    SAFE_ALLOCATE(isrtrqi, (gvec%ng))
    if (sig%elph) call die('Electron phonon not set up for sig%freq_dep.eq.0.or.sig%exact_ch.eq.1')
  endif
  ! ZL: ekin is of the length gvec%ng total G
  SAFE_ALLOCATE(ekin, (gvec%ng))
  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
    ! ZL: this includes GPP
    ! isrtq is of the length of all gvec pool
    SAFE_ALLOCATE(isrtq, (gvec%ng))
    SAFE_ALLOCATE(isrtqi, (gvec%ng))
    SAFE_ALLOCATE(ind, (gvec%ng))
    SAFE_ALLOCATE(indinv, (gvec%ng))
    SAFE_ALLOCATE(ph, (gvec%ng))
    if (peinf%inode .eq. 0) then
      SAFE_ALLOCATE(epsmpi%isrtq, (gvec%ng,sig%nq))
      SAFE_ALLOCATE(epsmpi%isrtqi, (gvec%ng,sig%nq))
    endif
    SAFE_ALLOCATE(epsmpi%qk, (3,sig%nq))
    SAFE_ALLOCATE(epsmpi%nmtx, (sig%nq))
    SAFE_ALLOCATE(epsmpi%inv_igp_index, (epsmpi%ngpown_max))
  endif

  call timing%stop(timing%read_neps)
 
!----------------------------
! Read eps^-1 from eps0mat/epsmat
!
! JRD: The matrices are read in from eps0mat/epsmat files and writen
! to temporary INT_EPS files on unit iunit_eps. The q->0 matrix is not
! symmetrized. The wavevector q is also read from subroutine epscopy.

  call timing%start(timing%epscopy)
  epshead = ZERO
  if (sig%freq_dep/=-1) then
    call logit('Calling epscopy')
    ! FHJ: distribute columns of epsinv with standard ScaLAPACK 1d block-column
    ! layout. For convenience, the ScaLAPACK auxiliary functions are defined
    ! in BerkeleyGW even if you compile the code without ScaLAPACK.
    ! TODO: get rid of array inv_igp_index, it can be computed on the fly.
    epsmpi%inv_igp_index = 0
    do igp_loc = 1, epsmpi%ngpown
      igp = INDXL2G(igp_loc, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
      epsmpi%inv_igp_index(igp_loc) = igp
    enddo
    call epscopy(crys,gvec,sig,neps,epsmpi,epshead,iunit_eps,fne)
  endif

  call timing%stop(timing%epscopy)

  if (sig%subsample) then
    SAFE_ALLOCATE(sig%subweights, (sig%nq0))
    call open_file(666, file='subweights.dat', form='formatted', status='old')
    read(666,*) nq0
    if (nq0/=sig%nq0) then
      call die('Inconsistency between nq0 in subweights.dat and eps0mat.h5 files.')
    endif
    do irq=1,sig%nq0
      read(666,*) sig%subweights(irq)
    enddo
    call close_file(666)
    if (peinf%inode==0) write(6,'(/,1x,a,f0.6,/)') &
      'Sum of subweights before renormalization: ', sum(sig%subweights(:))
    sig%subweights(:) = sig%subweights(:) / sum(sig%subweights(:))
    ! FHJ: Choose the smallest cutoff that doesn`t include a lattice vector
    ! in a periodic direction. FIXME: there might be a combination, such as
    ! |Gx+Gy|^2, that gives a smaller cutoff than |Gx|^2 or |Gy|^2 alone!
    if (sig%icutv==TRUNC_WIRE) then
      write (0,'(/a/)') 'WANING: subsampling not tested for quasi-1D geometries'
      subcut = crys%bdot(3,3)
    elseif (sig%icutv==TRUNC_SLAB) then
      subcut = INF
      do ii=1,2
        subcut = min(subcut,crys%bdot(ii,ii)*(1d0-TOL_SMALL))
      enddo
    else
      write (0,'(/a/)') 'WANING: subsampling not tested for the specified truncation geometry'
      subcut = INF
      do ii=1,3
        if (sig%qgrid(ii)>1) subcut = min(subcut,crys%bdot(ii,ii)*(1d0-TOL_SMALL))
      enddo
    endif
    if (peinf%inode==0) write(6,'(1x,a,f0.9)') &
      'Cutoff for subsampled q-points: ', subcut
  endif

!----------------------------
! Generate full Brillouin zone from irreducible wedge q -> gr%f

  call timing%start(timing%fullbz)
  ! gr%nr: number in reduced zone
  ! gr%nf: number in full zone
  ! gr%r: points in reduced zone
  ! gr%f: points in full zone
  gr%nr = sig%nq
  SAFE_ALLOCATE(gr%r, (3, sig%nq))
  ! reduced q points 
  gr%r(1:3,1:sig%nq) = sig%qpt(1:3,1:sig%nq)
  ! wigner_seitz = .false. for Epsilon, Sigma, use usual "box" BZ
  !              = .true.  for BSE, use wigner seitz cell
  call fullbz(crys,syms,gr,syms%ntran,skip_checkbz,wigner_seitz=.false.,paranoid=.true.,nfix=sig%nq0)
  qshift(:)=0.0d0
  if (sig%freq_dep.eq.-1) then
    ! for Hartree-Fock, there is no epsmat/eps0mat file
    tmpfn="sigma.inp"
  else
    if (sig%igamma.ne.0) then
      tmpfn='eps0mat'
    else
      tmpfn='epsmat'
    endif
  endif
  if (.not. skip_checkbz) then
    !FHJ: TODO: ignore change checkbz to support nq0
    !call checkbz(gr%nf,gr%f,sig%qgrid,qshift,crys%bdot,tmpfn,'q',.false.,sig%freplacebz,sig%fwritebz,nfix=sig%nq0)
    call checkbz(gr%nf,gr%f,sig%qgrid,qshift,crys%bdot,tmpfn,'q',.false.,sig%freplacebz,sig%fwritebz)
  endif
  call timing%stop(timing%fullbz)
  if (peinf%inode==0) then
    write(6,'(1x,a,i0)') 'Number of k-points in WFN_inner: ', kp%nrk
    write(6,'(1x,a,i0)') 'Number of k-points in the full BZ of WFN_inner: ', gr%nf
!#BEGIN_INTERNAL_ONLY
    ! ZL add
    if (sig%elph) write(6,*) "number of kpoints in reduced, full grids =", gr%nr, gr%nf
!#END_INTERNAL_ONLY
  endif

! ZL: above are all preparing procedures for sigma and EP perturbed self-energy calculations
  

!-------- Start computation of sigma operator ----------------------------------
! ZL: here is general start
  
  fact = 1D0/(dble(gr%nf-sig%nq0+1)*crys%celvol)
  coulfact = 8D0*PI_D/(dble(gr%nf-sig%nq0+1)*crys%celvol)

!----------------------------
! Initialize distribution of epsilon

  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
    SAFE_ALLOCATE(eps, (neps,epsmpi%ngpown))
    SAFE_ALLOCATE(epstemp, (neps))
  endif
  if (sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
    if(sig%do_sigma_subspace) then
      !MDB only static epsinv will be calculated in the original PW basis
      nfreq_fixwings = 1
      SAFE_ALLOCATE(epsR, (neps,epsmpi%ngpown,1))
    else
      nfreq_fixwings = sig%nFreq
      SAFE_ALLOCATE(epsR, (neps,epsmpi%ngpown,sig%nFreq))
      if (sig%need_advanced) then
        SAFE_ALLOCATE(epsA, (neps,epsmpi%ngpown,sig%nFreq))
      endif
    end if
  endif

  if (kp%nspinor == 2) then
    call require_reference(REF_Barker2018)
  endif
!---------- Check grid for uniformity

! icutv: truncation type
  if(peinf%inode == 0) call checkgriduniformity(sig%qgrid, crys, sig%icutv)

!-------- Loop over kpoints rkn-------------------------------------------------
!#BEGIN_INTERNAL_ONLY
! ZL: Calculate sigma at each k
! Each k is independent
! In doing phonon perturbation, each "k + phonq, phonv" is a independent calculation
! phonq: electron phonon q; phonv: electron phonon v (nu)
!#END_INTERNAL_ONLY

! nkn: number of k-points on which to calculate Sigma (from sigma.inp)
  do ika=1,sig%nkn
!----------------------------
! Read wavefunctions for rkn (if sig%nkn.gt.1)
! (else it is already in wfnk)
! ZL: indkn: mapping of k-points from sigma.inp to those in kp%rk from WFN files
! get the kpoint index in the WFN file
!
    ! ==================================================================
    ! ZL: begin preparation of wfnk
    !   in ika for outer
    !   in ikn for inner 
    ! ZL: the reason for this lengthy generation of wfnk is that wfnkmpi is
    ! distributed in input.f90

    ! ika is ik for outer in sigma.inp, ikn is ik for inner in kp%rk(1:3, 1:kp%nrk)
    ikn = sig%indkn(ika)
    call timing%start(timing%wf_comm)

    SAFE_ALLOCATE(wfnkq%occ, (sig%ntband, sig%nspin))

    if(sig%nkn.gt.1) then
      ! ZL: ntband       !< number of bands in dynamical sigma summation
      !     nvband       !< number of bands in bare exchange
      nbandi=sig%ntband

      ! ZL: in input.f90:
      ! wfnkmpi stores in the order of sigma input, outer wfn
      ! wfnkqmpi stores wavefunction read in, inner wfn
      wfnk%nkpt=wfnkmpi%nkptotal(ika)
      !  ndiag_max=sig%ndiag/npools
      !  ndv = ndiag * ngk

      ! ZL: def in input.f90
      !    wfnkqmpi%nkptotal(irk) = kp%ngk(irk)
      !    wfnkqmpi%isort(1:kp%ngk(irk),irk) = isort(1:kp%ngk(irk))
      wfnk%ndv=peinf%ndiag_max*wfnk%nkpt
      ! ZL: only first ngk numbers are effective
      !     this is done in read_wavefunctions() in input.f90, with findvector
      wfnk%isrtk(1:wfnk%nkpt)=wfnkmpi%isort(1:wfnk%nkpt,ika)

      ! ZL: in input.f90:  SAFE_ALLOCATE(wfnkmpi%qk, (3,sig%nkn)), nkn: kpoints
      ! in outer wfn
      ! ZL: in read_wavefunctions():
      !     do irk=1,kp%nrk 
      !         qk(:)=kp%rk(:,irk)
      !     then others copy this qk
      !     therefore qk(1:3) stores coordinates of k-points
      qk(1:3)=wfnkmpi%qk(1:3,ika)
      ! ZL: in ../Common/typedefs.f90
      !      el(:,:,:) !< band energies (band, kpoint, spin)
      !      elda(:,:,:) !< band energies before eqp correction
      wfnk%ek(1:sig%ntband,1:sig%nspin) = wfnkmpi%el(1:sig%ntband,1:sig%nspin,ika)
      wfnk%elda(1:sig%ntband,1:sig%nspin) = wfnkmpi%elda(1:sig%ntband,1:sig%nspin,ika)
      SAFE_ALLOCATE(wfnk%zk, (wfnk%ndv,sig%nspin*kp%nspinor))
      wfnk%zk=ZERO

      ! ZL: here k loops over spin index (unbelievably!)
      do k=1,sig%nspin*kp%nspinor
#ifdef MPI
        tag=1024
        ! ZL: jj represents how many ndiag*ngk coeffs will be calculated for
        ! each pool
        if (mod(wfnk%ndv,peinf%npes/peinf%npools).eq.0) then
          jj=wfnk%ndv/(peinf%npes/peinf%npools)
        else
          jj=wfnk%ndv/(peinf%npes/peinf%npools)+1
        endif
        do j=1,peinf%npools
          dest=(j-1)*(peinf%npes/peinf%npools)
          if (peinf%inode.eq.dest) then
            g1=1
            g2=min(jj,wfnk%ndv)
            if (g2.ge.g1) then
              ! ZL: wfnk gets coeffs for each pool, stored in zk
              !     from outer wfn at ika
              wfnk%zk(g1:g2,k)=wfnkmpi%cg(g1:g2,k,ika)
            endif ! g2.ge.g1
          endif
          do i=2,peinf%npes/peinf%npools
            source=(i-1)+dest
            g1=1+(i-1)*jj
            g2=min(jj+(i-1)*jj,wfnk%ndv)
            if (g2.ge.g1) then
              if (peinf%inode.eq.source) &
                call MPI_Send(wfnkmpi%cg(1,k,ika),g2-g1+1,MPI_SCALAR,dest,tag,MPI_COMM_WORLD,mpierr)
              if (peinf%inode.eq.dest) &
                call MPI_Recv(wfnk%zk(g1,k),g2-g1+1,MPI_SCALAR,source,tag,MPI_COMM_WORLD,mpistatus,mpierr)
            endif ! g2.ge.g1
          enddo ! i=2,peinf%npes/peinf%npools
        enddo ! j=1,peinf%npools
#else
        wfnk%zk(1:wfnk%ndv,k)=wfnkmpi%cg(1:wfnk%ndv,k,ika)
        ! wfnk done with collection from wfnkmpi (outer), which is distributed
        ! ZL: outer wfn is stored in wfnk when in calculations
#endif
      enddo ! k=1,sig%nspin*kp%spinor   

#ifdef MPI
      call MPI_Bcast(wfnk%nkpt,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(wfnk%ndv,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(nbandi,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(wfnk%isrtk,gvec%ng,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(qk,3,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(wfnk%ek,nbandi*sig%nspin,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(wfnk%elda,nbandi*sig%nspin,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)

      if (peinf%npools.eq.1) then
        call MPI_Bcast(wfnk%zk,wfnk%ndv*sig%nspin*kp%nspinor,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
      else
        tag=1024
        do j=1,peinf%npools
          source=(j-1)*(peinf%npes/peinf%npools)
          do i=2,peinf%npes/peinf%npools
            dest=(i-1)+source
            if (peinf%inode.eq.source) &
              call MPI_Send(wfnk%zk,wfnk%ndv*sig%nspin*kp%nspinor,MPI_SCALAR,dest,tag,MPI_COMM_WORLD,mpierr)
            if (peinf%inode.eq.dest) &
              call MPI_Recv(wfnk%zk,wfnk%ndv*sig%nspin*kp%nspinor,MPI_SCALAR,source,tag,MPI_COMM_WORLD,mpistatus,mpierr)
          enddo
        enddo
      endif
#endif
    endif ! sig%nkn.gt.1
    call timing%stop(timing%wf_comm)
    
    if(peinf%inode.eq.0) then
      write(6,*)
      call print_dealing_with(ika, sig%nkn, kp%rk(:,ikn), 'k')
    endif
    
    ! ZL: end preparation of wfnk
    !   in ika for outer
    !   in ikn for inner 
    ! ==================================================================


!#BEGIN_INTERNAL_ONLY
    ! =====================================================================
    ! ====================== wfnk_phonq from wfnk_phonq_mpi=======================
    ! =====================================================================
    ! ZL: again, we are taking advantage of nphonq=1
    if(sig%elph .and. (sig%nphonq .ne. 1)) then
      call die("For electron phonon calculation, we only support ONE phonon_q point!", only_root_writes=.true.)
    endif

    if(sig%elph) then
      ik_phonq_idx = sig%indk_phonq(ika)
      call timing%start(timing%wf_comm)
    endif

    if(sig%elph .and. (sig%nkn.gt.1)) then
      ! ZL: ntband       !< number of bands in dynamical sigma summation
      !     nvband       !< number of bands in bare exchange
!      nbandi=sig%ntband

      wfnk_phonq%nkpt=wfnk_phonq_mpi%nkptotal(ika)
      wfnk_phonq%ndv=peinf%ndiag_max*wfnk_phonq%nkpt
      wfnk_phonq%isrtk(1:wfnk_phonq%nkpt)=wfnk_phonq_mpi%isort(1:wfnk_phonq%nkpt,ika)

      k_phonq_coord(1:3)=wfnk_phonq_mpi%qk(1:3,ika)
      wfnk_phonq%ek(1:sig%ntband,1:sig%nspin) = wfnk_phonq_mpi%el(1:sig%ntband,1:sig%nspin,ika)
      wfnk_phonq%elda(1:sig%ntband,1:sig%nspin) = wfnk_phonq_mpi%elda(1:sig%ntband,1:sig%nspin,ika)
      SAFE_ALLOCATE(wfnk_phonq%zk, (wfnk_phonq%ndv,sig%nspin*kp%nspinor))
      wfnk_phonq%zk=ZERO

      do k=1,sig%nspin*kp%nspinor
#ifdef MPI
        tag=1024
        if (mod(wfnk_phonq%ndv,peinf%npes/peinf%npools).eq.0) then
          jj=wfnk_phonq%ndv/(peinf%npes/peinf%npools)
        else
          jj=wfnk_phonq%ndv/(peinf%npes/peinf%npools)+1
        endif
        do j=1,peinf%npools
          dest=(j-1)*(peinf%npes/peinf%npools)
          if (peinf%inode.eq.dest) then
            g1=1
            g2=min(jj,wfnk_phonq%ndv)
            if (g2.ge.g1) then
              wfnk_phonq%zk(g1:g2,k)=wfnk_phonq_mpi%cg(g1:g2,k,ika)
            endif ! g2.ge.g1
          endif
          do i=2,peinf%npes/peinf%npools
            source=(i-1)+dest
            g1=1+(i-1)*jj
            g2=min(jj+(i-1)*jj,wfnk_phonq%ndv)
            if (g2.ge.g1) then
              if (peinf%inode.eq.source) &
                call MPI_Send(wfnk_phonq_mpi%cg(1,k,ika),g2-g1+1,MPI_SCALAR,dest,tag,MPI_COMM_WORLD,mpierr)
              if (peinf%inode.eq.dest) &
                call MPI_Recv(wfnk_phonq%zk(g1,k),g2-g1+1,MPI_SCALAR,source,tag,MPI_COMM_WORLD,mpistatus,mpierr)
            endif ! g2.ge.g1
          enddo ! i=2,peinf%npes/peinf%npools
        enddo ! j=1,peinf%npools
#else
        wfnk_phonq%zk(1:wfnk_phonq%ndv,k)=wfnk_phonq_mpi%cg(1:wfnk_phonq%ndv,k,ika)
#endif
      enddo ! k=1,sig%nspin*kp%spinor   

#ifdef MPI
      call MPI_Bcast(wfnk_phonq%nkpt,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(wfnk_phonq%ndv,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!      call MPI_Bcast(nbandi,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(wfnk_phonq%isrtk,gvec%ng,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(k_phonq_coord,3,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(wfnk_phonq%ek,nbandi*sig%nspin,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(wfnk_phonq%elda,nbandi*sig%nspin,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)

      if (peinf%npools.eq.1) then
        call MPI_Bcast(wfnk_phonq%zk,wfnk_phonq%ndv*sig%nspin*kp%nspinor,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
      else
        tag=1024
        do j=1,peinf%npools
          source=(j-1)*(peinf%npes/peinf%npools)
          do i=2,peinf%npes/peinf%npools
            dest=(i-1)+source
            if (peinf%inode.eq.source) &
              call MPI_Send(wfnk_phonq%zk,wfnk_phonq%ndv*sig%nspin*kp%nspinor,MPI_SCALAR,dest,tag,MPI_COMM_WORLD,mpierr)
            if (peinf%inode.eq.dest) &
              call MPI_Recv(wfnk_phonq%zk,wfnk_phonq%ndv*sig%nspin*kp%nspinor,MPI_SCALAR,source,tag,MPI_COMM_WORLD, & 
                            mpistatus,mpierr)
          enddo
        enddo
      endif
#endif
    endif ! (sig%elph .and. (sig%nkn.gt.1)) 
!#END_INTERNAL_ONLY
    call timing%stop(timing%wf_comm)
    
!#BEGIN_INTERNAL_ONLY
    if(sig%elph.and.(peinf%inode.eq.0)) then
      write(6,*)
      call print_dealing_with(ika, sig%nkn, kp%rk(:,ik_phonq_idx)+sig%indk_phonq_g0(:,ika), 'k+phonon_q')
    endif
!#END_INTERNAL_ONLY
    
    ! =====================================================================
    ! =================== done wfnk_phonq from wfnk_phonq_mpi=====================
    ! =====================================================================

    ! ZL: for DEBUG use ONLY!
    ! wfnk_phonq passes test
!    if(sig%elph) ep_debug = .true.
!    if(ep_debug) then
!      ikn = ik_phonq_idx
!      wfnk = wfnk_phonq
!    endif

!----------------------------
! Initialize Matrix Elements
! ZL: here is the start of the calculation of sigma matrix elements

    alda=0.0d0
    ! alda = vxc
    ! ZL: initialization for kih
    akih=0.0d0
!#BEGIN_INTERNAL_ONLY
    if (sig%coul_mod_flag) then
      alda2=0.0D0
    endif
!#END_INTERNAL_ONLY
    ax=ZERO
!#BEGIN_INTERNAL_ONLY
    if(sig%elph) then
      ax_ep_one=ZERO
      ax_ep_two=ZERO
      ax_ep=ZERO
      alda_ep=0.0d0
    endif
!#END_INTERNAL_ONLY
    achcor(:,:)=(0.0d0,0.0d0)
    asig_imag(:,:)=(0.0d0,0.0d0)
!#BEGIN_INTERNAL_ONLY
    if(sig%elph) then
      achcor_ep_one(:,:)=(0.0d0,0.0d0)
      asig_imag_ep_one(:,:)=(0.0d0,0.0d0)
      achcor_ep_two(:,:)=(0.0d0,0.0d0)
      asig_imag_ep_two(:,:)=(0.0d0,0.0d0)
      achcor_ep(:,:)=(0.0d0,0.0d0)
      asig_imag_ep(:,:)=(0.0d0,0.0d0)
    endif
!#END_INTERNAL_ONLY
    if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
      achtcor_n1(:) = ZERO
!#BEGIN_INTERNAL_ONLY
      if (sig%elph) then
        achtcor_n1_ep_one(:) = ZERO
        achtcor_n1_ep_two(:) = ZERO
      endif
!#END_INTERNAL_ONLY
    endif
    achcor_n1(:,:,:) = ZERO
!#BEGIN_INTERNAL_ONLY
    if (sig%elph) then
      achcor_n1_ep_one(:,:,:) = ZERO
      achcor_n1_ep_two(:,:,:) = ZERO
      achcor_n1_ep(:,:,:) = ZERO
    endif
!#END_INTERNAL_ONLY
    if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
      ! array size: SAFE_ALLOCATE(asx, (3,sig%ndiag+sig%noffdiag,sig%nspin))
      asx(:,:,:)=ZERO
      ach(:,:,:)=ZERO
      asig(:,:)=ZERO
      ach_n1(:,:,:)=ZERO
!#BEGIN_INTERNAL_ONLY
      if (sig%elph) then
        asx_ep_one(:,:,:)=ZERO
        ach_ep_one(:,:,:)=ZERO
        asig_ep_one(:,:)=ZERO
        ach_n1_ep_one(:,:,:)=ZERO

        asx_ep_two(:,:,:)=ZERO
        ach_ep_two(:,:,:)=ZERO
        asig_ep_two(:,:)=ZERO
        ach_n1_ep_two(:,:,:)=ZERO

        asx_ep(:,:,:)=ZERO
        ach_ep(:,:,:)=ZERO
        asig_ep(:,:)=ZERO
        ach_n1_ep(:,:,:)=ZERO
      endif
!#END_INTERNAL_ONLY
    endif
    if (sig%freq_dep.eq.2) then
      asxDyn(:,:,:)=(0.0d0,0.0d0)
      achDyn(:,:,:)=(0.0d0,0.0d0)
      achDyn_cor(:,:,:)=(0.0d0,0.0d0)
      achDyn_corb(:,:,:)=(0.0d0,0.0d0)
      ach2Dyn(:,:,:)=(0.0d0,0.0d0)
      asigDyn(:,:)=(0.0d0,0.0d0)
      achD_n1(:,:,:)=(0.0d0,0.0d0)
    endif
    
!----------------------------
! Read matrix elements of Vxc from file vxc.dat
! or compute them on the fly from Vxc potential

    call timing%start(timing%vxc)
    if(sig%use_vxcdat .and. .not. sig%is_EXX) then
      if(peinf%inode == 0) write(6,*) 'Reading vxc.dat'
      call open_file(120,file='vxc.dat',form='formatted',status='old')
      qk(:)=INF
      ierr=0
      do while (ierr.eq.0)
        call read_matrix_elements_type(120, ierr, qk, sig, alda)
        if (all(abs(kp%rk(1:3,ikn)-qk(1:3)) .lt. TOL_Small)) exit
      enddo
      call close_file(120)
          
! Check k-point

      if(any(abs(kp%rk(1:3,ikn)-qk(1:3)) .ge. TOL_Small)) then
        call die('cannot find k-point in vxc.dat', only_root_writes = .true.)
      endif

! ZL: IMPORTANT: here is where the units of diag and offdiag become different
! Divide by ryd for diag
! this will be undone by shift_energy routines later

      do s2=1,sig%nspin
        do in=1,sig%ndiag
        alda(in,s2) = alda(in,s2)/ryd  ! ZL: if directly read from vxc.dat
        enddo
      enddo

    elseif (.not.sig%sigma_correction .and. .not. sig%is_EXX .and. .not.sig%use_kihdat) then ! not using vxc.dat
 
      ! ZL: here we calculate mtxel_vxc     
      call logit('Calling mtxel_vxc')
      call mtxel_vxc(kp,gvec,sig,wfnk,wfnkoff,alda,1) ! ZL: 1 for VXC

    elseif (sig%use_kihdat) then
      ! ZL: conflict check between kihdat and vxcdat has been done in inread.f90
      if(peinf%inode == 0) write(6,*) 'Reading kih.dat'
      call open_file(876,file='kih.dat',form='formatted',status='old')
      qk(:)=INF
      ierr=0
      do while (ierr.eq.0)
        call read_matrix_elements_type(876, ierr, qk, sig, akih)
        if (all(abs(kp%rk(1:3,ikn)-qk(1:3)) .lt. TOL_Small)) exit
      enddo
      call close_file(876)
          
! Check k-point

      if(any(abs(kp%rk(1:3,ikn)-qk(1:3)) .ge. TOL_Small)) then
        call die('cannot find k-point in vxc.dat', only_root_writes = .true.)
      endif

! ZL: IMPORTANT: again the units of diag and offdiag become different

      do s2=1,sig%nspin
        do in=1,sig%ndiag
        akih(in,s2) = akih(in,s2)/ryd  ! ZL: if directly read from kih.dat
        enddo
      enddo

    endif

!#BEGIN_INTERNAL_ONLY
    if (sig%elph) then
      call logit('Calling mtxel_vxc for dVXC')
      call mtxel_vxc(kp,gvec,sig,wfnk,wfnkoff,alda_ep,3,wfnk_phonq,wfnk_phonq_off) ! ZL: 3 for dVXC, with two more arguments
    endif
!#END_INTERNAL_ONLY

!#BEGIN_INTERNAL_ONLY
    ! This if for the hybrid functional calculations (one shot)
    ! If the matrix elements are already there then just read them 
    ! in -- why the hell are you doing this calculation then???
    if(sig%coul_mod_flag) then
      if (sig%elph) then
        call die("Electron-phonon does not support coul_mod_flag.")
      endif
      if(sig%use_vxc2dat) then
        if(peinf%inode == 0) write(6,*) 'Reading vxc2.dat'
        call open_file(121,file='vxc2.dat',form='formatted',status='old')
        qk(:)=INF
        ierr=0
        do while (ierr.eq.0)
          call read_matrix_elements_type(121, ierr, qk, sig, ax)
          if (all(abs(kp%rk(1:3,ikn)-qk(1:3)) .lt. TOL_Small)) exit
        enddo
        call close_file(121)
          
! Check k-point

        if(any(abs(kp%rk(1:3,ikn)-qk(1:3)) .ge. TOL_Small)) then
          call die('cannot find k-point in vxc2.dat', only_root_writes = .true.)
        endif

! Divide by ryd for diag
! this will be undone by shift_energy routines later

        do s2=1,sig%nspin
          do in=1,sig%ndiag
            ax(in,s2) = ax(in,s2)/ryd
          enddo
        enddo

      else ! not using vxc2.dat
      
        call logit('Calling mtxel_vxc for VXC2')
        call mtxel_vxc(kp,gvec,sig,wfnk,wfnkoff,alda2,2) ! ZL: 2 for VXC2

      endif

    endif
!#END_INTERNAL_ONLY

    call timing%stop(timing%vxc)
    
#ifdef CPLX
    if (any(abs(IMAG(alda(1:sig%ndiag,:)))*ryd > 1.0d-4)) imagvxcflag=.true.
    ! ZL: for kih, everything is already complex and imagkihflag has been set true
#endif

!----------------------------
! Read ax from existing data
! ZL: generally it is directly calculated

    if(sig%use_xdat .and. xflag .and. (.not. sig%coul_mod_flag)) then
      if(sig%elph) then
        call die('EP calculation cannot read x.dat')
      endif

      if(peinf%inode == 0) write(6,*) 'Reading x.dat'
      call open_file(119,file='x.dat',form='formatted',status='old')
      qk(:)=INF
      ierr=0
      do while (ierr.eq.0)
        call read_matrix_elements_type(119, ierr, qk, sig, ax)
        if (all(abs(kp%rk(1:3,ikn)-qk(1:3)) .lt. TOL_Small)) exit
      enddo
      ax(:,:) = ax(:,:) * sig%xfrac
      call close_file(119)
      
! Check k-point

      if(any(abs(kp%rk(1:3,ikn)-qk(1:3)) .ge. TOL_Small)) then
        call die('cannot find k-point in x.dat', only_root_writes = .true.)
      endif

! Divide by ryd for diag

      do s2=1,sig%nspin
        do in=1,sig%ndiag
          ax(in,s2) = ax(in,s2)/ryd
        enddo
      enddo

#ifdef CPLX
      if (any(abs(IMAG(ax(1:sig%ndiag,:)))*ryd > 1.0d-4)) imagxflag=.true.
#endif

    endif ! using x.dat

!----------------------------
! Find subgroup which leaves kn invariant
! Indices of group operations in subgroup stored in array indsub
! stored in structure syms

! ZL: for phonon purpose, there are k, q, phonq points, may ignore this for now
! we now disable this for EP
!#BEGIN_INTERNAL_ONLY
    if (sig%elph) sig%qgridsym = .false. ! ZL added
!#END_INTERNAL_ONLY
    ! ZL: note that this has already been dealt with in input.f90
    ! ZL: here get qk(:), and qk(:) is the coord of ikn, inner wfn kpoint
    qk(:) = kp%rk(:,ikn)
    call timing%start(timing%subgrp)
    if (sig%qgridsym) then
      call subgrp(qk,syms)
    else
      syms%ntranq=1
      syms%indsub(1)=1
      syms%kgzero(1:3,1)=0
    endif
    call timing%stop(timing%subgrp)

!----------------------------
! Reduce qpoints with respect to group of kn
! Keep track of number of q-points equivalent
! to give qpoint in irr-bz neq(irq)
!
! Keep track of q-point in the set given in epsmat
! RQ = R(Q) + G0 with transformation R and umklapp G0
! In order to unfold the inverse dielectric matrix from Q to RQ
!
! The q-points are the points for which we have epsilon
! (should be the unshifted irreducible grid)

    ! ZL: grid%nr  !< number in reduced zone
    !     grid%nf  !< number in full zone
    SAFE_ALLOCATE(indrq, (gr%nf))  ! full zone
    SAFE_ALLOCATE(neq, (gr%nf))  ! full zone
    SAFE_ALLOCATE(itnrq, (gr%nf))  ! full zone
    SAFE_ALLOCATE(rq, (3,gr%nf))  ! full zone
    SAFE_ALLOCATE(kg0, (3,gr%nf))  ! full zone
    call timing%start(timing%irrbz)
    call irrbz(syms,gr%nf,gr%f,nrq,neq,indrq,rq,sig%nq,sig%qpt,itnrq,kg0,nfix=sig%nq0)
    call timing%stop(timing%irrbz)
    

   if ( mtxel_algo /= CPU_ALGO ) then
      call setup_FFT_sizes(gvec%FFTgrid,Nfft,dummy_scale)
      accel_mtxel_sig%mtxel_band_block_size = sig%accel_mtxel_band_block_size
      accel_mtxel_sig%mtxel_band_block_size = MIN(accel_mtxel_sig%mtxel_band_block_size, peinf%ntband_node)
      accel_mtxel_sig%mtxel_band_block_size = MAX(1, accel_mtxel_sig%mtxel_band_block_size)
      call allocate_accel_mtxel_sig ( Nfft, kp%ngkmax, kp%ngkmax, gvec%ng, kp%ngkmax, mtxel_algo )
      if (peinf%inode == 0) then
        write(6,'(1x,A)')
        write(6,'(1x,A)') 'Using GPU support via OpenACC/OMP-Target implementation'
        write(6,'(1x,A,I6)')   'Band block size in MTXEL:', accel_mtxel_sig%mtxel_band_block_size
        if ( sigma_gpp_algo == OMP_TARGET_ALGO .and. (sig%freq_dep==1 .or. sig%freq_dep==3) ) then
          write(6,'(1x,A)') 'Using GPU support via OMP-Target for Sigma-GPP'
          write(6,'(1x,A,I6)') 'Band block size:',  sig%gpp_band_block_size
          write(6,'(1x,A,I6)') 'G-vec block size:', sig%gpp_ig_block_size
        end if
      end if
   end if

!!---------- Loop over k-points in irr bz with respect to kn (rq) --------------

#ifdef MPI
   call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif
   call timing%start(timing%eqp_loop_total)

    irq_min = 1
! FHJ: For subsampling calculations, we calculate the self energy for the first
! q-point (irq_==0, irq==1) using the input epsilon cutoff (ecuts), but only a
! correction to the self energy for the other q-points, which is calculated
! using a smaller cuttof. The way this is implemented is by looping twice over
! the first q-point: first (irq_==0, irq==1) using a weight of 1, and then
! (irq_==1, irq==1) using a weight of subweight(1)-1. All other subsampled
! q-points are calculated with the weight subweight(irq).
    if (sig%subsample) irq_min = 0
    call progress_init(prog_info, 'calculating Sigma', 'block', &
      (nrq-irq_min+1)*(peinf%ndiag_max+peinf%noffdiag_max)*sig%nspin)

    !======================================================================
    ! ZL: here is the start loop over q points
    do irq_ = irq_min, nrq
      irq = irq_
      is_subq = .false.
      ! FHJ: is this a q-point used in the subsampling of the voronoi cell for q==0?
      if (sig%subsample.and.irq<=sig%nq0) then
        is_subq = .true.
        if (irq==0) irq=1
      endif

      if (peinf%verb_debug .and. peinf%inode.eq.0) then
        write(6,60) irq,nrq
        write(6,*)
      endif
60    format(/,3x,'qpoint',i5,' out of ',i5)

! Compute energies: |q+g|**2
! ZL: stored as ekin(ig)
      if (is_subq) then
        call kinetic_energies(gvec, crys%bdot, ekin)
      else
        call kinetic_energies(gvec, crys%bdot, ekin, qvec = rq(1:3, irq))
      endif

! Sort ekin in ascending order for this q
! The indices are placed in array isrtrq
      !ZL: Note that ekin of |q+g|**2 is used
      call sortrx(gvec%ng, ekin, isrtrq, gvec = gvec%components)
      if ((sig%freq_dep.eq.0.or.sig%exact_ch.eq.1).and.irq_==irq_min) then
        isrtrqi=0
        do j=1,gvec%ng
          if (isrtrq(j).ge.1.and.isrtrq(j).le.gvec%ng) &
            isrtrqi(isrtrq(j))=j
        enddo
      endif

      !=========================================================
      ! ZL: calculate cutoffs, unchanged
! Compute cutoff in sums over G,G`
! definition:   real(DP) :: ecutb    !< energy cutoff of bare coulomb interaction in Ry
! definition:   real(DP) :: ecuts    !< energy cutoff of screened coulomb interaction in Ry
      ncoulb = gcutoff(gvec%ng, ekin, isrtrq, sig%ecutb)
      ncouls = gcutoff(gvec%ng, ekin, isrtrq, sig%ecuts)
      if (is_subq.and.irq_>0) then !FHJ: see comment before loop over irq_
        ncouls = min(ncouls, gcutoff(gvec%ng, ekin, isrtrq, subcut))
      endif
      ! ncoul is the maximum possible value
      ncoul = max(ncouls,ncoulb)

      if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then
        if(irq_==irq_min) ncoulch = ncoul
      else
        ncoulch = 0
      endif

      ! ZL: important
      ! this condition is assumed later in the code (e.g. wpeff array sizes), so we must check it
      if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
        if (ncouls.gt.neps) then
          write(tmpstr,'(a,i6,a,i6)')&
            "screened Coulomb cutoff is bigger than epsilon cutoff"//CHAR(10)&
            //' ncouls = ',ncouls,' neps = ',neps
          ! NB: CHAR(10) is the carriage return.
          call die(tmpstr, only_root_writes = .true.)
        endif
      endif

      ! ZL: done cutoffs, unchanged
      !=========================================================
!----------------------------
! Allocate arrays for q-point rq(:,irq)
      ! ZL: V(G)
      SAFE_ALLOCATE(vcoul, (ncoul))
      ! ZL: aqs - matrix elements M(N_G, N_band) now for a given k and q point (in the loops)
      SAFE_ALLOCATE(aqs, (ncoul,peinf%ntband_max))
      if (sig%noffdiag.gt.0) then 
        SAFE_ALLOCATE(aqsaug, (ncoul,peinf%ntband_max,sig%ndiag,sig%nspin))
      endif
      ! ZL: allocating for EP
!#BEGIN_INTERNAL_ONLY
      if (sig%elph) then
        SAFE_ALLOCATE(aqs_ep_dq_onel, (ncoul,peinf%ntband_max))
        SAFE_ALLOCATE(aqs_phonq_twol, (ncoul,peinf%ntband_max))
        SAFE_ALLOCATE(aqs_ep_dmq_twor, (ncoul,peinf%ntband_max))
        if (sig%noffdiag.gt.0) then 
         SAFE_ALLOCATE(aqsaug_ep_dq_onel, (ncoul,peinf%ntband_max,sig%ndiag,sig%nspin))
         SAFE_ALLOCATE(aqsaug_phonq_twol, (ncoul,peinf%ntband_max,sig%ndiag,sig%nspin))
         SAFE_ALLOCATE(aqsaug_ep_dmq_twor, (ncoul,peinf%ntband_max,sig%ndiag,sig%nspin))
        endif
      endif
!#END_INTERNAL_ONLY
      ! ZL: for GPP, we dont consider this if-statement
      if ((sig%freq_dep.eq.0.or.sig%exact_ch.eq.1).and.irq_==irq_min) then
        ! ZL: EP implemented for GPP only for now
!#BEGIN_INTERNAL_ONLY
        if (sig%elph) then
          ! ZL: we do not even allocate aqsch_ep_one and aqsch_ep_two for EP for now
          call die('Electron phonon calculation implemented for GPP only, for now.')
        endif
!#END_INTERNAL_ONLY
        SAFE_ALLOCATE(aqsch, (ncoulch))
        if (nrq.gt.1) then
          SAFE_ALLOCATE(aqsaugchd, (ncoulch,peinf%ndiag_max,sig%nspin))
          if (sig%noffdiag.gt.0) then
            SAFE_ALLOCATE(aqsaugcho, (ncoulch,peinf%noffdiag_max,sig%nspin))
          end if
        endif
      endif
      ! ZL: nm is index for q
      nm = indrq(irq)


!!!------- Calculate Vcoul -----------------------------------------------------
      call timing%start(timing%vcoul)
      vq=rq(:,irq)
      qlen = sqrt(DOT_PRODUCT(vq,MATMUL(crys%bdot,vq)))
      if(sig%freq_dep /= -1) call check_screening_trunc(sig%icutv,sig%iscreen,sig%q0vec,crys%bdot)
      iparallel=1
      avgcut = sig%avgcut
      ! ZL: MC subsampling
      ! FHJ: Disable MC averages when we perform the subsampling of the voronoi cell
      if (is_subq) avgcut = 0d0
      if (.not. sig%coul_mod_flag) then
        call vcoul_generator(sig%icutv,sig%truncval,gvec,crys%bdot,crys%celvol, &
          gr%nf-sig%nq0+1,ncoul,isrtrq,sig%iscreen,vq,sig%q0vec,vcoul, &
          sig%iwritecoul,iparallel,avgcut,oneoverq,sig%qgrid,epshead, &
          work_scell,sig%averagew,sig%wcoul0)
!#BEGIN_INTERNAL_ONLY
      else
          ! for hybrid
        call vcoul_generator(sig%icutv,sig%truncval,gvec,crys%bdot,crys%celvol, &
          gr%nf-sig%nq0+1,ncoul,isrtrq,sig%iscreen,vq,sig%q0vec,vcoul, &
          sig%iwritecoul,iparallel,avgcut,oneoverq,sig%qgrid,epshead, &
          work_scell,sig%averagew,sig%wcoul0,coulomb_mod=sig%coulomb_mod)
!#END_INTERNAL_ONLY
      endif

      fact = 1D0/(dble(gr%nf-sig%nq0+1)*crys%celvol)
      coulfact = 8D0*PI_D/(dble(gr%nf-sig%nq0+1)*crys%celvol)
      ! FHJ: we won`t actually use sig%wcoul0 when sig%subsample==.true.
      ! Rescale vcoul according to the subsample weights
      if (is_subq) then
        weight = sig%subweights(nm)
        if (irq_==0) then ! FHJ: see comment before loop over irq_
          weight = 1d0
        elseif (irq_==1) then
          weight = weight - 1d0
        endif
        fact = fact*weight
        coulfact = coulfact*weight
        if (irq/=nm) then
          if (peinf%inode==0) write(0,'(a,i0,a,i0)') 'ERROR: irq=',irq,' nm=',nm
          call die('Inconsistent indices for a sub-sampled q-point', only_root_writes=.true.)
        endif
      endif
      ! ZL: now we have V(G)
      do ig = 1, ncoul
        vcoul(ig)=fact*vcoul(ig)
      enddo
      if (ika.eq.1.and.irq_==irq_min) then
        sig%wcoul0 = sig%wcoul0 * fact
      endif
      call timing%stop(timing%vcoul)
      

!!!------- Read inverse dielectric matrix for q-point rq(:,irq) ----------------
! ZL: for EP, epsinv is unchanged, which determines wtilde and Omega,
!     and they are for the same q point (or p point in the formalism)
      call timing%start(timing%epsread)
  
      if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
        
        q0len = sqrt(DOT_PRODUCT(sig%q0vec,MATMUL(crys%bdot,sig%q0vec)))
        
!----------------------------
! Processor 0 read eps^-1 for this q and broadcast to others
!
! Note: the total number of g-vectors used during
! computation of the inverse dielectric matrix
! may be different than in present calculation
! although the sets must coincide for small g

        ngq = gvec%ng
        nmtx = epsmpi%nmtx(nm)

        if (peinf%inode .eq. 0) then
          ! ZL: isrtq now represents the order of q-points in eps
          isrtq(:) = epsmpi%isrtq(:,nm)
          isrtqi(:) = epsmpi%isrtqi(:,nm)
        endif
#ifdef MPI
        call MPI_Bcast(isrtq, ngq, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(isrtqi, ngq, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
        qk(:) = epsmpi%qk(:,nm)

! Find g=0 in main gvec list and eps gvector list

        ! write explicitly to avoid possible warning about array temporary
        nullvec(1:3) = 0
        call findvector(iout,nullvec,gvec)
        iout = isrtqi(iout)
          
        if (peinf%verb_debug .and. peinf%inode==0) write(6,*) 'Reading Eps Back'

! JRD: XXX This is sort of a waste of memory... Can we use pointers for this sort of thing?
        
        if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
          eps(:,:)=epsmpi%eps(:,:,nm)
        else ! sig%freq_dep.eq.2 .or. 3
          epsR(:,:,:)=epsmpi%epsR(:,:,:,nm)
          if (sig%need_advanced) then
            epsA(:,:,:)=epsmpi%epsA(:,:,:,nm)
          endif
        endif
        
! CHP: By first distributing and then doing the wing fix, we can let
!      all the processors work together, thus, saving some time.

! ZL: TODO q0->0 behavior not determined yet, to be fixed or double checked
!-------------------------------------------------------------------------------
! Fix wing divergence for semiconductors and graphene

! This should really be done for all "|q+G| < avgcut" - but for now,
! it is done if "|q| < avgcut and G=0"
 
        ! make sure corrections for wings are zero at the beginning     
        if (sig%do_sigma_subspace) then
           sig%epssub%eps_wings_correction_cols(:,:) = (0.0d0, 0.0d0)
           sig%epssub%eps_wings_correction_rows(:,:) = (0.0d0, 0.0d0)
           ! save actual nm
           sig%epssub%actual_nm = nm
        end if

        if (.not.sig%subsample) then
        ngpown_q = NUMROC(nmtx, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)

          if (sig%do_sigma_subspace) then
            ! initialize wings corrections
            ipe_wing = INDXG2P( iout, epsmpi%nb, iproc_dum, 0, peinf%npes_pool)
            if (ipe_wing == peinf%pool_rank) then
              ! make a copy of the unmodified wing
              my_pos_loc_wing = INDXG2L(iout, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
              sig%epssub%eps_wings_correction_cols(1:nmtx,1:sig%nFreq) = &
              sig%epssub%eps_wings_cols(1:nmtx,1:sig%nFreq,nm)
            end if
            do igp_loc = 1, ngpown_q
               igp = INDXL2G(igp_loc, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
               sig%epssub%eps_wings_correction_rows(igp,1:sig%nFreq) = &
               sig%epssub%eps_wings_rows(igp,1:sig%nFreq,nm)
            end do
            ! MDB calculate correction factors for fixing wings (here we try 
            ! to avoid copying stuff from other modules even if less efficient)
            SAFE_ALLOCATE(epsR_corrections, (neps,epsmpi%ngpown,1))
            epsR_corrections = (1.0,0.0)
          end if

! fix wings here
        do igp_loc = 1, ngpown_q

          igp = INDXL2G(igp_loc, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
          if (nm .eq. 1) then
            q0flag=.true.
            if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
              call fixwings(vcoul(1),sig%wcoul0,eps(:,igp_loc),sig%icutv, &
                sig%iscreen,igp,nmtx,iout,q0len,oneoverq,fact,q0flag,sig%averagew,crys%bdot)
            else ! sig%freq_dep.eq.2 .or. 3
! JRD XXX bad locality
              call fixwings_dyn(vcoul(1),epsR(1:nmtx,igp_loc,:), &
                sig%icutv,sig%iscreen,igp,nfreq_fixwings,nmtx,iout,q0len,oneoverq,fact,q0flag,crys%bdot)
              if (sig%need_advanced) then
                call fixwings_dyn(vcoul(1),epsA(1:nmtx,igp_loc,:), &
                  sig%icutv,sig%iscreen,igp,nfreq_fixwings,nmtx,iout,q0len,oneoverq,fact,q0flag,crys%bdot)
              endif
              if (sig%do_sigma_subspace) then
                call fixwings_dyn(vcoul(1),epsR_corrections(1:nmtx,igp_loc,:), &
                  sig%icutv,sig%iscreen,igp,1,nmtx,iout,q0len,oneoverq,fact,q0flag,crys%bdot)
              end if
            endif
          else if (qlen**2 .lt. sig%avgcut) then
            q0flag=.false.
            if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
              call fixwings(vcoul(1),sig%wcoul0,eps(:,igp_loc),sig%icutv, &
                sig%iscreen,igp,nmtx,iout,qlen,oneoverq,fact,q0flag,sig%averagew,crys%bdot)
            else ! sig%freq_dep.eq.2 .or. 3
! JRD XXX bad locality
              call fixwings_dyn(vcoul(1),epsR(1:nmtx,igp_loc,:), &
                sig%icutv,sig%iscreen,igp,nfreq_fixwings,nmtx,iout,qlen,oneoverq,fact,q0flag,crys%bdot)
              if (sig%need_advanced) then
                call fixwings_dyn(vcoul(1),epsA(1:nmtx,igp_loc,:), &
                  sig%icutv,sig%iscreen,igp,nfreq_fixwings,nmtx,iout,qlen,oneoverq,fact,q0flag,crys%bdot)
              endif
              if (sig%do_sigma_subspace) then
                 call fixwings_dyn(vcoul(1),epsR_corrections(1:nmtx,igp_loc,:), &
                   sig%icutv,sig%iscreen,igp,1,nmtx,iout,qlen,oneoverq,fact,q0flag,crys%bdot)
              end if
            endif
          endif

        enddo ! igp

          if (sig%do_sigma_subspace) then
            ! calculte corrections
            ipe_wing = INDXG2P( iout, epsmpi%nb, iproc_dum, 0, peinf%npes_pool)
            if (ipe_wing == peinf%pool_rank) then
              my_pos_loc_wing = INDXG2L(iout, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
              do ifreq = 1, sig%nFreq
                sig%epssub%eps_wings_correction_cols(1:nmtx,ifreq) = &
                sig%epssub%eps_wings_correction_cols(1:nmtx,ifreq) * &
                (epsR_corrections(1:nmtx,my_pos_loc_wing,1) - 1.0d0)
              end do
              ! broadcast to your fellow
#ifdef MPI
              if (peinf%my_pool/=-1) then
                call MPI_Bcast(sig%epssub%eps_wings_correction_cols(1:Neps,1:sig%nFreq),&
                               Neps * sig%nFreq, MPI_COMPLEX_DPC, ipe_wing, peinf%pool_comm, mpierr)
              end if
#endif
            else
              ! receive the wing / col
#ifdef MPI
              if (peinf%my_pool/=-1) then
                call MPI_Bcast(sig%epssub%eps_wings_correction_cols(1:Neps,1:sig%nFreq),&
                               Neps * sig%nFreq, MPI_COMPLEX_DPC, ipe_wing, peinf%pool_comm, mpierr)
              end if
#endif
            end if
            ! wing row
            do igp_loc = 1, ngpown_q
              igp = INDXL2G(igp_loc, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
              sig%epssub%eps_wings_correction_rows(igp,1:sig%nFreq) = &
              sig%epssub%eps_wings_correction_rows(igp,1:sig%nFreq)*(epsR_corrections(iout,igp_loc,1) - 1.0d0)
              if(igp == iout) then
                ! here compensate for the 1.0 on diagonal 
                sig%epssub%eps_wings_correction_rows(igp,1:sig%nFreq) = &
                sig%epssub%eps_wings_correction_rows(igp,1:sig%nFreq) + epsR_corrections(iout,igp_loc,1)
              end if
            end do
            ! sum up
#ifdef MPI
            if (peinf%my_pool/=-1) then
              call MPI_Allreduce(MPI_IN_PLACE, sig%epssub%eps_wings_correction_rows(1:Neps,1:sig%nFreq), &
                                 Neps * sig%nFreq, MPI_COMPLEX_DPC, MPI_SUM, peinf%pool_comm, mpierr)
            end if
#endif
            SAFE_DEALLOCATE_P(epsR_corrections)
          end if ! subspace

        endif !.not.sig%subsample

        call logit('Read eps from memory')
       
      endif  ! sig%freq_dep
      
      if (sig%freq_dep.eq.-1) then
        ngq=0
        nmtx=0
        qk(:)=sig%qpt(:,nm)
      endif
      
      call timing%stop(timing%epsread)
      
      if (peinf%verb_debug .and. peinf%inode==0) then
        write(6,'(3(a,i10))') 'nmtx =',nmtx,' ncouls =',ncouls
        write(6,*)
      endif
      if(nmtx.gt.neps) then
        call die('nmtx.gt.neps')
      endif
        
! Check q-vector

      if(any(abs(sig%qpt(1:3, nm) - qk(1:3)) .gt. TOL_SMALL)) then
        write(0,*) peinf%inode,sig%qpt(:,nm),qk(:)
        call die('q-vector check wrong')
      endif
        
      if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
        if(ncouls .gt. nmtx) ncouls = nmtx
      endif

      if (peinf%verb_debug .and. peinf%inode.eq.0) then
        if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
          write(6,100) (rq(i,irq),i=1,3),ncouls,dble(eps(1,1))
        endif
        if (sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
          write(6,100) (rq(i,irq),i=1,3),ncouls,dble(epsR(1,1,1))
        endif
      endif
100   format(3x,'q=',3f8.5,2x,'n=',i6,2x,'head of epsilon inverse =',f12.6,/)

! Map g-vectors required for eps**(-1)(r(q)) to those
! for the known matrix eps**(-1)(q) and calculate phases

      itran = itnrq(irq)
      kg(:) = kg0(:,irq)
      call timing%start(timing%gmap)
      call logit('Calling gmap')
      if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
        ind(:)=0
        indinv(:)=0
        ph(:)=ZERO
        call gmap(gvec,syms,ncouls,itran,kg,isrtrq,isrtqi,ind,ph, &
         sig%die_outside_sphere) ! TROUBLE
        do ig = 1, gvec%ng
          if (ind(ig) .gt. 0 .and. ind(ig) .le. gvec%ng) then
            indinv(ind(ig))=ig
          endif
        enddo
      endif
      call timing%stop(timing%gmap)
      

!!!------- Done reading inverse dielectric matrix for q-point rq(:,irq) --------


!=================================================================
!ZL: generate wfnkq from wfnkqmpi (WFN read in) for rkq = rk - rq
!    the format of wfnkq is wfnkqinfo

!--------------------
! Generate needed wavefunctions for rkq = rkn - rq
! stored in derived type wfnkq

      rkq(1:3) = kp%rk(1:3, ikn) - rq(1:3, irq)
      ! FHJ: when we perform subsampling, we use the actual q-vector for the
      ! Coulomb potential, but we pretend we have q=0 for the matrix elements.
      if (is_subq) rkq(1:3) = kp%rk(1:3, ikn)
      call timing%start(timing%genwf)
      call logit('Calling genwf')
      if (.not.is_subq.or.irq_==0) then ! FHJ: subsampling uses the same WFNs and matrix elements
        call genwf_mpi(rkq, syms, gvec, crys, kp, sig, wfnkq, wfnkqmpi)
      endif
      wfnkq%occ(:,:) = sig%occ(:,wfnkq%kmq_pt,:,1)
      call timing%stop(timing%genwf)
!ZL: done generating wfnkq
!=================================================================

!=================================================================
!#BEGIN_INTERNAL_ONLY
!ZL: generate ep_dq_wfnkq from ep_dq_wfnkqmpi (dWFN read in) for rkq_ep = rk - rq
!    the format of ep_dq_wfnkq is wfnkqinfo
!ZL: Also, generate wfnkq_phonq, ep_dmq_wfnkq_phonq
!--------------------
! Generate needed wavefunctions for rkq_ep = rkn - rq
! stored in derived type ep_dq_wfnkq
!--------------------
! ZL: for EP, this is where we should modify
!     we need:
!       Delta_phonq psi(rkn - rq) for the first term
!       Delta_{-phonq} psi(rkn - rq + phonq) for the second term
!#END_INTERNAL_ONLY

!#BEGIN_INTERNAL_ONLY
      if (sig%elph) then
        ! First term
        ! ep_dq_wfnkq
        rkq_ep(1:3) = kp%rk(1:3, ikn) - rq(1:3, irq)
        ! ZL: not sure if EP is compatible with subsampling yet
        if (is_subq) call die('Electron phonon does not support subsampling for now.')
        call logit('Calling genwf for ep_dq_wfnkq')
        if (.not.is_subq.or.irq_==0) then ! FHJ: subsampling uses the same WFNs and matrix elements
          ! ZL: use the same sig, but different syms
          ! We dont check norms for dWFN
          ! Also we dont use symmetries for dWFN, because it carries k+phonq.
          ! Unfolding is very different. So we completely do not use it.
          ! This requires ep_syms%ntran = 1, which we have checked already
          ! Actually, we dont use any rotations, but we need extra phase exp(iGr)
          check_norms_save = peinf%check_norms
          peinf%check_norms = .false.
          call genwf_mpi(rkq_ep, ep_syms, ep_gvec, ep_crys, ep_kp, sig, &
            ep_dq_wfnkq, ep_dq_wfnkqmpi)
          peinf%check_norms = check_norms_save
        endif

        ! Second term
        ! (1) wfnkq_phonq
        ! need wfnkq(-k-phonq+p), then c.c. the coeffs and gvecs
        ! NOTE: system MUST have time-reversal symmetry (TRS)
        rk_phonq(1:3) = kp%rk(1:3, ikn) + sig%phonq(1:3)
        r_mkq_mphonq(1:3) = rq(1:3, irq) - rk_phonq(1:3)
        call logit('Calling genwf for wfnkq_phonq')
        if (.not.is_subq.or.irq_==0) then
          ! TRS gauge fixing, to match ep_dmq_wfnkq_phonq
          ! TRS will generate arbitrary phase, but this can be fixed by using TRS for both inner wfn/dwfn
          call genwf_mpi(r_mkq_mphonq, syms, gvec, crys, kp, sig, wfnkq_phonq, wfnkqmpi)
          call time_reversal(wfnkq_phonq, gvec)
        endif
        ! (2) ep_dmq_wfnkq_phonq
        ! need ep_dq_wfnkq(-k-phonq+p), then c.c. the coeffs and gvecs
        call logit('Calling genwf for ep_dmq_wfnkq_phonq')
        if (.not.is_subq.or.irq_==0) then
          check_norms_save = peinf%check_norms
          peinf%check_norms = .false.
          call genwf_mpi(r_mkq_mphonq, ep_syms, ep_gvec, ep_crys, ep_kp, sig, &
            ep_dmq_wfnkq_phonq, ep_dq_wfnkqmpi)
          peinf%check_norms = check_norms_save
          call time_reversal(ep_dmq_wfnkq_phonq, ep_gvec)
        endif
      endif ! sig%elph for genwf_mpi() for ep_dq_wfnkq
!ZL: done generating wfnkq
!#END_INTERNAL_ONLY

!#BEGIN_INTERNAL_ONLY
!==========================================================================
! ZL: for EP, at this step, wfnk_phonq and ep_dq_wfnkq should be ready for
! matrix elements aqs_ep_dq_onel calculation
!==========================================================================
!#END_INTERNAL_ONLY

!!-------- Loop over spins for which Sigma is computed -------------------------
! ZL: Now get into core part of sigma
      do ispin=1,sig%nspin

!!-------- Loop over bands for which diag Sigma is computed --------------------
! Bands are relabelled according to sig%diag(1:sig%ndiag)
! Loop is parallelized according to peinf%index_diag(1:peinf%ndiag_max)

! ZL: Loop over diagonal bands, in parallel
        do in=1,peinf%ndiag_max
          
          call progress_step(prog_info)
          if (peinf%verb_debug .and. peinf%inode.eq.0) then
            if (peinf%npools.eq.1) then
              write(6,999) peinf%index_diag(in)+(ispin-1)*sig%ndiag, &
                peinf%ndiag_max*peinf%npools*sig%nspin
            else
              write(6,997) peinf%index_diag(in)+(ispin-1)*sig%ndiag, &
                peinf%index_diag(in)+(ispin-1)*sig%ndiag+peinf%npools-1, &
                peinf%ndiag_max*peinf%npools*sig%nspin
            endif
          endif
999       format(1x,"Computing Sigma diag",i4,1x,"of",i4)
997       format(1x,"Computing Sigma diag",i4,1x,"to",i4,1x,"of",i4)
          write(tmpstr,*) 'Working on band ', sig%diag(peinf%index_diag(in)), ' 1st pool'
          call logit(tmpstr)

!---------------------
! Compute planewave matrix elements of g <n,k|exp{i(q+g).r}|n1,k-q>
! Note: wfnk%zk array keeps only the bands specified in sig%diag(:)
! Must keep track of the right band label

          call logit('Calling mtxel')
          call timing%start(timing%mtxel)
          ! ZL: role of isrtrq
          ! index array for g-vectors in <nk|exp(i(q+g).r)|n1k-q>
          ! sorted with |q+g|**2
          ! in ../Common/fftw_inc.f90
          !   isrtrq is only used to find where the ncoul cutoff is
          !   when getting aqs data from already-calculated FFT
          !
          ! ZL: wfnk, wfnkq will use their own isrt in putting themselves in FFT
          !     box, and that will be cutoff-ed as well
          !
          ! ZL: aqs = <nk|exp(i(q+g).r)|n1k-q>
          !     aqs(nG,num of n1 distributed on each process)
          !     aqs is distributed over bands peinf%ntband_max

          call mtxel(in,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)

          ! ZL: for EP: gvec passed to mtxel provides ng, FFTgrid, components,
          ! which are the same for WFN_inner and dWFN
          ! ZL: in kp and ep_kp, ngk are not the same

!#BEGIN_INTERNAL_ONLY
          if (sig%elph) then
            ! ZL: Calculate <k+phonq|e^{i(q+G)|d_{phonq,v}k-q> matrix elements
            !               <k+phonq|e^{i(q+G)|k+phonq-q>
            !               <k|e^{i(q+G)|d_{-phonq,v}k+phonq-q>
            ! ZL: CAUTION: kp in mtxel, just provides nspinor, so we use regular kp
            !     Interestingly, kp in mtxel is intent inout, even though it is
            !     used as in only
            call mtxel(in, gvec, ep_dq_wfnkq, wfnk_phonq, ncoul, isrtrq, aqs_ep_dq_onel, ispin, kp)
            call mtxel(in, gvec, wfnkq_phonq, wfnk_phonq, ncoul, isrtrq, aqs_phonq_twol, ispin, kp)
            call mtxel(in, gvec, ep_dmq_wfnkq_phonq, wfnk, ncoul, isrtrq, aqs_ep_dmq_twor, ispin, kp)
          endif
!#END_INTERNAL_ONLY

          call timing%stop(timing%mtxel)
          if (sig%noffdiag.gt.0.and.peinf%flag_diag(in)) then
            do j=1,peinf%ntband_node
              do i=1,ncoul
                ! ZL: aqsaug accumulates all aqs
                !     this is used for off-diag, because even in off-diag
                !     all matrix elements have already been calculated
                aqsaug(i,j,peinf%index_diag(in),ispin)=aqs(i,j)
                ! ZL: for EP with off-diag in band index
!#BEGIN_INTERNAL_ONLY
                if (sig%elph) then
                  aqsaug_ep_dq_onel(i,j,peinf%index_diag(in),ispin)=aqs_ep_dq_onel(i,j)
                  aqsaug_phonq_twol(i,j,peinf%index_diag(in),ispin)=aqs_phonq_twol(i,j)
                  aqsaug_ep_dmq_twor(i,j,peinf%index_diag(in),ispin)=aqs_ep_dmq_twor(i,j)
                endif
!#END_INTERNAL_ONLY
              enddo
            enddo
          endif
          if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then
            ! ZL: EP does not include this situation for now
!#BEGIN_INTERNAL_ONLY
            if (sig%elph) then
              call die('EP does not support (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) for now.')
            endif
!#END_INTERNAL_ONLY
            if (irq_==irq_min) then
              call logit('Calling mtxel_occ')
              call timing%start(timing%mtxel_ch)
! JRD ALL PROCS DO THIS NOW. 
              call mtxel_occ(in,in,gvec,wfnk,ncoulch,isrtrq,aqsch,ispin,kp)
              call timing%stop(timing%mtxel_ch)
              if (nrq.gt.1) aqsaugchd(:,in,ispin)=aqsch(:)
            else
              aqsch(:)=aqsaugchd(:,in,ispin)
            endif
          endif

!---------------------
! Compute diag SX and CH

          if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
            call logit('Calling mtxel_cor for diagonal matrix elements')
            call timing%start(timing%mtxel_cor_tot)
            ! ZL: actually there are two M - aqsn and aqsm, following
            !     < psi_n ( k ) | Sigma_cor ( E ) | psi_m ( k ) >
            !     aqsn is the regular M directly calculated e^{i(q+G).r}
            !     aqsm is calculated in the same way, so it need to be
            !     complex-conjugated when being used
            !     so be careful when first calculating aqsm
            !
            ! ZL: in mtxel_cor(), wfnk and wfnkq only provide band energies
            !   - wfnk%ek is used as the energy argument E in Sigma(E), in
            !     general with +/- dE in GPP, for interpolation
            !   - wfnkq%eqk is used, as explicitly shown in equations, band
            !     enregy at e(k-q)
            ! ZL: when dealing with EP, just need to put the correct energy
            !     (i.e. the corresponding wfn variables) into mtxel_cor
            !     Note: we may need to average e(k+q) and e(k) for off-diag
            !           nature of EP, for E.
            ! ZL: NOTE: TODO: if regular GW is done already, maybe we can use
            !     quasiparticle energies instead of LDA energies
            !
            ! ZL Note:
            ! ind, indinv: for eps q+G gvec sort
            ! isrtrq, isrtrqi: sorting of |q+G|^2, only used to determin cutoff
            ! aqsn, aqsm: two M in such way M_n * M_m^*
            ! epsR, epsA, and *D*: are related to full frequency
            ! aqsch: for sig%freq_dep==0 or sig%exact_ch==1, not considered for EP for now
            ! achtcor: for exact_ch, not considered for EP now
            !
            ! Useful variables with out intent: 
            ! asigt_imag, acht_n1, asxt, acht
            ! We also allocate for EP of the following as place holder:
            ! achtcor, achtcor_n1
            !
            ! aqsch: not even place holder for EP. Not even allocated for GPP
            !
            call mtxel_cor(peinf%index_diag(in), &
             peinf%index_diag(in),ispin,ncouls,neps,gvec,eps,ph, &
             ind,indinv,isrtrqi,isrtrq,vcoul,crys,sig,wpg,wfnk,wfnkq,ncoulch, &
             aqs,aqs,aqsch,asigt_imag,acht_n1,asxt,acht,achtcor,achtcor_n1, &
             kp%nspin,rq(:,irq),coulfact, &
             epsmpi%inv_igp_index,epsmpi%ngpown, &
             epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,1, &
             .false.)  ! ZL: set ep_mtxel to be .false. because this is standard
                       !     diagonal quasiparticle energy calculation

!#BEGIN_INTERNAL_ONLY
            if (sig%elph) then
              ! ZL: EP, icalc in mtxel_cor MUST BE 2
              if(sig%ep_bracket.eq.1) then ! ZL: evaluate dSig(E) at |ket>=wfnk band energy
                ! ZL: term one
                call mtxel_cor(peinf%index_diag(in), &
                  peinf%index_diag(in),ispin,ncouls,neps,gvec,eps,ph, &
                  ind,indinv,isrtrqi,isrtrq,vcoul,crys,sig,wpg,wfnk,wfnkq,ncoulch, &
                  aqs_ep_dq_onel,aqs,aqsch,asigt_imag_ep_one,acht_n1_ep_one,asxt_ep_one, &
                  acht_ep_one,achtcor_ep_one,achtcor_n1_ep_one, &
                  kp%nspin,rq(:,irq),coulfact, &
                  epsmpi%inv_igp_index,epsmpi%ngpown, &
                  epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,2, &
                  .false.)  ! ZL: set ep_mtxel to be .false. here because this
                            !     is just diagonal, for corss-check purpose

                ! ZL: term two
                !     term two uses wfnkq_phonq for inner band energy
                call mtxel_cor(peinf%index_diag(in), &
                  peinf%index_diag(in),ispin,ncouls,neps,gvec,eps,ph, &
                  ind,indinv,isrtrqi,isrtrq,vcoul,crys,sig,wpg,wfnk,wfnkq_phonq,ncoulch, &
                  aqs_phonq_twol,aqs_ep_dmq_twor,aqsch,asigt_imag_ep_two,acht_n1_ep_two,asxt_ep_two, &
                  acht_ep_two,achtcor_ep_two,achtcor_n1_ep_two, &
                  kp%nspin,rq(:,irq),coulfact, &
                  epsmpi%inv_igp_index,epsmpi%ngpown, &
                  epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,2, &
                  .false.)  ! ZL: set ep_mtxel to be .false. here because this
                            !     is just diagonal, for corss-check purpose

              elseif(sig%ep_bracket.eq.-1) then ! ZL: evaluate dSig(E) at <bra|=wfnk_phonq band energy
                ! ZL: term one
                call mtxel_cor(peinf%index_diag(in), &
                  peinf%index_diag(in),ispin,ncouls,neps,gvec,eps,ph, &
                  ind,indinv,isrtrqi,isrtrq,vcoul,crys,sig,wpg,wfnk_phonq,wfnkq,ncoulch, &
                  aqs_ep_dq_onel,aqs,aqsch,asigt_imag_ep_one,acht_n1_ep_one,asxt_ep_one, &
                  acht_ep_one,achtcor_ep_one,achtcor_n1_ep_one, &
                  kp%nspin,rq(:,irq),coulfact, &
                  epsmpi%inv_igp_index,epsmpi%ngpown, &
                  epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,2, &
                  .false.)  ! ZL: set ep_mtxel to be .false. here because this
                            !     is just diagonal, for corss-check purpose

                ! ZL: term two
                !     term two uses wfnkq_phonq for inner band energy
                call mtxel_cor(peinf%index_diag(in), &
                  peinf%index_diag(in),ispin,ncouls,neps,gvec,eps,ph, &
                  ind,indinv,isrtrqi,isrtrq,vcoul,crys,sig,wpg,wfnk_phonq,wfnkq_phonq,ncoulch, &
                  aqs_phonq_twol,aqs_ep_dmq_twor,aqsch,asigt_imag_ep_two,acht_n1_ep_two,asxt_ep_two, &
                  acht_ep_two,achtcor_ep_two,achtcor_n1_ep_two, &
                  kp%nspin,rq(:,irq),coulfact, &
                  epsmpi%inv_igp_index,epsmpi%ngpown, &
                  epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,2, &
                  .false.)  ! ZL: set ep_mtxel to be .false. here because this
                            !     is just diagonal, for corss-check purpose
              endif ! sig%ep_bracket
            endif
!#END_INTERNAL_ONLY
            call timing%stop(timing%mtxel_cor_tot)
          else
            achtcor = ZERO
            if (sig%freq_dep/=0 .and. sig%exact_ch==1) achtcor_n1 = ZERO
            asigt_imag = ZERO
!#BEGIN_INTERNAL_ONLY
            if(sig%elph) then
              achtcor_ep_one = ZERO
              if (sig%freq_dep/=0 .and. sig%exact_ch==1) achtcor_n1_ep_one = ZERO
              asigt_imag_ep_one = ZERO

              achtcor_ep_two = ZERO
              if (sig%freq_dep/=0 .and. sig%exact_ch==1) achtcor_n1_ep_two = ZERO
              asigt_imag_ep_two = ZERO
              ! ZL: TODO: check if compatible. For now, we disable
              call die('Electron-phonon does not support this freq_dep!')
            endif
!#END_INTERNAL_ONLY
          endif

          if (peinf%flag_diag(in)) then
            ! FHJ: Store temporary arrays in other variables, and put in
            ! degeneracy factor due to star of q-points.
            jb = peinf%index_diag(in)
            achcor(jb,ispin) = achcor(jb,ispin) + neq(irq)*achtcor
            asig_imag(jb,ispin) = asig_imag(jb,ispin) + neq(irq)*asigt_imag
!#BEGIN_INTERNAL_ONLY
            if(sig%elph) then
              achcor_ep_one(jb,ispin) = achcor_ep_one(jb,ispin) + neq(irq)*achtcor_ep_one
              asig_imag_ep_one(jb,ispin) = asig_imag_ep_one(jb,ispin) + neq(irq)*asigt_imag_ep_one
              achcor_ep_two(jb,ispin) = achcor_ep_two(jb,ispin) + neq(irq)*achtcor_ep_two
              asig_imag_ep_two(jb,ispin) = asig_imag_ep_two(jb,ispin) + neq(irq)*asigt_imag_ep_two
            endif
!#END_INTERNAL_ONLY
            if (sig%freq_dep==2) then
              asxDyn(:,jb,ispin) = asxDyn(:,jb,ispin) + neq(irq)*asxtDyn(:)
              achDyn(:,jb,ispin) = achDyn(:,jb,ispin) + neq(irq)*achtDyn(:)
              achDyn_cor(:,jb,ispin) = achDyn_cor(:,jb,ispin) + neq(irq)*achtDyn_cor(:)
              achDyn_corb(:,jb,ispin) = achDyn_corb(:,jb,ispin) + neq(irq)*achtDyn_corb(:)
              ach2Dyn(:,jb,ispin) = ach2Dyn(:,jb,ispin) + neq(irq)*ach2tDyn(:)
              achD_n1(:,jb,ispin) = achD_n1(:,jb,ispin) + neq(irq)*achtD_n1(:)
            elseif (sig%freq_dep/=-1) then
              asx(:,jb,ispin) = asx(:,jb,ispin) + neq(irq)*asxt(:)
              ach(:,jb,ispin) = ach(:,jb,ispin) + neq(irq)*acht(:)
              ach_n1(:,jb,ispin) = ach_n1(:,jb,ispin) + neq(irq)*acht_n1(:)
!#BEGIN_INTERNAL_ONLY
              if (sig%elph) then
                asx_ep_one(:,jb,ispin) = asx_ep_one(:,jb,ispin) + neq(irq)*asxt_ep_one(:)
                ach_ep_one(:,jb,ispin) = ach_ep_one(:,jb,ispin) + neq(irq)*acht_ep_one(:)
                ach_n1_ep_one(:,jb,ispin) = ach_n1_ep_one(:,jb,ispin) + neq(irq)*acht_n1_ep_one(:)

                asx_ep_two(:,jb,ispin) = asx_ep_two(:,jb,ispin) + neq(irq)*asxt_ep_two(:)
                ach_ep_two(:,jb,ispin) = ach_ep_two(:,jb,ispin) + neq(irq)*acht_ep_two(:)
                ach_n1_ep_two(:,jb,ispin) = ach_n1_ep_two(:,jb,ispin) + neq(irq)*acht_n1_ep_two(:)
              endif
!#END_INTERNAL_ONLY
            endif
            if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
              achcor_n1(:,jb,ispin) = achcor_n1(:,jb,ispin) + neq(irq)*achtcor_n1(:)
!#BEGIN_INTERNAL_ONLY
              if (sig%elph) then
                achcor_n1_ep_one(:,jb,ispin) = achcor_n1_ep_one(:,jb,ispin) + neq(irq)*achtcor_n1_ep_one(:)
                achcor_n1_ep_two(:,jb,ispin) = achcor_n1_ep_two(:,jb,ispin) + neq(irq)*achtcor_n1_ep_two(:)
              endif
!#END_INTERNAL_ONLY
            endif
          endif
          
!---------------------
! Compute diag bare exchange (SKIP this computation if you already know it)

          call timing%start(timing%bare_x)
          if ((.not. sig%use_xdat .and. xflag .and. (.not. sig%coul_mod_flag)) .or. (sig%coul_mod_flag &
            .and. (.not.  sig%use_vxc2dat))) then
            call logit('Computing bare X')
            axt=0.0d0
            ! ZL: EP, bare perturbed exchange part
!#BEGIN_INTERNAL_ONLY
            if(sig%elph) then
              axt_ep_one = 0.0d0
              axt_ep_two = 0.0d0
            endif
!#END_INTERNAL_ONLY
! XXX THREAD?  
            do n1=1,peinf%nvband_node
              tempval = wfnkq%ekq(peinf%indext(n1),ispin) - sig%efermi
              if (.not. sig%is_EXX) then
                if (tempval < sig%tol) then ! ZL: sig%tol is positive, any negative values are smaller than sig%tol
                  if(abs(tempval) < sig%tol) then
                    occ=0.5  ! Fermi-Dirac distribution = 1/2 at Fermi level
                  else
                    occ=1D0
                  endif
                  do ig=1,ncoulb
                    ! ZL: here it calculates bare X term: GV
                    axt = axt + abs(aqs(ig,n1))**2 * occ * vcoul(ig)
                    ! ZL: for EP, involve more terms
!#BEGIN_INTERNAL_ONLY
                    if(sig%elph) then ! EP term one use the same occ(k-p), term two is dealt with in the next
                      axt_ep_one = axt_ep_one + aqs_ep_dq_onel(ig,n1)*MYCONJG(aqs(ig,n1)) * occ * vcoul(ig)
                    endif
!#END_INTERNAL_ONLY
                  enddo
                endif
                ! sig%ncrit = 0 and tempval > sig%tol should never happen!
              else ! sig%is_EXX
                ib = peinf%indext(n1)
                occ = wfnkq%occ(ib+sig%ncore_excl, ispin)
                do ig=1,ncoulb
                  axt = axt + abs(aqs(ig,n1))**2 * occ * vcoul(ig)
                enddo
              endif
!#BEGIN_INTERNAL_ONLY
              if(sig%elph) then ! ZL: for EP term two, use occ(k-p+phonq)
                tempval = wfnkq_phonq%ekq(peinf%indext(n1),ispin) - sig%efermi
                if (tempval < sig%tol) then
                  if(abs(tempval) < sig%tol) then
                    occ=0.5  ! Fermi-Dirac distribution = 1/2 at Fermi level
                  else
                    occ=1D0
                  endif
                  do ig=1,ncoulb ! term two
                    axt_ep_two = axt_ep_two + aqs_phonq_twol(ig,n1)*MYCONJG(aqs_ep_dmq_twor(ig,n1)) * occ * vcoul(ig)
                  enddo
                endif
              endif ! sig%elph
!#END_INTERNAL_ONLY
            enddo
            ! if(peinf%inode.eq.0) write(6,*) "axt, axt_1, axt_2", axt, axt_ep_one, axt_ep_two

            if (peinf%flag_diag(in)) then
              ib = peinf%index_diag(in)
              ax(ib,ispin) = ax(ib,ispin) - neq(irq)*axt*sig%xfrac
!#BEGIN_INTERNAL_ONLY
              if(sig%elph) then
                ax_ep_one(ib,ispin) = ax_ep_one(ib,ispin) - neq(irq)*axt_ep_one*sig%xfrac
                ax_ep_two(ib,ispin) = ax_ep_two(ib,ispin) - neq(irq)*axt_ep_two*sig%xfrac
              endif
!#END_INTERNAL_ONLY
            endif
          endif ! not use x.dat
          call timing%stop(timing%bare_x)
          
        enddo ! in (loop over bands for which we need diag Sigma)
        
        if (ispin.eq.sig%nspin) then
          if (.not.is_subq.or.irq_==sig%nq0) then
            SAFE_DEALLOCATE_P(wfnkq%zkq)
            ! ZL: EP
!#BEGIN_INTERNAL_ONLY
            if(sig%elph) then
              SAFE_DEALLOCATE_P(ep_dq_wfnkq%zkq)
              SAFE_DEALLOCATE_P(wfnkq_phonq%zkq)
              SAFE_DEALLOCATE_P(ep_dmq_wfnkq_phonq%zkq)
            endif
!#END_INTERNAL_ONLY
          endif
          SAFE_DEALLOCATE(aqs)
!#BEGIN_INTERNAL_ONLY
          ! ZL: EP
          if (sig%elph) then
            SAFE_DEALLOCATE(aqs_ep_dq_onel)
            SAFE_DEALLOCATE(aqs_phonq_twol)
            SAFE_DEALLOCATE(aqs_ep_dmq_twor)
          endif
!#END_INTERNAL_ONLY
          if ((sig%freq_dep.eq.0.or.sig%exact_ch.eq.1).and.irq_==nrq) then
            if (nrq.gt.1) then
              SAFE_DEALLOCATE_P(aqsaugchd)
            end if
          endif
        endif
        
!!-------- End diag band loop --------------------------------------------------

! (gsm) begin distributing aqsaug matrix elements for offdiag calculation

! $$$ inefficient communication, this should be rewritten $$$

        call timing%start(timing%mtxel_comm)
#ifdef MPI
        if (sig%noffdiag.gt.0.and.peinf%npools.gt.1) then
          tag=1024
          do in=1,sig%ndiag
            jj=mod(in-1,peinf%npools)+1
            do i=1,peinf%npes/peinf%npools
              source=(i-1)+(jj-1)*(peinf%npes/peinf%npools)
              do j=1,peinf%npools
                dest=(i-1)+(j-1)*(peinf%npes/peinf%npools)
                if (j.ne.jj) then
                  if (peinf%inode.eq.source) &
                    call MPI_Send(aqsaug(:,:,in,ispin),peinf%ntband_max*ncoul,MPI_SCALAR,dest,tag,MPI_COMM_WORLD,mpierr)
                  if (peinf%inode.eq.dest) &
                    call MPI_Recv(aqsaug(:,:,in,ispin),peinf%ntband_max*ncoul,MPI_SCALAR, &
                      source,tag,MPI_COMM_WORLD,mpistatus,mpierr)

!#BEGIN_INTERNAL_ONLY
                  if (sig%elph) then
                    if (peinf%inode.eq.source) &
                      call MPI_Send(aqsaug_ep_dq_onel(:,:,in,ispin),peinf%ntband_max*ncoul,MPI_SCALAR, &
                        dest,tag,MPI_COMM_WORLD,mpierr)
                    if (peinf%inode.eq.dest) &
                      call MPI_Recv(aqsaug_ep_dq_onel(:,:,in,ispin),peinf%ntband_max*ncoul,MPI_SCALAR, &
                        source,tag,MPI_COMM_WORLD,mpistatus,mpierr)

                    if (peinf%inode.eq.source) &
                      call MPI_Send(aqsaug_phonq_twol(:,:,in,ispin),peinf%ntband_max*ncoul,MPI_SCALAR, &
                        dest,tag,MPI_COMM_WORLD,mpierr)
                    if (peinf%inode.eq.dest) &
                      call MPI_Recv(aqsaug_phonq_twol(:,:,in,ispin),peinf%ntband_max*ncoul,MPI_SCALAR, &
                        source,tag,MPI_COMM_WORLD,mpistatus,mpierr)

                    if (peinf%inode.eq.source) &
                      call MPI_Send(aqsaug_ep_dmq_twor(:,:,in,ispin),peinf%ntband_max*ncoul,MPI_SCALAR, &
                        dest,tag,MPI_COMM_WORLD,mpierr)
                    if (peinf%inode.eq.dest) &
                      call MPI_Recv(aqsaug_ep_dmq_twor(:,:,in,ispin),peinf%ntband_max*ncoul,MPI_SCALAR, &
                        source,tag,MPI_COMM_WORLD,mpistatus,mpierr)
                  endif
!#END_INTERNAL_ONLY
                endif
              enddo
            enddo
          enddo
        endif ! sig%noffdiag.gt.0.and.peinf%npools.gt.1
#endif
        call timing%stop(timing%mtxel_comm)
        
! (gsm) end distributing aqsaug matrix elements for offdiag calculation

!!-------- Loop over bands for which offdiag Sigma is computed -----------------

! Bands are relabelled according to sig%off*(1:sig%noffdiag)
! Loop is parallelized according to peinf%index_offdiag(1:peinf%noffdiag_max)

        do ioff=1,peinf%noffdiag_max
          
          call progress_step(prog_info)
          if (peinf%verb_debug .and. peinf%inode==0) then
            if (peinf%npools.eq.1) then
              write(6,998) peinf%index_offdiag(ioff)+(ispin-1)*sig%noffdiag, &
                peinf%noffdiag_max*peinf%npools*sig%nspin
            else
              write(6,996) peinf%index_offdiag(ioff)+(ispin-1)*sig%noffdiag, &
                peinf%index_offdiag(ioff)+(ispin-1)*sig%noffdiag+peinf%npools-1, &
                peinf%noffdiag_max*peinf%npools*sig%nspin
            endif
          endif
998       format(1x,"Computing Sigma offdiag",i4,1x,"of",i4)
996       format(1x,"Computing Sigma offdiag",i4,1x,"to",i4,1x,"of",i4)
          write(tmpstr,'(a,2i6,a)') 'Working on bands ', sig%off1(peinf%index_offdiag(ioff)), &
            sig%off2(peinf%index_offdiag(ioff)), ' 1st pool'
          call logit(tmpstr)

! <n|Sigma|m> = 0 if n and m belong to different irreducible representations
! Even without assigning representations, we can tell they are different if
! the size of the degenerate subspace is different for n and m.
! This saves time and helps enforce symmetry. -- DAS
          bExactlyZero = .false.
          if(.not. sig%wfn_outer_present .and. sig%offdiagsym .and. &
             kp%degeneracy(sig%off1(peinf%index_offdiag(ioff))-sig%ncore_excl, ikn, ispin) /= &
             kp%degeneracy(sig%off2(peinf%index_offdiag(ioff))-sig%ncore_excl, ikn, ispin)) then
            ! ZL: comment out the following line due to messy output
            !if (peinf%inode.eq.0) write(6,'(a)') 'Zero by symmetry -- not computing.'
            ! the matrix elements are zeroed at the beginning, and at the end of the loop
            ! so we can just leave them as they are
            ! JRD - We cannot cycle here because other pools may have legitimate work to do
            ! and we need to be around for the communication.  Particularly inside the subroutine
            ! mtxel_cor. Thus we can`t really save time here. We could still zero out elements 
            ! at bottom of loop if want. I leave it up to DAS.
            !cycle
            bExactlyZero=.true.
          endif
          ! ZL: for EP, we always calculate everything because phonq may break symmetry
          if(sig%elph) bExactlyZero = .false.
            
!---------------------
! Compute planewave matrix elements for exact static CH
! ZL: wfnkoff is only used for this purpose

          if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then
            if (irq_==irq_min) then
              call timing%start(timing%wf_ch_comm)
              wfnkoff%nkpt=wfnk%nkpt
              if (mod(peinf%inode,peinf%npes/peinf%npools).eq.0) then
                SAFE_ALLOCATE(wfnkoff%isrtk, (gvec%ng))
                wfnkoff%isrtk=wfnk%isrtk
                SAFE_ALLOCATE(wfnkoff%zk, (2*wfnk%nkpt,kp%nspinor))
              endif
                
! (gsm) begin gathering wavefunctions over pools

! $$$ inefficient communication, this should be rewritten $$$

! BAB note: loop over spinors

              do jj=1,peinf%npools
                dest=(jj-1)*(peinf%npes/peinf%npools)
                do ii=1,2
                  i=sig%offmap(peinf%index_offdiag(ioff),ii)
#ifdef MPI
                  call MPI_Bcast(i,1,MPI_INTEGER,dest,MPI_COMM_WORLD,mpierr)
#endif
                  j=(i-1)/peinf%npools+1
                  source=mod(i-1,peinf%npools)*(peinf%npes/peinf%npools)
                  if (peinf%inode.eq.source.and.peinf%inode.eq.dest) then
                    do jsp=1,kp%nspinor
                      do k=1,wfnk%nkpt
                        wfnkoff%zk((ii-1)*wfnk%nkpt+k,jsp) = wfnk%zk((j-1)*wfnk%nkpt+k,ispin*jsp)
                      enddo
                    enddo
                  else
#ifdef MPI
                    tag=1024
                    do jsp=1,kp%nspinor
                      if (peinf%inode.eq.source) &
                        call MPI_Send(wfnk%zk((j-1)*wfnk%nkpt+1,ispin*jsp),wfnk%nkpt, &
                          MPI_SCALAR,dest,tag,MPI_COMM_WORLD,mpierr)
                      if (peinf%inode.eq.dest) &
                        call MPI_Recv(wfnkoff%zk((ii-1)*wfnk%nkpt+1,jsp),wfnk%nkpt, &
                          MPI_SCALAR,source,tag,MPI_COMM_WORLD,mpistatus,mpierr)
                    enddo
#else
                    do jsp=1,kp%nspinor
                      do k=1,wfnk%nkpt
                        wfnkoff%zk((ii-1)*wfnk%nkpt+k,jsp) = wfnk%zk((j-1)*wfnk%nkpt+k,ispin*jsp)
                      enddo
                    enddo
#endif
                  endif
                enddo
              enddo
              
! (gsm) end gathering wavefunctions over pools

              call timing%stop(timing%wf_ch_comm)
              call logit('Calling mtxel_occ')
              call timing%start(timing%mtxel_ch)
! JRD Everyone does this now YYYY
              if (mod(peinf%inode,peinf%npes/peinf%npools).eq.0) then
                call mtxel_occ(1,2,gvec,wfnkoff,ncoulch,isrtrq,aqsch,1,kp)
              endif
#ifdef MPI
              call MPI_Bcast(aqsch,ncoulch,MPI_SCALAR,0,peinf%pool_comm,mpierr)
#endif
              call timing%stop(timing%mtxel_ch)
              if (mod(peinf%inode,peinf%npes/peinf%npools).eq.0) then
                SAFE_DEALLOCATE_P(wfnkoff%isrtk)
                SAFE_DEALLOCATE_P(wfnkoff%zk)
              endif
              if (nrq.gt.1) aqsaugcho(:,ioff,ispin)=aqsch(:)
            else
              aqsch(:)=aqsaugcho(:,ioff,ispin)
            endif
          endif ! sig%freq_dep.eq.0.or.sig%exact_ch.eq.1

!---------------------
! Compute offdiag SX and CH
! ZL: EP related start from here

          if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
            call logit('Calling mtxel_cor for offdiagonal matrix elements')
            call timing%start(timing%mtxel_cor_tot)
            call mtxel_cor(sig%offmap(peinf%index_offdiag(ioff),1), &
              sig%offmap(peinf%index_offdiag(ioff),3), &
              ispin,ncouls,neps,gvec,eps,ph,ind,indinv,isrtrqi, &
              isrtrq,vcoul,crys,sig,wpg,wfnk,wfnkq,ncoulch, &
              aqsaug(:,:,sig%offmap(peinf%index_offdiag(ioff),1),ispin), &
              aqsaug(:,:,sig%offmap(peinf%index_offdiag(ioff),2),ispin), &
              aqsch,asigt_imag,acht_n1,asxt,acht,achtcor,achtcor_n1,kp%nspin, &
              rq(:,irq),coulfact, &
              epsmpi%inv_igp_index,epsmpi%ngpown, &
              epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,2, &
              .false.)  ! ZL: set ep_mtxel to be .false. because this is just
                        !     standard off-diagonal quasiparticle Sigma
                        !     calculation
            
             ! ZL: The following piece calculates off-diag (in band index) for
             !     perturbed Sigma operator matrix element
!#BEGIN_INTERNAL_ONLY
            if (sig%elph) then
              if(sig%ep_bracket.eq.1) then ! ZL: evaluate dSig(E) at |ket>=wfnk band energy
                ! ZL: term one off-diag in bands
                call mtxel_cor(sig%offmap(peinf%index_offdiag(ioff),1), &
                  sig%offmap(peinf%index_offdiag(ioff),3), &
                  ispin,ncouls,neps,gvec,eps,ph,ind,indinv,isrtrqi, &
                  isrtrq,vcoul,crys,sig,wpg,wfnk,wfnkq,ncoulch, &
                  aqsaug_ep_dq_onel(:,:,sig%offmap(peinf%index_offdiag(ioff),1),ispin), &
                  aqsaug(:,:,sig%offmap(peinf%index_offdiag(ioff),2),ispin), &
                  aqsch,asigt_imag_ep_one,acht_n1_ep_one,asxt_ep_one,acht_ep_one,achtcor_ep_one,achtcor_n1_ep_one,kp%nspin, &
                  rq(:,irq),coulfact, &
                  epsmpi%inv_igp_index,epsmpi%ngpown, &
                  epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,2, &
                  .true., wfnk_phonq, sig%offmap_ep(peinf%index_offdiag(ioff)))  ! ZL: set ep_mtxel to be .true.
                ! ZL: term two off-diag in bands
                !     term two uses wfnkq_phonq for inner band energy
                call mtxel_cor(sig%offmap(peinf%index_offdiag(ioff),1), &
                  sig%offmap(peinf%index_offdiag(ioff),3), &
                  ispin,ncouls,neps,gvec,eps,ph,ind,indinv,isrtrqi, &
                  isrtrq,vcoul,crys,sig,wpg,wfnk,wfnkq_phonq,ncoulch, &
                  aqsaug_phonq_twol(:,:,sig%offmap(peinf%index_offdiag(ioff),1),ispin), &
                  aqsaug_ep_dmq_twor(:,:,sig%offmap(peinf%index_offdiag(ioff),2),ispin), &
                  aqsch,asigt_imag_ep_two,acht_n1_ep_two,asxt_ep_two,acht_ep_two,achtcor_ep_two,achtcor_n1_ep_two,kp%nspin, &
                  rq(:,irq),coulfact, &
                  epsmpi%inv_igp_index,epsmpi%ngpown, &
                  epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,2, &
                  .true., wfnk_phonq, sig%offmap_ep(peinf%index_offdiag(ioff)))  ! ZL: set ep_mtxel to be .true.

              elseif(sig%ep_bracket.eq.-1) then ! ZL: evaluate dSig(E) at <bra|=wfnk_phonq band energy
                ! ZL: term one off-diag in bands
                call mtxel_cor(sig%offmap(peinf%index_offdiag(ioff),1), &
                  sig%offmap(peinf%index_offdiag(ioff),3), &
                  ispin,ncouls,neps,gvec,eps,ph,ind,indinv,isrtrqi, &
                  isrtrq,vcoul,crys,sig,wpg,wfnk_phonq,wfnkq,ncoulch, &
                  aqsaug_ep_dq_onel(:,:,sig%offmap(peinf%index_offdiag(ioff),1),ispin), &
                  aqsaug(:,:,sig%offmap(peinf%index_offdiag(ioff),2),ispin), &
                  aqsch,asigt_imag_ep_one,acht_n1_ep_one,asxt_ep_one,acht_ep_one,achtcor_ep_one,achtcor_n1_ep_one,kp%nspin, &
                  rq(:,irq),coulfact, &
                  epsmpi%inv_igp_index,epsmpi%ngpown, &
                  epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,2, &
                  .true., wfnk, sig%offmap_ep(peinf%index_offdiag(ioff)))  ! ZL: set ep_mtxel to be .true.
                ! ZL: term two off-diag in bands
                !     term two uses wfnkq_phonq for inner band energy
                call mtxel_cor(sig%offmap(peinf%index_offdiag(ioff),1), &
                  sig%offmap(peinf%index_offdiag(ioff),3), &
                  ispin,ncouls,neps,gvec,eps,ph,ind,indinv,isrtrqi, &
                  isrtrq,vcoul,crys,sig,wpg,wfnk_phonq,wfnkq_phonq,ncoulch, &
                  aqsaug_phonq_twol(:,:,sig%offmap(peinf%index_offdiag(ioff),1),ispin), &
                  aqsaug_ep_dmq_twor(:,:,sig%offmap(peinf%index_offdiag(ioff),2),ispin), &
                  aqsch,asigt_imag_ep_two,acht_n1_ep_two,asxt_ep_two,acht_ep_two,achtcor_ep_two,achtcor_n1_ep_two,kp%nspin, &
                  rq(:,irq),coulfact, &
                  epsmpi%inv_igp_index,epsmpi%ngpown, &
                  epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,2, &
                  .true., wfnk, sig%offmap_ep(peinf%index_offdiag(ioff)))  ! ZL: set ep_mtxel to be .true.
              endif ! sig%ep_bracket
            endif
!#END_INTERNAL_ONLY
            call timing%stop(timing%mtxel_cor_tot)
          endif
          
          if (bExactlyZero) cycle
         
          if (peinf%flag_offdiag(ioff)) then
            ! FHJ: Store temporary arrays in other variables, and put in
            ! degeneracy factor due to star of q-points.
            jb = peinf%index_offdiag(ioff) + sig%ndiag
            achcor(jb,ispin) = achcor(jb,ispin) + ryd*neq(irq)*achtcor
!#BEGIN_INTERNAL_ONLY
            if (sig%elph) then
              achcor_ep_one(jb,ispin) = achcor_ep_one(jb,ispin) + ryd*neq(irq)*achtcor_ep_one
              achcor_ep_two(jb,ispin) = achcor_ep_two(jb,ispin) + ryd*neq(irq)*achtcor_ep_two
            endif
!#END_INTERNAL_ONLY
            if (sig%freq_dep==2) then
              asxDyn(:,jb,ispin) = asxDyn(:,jb,ispin) + ryd*neq(irq)*asxtDyn(:)
              achDyn(:,jb,ispin) = achDyn(:,jb,ispin) + ryd*neq(irq)*achtDyn(:)
              achDyn_cor(:,jb,ispin) = achDyn_cor(:,jb,ispin) + ryd*neq(irq)*achtDyn_cor(:)
              achDyn_corb(:,jb,ispin) = achDyn_corb(:,jb,ispin) + ryd*neq(irq)*achtDyn_corb(:)
              ach2Dyn(:,jb,ispin) = ach2Dyn(:,jb,ispin) + ryd*neq(irq)*ach2tDyn(:)
              achD_n1(:,jb,ispin) = achD_n1(:,jb,ispin) + neq(irq)*achtD_n1(:)
            elseif (sig%freq_dep/=-1) then
              asig_imag(jb,ispin) = asig_imag(jb,ispin) + ryd*neq(irq)*asigt_imag
              asx(:,jb,ispin) = asx(:,jb,ispin) + ryd*neq(irq)*asxt(:)
              ach(:,jb,ispin) = ach(:,jb,ispin) + ryd*neq(irq)*acht(:)
              ach_n1(:,jb,ispin) = ach_n1(:,jb,ispin) + neq(irq)*acht_n1(:)
!#BEGIN_INTERNAL_ONLY
              if (sig%elph) then
                asig_imag_ep_one(jb,ispin) = asig_imag_ep_one(jb,ispin) + ryd*neq(irq)*asigt_imag_ep_one
                asx_ep_one(:,jb,ispin) = asx_ep_one(:,jb,ispin) + ryd*neq(irq)*asxt_ep_one(:)
                ach_ep_one(:,jb,ispin) = ach_ep_one(:,jb,ispin) + ryd*neq(irq)*acht_ep_one(:)
                ach_n1_ep_one(:,jb,ispin) = ach_n1_ep_one(:,jb,ispin) + neq(irq)*acht_n1_ep_one(:)

                asig_imag_ep_two(jb,ispin) = asig_imag_ep_two(jb,ispin) + ryd*neq(irq)*asigt_imag_ep_two
                asx_ep_two(:,jb,ispin) = asx_ep_two(:,jb,ispin) + ryd*neq(irq)*asxt_ep_two(:)
                ach_ep_two(:,jb,ispin) = ach_ep_two(:,jb,ispin) + ryd*neq(irq)*acht_ep_two(:)
                ach_n1_ep_two(:,jb,ispin) = ach_n1_ep_two(:,jb,ispin) + neq(irq)*acht_n1_ep_two(:)
              endif
!#END_INTERNAL_ONLY
            endif
            if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
              achcor_n1(:,jb,ispin) = achcor_n1(:,jb,ispin) + neq(irq)*achtcor_n1(:)
!#BEGIN_INTERNAL_ONLY
              if (sig%elph) then
                achcor_n1_ep_one(:,jb,ispin) = achcor_n1_ep_one(:,jb,ispin) + neq(irq)*achtcor_n1_ep_one(:)
                achcor_n1_ep_two(:,jb,ispin) = achcor_n1_ep_two(:,jb,ispin) + neq(irq)*achtcor_n1_ep_two(:)
              endif
!#END_INTERNAL_ONLY
            endif
          endif
          
!---------------------
! Compute offdiag bare exchange (SKIP this computation if you already know it)

          call timing%start(timing%bare_x)
          if ((.not. sig%use_xdat .and. xflag .and. (.not. sig%coul_mod_flag)) .or. (sig%coul_mod_flag &
            .and. (.not.  sig%use_vxc2dat))) then
            axt=0.0d0
!#BEGIN_INTERNAL_ONLY
            if(sig%elph) then
              axt_ep_one = 0.0d0
              axt_ep_two = 0.0d0
            endif
!#END_INTERNAL_ONLY
! XXX THREAD?
            do n1=1,peinf%nvband_node
              tempval = wfnkq%ekq(peinf%indext(n1),ispin) - sig%efermi
              if (tempval < sig%tol) then
                if(abs(tempval) < sig%tol) then
                  occ=0.5  ! Fermi-Dirac distribution = 1/2 at Fermi level
                else
                  occ=1D0
                endif
                do ig=1,ncoulb
                  axt=axt+aqsaug(ig,n1,sig%offmap(peinf%index_offdiag(ioff),1),ispin) &
                    *MYCONJG(aqsaug(ig,n1,sig%offmap(peinf%index_offdiag(ioff),2),ispin))*occ &
                    *vcoul(ig)
!#BEGIN_INTERNAL_ONLY
                  ! ZL: for EP term one
                  if(sig%elph) then
                    axt_ep_one=axt_ep_one &
                      +aqsaug_ep_dq_onel(ig,n1,sig%offmap(peinf%index_offdiag(ioff),1),ispin) &
                      *MYCONJG(aqsaug(ig,n1,sig%offmap(peinf%index_offdiag(ioff),2),ispin))*occ &
                      *vcoul(ig)
                  endif
!#END_INTERNAL_ONLY
                enddo
              endif
!#BEGIN_INTERNAL_ONLY
              if(sig%elph) then ! ZL: for EP term two
                tempval = wfnkq_phonq%ekq(peinf%indext(n1),ispin) - sig%efermi
                if (tempval < sig%tol) then
                  if(abs(tempval) < sig%tol) then
                    occ=0.5  ! Fermi-Dirac distribution = 1/2 at Fermi level
                  else
                    occ=1D0
                  endif
                  do ig=1,ncoulb
                    axt_ep_two=axt_ep_two &
                      +aqsaug_phonq_twol(ig,n1,sig%offmap(peinf%index_offdiag(ioff),1),ispin) &
                      *MYCONJG(aqsaug_ep_dmq_twor(ig,n1,sig%offmap(peinf%index_offdiag(ioff),2),ispin))*occ &
                      *vcoul(ig)
                  enddo
                endif
              endif ! sig%elph
!#END_INTERNAL_ONLY
            enddo
            if (peinf%flag_offdiag(ioff)) then
              ib = peinf%index_offdiag(ioff) + sig%ndiag
              ax(ib,ispin) = ax(ib,ispin) - neq(irq)*axt*ryd*sig%xfrac
!#BEGIN_INTERNAL_ONLY
              if(sig%elph) then
                ! ZL: unit is unified for EP as well
                ax_ep_one(ib,ispin) = ax_ep_one(ib,ispin) - neq(irq)*axt_ep_one*ryd*sig%xfrac
                ax_ep_two(ib,ispin) = ax_ep_two(ib,ispin) - neq(irq)*axt_ep_two*ryd*sig%xfrac
              endif
!#END_INTERNAL_ONLY
            endif
          endif ! not using x.dat
          call timing%stop(timing%bare_x)
          
        enddo ! ioff (loop over bands for which we need offdiag Sigma)
        
        if (ispin.eq.sig%nspin) then
          SAFE_DEALLOCATE(vcoul)
          if (sig%noffdiag.gt.0) then
            SAFE_DEALLOCATE(aqsaug)
!#BEGIN_INTERNAL_ONLY
            if(sig%elph) then
              SAFE_DEALLOCATE(aqsaug_ep_dq_onel)
              SAFE_DEALLOCATE(aqsaug_phonq_twol)
              SAFE_DEALLOCATE(aqsaug_ep_dmq_twor)
            endif
!#END_INTERNAL_ONLY
          end if
          if ((sig%freq_dep.eq.0.or.sig%exact_ch.eq.1).and.irq_==nrq) then
            SAFE_DEALLOCATE_P(aqsch)
            if (nrq.gt.1.and.sig%noffdiag.gt.0) then
              SAFE_DEALLOCATE_P(aqsaugcho)
            end if
          endif
        endif
!!-------- End offdiag band loop -----------------------------------------------
      enddo ! ispin
!!-------- End spin loop -------------------------------------------------------
    enddo ! irq (loop over rq point in BZ summation)
    call progress_free(prog_info)
!!-------- End loop over k-points in irr bz with respect to kn (rq) ------------

    call timing%stop(timing%eqp_loop_total)
#ifdef MPI
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

    ! deallocate GPU stuff
    if ( mtxel_algo /= CPU_ALGO ) then
      call deallocate_accel_mtxel_sig(mtxel_algo)
    end if

#ifdef MPI
    call logit('Reducing arrays from all PEs')
    if (.not.sig%use_xdat) call xreduce(ax, size(ax))
    call zreduce(achcor, size(achcor))
!#BEGIN_INTERNAL_ONLY
    if (sig%elph) then
      call xreduce(ax_ep_one, size(ax_ep_one))
      call xreduce(ax_ep_two, size(ax_ep_two))
      call xreduce(ax_ep, size(ax_ep))
      call zreduce(achcor_ep_one, size(achcor_ep_one))
      call zreduce(achcor_ep_two, size(achcor_ep_two))
      call zreduce(achcor_ep, size(achcor_ep))
    endif
!#END_INTERNAL_ONLY
    if (sig%freq_dep==2) then
      call zreduce(asxDyn, size(asxDyn))
      call zreduce(achDyn, size(achDyn))
      call zreduce(achDyn_cor, size(achDyn_cor))
      call zreduce(achDyn_corb, size(achDyn_corb))
      call zreduce(ach2Dyn, size(ach2Dyn))
      call zreduce(achD_n1, size(achD_n1))
    else
      !FHJ: FIXME - Why do we store sx, ch, etc. for HF calculations??
      call xreduce(asx, size(asx))
      call xreduce(ach, size(ach))
      call zreduce(asig_imag, size(asig_imag))
      call xreduce(ach_n1, size(ach_n1))
!#BEGIN_INTERNAL_ONLY
      if (sig%elph) then
        call xreduce(asx_ep_one, size(asx_ep_one))
        call xreduce(ach_ep_one, size(ach_ep_one))
        call zreduce(asig_imag_ep_one, size(asig_imag_ep_one))
        call xreduce(ach_n1_ep_one, size(ach_n1_ep_one))

        call xreduce(asx_ep_two, size(asx_ep_two))
        call xreduce(ach_ep_two, size(ach_ep_two))
        call zreduce(asig_imag_ep_two, size(asig_imag_ep_two))
        call xreduce(ach_n1_ep_two, size(ach_n1_ep_two))

        call xreduce(asx_ep, size(asx_ep))
        call xreduce(ach_ep, size(ach_ep))
        call zreduce(asig_imag_ep, size(asig_imag_ep))
        call xreduce(ach_n1_ep, size(ach_n1_ep))
      endif
!#END_INTERNAL_ONLY
    endif
    if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
      call xreduce(achcor_n1, size(achcor_n1))
!#BEGIN_INTERNAL_ONLY
      if (sig%elph) then
        call xreduce(achcor_n1_ep_one, size(achcor_n1_ep_one))
        call xreduce(achcor_n1_ep_two, size(achcor_n1_ep_two))
        call xreduce(achcor_n1_ep, size(achcor_n1_ep))
      endif
!#END_INTERNAL_ONLY
    endif
    call logit('Done with reductions')
#endif

!#BEGIN_INTERNAL_ONLY
    ! ZL: EP terms
    if (sig%elph) then
      ax_ep = ax_ep_one + ax_ep_two
      achcor_ep = achcor_ep_one + achcor_ep_two
      asx_ep = asx_ep_one + asx_ep_two
      ach_ep = ach_ep_one + ach_ep_two
      asig_imag_ep = asig_imag_ep_one + asig_imag_ep_two
      ach_n1_ep = ach_n1_ep_one + ach_n1_ep_two
      achcor_n1_ep = achcor_n1_ep_one + achcor_n1_ep_two
    endif
!#END_INTERNAL_ONLY

!#BEGIN_INTERNAL_ONLY
    ! If this is modified hybrid functional calculation
    if (sig%coul_mod_flag .and. (.not.sig%use_vxc2dat)) then
      ! ZL: EP does not consider this situation
      if (sig%elph) then
        call die("Electron-phonon does not consider Hybrid case.")
      endif
      ax(:,:) = ax(:,:) + alda2(:,:)
    endif
!#END_INTERNAL_ONLY
    
!----------------------------
! Output unsymmetrized values of X,SX,CH
! Symmetrize X,SX,CH matrix elements over degenerate states
! Convert to eV and output symmetrized values of X,SX,CH

    if (peinf%inode.eq.0) then

      ! This beautiful piece of code should never be executed and here only because it somehow
      ! prevents gfortran from an incorrect optimization of this routine that will produce NaNs.
      ! It has some mysterious relation to the igp loop above reading from iunit_eps.
      ! Insanely, the presence of function 'isNaN' appears to be crucial here. --DAS
#ifdef GNU
      if(ikn == -1) then
        call die("BUG: ikn = -1")
        write(0,*) isNaN(fact)
      endif
#endif

! DVF : sig%ncore_excl has to be substracted here because wfnk%ek/elda are defined in the 
! read_wavefunction subroutine in input.f90 to be referenced to the case with 
! no core states. Same goes for all subsequent occurences below. 

      do ispin=1,sig%nspin
        write(6,989)ikn,sig%spin_index(ispin)
        if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
          write(6,*)
          ! ZL: either print Vxc or KIH
          if(.not. sig%use_kihdat) then ! standard Vxc output
            write(6,900) &
              "n", "Emf", "Eo", "Vxc", "X", "SX-X", "CH", "Cor", "Sig"
            do i=1,sig%ndiag
              write(6,901) sig%diag(i), wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), &
                wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin), dble(alda(i,ispin))*ryd, &
                dble(ax(i,ispin))*ryd, dble(asx(2,i,ispin))*ryd, &
                dble(ach(2,i,ispin)+achcor(i,ispin))*ryd, &
                dble(asx(2,i,ispin)+ach(2,i,ispin)+achcor(i,ispin))*ryd, &
                dble(ax(i,ispin)+asx(2,i,ispin)+ach(2,i,ispin)+achcor(i,ispin))*ryd
            enddo
          else ! KIH output
            write(6,900) &
              "n", "Emf", "Eo", "KIH", "X", "SX-X", "CH", "Cor", "Sig"
            do i=1,sig%ndiag
              write(6,901) sig%diag(i), wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), &
                wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin), dble(akih(i,ispin))*ryd, &
                dble(ax(i,ispin))*ryd, dble(asx(2,i,ispin))*ryd, &
                dble(ach(2,i,ispin)+achcor(i,ispin))*ryd, &
                dble(asx(2,i,ispin)+ach(2,i,ispin)+achcor(i,ispin))*ryd, &
                dble(ax(i,ispin)+asx(2,i,ispin)+ach(2,i,ispin)+achcor(i,ispin))*ryd
            enddo
          endif

!#BEGIN_INTERNAL_ONLY
          if (sig%elph) then
            write(6,*)
            write(6,*) "================ Electron-phonon perturbed matrix elements ================"
            write(6,*)
            write(6,*) "Electron-phonon band index diagonal (real part):"
            if(sig%ep_bracket.eq.1) then
              write(6,*) "dSig(E) evaluated at |ket> = |k> band energy"
            elseif(sig%ep_bracket.eq.-1) then
              write(6,*) "dSig(E) evaluated at <bra| = <k+phonq| band energy"
            endif
            write(6,*)
            write(6,905) &
              "n", "Emf", "Eo", "dVxc", "dX", "dSX-dX", "dCH", "dCor", "dSig", "dSig-dVxc"
            do i=1,sig%ndiag
              if(sig%ep_bracket.eq.1) then
                Eo_save = wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin)
              elseif(sig%ep_bracket.eq.-1) then
                Eo_save = wfnk_phonq%ek(sig%diag(i)-sig%ncore_excl,ispin)
              endif
              write(6,906) sig%diag(i), wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), &
                Eo_save, &
                dble(alda_ep(i,ispin))*ryd, &
                dble(ax_ep(i,ispin))*ryd, &
                dble(asx_ep(2,i,ispin))*ryd, &
                dble(ach_ep(2,i,ispin)+achcor_ep(i,ispin))*ryd, &
                dble(asx_ep(2,i,ispin)+ach_ep(2,i,ispin)+achcor_ep(i,ispin))*ryd, &
                dble(ax_ep(i,ispin)+asx_ep(2,i,ispin)+ach_ep(2,i,ispin)+achcor_ep(i,ispin))*ryd, &
                dble(ax_ep(i,ispin)+asx_ep(2,i,ispin)+ach_ep(2,i,ispin)+achcor_ep(i,ispin)-alda_ep(i,ispin))*ryd
            enddo
            ! ZL: band index offdiag
            if (sig%noffdiag .gt. 0) then
              write(6,*)
              write(6,*) "Electron-phonon band index off-diagonal:"
              if(sig%ep_bracket.eq.1) then
                write(6,*) "dSig(E) evaluated at |ket> = |k> band energy"
              elseif(sig%ep_bracket.eq.-1) then
                write(6,*) "dSig(E) evaluated at <bra| = <k+phonq| band energy"
              endif
              write(6,'(/,3x,"n",3x,"m",3x,"l",2x,"lp",6x,6x,"dVxc",8x,"dX",4x,"dSX-dX",7x,"dCH",6x,"dCor",3x,"dSig(l)",2x,"dSig-dVxc",2x,"dSig1(lp)",2x,"dSig3([l+lp]/2)")')
              do i=sig%ndiag+1,sig%ndiag+sig%noffdiag
                write(6,'(4i4,1x,"real",1x,9f10.4)') sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
                  sig%off3(i-sig%ndiag), &
                  sig%off_ep(i-sig%ndiag), &
                  dble(alda_ep(i,ispin)), &
                  dble(ax_ep(i,ispin)), &
                  dble(asx_ep(2,i,ispin)), &
                  dble(ach_ep(2,i,ispin)+achcor_ep(i,ispin)), &
                  dble(asx_ep(2,i,ispin)+ach_ep(2,i,ispin)+achcor_ep(i,ispin)), &
                  dble(ax_ep(i,ispin)+asx_ep(2,i,ispin)+ach_ep(2,i,ispin)+achcor_ep(i,ispin)), &
                  dble(ax_ep(i,ispin)+asx_ep(2,i,ispin)+ach_ep(2,i,ispin)+achcor_ep(i,ispin)-alda_ep(i,ispin)), &
                  dble(ax_ep(i,ispin)+asx_ep(1,i,ispin)+ach_ep(1,i,ispin)+achcor_ep(i,ispin)), &
                  dble(ax_ep(i,ispin)+asx_ep(3,i,ispin)+ach_ep(3,i,ispin)+achcor_ep(i,ispin))
#ifdef CPLX
                ! ZL: EP only support complex version
                write(6,'(4i4,1x,"imag",1x,9f10.4)') sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
                  sig%off3(i-sig%ndiag), &
                  sig%off_ep(i-sig%ndiag), &
                  IMAG(alda_ep(i,ispin)), &
                  IMAG(ax_ep(i,ispin)), &
                  IMAG(asx_ep(2,i,ispin)), &
                  IMAG(ach_ep(2,i,ispin)+achcor_ep(i,ispin)), &
                  IMAG(asx_ep(2,i,ispin)+ach_ep(2,i,ispin)+achcor_ep(i,ispin)), &
                  IMAG(ax_ep(i,ispin)+asx_ep(2,i,ispin)+ach_ep(2,i,ispin)+achcor_ep(i,ispin)), &
                  IMAG(ax_ep(i,ispin)+asx_ep(2,i,ispin)+ach_ep(2,i,ispin)+achcor_ep(i,ispin)-alda_ep(i,ispin)), &
                  IMAG(ax_ep(i,ispin)+asx_ep(1,i,ispin)+ach_ep(1,i,ispin)+achcor_ep(i,ispin)), &
                  IMAG(ax_ep(i,ispin)+asx_ep(3,i,ispin)+ach_ep(3,i,ispin)+achcor_ep(i,ispin))
#endif
              enddo
            endif
            write(6,*)
            write(6,*) "================================ EP done! ================================="
          endif
!#END_INTERNAL_ONLY
        endif
          
        ! ZL: this is full frequency, EP does not support yet
        if (sig%freq_dep.eq.2) then
          write(6,*)
          if (sig%freq_dep_method==2) then
            ! ZL: Vxc vs. KIH
            if(.not. sig%use_kihdat) then
              write(6,900) &
                "n", "Emf", "Eo", "X", "Vxc", "Re Res", "Re Int", "Re Sig"
              write(6,900) &
                '', '', '', '', '', "Im Res", "Im Int", "Im Sig"
            else ! KIH
              write(6,900) &
                "n", "Emf", "Eo", "X", "KIH", "Re Res", "Re Int", "Re Sig"
              write(6,900) &
                '', '', '', '', '', "Im Res", "Im Int", "Im Sig"
            endif
          else
            ! ZL: Vxc vs. KIH
            if(.not. sig%use_kihdat) then
              write(6,900) &
                "n", "Emf", "Eo", "X", "Vxc", "Re SX-X", "Re CH", "Re Sig"
              write(6,900) &
                '', '', '', '', '', "Im SX-X", "Im CH", "Im Sig"
            else ! KIH
              write(6,900) &
                "n", "Emf", "Eo", "X", "KIH", "Re SX-X", "Re CH", "Re Sig"
              write(6,900) &
                '', '', '', '', '', "Im SX-X", "Im CH", "Im Sig"
            endif
          endif
          do i=1,sig%ndiag
            
! JRD: Find iw closest to e_lk
            diffmin = INF
            e_lk = wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin)
            ! FHJ: Figure out starting frequency for freq. grid
            if (sig%freq_grid_shift<2) then
              freq0 = sig%freqevalmin
            else
              freq0 = e_lk - sig%freqevalstep*(sig%nfreqeval-1)/2
            endif
            do iw=1,sig%nfreqeval
              diff = abs(freq0 + (iw-1)*sig%freqevalstep - e_lk)
              if (diff .lt. diffmin) then
                diffmin=diff
                iwlda=iw
              endif
            enddo
            ! ZL: print Vxc or KIH
            if(.not. sig%use_kihdat) then
              write(6,901) sig%diag(i), wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), &
                wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin), dble(ax(i,ispin))*ryd, &
                dble(alda(i,ispin))*ryd, dble(asxDyn(iwlda,i,ispin))*ryd, &
                dble(achDyn(iwlda,i,ispin)+achcor(i,ispin))*ryd, &
                dble(ax(i,ispin)+asxDyn(iwlda,i,ispin)+achDyn(iwlda,i,ispin)+achcor(i,ispin))*ryd, &
                dble(alda(i,ispin))*ryd
              write(6,902) IMAG(asxDyn(iwlda,i,ispin))*ryd, IMAG(achDyn(iwlda,i,ispin))*ryd, &
                IMAG(ax(i,ispin)+asxDyn(iwlda,i,ispin)+achDyn(iwlda,i,ispin))*ryd, &
                IMAG(ax(i,ispin)+asxDyn(iwlda,i,ispin)+ach2Dyn(iwlda,i,ispin))*ryd
            else ! KIH
              write(6,901) sig%diag(i), wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), &
                wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin), dble(ax(i,ispin))*ryd, &
                dble(akih(i,ispin))*ryd, dble(asxDyn(iwlda,i,ispin))*ryd, &
                dble(achDyn(iwlda,i,ispin)+achcor(i,ispin))*ryd, &
                dble(ax(i,ispin)+asxDyn(iwlda,i,ispin)+achDyn(iwlda,i,ispin)+achcor(i,ispin))*ryd, &
                dble(akih(i,ispin))*ryd
              ! ZL: turns out that following format 901, the last redundant akih will not be printed anyways
              !     but we keep it the same as the original code
              write(6,902) IMAG(asxDyn(iwlda,i,ispin))*ryd, IMAG(achDyn(iwlda,i,ispin))*ryd, &
                IMAG(ax(i,ispin)+asxDyn(iwlda,i,ispin)+achDyn(iwlda,i,ispin))*ryd, &
                IMAG(ax(i,ispin)+asxDyn(iwlda,i,ispin)+ach2Dyn(iwlda,i,ispin))*ryd
            endif
          enddo
        endif
      enddo
    endif ! node 0
989 format(1x,"Unsymmetrized values for ik =",i4,1x,"spin =",i2)
900 format(a6,8(a9))
901 format(i6,8(f9.3))
902 format(6x,4(9x),4(f9.3))
905 format(a6,8(a9),(a11))
906 format(i6,9(f9.3))

!----------------------------
! Symmetrize matrix elements

    if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
      call shiftenergy(sig, wfnk, alda, asx, ach, achcor, ach_n1, achcor_n1, &
        ax, efsto, asig, enew, zrenorm, nfreqgpp, sig%ncore_excl, akih) ! ZL: akih is optional argument
    endif
    if (sig%freq_dep.eq.2) then
      call shiftenergy_dyn(sig, wfnk, alda, asxDyn, achDyn, achDyn_cor, &
        achDyn_corb, ach2Dyn, achcor, achD_n1, achcor_n1, &
        ax, efstoDyn, asigDyn, enewDyn, enewDyn_nosr, &
        neqp1, neqp1_nosr, ikn, kp,sig%ncore_excl, akih) ! ZL: akih is optional argument
    endif
    
!----------------------------
! Write out matrix elements

    if (peinf%inode.eq.0) then
      write(6,'(a)')
      write(6,'(a)') ' Symmetrized values from band-averaging:'

      ! FHJ: write eqp0.dat and eqp1.dat to units 30 and 31.
      ! ZL: TODO - write separate files for EP
      write(30,'(3(f13.9),i8)') kp%rk(:,ikn), sig%nspin*sig%ndiag
      write(31,'(3(f13.9),i8)') kp%rk(:,ikn), sig%nspin*sig%ndiag
      if (sig%freq_dep==2) then
        do ispin = 1, sig%nspin
          do i = 1, sig%ndiag
            write(30,'(2i8,3f15.9)') ispin, sig%diag(i), &
              wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), efstoDyn(i,ispin)+achcor(i,ispin)
            write(31,'(2i8,3f15.9)') ispin, sig%diag(i), &
              wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), enewDyn(i,ispin)
          enddo
        enddo
      else
        do ispin = 1, sig%nspin
          do i = 1, sig%ndiag
            write(30,'(2i8,3f15.9)') ispin, sig%diag(i), &
              wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), efsto(i,ispin)+dble(achcor(i,ispin))
            write(31,'(2i8,3f15.9)') ispin, sig%diag(i), &
              wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), enew(i,ispin)+dble(achcor(i,ispin))*zrenorm(i,ispin)
          enddo
        enddo
      endif

      if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
        call write_result(kp, wfnk, sig, ax, asx, ach, &
          achcor, asig, alda, efsto, enew, zrenorm, ikn, sig%ncore_excl, akih)
        call write_result_hp(kp, wfnk, sig, ax, asx, ach, achcor, asig, alda, &
          efsto, enew, zrenorm, ikn, sig%ncore_excl, akih)
        if (sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
          do ispin=1,sig%nspin
            do i=1,sig%ndiag
              if (abs(efsto(i,ispin)+SCALARIFY(achcor(i,ispin)) - wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin)).gt.sig%dw) then
                eqp1_warns(1) = .true.
              endif
            enddo
          enddo
        endif
        call ch_converge(kp, sig, COMPLEXIFY(ach_n1), COMPLEXIFY(achcor_n1), ikn)
      endif
      
      if (sig%freq_dep==2) then
        call write_result_dyn(kp, wfnk, sig, ax, asxDyn, &
          achDyn, achDyn_corb, achcor, asigDyn, alda, efstoDyn, enewDyn, &
          enewDyn_nosr, neqp1, neqp1_nosr, ikn,sig%ncore_excl, akih)
        call write_result_dyn_hp(kp, wfnk, sig, ax, asxDyn, achDyn, &
          achDyn_corb, achcor, asigDyn, alda, efstoDyn, enewDyn, enewDyn_nosr, &
          neqp1, neqp1_nosr, ikn,sig%ncore_excl, akih)
        if (any(neqp1<0)) eqp1_warns(2) = .true. !FF extrap
        if (any(neqp1>1)) eqp1_warns(3) = .true. !FF mult solns
        if (any(neqp1==0)) eqp1_warns(4) = .true. !FF no solns
        call ch_converge(kp, sig, achD_n1, COMPLEXIFY(achcor_n1), ikn)
      endif
    endif

    if(peinf%inode == 0) then

!----------------------------
! If not using vxc.dat, create and write Vxc in it

      call timing%start(timing%vxc)
      if(.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX) then
        if(.not.sig%use_kihdat) then ! ZL: only if not using KIH
          call write_matrix_elements_type(120, kp%rk(:, ikn), sig, COMPLEXIFY(alda(:,:)))
        endif
      endif
      call timing%stop(timing%vxc)
      
!----------------------------
! If not using x.dat, create and write X in it

      call timing%start(timing%bare_x)
      if (sig%coul_mod_flag .and. (.not. sig%use_vxc2dat)) then
        call write_matrix_elements_type(121, kp%rk(:, ikn), sig, COMPLEXIFY(ax(:,:)))
      else if((.not. sig%use_xdat) .and. xflag .and. (.not. sig%coul_mod_flag)) then
        ax(:,:) = ax(:,:) / sig%xfrac
        call write_matrix_elements_type(119, kp%rk(:, ikn), sig, COMPLEXIFY(ax(:,:)))
      endif
      call timing%stop(timing%bare_x)
      
    endif
    
    SAFE_DEALLOCATE(indrq)
    SAFE_DEALLOCATE(neq)
    SAFE_DEALLOCATE(itnrq)
    SAFE_DEALLOCATE(rq)
    SAFE_DEALLOCATE(kg0)
    SAFE_DEALLOCATE_P(wfnk%zk)
!#BEGIN_INTERNAL_ONLY
    if(sig%elph) then
      SAFE_DEALLOCATE_P(wfnk_phonq%zk)
    endif
!#END_INTERNAL_ONLY

    if (peinf%inode==0) then
      !FHJ: Flush a couple of units so that the user can retrieve some data if
      !the calculation crashes
      if (sig%freq_dep==2 .or. (sig%fdf==-3 .and. sig%freq_dep==1)) then
        FLUSH(8000)
      endif
      FLUSH(0)
      FLUSH(6)
      FLUSH(8)
    endif
    
  SAFE_DEALLOCATE(wfnkq%occ)
  enddo ! ika (loop over k-points sig%nkn)

  call dealloc_grid(gr)
  call destroy_qran() ! from vcoul_generator

!--------- End loop over rkn points (for which we calculate sigma) -------------

  
  if (peinf%inode==0) then
    do iunit=6,8,2
      write(iunit,'()')
      write(iunit,'(a)') repeat('=', 80)
      write(iunit,'()')
      write(iunit,'(1x,a)') '   n = band index.'
      write(iunit,'(1x,a)') ' Emf = "inner" mean-field energy eigenvalue used to construct Sigma(E),'
      write(iunit,'(1x,a)') '       read from WFN_inner.'
      write(iunit,'(1x,a)') '  Eo = "outer" mean-field energy eigenvalue where we center the evaluation'
      write(iunit,'(1x,a)') '       frequency grid {E} of Sigma(E). Defaults to Emf, unless'
      write(iunit,'(1x,a)') '       you use WFN_outer and eqp_outer.dat / scissors_outer.'
      if (sig%freq_dep==2.and.sig%freq_grid_shift/=2) then
        write(iunit,'(1x,a)') '       (Note: your freq. grid {E} does not contain Eo, so we actually report'
        write(iunit,'(1x,a)') '        below the self energy evaluated at the freq. Eo` closest to Eo.)'
      endif
      write(iunit,'(1x,a)') ' Vxc = exchange-correlation pot., calculated from VXC or read from vxc.dat.'
      ! ZL: for KIH
      write(iunit,'(1x,a)') ' KIH = Kinetic energy + Ionic potential + Hartree, read from kih.dat.'
      write(iunit,'(1x,a)') '   X = bare exchange.'
      if (sig%freq_dep==2.and.sig%freq_dep_method==2) then
        write(iunit,'(1x,a)') ' Res = residue contrib. to Sigma(E) at energy E=Eo.'
        write(iunit,'(1x,a)') ' Int = contrib. to Sigma(E) at energy E=Eo due to integral over imag. freqs.'
        write(iunit,'(1x,a)') ' Cor = Res + Int = correlation portion of Sigma(E) at energy E=Eo.'
        write(iunit,'(1x,a)') ' Sig = X + Cor = self energy Sigma(E) at energy E=Eo.'
      elseif(sig%freq_dep==2) then
        write(iunit,'(1x,a)') '  SX = screened exchange contrib. to Sigma(E) at energy E=Eo.'
        write(iunit,'(1x,a)') '  CH = Coulomb hole contrib. to Sigma(E) at energy E=Eo.'
        write(iunit,'(1x,a)') ' Cor = correlation portion of Sigma(E) at energy E=Eo'
        write(iunit,'(1x,a)') ' Sig = SX + CH = X + Cor = self energy, Sigma(E), at energy E=Eo.'
      else
        write(iunit,'(1x,a)') '  SX = screened exchange contrib. to Sigma(E) at energy E=Eo'
        write(iunit,'(1x,a)') '  CH = Coulomb hole contrib. to Sigma(E) at energy E=Eo'
        write(iunit,'(1x,a)') ' Cor = SX-X + CH = correlation portion of Sigma(E) at energy E=Eo.'
        write(iunit,'(1x,a)') ' Sig = X + Cor = self energy, Sigma(E), at energy E=Eo.'
      endif
      write(iunit,'(1x,a)') 'Eqp0 = on-shell QP energy = Emf - Vxc + Sig(Eo)'
      write(iunit,'(1x,a)') '       Eqp0 is *not* the recommended quantity to use for QP properties.'
      write(iunit,'(1x,a)') 'Eqp1 = off-shell solution to the linearized  Dyson`s equation'
      write(iunit,'(1x,a)') '     = Eqp0 + (dSig/dE) / (1 - dSig/dE) * (Eqp0 - Eo),'
      write(iunit,'(1x,a)') '       or a full linear interpolation if more freq. points where computed.'
      write(iunit,'(1x,a)') '       Eqp1 is the recommended quantity to use for QP properties.'
      if (sig%freq_dep==2) then
        write(iunit,'(1x,a)') 'Soln = Whether we found a solution to Dyson`s equation for Eqp1 (see notes).'
      endif
      write(iunit,'(1x,a)') ' Znk = quasiparticle renormalization factor'
      if (sig%exact_ch==1) then
        write(iunit,'()')
        write(iunit,'(1x,a)') 'Notes on the static remainder:'
        write(iunit,'(1x,a)') '- Unprimed values, such as Eqp0, are evaluated WITH the static remainder'
        write(iunit,'(1x,a)') '- Primed values, such as Eqp0`, are evaluated WITHOUT the static remainder'
      endif
      if (sig%freq_dep==2) then
        write(iunit,'()')
        write(iunit,'(1x,a)') 'Notes on the solutions to Dyson`s equation (Soln):'
        write(iunit,'(1x,a)') '- When we solve Dyson`s eqn. for Eqp1, we might find several or no roots.'
        write(iunit,'(1x,a)') '- If Soln=unique, there is only one solution and the answer is unambiguous'
        write(iunit,'(1x,a)') '- If Soln=extrap{-/+}, we found a solution by extrapolating Sigma(E) to values of'
        write(iunit,'(1x,a)') '  E smaller/greater than our frequency grid. PROCEED WITH CAUTION'
        write(iunit,'(1x,a)') '- If Soln=MULT:, there are multiple solutions to Dyson`s equation. We'
        write(iunit,'(1x,a)') '  report the solution closer to Eqp0. PROCEED WITH CAUTION'
        write(iunit,'(1x,a)') '- If Soln=NO_SOLN!, there is no solution and extrap. is not possible. This'
        write(iunit,'(1x,a)') '  can happen if your frequency grid for Sigma is too small, and Sigma(E) is'
        write(iunit,'(1x,a)') '  ill-behaved. We simply copy Eqp0 to Eqp1, so DON`T TRUST THE ANSWER.'
        write(iunit,'(1x,a)') '- The number of sols. in full-freq. calcs. is affected by self-consistently'
        write(iunit,'(1x,a)') '  updating the eigenvalues in WFN_inner. Consider reruning your calculation'
        write(iunit,'(1x,a)') '  with eqp or scissors corrections if you have many unstable solutions.'
      elseif (sig%freq_dep==1.or.sig%freq_dep==3) then
        write(iunit,'()')
        write(iunit,'(1x,a)') 'Notes on the finite_difference_form from sigma.inp file:'
        write(iunit,'(1x,a)') '  none    : -2 => dSig/dE = 0 (skip the expansion)'
        write(iunit,'(1x,a)') '  backward: -1 => dSig/dE = (Sig(Eo) - Sig(Eo-dE)) / dE'
        write(iunit,'(1x,a)') '  central :  0 => dSig/dE = (Sig(Eo+dE) - Sig(Eo-dE)) / (2*dE)'
        write(iunit,'(1x,a)') '  forward :  1 => dSig/dE = (Sig(Eo+dE) - Sig(Eo)) / dE'
        write(iunit,'(1x,a)') '  default :  2 => forward for diagonal and none for off-diagonal'
        write(iunit,'(1x,a)') '  dE is finite_difference_spacing from Sigma.inp file.'
        write(iunit,'(1x,a,i0,a,f0.3,a)') '  We are using the form #', sig%fdf ,' with dE = ', sig%dw, ' eV.'
      endif
      write(iunit,'()')
      write(iunit,'(1x,a)') 'General notes:'
      write(iunit,'(1x,a)') '- All energies are reported here in eV.'
      write(iunit,'(1x,a)') '- Both Emf and Vxc contain the average pot. Vxc0, so Vxc0 doesn`t affect Sigma.'
      write(iunit,'(1x,a)') '- Eqp1 and Eqp0 are Eqs. (36-37) from Hybertsen & Louie PRB 34 5390.'
      write(iunit,'(1x,a)') '- We recommend you use Eqp1 for QP properties of materials.'
      write(iunit,'()')
      write(iunit,'(a)') repeat('=', 80)
      write(iunit,'()')
    enddo
    if (sig%noffdiag/=0) then
      write(6,776)
      write(8,776)
    endif
    ! FHJ: Flush these units otherwise the warnings that follow will break the
    ! formating of the messages we want to write to unit 6
    FLUSH(0)
    FLUSH(6)
    FLUSH(8)
  endif
776 format( &
      4x,"n = band index of bra wavefunction", &
      /,4x,"m = band index of ket wavefunction", &
      /,4x,"l = band index of energy eigenvalue",/, &
      /,1x,"< psi_n(k) |      X      | psi_m(k) >" &
      /,1x,"< psi_n(k) | SX(Eo_l(k)) | psi_m(k) >" &
      /,1x,"< psi_n(k) | CH(Eo_l(k)) | psi_m(k) >" &
      /,1x,"< psi_n(k) |     Vxc     | psi_m(k) >" &
      /,1x,"< psi_n(k) | T + V_ion + V_H | psi_m(k) >",/)
  
  if (peinf%inode==0) then
    FLUSH(0)
    FLUSH(6)
    FLUSH(8)
    if (eqp1_warns(1)) then
      write(0,'(/a)') "WARNING: |Eqp0 - Eo| > finite_difference_spacing. Linear extrapolation for eqp1"
      write(0,'(a)') "may be inaccurate. You should test the validity of eqp1 by rerunning the"
      write(0,'(a)') "calculation with the self energy evaluated at the eqp0 energies. For that,"
      write(0,'(a)') "use the eqp_outer.dat file, created with eqp.py script and point WFN_outer to"
      write(0,'(a/)') "WFN_inner, if you were not already using WFN_outer."
    endif
    if (eqp1_warns(2)) then
      write(0,'(/a)') "WARNING: Some solutions to Dyson`s equation for Eqp1 were extrapolated."
      write(0,'(a/)') "         PROCEED WITH CAUTION!"
    endif
    if (eqp1_warns(3)) then
      write(0,'(/a)') "WARNING: There are multiple solutions to Dyson`s equation for some values of Eqp1,"
      write(0,'(a)')  "         so we keep those closer to Eqp0."
      write(0,'(a/)') "         PROCEED WITH **EXTREME** CAUTION!"
    endif
    if (eqp1_warns(4)) then
      write(0,'(/a)') "WARNING: Could not find any solution to Dyson`s equation for some values of Eqp1,"
      write(0,'(a)')  "         so we are using Eqp0 for those states."
      write(0,'(a/)') "         PROCEED WITH **EXTREME** CAUTION! SOME VALUES FOR EQP1 ARE JUST *WRONG*!"
    endif
  endif
  if (peinf%inode.eq.0) then
    if (imagvxcflag) write(0,677) 'Vxc'
    if (imagxflag) write(0,677) 'exchange'
  endif
677 format( &
      1x,"WARNING: ",a," diagonal matrix elements have large imaginary part.",/, &
      3x,"This may indicate a problem with your calculation.",/)
  
    
!--------- Clean Up and Finish -------------------------------------------------

  call destroy_fftw_plans()
  call kp%free()

  if (sig%spl_tck%n>0) then
    SAFE_DEALLOCATE_P(sig%spl_tck%t)
    SAFE_DEALLOCATE_P(sig%spl_tck%c)
  endif
  if (sig%spl_tck_outer%n>0) then
    SAFE_DEALLOCATE_P(sig%spl_tck_outer%t)
    SAFE_DEALLOCATE_P(sig%spl_tck_outer%c)
  endif

  SAFE_DEALLOCATE_P(wfnk%isrtk)
  SAFE_DEALLOCATE_P(wfnk%ek)
  SAFE_DEALLOCATE_P(wfnk%elda)
!#BEGIN_INTERNAL_ONLY
  if (sig%elph) then
    SAFE_DEALLOCATE(ep_kp%rk)
    SAFE_DEALLOCATE_P(wfnk_phonq%isrtk)
    SAFE_DEALLOCATE_P(wfnk_phonq%ek)
    SAFE_DEALLOCATE_P(wfnk_phonq%elda)
  endif
!#END_INTERNAL_ONLY
  SAFE_DEALLOCATE_P(sig%qpt)
  if (sig%ndiag.ne.0) then
    SAFE_DEALLOCATE_P(sig%diag)
  end if
  if (sig%noffdiag.gt.0) then
    SAFE_DEALLOCATE_P(sig%off1)
    SAFE_DEALLOCATE_P(sig%off2)
    SAFE_DEALLOCATE_P(sig%off3)
    SAFE_DEALLOCATE_P(sig%offmap)
    if(sig%elph) then
      SAFE_DEALLOCATE_P(sig%off_ep)
      SAFE_DEALLOCATE_P(sig%offmap_ep)
    endif
  endif
  SAFE_DEALLOCATE_P(sig%indkn)
!#BEGIN_INTERNAL_ONLY
  if (sig%elph) then
    SAFE_DEALLOCATE_P(sig%indk_phonq)
    SAFE_DEALLOCATE_P(sig%indk_phonq_g0)
  endif
!#END_INTERNAL_ONLY
  if(.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX) then
    SAFE_DEALLOCATE_P(sig%vxc)
  end if
!#BEGIN_INTERNAL_ONLY
  if (sig%elph) then
    SAFE_DEALLOCATE_P(sig%dvxc)
  endif
!#END_INTERNAL_ONLY
!#BEGIN_INTERNAL_ONLY
  if(.not. sig%use_vxc2dat .and. sig%coul_mod_flag) then
    SAFE_DEALLOCATE_P(sig%vxc2)
  end if
!#END_INTERNAL_ONLY
  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE_P(gvec%index_vec)
!#BEGIN_INTERNAL_ONLY
  if(sig%elph) then
    SAFE_DEALLOCATE_P(ep_gvec%components)
    SAFE_DEALLOCATE_P(ep_gvec%index_vec)
  endif
!#END_INTERNAL_ONLY
  if(sig%freq_dep.eq.1) then
    SAFE_DEALLOCATE_P(wpg%rho)
!#BEGIN_INTERNAL_ONLY
    if(sig%elph) then
      SAFE_DEALLOCATE_P(ep_wpg%rho)
    endif
!#END_INTERNAL_ONLY
  endif
  SAFE_DEALLOCATE_P(wfnkq%isrtkq)
  SAFE_DEALLOCATE_P(wfnkq%ekq)
!#BEGIN_INTERNAL_ONLY
  if(sig%elph) then
    SAFE_DEALLOCATE_P(ep_dq_wfnkq%isrtkq)
    SAFE_DEALLOCATE_P(ep_dq_wfnkq%ekq)
    SAFE_DEALLOCATE_P(wfnkq_phonq%isrtkq)
    SAFE_DEALLOCATE_P(wfnkq_phonq%ekq)
    SAFE_DEALLOCATE_P(ep_dmq_wfnkq_phonq%isrtkq)
    SAFE_DEALLOCATE_P(ep_dmq_wfnkq_phonq%ekq)
  endif
!#END_INTERNAL_ONLY
  SAFE_DEALLOCATE_P(peinf%indext)
  SAFE_DEALLOCATE_P(peinf%indext_dist)
  SAFE_DEALLOCATE_P(peinf%ntband_dist)
  SAFE_DEALLOCATE(alda)
  SAFE_DEALLOCATE(akih)  ! ZL: for KIH
  SAFE_DEALLOCATE(ax)
!#BEGIN_INTERNAL_ONLY
  if(sig%elph) then
    SAFE_DEALLOCATE(ax_ep_one)
    SAFE_DEALLOCATE(ax_ep_two)
    SAFE_DEALLOCATE(ax_ep)
    SAFE_DEALLOCATE(alda_ep)
  endif
!#END_INTERNAL_ONLY
  SAFE_DEALLOCATE(achcor)
  SAFE_DEALLOCATE(asig_imag)
!#BEGIN_INTERNAL_ONLY
  if(sig%elph) then
    SAFE_DEALLOCATE(achcor_ep_one)
    SAFE_DEALLOCATE(asig_imag_ep_one)
    SAFE_DEALLOCATE(achcor_ep_two)
    SAFE_DEALLOCATE(asig_imag_ep_two)
    SAFE_DEALLOCATE(achcor_ep)
    SAFE_DEALLOCATE(asig_imag_ep)
  endif
!#END_INTERNAL_ONLY
  if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
    SAFE_DEALLOCATE_P(achtcor_n1)
!#BEGIN_INTERNAL_ONLY
    if (sig%elph) then
      SAFE_DEALLOCATE_P(achtcor_n1_ep_one)
      SAFE_DEALLOCATE_P(achtcor_n1_ep_two)
    endif
!#END_INTERNAL_ONLY
  endif
  SAFE_DEALLOCATE(achcor_n1)
!#BEGIN_INTERNAL_ONLY
  if (sig%elph) then
    SAFE_DEALLOCATE(achcor_n1_ep_one)
    SAFE_DEALLOCATE(achcor_n1_ep_two)
    SAFE_DEALLOCATE(achcor_n1_ep)
  endif
!#END_INTERNAL_ONLY
  if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
    SAFE_DEALLOCATE(asx)
    SAFE_DEALLOCATE(ach)
    SAFE_DEALLOCATE(asig)
    SAFE_DEALLOCATE_P(acht_n1)
    SAFE_DEALLOCATE(ach_n1)
    SAFE_DEALLOCATE(enew)
    SAFE_DEALLOCATE(efsto)
!#BEGIN_INTERNAL_ONLY
    if(sig%elph) then
      SAFE_DEALLOCATE(asx_ep_one)
      SAFE_DEALLOCATE(ach_ep_one)
      SAFE_DEALLOCATE(asig_ep_one)
      SAFE_DEALLOCATE_P(acht_n1_ep_one)
      SAFE_DEALLOCATE(ach_n1_ep_one)

      SAFE_DEALLOCATE(asx_ep_two)
      SAFE_DEALLOCATE(ach_ep_two)
      SAFE_DEALLOCATE(asig_ep_two)
      SAFE_DEALLOCATE_P(acht_n1_ep_two)
      SAFE_DEALLOCATE(ach_n1_ep_two)

      SAFE_DEALLOCATE(asx_ep)
      SAFE_DEALLOCATE(ach_ep)
      SAFE_DEALLOCATE(asig_ep)
      SAFE_DEALLOCATE(ach_n1_ep)
    endif
!#END_INTERNAL_ONLY
  endif
  if (sig%freq_dep.eq.2) then
    ! ZL: EP not implemented for this frequency dependence
    SAFE_DEALLOCATE(asxDyn)
    SAFE_DEALLOCATE(achDyn)
    SAFE_DEALLOCATE(achDyn_cor)
    SAFE_DEALLOCATE(achDyn_corb)
    SAFE_DEALLOCATE(ach2Dyn)
    SAFE_DEALLOCATE(asigDyn)
    SAFE_DEALLOCATE_P(achtD_n1)
    SAFE_DEALLOCATE(achD_n1)
    SAFE_DEALLOCATE(efstoDyn)
    SAFE_DEALLOCATE(enewDyn)
    SAFE_DEALLOCATE(enewDyn_nosr)
    SAFE_DEALLOCATE(neqp1)
    SAFE_DEALLOCATE(neqp1_nosr)
    SAFE_DEALLOCATE_P(asxtDyn)
    SAFE_DEALLOCATE_P(achtDyn)
    SAFE_DEALLOCATE_P(achtDyn_cor)
    SAFE_DEALLOCATE_P(achtDyn_corb)
    SAFE_DEALLOCATE_P(ach2tDyn)
  endif
  if (sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
    SAFE_DEALLOCATE_P(acht)
    SAFE_DEALLOCATE_P(asxt)
!#BEGIN_INTERNAL_ONLY
    if (sig%elph) then
      SAFE_DEALLOCATE_P(acht_ep_one)
      SAFE_DEALLOCATE_P(asxt_ep_one)
      SAFE_DEALLOCATE_P(acht_ep_two)
      SAFE_DEALLOCATE_P(asxt_ep_two)
    endif
!#END_INTERNAL_ONLY
  endif
  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
    SAFE_DEALLOCATE_P(eps)
  endif
  if (sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
    SAFE_DEALLOCATE_P(epsR)
    if (sig%need_advanced) then
      SAFE_DEALLOCATE_P(epsA)
    endif
  endif
  SAFE_DEALLOCATE(isrtrq)
  if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then
    SAFE_DEALLOCATE_P(isrtrqi)
  endif
  SAFE_DEALLOCATE(ekin)
  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
    SAFE_DEALLOCATE(isrtq)
    SAFE_DEALLOCATE(isrtqi)
    SAFE_DEALLOCATE(ind)
    SAFE_DEALLOCATE(indinv)
    SAFE_DEALLOCATE(ph)
    if (peinf%inode .eq. 0) then
      SAFE_DEALLOCATE_P(epsmpi%isrtq)
      SAFE_DEALLOCATE_P(epsmpi%isrtqi)
    endif
    SAFE_DEALLOCATE_P(epsmpi%qk)
    SAFE_DEALLOCATE_P(epsmpi%nmtx)
    if(sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
      SAFE_DEALLOCATE_P(epsmpi%eps)
    else
      SAFE_DEALLOCATE_P(epsmpi%epsR)
      if (sig%need_advanced) then
        SAFE_DEALLOCATE_P(epsmpi%epsA)
      endif
      if (sig%do_sigma_subspace) then
       SAFE_DEALLOCATE_P(sig%epssub%eps_sub_info)
       SAFE_DEALLOCATE_P(sig%epssub%eigenvec_sub)
       SAFE_DEALLOCATE_P(sig%epssub%eps_sub)
       SAFE_DEALLOCATE_P(sig%epssub%eps_wings_rows)
       SAFE_DEALLOCATE_P(sig%epssub%eps_wings_cols)
       !XXX SAFE_DEALLOCATE_P(sig%epssub%eps_wings_correction_rows)
       !XXX SAFE_DEALLOCATE_P(sig%epssub%eps_wings_correction_cols)
       !MDB name too long, workaround
       if(associated(sig%epssub%eps_wings_correction_rows))then 
         deallocate(sig%epssub%eps_wings_correction_rows)
         nullify(sig%epssub%eps_wings_correction_rows)
       endif
       if(associated(sig%epssub%eps_wings_correction_cols))then 
         deallocate(sig%epssub%eps_wings_correction_cols)
         nullify(sig%epssub%eps_wings_correction_cols)
       endif
       SAFE_DEALLOCATE_P(sig%epssub%vcoul_sub)
      end if
    endif
  endif
  SAFE_DEALLOCATE_P(wfnkqmpi%nkptotal)
  SAFE_DEALLOCATE_P(wfnkqmpi%isort)
  SAFE_DEALLOCATE_P(wfnkqmpi%el)
  SAFE_DEALLOCATE_P(wfnkqmpi%qk)
  SAFE_DEALLOCATE_P(wfnkqmpi%band_index)
  SAFE_DEALLOCATE_P(wfnkqmpi%cg)
!#BEGIN_INTERNAL_ONLY
  if (sig%elph) then
    SAFE_DEALLOCATE_P(ep_dq_wfnkqmpi%nkptotal)
    SAFE_DEALLOCATE_P(ep_dq_wfnkqmpi%isort)
    SAFE_DEALLOCATE_P(ep_dq_wfnkqmpi%el)
    SAFE_DEALLOCATE_P(ep_dq_wfnkqmpi%qk)
    SAFE_DEALLOCATE_P(ep_dq_wfnkqmpi%band_index)
    SAFE_DEALLOCATE_P(ep_dq_wfnkqmpi%cg)
  endif
!#END_INTERNAL_ONLY
  if(sig%nkn.gt.1) then
    SAFE_DEALLOCATE_P(wfnkmpi%nkptotal)
    SAFE_DEALLOCATE_P(wfnkmpi%isort)
    SAFE_DEALLOCATE_P(wfnkmpi%qk)
    SAFE_DEALLOCATE_P(wfnkmpi%el)
    SAFE_DEALLOCATE_P(wfnkmpi%elda)
    SAFE_DEALLOCATE_P(wfnkmpi%cg)
!#BEGIN_INTERNAL_ONLY
    if (sig%elph) then
      SAFE_DEALLOCATE_P(wfnk_phonq_mpi%nkptotal)
      SAFE_DEALLOCATE_P(wfnk_phonq_mpi%isort)
      SAFE_DEALLOCATE_P(wfnk_phonq_mpi%qk)
      SAFE_DEALLOCATE_P(wfnk_phonq_mpi%el)
      SAFE_DEALLOCATE_P(wfnk_phonq_mpi%elda)
      SAFE_DEALLOCATE_P(wfnk_phonq_mpi%cg)
    endif
!#END_INTERNAL_ONLY
  endif
  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
    SAFE_DEALLOCATE_P(epsmpi%inv_igp_index)
  endif
!#BEGIN_INTERNAL_ONLY
  if (sig%subsample) then
    SAFE_DEALLOCATE_P(sig%subweights)
  endif
!#END_INTERNAL_ONLY

  call show_references()

!----------------------------
! Time Accounting

  call timing%stop(timing%total)
  call timing%print(common_timing)

!----------------------------
! Close files and finish

  call close_file(55) ! file sigma.inp
  if(peinf%inode == 0) then
    call close_file(8) ! file sigma_hp.log
    call close_file(30)
    call close_file(31)
    if (sig%coul_mod_flag .and. (.not. sig%use_vxc2dat)) then
      call close_file(121) ! file vxc2.dat
    elseif ((.not. sig%use_xdat) .and. xflag .and. (.not. sig%coul_mod_flag)) then
      call close_file(119) ! file x.dat
    endif
    if(.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX) then
      if(.not.sig%use_kihdat) then  ! ZL: only if not using KIH
        call close_file(120) ! file vxc.dat
      endif
    end if
    if (.not.(sig%freq_dep .eq. 0 .and. sig%exact_ch .eq. 1)  .and. .not. (sig%freq_dep == -1)) then
      call close_file(127) ! file ch_converge.dat
    endif
    if (sig%iwritecoul .eq. 1) then
      call close_file(19) ! file vcoul
    endif
  endif  
  if (peinf%inode.eq.0 .and. (sig%freq_dep.eq.2 .or. (sig%fdf.eq.-3 .and. sig%freq_dep.eq.1))) then
    write(8000, '(/,a)')'# Please refer to Sigma/README for more information about this file.'
    call close_file(8000) ! file spectrum.dat
  endif
  
  call write_memory_usage()

#ifdef HDF5

  if(sig%use_hdf5 .or. sig%wfn_hdf5) call h5close_f(error)

#endif

#ifdef MPI
  call MPI_Finalize(mpierr)

contains

  subroutine xreduce(ar, siz)
    SCALAR, intent(inout) :: ar(*)
    integer, intent(in) :: siz

    PUSH_SUB(sigma.xreduce)
    call MPI_Allreduce(MPI_IN_PLACE, ar, siz, MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD, mpierr)
    POP_SUB(sigma.xreduce)

  end subroutine xreduce

  subroutine zreduce(ar, siz)
    complex(DPC), intent(inout) :: ar(*)
    integer, intent(in) :: siz

    PUSH_SUB(sigma.zreduce)
    call MPI_Allreduce(MPI_IN_PLACE, ar, siz, MPI_COMPLEX_DPC, MPI_SUM, MPI_COMM_WORLD, mpierr)
    POP_SUB(sigma.zreduce)

  end subroutine zreduce
#endif

end program sigma

