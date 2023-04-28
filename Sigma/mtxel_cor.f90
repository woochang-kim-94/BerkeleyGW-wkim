!================================================================================
!
! Module:
!
! (1) mtxel_cor()          Originally By ?       Last Modified 8/15/2015 (FHJ)
!
!     Computes the correlation contribution to the QP self-energy with
!     different approaches (SX+CH, or COR directly), and with different levels
!     of approximations to the dynamical effects of W (static, GPP, real-axis
!     integration and contour deformation).
!
!     More specifically, this routine computes the contribution from the
!     current processor to the correlation part of the self-energy operator
!     between bands `n` (sig%diag(in)) and `m` (sig%diag(im)) at the frequency
!     of the band `l` (sig%diag(il)).
!
!     Screened exchange is put into asxt(iw), Coulomb hole into acht(iw),
!     iw goes from 1 to 3 (see below for what it is).
!
!     On entry, asxt(:) and acht(:) are zeroed.
!
!     Physically, this routine will compute:
!
!     < psi_n ( k ) | Sigma_cor ( E ) | psi_m ( k ) >
!
!     where Sigma_cor(E) = Sigma(E) - Sigma_X(E). We compute the retarded
!     self-energy operator, except for contour-deformation calculations,
!     in which case we compute the time-ordered Sigma.
!
!     Sigma_cor(E) is computed within different approximations:
!     - static COHSEX approximation, both exact and partial sum CH (gsm)
!     - Generalized Plasmon Pole model (?)
!     - full frequency dependent inverse dielectric matrix (CDS/CHP/gsm)
!
!     ch2 variables are used to calculate the CH part of the imaginary part
!     of the self energy assuming a zero energy broadening in the energy
!     denominator for the evaluation of the self energy (CHP).
!
!================================================================================

#include "f_defs.h"

module mtxel_cor_m

  use accel_linalg_m
  use accel_memory_m
  use algos_sigma_m
  use global_m
  use misc_m
  use blas_m
  use timing_m, only: timing => sigma_timing
  use wpeff_m

#if defined (OPENACC)
  ! WPH: Needed for cuda_count_kind type
  use cudafor
#endif

  implicit none

  private

  public :: mtxel_cor

contains

subroutine mtxel_cor(in,il,ispin,ncouls,neps,gvec,eps,ph,ind, &
  indinv,isrtrqi,isrtrq,vcoul,crys,sig,wpg,wfnk,wfnkq,ncoulch, &
  aqsn,aqsm,aqsch,asigt_imag,acht_n1,asxt,acht,achtcor,achtcor_n1,nspin,qk, &
  coulfact,inv_igp_index,ngpown,epsR,epsA, &
  achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,icalc, &
  ep_mtxel, wfnk_ep_braket, ilp)  ! ZL: add the last three variables to deal with EP evaluation at bra and ket

  integer, intent(in) :: in,il,ispin,ncouls,neps
  type (gspace), intent(in) :: gvec
  !> (neps,ngpown) Uninitialized unless we are running the complex version
  SCALAR, pointer, intent(in) :: eps(:,:)
  SCALAR, intent(in) :: ph(:) !< (gvec%ng)
  integer, intent(in) :: ind(:), indinv(:), isrtrq(:) !< (gvec%ng)
  !> (gvec%ng) Uninitialized unless freq_dep=0 .or. exact_ch=1
  integer, pointer, intent(in) :: isrtrqi(:)
  real(DP), intent(in) :: vcoul(:) !< (ncoul)
  type (crystal), intent(in) :: crys
  type (siginfo), intent(in) :: sig
  type (wpgen), intent(in) :: wpg
  type (wfnkstates), intent(in) :: wfnk
  type (wfnkqstates), intent(in) :: wfnkq
  integer, intent(in) :: ncoulch
  SCALAR, intent(in) :: aqsn(:,:), aqsm(:,:) !< (peinf%ntband_max,ncoul)
  !> (ncoulch) Uninitialized unless freq_dep=0 .or. exact_ch=1
  SCALAR, pointer, intent(in) :: aqsch(:)
  ! FHJ: FIXME: make the following pointers allocatable arrays, so that
  ! they can be intent(out) instead of intent(inout).
  !> (sig%ntband) Uninitialized in FF calculations
  SCALAR, pointer, intent(inout) :: acht_n1(:)
  !> (3 or sig%nfreqeval) Uninitialized unless freq_dep=1
  SCALAR, pointer, intent(inout) :: asxt(:), acht(:)
! FHJ: static remainder, resolved over bands
  SCALAR, pointer, intent(inout) :: achtcor_n1(:)
  SCALAR, intent(out) :: achtcor
  SCALAR, intent(out) :: asigt_imag
  integer, intent(in) :: nspin
  real(DP), intent(in) :: qk(:) !< (3)
  real(DP), intent(in) :: coulfact
  integer, intent(in) :: inv_igp_index(:) !< (neps)
  integer, intent(in) :: ngpown
  !> The following pointers are uninitialized unless we are running the complex version
  complex(DPC), pointer, intent(in) :: epsR(:,:,:),epsA(:,:,:) !< (neps,ngpown,sig%nFreq)
  complex(DPC), pointer, intent(inout) :: achtD_n1(:) !< (sig%ntband)
  complex(DPC), pointer, intent(inout) :: asxtDyn(:), achtDyn(:), achtDyn_cor(:), ach2tDyn(:) !< (sig%nfreqeval)
  complex(DPC), pointer, intent(inout) :: achtDyn_corb(:) !< (sig%nfreqeval)
  integer, intent(in) :: icalc
  ! ZL: add for EP braket, it is a optional argument
  !     wfnk_ep_braket must be passed in if ep_mtxel is true
  ! NOTE: even if sig%elph is true, it does NOT mean that the call of mtxel_cor
  !       is for EP, because even sig%elph is true, there are still many calls
  !       that are standard calls for quasiparticle energies, etc.
  !       Therefore, we need to introduce the NEW MANDATORY control flag, ep_mtxel
  logical, intent(in) :: ep_mtxel  ! Flag to tell mtxel_cor if this is for EP offdiag matrix elements
  type (wfnkstates), intent(in), optional :: wfnk_ep_braket
  integer, intent(in), optional :: ilp

  SCALAR, allocatable :: epstemp(:), aqsntemp(:,:),aqsmtemp(:,:)
  complex(DPC), allocatable :: epsRtemp(:,:),epsAtemp(:,:)
  SCALAR, allocatable:: acht_n1_loc(:)

  integer :: my_igp, indigp, ipe, j, iout, my_id
  real(DP) :: diff, diffmin
  logical :: flag_occ
  integer :: ig,igp,igpp,igpp2,iw,iwlda,n1,iband,n1true,nstart,nend,gpp(3)
  real(DP), allocatable :: wx_array(:)
! chs - partial sum static CH, chx - exact static CH
  SCALAR :: achstemp, achxtemp, schx, epsggp, I_epsggp, achxtemp_gp, Omega2, wtilde2
  SCALAR, allocatable :: asxtemp(:),achtemp(:)
  SCALAR, allocatable :: wpmtx(:),I_eps_array(:,:)
  real(DP) :: e_lk, e_n1kq, lambda, phi, freq0, wx, occ
  complex(DPC) :: wtilde, wtilde2_temp, ctemp
  complex(DPC), allocatable :: wtilde_array(:,:)
  ! ZL: add for EP braket
  real(DP) :: e_lk_ep_braket

! full-frequency
  integer :: ifreq, nfreq_real
  real(DP) :: E_max, pref_zb, freqStart, freqEnd
  real(DP), allocatable :: pref(:), wxi(:)
  complex(DPC), allocatable :: asxDtemp(:), achDtemp(:), ach2Dtemp(:), achDtemp_cor(:)
  complex(DPC), allocatable :: achDtemp_corb(:), sch2Di(:), ssxDi(:)
  complex(DPC), allocatable :: schDi(:), schDi_cor(:), schDi_corb(:)
  complex(DPC) :: ssxDit,ssxDitt,ssxDittt,schDt,schDtt
  complex(DPC), allocatable :: epsRggp(:),I_epsRggp(:)
  complex(DPC), allocatable :: I_epsR_array(:,:,:),I_epsA_array(:,:,:)
  SCALAR, allocatable :: I_eps_imag(:,:,:)
  SCALAR :: mcph

! subspace global variables
  logical :: do_sigma_subspace
  integer :: actual_nm
  integer :: ngpown_sub_max, Nbas_own_max
  integer :: ngpown_sub,  Nbas_own, wing_pos
  integer, pointer :: eps_sub_info(:,:)
  complex(DPC), pointer :: eigenvec_sub(:,:)
  complex(DPC), pointer :: eps_sub(:,:,:)
  complex(DPC), pointer :: eps_wings_correction_rows(:,:)
  complex(DPC), pointer :: eps_wings_correction_cols(:,:)
  real(DP), pointer :: vcoul_sub(:)
! subspace local variables
  logical :: my_do_subspace
  logical :: fix_wings, fix_wings_res
  integer :: my_G_start, my_G_end, my_G
  integer :: my_Nb_start, my_Nb_end, my_Nb, Nb_tot, my_ib
  integer :: wing_pos_igp
  integer, allocatable :: ipe_2_LocSubSize(:,:)
  complex(DPC), allocatable :: wings_Re(:,:,:)  !(:,:,1/2) 1=row, 2=col
  complex(DPC), allocatable :: Re_eps_sub(:,:,:)
  complex(DPC) :: wings_Im_elem
  SCALAR, allocatable :: wings_Im(:,:,:)
  SCALAR, allocatable :: Im_eps_sub(:,:,:)
  SCALAR, allocatable :: Msub_m(:,:), Msub_n(:,:)
  SCALAR, allocatable :: Caux(:,:), Caux_send(:,:), Caux_rec(:,:), Maux(:,:)
  SCALAR, allocatable :: n_q_order(:,:), m_q_order(:,:)
  SCALAR, allocatable :: Msub_m_temp(:,:), Msub_n_temp(:,:)
  ! non-blocking cyclic communication
  integer :: ipe_inx
  integer :: isend_static, irec_static
  integer :: actual_send, actual_rec
  integer :: req_send_n, tag_send_n, req_rec_n, tag_rec_n
  integer :: req_send_m, tag_send_m, req_rec_m, tag_rec_m
  integer :: ib_start, ib_end, ib_size, my_num_band
#ifdef MPI
  integer :: stat(MPI_STATUS_SIZE)
#endif
#ifdef OMP
  integer, external :: omp_get_thread_num
#endif

  complex(DPC), allocatable :: mat_1(:,:), mat_2(:,:)

  PUSH_SUB(mtxel_cor)

  ! ZL: check ep_mtxel consistency
  if(ep_mtxel) then
    if(.not. sig%elph) then
      call die('ep_mtxel conflicts sig%elph in mtxel_cor!', only_root_writes=.true.)
    endif

    if(.not. present(wfnk_ep_braket)) then
      call die('To enable ep_mtxel, wfnk_ep_braket MUST be passed into mtxel_cor!', only_root_writes=.true.)
    endif

    if(.not. present(ilp)) then
      call die('To enable ep_mtxel, ilp MUST be passed into mtxel_cor!', only_root_writes=.true.)
    endif

  else
    ! ZL: note that no need to check sig%elph when ep_mtxel is false because even
    ! under sig$elph, there are many standard calls (not for ep_mtxel) such as
    ! quasiparticle energies.
    ! Only need to check consistency with wfnk_ep_braket. Note that if ep_mtxel
    ! if false, wfnk_ep_braket is not used even passed in. But to keep the code
    ! clean and without any confusions, we DO NOT pass in wfnk_ep_braket if
    ! ep_mtxel is false.  So is ilp.
    if(present(wfnk_ep_braket)) then
      call die('When ep_mtxel is false, do NOT pass wfnk_ep_braket to avoid confusion!', only_root_writes=.true.)
    endif

    if(present(ilp)) then
      call die('When ep_mtxel is false, do NOT pass ilp to avoid confusion!', only_root_writes=.true.)
    endif
  endif ! ep_mtxel

  my_do_subspace = .false.
  if(sig%do_sigma_subspace) then
    ! check if all pointers are assosiated
    if(associated(sig%epssub%eps_sub_info) .and. &
       associated(sig%epssub%eigenvec_sub) .and. &
       associated(sig%epssub%eps_sub) .and. &
       associated(sig%epssub%eps_wings_correction_rows) .and. &
       associated(sig%epssub%eps_wings_correction_cols) .and. &
       associated(sig%epssub%vcoul_sub)) then
      my_do_subspace = sig%do_sigma_subspace
    end if
  end if

  if(my_do_subspace) then
    ! 
    actual_nm = sig%epssub%actual_nm
    nullify(eps_sub_info, eigenvec_sub, eps_sub, &
            eps_wings_correction_rows, eps_wings_correction_cols, &
            vcoul_sub)
    eps_sub_info               => sig%epssub%eps_sub_info(:,:,actual_nm)
    eigenvec_sub               => sig%epssub%eigenvec_sub(:,:,actual_nm)
    eps_sub                    => sig%epssub%eps_sub(:,:,:,actual_nm)
    eps_wings_correction_rows  => sig%epssub%eps_wings_correction_rows(:,:)
    eps_wings_correction_cols  => sig%epssub%eps_wings_correction_cols(:,:)
    vcoul_sub                  => sig%epssub%vcoul_sub(:,actual_nm)
    ! 
    wing_pos       = sig%epssub%wing_pos(actual_nm)
    !
    ngpown_sub_max = sig%epssub%ngpown_sub_max
    Nbas_own_max   = sig%epssub%Nbas_own_max  
    ngpown_sub     = sig%epssub%ngpown_sub    
    Nbas_own       = sig%epssub%Nbas_own      
    !
    my_Nb_start = eps_sub_info(1,1)
    my_Nb_end   = eps_sub_info(2,1)
    Nb_tot      = eps_sub_info(3,1)
    my_Nb = my_Nb_end - my_Nb_start + 1 
    !
    my_G_start  = eps_sub_info(1,2)
    my_G_end    = eps_sub_info(2,2)
    my_G = my_G_end - my_G_start + 1

    SAFE_ALLOCATE(ipe_2_LocSubSize, (3, peinf%npes_pool))
    ipe_2_LocSubSize = 0
    ipe_2_LocSubSize(1,peinf%pool_rank+1) = my_Nb_start
    ipe_2_LocSubSize(2,peinf%pool_rank+1) = my_Nb_end
    ipe_2_LocSubSize(3,peinf%pool_rank+1) = my_Nb
#ifdef MPI
    call MPI_Allreduce(MPI_IN_PLACE, ipe_2_LocSubSize(:,:), &
                       3 * peinf%npes_pool, MPI_INTEGER, MPI_SUM, &
                       peinf%pool_comm, mpierr)
#endif

    if( (my_Nb_start .ne. my_G_start) .or. &
        (my_Nb_end   .ne. my_G_end) .or. &
        (my_Nb       .ne. my_G) ) then
      call die("BUG: Subspace Epsilon and eigenvector distribution don't match!")
    end if


    !XX SAFE_ALLOCATE(Caux, (ncouls,Nb_tot))
    !XX !XXXX replicate
    !XX Caux = ZERO
    !XX if(my_G_start > 0) then
    !XX   ! Caux(my_G_start:my_G_end,1:Nb_tot) = eigenvec_sub(1:my_G,1:Nb_tot)
    !XX   Caux(1:ncouls,my_G_start:my_G_end) = eigenvec_sub(1:ncouls,1:my_G)
    !XX end if
    !XX call MPI_Allreduce(MPI_IN_PLACE, Caux(:,:), &
    !XX                    ncouls * Nb_tot, MPI_COMPLEX_DPC, MPI_SUM, &
    !XX                    peinf%pool_comm, mpierr)

    ! matrix elements in subspace
    SAFE_ALLOCATE(Msub_m, (Nb_tot,peinf%ntband_max))
    SAFE_ALLOCATE(Msub_n, (Nb_tot,peinf%ntband_max))
  end if

  if(sig%freq_dep .eq. -1) then
    call die("BUG: cannot call mtxel_cor for Hartree-Fock!")
  endif
  nfreq_real = sig%nFreq - sig%nfreq_imag

  SAFE_ALLOCATE(acht_n1_loc, (sig%ntband))
  if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
    achtcor_n1 = ZERO
  endif

! Initialize Output Arrays
! SIB: Zero contribution to asx(n) and ach(n) for this irq

  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
    asxt(:) = ZERO
    acht(:) = ZERO
    acht_n1(:) = ZERO
    acht_n1_loc(:) = ZERO
  elseif (sig%freq_dep.eq.2) then
    asxtDyn(:) = (0.0d0,0.0d0)
    achtDyn(:) = (0.0d0,0.0d0)
    achtDyn_cor(:) = (0.0d0,0.0d0)
    achtDyn_corb(:) = (0.0d0,0.0d0)
    ach2tDyn(:) = (0.0d0,0.0d0)
    achtD_n1(:) = (0.0d0,0.0d0)
  endif
  achtcor = ZERO
  asigt_imag = ZERO

  if (peinf%my_pool .eq. -1) then
    POP_SUB(mtxel_cor)
    return
  endif


!-----------------------
! Initialization for full-frequency CH integral
  call timing%start(timing%m_cor_init)
  if (sig%freq_dep.eq.2) then
    SAFE_ALLOCATE(pref, (sig%nFreq))
#ifdef CPLX
    pref_zb = 0.5D0 / PI_D
#else
    pref_zb = 1D0 / PI_D
#endif
    do ifreq=1,sig%nFreq
      if (ifreq .lt. sig%nFreq) then
        pref(ifreq)=(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))/PI_D
      else
        pref(ifreq)=pref(ifreq-1)
      endif
    enddo
    pref(1)=pref(1)*0.5d0
    pref(sig%nFreq)=pref(sig%nFreq)*0.5d0
#ifdef CPLX
    do ifreq=1,sig%nFreq
      pref(ifreq)=pref(ifreq)*0.5d0
    enddo
#endif
    E_max=sig%dFreqGrid(sig%nFreq)
    SAFE_ALLOCATE(asxDtemp, (sig%nfreqeval))
    asxDtemp = 0D0
    SAFE_ALLOCATE(achDtemp, (sig%nfreqeval))
    achDtemp = 0D0
    SAFE_ALLOCATE(achDtemp_cor, (sig%nfreqeval))
    achDtemp_cor = 0D0
    SAFE_ALLOCATE(achDtemp_corb, (sig%nfreqeval))
    achDtemp_corb = 0D0
    SAFE_ALLOCATE(ach2Dtemp, (sig%nfreqeval))
    ach2Dtemp = 0D0
    SAFE_ALLOCATE(schDi, (sig%nfreqeval))
    schDi=0D0
    SAFE_ALLOCATE(schDi_cor, (sig%nfreqeval))
    schDi_cor=0D0
    SAFE_ALLOCATE(schDi_corb, (sig%nfreqeval))
    schDi_corb=0D0
    SAFE_ALLOCATE(sch2Di, (sig%nfreqeval))
    sch2Di=0D0
    SAFE_ALLOCATE(ssxDi, (sig%nfreqeval))
    ssxDi=0D0
    SAFE_ALLOCATE(wxi, (sig%nfreqeval))
    wxi=0D0
    if(.not. (my_do_subspace)) then
      ! in the subspace case we allocate this later
      SAFE_ALLOCATE(I_epsR_array, (ncouls,ngpown,nfreq_real))
      if (sig%need_advanced) then
        SAFE_ALLOCATE(I_epsA_array, (ncouls,ngpown,nfreq_real))
      endif
      if (sig%freq_dep_method==2) then
        SAFE_ALLOCATE(I_eps_imag, (ncouls,ngpown,sig%nfreq_imag))
      endif
    end if ! .not.(my_do_subspace)
  else if (sig%freq_dep.eq.1) then
    SAFE_ALLOCATE(I_eps_array, (ncouls,ngpown))
    SAFE_ALLOCATE(wtilde_array, (ncouls,ngpown))
    ! WPH:  The following OpenMP directives are causing segfaults on Nvidia
    !       HPC SDK builds in the buildbot
#ifndef PGI
    !$OMP PARALLEL DO private(j)
#endif
    do j = 1, ngpown
      I_eps_array(:,j)=ZERO
      wtilde_array(:,j)=(0D0,0D0)
    enddo
#ifndef PGI
   !$OMP END PARALLEL DO
#endif
  else if (sig%freq_dep.eq.3) then
    SAFE_ALLOCATE(I_eps_array, (ncouls,ngpown))
    SAFE_ALLOCATE(wtilde_array, (ncouls,ngpown))
    I_eps_array(:,:)=ZERO
    wtilde_array(:,:)=(0D0,0D0)
  else
    SAFE_ALLOCATE(I_eps_array, (ncouls,ngpown))
  endif

  if(my_do_subspace) then
     !XXX ! a bit stupid, just deallocate and reallocate
     !XXX SAFE_DEALLOCATE(I_epsR_array)
     !XXX SAFE_DEALLOCATE(I_eps_imag)
     SAFE_ALLOCATE(I_epsR_array,  (ncouls,ngpown,1))
     SAFE_ALLOCATE(I_eps_imag, (1,1,1))
     ! 
     SAFE_ALLOCATE(Im_eps_sub, (Nb_tot,my_Nb,sig%nfreq_imag))
     SAFE_ALLOCATE(Re_eps_sub, (Nb_tot,my_Nb,nfreq_real))
     Im_eps_sub = ZERO
     Re_eps_sub = (0D0,0D0)
     SAFE_ALLOCATE(wings_Im, (ncouls, sig%nfreq_imag, 2))
     SAFE_ALLOCATE(wings_Re, (ncouls, nfreq_real, 2))
     wings_Im = ZERO
     wings_Re = (0D0,0D0)
     ! reordered matrix elements 
     SAFE_ALLOCATE(n_q_order, (ncouls,peinf%ntband_max))
     SAFE_ALLOCATE(m_q_order, (ncouls,peinf%ntband_max))
     n_q_order = ZERO
     m_q_order = ZERO
  end if

! SIB: the array wfnk%ek has bands indexed by sig%diag
! take the right band index
! DVF : sig%ncore_excl has to be substracted here because wfnk%ek is defined in the 
! read_wavefunction subroutine in input.f90 to be referenced to the case with 
! no core states. 

  iband = sig%diag(il)
  e_lk = wfnk%ek(iband-sig%ncore_excl,ispin)
  ! ZL: for EP, only relevant to GPP
  if(ep_mtxel) then
    ! ZL: use ilp to identify which band
    e_lk_ep_braket = wfnk_ep_braket%ek(sig%diag(ilp)-sig%ncore_excl,ispin)
  endif
  ! FHJ: Figure out starting frequency for freq. grid
  if (sig%freq_grid_shift<2) then
    freq0 = sig%freqevalmin
  else
    freq0 = e_lk - sig%freqevalstep*(sig%nfreqeval-1)/2
    if (sig%freq_dep==2) then
      !FHJ: Avoid hitting a pole
      freq0 = freq0 + TOL_SMALL
    endif
  endif

! JRD: Initial frequencies for plasmon pole case.
! Note: below assumes that we will use iw=1,2,3 in
! subroutine shiftenergy depending on value of sig%fdf

  if (sig%freq_dep.eq.1 .or. sig%freq_dep.eq.0.or.sig%freq_dep.eq.3) then
    if (sig%fdf.eq.-3) then
      SAFE_ALLOCATE(asxtemp,(sig%nfreqeval))
      SAFE_ALLOCATE(achtemp,(sig%nfreqeval))
    else
      SAFE_ALLOCATE(asxtemp,(3))
      SAFE_ALLOCATE(achtemp,(3))
    endif
    asxtemp(:) = ZERO
    achtemp(:) = ZERO
  endif

  if (sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
    if (sig%fdf.eq.-3) then
      nstart=1
      nend=sig%nfreqeval
      SAFE_ALLOCATE(wxi, (sig%nfreqeval))
      SAFE_ALLOCATE(wx_array, (sig%nfreqeval))
      do iw=1,sig%nfreqeval
        wx = freq0 + (iw-1)*sig%freqevalstep
        wxi(iw) = wx
      enddo
    elseif (sig%fdf.eq.-1) then
      nstart = 1
      nend = 2
      SAFE_ALLOCATE(wx_array, (3))
    elseif (sig%fdf.eq.0) then
      nstart = 1
      nend = 3
      SAFE_ALLOCATE(wx_array, (3))
    elseif (sig%fdf.eq.1.or.(sig%fdf.eq.2.and.icalc.eq.1)) then
      nstart = 2
      nend = 3
      SAFE_ALLOCATE(wx_array, (3))
    else
      nstart = 2
      nend = 2
      SAFE_ALLOCATE(wx_array, (3))
    endif
  endif

  ! ZL: consistency check for ep_mtxel on fdf
  if(ep_mtxel) then
    if(sig%fdf .ne. 0) then
      call die('ep_mtxel requires fdf = 0. Add [ finite_difference_form 0 ] into sigma.inp.', only_root_writes=.true.)
    endif
  endif

  call timing%stop(timing%m_cor_init)
!-------------- Main loop over G` (igp) ---------------------------------------
  call timing%start(timing%m_cor_epsinit)
  achxtemp = ZERO

! In below OMP region:
! I_eps_array, I_epsR_array, I_epsA_array, I_eps_imag, wtilde_array are shared

!$OMP PARALLEL private (ig,epsggp,I_epsggp,gpp,iout,schx,igpp,igpp2,achxtemp_gp,igp,indigp, &
!$OMP                      epstemp,epsRtemp,I_epsRggp,epsRggp, &
#ifdef CPLX
!$OMP                      epsAtemp,wtilde2_temp,lambda,phi, &
#endif
!$OMP                      wpmtx,wtilde,wtilde2,Omega2,iw,mcph,ctemp,wings_Im_elem)

#ifdef OMP
  my_id = omp_get_thread_num()
#else
  my_id = 0
#endif

! Allocate Temporary Arrays

  select case(sig%freq_dep)
  case(0)
    SAFE_ALLOCATE(epstemp, (neps))
  case(1)
    SAFE_ALLOCATE(epstemp, (neps))
    SAFE_ALLOCATE(wpmtx, (neps))
  case(2)
    SAFE_ALLOCATE(epstemp, (neps))
  case(3)
    SAFE_ALLOCATE(epstemp, (neps))
    SAFE_ALLOCATE(epsRtemp, (neps,sig%nFreq))
    SAFE_ALLOCATE(epsRggp, (sig%nFreq))
    SAFE_ALLOCATE(I_epsRggp, (sig%nFreq))
  end select

!$OMP DO reduction(+:achxtemp)
  do my_igp=1,ngpown

    indigp = inv_igp_index(my_igp)
    igp = indinv(indigp)

    if (igp .gt. ncouls) write(6,*) "CATASTROPHE", peinf%inode, my_igp, igp
    if (igp .gt. ncouls .or. igp .le. 0) cycle

!!------------- Initialize eps^-1 for this G` ---------------------------------

    select case(sig%freq_dep)
    case(0,1)
      epstemp(:)=eps(:,my_igp)
    case(2)
      epstemp(:)=SCALARIFY(epsR(:,my_igp,1))
    case(3)
      epstemp(:)=SCALARIFY(epsR(:,my_igp,1))
      epsRtemp(:,:)=epsR(:,my_igp,:)

    end select

!------------------------------------------------------------------------------
! (gsm) Below you`ll find the code to compute SX & CH in the static COHSEX
!       approximation (both the exact and partial sum expressions for CH),
!       within the GPP model, and using the full frequency dependent RPA
!       inverse dielectric matrix. The GPP section of the code is well
!       commented (thanks to Sohrab), while the COHSEX and RPA sections
!       are not... But most of the variables are the same anyways,
!       so just look at the comments in the GPP section.
!------------------------------------------------------------------------------

! (gsm) <<<<<< static COHSEX approximation - exact static CH >>>>>>

! < n k | \Sigma_{CH} (r, r`; 0) | m k > =
! \frac{1}{2} \sum_{q G G`}
! < n k | e^{i (G - G`) \cdot r} | m k >
! [\eps_{G G`}^{-1} (q; 0) - \delta_{G G`}] v (q + G`)

    if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then

! Only Computed on One Processor Per Pool

      ! JRD: Since, we are distributed over igp now, all procs need to do this

      achxtemp_gp = ZERO

      do ig=1,ncouls

        epsggp = ph(ig)*MYCONJG(ph(igp))*SCALARIFY(epstemp(ind(ig)))

        if (ig.eq.igp) then
          I_epsggp = epsggp - 1.0d0
        else
          I_epsggp = epsggp
        endif
        if (abs(I_epsggp).lt.TOL_Small) cycle
        if (ig.ne.igp) then
          gpp(:)=gvec%components(:,isrtrq(ig))-gvec%components(:,isrtrq(igp))
          call findvector(iout,gpp,gvec)
          if (iout.eq.0) cycle
          igpp=isrtrqi(iout)
          if (igpp.lt.1.or.igpp.gt.ncoulch) cycle
          gpp(:)=gvec%components(:,isrtrq(igp))-gvec%components(:,isrtrq(ig))
          call findvector(iout,gpp,gvec)
          if (iout.eq.0) cycle
          igpp2=isrtrqi(iout)
          if (igpp2.lt.1.or.igpp2.gt.ncoulch) cycle
        else
          gpp(:)=0
          call findvector(iout,gpp,gvec)
          if (iout.eq.0) cycle
          igpp=isrtrqi(iout)
          if (igpp.lt.1.or.igpp.gt.ncoulch) cycle
        endif
        schx = aqsch(igpp) * I_epsggp
        achxtemp_gp = achxtemp_gp + schx
      enddo ! over G (ig)

      achxtemp = achxtemp + achxtemp_gp * vcoul(igp) * 0.5d0

    endif ! sig%freq_dep.eq.0.or.sig%exact_ch.eq.1

!!-----------------------------------------------------------------------------------------
! (gsm) <<<<<< static COHSEX approximation - CH as a partial sum over empty bands >>>>>>

! < n k | \Sigma_{SX} (r, r`; 0) | m k > =
! - \sum_{n_1}^{occ} \sum_{q G G`}
! < n k | e^{i (q + G) \cdot r} | n_1 k - q >
! < n_1 k - q | e^{- i (q + G`) \cdot r`} | m k >
! [\eps_{G G`}^{-1} (q; 0) - \delta_{G G`}] v (q + G`)

! < n k | \Sigma_{CH} (r, r`; 0) | m k > =
! \frac{1}{2} \sum_{n_1} \sum_{q G G`}
! < n k | e^{i (q + G) \cdot r} | n_1 k - q >
! < n_1 k - q | e^{- i (q + G`) \cdot r`} | m k >
! [\eps_{G G`}^{-1} (q; 0) - \delta_{G G`}] v (q + G`)


    if (sig%freq_dep.eq.0) then

      do ig=1,ncouls
        epsggp = ph(ig)*MYCONJG(ph(igp))*epstemp(ind(ig))
        if (ig.eq.igp) then
          I_epsggp = epsggp - 1.0d0
        else
          I_epsggp = epsggp
        endif

        I_eps_array(ig,my_igp) = I_epsggp

      enddo ! over G (ig)

! (gsm) <<<<<< Generalized Plasmon Pole model >>>>>>

    elseif (sig%freq_dep.eq.1) then

! Zero out temporary accumulation variables

!----------------------
! Calculate Plasma Frequencies
!
! SIB: Here we get the plasmon-pole effective plasma frequecies
! Omega(G,G`)^2 (equation (31) of Hybersten & Louie, PRB 34, 1986, p 5396)
! which come in the vector wp(G) for current G`.

! SIB: We calculate wp(G) for a given G` (trade-off for memory)
! Note that wp(G) G=1,ncouls requires O(N) operations
! Even if we redo it for each band n, not so bad

! JRD: I changed this to use qk instead of qkk because we use vcoul at q=0
! (instead of small q) throughout


! given how many times this routine is called, timing it appears to take a non-negligible amount of time

      call wpeff(crys,gvec,wpg,sig,neps,isrtrq,igp,ncouls,wpmtx,nspin,qk,vcoul,coulfact)

!!------ For given G` (igp), loop over all G vectors (ig) in lower triangle ---

      call timing%start(timing%m_cor_pp_prep)

      do ig=1,ncouls

! Put epsilon(G,G`) into epsggp

        epsggp = ph(ig)*MYCONJG(ph(igp))*epstemp(ind(ig))

! I_epsggp = Kronecker(G,G`) - eps(G,G`)

        if (ig.eq.igp) then
          I_epsggp = ONE - epsggp
        else
          I_epsggp = ZERO - epsggp
        endif

! If I_epsggp is too small, then we skip this (G,G`) entry
! This only happens when eps is 1 on diagonal or 0 off diagonal
! but, this means no screening correct and is already treated properly in bare
! exchange term

        if (abs(I_epsggp).lt.TOL_Small) cycle

        I_eps_array(ig,my_igp) = I_epsggp

! Omega2 = Omega(G,G`)^2 = effective plasma freq^2

        Omega2 = wpmtx(ig)

! If Omega2 is too small, then we skip this (G,G`) entry
! JRD: I am not sure why we have to cycle here... :/ Probably not needed
! FHJ: If Omega2->0, both the SX and the CH terms vanish. But the CH term
! goes to zero as Omega2/wtilde ~ Omega2/sqrt(Omega2). Skipping this term is
! probably better than risking a 0/sqrt(0) division.


        if (abs(Omega2).lt.TOL_Small) cycle

#ifdef CPLX

! <<<<<<<<<<<< COMPLEX GPP >>>>>>>>>>>>

! (gsm) equations (17), (20), (21) from [PRB 40, 3162 (1989)]

        wtilde2_temp = Omega2 / I_epsggp

        lambda = abs(wtilde2_temp)
        if (lambda .lt. TOL_Small) cycle
        phi = atan2(IMAG(wtilde2_temp), dble(wtilde2_temp))
        if (abs(cos(phi)) .lt. TOL_Small) cycle
        wtilde2 = lambda / cos(phi)
! this is not needed because we recalculate Omega2 below
!        Omega2 = Omega2 * CMPLX(1.0d0, -tan(phi))

#else

! <<<<<<<<<<<< REAL GPP >>>>>>>>>>>>

! (gsm) equation (30) from [PRB 34, 5390 (1986)]

        wtilde2 = Omega2 / I_epsggp

        if (abs(wtilde2) .lt. TOL_Small) cycle

#endif

        !FHJ: What do we do if we find an invalid mode with wtilde2<0?
        if (dble(wtilde2)<0) then
          select case (sig%invalid_gpp_mode)
            case (0) ! Skip invalid GPP mode and ignore its contribution to the self energy.
              wtilde = (0d0,0d0)
            case (1) ! "Find" a purely imaginary mode frequency (GSM & JRD).
              wtilde = sqrt(COMPLEXIFY(wtilde2))
            case (2) ! Set the GPP mode frequency to a fixed value of 2 Ry.
              wtilde = CMPLX(2d0*ryd, 0d0)
            case (3) ! Treat that mode within the static COHSEX approximation.
              wtilde = CMPLX(1d0/TOL_ZERO, 0d0)
            case default
              wtilde = CMPLX(1d0/TOL_ZERO, 0d0)
          endselect
        else
          wtilde = CMPLX(sqrt(dble(wtilde2)), 0d0)
        endif
        wtilde_array(ig,my_igp) = wtilde

      enddo ! G Loop for GPP Setup

      call timing%stop(timing%m_cor_pp_prep)

    elseif (sig%freq_dep.eq.3) then

      call timing%start(timing%m_cor_pp_prep)

      do ig=1,ncouls

! Put epsilon(G,G`) into epsRggp
! JRD XXX This is bad. But it is FJR code, so who cares...
        epsRggp(:) = ph(ig)*MYCONJG(ph(igp))*epsRtemp(ind(ig),:)

! I_epsRggp = Kronecker(G,G`) - eps(G,G`)
        if (ig.eq.igp) then
          I_epsRggp(:) = ONE - epsRggp(:)
        else
          I_epsRggp(:) = ZERO - epsRggp(:)
        endif

! If I_epsggp is too small, then we skip this (G,G`) entry
! This only happens when eps is 1 on diagonal or 0 off diagonal
! but, this means no screening correct and is already treated properly in bare
! exchange term

        if (all(abs(I_epsRggp(:)).lt.TOL_Small)) cycle

        I_eps_array(ig,my_igp) = SCALARIFY(I_epsRggp(1))
        wtilde2 = abs(sig%dFreqBrd(2))**2 * I_epsRggp(2)/ ( I_epsRggp(1) - I_epsRggp(2) )

        !FHJ: What do we do if we find an invalid mode with wtilde2<0?
        if (dble(wtilde2)<0) then
          select case (sig%invalid_gpp_mode)
            case (0) ! Skip invalid mode and ignore its contribution to the self energy.
              wtilde = (0d0,0d0)
            case (1) ! "Find" a purely complex mode frequency (GSM & JRD).
              wtilde = sqrt(COMPLEXIFY(wtilde2))
            case (2) ! Set the mode frequency to a fixed value of 2 Ry.
              wtilde = CMPLX(2d0*ryd, 0d0)
            case (3) ! Treat that mode within the static COHSEX approximation
              wtilde = CMPLX(1d0/TOL_ZERO, 0d0)
            case default
              wtilde = CMPLX(1d0/TOL_ZERO, 0d0)
          endselect
        else
          wtilde = sqrt(dble(wtilde2))
        endif
        wtilde_array(ig,my_igp) = wtilde

      enddo ! G Loop for GPP Setup

      call timing%stop(timing%m_cor_pp_prep)

!!--------------------------------------------------------------------
! (gsm) <<<<<< full frequency dependent inverse dielectric matrix >>>>>>

! The code below makes use of the following relations:
!
! {eps_{G G`}^r}^{-1}(q, -E) = {eps_{G G`}^a}^{-1}(q, E)
! for general systems (both real and complex versions of the code)
!
! {eps_{G G`}^a}^{-1}(q, E) = {{eps_{G G`}^r}^{-1}}^{*}(q, E)
! for systems with inversion symmetry (the real version of the code)
! since plane-wave matrix elements are real
!
! {eps_{G G}^a}^{-1}(q, E) = {{eps_{G G}^r}^{-1}}^{*}(q, E)
! for general systems, the diagonal of the matrix (G` = G)
! since plane-wave matrix elements are complex conjugates of each other
!
! The complex version of the code uses
! {eps_{G G`}^r}^{-1}(q, E) and {eps_{G G`}^a}^{-1}(q, E) for E >= 0
!
! The real version of the code uses
! {eps_{G G`}^r}^{-1}(q, E) for E >= 0

! CHP: full frequency routine - the case for sig%ggpsum == 1
!
! On top of the above relations, we need the following:
! Assuming that W_{G,G`}={eps_{G,G`}}^{-1} v(q+G`),
!
! W_{G`,G}^r(E) = {W_{G,G`}^a}^{*}(E)
!               = {W_{G,G`}^r}^{*}(-E) (this form is used if E<0)
! for general systems (both real and complex version of the code)
!
! W_{G`,G}^a(E) = {W_{G,G`}^r}^{*}(E)
!               = {W_{G,G`}^a}^{*}(-E) (this form is used if E<0)
! for general systems (both real and complex version of the code)
!
! W_{G`,G}^r(E) = W_{G,G`}^r(E)
! for systems with inversion symmetry
!
! W_{G`,G}^a(E) = W_{G,G`}^a(E)
! for systems with inversion symmetry
!
! Note that eps^{-1} does not have these symmetries.

    endif
  enddo
!$OMP END DO

  if (sig%freq_dep.eq.2) then
    if (sig%freq_dep_method==0) then

!$OMP DO
      do iw = 1, sig%nFreq
        do my_igp=1,ngpown

          indigp = inv_igp_index(my_igp)
          igp = indinv(indigp)

          if (igp .gt. ncouls .or. igp .le. 0) cycle

          mcph = MYCONJG(ph(igp))
! JRD: XXX ind() may kill us here
          do ig = 1, ncouls
            I_epsR_array(ig,my_igp,iw) = - mcph * ph(ig) * epsR(ind(ig),my_igp,iw)
          enddo
          I_epsR_array(igp,my_igp,iw) = I_epsR_array(igp,my_igp,iw) + 1.0d0
        enddo
      enddo
!$OMP END DO

    else

! JRD XXX
! May want to add schedule dynamic here since iw > nfreq_real has more work. Or move OMP in one level
!$OMP DO
      do iw = 1, sig%nFreq
        if (iw .le. nfreq_real) then

          if(my_do_subspace) then
            ! indigp = inv_igp_index(my_igp) -> local to global index for full Eps mat
            ! indinv() -> from the global Epsilon index to the k+q indexing
            ! ind() -> from the global k+q indexing go th the Epsilon indexing
            ! ph -> phase in the k+q indexing
            if(iw == 1) then
              ! copy full epsilon only for the static case
              I_epsR_array(:,:,1) = (0.0d0, 0.0d0)
              do my_igp=1,ngpown
                indigp = inv_igp_index(my_igp)
                igp = indinv(indigp)

                if (igp .gt. ncouls .or. igp .le. 0) cycle
                !XX if (indigp .gt. ncouls .or. indigp .le. 0) cycle

                mcph = MYCONJG(ph(igp))
                do ig = 1, ncouls
                  I_epsR_array(ig,my_igp,iw) = - mcph * ph(indinv(ig)) * epsR(ig, my_igp, iw)
                enddo
                I_epsR_array(indigp,my_igp,iw) = I_epsR_array(indigp,my_igp,iw) + 1.0d0

              end do
            end if

            ! fix wings, here we go with the ordering driven by the epsilon matrix
            ! each MPI task own the full wings correction, both row and col  
            wing_pos_igp = indinv(wing_pos)
            do ig = 1, ncouls
              ! row
              ! transform eps -> k+q index
              igp = indinv(ig)
              wings_Re(ig,iw,1) = - eps_wings_correction_rows(ig,iw) * MYCONJG(ph(igp)) * ph(wing_pos_igp)
              ! wings_Re(ig,iw,1) = wings_Re(ig,iw,1) * vcoul_sub(ig) / vcoul_sub(wing_pos)
            end do
            do ig = 1, ncouls
              ! column
              igp = indinv(ig)
              wings_Re(ig,iw,2) = - eps_wings_correction_cols(ig,iw) * ph(igp) * MYCONJG(ph(wing_pos_igp))
              ! wings_Re(ig,iw,2) = wings_Re(ig,iw,2) * vcoul_sub(wing_pos) / vcoul_sub(ig)
            end do
            ! and here the famous factor 1.0
            wings_Re(wing_pos,iw,1) = wings_Re(wing_pos,iw,1) + 1.0d0
            ! copy subspace epsilon
            Re_eps_sub(1:Nb_tot,1:my_Nb,iw) = - eps_sub(1:Nb_tot,1:my_Nb,iw)
          
          else ! my_do_subspace

            do my_igp=1,ngpown
  
              indigp = inv_igp_index(my_igp)
              igp = indinv(indigp)
  
              if (igp .gt. ncouls .or. igp .le. 0) cycle
  
              mcph = MYCONJG(ph(igp))
! JRD: XXX ind() may kill us here
              do ig = 1, ncouls
                I_epsR_array(ig,my_igp,iw) = - mcph * ph(ig) * epsR(ind(ig),my_igp,iw)
              enddo
              I_epsR_array(igp,my_igp,iw) = I_epsR_array(igp,my_igp,iw) + 1.0d0
            enddo

          end if ! my_do_subspace 

        endif
        if (iw .gt. nfreq_real) then

          if(my_do_subspace) then
            ! fix wings            
            wing_pos_igp = indinv(wing_pos)
            do ig = 1, ncouls
              ! row
              ! transform eps -> k+q index
              igp = indinv(ig)
              wings_Im_elem = - eps_wings_correction_rows(ig,iw) * MYCONJG(ph(igp)) * ph(wing_pos_igp)
              wings_Im(ig,iw-nfreq_real,1) = SCALARIFY(wings_Im_elem)
              !XXX wings_Im(ig,iw-nfreq_real,1) = wings_Im(ig,iw-nfreq_real,1) * vcoul_sub(ig) / vcoul_sub(wing_pos)
            end do
            do ig = 1, ncouls
              ! column
              igp = indinv(ig)
              wings_Im_elem = - eps_wings_correction_cols(ig,iw) * ph(igp) * MYCONJG(ph(wing_pos_igp))
              wings_Im(ig,iw-nfreq_real,2) = SCALARIFY(wings_Im_elem)
              !XXX wings_Im(ig,iw-nfreq_real,2) = wings_Im(ig,iw-nfreq_real,2) * vcoul_sub(wing_pos) / vcoul_sub(ig)
            end do
            ! and here the famous factor 1.0
            wings_Im(wing_pos,iw-nfreq_real,1) = wings_Im(wing_pos,iw-nfreq_real,1) + 1.0d0
            ! copy subspace epsilon
            Im_eps_sub(1:Nb_tot,1:my_Nb,iw-nfreq_real) = - SCALARIFY(eps_sub(1:Nb_tot,1:my_Nb,iw))
          else ! my_do_subspace

            do my_igp=1,ngpown
  
              indigp = inv_igp_index(my_igp)
              igp = indinv(indigp)
  
              if (igp .gt. ncouls .or. igp .le. 0) cycle
  
              mcph = MYCONJG(ph(igp))
! JRD: XXX ind() may kill us here
              do ig = 1, ncouls
                ctemp = - mcph * ph(ig) * epsR(ind(ig),my_igp,iw)
                I_eps_imag(ig,my_igp,iw-nfreq_real) = SCALARIFY(ctemp)
              enddo
              I_eps_imag(igp,my_igp,iw-nfreq_real) = I_eps_imag(igp,my_igp,iw-nfreq_real) + 1.0d0
            enddo

          end if ! my_do_subspace

        endif
      enddo
!OMP END DO

    endif

    if (sig%need_advanced) then

!$OMP DO
      do iw = 1, sig%nFreq
        do my_igp=1,ngpown

          indigp = inv_igp_index(my_igp)
          igp = indinv(indigp)

          if (igp .gt. ncouls .or. igp .le. 0) cycle

          mcph = MYCONJG(ph(igp))
! JRD: XXX ind() may kill us here
          do ig = 1, ncouls
            I_epsA_array(ig,my_igp,iw) = - mcph * ph(ig) * epsA(ind(ig),my_igp,iw)
          enddo
          I_epsA_array(igp,my_igp,iw) = I_epsA_array(igp,my_igp,iw) + 1.0d0
        enddo
      enddo
!$OMP END DO

    endif
  endif

  if (sig%freq_dep.eq.0) then
    SAFE_DEALLOCATE(epstemp)
  endif
  if (sig%freq_dep.eq.1) then
    SAFE_DEALLOCATE(wpmtx)
    SAFE_DEALLOCATE(epstemp)
  endif
  if (sig%freq_dep.eq.2) then
    SAFE_DEALLOCATE(epstemp)
  endif
  if (sig%freq_dep.eq.3) then
    SAFE_DEALLOCATE(epstemp)
    SAFE_DEALLOCATE(epsRtemp)
    SAFE_DEALLOCATE(epsRggp)
    SAFE_DEALLOCATE(I_epsRggp)
  endif

!$OMP END PARALLEL
  call timing%stop(timing%m_cor_epsinit)

  ! transform basis, moved outside OMP parallel region
  call timing%start(timing%sub_transf_tot)
  if(my_do_subspace) then
    ! transform basis function
    ! reorder wfn matrix elements 
    n_q_order = ZERO
    m_q_order = ZERO
    do ig = 1, ncouls
      igp = indinv(ig)
      if (igp .gt. ncouls .or. igp .le. 0) cycle
      n_q_order(ig,:) = aqsn(igp,:)
      m_q_order(ig,:) = aqsm(igp,:)
    end do

    ! copy matrix elements and scale them with the coulomb operator
    SAFE_ALLOCATE(aqsntemp,(ncouls,peinf%ntband_max))
    SAFE_ALLOCATE(aqsmtemp,(ncouls,peinf%ntband_max))
    aqsntemp = n_q_order
    aqsmtemp = m_q_order
    do ig = 1, ncouls
      igp = indinv(ig)
      if (igp .gt. ncouls .or. igp .le. 0) cycle
      aqsntemp(ig,:) = aqsntemp(ig,:) * ph(igp) * SQRT(vcoul_sub(ig))
      aqsmtemp(ig,:) = aqsmtemp(ig,:) * ph(igp) * vcoul(igp) / SQRT(vcoul_sub(ig))
    end do

    SAFE_ALLOCATE(Caux, (ncouls, Nbas_own_max))
    SAFE_ALLOCATE(Caux_send, (ncouls, Nbas_own_max))
    SAFE_ALLOCATE(Caux_rec, (ncouls, Nbas_own_max))
    SAFE_ALLOCATE(Maux, (Nbas_own_max, peinf%ntband_max))

    Caux      = ZERO  ! (0D0,0D0)
    Caux_send = ZERO  ! (0D0,0D0)
    Caux_rec  = ZERO  ! (0D0,0D0)
    Maux      = ZERO  ! (0D0,0D0)

    ! get ready for the first cycle
    isend_static = MOD(peinf%pool_rank + 1 + peinf%npes_pool, peinf%npes_pool)
    irec_static  = MOD(peinf%pool_rank - 1 + peinf%npes_pool, peinf%npes_pool)
    if(my_Nb_start > 0) then
      Caux(1:ncouls, 1:Nbas_own) = SCALARIFY(eigenvec_sub(1:ncouls, 1:Nbas_own))
    end if
    Msub_n = ZERO  ! (0D0,0D0)
    Msub_m = ZERO  ! (0D0,0D0)

    my_num_band = peinf%ntband_dist(peinf%pool_rank+1)

    ipe = peinf%pool_rank + 1
    do ipe_inx = 1, peinf%npes_pool
      actual_send = MOD(peinf%pool_rank + ipe_inx + peinf%npes_pool, peinf%npes_pool)
      actual_rec  = MOD(peinf%pool_rank - ipe_inx + peinf%npes_pool, peinf%npes_pool)
#ifdef MPI
      call timing%start(timing%sub_transf_com)
      ! post receiving mex
      tag_rec_n = 1
      req_rec_n = MPI_REQUEST_NULL
      CALL MPI_Irecv(Caux_rec, ncouls*Nbas_own_max, MPI_SCALAR, irec_static,&
                     tag_rec_n, peinf%pool_comm, req_rec_n, mpierr)
      ! post send mex
      tag_send_n = 1
      req_send_n = MPI_REQUEST_NULL
      Caux_send  = Caux
      CALL MPI_Isend(Caux_send, ncouls*Nbas_own_max, MPI_SCALAR, isend_static,&
                     tag_send_n, peinf%pool_comm, req_send_n, mpierr)
      call timing%stop(timing%sub_transf_com)
#endif
      ! go with zgemm
      ib_start = ipe_2_LocSubSize(1,ipe) 
      ib_end   = ipe_2_LocSubSize(2,ipe) 
      ib_size  = ipe_2_LocSubSize(3,ipe)
      call timing%start(timing%sub_transf_gemm)
      if(ib_start > 0) then
#ifdef CPLX
        call zgemm('T','N', Nbas_own_max, peinf%ntband_max, ncouls, &
                   (1D0,0D0), Caux(:,:), ncouls, &
                          aqsntemp(:,:), ncouls, &
                   (0D0,0D0), Maux, Nbas_own_max)
#else
        call dgemm('T','N', Nbas_own_max, peinf%ntband_max, ncouls, &
                    1.0D0, Caux(:,:), ncouls, &
                          aqsntemp(:,:), ncouls, &
                    0.0D0, Maux, Nbas_own_max)
#endif
        Msub_n(ib_start:ib_end, 1:my_num_band) = Msub_n(ib_start:ib_end, 1:my_num_band) + &
                                                 Maux(1:ib_size,1:my_num_band)

#ifdef CPLX
        call zgemm('T','N', Nbas_own_max, peinf%ntband_max, ncouls, &
                   (1D0,0D0), Caux(:,:), ncouls, &
                          aqsmtemp(:,:), ncouls, &
                   (0D0,0D0), Maux, Nbas_own_max)
#else
        call dgemm('T','N', Nbas_own_max, peinf%ntband_max, ncouls, &
                   1.0D0, Caux(:,:), ncouls, &
                          aqsmtemp(:,:), ncouls, &
                   0.0D0, Maux, Nbas_own_max)
#endif
        Msub_m(ib_start:ib_end, 1:my_num_band) = Msub_m(ib_start:ib_end, 1:my_num_band) + & 
                                                 Maux(1:ib_size,1:my_num_band)
      end if
      call timing%stop(timing%sub_transf_gemm)
#ifdef MPI
      call timing%start(timing%sub_transf_com)
      CALL MPI_Wait(req_rec_n,stat,mpierr)
      CALL MPI_Wait(req_send_n,stat,mpierr)
      call timing%stop(timing%sub_transf_com)
#endif
      ! get ready for next cycle
      Caux = Caux_rec
      ipe = actual_rec + 1
    end do
    SAFE_DEALLOCATE(aqsntemp)
    SAFE_DEALLOCATE(aqsmtemp)
  end if  ! my_do_subspace
  ! finish with basis transformation
  call timing%stop(timing%sub_transf_tot)

!-----------------------------------------------------------------------------
! End of setup. Begin computation of Sigma.
!-----------------------------------------------------------------------------

  if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
    ! Initialize the static part of the static remainder correction, resolved
    ! over bands. We put the static term into band #1, and later on
    ! substract the CH part. Note the famous factor of 1/2 here.
    achtcor_n1(1) = 0.5d0 * achxtemp
  endif
  SAFE_ALLOCATE(aqsntemp,(ncouls,peinf%ntband_max))
  SAFE_ALLOCATE(aqsmtemp,(ncouls,peinf%ntband_max))

  if (my_do_subspace) then
    ! create buffers for communication
    SAFE_ALLOCATE(Msub_n_temp, (Nb_tot,peinf%ntband_max))
    SAFE_ALLOCATE(Msub_m_temp, (Nb_tot,peinf%ntband_max))
    Msub_n_temp = ZERO
    Msub_m_temp = ZERO
    ! set the flag for fixing the wings
    fix_wings = .true.
    fix_wings_res = .true.
    !XXXX
    if(sig%subsample) then
      ! in the subsample case no fix is done, set flag to false
      fix_wings = .false.
      fix_wings_res = .false.
    end if
    !XXXX
    ! dry run just to fix wings, we pass Msub_m_temp and Msub_n_temp with 
    ! zeros and we pretend to run for our actual proc.
    call timing%start(timing%m_cor_sub_wings)
    aqsntemp = ZERO
    aqsmtemp = ZERO
    aqsntemp(1:ncouls,1:my_num_band) = n_q_order(1:ncouls,1:my_num_band)
    aqsmtemp(1:ncouls,1:my_num_band) = m_q_order(1:ncouls,1:my_num_band)
    ipe = peinf%pool_rank + 1 
    call sigma_cd_subspace()
    call timing%stop(timing%m_cor_sub_wings)
    !XXXX
  end if

  if (allocated(I_eps_array)) then
    call accel_enter_data_map_to(I_eps_array, sigma_gpp_algo)
  end if

  ! FHJ: Loop over each processor/pool, and over each band that ipe owns.
  ! Then, calculate contribution to Sigma due to that band, stored in the
  ! aqs*temp arrays.
  do ipe = 1, peinf%npes_pool

    call timing%start(timing%m_cor_comm)
    if (peinf%pool_rank .eq. ipe-1) then
      if(my_do_subspace) then
        !XXX aqsntemp(:,:) = n_q_order(:,:)
        !XXX aqsmtemp(:,:) = m_q_order(:,:)
        Msub_n_temp(:,:) = Msub_n(:,:)
        Msub_m_temp(:,:) = Msub_m(:,:)
      else
        aqsntemp(:,:) = aqsn(1:ncouls,:)
        aqsmtemp(:,:) = aqsm(1:ncouls,:)
      end if
    endif

#ifdef MPI
    if (peinf%my_pool/=-1) then
      if (my_do_subspace) then
        call MPI_Bcast(Msub_n_temp, Nb_tot * peinf%ntband_max, MPI_SCALAR, ipe-1, &
                       peinf%pool_comm, mpierr)
        call MPI_Bcast(Msub_m_temp, Nb_tot * peinf%ntband_max, MPI_SCALAR, ipe-1, &
                       peinf%pool_comm, mpierr)
      else

        call MPI_Bcast(aqsntemp, peinf%ntband_max*ncouls, MPI_SCALAR, ipe-1, &
          peinf%pool_comm, mpierr)
        ! FHJ: Only bother to BCast aqsmtemp if this is an off-diag. calculation
!#BEGIN_INTERNAL_ONLY
        ! ZL: EP diag also belongs to icalc=2
!#END_INTERNAL_ONLY
        if (icalc==2) then
          call MPI_Bcast(aqsmtemp, peinf%ntband_max*ncouls, MPI_SCALAR, ipe-1, &
            peinf%pool_comm, mpierr)
        else
          aqsmtemp = aqsntemp
        endif

      end if ! my_do_subspace
    endif
#endif
    call timing%stop(timing%m_cor_comm)

    ! here we want to move the band index inside the actual routine such
    ! that we can calculate sigma by calling ZGEMM
    if(my_do_subspace) then
      
      call sigma_cd_subspace()

    else ! my_do_subspace

      ! The accelerated versions of the GPP kernel are done in bulk,
      ! hence the difference in handling.
      if (sigma_gpp_algo /= CPU_ALGO .and. &
          (sig%freq_dep==1 .or. sig%freq_dep==3)) then
        call sigma_gpp()
      else ! sigma_gpp_algo /= CPU_ALGO

        do n1 = 1, peinf%ntband_dist(ipe)

          ! n1true = "True" band index of the band n1 w.r.t. all bands
          n1true = peinf%indext_dist(n1,ipe)
          ! energy of the |n1,k-q> state
          e_n1kq = wfnkq%ekq(n1true,ispin)
          ! occupation of the |n1,k-q> state
          flag_occ = (n1true<=(sig%nvband+sig%ncrit)) &
            .and.((sig%ncrit==0).or.(e_n1kq<=(sig%efermi+sig%tol)))
          if (abs(e_n1kq-sig%efermi)<sig%tol) then
            occ = 0.5d0 ! Fermi-Dirac distribution = 1/2 at Fermi level
          else
            occ = 1.0d0
          endif

          ! FHJ: CH term used for static remainder for current band "n1true".
          achstemp = ZERO

          !FHJ: Call specialized code to calculate contribution to <nk|Sigma^c|mk>
          !due to current band n1
          if (sig%freq_dep==0) then
            call sigma_cohsex()
          else if (sig%freq_dep==1 .or. sig%freq_dep==3) then
            call sigma_gpp()
          else if (sig%freq_dep==2 .and. sig%freq_dep_method==0) then
            call sigma_ra()
          else if (sig%freq_dep==2 .and. sig%freq_dep_method==2) then
            call sigma_cd()
          endif

          if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
            ! Compute the static remainder resolved over bands (just the CH part here).
            ! Note the famous factor of 1/2 here.
            achtcor_n1(n1true) = achtcor_n1(n1true) - 0.5d0*achstemp
          endif
        enddo ! over bands (n1)

      end if ! sigma_gpp_algo

    end if! my_do_subspace

  enddo ! over bands (ipe)

  if (allocated(I_eps_array)) then
    call accel_exit_data_map_delete(I_eps_array, sigma_gpp_algo)
  end if
#ifdef OPENACC
  if (sigma_gpp_algo == OPENACC_ALGO) then
    call acc_clear_freelists()
  end if
#endif

  ! FHJ: Sum contribution to Sigma due to different bands into a single
  ! contribution due to the portion of the dielectric matrix that I own. We`ll
  ! reduce these contributions in sigma_main.f90 after we loop over all q points.

  if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
    ! Compute the static remainder integrated over bands.
    ! Note the famous factor of 1/2 here.
    achtcor = sum(achtcor_n1)
  endif

  if (sig%freq_dep==0) then

    if (sig%exact_ch==1) then
      achtemp(2) = achxtemp
      achxtemp = ZERO
    endif
    asxt(:) = asxt(:) + asxtemp(2)
    acht(:) = acht(:) + achtemp(2)
    if (sig%exact_ch==0) then
      achtcor = achtcor + (achxtemp - achtemp(2))
    endif

  elseif (sig%freq_dep==1 .or. sig%freq_dep==3) then

    do iw=nstart,nend
      asxt(iw) = asxt(iw) + asxtemp(iw)
      acht(iw) = acht(iw) + achtemp(iw)
    enddo

  elseif (sig%freq_dep==2 .and. sig%freq_dep_method==0) then

    do iw = 1, sig%nfreqeval
      asxtDyn(iw) = asxtDyn(iw) + asxDtemp(iw)
      achtDyn(iw) = achtDyn(iw) + achDtemp(iw)
      achtDyn_cor(iw) = achtDyn_cor(iw) + achDtemp_cor(iw)
      achtDyn_corb(iw) = achtDyn_corb(iw) + achDtemp_corb(iw)
      ach2tDyn(iw) = ach2tDyn(iw) + ach2Dtemp(iw)
    enddo

  elseif (sig%freq_dep==2 .and. sig%freq_dep_method==2) then

    do iw = 1, sig%nfreqeval
      asxtDyn(iw) = asxtDyn(iw) + asxDtemp(iw)
      achtDyn(iw) = achtDyn(iw) + achDtemp(iw)
      achtDyn_cor(iw) = achtDyn_cor(iw) + achDtemp(iw)
      achtDyn_corb(iw) = achtDyn_corb(iw) + achDtemp(iw) + asxDtemp(iw)
    enddo

  endif

  SAFE_DEALLOCATE(aqsntemp)
  SAFE_DEALLOCATE(aqsmtemp)


!-----------------------------------------------------------------------------
! Done calculating Sigma! Deallocate and finish.
!-----------------------------------------------------------------------------

  if (sig%freq_dep.eq.1 .or.sig%freq_dep.eq.3) then
    SAFE_DEALLOCATE(wx_array)
    SAFE_DEALLOCATE(wtilde_array)
    if (sig%fdf.eq.-3) then
      SAFE_DEALLOCATE(wxi)
    endif
  endif

  if (sig%freq_dep.eq.1 .or. sig%freq_dep.eq.0 .or. sig%freq_dep.eq.3)  then
    SAFE_DEALLOCATE(asxtemp)
    SAFE_DEALLOCATE(achtemp)
    SAFE_DEALLOCATE(I_eps_array)
    acht_n1(1:sig%ntband) = acht_n1_loc(1:sig%ntband)
  endif

  if (sig%freq_dep.eq.2) then
    SAFE_DEALLOCATE(pref)
    SAFE_DEALLOCATE(asxDtemp)
    SAFE_DEALLOCATE(achDtemp)
    SAFE_DEALLOCATE(achDtemp_cor)
    SAFE_DEALLOCATE(achDtemp_corb)
    SAFE_DEALLOCATE(ach2Dtemp)
    SAFE_DEALLOCATE(schDi)
    SAFE_DEALLOCATE(schDi_cor)
    SAFE_DEALLOCATE(schDi_corb)
    SAFE_DEALLOCATE(sch2Di)
    SAFE_DEALLOCATE(ssxDi)
    SAFE_DEALLOCATE(wxi)
    SAFE_DEALLOCATE(I_epsR_array)
    if (sig%need_advanced) then
      SAFE_DEALLOCATE(I_epsA_array)
    endif
    if (sig%freq_dep_method==2) then
      SAFE_DEALLOCATE(I_eps_imag)
    endif
  endif
  SAFE_DEALLOCATE(acht_n1_loc)

  POP_SUB(mtxel_cor)
  return

contains

!==============================================================================
!                                   COHSEX
!==============================================================================
  !> Calculate Sigma matrix elements, COHSEX formalism
  subroutine sigma_cohsex()
    SCALAR :: aqsn_Ieps, asxtemp_loc, achtemp_loc

    PUSH_SUB(mtxel_cor.sigma_cohsex)

    call timing%start(timing%m_cor_sx_ch)
!$OMP PARALLEL DO private (my_igp,igp,indigp,ig,aqsn_Ieps,achtemp_loc, &
!$OMP                       asxtemp_loc) reduction(+:achtemp,asxtemp,acht_n1_loc)
    do my_igp = 1, ngpown
      indigp = inv_igp_index(my_igp)
      igp = indinv(indigp)
      if (igp>ncouls .or. igp<=0) cycle
      achtemp_loc = ZERO
      asxtemp_loc = ZERO

      do ig = 1, ncouls
        aqsn_Ieps = aqsntemp(ig,n1) * I_eps_array(ig,my_igp)
        if (sig%exact_ch==0) then
          achtemp_loc = achtemp_loc + 0.5d0*aqsn_Ieps
        endif
        if (flag_occ) asxtemp_loc = asxtemp_loc - aqsn_Ieps
      enddo ! ig

      asxtemp(2) = asxtemp(2) + asxtemp_loc*occ*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
      if (sig%exact_ch==0) then
        achtemp(2) = achtemp(2) + achtemp_loc*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
        acht_n1_loc(n1true) = acht_n1_loc(n1true) + achtemp_loc*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
      endif
    enddo ! my_igp
!$OMP END PARALLEL DO
    call timing%stop(timing%m_cor_sx_ch)

    POP_SUB(mtxel_cor.sigma_cohsex)

  end subroutine sigma_cohsex


!==============================================================================
!                                     GPP
!==============================================================================
  !> Calculate Sigma matrix elements, GPP (HL and GN) formalism
  subroutine sigma_gpp()
    PUSH_SUB(mtxel_cor.sigma_gpp)

    select case (sigma_gpp_algo)
    case (CPU_ALGO)
      call sigma_gpp_cpu()
    case (OPENACC_ALGO)
      call sigma_gpp_openacc()
    case (OMP_TARGET_ALGO)
      call sigma_gpp_omp_target()
    case default
      call die("Invald algorithm for sigma_gpp", only_root_writes = .true.)
    end select

    POP_SUB(mtxel_cor.sigma_gpp)

    return
  end subroutine sigma_gpp

  subroutine sigma_gpp_cpu()
    SCALAR, allocatable :: ssx_array(:), sch_array(:)
    real(DP) :: delwr, delw2, wdiffr, rden, wxt, ssxcutoff, limitone, limittwo
    complex(DPC) :: halfinvwtilde, delw, wdiff, cden
    integer :: igbeg, igend, igblk
    SCALAR :: ssx, sch, ssxt, scht, schtt

    PUSH_SUB(mtxel_cor.sigma_gpp_cpu)

    igblk = 512

    ! Some constants used in the loop below, computed here to save
    ! floating point operations
    limitone=1D0/(TOL_Small*4D0)
    limittwo=sig%gamma**2

    ! GSM: compute the static CH for the static remainder
    if (sig%exact_ch.eq.1) then
      call acc_static_ch(ngpown, ncouls, inv_igp_index, indinv, vcoul, &
        aqsntemp(:,n1), aqsmtemp(:,n1), achstemp, eps_scalar=I_eps_array)
    endif

!!!!--- Loop over three energy values which we compute Sigma -----------
!
! SIB: In terms of formulae, the two terms we will compute are given in
! formulae (34a,34b) on p. 5397 of above reference.
!
!                                         Omega^2(G,G`)
!    SX(E) = M(n,G)*conj(M(m,G`)) * ------------------------ * Vcoul(G`)
!                                   (E-E_n1(k-q))^2-wtilde^2
!
!                                           Omega^2(G,G`)
!    CH(E) = M(n,G)*conj(M(m,G`)) * ----------------------------- * Vcoul(G`)
!                                   2*wtilde*[E-E_n1(k-q)-wtilde]
!
! and we are evaluating both at E = E_l(k) and E = E_l(k) +/- dE.
! SX only gets contributions for the band n1 being an occupied state.
!
! For diagonal elements, we need matrix elements at various
! energies (to get quasi-particle energies), i.e. for iw=1,2,3;
! but for off-diagonal elements, we only need them at iw=2
! (the actual energy) so we won`t even bother calculating
! them at other energies

    call timing%start(timing%m_cor_sx_ch)
    do iw=nstart,nend
! wx = E_l(k) - E_n1(k-q) + dE = difference in energies
! of the two states appearing as E - E_n1(k-q) above.
      if (sig%fdf .eq. -3) then
        wx_array(iw) = wxi(iw) - e_n1kq
      else
        ! ZL: change the 3 GPP wx_array elements for EP
        if(ep_mtxel) then
          ! since sig%elph enforces GPP, we are certain that there are 3 values
          if(iw .eq. 2) wx_array(iw) = e_lk - e_n1kq                          ! evaluate at selected bra or ket
          if(iw .eq. 1) wx_array(iw) = e_lk_ep_braket - e_n1kq                ! evaluate at the other (unselected) ket or bra
          if(iw .eq. 3) wx_array(iw) = (e_lk + e_lk_ep_braket)/2.0 - e_n1kq   ! evaluate at the average energy of bra and ket
        else          
          wx_array(iw) = e_lk + sig%dw*(iw-2) - e_n1kq
        endif ! ep_mtxel
      endif
      if (abs(wx_array(iw)) .lt. TOL_Zero) wx_array(iw) = TOL_Zero
    enddo

! JRD: This Loop is Performance critical. Make Sure you don`t mess it up

!$OMP PARALLEL private (my_igp,igp,indigp,ssx_array,sch_array,ig, &
!$OMP                      wtilde,wtilde2,halfinvwtilde,ssxcutoff,sch,ssx, &
!$OMP                      iw,delw,delw2,Omega2,scht,schtt,ssxt,wxt, &
!$OMP                      rden,cden,delwr,wdiffr,wdiff,igbeg,igend)

    if (sig%fdf.eq.-3) then
      SAFE_ALLOCATE(ssx_array,(sig%nfreqeval))
      SAFE_ALLOCATE(sch_array,(sig%nfreqeval))
    else
      SAFE_ALLOCATE(ssx_array,(3))
      SAFE_ALLOCATE(sch_array,(3))
    endif

!$OMP DO reduction(+:asxtemp,acht_n1_loc,achtemp)
    do my_igp = 1, ngpown

      indigp = inv_igp_index(my_igp)
      igp = indinv(indigp)

      if (igp .gt. ncouls .or. igp .le. 0) cycle

      ssx_array = ZERO
      sch_array = ZERO

! delw measures how close to zero the difference
! wx - wtilde = E - E_n1(k-q) - wtilde is relative to wtilde:
! delw = (E - E_n1(k-q) - wtilde) / (2 * wtilde)
! delw2 is the squared absolute value of delw

! If delw is small, both SX and CH blow up, but their sum (for
! an occupied state n1) is finite.  In such a case, their sum
! is Omega^2 / (4 * wtilde2) / (1 + delw). Then the sum is
! assigned to ssx and sch is set to zero.

! If ssx is too large (which can happen due to the pole at
! wx + wtilde = 0 of the SX term), then we should drop this term.
! See the discussion at the bottom of p. 5411-5412 of Hybertsen & Louie.

! If G.neq.G`, then since we sum over only lower triangle,
! we include the contribution we would have had from (G`,G).

      if (flag_occ) then

        do iw=nstart,nend

          scht=0D0
          ssxt=0D0
          wxt = wx_array(iw)

          do ig = 1, ncouls

! Here we recompute Omega2 = wtilde2 * I_eps_array. This contains
! the factor of (1 - i tan phi) from Eqs. 21 & 22 of arXiv paper.

!FIXME: Here we use temporary variables wtilde, wtilde2, Omega2 while
! in the following sections we use wtilde_array and I_eps_array directly.
! JRD, please write a comment here explaining whether this is critical
! for performance or it doesn`t matter.

            wtilde = wtilde_array(ig,my_igp)
            wtilde2 = wtilde**2
            Omega2 = wtilde2 * I_eps_array(ig,my_igp)

! Cycle bad for vectorization. Not needed wtilde is zero
!                  if (abs(Omega2) .lt. TOL_Zero) cycle

            wdiff = wxt - wtilde

            cden = wdiff
            rden = cden * CONJG(cden)
            rden = 1D0 / rden
            delw = wtilde * CONJG(cden) * rden
            delwr = delw*CONJG(delw)
            wdiffr = wdiff*CONJG(wdiff)

! This Practice is bad for vectorization and understanding of the output.
! JRD: Complex division is hard to vectorize. So, we help the compiler.
            if (wdiffr.gt.limittwo .and. delwr.lt.limitone) then
              sch = delw * I_eps_array(ig,my_igp)
              cden = wxt**2 - wtilde2
              rden = cden*CONJG(cden)
              rden = 1D0 / rden
              ssx = Omega2 * CONJG(cden) * rden
            else if ( delwr .gt. TOL_Zero) then
              sch = 0.0d0
              cden = (4.0d0 * wtilde2 * (delw + 0.5D0 ))
              rden = cden*MYCONJG(cden)
              rden = 1D0 / rden
              ssx = -Omega2 * MYCONJG(cden) * rden * delw
            else
              sch = 0.0d0
              ssx = 0.0d0
            endif

! JRD: Breaks vectorization. But, I will have to fix later because
! leaving it out breaks GSM example.
            ssxcutoff = sig%sexcut*abs(I_eps_array(ig,my_igp))
            if (abs(ssx) .gt. ssxcutoff .and. wxt .lt. 0.0d0) ssx=0.0d0

            ssxt = ssxt + aqsntemp(ig,n1)*ssx
            scht = scht + aqsntemp(ig,n1)*sch

          enddo ! loop over g
          ssx_array(iw) = ssx_array(iw) + ssxt*MYCONJG(aqsmtemp(igp,n1))
          sch_array(iw) = sch_array(iw) + 0.5D0*scht*MYCONJG(aqsmtemp(igp,n1))
        enddo

      else

        do igbeg = 1,ncouls,igblk
        igend = min(igbeg+igblk-1,ncouls)
        do iw=nstart,nend

          scht=0D0
          ssxt=0D0
          wxt = wx_array(iw)

!dir$ no unroll
          do ig = igbeg, igend

! Here we recompute Omega2 = wtilde2 * I_eps_array. This contains
! the factor of (1 - i tan phi) from Eqs. 21 & 22 of arXiv paper.

!FIXME: Here we use wtilde_array and I_eps_array directly while in the
! previous sections we use temporary variables wtilde, wtilde2, Omega2.
! JRD, please write a comment here explaining whether this is critical
! for performance or it doesn`t matter.

            wdiff = wxt - wtilde_array(ig,my_igp)

            cden = wdiff
            rden = cden * CONJG(cden)
            rden = 1D0 / rden
            delw = wtilde_array(ig,my_igp) * CONJG(cden) * rden
            delwr = delw*CONJG(delw)
            wdiffr = wdiff*CONJG(wdiff)

            schtt = aqsntemp(ig,n1) * delw * I_eps_array(ig,my_igp)

! JRD: This if is OK for vectorization
            if (wdiffr.gt.limittwo .and. delwr.lt.limitone) then
              scht = scht + schtt
            endif

          enddo ! loop over g

          sch_array(iw) = sch_array(iw) + 0.5D0*scht*MYCONJG(aqsmtemp(igp,n1))

        enddo
        enddo

      endif

! If a valence band, then accumulate SX contribution.

      if (flag_occ) then
        do iw=nstart,nend
          asxtemp(iw) = asxtemp(iw) - ssx_array(iw) * occ * vcoul(igp)
        enddo
      endif

! Accumulate CH contribution.

      do iw=nstart,nend
        achtemp(iw) = achtemp(iw) + sch_array(iw) * vcoul(igp)
      enddo

! Logging CH convergence.

      acht_n1_loc(n1true) = acht_n1_loc(n1true) + sch_array(2) * vcoul(igp)

    enddo ! igp
!$OMP END DO
    SAFE_DEALLOCATE(ssx_array)
    SAFE_DEALLOCATE(sch_array)
!$OMP END PARALLEL
    call timing%stop(timing%m_cor_sx_ch)

    POP_SUB(mtxel_cor.sigma_gpp_cpu)

  end subroutine sigma_gpp_cpu

  subroutine sigma_gpp_openacc()
#ifdef OPENACC
    SCALAR, allocatable :: ssx_array(:), sch_array(:)
    SCALAR, allocatable :: aqsmtemp_local(:,:)
    SCALAR, allocatable :: mat_stat_local(:,:)
    real(DP) :: delwr, delw2, wdiffr, rden, ssxcutoff, limitone, limittwo
    complex(DPC) :: halfinvwtilde, delw, wdiff, cden
    real(DP) :: cden_mod_square
    SCALAR :: ssx_array_1, ssx_array_2, ssx_array_3, sch_array_1, sch_array_2, sch_array_3
    integer :: igbeg, igend, igblk
    SCALAR :: ssx, sch, schtt, matngmatmgp, epsa, aqsmconj
    integer :: dummy(32)
    integer :: ntband_dist
    real(DP), allocatable :: occ_array(:)
    real(DP), allocatable :: wx_array_t(:,:)
    integer, allocatable :: n1true_array(:)
    real(DP) :: nvband, ncrit, efermi, tol, sexcut
    real(DP) :: wxt, vcoulx, vcoulxocc, e_n1kq
    SCALAR :: upd1, upd2
    real(DP), allocatable :: vcoul_loc(:)
    integer(kind=8) ::  my_igp, n1_loc, ig
    integer :: cuErr
    integer(kind=cuda_count_kind) :: freemem, totmem
    integer   :: stat
    integer wta_blk, wta_blks, wta_blke, wta_blksize, wta_blknum, wta_blksize_last
    integer :: ig_blksize, ig_blk
    integer :: n1loc_blksize, n1loc_blk

    PUSH_SUB(mtxel_cor.sigma_gpp_openacc)

    call accel_enter_data_map_to(aqsntemp, sigma_gpp_algo)

    ! Some constants used in the loop below, computed here to save
    ! floating point operations
    limitone=1D0/(TOL_Small*4D0)
    limittwo=sig%gamma**2

!Static
    ! GSM: compute the static CH for the static remainder
    if ( sig%exact_ch.eq.1 ) then

      call timing%start(timing%m_cor_remain)

      SAFE_ALLOCATE(mat_stat_local, (ngpown,peinf%ntband_max))
      ntband_dist = peinf%ntband_dist(ipe)

      call accel_enter_data_map_alloc(mat_stat_local, sigma_gpp_algo)
      call accel_xgemm('T', 'N', ngpown, ntband_dist, ncouls, &
                     ONE, I_eps_array, ncouls, aqsntemp, ncouls,&
                     ZERO, mat_stat_local, ngpown, &
                     sigma_gpp_algo)
      call accel_update_from(mat_stat_local, sigma_gpp_algo)

      !$omp parallel do private(n1_loc, my_igp, n1true, indigp, igp)
      do n1_loc = 1, ntband_dist
        n1true = peinf%indext_dist(n1_loc,ipe)
        do my_igp = 1, ngpown
          indigp = inv_igp_index(my_igp)
          igp = indinv(indigp)
          if (igp>ncouls .or. igp<=0) cycle
          achtcor_n1(n1true) = achtcor_n1(n1true) + MYCONJG(aqsmtemp(igp,n1_loc))* mat_stat_local(my_igp,n1_loc) *vcoul(igp)*0.25D0
        end do
      end do
      !$omp end parallel do

      call accel_exit_data_map_delete(mat_stat_local, sigma_gpp_algo)
      call acc_clear_freelists()

      SAFE_DEALLOCATE(mat_stat_local)

      call timing%stop(timing%m_cor_remain)
    end if  ! static reminder

    call timing%start(timing%m_cor_sx_ch)

    ssx_array_1 = ZERO
    sch_array_1 = ZERO
    ssx_array_2 = ZERO
    sch_array_2 = ZERO
    ssx_array_3 = ZERO
    sch_array_3 = ZERO

    SAFE_ALLOCATE(wx_array_t,(peinf%ntband_dist(ipe),3))

    SAFE_ALLOCATE(aqsmtemp_local, (peinf%ntband_dist(ipe),ngpown))
    aqsmtemp_local = ZERO

    do n1_loc = 1, peinf%ntband_dist(ipe)
      ! n1true = "True" band index of the band n1 w.r.t. all bands
      n1true = peinf%indext_dist(n1_loc,ipe) ! changed to input
      ! energy of the |n1,k-q> state
      e_n1kq = wfnkq%ekq(n1true,ispin)
      do iw=nstart,nend
        if (sig%fdf .eq. -3) then
          wx_array_t(n1_loc,iw) = wxi(iw) - e_n1kq
        else
          if(ep_mtxel) then
            ! since sig%elph enforces GPP, we are certain that there are 3 values
            if(iw .eq. 2) wx_array_t(n1_loc,iw) = e_lk - e_n1kq                          ! evaluate at selected bra or ket
            if(iw .eq. 1) wx_array_t(n1_loc,iw) = e_lk_ep_braket - e_n1kq                ! evaluate at the other (unselected) ket or bra
            if(iw .eq. 3) wx_array_t(n1_loc,iw) = (e_lk + e_lk_ep_braket)/2.0 - e_n1kq   ! evaluate at the average energy of bra and ket
          else
            wx_array_t(n1_loc,iw) = e_lk + sig%dw*(iw-2) - e_n1kq
          endif ! ep_mtxel
        endif
        if (abs(wx_array_t(n1_loc,iw)) .lt. TOL_Zero) wx_array_t(n1_loc,iw) = TOL_Zero
      enddo

      ! fill the aqsmtemp_local array
      do my_igp = 1,  ngpown
        indigp = inv_igp_index(my_igp)
        igp = indinv(indigp)
        if (igp .le. ncouls .and. igp .gt. 0) then
          aqsmtemp_local(n1_loc,my_igp) = aqsmtemp(igp,n1_loc)
        end if
      end do
    enddo

    ntband_dist = peinf%ntband_dist(ipe)

    SAFE_ALLOCATE(n1true_array, (ntband_dist))
    SAFE_ALLOCATE(occ_array, (ntband_dist))

    ! variables
    nvband = sig%nvband
    ncrit  = sig%ncrit
    efermi = sig%efermi
    tol    = sig%tol
    sexcut = sig%sexcut

    do n1_loc = 1, ntband_dist
      n1true_array(n1_loc) = peinf%indext_dist(n1_loc,ipe)
      e_n1kq = wfnkq%ekq(n1true_array(n1_loc),ispin)
      occ_array(n1_loc) = 0.0d0
      !XXXX if (peinf%indext_dist(n1_loc,ipe) .le. nvband) then
      !XXXX This is the fix for metals (partial occ), in case there are problems, start investigating here.
      if ( (peinf%indext_dist(n1_loc,ipe) .le. (sig%nvband+sig%ncrit)) &
        .and.( (sig%ncrit==0) .or. (e_n1kq<=(sig%efermi+sig%tol)) ) ) then
        if (abs(e_n1kq-efermi)<tol) then
          occ_array(n1_loc) = 0.5d0 ! Fermi-Dirac distribution = 1/2 at Fermi level
        else
          occ_array(n1_loc) = 1.0d0
        endif
      endif
    enddo

    SAFE_ALLOCATE(vcoul_loc, (ngpown))
    vcoul_loc = 0.0_dp
    do my_igp = 1,  ngpown
      indigp = inv_igp_index(my_igp)
      igp = indinv(indigp)
      if (igp .le. ncouls .and. igp .gt. 0) then
        vcoul_loc(my_igp) = vcoul(igp)
      end if
    end do

! This is the Outter loop

    stat = cudaMemGetInfo(freemem, totmem)
    wta_blksize = (freemem * 0.7 - ntband_dist * (sizeof(n1true_array(1)) + DP*2 + 3*DP +size_of_SCALAR*ngpown)-ngpown*DP) / (ncouls * DPC * 2)

    if (wta_blksize .gt. ngpown) wta_blksize = ngpown
    if (wta_blksize .le. 0) wta_blksize = 1
    wta_blknum = ngpown / wta_blksize
    wta_blksize_last = mod(ngpown,wta_blksize)
    if (wta_blksize_last .gt.0) wta_blknum=wta_blknum+1

    ig_blksize = 256
    n1loc_blksize = 16
    !$acc data copyin(n1true_array, occ_array, wx_array_t, aqsmtemp_local, vcoul_loc)
    do wta_blk = 1, wta_blknum
      wta_blks=(wta_blk-1)*wta_blksize+1
      wta_blke=min(wta_blk*wta_blksize,ngpown)

      !$acc data copyin(wtilde_array(:,wta_blks:wta_blke))
      do iw=nstart,nend
        !$acc parallel present(I_eps_array, aqsntemp) vector_length(512)
        !$acc loop gang vector reduction(+: ssx_array_3, sch_array_3) collapse(3)
        do n1loc_blk = 1, n1loc_blksize
          do my_igp = wta_blks, wta_blke
            do ig_blk = 1, ig_blksize
              !$acc loop seq
              do ig = ig_blk, ncouls, ig_blksize
                wtilde = wtilde_array(ig,my_igp)
                wtilde2 = wtilde**2
                epsa    = I_eps_array(ig,my_igp)
                Omega2 = wtilde2 * epsa
                ssxcutoff = sexcut**2 * epsa * MYCONJG(epsa)
                vcoulx = vcoul_loc(my_igp)

                !$acc loop seq
                do n1_loc = n1loc_blk, ntband_dist, n1loc_blksize
                  aqsmconj = MYCONJG(aqsmtemp_local(n1_loc,my_igp))
                  matngmatmgp =  aqsmconj * aqsntemp(ig,n1_loc)
                  vcoulxocc = vcoulx * occ_array(n1_loc)

                  wxt = wx_array_t(n1_loc,iw)
                  wdiff = wxt - wtilde
                  wdiffr = wdiff*CONJG(wdiff)
                  rden = 1d0 / wdiffr
                  delw = wtilde*CONJG(wdiff)*rden
                  delwr = delw*CONJG(delw)

                  sch = 0.0d0
                  ssx = 0.0d0
                  if (wdiffr.gt.limittwo .and. delwr.lt.limitone) then
                    sch = delw * epsa
                    cden = wxt**2 - wtilde2
                    rden = cden*CONJG(cden)
                    rden = 1D0 / rden
                    ssx = Omega2 * CONJG(cden) * rden
                  else if ( delwr .gt. TOL_Zero) then
                    cden = (4.0d0 * wtilde2 * (delw + 0.5D0 ))
                    rden = cden*MYCONJG(cden)
                    rden = 1D0 / rden
                    ! ssx = -Omega2 * MYCONJG(cden) * rden * delw
                    wdiff = -Omega2 * CONJG(cden)
                    ssx =  rden * delw * wdiff
                  endif
                  upd1 = 0.0d0
                  rden = ssx * MYCONJG(ssx)
                  if (rden .le. ssxcutoff .or. wxt .ge. 0.0d0) then
                    upd1 = vcoulxocc * ssx * matngmatmgp
                  endif
                  upd2 = vcoulx * sch * matngmatmgp * 0.5d0

                  ssx_array_3 = ssx_array_3 + upd1
                  sch_array_3 = sch_array_3 + upd2
                end do ! loop over n1_loc
              end do ! loop over ig
            end do ! loop over ig_blk
          end do ! loop over igp
          ! NEED FIX
          ! acht_n1_loc(n1true_array(n1_loc)) = acht_n1_loc(n1true_array(n1_loc)) + sch_array_1 * vcoul_loc(my_igp)
        end do ! n1loc_blk
        !$acc end parallel

        if (iw.eq.1) then
          ssx_array_1 = ssx_array_3
          sch_array_1 = sch_array_3
          ssx_array_3 = ZERO
          sch_array_3 = ZERO
        else if (iw.eq.2) then
          ssx_array_2 = ssx_array_3
          sch_array_2 = sch_array_3
          ssx_array_3 = ZERO
          sch_array_3 = ZERO
        end if
      end do ! iw
      !$acc end data

      ! WPH: For streams, do not delete!
      !do n1_loc = 1, ntband_dist
      !  !$acc wait(n1_loc)
      !end do
    end do ! wta_blk
    !$acc end data

    SAFE_DEALLOCATE(wx_array_t)
    SAFE_DEALLOCATE(n1true_array)
    SAFE_DEALLOCATE(vcoul_loc)
    SAFE_DEALLOCATE(aqsmtemp_local)

    asxtemp(1) = asxtemp(1) - ssx_array_1
    asxtemp(2) = asxtemp(2) - ssx_array_2
    asxtemp(3) = asxtemp(3) - ssx_array_3

    achtemp(1) = achtemp(1) + sch_array_1
    achtemp(2) = achtemp(2) + sch_array_2
    achtemp(3) = achtemp(3) + sch_array_3

    call timing%stop(timing%m_cor_sx_ch)

    call accel_exit_data_map_delete(aqsntemp, sigma_gpp_algo)
#else
    PUSH_SUB(mtxel_cor.sigma_gpp_openacc)

    call die("OpenACC version of sigma_gpp requested, but OpenACC not compiled"&
            &" into this executable", only_root_writes = .true.)
#endif

    POP_SUB(mtxel_cor.sigma_gpp_openacc)
  end subroutine sigma_gpp_openacc

  subroutine sigma_gpp_omp_target()
#ifdef OMP_TARGET
    SCALAR, allocatable :: ssx_array(:), sch_array(:)
    SCALAR, allocatable :: aqsmtemp_local(:,:)
    SCALAR, allocatable :: mat_stat_local(:,:)
    real(DP) :: delwr, delw2, wdiffr, rden, ssxcutoff, limitone, limittwo
    complex(DPC) :: halfinvwtilde, delw, wdiff, cden
    real(DP) :: cden_mod_square
    SCALAR :: ssx_array_1, ssx_array_2, ssx_array_3, sch_array_1, sch_array_2, sch_array_3
    integer :: igbeg, igend, igblk
    SCALAR :: ssx, sch, schtt, matngmatmgp, epsa, aqsmconj
    integer :: dummy(32)
    integer :: ntband_dist
    real(DP), allocatable :: occ_array(:)
    real(DP), allocatable :: wx_array_t(:,:)
    integer, allocatable :: n1true_array(:)
    real(DP) :: nvband, ncrit, efermi, tol, sexcut
    real(DP) :: wxt, vcoulx, vcoulxocc, e_n1kq
    SCALAR :: upd1, upd2
    real(DP), allocatable :: vcoul_loc(:)
    integer(kind=8) ::  my_igp, n1_loc, ig
    integer   :: stat
    integer wta_blk, wta_blks, wta_blke, wta_blksize, wta_blknum, wta_blksize_last
    integer :: ig_blksize, ig_blk
    integer :: n1loc_blksize, n1loc_blk

    real(DP) :: sexcut_sq, ssx_sq
    complex(DPC) :: ssx_array_temp, sch_array_temp, ssx_array_loc, sch_array_loc

    PUSH_SUB(mtxel_cor.sigma_gpp_omp_target)

    call accel_enter_data_map_to(aqsntemp, sigma_gpp_algo)

    ! Some constants used in the loop below, computed here to save
    ! floating point operations
    limitone=1D0/(TOL_Small*4D0)
    limittwo=sig%gamma**2

!Static
    ! GSM: compute the static CH for the static remainder
    if ( sig%exact_ch.eq.1 ) then

      call timing%start(timing%m_cor_remain)

      SAFE_ALLOCATE(mat_stat_local, (ngpown,peinf%ntband_max))
      ntband_dist = peinf%ntband_dist(ipe)

      call accel_enter_data_map_alloc(mat_stat_local, sigma_gpp_algo)
      call accel_xgemm('T', 'N', ngpown, ntband_dist, ncouls, &
                     ONE, I_eps_array, ncouls, aqsntemp, ncouls,&
                     ZERO, mat_stat_local, ngpown, &
                     sigma_gpp_algo)
      call accel_update_from(mat_stat_local, sigma_gpp_algo)

      do n1_loc = 1, ntband_dist
        n1true = peinf%indext_dist(n1_loc,ipe)
        do my_igp = 1, ngpown
          indigp = inv_igp_index(my_igp)
          igp = indinv(indigp)
          if (igp>ncouls .or. igp<=0) cycle
          achtcor_n1(n1true) = achtcor_n1(n1true) + MYCONJG(aqsmtemp(igp,n1_loc))* mat_stat_local(my_igp,n1_loc) *vcoul(igp)*0.25D0
        end do
      end do
      call accel_exit_data_map_delete(mat_stat_local, sigma_gpp_algo)

      SAFE_DEALLOCATE(mat_stat_local)

      call timing%stop(timing%m_cor_remain)
    end if  ! static reminder

    call timing%start(timing%m_cor_sx_ch)


    if ( use_gpu ) then
      SAFE_ALLOCATE(wx_array_t,(peinf%ntband_dist(ipe),3))
      SAFE_ALLOCATE(aqsmtemp_local, (peinf%ntband_dist(ipe),ngpown))
      aqsmtemp_local = ZERO
    else
      SAFE_ALLOCATE(wx_array_t,(3,peinf%ntband_dist(ipe)))
      SAFE_ALLOCATE(aqsmtemp_local, (ngpown,peinf%ntband_dist(ipe)))
      aqsmtemp_local = ZERO
    end if

    do n1_loc = 1, peinf%ntband_dist(ipe)
      ! n1true = "True" band index of the band n1 w.r.t. all bands
      n1true = peinf%indext_dist(n1_loc,ipe) ! changed to input
      ! energy of the |n1,k-q> state
      e_n1kq = wfnkq%ekq(n1true,ispin)

      if ( use_gpu ) then
        do iw=nstart,nend
          if (sig%fdf .eq. -3) then
            wx_array_t(n1_loc,iw) = wxi(iw) - e_n1kq
          else
            if(ep_mtxel) then
              ! since sig%elph enforces GPP, we are certain that there are 3 values
              if(iw .eq. 2) wx_array_t(n1_loc,iw) = e_lk - e_n1kq                          ! evaluate at selected bra or ket
              if(iw .eq. 1) wx_array_t(n1_loc,iw) = e_lk_ep_braket - e_n1kq                ! evaluate at the other (unselected) ket or bra
              if(iw .eq. 3) wx_array_t(n1_loc,iw) = (e_lk + e_lk_ep_braket)/2.0 - e_n1kq   ! evaluate at the average energy of bra and ket
            else
              wx_array_t(n1_loc,iw) = e_lk + sig%dw*(iw-2) - e_n1kq
            endif ! ep_mtxel
          endif
          if (abs(wx_array_t(n1_loc,iw)) .lt. TOL_Zero) wx_array_t(n1_loc,iw) = TOL_Zero
        end do
      else
        do iw=nstart,nend
          if (sig%fdf .eq. -3) then
            wx_array_t(iw,n1_loc) = wxi(iw) - e_n1kq
          else
            if(ep_mtxel) then
              ! since sig%elph enforces GPP, we are certain that there are 3 values
              if(iw .eq. 2) wx_array_t(iw,n1_loc) = e_lk - e_n1kq                          ! evaluate at selected bra or ket
              if(iw .eq. 1) wx_array_t(iw,n1_loc) = e_lk_ep_braket - e_n1kq                ! evaluate at the other (unselected) ket or bra
              if(iw .eq. 3) wx_array_t(iw,n1_loc) = (e_lk + e_lk_ep_braket)/2.0 - e_n1kq   ! evaluate at the average energy of bra and ket
            else
              wx_array_t(iw,n1_loc) = e_lk + sig%dw*(iw-2) - e_n1kq
            endif ! ep_mtxel
          endif
          if (abs(wx_array_t(iw,n1_loc)) .lt. TOL_Zero) wx_array_t(iw,n1_loc) = TOL_Zero
        enddo
      end if

      ! fill the aqsmtemp_local array
      do my_igp = 1,  ngpown
        indigp = inv_igp_index(my_igp)
        igp = indinv(indigp)
        if (igp .le. ncouls .and. igp .gt. 0) then
          if ( use_gpu ) then
            aqsmtemp_local(n1_loc,my_igp) = aqsmtemp(igp,n1_loc)
          else
            aqsmtemp_local(my_igp,n1_loc) = aqsmtemp(igp,n1_loc)
          end if
        end if
      end do
    enddo

    ntband_dist = peinf%ntband_dist(ipe)

    SAFE_ALLOCATE(n1true_array, (ntband_dist))
    SAFE_ALLOCATE(occ_array, (ntband_dist))

    ! variables
    nvband = sig%nvband
    ncrit  = sig%ncrit
    efermi = sig%efermi
    tol    = sig%tol
    sexcut = sig%sexcut
    sexcut_sq = sexcut * sexcut

    do n1_loc = 1, ntband_dist
      n1true_array(n1_loc) = peinf%indext_dist(n1_loc,ipe)
      e_n1kq = wfnkq%ekq(n1true_array(n1_loc),ispin)
      occ_array(n1_loc) = 0.0d0
      !XXXX if (peinf%indext_dist(n1_loc,ipe) .le. nvband) then
      !XXXX This is the fix for metals (partial occ), in case there are problems, start investigating here.
      if ( (peinf%indext_dist(n1_loc,ipe) .le. (sig%nvband+sig%ncrit)) &
        .and.( (sig%ncrit==0) .or. (e_n1kq<=(sig%efermi+sig%tol)) ) ) then
        if (abs(e_n1kq-efermi)<tol) then
          occ_array(n1_loc) = 0.5d0 ! Fermi-Dirac distribution = 1/2 at Fermi level
        else
          occ_array(n1_loc) = 1.0d0
        endif
      endif
    enddo

    SAFE_ALLOCATE(vcoul_loc, (ngpown))
    vcoul_loc = 0.0_dp
    do my_igp = 1,  ngpown
      indigp = inv_igp_index(my_igp)
      igp = indinv(indigp)
      if (igp .le. ncouls .and. igp .gt. 0) then
        vcoul_loc(my_igp) = vcoul(igp)
      end if
    end do

    if (sig%fdf.eq.-3) then
      SAFE_ALLOCATE(ssx_array,(sig%nfreqeval))
      SAFE_ALLOCATE(sch_array,(sig%nfreqeval))
    else
      SAFE_ALLOCATE(ssx_array,(3))
      SAFE_ALLOCATE(sch_array,(3))
    endif

    ssx_array = ZERO
    sch_array = ZERO

    ! This is the Outter loop
    if( use_gpu ) then
      n1loc_blksize = sig%gpp_band_block_size
      n1loc_blksize = MIN(ntband_dist,n1loc_blksize)
      ig_blksize = sig%gpp_ig_block_size
      ig_blksize = MIN(ig_blksize, ncouls)

      !$omp target data map(to: wtilde_array, wx_array_t, aqsmtemp_local, vcoul_loc, occ_array )
      do iw=nstart,nend
        ssx_array_temp = ZERO
        sch_array_temp = ZERO

        !$omp target teams loop thread_limit(128) default(none) &
        !$omp&         shared(n1loc_blksize, ntband_dist, ngpown, ig_blksize, ncouls, iw, nstart, nend, &
        !$omp&               wtilde_array, I_eps_array, wx_array_t, &
        !$omp&               limittwo, limitone, sexcut, sexcut_sq, &
        !$omp&               aqsmtemp_local, aqsntemp, vcoul_loc, occ_array, TOL_Zero ) &
        !$omp&         private(n1loc_blk, n1_loc, my_igp, ig_blk, ig, wtilde, wtilde2, Omega2, wdiff, delw, delwr, rden, wdiffr, sch, cden, &
        !$omp&                 ssx, ssx_sq, matngmatmgp, ssxcutoff, ssx_array_loc, sch_array_loc ) &
        !$omp&         map(tofrom: ssx_array_temp, sch_array_temp) &
        !$omp& collapse(3) reduction(+:sch_array_temp, ssx_array_temp)
        do n1loc_blk = 1, n1loc_blksize
          do my_igp = 1,  ngpown
            do ig_blk = 1, ig_blksize
              ssx_array_loc = ZERO
              sch_array_loc = ZERO

              !$OMP loop bind(thread)
              do ig = ig_blk, ncouls, ig_blksize
                wtilde = wtilde_array(ig,my_igp)
                wtilde2 = wtilde**2
                Omega2 = wtilde2 * I_eps_array(ig,my_igp)
                ssxcutoff = sexcut_sq * I_eps_array(ig,my_igp) * MYCONJG(I_eps_array(ig,my_igp))

                !$OMP loop bind(thread)
                do n1_loc = n1loc_blk, ntband_dist, n1loc_blksize
                  wdiff = wx_array_t(n1_loc,iw) - wtilde

                  !XXX delw = wtilde / wdiff ! replace with reciprocal
                  wdiffr = wdiff*CONJG(wdiff)
                  rden = 1.0D+00 / wdiffr
                  delw = wtilde*CONJG(wdiff)*rden

                  delwr = delw*CONJG(delw)
                  wdiffr = wdiff*CONJG(wdiff)

                  if (wdiffr.gt.limittwo .and. delwr.lt.limitone) then
                    sch = delw * I_eps_array(ig,my_igp)
                    cden = wx_array_t(n1_loc,iw)**2 - wtilde2
                    rden = cden * CONJG(cden)
                    rden = 1.0D+00 / rden
                    ssx  = Omega2 * CONJG(cden) * rden
                  else if ( delwr .gt. TOL_Zero) then
                    sch = ZERO
                    cden = (4.0d0 * wtilde2 * (delw + 0.5D0 ))
                    rden = cden * CONJG(cden)
                    rden = 1.0D+00 / rden
                    ssx  = -Omega2 * delw * CONJG(cden) * rden
                  else
                    sch = ZERO
                    ssx = ZERO
                  endif

                  matngmatmgp =  MYCONJG(aqsmtemp_local(n1_loc,my_igp)) * aqsntemp(ig,n1_loc)

                  ssx_sq    = ssx * MYCONJG(ssx)
                  if ( ssx_sq .gt. ssxcutoff .and. wx_array_t(n1_loc,iw) .lt. 0.0d0) ssx=ZERO

                  ssx_array_loc = ssx_array_loc + vcoul_loc(my_igp)* occ_array(n1_loc) *ssx * matngmatmgp
                  sch_array_loc = sch_array_loc + vcoul_loc(my_igp)*sch * matngmatmgp *0.5d0
                end do ! n1_loc
              end do ! ig

              ssx_array_temp = ssx_array_temp + ssx_array_loc
              sch_array_temp = sch_array_temp + sch_array_loc
              !XXX need to be fixed for reduction
              !XXX n1true = n1true_array(n1_loc)
              !XXX acht_n1_loc(n1true) = acht_n1_loc(n1true) + sch_array(2) ! * v1828coul_loc(my_igp)
            end do ! ig_blk
          end do ! my_igp
        end do  ! n1loc_blk
        !$omp end target teams loop

        ssx_array(iw) = ssx_array(iw)+ssx_array_temp
        sch_array(iw) = sch_array(iw)+sch_array_temp
      end do ! iw
      !$omp end target data

    else ! .not. use_gpu
      ! CPU version
      !$omp parallel do default(private) &
      !$omp         shared(ntband_dist, ngpown, ncouls, nstart, nend, &
      !$omp               wtilde_array, I_eps_array, wx_array_t, &
      !$omp               limittwo, limitone, n1true_array, &
      !$omp               aqsmtemp_local, aqsntemp, vcoul_loc, occ_array ) &
      !$omp collapse(3) reduction(+:acht_n1_loc,ssx_array,sch_array)
      do n1_loc = 1, ntband_dist
        do my_igp = 1,  ngpown
          do ig = 1, ncouls
            do iw=nstart,nend
              wtilde = wtilde_array(ig,my_igp)
              wtilde2 = wtilde**2
              Omega2 = wtilde2 * I_eps_array(ig,my_igp)

              wdiff = wx_array_t(iw,n1_loc) - wtilde

              delw = wtilde / wdiff
              delwr = delw*CONJG(delw)
              wdiffr = wdiff*CONJG(wdiff)

              if (wdiffr.gt.limittwo .and. delwr.lt.limitone) then
                sch = delw * I_eps_array(ig,my_igp)
                cden = wx_array_t(iw,n1_loc)**2 - wtilde2
                ssx = Omega2 / cden
              else if ( delwr .gt. TOL_Zero) then
                sch = 0.0d0
                cden = (4.0d0 * wtilde2 * (delw + 0.5D0 ))
                ssx = -Omega2 * delw / cden
              else
                sch = 0.0d0
                ssx = 0.0d0
              endif

              matngmatmgp =  MYCONJG(aqsmtemp_local(my_igp,n1_loc)) * aqsntemp(ig,n1_loc)

              ssxcutoff = sexcut*abs(I_eps_array(ig,my_igp))
              if (abs(ssx) .gt. ssxcutoff .and. wx_array_t(iw,n1_loc) .lt. 0.0d0) ssx=0.0d0

              ssx_array(iw) = ssx_array(iw)+vcoul_loc(my_igp)* occ_array(n1_loc) *ssx * matngmatmgp
              sch_array(iw) = sch_array(iw)+vcoul_loc(my_igp)*sch * matngmatmgp *0.5d0
            end do  ! iw

            n1true = n1true_array(n1_loc)
            acht_n1_loc(n1true) = acht_n1_loc(n1true) + sch_array(2) ! * v1828coul_loc(my_igp)
          end do ! ig
        end do ! my_igp
      end do  ! n1_loc
      !$omp end parallel do
    end if

    do iw=nstart,nend
      achtemp(iw) = achtemp(iw) + sch_array(iw)
      asxtemp(iw) = asxtemp(iw) - ssx_array(iw)
    enddo

    SAFE_DEALLOCATE(ssx_array)
    SAFE_DEALLOCATE(sch_array)

    SAFE_DEALLOCATE(wx_array_t)
    SAFE_DEALLOCATE(n1true_array)
    SAFE_DEALLOCATE(vcoul_loc)
    SAFE_DEALLOCATE(aqsmtemp_local)

    call accel_exit_data_map_delete(aqsntemp, sigma_gpp_algo)
#else
    PUSH_SUB(mtxel_cor.sigma_gpp_openacc)

    call die("OpenMP Target version of sigma_gpp requested, but OpenMP Target&
             & not compiled into this executable", only_root_writes = .true.)
#endif

    POP_SUB(mtxel_cor.sigma_gpp_omp_target)
  end subroutine sigma_gpp_omp_target


!==============================================================================
!                                  REAL AXIS
!==============================================================================
  !> Calculate Sigma matrix elements, full-freq / REAL-AXIS formalism
  subroutine sigma_ra()
    complex(DPC) :: I_epsRggp_int, I_epsAggp_int
    real(DP) :: cedifft_zb,intfact,cedifft_zb_left,cedifft_zb_right
    complex(DPC) :: cedifft_coh, cedifft_cor
    complex(DPC) :: schDt_avg, schDt_right, schDt_left, schDt_lin, schDt_lin2, schDt_lin3
    complex(DPC) :: schDttt, schDttt_cor, schD, sch2Dt, sch2Dtt
    complex(DPC), allocatable :: schDt_array(:)
    real(DP) :: fact1, fact2
    integer :: ijk
    logical, save :: warned=.false.

    PUSH_SUB(mtxel_cor.sigma_ra)

    ! JRD: Find iw closest to e_lk
    diffmin = INF
    do iw=1,sig%nfreqeval
      diff = abs(freq0 + (iw-1)*sig%freqevalstep - e_lk)
      if (diff .lt. diffmin) then
        diffmin=diff
        iwlda=iw
      endif
    enddo

    do iw=1,sig%nfreqeval
      wx = freq0 - e_n1kq + (iw-1)*sig%freqevalstep
      wxi(iw) = wx
    enddo

    ! GSM: compute the static CH for the static remainder
    if (sig%exact_ch.eq.1) then
      call acc_static_ch(ngpown, ncouls, inv_igp_index, indinv, vcoul, &
        aqsntemp(:,n1), aqsmtemp(:,n1), achstemp, eps_cplx=I_epsR_array(:,:,1))
    endif

    ssxDi = (0D0,0D0)
    schDi = (0D0,0D0)
    schDi_cor = (0D0,0D0)
    schDi_corb = (0D0,0D0)
    sch2Di = (0D0,0D0)

    call timing%start(timing%m_cor_ra_sx)
    ! JRD: Don`t want to thread here, nfreqeval could be small
    do iw=1,sig%nfreqeval
      wx = wxi(iw)
      ! SX and CH terms: equation (1.42) of Catalin`s thesis
      ! Note the negative sign in I_epsRggp and I_epsAggp

      if (flag_occ) then

        if(wx.ge.0.0d0) then
          ifreq=0
          do ijk = 1, sig%nFreq-1
            if (wx .ge. sig%dFreqGrid(ijk) .and. wx .lt. sig%dFreqGrid(ijk+1)) then
              ifreq=ijk
            endif
          enddo
          if (ifreq .eq. 0) then
            ifreq = sig%nFreq+3 ! three is for luck
          endif
        else
          ifreq=0
          do ijk = 1, sig%nFreq-1
            if (-wx .ge. sig%dFreqGrid(ijk) .and. -wx .lt. sig%dFreqGrid(ijk+1)) then
              ifreq=ijk
            endif
          enddo
          if (ifreq .eq. 0) then
            ifreq = sig%nFreq+3 ! three is for luck
          endif
        endif

        if(ifreq.ge.sig%nFreq) then
          if (.not.warned .and. peinf%pool_rank==0) then
            write(0,777) iband, n1true, e_lk, e_n1kq, wx, E_max
          endif
          warned = .true.
          ifreq=sig%nFreq-1
        endif

#ifdef CPLX
        if(wx.ge.0.d0) then

          fact1 = (sig%dFreqGrid(ifreq+1)-wx)/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))
          fact2 = (wx-sig%dFreqGrid(ifreq))/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))

          ssxDittt = 0D0

!$OMP PARALLEL do private (my_igp,igp,indigp,ssxDitt,ig, &
!$OMP                       ssxDit) reduction(+:ssxDittt)
          do my_igp = 1, ngpown
            indigp = inv_igp_index(my_igp)
            igp = indinv(indigp)

            if (igp .gt. ncouls .or. igp .le. 0) cycle

            ssxDitt = (0D0,0D0)
            do ig = 1, ncouls
              ssxDit=I_epsR_array(ig,my_igp,ifreq)*fact1 + &
              I_epsR_array(ig,my_igp,ifreq+1)*fact2
              ssxDitt = ssxDitt + aqsntemp(ig,n1)*ssxDit
            enddo
            ssxDittt = ssxDittt + ssxDitt*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
          enddo
!$OMP END PARALLEL DO

          ssxDi(iw) = ssxDi(iw) + ssxDittt

        else

          fact1 = (sig%dFreqGrid(ifreq+1)+wx)/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))
          fact2 = (-sig%dFreqGrid(ifreq)-wx)/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))

          ssxDittt = 0D0

!$OMP PARALLEL do private (my_igp,igp,indigp,ssxDitt,ig, &
!$OMP                       ssxDit) reduction(+:ssxDittt)
          do my_igp = 1, ngpown
            indigp = inv_igp_index(my_igp)
            igp = indinv(indigp)

            if (igp .gt. ncouls .or. igp .le. 0) cycle

            ssxDitt = (0D0,0D0)
            do ig = 1, ncouls
              ssxDit=I_epsA_array(ig,my_igp,ifreq)*fact1+ &
                I_epsA_array(ig,my_igp,ifreq+1)*fact2
              ssxDitt = ssxDitt + aqsntemp(ig,n1)*ssxDit
            enddo
            ssxDittt = ssxDittt + ssxDitt*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
          enddo
!$OMP END PARALLEL DO

          ssxDi(iw) = ssxDi(iw) + ssxDittt

        endif
#else
        if(wx.ge.0.d0) then

          fact1 = (sig%dFreqGrid(ifreq+1)-wx)/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))
          fact2 = (wx-sig%dFreqGrid(ifreq))/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))

          ssxDittt = 0D0

!$OMP PARALLEL do private (my_igp,igp,indigp,ssxDitt,ig, &
!$OMP                       ssxDit) reduction(+:ssxDittt)
          do my_igp = 1, ngpown
            indigp = inv_igp_index(my_igp)
            igp = indinv(indigp)

            if (igp .gt. ncouls .or. igp .le. 0) cycle

            ssxDitt = (0D0,0D0)
            do ig = 1, ncouls
              ssxDit=I_epsR_array(ig,my_igp,ifreq)*fact1+ &
                I_epsR_array(ig,my_igp,ifreq+1)*fact2
              ssxDitt = ssxDitt + aqsntemp(ig,n1)*ssxDit
            enddo
            ssxDittt = ssxDittt + ssxDitt*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
          enddo
!$OMP END PARALLEL DO

          ssxDi(iw) = ssxDi(iw) + ssxDittt

        else

          fact1 = (sig%dFreqGrid(ifreq+1)+wx)/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))
          fact2 = (-sig%dFreqGrid(ifreq)-wx)/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))

          ssxDittt = 0D0

!$OMP PARALLEL do private (my_igp,igp,indigp,ssxDitt,ig, &
!$OMP                       ssxDit) reduction(+:ssxDittt)
          do my_igp = 1, ngpown
            indigp = inv_igp_index(my_igp)
            igp = indinv(indigp)

            if (igp .gt. ncouls .or. igp .le. 0) cycle

            ssxDitt = (0D0,0D0)
            do ig = 1, ncouls
              ssxDit=CONJG(I_epsR_array(ig,my_igp,ifreq))*fact1 + &
                CONJG(I_epsR_array(ig,my_igp,ifreq+1))*fact2
              ssxDitt = ssxDitt + aqsntemp(ig,n1)*ssxDit
            enddo
            ssxDittt = ssxDittt + ssxDitt*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
          enddo
!$OMP END PARALLEL DO

          ssxDi(iw) = ssxDi(iw) + ssxDittt

        endif
#endif
      endif
    enddo
    call timing%stop(timing%m_cor_ra_sx)

    ! JRD: Now do CH term
    call timing%start(timing%m_cor_ra_ch)
    SAFE_ALLOCATE(schDt_array,(sig%Nfreq))
    schDt_array(:) = 0D0
!$OMP PARALLEL do private (my_igp,igp,indigp,ig,schDtt,I_epsRggp_int, &
!$OMP                      I_epsAggp_int,schD,schDt)
    do ifreq=1,sig%Nfreq
        schDt = (0D0,0D0)

        do my_igp = 1, ngpown
          indigp = inv_igp_index(my_igp)
          igp = indinv(indigp)

          if (igp .gt. ncouls .or. igp .le. 0) cycle

! JRD: The below loop is performance critical

          schDtt = (0D0,0D0)
          do ig = 1, ncouls

            I_epsRggp_int = I_epsR_array(ig,my_igp,ifreq)

#ifdef CPLX
            I_epsAggp_int = I_epsA_array(ig,my_igp,ifreq)

            ! for G,G` components
            schD=I_epsRggp_int-I_epsAggp_int

            ! for G`,G components
            schDtt = schDtt + aqsntemp(ig,n1)*schD
#else
            schD= CMPLX(IMAG(I_epsRggp_int),0.0d0)
            schDtt = schDtt + aqsntemp(ig,n1)*schD
#endif
          enddo
          schDt = schDt + schDtt*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
        enddo

! XXX: Threads could be stomping on each-other`s cache over this... try reduction?
        schDt_array(ifreq) = schDt

    enddo
!$OMP END PARALLEL DO
    call timing%stop(timing%m_cor_ra_ch)

    call timing%start(timing%m_cor_ra_sum)
!$OMP PARALLEL do private (ifreq,schDt,cedifft_zb,cedifft_coh,cedifft_cor, &
!$OMP                      cedifft_zb_right,cedifft_zb_left,schDt_right,schDt_left, &
!$OMP                      schDt_avg,schDt_lin,schDt_lin2,intfact,iw, &
!$OMP                      schDt_lin3) reduction(+:schDi,schDi_corb,schDi_cor,sch2Di)
    do ifreq=1,sig%Nfreq

        schDt = schDt_array(ifreq)

        cedifft_zb = sig%dFreqGrid(ifreq)
        cedifft_coh = CMPLX(cedifft_zb,0D0)- sig%dFreqBrd(ifreq)

        if( flag_occ) then
          cedifft_cor = -1.0d0*CMPLX(cedifft_zb,0D0) - sig%dFreqBrd(ifreq)
        else
          cedifft_cor = CMPLX(cedifft_zb,0D0) - sig%dFreqBrd(ifreq)
        endif

        if (ifreq .ne. 1) then
          cedifft_zb_right = cedifft_zb
          cedifft_zb_left = sig%dFreqGrid(ifreq-1)
          schDt_right = schDt
          schDt_left = schDt_array(ifreq-1)
          schDt_avg = 0.5D0 * ( schDt_right + schDt_left )
          schDt_lin = schDt_right - schDt_left
          schDt_lin2 = schDt_lin/(cedifft_zb_right-cedifft_zb_left)
        endif

#ifdef CPLX
! The below two lines are for sigma1 and sigma3
        if (ifreq .ne. sig%Nfreq) then
          schDi(:) = schDi(:) - CMPLX(0.d0,pref(ifreq)) * schDt / ( wxi(:)-cedifft_coh)
          schDi_corb(:) = schDi_corb(:) - CMPLX(0.d0,pref(ifreq)) * schDt / ( wxi(:)-cedifft_cor)
        endif
        if (ifreq .ne. 1) then
          do iw = 1, sig%nfreqeval
!These lines are for sigma2
            intfact=abs((wxi(iw)-cedifft_zb_right)/(wxi(iw)-cedifft_zb_left))
            if (intfact .lt. 1d-4) intfact = 1d-4
            if (intfact .gt. 1d4) intfact = 1d4
            intfact = -log(intfact)
            sch2Di(iw) = sch2Di(iw) - CMPLX(0.d0,pref_zb) * schDt_avg * intfact
!These lines are for sigma4
            if (flag_occ) then
              intfact=abs((wxi(iw)+cedifft_zb_right)/(wxi(iw)+cedifft_zb_left))
              if (intfact .lt. 1d-4) intfact = 1d-4
              if (intfact .gt. 1d4) intfact = 1d4
              intfact = log(intfact)
              schDt_lin3 = (schDt_left + schDt_lin2*(-wxi(iw)-cedifft_zb_left))*intfact
            else
              schDt_lin3 = (schDt_left + schDt_lin2*(wxi(iw)-cedifft_zb_left))*intfact
            endif
            schDt_lin3 = schDt_lin3 + schDt_lin
            schDi_cor(iw) = schDi_cor(iw) - CMPLX(0.d0,pref_zb) * schDt_lin3
          enddo
        endif
#else
! The below two lines are for sigma1 and sigma3
        if (ifreq .ne. sig%Nfreq) then
          schDi(:) = schDi(:) + pref(ifreq)*schDt / ( wxi(:)-cedifft_coh)
          schDi_corb(:) = schDi_corb(:) + pref(ifreq)*schDt / ( wxi(:)-cedifft_cor)
        endif
        if (ifreq .ne. 1) then
          do iw = 1, sig%nfreqeval
!These lines are for sigma2
            intfact=abs((wxi(iw)-cedifft_zb_right)/(wxi(iw)-cedifft_zb_left))
            if (intfact .lt. 1d-4) intfact = 1d-4
            if (intfact .gt. 1d4) intfact = 1d4
            intfact = -log(intfact)
            sch2Di(iw) = sch2Di(iw) + pref_zb * schDt_avg * intfact
!These lines are for sigma4
            if (flag_occ) then
              intfact=abs((wxi(iw)+cedifft_zb_right)/(wxi(iw)+cedifft_zb_left))
              if (intfact .lt. 1d-4) intfact = 1d-4
              if (intfact .gt. 1d4) intfact = 1d4
              intfact = log(intfact)
              schDt_lin3 = (schDt_left + schDt_lin2*(-wxi(iw)-cedifft_zb_left))*intfact
            else
              schDt_lin3 = (schDt_left + schDt_lin2*(wxi(iw)-cedifft_zb_left))*intfact
            endif
            schDt_lin3 = schDt_lin3 + schDt_lin
            schDi_cor(iw) = schDi_cor(iw) + pref_zb * schDt_lin3
          enddo
        endif
#endif
    enddo
!$OMP END PARALLEL DO
    SAFE_DEALLOCATE(schDt_array)
    call timing%stop(timing%m_cor_ra_sum)

    ! JRD: Compute Sigma2 and Sigma4 delta function contributions
    call timing%start(timing%m_cor_ra_ch2)
    do iw = 1, sig%nfreqeval
      wx = wxi(iw)
      if(wx .ge. 0.0d0) then
        ifreq=0
        do ijk = 1, sig%nFreq-1
          if (wx .ge. sig%dFreqGrid(ijk) .and. wx .lt. sig%dFreqGrid(ijk+1)) then
            ifreq=ijk
          endif
        enddo
        if (ifreq .eq. 0) then
          ifreq=sig%nFreq-1
        endif

        fact1=(sig%dFreqGrid(ifreq+1)-wx)/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))
        fact2=(wx-sig%dFreqGrid(ifreq))/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))

        schDttt = 0D0
        schDttt_cor = 0D0

!$OMP PARALLEL do private (my_igp,igp,indigp,ig, &
!$OMP                      sch2Dt,sch2Dtt) reduction(+:schDttt,schDttt_cor)
        do my_igp = 1, ngpown
          indigp = inv_igp_index(my_igp)
          igp = indinv(indigp)

          if (igp .gt. ncouls .or. igp .le. 0) cycle

          sch2Dtt = (0D0,0D0)
          do ig = 1, ncouls
#ifdef CPLX
            sch2Dt=-0.5D0*((I_epsR_array(ig,my_igp,ifreq)-I_epsA_array(ig,my_igp,ifreq))*fact1 + &
                   (I_epsR_array(ig,my_igp,ifreq+1)-I_epsA_array(ig,my_igp,ifreq+1))*fact2)
#else
            sch2Dt = CMPLX(0D0,-1D0)* &
              (IMAG(I_epsR_array(ig,my_igp,ifreq))*fact1 + IMAG(I_epsR_array(ig,my_igp,ifreq+1))*fact2)
#endif
            sch2Dtt = sch2Dtt + aqsntemp(ig,n1)*sch2Dt
          enddo
          schDttt = schDttt + sch2Dtt*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
          if (.not.flag_occ) then
            schDttt_cor = schDttt_cor + sch2Dtt*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
          endif
        enddo
!$OMP END PARALLEL DO

        sch2Di(iw) = sch2Di(iw) + schDttt
        schDi_cor(iw) = schDi_cor(iw) + schDttt_cor
      else if (flag_occ) then
        wx=-wx
        ifreq=0
        do ijk = 1, sig%nFreq-1
          if (wx .ge. sig%dFreqGrid(ijk) .and. wx .lt. sig%dFreqGrid(ijk+1)) then
            ifreq=ijk
          endif
        enddo
        if (ifreq .eq. 0) then
          ifreq=sig%nFreq-1
        endif

        fact1=(sig%dFreqGrid(ifreq+1)-wx)/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))
        fact2=(wx-sig%dFreqGrid(ifreq))/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))

        schDttt_cor = 0D0

!$OMP PARALLEL do private (my_igp,igp,indigp,ig, &
!$OMP                      sch2Dt,sch2Dtt) reduction(+:schDttt_cor)
        do my_igp = 1, ngpown
          indigp = inv_igp_index(my_igp)
          igp = indinv(indigp)

          if (igp .gt. ncouls .or. igp .le. 0) cycle

          sch2Dtt = (0D0,0D0)
          do ig = 1, ncouls
#ifdef CPLX
            sch2Dt=-0.5D0*((I_epsR_array(ig,my_igp,ifreq)-I_epsA_array(ig,my_igp,ifreq))*fact1 + &
                   (I_epsR_array(ig,my_igp,ifreq+1)-I_epsA_array(ig,my_igp,ifreq+1))*fact2)
#else
            sch2Dt = CMPLX(0D0,-1D0)* &
              (IMAG(I_epsR_array(ig,my_igp,ifreq))*fact1 + IMAG(I_epsR_array(ig,my_igp,ifreq+1))*fact2)
#endif
            sch2Dtt = sch2Dtt + aqsntemp(ig,n1)*sch2Dt
          enddo
          schDttt_cor = schDttt_cor + sch2Dtt*vcoul(igp)*MYCONJG(aqsmtemp(igp,n1))
        enddo
!$OMP END PARALLEL DO
        schDi_cor(iw) = schDi_cor(iw) + schDttt_cor
      endif
    enddo
    call timing%stop(timing%m_cor_ra_ch2)

    do iw = 1, sig%nfreqeval
      if (flag_occ) then
        asxDtemp(iw) = asxDtemp(iw) + ssxDi(iw)*occ
      endif
      achDtemp(iw) = achDtemp(iw) + schDi(iw)
      achDtemp_cor(iw) = achDtemp_cor(iw) + schDi_cor(iw)
      achDtemp_corb(iw) = achDtemp_corb(iw) + schDi_corb(iw)
      ach2Dtemp(iw) = ach2Dtemp(iw) + sch2Di(iw)
    ! JRD: This is now close to LDA
      if (iw.eq.iwlda) achtD_n1(n1true) = &
        achtD_n1(n1true) + schDi(iw)
    enddo ! over iw
777         format(1x,"WARNING: The real frequency range is too small." &
              ,/,3x,"l =",i3,1x,"n1 =",i5,1x,"E_l =",f8.3,1x,"E_n1" &
              ,1x,"=",f8.3,1x,"wx =",f8.3,1x,"E_max =",f8.3)

    POP_SUB(mtxel_cor.sigma_ra)

  end subroutine sigma_ra


!==============================================================================
!                             CONTOUR DEFORMATION
!==============================================================================
  !> Calculate Sigma matrix elements, full-freq / CONTOUR-DEFORMATION formalism
  subroutine sigma_cd()
    real(DP) :: imag_freqs(sig%nfreq_imag)
    SCALAR :: c0, c1, c2, c3, y1, y2, y3, m1, m2, mm(sig%nfreq_imag)
    real(DP) :: x1, x2, x3, dx
    !XXX SCALAR :: sW_imag_freqs(sig%nfreq_imag), sint(sig%nfreqeval)
    SCALAR, allocatable :: sW_imag_freqs(:), sint(:)
    SCALAR :: sW_iomega, sW_iomega_acc
    complex(DPC) :: sres(sig%nfreqeval), sres_omega, sres_omega_acc
    integer :: occ_sign, ijk, ifreq_my_igp
    real(DP) :: fact1, fact2
    logical, save :: warned=.false.

    ! FHJ: WARNING: we calc. the TO Sigma here, while RA/FF calculates Ret Sigma.
    ! Everything we implement here is Eqn. (14) from PRB 67, 155208 (2003).
    ! Note the apparent sign flip; that`s because we have (1-epsinv)*v instead
    ! of (epsinv-1)*v = W^c.
    !
    ! FHJ: Definitions (following PRB 67, 155208 notation):
    ! sW_imag_freqs(i \omega) = - \sum_{G, G`} M^*_{G} W^c_{G G`}(i \omega) M(G`)
    ! sint(\omega) = \frac{1}{\pi} \int_0^{\inf} d{\omega`}
    !   sW_imag_freqs(i \omega`) (\omega-Elk)/((\omega-Elk)^2 + (\omega`)^2)
    !
    ! Note: we actually evaluate sint(\omega) using the quadrature scheme
    ! from Fabien`s thesis.

    PUSH_SUB(mtxel_cor.sigma_cd)

    !XXX This allocation seems to be necessary with PGI compiler for the OMP reduction with -O flag X>1
    SAFE_ALLOCATE(sW_imag_freqs, ( sig%nfreq_imag ) )
    SAFE_ALLOCATE(sint, (sig%nfreqeval) )

    ! JRD: Find iw closest to e_lk
    diffmin = INF
    do iw=1,sig%nfreqeval
      diff = abs(freq0 + (iw-1)*sig%freqevalstep - e_lk)
      if (diff .lt. diffmin) then
        diffmin=diff
        iwlda=iw
      endif
    enddo

    ! wxi = omega - e_n``k
    do iw=1,sig%nfreqeval
      wx = freq0 - e_n1kq + (iw-1)*sig%freqevalstep
      wxi(iw) = wx
    enddo

    ! GSM: compute the static CH for the static remainder
    if (sig%exact_ch.eq.1) then
      call acc_static_ch(ngpown, ncouls, inv_igp_index, indinv, vcoul, &
        aqsntemp(:,n1), aqsmtemp(:,n1), achstemp, eps_cplx=I_epsR_array(:,:,1))
    endif

    sres(:) = (0D0,0D0) ! Residue contrib.
    sint(:) = ZERO ! Integral contrib.

    ! FHJ: residue contribution to the correlation self energy
    ! JRD: Don`t want to thread here, nfreqeval could be small
    call timing%start(timing%m_cor_cd_res)
    do iw=1,sig%nfreqeval

      wx = wxi(iw)
      ! FHJ: need either wx>=0 and empty or wx<0 and occupied
      if ( (wx>=0.0d0) .eqv. flag_occ) cycle
      occ_sign = idint(sign(1d0, wx))
      ! FHJ: and from now on we can work on | omega - e_n``k |
      wx = abs(wx)

      ifreq = 0
      do ijk = 1, nfreq_real-1
        if (wx>=sig%dFreqGrid(ijk) .and. wx<sig%dFreqGrid(ijk+1)) then
          ifreq = ijk
        endif
      enddo
      if (ifreq==0) then
        ifreq = nfreq_real+3 ! three is for luck
      endif
      if(ifreq>=nfreq_real) then
        if (.not.warned .and. peinf%pool_rank==0) then
          write(0,777) iband, n1true, e_lk, e_n1kq, wx, E_max
        endif
        warned = .true.
        ifreq = nfreq_real-1
      endif

      sres_omega = 0D0
      if (nfreq_real>1) then
        fact1 = (sig%dFreqGrid(ifreq+1)-wx)/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))
        fact2 = (wx-sig%dFreqGrid(ifreq))/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))
      endif

!$OMP PARALLEL do private (my_igp,igp,indigp,sres_omega_acc,ig) reduction(+:sres_omega)
      do my_igp = 1, ngpown
        indigp = inv_igp_index(my_igp)
        igp = indinv(indigp)
        if (igp>ncouls .or. igp<=0) cycle

        sres_omega_acc = (0D0,0D0)
        ! JRD: The below loop is performance critical
        if (nfreq_real>1) then
          do ig = 1, ncouls
            sres_omega_acc = sres_omega_acc + aqsntemp(ig,n1)*( &
              I_epsR_array(ig,my_igp,ifreq)*fact1 + &
              I_epsR_array(ig,my_igp,ifreq+1)*fact2)
          enddo
        else
          do ig = 1, ncouls
            sres_omega_acc = sres_omega_acc + aqsntemp(ig,n1)*I_epsR_array(ig,my_igp,1)
          enddo
        endif
        sres_omega = sres_omega + MYCONJG(aqsmtemp(igp,n1))*sres_omega_acc*vcoul(igp)
      enddo
!$OMP END PARALLEL DO

      sres(iw) = sres(iw) - occ_sign*sres_omega

    enddo ! iw loop
    call timing%stop(timing%m_cor_cd_res)

    ! JRD: Now do Integral contribution
    call timing%start(timing%m_cor_cd_int)
    sW_imag_freqs(:) = ZERO
!$OMP PARALLEL DO PRIVATE(igp,indigp,ig,sW_iomega_acc,ifreq,my_igp,ifreq_my_igp) &
!$OMP REDUCTION(+:sW_imag_freqs)
    do ifreq_my_igp = 0, sig%nfreq_imag*ngpown - 1
      ifreq = ifreq_my_igp/ngpown + 1
      my_igp = mod(ifreq_my_igp, ngpown) + 1
      indigp = inv_igp_index(my_igp)
      igp = indinv(indigp)
      if (igp>ncouls .or. igp<=0) cycle

      sW_iomega_acc = ZERO
      ! JRD: The below loop is performance critical
      do ig = 1, ncouls
        sW_iomega_acc = sW_iomega_acc + aqsntemp(ig,n1)*I_eps_imag(ig,my_igp,ifreq)
      enddo
      sW_imag_freqs(ifreq) = sW_imag_freqs(ifreq) + &
        MYCONJG(aqsmtemp(igp,n1))*sW_iomega_acc*vcoul(igp)
    enddo
!$OMP END PARALLEL DO
    call timing%stop(timing%m_cor_cd_int)

    call timing%start(timing%m_cor_cd_sum)
    imag_freqs(:) = IMAG(sig%dFreqGrid(nfreq_real+1:)+sig%dFreqBrd(nfreq_real+1:))
    select case (sig%cd_int_method)

      case (0)
        ! FHJ: Integration from Fabien`s thesis; assumes that the matrix elements
        ! are piecewise constant functions, and perform the integral over each
        ! interval centered are the frequency evaluation point.

        ! JRD: XXX Should center this integral....
        ! FHJ: What do you mean?
        ifreq = 1
        sW_iomega = sW_imag_freqs(ifreq)
        freqStart = imag_freqs(ifreq)
        freqEnd = (imag_freqs(ifreq) + imag_Freqs(ifreq+1)) * 0.5d0
        sint(:) = sint(:) + sW_iomega * ( atan(freqEnd/wxi(:)) - atan(freqStart/wxi(:)) )

        !FHJ: no need to do OMP here, each trip is very short, and OMP can
        !hurt vectorization.
        do ifreq = 2, sig%nfreq_imag-1
            sW_iomega = sW_imag_freqs(ifreq)
            freqStart = (imag_freqs(ifreq-1) + imag_Freqs(ifreq)) * 0.5d0
            freqEnd = (imag_freqs(ifreq) + imag_Freqs(ifreq+1)) * 0.5d0
            sint(:) = sint(:) + sW_iomega * ( atan(freqEnd/wxi(:)) - atan(freqStart/wxi(:)) )
        enddo

        sW_iomega = sW_imag_freqs(sig%nfreq_imag)
        freqStart = (imag_freqs(sig%nfreq_imag-1) + imag_Freqs(sig%nfreq_imag)) * 0.5d0
        freqEnd = (-imag_freqs(sig%nfreq_imag-1) + 3d0*imag_Freqs(sig%nfreq_imag)) * 0.5d0
        sint(:) = sint(:) + sW_iomega * ( atan(freqEnd/wxi(:)) - atan(freqStart/wxi(:)) )


      case (2)
        ! FHJ: Integration scheme: assume piecewise quadratic matrix elements
        ! on each integration segment

        ! FHJ: First segment is special. TODO: assert that freqStart==0d0
        freqStart = imag_freqs(1)
        freqEnd = imag_freqs(2)
        c0 = sW_imag_freqs(1)
        c1 = 0d0
        c2 = (sW_imag_freqs(2) - c0) / (freqEnd-freqStart)**2
        sint(:) = sint(:) + &
          c2 * wxi(:) * (freqEnd-freqStart) + &
          (c0 - c2*wxi(:)**2) * (atan(freqEnd/wxi(:)) - atan(freqStart/wxi(:))) + &
          c1*wxi(:)*0.5d0*log((wxi(:)**2 + freqEnd**2)/(wxi(:)**2 + freqStart**2))

        ! FHJ: other segments: find c0, c1 and c2 by using each previous point.
        ! This also works for complex y`s. This sum is over the starting points.
        do ifreq = 2, sig%nfreq_imag-1
          x1 = imag_freqs(ifreq-1)
          x2 = imag_freqs(ifreq)
          x3 = imag_freqs(ifreq+1)
          y1 = sW_imag_freqs(ifreq-1)
          y2 = sW_imag_freqs(ifreq)
          y3 = sW_imag_freqs(ifreq+1)
          c2 = ((y2-y1)*(x1-x3) + (y3-y1)*(x2-x1)) / &
               ((x1-x3)*(x2**2-x1**2) + (x2-x1)*(x3**2-x1**2))
          c1 = ((y2-y1) - c2*(x2**2-x1**2)) / (x2-x1)
          c0 = y1 - c2*x1**2 - c1*x1
          freqEnd = x3
          freqStart = x2
          sint(:) = sint(:) + &
            c2 * wxi(:) * (freqEnd-freqStart) + &
            (c0 - c2*wxi(:)**2) * (atan(freqEnd/wxi(:)) - atan(freqStart/wxi(:))) + &
            c1*wxi(:)*0.5d0*log((wxi(:)**2 + freqEnd**2)/(wxi(:)**2 + freqStart**2))
        enddo

        ! FHJ: add Lorentzian tail? This is only a good idea is the last freq. is really large
        !freqStart = x3
        !AA = y3*x3**2
        !sint(:) = sint(:) + AA/(2d0*wxi(:)**2*freqStart)*( &
        !  2d0*wxi(:) - sign(1d0,wxi(:))*PI_D*freqStart + 2d0*freqStart*atan(freqStart/wxi(:)) )


      case (3)
        ! FHJ: Integration scheme: assume piecewise cubic matrix elements
        ! on each integration segment. Find cubic Hermite spline representation
        ! using a finite difference approximations for the derivatives

        ! Estimate derivatives
        mm(1) = ZERO
        mm(2:sig%nfreq_imag-1) = &
          0.5d0*(sW_imag_freqs(3:sig%nfreq_imag) - sW_imag_freqs(2:sig%nfreq_imag-1)) / &
          ( imag_freqs(3:sig%nfreq_imag) - imag_freqs(2:sig%nfreq_imag-1) ) + &
          0.5d0*(sW_imag_freqs(2:sig%nfreq_imag-1) - sW_imag_freqs(1:sig%nfreq_imag-2)) / &
          ( imag_freqs(2:sig%nfreq_imag-1) - imag_freqs(1:sig%nfreq_imag-2) )
        mm(sig%nfreq_imag) = mm(sig%nfreq_imag-1) + &
          (imag_freqs(sig%nfreq_imag)-imag_freqs(sig%nfreq_imag-1)) * &
          (mm(sig%nfreq_imag-1)-mm(sig%nfreq_imag-2)) / &
          (imag_freqs(sig%nfreq_imag-1)-imag_freqs(sig%nfreq_imag-2))

        ! FHJ: other segments: find c0, c1 and c2 by using each previous point.
        ! This also works for complex y`s. This sum is over the starting points.
        do ifreq = 1, sig%nfreq_imag-1
          x1 = imag_freqs(ifreq)
          x2 = imag_freqs(ifreq+1)
          y1 = sW_imag_freqs(ifreq)
          y2 = sW_imag_freqs(ifreq+1)
          dx = x2 - x1

          c0 = -(m1*x1) - (2*m1*x1**2)/dx - (m2*x1**2)/dx - (m1*x1**3)/dx**2 - &
                (m2*x1**3)/dx**2 + y1 - (3*x1**2*y1)/dx**2 - (2*x1**3*y1)/dx**3 + &
                (3*x1**2*y2)/dx**2 + (2*x1**3*y2)/dx**3
          c1 =  m1 + (4*m1*x1)/dx + (2*m2*x1)/dx + (3*m1*x1**2)/dx**2 + &
                (3*m2*x1**2)/dx**2 + (6*x1*y1)/dx**2 + (6*x1**2*y1)/dx**3 - &
                (6*x1*y2)/dx**2 - (6*x1**2*y2)/dx**3
          c2 =  (-2*m1)/dx - m2/dx - (3*m1*x1)/dx**2 - (3*m2*x1)/dx**2 - &
                (3*y1)/dx**2 - (6*x1*y1)/dx**3 + (3*y2)/dx**2 + (6*x1*y2)/dx**3
          c3 = m1/dx**2 + m2/dx**2 + (2*y1)/dx**3 - (2*y2)/dx**3

          sint(:) = sint(:) + ((wxi*dx)*(2*c2 + c3*(x1+x2)) - &
            2d0*(c0 - wxi**2*c2)*(atan(x1/wxi) - atan(x2/wxi)) + &
            wxi*(-c1 + wxi**2*c3)*log((wxi**2 + x1**2)/(wxi**2 + x2**2)))/2d0
        enddo

      case default
        call die('Invalid integration method')

    endselect !sig%cd_int_method
    call timing%stop(timing%m_cor_cd_sum)

    do iw = 1, sig%nfreqeval
      !FHJ: the output accumulated arrays "asxtDyn" and "asxtDyn" actually store
      !the residue and integral contribution to the GW self energy for CD
      !calculations, respectively.
      asxDtemp(iw) = asxDtemp(iw) + sres(iw)
      achDtemp(iw) = achDtemp(iw) + sint(iw)/PI_D
      ! JRD: This is now close to LDA
      if (iw==iwlda) achtD_n1(n1true) = achtD_n1(n1true) + sint(iw)/PI_D
    enddo ! over iw

777         format(1x,"WARNING: The real frequency range is too small." &
              ,/,3x,"l =",i3,1x,"n1 =",i5,1x,"E_l =",f8.3,1x,"E_n1" &
              ,1x,"=",f8.3,1x,"wx =",f8.3,1x,"E_max =",f8.3)

    SAFE_DEALLOCATE( sW_imag_freqs )
    SAFE_DEALLOCATE( sint )

    POP_SUB(mtxel_cor.sigma_cd)

  end subroutine sigma_cd

!==============================================================================
!                        SUBSPACE CONTOUR DEFORMATION
!==============================================================================
  !> Calculate Subspace Sigma matrix elements, full-freq / CONTOUR-DEFORMATION formalism
  subroutine sigma_cd_subspace()
    real(DP) :: imag_freqs(sig%nfreq_imag)
    SCALAR :: c0, c1, c2, c3, y1, y2, y3, m1, m2, mm(sig%nfreq_imag)
    real(DP) :: x1, x2, x3, dx
    !XXX SCALAR :: sW_imag_freqs(sig%nfreq_imag, peinf%ntband_max), sint(sig%nfreqeval)
    SCALAR, allocatable :: sW_imag_freqs(:,:), sint(:)
    SCALAR :: sW_iomega, sW_iomega_acc
    complex(DPC) :: sres(sig%nfreqeval), sres_omega, sres_omega_acc
    integer :: occ_sign, ijk, ifreq_my_igp
    real(DP) :: fact1, fact2
    logical, save :: warned=.false.

    integer :: n1
    SCALAR, allocatable :: tempz(:,:,:), achstemp_array(:)
    ! variables for static reminder
    SCALAR, ALLOCATABLE :: schs_array(:,:), eps_scalar(:,:), achs_n1(:)
    logical :: use_zgemm
    integer :: ipe_block, ipe_real, my_ib
    integer :: n1_start, n1_end, n1_size_actual
    integer :: n1_ipe

    PUSH_SUB(mtxel_cor.sigma_cd_subspace)

    !XXX This allocation seems to be necessary with PGI compiler for the OMP reduction with -O flag X>1
    SAFE_ALLOCATE( sW_imag_freqs, (sig%nfreq_imag, peinf%ntband_max) )
    SAFE_ALLOCATE( sint, (sig%nfreqeval) )

    ! check if we need to calculate the static reminder
    if (sig%exact_ch.eq.1) then
      call timing%start(timing%m_cor_remain)
      SAFE_ALLOCATE(schs_array, (my_Nb, peinf%ntband_dist(ipe)))
      SAFE_ALLOCATE(achs_n1, (peinf%ntband_dist(ipe)))

!$OMP PARALLEL do private (n1)
      do n1 = 1, peinf%ntband_dist(ipe)
        achs_n1(n1) = ZERO
      end do
!$OMP END PARALLEL DO

      if (my_Nb_start > 0) then
#ifdef CPLX
        call zgemm('t','n', my_Nb,  peinf%ntband_dist(ipe), Nb_tot, (-1D0,0D0), &
                   Re_eps_sub(:,:,1), Nb_tot, &
                   Msub_n_temp, Nb_tot, (0D0,0D0), &
                   schs_array, my_Nb)
#else
        SAFE_ALLOCATE(eps_scalar, (Nb_tot,my_Nb))
        eps_scalar(:,:) = SCALARIFY(Re_eps_sub(:,:,1))
        CALL dgemm('t','n',  my_Nb,  peinf%ntband_dist(ipe), Nb_tot, -1.0D+00, &
                   eps_scalar(:,:), Nb_tot, &
                   Msub_n_temp, Nb_tot, 0.0D+00, &
                   schs_array, my_Nb)
        SAFE_DEALLOCATE(eps_scalar)
#endif

!$OMP PARALLEL do private (n1,my_ib) reduction(+:achs_n1)
        do n1 = 1, peinf%ntband_dist(ipe)
          do my_ib = 1, my_Nb
            achs_n1(n1) = achs_n1(n1) + MYCONJG(Msub_m_temp(my_ib+my_Nb_start-1,n1)) * &
                                        schs_array(my_ib,n1) * 0.5D0
          enddo
        end do
!$OMP END PARALLEL DO
      end if ! my_Nb_start > 0

      ! correct wings, same flag as residual part (computed in the same loop later)
      ! note the sign is changed and the famous factor 0.5
      if(ipe-1 == peinf%pool_rank) then
        if(fix_wings_res) then
          do n1 = 1, peinf%ntband_dist(ipe)
            achs_n1(n1) = achs_n1(n1) - MYCONJG(aqsmtemp(wing_pos,n1)) * &
                                        SUM(wings_Re(1:ncouls,1,2) * aqsntemp(1:ncouls,n1)) * &
                                            vcoul(indinv(wing_pos)) * 0.5d0
            do ig = 1, ncouls
              achs_n1(n1) = achs_n1(n1) - MYCONJG(aqsmtemp(ig,n1)) * &
                                              wings_Re(ig,1,1) * &
                                              aqsntemp(wing_pos,n1) * &
                                              vcoul(indinv(ig)) * 0.5d0
            end do                                    
          end do
        end if
      end if ! wings

      do n1 = 1, peinf%ntband_dist(ipe)
        n1true = peinf%indext_dist(n1,ipe)
        achtcor_n1(n1true) = achtcor_n1(n1true) - 0.5d0 * achs_n1(n1)
      end do

      SAFE_DEALLOCATE(schs_array)
      SAFE_DEALLOCATE(achs_n1)

      call timing%stop(timing%m_cor_remain)
    end if

    ! residual contribution
    do n1 = 1, peinf%ntband_dist(ipe)

      ! n1true = "True" band index of the band n1 w.r.t. all bands
      n1true = peinf%indext_dist(n1,ipe)
      ! energy of the |n1,k-q> state
      e_n1kq = wfnkq%ekq(n1true,ispin)
      ! occupation of the |n1,k-q> state
      flag_occ = (n1true<=(sig%nvband+sig%ncrit)) &
        .and.((sig%ncrit==0).or.(e_n1kq<=(sig%efermi+sig%tol)))
      if (abs(e_n1kq-sig%efermi)<sig%tol) then
        occ = 0.5d0 ! Fermi-Dirac distribution = 1/2 at Fermi level
      else
        occ = 1.0d0
      endif

      ! JRD: Find iw closest to e_lk
      diffmin = INF
      do iw=1,sig%nfreqeval
        diff = abs(freq0 + (iw-1)*sig%freqevalstep - e_lk)
        if (diff .lt. diffmin) then
          diffmin=diff
          iwlda=iw
        endif
      enddo

      ! wxi = omega - e_n``k
      do iw=1,sig%nfreqeval
        wx = freq0 - e_n1kq + (iw-1)*sig%freqevalstep
        wxi(iw) = wx
      enddo

      ! FHJ: residue contribution to the correlation self energy
      ! JRD: Don`t want to thread here, nfreqeval could be small
      call timing%start(timing%m_cor_cd_res)
      sres(:) = (0D0,0D0) ! Residue contrib.
      do iw=1,sig%nfreqeval

        wx = wxi(iw)
        ! FHJ: need either wx>=0 and empty or wx<0 and occupied
        if ( (wx>=0.0d0) .eqv. flag_occ) cycle
        occ_sign = idint(sign(1d0, wx))
        ! FHJ: and from now on we can work on | omega - e_n``k |
        wx = abs(wx)

        ifreq = 0
        do ijk = 1, nfreq_real-1
          if (wx>=sig%dFreqGrid(ijk) .and. wx<sig%dFreqGrid(ijk+1)) then
            ifreq = ijk
          endif
        enddo
        if (ifreq==0) then
          ifreq = nfreq_real+3 ! three is for luck
        endif
        if(ifreq>=nfreq_real) then
          if (.not.warned .and. peinf%pool_rank==0) then
            write(0,777) iband, n1true, e_lk, e_n1kq, wx, E_max
          endif
          warned = .true.
          ifreq = nfreq_real-1
        endif

        sres_omega = 0D0
        if (nfreq_real>1) then
          fact1 = (sig%dFreqGrid(ifreq+1)-wx)/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))
          fact2 = (wx-sig%dFreqGrid(ifreq))/(sig%dFreqGrid(ifreq+1)-sig%dFreqGrid(ifreq))
        endif

        if (my_Nb_start > 0) then

!$OMP PARALLEL do private (my_ib, sres_omega_acc, ig) &
!$OMP          reduction(+:sres_omega)
          do my_ib = 1, my_Nb
            sres_omega_acc = (0D0,0D0)

            ! JRD: The below loop is performance critical
            if (nfreq_real>1) then
              do ig = 1, Nb_tot
                sres_omega_acc = sres_omega_acc + Msub_n_temp(ig,n1)*(&
                     Re_eps_sub(ig,my_ib,ifreq)*fact1 + &
                     Re_eps_sub(ig,my_ib,ifreq+1)*fact2)
              end do
            else
              do ig = 1, Nb_tot
                sres_omega_acc = sres_omega_acc + Msub_n_temp(ig,n1)*Re_eps_sub(ig,my_ib,1)
              end do
            end if

            sres_omega = sres_omega + MYCONJG(Msub_m_temp(my_ib+my_Nb_start-1,n1))*sres_omega_acc
          end do
!$OMP END PARALLEL DO

        end if ! my_Nb_start > 0

        if(ipe-1 == peinf%pool_rank) then
          if(fix_wings_res) then
            !XXXX if (nfreq_real>1) then
            !XXXX   sres_omega = sres_omega + MYCONJG(aqsmtemp(wing_pos,n1)) * &
            !XXXX                             SUM((wings_Re(1:ncouls,ifreq,2)*fact1 + &
            !XXXX                                  wings_Re(1:ncouls,ifreq+1,2)*fact2) * aqsntemp(1:ncouls,n1)) * &
            !XXXX                             vcoul(indinv(wing_pos))
            !XXXX   do ig = 1, ncouls
            !XXXX     sres_omega = sres_omega + MYCONJG(aqsmtemp(ig,n1)) * &
            !XXXX                               (wings_Re(ig,ifreq,1)*fact1 + & 
            !XXXX                                wings_Re(ig,ifreq+1,1)*fact2) * aqsntemp(wing_pos,n1) * &
            !XXXX                               vcoul(indinv(ig))
            !XXXX   end do
            !XXXX else
            !XXXX   sres_omega = sres_omega + MYCONJG(aqsmtemp(wing_pos,n1)) * &
            !XXXX                             SUM(wings_Re(1:ncouls,1,2) * aqsntemp(1:ncouls,n1)) * &
            !XXXX                             vcoul(indinv(wing_pos))
            !XXXX   do ig = 1, ncouls
            !XXXX     sres_omega = sres_omega + MYCONJG(aqsmtemp(ig,n1)) * &
            !XXXX                               wings_Re(ig,1,1) * &
            !XXXX                               aqsntemp(wing_pos,n1) * &
            !XXXX                               vcoul(indinv(ig))
            !XXXX   end do
            !XXXX end if
            sres_omega_acc = (0D0,0D0)
            if (nfreq_real>1) then
!$OMP PARALLEL do private (ig) &
!$OMP          reduction(+:sres_omega_acc)
               do ig = 1, ncouls
                 sres_omega_acc = sres_omega_acc + MYCONJG(aqsmtemp(wing_pos,n1)) * &
                                           (wings_Re(ig,ifreq,2)*fact1 + &
                                            wings_Re(ig,ifreq+1,2)*fact2) * aqsntemp(ig,n1) * &
                                            vcoul(indinv(wing_pos))
                 sres_omega_acc = sres_omega_acc + MYCONJG(aqsmtemp(ig,n1)) * &
                                           (wings_Re(ig,ifreq,1)*fact1 + & 
                                            wings_Re(ig,ifreq+1,1)*fact2) * aqsntemp(wing_pos,n1) * &
                                            vcoul(indinv(ig))
               end do
!$END OMP PARALLEL DO
            else
!$OMP PARALLEL do private (ig) &
!$OMP          reduction(+:sres_omega_acc)
               do ig = 1, ncouls
                 sres_omega_acc = sres_omega_acc + MYCONJG(aqsmtemp(wing_pos,n1)) * &
                                            wings_Re(ig,1,2) * aqsntemp(ig,n1) * &
                                            vcoul(indinv(wing_pos))
                 sres_omega_acc = sres_omega_acc + MYCONJG(aqsmtemp(ig,n1)) * &
                                            wings_Re(ig,1,1) * aqsntemp(wing_pos,n1) * &
                                            vcoul(indinv(ig))
               end do
!$END OMP PARALLEL DO
            end if
            sres_omega = sres_omega + sres_omega_acc
          end if
        end if ! fix tail-wings

        sres(iw) = sres(iw) - occ_sign*sres_omega

      enddo ! iw loop
      call timing%stop(timing%m_cor_cd_res)

      do iw = 1, sig%nfreqeval
        !FHJ: the output accumulated arrays "asxtDyn" and "asxtDyn" actually store
        !the residue and integral contribution to the GW self energy for CD
        !calculations, respectively.
        asxDtemp(iw) = asxDtemp(iw) + sres(iw)
      enddo ! over iw

    end do ! n1 RESIDUAL PART
    if(ipe-1 == peinf%pool_rank) fix_wings_res = .false.

    ! integral part
!$OMP PARALLEL do private (n1,ifreq)
      do n1 = 1, peinf%ntband_max
        do ifreq = 1, sig%nfreq_imag
          sW_imag_freqs(ifreq,n1) = ZERO
        end do
      end do
!$OMP END PARALLEL DO

    call timing%start(timing%m_cor_cd_int)

    if (my_Nb_start > 0) then
      ! go with ZGEMM
      SAFE_ALLOCATE(tempz, (my_Nb, sig%nfreq_imag, peinf%ntband_dist(ipe)))
      !
      call timing%start(timing%m_cor_cd_gemm)
#ifdef CPLX
      call zgemm('t','n', my_Nb*sig%nfreq_imag,  peinf%ntband_dist(ipe), Nb_tot, &
                 (1D0,0D0), Im_eps_sub, Nb_tot, &
                            Msub_n_temp, Nb_tot,&
                 (0D0,0D0), tempz, my_Nb*sig%nfreq_imag)
#else
      call dgemm('t','n', my_Nb*sig%nfreq_imag,  peinf%ntband_dist(ipe), Nb_tot, &
                 1.0D0, Im_eps_sub, Nb_tot, &
                        Msub_n_temp, Nb_tot,&
                 0.0D0, tempz, my_Nb*sig%nfreq_imag)

#endif
     call timing%stop(timing%m_cor_cd_gemm)
  
!$OMP PARALLEL do private (n1,my_ib,my_igp,igp,indigp,sW_iomega_acc,ifreq) &
!$OMP          reduction(+:sW_imag_freqs)
      do n1 = 1, peinf%ntband_dist(ipe)
        do ifreq = 1, sig%nfreq_imag
          sW_iomega_acc = ZERO
          do my_ib = 1, my_Nb
            sW_iomega_acc = sW_iomega_acc + &
                MYCONJG(Msub_m_temp(my_Nb_start+my_ib-1,n1))*tempz(my_ib,ifreq,n1)
          end do
          sW_imag_freqs(ifreq,n1) = sW_imag_freqs(ifreq,n1) + sW_iomega_acc
        end do
      end do
!$OMP END PARALLEL DO

      SAFE_DEALLOCATE(tempz)
    end if ! my_Nb_start > 0
  
    ! fix wings
    if (ipe-1 == peinf%pool_rank) then
      if(fix_wings) then
!$OMP PARALLEL do private (n1,ifreq)  collapse(2)
        do n1 = 1, peinf%ntband_dist(ipe)
          do ifreq = 1, sig%nfreq_imag
            sW_imag_freqs(ifreq,n1) = sW_imag_freqs(ifreq,n1) + &
                                      MYCONJG(aqsmtemp(wing_pos,n1)) * &
                                      SUM(wings_Im(1:ncouls,ifreq,2) * aqsntemp(1:ncouls,n1)) * &
                                      vcoul(indinv(wing_pos))
          end do 
        end do
!$OMP END PARALLEL DO  

!$OMP PARALLEL do private (n1,ifreq,ig)  collapse(2)
        do n1 = 1, peinf%ntband_dist(ipe)
          do ifreq = 1, sig%nfreq_imag
            do ig = 1, ncouls
            sW_imag_freqs(ifreq,n1) = sW_imag_freqs(ifreq,n1) + &
                                      MYCONJG(aqsmtemp(ig,n1)) * &
                                      wings_Im(ig,ifreq,1) * aqsntemp(wing_pos,n1) * &
                                      vcoul(indinv(ig))
            end do
          end do
        end do
!$OMP END PARALLEL DO

        fix_wings = .false.
      end if
    end if
    call timing%stop(timing%m_cor_cd_int)

    call timing%start(timing%m_cor_cd_sum)
    ! calculate integral contribution
    do n1 = 1, peinf%ntband_dist(ipe)
  
      ! n1true = "True" band index of the band n1 w.r.t. all bands
      n1true = peinf%indext_dist(n1,ipe)
      ! energy of the |n1,k-q> state
      e_n1kq = wfnkq%ekq(n1true,ispin)
      ! occupation of the |n1,k-q> state
      flag_occ = (n1true<=(sig%nvband+sig%ncrit)) &
        .and.((sig%ncrit==0).or.(e_n1kq<=(sig%efermi+sig%tol)))
      if (abs(e_n1kq-sig%efermi)<sig%tol) then
        occ = 0.5d0 ! Fermi-Dirac distribution = 1/2 at Fermi level
      else
        occ = 1.0d0
      endif
  
      ! FHJ: CH term used for static remainder for current band "n1true".
      achstemp = ZERO
  
      ! JRD: Find iw closest to e_lk
      diffmin = INF
      do iw=1,sig%nfreqeval
        diff = abs(freq0 + (iw-1)*sig%freqevalstep - e_lk)
        if (diff .lt. diffmin) then
          diffmin=diff
          iwlda=iw
        endif
      enddo
  
      ! wxi = omega - e_n``k
      do iw=1,sig%nfreqeval
        wx = freq0 - e_n1kq + (iw-1)*sig%freqevalstep
        wxi(iw) = wx
      enddo
  
      sint(:) = ZERO ! Integral contrib.
      imag_freqs(:) = IMAG(sig%dFreqGrid(nfreq_real+1:)+sig%dFreqBrd(nfreq_real+1:))
      select case (sig%cd_int_method)
  
        case (0)
          ! FHJ: Integration from Fabien`s thesis; assumes that the matrix elements
          ! are piecewise constant functions, and perform the integral over each
          ! interval centered are the frequency evaluation point.
  
          ! JRD: XXX Should center this integral....
          ! FHJ: What do you mean?
          ! MDB: explicit OMP parallelization over the sigma grid index

!$OMP PARALLEL DO private (iw,ifreq,sW_iomega,freqStart,freqEnd) &
!$OMP             reduction(+:sint)
          do iw=1,sig%nfreqeval
            ifreq = 1
            sW_iomega = sW_imag_freqs(ifreq,n1)
            freqStart = imag_freqs(ifreq)
            freqEnd = (imag_freqs(ifreq) + imag_Freqs(ifreq+1)) * 0.5d0

            !XXX sint(:) = sint(:) + sW_iomega * ( atan(freqEnd/wxi(:)) - atan(freqStart/wxi(:)) )
            sint(iw) = sint(iw) + sW_iomega * ( atan(freqEnd/wxi(iw)) - atan(freqStart/wxi(iw)) )

            !FHJ: no need to do OMP here, each trip is very short, and OMP can
            !hurt vectorization.
            do ifreq = 2, sig%nfreq_imag-1
                sW_iomega = sW_imag_freqs(ifreq,n1)
                freqStart = (imag_freqs(ifreq-1) + imag_Freqs(ifreq)) * 0.5d0
                freqEnd = (imag_freqs(ifreq) + imag_Freqs(ifreq+1)) * 0.5d0

                !XXXX sint(:) = sint(:) + sW_iomega * ( atan(freqEnd/wxi(:)) - atan(freqStart/wxi(:)) )
                sint(iw) = sint(iw) + sW_iomega * ( atan(freqEnd/wxi(iw)) - atan(freqStart/wxi(iw)) )

            enddo
  
            sW_iomega = sW_imag_freqs(sig%nfreq_imag,n1)
            freqStart = (imag_freqs(sig%nfreq_imag-1) + imag_Freqs(sig%nfreq_imag)) * 0.5d0
            freqEnd = (-imag_freqs(sig%nfreq_imag-1) + 3d0*imag_Freqs(sig%nfreq_imag)) * 0.5d0
  
            !XXXX sint(:) = sint(:) + sW_iomega * ( atan(freqEnd/wxi(:)) - atan(freqStart/wxi(:)) )
            sint(iw) = sint(iw) + sW_iomega * ( atan(freqEnd/wxi(iw)) - atan(freqStart/wxi(iw)) )
 
          end do  ! iw
!$OMP END PARALLEL DO

        case (2)
          ! FHJ: Integration scheme: assume piecewise quadratic matrix elements
          ! on each integration segment
  
          ! FHJ: First segment is special. TODO: assert that freqStart==0d0
          freqStart = imag_freqs(1)
          freqEnd = imag_freqs(2)
          c0 = sW_imag_freqs(1,n1)
          c1 = 0d0
          c2 = (sW_imag_freqs(2,n1) - c0) / (freqEnd-freqStart)**2
          sint(:) = sint(:) + &
            c2 * wxi(:) * (freqEnd-freqStart) + &
            (c0 - c2*wxi(:)**2) * (atan(freqEnd/wxi(:)) - atan(freqStart/wxi(:))) + &
            c1*wxi(:)*0.5d0*log((wxi(:)**2 + freqEnd**2)/(wxi(:)**2 + freqStart**2))
  
          ! FHJ: other segments: find c0, c1 and c2 by using each previous point.
          ! This also works for complex y`s. This sum is over the starting points.
          do ifreq = 2, sig%nfreq_imag-1
            x1 = imag_freqs(ifreq-1)
            x2 = imag_freqs(ifreq)
            x3 = imag_freqs(ifreq+1)
            y1 = sW_imag_freqs(ifreq-1,n1)
            y2 = sW_imag_freqs(ifreq,n1)
            y3 = sW_imag_freqs(ifreq+1,n1)
            c2 = ((y2-y1)*(x1-x3) + (y3-y1)*(x2-x1)) / &
                 ((x1-x3)*(x2**2-x1**2) + (x2-x1)*(x3**2-x1**2))
            c1 = ((y2-y1) - c2*(x2**2-x1**2)) / (x2-x1)
            c0 = y1 - c2*x1**2 - c1*x1
            freqEnd = x3
            freqStart = x2
            sint(:) = sint(:) + &
              c2 * wxi(:) * (freqEnd-freqStart) + &
              (c0 - c2*wxi(:)**2) * (atan(freqEnd/wxi(:)) - atan(freqStart/wxi(:))) + &
              c1*wxi(:)*0.5d0*log((wxi(:)**2 + freqEnd**2)/(wxi(:)**2 + freqStart**2))
          enddo
  
          ! FHJ: add Lorentzian tail? This is only a good idea is the last freq. is really large
          !freqStart = x3
          !AA = y3*x3**2
          !sint(:) = sint(:) + AA/(2d0*wxi(:)**2*freqStart)*( &
          !  2d0*wxi(:) - sign(1d0,wxi(:))*PI_D*freqStart + 2d0*freqStart*atan(freqStart/wxi(:)) )
  
  
        case (3)
          ! FHJ: Integration scheme: assume piecewise cubic matrix elements
          ! on each integration segment. Find cubic Hermite spline representation
          ! using a finite difference approximations for the derivatives
  
          ! Estimate derivatives
          mm(1) = ZERO
          mm(2:sig%nfreq_imag-1) = &
            0.5d0*(sW_imag_freqs(3:sig%nfreq_imag,n1) - sW_imag_freqs(2:sig%nfreq_imag-1,n1)) / &
            ( imag_freqs(3:sig%nfreq_imag) - imag_freqs(2:sig%nfreq_imag-1) ) + &
            0.5d0*(sW_imag_freqs(2:sig%nfreq_imag-1,n1) - sW_imag_freqs(1:sig%nfreq_imag-2,n1)) / &
            ( imag_freqs(2:sig%nfreq_imag-1) - imag_freqs(1:sig%nfreq_imag-2) )
          mm(sig%nfreq_imag) = mm(sig%nfreq_imag-1) + &
            (imag_freqs(sig%nfreq_imag)-imag_freqs(sig%nfreq_imag-1)) * &
            (mm(sig%nfreq_imag-1)-mm(sig%nfreq_imag-2)) / &
            (imag_freqs(sig%nfreq_imag-1)-imag_freqs(sig%nfreq_imag-2))
  
          ! FHJ: other segments: find c0, c1 and c2 by using each previous point.
          ! This also works for complex y`s. This sum is over the starting points.
          do ifreq = 1, sig%nfreq_imag-1
            x1 = imag_freqs(ifreq)
            x2 = imag_freqs(ifreq+1)
            y1 = sW_imag_freqs(ifreq,n1)
            y2 = sW_imag_freqs(ifreq+1,n1)
            dx = x2 - x1
  
            c0 = -(m1*x1) - (2*m1*x1**2)/dx - (m2*x1**2)/dx - (m1*x1**3)/dx**2 - &
                  (m2*x1**3)/dx**2 + y1 - (3*x1**2*y1)/dx**2 - (2*x1**3*y1)/dx**3 + &
                  (3*x1**2*y2)/dx**2 + (2*x1**3*y2)/dx**3
            c1 =  m1 + (4*m1*x1)/dx + (2*m2*x1)/dx + (3*m1*x1**2)/dx**2 + &
                  (3*m2*x1**2)/dx**2 + (6*x1*y1)/dx**2 + (6*x1**2*y1)/dx**3 - &
                  (6*x1*y2)/dx**2 - (6*x1**2*y2)/dx**3
            c2 =  (-2*m1)/dx - m2/dx - (3*m1*x1)/dx**2 - (3*m2*x1)/dx**2 - &
                  (3*y1)/dx**2 - (6*x1*y1)/dx**3 + (3*y2)/dx**2 + (6*x1*y2)/dx**3
            c3 = m1/dx**2 + m2/dx**2 + (2*y1)/dx**3 - (2*y2)/dx**3
  
            sint(:) = sint(:) + ((wxi*dx)*(2*c2 + c3*(x1+x2)) - &
              2d0*(c0 - wxi**2*c2)*(atan(x1/wxi) - atan(x2/wxi)) + &
              wxi*(-c1 + wxi**2*c3)*log((wxi**2 + x1**2)/(wxi**2 + x2**2)))/2d0
          enddo
  
        case default
          call die('Invalid integration method')
  
      endselect !sig%cd_int_method
  
      do iw = 1, sig%nfreqeval
        !FHJ: the output accumulated arrays "asxtDyn" and "asxtDyn" actually store
        !the residue and integral contribution to the GW self energy for CD
        !calculations, respectively.
        achDtemp(iw) = achDtemp(iw) + sint(iw)/PI_D
        ! JRD: This is now close to LDA
        if (iw==iwlda) achtD_n1(n1true) = achtD_n1(n1true) + sint(iw)/PI_D
      enddo ! over iw
  
    end do ! n1
    call timing%stop(timing%m_cor_cd_sum)

777         format(1x,"WARNING: The real frequency range is too small." &
              ,/,3x,"l =",i3,1x,"n1 =",i5,1x,"E_l =",f8.3,1x,"E_n1" &
              ,1x,"=",f8.3,1x,"wx =",f8.3,1x,"E_max =",f8.3)

    SAFE_DEALLOCATE( sW_imag_freqs )
    SAFE_DEALLOCATE( sint )

    POP_SUB(mtxel_cor.sigma_cd_subspace)

  end subroutine sigma_cd_subspace

end subroutine mtxel_cor

!> Accumulates the partial static CH as a sum over all bands. We use this to
!! calculate the static remainder later on as 1/2 of the difference between the
!! exact static CH and the partial static CH obtained here.
subroutine acc_static_ch(ngpown, ncouls, inv_igp_index, indinv, vcoul, &
  aqsn_n1, aqsm_n1, achs, eps_scalar, eps_cplx)
  integer, intent(in) :: ngpown, ncouls, inv_igp_index(:), indinv(:)
  real(DP), intent(in) :: vcoul(:)
  SCALAR, intent(in) :: aqsn_n1(:), aqsm_n1(:)
  SCALAR, intent(inout) :: achs
  SCALAR, intent(in), optional :: eps_scalar(:,:)
  complex(DPC), intent(in), optional :: eps_cplx(:,:)

  SCALAR :: schs
  integer :: my_igp, igp, indigp, ig

  PUSH_SUB(acc_static_ch)

  if (present(eps_scalar).eqv.present(eps_cplx)) &
    call die('Internal error calling acc_static_ch', only_root_writes=.true.)

  call timing%start(timing%m_cor_remain)
!$OMP PARALLEL do private (my_igp,indigp,igp,ig,schs) reduction(+:achs)
  do my_igp = 1, ngpown
    indigp = inv_igp_index(my_igp)
    igp = indinv(indigp)
    if (igp>ncouls .or. igp<=0) cycle

    schs = ZERO
    if (present(eps_scalar)) then
      do ig = 1, ncouls
        schs = schs - aqsn_n1(ig)*eps_scalar(ig,my_igp)
      enddo
    else
      do ig = 1, ncouls
        schs = schs - aqsn_n1(ig)*SCALARIFY(eps_cplx(ig,my_igp))
      enddo
    endif
    achs = achs + MYCONJG(aqsm_n1(igp))*schs*vcoul(igp)*0.5D0
  enddo
!$OMP END PARALLEL DO
  call timing%stop(timing%m_cor_remain)

  POP_SUB(acc_static_ch)

end subroutine acc_static_ch

end module mtxel_cor_m

