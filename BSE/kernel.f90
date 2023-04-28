!===================================================================================
!
! Routines:
!
! (1) kernel (main)      Originally By MLT       Last Modified: 5/5/2008 (JRD)
!
!     See README_kernel for more information.
!
!     Calculates the kernel, direct and exchange parts, of the Bethe-Salpeter
!     equation. The direct part is decomposed in head, wing, body, and the
!     exchange part involves only the "proper" part of the Coulomb
!     interaction. Spin-polarized case implemented.
!
!     For more details, see:
!     Rohlfing & Louie, PRB 62:(8), 4927 (2000)
!     G. Strinati, Rivista del Nuovo Cimento 11:(12), 1 (1988)
!
!     Code originally written by Gian-Marco Rignanese, Eric K Chang.
!
!     All cited equations refer to Rohlfing & Louie (PRB (62):4927, 2000)
!     unless specified otherwise.
!
!===================================================================================

#include "f_defs.h"

program kernel

#ifdef HDF5
  use hdf5
#endif
  use global_m
  use accel_memory_m
  use accel_fft_m
  use algos_kernel_m
  use check_screening_m
  use fftw_m
  use fullbz_m
  use io_utils_m
  use mtxel_kernel_m
!#BEGIN_INTERNAL_ONLY
  use mtxel_tddft_m
!#END_INTERNAL_ONLY
  use vcoul_generator_m
  use bsewrite_m
  use epscopy_m
  use genwf_kernel_m
  use input_kernel_m
  use inread_kernel_m
  use references_m
  use sortbyq_m
  use write_program_header_m
  use references_m
  implicit none

  type (crystal) :: crys
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (xctinfo) :: xct
  type (grid) :: kg,qg,kgq
  type (kpoints) :: kp
  type (wavefunction) :: wfnc,wfncp
  type (wavefunction) :: wfnv,wfnvp
  type (int_wavefunction) :: intwfnv,intwfnc
  type (twork_scell) :: work_scell

  character :: filename*20
  integer :: ii,iparallel,iownsize,ijk
  integer :: ik,ikp
  integer :: ncount,ntim,flagbz
  integer :: ic,iv,icp,ivp,error
  integer :: ig
  real(DP) :: qpg_len, gpq0_len, qpg(3)
  integer :: ifqa_dummy,irqa_dummy,g0a_dummy(3)
  real(DP) :: vq0, avgcut, oneoverq, q0len
  real(DP) :: tsec(2),tmin(2),tmax(2),vq(3)
  real(DP) :: qqa_dummy,fqa_dummy(3),vcoul0(1)
  SCALAR :: epsheaddummy, wcoul0dummy

  character*16, allocatable :: routnam(:)
  integer, allocatable :: indexq(:),irqa(:),ifqa(:),g0a(:,:)
  integer, allocatable :: isrtq(:)
  real(DP), allocatable :: vcoul(:),vcoularray(:,:),fqa(:,:),qqa(:),vcoularray_mod(:,:)
  SCALAR, allocatable :: &
    bsedbody(:,:,:),bsedhead(:,:,:), &
    bsedwing(:,:,:),bsex(:,:,:),bset(:,:,:)

  type(progress_info) :: prog_info !< a user-friendly progress report

  ! Buffer sizes for GPU allocation
  integer :: n_bands_max
  real(DP) :: dummy_scale
  integer :: Nfft(3)

!----------------- Begin Program ----------------------------------------------

  call peinfo_init()

#ifdef PGI
  call setvbuf3f(6,2,0)
#endif

!----------------------
! Dont create random numbers

  peinf%jobtypeeval = 0

!----------------------
! Initialize timer

  call timacc(0,0)
  call timacc(1,1)
  
  call write_program_header('BSE/Kernel', .false.)
  
!------------------------
! Read kernel.inp

  call logit('Calling inread_kernel')
  call open_file(8,file='kernel.inp',form='formatted',status='old')
  call inread_kernel(xct,flagbz,qg)
  call close_file(8)

  if (xct%iwritecoul .eq. 1 .and. peinf%inode .eq. 0) then
    call open_file(19,file='vcoul',form='formatted',status='replace')
  endif

!----------------------

#ifdef HDF5
  if(xct%use_hdf5 .or. xct%use_wfn_hdf5) call h5open_f(error)
#endif

!-------------------------
! Read WFN_co and WFNq_co if using finite momentum

  call timacc(2,1)
  call logit('Calling input_kernel')

  call input_kernel(crys,gvec,kg,kgq,kp,syms,xct,flagbz,intwfnv,intwfnc)
  SAFE_ALLOCATE(indexq, (kg%nf))
  indexq(:)=xct%indexq(:)

  if (xct%nspinor == 2) then
    call require_reference(REF_Wu2020)
  endif
  
  call timacc(2,2)

!      if(peinf%inode.eq.0) write(6,*) 'Exit input_kernel'

! JRD: Write some info about our calculation

  if (peinf%inode.eq.0) then
    write(6,*)
    call output_algos()

    write(6,'(/1x,a)') 'Calculation parameters:'
    if (xct%theory.eq.0) then
      write(6,'(1x,a)') '- This is a Bethe-Salpeter-equation calculation.'
    elseif (xct%theory.eq.1) then
      write(6,'(1x,a)') '- This is a time-dependent density-functional-theory calculation.'
    endif
    write(6,'(1x,a,f0.2)') '- Cutoff of the bare Coulomb interaction (Ry): ', xct%ecutg
    if (xct%theory==0) &
      write(6,'(1x,a,f0.2)') '- Cutoff of the screened Coulomb interaction (Ry): ', xct%ecute
    write(6,'(1x,a,i0)') '- Number of G-vectors up to the bare int. cutoff: ', xct%ng
    if (xct%theory==0) &
      write(6,'(1x,a,i0)') '- Number of G-vectors up to the screened int. cutoff: ', xct%neps
    write(6,'(1x,a,i0)') '- Number of valence bands: ', xct%nvb_co
    write(6,'(1x,a,i0)') '- Number of conduction bands: ', xct%ncb_co
    write(6,'(1x,a,i0)') '- Number of spins: ', xct%nspin
    write(6,'(1x,a,i0)') '- Number of spinors: ', xct%nspinor   ! ZL: add
    write(6,*)
  endif

!------------------------
! Initialize GPU
  if (w_sum_algo == OPENACC_ALGO .or. g_sum_algo == OPENACC_ALGO) then
    call initialize_gpu(OPENACC_ALGO)
  end if
  if (w_sum_algo == OMP_TARGET_ALGO .or. g_sum_algo == OMP_TARGET_ALGO) then
    call initialize_gpu(OMP_TARGET_ALGO)
  end if

!----------------------------------
! Read eps0mat and epsmat
! FHJ: TODO - write epsilon header into bsemat.h5 file
  if (xct%theory.eq.0) then
    call logit('Calling epscopy')
    call timacc(3,1)
    call epscopy(crys,gvec,syms,qg,xct,.false.)
!#BEGIN_INTERNAL_ONLY
    if (xct%screen_exchange) then
      call epscopy(crys,gvec,syms,qg,xct,.true.)
    endif
!#END_INTERNAL_ONLY
    call timacc(3,2)
  elseif (xct%theory.eq.1) then
    call logit('Calling tddft_bz_gen')
    call timacc(3,1)
    call tddft_bz_gen(crys,syms,qg,xct)
    call timacc(3,2)
  endif

!------------------ Initialize BSE Arrays ----------------------------------------------
  if ( xct%ivpar .eq. 1) then
    iownsize=1
  else if ( xct%icpar .eq. 1) then
    iownsize=(xct%nvb_co)**2
  else
    iownsize=(xct%n1b_co*xct%n2b_co)**2
  endif

  if (xct%theory .eq. 0) then
    if (g_sum_algo /= CPU_ALGO) then
      call allocate_accel_kernel(bsedbody, bsedhead, bsedwing, &
                               peinf%myown * iownsize, xct%nspin)
    else
      SAFE_ALLOCATE(bsedbody, (peinf%myown*iownsize,xct%nspin,xct%nspin))
      SAFE_ALLOCATE(bsedhead, (peinf%myown*iownsize,xct%nspin,xct%nspin))
      SAFE_ALLOCATE(bsedwing, (peinf%myown*iownsize,xct%nspin,xct%nspin))
      bsedbody(:,:,:) = ZERO
      bsedhead(:,:,:) = ZERO
      bsedwing(:,:,:) = ZERO
    end if
  else if (xct%theory .eq. 1) then
    SAFE_ALLOCATE(bsedbody, (peinf%myown*iownsize,xct%nspin,xct%nspin))
    SAFE_ALLOCATE(bsedhead, (peinf%myown*iownsize,xct%nspin,xct%nspin))
    SAFE_ALLOCATE(bset, (peinf%myown*iownsize,xct%nspin,xct%nspin))
    bsedbody(:,:,:) = ZERO
    bsedhead(:,:,:) = ZERO
    bset(:,:,:) = ZERO
  endif
  ! We always need exchange
  SAFE_ALLOCATE(bsex, (peinf%myown*iownsize,xct%nspin,xct%nspin))
  bsex(:,:,:) = ZERO

!--------- Calculate Needed Coulomb Interaction -------------------------------

  SAFE_ALLOCATE(vcoul, (xct%ng))
  if (xct%qflag.eq.1) then
    SAFE_ALLOCATE(vcoularray, (xct%ng,qg%nf))
  else
    ! DYQ: if finite Q is used, we will need vcoul at q=-Q.
    SAFE_ALLOCATE(vcoularray, (xct%ng,qg%nf+1))
  endif
  vcoularray=0d0

  ! This array is for the direct part of the TDDFT and
  ! vcoularray will be for the exchange. This part will
  ! get modified if one is doing TD-Hybrids
  if (xct%theory .eq. 1) then
    SAFE_ALLOCATE(vcoularray_mod, (xct%ng,qg%nf))
    vcoularray_mod=0d0
  endif
  SAFE_ALLOCATE(isrtq, (xct%ng))
  do ijk=1,xct%ng
    isrtq(ijk) = ijk
  enddo
  
  if (peinf%inode==0) write(6,'(/1x,a)') 'Calculating Coulomb potential.'
  avgcut=TOL_ZERO
  q0len = sqrt(DOT_PRODUCT(xct%q0vec,MATMUL(crys%bdot,xct%q0vec)))
  iparallel=1

  xct%qpg0_ind = 1 !< store index for smallest G vector here.
  ! DYQ: if using finite Q, store v(-Q+G) as the last element of vcoularray
  if (xct%qflag.ne.1) then
    ! JBH: Remember that Q = -xct%finiteQ, but we want v(-Q+G)
    call vcoul_generator(xct%icutv,xct%truncval,gvec, &
      crys%bdot,crys%celvol,kg%nf,xct%ng,isrtq,xct%iscreen,xct%finiteQ,xct%q0vec, &
      vcoul,xct%iwritecoul,iparallel,avgcut,oneoverq, &
      kp%kgrid,epsheaddummy,work_scell,.false.,wcoul0dummy)
    vcoularray(:,qg%nf+1)=vcoul(:)

    ! JBH: if using finite Q, find index of G vector corresponding to smallest |-Q+G|^2
    ! so we can zero out later if we choose to run the finite Q code without the energy_loss flag.
    ! DYQ: Note, we should NEVER use finite Q without the energy_loss flag. Results without energy
    ! loss flag are not physical.
    gpq0_len = 1.0 / TOL_ZERO
    do ig=1,xct%ng
      qpg = gvec%components(:,ig) + xct%finiteq
      qpg_len = DOT_PRODUCT(qpg,MATMUL(crys%bdot,qpg))
      if (qpg_len < gpq0_len) then
          gpq0_len = qpg_len
          xct%qpg0_ind = ig !< store index for smallest G vector here.
      endif
    enddo
  endif

  do ik=1,qg%nf
    vq(:)=qg%f(:,ik)
    vq0 = DOT_PRODUCT(vq,MATMUL(crys%bdot,vq))
    
    if (peinf%verb_debug .and. peinf%inode==0) then
      write(6,'(1x,a,2i0)') 'Calculating Vcoul', ik, qg%nf
    endif

    call vcoul_generator(xct%icutv,xct%truncval,gvec, &
      crys%bdot,crys%celvol,kg%nf,xct%ng,isrtq,xct%iscreen,vq,xct%q0vec, &
      vcoul,xct%iwritecoul,iparallel,avgcut,oneoverq, &
      kp%kgrid,epsheaddummy,work_scell,.false.,wcoul0dummy)
    vcoularray(:,ik)=vcoul(:)

    if (xct%theory .eq. 1) then
      if (xct%coul_mod_flag) then
        call vcoul_generator(xct%icutv,xct%truncval,gvec, &
          crys%bdot,crys%celvol,kg%nf,xct%ng,isrtq,xct%iscreen,vq,xct%q0vec, &
          vcoul,xct%iwritecoul,iparallel,avgcut,oneoverq, &
          kp%kgrid,epsheaddummy,work_scell,.false.,wcoul0dummy,coulomb_mod=xct%coulomb_mod)
        vcoularray_mod(:,ik)=vcoul(:)
      else
        vcoularray_mod(:,ik)=vcoul(:)
      endif
    endif

    if (vq0 .lt. TOL_Zero) then
      if (peinf%inode .eq. 0) then
        write(6,801) xct%q0vec
      endif
      vq(:) = xct%q0vec(:)
801 format(1x,'Note: for G=0: setting q0 =',3f10.6)

      if (xct%theory .ne. 1) call check_screening_trunc(xct%icutv,xct%iscreen,xct%q0vec,crys%bdot)
    
      vcoul0(1)=0d0
      call vcoul_generator(xct%icutv,xct%truncval,gvec, &
      crys%bdot,crys%celvol,kg%nf,1,isrtq,xct%iscreen,vq,xct%q0vec, &
        vcoul0(:),xct%iwritecoul,iparallel,avgcut,oneoverq, &
        kp%kgrid,epsheaddummy,work_scell,.false.,wcoul0dummy)

      vcoularray(1,ik)=vcoul0(1)

      if (xct%theory .eq. 1) then
        if (xct%coul_mod_flag) then
          vcoul0(1)=0d0
          call vcoul_generator(xct%icutv,xct%truncval,gvec, &
            crys%bdot,crys%celvol,kg%nf,1,isrtq,xct%iscreen,vq,xct%q0vec, &
            vcoul0(:),xct%iwritecoul,iparallel,avgcut,oneoverq, &
            kp%kgrid,epsheaddummy,work_scell,.false.,wcoul0dummy,coulomb_mod=xct%coulomb_mod)

          vcoularray_mod(1,ik)=vcoul0(1)
        else
          vcoularray_mod(1,ik)=vcoul0(1)
        endif
      endif
    endif

  enddo ! ik

  SAFE_DEALLOCATE(isrtq)
  call destroy_qran()

  vcoularray=vcoularray/(8d0*PI_D)
  if (xct%theory .eq. 1) then
    vcoularray_mod=vcoularray_mod/(8d0*PI_D)
  endif
  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,'(1x,a,2i0)') 'Finished Vcoul', ik, qg%nf
  endif

  SAFE_DEALLOCATE(vcoul)

  ! Do a dummy FFT setup to get the size for the FFT box buffer
  call setup_FFT_sizes(gvec%FFTgrid, Nfft, dummy_scale)
  n_bands_max = max(xct%nvb_co, xct%ncb_co)
  call allocate_accel_mtxel_ker(Nfft, gvec%components, kp%ngkmax, kp%ngkmax, &
                              gvec%ng, xct%ng, n_bands_max, n_bands_max, &
                              xct%n_ffts_per_batch, mtxel_algo)

!--------- Start the computation. ---------------------------------------------

  call logit('Starting main kernel loop')
  
  SAFE_ALLOCATE(qqa, (peinf%myown))
  SAFE_ALLOCATE(fqa, (3,peinf%myown))
  SAFE_ALLOCATE(g0a, (3,peinf%myown))
  SAFE_ALLOCATE(irqa, (peinf%myown))
  SAFE_ALLOCATE(ifqa, (peinf%myown))
  
  call sortbyq(fqa,qqa,g0a,ifqa,irqa,qg,kg,crys)

! JRD: The below may be useful for debugging
!
!  call mpi_barrier(mpi_comm_world,mpierr)
!  write(6,*) 'myown',peinf%inode,peinf%myown
!  call mpi_barrier(mpi_comm_world,mpierr)

  ! FHJ: this is to generate nice output / time estimate
  call progress_init(prog_info, 'calculation of matrix elements', 'blocked transition', peinf%nckpe)

  do ii=1,peinf%nckpe

    ! FHJ : friendly output / running time estimate
    call progress_step(prog_info, ii)

    call logitint('   Main loop:  ii=',ii)
    if (ii .le. peinf%myown) then
      ik=peinf%ik(peinf%inode+1,ii)
      ikp=peinf%ikp(peinf%inode+1,ii)
      
      ! DYQ: If using patched sampling and finite Q, skip k+Q if it falls outside
      ! the patch
      if (xct%patched_sampling_co .and. xct%qflag.eq.2) then
        if (xct%indexq(ik).eq.0 .or. xct%indexq(ikp).eq.0) then 
          if (peinf%inode.eq.0) write(6,*) "Skipping: ik,ikp",ik,ikp
          cycle
        endif
      endif

      if (xct%icpar .eq. 1) then
        ic=peinf%ic(peinf%inode+1,ii)
        icp=peinf%icp(peinf%inode+1,ii)
      else
        ic=1
        icp=1
      endif
      
      if (xct%ivpar .eq. 1) then
        iv=peinf%iv(peinf%inode+1,ii)
        ivp=peinf%ivp(peinf%inode+1,ii)
      else
        iv=1 
        ivp=1 
      endif
      
      call timacc(4,1)
      call logit('   Calling genwf_kernel')
      
      !write(6,*) peinf%inode,'calling gw 1', ii, ik
      call genwf_kernel(crys,gvec,kg,kgq,syms,wfnc, &
        wfnv,xct%nspin,ik,ic,iv,indexq,xct,intwfnv,intwfnc)
      
      !write(6,*) peinf%inode,'calling gw 2', ii, ikp

      ! returns wfnv at indexq(ikp) and wfnc at ikp
      call genwf_kernel(crys,gvec,kg,kgq,syms,wfncp, &
        wfnvp,xct%nspin,ikp,icp,ivp,indexq,xct,intwfnv,intwfnc)

      call timacc(4,2)
      
    else
      ik=-1
      ikp=-1
      ic=-1
      icp=-1
      iv=-1
      ivp=-1
    endif
    if (peinf%verb_debug .and. peinf%inode.eq.0) then
      write(6,772) ii
      write(6,773) ik,ikp,ic,icp,iv,ivp
772 format(1x,"PE # 0 dealing with block",i6)
773 format(1x,"ik =",i6,1x,"ikp =",i6,1x,"ic =",i6,1x,"icp =",i6,1x,"iv =",i6,1x,"ivp =",i6)
    endif

    call logit('      Calling mtxel_kernel')
    call timacc(6,1)

    if (ii .le. peinf%myown) then
      if (xct%theory .eq. 0) then
        call mtxel_kernel(crys,gvec,syms,qg,wfnc,wfncp,wfnvp, &
          wfnv,xct,peinf%myown*iownsize,bsedbody,bsedhead,bsedwing,bsex,ii,ik,ikp, &
          ic,icp,iv,ivp, &
          vcoularray,fqa(:,ii),qqa(ii),g0a(:,ii),ifqa(ii),irqa(ii),q0len)
      else if (xct%theory .eq. 1) then
!#BEGIN_INTERNAL_ONLY
        call mtxel_tddft(crys,gvec,syms,qg,wfnc,wfncp,wfnvp, &
          wfnv,xct,peinf%myown*iownsize,bsedbody,bsedhead,bsex,bset,ii,ik,ikp, &
          ic,icp,iv,ivp,vcoularray_mod, &
          vcoularray,fqa(:,ii),qqa(ii),g0a(:,ii),ifqa(ii),irqa(ii),q0len)
!#END_INTERNAL_ONLY
      endif

      SAFE_DEALLOCATE_P(wfncp%cg)
      SAFE_DEALLOCATE_P(wfncp%isort)
      SAFE_DEALLOCATE_P(wfnvp%cg)
      SAFE_DEALLOCATE_P(wfnvp%isort)
      SAFE_DEALLOCATE_P(wfnc%cg)
      SAFE_DEALLOCATE_P(wfnc%isort)
      SAFE_DEALLOCATE_P(wfnv%cg)
      SAFE_DEALLOCATE_P(wfnv%isort)
    else
      !            write(6,*) peinf%inode, 'Calling mtxel without task'
      fqa_dummy(:)=0.0d0
      qqa_dummy=0.0d0
      g0a_dummy(:)=0
      ifqa_dummy=0
      irqa_dummy=0
      if (xct%theory .eq. 0) then
        call mtxel_kernel(crys,gvec,syms,qg,wfnc,wfncp,wfnvp,wfnv, &
          xct,peinf%myown*iownsize,bsedbody,bsedhead,bsedwing,bsex,ii,ik,ikp,ic,icp, &
          iv,ivp, &
          vcoularray,fqa_dummy,qqa_dummy,g0a_dummy,ifqa_dummy,irqa_dummy,q0len)
      else if (xct%theory .eq. 1) then
!#BEGIN_INTERNAL_ONLY
        call mtxel_tddft(crys,gvec,syms,qg,wfnc,wfncp,wfnvp,wfnv, &
          xct,peinf%myown*iownsize,bsedbody,bsedhead,bsex,bset,ii,ik,ikp,ic,icp, &
          iv,ivp,vcoularray_mod, &
          vcoularray,fqa_dummy,qqa_dummy,g0a_dummy,ifqa_dummy,irqa_dummy,q0len)
!#END_INTERNAL_ONLY
      endif
    endif
    
    call timacc(6,2)

  enddo !ii

  if (xct%qflag .eq. 0) call dealloc_grid(kgq)
  call dealloc_grid(qg)
  call progress_free(prog_info)

  SAFE_DEALLOCATE_P(intwfnv%cg)
  SAFE_DEALLOCATE_P(intwfnc%cg)
  SAFE_DEALLOCATE_P(intwfnv%isort)
  SAFE_DEALLOCATE_P(intwfnc%isort)
  SAFE_DEALLOCATE_P(intwfnv%ng)
  SAFE_DEALLOCATE_P(intwfnc%ng)
  
  SAFE_DEALLOCATE(qqa)
  SAFE_DEALLOCATE(fqa)
  SAFE_DEALLOCATE(g0a)
  SAFE_DEALLOCATE(irqa)
  SAFE_DEALLOCATE(ifqa)
  if (xct%theory .eq. 0) then
    SAFE_DEALLOCATE_P(xct%epsdiag)
  endif
  if (peinf%inode .eq. 0 .or. xct%bLowComm) then
    SAFE_DEALLOCATE_P(xct%isrtqi)
  endif
  if (xct%theory .eq. 0) then
    SAFE_DEALLOCATE_P(xct%epscol)
  endif
  SAFE_DEALLOCATE_P(peinf%nxqown)
  SAFE_DEALLOCATE_P(peinf%nxqi)
  

!--------------- Write BSE matrices ----------------------------------------------

#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif
  call timacc(8,1)
  if (xct%theory .eq. 0) then
    call accel_update_from(bsedhead, g_sum_algo)
    call accel_update_from(bsedbody, g_sum_algo)
    call accel_update_from(bsedwing, g_sum_algo)
    call bsewrite(xct,iownsize,bsedbody,bsedhead,bsex,kg,kp,gvec,syms,crys, &
         bsedwing=bsedwing)
  else if (xct%theory .eq. 1) then
    call bsewrite(xct,iownsize,bsedbody,bsedhead,bsex,kg,kp,gvec,syms,crys, &
      bset=bset)
  endif
  call timacc(8,2)

  call dealloc_grid(kg)

  call destroy_fftw_plans()

  if (xct%theory .eq. 0) then
    if (g_sum_algo /= CPU_ALGO) then
      call deallocate_accel_kernel(bsedbody, bsedhead, bsedwing)
    else
      SAFE_DEALLOCATE(bsedhead)
      SAFE_DEALLOCATE(bsedwing)
      SAFE_DEALLOCATE(bsedbody)
    end if
  else if (xct%theory .eq. 1) then
    SAFE_DEALLOCATE(bsedhead)
    SAFE_DEALLOCATE(bsedbody)
    SAFE_DEALLOCATE(bset)
  endif
  SAFE_DEALLOCATE(bsex)
  SAFE_DEALLOCATE(indexq)

  if (xct%iwritecoul .eq. 1) then
    if (peinf%inode .eq. 0) then
      call close_file(19) ! file vcoul
    endif
  endif

  call show_references()

!---------------- Time accounting -----------------------------------------------


  ntim=8
  SAFE_ALLOCATE(routnam, (76))
  routnam(1)='TOTAL:'
  routnam(2)='INPUT:'
  routnam(3)='EPSCOPY:'
  routnam(4)='GENWF:'
  routnam(5)='EXCWF:'
  routnam(6)='MTXEL:'
  routnam(7)='FULLBZ:'
  routnam(8)='BSEWRITE:'
  
  call timacc(1,2)
  if(peinf%inode.eq.0) then
    write(6,*)
    write(6,9000) 'CPU (s)','WALL (s)','#'
    write(6,*)
  endif
  
  do ii=2,ntim
    call timacc(ii,3,tsec,ncount)
#ifdef MPI
    call MPI_ALLREDUCE(tsec,tmin,2,MPI_REAL_DP,MPI_MIN,MPI_COMM_WORLD,mpierr)
    call MPI_ALLREDUCE(tsec,tmax,2,MPI_REAL_DP,MPI_MAX,MPI_COMM_WORLD,mpierr)
#else
    tmin = tsec
    tmax = tsec
#endif
    if(peinf%inode.eq.0) then
      write(6,9001) routnam(ii),tmin(1),tmin(2),ncount
      write(6,9002) tsec(1),tsec(2)
      write(6,9003) tmax(1),tmax(2)
    endif
  enddo

! JRD More Time Accounting for mtxel_kernel

  routnam(61)='MTXEL Setup:'
  routnam(62)='MTXEL Vcoul:'
  routnam(63)='MTXEL W:'
  routnam(64)='MTXEL W-Sum:'
  routnam(65)='MTXEL FFT Dir:'
  routnam(66)='MTXEL High G:'
  routnam(67)='MTXEL FFT X:'
  routnam(68)='MTXEL BSEX:'
  routnam(69)='MTXEL INDEX:'
  routnam(70)='MTXEL EPSHEAD:'
  routnam(71)='MTXEL EPSREAD:'
  routnam(72)='MTXEL EPSOPEN:'
  routnam(73)='MTXEL G SUM:'
  routnam(74)='MTXEL COMM:'
  routnam(75)='BSEWRITE COMM:'
  routnam(76)='BSEWRITE I/O:'
  
  do ii=61,76
    call timacc(ii,3,tsec,ncount)
#ifdef MPI
    call MPI_Reduce(tsec,tmin,2,MPI_REAL_DP,MPI_MIN,0,MPI_COMM_WORLD,mpierr)
    call MPI_Reduce(tsec,tmax,2,MPI_REAL_DP,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#else
    tmin=tsec
    tmax=tsec
#endif
    if(peinf%inode.eq.0) then
      write(6,9001) routnam(ii),tmin(1),tmin(2),ncount
      write(6,9002) tsec(1),tsec(2)
      write(6,9003) tmax(1),tmax(2)
    endif
  enddo

  routnam(30)='MTXEL FFT CC:'
  routnam(31)='MTXEL FFT VV:'
  routnam(32)='MTXEL FIND:'
  routnam(33)='MTXEL ALLOC:'
  routnam(34)='MTXEL FFTBOX:'
  routnam(35)='MTXEL FFT VC:'

  do ii=30,35
    call timacc(ii,3,tsec,ncount)
#ifdef MPI
    call MPI_Reduce(tsec,tmin,2,MPI_REAL_DP,MPI_MIN,0,MPI_COMM_WORLD,mpierr)
    call MPI_Reduce(tsec,tmax,2,MPI_REAL_DP,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#else
    tmin=tsec
    tmax=tsec
#endif
    if(peinf%inode.eq.0) then
      write(6,9001) routnam(ii),tmin(1),tmin(2),ncount
      write(6,9002) tsec(1),tsec(2)
      write(6,9003) tmax(1),tmax(2)
    endif
  enddo
  
  if(peinf%inode.eq.0) write(6,*)
  
  call timacc(1,3,tsec,ncount)
#ifdef MPI
  call MPI_Reduce(tsec,tmin,2,MPI_REAL_DP,MPI_MIN,0,MPI_COMM_WORLD,mpierr)
  call MPI_Reduce(tsec,tmax,2,MPI_REAL_DP,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#else
  tmin=tsec
  tmax=tsec
#endif
  
  if(peinf%inode.eq.0) then
    write(6,9004) routnam(1),tmin(1),tmin(2)
    write(6,9002) tsec(1),tsec(2)
    write(6,9003) tmax(1),tmax(2)
    write(6,*)
  endif
  
9000 format(23x,a13,3x,a13,3x,a8)
9001 format(1x,a16,'(min.)',f13.3,3x,f13.3,3x,i8)
9002 format(   17x,'(PE 0)',f13.3,3x,f13.3)
9003 format(   17x,'(max.)',f13.3,3x,f13.3)
9004 format(1x,a16,'(min.)',f13.3,3x,f13.3)
  
  call write_memory_usage()

  call deallocate_accel_mtxel_ker(mtxel_algo)

#ifdef MPI
  call MPI_FINALIZE(mpierr)
#endif

end program kernel
