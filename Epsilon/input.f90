!================================================================================
!
! Routines:
!
! (1) input()   Originally By ?                 Last Modified 7/7/2008
!
! Read in and setup up various data structures,
! read in and distributes wave functions, do all sorts of setup, etc.
!
!===============================================================================

#include "f_defs.h"

module input_m

  use global_m
  use checkbz_m
  use createpools_m
  use eqpcor_m
  use fftw_m
  use fullbz_m
#ifdef HDF5
  use hdf5
#endif
  use hdf5_io_m
  use wfn_io_hdf5_m
  use epswrite_hdf5_m
  use input_utils_m
  use io_utils_m
  use irrbz_m
  use misc_m
  use scissors_m
  use sort_m
  use wfn_rho_vxc_io_m

  use timing_m, only: timing => epsilon_timing
  use inread_m
  use subgrp_m

  implicit none

  private
  public :: input

contains

  subroutine input(kp,crys,syms,gvec,pol,cwfn,vwfn,intwfnv,intwfnc,omega_plasma, gr)
    type (kpoints), intent(out) :: kp
    type (crystal), intent(out) :: crys
    type (symmetry), intent(out) :: syms
    type (gspace), intent(out) :: gvec
    type (polarizability), intent(out) :: pol
    type (conduction_wfns), intent(out) :: cwfn
    type (valence_wfns), intent(out) :: vwfn
    type (int_wavefunction), intent(out) :: intwfnv
    type (int_wavefunction), intent(out) :: intwfnc
    real(DP), intent(out) :: omega_plasma
    type(grid), intent(out) :: gr

    integer, allocatable :: global_nvown_temp(:),global_ncown_temp(:)
    integer, allocatable :: global_pairowner_temp(:,:),global_indexv_temp(:,:)

    character :: aheadinput*60, ajname*6, adate*11, atime*14
    character :: filenameh5*80
    character :: tmpfn*16
    character :: fncor*32
    character :: tmpstr1*16,tmpstr2*16,tmpstr3*16,tmpstr*120
    integer :: ii,ig,jj,itran,ib,ispin,ik,nband_temp,ncrit_band_temp,ncrit_band,ncrit_temp
    integer :: error
    integer :: nrq,nrqmax,iq,ic,npools,mypool,mypoolrank
    integer :: iv,nrkq,myipe,ipool
    integer :: ipe,iw
    integer :: nmtx,nmtx_l,nmtx_t,nmtx0,npairs
    integer :: Nrod,Nplane,Nfft(3),dNfft(3),dkmax(3),nmpinode
    real(DP) :: occ, occq
    real(DP) :: qk(3),vcell,qtot,qshift(3),norm
    real(DP) :: mem,fmem,rmem,rmem2,smem,scale,dscale
    real(DP), allocatable :: energies(:,:,:)
    type(grid) :: qgr

    character(len=3) :: sheader
    integer :: iflavor, matrix_type
    type(gspace) :: gvecq
    type(kpoints) :: kpq
    type(crystal) :: crysq
    type(symmetry) :: symsq

    logical :: skip_checkbz

    PUSH_SUB(input)

!-------------------------------
! SIB: Read the input file

    call inread(pol,vwfn,cwfn)

! DVF: Set up the variables needed for parallel frequencies, or set the variables
! to the corresponding values to do non-parallel frequency calculations. This
! subroutine REDEFINES peinf%npes to be peinf%npes/pol%nfreq_group*pol%nfreq_group
! and defines peinf%npes_orig the original peinf%npes, for the few instances
! where we still need the original number of processors. There is no effect on
! peinf%npes for static calculations, calculations without parallel frequencies,
! and calculations for which pol%nfreq_group divides peinf%npes.

    call setup_parallel_freqs(pol)

!-------------------------------
! JRD: Initial hdf interface
#ifdef HDF5
    if(pol%use_hdf5 .or. pol%wfn_hdf5) call h5open_f(error)
#endif

! FHJ: Read header of WFN file (and WFNq, if requred), perform consistency
! checks, and figure out the number of bands and ncrit.
    sheader = 'WFN'
    iflavor = 0
    call timing%start(timing%io_total)
    if(pol%wfn_hdf5) then
#ifdef HDF5
      if (peinf%inode==0) write(6,'(a)') ' Reading header of WFN.h5'
      call read_hdf5_header_type('WFN.h5', sheader, iflavor, kp, gvec, syms, crys)
#endif
    else
      if(peinf%inode == 0) then
        write(6,'(a)') ' Reading header of WFN'
        call open_file(25,file='WFN',form='unformatted',status='old')
      endif
      !call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys, &
        !dont_warn_kgrid=pol%subsample)
      call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys, &
        dont_warn_kgrid=.true.)   ! WK: Off the warning for tilted grid.
    endif
    call timing%stop(timing%io_total)
    call check_trunc_kpts(pol%icutv, kp)
    if (cwfn%nband==0) cwfn%nband = kp%mnband - 1
! DVF: Even though it is included in the conduction wavefunction structure, the nband
! in cwfn is the total number of bands. When we exclude deep core states, we need to
! redefine this number.
    cwfn%nband=cwfn%nband-vwfn%ncore_excl

    call scissors_shift(kp, pol%scis)
    if (pol%eqp_corrections) then
      fncor = 'eqp.dat'
      call eqpcor(fncor, peinf%inode, peinf%npes, kp, 1+vwfn%ncore_excl, cwfn%nband+vwfn%ncore_excl, 0, 0, &
        kp%el, kp%el, kp%el, 1, 0)
    endif
    call find_efermi(pol%rfermi, pol%efermi, pol%efermi_input, kp, cwfn%nband+vwfn%ncore_excl, 1+vwfn%ncore_excl, &
      "unshifted grid", should_search = .true., should_update = .true., write7 = .true.)
  
    call find_occupations(kp%el, kp%nspin, kp%nrk, kp%mnband, pol%efermi, pol%occ_broadening, pol%occ, pol%do_smearing)

    if (pol%need_WFNq) then
      call timing%start(timing%io_total)
      if(pol%wfn_hdf5) then
#ifdef HDF5
        if (peinf%inode==0) write(6,'(a)') ' Reading header of WFNq.h5'
        call read_hdf5_header_type('WFNq.h5', sheader, iflavor, kpq, gvecq, symsq, crysq)
#endif
      else
        if(peinf%inode == 0) then
          write(6,'(a)') ' Reading header of WFNq'
          call open_file(26,file='WFNq',form='unformatted',status='old')
        endif
        !call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, symsq, crysq, &
        !  dont_warn_kgrid=pol%subsample, warn=.false.)

        call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, symsq, crysq, &
          dont_warn_kgrid=.true., warn=.false.)

        if(peinf%inode == 0) call close_file(26)
      endif
      call timing%stop(timing%io_total)

      nrkq = kpq%nrk
      call check_trunc_kpts(pol%icutv, kpq)
      call check_header('WFN', kp, gvec, syms, crys, 'WFNq', kpq, gvecq, symsq, crysq, is_wfn=.true.)
      if (any(pol%qgrid==0) .and. any(kp%kgrid(1:3)/=kpq%kgrid(1:3)) .and. .not.pol%subsample) then
        if(peinf%inode == 0) then
          write(0,*) 'ERROR: kgrids for WFN and WFNq are not the same'
          write(0,*) 'WFN  kgrid = ', kp%kgrid(1:3)
          write(0,*) 'WFNq kgrid = ', kpq%kgrid(1:3)
          write(0,*) 'Note: you can use the `qgrid` keyword if this behavior is really desired.'
        endif
        call die('kgrids for WFN and WFNq must be the same', only_root_writes=.true.)
      endif
      call scissors_shift(kpq, pol%scis)
      if (peinf%inode==0) then
        if(pol%eqp_corrections) then
          fncor='eqp_q.dat'
          ! FHJ: We should technically ask for vwfn%nband+pol%ncrit bands here.
          ! But we don`t know what`s vwfn%nband+pol%ncrit. However, unless you use
          ! patched_sampling or intraband_flag, we can determine vwfn%nband+pol%ncrit
          ! from WFN alone, which is maxval(kp%ifmax).
          call eqpcor(fncor,0,1,kpq,1+vwfn%ncore_excl,maxval(kp%ifmax),0,0,kpq%el,kpq%el,kpq%el,1,0)
        endif
        call find_efermi(pol%rfermi, pol%efermi, pol%efermi_input, kpq, kpq%mnband, 1+vwfn%ncore_excl, &
        "shifted grid", should_search = .false., should_update = .false., write7 = .true.)
        if (.not.pol%patched_sampling .and. pol%intraband_flag==0) then
          if (minval(kpq%ifmin)/=minval(kp%ifmin) .or. maxval(kpq%ifmin)/=maxval(kp%ifmin)) then
            write(0,*) 'ERROR: occupations (ifmin) inconsistent between WFN and WFNq files.'
            write(0,*) 'Remember that you should NOT use WFNq for metals and graphene.'
            call die('occupations (ifmin field) inconsistent between WFN and WFNq files.', &
              only_root_writes=.true.)
          endif
          if (minval(kpq%ifmax)/=minval(kp%ifmax) .or. maxval(kpq%ifmax)/=maxval(kp%ifmax)) then
            write(0,*) 'ERROR: occupations (ifmax) inconsistent between WFN and WFNq files.'
            write(0,*) 'Remember that you should NOT use WFNq for metals and graphene.'
            call die('occupations (ifmax field) inconsistent between WFN and WFNq files.', &
              only_root_writes=.true.)
          endif
        endif
      endif
      call find_occupations(kpq%el, kpq%nspin, kpq%nrk, kpq%mnband, pol%efermi, pol%occ_broadening, pol%occq, pol%do_smearing)
      if (vwfn%nband==0) then
        if (pol%do_smearing) then
          call get_nband_ncrit(kp%nspin, kp%nrk, kp%mnband, vwfn%nband, pol%ncrit, pol%occ, kpq=kpq, occqs=pol%occq)
        else ! no smearing
          vwfn%nband = min(minval(kpq%ifmax), minval(kp%ifmax))
          pol%ncrit = max(maxval(kpq%ifmax), maxval(kp%ifmax)) - vwfn%nband
        endif
      endif
    call dealloc_header_type(sheader, crysq, kpq)
    else ! FHJ: No WFNq case below
      if (vwfn%nband==0) then
        if (pol%do_smearing) then
          call get_nband_ncrit(kp%nspin, kp%nrk, kp%mnband, vwfn%nband, pol%ncrit, pol%occ)
        else ! No smearing
          vwfn%nband = min(minval(kp%ifmax), minval(kp%ifmax))
          pol%ncrit = max(maxval(kp%ifmax), maxval(kp%ifmax)) - vwfn%nband
        endif
      endif
    endif
! DVF: Here we re-define the number of valence bands to exclude the deep core
! states. This ensures that the pools will be set up properly. This achieves
! most of what is needed to exclude core states. There is a little more work
! in the i/o routines to make sure you`re getting the higher valence states.
    vwfn%nband=vwfn%nband-vwfn%ncore_excl

    if (any(pol%qgrid==0)) pol%qgrid = kp%kgrid
    if (peinf%inode==0) then
      write(6,'(/1x,a)') 'Calculation parameters:'
      write(6,'(1x,a,f0.2)') '- Cutoff of the dielectric matrix (Ry): ', pol%ecuts
      write(6,'(1x,a,i0)') '- Total number of bands in the calculation: ', cwfn%nband
      write(6,'(1x,a,i0)') '- Number of fully occupied valence bands: ', vwfn%nband
      write(6,'(1x,a,i0)') '- Number of partially occ. conduction bands: ', pol%ncrit
      write(6,'(1x,a,3(1x,i0))') '- Monkhorst-Pack q-grid for epsilon(q):', pol%qgrid
      write(6,*)
    endif
#ifdef MPI
    call MPI_Bcast(vwfn%nband, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(pol%ncrit, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
    call distribution()

!---------------------------------
! Determine the available memory
    if(pol%nfreq_group .eq. 1) then
      call procmem(mem,nmpinode)
    else
      call procmem(mem,nmpinode,pol%nfreq_group)
    endif

    if(peinf%inode.eq.0) then
      write(6,998) mem / 1024**2
      write(7,998) mem / 1024**2
    endif
998 format(1x,'Memory available: ',f0.1,' MB per PE')
    fmem = mem / 8


    SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
    call timing%start(timing%io_total)
    if(pol%wfn_hdf5) then
#ifdef HDF5
      call read_hdf5_gvectors('WFN.h5', gvec%ng, gvec%components)
#endif
    else
      call read_binary_gvectors(25, gvec%ng, gvec%ng, gvec%components)
    endif
    call timing%stop(timing%io_total)

!-------------------------------
! (gsm) Estimate the required memory

    nmtx=0

    gr%nr = kp%nrk
    SAFE_ALLOCATE(gr%r, (3, gr%nr))
    gr%r = kp%rk
    call timing%start(timing%fullbz)
    if (pol%patched_sampling) then
      call fullbz(crys,syms,gr,1,skip_checkbz,wigner_seitz=.false.,paranoid=.true.)
    else
      call fullbz(crys,syms,gr,syms%ntran,skip_checkbz,wigner_seitz=.false.,paranoid=.true.)
    endif
    tmpfn='WFN'
    if (.not.skip_checkbz.and..not.pol%patched_sampling) then
      call checkbz(gr%nf,gr%f,kp%kgrid,kp%shift,crys%bdot, &
        tmpfn,'k',.false.,pol%freplacebz,pol%fwritebz)
    endif

    ! do not worry if this is the separate q->0 calculation for metals
    if(pol%nq > 1 .or. pol%valueq0 /= 2) then
      qgr%nr = pol%nq
      SAFE_ALLOCATE(qgr%r, (3, pol%nq))
      qgr%r(1:3,1:pol%nq) = pol%qpt(1:3,1:pol%nq)
      if(pol%nq0>0) qgr%r(1:3,1) = 0d0
      if (pol%patched_sampling) then
        call fullbz(crys,syms,qgr,1,skip_checkbz,wigner_seitz=.false.,paranoid=.true.)
      else
        call fullbz(crys,syms,qgr,syms%ntran,skip_checkbz,wigner_seitz=.false.,paranoid=.true.)
      endif
      qshift(:)=0.0d0
      if (.not.skip_checkbz.and..not.pol%subsample.and..not.pol%patched_sampling) then
        call checkbz(qgr%nf,qgr%f,kp%kgrid,qshift,crys%bdot,'epsilon.inp','q',.false.,pol%freplacebz,pol%fwritebz)
      endif
      call dealloc_grid(qgr)
    endif
    call timing%stop(timing%fullbz)

! XXX - Do this in parallel?

! compute max for nmtx
    if(peinf%inode == 0) then
      SAFE_ALLOCATE(pol%nmtx_of_q,(pol%nq))

      nmtx0=-1

! XXX - THREAD?
      do iq = 1, pol%nq
        SAFE_ALLOCATE(gvec%ekin, (gvec%ng))
        if (iq<=pol%nq0) then
          call kinetic_energies(gvec, crys%bdot, gvec%ekin)
        else
          call kinetic_energies(gvec, crys%bdot, gvec%ekin, qvec = pol%qpt(:,iq))
        endif
        SAFE_ALLOCATE(pol%isrtx, (gvec%ng))
        call sortrx(gvec%ng, gvec%ekin, pol%isrtx, gvec = gvec%components)
        nmtx_t = gcutoff(gvec%ng, gvec%ekin, pol%isrtx, pol%ecuts)
        SAFE_DEALLOCATE_P(pol%isrtx)
        SAFE_DEALLOCATE_P(gvec%ekin)
        pol%nmtx_of_q(iq) = nmtx_t
        if (iq<=pol%nq0) then
          nmtx0 = nmtx_t
        endif
        if (nmtx_t .gt. nmtx) nmtx = nmtx_t
      enddo

      nmtx_l = nmtx
      nmtx_l = int(sqrt(dble(nmtx_l)**2 / peinf%npes_freqgrp))
      ! JRD: DUMB Debugging write(6,*) 'computed nmtx_l', nmtx, peinf%npes, nmtx_l, pol%nq

      nrqmax=0
      do iq=1,pol%nq
        nrq=0
        call subgrp(pol%qpt(:,iq),syms)
        call irrbz(syms,gr%nf,gr%f,nrq)
        if (nrq .gt. nrqmax) nrqmax=nrq
      enddo

! required memory
      rmem=0.0d0
! array pol%chi in program chi_summation (and chilocal, chilocal2 (if gcomm=-1))
      if (pol%freq_dep.eq.0) then
        rmem = rmem + 2 * dble(kp%nspin) * dble(nmtx_l)**2
        if (pol%gcomm .eq. -1) then
          rmem = rmem + dble(kp%nspin) * dble(nmtx_l)**2
        endif
      endif
! DVF : Arrays pol%chiRDyn, chilocalRDyn in program chi_summation (their advanced counterparts
! are taken care of via multiplication by sizeof_scalar below). These arrays are present
! for both gcomm=-1 (matrix comm) and 0 (elements comm). The gcomm=-1, you also have
! chilocal2RDyn and its advanced counterpart. The other factor of two comes from the
! fact that these arrays are complex. Note : I think the multiplication of the some
! of the memory by sizeof_scalar is not proper because many of these arrays are
! always complex, so the factor or 2 between real and complex versions of the code
! makes no sense. I guess this is just an estimate, but I think it may be a decent
! bit off because of this.
      if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
        rmem = rmem + dble(kp%nspin) * dble(nmtx_l)**2 * dble(pol%nfreq_in_group) *4
        if (pol%gcomm .eq. -1) then
          rmem = rmem + dble(kp%nspin) * dble(nmtx_l)**2 * dble(pol%nfreq_in_group) * 2
        endif
      endif

      if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
        !pol%chiRDyn
        rmem = rmem + dble(kp%nspin) * dble(nmtx_l)**2 * dble(pol%nfreq_in_group)
        !chilocalRDyn
        rmem = rmem + dble(nmtx_l)**2 * dble(pol%os_nsfreq_para)
        !pol%chiTDyn
        rmem = rmem + dble(kp%nspin) * dble(nmtx_l)**2 * dble(pol%os_nsfreq_para)
        !chilocal2RDyn
        if (pol%gcomm .eq. -1) then
          rmem = rmem + dble(kp%nspin) * dble(nmtx_l)**2 * dble(pol%os_nsfreq_para)
        endif
      endif

      if (pol%freq_dep .eq. 3) then
        rmem = rmem + dble(kp%nspin) * dble(nmtx_l)**2 * dble(pol%nfreq_in_group) *4
        if (pol%gcomm .eq. -1) then
          rmem = rmem + dble(kp%nspin) * dble(nmtx_l)**2 * dble(pol%nfreq_in_group) * 2
        endif
      endif

! arrays gmetempr and gmetempc in program chi_summation
      if (pol%freq_dep.eq.0) then
        if (kp%nrk .eq. 1) then
          rmem = rmem + 2 * dble(nmtx) * dble(vwfn%nband + pol%ncrit) / sqrt(dble(peinf%npes))
        else
          rmem = rmem + 2 * dble(nmtx) * dble(vwfn%nband + pol%ncrit) * dble(syms%ntran) / sqrt(dble(peinf%npes))
        endif
      endif
! arrays gmeR(A)Dyntempr2 and gmeR(A)Dyntempc in program chi_summation
! DVF : I think this is an overestimate because gmeRDyntempc does not
! have frequency dependence
      if ( ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) .or. pol%freq_dep .eq. 3) then
        if (pol%gcomm .eq. 0) then
          if (kp%nrk .eq. 1) then
            rmem = rmem + dble(nmtx) * dble(pol%nfreq + 1) * dble(vwfn%nband + pol%ncrit) / sqrt(dble(peinf%npes))
          else
            rmem = rmem + dble(nmtx) * dble(pol%nfreq + 2) * &
              dble(vwfn%nband + pol%ncrit) * dble(syms%ntran) / sqrt(dble(peinf%npes))
          endif
#ifdef CPLX
          if (pol%freq_dep .ne. 3) then
            if (kp%nrk .eq. 1) then
              rmem = rmem + dble(nmtx) * dble(pol%nfreq + 2) * dble(vwfn%nband + pol%ncrit) / sqrt(dble(peinf%npes))
            else
              rmem = rmem + dble(nmtx) * dble(pol%nfreq + 2) * dble(vwfn%nband + pol%ncrit)&
                * dble(syms%ntran) / sqrt(dble(peinf%npes))
            endif
          endif
#endif
        endif
        npairs=(vwfn%nband+pol%ncrit)*(cwfn%nband-vwfn%nband)
        if (pol%gcomm .eq. -1) then
          if (kp%nrk .eq. 1) then
            rmem = rmem + 2 * dble(nmtx) * dble(npairs) * dble(pol%nfreq_in_group) / peinf%npes_freqgrp**1.5
          else
            rmem = rmem + 2 * dble(nmtx) * dble(kp%nrk) * dble(npairs) * dble(syms%ntran) * dble(pol%nfreq_in_group) &
              / peinf%npes_freqgrp**1.5
          endif
#ifdef CPLX
          if (pol%freq_dep .ne. 3) then
            if (kp%nrk .eq. 1) then
              rmem = rmem + 2 * dble(nmtx) * dble(pol%nfreq) * dble(npairs) / peinf%npes_freqgrp**1.5
            else
              rmem = rmem + 2 * dble(nmtx) * dble(kp%nrk) * dble(npairs) * dble(syms%ntran) * dble(pol%nfreq) &
                / peinf%npes_freqgrp**1.5
            endif
          endif
#endif
        endif
      endif

      if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
        if (pol%gcomm .eq. 0) then
          !gmeRDyntempr2, gmeRDyntempcs, and gmeRDyntempc
          !It is an overestimate because max_nv is less than nv.
          if (kp%nrk .eq. 1) then
            rmem = rmem + dble(nmtx) * dble(pol%nsfreq + 2)&
              * dble(vwfn%nband + pol%ncrit) / sqrt(dble(peinf%npes))
          else
            rmem = rmem + dble(nmtx) * dble(pol%nsfreq + 3)&
              * dble(vwfn%nband + pol%ncrit) * dble(syms%ntran) / sqrt(dble(peinf%npes))
          endif
        endif
        npairs=(vwfn%nband+pol%ncrit)*(cwfn%nband-vwfn%nband)
        if (pol%gcomm .eq. -1) then
          if (kp%nrk .eq. 1) then
            rmem = rmem + 2 * dble(nmtx) * dble(npairs)&
              * dble(pol%os_nsfreq_para) / peinf%npes_freqgrp**1.5
          else
            rmem = rmem + 2 * dble(nmtx) * dble(kp%nrk)&
              * dble(npairs) * dble(syms%ntran) * dble(pol%os_nsfreq_para) &
              / peinf%npes_freqgrp**1.5
          endif
        endif
      endif
! array gmetempn in program chi_summation
      rmem = rmem + dble(nmtx)
! array tmparray in program mtxel
      rmem = rmem + dble(nmtx)
! array pol%gme in program epsilon_main
      rmem = rmem + dble(kp%nspin) * dble(nmtx) * dble(peinf%ncownmax) * &
        dble(peinf%nvownmax) * dble(nrqmax)*dble(pol%nfreq_group)
! arrays zin, vwfn%zv, cwfn%zc, zinc and zinc_old in subroutine genwf
      rmem = rmem + dble(kp%nspin) * dble(2 * 1 + 2 * peinf%ncownmax + 1) * dble(kp%ngkmax)
! intwfn_files
      rmem = rmem + dble(peinf%nvownmax + peinf%ncownmax) * dble(kp%ngkmax) * dble(kp%nrk) * dble(kp%nspin)
      if (pol%need_WFNq) then
        rmem = rmem + dble(peinf%nvownmax) * dble(kp%ngkmax) * dble(nrkq) * dble(kp%nspin)
      endif

      rmem = rmem * sizeof_scalar()

! memory for pol%eden or pol%edenDyn
      if (pol%freq_dep .eq. 0) then
        rmem = rmem + dble(kp%nspin) * dble(vwfn%nband + pol%ncrit) * dble(cwfn%nband) * 8
      endif
      if (pol%freq_dep .eq. 2 .or. pol%freq_dep .eq. 3) then
        rmem = rmem + dble(kp%nspin) * dble(peinf%nvownmax) * dble(peinf%ncownmax) * dble(nrqmax) * dble(pol%nfreq_group)*8
      endif
! array gvec%index_vec in input
      rmem = rmem + dble(gvec%nFFTgridpts) * 4
! arrays fftbox1 and fftbox2 in subroutines mtxel
      call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)
      rmem = rmem + dble(product(Nfft(1:3))) * 32
! arrays cwfn%wfn_fft in subroutines mtxel
      if (pol%os_opt_ffts.ge.1) then
        rmem = rmem + dble(product(Nfft(1:3))) * dble(peinf%ncownmax) * 16
      endif
      if (pol%os_opt_ffts.eq.2) then
        rmem = rmem + dble(product(Nfft(1:3))) * dble(peinf%nvownmax) * 16
      endif
      write(6,989) rmem / 1024**2
      write(7,989) rmem / 1024**2

! store required memory in variable smem
      smem=rmem

! random numbers
      rmem=0.0D0
! various truncation schemes
      rmem2=0.0d0
! cell wire truncation
      if (pol%icutv==TRUNC_WIRE) then
        dkmax(1) = gvec%FFTgrid(1) * n_in_wire
        dkmax(2) = gvec%FFTgrid(2) * n_in_wire
        dkmax(3) = 1
        call setup_FFT_sizes(dkmax,dNfft,dscale)
! array fftbox_2D
        rmem2 = rmem2 + dble(dNfft(1) * dNfft(2)) * 16
! array inv_indx
        rmem2 = rmem2 + dble(product(Nfft(1:3))) * 4
! array qran
        rmem2 = rmem2 + 3 * dble(nmc) * 8
      endif
! cell box truncation (parallel version only)
      if (pol%icutv==TRUNC_BOX) then
        dkmax(1:3) = gvec%FFTgrid(1:3) * n_in_box
        call setup_FFT_sizes(dkmax,dNfft,dscale)
        if (mod(dNfft(3),peinf%npes) == 0) then
          Nplane = dNfft(3)/peinf%npes
        else
          Nplane = dNfft(3)/peinf%npes+1
        endif
        if (mod(dNfft(1)*dNfft(2),peinf%npes) == 0) then
          Nrod = (dNfft(1)*dNfft(2))/peinf%npes
        else
          Nrod = (dNfft(1)*dNfft(2))/peinf%npes+1
        endif
! array fftbox_2D
        rmem2 = rmem2 + dble(dNfft(1)) * dble(dNfft(2)) * dble(Nplane) * 16
! array fftbox_1D
        rmem2 = rmem2 + dble(dNfft(3)) * dble(Nrod) * 16
! array dummy
!          rmem2=rmem2+dble(dNfft(1))*dble(dNfft(2))*16.0d0
! arrays dummy1 and dummy2
        rmem2 = rmem2 + dble(Nrod) * (peinf%npes + 1) * 16
! array inv_indx
        rmem2 = rmem2 + dble(product(Nfft(1:3))) * 4
      endif
      if (rmem2 .gt. rmem) rmem = rmem2
      write(6,988) rmem / 1024**2
      write(7,988) rmem / 1024**2
      if (smem .gt. mem) write(0,777)
    endif
989 format(1x,'Memory required for execution: ',f0.1,' MB per PE')
988 format(1x,'Memory required for vcoul: ',f0.1,' MB per PE',/)
777 format(1x,'WARNING: Required memory greater than memory available.',/, &
      3x,'This calculation will likely die. To reduce memory,',/, &
      3x,'first make sure you are running on an ideal # of PEs.',/, &
      3x,'Then, try increasing the number of PEs or try adding',/, &
      3x,'gcomm_elements to epsilon.inp',/)

!--------------------------------------------
! SIB: Read in crys%bdot and crys%celvol from file #25
! Complain if crys%celvol and (2pi)^3/sqrt[|det(b)|] are different
! bdot is read in row by row.  get_volume() assumes that bdot is symmetric.
! bdot should be dot products of the reciprocal lattice vectors.
! Therefore vol is the real space volume of the cell.

    if(peinf%inode.eq.0) then

      write(7,'(1x,"Cell Volume =",e16.9,/)') crys%celvol

!-------------------------
! Compute Cell Volume

      call get_volume(vcell,crys%bdot)
      if (abs(crys%celvol-vcell).gt.TOL_Small) then
        call die('volume mismatch')
      endif

      if (peinf%verb_medium) then
        write(6,'()')
        write(6,'(1x,a)') 'Symmetry information:'
        write(6,'(1x,a,i0)') '- Number of symmetry operations: ', syms%ntran
        write(6,'(1x,a)') '- Symmetry matrices and fractional translations:'
        do itran=1,syms%ntran
          write(6,'(3x,a,i0,a,3(1x,f10.6))') &
            'Matrix # ', itran, ' , frac. transl.:',syms%tnp(:,itran)
          write(6,'(1(2x,3(1x,i2)))') (syms%mtrx(ii,:,itran), ii=1,3)
        enddo
        write(6,'()')
      endif
    endif

    if(cwfn%nband.gt.kp%mnband) then
      if(peinf%inode == 0) then
        write(tmpstr1,660) cwfn%nband
        write(tmpstr2,660) kp%mnband
        write(tmpstr,666) TRUNC(tmpstr1), TRUNC(tmpstr2)
      endif
660   format(i16)
666   format(1x,'The number of bands (',a,&
        ') specified in epsilon.inp is larger than the number of bands (',a,') available in WFN.')
      call die(tmpstr, only_root_writes = .true.)
    endif
    !if(cwfn%nband.eq.kp%mnband) then WK
    !  call die("You must provide one more band in WFN than used in epsilon.inp in order to assess degeneracy.")
    !endif


!----------------------------------------------------------------
! If quasi-particle corrections were requested, read the corrected
! quasiparticle energies from file (in eV)

    if(peinf%inode == 0) then
      ! WK
      !if(any(abs(kp%el(cwfn%nband+vwfn%ncore_excl, 1:kp%nrk, 1:kp%nspin) - & 
      !           kp%el(cwfn%nband + 1 +vwfn%ncore_excl, 1:kp%nrk, 1:kp%nspin)) .lt. TOL_Degeneracy)) then
      !  if(pol%degeneracy_check_override) then
      !    write(0,'(a)') &
      !      "WARNING: Selected number of bands breaks degenerate subspace. " // &
      !      "Run degeneracy_check.x for allowable numbers."
      !    write(0,*)
      !  else
      !    write(0,'(a)') &
      !      "Run degeneracy_check.x for allowable numbers, or use keyword " // &
      !      "degeneracy_check_override to run anyway (at your peril!)."
      !    call die("Selected number of bands breaks degenerate subspace.")
      !  endif
      !endif

      if(any (kp%ifmax(:,:) < vwfn%nband+vwfn%ncore_excl .or. kp%ifmax(:,:) > vwfn%nband+vwfn%ncore_excl + pol%ncrit)) then
        write(0,'(a,i6,a,i6,a)') 'epsilon.inp says there are ', vwfn%nband, ' fully occupied bands and ', &
          pol%ncrit, ' partially occupied.'
        write(0,'(a,2i6)') 'This is inconsistent with highest bands in WFN file; min, max = ', minval(kp%ifmax), maxval(kp%ifmax)
        call die("band_occupation, number_partial_occup, and WFN inconsistent.")
      endif

      if(maxval(kp%ifmax) - minval(kp%ifmax) > pol%ncrit) then
        write(0,'(a,i6,a)') 'epsilon.inp says there are ', pol%ncrit, ' partially occupied bands.'
        write(0,'(a,i6)') 'This is less than the number partially occupied in WFN file: ', maxval(kp%ifmax) - minval(kp%ifmax)
        call die("number_partial_occup and WFN inconsistent.")
      endif

      call assess_degeneracies(kp, kp%el(cwfn%nband + 1+vwfn%ncore_excl, :, :), cwfn%nband, pol%efermi, TOL_Degeneracy, &
                               ncore_excl=vwfn%ncore_excl)

      call calc_qtot(kp, crys%celvol, pol%efermi, qtot, omega_plasma, write7 = .true.)

      if (mod(peinf%npes,npools) .ne.0) then
        write(0,'(/1x,a)') 'WARNING: The number of cpus does not divide evenly in the optimal number of pools.'
        write(0,'(1x,i0,a/)') mod(peinf%npes,npools), 'cpus are doing no work'
      endif
      if (peinf%nvownactual .ne. peinf%nvownmax) then
        write(0,'(/1x,a)') ' WARNING: Your valence bands are not equally distributed among the pools'
        write(0,'(1x,a,i0)') 'Max valence bands per pool is', peinf%nvownmax
        write(0,'(1x,a,i0/)') 'Min valence bands per pool is', peinf%nvownactual
      endif
    endif

#ifdef MPI
    call MPI_BCAST(omega_plasma, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
    if(peinf%inode.ne.0) then
      SAFE_ALLOCATE(kp%degeneracy, (cwfn%nband, kp%nrk, kp%nspin))
    endif
    call MPI_Bcast(kp%degeneracy(1,1,1), cwfn%nband * kp%nrk * kp%nspin, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif

    call gvec_index(gvec)

    call timing%start(timing%io_total)
    if (.not. pol%skip_chi) then
      if (pol%wfn_hdf5) then
#ifdef HDF5
        call read_hdf5_wavefunctions(kp, gvec, pol, cwfn, vwfn, intwfnv, intwfnc)
#endif
      else
        call read_wavefunctions(kp, gvec, pol, cwfn, vwfn, intwfnv, intwfnc)
        if (peinf%inode.eq.0) then
          call close_file(25)
        endif
      endif
    endif
    call timing%stop(timing%io_total)

    if (pol%eqp_corrections) call scissors_zero(pol%scis)

    if(peinf%inode == 0) then
      write(aheadinput,'(60x)')
      write(ajname,'("chiGG0")')
      call date_time(adate,atime)
    endif

    call eps_setup_sizes(pol, SCALARSIZE, kp%nspin)

    if(pol%subspace) then
      ! just allocate this for all
      SAFE_ALLOCATE(pol%neigen_of_q,(pol%nq))
      pol%neigen_of_q = 0
    end if

    !XXXX
    ! make sure we have a derfined size for the subspace basis
#ifdef HDF5
    if (pol%use_hdf5) then
      if(pol%subspace .and. pol%matrix_in_subspace_basis) then
        if(pol%neig_sub_input < 0) then
          if(peinf%inode.eq.0) then
            pol%neig_sub_input = INT(nmtx * 0.20d0)
          end if
#ifdef MPI
            call MPI_BCAST(pol%neig_sub_input, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
        end if
      end if
    end if
#endif
    !XXXX

    call timing%start(timing%io_total)

    if (peinf%inode.eq.0) then

      if (pol%nq0>0) then

#ifdef HDF5
        if (pol%use_hdf5) then
          if (pol%skip_epsilon) then
            filenameh5 = 'chi0mat.h5'
          else
            filenameh5 = 'eps0mat.h5'
          endif
          call eps_hdf5_setup(kp, gvec, syms, crys, pol, pol%qgrid, &
            pol%nq0, pol%qpt(:,1:pol%nq0), pol%nmtx_of_q(1:pol%nq0), nmtx0, &
            cwfn%nband, trim(filenameh5), restart=pol%restart)
        else
#endif
!-------------------------------------------
! Initialize chi0mat and eps0mat:
! SIB: only printed if pol%nq0>0
! All sorts of stuff is printed, such as the epsilon cutoff, number of bands,
! gvectors, bdot matrix, the addresses (gvec%index_vec), and pol%qpt (q-points)
! to chi0mat, and eps0mat has some of it.

          if (pol%skip_epsilon) then
            write(6,'(a)') 'File chi0mat will be written.'
            call open_file(10,file='chi0mat',form='unformatted',status='replace')
            write(10) aheadinput,ajname,adate
            write(10) pol%freq_dep,pol%nFreq
            write(10) pol%qgrid(1:3)
            if (pol%freq_dep .ne. 0) then
              write(10) (pol%dFreqGrid(ii),ii=1,pol%nFreq),(pol%dFreqBrd(ii),ii=1,pol%nFreq)
            else
              write(10)
            endif
            write(10)
            write(10)
            write(10) pol%ecuts,cwfn%nband
            write(10) kp%nrk, 1 ! invflag
            write(10) gvec%ng,gvec%nFFTgridpts,gvec%FFTgrid(1:3), &
              ((gvec%components(ii,ig),ii=1,3),ig=1,gvec%ng), &
              ((crys%bdot(ii,jj),jj=1,3),ii=1,3),(gvec%index_vec(ig),ig=1,gvec%nFFTgridpts)
            write(10) ((pol%qpt(ii,iq),ii=1,3),iq=1,pol%nq0)
          else
            call open_file(12,file='eps0mat',form='unformatted',status='replace')
            write(12) ajname,adate
            write(12) pol%freq_dep,pol%nFreq, pol%subspace, pol%matrix_in_subspace_basis, pol%keep_full_eps_static
            write(12) pol%qgrid(1:3)
            if (pol%freq_dep .ne. 0) then
              write(12) (pol%dFreqGrid(ii),ii=1,pol%nFreq),(pol%dFreqBrd(ii),ii=1,pol%nFreq)
            else
              write(12)
            endif
            write(12)
            write(12)
            write(12) pol%ecuts
            write(12) pol%nq0,((pol%qpt(ii,iq),ii=1,3),iq=1,pol%nq0)
            write(12) gvec%ng, ((gvec%components(jj,ig),jj=1,3),ig=1,gvec%ng)
          endif
#ifdef HDF5
        endif
#endif

      endif

!----------------------------------------
! Initialize chimat and epsmat:
! SIB: chimat and epsmat are printed to for all non-zero q vectors.
! Same sort of information as above.

      if (pol%nq1>0) then

#ifdef HDF5
        if (pol%use_hdf5) then
          if (pol%skip_epsilon) then
            filenameh5 = 'chimat.h5'
          else
            filenameh5 = 'epsmat.h5'
          endif
          call eps_hdf5_setup(kp, gvec, syms, crys, pol, pol%qgrid, &
            pol%nq1, pol%qpt(:,pol%nq0+1:pol%nq), pol%nmtx_of_q(pol%nq0+1:pol%nq), nmtx, &
            cwfn%nband, trim(filenameh5), restart=pol%restart)
        else
#endif
          if (pol%skip_epsilon) then
            write(6,'(a)') 'File chimat will be written.'
            call open_file(11,file='chimat',form='unformatted',status='replace')
            write(11) aheadinput,ajname,adate
            write(11) pol%freq_dep,pol%nFreq
            write(11) pol%qgrid(1:3)
            if (pol%freq_dep .ne. 0) then
              write(11) (pol%dFreqGrid(ii),ii=1,pol%nFreq),(pol%dFreqBrd(ii),ii=1,pol%nFreq)
            else
              write(11)
            endif
            write(11)
            write(11)
            write(11) pol%ecuts,cwfn%nband
            write(11) kp%nrk, 1 ! invflag
            write(11) gvec%ng,gvec%nFFTgridpts,gvec%FFTgrid(1:3), &
              ((gvec%components(jj,ig),jj=1,3),ig=1,gvec%ng), &
              ((crys%bdot(ii,jj),jj=1,3),ii=1,3),(gvec%index_vec(ig),ig=1,gvec%nFFTgridpts)
            write(11) pol%nq,pol%nq0,((pol%qpt(jj,iq),jj=1,3),iq=1,pol%nq) !??????
          else
            call open_file(13,file='epsmat',form='unformatted',status='replace')
            write(13) ajname,adate
            write(13) pol%freq_dep,pol%nFreq, pol%subspace, pol%matrix_in_subspace_basis, pol%keep_full_eps_static
            write(13) pol%qgrid(1:3)
            if (pol%freq_dep .ne. 0) then
              write(13) (pol%dFreqGrid(ii),ii=1,pol%nFreq),(pol%dFreqBrd(ii),ii=1,pol%nFreq)
            else
              write(13)
            endif
            write(13)
            write(13)
            write(13) pol%ecuts
            write(13) pol%nq1, ((pol%qpt(jj,iq),jj=1,3),iq=pol%nq0+1,pol%nq)
            write(13) gvec%ng, ((gvec%components(jj,ig),jj=1,3),ig=1,gvec%ng)
          endif
#ifdef HDF5
        endif
#endif

      endif

!-------------------------------
! SIB: print to stdout and unit=7 (epsilon.log) some information:  cutoffs,
! number of bands, number of q points, and list of q-points.
      write(6,120) pol%ecuts,cwfn%nband,pol%nq
      write(7,120) pol%ecuts,cwfn%nband,pol%nq
120   format(1x,'- Screened Coulomb cutoff: ',f0.3,' Ry'/,&
        1x,'- Total number of bands: ',i0,/,&
        1x,'- Number of q-points: ',i0)

! SIB: report on scissors operator parameters to units 6 and 7
      call scissors_write(6, pol%scis)
      call scissors_write(7, pol%scis)

      write(6,'()')
      do iq=1,pol%nq
        norm=sqrt(DOT_PRODUCT(MATMUL(crys%bdot, pol%qpt(:,iq)),pol%qpt(:,iq)))
        if (iq<=pol%nq0) then
          write(6,925) (pol%qpt(jj,iq),jj=1,3),norm
        else
          write(6,926) (pol%qpt(jj,iq),jj=1,3),norm
        endif
      enddo
      write(6,'()')
925   format(1x,'Q0 and |Q0| =',4f10.6)
926   format(1x,'Q  and |Q|  =',4f10.6)
    endif

    call timing%stop(timing%io_total)

#ifdef MPI
    ! FHJ: the pol%restart flag might have been changed by rank 0; need to bcast
    call MPI_Bcast(pol%restart, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
#endif

    POP_SUB(input)

    return

  contains

    ! Distribute Val/Cond Bands over Processors
    subroutine distribution()

      PUSH_SUB(input.distribution)

      ! Create pools if not read from epsilon.inp

      if (peinf%npools .le. 0 .or. peinf%npools .gt. peinf%npes) then
        call createpools(vwfn%nband+pol%ncrit,cwfn%nband-vwfn%nband,peinf%npes,npools,peinf%nvownmax,peinf%ncownmax)
        peinf%npools = npools
      else
        npools = peinf%npools

        if (mod((vwfn%nband+pol%ncrit),npools) .eq. 0) then
          peinf%nvownmax = (vwfn%nband+pol%ncrit) / npools
        else
          peinf%nvownmax = ((vwfn%nband+pol%ncrit) / npools) + 1
        endif

        if (mod((cwfn%nband-vwfn%nband),(peinf%npes/npools)) .eq. 0) then
          peinf%ncownmax = (cwfn%nband-vwfn%nband) / (peinf%npes / npools)
        else
          peinf%ncownmax = (cwfn%nband-vwfn%nband) / (peinf%npes / npools) + 1
        endif
      endif

      if (peinf%inode .eq. 0) then
        write(tmpstr1,440) npools
        write(tmpstr2,440) peinf%ncownmax
        write(tmpstr3,440) peinf%nvownmax
        write(6,444) TRUNC(tmpstr1),TRUNC(tmpstr2),TRUNC(tmpstr3)
        write(7,444) TRUNC(tmpstr1),TRUNC(tmpstr2),TRUNC(tmpstr3)
      endif
440   format(i16)
444   format(1x,"Running with",1x,a,1x,"valence pools",/, &
        1x,"Number of conduction bands per processor:",1x,a,/, &
        1x,"Number of valence bands per processor:",1x,a,/)

      SAFE_ALLOCATE(peinf%global_pairowner,((vwfn%nband+pol%ncrit),(cwfn%nband-vwfn%nband)))
      peinf%global_pairowner=0
      SAFE_ALLOCATE(peinf%doiownv,(vwfn%nband+pol%ncrit))
      peinf%doiownv=.false.
      SAFE_ALLOCATE(peinf%doiownc,(cwfn%nband-vwfn%nband))
      peinf%doiownc=.false.
      SAFE_ALLOCATE(peinf%does_it_ownv,(vwfn%nband+pol%ncrit,peinf%npes))
      peinf%does_it_ownv=.false.
      SAFE_ALLOCATE(peinf%does_it_ownc,(cwfn%nband-vwfn%nband,peinf%npes))
      peinf%does_it_ownc=.false.
      SAFE_ALLOCATE(peinf%global_nvown,(peinf%npes))
      peinf%global_nvown=0
      SAFE_ALLOCATE(peinf%global_ncown,(peinf%npes))
      peinf%global_ncown=0
      SAFE_ALLOCATE(peinf%indexv,(vwfn%nband+pol%ncrit))
      peinf%indexv=0
      SAFE_ALLOCATE(peinf%global_indexv,(vwfn%nband+pol%ncrit,peinf%npes))
      peinf%global_indexv=0
      SAFE_ALLOCATE(peinf%indexc,(cwfn%nband-vwfn%nband))
      peinf%indexc=0
      SAFE_ALLOCATE(peinf%invindexv,(peinf%nvownmax))
      peinf%invindexv=0
      SAFE_ALLOCATE(peinf%invindexc,(peinf%ncownmax))
      peinf%invindexc=0

      peinf%nvownactual=0
      peinf%ncownactual=0

      ! FHJ: each pool consist of (peinf%npes/npools) (rounded down) sequential
      !      MPI processes.
      !      Except for iv and ic all indices here start at 0.
      mypool = (peinf%inode/(peinf%npes/npools))
      mypoolrank = mod(peinf%inode,(peinf%npes/npools)) ! my rank within the pool
      myipe = peinf%inode + 1

      do iv = 1,vwfn%nband+pol%ncrit
        ! FHJ: The valence bands are assigned to the pools in a sequential fashion.
        !      This is more efficient for the parallel IO.  Eg:
        !       for 2 (val. bands)/pool, v1->p0, v2->p0, v3->p1, etc.
        ipool = (iv-1)/peinf%nvownmax

        if (mypool .eq. ipool .and. peinf%inode .lt. peinf%npes) then

          peinf%nvownactual=peinf%nvownactual+1
          peinf%global_nvown(myipe)=peinf%nvownactual
          peinf%indexv(iv)=peinf%nvownactual
          peinf%global_indexv(iv,myipe)=peinf%indexv(iv)
          peinf%invindexv(peinf%nvownactual)=iv
          peinf%doiownv(iv)=.true.

          do ic = 1, cwfn%nband-vwfn%nband
            ! FHJ: Conduction bands are also assigned sequentially. Eg:
            !       c1..c4 -> pr0, c5..c8 -> pr1, etc.
            if ( (ic-1)/peinf%ncownmax == mypoolrank ) then
              ! We only have to create the local list of conduction bands once
              if (peinf%nvownactual .eq. 1) then
                peinf%ncownactual=peinf%ncownactual+1
                peinf%global_ncown(myipe) = peinf%ncownactual
                peinf%invindexc(peinf%global_ncown(myipe))=ic
                peinf%indexc(ic)=peinf%ncownactual
                peinf%doiownc(ic)=.true.
              endif
              peinf%global_pairowner(iv,ic)=myipe
            endif
          enddo

        endif
      enddo

      if (peinf%inode .lt. peinf%npes) then
        peinf%does_it_ownv(:,peinf%inode+1) = peinf%doiownv(:)
        peinf%does_it_ownc(:,peinf%inode+1) = peinf%doiownc(:)
      endif
#ifdef MPI
      do ipe=0, peinf%npes-1
        call MPI_Bcast(peinf%does_it_ownv(1,ipe+1), vwfn%nband+pol%ncrit,&
          MPI_LOGICAL, ipe, MPI_COMM_WORLD, mpierr)
        call MPI_Bcast(peinf%does_it_ownc(1,ipe+1), cwfn%nband-vwfn%nband,&
          MPI_LOGICAL, ipe, MPI_COMM_WORLD, mpierr)
      enddo

      SAFE_ALLOCATE(global_pairowner_temp,((vwfn%nband+pol%ncrit),(cwfn%nband-vwfn%nband)))
      call MPI_ALLREDUCE(peinf%global_pairowner(1,1),global_pairowner_temp(1,1), &
        (vwfn%nband+pol%ncrit)*(cwfn%nband-vwfn%nband),MPI_INTEGER,MPI_SUM, &
        MPI_COMM_WORLD,mpierr)
      peinf%global_pairowner=global_pairowner_temp
      SAFE_DEALLOCATE(global_pairowner_temp)

      SAFE_ALLOCATE(global_nvown_temp,(peinf%npes))
      call MPI_ALLREDUCE(peinf%global_nvown,global_nvown_temp,peinf%npes,MPI_INTEGER,MPI_SUM, &
        MPI_COMM_WORLD,mpierr)
      peinf%global_nvown=global_nvown_temp
      SAFE_DEALLOCATE(global_nvown_temp)

      SAFE_ALLOCATE(global_ncown_temp,(peinf%npes))
      call MPI_ALLREDUCE(peinf%global_ncown,global_ncown_temp,peinf%npes,MPI_INTEGER,MPI_SUM, &
        MPI_COMM_WORLD,mpierr)
      peinf%global_ncown=global_ncown_temp
      SAFE_DEALLOCATE(global_ncown_temp)

! The below command works for MPI-1 and MPI-2, so we will use it for now. Soon we will have the in-place command
! that is commented out below to use, since it is more efficient with memory and most machines can run with MPI-2 nowadays

      SAFE_ALLOCATE(global_indexv_temp,((vwfn%nband+pol%ncrit),peinf%npes))
      call MPI_ALLREDUCE(peinf%global_indexv(1,1),global_indexv_temp(1,1),peinf%npes*(vwfn%nband+pol%ncrit),MPI_INTEGER,MPI_SUM, &
        MPI_COMM_WORLD,mpierr)
      peinf%global_indexv=global_indexv_temp
      SAFE_DEALLOCATE(global_indexv_temp)

! Here is the command for MPI-2. It will be commented back in when we have appropriate compiler flags set up. All other
! mpi allreduces should be done in the same way once things are set up
!  call MPI_ALLREDUCE(MPI_IN_PLACE,peinf%global_indexv(1,1),peinf%npes*(vwfn%nband+pol%ncrit),MPI_INTEGER,MPI_SUM, &
!    MPI_COMM_WORLD,mpierr)

#endif

      POP_SUB(input.distribution)

    end subroutine distribution

  end subroutine input

#ifdef HDF5
!================================================================================
!
! read_hdf5_wavefunctions()   Originally By JIM                 Last Modified 5/1/2012
!
! Read wavefunctions from WFN.h5 and distribute to processors
! depending on the bands that they own
!
!===============================================================================
  subroutine read_hdf5_wavefunctions(kp, gvec, pol, cwfn, vwfn, intwfnv, intwfnc)
    type (kpoints), intent(in) :: kp
    type (gspace), intent(in) :: gvec
    type (polarizability), intent(in) :: pol
    type (conduction_wfns), intent(in) :: cwfn
    type (valence_wfns), intent(in) :: vwfn
    type (int_wavefunction), intent(out) :: intwfnv
    type (int_wavefunction), intent(out) :: intwfnc

    integer, allocatable :: isort(:)
    integer, allocatable :: components(:,:)
    SCALAR, allocatable :: wfns(:,:,:)

    integer :: ii,ig,ik,iiii,ib,is
    real(DP) :: qk(3)
    type(gspace) :: gvec_kpt

    integer(HID_T) :: file_id
    integer(HID_T) :: plist_id
    integer :: error
    integer :: ngktot
    integer :: istart, ib_first

    PUSH_SUB(read_hdf5_wavefunctions)

    call logit('Inside read hdf5 wfns')
    SAFE_ALLOCATE(intwfnv%ng, (kp%nrk))
    SAFE_ALLOCATE(intwfnv%isort, (kp%ngkmax,kp%nrk))
    SAFE_ALLOCATE(intwfnv%cg, (kp%ngkmax,kp%nrk*peinf%nvownactual,kp%nspin*kp%nspinor))
    SAFE_ALLOCATE(intwfnv%qk, (3,kp%nrk))
    SAFE_ALLOCATE(intwfnc%ng, (kp%nrk))
    SAFE_ALLOCATE(intwfnc%isort, (kp%ngkmax,kp%nrk))
    SAFE_ALLOCATE(intwfnc%cg, (kp%ngkmax,kp%nrk*peinf%ncownactual,kp%nspin*kp%nspinor))
    SAFE_ALLOCATE(intwfnc%cbi, (kp%nrk*peinf%ncownactual))
    SAFE_ALLOCATE(intwfnc%qk, (3,kp%nrk))

    ngktot = SUM(kp%ngk)
    SAFE_ALLOCATE(gvec_kpt%components, (3, ngktot))

    call logit('Reading G-vectors')
    call read_hdf5_wfn_gvectors('WFN.h5', gvec_kpt%components, ngktot)

    istart=1
    do ik= 1, kp%nrk
      qk(:)=kp%rk(:,ik)

      SAFE_ALLOCATE(isort, (kp%ngk(ik)))
      SAFE_ALLOCATE(components, (3,kp%ngk(ik)))

      components(:,:) = gvec_kpt%components(:,istart:istart+kp%ngk(ik)-1)
      istart=istart+kp%ngk(ik)

      do ig = 1, kp%ngk(ik)
        call findvector(isort(ig), components(:, ig), gvec)
        if (isort(ig) ==  0)  call die('input: could not find gvec')
      enddo

      intwfnv%ng(ik)=kp%ngk(ik)
      intwfnv%isort(1:kp%ngk(ik),ik)=isort(1:kp%ngk(ik))
      intwfnv%qk(:,ik)=qk(:)

      intwfnc%ng(ik)=kp%ngk(ik)
      intwfnc%isort(1:kp%ngk(ik),ik)=isort(1:kp%ngk(ik))
      intwfnc%qk(:,ik)=qk(:)

      SAFE_DEALLOCATE(isort)
      SAFE_DEALLOCATE(components)
    end do

    SAFE_DEALLOCATE_P(gvec_kpt%components)

    call logit("Before interface setup")
    ! HDF5 interface setup
    call hdf5_open_file('WFN.h5', 'r', file_id, parallel_io=.true.)

    call logit('Reading valence bands')
    ib_first = 0
    if (peinf%nvownactual>0) ib_first = peinf%invindexv(1)
    SAFE_ALLOCATE(wfns, (ngktot,kp%nspin*kp%nspinor,peinf%nvownactual))
    call read_hdf5_bands_block(file_id, kp, peinf%nvownmax, peinf%nvownactual, &
      peinf%does_it_ownv, ib_first, wfns, ioffset=vwfn%ncore_excl)

! DVF : here we flip from hdf5 wfn ordering of indices to the `traditional` BGW ordering of indices, with
! respect to spin. The traditional BGW ordering should change soon to reflect the newer, better hdf5 setup
! Note that the sum for switching the indices is over kp%nspin*kp%nspinor, while the sum for check norm
! is over kp%nspin but you pass kp%nspinor. This is because the norm requires special consideration when you
! have spinors. So, you can`t have a simple sum over kp%nspin*kp%nspinor.

    call logit('Checking norms')
    ! perform checks and write to valence file
    do ib = 1, peinf%nvownactual
      istart=1
      do ik = 1, kp%nrk
        if (peinf%verb_debug .and. peinf%inode==0) then
          write(6,'(1x,a,i0,", ",i0)') "val: ib, ik = ", ib, ik
        endif
        iiii=peinf%indexv(peinf%invindexv(ib))+(ik-1)*peinf%nvownactual
        do is=1,kp%nspin*kp%nspinor
          intwfnv%cg(1:kp%ngk(ik),iiii,is)=wfns(istart:istart+kp%ngk(ik)-1,is,ib)
        enddo
        do is = 1, kp%nspin
          call checknorm('WFN',peinf%invindexv(ib),ik,kp%ngk(ik),is,kp%nspinor,intwfnv%cg(1:kp%ngk(ik),iiii,:))
        enddo
        ! for metals
        ! FHJ: FIXME: isn`t the following block redundant?
        if(peinf%invindexv(ib).le.vwfn%nband+pol%ncrit) then
          iiii=peinf%indexv(peinf%invindexv(ib))+(ik-1)*peinf%nvownactual
          do is=1,kp%nspin*kp%nspinor
            intwfnv%cg(1:kp%ngk(ik),iiii,is)=wfns(istart:istart+kp%ngk(ik)-1,is,ib)
          enddo
        end if
        istart=istart+kp%ngk(ik)
      enddo
    enddo
    SAFE_DEALLOCATE(wfns)

    call logit('Reading conduction bands')
    ib_first = 0
    if (peinf%ncownactual>0) ib_first = peinf%invindexc(1)
    SAFE_ALLOCATE(wfns, (ngktot,kp%nspin*kp%nspinor,peinf%ncownactual))
    call read_hdf5_bands_block(file_id, kp, peinf%ncownmax, peinf%ncownactual, &
      peinf%does_it_ownc, ib_first, wfns, ioffset = vwfn%nband+vwfn%ncore_excl)

    call logit('Checking norms')
    ! write to conduction file
    do ib = 1, peinf%ncownactual
      istart=1
      do ik = 1, kp%nrk
        if (peinf%verb_debug .and. peinf%inode==0) then
          write(6,'(1x,a,i0,", ",i0)') "cond: ib, ik = ", ib, ik
        endif
        iiii=ib+(ik-1)*peinf%ncownactual
        intwfnc%cbi(iiii)=peinf%invindexc(ib)
        do is=1,kp%nspin*kp%nspinor
          intwfnc%cg(1:kp%ngk(ik),iiii,is)=wfns(istart:istart+kp%ngk(ik)-1,is,ib)
        enddo
        do is = 1, kp%nspin
          call checknorm('WFN',peinf%invindexc(ib)+vwfn%nband,ik,kp%ngk(ik),is,kp%nspinor,intwfnc%cg(1:kp%ngk(ik),iiii,:))
        enddo
        istart=istart+kp%ngk(ik)
      enddo
    enddo
    SAFE_DEALLOCATE(wfns)

    call hdf5_close_file(file_id)
    call logit('Finished reading HDF5 wavefunctions')

    POP_SUB(read_hdf5_wavefunctions)
    return
  end subroutine read_hdf5_wavefunctions
#endif

  subroutine read_wavefunctions(kp, gvec, pol, cwfn, vwfn, intwfnv, intwfnc)
    type (kpoints), intent(in) :: kp
    type (gspace), intent(in) :: gvec
    type (polarizability), intent(in) :: pol
    type (conduction_wfns), intent(in) :: cwfn
    type (valence_wfns), intent(in) :: vwfn
    type (int_wavefunction), intent(out) :: intwfnv
    type (int_wavefunction), intent(out) :: intwfnc

    integer, allocatable :: isort(:)
    SCALAR, allocatable :: zc(:,:)

    character :: filename*20
    character :: filenamev*20
    integer :: ii,ig,ik,iiii,ib,ib2,is
    integer :: dest
    integer :: iwritecb
    integer :: iunit_v,iunit_c
    real(DP) :: qk(3)
    logical :: dont_read

    type(gspace) :: gvec_kpt
    type(progress_info) :: prog_info !< a user-friendly progress report

    PUSH_SUB(read_wavefunctions)

!------------------------------------
! Loop over kpoints and read in wavefunctions

    SAFE_ALLOCATE(intwfnv%ng, (kp%nrk))
    SAFE_ALLOCATE(intwfnv%isort, (kp%ngkmax,kp%nrk))
    SAFE_ALLOCATE(intwfnv%cg, (kp%ngkmax,kp%nrk*peinf%nvownactual,kp%nspin*kp%nspinor))
    SAFE_ALLOCATE(intwfnv%qk, (3,kp%nrk))
    SAFE_ALLOCATE(intwfnc%ng, (kp%nrk))
    SAFE_ALLOCATE(intwfnc%isort, (kp%ngkmax,kp%nrk))
    SAFE_ALLOCATE(intwfnc%cg, (kp%ngkmax,kp%nrk*peinf%ncownactual,kp%nspin*kp%nspinor))
    SAFE_ALLOCATE(intwfnc%cbi, (kp%nrk*peinf%ncownactual))
    SAFE_ALLOCATE(intwfnc%qk, (3,kp%nrk))

!----------------------------------------------------------------------------
! Beginning k-point loop

    call progress_init(prog_info, 'reading wavefunctions (WFN)', 'state', kp%nrk*cwfn%nband)
    do ik=1,kp%nrk
      qk(:)=kp%rk(:,ik)

! For each G-vector read, tries to find its match in gvec%components(:,:).
! If not found, aborts.  Otherwise stores its index in isort(i)
! where ik is the current k-point we are considering (loop we are in)
! and i is the index of the kx,ky,kz just read.

      SAFE_ALLOCATE(gvec_kpt%components, (3, kp%ngk(ik)))

      call read_binary_gvectors(25, kp%ngk(ik), kp%ngk(ik), gvec_kpt%components)

      SAFE_ALLOCATE(isort, (kp%ngk(ik)))
      do ig = 1, kp%ngk(ik)
        call findvector(isort(ig), gvec_kpt%components(:, ig), gvec)
        if (isort(ig) ==  0)  call die('input: could not find gvec')
      enddo
      SAFE_DEALLOCATE_P(gvec_kpt%components)

      intwfnv%ng(ik)=kp%ngk(ik)
      intwfnv%isort(1:kp%ngk(ik),ik)=isort(1:kp%ngk(ik))
      intwfnv%qk(:,ik)=qk(:)

! SIB: each proc writes to unit iunit_c qk, kp%el, kp%ngk(ik), and isort
! for current kpoint (ik) and all bands.

      intwfnc%ng(ik)=kp%ngk(ik)
      intwfnc%isort(1:kp%ngk(ik),ik)=isort(1:kp%ngk(ik))
      intwfnc%qk(:,ik)=qk(:)

! SIB: loop over i for kp%mnbands (number of bands) and do the following:
!   - reads the info from unit=25 and checks norm of wavefunction.
!   - if a valence band, proc 0 writes it to unit iunit_v
!   - if a conduction band, then it is sent to/received by the correct
!     destination node (which is proc # i-nvalence-1 mod npes)
!   - the destination processor then writes the data to unit iunit_c

      SAFE_ALLOCATE(zc, (kp%ngk(ik), kp%nspinor*kp%nspin))
      do ib=1,kp%mnband

! JRD: Dumb debugging
        if (peinf%verb_debug .and. peinf%inode==0) then
          write(6,'("Reading Wavefunction: ik = ",i6," n = ",i6)') ik, ib
        endif

        dont_read = (ib > cwfn%nband + vwfn%ncore_excl .or. ib <= vwfn%ncore_excl)
        call read_binary_data(25, kp%ngk(ik), kp%ngk(ik), kp%nspin*kp%nspinor, zc, dont_read = dont_read)

        ! FHJ: the following lines were introduced in r6294 and are supposed to
        ! be a shortcut if we are past the last band of the last k-point. However,
        ! in light of a previous bug (#223), this feature is commented out for now.
        !! FHJ: shortcut if this is past the last band of the last k-point
        !if (dont_read .and. ik==kp%nrk) exit
        !if (dont_read) cycle
        if(ib > cwfn%nband + vwfn%ncore_excl .or. ib <= vwfn%ncore_excl) cycle
        call progress_step(prog_info)

        if(peinf%inode == 0) then
          do is = 1, kp%nspin
            call checknorm('WFN',ib,ik,kp%ngk(ik),is,kp%nspinor,zc(:,:))
          enddo
        endif

! DVF: recall that we redefined the number of valence bands to exclude the
! core states. So, we have to subtract ncore_excl right here because ib is
! referenced to the full wavefunction file including all the core states.
! The arrays for the val/cond pools are setup with the core states excluded, while ib
! includes the core states. To avoid going past array bounds, we must subtract vwfn%ncore_excl
! from ib here. If we had ib in this if statement instead of ib2, then this line and the line
! above where we cycle if ib < vwfn%ncore_excl would mean that we would not enter this loop
! and store the wavefunctions for the highest nv-ncore_excl valence states, where nv is the
! number of valence states kept in the calculation. In the worst case scenario, where
! ncore_excl > nv, we would store no wavefunctions at all.
        ib2=ib-vwfn%ncore_excl
        if(ib2.le. vwfn%nband) then
! Write to valence file

          if (peinf%doiownv(ib2)) then
            iiii=peinf%indexv(ib2)+(ik-1)*peinf%nvownactual
            intwfnv%cg(1:kp%ngk(ik),iiii,1:kp%nspin*kp%nspinor)=zc(1:kp%ngk(ik),1:kp%nspinor*kp%nspin)
          endif

        else

! JRD/JBN: For Metals

          if(ib2.le.vwfn%nband+pol%ncrit) then
            if (peinf%doiownv(ib2)) then
              iiii=peinf%indexv(ib2)+(ik-1)*peinf%nvownactual
              intwfnv%cg(1:kp%ngk(ik),iiii,1:kp%nspin*kp%nspinor)=zc(1:kp%ngk(ik),1:kp%nspinor*kp%nspin)
            endif
          endif

! Write to conduction file
          iwritecb=0
          dest=ib2-vwfn%nband
          if (peinf%doiownc(dest)) then
            iwritecb=1
          endif
          if(iwritecb .eq. 1) then
            iiii=peinf%indexc(ib2-vwfn%nband)+(ik-1)*peinf%ncownactual
            intwfnc%cbi(iiii)=ib2-vwfn%nband
            intwfnc%cg(1:kp%ngk(ik),iiii,1:kp%nspin*kp%nspinor)=zc(1:kp%ngk(ik),1:kp%nspinor*kp%nspin)
          endif
        endif
      enddo  ! of loop i over kp%mnbands

      SAFE_DEALLOCATE(isort)
      SAFE_DEALLOCATE(zc)

    enddo ! of loop over k-points
    call progress_free(prog_info)

    POP_SUB(read_wavefunctions)
    return
  end subroutine read_wavefunctions

!================================================================================
!
! subroutine setup_parallel_freqs   Originally By DVF                 Last Modified 05/08/2015
!
! Setup the global communicator for epsilon. If using paralel frequencies,
! also setup the MPI groups needed for that.
!
!===============================================================================

  subroutine setup_parallel_freqs(pol)
    type (polarizability), intent(inout) :: pol

    integer :: orig_group,ii,ifreq,group_size

    PUSH_SUB(setup_parallel_freqs)

! Get number of processors in your frequency eval group
    peinf%npes_freqgrp = peinf%npes/pol%nfreq_group

    if(pol%nfreq_group .gt. 1) then
! Get the handle from the original group
#ifdef MPI
      call MPI_Comm_Group(MPI_COMM_WORLD,orig_group,mpierr)
! NOTE : we redefine npes here!!!
! Through this one simple command we select the largest subset of processors such that
! the number of processors is divisible by pol%nfreq_group. This one variable determines
! which processors get included in the band distribution and ScalaPACK. A few lines, with
! a very large effect.

      peinf%npes_orig=peinf%npes
      peinf%npes=peinf%npes_freqgrp*pol%nfreq_group

! DVF: if we are doing frequencies in parallel, we need to create the needed MPI groups.
! Get frequency evaluation group number and rank, then determine array of ranks
! that tells you what subset of processors are in your frequency evaluation group.
! The create the mpi group and group communicator. For the processors that are not used
! in the calculation, we also create groups for them. THIS MUST BE DONE because MPI
! only allows child communicators to be created from parent communicators when all the
! processors in the parent communicators get mapped to one of the child communicators, i.e.
! you can`t just leave some processors out of MPI_WORLD_COMM idle. They have to have their
! own garbage group.

      peinf%igroup_f=peinf%inode/peinf%npes_freqgrp
      peinf%rank_f = mod(peinf%inode,peinf%npes_freqgrp)
      if(peinf%inode .lt. peinf%npes) then
        group_size=peinf%npes_freqgrp
        SAFE_ALLOCATE(peinf%ranks,(group_size))
        do ii=1,group_size
          peinf%ranks(ii)=peinf%igroup_f*peinf%npes_freqgrp+ii-1
        enddo
      else
        group_size=peinf%npes_orig-peinf%npes
        SAFE_ALLOCATE(peinf%ranks,(group_size))
        do ii=1,group_size
          peinf%ranks(ii)=peinf%npes+ii-1
        enddo
      endif
      call create_mpi_group(orig_group,group_size,peinf%ranks,peinf%freq_group,peinf%freq_comm)
      SAFE_DEALLOCATE_P(peinf%ranks)
! Same thing, but for matrix element communication group
      peinf%igroup_mtxel= mod(peinf%inode,peinf%npes_freqgrp)
      peinf%rank_mtxel = peinf%inode/peinf%npes_freqgrp
      if(peinf%inode .lt. peinf%npes) then
        group_size=pol%nfreq_group
        SAFE_ALLOCATE(peinf%ranks,(group_size))
        do ii=1,group_size
          peinf%ranks(ii)=mod(peinf%igroup_mtxel*pol%nfreq_group+ii-1,pol%nfreq_group)*peinf%npes_freqgrp+&
            (peinf%igroup_mtxel*pol%nfreq_group+ii-1)/pol%nfreq_group
        enddo
      else
        group_size=peinf%npes_orig-peinf%npes
        SAFE_ALLOCATE(peinf%ranks,(group_size))
        do ii=1,group_size
          peinf%ranks(ii)=peinf%npes+ii-1
        enddo
      endif
      call create_mpi_group(orig_group,group_size,peinf%ranks,peinf%mtxel_group,peinf%mtxel_comm)
      SAFE_DEALLOCATE_P(peinf%ranks)
#endif
    else
      peinf%rank_mtxel = 0
      peinf%rank_f = peinf%inode
      peinf%igroup_f=0
      peinf%igroup_mtxel=peinf%inode
      peinf%npes_orig = peinf%npes
    endif

! DVF : number of frequencies owned by a given processor. Each frequency group
! owns frequencies igroup_f, igroup_f+pol%nfreq_group, igroup_f+2*pol%nfreq_group, etc.
! Note that since not all processors own the same number of frequencies (some may have one
! more than others), some processors will be idle during the computation of one of the
! frequencies. To have no processors idle, you would have to have pol%nfreq_group divide
! pol%nfreq which is in different general to do since you don`t know pol%nfreq in advance (it`s
! of course not that hard to figure out). We don`t want to make life harder on the user, so
! we just let the processors be idle by default. When calculating many frequencies the performance
! hit should be quite minimal.

    pol%nfreq_in_group=0
    do ifreq=1,pol%nFreq
      if( mod(ifreq-1,pol%nfreq_group) .eq. peinf%igroup_f) then
        pol%nfreq_in_group=pol%nfreq_in_group+1
      endif
    enddo

    pol%os_nsfreq_para=0
    do ifreq=1,pol%nsFreq
      if( mod(ifreq-1,pol%nfreq_group) .eq. peinf%igroup_f) then
        pol%os_nsfreq_para=pol%os_nsfreq_para+1
      endif
    enddo

    POP_SUB(setup_parallel_freqs)
    return
  end subroutine setup_parallel_freqs
end module input_m
