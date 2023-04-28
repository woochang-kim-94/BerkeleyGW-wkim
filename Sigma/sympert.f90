!=============================================================================
! Symmetry unfolding program for electron-phonon matrix elements from
! GWPT and DFPT.
!
! Input calculation of gkq matrix elements should be performed with a
! fixed set of wavefunctions on full k-grid Brillouin zone (BZ), and
! with symmetry reduced q-grid wedge and full atom-displacement perturbations
! at each q. The program  will unfold the q-grid to full BZ.
!
! The KEY feature is to include the GAUGE-RECOVERING technique, such
! that the unfoled gkq matrix elements on full BZ were calculated AS IF
! they were directly calculated with no symmetry used at all. This
! gauge-recovering technique MUST BE USED for Wannier interpolation of
! gkq later, because the set of Wannier functions were derived based on
! the original set of wavefunctions with the specific gauge there.
!
! This program is a serial program, with two major steps:
! 1. Find out the unitary transformation between the originally calculated
!    wavefunctions at k`, and from the the symmetry-rotated wavefunctions
!    from k, where k is in the irreducible wedge.
! 2. Symmetry unfolding phonon atom-displacement perturbations.
!
!
! Originally by Zhenglu Li, at LBNL and UC Berkeley, 2020.
! Distributed as a part of the BerkeleyGW package.
! Email: lzlwell@berkeley.edu
!
! Complex flavor only.
!
!=============================================================================

#include "f_defs.h"

program sympert

#ifdef CPLX

  use global_m
  use blas_m
  use find_kpt_match_m
  use gmap_m
  use input_utils_m
  use io_utils_m
  use misc_m
  use sort_m
  use wfn_rho_vxc_io_m
  use sympert_utils_m
  use sympert_phonon_m

  implicit none

  type(crystal) :: crys, crys_dummy
  type(symmetry) :: syms, syms_dummy
  type(kpoints) :: kp, kp_dummy
  type(gspace) :: gvec, gvec_dummy
  character(len=3) :: sheader
  character(len=256) :: wfnfile
  character(len=256) :: rhofile
  character(len=256) :: gkqfile
  character(len=256) :: gkqfile_output
  integer :: ik, ib, is, iflavor, ig
  integer :: isym
  !
  ! variables for I/O
  SCALAR, allocatable :: zwfn(:,:), wfn_coef_fbz(:,:,:,:)
  integer :: ioffsetk
  integer, allocatable :: wfn_gvec_local(:,:), gvlocal(:,:) !< G-vectors corresponding to each coefficient, local to k
  integer, allocatable :: gindex(:) !< G-index map from local list to master list
  integer :: narg
  logical :: is_fmt
  type(progress_info) :: prog_info !< a user-friendly progress report
  !
  ! variables for wavefunctions gauge recovering
  type(wavefunctions) :: wfn_input, wfn_sym  ! input wavefunctions and symmetry-generated wavefunctions
  real(DP), allocatable :: ekin(:) ! kinetic energy array
  integer, allocatable :: gvmaster(:,:) ! G-vectors master list for temporary use
  integer, allocatable :: isort(:) ! sorting order index
  real(DP) :: eb1, eb2, eb3 ! band energies for degenercy check
  integer, allocatable :: ngk_original(:), ksave_map(:) ! for saving ngk before resizing kp
  integer :: nktot, iktot ! original total nk and its looper before resizing kp
  !
  ! variables for phonon perturbations
  type(perturbations) :: phonpert_input, phonpert_sym
  !
  ! variables for e-ph matrix elements
  type(gkqmatrix) :: gkqmat_input, gkqmat_output
  integer, allocatable :: sym_to_invsym(:), invsym_to_sym(:) ! maps between S and S-inverse
  integer :: iinvsym
  real(DP), allocatable :: qlist_fbz(:,:), qlist_irr(:,:)
  integer, allocatable :: symq_map(:,:) ! size (2, nqfbz), for each iqfbz, first index saves iqirr, second saves isym
  integer :: nqfbz, nqirr ! number of qpoints in full and irreducible BZ
  integer :: iqfbz, iqirr
  integer :: iatom, idir, jdir ! labels atom index and direction index
  real(DP) :: sym_inv_mtrx_k(3,3), sym_inv_mtrx_r(3,3) ! inverse operation symmetry atrix in k and r spaces
  real(DP) :: mtrx_aux(3,3) ! auxiliary matrix for general purpose
  integer, allocatable :: iqcomputed(:)
  integer :: ib_k, ib_kq ! band index for k and k+q point
  integer :: ibn, ibm ! band index for unitary transformation
  real(DP) :: phi 
  SCALAR :: phase ! need to be complex
  integer :: idx_k, idx_kq ! index for k and k+q in wfn_sym
  real(DP) :: kpt(3), kqpt(3) ! k and k+q points
  SCALAR :: dk, dkq ! coefficient for D^k and D^{k+q}
  integer :: ikq ! index for k+q in wfn_input
  integer :: iinvsk ! index for S^{-1}k in wfn_input
  real(DP) :: invskpt(3) ! S^{-1}k point 
  real(DP), allocatable :: apos_frac(:,:) ! fractional coordinates of atom positions for phase calculation
  real(DP) :: sqirr(3) ! S^qirr point

!=========================================================================
!> check complex compilation

#ifndef CPLX
  call die('sympert.x MUST be compiled with complex flavor.')
  ! will not reach here. keep it for clarity for developers
#endif

!#BEGIN_INTERNAL_ONLY

!=========================================================================
!> reading input files

  narg = command_argument_count()
  if (narg.ne.3) then
    call die('Usage: sympert.x rho.sym.cplx wfn.fullBZ.cplx gkq.mat')
  endif
  call get_command_argument(1, rhofile)
  call get_command_argument(2, wfnfile)
  call get_command_argument(3, gkqfile)

  ! always dealing with binary rho and wavefunctions
  is_fmt = .false.

  ! we read symmetries from charge density file which usually includes symmetry in scf calculation
  ! other header information are all read from wavefunction, which must be generated with no symmetry for GWPT
  call open_file(unit=95, file=trim(rhofile), form='unformatted', status='old')

  write(*,*) 'Reading symmetries from ', trim(rhofile)
  write(*,*)

  sheader = 'GET'
  iflavor = -1
  call read_header_type(95, is_fmt, sheader, iflavor, kp_dummy, gvec_dummy, syms, crys_dummy, dont_warn_kgrid=.true.)
  if (sheader .ne. 'RHO') call die('Please provide the correct RHO file.')

  call close_file(95)

  ! reading wavefunction file
  call open_file(unit=7, file=trim(wfnfile), form='unformatted', status='old')

  write(*,*) 'Reading header and gvectors from ', trim(wfnfile)
  write(*,*)

  sheader = 'GET'
  iflavor = -1
  call read_header_type(7, is_fmt, sheader, iflavor, kp, gvec, syms_dummy, crys, dont_warn_kgrid=.true.)
  if (sheader .ne. 'WFN') call die('Please provide the correct WFN file.')

  ! print header information
  call print_header(kp, gvec, syms, crys)

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  SAFE_ALLOCATE(wfn_gvec_local, (3,sum(kp%ngk)))
  SAFE_ALLOCATE(gvlocal, (3,kp%ngkmax))

  ! read in charge density G-vectors, i.e. the master list of G-pool
  ! stored in gvec%components
  call read_gvectors(7, is_fmt, gvec%ng, gvec%ng, gvec%components)

  if (iflavor .eq. 1) then
    call die('sympert.x should only be used with complex flavor of wavefunctions.')
  endif

  SAFE_ALLOCATE(zwfn, (kp%ngkmax,kp%nspin*kp%nspinor))
  SAFE_ALLOCATE(wfn_coef_fbz, (kp%ngkmax,kp%mnband,kp%nspin*kp%nspinor,kp%nrk))
  wfn_coef_fbz(:,:,:,:) = ZERO

  if (kp%nspin*kp%nspinor .ne. 1) then
    call die('sympert.x only supports nspin = 1 now.')
  endif

  write(*,*) 'Reading wavefunctions from ', trim(wfnfile)

  ioffsetk = 0
  call progress_init(prog_info, 'reading wavefunctions', 'k-point', kp%nrk)
  do ik = 1, kp%nrk
    call progress_step(prog_info)
    ! gvlocal will be initialized as zeros in the subroutine call
    call read_gvectors(7, is_fmt, kp%ngk(ik), kp%ngkmax, gvlocal)
    wfn_gvec_local(:, ioffsetk+1:ioffsetk+kp%ngk(ik))=gvlocal(:,1:kp%ngk(ik))
    do ib = 1, kp%mnband
      !write(*,'(8X,A,I6,A,I6)') 'kpoint ', ik, '     band ', ib 
      ! read in wavefunction of one band ib at one kpoint ik, stored in zwfn
      ! zwfn will be initilzed as zeros in the subroutine call
      call read_complex_data(7, is_fmt, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, zwfn)
      wfn_coef_fbz(1:kp%ngk(ik),ib,:,ik) = zwfn(1:kp%ngk(ik),:)
    enddo
    ioffsetk = ioffsetk + kp%ngk(ik)
  enddo
  call progress_free(prog_info)

  ! end of reading wavefunctions
  call close_file(7)

!=========================================================================
!> reading gkq.mat

  ! read in gkq.mat and store in gkqmat_input
  call init_read_free_gkqmatrix(1, gkqmat_input, -1, crys=crys, gkqfile=gkqfile)

!=========================================================================
!> sort global gvec%components 

  ! Note: This step of soriting global G components is not necessary for sympert purpose. 
  !       We do this just to be consistent with input.f90 for sigma.  We will have to do 
  !       the sorting for gmap.  This is because gmap uses a logic that it asks the new rotated k 
  !       to provide all its lowest energy G-components (sorted by ekin), and then go back to 
  !       the original k point to look for the corresponding G vector.
  !       This sorting procedure will also make the master G-list more intuitive.
  SAFE_ALLOCATE(isort, (gvec%ng))
  SAFE_ALLOCATE(ekin, (gvec%ng))

  call kinetic_energies(gvec, crys%bdot, ekin)
  call sortrx(gvec%ng, ekin, isort, gvec=gvec%components)

  SAFE_ALLOCATE(gvmaster, (3,gvec%ng))
  gvmaster(:,:) = gvec%components(:,:)

  do ig = 1, gvec%ng
    gvec%components(:, ig) = gvmaster(:, isort(ig))
  enddo

  SAFE_DEALLOCATE(isort)
  SAFE_DEALLOCATE(ekin)
  SAFE_DEALLOCATE(gvmaster)

  ! write out master G-list
  call write_master_g_list(gvec)

  ! initialize gvec%index_vec to map 3D FFT index to 1D array, to enable findvector
  ! this step must be done after sorting G-vectors, because the sorted new G-list is
  ! now the master G-list
  call gvec_index(gvec)
  ! find gindex from local G-list to master G-list
  SAFE_ALLOCATE(gindex, (sum(kp%ngk)))
  do ig = 1, sum(kp%ngk)
    call findvector(gindex(ig), wfn_gvec_local(:,ig), gvec)
  enddo

  ! write out local G-list with map to master G-list
  call write_local_g_list(kp, wfn_gvec_local, gindex)

  ! free up space
  SAFE_DEALLOCATE(wfn_gvec_local)
  SAFE_DEALLOCATE(gvlocal)
  SAFE_DEALLOCATE(zwfn)

!=========================================================================
!> resize kp object, based on gkqmat_inpt
  ! after I/O and store wavefunction is done, we can resize the kp object
  ! because the file WFN may contain more bands/kpoints than needed that
  ! gkq matrix elements were calculated in DFPT and GWPT

  SAFE_ALLOCATE(ngk_original, (kp%nrk))
  SAFE_ALLOCATE(ksave_map, (kp%nrk))

  nktot = kp%nrk
  ngk_original(:) = kp%ngk(:)
  ksave_map(:) = 0

  call resize_kpoints_object(kp, gkqmat_input, ksave_map)

!=========================================================================
!> wavefunctions objects for symmetry rotation and gauge recovering

  call init_free_wavefunctions(1, wfn_input, kp)

  ! pass input data to wfn_input
  wfn_input%ngk(:) = kp%ngk(:)
  wfn_input%klist(:,:) = kp%rk(:,:)
  wfn_input%isym = 1  ! for wfn_input, it is related to global kp%rk via identity, i.e. isym = 1 

  ioffsetk = 0
  do iktot = 1, nktot
    ! gindex is saved as continuous order with original kpoints (before resizing kp)
    ! so we need the main loop to be over nktot (before resizing), and ksave_map designed correspondingly
    ik = ksave_map(iktot)

    if (ik .gt. 0) then
      wfn_input%kmap(ik) = ik ! for wfn_input, each ik points to the corresponding ik in kp%rk
      wfn_input%kg0(:,ik) = 0 ! for wfn_input, all G0 umklapp are zero
  
      wfn_input%gindex(1:kp%ngk(ik), ik) = gindex(ioffsetk+1:ioffsetk+ngk_original(iktot))
  
      do is = 1, kp%nspin*kp%nspinor
        do ib = 1, kp%mnband
          wfn_input%cg(1:kp%ngk(ik), ib, is, ik) = wfn_coef_fbz(1:kp%ngk(ik), ib, is, iktot)
        enddo
      enddo
    endif

    ioffsetk = ioffsetk + ngk_original(iktot)
  enddo

  ! free up space
  SAFE_DEALLOCATE(gindex)
  SAFE_DEALLOCATE(wfn_coef_fbz)
  SAFE_DEALLOCATE(ngk_original)
  SAFE_DEALLOCATE(ksave_map)

  ! initialize wfn_input%umat_gauge for wfn_input (not used, only for data consistency)
  wfn_input%umat_gauge(:,:,:) = ONE

  ! find band degeneracy
  do ik = 1, kp%nrk
    do is = 1, kp%nspin*kp%nspinor
      if (kp%mnband .eq. 1) then
        wfn_input%degeneracy(:,is,ik) = .false.
      elseif (kp%mnband .eq. 2) then
        wfn_input%degeneracy(:,is,ik) = .false.
        eb1 = kp%el(1, ik, is)  ! in Ry unit
        eb2 = kp%el(2, ik, is)
        if (abs(eb1-eb2) .lt. TOL_Small) then
          wfn_input%degeneracy(:,is,ik) = .true.  ! both bands belong to degenerate set
        endif
      else ! for all cases that kp%mnband >= 3
        ! kp%el is ordered with energies from lowest to highest
        wfn_input%degeneracy(:,is,ik) = .false.

        ! check the first band
        eb1 = kp%el(1, ik, is)
        eb2 = kp%el(2, ik, is)
        if (abs(eb1-eb2) .lt. TOL_Small) then
          wfn_input%degeneracy(1,is,ik) = .true.  ! label the first band as degenerate
        endif

        ! check 2 <= ib <= kp%mnband-1
        do ib = 2, kp%mnband-1
          eb1 = kp%el(ib-1, ik, is)
          eb2 = kp%el(ib, ik, is)
          eb3 = kp%el(ib+1, ik, is)
          if ((abs(eb1-eb2) .lt. TOL_Small) .or. (abs(eb3-eb2) .lt. TOL_Small)) then
            wfn_input%degeneracy(ib,is,ik) = .true.
          endif
        enddo

        ! check the last band
        eb1 = kp%el(kp%mnband-1, ik, is)
        eb2 = kp%el(kp%mnband, ik, is) 
        if (abs(eb1-eb2) .lt. TOL_Small) then
          wfn_input%degeneracy(kp%mnband,is,ik) = .true.  ! label the last band as degenerate
        endif
      endif ! different kp%mnband cases
    enddo
  enddo

  ! print band energies and degeneracies information of the input wavefunction
  ! kp type itself has a degeneracy attribute, but we do not use it here
  write(*,*)
  write(*,*) 'Band energies (eV) and the degeneracy labels (T/F) of the input wavefunction'
  do ik = 1, kp%nrk
    write(*,*)
    write(*,'(3x,a,3x,a,5x,a,3x,a,3x,a)') 'ik', 'is', 'ib', 'energy (eV)', 'degenerate'
    do is = 1, kp%nspin*kp%nspinor
      do ib = 1, kp%mnband
        write(*,'(i5,2x,i3,1x,i6,2x,f12.6,3x,l3)') ik, is, ib, kp%el(ib, ik, is)*ryd, &
                                                   wfn_input%degeneracy(ib, is, ik)
      enddo
    enddo
  enddo

!=========================================================================
!> find maps between symmetries and their inverse symmetries which are all in one group
  SAFE_ALLOCATE(sym_to_invsym, (syms%ntran))
  SAFE_ALLOCATE(invsym_to_sym, (syms%ntran))

  call find_maps_sym_invsym(syms, sym_to_invsym, invsym_to_sym)

!=========================================================================
!> for eqch iqfbz of qpoints in full BZ, find the unique pair of (iqirr, iqsym) to generate this iqfbz
  nqfbz = gkqmat_input%qmesh(1) * gkqmat_input%qmesh(2) * gkqmat_input%qmesh(3)
  nqirr = gkqmat_input%nq

  SAFE_ALLOCATE(qlist_fbz, (3, nqfbz))
  SAFE_ALLOCATE(symq_map, (2, nqfbz))

  call gen_unfold_q_sym(gkqmat_input, syms, qlist_fbz, symq_map)

  SAFE_ALLOCATE(qlist_irr, (3, nqirr))
  qlist_irr(:,:) = gkqmat_input%qlist(:,:)

!
!=========================================================================
!>                              MAIN LOOP
  ! entering loop over phonon q points on full q-BZ. Each iq on full BZ is associated with
  ! a pair (irq, isym) of index of reduced q-list and the i-th symmetry that brings irq to iq
  !

  write(*,*)
  write(*,*)
  write(*,*) '==================== Entering main loop ===================='
  write(*,*)

  ! prepare gkqmat_output
  call init_read_free_gkqmatrix(1, gkqmat_output, nqfbz, crys=crys, kp=kp)
  gkqmat_output%qmesh(:) = gkqmat_input%qmesh(:)
  gkqmat_output%qlist(:,:) = qlist_fbz(:,:)

  ! record computed iqfbz points, at the end, all elements should be 1 (calculated and only once)
  SAFE_ALLOCATE(iqcomputed, (nqfbz))
  iqcomputed(:) = 0

  ! convert atom positions into fractional coordinates
  SAFE_ALLOCATE(apos_frac, (3,crys%nat))
  do iatom = 1, crys%nat
    apos_frac(:, iatom) = matmul(transpose(crys%bvec), crys%apos(:,iatom))
  enddo

  ! we make the out-most loop as over symmetry, which is usually few tens, no more than 48
  ! because under each symmetry, the rotated wavefunction and computed umat_gauge is the same
  do isym = 1, syms%ntran
    ! under each isym, we check if we need to compute this iqfbz using this isym
    ! each iqfbz must be calculated once and only once, we save checks in iqcomputed

    write(*,'(1x,a,i3)') '======> Symmetry loop, working on isym =', isym

    !======> electron wavefunction
    ! for each isym/iq, initialize wfn_sym
    call init_free_wavefunctions(1, wfn_sym, kp)

    ! generate wfn_sym using wfn_input and isym
    call gen_wfn_sym(wfn_input, kp, gvec, syms, crys, isym, wfn_sym)

    ! calculate the unitary gauge transformation matrix
    call calc_gauge_transform(kp, gvec, wfn_input, wfn_sym)

    ! write gauge transformation matrix to file for each wfn_sym
    ! if different iqfbz uses the same isym, the file will be overwritten with exactly the same content
    ! because umat only depends on which isym to be used
    call write_umat_gauge(kp, wfn_sym)

    do iqfbz = 1, nqfbz
      ! we check whether this iqfbz uses this isym
      if (symq_map(2, iqfbz) .ne. isym) cycle

      ! after the above if-statement, if we need to calculate this iqfbz, we record it
      iqcomputed(iqfbz) = iqcomputed(iqfbz) + 1

      ! for the given iqfbz, we use the specific iqirr associted with the specific symmetry isym
      iqirr = symq_map(1, iqfbz)
      iinvsym = sym_to_invsym(isym)
 
      write(*,*)
      write(*,'(1x,a,i3,1x,a,i4)') '---> Progress: ', sum(iqcomputed(:)), '/', nqfbz
      write(*,'(1x,a,i3,a,i3,a,i3)') '---> Within isym =', isym, ', unfolding iqfbz =', iqfbz, ' in qmesh from irreducible iqirr =', iqirr
      write(*,*)

      !======> phonon perturbation symmetry
      ! Irregardless of how DFT codes (QE, Abinit, etc.) treat and record symmetry,
      ! the symmetry matrices recorded in BerkeleyGW WFN are in the reciprocal basis
      ! for direct rotation of wavevectors k/q and planewave components G.
      !
      ! Here, we directly take the symmetry operation from WFN. For a specific symmetry
      ! operation isym, we need to convert the rotation matrix to lattice-vector basis
      ! for the atomic position rotation. For wavevector q, we use the original matrices
      ! in reciprocal basis.
      !
      ! The BGW syms object saves the fractional translation as ft*2*Pi, where ft is
      ! in the lattice-vector basis.
    
      ! for each isym, we find out the equivalent site of each atom.
      call init_free_perturbations(1, phonpert_sym, crys)
  
      ! find equivalent site for this perturbation that is related to original perturbation by isym
      call find_equiv_site(crys, syms, isym, phonpert_sym)
  
      ! transfer data to perturbations object
      phonpert_sym%iqfbz = iqfbz
      phonpert_sym%qpoint(:) = qlist_fbz(:, iqfbz)
      phonpert_sym%iqirr = iqirr
      phonpert_sym%isym = isym
  
      ! convert symmetry in k to symmetry in crystal lattice, for this iinvsym
      ! sym_inv_mtrx_r will be used for phonon perturbations, which are defined in crystal lattice
      sym_inv_mtrx_k(:,:) = real(syms%mtrx(:, :, iinvsym), DP)
      mtrx_aux(:,:) = transpose(sym_inv_mtrx_k)
      call invert_matrix(mtrx_aux, sym_inv_mtrx_r)
 
      !======> core part of gkqmat unfolding
      call progress_init(prog_info, 'unfolding electron-phonon matrix elements', 'band', kp%mnband)

      do ib_kq = 1, kp%mnband
        call progress_step(prog_info)

        do iatom = 1, crys%nat
          ! idir is related to symmetry rotation of perturbed potential
          do idir = 1, 3
            do ik = 1, kp%nrk
              ! find indices for k+Sq and k in wfn_sym%klist
              kpt(:) = gkqmat_output%klist(:,ik)
              kqpt(:) = kpt(:) + qlist_fbz(:,iqfbz)

              idx_k = kindex(kpt, wfn_sym%klist, kp%nrk)  ! index in wfn_sym
              idx_kq = kindex(kqpt, wfn_sym%klist, kp%nrk)  ! index in wfn_sym

              ! find k+Sq in wfn_input%klist for degeneracy information
              ! (although we can use wfn_sym for such info, but using wfn_input is more intuitive)
              ikq = kindex(kqpt, wfn_input%klist, kp%nrk)  ! index in wfn_input

              ! find S^{-1}k index in wfn_input
              invskpt(:) = matmul(sym_inv_mtrx_k(:,:), kpt(:))
              iinvsk = kindex(invskpt, wfn_input%klist, kp%nrk)

              ! compute the phase e^{i q \cdot x(kappa) - i (S^k q) \cdot x(KAPPA) + i (S^k q) \cdot v  } 
              ! using the formalism based on the inclusion of e^{-iSk.v} phase
              ! we keep this part under ik loop because logically the last e^{i (S^k q) \cdot v} is 
              ! fundamentally from the difference between kpoint difference between k+Sq and k
              sqirr(:) = matmul( real(syms%mtrx(:,:,isym), DP), qlist_irr(:, iqirr) )
              phi = dot_product( qlist_irr(:, iqirr), apos_frac(:, phonpert_sym%apos_map(iatom)) )
              phi = phi - dot_product( sqirr(:), apos_frac(:, iatom) )
              phi = phi * 2.0d0 * PI_D  ! note the unit and 2pi factor
              phi = phi + dot_product( sqirr(:), syms%tnp(:, isym) )
              phase = CMPLX(cos(phi), sin(phi))

              do ib_k = 1, kp%mnband
                ! initialize to be zero, for accumulation of value in the following loops
                gkqmat_output%gkq(ib_k, ik, idir, iatom, ib_kq, iqfbz) = ZERO

                ! nspin = 1 only
                is = 1

                ! the following two loops are for summation in unitary transformation
                do ibm = 1, kp%mnband
                  ! degeneracy check to speed up
                  if ((.not. wfn_input%degeneracy(ib_kq, is, ikq)) .and. (ibm .ne. ib_kq)) cycle
                  !
                  do ibn = 1, kp%mnband
                    ! degeneracy check to speed up
                    if ((.not. wfn_input%degeneracy(ib_k, is, ik)) .and. (ibn .ne. ib_k)) cycle
                    !
                    dk = wfn_sym%umat_gauge(ibn, ib_k, idx_k)
                    dkq = wfn_sym%umat_gauge(ibm, ib_kq, idx_kq)
                    !
                    ! jdir is used for symmetry rotation of perturbed potential
                    do jdir = 1, 3
                      gkqmat_output%gkq(ib_k, ik, idir, iatom, ib_kq, iqfbz) = &
                        gkqmat_output%gkq(ib_k, ik, idir, iatom, ib_kq, iqfbz) + &
                        phase * MYCONJG(dkq) * dk * &
                        gkqmat_input%gkq(ibn, iinvsk, jdir, phonpert_sym%apos_map(iatom), ibm, iqirr) * &
                        sym_inv_mtrx_r(jdir, idir)
                    enddo ! jdir
                  enddo ! ibn
                enddo ! ibm
              enddo ! ib_k
            enddo ! ik
          enddo ! idir
        enddo ! iatom
      enddo ! ib_kq
      call progress_free(prog_info)
  
      ! free up space at the end of each iteration in the loop iqfbz
      call init_free_perturbations(0, phonpert_sym)
  
      ! loop over output atom perturbations
    enddo ! iqfbz

    ! free up space at the end of each iteration in the loop isym
    call init_free_wavefunctions(0, wfn_sym)
  enddo ! isym

  ! final check computed iqfbz
  write(*,*) 'check iqcomputed'
  write(*,*) iqcomputed

  do iqfbz = 1, nqfbz
    if (iqcomputed(iqfbz) .lt. 1) then
      call die('Some qpoint has not been computed.')
    elseif (iqcomputed(iqfbz) .gt. 1) then
      call die('Some qpoint has been calculated more than once.')
    endif
  enddo

  ! write output file
  gkqfile_output = trim(gkqfile) // '.unfolded'
  call write_gkqmatrix(gkqmat_output, crys, gkqfile_output)

!=========================================================================
!> end of program, deallocating
  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE(sym_to_invsym)
  SAFE_DEALLOCATE(invsym_to_sym)
  SAFE_DEALLOCATE(qlist_fbz)
  SAFE_DEALLOCATE(qlist_irr)
  SAFE_DEALLOCATE(symq_map)
  SAFE_DEALLOCATE(iqcomputed)
  SAFE_DEALLOCATE(apos_frac)

  call init_free_wavefunctions(0, wfn_input)
  call init_read_free_gkqmatrix(0, gkqmat_input)
  call init_read_free_gkqmatrix(0, gkqmat_output)

  call kp%free()

!#END_INTERNAL_ONLY

  write(*,*) 'Program sympert.x finishes.'

#endif

end program sympert

