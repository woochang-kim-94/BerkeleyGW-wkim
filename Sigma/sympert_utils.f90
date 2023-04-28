!=============================================================================
! Utility subroutines for sympert.f90
!
! Originally by Zhenglu Li
!=============================================================================

#include "f_defs.h"

module sympert_utils_m

#ifdef CPLX

  use global_m
  use misc_m
  use gmap_m
  use input_utils_m
  use sort_m
  use io_utils_m
  use blas_m
  use sympert_phonon_m, only: gkqmatrix

  implicit none

  private

  public :: &
    wavefunctions, &
    print_header, &
    write_master_g_list, &
    write_local_g_list, &
    resize_kpoints_object, &
    init_free_wavefunctions, &
    gen_wfn_sym, &
    calc_gauge_transform, &
    write_umat_gauge, &
    find_maps_sym_invsym, &
    kindex

  type wavefunctions
    integer, pointer :: ngk(:) ! ngk
    integer, pointer :: gindex(:,:) ! local G to master G-list
    real(DP), pointer :: klist(:,:)
    integer, pointer :: kmap(:) ! mapping relation from local klist to global kp%rk
    integer, pointer :: kg0(:,:) ! umklapp vectors combined with kmap
    SCALAR, pointer :: cg(:,:,:,:) ! wavefunction coefficient (ngmax,nband,nspin*nspinor,nk)
    integer :: isym ! i-th symmetry operation in syms used to generate this klist
    SCALAR, pointer :: umat_gauge(:,:,:) ! unitary matrix of <wfn_sym|wfn_input>, (nband,nband,nk), for spinless
    logical, pointer :: degeneracy(:,:,:) ! degeneracy labels for (nband,nspin*nspinor,nk)
  end type wavefunctions

contains

!#BEGIN_INTERNAL_ONLY

!=========================================================================
!> print detailed header information
subroutine print_header(kp, gvec, syms, crys)

  type(kpoints), intent(in) :: kp
  type(gspace), intent(in) :: gvec
  type(symmetry), intent(in) :: syms
  type(crystal), intent(in) :: crys
  !
  integer :: isym, i

  PUSH_SUB(print_header)

  write(*,*)
  write(*,*) 'Printing header information for input charge density and wavefunction files'
  write(*,*)

  write(*,*) 'kp%nspin = ', kp%nspin
  write(*,*) 'kp%nspinor = ', kp%nspinor
  write(*,*) 'kp%ecutwfc = ', kp%ecutwfc
  write(*,*) 'kp%ngkmax = ', kp%ngkmax
  write(*,*) 'kp%nrk = ', kp%nrk
  write(*,*) 'kp%mnband = ', kp%mnband
  write(*,*)

  write(*,*) 'gvec%ng = ', gvec%ng
  write(*,*) 'gvec%ecutrho = ', gvec%ecutrho
  write(*,'(1x,a,i6,i6,i6)') 'gvec%FFTgrid = ', gvec%FFTgrid
  write(*,*)

  write(*,*) 'syms%ntran = ', syms%ntran
  write(*,*) 'syms%cell_symmetry = ', syms%cell_symmetry
  do isym = 1, syms%ntran
    write(*,*) 'Symmetry matrix #', isym
    do i = 1, 3
      write(*,'(10x,i6,i6,i6)') syms%mtrx(i, 1:3, isym)
    enddo
    write(*,*) '    -- associated fractional translation --'
    write(*,'(3(f15.8))') syms%tnp(1:3, isym)
  enddo
  write(*,*)

  write(*,*) 'crys%nat = ', crys%nat
  write(*,*) 'crys%atyp = ', crys%atyp
  do i = 1, crys%nat
    write(*,'(1x,a,i3,a,3(f15.8))') 'crys%apos(:, ', i, ') = ', crys%apos(1:3, i)
  enddo
  write(*,*) 'crys%alat = ', crys%alat
  write(*,*) 'crys%blat = ', crys%blat
  write(*,*) 'crys%celvol = ', crys%celvol
  write(*,*) 'crys%recvol = ', crys%recvol
  write(*,*) 'avec'
  do i = 1, 3
    write(*,'(1x,a,i1,1x,a,1x,3(f15.8))') 'a', i , '=', crys%avec(1:3, i)
  enddo
  write(*,*) 'bvec'
  do i = 1, 3
    write(*,'(1x,a,i1,1x,a,1x,3(f15.8))') 'b', i , '=', crys%bvec(1:3, i)
  enddo
  write(*,*)

  POP_SUB(print_header)

  return

end subroutine print_header

!=========================================================================
!> write the master G-list to file
subroutine write_master_g_list(gvec)

  type(gspace), intent(in) :: gvec
  !
  character(len=256) :: gfile
  integer :: ig

  PUSH_SUB(write_master_g_list)

  gfile = 'master_G.list'

  call open_file(unit=51, file=trim(gfile), form='formatted', status='replace')

  write(*,*) 'Writing master G-list to file ', trim(gfile)
  write(*,*)

  do ig = 1, gvec%ng
    write(51,'(a,i8,5x,a,3(i6),a)') 'iG = ', ig, 'G = ( ', gvec%components(1:3, ig), ' )'
  enddo

  call close_file(51)

  POP_SUB(write_master_g_list)

  return

end subroutine write_master_g_list


!=========================================================================
!> write local G-list to file, local to each k, with map to master G-list
subroutine write_local_g_list(kp, wfn_gvec_local, gindex)

  type(kpoints), intent(in) :: kp
  integer, allocatable, intent(in) :: wfn_gvec_local(:,:) 
  integer, allocatable, intent(in) :: gindex(:) 
  !
  character(len=256) :: gfile
  integer :: ik, ioffsetk
  integer :: igloc

  PUSH_SUB(write_local_g_list)

  gfile = 'local_G.list'

  call open_file(unit=52, file=trim(gfile), form='formatted', status='replace')

  write(*,*) 'Writing local G-list to file ', trim(gfile)
  write(*,*) 'Sum of ngk is', sum(kp%ngk)
  write(*,*)

  ioffsetk = 0
  do ik = 1, kp%nrk
    do igloc = 1, kp%ngk(ik)
      write(52,'(a,i5,a,i6,3x,a,i6,3x,a,3(i5),a,3x,a,i8)') &
        'ik = ',ik, '  ngk = ', kp%ngk(ik), 'iGloc = ', igloc, 'G = (', &
        wfn_gvec_local(1:3, ioffsetk+igloc), ' )', &
        'map to master iG = ', gindex(ioffsetk+igloc)
    enddo
    ioffsetk = ioffsetk + kp%ngk(ik)
  enddo

  call close_file(52)

  POP_SUB(write_local_g_list)

  return

end subroutine write_local_g_list


!=========================================================================
!> resize the kpoints object with reduced number of kpoints/bands based on gkqmat
subroutine resize_kpoints_object(kp, gkqmat, ksave_map)

  type(kpoints), intent(inout) :: kp
  type(gkqmatrix), intent(in) :: gkqmat
  integer, allocatable, intent(inout) :: ksave_map(:)
  !
  integer :: nktot, iktot
  integer :: nksave, iksave
  integer :: nbtot, nbsave  ! bands saved always start from lowest band, no map needed
  integer :: nfound
  real(DP) :: kdiff(3), knorm
  ! kp variables original values, keep those only used by symperta
  integer :: nspinor, nspin, is
  real(DP) :: ecutwfc
  integer :: ngkmax
  integer, allocatable :: ngk_original(:)
  real(DP), allocatable :: rk_original(:,:)
  real(DP), allocatable :: el_original(:,:,:)

  PUSH_SUB(resize_kpoints_object)

  nspinor = kp%nspinor
  nspin = kp%nspin

  if ((nspinor .ne. 1) .or. (nspin .ne. 1)) then
    call die('resize_kpoints_object: current sympert only supports spinless case.')
  endif

  SAFE_ALLOCATE(ngk_original, (kp%nrk))
  SAFE_ALLOCATE(rk_original, (3, kp%nrk))
  SAFE_ALLOCATE(el_original, (kp%mnband, kp%nrk, kp%nspin))

  nktot = kp%nrk
  nbtot = kp%mnband
  ecutwfc = kp%ecutwfc
  ngkmax = kp%ngkmax  ! we do not need to redefine ngkmax, which is used as an upper bound for safety
  ngk_original(:) = kp%ngk(:)
  rk_original(:,:) = kp%rk(:,:)
  el_original(:,:,:) = kp%el(:,:,:)

  ! after saving needed data, free kp
  call kp%free()

  ! get target size from gkqmat
  nksave = gkqmat%nk
  nbsave = gkqmat%nband

  write(*,*) 'Resizing kpoints and bands:'
  write(*,'(6x,a,i5,2x,a,i5)') 'number of kpoints is reduced:', nktot, '->', nksave
  write(*,'(6x,a,i5,2x,a,i5)') 'number of  bands  is reduced:', nbtot, '->', nbsave
  write(*,*)

  if ((nksave .gt. nktot) .or. (nbsave .gt. nbtot)) then
    call die('Cannot resize kp because requested size is larger than input wavefunction.')
  endif

  ! update kp with reduced size
  kp%nspinor = nspinor
  kp%nspin = nspin
  kp%nrk = nksave ! new size
  kp%mnband = nbsave ! new size
  kp%ecutwfc = ecutwfc
  kp%ngkmax = ngkmax

  SAFE_ALLOCATE(kp%ngk, (kp%nrk))
  SAFE_ALLOCATE(kp%rk, (3, kp%nrk))
  SAFE_ALLOCATE(kp%el, (kp%mnband, kp%nrk, kp%nspin))

  ! generate ksave_map
  ! kpoints are exactly matched, no umklapp involved
  ksave_map(:) = 0
  nfound = 0
  do iktot = 1, nktot
    do iksave = 1, nksave
      kdiff(:) = rk_original(:, iktot) - gkqmat%klist(:, iksave)
      knorm = sqrt(dot_product(kdiff, kdiff))
      
      if(abs(knorm) .lt. TOL_Small) then
        nfound = nfound + 1
        ksave_map(iktot) = iksave
      endif
    enddo
  enddo

  if (nfound .ne. nksave) then
    call die('kpoints in gkqmat is not a subset of kpoints in the input wavefunction.')
  endif

  ! now kpsave_map is successfully computed, and saved as intent out
  ! we put resized data to kp
  is = 1 ! spinless
  do iktot = 1, nktot
    iksave = ksave_map(iktot)

    if (iksave .gt. 0) then
      kp%ngk(iksave) = ngk_original(iktot)
      kp%rk(:,iksave) = rk_original(:,iktot)
      kp%el(1:nbsave,iksave,is) = el_original(1:nbsave,iktot,is) ! reduce bands here
    endif
  enddo

  ! print information
  write(*,*) 'kpoints list used by sympert (after resizing kp)'
  do iksave = 1, kp%nrk
    write(*,'(3x,a,i5,1x,a,3(f9.5))') 'k(', iksave, ') = ', kp%rk(:,iksave)
  enddo
  write(*,*)

  SAFE_DEALLOCATE(ngk_original)
  SAFE_DEALLOCATE(rk_original)
  SAFE_DEALLOCATE(el_original)

  POP_SUB(resize_kpoints_object)

  return

end subroutine resize_kpoints_object


!=========================================================================
!> initialize and free wavefunctions object
subroutine init_free_wavefunctions(switch, wfn, kp)

  integer, intent(in) :: switch
  type(wavefunctions), intent(inout) :: wfn
  type(kpoints), optional, intent(in) :: kp

  PUSH_SUB(init_free_wavefunctions)

  if (switch .eq. 1) then
    if (.not. present(kp)) then
      call die('You need to pass kp into init_free_wavefunctions to initialize wfn object.')
    endif
    SAFE_ALLOCATE(wfn%ngk, (kp%nrk))
    SAFE_ALLOCATE(wfn%gindex, (kp%ngkmax, kp%nrk))
    SAFE_ALLOCATE(wfn%klist, (3, kp%nrk))
    SAFE_ALLOCATE(wfn%kmap, (kp%nrk))
    SAFE_ALLOCATE(wfn%kg0, (3, kp%nrk))
    SAFE_ALLOCATE(wfn%cg, (kp%ngkmax, kp%mnband, kp%nspin*kp%nspinor, kp%nrk))
    SAFE_ALLOCATE(wfn%umat_gauge, (kp%mnband, kp%mnband, kp%nrk))
    SAFE_ALLOCATE(wfn%degeneracy, (kp%mnband, kp%nspin*kp%nspinor, kp%nrk))

    wfn%ngk(:) = 0
    wfn%gindex(:,:) = 0
    wfn%klist(:,:) = 0.0d0
    wfn%kmap(:) = 0
    wfn%kg0(:,:) = 0
    wfn%cg(:,:,:,:) = ZERO
    wfn%isym = 0
    wfn%umat_gauge(:,:,:) = ZERO
    wfn%degeneracy(:,:,:) = .false.
  elseif (switch .eq. 0) then
    SAFE_DEALLOCATE_P(wfn%ngk)
    SAFE_DEALLOCATE_P(wfn%gindex)
    SAFE_DEALLOCATE_P(wfn%klist)
    SAFE_DEALLOCATE_P(wfn%kmap)
    SAFE_DEALLOCATE_P(wfn%kg0)
    SAFE_DEALLOCATE_P(wfn%cg)
    SAFE_DEALLOCATE_P(wfn%umat_gauge)
    wfn%isym = 0
  else
    call die('switch in init_free_wavefunctions can only be 1 (init) or 0 (free).')
  endif

  POP_SUB(init_free_wavefunctions)

  return

end subroutine init_free_wavefunctions


!=========================================================================
!> generate wfn_sym on full k-grid using wfn_input and isym
subroutine gen_wfn_sym(wfn_input, kp, gvec, syms, crys, isym, wfn_sym)

  type(wavefunctions), intent(in) :: wfn_input
  type(kpoints), intent(in) :: kp
  type(gspace), intent(in) :: gvec
  type(symmetry), intent(in) :: syms
  type(crystal), intent(in) :: crys
  integer, intent(in) :: isym
  type(wavefunctions), intent(inout) :: wfn_sym
  !
  integer :: ik_input, ik_sym
  integer :: i, ig, is, ib
  integer :: g0(3)
  real(DP) :: del
  logical :: found
  real(DP), allocatable :: ekin(:)
  integer, allocatable :: isort_sym(:), isortinv_inp(:) ! isort for wfn_sym, and isortinv for wfn_input
  integer, allocatable :: ind(:) ! output of gmap, index map of G vectors
  SCALAR, allocatable :: ph(:) ! output of gmap, phase for each G vector, due to fractional translation
  real(DP) :: skv, sk(3)
  SCALAR :: skv_phase
  type(progress_info) :: prog_info

  PUSH_SUB(gen_wfn_sym)

  wfn_sym%isym = isym  ! this wfn_sym is generated with isym-th symmetry operation

  do ik_sym = 1, kp%nrk
    do i =1, 3
      ! generate wfn_sym%klist by applying the specific isym operation on all wfn_input%klist
      ! generating k_sym(ik) = S * k_input(ik) 
      ! we always move S*k back to the same first BZ (as in kp%rk)
      wfn_sym%klist(i, ik_sym) = dot_product(dble(syms%mtrx(i,:,isym)), kp%rk(:,ik_sym))
    enddo
    wfn_sym%ngk(ik_sym) = wfn_input%ngk(ik_sym) ! they are symmetry related, by construction
    wfn_sym%degeneracy(:,:,ik_sym) = wfn_input%degeneracy(:,:,ik_sym) ! sym-related kpionts have the same band energies
  enddo

  ! find kmap and kg0
  ! mapping between k_sym(ik) and k_input(kmap(ik)):
  !      k_sym_old(ik) + kg0 = k_input(kmap(ik)) 
  ! such that k_sym(ik) is moved back to 1st BZ, corresponding to k_input(kmap(ik)) point
  ! then we udpate
  !      k_sym_new(ik) = k_sym_old(ik) + kg0 = k_input(kmap(ik))
  do ik_sym = 1, kp%nrk
    found = .false.
    ik_input_loop: do ik_input = 1, kp%nrk
      do i = 1, 3
        ! note the definition of del, it must be compatible with gmap
        ! del = A - B
        ! A: wfn at the k to be generated, B: wfn used as input
        ! Here, the rotation in B has been included above, with isym
        del =  kp%rk(i,ik_input) - wfn_sym%klist(i, ik_sym)
        if (del .ge. 0.0d0) g0(i) = del + TOL_Small
        if (del .lt. 0.0d0) g0(i) = del - TOL_Small
        if (abs(del-g0(i)) .gt. TOL_Small) cycle ik_input_loop
      enddo
      found = .true.
      wfn_sym%kmap(ik_sym) = ik_input
      wfn_sym%kg0(:, ik_sym) = g0(:)
      ! now we update klist, to move S*k back to 1st BZ
      ! it is for the convenience for future dot_product to get the gauge transform matrix
      ! because the periodic part of Bloch functions are of the same k
      ! S*k_inp(ik_sym) + G0 = k_inp(ik_input)
      wfn_sym%klist(:, ik_sym) = wfn_sym%klist(:, ik_sym) + g0(:)
      exit ik_input_loop
    enddo ik_input_loop
    if (.not. found) call die('gen_wfn_sym: kpoints mismatch.')  ! it should never happen
  enddo

  ! print information of kmap and kg0 map
  write(*,*)
  write(*,'(1x,a,i2,a)') 'Generating new wavefunctions using symmetry isym = ', isym, ' on input wavefunctions'
  write(*,*)
  write(*,'(1x,a)') 'Mapping relation: k_sym(ik_sym) = k_input(kmap(ik_sym)) = S(isym) * k_input(ik_sym) + G0'
  write(*,'(1x,a,11x,a,11x,a,6x,a,9x,a)') 'ik_sym', 'k_sym(ik_sym)', 'G0 umklapp', 'kmap(ik_sym)->ik_input', 'k_input(ik_sym)'
  do ik_sym = 1, kp%nrk
    write(*,'(1x,i4,3x,a,3(f9.5),1x,a,3x,a,3(i3),1x,a,10x,i6,12x,a,3(f9.5),1x,a)') ik_sym, '(', wfn_sym%klist(:, ik_sym), ')', &
      '(', wfn_sym%kg0(:,ik_sym), ')', wfn_sym%kmap(ik_sym), '(', kp%rk(:,ik_sym),')'
  enddo

  ! generate rotated wavefunctions, following genwf_mpi.f90
  SAFE_ALLOCATE(ekin, (gvec%ng))
  SAFE_ALLOCATE(isort_sym, (gvec%ng))
  SAFE_ALLOCATE(isortinv_inp, (gvec%ng))
  SAFE_ALLOCATE(ind, (gvec%ng))
  SAFE_ALLOCATE(ph, (gvec%ng))

  ! Note:
  !   We are generating a wavefunction at wfn_sym%klist(ik_sym) using wfn_input%klist(ik_sym), because
  !
  !      S * k_input(ik_sym) + G0(ik_sym) = k_inp(ik_inp), with ik_inp = wfn_sym%kmap(ik_sym)
  !
  !   the left-hand side is just wfn_sym%klist(ik_sym), which can be derived using wfn_input%klist(ik_sym)
  !   coefficients with isym and umklapp G0
  !
  call progress_init(prog_info, 'generating wavefunctions using symmetry', 'k-point', kp%nrk)

  do ik_sym = 1, kp%nrk
    call progress_step(prog_info)

    ! we use wfn_input(ik_sym) [NOT ik_input = kmap(ik_sym) !!] combined with rotation + unklapp to generate wfn_sym(ik_sym)
    isortinv_inp(:) = 0
    do ig = 1, wfn_input%ngk(ik_sym)
      isortinv_inp(wfn_input%gindex(ig, ik_sym)) = ig
    enddo

    ! now for k_sym, we sort all (G+k_sym)**2 to find its lowest G vectors up to wfn_sym%ngk(ik_sym), 
    ! or equivalently, wfn_input%ngk(ik_input), where ik_input = wfn_sym%kmap(ik_sym)

    ! first get kinetic energies
    ekin(:) = 0.0d0
    call kinetic_energies(gvec, crys%bdot, ekin, qvec = wfn_sym%klist(:,ik_sym))

    ! then sort and get the sorting order, pointing to global master G-list
    isort_sym(:) = 0
    call sortrx(gvec%ng, ekin, isort_sym, gvec = gvec%components)

    ! save gindex, according to the use in gmap
    wfn_sym%gindex(1:wfn_sym%ngk(ik_sym), ik_sym) = isort_sym(1:wfn_sym%ngk(ik_sym))

    ! call gmap to find out ind and ph, using G vectors in new wfn_sym to find the corresponding (sym-related) G vectors in wfn_input
    ind(:) = 0
    ph(:) = ZERO
    call gmap(gvec, syms, wfn_sym%ngk(ik_sym), isym, wfn_sym%kg0(:,ik_sym), isort_sym, isortinv_inp, ind, ph, .true.)
    ! we MUST include the e^{-iSk.v} phase for gkq unfolding
    sk(:) = matmul(real(syms%mtrx(:,:,isym), DP), wfn_input%klist(:,ik_sym))
    skv = dot_product(sk, syms%tnp(:,isym))
    skv_phase = CMPLX(cos(skv), -sin(skv))
    ph(:) = ph(:) * skv_phase

    ! now we have coefficients for rotated wfn_sym
    do is = 1, kp%nspin*kp%nspinor
      do ib = 1, kp%mnband
        do ig = 1, wfn_input%ngk(ik_sym)
          if (ind(ig) .gt. 0) then
            wfn_sym%cg(ig, ib, is, ik_sym) = wfn_input%cg(ind(ig), ib, is, ik_sym) * ph(ig)
          else
            call die('gen_wfn_sym: some gvectors are not found in symmetry equivalent wavefunctions.')
          endif
        enddo
      enddo
    enddo 

  enddo ! ik_sym main loop

  call progress_free(prog_info)

  POP_SUB(gen_wfn_sym)

  return
 
end subroutine gen_wfn_sym

!=========================================================================
!> calculate the gauge transformation unitary matrix <wfn_sym|wfn_input>, Dmn
subroutine calc_gauge_transform(kp, gvec, wfn_input, wfn_sym)

  type(kpoints), intent(in) :: kp
  type(gspace), intent(in) :: gvec
  type(wavefunctions), intent(in) :: wfn_input
  type(wavefunctions), intent(inout) :: wfn_sym ! one should not pass wfn_input here, to avoid conflicting aliasing
  !
  integer :: ik_sym, ik_input, is, ibm, ibn, ig
  SCALAR, allocatable :: evinp(:), evsym(:) ! single band wavefunction eigenvector, in global master G-list basis
  real(DP), allocatable :: projection(:) ! save summed projection to check overlap
  real(DP) :: max_err ! to report the maxium error deviation from 1.0 in projection
  type(progress_info) :: prog_info

  PUSH_SUB(calc_gauge_transform)

  SAFE_ALLOCATE(evinp, (gvec%ng))
  SAFE_ALLOCATE(evsym, (gvec%ng))

  ! calculating overlap
  call progress_init(prog_info, 'calclating the unitary gauge transformation matrix', 'k-point', kp%nrk)

  is = 1  ! spinless case only
  do ik_sym = 1, kp%nrk
    call progress_step(prog_info)

    ! we need to use kmap to find the corresponding k in the input list
    ik_input = wfn_sym%kmap(ik_sym)

    ! ibn band index for wfn_input
    do ibn = 1, kp%mnband
      ! map wfn_input%cg to evinp in global master G-list basis for wfn dot product
      ! this mapping is needed because wfn_input cg g-vectors have not be sorted,
      ! even if both wfn_input and wfn_sym have the same ngk
      evinp(:) = ZERO
      do ig = 1, wfn_input%ngk(ik_input)
        evinp(wfn_input%gindex(ig, ik_input)) = wfn_input%cg(ig, ibn, is, ik_input)
      enddo
      !
      ! ibm band index for wfn_sym
      do ibm = 1, kp%mnband
        ! make use of the degeneracy information
        ! note that wfn_sym%klist(ik_sym) = wfn_input%klist(ik_input)
        ! they have the same band energies thus degeneracies
        !
        ! if not degenerate, we only need to calculate the diagonal part of the umat_gauge, i.e. ibm=ibn
        ! other off-diagonals are zeros, and are initialized in init_free_wavefunctions
        if ((.not. wfn_input%degeneracy(ibn, is, ik_input)) .and. (ibm .ne. ibn) ) cycle
 
        ! if a band is degenerate, in principle we just need to consider the degenerate subspace
        ! for simplicity, we just loop over all bands for this case
        ! map wfn_sym%cg to evsym in global master G-list basis for wfn dot product
        evsym(:) = ZERO
        do ig = 1, wfn_sym%ngk(ik_sym)
          evsym(wfn_sym%gindex(ig, ik_sym)) = wfn_sym%cg(ig, ibm, is, ik_sym)
        enddo
  
        ! now, we do dot product <wfn_sym(ibm, ik_sym)|wfn_input(ibn, ik_input)
        ! store results in wfn_sym%umat_gauge
        wfn_sym%umat_gauge(ibm, ibn, ik_sym) = zdotc(gvec%ng, evsym(:), 1, evinp(:), 1)
      enddo
    enddo
  enddo

  call progress_free(prog_info)

  ! check completenss relation for overlap
  SAFE_ALLOCATE(projection, (kp%mnband))

  max_err = 0.0d0
  write(*,*) 'Checking projections of wfn_sym, i.e. sum_n(|Dmn|^2), all bands at every k'
  do ik_sym = 1, kp%nrk
    do ibm = 1, kp%mnband
      projection(ibm) = 0.0d0
      do ibn = 1, kp%mnband
        projection(ibm) = projection(ibm) + dble(wfn_sym%umat_gauge(ibm,ibn,ik_sym))**2 & 
                                         + aimag(wfn_sym%umat_gauge(ibm,ibn,ik_sym))**2
      enddo
    enddo
    write(*,'(1x,a,i3,2x,a,*(f6.3))') 'ik_sym =', ik_sym, 'projection = ', projection(:)
    ! compute maximum error
    projection(:) = abs(projection(:) - 1.0d0)
    if (maxval(projection) .gt. max_err) then
      max_err = maxval(projection)
    endif
  enddo
  write(*,'(1x,a,1x,e20.10)') 'The maximum error in the projection (deviation from 1.0) is', max_err
  if (max_err .gt. TOL_Small) then
    write(*,*) 'WARNING: Projection error too large.'
    write(*,*) 'Your choice of nband_use may have break band degeneracy.'
    write(*,*) 'Or your wavefunction quality is low, especially for the highest few bands.'
    write(*,*) 'Try to adjust nband_use, or to compute more bands in WFN for high quality.'
  endif
  write(*,*)

  max_err = 0.0d0
  write(*,*) 'Checking projections of wfn_input, i.e. sum_m(|Dmn|^2), all bands at every k'
  do ik_sym = 1, kp%nrk
    do ibn = 1, kp%mnband
      projection(ibn) = 0.0d0
      do ibm = 1, kp%mnband
        projection(ibn) = projection(ibn) + dble(wfn_sym%umat_gauge(ibm,ibn,ik_sym))**2 & 
                                         + aimag(wfn_sym%umat_gauge(ibm,ibn,ik_sym))**2
      enddo
    enddo
    write(*,'(1x,a,i3,2x,a,*(f6.3))') 'ik_inp =', wfn_sym%kmap(ik_sym), 'projection = ', projection(:)
    ! compute maximum error
    projection(:) = abs(projection(:) - 1.0d0)
    if (maxval(projection) .gt. max_err) then
      max_err = maxval(projection)
    endif
  enddo
  write(*,'(1x,a,1x,e20.10)') 'The maximum error in the projection (deviation from 1.0) is', max_err
  if (max_err .gt. TOL_Small) then
    write(*,*) 'WARNING: Projection error too large.'
    write(*,*) 'Your choice of nband_use may have break band degeneracy.'
    write(*,*) 'Or your wavefunction quality is low, especially for the highest few bands.'
    write(*,*) 'Try to adjust nband_use, or to compute more bands in WFN for high quality.'
  endif
  write(*,*)

  SAFE_DEALLOCATE(evinp)
  SAFE_DEALLOCATE(evsym)
  SAFE_DEALLOCATE(projection)

  POP_SUB(calc_gauge_transform)

  return

end subroutine calc_gauge_transform

!=========================================================================
!> write the gauge transformation unitary matrix <wfn_sym|wfn_input>, Dmn, to file
subroutine write_umat_gauge(kp, wfn)

  type(kpoints), intent(in) :: kp
  type(wavefunctions), intent(in) :: wfn
  !
  integer :: ik, ibm, ibn, nblock, iblock, i, ibsave(10), nremain, ib_end
  character(len=256) :: ufile, fmtline

  PUSH_SUB(write_umat_gauge)

  write(ufile, '(a,i0,a)') 'umat_gauge.isym', wfn%isym,'.dat'

  call open_file(unit=53, file=trim(ufile), form='formatted', status='replace')

  write(*,'(1x,a,1x,a)') 'Writing umat_gauge to file', trim(ufile)
  write(*,*)

  ! write umat_gauge matrix as 10 columns per block
  nblock = kp%mnband / 10
  nremain = mod(kp%mnband, 10)
  if (nremain .gt. 0) then
    nblock = nblock + 1  ! there are some remainder columns
  endif

  do ik = 1, kp%nrk
    write(53, '(a,i5,a,3(f10.6),a)') 'ik = ', ik, ' , kpoint = ( ', wfn%klist(:,ik), ' )'

    ! write umat_gauge in blocks
    do iblock = 1, nblock
      if ((iblock .eq. nblock) .and. (nremain .gt. 0)) then
        ib_end = nremain
      else
        ib_end = 10
      endif

      ibsave(:) = 0
      do i = 1, ib_end
        ibn = (iblock-1) * 10 + i
        ibsave(i) = ibn
      enddo

      write(53, '(2x,a,10i16)') 'ibm\ibn', ibsave(1:ib_end)

      ! each block writes full row range of ibm
      do ibm = 1, kp%mnband
        ! complex number need two data descriptors
        ! write the format line first
        write(fmtline, '(a,i0,a)') "(i5,12x,", ib_end, "('('f6.3,',',f6.3,') '))"
        write(53, fmtline) ibm, wfn%umat_gauge(ibm, ibsave(1):ibsave(ib_end) ,ik)
      enddo
    enddo ! iblock
    write(53,*)
  enddo

  call close_file(53)

  POP_SUB(write_umat_gauge)

  return

end subroutine write_umat_gauge

!=========================================================================
!> find maps between symmetries and their inverse symmetries
subroutine find_maps_sym_invsym(syms, sym_to_invsym, invsym_to_sym)

  type(symmetry), intent(in) :: syms
  integer, allocatable, intent(inout) :: sym_to_invsym(:)
  integer, allocatable, intent(inout) :: invsym_to_sym(:)
  !
  integer :: isym, iinvsym
  real(DP) :: sym_mtrx_k(3,3), sym_inv_mtrx_r(3,3) ! symmetry matrix in k and r spaces
  real(DP) :: mtrx_aux(3,3) ! auxiliary matrix for general purpose
  real(DP) :: ftrans(3), invftrans(3) ! fractional translations in crystal lattice basis
  real(DP), allocatable :: syms_inv_mtrx(:,:,:) ! store all inverse matrices of syms%mtrx
  real(DP), allocatable :: syms_inv_tnp(:,:) ! store all inverse ftrans of syms%tnp
  real(DP) :: symdiff
  integer :: i, j

  PUSH_SUB(find_maps_sym_invsym)

  SAFE_ALLOCATE(syms_inv_mtrx, (3,3,syms%ntran))
  SAFE_ALLOCATE(syms_inv_tnp, (3,syms%ntran))

  ! first compute and store the inverse symmetry operation for each isym
  do isym = 1, syms%ntran
    sym_mtrx_k(:,:) = real(syms%mtrx(:, :, isym), DP)
    ! first invert the matrix to get rotation inverse S^{-1}, directly, in k-basis
    call invert_matrix(sym_mtrx_k(:,:), syms_inv_mtrx(:,:,isym))

    ! then we find inv_t = -S^{-1} t, we need S^{-1} in real space to do this
    ! convert syminv matrix in k space to r space in crystal lattices basis: related as transpose + inverse
    mtrx_aux(:,:) = transpose(syms_inv_mtrx(:,:,isym))
    call invert_matrix(mtrx_aux, sym_inv_mtrx_r)
  
    ! convert BGW fractional translation to r space in crystal lattice basis
    ftrans(:) = syms%tnp(:, isym) / 2.0d0 / PI_D
    invftrans(:) = -matmul(sym_inv_mtrx_r(:,:), ftrans(:))

    syms_inv_tnp(:,isym) = invftrans(:) * 2.0d0 * PI_D
  enddo

  sym_to_invsym(:) = 0
  invsym_to_sym(:) = 0

  ! now we find the maps
  do isym = 1, syms%ntran
    ! for each isym, find the index in sym%ntran that which iinvsym is isym`s inverse correspondence
    iinvsym_loop: do iinvsym = 1, syms%ntran
      symdiff = 0.0d0

      sym_mtrx_k(:,:) = real(syms%mtrx(:, :, iinvsym), DP)
      mtrx_aux(:,:) = syms_inv_mtrx(:,:,isym) - sym_mtrx_k(:,:)

      do i = 1, 3
        do j = 1, 3
          symdiff = symdiff + abs(mtrx_aux(i,j))
        enddo
      enddo

      ftrans(:) = (syms_inv_tnp(:,isym) - syms%tnp(:,iinvsym)) / 2.0d0 / PI_D
      do i = 1, 3
        symdiff = symdiff + abs(ftrans(i) - nint(ftrans(i)))
      enddo

      !write(*,*) 'DEBUG:', isym, iinvsym, symdiff

      if (symdiff .lt. TOL_Small) then
        sym_to_invsym(isym) = iinvsym
        invsym_to_sym(iinvsym) = isym
        exit iinvsym_loop
      endif
    enddo iinvsym_loop
  enddo

  write(*,*)
  write(*,*) 'Maps from symmetries to inverse symmetries'
  do isym = 1, syms%ntran
    write(*,'(3x,a,i3,a,i3,a)') 'The', isym, '-th symmetry has its inverse operation as the', sym_to_invsym(isym), '-th symmetry'
  enddo
  write(*,*)
  write(*,*) 'Maps from inverse symmetries to symmetries'
  do iinvsym = 1, syms%ntran
    write(*,'(3x,a,i3,a,i3,a)') 'The', iinvsym, '-th symmetry is the inverse operation of the', invsym_to_sym(iinvsym), '-th symmetry'
  enddo
  write(*,*)

  SAFE_DEALLOCATE(syms_inv_mtrx)
  SAFE_DEALLOCATE(syms_inv_tnp)

  if (any(sym_to_invsym(:) .eq. 0) .or. any(invsym_to_sym(:) .eq. 0)) then
    call die('Symmetries do not form a group. Inverse symmetries are not properly identified.')
  endif

  POP_SUB(find_maps_sym_invsym)

  return

end subroutine find_maps_sym_invsym

!=========================================================================
!> given a kpoint coordinate and a list of k, return its index in the list
!> if two kpoints are differed by an umklapp G vector, they are considered the same kpoint
integer function kindex(kpt, klist, nk)
  real(DP), intent(in) :: kpt(3)
  real(DP), pointer, intent(in) :: klist(:,:)
  integer, intent(in) :: nk
  !
  integer :: ik
  real(DP) :: kdiff(3), knorm

  ! no PUSH/POP because called frequently

  kindex = 0  ! initialize with an invalid value

  do ik = 1, nk
    kdiff(:) = kpt(:) - klist(:, ik)
    kdiff(:) = kdiff(:) - nint(kdiff(:))
    knorm = sqrt(dot_product(kdiff, kdiff))

    if (knorm .lt. TOL_Small) then
      kindex = ik
      exit
    endif
  enddo

  if (kindex .eq. 0) then
    call die('kindex: kpt not found in provided klist.')
  endif

  return

end function kindex

!#END_INTERNAL_ONLY

#endif

end module sympert_utils_m
