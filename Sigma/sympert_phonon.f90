!=============================================================================
! Phonon relevant subroutines for sympert.f90
!
! Originally by Zhenglu Li
!=============================================================================

#include "f_defs.h"

module sympert_phonon_m
 
#ifdef CPLX

  use global_m
  use misc_m
  use input_utils_m
  use io_utils_m

  implicit none

  private

  public :: &
    perturbations, &
    gkqmatrix, &
    init_free_perturbations, &
    init_read_free_gkqmatrix, &
    find_equiv_site, &
    gen_unfold_q_sym, &
    write_gkqmatrix

  type perturbations
    integer :: iqfbz
    real(DP) :: qpoint(3)
    integer :: iqirr ! sym-equiv. iqirr using isym
    integer :: isym ! i-th symmetry operation in syms
    integer, pointer :: apos_map(:)
  end type perturbations

  type gkqmatrix
    integer :: nq
    integer :: qmesh(3) ! define the full mesh to be unfoled to
    real(DP), pointer :: qlist(:,:)
    integer :: nk
    real(DP), pointer :: klist(:,:) ! kpoints where gkq are calculated
    integer :: nband ! number of bands for gkq
    SCALAR, pointer :: gkq(:,:,:,:,:,:) ! e-ph matrix elements, size (nband_k, nk, ndir, natom, nband_kq, nq)
  end type gkqmatrix

contains

!#BEGIN_INTERNAL_ONLY

!=========================================================================
!> initialize or free perturbations object
subroutine init_free_perturbations(switch, phonpert, crys)

  integer, intent(in) :: switch
  type(perturbations), intent(inout) :: phonpert
  type(crystal), optional, intent(in) :: crys

  PUSH_SUB(init_free_perturbations)

  if (switch .eq. 1) then
    if (.not. present(crys)) then
      call die('You need to pass crys into init_free_perturbations to initialize phonpert object.')
    endif
    SAFE_ALLOCATE(phonpert%apos_map, (crys%nat))

    phonpert%iqfbz = 0
    phonpert%qpoint(:) = 0.0d0
    phonpert%iqirr = 0
    phonpert%isym = 0
    phonpert%apos_map(:) = 0
  elseif (switch .eq. 0) then
    SAFE_DEALLOCATE_P(phonpert%apos_map)
    phonpert%isym = 0
  else
    call die('switch in init_free_perturbations can only be 1 (init) or 0 (free).')
  endif

  POP_SUB(init_free_perturbations)

  return

end subroutine init_free_perturbations

!=========================================================================
!> initialize, read or free gkqmat object
subroutine init_read_free_gkqmatrix(switch, gkqmat, nq_flag, crys, kp, gkqfile)

  integer, intent(in) :: switch  ! 1: initialize, 0: free
  type(gkqmatrix), intent(inout) :: gkqmat
  integer, optional, intent(in) :: nq_flag  ! -1, read from gkq.mat, >0, initialize with zero
  type(crystal), optional, intent(in) :: crys
  type(kpoints), optional, intent(in) :: kp
  character(len=256), optional, intent(in) :: gkqfile
  !
  integer :: nq, iq, ik
  integer :: nb_k, nk, ndir, natom, nb_kq! dimension of matrix for each iq
  integer :: qmesh(3)
  type(progress_info) :: prog_info

  PUSH_SUB(init_read_free_gkqmatrix)

  if (switch .eq. 1) then
    ! nq_flag = -1: read from file gkq.mat
    ! nq_flag > 0 : number of q points to be initialized, with zero
    if (nq_flag .eq. -1) then
      if ((.not. present(nq_flag)) .or. (.not. present(crys)) .or. (.not. present(gkqfile))) then
        call die('You need to pass nq_flag, crys and gkqfile to init_free_gkqmatrix to initialize gkq object from file.')
      endif

      ! read from gkq.mat file
      write(*,*) 'Reading gkqmat from file ', trim(gkqfile)
  
      call open_file(unit=55, file=trim(gkqfile), form='formatted', status='old')
  
      read(55,*) nq
      read(55,*) qmesh(:)
      read(55,*) nb_k, nk, ndir, natom, nb_kq
  
      write(*,*)
      write(*,'(1x,a,i6)') 'number of q-points in file:', nq
      write(*,'(1x,a,3i5)') 'uniform q-mesh grid:', qmesh
      write(*,'(1x,a,1x,5i6)') 'gkq matrix (of each q) dimension: nb_k, nk, ndir, natom, nb_kq     <=>', nb_k, nk, ndir, natom, nb_kq

      if ((ndir.ne.3) .or. (natom.ne.crys%nat)) then
        call die('init_read_free_gkqmatrix: wrong atom perturbation dimension in gkqfile.')
      endif

      if (nb_k .ne. nb_kq) then
        call die('Number of bands on k and k+q must be equal.')
      endif

      if (nq .gt. qmesh(1)*qmesh(2)*qmesh(3)) then
        call die('Number of listed q-points is larger than number of points in q-mesh.')
      endif

      gkqmat%nq = nq
      gkqmat%qmesh(:) = qmesh(:)
      gkqmat%nk = nk
      gkqmat%nband = nb_k  ! we have checked nb_k .eq. nb_kq

      ! read kpoints list on which gkq is calculated
      SAFE_ALLOCATE(gkqmat%klist, (3, nk))
      do ik = 1, nk
        read(55,*) gkqmat%klist(:, ik)
      enddo

      ! read qlist and gkq matrix
      SAFE_ALLOCATE(gkqmat%qlist, (3, nq))
      SAFE_ALLOCATE(gkqmat%gkq, (nb_k, nk, 3, natom, nb_kq, nq))
  
      call progress_init(prog_info, 'reading e-ph matrix elements', 'q-point', nq)
  
      do iq = 1, nq
        call progress_step(prog_info)
  
        read(55,*) gkqmat%qlist(:, iq)
        read(55,*) gkqmat%gkq(:, :, :, :, :, iq)
      enddo
      
      call progress_free(prog_info)
  
      call close_file(55)
  
      write(*,*) trim(gkqfile), ' has been successfully read in.'
      write(*,*)

      write(*,*) 'qpoints list (irreducible wedge) of input gkq e-ph matrix elements:'
      do iq = 1, gkqmat%nq
        write(*,'(3x,a,i5,1x,a,3(f9.5))') 'q(', iq, ') = ', gkqmat%qlist(:,iq)
      enddo
      write(*,*)

    elseif (nq_flag .gt. 0) then
      if ((.not. present(nq_flag)) .or. (.not. present(crys)) .or. (.not. present(kp))) then
        call die('You need to pass nq_flag, crys and kp to init_free_gkqmatrix to initialize gkq object with zero.')
      endif
  
      nq = nq_flag
      SAFE_ALLOCATE(gkqmat%qlist, (3, nq))
      SAFE_ALLOCATE(gkqmat%klist, (3, kp%nrk))
      SAFE_ALLOCATE(gkqmat%gkq, (kp%mnband, kp%nrk, 3, crys%nat, kp%mnband, nq))
  
      gkqmat%nq = nq
      gkqmat%qmesh(:) = 0
      gkqmat%qlist(:,:) = 0.0d0
      gkqmat%nk = kp%nrk
      gkqmat%klist(:,:) = kp%rk(:,:)
      gkqmat%nband = kp%mnband
      gkqmat%gkq(:,:,:,:,:,:) = ZERO
    else
      call die('init_read_free_gkqmatrix: nq_flag can only be -1 (read from file) and positive (initialize with zero).')
    endif ! nq_flag
  elseif (switch .eq. 0) then
    SAFE_DEALLOCATE_P(gkqmat%qlist)
    SAFE_DEALLOCATE_P(gkqmat%klist)
    SAFE_DEALLOCATE_P(gkqmat%gkq)

    gkqmat%nq = 0
    gkqmat%qmesh(:) = 0
    gkqmat%nk = 0
    gkqmat%nband = 0
  else
    call die('switch in init_read_free_gkqmatrix can only be 1 (init) or 0 (free).')
  endif

  POP_SUB(init_read_free_gkqmatrix)

  return

end subroutine init_read_free_gkqmatrix

!=========================================================================
!> find equivalent atom sites for a given symmetry isym
subroutine find_equiv_site(crys, syms, isym, phonpert)

  type(symmetry), intent(in) :: syms
  type(crystal), intent(in) :: crys
  integer, intent(in) :: isym
  type(perturbations), intent(inout) :: phonpert
  !
  real(DP) :: sym_mtrx_k(3,3), sym_mtrx_r(3,3) ! symmetry matrix in k and r spaces
  real(DP) :: mtrx_aux(3,3) ! auxiliary matrix for general purpose
  real(DP) :: ftrans(3) ! fractional translation in crystal lattice basis
  real(DP), allocatable :: apos_inp(:,:) ! atom positions from input
  real(DP), allocatable :: apos_sym(:,:) ! atom positions derived from symmetry operation isym
  real(DP) :: apos_diff(3)
  integer :: ia_inp, ia_sym, i
  logical :: found

  PUSH_SUB(find_equiv_site)

  write(*,'(1x,a,1x,i3)') 'Looking for equivalent atom sites in crystal using isym =', isym

  ! convert BGW sym matrix in k space to r space in crystal lattices basis
  sym_mtrx_k(:,:) = real(syms%mtrx(:, :, isym), DP)
  mtrx_aux(:,:) = transpose(sym_mtrx_k)
  call invert_matrix(mtrx_aux, sym_mtrx_r)

  ! convert BGW fractional translation to r space in crystal lattice basis
  ftrans(:) = syms%tnp(:, isym) / 2.0d0 / PI_D

  phonpert%apos_map(:) = 0

  SAFE_ALLOCATE(apos_inp, (3, crys%nat))
  SAFE_ALLOCATE(apos_sym, (3, crys%nat))

  do ia_inp = 1, crys%nat
    ! convert apos_inp from cartesian to fractional coordinate, and save all positions
    apos_inp(:, ia_inp) = matmul(transpose(crys%bvec), crys%apos(:,ia_inp))
    ! find all positions generated from apos_inp by isym
    ! the ia_inp-th atom coordinate is rotated and fractional translated to generate apos_sym
    apos_sym(:, ia_inp) = matmul(sym_mtrx_r(:,:), apos_inp(:, ia_inp))
    apos_sym(:, ia_inp) = apos_sym(:, ia_inp) + ftrans(:)
  enddo

  ! it is very important to define apos_map properly
  ! the definition of transformation of atom coordinates directly determins how perturbed potentials transform
  ! the relation we use is,
  !  x_KAPPA = S * x_kappa + v, up to real-space crystal lattice vectors
  ! where x represents atom position, KAPPA labels the new atom position derived from symmetry operation on kappa
  !
  ! here we define the apos_map as: for the same atom coordinates as input,
  ! the ia_sym-th atom (x_KAPPA) can be GENERATED from apos_map(ia_sym)-th atom (x_kappa) using isym-th symmetry
  !
  ! this way it is straightforward to generate the output gkqmat_output, that we loop over KAPPA, and look for 
  ! equivalent sites to generate the unknown gkqmat
  !
  do ia_inp = 1, crys%nat
    found = .false.
    ia_sym_loop: do ia_sym = 1, crys%nat
      apos_diff(:) = apos_inp(:,ia_inp) - apos_sym(:,ia_sym)

      if ((abs(apos_diff(1) - nint(apos_diff(1))) < TOL_Small) .and. &
          (abs(apos_diff(2) - nint(apos_diff(2))) < TOL_Small) .and. &
          (abs(apos_diff(3) - nint(apos_diff(3))) < TOL_Small)) then
        found = .true.
        phonpert%apos_map(ia_inp) = ia_sym
        exit ia_sym_loop
      endif
    enddo ia_sym_loop

    if (.not. found) call die('find_equiv_site: no equivalent atom site found. You may need to increase structure accuracy.')
  enddo

  ! print apos_map information
  write(*,*) 'Symmetry equivalent site information:'
  do ia_inp = 1, crys%nat
    write(*,'(3x,a,i3,1x,a,i3,1x,a)') 'apos_inp(', ia_inp, ')  <=>  S * apos_inp(', phonpert%apos_map(ia_inp), ') + ftrans'
  enddo
  write(*,*)

  SAFE_DEALLOCATE(apos_inp)

  POP_SUB(find_equiv_site)

  return

end subroutine find_equiv_site

!=========================================================================
!> generate full qlist in full BZ, and unique symmetry isym for each 
!> irredicible iq to unfold to full BZ
subroutine gen_unfold_q_sym(gkqmat, syms, qlist_fbz, symq_map)

  type(gkqmatrix), intent(in) :: gkqmat
  type(symmetry), intent(in) :: syms
  real(DP), allocatable, intent(inout) :: qlist_fbz(:,:)
  integer, allocatable, intent(inout) :: symq_map(:,:)
  !
  integer :: nqirr, nqfbz  ! number of q points in irreducible and full BZ
  integer :: iqirr, iqfbz, isym, iq1, iq2, iq3
  integer :: qmesh(3)
  real(DP) :: qdiff(3), qnorm
  logical :: found

  PUSH_SUB(gen_unfold_q_sym)

  nqirr = gkqmat%nq
  qmesh(:) = gkqmat%qmesh(:)
  nqfbz = qmesh(1) * qmesh(2) * qmesh(3)

  ! generate q points in full BZ, (-0.5, 0.5]
  iqfbz = 0
  do iq1 = 1, qmesh(1)
    do iq2 = 1, qmesh(2)
      do iq3 = 1, qmesh(3)
        iqfbz = iqfbz + 1
        qlist_fbz(1,iqfbz) = real((iq1-1), DP)/qmesh(1)
        qlist_fbz(2,iqfbz) = real((iq2-1), DP)/qmesh(2)
        qlist_fbz(3,iqfbz) = real((iq3-1), DP)/qmesh(3)

        if (qlist_fbz(1,iqfbz) .gt. (0.5d0+TOL_Small)) qlist_fbz(1,iqfbz) = qlist_fbz(1,iqfbz) - 1.0d0
        if (qlist_fbz(2,iqfbz) .gt. (0.5d0+TOL_Small)) qlist_fbz(2,iqfbz) = qlist_fbz(2,iqfbz) - 1.0d0
        if (qlist_fbz(3,iqfbz) .gt. (0.5d0+TOL_Small)) qlist_fbz(3,iqfbz) = qlist_fbz(3,iqfbz) - 1.0d0
      enddo
    enddo
  enddo

  ! loop over symmetries to find a unique (iqirr, isym) pair to generate an iqfbz
  ! we do not need to track the umklapp vector because:
  !   d_q V(r_{BvK}) = d_{q+G} V(r_{BvK}), e^{i G \cdot R} = 1
  do iqfbz = 1, nqfbz
    found = .false.

    iqirr_loop: do iqirr = 1, nqirr
      isym_loop: do isym = 1, syms%ntran
        qdiff(:) = qlist_fbz(:, iqfbz) - matmul(real(syms%mtrx(:,:,isym), DP), gkqmat%qlist(:, iqirr))
        qdiff(:) = qdiff(:) - nint(qdiff(:))
        qnorm = sqrt(dot_product(qdiff, qdiff))

        !write(*,*) 'DEBUG:', iqfbz, iqirr, isym
        !write(*,*) 'q0', gkqmat%qlist(:, iqirr)
        !write(*,*) 'Sq0', matmul(real(syms%mtrx(:,:,isym), DP), gkqmat%qlist(:, iqirr))
        !write(*,*) qdiff, qnorm

        if (qnorm .lt. TOL_Small) then
          found = .true.

          symq_map(1, iqfbz) = iqirr ! a pair of (iqirr, isym) is found
          symq_map(2, iqfbz) = isym  ! isym acting on iqirr qpoint will generate iqfbz qpoint

          exit iqirr_loop ! then we exit the search loops for iqfbz to make this pair unique
        endif
      enddo isym_loop
    enddo iqirr_loop

    if (.not. found) then
      call die('Input irreducible qpoint cannot be mapped to full qmesh with any symmetries.')
    endif
  enddo

  write(*,*)
  write(*,*) 'Generation of full q-points mesh BZ using irreducible q-wedge and symmetries'
  do iqfbz = 1, nqfbz
    write(*,'(1x,a,i3,a,3(f9.5),a,i3,1x,a,i3)') 'q_full_BZ(', iqfbz, ') = (', qlist_fbz(:, iqfbz), &
                ') is generated from irreducible iqirr =', symq_map(1, iqfbz), 'and isym =', symq_map(2, iqfbz)
  enddo

  POP_SUB(gen_unfold_q_sym)

  return
end subroutine gen_unfold_q_sym

!=========================================================================
!> write gkqmat object to output file
subroutine write_gkqmatrix(gkqmat, crys, gkqfile)

  type(gkqmatrix), intent(in) :: gkqmat
  type(crystal), intent(in) :: crys
  character(len=256), intent(in) :: gkqfile
  !
  integer :: iq, ik, ib_k, ib_kq, idir, iatom
  type(progress_info) :: prog_info

  PUSH_SUB(write_gkqmatrix)

  ! write gkq.mat to file
  write(*,*) 'Writing gkqmat to file ', trim(gkqfile)

  call open_file(unit=56, file=trim(gkqfile), form='formatted', status='replace')

  ! write dimensions
  write(56,'(i0)') gkqmat%nq
  write(56,'(3i6)') gkqmat%qmesh(:)
  write(56,'(5i6)') gkqmat%nband, gkqmat%nk, 3, crys%nat, gkqmat%nband

  write(*,*)
  write(*,'(1x,a,i5,1x,a)') 'Writing unfolded gkq matrix elements with', gkqmat%nq, 'q-points.'

  ! write klist
  do ik = 1, gkqmat%nk
    write(56,'(3(f12.8))') gkqmat%klist(:, ik)
  enddo

  ! progress report
  call progress_init(prog_info, 'writing e-ph matrix elements', 'q-point', gkqmat%nq)

  ! write main iq loop
  do iq = 1, gkqmat%nq
    call progress_step(prog_info)

    ! write q point
    write(56,'(3(f12.8))') gkqmat%qlist(:, iq)

    ! write gkq matrix

    do ib_kq = 1, gkqmat%nband
      do iatom = 1, crys%nat
        do idir = 1, 3
          do ik = 1, gkqmat%nk
            do ib_k = 1, gkqmat%nband
              write(56,*) gkqmat%gkq(ib_k, ik, idir, iatom, ib_kq, iq)  ! default Fortran complex number format
            enddo ! ib_k
          enddo ! ik
        enddo ! idir
      enddo ! iatom
    enddo ! ib_kq
  enddo ! iq

  call progress_free(prog_info)

  call close_file(56)

  POP_SUB(write_gkqmatrix)

  return

end subroutine write_gkqmatrix

!#END_INTERNAL_ONLY

#endif

end module sympert_phonon_m
