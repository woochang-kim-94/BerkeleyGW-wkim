!===============================================================================
!
! Program:
!
! (1) compare_wfns          Originally By FHJ      Last Modified May/2021 (FHJ)
!
!===============================================================================

#include "f_defs.h"

program compare_wfns
#ifdef HDF5
  use hdf5
  use hdf5_io_m
  use input_utils_m,    only: gvec_index
  use misc_m,           only: findvector
  implicit none

  ! Input
  character(len=256) :: sarg, fname_ref, fname_new
  integer :: nargs, iarg, nb, err, ios
  real(DP) :: tol_deg
  logical :: terse, no_deg
  ! General
  integer, allocatable :: dims(:), dims2(:)
  integer(hid_t) :: file_id_ref, file_id_new
  real(DP) :: err_evals, err_ortho1, err_ortho2, err_cross, err_ham
  ! Headers
  integer :: nk, nk2, ns, ns2, nb1, nb2
  integer, allocatable :: ngk(:), ngk2(:)
  ! Energies
  real(DP), allocatable :: en1(:,:,:), en2(:,:,:) !(nb, rk, nspin)
  ! Wavefunctions
  SCALAR, allocatable :: wfns1(:,:,:), wfns2(:,:,:) !(ngktot, nspin*nspinor, nb)


  !------------
  ! Parse input
  fname_ref = ''
  fname_new = ''
  terse = .false.
  nb = -1
  tol_deg = 1d-10
  no_deg = .false.

  nargs = command_argument_count()
  iarg = 1
  do while (iarg<=nargs)
    call get_command_argument(iarg, sarg, status=err)
    if (err/=0) call die_usage()

    select case(trim(sarg))
      case ('--terse')
        terse = .true.
        continue
      case ('--no_deg')
        no_deg = .true.
        continue
      case ('--nb')
        iarg = iarg + 1
        call get_command_argument(iarg, sarg, status=err)
        if (err/=0) call die_usage('nb')
        read(sarg,*,iostat=ios) nb
        if (ios/=0) call die_usage('nb')
      case ('--tol_deg')
        iarg = iarg + 1
        call get_command_argument(iarg, sarg, status=err)
        if (err/=0) call die_usage('tol_deg')
        read(sarg,*,iostat=ios) tol_deg
        if (ios/=0) call die_usage('tol_deg')
      case default
        if (len_trim(fname_ref)==0) then
          fname_ref = trim(sarg)
        elseif (len_trim(fname_new)==0) then
          fname_new = trim(sarg)
        else
          call die_usage()
        endif
    endselect
    iarg = iarg + 1
  enddo

  if (len_trim(fname_ref)==0 .and. len_trim(fname_new)==0) call die_usage()

  call h5open_f(err)

  call hdf5_open_file(trim(fname_ref), 'r', file_id_ref)
  call hdf5_open_file(trim(fname_new), 'r', file_id_new)

  !--------------
  ! Check headers
  call hdf5_read_int(file_id_ref, 'mf_header/kpoints/nrk', nk)
  call hdf5_read_int(file_id_new, 'mf_header/kpoints/nrk', nk2)
  if (nk/=nk2) call die('number of k-point mismatch')

  call hdf5_read_int(file_id_ref, 'mf_header/kpoints/nspin', ns)
  call hdf5_read_int(file_id_new, 'mf_header/kpoints/nspin', ns2)
  if (nk/=nk2) call die('number of spin components mismatch')

  SAFE_ALLOCATE(ngk, (nk))
  SAFE_ALLOCATE(ngk2, (nk))
  call hdf5_read_int_array(file_id_ref, 'mf_header/kpoints/ngk', [nk], ngk)
  call hdf5_read_int_array(file_id_new, 'mf_header/kpoints/ngk', [nk], ngk2)
  if (any(ngk/=ngk2)) call die('number of G-vectors per k-point mismatch')

  call hdf5_read_int(file_id_ref, 'mf_header/kpoints/mnband', nb1)
  call hdf5_read_int(file_id_new, 'mf_header/kpoints/mnband', nb2)
  if (nb == -1) then
    nb = min(nb1, nb2)
  else
    if (min(nb1,nb2) < nb) call die('Input parameter nb too large')
  endif


  !--------------------
  ! Compare eigenvalues
  ! 1) Check that the eigenvalues match

  dims = get_dims(file_id_ref, 'mf_header/kpoints/el')
  dims(1) = nb
  SAFE_ALLOCATE(en1, (dims(1),dims(2),dims(3)))
  call hdf5_read_double_hyperslab(file_id_ref, 'mf_header/kpoints/el', dims, [0,0,0], en1)

  dims = get_dims(file_id_new, 'mf_header/kpoints/el')
  dims(1) = nb
  SAFE_ALLOCATE(en2, (dims(1),dims(2),dims(3)))
  call hdf5_read_double_hyperslab(file_id_new, 'mf_header/kpoints/el', dims, [0,0,0], en2)

  if (any(shape(en1)/=shape(en2))) call die('eigenvalue shape mismatch')
  err_evals = maxval(abs(en1 - en2))


  !-------------------
  ! Read wavefunctions
  dims2 = get_dims(file_id_ref, 'wfns/coeffs')
  dims = dims2(2:) ! Ignore the real/complex part
  dims(3) = nb
  SAFE_ALLOCATE(wfns1, (dims(1),dims(2),dims(3)))
  call hdf5_read_scalar_hyperslab(file_id_ref, 'wfns/coeffs', dims, [0,0,0], wfns1)

  dims2 = get_dims(file_id_new, 'wfns/coeffs')
  dims = dims2(2:) ! Ignore the real/complex part
  dims(3) = nb
  SAFE_ALLOCATE(wfns2, (dims(1),dims(2),dims(3)))
  call hdf5_read_scalar_hyperslab(file_id_new, 'wfns/coeffs', dims, [0,0,0], wfns2)

  if (any(shape(wfns1)/=shape(wfns1))) call die('wavefunction shape mismatch')


  !----------------------
  ! Compare wavefunctions
  call check_wfns()

  !-------
  ! Output

  if (terse) then
    write(6,'(es10.3)') err_evals
    write(6,'(es10.3)') err_ortho1
    write(6,'(es10.3)') err_ortho2
    write(6,'(es10.3)') err_cross
    write(6,'(es10.3)') err_ham
  else
    write(6,'(a)') 'Maximum error in:'
    write(6,'(a,es10.3,a)') '- eigenvalues |E_new - E_ref|: ', err_evals, ' Ry'
    write(6,'(a,es10.3)') '- ref. wavefunctions orthonormality |C_ref^H C_ref - 1|: ', err_ortho1
    write(6,'(a,es10.3)') '- new wavefunctions orthonormality |C_new^H C_new - 1|: ', err_ortho2
    write(6,'(a,es10.3)') '- cross wavefunctions orthonormality |C_new^H C_new - C_ref^H C_ref|: ', err_cross
    write(6,'(a,es10.3)') '- reassembled Hamiltonian |H_ref - H_new|/(|H_ref||H_new|): ', err_ham
  endif

  call h5close_f(err)

contains

  subroutine die_usage(sarg)
    character(len=*), optional, intent(in) :: sarg

    if (present(sarg)) then
      write(0,'(a,a/)') 'Error parsing parameter for flag ', trim(sarg)
    endif

    write(0,'(100(a/))') &
'usage: compare_wfns.x [--terse] [--nb NB] [--tol_deg TOL_DEG] [--no_deg] fname_ref fname_new', &
'',&
'positional arguments:',&
'  fname_ref          Reference WFN.h5 file.',&
'  fname_new          New WFN.h5 file to check against reference.',&
'',&
'optional arguments:',&
'  --terse            Return raw numbers without any additional explanation.',&
'  --nb NB            Number of bands to compare. Defaults to all.',&
'  --tol_deg TOL_DEG  Tolerance for degenerate states, in Ry.',&
'  --no_deg           Do not try to find a non-degenerate subspace. Only use this if the WFN does not slice a degenerate subspace'
    stop

  end subroutine die_usage

  function get_dims(file_id, path) result(dims)
    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: path
    integer, allocatable :: dims(:)

    integer(hid_t) :: dset_id
    integer(hid_t) :: dataspace_id
    integer(hsize_t), allocatable :: dims_(:), maxdims_(:)
    integer :: rank, err

    call safe_h5dopen(file_id, path, dset_id)
    call safe_h5dget_space(dset_id, dataspace_id)
    call h5sget_simple_extent_ndims_f(dataspace_id, rank, err)
    if (err/=0) call die('Error calling h5sget_simple_extent_ndims_f')
    SAFE_ALLOCATE(dims_, (rank))
    SAFE_ALLOCATE(maxdims_, (rank))
    call h5sget_simple_extent_dims_f(dataspace_id, dims_, maxdims_, err)
    if (err<0) call die('Error calling h5sget_simple_extent_dims_f')
    SAFE_ALLOCATE(dims, (rank))
    dims(:) = int(dims_(:))
    call safe_h5sclose(dataspace_id)
    call safe_h5dclose(dset_id)

  end function get_dims

  function get_norm(x) result(norm)
    SCALAR, intent(in) :: x(:,:)
    real(DP) :: norm

    norm = sqrt(sum(MYABS2(x)))

  end function get_norm

  ! We make a few consistency checks, where C = C(ig,ib), i.e. Fortran convention:
  ! 2-3) Check that each WFN file is orthonormal to within itself, C* C ~ 1,
  ! 4) Check that both WFN file spawn the same identity, C2* C2 - C1* C1.
  !    This is useful for pseudobands.
  ! 5) Check that the reassembled Hamiltonian is the same, where the
  !    reassembled Hamiltonian is constructed as H = C E C*, on a non-degenerate
  !    subspace.
  subroutine check_wfns()
    integer :: ik, igstart, igend, ngk_ik
    integer :: gvecs1(3,maxval(ngk)), gvecs2(3,maxval(ngk))
    integer :: isort1(maxval(ngk)), isort2(maxval(ngk))
    integer, allocatable :: isorti1(:)
    SCALAR, allocatable :: wfns1_ks(:,:), wfns2_ks(:,:), wfns1_ks_tmp(:,:), wfns2_ks_tmp(:,:)
    SCALAR, allocatable :: ham1(:,:), ham2(:,:), ham_diff(:,:)
    real(DP) :: err_ham_tmp
    real(DP), allocatable :: err1(:,:), err2(:,:)
    type(gspace) :: gvec
    real(DP) :: e1(nb), e2(nb)
    integer :: nb_, ib_found, ig, ib, ispin

    err_ortho1 = 0d0
    err_ortho2 = 0d0
    err_cross = 0d0
    err_ham = 0d0

    call hdf5_read_int(file_id_ref, 'mf_header/gspace/ng', gvec%ng)
    call hdf5_read_double(file_id_ref, 'mf_header/gspace/ecutrho', gvec%ecutrho)
    call hdf5_read_int_array(file_id_ref, 'mf_header/gspace/FFTgrid', [3], gvec%FFTgrid)
    SAFE_ALLOCATE(gvec%components, (3,gvec%ng))
    call hdf5_read_int_array(file_id_ref, 'mf_header/gspace/components', [3,gvec%ng], gvec%components)
    call gvec_index(gvec)
    SAFE_ALLOCATE(isorti1, (gvec%ng))

    igstart = 1
    do ik = 1, nk
      ngk_ik = ngk(ik)
      igend = igstart + ngk_ik - 1
      ! For each k-point, find mapping array between G-vectors in original file,
      ! f1, and the new file, f2.

      call hdf5_read_int_hyperslab(file_id_ref, 'wfns/gvecs', [3,ngk_ik], [0,igstart-1], gvecs1)
      call hdf5_read_int_hyperslab(file_id_new, 'wfns/gvecs', [3,ngk_ik], [0,igstart-1], gvecs2)
      isorti1(:) = 0
      do ig = 1, ngk_ik
        call findvector(isort1(ig), gvecs1(:,ig), gvec)
        isorti1(isort1(ig)) = ig
        call findvector(isort2(ig), gvecs2(:,ig), gvec)
      enddo
      do ig = 1, ngk_ik
        if (isorti1(isort2(ig))<1) call die('Error mapping G-vectors')
      enddo

      do ispin = 1, ns
        wfns1_ks = wfns1(igstart:igend,ispin,1:nb)
        wfns1_ks(:,:) = wfns1_ks(isorti1(isort2(1:ngk_ik)),:)
        wfns2_ks = wfns2(igstart:igend,ispin,1:nb)

        ! Check 2-3: WFNs orthonormal to within thenselves
        err1 = abs(matmul(transpose(wfns1_ks), MYCONJG(wfns1_ks)))
        forall(ib=1:nb) err1(ib,ib) = err1(ib,ib) - 1
        err_ortho1 = max(err_ortho1, maxval(err1))

        err2 = abs(matmul(transpose(wfns2_ks), MYCONJG(wfns2_ks)))
        forall(ib=1:nb) err2(ib,ib) = err2(ib,ib) - 1
        err_ortho2 = max(err_ortho2, maxval(err2))

        ! Check 4: cross WFN orthonormality
        !err12 = matmul(transpose(wfns2_ks),MYCONJG(wfns2_ks)) - matmul(transpose(wfns2_ks),MYCONJG(wfns2_ks))

        ! Check 5: reassembled Hamiltonian
        ! First, find the last non-degenerate band. This test is very picky
        ! with degenerancies!
        e1 = en1(:,ik,ispin)
        e2 = en2(:,ik,ispin)
        if (no_deg) then
          nb_ = nb
        else
          ib_found = -1
          do ib = nb-1, 1, -1
            if ((abs(e1(ib)-e1(ib+1)) > tol_deg) .and. &
                (abs(e2(ib)-e2(ib+1)) > tol_deg)) then
                ib_found = ib
                exit
            endif
          enddo
          if (ib_found==-1) call die('All bands are degenerate')
          nb_ = ib_found
        endif

        allocate(wfns1_ks_tmp, mold=wfns1_ks)
        allocate(wfns2_ks_tmp, mold=wfns2_ks)
        forall (ib=1:nb_) wfns1_ks_tmp(:,ib) = wfns1_ks(:,ib) * e1(ib)
        forall (ib=1:nb_) wfns2_ks_tmp(:,ib) = wfns2_ks(:,ib) * e2(ib)

        ham1 = matmul(MYCONJG(wfns1_ks(:,:nb_)), transpose(wfns1_ks_tmp(:,:nb_)))
        ham2 = matmul(MYCONJG(wfns2_ks(:,:nb_)), transpose(wfns2_ks_tmp(:,:nb_)))
        ham_diff = ham1 - ham2

        err_ham_tmp = get_norm(ham_diff) / (get_norm(ham1) * get_norm(ham2))
        err_ham = max(err_ham, err_ham_tmp)
        deallocate(wfns1_ks_tmp,wfns2_ks_tmp)

      enddo !ispin

      igstart = igend + 1
    enddo !ik


  end subroutine check_wfns

#endif
end program compare_wfns
