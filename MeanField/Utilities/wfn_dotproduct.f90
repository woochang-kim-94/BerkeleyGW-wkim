!=============================================================
! Utility:
! 
! wfn_dotproduct (flavorless). Originally by DAS.
!
! Form the overlap between the bands in two wavefunction files.
! If they are copies of the same file, you can check orthonormality.
! Only bands at corresponding k-points are considered, since
! the overlap is zero by Bloch`s theorem if the k-points differ.
!
! Warning: a known issue is possible errors from gmap when relating
! k-points by symmetry operations.
!
! 11/06/2014 - FHJ: Added HDF5 support. Note that, while probably faster the
! the native Fortran binary format, the HDF5 support is far from optimal
! because we perform the loops in the same order as the binary WFN files
! (band is the fast index), and not in the native order of the HDF5 files
! (kpt it the fast index).
!
!=============================================================

#include "f_defs.h"

program wfn_dotproduct

  use global_m
  use blas_m
  use find_kpt_match_m
  use gmap_m
  use input_utils_m
  use misc_m
  use sort_m
  use wfn_rho_vxc_io_m
  use wfn_io_hdf5_m
#ifdef HDF5
  use hdf5
#endif
  implicit none

  type(mf_header_t) :: mf1, mf2
  type(gspace) :: gvec_kpt
  character*11 :: inform
  character(len=16) :: sdq
  character*256 :: sformat, file_1, file_2, usage
  character :: sform
  integer :: ik1, ib1, ik2, ik3, ib2, ngkmax, ig, ispin, ispinor
  integer :: nargs, itqq, kgqq(3)
  real(DP), allocatable :: dwfn1(:,:), dwfn2(:,:,:), dph(:)
  complex(DPC), allocatable :: zwfn1(:,:), zwfn2(:,:,:), zph(:)
  integer, allocatable :: gindex1(:), gindex2(:), ind(:), isrt(:), ginv(:)
  logical :: informat
  real(DP) :: doverlap, q_shift(3)
  complex(DPC) :: zoverlap
#ifdef HDF5
  integer, allocatable :: all_gvecs1(:,:), all_gvecs2(:,:)
  real(DP), allocatable :: dwfn(:,:)
  complex(DPC), allocatable :: zwfn(:,:)
  integer :: ierr, ioffsetk, ngk, ngktot
#endif

  usage = 'Usage: wfn_dotproduct.x A|B WFN_1 WFN_2'
#ifdef HDF5
  usage = 'Usage: wfn_dotproduct.x A|B|H WFN_1 WFN_2'
#endif

  ! Get file names from command-line arguments

  nargs = command_argument_count()

  if (nargs < 3) then
    call die(usage)
  endif

  ! mpiexec from mpich1 may add 4 extra arguments to the list
  if (nargs > 3) then
    write(0,'(a,i3,a)') 'WARNING: ', nargs, ' arguments found, only first 3 being used.'
  endif

  call get_command_argument(1, sformat)
  call get_command_argument(2, file_1)
  call get_command_argument(3, file_2)

  ! Open units
  if (len(TRUNC(sformat))/=1) then
    call die(usage)
  endif
  sform = sformat(1:1)

  if (sform=='A') then
    informat = .true.
    inform = 'formatted'
  elseif (sform=='B') then
    informat = .false.
    inform = 'unformatted'
#ifdef HDF5
  elseif (sform=='H') then
    informat = .false.
    inform = 'HDF5'
#endif
  else  
    call die(usage)
  endif

  ! Read headers

  write(6,'(3x,a)') 'Reading ' // TRUNC(inform) // ' file ' // TRUNC(file_1)
  if (sform=='H') then
#ifdef HDF5
    call h5open_f(ierr)
    call read_hdf5_mf_header(TRUNC(file_1), mf1, sheader='WFN', &
      warn=.false., dont_warn_kgrid=.true.)
#endif
  else
    call open_file(unit=7, file=TRUNC(file_1), form=inform, status='old')
    call read_mf_header(7, mf1, sheader='WFN', warn=.false., dont_warn_kgrid=.true.)
  endif

  ! must have same flavor as first one
  write(6,'(3x,a)') 'Reading ' // TRUNC(inform) // ' file ' // TRUNC(file_2)
  if (sform=='H') then
#ifdef HDF5
    call read_hdf5_mf_header(TRUNC(file_2), mf2, sheader=mf1%sheader, &
      iflavor=mf1%iflavor, warn=.false., dont_warn_kgrid=.true.)
#endif
  else
    call open_file(unit=8, file=TRUNC(file_2), form=inform, status='old')
    call read_mf_header(8, mf2, iflavor=mf1%iflavor, sheader=mf1%sheader, warn=.false., dont_warn_kgrid=.true.)
  endif

  call check_header(TRUNC(file_1), mf1%kp, mf1%gvec, mf1%syms, mf1%crys, &
    TRUNC(file_2), mf2%kp, mf2%gvec, mf2%syms, mf2%crys, is_wfn = .true., tolerant = .true.)

  write(6,'(2(a,a,i8,a))') TRUNC(file_1), ' has ', mf1%kp%nrk, ' k-points; ', &
    TRUNC(file_2), ' has ', mf2%kp%nrk, ' k-points'
  if(mf2%kp%nrk > mf1%kp%nrk) then
    write(6,'(a)') "Note: only the second file's k-points will have symmetries applied. Switch the order on the command-line"
    write(6,'(a)') "if you want matches where a symmetry is applied to the first file's k-points."
  endif

  if(all(mf1%kp%kgrid(:) > 0 .and. mf2%kp%kgrid(:) > 0)) then
    q_shift(1:3) = dble(mf2%kp%shift(1:3))/mf2%kp%kgrid(1:3) - dble(mf1%kp%shift(1:3))/mf1%kp%kgrid(1:3)
  else
    q_shift(1:3) = 0d0
  endif
  write(6,'(a,3f12.6)') 'Using q-shift = ', q_shift(1:3)

  ! Read charge density gvectors
  SAFE_ALLOCATE(mf1%gvec%components, (3, mf1%gvec%ng))
  SAFE_ALLOCATE(mf2%gvec%components, (3, mf2%gvec%ng))

  ! use mf1%gvec as the master list of G-vectors
  if (sform=='H') then
#ifdef HDF5
    call read_hdf5_gvectors(TRUNC(file_1), mf1%gvec%ng, mf1%gvec%components)
#endif
  else
    call read_gvectors(7, informat, mf1%gvec%ng, mf1%gvec%ng, mf1%gvec%components)
  endif
  ! needed to be able to use findvector
  call gvec_index(mf1%gvec)

  ngkmax = max(mf1%kp%ngkmax, mf2%kp%ngkmax)
  if (mf1%iflavor .eq. 1) then
#ifdef HDF5
    if (sform=='H') then
      SAFE_ALLOCATE(dwfn, (mf1%gvec%ng, mf1%kp%nspin))
    endif
#endif
    SAFE_ALLOCATE(dwfn1, (mf1%gvec%ng, mf1%kp%nspin))
    SAFE_ALLOCATE(dwfn2, (mf1%gvec%ng, mf2%kp%nspin, mf2%kp%mnband))
    SAFE_ALLOCATE(dph, (mf1%gvec%ng))
  else
#ifdef HDF5
    if (sform=='H') then
      SAFE_ALLOCATE(zwfn, (mf2%gvec%ng, mf1%kp%nspin*mf1%kp%nspinor))
    endif
#endif
    SAFE_ALLOCATE(zwfn1, (mf2%gvec%ng, mf1%kp%nspin*mf1%kp%nspinor))
    SAFE_ALLOCATE(zwfn2, (mf2%gvec%ng, mf2%kp%nspin*mf2%kp%nspinor, mf2%kp%mnband))
    SAFE_ALLOCATE(zph, (mf1%gvec%ng))
  end if
  SAFE_ALLOCATE(gindex1, (ngkmax))
  SAFE_ALLOCATE(gindex2, (ngkmax))
  SAFE_ALLOCATE(ind, (ngkmax))
  SAFE_ALLOCATE(gvec_kpt%components, (3, ngkmax))
  SAFE_ALLOCATE(ginv, (mf1%gvec%ng))
  SAFE_ALLOCATE(mf1%gvec%ekin, (mf1%gvec%ng))
  SAFE_ALLOCATE(isrt, (mf1%gvec%ng))

  if (sform=='H') then
#ifdef HDF5
    ngktot = SUM(mf1%kp%ngk)
    SAFE_ALLOCATE(all_gvecs1, (3, ngktot))
    call read_hdf5_wfn_gvectors(TRUNC(file_1), all_gvecs1, ngktot)
    ngktot = SUM(mf2%kp%ngk)
    SAFE_ALLOCATE(all_gvecs2, (3, ngktot))
    call read_hdf5_wfn_gvectors(TRUNC(file_2), all_gvecs2, ngktot)
#endif
  endif

  do ik1 = 1, mf1%kp%nrk
    write(6,'(/,a,i6,a,3f12.6)') 'k-point ', ik1, ': ', mf1%kp%rk(1:3, ik1)
    call find_kpt_match(mf2%kp, mf2%syms, mf1%kp%rk(1:3, ik1) + q_shift(1:3), ik2, itqq, kgqq)
    if(ik2 == 0) then
      write(6,'(a)') 'No match in ' // TRUNC(file_2)
      if (sform/='H') then
        call read_gvectors(7, informat, mf1%kp%ngk(ik1), ngkmax, gvec_kpt%components, dont_read=.not.informat)
      endif
    else
      write(6,'(a,i6,3a,3f12.6)') 'Matches k-point ', ik2, ' in ', TRUNC(file_2), ':', mf2%kp%rk(1:3, ik2)
      write(6,'(a,i3,a,3i3)') 'with symmetry operation ', itqq, ' and Umklapp ', kgqq
      if (sform=='H') then
#ifdef HDF5
        ngk = mf1%kp%ngk(ik1)
        ioffsetk = sum(mf1%kp%ngk(1:ik1-1))
        gvec_kpt%components(1:3, 1:ngk) = all_gvecs1(1:3, ioffsetk+1:ioffsetk+ngk)
#endif
      else
        call read_gvectors(7, informat, mf1%kp%ngk(ik1), ngkmax, gvec_kpt%components)
      endif
      write(6,'(2(a10,a16,4x),a10,6x,a20)',advance='no') 'band 1', 'energy 1 ', 'band 2', 'energy 2 ', 'spin ', 'Re overlap '
      if(mf1%iflavor == 2) write(6,'(a20)',advance='no') 'Im overlap '
      write(6,'(a20)') '|overlap|^2 '

      do ig = 1, mf1%kp%ngk(ik1)
        call findvector(gindex1(ig), gvec_kpt%components(:, ig), mf1%gvec)
        if (gindex1(ig) ==  0)  call die('could not match G-vector')
      enddo
    endif

    if(ik2 > 0) then
      if(mf1%kp%ngk(ik1) /= mf2%kp%ngk(ik2)) then
        call die("Internal error: ngk mismatch")
      endif

      if (sform=='H') then
#ifdef HDF5
        ngk = mf2%kp%ngk(ik2)
        ioffsetk = sum(mf2%kp%ngk(1:ik2-1))
        gvec_kpt%components(1:3, 1:ngk) = all_gvecs2(1:3, ioffsetk+1:ioffsetk+ngk)
#endif
      else
        rewind(8)
        ! FHJ: Note: read_header_type has intent(out) for gvec, so we need to reallocate
        ! gvec%components
        SAFE_DEALLOCATE_P(mf2%gvec%components)
        call dealloc_header_type(mf2%sheader, mf2%crys, mf2%kp)
        call read_header_type(8, informat, mf2%sheader, mf1%iflavor, mf2%kp, mf2%gvec, mf2%syms, mf2%crys, &
          warn=.false., dont_warn_kgrid=.true.)
        SAFE_ALLOCATE(mf2%gvec%components, (3, mf2%gvec%ng))
        call check_header(TRUNC(file_1), mf1%kp, mf1%gvec, mf1%syms, mf1%crys, &
          TRUNC(file_2), mf2%kp, mf2%gvec, mf2%syms, mf2%crys, is_wfn = .true., tolerant = .true.)
        call read_gvectors(8, informat, mf2%gvec%ng, mf2%gvec%ng, mf2%gvec%components, dont_read = .not. informat)

        do ik3 = 1, ik2 - 1
          call read_gvectors(8, informat, mf2%kp%ngk(ik3), ngkmax, gvec_kpt%components, dont_read = .not. informat)
          do ib2 = 1, mf2%kp%mnband
            if(mf1%iflavor == 1) then
              call read_real_data(8, informat, mf2%kp%ngk(ik3), mf2%gvec%ng, mf2%kp%nspin, &
                dwfn2(:,:,ib2), dont_read =.not.informat)
            else
              call read_complex_data(8, informat, mf2%kp%ngk(ik3), mf2%gvec%ng, mf2%kp%nspin*mf2%kp%nspinor, &
                zwfn2(:,:,ib2), dont_read =.not.informat)
            endif
          enddo
        enddo
        call read_gvectors(8, informat, mf2%kp%ngk(ik2), ngkmax, gvec_kpt%components)
        ! now we are at the right point in the file
      endif

      ginv = 0
      do ig = 1, mf2%kp%ngk(ik2)
        call findvector(gindex2(ig), gvec_kpt%components(:, ig), mf1%gvec)
        if (gindex2(ig) ==  0)  call die('could not match G-vector')
        ginv(gindex2(ig)) = ig
      enddo

      call kinetic_energies(mf1%gvec, mf1%crys%bdot, mf1%gvec%ekin, qvec=mf2%kp%rk(:, ik2))
      call sortrx(mf1%gvec%ng, mf1%gvec%ekin, isrt, gvec=mf1%gvec%components)

      ! FHJ: Keep the following line. It prevents a compiler bug with sunf90
      ! which changes the value of mf1%iflavor.
      mf1%iflavor = mf1%iflavor
      if(mf1%iflavor == 1) then
        call gmap(mf1%gvec, mf1%syms, mf2%kp%ngk(ik2), itqq, kgqq, isrt, ginv, ind, dph, .true.)
      else
        call gmap(mf1%gvec, mf1%syms, mf2%kp%ngk(ik2), itqq, kgqq, isrt, ginv, ind, zph, .true.)
      endif

      if(any(ind(1:mf2%kp%ngk(ik2)) == 0)) call die("ind array from gmap has a zero")

      ! Read "second" WFN
      if (sform=='H') then
#ifdef HDF5
        ngk = mf2%kp%ngk(ik2)
        ioffsetk = sum(mf2%kp%ngk(1:ik2-1))
        do ib2 = 1, mf2%kp%mnband
          if(mf1%iflavor == 1) then
            call read_hdf5_band(TRUNC(file_2), dwfn(1:ngk,:), ngk, &
              mf2%kp%nspin*mf2%kp%nspinor, ioffsetk, ib2-1)
            dwfn2(:,:,ib2) = 0d0
            dwfn2(gindex2(1:ngk),:,ib2) = dwfn(1:ngk,:)
          else
            call read_hdf5_band(TRUNC(file_2), zwfn(1:ngk,:), ngk, &
              mf2%kp%nspin*mf2%kp%nspinor, ioffsetk, ib2-1)
            zwfn2(:,:,ib2) = ZERO
            zwfn2(gindex2(1:ngk),:,ib2) = zwfn(1:ngk,:)
          endif

        enddo
#endif
      else
        do ib2 = 1, mf2%kp%mnband
          if(mf1%iflavor == 1) then
            call read_real_data(8, informat, mf2%kp%ngk(ik2), mf2%gvec%ng, mf2%kp%nspin, &
              dwfn2(:,:,ib2), gindex=gindex2)
          else
            call read_complex_data(8, informat, mf2%kp%ngk(ik2), mf2%gvec%ng, mf2%kp%nspin*mf2%kp%nspinor, &
              zwfn2(:,:,ib2), gindex=gindex2)
          endif
        enddo
      endif

      if(mf1%iflavor == 1) then
        do ig = 1, mf2%kp%ngk(ik2)
          dwfn2(ind(ig), :, :) = dwfn2(ind(ig), :, :) * dph(ig)
        enddo
      else
        do ig = 1, mf2%kp%ngk(ik2)
          zwfn2(ind(ig), :, :) = zwfn2(ind(ig), :, :) * zph(ig)
        enddo
      endif
    endif !ik2/=0

    ! Read "first" WFN
    do ib1 = 1, mf1%kp%mnband
      if (sform=='H') then
#ifdef HDF5
        if(ik2 == 0) cycle
        ngk = mf1%kp%ngk(ik1)
        ioffsetk = sum(mf1%kp%ngk(1:ik1-1))
        if(mf1%iflavor == 1) then
          call read_hdf5_band(TRUNC(file_1), dwfn(1:ngk,:), ngk, &
            mf1%kp%nspin*mf1%kp%nspinor, ioffsetk, ib1-1)
          dwfn1(:,:) = 0d0
          dwfn1(gindex1(1:ngk),:) = dwfn(1:ngk,:)
        else
          call read_hdf5_band(TRUNC(file_1), zwfn(1:ngk,:), ngk, &
            mf1%kp%nspin*mf1%kp%nspinor, ioffsetk, ib1-1)
          zwfn1(:,:) = ZERO
          zwfn1(gindex1(1:ngk),:) = zwfn(1:ngk,:)
        endif
#endif
      else
        if(mf1%iflavor == 1) then
          call read_real_data(7, informat, mf1%kp%ngk(ik1), mf1%gvec%ng, mf1%kp%nspin, dwfn1, &
            dont_read = (ik2 == 0) .and. .not. informat, gindex = gindex1)
        else
          call read_complex_data(7, informat, mf1%kp%ngk(ik1), mf1%gvec%ng, mf1%kp%nspin*mf1%kp%nspinor, zwfn1, &
            dont_read = (ik2 == 0) .and. .not. informat, gindex = gindex1)
        endif
        if(ik2 == 0) cycle
      endif

      do ib2 = 1, mf2%kp%mnband
        if(mf1%iflavor == 1) then
          do ispin = 1, mf1%kp%nspin
            doverlap = ddot(mf1%gvec%ng, dwfn1(:, ispin), 1, dwfn2(:, ispin, ib2), 1)
            write(6,'(2(i8,2x,g20.8),i8,8x,2g20.10)') ib1, mf1%kp%el(ib1, ik1, ispin), ib2, mf2%kp%el(ib2, ik2, ispin), &
              ispin, doverlap, dble(doverlap)**2
          enddo
        else
          do ispin = 1, mf1%kp%nspin
            if(mf1%kp%nspinor.eq.2) then
              zoverlap = 0.0D0
              do ispinor = 1, mf1%kp%nspinor 
                zoverlap = zoverlap + zdotc(mf1%gvec%ng, zwfn1(:, ispin*ispinor), 1, zwfn2(:, ispin*ispinor, ib2), 1)
              enddo
            else
              zoverlap = zdotc(mf1%gvec%ng, zwfn1(:, ispin), 1, zwfn2(:, ispin, ib2), 1)
            endif
            write(6,'(2(i8,2x,g20.8),i8,8x,3g20.10)') ib1, mf1%kp%el(ib1, ik1, ispin), ib2, mf2%kp%el(ib2, ik2, ispin), &
              ispin, zoverlap, dble(zoverlap)**2 + aimag(zoverlap)**2
          enddo
        endif
      enddo
    enddo
  enddo

  ! Close files

  if (sform=='H') then
#ifdef HDF5
  call h5close_f(ierr)
#endif
  else
    call close_file(7)
    call close_file(8)
  endif

  ! Deallocate arrays

  if (mf1%iflavor .eq. 1) then
#ifdef HDF5
    if (sform=='H') then
      SAFE_DEALLOCATE(dwfn)
    endif
#endif
    SAFE_DEALLOCATE(dwfn1)
    SAFE_DEALLOCATE(dwfn2)
    SAFE_DEALLOCATE(dph)
  else
#ifdef HDF5
    if (sform=='H') then
      SAFE_DEALLOCATE(zwfn)
    endif
#endif
    SAFE_DEALLOCATE(zwfn1)
    SAFE_DEALLOCATE(zwfn2)
    SAFE_DEALLOCATE(zph)
  endif

  SAFE_DEALLOCATE_P(mf1%gvec%components)
  SAFE_DEALLOCATE_P(mf2%gvec%components)
  SAFE_DEALLOCATE_P(mf1%gvec%index_vec)
  SAFE_DEALLOCATE_P(gvec_kpt%components)
  SAFE_DEALLOCATE_P(mf1%gvec%ekin)
  SAFE_DEALLOCATE(isrt)
  SAFE_DEALLOCATE(ind)
  SAFE_DEALLOCATE(ginv)
  SAFE_DEALLOCATE(gindex1)
  SAFE_DEALLOCATE(gindex2)
  call dealloc_header_type(mf1%sheader, mf1%crys, mf1%kp)
  call dealloc_header_type(mf2%sheader, mf2%crys, mf2%kp)

end program wfn_dotproduct
