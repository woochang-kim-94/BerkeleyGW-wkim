!=============================================================
! Utility:
! 
! wfn_time_reversal (flavorless). Originally by DAS.
!
!  Unfold the k-points in a WFN file according to time-reversal.
!  This allows use of the 'automatic' k-points from ESPRESSO
!  for an unshifted or half-shifted Monkhorst-Pack grid.
!  u_{-k}(G) = u_{k}(-G)* (spinless)
!  For spin-polarized or spinors:
!  u_{-k}(G,up) = -u_{k}(-G,dn)*
!  u_{-k}(G,dn) =  u_{k}(-G,up)*
!
!  You can use this utility for real WFNs but it is probably
!  useless: time-reversal has the same effect on k-points as
!  inversion, which needs to present for real WFNs to be allowed.
!
!=============================================================
! In this section we study another discrete symmetry operator,
! called 'time reversal'. This is a difficult topic for the novice,
! partly because the term 'time reversal' is a misnomer; it reminds
! us of science fiction. What we do in this section can be more
! appropriately characterized by the term 'reversal of motion.'
!           - J. J. Sakurai, Modern Quantum Mechanics
!=============================================================

#include "f_defs.h"

program wfn_time_reversal

  use global_m
  use find_kpt_match_m
  use misc_m
  use wfn_rho_vxc_io_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp, kp2
  type(gspace) :: gvec, gvec_kpt
  character*3 :: sheader
  character*11 :: inform
  character*256 :: sformat, file_in, file_out, usage
  integer :: ik, ib, iflavor, ig, nargs, sym_umklapp(3), isym, ik_match
  real(DP) :: new_kpt(3)
  real(DP), allocatable :: dwfn(:,:), dwfn_spin(:)
  complex(DPC), allocatable :: zwfn(:,:), zwfn_spin(:)
  logical :: informat
  logical, allocatable :: time_reversible(:)
  integer, allocatable :: umklapp(:,:), new_gvecs(:,:)

  usage = 'Usage: wfn_time_reversal.x A|B WFN_in WFN_out'

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
  call get_command_argument(2, file_in)
  call get_command_argument(3, file_out)

  ! Open units

  if (TRUNC(sformat) .eq. 'A') then
    informat = .true.
    inform = 'formatted'
  elseif (TRUNC(sformat) .eq. 'B') then
    informat = .false.
    inform = 'unformatted'
  else  
    call die(usage)
  endif

  ! Read header

  sheader='WFN'
  iflavor=-1
  write(6,'(a)') 'Reading ' // TRUNC(inform) // ' file ' // TRUNC(file_in)
  call open_file(unit=7,file=TRUNC(file_in),form=inform,status='old')
  call read_header_type(7, informat, sheader, iflavor, kp, gvec, syms, crys, warn = .false., dont_warn_kgrid = .true.)

  kp2%nspinor  = kp%nspinor
  kp2%nspin    = kp%nspin
  kp2%mnband   = kp%mnband
  kp2%kgrid(:) = kp%kgrid
  kp2%shift(:) = kp%shift(:)
  kp2%ecutwfc  = kp%ecutwfc
  kp2%ngkmax   = kp%ngkmax

  SAFE_ALLOCATE(kp2%ngk, (2*kp%nrk))
  kp2%ngk(1:kp%nrk) = kp%ngk(1:kp%nrk)
  SAFE_ALLOCATE(kp2%rk, (3, 2*kp%nrk))
  kp2%rk(1:3, 1:kp%nrk) = kp%rk(1:3, 1:kp%nrk)
  SAFE_ALLOCATE(kp2%ifmin, (2*kp%nrk, kp%nspin))
  kp2%ifmin(1:kp%nrk, 1:kp%nspin) = kp%ifmin(1:kp%nrk, 1:kp%nspin)
  SAFE_ALLOCATE(kp2%ifmax, (2*kp%nrk, kp%nspin))
  kp2%ifmax(1:kp%nrk, 1:kp%nspin) = kp%ifmax(1:kp%nrk, 1:kp%nspin)
  SAFE_ALLOCATE(kp2%w, (2*kp%nrk))
  kp2%w(1:kp%nrk) = kp%w(1:kp%nrk)
  SAFE_ALLOCATE(kp2%el, (kp%mnband, 2*kp%nrk, kp%nspin))
  kp2%el(1:kp%mnband, 1:kp%nrk, 1:kp%nspin) = kp%el(1:kp%mnband, 1:kp%nrk, 1:kp%nspin)
  SAFE_ALLOCATE(kp2%occ, (kp%mnband, 2*kp%nrk, kp%nspin))
  kp2%occ(1:kp%mnband, 1:kp%nrk, 1:kp%nspin) = kp%occ(1:kp%mnband, 1:kp%nrk, 1:kp%nspin)

  kp2%nrk = kp%nrk
  SAFE_ALLOCATE(time_reversible, (1:kp%nrk))
  SAFE_ALLOCATE(umklapp, (1:3, 1:kp%nrk))

  do ik = 1, kp%nrk
    write(6,'(a,i6,a,3f12.6)') 'Old k-point ', ik, ': ', kp%rk(1:3, ik)

    new_kpt(:) = -kp%rk(:, ik)
    call k_range(new_kpt(:), umklapp(:, ik), TOL_Zero)
    write(6,'(a,3f12.6)') '  Time-reverses to : ', new_kpt(1:3)

    ! Make sure this is not related by symmetry/Umklapp to a k-point already in the set
    ! k-coordinates 0 or 0.5 are mapped to themselves by time-reversal (and possible Umklapp)
    call find_kpt_match(kp, syms, new_kpt(:), ik_match, isym, sym_umklapp)
    time_reversible(ik) = (ik_match == 0)

    if(time_reversible(ik)) then
      kp2%nrk = kp2%nrk + 1
      kp2%rk(:, kp2%nrk) = new_kpt(:)

      ! divide weight between new and old k-points
      kp2%w(ik)      = kp%w(ik) * 0.5d0
      kp2%w(kp2%nrk) = kp%w(ik) * 0.5d0

      kp2%ngk(kp2%nrk) = kp%ngk(ik)
      kp2%ifmin(kp2%nrk, 1:kp%nspin) = kp%ifmin(ik, 1:kp%nspin)
      kp2%ifmax(kp2%nrk, 1:kp%nspin) = kp%ifmax(ik, 1:kp%nspin)
      kp2%el (1:kp%mnband, kp2%nrk, 1:kp%nspin) = kp%el (1:kp%mnband, ik, 1:kp%nspin)
      kp2%occ(1:kp%mnband, kp2%nrk, 1:kp%nspin) = kp%occ(1:kp%mnband, ik, 1:kp%nspin)
    else
      write(6,'(a,i6,a,i2,a,3i3)') '  But matches old k-point ', ik_match, &
        ' with symmetry number ', isym, ' and Umklapp ', sym_umklapp(:)
      kp2%w(ik) = kp%w(ik)
    endif
  enddo

  if(.not. any(time_reversible(1:kp%nrk))) then
    write(6,'(a)') 'No new unique k-points are generated by time-reversal. Not writing a new WFN file.'
    stop
  endif

  write(6,'(a)') 'Writing ' // TRUNC(inform) // ' file ' // TRUNC(file_out)
  call open_file(unit=8,file=TRUNC(file_out),form=inform,status='replace')
  call write_header_type(8, informat, sheader, iflavor, kp2, gvec, syms, crys, warn = .false.)

  ! charge density gvectors
  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  call read_gvectors (7, informat, gvec%ng, gvec%ng, gvec%components)
  call write_gvectors(8, informat, gvec%ng, gvec%ng, gvec%components)

  if (iflavor .eq. 1) then
    SAFE_ALLOCATE(dwfn, (gvec%ng, kp%nspin))
    if(kp%nspin == 2) then
      SAFE_ALLOCATE(dwfn_spin, (gvec%ng))
    endif
  else
    SAFE_ALLOCATE(zwfn, (gvec%ng, kp%nspin))
    if(kp%nspin*kp%nspinor == 2) then
      SAFE_ALLOCATE(zwfn_spin, (gvec%ng))
    endif
  end if

  SAFE_ALLOCATE(gvec_kpt%components, (3, kp%ngkmax))
  SAFE_ALLOCATE(new_gvecs, (3, kp%ngkmax))

  ! Re-write original k-points
  do ik = 1, kp%nrk
    call read_gvectors (7, informat, kp%ngk(ik), kp%ngkmax, gvec_kpt%components)
    call write_gvectors(8, informat, kp%ngk(ik), kp%ngkmax, gvec_kpt%components)

    do ib = 1, kp%mnband
      if(iflavor == 1) then
        call read_real_data    (7, informat, kp%ngk(ik), gvec%ng, kp%nspin, dwfn(:,:))
        call write_real_data   (8, informat, kp%ngk(ik), gvec%ng, kp%nspin, dwfn(:,:))
      else
        call read_complex_data (7, informat, kp%ngk(ik), gvec%ng, kp%nspin*kp%nspinor, zwfn(:,:))
        call write_complex_data(8, informat, kp%ngk(ik), gvec%ng, kp%nspin*kp%nspinor, zwfn(:,:))
      endif
    enddo
  enddo

  call close_file(7)
  call dealloc_header_type(sheader, crys, kp)
  call open_file(unit=7,file=TRUNC(file_in),form=inform,status='old')
  ! FHJ: Note: read_header_type has intent(out) for gvec, so we need to reallocate
  ! gvec%components
  SAFE_DEALLOCATE_P(gvec%components)
  call read_header_type(7, informat, sheader, iflavor, kp, gvec, syms, crys, warn = .false., dont_warn_kgrid = .true.)
  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  call read_gvectors(7, informat, gvec%ng, gvec%ng, gvec%components)

  ! Write new k-points
  do ik = 1, kp%nrk
    call read_gvectors (7, informat, kp%ngk(ik), kp%ngkmax,  gvec_kpt%components, dont_read = .not. time_reversible(ik))

    if(time_reversible(ik)) then
      do ig = 1, kp%ngk(ik)
        new_gvecs(1:3, ig) = -umklapp(1:3, ik) - gvec_kpt%components(1:3, ig)
      enddo
      call write_gvectors(8, informat, kp%ngk(ik), kp%ngkmax, new_gvecs)
    endif

    ! Even if we will not generate a new k-point, we must read through the data to get to the next k-point.
    do ib = 1, kp%mnband
      if(iflavor == 1) then
        call read_real_data(7, informat, kp%ngk(ik), gvec%ng, kp%nspin, dwfn(:,:), &
          dont_read = .not. time_reversible(ik) .and. .not. informat)
        if(time_reversible(ik)) then
          if(kp%nspin == 2) then
            dwfn_spin(:) = dwfn(:, 2)
            dwfn(:, 2) = dwfn(:, 1)
            dwfn(:, 1) = -dwfn_spin(:)
          endif
          call write_real_data   (8, informat, kp%ngk(ik), gvec%ng, kp%nspin, dwfn(:,:))
        endif
      else
        call read_complex_data(7, informat, kp%ngk(ik), gvec%ng, kp%nspin*kp%nspinor, zwfn(:,:), &
          dont_read = .not. time_reversible(ik) .and. .not. informat)
        if(time_reversible(ik)) then
          if(kp%nspin*kp%nspinor == 1) then
            zwfn(:,:) = conjg(zwfn(:,:))
          else
            zwfn_spin(:) = conjg(zwfn(:, 2))
            zwfn(:, 2) = conjg(zwfn(:, 1))
            zwfn(:, 1) = -conjg(zwfn_spin(:))
          endif
          call write_complex_data(8, informat, kp%ngk(ik), gvec%ng, kp%nspin*kp%nspinor, zwfn(:,:))
        endif
      endif
    enddo
  enddo

  ! Close files

  call close_file(7)
  call close_file(8)

  ! Deallocate

  if (iflavor .eq. 1) then
    SAFE_DEALLOCATE(dwfn)
    if(kp%nspin == 2) then
      SAFE_DEALLOCATE(dwfn_spin)
    endif
  else
    SAFE_DEALLOCATE(zwfn)
    if(kp%nspin*kp%nspinor == 2) then
      SAFE_DEALLOCATE(zwfn_spin)
    endif
  endif

  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE_P(gvec_kpt%components)
  SAFE_DEALLOCATE(time_reversible)
  SAFE_DEALLOCATE(umklapp)
  SAFE_DEALLOCATE(new_gvecs)
  call dealloc_header_type(sheader, crys, kp)
  call dealloc_kp(kp2)

end program wfn_time_reversal
