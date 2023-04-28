!===============================================================================
!
! Program:
!
! (1) forces         Originally by DAS       Last Edited: 2/3/2012 (DAS)
!
!     Calculates excited-state forces for excitons calculated by BSE diagonalization.
!
! A reformulation of the approach described in
! Sohrab Ismail-Beigi and Steven G. Louie, Excited-State Forces within a First-Principles
! Green`s Function Formalism, Phys. Rev. Lett. 90, 076401 (2003)
!
! Needed files: WFN_fi, eigenvectors, eqp.dat, dfpt.dat_mode*, displacements.dat
!     
!================================================================================

#include "f_defs.h"

program forces

  use global_m
  use eqpcor_m
  use wfn_rho_vxc_io_m
  use write_program_header_m
  implicit none

  character :: sheader*3, filename*20, fncor*32
  integer :: iflavor, iunit, ik, ierr, ndiag, noffdiag, spin_index(2), icp, ivp, ncband_max, nvband_max, ii, iat
  integer :: nkpt, ncband, nvband, nspin, ifmax, iexc, jexc, nexc, imode, nmode, mm, is, ic, iv, ib, ib2, ioffdiag
  integer, allocatable :: vband(:), cband(:), diag(:), offdiag1(:), offdiag2(:), &
    invoffdiagv(:,:), invoffdiagc(:,:), invdiagv(:), invdiagc(:)
  type(kpoints), target :: kp
  type(gspace) :: gvec
  type(symmetry) :: syms
  type(crystal) :: crys
  real(DP) :: kk(3), energy
  logical :: believe_sohrab = .true.
  real(DP), pointer :: el_denom(:,:,:) !< energies to use in denominators
  SCALAR, allocatable :: force_diag(:), force_offdiag(:), mtxel_sum(:,:,:), force_mode(:), force_atoms(:,:)
  SCALAR, allocatable :: mtxel(:,:,:) !< idiag/ioffdiag, ispin, ik
  SCALAR, allocatable :: asvck(:,:,:,:,:) !< iexc, ispin, iv, ic, ik
  SCALAR, allocatable :: aread(:) !< iexc
  complex(DPC), allocatable :: disp(:,:,:) !< idir, iat, imode 

!  call peinfo_init()

!---------------------------
! Write header

  call write_program_header('BSE/Forces', .false.)

!---------------------------
! Read input file

!  call open_file(8,file='forces.inp',form='formatted',status='old')
  ! what do we need to know? not much I thinks

! nvband_max, ncband_max, iexc, jexc (range perhaps)

!  call close_file(8)

!-----------------------------------------------------------------------
! Read info for crystal from WFN_fi

  if(peinf%inode == 0) call open_file(25,file='WFN_fi',form='unformatted',status='old')

  sheader = 'WFN'
  iflavor = 0
  call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys, dont_warn_kgrid = .true.)
  call close_file(25)

  ifmax = kp%ifmax(1,1)
  if(any(kp%ifmax(:,:) /= ifmax)) call die("not implemented for non-constant ifmax")

  ! only implemented for momentum operator so far.

!-----------------------------------------------------------------------
! Read mode displacements

  nmode = 3 * crys%nat
  call open_file(unit = 16, file = 'displacements.dat', form = 'formatted', status = 'old')
  SAFE_ALLOCATE(disp, (3, crys%nat, nmode))
  do imode = 1, nmode
    read(16,*) ((disp(ii, iat, imode), ii = 1, 3), iat = 1, crys%nat)
  enddo
  call close_file(16)

!-----------------------------------------------------------------------
! Read A_svck

  iunit = 77
  call open_file(unit=iunit,file='eigenvectors',form='unformatted',status='old')

  read(iunit) nspin
  read(iunit) nvband
  read(iunit) ncband
  read(iunit) nkpt

  write(6,*) 'nspin = ', nspin
  write(6,*) 'nvband = ', nvband
  write(6,*) 'ncband = ', ncband
  write(6,*) 'nkpt = ', nkpt
  write(6,*) 'ifmax = ', ifmax

  if(nkpt /= kp%nrk) write(0,'(a,i6,a,i5)') 'WARNING: k-point mismatch. eigenvectors has ', nkpt, ' but WFN_fi has ', kp%nrk
  if(nspin /= kp%nspin) write(0,'(a,i6,a,i5)') 'WARNING: nspin mismatch. eigenvectors has ', nspin, ' but WFN_fi has ', kp%nspin
  if(nvband > minval(kp%ifmax(:,:))) call die("Too many valence bands in eigenvectors.")
  if(ncband > kp%mnband - maxval(kp%ifmax(:,:))) call die("Too many conduction bands in eigenvectors.")

  read(iunit) ! k-points coords. We do not care.

  nexc = nspin * nvband * ncband * nkpt
  write(6,*) 'nexc = ', nexc
  SAFE_ALLOCATE(asvck, (nexc, nspin, nvband, ncband, nkpt))
  SAFE_ALLOCATE(aread, (nexc))

  do iexc = 1, nexc

    read(iunit) energy
!    write(6,*) 'energy = ', energy
    read(iunit) aread(:)
!    if (nmown(ii) .eq. peinf%inode)  then
    mm = 0
    do ik = 1, nkpt
      do ic = 1, ncband
        do iv = 1, nvband
          do is = 1, nspin
            mm = mm + 1
!              A(myown(ii),is,iv,ic,ik) = Aread(m)
            asvck(iexc, is, iv, ic, ik) = aread(mm)
          enddo
        enddo
      enddo
    enddo
!    endif

  enddo

  SAFE_DEALLOCATE(aread)
  call close_file(iunit)

  do iexc = 1, nexc
    if(abs(sum(abs(asvck(iexc,:,:,:,:))**2) - 1d0) > TOL_Small) call die("eigenvectors are not normalized.")
  enddo

  write(*,*) 'asvck 1: ', asvck(1,1,1,1,1), asvck(1,1,1,2,1)
  write(*,*) 'asvck 2: ', asvck(2,1,1,1,1), asvck(2,1,1,2,1)
  
!-----------------------------------------------------------------------
! Read quasiparticle energies from eqp.dat

  SAFE_ALLOCATE(kp%elda, (kp%mnband, kp%nrk, kp%nspin))
  kp%elda(:,:,:) = kp%el(:,:,:)

! (what about scissors?) intwfn may write eqp.dat correctly if interpolation is done.
!  call scissors_shift(kp, eqp%scis, eqp%spl_tck)

  fncor='eqp.dat'
  ! FIXME: for metals this is asking for a few more bands than actually needed on some k-points
  call eqpcor(fncor,peinf%inode,peinf%npes,kp,minval(kp%ifmax(:,:)-nvband+1),&
    maxval(kp%ifmax(:,:)+ncband),0,0,kp%el,kp%el,kp%el,1,0)

  if(believe_sohrab) then
    el_denom => kp%el
  else
    el_denom => kp%elda
  endif

  write(*,*) 'elda = ', kp%elda(1:nvband+ncband,:,:)*ryd
  write(*,*) 'eqp  = ', kp%el(1:nvband+ncband,:,:)*ryd
  write(*,*) 'denom = ', el_denom(1:nvband+ncband,:,:)*ryd

!-----------------------------------------------------------------------
! Read electron-phonon matrix elements

  if(peinf%inode == 0) write(6,'(a)') 'Reading dfpt.dat files'

  SAFE_ALLOCATE(vband, (nvband))
  SAFE_ALLOCATE(cband, (ncband))

  do iv = 1, nvband
    vband(iv) = ifmax + 1 - iv
  enddo
  do ic = 1, ncband
    cband(ic) = ifmax + ic
  enddo

!  write(*,*) 'vband = ', vband(:)
!  write(*,*) 'cband = ', cband(:)

  ndiag = nvband + ncband
  noffdiag = ndiag**2
  SAFE_ALLOCATE(diag, (ndiag))
  SAFE_ALLOCATE(offdiag1, (noffdiag))
  SAFE_ALLOCATE(offdiag2, (noffdiag))
  SAFE_ALLOCATE(invdiagv, (nvband))               ! map iv to mtxel indexing
  SAFE_ALLOCATE(invdiagc, (ncband))               ! map ic to mtxel indexing
  SAFE_ALLOCATE(invoffdiagv, (nvband, nvband))    ! map iv,ivp to mtxel indexing
  SAFE_ALLOCATE(invoffdiagc, (ncband, ncband))    ! map ic,icp to mtxel indexing

  do ib = 1, nvband + ncband
    diag(ib) = ifmax - nvband + ib
    if(ib <= nvband) invdiagv(ib)          = ifmax - ib + 1
    if(ib >  nvband) invdiagc(ib - nvband) = ib
  enddo

!  write(*,*) 'invdiagv = ', invdiagv(:)
!  write(*,*) 'invdiagc = ', invdiagc(:)
  
  invoffdiagv(:,:) = ZERO
  invoffdiagc(:,:) = ZERO

  ioffdiag = 1
  do ib2 = 1, ndiag
    do ib = 1, ndiag
      offdiag1(ioffdiag) = diag(ib)
      offdiag2(ioffdiag) = diag(ib2)
      if(ib <= nvband .and. ib2 <= nvband) then
        invoffdiagv(nvband - ib + 1, nvband - ib2 + 1) = ioffdiag + ndiag  ! this is probably wrong!!!!!
!        write(6,*) 'invoffdiagv(', ib, ib2, ') = ', ioffdiag + ndiag
      endif
      if(ib >  nvband .and. ib2 >  nvband) then
        invoffdiagc(ib - nvband, ib2 - nvband) = ioffdiag + ndiag
!        write(6,*) 'invoffdiagc(', ib, ib2, ') = ', ioffdiag + ndiag
      endif
      ioffdiag = ioffdiag + 1
    enddo
  enddo

!  write(*,*) 'invoffdiagv = ', invoffdiagv(:,:)
!  write(*,*) 'invoffdiagc = ', invoffdiagc(:,:)

  spin_index(1) = 1
  if(kp%nspin > 1) spin_index(2) = 2
  SAFE_ALLOCATE(mtxel, (ndiag + noffdiag, nspin, nkpt))
  SAFE_ALLOCATE(mtxel_sum, (ndiag + noffdiag, nspin, nkpt)) ! for testing
  mtxel_sum(:,:,:) = ZERO

! Let us actually calculate some stuff

  iexc = 1
  jexc = 1
  nvband_max = nvband
  ncband_max = ncband
!  nvband_max = 1
!  ncband_max = 2

  write(6,*) 'Computing <', iexc, '|dH^BSE|', jexc, '>'
  write(6,*) 'nvband_max = ', nvband_max, ' ncband_max = ', ncband_max

!  <Ai|H^BSE|Aj>

  SAFE_ALLOCATE(force_diag, (nmode))
  SAFE_ALLOCATE(force_offdiag, (nmode))
  SAFE_ALLOCATE(force_mode, (nmode))
  SAFE_ALLOCATE(force_atoms, (3, crys%nat))
  force_diag(:) = ZERO
  force_offdiag(:) = ZERO
  force_mode(:) = ZERO
  force_atoms(:,:) = ZERO

  if(peinf%inode == 0) write(6,'(a)') 'Beginning force calculation'

  mode_loop: do imode = 1, nmode
    write(filename, '(a,i6.6)') 'dfpt.dat_mode', imode
    if(peinf%inode == 0) write(6,'(a)') 'Reading matrix elements from ' // TRUNC(filename)
    do ik = 1, nkpt
      kk(:)=INF
      ierr=0
      call open_file(120,file=trim(filename),form='formatted',status='old',iostat=ierr)
      if(ierr /= 0) then
        write(0,'(a)') 'WARNING: mtxel file not found. Skipping.'
        cycle mode_loop
      endif
      do while (ierr.eq.0)
        call read_matrix_elements(120, ierr, kk, kp%nspin, ndiag, noffdiag, spin_index, diag, offdiag1, offdiag2, mtxel(:, :, ik))
        if (all(abs(kp%rk(1:3,ik)-kk(1:3)) .lt. TOL_Small)) exit
      enddo
      call close_file(120)

      ! Check k-point
      if(any(abs(kp%rk(1:3,ik)-kk(1:3)) .ge. TOL_Small)) then
        call die('cannot find k-point in ' // trim(filename), only_root_writes = .true.)
      endif
    enddo

!    do ib2 = 1, ndiag
!      do ib = 1, ndiag       
!        if(ib <= nvband .and. ib2 <= nvband) then
!          write(*,*) 'mtxelv(',ib,ib2,') = ', mtxel(invoffdiagv(ib,ib2),:,:)
!          write(*,*) 'invoffdiagv = ', invoffdiagv(ib,ib2)
!        endif
!        if(ib >  nvband .and. ib2 >  nvband) then
!          write(*,*) 'mtxelc(',ib,ib2,') = ', mtxel(invoffdiagc(ib-nvband,ib2-nvband),:,:)
!          write(*,*) 'invoffdiagc = ', invoffdiagc(ib-nvband,ib2-nvband)
!        endif
!      enddo
!    enddo

    mtxel_sum(:,:,:) = mtxel_sum(:,:,:) + mtxel(:,:,:)
    mtxel(:,:,:) = mtxel(:,:,:) / RYD ! convert from eV

    if(peinf%inode == 0) write(6,'(a)') 'Read matrix elements from ' // TRUNC(filename)

    ! could use some fancy matrix multiplication scheme here to speed up
    do ik = 1, nkpt
      do is = 1, nspin
        do ic = 1, ncband_max
          do iv = 1, nvband_max
            force_diag(imode) = force_diag(imode) + &
              (mtxel(invdiagc(ic), is, ik) - mtxel(invdiagv(iv), is, ik)) &
              * MYCONJG(asvck(iexc, is, iv, ic, ik)) * asvck(jexc, is, iv, ic, ik)
            
!            write(*,*) 'mtxelc = ', mtxel(invdiagc(ic), is, ik) * RYD, invdiagc(ic), 'ic = ', ic
!            write(*,*) 'mtxelv = ', mtxel(invdiagv(iv), is, ik) * RYD, invdiagv(iv), 'iv = ', iv
!            write(*,*) 'asvck = ', asvck(iexc, is, iv, ic, ik), asvck(jexc, is, iv, ic, ik)
!            write(*,*) 'force_diag = ', force_diag(imode)
            
!            icp = ic
            do ivp = 1, nvband_max
              if(ivp == iv) cycle
              ! we must avoid dividing by zero. at least for iexc = jexc, such terms should not contribute.
              if(abs(el_denom(vband(iv), ik, is) - el_denom(vband(ivp), ik, is)) < TOL_Degeneracy) cycle
              
              force_offdiag(imode) = force_offdiag(imode) + &
                MYCONJG(asvck(iexc, is, ivp, ic, ik)) * asvck(jexc, is, iv, ic, ik) * &
                (kp%el(vband(ivp), ik, is) - kp%el(vband(iv), ik, is)) * &
                MYCONJG(mtxel(invoffdiagv(ivp, iv), is, ik)) / (el_denom(vband(iv), ik, is) - el_denom(vband(ivp), ik, is))
!              write(*,*) 'v', ivp, el_denom(vband(iv), ik, is)*ryd, el_denom(vband(ivp), ik, is)*ryd, &
              !mtxel(invoffdiagv(ivp, iv), is, ik)*ryd, force_offdiag(imode) * RYD / BOHR
!                (kp%el(cband(ic), ik, is) - kp%el(vband(iv), ik, is) &
!                - kp%el(cband(ic), ik, is) + kp%el(vband(ivp), ik, is)) * &
            enddo

!            ivp = iv
            do icp = 1, ncband_max
              if(icp == ic) cycle
              if(abs(el_denom(cband(ic), ik, is) - el_denom(cband(icp), ik, is)) < TOL_Degeneracy) cycle
              
              force_offdiag(imode) = force_offdiag(imode) + &
                MYCONJG(asvck(iexc, is, iv, icp, ik)) * asvck(jexc, is, iv, ic, ik) * &
                (kp%el(cband(ic), ik, is) - kp%el(cband(icp), ik, is)) * &
                mtxel(invoffdiagc(icp, ic), is, ik) / (el_denom(cband(ic), ik, is) - el_denom(cband(icp), ik, is))
 !             write(*,*) 'c', icp, el_denom(cband(ic), ik, is)*ryd, el_denom(cband(icp), ik, is)*ryd, &
              !mtxel(invoffdiagc(icp, ic), is, ik)*ryd, force_offdiag(imode) * RYD / BOHR
            enddo

          enddo
        enddo 
      enddo ! is
    enddo ! ik

    ! convert to eV/Ang. minus sign is from F = -dE/dR
    force_mode(imode) = -(force_diag(imode) + force_offdiag(imode))
    write(6,'(4f12.6)') -force_diag(imode) * RYD / BOHR, -force_offdiag(imode) * RYD / BOHR
    write(6,*) 'Total force = ', force_mode(imode) * RYD / BOHR
  enddo mode_loop ! imode

  write(6,*)
  write(6,*) 'Atomic forces:'
  do iat = 1, crys%nat
    do ii = 1, 3
      force_atoms(ii, iat) = sum(disp(ii, iat, 1:nmode) * force_mode(1:nmode))
    enddo
    write(6,'(i4,999f12.6)') iat, dble(force_atoms(1:3, iat) * RYD/BOHR)
  enddo

  ! according to the acoustic sum rule, these should be zero
  call open_file(120,file='dfpt.dat_sum',form='formatted',status='replace')

  call write_matrix_elements(120, kk, kp%nspin, ndiag, noffdiag, spin_index, diag, &
    offdiag1, offdiag2, COMPLEXIFY(mtxel_sum(:, :, 1)))
  call close_file(120)

  call dealloc_header_type(sheader, crys, kp)
  call kp%free()
  SAFE_DEALLOCATE(mtxel)
  SAFE_DEALLOCATE(mtxel_sum)
  SAFE_DEALLOCATE(asvck)
  SAFE_DEALLOCATE(vband)
  SAFE_DEALLOCATE(cband)
  SAFE_DEALLOCATE(force_diag)
  SAFE_DEALLOCATE(force_offdiag)
  SAFE_DEALLOCATE(force_mode)
  SAFE_DEALLOCATE(force_atoms)
  SAFE_DEALLOCATE(diag)
  SAFE_DEALLOCATE(offdiag1)
  SAFE_DEALLOCATE(offdiag2)
  SAFE_DEALLOCATE(invdiagv)
  SAFE_DEALLOCATE(invdiagc)
  SAFE_DEALLOCATE(invoffdiagv)
  SAFE_DEALLOCATE(invoffdiagc)

end program forces
