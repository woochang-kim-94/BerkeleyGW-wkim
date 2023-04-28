!
! SIB: This routine looks the same as input().
! Except it reads from WFNq and
! writes to iunit_v='INT_VWFQ' and only valence bands.
! And k-point information is read into kpq.
!
!     SUBROUTINE READS CRYSTAL DATA AND WAVEFUNCTIONS FROM TAPE26
!     AND PARAMETERS FOR POLARIZABILITY CALCULATION FROM TAPE5
!     TAPE10 (OUTPUT TAPE) IS INITIALIZED

#include "f_defs.h"

module input_q_m

  use global_m
  use eqpcor_m
  use input_utils_m
  use misc_m
  use scissors_m
  use wfn_rho_vxc_io_m
  use io_utils_m
#ifdef HDF5
  use hdf5
#endif
  use hdf5_io_m
  use wfn_io_hdf5_m
  use timing_m, only: timing => epsilon_timing

  implicit none

  private
  public :: input_q

contains

subroutine input_q(gvec,kpq,cwfn,vwfn,pol,intwfnvq)
  type (gspace), intent(in) :: gvec
  type (kpoints), intent(out) :: kpq
  type (conduction_wfns), intent(in) :: cwfn
  type (valence_wfns), intent(in) :: vwfn
  type (polarizability), intent(inout) :: pol
  type (int_wavefunction), intent(out) :: intwfnvq

  type (crystal) :: crys
  type (symmetry) :: syms

  real(DP) :: vcell
  integer :: dummygvec(1, 1)
  character :: fncor*32

  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvecq

  PUSH_SUB(input_q)

  sheader = 'WFN'
  iflavor = 0
  call timing%start(timing%io_total)
  if(pol%wfn_hdf5) then
#ifdef HDF5
    if (peinf%inode==0) write(6,'(a)') ' Reading header of WFNq.h5'
    call read_hdf5_header_type('WFNq.h5', sheader, iflavor, kpq, gvecq, syms, crys)
#endif
  else
    if(peinf%inode == 0) then
      write(6,'(a)') ' Reading header of WFNq'
      call open_file(26,file='WFNq',form='unformatted',status='old')
    endif
    !call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, syms, crys, &
    !  dont_warn_kgrid=pol%subsample, warn=.false.)
    call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, syms, crys, &
      dont_warn_kgrid=.true., warn=.false.)

    call read_binary_gvectors(26, gvecq%ng, gvecq%ng, dummygvec, dont_read = .true.)
  endif
  call timing%stop(timing%io_total)
  call check_trunc_kpts(pol%icutv, kpq)

  call scissors_shift(kpq, pol%scis)
  call get_volume(vcell,crys%bdot)
  if (abs(crys%celvol-vcell).gt.TOL_Small) then
    call die('volume mismatch')
  endif

!-----------------------------------------------------------------
! If a quasi-particle correction file exists, read the corrected
! quasiparticle energies from file (in eV)

  if (pol%eqp_corrections) then
    fncor='eqp_q.dat'
    call eqpcor(fncor,peinf%inode,peinf%npes,kpq,1+vwfn%ncore_excl,vwfn%nband+vwfn%ncore_excl+pol%ncrit,0,0,kpq%el,kpq%el,&
                 kpq%el,1,0)
  endif

  call find_efermi(pol%rfermi, pol%efermi, pol%efermi_input, kpq, kpq%mnband, 1+vwfn%ncore_excl, &
    "shifted grid", should_search = .false., should_update = .false., write7 = .true.)

  if (peinf%inode==0) then
    if(any (kpq%ifmax(:,:) < vwfn%nband+vwfn%ncore_excl .or. kpq%ifmax(:,:) > vwfn%nband +vwfn%ncore_excl+ pol%ncrit)) then
      write(0,'(a,i6,a,i6,a)') 'epsilon.inp says there are ', vwfn%nband, ' fully occupied bands and ', &
        pol%ncrit, ' partially occupied.'
      write(0,'(a,2i6)') 'This is inconsistent with highest bands in WFNq file; min, max = ', minval(kpq%ifmax), maxval(kpq%ifmax)
      call die("band_occupation, number_partial_occup, and WFNq inconsistent.")
    endif

    if(maxval(kpq%ifmax) - minval(kpq%ifmax) > pol%ncrit) then
      write(0,'(a,i6,a)') 'epsilon.inp says there are ', pol%ncrit, ' partially occupied bands.'
      write(0,'(a,i6)') 'This is less than the number partially occupied in WFNq file: ', maxval(kpq%ifmax) - minval(kpq%ifmax)
      call die("number_partial_occup and WFNq inconsistent.")
    endif
  endif

  call timing%start(timing%io_total)
  if (.not.pol%skip_chi) then
    if (pol%wfn_hdf5) then
#ifdef HDF5
      call read_hdf5_wavefunctions(kpq, gvec, pol, cwfn, vwfn, intwfnvq)
#endif
    else
      call read_wavefunctions(kpq, gvec, pol, cwfn, vwfn, intwfnvq)
      if (peinf%inode.eq.0) then
        call close_file(26)
      endif
    endif
  endif
  call timing%stop(timing%io_total)
  POP_SUB(input_q)

  return
end subroutine input_q

#ifdef HDF5
subroutine read_hdf5_wavefunctions(kpq, gvec, pol, cwfn, vwfn, intwfnvq)
  type (kpoints), intent(in) :: kpq
  type (gspace), intent(in) :: gvec
  type (polarizability), intent(in) :: pol
  type (conduction_wfns), intent(in) :: cwfn
  type (valence_wfns), intent(in) :: vwfn
  type (int_wavefunction), intent(out) :: intwfnvq

  integer, allocatable :: isort(:)
  integer, allocatable :: components(:,:)
  SCALAR, allocatable :: wfns(:,:,:)

  integer :: ig,ik,iiii,ib,is
  real(DP) :: qk(3)
  type(gspace) :: gvec_kpt

  integer(HID_T) :: file_id
  integer(HID_T) :: plist_id
  integer :: error
  integer :: ngktot
  integer :: istart, ib_first

  PUSH_SUB(read_hdf5_wavefunctions)

  SAFE_ALLOCATE(intwfnvq%ng, (kpq%nrk))
  SAFE_ALLOCATE(intwfnvq%isort, (kpq%ngkmax,kpq%nrk))
  SAFE_ALLOCATE(intwfnvq%cg, (kpq%ngkmax,kpq%nrk*peinf%nvownactual,kpq%nspin*kpq%nspinor))
  SAFE_ALLOCATE(intwfnvq%qk, (3,kpq%nrk))

  ngktot = SUM(kpq%ngk)
  SAFE_ALLOCATE(gvec_kpt%components, (3, ngktot))

  call logit('Reading G-vectors')
  call read_hdf5_wfn_gvectors('WFNq.h5', gvec_kpt%components, ngktot)

  istart=1
  do ik= 1, kpq%nrk
    qk(:)=kpq%rk(:,ik)

    SAFE_ALLOCATE(isort, (kpq%ngk(ik)))
    SAFE_ALLOCATE(components, (3,kpq%ngk(ik)))

    components(:,:) = gvec_kpt%components(:,istart:istart+kpq%ngk(ik)-1)
    istart=istart+kpq%ngk(ik)

    do ig = 1, kpq%ngk(ik)
      call findvector(isort(ig), components(:, ig), gvec)
      if (isort(ig) ==  0)  call die('input: could not find gvec')
    enddo

    intwfnvq%ng(ik)=kpq%ngk(ik)
    intwfnvq%isort(1:kpq%ngk(ik),ik)=isort(1:kpq%ngk(ik))
    intwfnvq%qk(:,ik)=qk(:)

    SAFE_DEALLOCATE(isort)
    SAFE_DEALLOCATE(components)
  end do

  SAFE_DEALLOCATE_P(gvec_kpt%components)

  ! HDF5 inteface setup
  call hdf5_open_file('WFNq.h5', 'r', file_id, parallel_io=.true.)

  call logit('Reading shifted valence bands')
  ib_first = 0
  if (peinf%nvownactual>0) ib_first = peinf%invindexv(1)
  SAFE_ALLOCATE(wfns, (ngktot,kpq%nspin*kpq%nspinor,peinf%nvownactual))
  call read_hdf5_bands_block(file_id, kpq, peinf%nvownmax, peinf%nvownactual, &
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
    do ik = 1, kpq%nrk
      iiii=peinf%indexv(peinf%invindexv(ib))+(ik-1)*peinf%nvownactual
      do is=1,kpq%nspin*kpq%nspinor
        intwfnvq%cg(1:kpq%ngk(ik),iiii,is)=wfns(istart:istart+kpq%ngk(ik)-1,is,ib)
      enddo
      do is = 1, kpq%nspin
        call checknorm('WFNq',peinf%invindexv(ib),ik,kpq%ngk(ik),is,kpq%nspinor,intwfnvq%cg(1:kpq%ngk(ik),iiii,:))
      enddo
      ! for metals
      ! FHJ: FIXME: isn`t the following block redundant?
      if(peinf%invindexv(ib).le.vwfn%nband+pol%ncrit) then
        iiii=peinf%indexv(peinf%invindexv(ib))+(ik-1)*peinf%nvownactual
        do is=1,kpq%nspin*kpq%nspinor
          intwfnvq%cg(1:kpq%ngk(ik),iiii,is)=wfns(istart:istart+kpq%ngk(ik)-1,is,ib)
        enddo
      end if
      istart=istart+kpq%ngk(ik)
    enddo
  enddo
  SAFE_DEALLOCATE(wfns)

  call hdf5_close_file(file_id)
  call logit('Finished reading HDF5 wavefunctions')

  POP_SUB(read_hdf5_wavefunctions)
  return
end subroutine read_hdf5_wavefunctions
#endif

subroutine read_wavefunctions(kpq, gvec, pol, cwfn, vwfn, intwfnvq)
  type (kpoints), intent(in) :: kpq
  type (gspace), intent(in) :: gvec
  type (polarizability), intent(in) :: pol
  type (conduction_wfns), intent(in) :: cwfn
  type (valence_wfns), intent(in) :: vwfn
  type (int_wavefunction), intent(out) :: intwfnvq

  integer, allocatable :: isort(:)
  SCALAR, allocatable :: zc(:,:)

  character :: filenamevq*20
  integer :: i,i2,j,k,ik,iiii,is
  integer :: iunit_v
  real(DP) :: qk(3)
  logical :: dont_read

  type(gspace) :: gvec_kpt
  type(progress_info) :: prog_info !< a user-friendly progress report

  PUSH_SUB(read_wavefunctions)

  SAFE_ALLOCATE(intwfnvq%ng, (kpq%nrk))
  SAFE_ALLOCATE(intwfnvq%isort, (kpq%ngkmax,kpq%nrk))
  SAFE_ALLOCATE(intwfnvq%cg, (kpq%ngkmax,kpq%nrk*peinf%nvownactual,kpq%nspin*kpq%nspinor))
  SAFE_ALLOCATE(intwfnvq%qk, (3,kpq%nrk))

  call progress_init(prog_info, 'reading wavefunctions (WFNq)', 'state', kpq%nrk*(vwfn%nband+pol%ncrit))
  do ik=1,kpq%nrk
    qk(1:3) = kpq%rk(1:3, ik)
    SAFE_ALLOCATE(gvec_kpt%components, (3, kpq%ngk(ik)))

    call read_binary_gvectors(26, kpq%ngk(ik), kpq%ngk(ik), gvec_kpt%components)

    SAFE_ALLOCATE(isort, (kpq%ngk(ik)))
    do i = 1, kpq%ngk(ik)
      call findvector(isort(i), gvec_kpt%components(:, i), gvec)
      if (isort(i) == 0) then
        if(peinf%inode == 0) write(0,*) 'ik = ', ik, 'ig = ', i, 'gvec = ', gvec_kpt%components(:, i)
        call die('input_q: could not find gvec')
      endif
    enddo
    SAFE_DEALLOCATE_P(gvec_kpt%components)

    intwfnvq%ng(ik)=kpq%ngk(ik)
    intwfnvq%isort(1:kpq%ngk(ik),ik)=isort(1:kpq%ngk(ik))
    intwfnvq%qk(:,ik)=qk(:)
!
! SIB:  loop on max number of bands, and proc 0 reads the wave function
! from unit 26, checks normalization, and if the band is less than
! cwfn%nband (# of bands in total) **AND** is a valence band (so
! its index is <= vwfn%nband), then it is written to iunit_v.
!
    SAFE_ALLOCATE(zc, (kpq%ngk(ik), kpq%nspinor*kpq%nspin))

    do i=1,kpq%mnband

      dont_read = (i > cwfn%nband .or. i <= vwfn%ncore_excl)
      if(.not. dont_read) dont_read = i > vwfn%nband+vwfn%ncore_excl+pol%ncrit
      call read_binary_data(26, kpq%ngk(ik), kpq%ngk(ik), kpq%nspin*kpq%nspinor, zc, dont_read = dont_read)

      ! FHJ: the following lines were introduced in r6294 and are supposed to
      ! be a shortcut if we are past the last band of the last k-point. However,
      ! in light of a previous bug (#223), this feature is commented out for now.
      !! FHJ: shortcut if this is past the last band of the last k-point
      !if (dont_read .and. ik==kpq%nrk) exit
      if (.not.dont_read) then
        ! DVF: recall that we redefined the number of valence bands to exclude the
        ! core states. So, we have to subtract ncore_excl right here because ib is
        ! referenced to the full wavefunction file including all the core states.
        i2=i-vwfn%ncore_excl
        call progress_step(prog_info, (ik-1)*(vwfn%nband+pol%ncrit) + i)
        if (peinf%inode == 0) then
          do is = 1, kpq%nspin
            call checknorm('WFNq',i,ik,kpq%ngk(ik),is,kpq%nspinor,zc(:,:))
          enddo
        endif

        if (peinf%doiownv(i2)) then
          iiii=peinf%indexv(i2)+(ik-1)*peinf%nvownactual
          intwfnvq%cg(1:kpq%ngk(ik),iiii,1:kpq%nspin*kpq%nspinor)=zc(1:kpq%ngk(ik),1:kpq%nspinor*kpq%nspin)
        endif
      endif
    enddo
    SAFE_DEALLOCATE(isort)
    SAFE_DEALLOCATE(zc)
  enddo                     ! end loop over k+q points
  call progress_free(prog_info)

  POP_SUB(read_wavefunctions)
  return

end subroutine read_wavefunctions

end module input_q_m
