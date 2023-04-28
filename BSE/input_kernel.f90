#include "f_defs.h"

module input_kernel_m

  use checkbz_m
  use distrib_kernel_m
  use epsread_hdf5_m
  use fullbz_m
  use global_m
#ifdef HDF5
  use hdf5
#endif
  use io_utils_m
  use input_utils_m
  use misc_m
  use read_rho_vxc_m
  use sort_m
  use wfn_rho_vxc_io_m
  use wfn_io_hdf5_m
  implicit none

  private

  public :: &
    input_kernel

contains

!-----------------------------------------------------------------------
subroutine input_kernel(crys,gvec,kg,kgq,kp,syms,xct,flagbz, &
  intwfnv,intwfnc)
!-----------------------------------------------------------------------
!
!     Read parameters from file WFN_co
!     Initialize k-points sampling, kg type
!     Initialize G-space, gvec
!
!     input: xct type
!
!     output: crys,gvec,syms,kg types
!             INT_VWFN_* and INT_CWFN_* files
!
  type (crystal), intent(out) :: crys
  type (gspace), intent(out) :: gvec
  type (grid), intent(out) :: kg,kgq
  type (kpoints), intent(out) :: kp
  type (symmetry), intent(out) :: syms
  type (xctinfo), intent(inout) :: xct
  integer, intent(in) :: flagbz
  type (int_wavefunction), intent(out) :: intwfnc,intwfnv


  type (wavefunction) :: wfnc,wfnv
  character :: tmpfn*16
  integer :: iwritev,iwritec,iwritek
  integer :: ickmem,irk
  integer :: ii,jj,is,isp,ig,ikq,ik,umk
  integer :: irks,ivband,icband

  real(DP) :: diffvol,vcell,kt(3),div,tol,delta,qq_temp(3)
  
  real(DP), allocatable :: ek_tmp(:)
  integer, allocatable :: index(:),indxk(:),indxkq(:),k_tmp(:,:)
  integer, allocatable :: isrti(:)
  SCALAR, allocatable :: cg(:,:)

  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvec_kpt
  logical :: skip_checkbz, need_band
  type(progress_info) :: prog_info

! DYQ: Variables needed for WFNq_co
  type (kpoints) :: kpq
  type (crystal) :: crysq
  type (gspace) :: gvecq
  type (symmetry) :: symsq

  ! WFN HDF5 stuff
#ifdef HDF5
  integer(HID_T) :: file_id
  integer(HID_T) :: plist_id
#endif
  integer, allocatable :: gvec_kpt_components_all(:,:)
  integer, allocatable :: ib_size_array(:)
  integer :: istart, ngktot
  integer :: iband_min, iband_max, tot_nbands, band_block, min_band_block
  integer :: ib_first, ib_last, ib_size, ib_size_max
  integer :: ipe, npes_hdf5, ii_start, ii_end, ii_loc
  SCALAR, allocatable :: my_cg_all(:,:,:), ipe_cg_all(:,:,:)
  integer :: error

  PUSH_SUB(input_kernel)

  call logit('input_kernel:  reading WFN_co')
  sheader = 'WFN'
  iflavor = 0
  if ( xct%use_wfn_hdf5 ) then
    if (peinf%inode==0) write(6,'(a)') ' Reading header of WFN_co.h5'
#ifdef HDF5
    call read_hdf5_header_type('WFN_co.h5', sheader, iflavor, kp, gvec, syms, crys)
#endif
  else
    if (peinf%inode == 0) call open_file(25,file='WFN_co',form='unformatted',status='old')
    call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys, dont_warn_kgrid=.true.)
  end if
  call check_trunc_kpts(xct%icutv, kp)

  if(any(kp%shift > TOL_Zero) .and. peinf%inode == 0) then
    write(0,'(a)') "WARNING: WFN_co has a shift. This is not recommended."
  endif

  call logit('input_kernel:  reading gvec info')
  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  if( xct%use_wfn_hdf5 ) then
#ifdef HDF5
    call read_hdf5_gvectors('WFN_co.h5', gvec%ng, gvec%components)
#endif
  else
    call read_binary_gvectors(25, gvec%ng, gvec%ng, gvec%components)
  end if

  call get_volume(vcell,crys%bdot)
  diffvol=abs(crys%celvol-vcell)
  if (diffvol.gt.TOL_Small) then
    call die('volume mismatch', only_root_writes = .true.)
  endif

  ! there is no fine grid to set the Fermi level in kernel
  call find_efermi(xct%rfermi, xct%efermi, xct%efermi_input, kp, kp%mnband, 1, &
    "coarse grid", should_search = .true., should_update = .true., write7 = .false.)

  call assess_degeneracies(kp, kp%el(kp%mnband, :, :), kp%mnband - 1, xct%efermi, TOL_Degeneracy)

  xct%nspin = kp%nspin
  xct%nspinor = kp%nspinor

  if(any(kp%ifmax(:,:) == 0)) &
    call die("BSE codes cannot handle a system where some k-points have no occupied bands.", only_root_writes = .true.)

  kp%nvband=minval(kp%ifmax(:,:)-kp%ifmin(:,:))+1
  kp%ncband=kp%mnband-maxval(kp%ifmax(:,:))

!----------------------------------------------------------------
! (gsm) check whether the requested number of bands
!       is available in the wavefunction file

  if(xct%nvb_co.gt.kp%nvband) then
    call die("The requested number of valence bands is not available in WFN_co.", &
             only_root_writes=.true.)
  endif
  if(xct%ncb_co.gt.kp%ncband) then
    call die("The requested number of conduction bands is not available in WFN_co.", &
             only_root_writes=.true.)
  endif

! DAS: degenerate subspace check

  if (peinf%inode.eq.0) then
      ! MHN
    !if(xct%ncb_co.eq.kp%ncband) then
    !  call die("You must provide one more conduction band in WFN_co in order to assess degeneracy.", &
    !           only_root_writes=.true.)
    !endif
    do jj = 1, kp%nspin
      do ii = 1, kp%nrk
        if(kp%ifmax(ii, jj) - xct%nvb_co > 0) then
          ! no need to compare against band 0 if all valence are included
          if(abs(kp%el(kp%ifmax(ii, jj) - xct%nvb_co + 1, ii, jj) &
            - kp%el(kp%ifmax(ii, jj) - xct%nvb_co, ii, jj)) .lt. TOL_Degeneracy) then
            if(xct%degeneracy_check_override) then
              write(0,'(a)') &
                "WARNING: Selected number of valence bands breaks degenerate subspace in WFN_co. " // &
                "Run degeneracy_check.x for allowable numbers."
              write(0,*)
            else
              write(0,'(a)') &
                "Run degeneracy_check.x for allowable numbers, or use keyword " // &
                "degeneracy_check_override to run anyway (at your peril!)."
              call die("Selected number of valence bands breaks degenerate subspace in WFN_co.", &
                       only_root_writes=.true.)
            endif
          endif
        endif
        if(abs(kp%el(kp%ifmax(ii, jj) + xct%ncb_co, ii, jj) &
          - kp%el(kp%ifmax(ii, jj) + xct%ncb_co + 1, ii, jj)) .lt. TOL_Degeneracy) then 
          if(xct%degeneracy_check_override) then
            write(0,'(a)') &
              "WARNING: Selected number of conduction bands breaks degenerate subspace in WFN_co. " // &
              "Run degeneracy_check.x for allowable numbers."
            write(0,*)
          else
            write(0,'(a)') &
              "Run degeneracy_check.x for allowable numbers, or use keyword " // &
              "degeneracy_check_override to run anyway (at your peril!)."
            call die("Selected number of conduction bands breaks degenerate subspace in WFN_co.",&
                     only_root_writes=.true.)
          endif
        endif
      enddo
    enddo
  endif

!-----------------------------------------------------------------------
!     Read the k-point sampling from kpoints (if it exists) or from
!     WFN_co

  if (xct%read_kpoints) then
    if (peinf%inode.eq.0) then
      call open_file(9,file='kpoints_co',form='formatted',status='old')
      read(9,*) kg%nr
      SAFE_ALLOCATE(kg%r, (3,kg%nr))
      do ii=1,kg%nr
        read(9,*) (kg%r(jj,ii),jj=1,3),div
        kg%r(:,ii) = kg%r(:,ii)/div
      enddo
      call close_file(9)
    endif
#ifdef MPI
    call MPI_BCAST(kg%nr,   1,     MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if(peinf%inode.ne.0) then
      SAFE_ALLOCATE(kg%r, (3,kg%nr))
    endif
    call MPI_BCAST(kg%r,    3*kg%nr,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif

!----------------------------------------------------------------
!     indxk : stores the correspondence between k-points kg%r and kp%rk
!     (it is used to select the set of wavefunctions to be stored)
!     tol : tolerance in the coordinates of k-points

    tol = 1.d-4
    SAFE_ALLOCATE(indxk, (kg%nr))
    indxk=0
    do jj=1,kg%nr
      do ii=1,kp%nrk
        kt(:) = kg%r(:,jj) - kp%rk(:,ii)
        if ((abs(kt(1)).lt.tol).and.(abs(kt(2)).lt.tol) &
          .and.(abs(kt(3)).lt.tol)) then
          if (indxk(jj).ne.0) write(0,*) 'WARNING: multiple definition of k-point',jj,indxk(jj),kg%r(:,jj)
          indxk(jj)=ii
        endif
      enddo
      if (indxk(jj).eq.0) write(0,*) 'WARNING: could not find vector ',kg%r(:,jj),' in WFN_co'
!
!     no need to stop here; if indxk.eq.0, the job will stop in genwf
!
    enddo
  else
    kg%nr=kp%nrk
    SAFE_ALLOCATE(kg%r, (3,kg%nr))
    kg%r(1:3,1:kg%nr)=kp%rk(1:3,1:kp%nrk)
    SAFE_ALLOCATE(indxk, (kg%nr))
    do ii=1,kg%nr
      indxk(ii) = ii
    enddo
  endif

!-----------------------------------------------------------------------
!     Order g-vectors with respect to their kinetic energy
!

  call logit('input_kernel:  reordering gvecs')

  ! FHJ: Figure out ecute and ecutg
  call get_ecut()
  
  SAFE_ALLOCATE(index, (gvec%ng))
  SAFE_ALLOCATE(gvec%ekin, (gvec%ng))
  call kinetic_energies(gvec, crys%bdot, gvec%ekin)
  call sortrx(gvec%ng, gvec%ekin, index, gvec = gvec%components)

  xct%ng   = gcutoff(gvec%ng, gvec%ekin, index, xct%ecutg)
  if (xct%theory.eq.0) &
  xct%neps = gcutoff(gvec%ng, gvec%ekin, index, xct%ecute) ! cuteness energy??

  if(xct%ecute > xct%ecutg .and. xct%theory.eq.0) then
    write(0,*) 'ecute = ', xct%ecute, ' ecutg = ', xct%ecutg
    call die("The screened_coulomb_cutoff cannot be greater than the bare_coulomb_cutoff.", only_root_writes = .true.)
  endif

  SAFE_ALLOCATE(ek_tmp, (gvec%ng))
  ek_tmp = gvec%ekin
  SAFE_ALLOCATE(k_tmp, (3,gvec%ng))
  k_tmp = gvec%components
  do ii=1,gvec%ng
    gvec%ekin(ii) = ek_tmp(index(ii))
    gvec%components(:,ii) = k_tmp(:,index(ii))
  enddo

  call gvec_index(gvec)

  ! If we are not doing just purely TDHF then
  ! read the charge density/fxc
  if ((xct%theory == 1) .and. ((1.0d0 - xct%coulomb_mod%long_range_frac_fock > TOL_SMALL) .or. &
     (1.0d0 - xct%coulomb_mod%short_range_frac_fock > TOL_SMALL))) then
    xct%coul_mod_flag=.true.
    SAFE_ALLOCATE(isrti, (gvec%ng))
    do ii=1,gvec%ng
      isrti(index(ii)) = ii
    enddo
    call read_rho(xct%wpg, gvec, kp, syms, crys, isrti, index, 'WFN_co')
    SAFE_DEALLOCATE(isrti)
  endif

  SAFE_DEALLOCATE(index)
  SAFE_DEALLOCATE(ek_tmp)
  SAFE_DEALLOCATE(k_tmp)

  
!-----------------------------------------------------------------------
!     Generate full brillouin zone from irreducible wedge, rk -> fk
!
!     If flagbz.eq.1, only Identity will be used as
!     symmetry operation. In this case, kg%r (irreducible BZ) and kg%f
!     (full BZ) will be identical.
!
  if (flagbz.eq.0.and.peinf%inode.eq.0) write(6,801)
  if (flagbz.eq.1.and.peinf%inode.eq.0) write(6,802)
801 format(1x,'Using symmetries to expand the coarse grid sampling')
802 format(1x,'No symmetries used in the coarse grid sampling')
!
  call timacc(7,1)
  if (flagbz.eq.1) then
    call fullbz(crys,syms,kg,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  else
    call fullbz(crys,syms,kg,syms%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  endif
  tmpfn='WFN_co'
  if (.not. skip_checkbz.and..not.xct%patched_sampling_co) then
    call checkbz(kg%nf,kg%f,kp%kgrid,kp%shift,crys%bdot, &
      tmpfn,'k',.true.,xct%freplacebz,xct%fwritebz)
  endif
  call timacc(7,2)
  xct%nkpt_co=kg%nf

  if (peinf%verb_high .and. peinf%inode==0) then
    write(6,'(/1x,a6,14x,a7,12x,2(1x,a6),3x,a3)') 'i', 'k-point', 'indr', 'itran', 'kg0'
    write(6,'(1x,6("-"),1x,32("-"),2(1x,6("-")),1x,8("-"))')
    do ii=1,kg%nf
      write(6,'(1x,i6,3(1x,f10.6),2(1x,i6),3(1x,i2))') &
        ii, kg%f(:,ii), kg%indr(ii), kg%itran(ii), kg%kg0(:,ii)
    enddo
  endif

!------------------------------------------------------------------------
! If there is a finite center-of-mass momentum, Q, find mapping between k
! and k+Q

  SAFE_ALLOCATE(xct%indexq,(kg%nf))
  if (xct%qflag.eq.1) then
    do ik=1,kg%nf
      xct%indexq(ik) = ik
    enddo
  endif

  SAFE_ALLOCATE(kgq%f,(3,kg%nf))
  SAFE_ALLOCATE(kgq%kg0,(3,kg%nf))
  SAFE_ALLOCATE(kgq%indr,(kg%nf))
  ! DYQ: Use a finite momentum that corresponds to a small shift, when qflag=1.
  ! Shifted wave functions will be read from WFNq_co.
  if (xct%qflag.eq.0) then

    sheader = 'WFN'
    iflavor = 0
    ! Read the header of WFNq_co
    if ( xct%use_wfn_hdf5 ) then
      if (peinf%inode==0) write(6,'(a)') ' Reading header of WFNq_co.h5'
#ifdef HDF5
      call read_hdf5_header_type('WFNq_co.h5', sheader, iflavor, kpq, gvecq, symsq, crysq)
#endif
    else
      if (peinf%inode == 0) call open_file(26,file='WFNq_co',form='unformatted',status='old')
      call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, symsq, crysq, dont_warn_kgrid=.true.)
    end if
    call check_header('WFN_co', kp, gvec, syms, crys, 'WFNq_co', kpq, gvecq, symsq, crysq, is_wfn = .true.)

    if(any(kp%kgrid(1:3) /= kpq%kgrid(1:3))) then
      if(peinf%inode == 0) then
        write(0,*) 'WFN_co  kgrid = ', kp%kgrid(1:3)
        write(0,*) 'WFNq_co kgrid = ', kpq%kgrid(1:3)
      endif
      call die('kgrids for WFN_co and WFNq_co must be the same', only_root_writes = .true.)
    endif

    SAFE_ALLOCATE(gvecq%components, (3, gvecq%ng))
    if( xct%use_wfn_hdf5 ) then
#ifdef HDF5
      call read_hdf5_gvectors('WFNq_co.h5', gvecq%ng, gvecq%components)
#endif
    else
      call read_binary_gvectors(26, gvecq%ng, gvecq%ng, gvecq%components)
    end if

    call get_volume(vcell,crysq%bdot)
    diffvol=abs(crysq%celvol-vcell)
    if (diffvol.gt.TOL_Small) then
      call die('volume mismatch', only_root_writes = .true.)
    endif

    call assess_degeneracies(kpq, kpq%el(kpq%mnband, :, :), kpq%mnband - 1, xct%efermi, TOL_Degeneracy)

    if(any(kpq%ifmax(:,:) == 0)) &
      call die("BSE codes cannot handle a system where some k-points have no occupied bands.", only_root_writes = .true.)

    kpq%nvband=minval(kpq%ifmax(:,:)-kpq%ifmin(:,:))+1
    kpq%ncband=kpq%mnband-maxval(kpq%ifmax(:,:))

    kgq%nr=kpq%nrk
    SAFE_ALLOCATE(kgq%r, (3,kgq%nr))
    kgq%r(1:3,1:kgq%nr)=kpq%rk(1:3,1:kpq%nrk)
    SAFE_ALLOCATE(indxkq, (kgq%nr))
    do ii=1,kgq%nr
      indxkq(ii) = ii
    enddo

    ! Check that finite Q is the same as the difference between kgq and kg

    ! DYQ TODO: fix symmetries. For now, symmetries are not used for WFNq_co
    if (.true.) then
      call fullbz(crysq,symsq,kgq,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
    else
      call fullbz(crysq,symsq,kgq,symsq%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
    endif
    
#ifdef DEBUG
  if (peinf%inode .eq. 0) then
    write(6,*) 'kgq%nf',kgq%nf
  endif
#endif
    ! Find mapping between kgq and kg
    do ik=1,kg%nf
      ikq = 0
      delta = 0.1d0
      do while ((delta .gt. TOL_Small) .and. (ikq.lt.kgq%nf))
        ikq = ikq+1
        qq_temp(:) = kgq%f(:,ikq) - kg%f(:,ik) - xct%finiteq(:) 
        do jj=1,3
          qq_temp(jj) = qq_temp(jj) - anint( qq_temp(jj) )
        enddo
        delta=sqrt((qq_temp(1))**2+(qq_temp(2))**2+(qq_temp(3))**2)
      enddo
      if (delta .gt. TOL_Small) then 
        if(peinf%inode.eq.0) then
          write(0,*) '  Could not find point equivalent to ', (kg%f(ii,ik),ii=1,3)
        endif
        call die('Finite momentum not commensurate with kgrid of WFN_co',only_root_writes = .true.)
      endif
      xct%indexq(ik)=ikq
      ! kg%f and (kgq%f-xct%finiteq) may differ by a lattice vector near the zone edge 
      do jj=1,3
        umk = nint(kg%f(jj,ik) - kgq%f(jj,ikq) + xct%finiteq(jj) )
        kgq%f(jj,ikq) = kgq%f(jj,ikq) + dble(umk) !kgq(indexq(ik)) = kg(ikq) 
        kgq%kg0(jj,ikq) = kgq%kg0(jj,ikq) + umk
      enddo
    enddo
  endif

  ! DYQ: Use a center-of-mass momentum commensurate with the kgrid of WFN_co, when qflag=2 
  ! Find the index of the point in kg corresponding to kg%f+xct%finiteq 
  ! and store the index mapping in indexq. kgq is the grid kg shifted by Q.
  if (xct%qflag .eq. 2) then
    do ik=1,kg%nf
      kgq%f(:,ik) = kg%f(:,ik) + xct%finiteq(:) !k+Q
      ikq = 0
      delta=0.1d0
      do while ((delta .gt. TOL_Small) .and. (ikq.lt.kg%nf))
        ikq = ikq+1
        qq_temp(:) = kgq%f(:,ik) - kg%f(:,ikq)
        do jj=1,3
          qq_temp(jj) = qq_temp(jj) - anint( qq_temp(jj) )
        enddo
        delta=sqrt((qq_temp(1))**2+(qq_temp(2))**2+(qq_temp(3))**2)
      enddo
      ! With patched sampling, k+Q might fall outside the patch
      if (delta.gt.TOL_Small .and. xct%patched_sampling_co) then
        if(peinf%inode.eq.0) then
          write(6,*) '  Could not find point equivalent to ', (kg%f(ii,ik),ii=1,3)
          write(6,*) '  Skipping point.'
        endif
        xct%indexq(ik)= 0
        kgq%indr(ik) = ik
      elseif (delta.gt.TOL_Small) then
        if(peinf%inode.eq.0) then
          write(0,*) '  Could not find point equivalent to ', (kg%f(ii,ik),ii=1,3)
        endif
        call die('Finite momentum not commensurate with kgrid of WFN_co',only_root_writes = .true.)
      else
        xct%indexq(ik)=ikq
        ! kg%f and kgq%f may differ by a lattice vector near the zone edge     
        do jj=1,3
          umk = nint( kgq%f(jj,ik)-kg%f(jj,ikq))
          !kgq%f(jj,ik) = kgq%f(jj,ik) - dble(umk) !kgq(ik) = kg(indexq(ik))   
          kgq%kg0(jj,ik) = kg%kg0(jj,ikq) + umk
        enddo
        kgq%indr(ik)=kg%indr(ikq)
      endif
    enddo
  endif    

!------------------------------------------------
! Distribute vcks-quadruplets among the PEs

  call logit('input_kernel:  calling distrib_kernel')
  if (xct%qflag.eq.1) then
    call distrib_kernel(xct,kp%ngkmax,kg,kg,gvec)
  elseif (xct%qflag.ne.1) then
    ! DYQ: a processor that owns ik cond. bands will also own ikq val. bands
    call distrib_kernel(xct,kp%ngkmax,kg,kgq,gvec)
  endif
  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,'(/1x,a,3(1x,i0)/)') "Allocating isort", peinf%iownwfv(peinf%inode+1), peinf%iownwfc(peinf%inode+1), gvec%ng
  endif
  if (xct%qflag.eq.1) then
    SAFE_ALLOCATE(intwfnv%cg, (kp%ngkmax,peinf%iownwfv(peinf%inode+1),kp%nspin*kp%nspinor))
  else
    SAFE_ALLOCATE(intwfnv%cg, (max(kp%ngkmax,kpq%ngkmax),peinf%iownwfv(peinf%inode+1),kp%nspin*kp%nspinor))
  endif
  SAFE_ALLOCATE(intwfnc%cg, (kp%ngkmax,peinf%iownwfc(peinf%inode+1),kp%nspin*kp%nspinor))
  if (xct%qflag.eq.1) then
    SAFE_ALLOCATE(intwfnv%isort, (gvec%ng,peinf%iownwfk(peinf%inode+1)))
    SAFE_ALLOCATE(intwfnv%ng, (peinf%iownwfk(peinf%inode+1)))
  else
    SAFE_ALLOCATE(intwfnv%isort, (gvec%ng,peinf%iownwfkq(peinf%inode+1)))
    SAFE_ALLOCATE(intwfnv%ng, (peinf%iownwfkq(peinf%inode+1)))
  endif
  SAFE_ALLOCATE(intwfnc%isort, (gvec%ng,peinf%iownwfk(peinf%inode+1)))
  SAFE_ALLOCATE(intwfnc%ng, (peinf%iownwfk(peinf%inode+1)))
  intwfnv%nspin=kp%nspin
  intwfnv%nspinor=kp%nspinor
  intwfnc%nspin=kp%nspin
  intwfnc%nspinor=kp%nspinor

#ifdef DEBUG
  if (peinf%inode.eq.0) then
    write(6,*) 'kp%ngkmax',kp%ngkmax 
  endif
#endif

!------------------------------------------------
! Begin loop that distributes wave functions

  SAFE_ALLOCATE(wfnv%isort, (gvec%ng))
  SAFE_ALLOCATE(wfnc%isort, (gvec%ng))
  wfnv%nspin=kp%nspin
  wfnv%nspinor=kp%nspinor
  wfnc%nspin=kp%nspin
  wfnc%nspinor=kp%nspinor

  if ( xct%use_wfn_hdf5  ) then
    if (peinf%inode==0) write(6,*)
    if (peinf%inode==0) write(6,'(a)') ' Reading HDF5 wavefuntion (WFN_co.h5)'
    ngktot = SUM(kp%ngk)
    SAFE_ALLOCATE( gvec_kpt_components_all, (3,ngktot) )
#ifdef HDF5
    call read_hdf5_wfn_gvectors('WFN_co.h5', gvec_kpt_components_all, ngktot)
#endif

   iband_min      = MINVAL(kp%ifmax(:,:)) - xct%nvb_co + 1
   iband_min      = MAX(iband_min, 1)
   iband_max      = MAXVAL(kp%ifmax(:,:)) + xct%ncb_co
   iband_max      = MIN(iband_max, kp%mnband)
   tot_nbands     = iband_max - iband_min + 1
   band_block     = (tot_nbands + peinf%npes - 1) / peinf%npes
   ! read at least 128 Mb per MPI task
   min_band_block = 128.0D+00 / ( dble(ngktot) * dble(kp%nspin*kp%nspinor) * 16.0D+00 /1024.0D+00/1024D+00 )
   min_band_block = MAX(min_band_block, 1)
   if ( xct%wfn_hdf5_min_band_block > 0 )  min_band_block = xct%wfn_hdf5_min_band_block
   if ( min_band_block > band_block ) then
     band_block = min_band_block
   end if

   ib_first = iband_min + band_block * peinf%inode
   ib_last  = min(ib_first + band_block - 1, iband_max)
   ib_size  = ib_last - ib_first + 1
   if ( ib_size < 1 ) then
     ! don`t read
     ib_first = -1
     ib_last  = -1
     ib_size  =  0
   end if
   SAFE_ALLOCATE( ib_size_array, (peinf%npes) )
   ib_size_array = 0
   ib_size_array(peinf%inode+1) = ib_size

   !XXX Debug
   !XXX write(*,*) peinf%inode, iband_min, iband_max, tot_nbands, band_block, min_band_block, ib_first, ib_last, ib_size

#ifdef HDF5
#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, mpierr)
    call h5fopen_f('WFN_co.h5', H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id,error)
    call MPI_AllReduce(MPI_IN_PLACE, ib_size_array, peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
    call h5fopen_f('WFN_co.h5', H5F_ACC_RDONLY_F, file_id, error)
#endif
#endif

    ib_size_max = MAXVAL( ib_size_array )
    SAFE_ALLOCATE( my_cg_all, (ngktot, kp%nspin*kp%nspinor, MAX(ib_size,1)) )
    SAFE_ALLOCATE( ipe_cg_all, (ngktot, kp%nspin*kp%nspinor, ib_size_max) )
    my_cg_all = ZERO

#ifdef HDF5
    ! dummy allocation
    SAFE_ALLOCATE( peinf%does_it_ownv, (1, 1) )
    peinf%does_it_ownv = .true.
    call read_hdf5_bands_block(file_id, kp, ib_size_max, ib_size, &
                               peinf%does_it_ownv, ib_first, my_cg_all, is_bse = .true.)

    SAFE_DEALLOCATE_P(peinf%does_it_ownv)

    call h5fclose_f(file_id, error)
#endif

  end if

  if (xct%qflag/=0) then
    call progress_init(prog_info, 'reading wavefunctions (WFN_co)', 'state', &
      kp%nrk*(xct%nvb_co+xct%ncb_co))
  else
    call progress_init(prog_info, 'reading wavefunctions (WFN_co)', 'state', &
      kp%nrk*(xct%ncb_co))
  endif
 
  npes_hdf5 = 0
  if ( xct%use_wfn_hdf5 ) npes_hdf5 = peinf%npes-1

do ipe = 0, npes_hdf5

  if ( xct%use_wfn_hdf5 ) then
    if ( ib_size_array(ipe+1) > 0) then
     ipe_cg_all = ZERO
     if ( ipe == peinf%inode ) then
       ipe_cg_all(:, :, 1:ib_size) = my_cg_all(:, :, 1:ib_size)
     end if
#ifdef MPI
     call MPI_Bcast(ipe_cg_all, ngktot * kp%nspin*kp%nspinor * ib_size_max, &
                    MPI_SCALAR, ipe, MPI_COMM_WORLD, mpierr)
#endif
    else
      ! no elements on this MPI task, cycle
      cycle
    end if
  end if

  istart = 1
  do irk=1,kp%nrk
    irks = 0
    do ii=1,kg%nr
      if (irk == indxk(ii)) then
        irks=ii
        exit
      endif
    enddo

    SAFE_ALLOCATE(gvec_kpt%components, (3, kp%ngk(irk)))

    if ( xct%use_wfn_hdf5 ) then
      gvec_kpt%components(:,:) = gvec_kpt_components_all(1:3, istart:istart+kp%ngk(irk)-1)
    else
      call read_binary_gvectors(25, kp%ngk(irk), kp%ngk(irk), gvec_kpt%components)
    end if

    SAFE_ALLOCATE(cg, (kp%ngk(irk),kp%nspin*kp%nspinor))
    if(irks > 0) then
      do ii = 1, kp%ngk(irk)
        call findvector(wfnv%isort(ii), gvec_kpt%components(:, ii), gvec)
        if(wfnv%isort(ii) == 0) call die('could not find gvec', only_root_writes=.true.)
      enddo

      wfnv%ng=kp%ngk(irk)
      wfnc%ng=kp%ngk(irk)
      wfnc%isort(1:gvec%ng)=wfnv%isort(1:gvec%ng)
      ickmem=0
    endif

!        write(6,*) peinf%inode, 'loop wfnup', irk

! Loop Over Bands

    ii_start = 1
    ii_end   = kp%mnband
    if ( xct%use_wfn_hdf5 ) then
      ii_start = iband_min + band_block * ipe
      ii_end   = min(ii_start+ band_block - 1, iband_max)
    end if
    ii_loc = 0
    do ii = ii_start, ii_end
      ii_loc = ii_loc + 1

      ! FHJ: do we care about the band (any spin)?
      need_band = .false.
      if (irks>0) then
        do is = 1, kp%nspin
          ! If ii is one of the selected valence band...
          ! (Skip current valence band if valence bands will be read from WFNq_co)
          if ((ii<=kp%ifmax(irk,is)).and. &
              (ii>=(kp%ifmax(irk,is)-xct%nvb_co+1)).and.&
              xct%qflag/=0) then
            need_band = .true.
          endif
          ! If ii is one of the selected conduction band...
          if ((ii>=(kp%ifmax(irk,is)+1)) .and. &
              (ii<=(kp%ifmax(irk,is)+xct%ncb_co))) then
            need_band = .true.
          endif
        enddo
      endif
      if (need_band) call progress_step(prog_info)

      if ( xct%use_wfn_hdf5 ) then
        cg = ZERO 
        cg (:,:) = ipe_cg_all(istart:istart+kp%ngk(irk)-1, 1:kp%nspin*kp%nspinor, ii_loc)
      else
        call read_binary_data(25, kp%ngk(irk), kp%ngk(irk), kp%nspin*kp%nspinor, &
          cg, bcast=.true., dont_read=.not.need_band)
      end if

      ! If we do not need this band, skip it...
      if (irks == 0) cycle

      do is = 1, kp%nspin
        ! If ii is one of the selected valence band...
        ! (Skip current valence band if valence bands will be read from WFNq_co)
        if ((ii<=kp%ifmax(irk,is)).and. &
            (ii>=(kp%ifmax(irk,is)-xct%nvb_co+1)).and.&
            xct%qflag/=0) then

          ivband = kp%ifmax(irk,is)-ii+1
          iwritev = peinf%ipev(peinf%inode+1,ivband,irk)
          if (xct%qflag==1) then
            iwritek = peinf%ipek(peinf%inode+1,irk)
          else if (xct%qflag==2) then
            iwritek = peinf%ipekq(peinf%inode+1,irk)
          endif

          if (iwritev/=0) then
            call checknorm('WFN_co',ii,irks,kp%ngk(irk),kp%nspin,kp%nspinor,cg(:,:))
            do isp = 1, kp%nspinor
              intwfnv%cg(1:wfnv%ng,iwritev,is*isp) = cg(1:kp%ngk(irk), is*isp)
            enddo
          endif
          if (iwritek/=0) then
            intwfnv%isort(1:gvec%ng,iwritek) = wfnv%isort(1:gvec%ng)
            intwfnv%ng(iwritek) = kp%ngk(irk)
          endif
        endif !ii is one of the selected valence band

        ! If ii is one of the selected conduction band...
        if ((ii>=(kp%ifmax(irk,is)+1)) .and. &
            (ii<=(kp%ifmax(irk,is)+xct%ncb_co))) then

          icband = ii-kp%ifmax(irk,is)
          iwritec = peinf%ipec(peinf%inode+1,icband,irk)
          iwritek = peinf%ipek(peinf%inode+1,irk)

          if (iwritec/=0) then
            call checknorm('WFN_co',ii,irks,kp%ngk(irk),kp%nspin,kp%nspinor,cg(:,:))
            do isp = 1, kp%nspinor
              intwfnc%cg(1:wfnc%ng,iwritec,is*isp) = cg(1:kp%ngk(irk), is*isp)
            enddo
          endif
          if (iwritek/=0) then
            intwfnc%isort(1:gvec%ng,iwritek) = wfnc%isort(1:gvec%ng)
            intwfnc%ng(iwritek) = kp%ngk(irk)
          endif
        endif ! ii is one of the selected conduction bands
      end do ! is
      if (ii>maxval(kp%ifmax)+xct%ncb_co .and. irk==kp%nrk) then
        exit
      endif
    enddo ! ii (loop on bands)

    SAFE_DEALLOCATE(cg)
   
    istart = istart + kp%ngk(irk) 
  enddo !end loop over k-points
end do ! ipe loop for HDF5 case
  call progress_free(prog_info)

  if ( xct%use_wfn_hdf5  ) then
    SAFE_DEALLOCATE( gvec_kpt_components_all )
    SAFE_DEALLOCATE( ib_size_array )
    SAFE_DEALLOCATE( my_cg_all )
    SAFE_DEALLOCATE( ipe_cg_all )
  end if
! Write out info about xtal

  if (peinf%inode.eq.0) then
    write(6,'(/1x,a)') 'Crystal wavefunctions read from file WFN_co:'
    write(6,'(1x,a,i0)') '- Number of k-points in WFN_co: ', kp%nrk
    write(6,'(1x,a,i0)') '- Number of k-points in the full BZ of WFN_co: ', kg%nf
    if (peinf%verb_high) then
      write(6,'(1x,a)') '- K-points:'
      write(6,'(1(2x,3(1x,f10.6)))') kg%r(1:3,1:kg%nr)
    endif
    if ( .not. xct%use_wfn_hdf5 ) call close_file(25)
  endif ! node 0

!-----------------------------------------------------------------
! DYQ: If doing a finite Q calculation with WFNq_co read and distribute valence wave functions here
!
  if (xct%qflag.eq.0) then

    if ( xct%use_wfn_hdf5  ) then
      if (peinf%inode==0) write(6,*)
      if (peinf%inode==0) write(6,'(a)') ' Reading HDF5 wavefuntion (WFNq_co.h5)'
      ngktot = SUM(kpq%ngk)
      SAFE_ALLOCATE( gvec_kpt_components_all, (3,ngktot) )
#ifdef HDF5
      call read_hdf5_wfn_gvectors('WFNq_co.h5', gvec_kpt_components_all, ngktot)
#endif

      ! We only need the valence states
      iband_min      = MINVAL(kpq%ifmax(:,:)) - xct%nvb_co + 1
      iband_min      = MAX(iband_min, 1)
      iband_max      = MAXVAL(kpq%ifmax(:,:))
      iband_max      = MIN(iband_max, kpq%mnband)
      tot_nbands     = iband_max - iband_min + 1
      band_block     = (tot_nbands + peinf%npes - 1) / peinf%npes
      ! read at least 128 Mb per MPI task
      min_band_block = 128.0D+00 / ( dble(ngktot) * dble(kpq%nspin*kpq%nspinor) * 16.0D+00 /1024.0D+00/1024D+00 )
      min_band_block = MAX(min_band_block, 1)
      if ( xct%wfn_hdf5_min_band_block > 0 )  min_band_block = xct%wfn_hdf5_min_band_block
      if ( min_band_block > band_block ) then
        band_block = min_band_block
      end if

      ib_first = iband_min + band_block * peinf%inode
      ib_last  = min(ib_first + band_block - 1, iband_max)
      ib_size  = ib_last - ib_first + 1
      if ( ib_size < 1 ) then
        ! don`t read
        ib_first = -1
        ib_last  = -1
        ib_size  =  0
      end if
      SAFE_ALLOCATE( ib_size_array, (peinf%npes) )
      ib_size_array = 0
      ib_size_array(peinf%inode+1) = ib_size

      !XXX Debug
      !XXX write(*,*) peinf%inode, iband_min, iband_max, tot_nbands, band_block, min_band_block, ib_first, ib_last, ib_size
#ifdef HDF5
#ifdef MPI
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, mpierr)
      call h5fopen_f('WFNq_co.h5', H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
      call h5pclose_f(plist_id,error)
      call MPI_AllReduce(MPI_IN_PLACE, ib_size_array, peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
      call h5fopen_f('WFNq_co.h5', H5F_ACC_RDONLY_F, file_id, error)
#endif
#endif

      ib_size_max = MAXVAL( ib_size_array )
      SAFE_ALLOCATE( my_cg_all, (ngktot, kpq%nspin*kpq%nspinor, MAX(ib_size,1)) )
      SAFE_ALLOCATE( ipe_cg_all, (ngktot, kpq%nspin*kpq%nspinor, ib_size_max) )
      my_cg_all = ZERO

#ifdef HDF5
      ! dummy allocation
      SAFE_ALLOCATE( peinf%does_it_ownv, (1, 1) )
      peinf%does_it_ownv = .true.
      call read_hdf5_bands_block(file_id, kpq, ib_size_max, ib_size, &
                                 peinf%does_it_ownv, ib_first, my_cg_all, is_bse = .true.)

      SAFE_DEALLOCATE_P(peinf%does_it_ownv)

      call h5fclose_f(file_id, error)
#endif

    end if

    call progress_init(prog_info, 'reading Q-shifted valence wavefunctions (WFNq_co)', 'state', &
      kpq%nrk*(xct%nvb_co))

  npes_hdf5 = 0
  if ( xct%use_wfn_hdf5 ) npes_hdf5 = peinf%npes-1

  ! loop over process only in case of wfn_hdf5
  do ipe = 0, npes_hdf5

    if ( xct%use_wfn_hdf5 ) then
      if ( ib_size_array(ipe+1) > 0) then
       ipe_cg_all = ZERO
       if ( ipe == peinf%inode ) then
         ipe_cg_all(:, :, 1:ib_size) = my_cg_all(:, :, 1:ib_size)
       end if
#ifdef MPI
       call MPI_Bcast(ipe_cg_all, ngktot * kpq%nspin*kpq%nspinor * ib_size_max, &
                      MPI_SCALAR, ipe, MPI_COMM_WORLD, mpierr)
#endif
      else
        ! no elements on this MPI task, cycle
        cycle
      end if
    end if

    istart = 1
    do irk=1,kpq%nrk
      irks = 0
      do ii=1,kgq%nr
        if (irk == indxkq(ii)) then
          irks=ii
          exit
        endif
      enddo

      SAFE_ALLOCATE(gvec_kpt%components, (3, kpq%ngk(irk)))
      if ( xct%use_wfn_hdf5 ) then
        gvec_kpt%components(:,:) = gvec_kpt_components_all(1:3, istart:istart+kpq%ngk(irk)-1)
      else
        call read_binary_gvectors(26, kpq%ngk(irk), kpq%ngk(irk), gvec_kpt%components)
      end if

      SAFE_ALLOCATE(cg, (kpq%ngk(irk),kpq%nspin*kpq%nspinor))
      if(irks > 0) then
        do ii = 1, kpq%ngk(irk)
          call findvector(wfnv%isort(ii), gvec_kpt%components(:, ii), gvec)
          if(wfnv%isort(ii) == 0) call die('could not find gvec', only_root_writes=.true.)
        enddo

        wfnv%ng=kpq%ngk(irk)
        ickmem=0
      endif
      
      ! Loop over bands
      ii_start = 1
      ii_end   = kpq%mnband
      if ( xct%use_wfn_hdf5 ) then
        ii_start = iband_min + band_block * ipe
        ii_end   = min(ii_start + band_block - 1, iband_max)
      end if
      ii_loc = 0
      do ii = ii_start, ii_end
        ii_loc = ii_loc + 1

        ! FHJ: do we care about the data?
        need_band = .false.
        if (irks>0) then
          do is = 1, kpq%nspin
            if ((ii>=kpq%ifmax(irk,is)-xct%nvb_co+1) .and. (ii<=kpq%ifmax(irk,is))) then
              need_band = .true.
            endif
          enddo
        endif
        if (need_band) call progress_step(prog_info)

        if ( xct%use_wfn_hdf5 ) then
          cg = ZERO
          cg (:,:) = ipe_cg_all(istart:istart+kpq%ngk(irk)-1, 1:kpq%nspin*kpq%nspinor, ii_loc)
        else
          call read_binary_data(26, kpq%ngk(irk), kpq%ngk(irk), kpq%nspin*kpq%nspinor, &
            cg, bcast=.true., dont_read=.not.need_band)
        end if

        ! If we do not need this band, skip it...
        if(irks == 0) cycle

        do is=1, kpq%nspin
          ! if current band is one of the selected valence bands
          if ((ii>=kpq%ifmax(irk,is)-xct%nvb_co+1) .and. (ii<=kpq%ifmax(irk,is))) then
            ivband = kpq%ifmax(irk,is)-ii+1
            iwritev = peinf%ipev(peinf%inode+1,ivband,irk)
            iwritek = peinf%ipekq(peinf%inode+1,irk)

            if (iwritev/=0) then
              call checknorm('WFNq_co',ii,irks,kpq%ngk(irk),kpq%nspin,kpq%nspinor,cg(:,:))
              do isp = 1, kp%nspinor
                intwfnv%cg(1:wfnv%ng,iwritev,is*isp) = cg(1:kpq%ngk(irk), is*isp)
              enddo
            endif
            if (iwritek/=0) then
              intwfnv%isort(1:gvec%ng,iwritek) = wfnv%isort(1:gvec%ng)
              intwfnv%ng(iwritek) = kpq%ngk(irk)
            endif
          endif !ii is one of the selected valence bands
        enddo ! loop over nspin
        if (ii>maxval(kpq%ifmax) .and. irk==kpq%nrk) then
          exit
        endif
      enddo ! loop over all bands

      SAFE_DEALLOCATE(cg)

      istart = istart + kpq%ngk(irk)
    enddo !end loop over k-points
  end do ! ipe loop for HDF5 case
    call progress_free(prog_info)

    if ( xct%use_wfn_hdf5  ) then
      SAFE_DEALLOCATE( gvec_kpt_components_all )
      SAFE_DEALLOCATE( ib_size_array )
      SAFE_DEALLOCATE( my_cg_all )
      SAFE_DEALLOCATE( ipe_cg_all )
    end if

    if (peinf%inode.eq.0) then
      write(6,'(/1x,a)') 'Crystal wavefunctions read from file WFNq_co:'
      write(6,'(1x,a,i0)') '- Number of k-points in WFNq_co: ', kgq%nr
      write(6,'(1x,a,i0)') '- Number of k-points in the full BZ of WFNq_co: ', kgq%nf
      if (peinf%verb_high) then
        write(6,'(1x,a)') '- K-points:'
        write(6,'(1(2x,3(1x,f10.6)))') kgq%r(1:3,1:kgq%nr)
      endif
      if ( .not. xct%use_wfn_hdf5 ) call close_file(26)
    endif ! node 0

  endif ! Finished distributing WFNq_co

  SAFE_DEALLOCATE_P(wfnv%isort)
  SAFE_DEALLOCATE_P(wfnc%isort)
  if (xct%qflag.eq.0) then
    SAFE_DEALLOCATE_P(gvecq%components)
  endif
  SAFE_DEALLOCATE_P(gvec_kpt%components)
  SAFE_DEALLOCATE(indxk)
  SAFE_DEALLOCATE(indxkq)
      
  POP_SUB(input_kernel)

  return

contains

  ! FHJ: Figure out xct%ecute and xct%ecutg, if the user didn`t specify them.
  subroutine get_ecut()
#ifdef HDF5
    integer :: ng, nq, nfreq, nfreq_imag, nmtx_max, qgrid(3), freq_dep
#endif
    real(DP) :: ecuts
    !real(DP) :: raux

    PUSH_SUB(input_kernel.get_ecut)

    if (xct%theory/=0) then
      xct%ecute = 0
      if (xct%ecutg < TOL_ZERO) xct%ecutg = kp%ecutwfc
      return
    endif
    if (peinf%inode==0) then
#ifdef HDF5
      if(xct%use_hdf5) then
        call read_eps_grid_sizes_hdf5(ng, nq, ecuts, nfreq, &
          nfreq_imag, nmtx_max, qgrid, freq_dep, 'eps0mat.h5')
      else
#endif
        call open_file(10, file='eps0mat', form='unformatted', status='old')
        read(10)
        read(10)
        read(10)
        read(10)
        read(10)
        read(10)
        read(10) ecuts
        call close_file(10)
      endif
#ifdef HDF5
    endif
#endif
#ifdef MPI
    call MPI_Bcast(ecuts, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
#endif
    if (xct%ecute < TOL_ZERO) xct%ecute = ecuts
    !raux = maxval(crys%bdot)
    !if (xct%ecutg < TOL_ZERO) xct%ecutg = xct%ecute + raux
    if (xct%ecutg < TOL_ZERO) xct%ecutg = kp%ecutwfc

    POP_SUB(input_kernel.get_ecut)

  end subroutine get_ecut

end subroutine input_kernel

end module input_kernel_m
