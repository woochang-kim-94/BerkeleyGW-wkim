!==========================================================================
!
! Routines:
!
! (1) genwf_mpi()  Originally By (JRD)             Last Modified 11/2009 (JRD)
!
!     On entry:
!     qq  = current q-vector
!     rk  = current k point in irr. zone
!     irk = its index
!
!     On exit:
!     vwfn%ev and vwfn%zv hold eigenvalues and wavefunctions (valence)
!     cwfn%ec and cwfn%zc hold eigenvalues and wavefunctions (conduction)
!
!     with proper phases and ordering for the k-point rk (given the
!     data for irr. zone)
!
!     subroutine generates valence-band wavefunctions for rk(irk)
!     and conduction-band wavefunctions for rk(irk) from the
!     wavefunctions available for rk in the irreducible wedge of the
!     bz
!
!     i   rk                 rk(irk)
!     o c ngv,...,ev       valence-band wavefunctions for rk+q
!     and associated data
!     o c ngc,...,ec       conduction-band wavefunctions for rk
!     and associated data
!
!==========================================================================

#include "f_defs.h"

module genwf_mpi_m

  use global_m
  use find_kpt_match_m
  use gmap_m
  use susymmetries_m
  use input_utils_m
  use misc_m
  use sort_m
  use timing_m, only: timing => epsilon_timing

  implicit none

  private 
    integer, allocatable :: GK(:)
    logical :: GK_allocated=.false.
  public :: genwf_mpi

contains  

subroutine genwf_mpi(syms,gvec,crys,kp,kpq,irk,rk,qq,vwfn,pol,cwfn,use_wfnq,intwfnv,intwfnvq,intwfnc,ivin)
  type (symmetry), intent(in) :: syms
  type (gspace), intent(in) :: gvec
  type (crystal), intent(in) :: crys
  type (kpoints), target, intent(in) :: kp
  type (kpoints), target, intent(in) :: kpq
  integer, intent(in) :: irk
  real(DP), intent(in) :: rk(3)
  real(DP), intent(in) :: qq(3)
  type (valence_wfns), intent(inout) :: vwfn
  type (polarizability), intent(in) :: pol
  type (conduction_wfns), intent(inout) :: cwfn
  logical, intent(in) :: use_wfnq
  type (int_wavefunction), intent(in) :: intwfnv
  type (int_wavefunction), intent(in) :: intwfnvq
  type (int_wavefunction), intent(in) :: intwfnc
  integer, intent(in) :: ivin

  integer :: itval,jband
  integer :: ng
  integer, allocatable :: isortc(:)
  real(DP), allocatable :: eig(:,:)
  SCALAR, allocatable :: zin(:,:),zinc(:,:)
  SCALAR, allocatable :: zintemp(:,:)
  type(kpoints), pointer :: kp_point

  character :: tmpstr*120
  integer :: itqq,i
  integer :: n,ig,ispin,iband
  integer :: naddc
  integer :: ikrkq,kgq(3),kgqq(3)
  integer, allocatable :: ind(:),isorti(:)
  real(DP), allocatable :: xnorm(:)
  real(DP) :: qk(3),rkq(3)
  real(DP), allocatable :: ekin(:)
  real(DP) :: rkmatch(3)
#ifdef CPLX
  SCALAR, allocatable :: zvtemp(:,:),zctemp(:,:)
  complex(DP) :: umtrx(2,2)
#endif

  SCALAR, allocatable :: ph(:)

  SCALAR, allocatable, save :: ph_old(:)

  integer, allocatable, save :: ind_old(:)
  integer, save :: irk_old=0

  PUSH_SUB(genwf_mpi)

  if (.not.GK_allocated) then
    SAFE_ALLOCATE(GK, (gvec%ng))
    call get_GK_array_from_gvecs(gvec%ng, gvec%components, GK)
    GK_allocated = .true.
  endif
  SAFE_ALLOCATE(isortc, (gvec%ng))
  SAFE_ALLOCATE(eig, (cwfn%nband,kp%nspin))
  SAFE_ALLOCATE(xnorm, (kp%nspin))

  call timing%start(timing%genwf_val)

! rkq = rk + qq   (i.e. rkq = current kpoint in irr. zone + qq)

  rkq(1:3) = rk(1:3) + qq(1:3)

! Should we look for the valence WFNs in WFNq or WFN?
  if (use_wfnq) then
    kp_point => kpq
  else
    kp_point => kp
  endif

! We used to assume WFNq contained the needed q->0 point,
! but it is better to unfold the shifted grid for more flexibility.

  call find_kpt_match(kp_point, syms, rkq, ikrkq, itqq, kgqq)

  if(ikrkq == 0) then
    if(peinf%inode == 0) write(0,'(a,3f12.6,/,a,3f12.6)') 'rk = ', rk(:), 'qq = ', qq(:)
    write(tmpstr,'(a,3f8.3,a)') 'genwf_mpi: No match for rkq point:',rkq,' in file '
    if(use_wfnq) then
      write(tmpstr,'(a,a)') TRUNC(tmpstr), ' WFNq'
    else
      write(tmpstr,'(a,a)') TRUNC(tmpstr), ' WFN'
    endif
    call die(tmpstr, only_root_writes = .true.)
  endif

  rkmatch(1:3) = kp_point%rk(1:3, ikrkq)
  vwfn%idx_kp = ikrkq

!      if(peinf%inode.eq.0) then
  itval=vwfn%nband+pol%ncrit
  if (.not.use_wfnq) then
    ng = intwfnv%ng(ikrkq)
    if (irk .ne. irk_old) then
      eig(1:itval,1:kp%nspin)=kp_point%el(1+vwfn%ncore_excl:itval+vwfn%ncore_excl,ikrkq,1:kp%nspin)
      isortc(1:ng)=intwfnv%isort(1:ng,ikrkq)
    endif
    qk(:)=intwfnv%qk(:,ikrkq)
  else
    ng = intwfnvq%ng(ikrkq)
    if (irk .ne. irk_old) then
      eig(1:itval,1:kp%nspin)=kp_point%el(1+vwfn%ncore_excl:itval+vwfn%ncore_excl,ikrkq,1:kp%nspin)
      isortc(1:ng)=intwfnvq%isort(1:ng,ikrkq)
    endif
    qk(:)=intwfnvq%qk(:,ikrkq)
  endif
!      endif

  SAFE_ALLOCATE(zintemp, (ng,kp%nspin*kp%nspinor))
  SAFE_ALLOCATE(zin, (ng,kp%nspin*kp%nspinor))
  
  jband = (ikrkq-1)*peinf%nvownactual + ivin
  zintemp=0D0
  if (use_wfnq) then
    zintemp(1:ng,1:kp%nspin*kp%nspinor)=intwfnvq%cg(1:ng,jband,1:kp%nspin*kp%nspinor)
  else
    zintemp(1:ng,1:kp%nspin*kp%nspinor)=intwfnv%cg(1:ng,jband,1:kp%nspin*kp%nspinor)
  endif

  zin(:,:)=zintemp(:,:) 
  SAFE_DEALLOCATE(zintemp)

! Check kpoint

  if(any(abs(rkmatch(1:3) - qk(1:3)) .gt. TOL_Small)) call die('genwf_mpi: rkmatch')

  if (irk .ne. irk_old) then

! Compute kinetic energies for rkq+g

    SAFE_ALLOCATE(ekin, (gvec%ng))
    call timing%start(timing%genwf_ekin)
    call kinetic_energies(gvec, crys%bdot, ekin, qvec = rkq)
    call timing%stop(timing%genwf_ekin)

! sort array ekin to ascending order, store indices in array vwfn%isort
! WARNING: one should not assume that in the case of
! q-->0 the orders as provided below and as read in from
! WFNq file is the same (sorting may change....)
! ===> need to call gmap also for q-->0 (see below)
! EKC: We initialize vwfn%isort to the appropriate array
! before reading in. this way we do not get zeroes in the array
! these are the valence wave-functions that do not need
! to be changed

    call timing%start(timing%genwf_sort)
    call sort_threaded(gvec%ng, ekin, vwfn%isort, GK)
    !call sortrx(gvec%ng,ekin,vwfn%isort,gvec=gvec%components)
    call timing%stop(timing%genwf_sort)

    SAFE_ALLOCATE(isorti, (gvec%ng))
    do i=1,gvec%ng
      isorti(vwfn%isort(i))=i
    enddo
    do i=1,ng
      isorti(isortc(i))=i
    enddo

    vwfn%ngv=ng

! SIB: put read eigenvalues into vwfn%ev(band,spin).

    SAFE_ALLOCATE(vwfn%ev, ((vwfn%nband+pol%ncrit),kp%nspin))
    vwfn%ev(1:(vwfn%nband+pol%ncrit),:) = eig(1:(vwfn%nband+pol%ncrit),:)

! JRD: Map planewave components for rk+q, to those of rk
! (even for q--> 0)
!
! SIB: get phases (ph) and indices (ind) for g-vectors
! gvec%components(:,vwfn%isort(1:vwfn%ngv))+kgqq

    SAFE_ALLOCATE(ind, (vwfn%ngv))
    SAFE_ALLOCATE(ph, (vwfn%ngv))
    call gmap(gvec,syms,vwfn%ngv,itqq,kgqq,vwfn%isort,isorti,ind,ph,.true.)
    if (irk_old .ne. 0) then
      SAFE_DEALLOCATE(ph_old)
      SAFE_DEALLOCATE(ind_old)
    endif
    SAFE_ALLOCATE(ph_old, (vwfn%ngv))
    SAFE_ALLOCATE(ind_old, (vwfn%ngv))
    ph_old(:) = ph(:)
    ind_old(:) = ind(:)
  else

    SAFE_ALLOCATE(ph, (vwfn%ngv))
    SAFE_ALLOCATE(ind, (vwfn%ngv))
    ph(:)=ph_old(:)
    ind(:)=ind_old(:)

  endif  ! irk = irk_old
  SAFE_DEALLOCATE(eig)

  SAFE_ALLOCATE(vwfn%zv, (vwfn%ngv,kp%nspin*kp%nspinor))

! XAV: vwfn%zv(ig) corresponds really to the
! vwfn%isort(ig) G-vector (between 1 and ng)
! The subroutine gmap assumes that, as read from WFNq or WFN,
! zin(ig) corresponds really to isortc(ig) G-vector !!!
!
! SIB:  checks that zin vectors have norm greater than 1e-6, and then
! normalizes them to have unit square modulus.

  do ispin=1,kp%nspin*kp%nspinor
    vwfn%zv(:,ispin)=zin(ind(:),ispin)*ph(:)
  enddo

! In spinor case, we must rotate spinors according to spinor rotation matrix umtrx
#ifdef CPLX
  if (kp%nspinor.eq.2) then
    SAFE_ALLOCATE(zvtemp, (vwfn%ngv,kp%nspin*kp%nspinor))
    call susymmetries(crys%bvec(:,:),syms%mtrx(:,:,itqq),umtrx,itqq)
    zvtemp = matmul(vwfn%zv,umtrx)
    vwfn%zv = zvtemp
    SAFE_DEALLOCATE(zvtemp)
  endif
#endif

! BAB:  we check if the norm differs appreciably from unity.
! there is no longer a need to further normalize the vector
  if (peinf%check_norms) then
    do ispin=1,kp%nspin
      call compute_norm(xnorm(ispin),ispin,vwfn%ngv,kp%nspinor,vwfn%zv(:,:))
      if(abs(xnorm(ispin) - 1d0) > TOL_Small) then
        write(0,*) 'Bad norm: sym op=',itqq,'ik=',ikrkq,'ispin=',ispin,'norm=',xnorm(ispin)
        call die('genwf_mpi: Bad norm')
      endif
    enddo
  endif

! End calculation of valence-band wavefunctions

  SAFE_DEALLOCATE(zin)
  SAFE_DEALLOCATE(ind)
  SAFE_DEALLOCATE(ph)
  SAFE_DEALLOCATE(xnorm)
  call timing%stop(timing%genwf_val)

  !FHJ: only generate conduction WFNs if ivin=1
  if(irk .ne. irk_old) then

!------------------------------------------------------------------------
! Generate conduction-band wavefunctions for rk
! find rk, r, g0 such that rk=r(rk)+g0

! SIB: This seems redundant, but find a k-point and symmetry so that
! rk = sym%mtrx(:,:,itqq)*kp%rk(irk_,:) + kgq, where kgq is integer 3-vec

    call timing%start(timing%genwf_cond)

    call find_kpt_match(kp, syms, rk, ikrkq, itqq, kgq)
    cwfn%idx_kp = ikrkq

    if(ikrkq == 0) call die('genwf_mpi: kgq mismatch')

! Write out rk, it and kgq

    if (peinf%verb_debug .and. peinf%inode==0) then
      write(6,7000) (rk(i),i=1,3),ikrkq,(kp%rk(i,ikrkq),i=1,3),itqq,(kgq(i),i=1,3)
7000  format(1x,'rk=',3f7.3,1x,'irk_=',i5,1x,' rk=',3f7.3,1x,'it=',i5,1x,'kg0=',3i3)
    endif

! SIB: if we already found this k-point last time, get its qk, ng,
! and isortc(:).  Otherwise, skip ikrkq-1 records, and read in information (qk,cwfn%ec,ng,isortc),

    SAFE_ALLOCATE(cwfn%ec, (cwfn%nband,kp%nspin))
    
    ng = intwfnc%ng(ikrkq)
    cwfn%ec(1:cwfn%nband,1:kp%nspin)=kp%el(1+vwfn%ncore_excl:cwfn%nband+vwfn%ncore_excl,ikrkq,1:kp%nspin)
    qk(:)=intwfnc%qk(:,ikrkq)
    isortc(1:ng)=intwfnc%isort(1:ng,ikrkq)
    
! Check kpoint (again ...  boring...)
! Check that kp%rk(:,ikrkq) = qk  (why it wouldn`t is a mystery!)

    if(any(abs(kp%rk(1:3, ikrkq) - qk(1:3)) .gt. TOL_Small)) &
      call die('genwf_mpi: qk mismatch')

    cwfn%ngc=ng

! Compute inverse to isort
! NOTE: isortc orders   |kp%rk+G|^2
! It is not necessarily the same order as |rk+G|^2
! (particularly if umklapp, ie kgq non zero)

    call timing%start(timing%genwf_ekin)
    call kinetic_energies(gvec, crys%bdot, ekin, qvec = rk)
    call timing%stop(timing%genwf_ekin)

! Sort array ekin to ascending order
! store indices in array isort

    call timing%start(timing%genwf_sort)
    call sort_threaded(gvec%ng, ekin, cwfn%isort, GK)
    !call sortrx(gvec%ng,ekin,cwfn%isort,gvec=gvec%components)
    call timing%stop(timing%genwf_sort)

    do i=1,ng
      isorti(isortc(i))=i
    enddo
    SAFE_DEALLOCATE(ekin)

! map planewave components for rk to those of rk
! compute phases
! We do not the isorti related to kp%rk BUT the cwfn%isort related to rk

    SAFE_ALLOCATE(ind, (cwfn%ngc))
    SAFE_ALLOCATE(ph, (cwfn%ngc))
    call gmap(gvec,syms,cwfn%ngc,itqq,kgq,cwfn%isort,isorti,ind,ph,.true.)
    SAFE_DEALLOCATE(isorti)

! generate conduction-band wavefunctions
! loop over wavefunctions
! read conduction-band from file one by one

    SAFE_ALLOCATE(cwfn%zc, (peinf%ncownactual*cwfn%ngc,kp%nspin*kp%nspinor))
    SAFE_ALLOCATE(zinc, (cwfn%ngc,kp%nspin*kp%nspinor))

    do n=1,peinf%ncownactual
      
      jband= (ikrkq-1)*peinf%ncownactual + n
      
      iband = intwfnc%cbi(jband)
      zinc(1:ng,1:kp%nspin*kp%nspinor)=intwfnc%cg(1:ng,jband,1:kp%nspin*kp%nspinor)
      
      if(iband .ne. peinf%invindexc(n)) call die('genwf_mpi: invindexc mismatch')

      naddc=(n-1)*cwfn%ngc

!
! Loop over components of wfns
! note that only conduction-band wfns are stored-
! they start in the 1st position with the state nvband+1
!

      do ig=1,cwfn%ngc
        cwfn%zc(naddc+ig,:)=ph(ig)*zinc(ind(ig),:)
      enddo

! In spinor case, we must rotate spinors according to spinor rotation matrix umtrx

#ifdef CPLX
      if (kp%nspinor.eq.2) then
        SAFE_ALLOCATE(zctemp, (cwfn%ngc,kp%nspin*kp%nspinor))
        call susymmetries(crys%bvec(:,:),syms%mtrx(:,:,itqq),umtrx,itqq)
        zctemp = matmul(cwfn%zc(naddc+1:naddc+cwfn%ngc,:),umtrx)
        cwfn%zc(naddc+1:naddc+cwfn%ngc,:) = zctemp
        SAFE_DEALLOCATE(zctemp)
      endif
#endif

    enddo  ! n (cond-bands per node [ncownactual] loop)

    SAFE_DEALLOCATE(zinc)
    SAFE_DEALLOCATE(ind)
    SAFE_DEALLOCATE(ph)

!     end generation of conduction-band wavefunctions for rk

    call timing%stop(timing%genwf_cond)

  endif ! (irk .ne. irk_old)

  irk_old = irk
  SAFE_DEALLOCATE(isortc)

  POP_SUB(genwf_mpi)

  return
end subroutine genwf_mpi

end module genwf_mpi_m
