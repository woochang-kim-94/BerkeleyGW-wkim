#include "f_defs.h"

module genwf_kernel_m

  use blas_m
  use global_m
  use gmap_m
  use input_utils_m
  use misc_m
  use sort_m
  use susymmetries_m
  implicit none

  private

  public :: &
    genwf_kernel

contains

!-----------------------------------------------------------------------
subroutine genwf_kernel(crys,gvec,kg,kgq,syms, &
  wfnc,wfnv,nspin,ik,ic,iv,indexq,xct,intwfnv,intwfnc)
!
!     input: crys, gvec, kg,  syms types
!            nspin  (number of spins)
!            ik     (index of k-point in full BZ)
!            ib     (index of ck block)
!     output: wfnc   (conduction states at point ik, block ib)
!             wfnv   (valence states at point ik)
!
!-----------------------------------------------------------------------
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (grid), intent(in) :: kg
  type (grid), intent(in) :: kgq
  type (symmetry), intent(in) :: syms
  type (wavefunction), intent(out) :: wfnc,wfnv
  integer, intent(in) :: nspin,ik,ic,iv
  integer, intent(in) :: indexq(kg%nf)
  type (int_wavefunction), intent(in) :: intwfnc,intwfnv
  type (xctinfo), intent(in) :: xct
      
  character :: errmsg*100
  integer :: cnb,cng,cns,cnsp,vng,vns,vnsp
  integer :: ig,ii,jj,kk,info
  integer :: ivown,icown,ikown,invband,incband,ijk
  real(DP) :: xnorm,qk(3)
  
  integer, save :: ikold=0
  integer, save :: ikold2=0
  integer, save :: vngold=0
  integer, save :: cngold=0
  integer, save :: ifirst=0
  
  integer, allocatable :: isorti(:)
  integer, allocatable :: ind(:),isort(:)
  integer, allocatable :: indq(:),isortq(:)
  integer, allocatable, save :: ind_old(:),isort_old(:)
  integer, allocatable, save :: indq_old(:),isortq_old(:)
  integer, allocatable, save :: ind_old2(:),isort_old2(:)
  integer, allocatable, save :: indq_old2(:),isortq_old2(:)
  real(DP), allocatable :: ekin(:)
  SCALAR, allocatable :: ccg(:,:,:),vcg(:,:,:),ph(:),phq(:)
  SCALAR, allocatable, save :: ph_old(:),phq_old(:)
  SCALAR, allocatable, save :: ph_old2(:),phq_old2(:)
#ifdef CPLX
    SCALAR, allocatable :: zvtemp(:,:), zctemp(:,:)
    complex(DP) :: umtrx(2,2)
#endif

  PUSH_SUB(genwf_kernel)

!-----------------------------------------------------------------------
! Deal with the valence wavefunctions

!      write(6,*) peinf%inode, 'in genwf',iv,ic,ik

  if (xct%qflag .eq. 0) then
    ivown = peinf%ipev(peinf%inode+1,iv,kgq%indr(indexq(ik)))
    ikown = peinf%ipekq(peinf%inode+1,kgq%indr(indexq(ik)))
  elseif (xct%qflag .eq. 2) then
    ! DYQ: For finite center-of-mass Q, the valence wavefunctions
    ! will be at k+Q, and the conduction band will be at k
    ivown = peinf%ipev(peinf%inode+1,iv,kg%indr(indexq(ik)))
    ikown = peinf%ipekq(peinf%inode+1,kg%indr(indexq(ik)))
#ifdef DEBUG
    if(ikown.eq.0) then
      write (6,*) '   ikown = ', ikown
      write (6,*) '   indexq(ik)= ', indexq(ik)
      write (6,*) '   ik  ', ik
      write (6,*) '   kg%indr(indexq(ik))  ', kg%indr(indexq(ik))
    endif
#endif
  else
    ivown = peinf%ipev(peinf%inode+1,iv,kg%indr(ik))
    ikown = peinf%ipek(peinf%inode+1,kg%indr(ik))
  endif
  
  vng = intwfnv%ng(ikown)
  vns = intwfnv%nspin
  vnsp = intwfnv%nspinor
  if (xct%ivpar .eq. 1) then
    SAFE_ALLOCATE(vcg, (vng,1,vns*vnsp))
  else
    SAFE_ALLOCATE(vcg, (vng,xct%nvb_co,vns*vnsp))
  endif
  
  SAFE_ALLOCATE(indq, (vng))
  SAFE_ALLOCATE(phq, (vng))
  SAFE_ALLOCATE(isortq, (gvec%ng))

  if (ik .ne. ikold .and. ik .ne. ikold2 .and. ifirst .ne. 0) then
    SAFE_DEALLOCATE(indq_old2)
    SAFE_DEALLOCATE(phq_old2)
    SAFE_DEALLOCATE(isortq_old2)
    SAFE_ALLOCATE(indq_old2, (vngold))
    SAFE_ALLOCATE(phq_old2, (vngold))
    SAFE_ALLOCATE(isortq_old2, (gvec%ng))

    indq_old2 = indq_old
    phq_old2 = phq_old
    isortq_old2 = isortq_old

    SAFE_DEALLOCATE(indq_old)
    SAFE_DEALLOCATE(phq_old)
    SAFE_DEALLOCATE(isortq_old)
    SAFE_ALLOCATE(indq_old, (vng))
    SAFE_ALLOCATE(phq_old, (vng))
    SAFE_ALLOCATE(isortq_old, (gvec%ng))
  endif
  if (ifirst .eq. 0) then
    SAFE_ALLOCATE(indq_old, (vng))
    SAFE_ALLOCATE(phq_old, (vng))
    SAFE_ALLOCATE(isortq_old, (gvec%ng))
    SAFE_ALLOCATE(indq_old2, (vng))
    SAFE_ALLOCATE(phq_old2, (vng))
    SAFE_ALLOCATE(isortq_old2, (gvec%ng))
  endif
  
! Initalize parameters in variable wfnv

  wfnv%ng=vng
  wfnv%nband=1
  wfnv%nspin=vns
  wfnv%nspinor=vnsp
  if (intwfnv%nspin.ne.nspin) then
    write(errmsg,'(2a,2i2)') 'spin number mismatch: ', nspin, vns
    call die(errmsg, only_root_writes = .true.)
  endif

  if (xct%ivpar .eq. 1) then
    SAFE_ALLOCATE(wfnv%cg, (wfnv%ng,1,wfnv%nspin*wfnv%nspinor))
    vcg(1:vng,1,:) = intwfnv%cg(1:vng,ivown,:)
  else
    SAFE_ALLOCATE(wfnv%cg, (wfnv%ng,xct%nvb_co,wfnv%nspin*wfnv%nspinor))        
    vcg(1:vng,:,:) = intwfnv%cg(1:vng,ivown:ivown+xct%nvb_co-1,:)
  endif
  SAFE_ALLOCATE(wfnv%isort, (gvec%ng))
  
  isortq(:) = intwfnv%isort(:,ikown)

! Compute inverse index array of Fourier components around rk-kpoint

  if (ik .ne. ikold .and. ik .ne. ikold2) then
    SAFE_ALLOCATE(isorti, (gvec%ng))
    isorti(:)=0
    do ii=1,wfnv%ng
      isorti(isortq(ii))=ii
    enddo
  endif

! Compute index array of Fourier components around fk-kpoint

  if (ik .ne. ikold .and. ik .ne. ikold2) then
    if (peinf%verb_debug) then
      write(6,*) 'ikneikold',peinf%inode,ik,ikold
    endif
    SAFE_ALLOCATE(ekin, (gvec%ng))
    if (xct%qflag .eq. 0) then
      qk(1:3)=kgq%f(1:3,indexq(ik))
    elseif (xct%qflag .eq. 2) then
      qk(1:3)=kg%f(1:3,indexq(ik))
    else
      qk(1:3)=kg%f(1:3,ik)
    endif
    call kinetic_energies(gvec, crys%bdot, ekin, qvec = qk)
    call sortrx(gvec%ng,ekin,isortq)
    SAFE_DEALLOCATE(ekin)
    isortq_old(:)=isortq(:)
  elseif (ik .eq. ikold) then
    isortq(:)=isortq_old(:)
  else 
    isortq(:)=isortq_old2(:)
  endif

! Find ind and ph relating wavefunctions in fk to rk-kpoint

  if (ik .ne. ikold .and. ik .ne. ikold2) then
    indq=0
    phq=ZERO
    if (xct%qflag .eq. 0) then
      call gmap(gvec,syms,wfnv%ng,kgq%itran(indexq(ik)), &
        kgq%kg0(:,indexq(ik)),isortq,isorti,indq,phq,.true.)
    elseif (xct%qflag .eq. 2) then
      call gmap(gvec,syms,wfnv%ng,kg%itran(indexq(ik)), &
        kgq%kg0(:,ik),isortq,isorti,indq,phq,.true.)
    else
      call gmap(gvec,syms,wfnv%ng,kg%itran(ik), &
        kg%kg0(:,ik),isortq,isorti,indq,phq,.true.)
    endif
    SAFE_DEALLOCATE(isorti)
    indq_old(:)=indq(:)
    phq_old(:)=phq(:)
  else if (ik .eq. ikold) then
    indq(:)=indq_old(:)
    phq(:)=phq_old(:)
  else
    indq(:)=indq_old2(:)
    phq(:)=phq_old2(:)
  endif

! Compute and renormalize valence wavefunctions in fk-kpoint

  if (xct%ivpar .eq. 1) invband = 1
  if (xct%ivpar .eq. 0) invband = xct%nvb_co
  
  do ijk = 1, invband
    do kk=1,wfnv%nspin
      do jj=1,wfnv%nspinor
        do ii=1,wfnv%ng
          if (indq(ii) .gt. 0) then
            wfnv%cg(ii,ijk,kk*jj)=phq(ii)*vcg(indq(ii),ijk,kk*jj)
          else
            wfnv%cg(ii,ijk,kk*jj) = ZERO
          endif
        enddo !ii
      enddo ! jj
      call checknorm('wfnv%cg',ijk,ik,wfnv%ng,kk,wfnv%nspinor,wfnv%cg(1:wfnv%ng,ijk,:))
    enddo ! kk
  enddo ! ijk
  vcg=wfnv%cg

  !   In spinor case, we must rotate spinors according to spinor rotation matrix umtrx
#ifdef CPLX
  if (vnsp.eq.2) then
    do ig=1,vng
      if (xct%ivpar .eq. 1) then
        SAFE_ALLOCATE(zvtemp, (1,vns*vnsp))
      else
        SAFE_ALLOCATE(zvtemp, (xct%nvb_co,vns*vnsp))
      endif
      call susymmetries(crys%bvec(:,:),syms%mtrx(:,:,kg%itran(ik)),umtrx,kg%itran(ik))
      zvtemp = matmul(vcg(ig,:,:),umtrx)
      vcg(ig,:,:) = zvtemp
      SAFE_DEALLOCATE(zvtemp)
    enddo
  endif
#endif

  wfnv%cg=vcg
  wfnv%isort=isortq
  
  
  SAFE_DEALLOCATE(vcg)
  SAFE_DEALLOCATE(indq)
  SAFE_DEALLOCATE(phq)
  SAFE_DEALLOCATE(isortq)
  
!-----------------------------------------------------------------------
! Deal with the conduction wavefunctions

  icown = peinf%ipec(peinf%inode+1,ic,kg%indr(ik))
  ikown = peinf%ipek(peinf%inode+1,kg%indr(ik))
  
  cng = intwfnc%ng(ikown)
  cnb = 1
  cns = intwfnc%nspin
  cnsp = intwfnc%nspinor
  
  if (xct%icpar .eq. 1) then
    SAFE_ALLOCATE(ccg, (cng,1,cns*cnsp))
  else
    SAFE_ALLOCATE(ccg, (cng,xct%ncb_co,cns*cnsp))
  endif
  
  SAFE_ALLOCATE(ind, (cng))
  SAFE_ALLOCATE(ph, (cng))
  SAFE_ALLOCATE(isort, (gvec%ng))

  if (ik .ne. ikold .and. ik .ne. ikold2 .and. ifirst .ne. 0) then
    SAFE_DEALLOCATE(ind_old2)
    SAFE_DEALLOCATE(ph_old2)
    SAFE_DEALLOCATE(isort_old2)
    SAFE_ALLOCATE(ind_old2, (cngold))
    SAFE_ALLOCATE(ph_old2, (cngold))
    SAFE_ALLOCATE(isort_old2, (gvec%ng))

    ind_old2 = ind_old
    ph_old2 = ph_old
    isort_old2 = isort_old

    SAFE_DEALLOCATE(ind_old)
    SAFE_DEALLOCATE(ph_old)
    SAFE_DEALLOCATE(isort_old)
    SAFE_ALLOCATE(ind_old, (cng))
    SAFE_ALLOCATE(ph_old, (cng))
    SAFE_ALLOCATE(isort_old, (gvec%ng))
  endif
  if (ifirst .eq. 0) then
    SAFE_ALLOCATE(ind_old, (cng))
    SAFE_ALLOCATE(ph_old, (cng))
    SAFE_ALLOCATE(isort_old, (gvec%ng))
    SAFE_ALLOCATE(ind_old2, (cng))
    SAFE_ALLOCATE(ph_old2, (cng))
    SAFE_ALLOCATE(isort_old2, (gvec%ng))
  endif
  
  wfnc%ng=cng
  wfnc%nband=1
  wfnc%nspin=cns
  wfnc%nspinor=cnsp
  
  if (cns.ne.nspin) then
    write(errmsg,'(2a,2i2)') 'spin number mismatch: ', nspin, cns
    call die(errmsg, only_root_writes = .true.)
  endif

  SAFE_ALLOCATE(wfnc%isort, (gvec%ng))
  
  isort(:)=intwfnc%isort(:,ikown)
  if (xct%icpar .eq. 1) then
    SAFE_ALLOCATE(wfnc%cg, (wfnc%ng,1,wfnc%nspin*wfnc%nspinor))
    ccg(1:cng,1,:) = intwfnc%cg(1:cng,icown,:)
  else
    SAFE_ALLOCATE(wfnc%cg, (wfnc%ng,xct%ncb_co,wfnc%nspin*wfnc%nspinor))
    ccg(1:cng,:,:) = intwfnc%cg(1:cng,icown:icown+xct%ncb_co-1,:)
  endif

! JRD: Below is now necessary because kg might be different from kgq
! Compute inverse index array of Fourier components around rk-kpoint

  if (ik .ne. ikold .and. ik .ne. ikold2) then
    SAFE_ALLOCATE(isorti, (gvec%ng))
    isorti(:)=0
    do ii=1,wfnc%ng
      isorti(isort(ii))=ii
    enddo
  endif
  
! Compute index array of Fourier components around fk-kpoint

  if (ik .ne. ikold .and. ik .ne. ikold2) then
    SAFE_ALLOCATE(ekin, (gvec%ng))
    call kinetic_energies(gvec, crys%bdot, ekin, qvec = kg%f(:, ik))
    call sortrx(gvec%ng,ekin,isort)
    SAFE_DEALLOCATE(ekin)
    isort_old(:) = isort(:)
  else if (ik .eq. ikold) then
    isort(:)=isort_old(:)
  else
    isort(:)=isort_old2(:)
  endif

! Find ind and ph relating wavefunctions in fk to rk-kpoint

  if (ik .ne. ikold .and. ik .ne. ikold2) then
    ind=0
    ph=ZERO
    call gmap(gvec,syms,wfnc%ng,kg%itran(ik), &
      kg%kg0(:,ik),isort,isorti,ind,ph,.true.)
    SAFE_DEALLOCATE(isorti)
    ind_old(:)=ind(:)
    ph_old(:)=ph(:)
  else if (ik .eq. ikold) then
    ind(:)=ind_old(:)
    ph(:)=ph_old(:)
  else
    ind(:)=ind_old2(:)
    ph(:)=ph_old2(:)
  endif

! Compute and renormalize conduction wavefunctions

  if (xct%icpar .eq. 1) incband =1
  if (xct%icpar .eq. 0) incband = xct%ncb_co
  
  do ijk =1, incband
    do kk=1,wfnc%nspin
      do jj=1,wfnc%nspinor
        do ii=1,wfnc%ng
          if (ind(ii) .gt. 0) then
            wfnc%cg(ii,ijk,kk*jj)=ph(ii)*ccg(ind(ii),ijk,kk*jj)
          else
            wfnc%cg(ii,ijk,kk*jj)=ZERO
          endif
        enddo
      enddo ! jj
      call checknorm('wfnc%cg',ijk,ik,wfnc%ng,kk,wfnc%nspinor,wfnc%cg(1:wfnc%ng,ijk,:))
    enddo ! kk
  enddo ! ijk
  ccg=wfnc%cg

  !   In spinor case, we must rotate spinors according to spinor rotation matrix umtrx
#ifdef CPLX
  if (cnsp.eq.2) then
    do ig=1,cng
      if (xct%icpar .eq. 1) then
        SAFE_ALLOCATE(zctemp, (1,cns*cnsp))
      else
        SAFE_ALLOCATE(zctemp, (xct%ncb_co,cns*cnsp))
      endif
      call susymmetries(crys%bvec(:,:),syms%mtrx(:,:,kg%itran(ik)),umtrx,kg%itran(ik))
      zctemp = matmul(ccg(ig,:,:),umtrx)
      ccg(ig,:,:) = zctemp
      SAFE_DEALLOCATE(zctemp)
    enddo
  endif
#endif

  wfnc%cg=ccg
  wfnc%isort=isort
  
  SAFE_DEALLOCATE(ccg)
  SAFE_DEALLOCATE(ind)
  SAFE_DEALLOCATE(ph)
  SAFE_DEALLOCATE(isort)
  
  if (ik .ne. ikold .and. ik .ne. ikold2) then
    ikold2=ikold
    ikold=ik
    cngold=cng
    vngold=vng
  endif

  ifirst=-1
  
  POP_SUB(genwf_kernel)
  
  return
end subroutine genwf_kernel

end module genwf_kernel_m
