!================================================================================
!
! Routines:
!
! (1) mtxel_md()        Originally By MLT               Last Modified 7/1/2008 (JRD)
!
!     Calculates matrix elements based on model.  Is currently not used in code.
!     Was replaced by ab initio matrix elements in kernel.
!
!================================================================================

#include "f_defs.h"

module mtxel_md_m

  use global_m
  implicit none

  private

  public :: &
    mtxel_md

contains

subroutine mtxel_md(crys,gvec,syms,kg,wfnc,wfncp, &
  wfnvp,wfnv,xct,bsedbody,bsedhead,bsedwing,bsex,ii,iip,ipe)
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (symmetry), intent(in) :: syms
  type (grid), intent(in) :: kg
  type (wavefunction), intent(in) :: wfnc,wfncp,wfnvp,wfnv
  type (xctinfo), intent(in) :: xct
  SCALAR, intent(inout) :: &
    bsedbody(xct%nvband,peinf%nck,xct%nvband,peinf%nckpe,xct%nspin,xct%nspin), &
    bsedhead(xct%nvband,peinf%nck,xct%nvband,peinf%nckpe,xct%nspin,xct%nspin), &
    bsedwing(xct%nvband,peinf%nck,xct%nvband,peinf%nckpe,xct%nspin,xct%nspin), &
    bsex(xct%nvband,peinf%nck,xct%nvband,peinf%nckpe,xct%nspin,xct%nspin)
  integer, intent(in) :: ii,iip,ipe

  integer :: jj,jjp,kk,i1,i2,i3,irq
  integer :: ic,ick,ickp,icp,ik,ikp,ifq,ig,ign,igp, &
    igadd,igaddn,igaddp,isc,iscp,isv,isvp,iv,ivp
  integer :: ngq,neps
  integer :: kxm,kxn,kym,kyn,kzm,kzn,kaddn,g0(3),ifqmin
  real(DP) :: ekin,wval,fq(3),qq,fqq(3),diffl,difmin,diff(3)
  real(DP) :: fqmin(3),lmin,length,epsinv
      
  integer, allocatable :: isorti(:)
  SCALAR :: bsed
  SCALAR, allocatable :: wpt(:),mcc(:,:,:,:),mvv(:,:,:,:),mvc(:,:,:,:),mvpcp(:,:,:,:)
  
! Compute maximum size of the cube

  PUSH_SUB(mtxel_md)

  kxm=gvec%FFTgrid(1)
  kym=gvec%FFTgrid(2)
  kzm=gvec%FFTgrid(3)
  ik=peinf%ik(peinf%inode+1,ii)
  ikp=peinf%ik(ipe+1,iip)

!-------------------------------------
! Read epsinv for q + g0 =kp-k

  fq(:)=kg%f(:,ik)-kg%f(:,ikp)
  fqq(:) = fq(:)

! Reduce to BZ

  jj=2
  lmin = 1.d10
  do i1=-jj,jj
    diff(1) = fq(1) - i1
    do i2=-jj,jj
      diff(2) = fq(2) - i2
      do i3=-jj,jj
        diff(3) = fq(3) - i3
        length = DOT_PRODUCT(diff,MATMUL(crys%bdot,diff))
        if (length.lt.lmin) then
          lmin = length
          fqmin(:) = diff(:)
        endif
      enddo
    enddo
  enddo
  fq(:) = fqmin(:)
  g0(:) = fqq(:) - fq(:)
      
! Model

  SAFE_ALLOCATE(wpt, (xct%neps))
  wpt=0.d0

!-------------------------------------
! Compute W(g,gp) for q + g0 =kp-k

! head for q=0 and q<>0, *no interpolation* !!
  wpt(1)= 1.d0
! body for q=0 and qq<>0
  do igp=2,xct%neps
    fqq(:) = fq(:) + gvec%components(:,igp)
    ekin= DOT_PRODUCT(fqq,MATMUL(crys%bdot,fqq))
    call epsmodel(crys,fqq,epsinv)
    wpt(igp)=epsinv/ekin
  enddo

!-------------------------------------
! Compute direct term matrices

  SAFE_ALLOCATE(mcc, (xct%ng,xct%ncband,xct%ncband,xct%nspin))
  SAFE_ALLOCATE(mvv, (xct%ng,xct%nvband,xct%nvband,xct%nspin))
  mcc=0.0d0
  mvv=0.0d0

! Compute inverse array to wfnc%isort

  SAFE_ALLOCATE(isorti, (gvec%ng))
  isorti=0
  do ig=1,gvec%ng
    isorti(wfnc%isort(ig))=ig
  enddo

! Loop G

  do ig=1,xct%ng
    igadd=ig

! Loop over G`

    do igp=1,wfncp%ng
      igaddp=wfncp%isort(igp)

! Compute address of Gn=Gp+G-G0:
! if the result is beyond planewaves included
! in wavefunction, skip

      kxn=gvec%components(1,igaddp)+gvec%components(1,igadd)-g0(1)+gvec%FFTgrid(1)/2+1
      if((kxn.lt.1).or.(kxn.gt.kxm)) cycle
      kyn=gvec%components(2,igaddp)+gvec%components(2,igadd)-g0(2)+gvec%FFTgrid(2)/2+1
      if((kyn.lt.1).or.(kyn.gt.kym)) cycle
      kzn=gvec%components(3,igaddp)+gvec%components(3,igadd)-g0(3)+gvec%FFTgrid(3)/2+1
      if((kzn.lt.1).or.(kzn.gt.kzm)) cycle
      kaddn=((kxn-1)*kym+kyn-1)*kzm+kzn
      igaddn=gvec%index_vec(kaddn) ! relate the cube to the sphere
      ign=isorti(igaddn)
      if(ign.gt.wfnc%ng) cycle

! Compute matrix element: <c`k`|exp(i(k`-k-G+G0).r)|ck>

      do jj=ii,ii+peinf%nckmem-1
        do jjp=iip,iip+peinf%nckmem-1
          ic=peinf%ic(peinf%inode+1,jj)
          icp=peinf%ic(ipe+1,jjp)
          ick=peinf%ick(ic,ik)
          ickp=peinf%ick(icp,ikp)
          do isc=1,xct%nspin
            mcc(ig,ic,icp,isc)=mcc(ig,ic,icp,isc)+ &
              MYCONJG(wfncp%cg(igp,icp,isc))*(wfnc%cg(ign,ic,isc))
          enddo
        enddo
      enddo

! Compute matrix element: <v`k`|exp(i(k`-k-G+G0).r)|vk>

      do iv=1,xct%nvband
!          do ivp=iv,xct%nvband
        do ivp=1,xct%nvband
          do isv=1,xct%nspin
            mvv(ig,iv,ivp,isv)=mvv(ig,iv,ivp,isv)+ &
              MYCONJG(wfnvp%cg(igp,ivp,isv))*wfnv%cg(ign,iv,isv)
          enddo
        enddo
      enddo

    enddo !end loop over gp

  enddo !end loop over g

  SAFE_DEALLOCATE(isorti)

  do ig=1,xct%neps
    wval=wpt(ig)
    if(abs(wval).lt.TOL_Zero) cycle
    do jj=ii,ii+peinf%nckmem-1
      do jjp=iip,iip+peinf%nckmem-1
        ic=peinf%ic(peinf%inode+1,jj)
        icp=peinf%ic(ipe+1,jjp)
        ick=peinf%ick(ic,ik)
        ickp=peinf%ick(icp,ikp)
        do iv=1,xct%nvband
          do ivp=1,xct%nvband
            do isc=1,xct%nspin
              do isv=1,xct%nspin
                bsed=wval*MYCONJG(mcc(ig,ic,icp,isc))*mvv(ig,iv,ivp,isv)
                if(ig.ne.1) then
                  bsedbody(ivp,ickp,iv,jj,isc,isv)= &
                    bsedbody(ivp,ickp,iv,jj,isc,isv)+bsed
                else
                  bsedhead(ivp,ickp,iv,jj,isc,isv)= &
                    bsedhead(ivp,ickp,iv,jj,isc,isv)+bsed
                endif
              enddo !isv
            enddo !isc
          enddo !ivp
        enddo !iv
      enddo !icp
    enddo !ic
  enddo !ig
  
  SAFE_DEALLOCATE(mcc)
  SAFE_DEALLOCATE(mvv)
  SAFE_DEALLOCATE(wpt)

!-------------------------------------
! Compute exchange term matrices

  SAFE_ALLOCATE(mvc, (xct%ng,xct%nvband,xct%ncband,xct%nspin))
  mvc=0.0d0

! Compute inverse array to wfnc%isort

  SAFE_ALLOCATE(isorti, (gvec%ng))
  isorti=0
  do ig=1,gvec%ng
    isorti(wfnc%isort(ig))=ig
  enddo

! Loop G

  do ig=2,xct%ng
    igadd=ig

! Loop over G`

    do igp=1,wfnv%ng
      igaddp=wfnv%isort(igp)

! Compute address of gn=gp-g:
! If the result is beyond planewaves included
! in wavefunction, skip

      kxn=gvec%components(1,igaddp)-gvec%components(1,igadd)+gvec%FFTgrid(1)/2+1
      if((kxn.lt.1).or.(kxn.gt.kxm)) cycle
      kyn=gvec%components(2,igaddp)-gvec%components(2,igadd)+gvec%FFTgrid(2)/2+1
      if((kyn.lt.1).or.(kyn.gt.kym)) cycle
      kzn=gvec%components(3,igaddp)-gvec%components(3,igadd)+gvec%FFTgrid(3)/2+1
      if((kzn.lt.1).or.(kzn.gt.kzm)) cycle
      kaddn=((kxn-1)*kym+kyn-1)*kzm+kzn
      igaddn=gvec%index_vec(kaddn) ! relate the cube to the sphere
      ign=isorti(igaddn)
      if(ign.gt.wfnc%ng) cycle

! Compute matrix elements: <vk|exp(iG.r)|ck>  ,  G.ne.0

      do jj=ii,ii+peinf%nckmem-1
        ic=peinf%ic(peinf%inode+1,jj)
        do iv=1,xct%nvband
          do isc=1,xct%nspin

! MLT: valence,conduction states should have the same spin!

!              isv=xct%nspin+1-isc
            isv=isc
            mvc(ig,iv,ic,isc)=mvc(ig,iv,ic,isc)+ &
              MYCONJG(wfnv%cg(igp,iv,isv))*wfnc%cg(ign,ic,isc)
          enddo
        enddo
      enddo
      
    enddo !end loop over gp
  enddo !end loop over g

  SAFE_DEALLOCATE(isorti)
  
  SAFE_ALLOCATE(mvpcp, (xct%ng,xct%nvband,xct%ncband,xct%nspin))
  mvpcp=0.0d0

! Compute inverse array to wfncp%isort

  SAFE_ALLOCATE(isorti, (gvec%ng))
  isorti=0
  do ig=1,gvec%ng
    isorti(wfncp%isort(ig))=ig
  enddo

! Loop G

  do ig=2,xct%ng
    igadd=ig

! Loop over G`

    do igp=1,wfnvp%ng
      igaddp=wfnvp%isort(igp)

! Compute address of gn=gp-g:
! if the result is beyond planewaves included
! in wavefunction, skip

      kxn=gvec%components(1,igaddp)-gvec%components(1,igadd)+gvec%FFTgrid(1)/2+1
      if((kxn.lt.1).or.(kxn.gt.kxm)) cycle
      kyn=gvec%components(2,igaddp)-gvec%components(2,igadd)+gvec%FFTgrid(2)/2+1
      if((kyn.lt.1).or.(kyn.gt.kym)) cycle
      kzn=gvec%components(3,igaddp)-gvec%components(3,igadd)+gvec%FFTgrid(3)/2+1
      if((kzn.lt.1).or.(kzn.gt.kzm)) cycle
      kaddn=((kxn-1)*kym+kyn-1)*kzm+kzn
      igaddn=gvec%index_vec(kaddn) ! relate the cube to the sphere
      ign=isorti(igaddn)
      if(ign.gt.wfncp%ng) cycle

! Compute matrix element:
!   mvpcp = <v`k`|exp(iG.r)|c`k`> ,   G.ne.0

      do jjp=iip,iip+peinf%nckmem-1
        icp=peinf%ic(ipe+1,jjp)
        do ivp=1,xct%nvband
          do iscp=1,xct%nspin

! MLT: valence,conduction states should have the same spin!

!              isvp=xct%nspin+1-iscp
            isvp=iscp
            mvpcp(ig,ivp,icp,iscp)=mvpcp(ig,ivp,icp,iscp)+ &
              MYCONJG(wfnvp%cg(igp,ivp,isvp))*wfncp%cg(ign,icp,iscp)
          enddo
        enddo
      enddo
      
    enddo !end loop over gp
  enddo !end loop over g

  SAFE_DEALLOCATE(isorti)
  
  do ig=2,xct%neps
    ekin=gvec%ekin(ig)
    do jj=ii,ii+peinf%nckmem-1
      do jjp=iip,iip+peinf%nckmem-1
        ic=peinf%ic(peinf%inode+1,jj)
        icp=peinf%ic(ipe+1,jjp)
        ick=peinf%ick(ic,ik)
        ickp=peinf%ick(icp,ikp)
        do iv=1,xct%nvband
          do ivp=1,xct%nvband
            do isc=1,xct%nspin
              do iscp=1,xct%nspin
                bsex(ivp,ickp,iv,jj,isc,iscp)= &
                  bsex(ivp,ickp,iv,jj,isc,iscp)+1.0d0/ekin* &
                  mvpcp(ig,ivp,icp,iscp)*MYCONJG(mvc(ig,iv,ic,isc))
              enddo !iscp
            enddo !isc
          enddo !ivp
        enddo !iv
      enddo !icp
    enddo !ic
  enddo !ig
  
  SAFE_DEALLOCATE(mvc)
  SAFE_DEALLOCATE(mvpcp)
  
  POP_SUB(mtxel_md)
  
  return
end subroutine mtxel_md

end module mtxel_md_m
