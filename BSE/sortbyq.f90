!===================================================================================
!
! Routines:
!
! (1) sortbyq      Originally By MLT       Last Modified: 9/27/2010 (JRD)
!
!===================================================================================

#include "f_defs.h"

module sortbyq_m

  use global_m
  implicit none

  private

  public :: &
    sortbyq

contains

subroutine sortbyq(fqa,qqa,g0a,ifqa,irqa,qg,kg,crys)
  real(DP), intent(out) :: fqa(3,peinf%myown), qqa(peinf%myown)
  integer, intent(out) :: g0a(3,peinf%myown)
  integer, intent(out) :: ifqa(peinf%myown)
  integer, intent(out) :: irqa(peinf%myown)
  type (grid), intent(in) :: qg, kg
  type (crystal), intent(in) :: crys
  
  integer :: sinv,ik,ic,iv,ikp,icp,ivp,ii,ifq,irq,ifqmin
  
  real(DP) :: fq(3),difmin,diff(3),length
  
  PUSH_SUB(sortbyq)
  
  do ii = 1, peinf%myown
    
    ik=peinf%ik(peinf%inode+1,ii)
    ic=peinf%ic(peinf%inode+1,ii)
    iv=peinf%iv(peinf%inode+1,ii)
    ikp=peinf%ikp(peinf%inode+1,ii)
    icp=peinf%icp(peinf%inode+1,ii)
    ivp=peinf%ivp(peinf%inode+1,ii)

! If ik.gt.ikp, define q + g0 = kp-k
! all factors of (G0-G) and (q) must get a minus sign,
! and epsinv must be replaced by its complex conjugate
! otherwise, work with q + g0 = k-kp
! sinv keeps track of signs (in the comments below, assume ik.le.ikp)

    sinv = 1
    if (ik.gt.ikp) sinv = -1

    fq(:)=dble(sinv)*( kg%f(:,ik)-kg%f(:,ikp) )

! For a given pair of k-points (ik,ikp), calculate the shortest
! q = k - kp + G_0 , with G_0 a reciprocal lattice vector

    difmin=1.d10
    ifqmin=0     
    do ifq=1,qg%nf
      diff(1:3) = fq(1:3) - qg%f(1:3, ifq)
      diff(1:3) = diff(1:3) - anint( diff(1:3) )
      
      length = DOT_PRODUCT(diff,MATMUL(crys%bdot,diff))
      
! Found a good matching point. Save it

      if (length <= TOL_Zero) then
        ifqmin = ifq
        difmin = 0.d0
        exit
      endif

! Check other points for the best matching point

      if (length.lt.difmin) then  
        difmin = length
        ifqmin = ifq
      endif
    enddo ! ifq
    
    ifq=ifqmin
    
    if (difmin.gt.TOL_Zero .or. ifq .eq. 0) then
      write(0,'(a,3f12.6,a,i6,a,i6,a)') 'epsmat does not contain q-point ', fq, ' (= k ', ik, ' - k ', ikp, ')'
      write(0,'(a,i6,a,f10.5)') 'Closest q-point is', ifq, ' of distance ', difmin
      call die("Needed q-point not found. This likely means your WFN_co file is shifted, which is not recommended.")
    endif
    
    fqa(:,ii) = qg%f(:,ifq)
    g0a(:,ii) = anint(fq(:) - fqa(:,ii))
    ifqa(ii)=ifq
    irqa(ii)=qg%indr(ifq)
    
    qqa(ii) = sqrt( DOT_PRODUCT(fqa(:,ii),MATMUL(crys%bdot,fqa(:,ii))) )
    
  enddo ! ii
  
  SAFE_ALLOCATE(peinf%nxqown, (qg%nr))
  SAFE_ALLOCATE(peinf%nxqi, (peinf%myown))
  peinf%nxqown = 0
  peinf%nxqi = 0
  
  do irq = 1, qg%nr
    peinf%nxqown(irq)=0
    do ii = 1, peinf%myown
      if (irqa(ii) .eq. irq) then
        peinf%nxqown(irq)=peinf%nxqown(irq)+1
        peinf%nxqi(peinf%nxqown(irq))=ii
      endif
    enddo
  enddo
  
  POP_SUB(sortbyq)
  
  return
end subroutine sortbyq

end module sortbyq_m
