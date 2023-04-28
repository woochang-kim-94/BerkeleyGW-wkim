!===============================================================================
!
! Routines:
!
! (1) epsmodel()        Originally By MLT               Last Modified 7/3/2008 (JRD)
!
!     Model dielectric function for GaAs
!     see G. Cappellini et al, PRB 47, 9892 (1993)
!
! (2) epsmodel_pac()    Orignally By MLT                Last Modified 7/3/2008 (JRD)
!
!     Model dielectric function for polyacetylene
!     see B. Wenzien et al, PRB 51, 14701 (1995)
!     omega_plasma = 2.0 (actual value)
!     q_ThomasFermi = 1.0
!     epsp : dielectric constant perpendicular to chain
!     epsc : dielectric constant along chain
!
! (3) epsmodel_cnt()    Orignally By Li Yang            Last Modified 7/3/2008 (JRD)
!
!     Model dielectric function for carbon nanotube
!     zcut=inf, rho=15.11, grid 1x1x64
!
!  These routines are not used in current version of the code.
!  Calculated dielectric function from Epsilon is used instead.
!
!===============================================================================

#include "f_defs.h"

module epsmodel_m

  use global_m
  implicit none

  public ::       &
    epsmodel,     &
    epsmodel_pac, &
    epsmodel_cnt

contains

subroutine epsmodel(bdot,fq,epsinv)
  real(DP), intent(in) :: bdot(3,3),fq(3)
  real(DP), intent(out) :: epsinv

  real(DP) :: eps0,alpha,beta,q2

  PUSH_SUB(epsmodel)

  q2= DOT_PRODUCT(fq,MATMUL(bdot,fq))
!      eps0 = 0.0747
!      alpha = 1.37
!      beta = 0.76
  eps0 = 0.0971
  alpha = 1.261
  beta = 0.67
  epsinv = 1.d0/(1.d0 + 1.d0/(eps0 + alpha*q2 + beta*q2*q2))

  POP_SUB(epsmodel)

  return
end subroutine epsmodel


!==============================================================================
subroutine epsmodel_pac(bdot,fq,epsinv)
  real(DP), intent(in) :: bdot(3,3),fq(3)
  real(DP), intent(out) :: epsinv
  ! epsinv was complex originally, but the calc here seems to give always a real result

  integer :: jj
  real(DP) :: epsp,epsc,alphac,alphap,sum,qc,qp,cost2,sint2,qq

  PUSH_SUB(epsmodel_pac)

  epsp=3.7
  epsc=20.0
  alphac=8.0
  alphap=1.0
  
  qq= sqrt( DOT_PRODUCT(fq,MATMUL(bdot,fq)) )
  if (abs(qq).lt.1.d-5) then
    epsinv=0.05
  else
    sum=0.0
    do jj=1,3
      sum = sum + fq(jj)*bdot(jj,3)
    enddo
    qc=sum/sqrt(bdot(3,3))
    qp=sqrt(abs(qq**2 - qc**2))
    cost2=(qc/qq)**2
    sint2=(qp/qq)**2
    epsinv = cost2/(1.d0+1.d0/(1.d0/(epsc-1)+alphac*qc*qc+ &
      qq**4/8.d0)) + sint2/(1.d0+1.d0/(1.d0/(epsp-1d0)+ &
      alphap*qp*qp+qq**4/8.d0))
  endif
  
  POP_SUB(epsmodel_pac)
  
  return
end subroutine epsmodel_pac


!==============================================================================
subroutine epsmodel_cnt(diff,epsinv)
  real(DP), intent(in) :: diff(3)
  real(DP), intent(out) :: epsinv

  real(DP) :: a0,a1,a2,x,length2,qz2,qp2,epsz,epsp,v
  
  PUSH_SUB(epsmodel_cnt)
  
  a0 = -248.187d0
  a1 = 921.851d0
  a2 = 347.231d0
  
  x=abs(diff(3))
  
  if(x.lt.1.0d-15) then
    epsinv=1.d0
  else
    epsinv = (1+a2*x*x)/(1+a0*x*x*log(x*0.7889*15.11)+a1*x*x)
  endif
  
  POP_SUB(epsmodel_cnt)
  
  return
end subroutine epsmodel_cnt

end module epsmodel_m
