!=======================================================================
!
! Routines:
!
! (1) minibzaverage_3d_oneoverq2() Originally by JRD/MJ Last Modified: 8/27/2009 (MJ/JRD)
!
! (2) minibzaverage_3d_oneoverq() Originally by JRD/MJ  Last Modified: 8/27/2009 (MJ/JRD)
!
! (3) minibzaverage_2d_oneoverq2() Originally by JRD/MJ Last Modified: 9/15/2009 (MJ/JRD)
!
! (4) minbizaverage_1d() Originally by JRD/MJ           Last Modified: 8/27/2009 (MJ/JRD)
!
!  Output: average of <V_q> on the mini-BZ for a 1-D system.
!  output units: units equivalent to 8Pi/q^2
!
!=======================================================================

#include "f_defs.h"

module minibzaverage_m

  use global_m
  use bessel_m
  use misc_m

  implicit none
  
  private

  public :: minibzaverage_3d_oneoverq2, minibzaverage_3d_oneoverq, &
    minibzaverage_2d_oneoverq2, minibzaverage_1d, minibzaverage_3d_oneoverq2_mod

contains

subroutine minibzaverage_3d_oneoverq2(nn,bdot,integral,qran,qk,averagew,epshead,wcoul0,q0sph2,celvol,nfk)
  integer, intent(in) :: nn
  real(DP), intent(in) :: bdot(3,3), qran(:,:) !< (3, nn)
  real(DP), intent(out) :: integral
  real(DP), intent(in) ::  qk(3)
  logical, intent(in)  :: averagew
  SCALAR, intent(in) :: epshead  
  real(DP), intent(in) ::  q0sph2
  real(DP), intent(in) :: celvol
  integer, intent(in) :: nfk
  real(DP) :: gkq(3), length,length_qk
  SCALAR, intent(inout)  :: wcoul0 
 
  integer :: ii, nn2

  PUSH_SUB(minibzaverage_3d_oneoverq2)

  integral = 0D0
  length_qk = DOT_PRODUCT(qk,MATMUL(bdot,qk))

  if( length_qk < TOL_Zero ) then
    nn2 = nn
    do ii = 1, nn2
      gkq(:) = qk(:) + qran(:,ii)
      length = DOT_PRODUCT(gkq,MATMUL(bdot,gkq))
      ! Skip the small value of q that will be integrated analytically later on
      if ( length > q0sph2 ) integral = integral + 1D0/length
    enddo
  else
    ! FHJ: for spherical integration regions, one can make the error per MC
    ! integration ~const. by choosing the number of points such that N ~ 1/ekinx.
    ! This is because, in 3D, error = sigma/N^{3/2}, and sigma ~ 1/ekinx^{3/2}
    ! If we fix the number of points such that N(ekinx=4*q0sph2) = nmc_coarse,
    nn2 = idnint(nmc_coarse * 4d0 * q0sph2 / length_qk)
    nn2 = max(1, min(nn2, nn))
    do ii = 1, nn2
      gkq(:) = qk(:) + qran(:,ii)
      length = DOT_PRODUCT(gkq,MATMUL(bdot,gkq))
      integral = integral + 1D0/length
    enddo
  endif

  integral = integral * 8D0 * PI_D / dble(nn2)

  if( length_qk < TOL_Zero ) then
    integral = integral + 32.0D0 * PI_D**2 * SQRT(q0sph2) / ( 8.0D0 * PI_D**3 / (celvol * dble(nfk)) )
  endif

  if (length_qk .lt. TOL_Zero .and. averagew) then
    wcoul0 = integral * epshead
  endif
  
  POP_SUB(minibzaverage_3d_oneoverq2)
  
  return
end subroutine minibzaverage_3d_oneoverq2


!========================================================================

! This is for Slab Truncation

subroutine minibzaverage_2d_oneoverq2(nn,bdot,integral,qran,qk,kz,zc,epshead,q0len,averagew,wcoul0)
  integer, intent(in) :: nn
  real(DP), intent(in) :: bdot(3,3), qran(:,:) !< (3,nn)
  real(DP), intent(out) :: integral
  real(DP), intent(in) ::  qk(3)
  logical, intent(in)  :: averagew
  SCALAR, intent(in) :: epshead  
  SCALAR, intent(inout)  :: wcoul0 
  real(DP), intent(in) ::  zc, q0len
  real(DP), intent(out) :: kz

  integer ::  ii
  real(DP) :: gkq(3), length, kxy, gkqxy(3),lengthqk
  real(DP) :: gkqz(3),epsmodel,gamma,alpha,vc,vc_qtozero
  SCALAR  :: integralW
  
  PUSH_SUB(minibzaverage_2d_oneoverq2)
  
!
! Sahar:
! Define Gamma parameter for model epsilon (see Sohrab, PRB 2006)
! Extract the quadratic dependence of 1/epsinv(00)
! 1/epsinv(q;0,0) = 1 + q^2*vc(q)*gamma

!get Vc
  vc_qtozero=((1.0d0 - exp(-q0len*zc))/q0len**2)
! Define Gamma
  gamma = (1.0d0/epshead-1.0d0)/((q0len**2)*vc_qtozero)
!
! Define alpha 
! Set to zero for now
  alpha = 0.0d0

! length of q + G
  lengthqk = sqrt(DOT_PRODUCT(qk,MATMUL(bdot,qk)))

  integral = 0D0
  integralW = 0D0

  do ii = 1, nn
    gkq(:) = qk(:)
    gkq(1:2) = gkq(1:2) + qran(1:2,ii)
    gkqxy(1:2) = gkq(1:2)
    gkqxy(3) = 0D0
    kxy=sqrt(DOT_PRODUCT(gkqxy,MATMUL(bdot,gkqxy)))
    length = DOT_PRODUCT(gkq,MATMUL(bdot,gkq))
    
    ! This is Temporary??
    gkqz(:)=gkq(:)
    gkqz(1)=0D0
    gkqz(2)=0D0
    kz=sqrt(DOT_PRODUCT(gkqz,MATMUL(bdot,gkqz)))
    
! First average v

    integral = integral + (1.0d0+exp(-kxy*zc)* &
      ((kz/kxy)*sin(kz*zc) -  cos(kz*zc))) &
      / length

! Do we also want to average W?

! This is a waste of time if we are not qk=0
    if (lengthqk.lt.TOL_zero.and.averagew) then

! Use model epsilon here 
! Normalize integral by head of epsilon
          
      vc = ((1.0d0 - exp(-kxy*zc))/kxy**2)
      epsmodel=1.0d0 + vc * kxy**2 * gamma*exp(-alpha*kxy)
      integralW = integralW + (vc/epsmodel) 

!      write(6,*)  'USING MODEL EPSILON FOR AVERAGING OF W'
!      write(6,*)  'gamma: ', gamma, 'alpha: ', alpha, 'qk', qk  
!      write(6,*)   'qk', qk  

! No model epsilon here
    endif
  enddo

! Convert integral to Ry

  integral = integral * 8D0 * PI_D / dble(nn) 
  if (lengthqk.lt.TOL_zero.and.averagew) then
    wcoul0 = integralW * 8D0 * PI_D / dble(nn) 
  endif
  if (peinf%inode .eq. 0) then
    write(6,*) "wcoul0, integral, epshead", wcoul0, integral, epshead
  endif
  POP_SUB(minibzaverage_2d_oneoverq2)
     
  return
end subroutine minibzaverage_2d_oneoverq2


!========================================================================


subroutine minibzaverage_3d_oneoverq(nn,bdot,integral,qran,qk)
  integer, intent(in) :: nn
  real(DP), intent(in) :: bdot(3,3), qran(:,:) !< (3,nn)
  real(DP), intent(out) :: integral
  real(DP), intent(in) ::  qk(3)

  integer :: ii
  real(DP) :: gkq(3), length
  
  PUSH_SUB(minibzaverage_3d_oneoverq)
  
  integral = 0D0
  
  do ii = 1, nn
    
    gkq(:) = qk(:) + qran(:,ii)
    length = DOT_PRODUCT(gkq,MATMUL(bdot,gkq))
    length = sqrt(length)
    integral = integral + 1D0/length
    
  enddo
  
  integral = integral * 8D0 * PI_D / dble(nn)
  
  POP_SUB(minibzaverage_3d_oneoverq)

  return
end subroutine minibzaverage_3d_oneoverq

!===========================================================================

subroutine minibzaverage_1d(gvec,N_k,bdot,integral,iparallel,qk,epshead,q0len,averagew,wcoul0)
  type (gspace), intent(in) :: gvec
  integer, intent(in) :: N_k ! number of k-points
  integer, intent(in) :: iparallel
  real(DP), intent(in) :: bdot(3,3),qk(3)
  real(DP), intent(out) :: integral
  logical, intent(in) :: averagew
  real(DP), intent(in) ::  q0len
  SCALAR, intent(in) ::  epshead
  SCALAR, intent(inout) ::  wcoul0

  real(DP) :: integralTemp
  SCALAR :: wcoul0temp
  real(DP) :: epsmodel,gamma,vc_qtozero

  logical :: first_minibz
  integer :: i, j, i1, i2, l1, l2, iline
  real(DP) :: sum_vt, adot(3,3), rr(3), tt(3), rx, ry
  real(DP) :: gpq_xy(2), gpq_z, r_len, t_len, scale, xline

  integer, parameter :: nline = 1000 ! Number of points in 1-D integral

  PUSH_SUB(minibzaverage_1d)

  integral = 0.0d0

  first_minibz = all(abs(qk(1:3)) .lt. Tol_Zero)

  rr = 0.0d0
  tt = 0.0d0

  call invert_matrix(bdot, adot)
  adot = adot * 4.d0 * PI_D * PI_D
  
  do i=1,2
    do j=1,2
      adot(i,j)=adot(i,j)/(dble(gvec%FFTgrid(i)) * dble(gvec%FFTgrid(j)))
    enddo
  enddo

  scale = adot(1,1)*adot(2,2) - adot(1,2)*adot(2,1)
  scale = 4.d0 * sqrt(scale)

! Compute parameters of epsilon model

  if (first_minibz .and. averagew) then
    vc_qtozero = 0D0

    do i2 = 1, gvec%FFTgrid(2)
      rr(2) = dble(i2-1) + trunc_shift(2)
      do i1 = 1, gvec%FFTgrid(1)
        rr(1) = dble(i1-1) + trunc_shift(1)

        r_len = INF

        do l2 = -ncell+1, ncell
          tt(2) = rr(2) - dble(l2 * gvec%FFTgrid(2))
          do l1 = -ncell+1, ncell
            tt(1) = rr(1) - dble(l1 * gvec%FFTgrid(1))   
            t_len = dot_product(tt,matmul(adot,tt))
            if (t_len < r_len) then
              r_len = t_len
            endif
          enddo ! l1
        enddo ! l2

        r_len = sqrt(r_len)
        vc_qtozero = vc_qtozero + dbesk0(q0len * r_len)

      enddo ! i1
    enddo ! i2

    vc_qtozero = vc_qtozero * scale
    gamma = ((1/epshead)-1.0D0) / (q0len**2 * vc_qtozero)
  endif

! Compute integral along z direction of minibz

  do iline = 1, nline
    if (iparallel .eq. 1 .and. mod(iline-1,peinf%npes) .ne. peinf%inode) cycle

    xline = ((dble(iline) - 0.5d0) / dble(nline) - 0.5d0) / dble(N_k)
    gpq_z = abs(qk(3)+xline)*sqrt(bdot(3,3))

    sum_vt = 0D0

    do i2 = 1, gvec%FFTgrid(2)
      rr(2) = dble(i2-1) + trunc_shift(2)

      do i1 = 1, gvec%FFTgrid(1)
        rr(1) = dble(i1-1) + trunc_shift(1)

        r_len = INF

        do l2 = -ncell+1, ncell
          tt(2) = rr(2) - dble(l2 * gvec%FFTgrid(2))
          do l1 = -ncell+1, ncell
            tt(1) = rr(1) - dble(l1 * gvec%FFTgrid(1))   
            t_len = dot_product(tt,matmul(adot,tt))
            if (t_len < r_len) then
              r_len = t_len
              rx = tt(1)
              ry = tt(2)
            endif
          enddo ! l1
        enddo ! l2

        r_len = sqrt(r_len)
        rx = rx/dble(gvec%FFTgrid(1))
        ry = ry/dble(gvec%FFTgrid(2))

        gpq_xy(1:2) = qk(1:2)
        sum_vt = sum_vt + dbesk0(gpq_z * r_len) * &
         cos(2.0d0 * PI_D * (gpq_xy(1)*rx + gpq_xy(2)*ry))

      enddo ! i1
    enddo ! i2

    sum_vt = sum_vt * scale
    integral = integral + sum_vt

    if (first_minibz .and. averagew) then
      epsmodel = 1.0D0 + gamma * gpq_z**2 * sum_vt
      wcoul0 = wcoul0 + (sum_vt / epsmodel) 
    endif

  enddo

  integral = integral / dble(nline)
  if (first_minibz .and. averagew) then
    wcoul0 = wcoul0 / dble(nline)
  endif

  if (iparallel .eq. 1) then
#ifdef MPI
    integralTemp=integral
    wcoul0Temp=wcoul0
    call MPI_Allreduce(integralTemp,integral,1,MPI_REAL_DP, &
     MPI_SUM,MPI_COMM_WORLD,mpierr)
    if (first_minibz .and. averagew) then
      call MPI_Allreduce(wcoul0Temp,wcoul0,1,MPI_SCALAR, &
       MPI_SUM,MPI_COMM_WORLD,mpierr)
    endif
#endif

    if (peinf%verb_debug .and. peinf%inode==0) then
      write(6,'(3x,"vcoul =",e20.12,1x,"wcoul0 =",e20.12)') integral, wcoul0
    endif
  endif

  POP_SUB(minibzaverage_1d)
  
  return
end subroutine minibzaverage_1d

! For modified coulomb interaction
subroutine minibzaverage_3d_oneoverq2_mod(nn,bdot,integral,qran,qk,coulomb_mod)
  integer, intent(in) :: nn
  real(DP), intent(in) :: bdot(3,3), qran(:,:) !< (3, nn)
  real(DP), intent(out) :: integral
  real(DP), intent(in) ::  qk(3)
  type(coulomb_modifier_t), intent(in)  :: coulomb_mod

  real(DP) :: gkq(3), screeninv, temp_exp, length 
  integer :: ii

  PUSH_SUB(minibzaverage_3d_oneoverq2_mod)
  
  integral = 0D0
  screeninv = 1.0D0/(4.0D0 * coulomb_mod%screening_length *coulomb_mod%screening_length)
  ! convert screening_length from A^{-1} to Bohr
  screeninv = screeninv/(BOHR*BOHR)
  
  do ii = 1, nn
    
    gkq(:) = qk(:) + qran(:,ii)
    length = DOT_PRODUCT(gkq,MATMUL(bdot,gkq))
    temp_exp = exp(-screeninv*length)
    integral = integral + 1D0/length*(temp_exp*coulomb_mod%long_range_frac_fock + &
                (1.0D0 - temp_exp)*coulomb_mod%short_range_frac_fock)
    
  enddo
  
  integral = integral * 8D0 * PI_D / dble(nn) 

  POP_SUB(minibzaverage_3d_oneoverq2_mod)
  
  return
end subroutine minibzaverage_3d_oneoverq2_mod

end module minibzaverage_m
