!==========================================================================
!
! lin_denominator_m   Originally By Liang Z. Tan         2/8/2013
!
!    computes the integral of
!
!      Theta(Ef-Ev(k))Theta(Ec(k+q)-Ef)/(Ev(k)-Ec(k+q)+w+i0+)
!
!    where the range of integration over k is the mini-BZ and 
!    Ec(k+q) and Ev(k) are approximated as linear functions of k 
!    (i.e. planes).
!    Currently implemented for 2D. 
!    z-direction is the out-of-plane direction.
!
!    input: k, energies and velocities at k, k+q, frequency, Ef.
!    output: (complex) value of the integral.
!
!==========================================================================



#include "f_defs.h"

module lin_denominator_m

  use global_m
  use sort_m
  implicit none

  !> define boundaries as a useful type
  !! x = slope * y + intercept
  type boundary 
    real(DP) :: slope, intercept
    real(DP) :: sgn 
  end type boundary

  private
  public :: integrate_mbz
   
contains
 
  !> this integrates log|ax+b| over x in [0,1].
  real(DP) function integratelog(a,b)
    real(DP), intent(in) :: a,b

    PUSH_SUB(integratelog)

    if(abs(a) > TOL_Zero .and. abs(b)>TOL_Zero) then
      if(abs(a+b) > TOL_Zero) then
        integratelog = -1.0d0 -(b/a)*log(abs(b))+((a+b)/a)*log(abs(a+b))
      else
        integratelog = -1.0d0 -(b/a)*log(abs(b))
      endif
    elseif(abs(a) > TOL_Zero .and. abs(b)<TOL_Zero) then
      if(abs(a+b) > TOL_Zero) then
        integratelog = -1.0d0 +((a+b)/a)*log(abs(a+b))
      else
        integratelog = -1.0d0 
      endif
    else
      integratelog = log(abs(b))
    endif
    !write(6,*) "integratelog: a= ",a," b= ",b," ans= ",integratelog

    POP_SUB(integratelog)

  end function integratelog

  !> evaluate a linear function at a point
  real(DP) function evalb(b,y)
    type(boundary), intent(in) :: b
    real(DP), intent(in) :: y

    PUSH_SUB(evalb)

    evalb = b%slope * y + b%intercept

    POP_SUB(evalb)
  end function evalb

  !> intersection point of two lines
  real(DP) function critfn(b1,b2)
    type(boundary), intent(in) :: b1, b2

    PUSH_SUB(critfn)
    if(b1%slope /= b2%slope) then
      critfn = (b2%intercept - b1%intercept)/(b1%slope - b2%slope)
    else
      critfn = 100
    end if

    POP_SUB(critfn)
  end function critfn


  !> imaginary part of the integral 1/(ux*kx+uy*ky+1+J*eta)
  !! over a trapezoid in 2D
  real(DP) function trapzimagpart(ux,uy,y0,y1,bdl,bdu)
    real(DP),intent(in) :: ux,uy,y0,y1
    type(boundary),intent(in) :: bdl, bdu

    real(DP) :: yl,yu,x1,x0,x0l,x0u,x1l,x1u,pi,yc
    type(boundary) :: bdsing

    PUSH_SUB(trapzimagpart)

    pi = 4.d0*datan(1.d0)

    if(evalb(bdl,(y0+y1)/2) > evalb(bdu,(y0+y1)/2)) then
      trapzimagpart = 0d0 
    elseif(abs(ux)<TOL_Zero .and. abs(uy)<TOL_Zero) then
      trapzimagpart = 0d0
    elseif (abs(ux)<TOL_Zero) then
      yc = -1.0d0/uy
      if( yc<y1 .and. yc>y0) then
        trapzimagpart = -pi*(evalb(bdu,yc)-evalb(bdl,yc))/abs(uy)
      else
        trapzimagpart = 0d0
      endif
    else
      bdsing%slope = -uy/ux
      bdsing%intercept = -1.0d0/ux
      yu = critfn(bdsing,bdu)
      yl = critfn(bdsing,bdl)
      x1 = evalb(bdsing,y1)
      x0 = evalb(bdsing,y0)
      x0l= evalb(bdl,y0)
      x0u= evalb(bdu,y0)
      x1l= evalb(bdl,y1)
      x1u= evalb(bdu,y1)
      if ((x0<x0l .and. x1<x1l) .or. (x0>x0u .and. x1>x1u)) then 
        trapzimagpart=0d0
      elseif(x0>x0l .and. x0<x0u .and. x1>x1l .and. x1<x1u) then
        trapzimagpart= -pi*abs(y1-y0)/abs(ux)
      elseif((x0<x0l .and. x1>x1u) .or. (x0>x0u .and. x1<x1l))then
        trapzimagpart = -pi*abs(yl-yu)/abs(ux)
      elseif(x0<x0l .and. x1>x1l .and. x1<x1u) then
        trapzimagpart = -pi*abs(y1-yl)/abs(ux)
      elseif(x0>x0u .and. x1>x1l .and. x1<x1u) then 
        trapzimagpart = -pi*abs(y1-yu)/abs(ux)
      elseif(x1<x1l .and. x0>x0l .and. x0<x0u) then 
        trapzimagpart = -pi*abs(yl-y0)/abs(ux)
      !elseif(x1>x1u .and. x0>x0l .and. x0<x0u) then
      else
        trapzimagpart = -pi*abs(yu-y0)/abs(ux)
      endif
    endif

    POP_SUB(trapzimagpart)
  end function trapzimagpart

  !> real part of the integral 1/(ux*kx+uy*ky+1+J*eta)
  !! over a trapezoid in 2D
  real(DP) function trapzrealpart(ux,uy,y0,y1,bdl,bdu)
    real(DP),intent(in) :: ux,uy,y0,y1
    type(boundary),intent(in) :: bdl, bdu

    real(DP) :: a,b,part1, part2

    PUSH_SUB(trapzrealpart)

    if(evalb(bdl,(y0+y1)/2) > evalb(bdu,(y0+y1)/2)) then
      trapzrealpart = 0d0 
    elseif( abs(ux)<TOL_Small .and. abs(uy)<TOL_Small) then
      trapzrealpart = 0.5*(y1-y0)&
      *(evalb(bdu,y0)-evalb(bdl,y0)+evalb(bdu,y1)-evalb(bdl,y1))
    elseif(abs(ux)<TOL_Small) then
      a = (bdu%slope-bdl%slope)/uy
      b = bdu%intercept-bdl%intercept - a
      part1 = (y1-y0)*a
      part2 = (b/uy)*(log(abs(uy*y1+1.0d0))-log(abs(uy*y0+1.0d0)))
      trapzrealpart = part1 + part2
    else
      a = (ux*bdu%slope+uy)*(y1-y0)
      b = (ux*bdu%slope+uy)*y0+ux*bdu%intercept+1.0d0
      part1 = (y1-y0)/ux * integratelog(a,b)
      a = (ux*bdl%slope+uy)*(y1-y0)
      b = (ux*bdl%slope+uy)*y0+ux*bdl%intercept+1.0d0
      part2 = (y1-y0)/ux * integratelog(a,b)
      trapzrealpart = part1 - part2
    endif

    POP_SUB(trapzrealpart)
  end function trapzrealpart

  !> picks the relevant boundaries from a list
  subroutine actb(bs,y,b0,b1)
    type(boundary), intent(in) :: bs(4)
    type(boundary), intent(out) :: b0,b1

    type(boundary) :: b,bs0(4), bs1(4)
    integer :: i,nbs0,nbs1, xidx
    real(DP), intent(in) :: y
    real(DP) :: x,xtry

    PUSH_SUB(actb)

    nbs0=0
    nbs1=0
    do i = 1,4
      b = bs(i)
      if(b%sgn>0) then
        nbs0 = nbs0 +1
        bs0(nbs0) = b
      else
        nbs1 = nbs1 +1
        bs1(nbs1) = b
      endif
    enddo

    x = -100
    do i = 1,nbs0
      xtry = evalb(bs0(i),y)
      if(xtry> x) then
        x = xtry
        xidx = i
      endif
    enddo
    b0 = bs0(xidx) 
   
    x = 100
    do i = 1,nbs1
      xtry = evalb(bs1(i),y)
      if(xtry<x) then
        x = xtry
        xidx = i
      endif
    enddo
    b1 = bs1(xidx)
    POP_SUB(actb) 
  end subroutine actb

  !> computes the interval (ymin,ymax) from a set of n linear inequalities
  !! that y satisfies.
  !! Each inequality is modelled as a boundary a*y + c > 0 or a*y + c < 0
  !! where a = slope, c = intercept, > : sgn=1.0, < : sgn=-1.0
  !! Implicit bounds (0,1) assumed.
  !! If there is no solution, a zero-length interval is returned.
  subroutine intervals(n,ineqs,ymin,ymax)
    integer,intent(in) :: n
    type(boundary),intent(in) :: ineqs(:)
    real(DP),intent(out) :: ymin,ymax

    integer :: ii
    real(DP) :: a,c,s,y2

    PUSH_SUB(intervals)

    ymin = 0.0d0
    ymax = 1.0d0

    do ii = 1,n
      a=ineqs(ii)%slope
      c=ineqs(ii)%intercept
      s=ineqs(ii)%sgn
      if(abs(a)<TOL_Small)then
        if(s*c<0.0d0) then
          ymax = ymin
          exit
        endif
      else
        if(s*a>0.0) then
          y2 = -c/a
          if(y2>ymin)then 
            ymin = y2
          endif
          if(ymin>ymax) then 
            ymin = ymax
            exit
          endif
        else
          y2 = -c/a
          if(y2<ymax)then 
            ymax = y2
          endif
          if(ymax<ymin) then 
            ymax = ymin
            exit
          endif
        endif
      endif
    enddo

    POP_SUB(intervals)
  end subroutine intervals

  !> exact value of the integral Theta(a.r+c)Theta(b.r+c)/(u.r+1) 
  !! over the unit square in 2D.
  complex(DP) function exactmodel1(ax,ay,c,bx,by,d,ux,uy)
    real(DP), intent(in) :: ax,ay,c,bx,by,d,ux,uy 

    type(boundary) :: bd0, bd1, bd2, bd3, bds(4)
    type(boundary) :: bl,bu,bls(7),bus(7),intv1,intv2
    real(DP) :: critpts0(5), y(7), yval, accu_real, accu_imag,ymin,ymax
    complex(DP) :: accu
    integer ::  critpts0_idx(5)
    integer :: num_crit
    integer :: i 

    PUSH_SUB(exactmodel1)

    !make boundaries
    bd0%slope = 0d0
    bd0%intercept = 0d0
    bd0%sgn = 1.0d0
    bd1%slope = 0d0
    bd1%intercept = 1.0d0
    bd1%sgn = -1.0d0

    if(abs(ax)>TOL_Small .and. abs(bx)>TOL_Small) then
      bd2%slope = -ay/ax
      bd2%intercept = -c/ax
      bd2%sgn = sign(1.0d0,ax)
      bd3%slope = -by/bx
      bd3%intercept = -d/bx
      bd3%sgn = sign(1.0d0,bx)

      ymin = 0.0d0
      ymax = 1.0d0
    elseif(abs(ax)<TOL_Small .and. abs(bx)>TOL_Small) then
      bd2%slope = 1.0d0 !< dummy boundary
      bd2%intercept = 100.0 !<dummy boundary
      bd2%sgn = -1.0d0 !<dummy boundary

      bd3%slope = -by/bx
      bd3%intercept = -d/bx
      bd3%sgn = sign(1.0d0,bx)

      intv1 = boundary(ay,c,1.0d0)
      call intervals(1,(/intv1/),ymin,ymax)
      !write(6,*) "ymin= ",ymin,"ymax= ",ymax

 
    elseif(abs(ax)>TOL_Small .and. abs(bx)<TOL_Small) then
      bd2%slope = -ay/ax
      bd2%intercept = -c/ax
      bd2%sgn = sign(1.0d0,ax)
      bd3%slope = 1.0d0 !< dummy boundary
      bd3%intercept = 100.0d0 !<dummy boundary
      bd3%sgn = -1.0d0 !<dummy boundary

      intv1 = boundary(by,d,1.0d0)
      call intervals(1,(/intv1/),ymin,ymax)
 
    else !both ax and bx small
      bd2%slope = 0.9d0 !< dummy boundary
      bd2%intercept = 101.0d0 !<dummy boundary
      bd2%sgn = -1.0d0 !<dummy boundary
      bd3%slope = 1.0d0 !< dummy boundary
      bd3%intercept = 100.0d0 !<dummy boundary
      bd3%sgn = -1.0d0 !<dummy boundary

      intv1 = boundary(ay,c,1.0d0)
      intv2 = boundary(by,d,1.0d0)
      call intervals(2,(/intv1,intv2/),ymin,ymax)
      !write(6,*) "ymin= ",ymin,"ymax= ",ymax
    endif

    bds(1) = bd0
    bds(2) = bd1
    bds(3) = bd2
    bds(4) = bd3
     
    !first split integration region into trapezoids
    !start by finding where the boundaries intersect
    critpts0(1) = critfn(bd0,bd2)
    critpts0(2) = critfn(bd0,bd3)
    critpts0(3) = critfn(bd1,bd2)
    critpts0(4) = critfn(bd1,bd3)
    critpts0(5) = critfn(bd2,bd3)
  
    call sortrx(5,critpts0,critpts0_idx)
    num_crit = 2
    y(1) = ymin
    do i = 1,5
      yval = critpts0(critpts0_idx(i))
      if(yval>ymin .and. yval<ymax) then
        y(num_crit) = yval
        num_crit = num_crit + 1
      endif
    enddo
    y(num_crit) = ymax
    

    !make the list of lower and upper boundaries
    do i = 1,num_crit-1
      call actb(bds,(y(i)+y(i+1))/2,bl,bu)
      bls(i) = bl
      bus(i) = bu
    enddo

    !sum over all the trapezoids
    accu = (0,0)
    do i = 1,num_crit-1
      !write(6,*) "ux= ",ux,"uy= ",uy,"yi= ",y(i),"yi1= ",y(i+1),&
      !"bls= ",bls(i),"bus= ",bus(i)
      accu_real = trapzrealpart(ux,uy,y(i),y(i+1),bls(i),bus(i))
      accu_imag = trapzimagpart(ux,uy,y(i),y(i+1),bls(i),bus(i))
      accu = accu + CMPLX(accu_real,accu_imag)
    enddo

    exactmodel1 = accu

    !write(6,*) "exactmodel1= ",exactmodel1

    POP_SUB(exactmodel1)
  end function exactmodel1

  !> compute integral in miniBZ from energies and velocities
  !! NOTE: currently for 2D only. z-direction is the out-of-plane direction.
  complex(DP) function integrate_mbz(kpts,k,ek,ekq,vk,vkq,w,ef)
    type(kpoints), intent(in) :: kpts
    real(DP), intent(in) :: k(3)
    real(DP), intent(in):: ek, ekq, vk(2), vkq(2), w,ef 

    real(DP) :: kx,ky,dkx,dky,kxmin,kymin,vkx,vky,vkqx,vkqy, &
                cc,ux,uy,ax,ay,c,bx,by,d
    complex(DP) :: accu
    integer :: ii, ix, iy

    PUSH_SUB(integrate_mbz)

    ! need to know which is the out of plane direction
    ix = 1
    iy = 2
    do ii = 1,3
      if(kpts%kgrid(ii) == 1) then
        if(ii <= ix) then
          ix = ix +1
        endif
        if(ii <= iy) then
          iy = iy +1
        endif
      endif 
    enddo

    dkx = 1.0d0/real(kpts%kgrid(ix)) 
    dky = 1.0d0/real(kpts%kgrid(iy)) 
    kx = k(ix)
    ky = k(iy)
    vkx = vk(1)
    vky = vk(2)
    vkqx = vkq(1)
    vkqy = vkq(2)
    kxmin = kx-dkx/2.0
    kymin = ky-dky/2.0

    cc = (vkx-vkqx)*kxmin + (vky-vkqy)*kymin + &
         (ek - ekq +(vkqx-vkx)*kx +(vkqy-vky)*ky + w)
    ux = (vkx-vkqx)*dkx/cc
    uy = (vky-vkqy)*dky/cc

    
    accu = (0,0)
    !k occ, kq unocc
    ax = vkx*dkx
    ay = vky*dky
    c = vkx*kxmin+vky*kymin-vkx*kx-vky*ky+ek-ef
    bx = -vkqx*dkx
    by = -vkqy*dky
    d = -vkqx*kxmin-vkqy*kymin+vkqx*kx+vkqy*ky+(ef-ekq)
    accu = accu -  exactmodel1(ax,ay,c,bx,by,d,ux,uy)/cc

    !if(dble(accu)> 0.0) then
    !  write(6,*) "accu1= ",accu,"vkx= ",vkx,"vky= ",vky,&
    !  "vkqx= ",vkqx,"vkqy= ",vkqy,&
    !  "ax",ax,"ay",ay,"bx",bx,"by",by,&
    !  "c= ",c,"d= ",d,"cc= ",cc
    !endif

    ! k occ, kq unocc
    ax = -vkx*dkx
    ay = -vky*dky
    c = -1.0d0*(vkx*kxmin+vky*kymin-vkx*kx-vky*ky+ek-ef)
    bx = vkqx*dkx
    by = vkqy*dky
    d = -1.0d0*(-vkqx*kxmin-vkqy*kymin+vkqx*kx+vkqy*ky+(ef-ekq))
    accu = accu +  exactmodel1(ax,ay,c,bx,by,d,ux,uy)/cc

    !if(dble(accu)> 0.0) then
    !  write(6,*) "accu2= ",accu,"vkx= ",vkx,"vky= ",vky,&
    !  "vkqx= ",vkqx,"vkqy= ",vkqy,&
    !  "ax",ax,"ay",ay,"bx",bx,"by",by,&
    !  "c= ",c,"d= ",d,"cc= ",cc
    !endif

    integrate_mbz = accu
    POP_SUB(integrate_mbz)
  end function integrate_mbz

end module lin_denominator_m
