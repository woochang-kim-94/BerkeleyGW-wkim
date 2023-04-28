!============================================================================
!
! Routines:
!
! (1) gmap          Originally by ?          Last Modified: 1/24/2011 (gsm)
!
!     Find the index array and phases needed to map G-vectors for a
!     wavefunction at one k-point to the corresponding wavefunction
!     at a symmetry-related k-point. This routine is flavorless.
!     In the real version, pass dphase. In complex, pass zphase.
!
!============================================================================

#include "f_defs.h"

module gmap_m

  use global_m
  use misc_m

  implicit none

  private

  public :: gmap

  interface gmap
    module procedure dgmap, zgmap
  end interface gmap

contains

subroutine dgmap(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, phase, die_outside_sphere)
  type (gspace), intent(in) :: gvec         !< uses gvec%components(:,isrtc(1:ngk)), and index_vec
  type (symmetry), intent(in) :: syms       !< uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
  integer, intent(in) :: ngk                !< number of g-vector entries in a wavefunction
  integer, intent(in) :: itran              !< index of transformation
  integer, intent(in) :: kgq(3)             !< an umklapp vector (i.e. integer 3-vector)
  integer, intent(in) :: isortc(:)          !< index array for R(q) (1:ngk)
  integer, intent(in) :: isorti(:)          !< inverse index array for q (1:gvec%ng)
  integer, intent(out) :: ind(:)            !< indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
  real(DP), intent(out) :: phase(:)         !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
  logical, intent(in) :: die_outside_sphere !< specifies whether to die if G-vectors are falling outside of the sphere

  PUSH_SUB(dgmap)

  call gmap_base(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, die_outside_sphere, dphase = phase)

  POP_SUB(dgmap)
  return
end subroutine dgmap

!===============================================================================

subroutine zgmap(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, phase, die_outside_sphere)
  type (gspace), intent(in) :: gvec         !< uses gvec%components(:,isrtc(1:ngk)), and index_vec
  type (symmetry), intent(in) :: syms       !< uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
  integer, intent(in) :: ngk                !< number of g-vector entries in a wavefunction
  integer, intent(in) :: itran              !< index of transformation
  integer, intent(in) :: kgq(3)             !< an umklapp vector (i.e. integer 3-vector)
  integer, intent(in) :: isortc(:)          !< index array for R(q) (1:ngk)
  integer, intent(in) :: isorti(:)          !< inverse index array for q (1:gvec%ng)
  integer, intent(out) :: ind(:)            !< indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
  complex(DPC), intent(out) :: phase(:)     !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
  logical, intent(in) :: die_outside_sphere !< specifies whether to die if G-vectors are falling outside of the sphere

  PUSH_SUB(zgmap)

  call gmap_base(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, die_outside_sphere, zphase = phase)

  POP_SUB(zgmap)
  return
end subroutine zgmap

!===============================================================================

subroutine gmap_base(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, die_outside_sphere, dphase, zphase)
  type (gspace), intent(in) :: gvec         !< uses gvec%components(:,isrtc(1:ngk)), and index_vec
  type (symmetry), intent(in) :: syms       !< uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
  integer, intent(in) :: ngk                !< number of g-vector entries in a wavefunction
  integer, intent(in) :: itran              !< index of transformation
  integer, intent(in) :: kgq(3)             !< an umklapp vector (i.e. integer 3-vector)
  integer, intent(in) :: isortc(:)          !< index array for R(q) (1:ngk)
  integer, intent(in) :: isorti(:)          !< inverse index array for q (1:gvec%ng)
  integer, intent(out) :: ind(:)            !< indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
  logical, intent(in) :: die_outside_sphere !< specifies whether to die if G-vectors are falling outside of the sphere
  real(DP),     optional, intent(out) :: dphase(:) !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
  complex(DPC), optional, intent(out) :: zphase(:) !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)

  integer :: ig, kd(3), kadd, kgrad, kgrad1
  integer :: kg(3), kgr(3), mtrxi(3,3), nout, nin
  real(DP) :: fi
  
  PUSH_SUB(gmap_base)

  if(present(dphase) .and. present(zphase)) then
    call die("gmap: cannot pass both dphase and zphase")
  else if(.not. present(dphase) .and. .not. present(zphase)) then
    call die("gmap: must pass either dphase or zphase")
  endif

  if(ngk > gvec%ng) call die("gmap: ngk (wfn cutoff) is greater than gvec%ng (rho cutoff)")

  if(ubound(isorti, 1) < gvec%ng)  call die("gmap: isorti size < gvec%ng")
!  if(any(isorti(1:gvec%ng) < 1))       call die("gmap: isorti cannot be less than 1.")
  if(any(isorti(1:gvec%ng) > gvec%ng)) call die("gmap: isorti cannot be greater than gvec%ng.")

  if(ubound(isortc, 1) < ngk)      call die("gmap: isortc size < ngk")
  if(any(isortc(1:ngk) < 1))       call die("gmap: isortc cannot be less than 1.")
  if(any(isortc(1:ngk) > gvec%ng)) call die("gmap: isortc cannot be greater than ng.")

  if(ubound(gvec%index_vec, 1) /= gvec%nFFTgridpts) call die("gmap: gvec%index_vec has wrong size")
  if(any(gvec%index_vec(1:gvec%nFFTgridpts) < 0)) call die("gmap: index_vec cannot be less than 0")
  if(any(gvec%index_vec(1:gvec%nFFTgridpts) > gvec%ng)) call die("gmap: index_vec cannot be greater than ng")

  if(present(dphase)) then
    if(ubound(dphase, 1) < ngk) call die("gmap: dphase size < ngk")
  else
    if(ubound(zphase, 1) < ngk) call die("gmap: zphase size < ngk")
  endif
  if(ubound(ind, 1) < ngk) call die("gmap: ind size < ngk")

! Invert rotation matrix that gives rq

  call invert_matrix_int(syms%mtrx(1:3, 1:3, itran), mtrxi(1:3, 1:3))

! JRD: Temporary Debugging

!      write(6,*) peinf%inode,'itran: ',itran
!      write(6,*) peinf%inode,'mtrxi: '
!      do i = 1, 3
!        write(6,*) peinf%inode,(mtrxi(i,j),j=1,3)
!      enddo
!      write(6,*) peinf%inode,'kgq',kgq

! Loop over g-vectors in function of r(q)

  nout = 0  ! number of waves outside sphere
  nin = 0  ! number of waves inside sphere

  do ig = 1, ngk

! kg = g(ig) + kgq

    kg(1:3) = gvec%components(1:3, isortc(ig)) + kgq(1:3)

! kgr = (r**-1) kg

    kgr(1:3) = MATMUL(mtrxi(1:3, 1:3), kg(1:3))
    
! Compute address of kgr -> kgrad

    kd(1:3) = kgr(1:3) + gvec%FFTgrid(1:3) / 2 + 1
    if (any(kd(1:3) .lt. 1 .or. kd(1:3) .gt. gvec%FFTgrid(1:3))) then
      call die('gmap: kd out of bounds')
    endif
    kadd = ((kd(1) - 1) * gvec%FFTgrid(2) + kd(2) - 1) * gvec%FFTgrid(3) + kd(3)
    kgrad1 = gvec%index_vec(kadd) ! index_vec relate cube and sphere
    if (kgrad1 .lt. 1 .or. kgrad1 .gt. gvec%ng) then
      write(0,*) 'itran = ', itran, 'ig = ', ig, ', kadd = ', kadd, ', kgrad1 = ', kgrad1
      call die('gmap: G-vectors falling outside of the charge-density G-sphere')
    endif
    kgrad = isorti(kgrad1)

! SIB: if kgr is outside the sphere, then increment out counter,
! set its phase to zero, and its ind() entry to the maximum.
    if (kgrad .gt. ngk) then ! outside sphere
      nout = nout + 1
      ind(ig) = ngk
      if(present(zphase)) then
        zphase(ig) = CMPLX(0d0, 0d0)
      else
        dphase(ig) = 0d0
      endif
!    else if (kgrad <  1) then ! inside sphere
!      nin = nin + 1
!      ind(ig) = 1
!      if(present(zphase)) then
!        zphase(ig) = CMPLX(0d0, 0d0)
!      else
!        dphase(ig) = 0d0
!      endif
    else
! SIB:  Otherwise, record the index of kgr (kgrad) into ind(ig)
! and compute the phase = exp(-i*kg.dot.syms%tnp(:,itran))
      ind(ig) = kgrad
      fi = dot_product(dble(kg(:)), syms%tnp(:,itran))

      if(present(zphase)) then
        zphase(ig) = CMPLX(cos(fi), -sin(fi))
      else

! DAS: The imaginary part can be thrown away because it is always zero
! if we have inversion and time-reversal symmetries, and so can use the real version.
! phase = +/- 1. Otherwise the wavefunction would not be normalized.
! c(G) -> c(G) e^iGt with fractional translation.
! c(G) e^iGt = c(-G)* e^-iGt by time-reversal symmetry
!  = c(G)* e^-iGt by inversion symmetry. c(G) = c(G)* since real.
! Therefore e^iGt = e^-iGt. e^iGt is real, and hence 1 or -1.
! Note there is also a global phase e^ikt, but it is just a convention
! and can be safely ignored here.

        dphase(ig) = cos(fi)

        if(abs(abs(dphase(ig)) - 1) .gt. TOL_Small) then
          write(0,'(a,i8,a,f12.8,a)') 'phase(', ig, ') = ', dphase(ig), ' != +/- 1'
          call die("Illegal non-unity phase in gmap, error in fractional translation.")
        endif
        if(abs(sin(fi)) .gt. TOL_Small) then
          write(0,'(a,i8,a,f12.8)') 'Im phase(', ig, ') = ', -sin(fi)
          call die("Illegal complex phase in gmap, error in fractional translation.")
        endif
      endif
    endif
  enddo !end loop over g-vectors (ig)

  if (die_outside_sphere .and. nout .gt. 0) then
    call die('G-vectors are falling outside of the sphere in gmap')
  endif

!  if (die_outside_sphere .and. nin .gt. 0) then
!    call die('G-vectors are falling inside of the sphere in gmap')
!  endif

  POP_SUB(gmap_base)

  return
  end subroutine gmap_base

!===============================================================================

end module gmap_m
