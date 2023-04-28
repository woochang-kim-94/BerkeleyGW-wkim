!=================================================================================
!
! Routines:
!
! (1) gx_sum()         Originally By JRD       Last Modified 4/1/2012 (JRD)
!
!  Sum over G vectors for a given (v,c,v',c') set for the exchange contribution.
!
!  This routine scales as N^5 with number of atoms.
!
!=================================================================================

#include "f_defs.h"

module gx_sum_m

  use global_m
  use blas_m
  implicit none

  public :: gx_sum_TDA, gx_sum_extended

  private

  contains

    !> Sum G vectors to get the exchange kernel block.
    !! Works only with TDA kernel calculations.
    !! Note that, while we should compute
    !! \sum_G M_cv(G) v(G+Q) [M_cpvp(G)]^*,
    !! for historical reasons we end up evaluating this via
    !! \sum_G [M_vc(G)]^* v(G-Q) M_cpvp(G)
    subroutine gx_sum_TDA(xct, invband, incband, vcoul, mvc, mvpcp, &
                     bsex, ivp_in, icp_in, ikp, iv_in, ic_in, ik)
      type (xctinfo), intent(in) :: xct
      integer, intent(in) :: invband,incband
      integer, intent(in) :: ivp_in,icp_in,ikp,iv_in,ic_in,ik
      real(DP), intent(in) :: vcoul(:) !< Modified Coulomb potential (vbar).
      SCALAR, intent(in) :: mvc(:,:,:,:) !< (ng, invband, incband, nspin)
      SCALAR, intent(inout) :: mvpcp(:,:,:,:) !< (ng, invband, incband, nspin)
      SCALAR, intent(inout) :: bsex(:,:,:) !< (leadingdim, nspin, nspin)

      SCALAR, allocatable :: outtemp(:,:,:,:)
      integer :: iscp, isc, iv, ivp, ic, icp, iit

      PUSH_SUB(gx_sum_TDA)

      SAFE_ALLOCATE(outtemp,(invband,incband,invband,incband))
      ! FHJ: Note: we are using the modified Coulomb potential here, where vbar(G=0)=0.
      do iscp=1,xct%nspin
        do icp=1,incband
          do ivp=1,invband
            mvpcp(:,ivp,icp,iscp) = vcoul(:) * mvpcp(:,ivp,icp,iscp)
          enddo
        enddo

        do isc=1,xct%nspin
          ! JRD: Sum over all G
          call X(gemm)('c','n',invband*incband,invband*incband,xct%ng, &
           ONE,mvc(:,:,:,isc),xct%ng,mvpcp(:,:,:,iscp),xct%ng, &
           ZERO,outtemp,invband*incband)

          do icp=1,incband
            do ivp=1,invband
              do ic=1,incband
                do iv=1,invband
                  if (xct%icpar .eq. 0) then
                    iit = peinf%wown(1,1,ikp,1,1,ik) + xct%nvb_co*xct%nvb_co*xct%ncb_co*(icp-1) &
                     + xct%nvb_co*xct%nvb_co*(ic-1) + xct%nvb_co*(ivp-1) + iv -1
                  else if (xct%ivpar .eq. 0) then
                    iit = peinf%wown(1,icp_in,ikp,1,ic_in,ik) + xct%nvb_co*(ivp-1) + iv -1
                  else
                    iit = peinf%wown(ivp_in,icp_in,ikp,iv_in,ic_in,ik)
                  endif
                  bsex(iit,isc,iscp) = bsex(iit,isc,iscp) + outtemp(iv,ic,ivp,icp)
                enddo !iv
              enddo !ic
            enddo !ivp
          enddo !icp
        enddo !isc
      enddo !iscp

      SAFE_DEALLOCATE(outtemp)

      POP_SUB(gx_sum_TDA)

    end subroutine gx_sum_TDA

    !> Sum G vectors to get the exchange kernel blocks.
    !! Works with TDA and extended kernel calculations.
    subroutine gx_sum_extended(xct, ofs1, ofs2, ofs1p, ofs2p, n1, n2, n1p, n2p, &
      m12, v_m1p2p, bsex, i1p_in, i2p_in, ikp, i1_in, i2_in, ik)
      type (xctinfo), intent(in) :: xct
      !> offset (i.e., add this number to map a local index to the global band index)
      integer, intent(in) :: ofs1, ofs2, ofs1p, ofs2p
      !> number of bands for each wfn
      integer, intent(in) :: n1, n2, n1p, n2p
      integer, intent(in) :: i1p_in, i2p_in, ikp, i1_in, i2_in, ik
      SCALAR, intent(in) :: m12(:,:,:,:) !< (ng, in1band, in2band, nspin)
      SCALAR, intent(in) :: v_m1p2p(:,:,:,:) !< (ng, in1band, in2band, nspin)
      SCALAR, intent(inout) :: bsex(:,:,:) !< (leadingdim, nspin, nspin)

      SCALAR, allocatable :: outtemp(:,:,:,:)
      integer :: iscp, isc, i1, i1p, i2, i2p, iit
      integer :: gi1, gi1p, gi2, gi2p

      PUSH_SUB(gx_sum_extended)

      SAFE_ALLOCATE(outtemp,(n1,n2,n1p,n2p))

      do iscp=1,xct%nspin
        do isc=1,xct%nspin
          ! JRD: Sum over all G
          call X(gemm)('c','n',n1*n2,n1p*n2p,xct%ng, &
           ONE,m12(:,:,:,isc),xct%ng,v_m1p2p(:,:,:,iscp),xct%ng, &
           ZERO,outtemp,n1*n2)

          do i2p=1,n2p
            gi2p = i2p + ofs2p
            do i1p=1,n1p
              gi1p = i1p + ofs1p
              do i2=1,n2
                gi2 = i2 + ofs2
                do i1=1,n1
                  gi1 = i1 + ofs1
                  if (xct%icpar .eq. 0) then
                    iit = peinf%wown(1,1,ikp,1,1,ik) + xct%n1b_co*xct%n1b_co*xct%n2b_co*(gi2p-1) &
                     + xct%n1b_co*xct%n1b_co*(gi2-1) + xct%n1b_co*(gi1p-1) + gi1 -1
                  else if (xct%ivpar .eq. 0) then
                    iit = peinf%wown(1,i2p_in,ikp,1,i2_in,ik) + xct%n1b_co*(gi1p-1) + gi1 -1
                  else
                    iit = peinf%wown(i1p_in,i2p_in,ikp,i1_in,i2_in,ik)
                  endif
                  bsex(iit,isc,iscp) = bsex(iit,isc,iscp) + outtemp(i1,i2,i1p,i2p)
                enddo !i1
              enddo !i2
            enddo !i1p
          enddo !i2p
        enddo !isc
      enddo !iscp

      SAFE_DEALLOCATE(outtemp)

      POP_SUB(gx_sum_extended)

    end subroutine gx_sum_extended

end module gx_sum_m
