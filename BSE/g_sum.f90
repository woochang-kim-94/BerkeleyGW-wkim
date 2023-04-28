!=================================================================================
!
! Routines:
!
! (1) g_sum()         Originally By JRD       Last Modified 2/11/2012 (FHJ)
!
!  Performs the sum over G for each (v,c,v',c') set.
!
!  This routine scales as N^5. This routine represents the worst 
!  scaling with number of atoms in Kernel code.
!
!  Bow down before it.
!
!=================================================================================

#include "f_defs.h"

module g_sum_m

  use accel_linalg_m
  use algos_kernel_m
  use blas_m
  use global_m

  implicit none

  public :: g_sum_TDA, g_sum_extended

  private

contains

  !> FHJ: Performs the summation over Gp to get the direct kernel term.
  !! This routine is specialized for TDA kernels.
  subroutine g_sum_TDA(xct,invband,incband,temph,tempb,mccp, &
    bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik,&
    tempw,bsedwing)
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: invband,incband,leading_dim
    integer, intent(in) :: ivp_in,icp_in,ikp,iv_in,ic_in,ik
    SCALAR, intent(in), target :: tempb(:,:,:,:), temph(:,:,:), mccp(:,:,:,:)
    SCALAR, intent(inout), target :: bsedhead(:,:,:), bsedbody(:,:,:)
    SCALAR, intent(in), optional, target :: tempw(:,:,:,:)
    SCALAR, intent(inout), optional, target :: bsedwing(:,:,:)

    PUSH_SUB(g_sum_TDA)

    if (present(bsedwing) .and. .not. present(tempw)) then
      call die("Internal error in g_sum_TDA: bsedwing present but tempw not present", &
               only_root_writes = .true.)
    end if

    select case (g_sum_algo)
    case (CPU_ALGO)
      if (present(bsedwing)) then
        call g_sum_TDA_cpu(xct, invband, incband, temph, tempb, mccp, &
                           bsedhead, bsedbody, leading_dim, ivp_in, icp_in, &
                           ikp, iv_in, ic_in, ik, tempw, bsedwing)
      else
        call g_sum_TDA_cpu(xct, invband, incband, temph, tempb, mccp, &
                           bsedhead, bsedbody, leading_dim, ivp_in, icp_in, &
                           ikp, iv_in, ic_in, ik)
      end if
    case (OPENACC_ALGO)
      if (present(bsedwing)) then
        call g_sum_TDA_openacc(xct, invband, incband, temph, tempb, mccp, &
                           bsedhead, bsedbody, leading_dim, ivp_in, icp_in, &
                           ikp, iv_in, ic_in, ik, tempw, bsedwing)
      else
        call g_sum_TDA_openacc(xct, invband, incband, temph, tempb, mccp, &
                           bsedhead, bsedbody, leading_dim, ivp_in, icp_in, &
                           ikp, iv_in, ic_in, ik)
      end if
    case (OMP_TARGET_ALGO)
      if (present(bsedwing)) then
        call g_sum_TDA_omp_target(xct, invband, incband, temph, tempb, mccp, &
                           bsedhead, bsedbody, leading_dim, ivp_in, icp_in, &
                           ikp, iv_in, ic_in, ik, tempw, bsedwing)
      else
        call g_sum_TDA_omp_target(xct, invband, incband, temph, tempb, mccp, &
                           bsedhead, bsedbody, leading_dim, ivp_in, icp_in, &
                           ikp, iv_in, ic_in, ik)
      end if
    case default
      call die("Invald algorithm for g_sum_TDA", &
               only_root_writes = .true.)
    end select

    POP_SUB(g_sum_TDA)

  end subroutine g_sum_TDA

  subroutine g_sum_TDA_cpu(xct,invband,incband,temph,tempb,mccp, &
    bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik,&
    tempw,bsedwing)
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: invband,incband,leading_dim
    integer, intent(in) :: ivp_in,icp_in,ikp,iv_in,ic_in,ik
    SCALAR, intent(in) :: tempb(:,:,:,:), temph(:,:,:), mccp(:,:,:,:)
    SCALAR, intent(inout) :: bsedhead(:,:,:), bsedbody(:,:,:)
    SCALAR, intent(in), optional :: tempw(:,:,:,:)
    SCALAR, intent(inout), optional :: bsedwing(:,:,:)

    integer :: isv, isc, iv, ivp, ic, icp, iit, iitend, iitbeG

    PUSH_SUB(g_sum_TDA_cpu)

    do isv=1,xct%nspin
      do isc=1,xct%nspin

        if (xct%icpar .eq. 0) then
          iit = peinf%wown(1,1,ikp,1,1,ik) - 1
        else if (xct%ivpar .eq. 0) then 
          iit = peinf%wown(1,icp_in,ikp,1,ic_in,ik) - 1
        else
          iit = peinf%wown(ivp_in,icp_in,ikp,iv_in,ic_in,ik) - 1
        endif

        iitbeg = iit + 1
        iitend = iit + invband*invband*incband*incband

        if (incband .gt. 1) then

          call X(gemm)('t','n',invband*invband,incband*incband,xct%ng, &
            ONE,tempb(:,:,:,isv),xct%ng,mccp(:,:,:,isc),xct%ng, &
            ONE,bsedbody(iitbeg:iitend,isc,isv),invband*invband)

          if (present(bsedwing)) &
            call X(gemm)('t','n',invband*invband,incband*incband,xct%ng, &
            ONE,tempw(:,:,:,isv),xct%ng,mccp(:,:,:,isc),xct%ng, &
            ONE,bsedwing(iitbeg:iitend,isc,isv),invband*invband)

        else

          call X(gemv)('t',xct%ng,invband*invband,ONE,tempb(:,:,:,isv),xct%ng,mccp(:,:,:,isc), &
            1,ZERO,bsedbody(iitbeg:iitend,isc,isv),1)

          if (present(bsedwing)) &
            call X(gemv)('t',xct%ng,invband*invband,ONE,tempw(:,:,:,isv),xct%ng,mccp(:,:,:,isc), &
            1,ZERO,bsedwing(iitbeg:iitend,isc,isv),1)

        endif

        if (iit + (incband*invband)**2 > leading_dim) then
          write(0,*) peinf%inode, ik, ikp
          write(0,*) peinf%wown(1,1,ikp,1,1,ik), iit
          write(0,*) peinf%nckpe, peinf%myown, leading_dim
          call die("Internal error in g_sum with array dimensions.")
        endif

        ! Head
        do icp = 1, incband
          do ic = 1, incband      
            do ivp = 1, invband
              do iv = 1, invband
                iit = iit + 1
                bsedhead(iit,isc,isv) = bsedhead(iit,isc,isv) + &
                  mccp(1,ic,icp,isc) * temph(iv,ivp,isv)
              enddo
            enddo
          enddo
        enddo

      enddo
    enddo

    POP_SUB(g_sum_TDA_cpu)

  end subroutine g_sum_TDA_cpu

  subroutine g_sum_TDA_openacc(xct,invband,incband,temph,tempb,mccp, &
    bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik,&
    tempw,bsedwing)
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: invband,incband,leading_dim
    integer, intent(in) :: ivp_in,icp_in,ikp,iv_in,ic_in,ik
    SCALAR, intent(in), target :: tempb(:,:,:,:), temph(:,:,:), mccp(:,:,:,:)
    SCALAR, intent(inout), target :: bsedhead(:,:,:), bsedbody(:,:,:)
    SCALAR, intent(in), optional, target :: tempw(:,:,:,:)
    SCALAR, intent(inout), optional, target :: bsedwing(:,:,:)
#ifdef OPENACC
    integer :: isv, isc, iv, ivp, ic, icp, iit, iitend, iitbeg, iitbeg_h, iblock

    PUSH_SUB(g_sum_TDA_openacc)

    do isv=1,xct%nspin
      do isc=1,xct%nspin

        if (xct%icpar .eq. 0) then
          iit = peinf%wown(1,1,ikp,1,1,ik) - 1
        else if (xct%ivpar .eq. 0) then
          iit = peinf%wown(1,icp_in,ikp,1,ic_in,ik) - 1
        else
          iit = peinf%wown(ivp_in,icp_in,ikp,iv_in,ic_in,ik) - 1
        endif

        if (iit + (incband*invband)**2 > leading_dim) then
          write(0,*) peinf%inode, ik, ikp
          write(0,*) peinf%wown(1,1,ikp,1,1,ik), iit
          write(0,*) peinf%nckpe, peinf%myown, leading_dim
          call die("Internal error in g_sum with array dimensions.")
        endif

        iblock = iit / (invband*invband*incband*incband) + 1

        iitbeg = iit + 1
        iitend = iit + invband*invband*incband*incband

        iitbeg_h = iit

        ! Body term
        !$acc data present(tempb_accel, tempw_accel, temph_accel, bsedbody, bsedwing, bsedhead) copyin(mccp)
        if (incband .gt. 1) then
          call accel_xgemm('t', 'n', &
                         invband*invband, incband*incband, xct%ng, &
                         ONE, &
                         tempb_accel(:,:,:,isv), xct%ng, &
                         mccp(:,:,:,isc), xct%ng, &
                         ONE, &
                         bsedbody(iitbeg:iitend,isc,isv), invband*invband, &
                         g_sum_algo)
        else
          call accel_xgemv('t', &
                         xct%ng, invband*invband, &
                         ONE, &
                         tempb_accel(:,:,:,isv), xct%ng, &
                         mccp(:,:,:,isc), 1, &
                         ZERO, &
                         bsedbody(iitbeg:iitend,isc,isv), 1, &
                         g_sum_algo)
        endif

        if (present(bsedwing)) then
          if (incband .gt. 1) then
            call accel_xgemm('t', 'n', &
                           invband*invband, incband*incband, xct%ng, &
                           ONE, &
                           tempw_accel(:,:,:,isv), xct%ng, &
                           mccp(:,:,:,isc), xct%ng, &
                           ONE, &
                           bsedwing(iitbeg:iitend,isc,isv), invband*invband, &
                           g_sum_algo)
          else
            call accel_xgemv('t', &
                           xct%ng, invband*invband, &
                           ONE, &
                           tempw_accel(:,:,:,isv), xct%ng, &
                           mccp(:,:,:,isc), 1, &
                           ZERO, &
                           bsedwing(iitbeg:iitend,isc,isv), 1, &
                           g_sum_algo)
          end if
        endif

        ! Head term
        !$acc parallel loop gang vector collapse(4)
        do icp = 1, incband
          do ic = 1, incband
            do ivp = 1, invband
              do iv = 1, invband
                ! "iit = iit + 1" rewritten to eliminate loop carried dependency
                iit = iitbeg_h+ iv + &
                               (ivp-1)*invband + &
                               (ic-1)*invband*invband + &
                               (icp-1)*incband*invband*invband
                bsedhead(iit,isc,isv) = bsedhead(iit,isc,isv) + &
                  mccp(1,ic,icp,isc) * temph_accel(iv,ivp,isv)
              enddo
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$acc end data
      enddo
    enddo
#else
    PUSH_SUB(g_sum_TDA_openacc)

    call die("OpenACC version of g_sum requested, but OpenACC not compiled "//&
             "into this executable", only_root_writes = .true.)
#endif

    POP_SUB(g_sum_TDA_openacc)
  end subroutine g_sum_TDA_openacc

  subroutine g_sum_TDA_omp_target(xct,invband,incband,temph,tempb,mccp, &
    bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik,&
    tempw,bsedwing)
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: invband,incband,leading_dim
    integer, intent(in) :: ivp_in,icp_in,ikp,iv_in,ic_in,ik
    SCALAR, intent(in), target :: tempb(:,:,:,:), temph(:,:,:), mccp(:,:,:,:)
    SCALAR, intent(inout), target :: bsedhead(:,:,:), bsedbody(:,:,:)
    SCALAR, intent(in), optional, target :: tempw(:,:,:,:)
    SCALAR, intent(inout), optional, target :: bsedwing(:,:,:)
#ifdef OMP_TARGET
    integer :: isv, isc, iv, ivp, ic, icp, iit, iitend, iitbeg, iitbeg_h, iblock

    PUSH_SUB(g_sum_TDA_omp_target)

    do isv=1,xct%nspin
      do isc=1,xct%nspin

        if (xct%icpar .eq. 0) then
          iit = peinf%wown(1,1,ikp,1,1,ik) - 1
        else if (xct%ivpar .eq. 0) then
          iit = peinf%wown(1,icp_in,ikp,1,ic_in,ik) - 1
        else
          iit = peinf%wown(ivp_in,icp_in,ikp,iv_in,ic_in,ik) - 1
        endif

        if (iit + (incband*invband)**2 > leading_dim) then
          write(0,*) peinf%inode, ik, ikp
          write(0,*) peinf%wown(1,1,ikp,1,1,ik), iit
          write(0,*) peinf%nckpe, peinf%myown, leading_dim
          call die("Internal error in g_sum with array dimensions.")
        endif

        iblock = iit / (invband*invband*incband*incband) + 1

        iitbeg = iit + 1
        iitend = iit + invband*invband*incband*incband

        iitbeg_h = iit

        ! Body term
        !$omp target data map(to: mccp)
        if (incband .gt. 1) then
          call accel_xgemm('t', 'n', &
                         invband*invband, incband*incband, xct%ng, &
                         ONE, &
                         tempb_accel(:,:,:,isv), xct%ng, &
                         mccp(:,:,:,isc), xct%ng, &
                         ONE, &
                         bsedbody(iitbeg:iitend,isc,isv), invband*invband, &
                         g_sum_algo)
        else
          call accel_xgemv('t', &
                         xct%ng, invband*invband, &
                         ONE, &
                         tempb_accel(:,:,:,isv), xct%ng, &
                         mccp(:,:,:,isc), 1, &
                         ZERO, &
                         bsedbody(iitbeg:iitend,isc,isv), 1, &
                         g_sum_algo)
        endif

        if (present(bsedwing)) then
          if (incband .gt. 1) then
            call accel_xgemm('t', 'n', &
                           invband*invband, incband*incband, xct%ng, &
                           ONE, &
                           tempw_accel(:,:,:,isv), xct%ng, &
                           mccp(:,:,:,isc), xct%ng, &
                           ONE, &
                           bsedwing(iitbeg:iitend,isc,isv), invband*invband, &
                           g_sum_algo)
          else
            call accel_xgemv('t', &
                           xct%ng, invband*invband, &
                           ONE, &
                           tempw_accel(:,:,:,isv), xct%ng, &
                           mccp(:,:,:,isc), 1, &
                           ZERO, &
                           bsedwing(iitbeg:iitend,isc,isv), 1, &
                           g_sum_algo)
          end if
        endif

        ! Head term
        !$omp target teams loop default(none) &
        !$omp&     shared(incband, invband, iitbeg_h, isv, isc) &
        !$omp&     private(icp, ic, ivp, iv, iit) &
        !$omp& collapse(4)
        do icp = 1, incband
          do ic = 1, incband
            do ivp = 1, invband
              do iv = 1, invband
                ! "iit = iit + 1" rewritten to eliminate loop carried dependency
                iit = iitbeg_h+ iv + &
                               (ivp-1)*invband + &
                               (ic-1)*invband*invband + &
                               (icp-1)*incband*invband*invband
                bsedhead(iit,isc,isv) = bsedhead(iit,isc,isv) + &
                  mccp(1,ic,icp,isc) * temph_accel(iv,ivp,isv)
              enddo
            enddo
          enddo
        enddo
        !$omp end target teams loop
        !$omp end target data
      enddo
    enddo
#else
    PUSH_SUB(g_sum_TDA_omp_target)

    call die("OpenMP Target version of g_sum requested, but OpenMP Target not compiled "//&
             "into this executable", only_root_writes = .true.)
#endif

    POP_SUB(g_sum_TDA_omp_target)
  end subroutine g_sum_TDA_omp_target

  !> FHJ: Performs the summation over Gp to get the direct kernel term.
  !! This routine is generalized for extended kernels calculations.
  subroutine g_sum_extended(xct,ofs2,ofs2p,n2,n2p,temph,tempb,m22p, &
    bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik, &
    tempw, bsedwing)
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: ofs2, ofs2p
    integer, intent(in) :: n2, n2p
    integer, intent(in) :: leading_dim
    integer, intent(in) :: ivp_in,icp_in,ikp,iv_in,ic_in,ik
    SCALAR, intent(in) :: tempb(:,:,:,:), temph(:,:,:), m22p(:,:,:,:)
    SCALAR, intent(inout) :: bsedhead(:,:,:), bsedbody(:,:,:)
    SCALAR, intent(in), optional :: tempw(:,:,:,:)
    SCALAR, intent(inout), optional :: bsedwing(:,:,:)

    integer :: isv, isc, i1, i1p, i2, i2p, it, itbeg, it_buf
    integer :: n_left, n_right, buf_sz
    SCALAR, pointer :: buf_b(:), buf_w(:)

    PUSH_SUB(g_sum_extended)

    ! FHJ: n_left/n_right is the total number of val/cond states that this PE
    ! deals with. The number will be nv+nc for extended kernel calculations.
    ! However, whenever we deal with m22p, we should actually use n2/n2p +
    ! offsets, b/c the matrix m22p is decomposed in blocks.
    n_left = xct%n1b_co
    if (xct%ivpar==1) n_left = 1
    n_right = xct%n2b_co
    if (xct%icpar==1) n_right = 1
    buf_sz = n_left**2 * n2 * n2p
    SAFE_ALLOCATE(buf_b, (buf_sz))
    if (present(bsedwing)) then
      SAFE_ALLOCATE(buf_w, (buf_sz))
    endif

    do isv=1,xct%nspin
      do isc=1,xct%nspin

        if (xct%icpar .eq. 0) then
          itbeg = peinf%wown(1,1,ikp,1,1,ik)
        else if (xct%ivpar .eq. 0) then 
          itbeg = peinf%wown(1,icp_in,ikp,1,ic_in,ik)
        else
          itbeg = peinf%wown(ivp_in,icp_in,ikp,iv_in,ic_in,ik)
        endif

        if (n2>1 .or. n2p>1) then
          call X(gemm)('t','n',n_left**2,n2*n2p,xct%ng, &
            ONE,tempb(:,:,:,isv),xct%ng,m22p(:,:,:,isc),xct%ng, &
            ZERO,buf_b(:),n_left**2)

          if (present(bsedwing)) &
            call X(gemm)('t','n',n_left**2,n2*n2p,xct%ng, &
            ONE,tempw(:,:,:,isv),xct%ng,m22p(:,:,:,isc),xct%ng, &
            ZERO,buf_w(:),n_left**2)
        else
          call X(gemv)('t',xct%ng,n_left**2,ONE,tempb(:,:,:,isv),xct%ng,m22p(:,:,:,isc), &
            1,ZERO,buf_b(:),1)

          if (present(bsedwing)) &
            call X(gemv)('t',xct%ng,n_left**2,ONE,tempw(:,:,:,isv),xct%ng,m22p(:,:,:,isc), &
            1,ZERO,buf_w(:),1)
        endif

        if (itbeg-1 + (n_left*n_right)**2 > leading_dim) then
          write(0,*) peinf%inode, ik, ikp
          write(0,*) peinf%wown(1,1,ikp,1,1,ik), itbeg, n_left, n_right
          write(0,*) peinf%nckpe, peinf%myown, leading_dim
          call die("Internal error in g_sum with array dimensions.")
        endif

        do i2p = 1, n2p
          do i2 = 1, n2
            ! FHJ: max(it_buf) = n_left**2 * n2*n2p
            it_buf = 1 + n_left*n_left*((i2-1) + n2*(i2p-1))
            ! FHJ: max(it) = n_left**2 * n_right**2
            it = itbeg + n_left*n_left*((ofs2+i2-1) + n_right*(ofs2p+i2p-1))
            do i1p = 1, n_left
              do i1 = 1, n_left
                bsedhead(it,isc,isv) = temph(i1,i1p,isv) * m22p(1,i2,i2p,isc)
                if (present(bsedwing)) bsedwing(it,isc,isv) = buf_w(it_buf)
                bsedbody(it,isc,isv) = buf_b(it_buf)
                it = it + 1
                it_buf = it_buf + 1
              enddo
            enddo
          enddo
        enddo

      enddo
    enddo

    SAFE_DEALLOCATE_P(buf_b)
    if (present(bsedwing)) then
      SAFE_DEALLOCATE_P(buf_w)
    endif

    POP_SUB(g_sum_extended)

  end subroutine g_sum_extended

end module g_sum_m
