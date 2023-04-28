!==================================================================================
!
! Routines:
!
! (1) iterate()         Originally By MLT               Last Modified 7/1/2008 (JRD)
!
!     Performs a pair of iterations using the moments (Haydock) algorithm
!     for the haydock code.
!
!     input:  mmts, xct types
!            nmat        number of exciton states
!            nit         order of iteration
!            hbse_a       Hamiltonian
!            s1          Haydock vector of order |nit-1>
!            s0          Haydock vector of order |nit>
!            mmts%an(1:nit-1)
!            mmts%bn(1:nit-1)
!
!     output: s1          Haydock vector of order |nit+1>
!             s0          Haydock vector of order |nit+2>
!             mmts%an(1:nit+1)
!             mmts%bn(1:nit+1)
!
! (2) local_s()         Originally By MLT               Last Modified 7/1/2008 (JRD)
!
!     Given a global state ss, defines ss_l with the same parallel
!     distribution as the second index of the Hamiltonian, so that
!     one can easily apply ss to the right of hbse_a by doing
!     matmul(hbse_a,ss_l).
!
!     input: peinf type
!            nmat        number of exciton states
!            ss          global state, size nmat
!
!     output: ss_l       local state, size peinf%nkpe*nmat
!
!====================================================================================

#include "f_defs.h"

module iterate_m

  use global_m
  use blas_m
  use scalapack_m
  use misc_m
  implicit none

  private

  public :: iterate, local_s

contains


subroutine iterate(mmts,xct,nmat,nit,hbse_a,s1,s0)
  type (mmtsinfo), intent(inout) :: mmts
  type (xctinfo), intent(in) :: xct
  integer, intent(in) :: nmat,nit
  SCALAR, intent(in) :: hbse_a(nmat,peinf%nblocks*peinf%block_sz)
  SCALAR, intent(inout) :: s1(nmat),s0(nmat)

  character :: filename*12
  integer :: ii,ik,ic,iv,is,ikcvs,iprint,itth,ith
  real(DP) :: sum_, sum1, dtemp
  SCALAR :: s1_l(peinf%nblocks*peinf%block_sz), xtemp, &
    s0_l(peinf%nblocks*peinf%block_sz),dum1(nmat),dum2(nmat), &
    xtemp_vec(nmat)

  PUSH_SUB(iterate)

!----------------------------
! Print out

  if (peinf%inode.eq.0) then
    iprint=0
    if (iprint.ne.0) then
      write(6,*) ' Output of Haydock vectors',nit-1,nit
      write(filename,'(a,i3.3)')  'HAYDOCKN_', nit

      itth=31
      call open_file(itth,file=filename,form='formatted',status='replace')
      do ik=1,xct%nkpt_fi
        do ic=1,xct%ncb_fi
          do iv=1,xct%nvb_fi
            do is=1,xct%nspin
              ikcvs = bse_index(ik, ic, iv, is, xct)
              write(itth,*) ik,ic,iv,is,ikcvs,s1(ikcvs),s0(ikcvs)
            enddo
          enddo
        enddo
      enddo
      call close_file(itth)
    endif
  endif

!---------------------------------
! Initialize local states : s1_l  --> local part of s1
!                           s0_l  --> local part of s0
! Whenever s1 and s0 change, their local parts must be updated

  call local_s(nmat,s1,s1_l)
  call local_s(nmat,s0,s0_l)

! Calculate an(nit) from vector s0 = |nit>

  ! FHJ: We have to use BLAS because Intel incorrectly optimizes this
  ! sum1 = dble( DOT_PRODUCT(s0,MATMUL(hbse_a,s0_l)) )
  call X(gemv)('N', nmat, peinf%nblocks*peinf%block_sz, ONE, &
    hbse_a, nmat, s0_l, 1, ZERO, xtemp_vec, 1)
  xtemp = blas_dot(nmat, s0, 1, xtemp_vec, 1)
  sum1 = MYREAL(xtemp)

#ifdef MPI
  call MPI_Allreduce(sum1,sum_,1,MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
  sum1 = sum_
#endif
  mmts%an(nit) = sum1

!----------------------------------
! Calculate vector s1 = |nit+1> and bn(nit)
! Note that bn corresponds to the b^2_n in Benedict & Shirley, eq. (22)

  ! FHJ: We have to use BLAS because Intel incorrectly optimizes this
  ! dum2 = MATMUL(hbse_a,s0_l)
  call X(gemv)('N', nmat, peinf%nblocks*peinf%block_sz, ONE, &
    hbse_a, nmat, s0_l, 1, ZERO, dum2, 1)

#ifdef MPI
  call MPI_Allreduce(dum2,dum1,nmat,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
  dum2 = dum1
#endif

  dum1(:) = dum2(:) - mmts%an(nit)*s0(:) - sqrt(mmts%bn(nit-1))*s1(:)
  s1(:) = dum1(:)

  ! FHJ: We have to use BLAS because Intel incorrectly optimizes this
  dtemp = blas_nrm2(nmat, s1, 1)
  mmts%bn(nit) = dtemp**2
  s1(:) = s1(:)/dtemp
  call local_s(nmat,s1,s1_l)
  
  if (peinf%inode.eq.0) write(6,120) nit,mmts%an(nit),mmts%bn(nit)
  
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

!--------------------------------
! Calculate an(nit+1) from vector s1 = |nit+1>

  ! FHJ: We have to use BLAS because Intel incorrectly optimizes this
  ! sum1 = dble( DOT_PRODUCT(s1,MATMUL(hbse_a,s1_l)) )
  call X(gemv)('N', nmat, peinf%nblocks*peinf%block_sz, ONE, &
    hbse_a, nmat, s1_l, 1, ZERO, xtemp_vec, 1)
  xtemp = blas_dot(nmat, s1, 1, xtemp_vec, 1)
  sum1 = MYREAL(xtemp)

#ifdef MPI
  call MPI_Allreduce(sum1,sum_,1,MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
  sum1 = sum_
#endif
  mmts%an(nit+1) = sum1
      
!--------------------------------
! Calculate vector s0 = |nit+2> and bn(nit+1)

  ! FHJ: We have to use BLAS because Intel incorrectly optimizes this
  ! dum2 = MATMUL(hbse_a,s1_l)
  call X(gemv)('N', nmat, peinf%nblocks*peinf%block_sz, ONE, &
    hbse_a, nmat, s1_l, 1, ZERO, dum2, 1)

#ifdef MPI
  call MPI_Allreduce(dum2,dum1,nmat,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
  dum2 = dum1
#endif

  dum1(:) = dum2(:) - mmts%an(nit+1)*s1(:) - sqrt(mmts%bn(nit))*s0(:)
  s0(:) = dum1(:)

  ! FHJ: We have to use BLAS because Intel incorrectly optimizes this
  dtemp = blas_nrm2(nmat, s0, 1)
  mmts%bn(nit+1) = dtemp**2
  s0(:) = s0(:)/dtemp
  call local_s(nmat,s0,s0_l)
  
  if (peinf%inode.eq.0) write(6,120) nit+1,mmts%an(nit+1),mmts%bn(nit+1)

! Save coefficients after every 10 iterations

  iprint=0
  if (peinf%inode.eq.0 .and. mod(nit,10).eq.0) iprint=1
  if (iprint.ne.0) then
    ith = 21
    call open_file(ith,file='eps2_moments',form='unformatted',status='replace')
    write(ith) nit+1,mmts%norm,mmts%vol,nmat,xct%nspin
    write(ith) (mmts%an(ii),ii=1,nit+1)
    write(ith) (mmts%bn(ii),ii=1,nit+1)
    write(ith) (s1(ii),ii=1,nmat)
    write(ith) (s0(ii),ii=1,nmat)
    call close_file(ith)
    write(6,*) 'Vectors written in file'
  endif
  iprint=0
  
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
  
120 format(3x,i4,3x,2e12.5)
  
  POP_SUB(iterate)
  
  return
end subroutine iterate


!===========================================================================


subroutine local_s(nmat,ss,ss_l)
  integer, intent(in) :: nmat
  SCALAR, intent(in) :: ss(nmat)
  SCALAR, intent(out) :: ss_l(peinf%nblocks*peinf%block_sz)
  
  integer :: ii

  PUSH_SUB(local_s)
  
  ss_l = 0.d0
  
  do ii=1,peinf%ibt(peinf%inode+1)*peinf%block_sz
    ss_l(ii) = ss(peinf%peig(peinf%inode+1,ii))
  enddo
  
  POP_SUB(local_s)
  
  return
end subroutine local_s

end module iterate_m
