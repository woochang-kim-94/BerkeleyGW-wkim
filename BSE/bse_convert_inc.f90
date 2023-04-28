!============================================================================
!
! included from bse_convert.F90
!
!============================================================================

subroutine NAME(bsemat)(iunit, ounit)
  integer, intent(in) :: iunit, ounit

  type(kernel_header_t) :: kernel
  SCALAR, allocatable :: bsemat(:,:,:,:,:)
  integer :: ik, ikp, imat, ibb, ic, iv

  PUSH_SUB(NAME(bsemat))

  call read_kernel_header(iunit, kernel)
  call write_kernel_header(ounit, kernel)
  SAFE_ALLOCATE(bsemat, (kernel%nk, kernel%n2b, kernel%n1b, kernel%ns, kernel%ns))
  do ik = 1, kernel%nk
    do imat = 1, kernel%nmat
      do ibb = 1, kernel%n1b*kernel%n2b
        READWRITE(ikp,ic,iv,bsemat(:,:,:,:,:))
      enddo
    end do
  end do
  
  POP_SUB(NAME(bsemat))

end subroutine NAME(bsemat)

!============================================================================
! compare read and write in intwfn.f90
subroutine NAME(dtmat)(iunit, ounit)
  integer, intent(in) :: iunit, ounit

  integer :: nk, nc, nv, nkf, ncf, nvf, ns, nmat, &
    ik, ic, jk, jc, is, ii, iv, jv, npts, ndims
  logical :: per(3)
  real(DP) :: kk(1:3)
  SCALAR :: dcc, dvv

  PUSH_SUB(NAME(dtmat))

  READWRITE(ndims, per(1:3), npts, nk)

  READWRITE(nk,nc,nv,nkf,ncf,nvf,ns)
  do ik=1, nk
    READWRITE(kk(1:3))
  enddo

  ! cc
  nmat=nkf*ncf*nc*ns
  do ii=1,nmat
    READWRITE(ik,ic,jk,jc,is,dcc)
  enddo

  ! vv
  nmat=nkf*nvf*nv*ns
  do ii=1,nmat
    READWRITE(ik,iv,jk,jv,is,dvv)
  enddo
  
  POP_SUB(NAME(dtmat))
  return
end subroutine NAME(dtmat)

!============================================================================
subroutine NAME(vmtxel)(iunit, ounit)
  integer, intent(in) :: iunit, ounit

  integer :: nkf, ncf, nvf, ns, ic, nmat, ii
  SCALAR, allocatable :: s1(:)

  PUSH_SUB(NAME(vmtxel))

  READWRITE(nkf,ncf,nvf,ns,ic)

  nmat=nkf*ncf*nvf*ns
  SAFE_ALLOCATE(s1, (nmat))
  
  READWRITE((s1(ii),ii=1,nmat))
  SAFE_DEALLOCATE(s1)

  PUSH_SUB(NAME(vmtxel))
  return
end subroutine NAME(vmtxel)

!============================================================================
subroutine NAME(eps2_moments)(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  
  integer :: nn, nmat, ii
  real(DP) :: tmp(2)
  SCALAR, allocatable :: s1(:)
  real(DP), allocatable :: array(:)

  PUSH_SUB(NAME(eps2_moments))

  READWRITE(nn,tmp(1:2),nmat,ii)

  SAFE_ALLOCATE(array, (nn))
  ! an
  READWRITE((array(ii),ii=1,nn))
  ! bn
  READWRITE((array(ii),ii=1,nn))
  SAFE_DEALLOCATE(array)
  
  SAFE_ALLOCATE(s1, (nmat))
  READWRITE((s1(ii),ii=1,nmat))
  READWRITE((s1(ii),ii=1,nmat))
  SAFE_DEALLOCATE(s1)
  
  PUSH_SUB(NAME(eps2_moments))
  return
end subroutine NAME(eps2_moments)
  
!============================================================================
subroutine NAME(eigenvectors)(iunit, ounit)
  integer, intent(in) :: iunit, ounit

  integer :: ns, nv, nc, nk, nmat, ii, ik, isvck
  real(DP), allocatable :: kg(:,:)
  real(DP) :: energy
  SCALAR, allocatable :: Asvck(:)
  
  PUSH_SUB(NAME(eigenvectors))

  READWRITE(ns)
  READWRITE(nv)
  READWRITE(nc)
  READWRITE(nk)

  SAFE_ALLOCATE(kg, (3,nk))
  READWRITE(((kg(ii,ik),ii=1,3),ik=1,nk))
  SAFE_DEALLOCATE(kg)

  nmat = ns*nv*nc*nk
  SAFE_ALLOCATE(Asvck, (nmat))

  ! FIXME: there can be between 1 and nmat eigenvectors here, actually
  ! The code will crash if it is less than nmat, though what is written
  ! will be correct and complete nonetheless.
  do isvck = 1, nmat
    READWRITE(energy)
    READWRITE((Asvck(ik),ik=1,nmat))
  enddo
    
  SAFE_DEALLOCATE(Asvck)
    
  POP_SUB(NAME(eigenvectors))
  return
end subroutine NAME(eigenvectors)

