
!===============================================================================
!
! Routines:
!
! 1. orthonormalize()      Originally By gsm      Last Modified 9/18/2010 (gsm)
!
!    Orthonormalizes distributed wavefunctions in G-space by applying
!    Gram-Schmidt process.
!
! 2. orderindex()          Originally By gsm      Last Modified 9/18/2010 (gsm)
!
!    Builds an index array for reordering wavefunctions, eigenvalues and
!    arrays for fast scalar products.
!
! 3. orderarray()          Originally By gsm      Last Modified 9/18/2010 (gsm)
!
!    Reorders wavefunctions, eigenvalues and arrays for fast scalar products
!    according to the index array.
!
!===============================================================================

#include "f_defs.h"

module ortho_m

  use global_m
  use blas_m
  use sort_m

  implicit none

  private

  public :: orthonormalize, orderindex, orderarray

contains

subroutine orthonormalize(ffastsp,fblas,fecor,nk,ns,nbstart, &
  nbend,nbmax,ngk_l,npwdeg,pw_num,pw_ind,en,wfn_d)

  logical, intent(in) :: ffastsp,fblas,fecor
  integer, intent(in) :: nk,ns,nbstart,nbend,nbmax,ngk_l,npwdeg
  integer, intent(in) :: pw_num(:,:,:) !< (nbmax,ns,nk)
  integer, intent(in) :: pw_ind(:,:,:,:,:) !< (2,npwdeg,nbmax,ns,nk)
  real(DP), intent(inout) :: en(:,:,:) !< (nbmax,nk,ns)
  SCALAR, intent(inout) :: wfn_d(:,:,:,:) !< (ngk_l,nbmax,ns,nk)
  
  integer :: i,ik,is,ib,jb,ig
  logical :: ffastsp_ib
  real(DP) :: norm2,normi,de,desp2,norm2_dummy
  SCALAR :: z
  character(len=16) :: s1,s2,s3
  real(DP), allocatable :: wfnsp2(:)
  SCALAR, allocatable :: wfnsp(:)
  SCALAR, allocatable :: wfnspdummy(:)
  
  PUSH_SUB(orthonormalize)
  
  SAFE_ALLOCATE(wfnsp2, (nbend-1))
  SAFE_ALLOCATE(wfnsp, (nbend-1))
#ifdef MPI
  SAFE_ALLOCATE(wfnspdummy, (nbend-1))
#endif

  do ik=1,nk
    do is=1,ns
      do ib=nbstart,nbend

        ! check whether to compute fast scalar products
        ! for the current wavefunction
        ! (if .not.ffastsp pw_num may be unallocated)
        if (ffastsp) then
          ffastsp_ib=pw_num(ib,is,ik).ne.0
        else ! ffastsp
          ffastsp_ib=.false.
        endif ! ffastsp
        
        ! compute scalar products of the current wavefunction
        ! with the previous wavefunctions
        wfnsp(:)=ZERO
        if (ffastsp_ib) then
          ! fast scalar product
          do jb=1,ib-1
            z=ZERO
            do i=1,pw_num(ib,is,ik)
              if (peinf%inode.eq.pw_ind(1,i,ib,is,ik)) &
                z=z+MYCONJG(wfn_d(pw_ind(2,i,ib,is,ik),jb,is,ik))* &
                wfn_d(pw_ind(2,i,ib,is,ik),ib,is,ik)
            enddo ! i
            wfnsp(jb)=z
          enddo ! jb
        else ! ffastsp_ib
          ! full scalar product
          if (fblas) then
            call X(gemv)('C',ngk_l,ib-1,ONE,wfn_d(:,:,is,ik), &
              ngk_l,wfn_d(:,ib,is,ik),1,ZERO,wfnsp(:),1)
          else ! fblas
            do jb=1,ib-1
              z=ZERO
              do ig=1,ngk_l
                z=z+MYCONJG(wfn_d(ig,jb,is,ik))* &
                  wfn_d(ig,ib,is,ik)
              enddo ! ig
              wfnsp(jb)=z
            enddo ! jb
          endif ! fblas
        endif ! ffastsp_ib
#ifdef MPI
        wfnspdummy(1:ib-1)=wfnsp(1:ib-1)
        call MPI_Allreduce(wfnspdummy,wfnsp,ib-1, &
          MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
        
        ! orthogonalize the current wavefunction
        ! to the previous wavefunctions
        if (fblas) then
          call X(gemv)('N',ngk_l,ib-1,-ONE,wfn_d(:,:,is,ik), &
            ngk_l,wfnsp(:),1,ONE,wfn_d(:,ib,is,ik),1)
        else ! fblas
          do jb=1,ib-1
            do ig=1,ngk_l
              wfn_d(ig,ib,is,ik)=wfn_d(ig,ib,is,ik)- &
                wfnsp(jb)*wfn_d(ig,jb,is,ik)
            enddo ! ig
          enddo ! jb
        endif ! fblas
        
        ! compute the squared scalar products
        wfnsp2(:)=0.0d0
        do jb=1,ib-1
          wfnsp2(jb)=abs(wfnsp(jb))**2
        enddo ! jb
        
        ! normalize the current wavefunction
        norm2 = sum(abs(wfn_d(:,ib,is,ik))**2)
#ifdef MPI
        norm2_dummy = norm2
        call MPI_Allreduce(norm2_dummy, norm2, 1, &
          MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD, mpierr)
#endif
        
        if (norm2.gt.TOL_Zero) then
          normi=1.0d0/sqrt(norm2)
          if (fblas) then
            call X(SCAL)(ngk_l,normi*ONE,wfn_d(:,ib,is,ik),1)
          else ! fblas
            do ig=1,ngk_l
              wfn_d(ig,ib,is,ik)=wfn_d(ig,ib,is,ik)*normi
            enddo ! ig
          endif ! fblas
        else ! norm2.gt.TOL_Zero
          write(s1,151)ib
          write(s2,151)is
          write(s3,151)ik
          ! aka, "ortho went bananas"
          call die('Failed to orthonormalize band ' // TRUNC(s1) // &
            ' for spin ' // TRUNC(s2) // ' and k-point ' // TRUNC(s3))
        endif ! norm2.gt.TOL_Zero
        
        ! correct the eigenvalue of the current
        ! wavefunction by perturbative approach
        if (fecor) then
          desp2=0.0d0
          do jb=1,ib-1
            de=en(ib,ik,is)-en(jb,ik,is)
            desp2=desp2+de*wfnsp2(jb)
          enddo ! jb
          en(ib,ik,is)=en(ib,ik,is)+desp2/norm2
        endif ! fecor
        
      enddo ! ib
    enddo ! is
  enddo ! ik
  
  SAFE_DEALLOCATE(wfnsp2)
  SAFE_DEALLOCATE(wfnsp)
#ifdef MPI
  SAFE_DEALLOCATE(wfnspdummy)
#endif
  
  POP_SUB(orthonormalize)
  
  return
  
151 format(i16)
  
end subroutine orthonormalize

!-------------------------------------------------------------------------------

subroutine orderindex(iblock,iorder,nk,ns,nbtot,nbinp,nbpw, &
  nbaux,nbmax,isrtort,en)
  
  integer, intent(in) :: iblock,iorder,nk,ns,nbtot,nbinp,nbpw, &
    nbaux,nbmax
  integer, intent(out) :: isrtort(:,:,:) !< (nbtot,ns,nk)
  real(DP), intent(in) :: en(:,:,:) !< (nbmax,nk,ns)
  
  integer :: ik,is,ib
  real(DP), allocatable :: ensrt(:)
  integer, allocatable :: isrt(:)
  
  PUSH_SUB(orderindex)
  
  isrtort(:,:,:)=0
  
  if (iblock.eq.2) then
    SAFE_ALLOCATE(ensrt, (nbtot-nbinp))
    SAFE_ALLOCATE(isrt, (nbtot-nbinp))
  endif ! iblock.eq.2
  
  do ik=1,nk
    do is=1,ns
      
      if (iblock.eq.0) then
        
        if (iorder.eq.0) then
          do ib=1,nbtot
            isrtort(ib,is,ik)=ib
          enddo ! ib
        elseif (iorder.eq.1) then
          do ib=1,nbinp
            isrtort(ib,is,ik)=nbinp+1-ib
          enddo ! ib
          do ib=nbinp+1,nbinp+nbpw
            isrtort(ib,is,ik)=nbinp+1+nbinp+nbpw-ib
          enddo ! ib
          do ib=nbinp+nbpw+1,nbtot
            isrtort(ib,is,ik)=nbinp+nbpw+1+nbtot-ib
          enddo ! ib
        endif ! iorder
        
      elseif (iblock.eq.1) then
        
        if (iorder.eq.0) then
          do ib=1,nbinp
            isrtort(ib,is,ik)=ib
          enddo ! ib
          do ib=nbinp+1,nbinp+nbaux
            isrtort(ib,is,ik)=ib+nbpw
          enddo ! ib
          do ib=nbinp+nbaux+1,nbtot
            isrtort(ib,is,ik)=ib-nbaux
          enddo ! ib
        elseif (iorder.eq.1) then
          do ib=1,nbinp
            isrtort(ib,is,ik)=nbinp+1-ib
          enddo ! ib
          do ib=nbinp+1,nbtot
            isrtort(ib,is,ik)=nbinp+1+nbtot-ib
          enddo ! ib
        endif ! iorder
        
      elseif (iblock.eq.2) then
        
        do ib=nbinp+1,nbtot
          ensrt(ib-nbinp)=en(ib,ik,is)
        enddo ! ib
        call sortrx(nbtot-nbinp,ensrt,isrt)
        
        if (iorder.eq.0) then
          do ib=1,nbinp
            isrtort(ib,is,ik)=ib
          enddo ! ib
          do ib=nbinp+1,nbtot
            isrtort(ib,is,ik)=isrt(ib-nbinp)+nbinp
          enddo ! ib
        elseif (iorder.eq.1) then
          do ib=1,nbinp
            isrtort(ib,is,ik)=nbinp+1-ib
          enddo ! ib
          do ib=nbinp+1,nbtot
            isrtort(ib,is,ik)=isrt(nbtot+1-ib)+nbinp
          enddo ! ib
        endif ! iorder
        
      endif ! iblock
      
    enddo ! is
  enddo ! ik
  
  if (iblock.eq.2) then
    SAFE_DEALLOCATE(ensrt)
    SAFE_DEALLOCATE(isrt)
  endif ! iblock.eq.2
  
  POP_SUB(orderindex)
  
  return
  
end subroutine orderindex

!-------------------------------------------------------------------------------

subroutine orderarray(ffastsp,idir,nk,ns,nbtot,nbmax,ngk_l, &
  npwdeg,isrtort,pw_num,pw_ind,en,wfn_d)
  
  logical, intent(in) :: ffastsp
  integer, intent(in) :: idir,nk,ns,nbtot,nbmax,ngk_l,npwdeg
  integer, intent(in) :: isrtort(:,:,:) !< (nbtot,ns,nk)
  integer, intent(inout) :: pw_num(:,:,:) !< (nbmax,ns,nk)
  integer, intent(inout) :: pw_ind(:,:,:,:,:) !< (2,npwdeg,nbmax,ns,nk)
  real(DP), intent(inout) :: en(:,:,:) !< (nbmax,nk,ns)
  SCALAR, intent(inout) :: wfn_d(:,:,:,:) !< (ngk_l,nbmax,ns,nk)
  
  integer :: ik,is,ib
  integer, allocatable :: pw_numsrt(:)
  integer, allocatable :: pw_indsrt(:,:,:)
  real(DP), allocatable :: ensrt(:)
  SCALAR, allocatable :: wfnsrt_d(:,:)
  
  PUSH_SUB(orderarray)
  
  if (ffastsp) then
    SAFE_ALLOCATE(pw_numsrt, (nbtot))
    SAFE_ALLOCATE(pw_indsrt, (2,npwdeg,nbtot))
  endif ! ffastsp
  SAFE_ALLOCATE(ensrt, (nbtot))
  SAFE_ALLOCATE(wfnsrt_d, (ngk_l,nbtot))
  
  if (idir.eq.1) then
    
    do ik=1,nk
      do is=1,ns
        if (ffastsp) then
          do ib=1,nbtot
            pw_numsrt(ib)=pw_num(isrtort(ib,is,ik),is,ik)
          enddo ! ib
          do ib=1,nbtot
            pw_num(ib,is,ik)=pw_numsrt(ib)
          enddo ! ib
          do ib=1,nbtot
            pw_indsrt(:,:,ib)=pw_ind(:,:,isrtort(ib,is,ik),is,ik)
          enddo ! ib
          do ib=1,nbtot
            pw_ind(:,:,ib,is,ik)=pw_indsrt(:,:,ib)
          enddo ! ib
        endif ! ffastsp
        do ib=1,nbtot
          ensrt(ib)=en(isrtort(ib,is,ik),ik,is)
        enddo ! ib
        do ib=1,nbtot
          en(ib,ik,is)=ensrt(ib)
        enddo ! ib
        do ib=1,nbtot
          wfnsrt_d(:,ib)=wfn_d(:,isrtort(ib,is,ik),is,ik)
        enddo ! ib
        do ib=1,nbtot
          wfn_d(:,ib,is,ik)=wfnsrt_d(:,ib)
        enddo ! ib
      enddo ! is
    enddo ! ik
    
  elseif (idir.eq.2) then
    
    do ik=1,nk
      do is=1,ns
        if (ffastsp) then
          do ib=1,nbtot
            pw_numsrt(isrtort(ib,is,ik))=pw_num(ib,is,ik)
          enddo ! ib
          do ib=1,nbtot
            pw_num(ib,is,ik)=pw_numsrt(ib)
          enddo ! ib
          do ib=1,nbtot
            pw_indsrt(:,:,isrtort(ib,is,ik))=pw_ind(:,:,ib,is,ik)
          enddo ! ib
          do ib=1,nbtot
            pw_ind(:,:,ib,is,ik)=pw_indsrt(:,:,ib)
          enddo ! ib
        endif ! ffastsp
        do ib=1,nbtot
          ensrt(isrtort(ib,is,ik))=en(ib,ik,is)
        enddo ! ib
        do ib=1,nbtot
          en(ib,ik,is)=ensrt(ib)
        enddo ! ib
        do ib=1,nbtot
          wfnsrt_d(:,isrtort(ib,is,ik))=wfn_d(:,ib,is,ik)
        enddo ! ib
        do ib=1,nbtot
          wfn_d(:,ib,is,ik)=wfnsrt_d(:,ib)
        enddo ! ib
      enddo ! is
    enddo ! ik
    
  endif ! idir
  
  if (ffastsp) then
    SAFE_DEALLOCATE(pw_numsrt)
    SAFE_DEALLOCATE(pw_indsrt)
  endif ! ffastsp
  SAFE_DEALLOCATE(ensrt)
  SAFE_DEALLOCATE(wfnsrt_d)
  
  POP_SUB(orderarray)
  
  return
  
end subroutine orderarray

end module ortho_m

