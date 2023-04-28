!============================================================================
!
! Routines:
!
! (1) distrib()         Originally By MLT       Last Modified (JRD) 5/29/08
!
!     Distributes k-points among processors.  It is best if number of kpoints
!     divides number of processor evenly, otherwise who knows...
!
!     nk = # of kpts
!     np = # of processors
!
!           peinf%nkpe = floor(nk/np) = # of kpts per proc
!           peinf%ikt(p) = # of kpoints belonging to proc p
!           peinf%ik(p,i) = kpt label (1..xct%nkpt_fi) of ith kpts on proc p
!
!
!     input: xct%nkpt_fi
!            xct%ncb_fi
!            xct%nvb_fi
!            xct%neps
!            xct%nspin
!
!=============================================================================

#include "f_defs.h"

module distrib_m

  use fftw_m
  use global_m
  use misc_m
  use scalapack_m
  implicit none

  private

  public :: &
    distrib

contains

subroutine distrib(xct,FFTgrid,kgrid,is_diag)
  type (xctinfo), intent(inout) :: xct
  integer, intent(in) :: FFTgrid(3),kgrid(3)
  logical, intent(in) :: is_diag

  integer :: ik,ipe,nmat,nmat_co,ic,iv,nblockstotal
  real(DP) :: mem,fmem,dmem,hmem,rmem,rmem2,scale,dscale,facdyn
  integer :: Nrod,Nplane,Nfft(3),dNfft(3),dkmax(3),nmpinode
  integer :: irank, rank_last, ibt_last
#ifdef MPI
  integer :: numbl_max,nnp,nn,np0,mq0,nq0,cls,neig,lrwork,lwork,nbr,nbc
  type(scalapack) :: scal
#endif

!---------------------------------
! Distribute work
! Each processor will own peinf%nblocks of size peinf%block_sz

  PUSH_SUB(distrib)

  peinf%nkpe=xct%nkpt_fi/peinf%npes
  if(xct%nkpt_fi-peinf%npes*peinf%nkpe.gt.0) peinf%nkpe=peinf%nkpe+1
  
  if (peinf%inode .eq. 0) then
    write(6,*)
  endif
  if (peinf%npes .le. xct%nkpt_fi) then
    peinf%nblocks=xct%nkpt_fi/peinf%npes
    xct%ipar = 1
    peinf%nc_block = xct%ncb_fi
    peinf%nv_block = xct%nvb_fi
    if(xct%nkpt_fi-peinf%npes*peinf%nblocks.gt.0) &
      peinf%nblocks=peinf%nblocks+1
    nblockstotal = xct%nkpt_fi
    if (peinf%inode .eq. 0) then
      write(6,*) 'You are doing a class 1 calculation - npes <= nk' 
    endif
  else if (peinf%npes .le. xct%nkpt_fi * xct%ncb_fi) then
    peinf%nblocks=xct%nkpt_fi*xct%ncb_fi/peinf%npes
    xct%ipar = 2
    peinf%nc_block = 1
    peinf%nv_block = xct%nvb_fi
    if(xct%nkpt_fi*xct%ncb_fi-peinf%npes*peinf%nblocks.gt.0) &
      peinf%nblocks=peinf%nblocks+1
    nblockstotal = xct%nkpt_fi*xct%ncb_fi
    if (peinf%inode .eq. 0) then
      write(6,*) 'You are doing a class 2 calculation - npes <= nk*nc' 
    endif
  else
    peinf%nblocks=xct%nkpt_fi*xct%nvb_fi*xct%ncb_fi/peinf%npes
    xct%ipar = 3
    peinf%nc_block = 1
    peinf%nv_block = 1
    if(xct%nkpt_fi*xct%nvb_fi*xct%ncb_fi-peinf%npes*peinf%nblocks.gt.0) &
      peinf%nblocks=peinf%nblocks+1
    nblockstotal = xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi
    if (peinf%inode .eq. 0) then
      write(6,*) 'You are doing a class 3 calculation - npes > nk*nc' 
    endif
  endif
  
  SAFE_ALLOCATE(peinf%ik, (peinf%npes,peinf%nkpe))
  SAFE_ALLOCATE(peinf%ikb, (peinf%nblocks))
  SAFE_ALLOCATE(peinf%ivb, (peinf%nblocks))
  SAFE_ALLOCATE(peinf%icb, (peinf%nblocks))
  SAFE_ALLOCATE(peinf%ikt, (peinf%npes))
  SAFE_ALLOCATE(peinf%ibt, (peinf%npes))

  peinf%ik=0
  peinf%ikt=0
  peinf%ikb=0
  peinf%icb=0
  peinf%ivb=0
  peinf%ibt=0

  ! FHJ: This is for the WFN distribution in intwfn
  ipe = 0
  do ik=1,xct%nkpt_fi
    ipe = mod(ipe,peinf%npes) + 1
    peinf%ikt(ipe)=peinf%ikt(ipe)+1
    peinf%ik(ipe,peinf%ikt(ipe))=ik
  enddo

  ! FHJ: This is for the hbse distribution in intkernel/diagonalize/etc.
  ipe = 0
  ! FHJ: Always loop over k because we always distribute over k-points
  ! (ie, implicitly assume that "nk_block==1")
  do ik=1,xct%nkpt_fi
    ! FHJ: Only loop over cond bands if we distribute them
    ! (i.e., only if peinf%nc_block==1)
    do ic=1,xct%ncb_fi - peinf%nc_block + 1
      ! FHJ: Only loop over val bands if we distribute them
      ! (i.e., only if peinf%nv_block==1)
      do iv=1,xct%nvb_fi - peinf%nv_block + 1
        ipe = mod(ipe,peinf%npes) + 1
        peinf%ibt(ipe) = peinf%ibt(ipe) + 1
        if (ipe==peinf%inode+1) then
          peinf%ikb(peinf%ibt(ipe)) = ik
          peinf%ivb(peinf%ibt(ipe)) = iv
          peinf%icb(peinf%ibt(ipe)) = ic
        endif
      enddo
    enddo
  enddo
  
!---------------------------------
! Determine the available memory

  call procmem(mem,nmpinode)
  
  if(peinf%inode.eq.0) then
    write(6,998) mem/1024.0d0**2
  endif
998 format(/,1x,'Memory available: ',f0.1,' MB per PE')

  fmem=mem/8.0d0

!---------------------------------------------------------
! (gsm) Determine the amount of memory required for vcoul

! random numbers
  rmem=0.0D0
  if (xct%icutv/=TRUNC_BOX) then
! arrays ran, qran, and qran2
! (ran is deallocated before qran2 is allocated)
    rmem=rmem+6.0D0*dble(nmc)*8.0D0
  endif
! various truncation schemes
  call setup_FFT_sizes(FFTgrid,Nfft,scale)
  rmem2=0.0d0
! cell wire truncation
  if (xct%icutv==TRUNC_WIRE) then
    dkmax(1) = FFTgrid(1) * n_in_wire
    dkmax(2) = FFTgrid(2) * n_in_wire
    dkmax(3) = 1
    call setup_FFT_sizes(dkmax,dNfft,dscale)
! array fftbox_2D
    rmem2=rmem2+dble(dNfft(1))*dble(dNfft(2))*16.0d0
! array inv_indx
    rmem2=rmem2+dble(Nfft(1))*dble(Nfft(2))*dble(Nfft(3))*4.0d0
! array qran
    rmem2=rmem2+3.0D0*dble(nmc)*8.0D0
  endif
! cell box (parallel version only) & supercell box truncation
  if (xct%icutv==TRUNC_BOX .or. xct%icutv==TRUNC_SUPERCELL) then
    dkmax(1) = FFTgrid(1) * n_in_box
    dkmax(2) = FFTgrid(2) * n_in_box
    dkmax(3) = FFTgrid(3) * n_in_box
    call setup_FFT_sizes(dkmax,dNfft,dscale)
    if (xct%icutv==TRUNC_SUPERCELL) then
      dNfft(1:3) = dNfft(1:3) * kgrid(1:3)
    endif
    if (mod(dNfft(3),peinf%npes) == 0) then
      Nplane = dNfft(3)/peinf%npes
    else
      Nplane = dNfft(3)/peinf%npes+1
    endif
    if (mod(dNfft(1)*dNfft(2),peinf%npes) == 0) then
      Nrod = (dNfft(1)*dNfft(2))/peinf%npes
    else
      Nrod = (dNfft(1)*dNfft(2))/peinf%npes+1
    endif
! array fftbox_2D
    rmem2=rmem2+dble(dNfft(1))*dble(dNfft(2))*dble(Nplane)*16.0d0
! array fftbox_1D
    rmem2=rmem2+dble(dNfft(3))*dble(Nrod)*16.0d0
! array dummy
!        rmem2=rmem2+dble(dNfft(1))*dble(dNfft(2))*16.0d0
! arrays dummy1 and dummy2
    rmem2=rmem2+dble(Nrod)*dble(peinf%npes+1)*16.0d0
! array inv_indx
    rmem2=rmem2+dble(Nfft(1))*dble(Nfft(2))*dble(Nfft(3))*4.0d0
  endif
  if (rmem2 .gt. rmem) rmem = rmem2
  if(peinf%inode.eq.0) then
    write(6,988) rmem/1024.0d0**2
  endif
988 format(1x,'Memory required for vcoul: ',f0.1,' MB per PE')
  
!---------------------------------
! Check how much memory is needed to allocate hamiltonian and
! eigenstates, in case of doing diagonalization. Only arrays with size
! that increase with ( # k-points )^2 are taken into account.
! Also include, bsedmatrix, bsedmt

! hmem : memory needed to store hmtrx and BSE arrays in intkernel
! dmem : memory needed to store all eigenvectors and the
! temporary arrays hbse_a_bl, evecs_r_bl (assumes that all eigenvectors
! are computed and stored)

  nmat= xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin
  nmat_co= xct%ncb_co*xct%nvb_co*xct%nspin

  facdyn = 1D0
  
! JRD: The factor of 3 here comes from bsedmattrix, bsedmatrix_loc and data in the HDF5 read area
#ifdef HDF5
  hmem = 3.0d0*8.0d0*nmat_co*nmat_co*xct%nkpt_co*facdyn
#else
  hmem = 2.0d0*8.0d0*nmat_co*nmat_co*xct%nkpt_co*facdyn
#endif  

  if (xct%ipar .eq. 1) then
    hmem = hmem + 8.0d0*nmat*peinf%nblocks*xct%nvb_fi*xct%ncb_fi*xct%nspin
    hmem = hmem + 8.0d0*(xct%nvb_fi*xct%ncb_fi*xct%nspin)**2
  else if (xct%ipar .eq. 2) then
    hmem = hmem + 8.0d0*nmat*peinf%nblocks*xct%nvb_fi*xct%nspin
    hmem = hmem + 8.0d0*(xct%nvb_fi*xct%ncb_fi*xct%nspin)*(xct%nvb_fi*xct%nspin)
  else
    hmem = hmem + 8.0d0*nmat*peinf%nblocks*xct%nspin
    hmem = hmem + 8.0d0*(xct%nvb_fi*xct%ncb_fi*xct%nspin)*xct%nspin
  endif

#ifdef CPLX
  hmem = hmem * sizeof_scalar() / 8.0d0
#endif

  if (is_diag) then 
    if (peinf%npes.eq.1) then
      dmem = 2D0*hmem + 4D0*6D0*nmat + 8D0*10D0*nmat
#ifdef CPLX
      dmem = dmem + 8D0*7D0*nmat
#endif

! MPI must be running if more than 1 processor is used...

#if defined MPI && !defined USESCALAPACK
#error ScaLAPACK is required for MPI builds.
#endif

#ifdef MPI
    else

! Estimate size of layout arrays in diagonalize
! assume mpi support

      call layout_scalapack(nmat, scal%nbl, peinf%npes, scal%nprow, scal%npcol)
      nbc = nmat/(scal%nbl*scal%npcol)
      if(mod(nmat,(scal%nbl*scal%npcol)).gt.scal%nbl*scal%npcol) nbc=nbc+1
      nbr = nmat/(scal%nbl*scal%nprow)
      if(mod(nmat,(scal%nbl*scal%nprow)).gt.scal%nbl*scal%nprow) nbr=nbr+1
      numbl_max = nbc*nbr
      
      cls = 50
      neig = nmat
      nn = max(nmat,scal%nbl,2)
      nnp = max(nmat,scal%nprow*scal%npcol+1,4)
      np0 = numroc(nn,scal%nbl,0,0,scal%nprow)
      mq0 = numroc(max(neig,scal%nbl,2),scal%nbl,0,0,scal%npcol)
#ifdef CPLX
!  use nq0 for lwork to be consistent with pcheevx, man pages wrong
      nq0 = numroc(nn,scal%nbl,0,0,scal%npcol)
      lwork = nmat + (np0+nq0+scal%nbl)*scal%nbl

      lrwork=4*nmat+max(5*nn,np0*mq0)+iceil(neig,scal%nprow*scal%npcol)*nn + &
        (cls-1)*nmat

#else
      lwork=5*nmat+max(5*nn,np0*mq0+2*scal%nbl*scal%nbl)+ &
        iceil(neig,scal%nprow*scal%npcol)*nn+(cls-1)*nmat
      lrwork=0
#endif
      dmem = 8D0*(nmat+lrwork+scal%nprow*scal%npcol) + sizeof_scalar() * (2D0*numbl_max*scal%nbl*scal%nbl+ &
       lwork+(dble(nmat)*dble(nmat)/dble(peinf%npes)))
#endif
    endif
    
  endif

  if(peinf%inode.eq.0) then
    write(6,996) hmem/0.1024d4**2
    if (is_diag) then
      write(6,997) dmem/0.1024d4**2
    endif
    write(6,*)
  endif
996 format(1x,'Memory needed to store the effective Ham. and intkernel arrays: ',&
      f0.1,' MB per PE')
997 format(1x,'Additional memory needed for evecs and diagonalization: ',f0.1,' MB per PE')

!      if ((dmem+hmem).gt.8.0*fmem) then
!        if(peinf%inode.eq.0) then
!          write(6,995)dmem/0.1024d4**2,hmem/0.1024d4**2, &
!           fmem/0.1024d4**2,fmem*8.0d0/0.1024d4**2
!        endif
!  995 format(1x,'WARNING: not enough memory to solve BSE equations.' &
!       ,/,3x,'The job may stop before completion.',/,3x,'dmem:' &
!       ,f8.1,1x,'hmem:',f8.1,1x,'fmem',f8.1,1x,'8*fmem',f8.1)

! If there is not enough memory the job may stop
! while allocating hmtrx or in some other allocation

!      endif

  if(peinf%inode.eq.0) then
    write(6,'(1x,a)') 'Distribution of kcv blocks among PE ranks:'
    rank_last = 0
    ibt_last = peinf%ibt(rank_last+1)
    do irank = 1, peinf%npes-1
      if (peinf%ibt(irank+1)/=ibt_last) then
        write(6,70) rank_last, irank-1, ibt_last
        rank_last = irank
        ibt_last = peinf%ibt(rank_last+1)
      endif
    enddo
    write(6,70) rank_last, peinf%npes-1, ibt_last
    if (peinf%nblocks*peinf%npes/=nblockstotal .and. peinf%verb_medium) then
      write(6,'(/1x,a)') 'Note: workload not evenly distributed among PEs.'
      write(6,'(1x,a/)') 'This job will not run with optimum load balance.'
    endif
  endif
70 format(3x,'- PEs ',i0,' to ',i0,' : ',i0,' blocks')

  
  POP_SUB(distrib)
  
  return
end subroutine distrib

end module distrib_m
