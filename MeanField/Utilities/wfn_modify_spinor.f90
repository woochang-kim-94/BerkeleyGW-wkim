!=============================================================
! Utilities:
! 
! wfnmerge  Originally by JLL. Modifications by JRD. Version compatible 
! with new wavefunction format by BDM.
!
! Merges many WFN files into one. It assumes that all input files have the 
! same number and same ordering of G-vectors for the charge density.
! The number and name of input files is read from "wfnmerge.inp", as well
! as the k-grid and the k-shift.
!
! FAQ:
! Q: Why would someone want to use this?
! A: For example, maybe to do a converged calculation one needs to include a
!     large number of unoccupied states. There may not be the resources (either
!     CPUs or wallclock) to do this all in one shot, so it may be beneficial to
!     split up a calculation of kpoints (e.g., 4x4x4 MP grid) into smaller
!     pieces and then add up the final wavefunctions.
! Q: What is the deal with kpoint weights?
! A: Quantum Espresso renormalizes the kpoint weights that you use in a 
!    calculation. If you are splitting up a MP grid into different groups of 
!    kpoints, and each group is normalized, the relative weights between
!    kpoints in different groups is no longer correct. BerkeleyGW does not 
!    make use of the kpoint weights (except for in determining whether the 
!    Fermi level is reasonable for metallic calculations). So for most uses
!    the fact that these weights are incorrect does not matter and you can
!    set the relative weights for each group of k-points to 1 in the input
!    file. If it matters to you, set the relative weights according to the
!    weights of k-points in the original MP grid and check the weights in
!    the final WFN file using wfn_rho_vxc_info.x.
!=============================================================

#include "f_defs.h"

program wfn_modify_spinor

  use global_m
  use wfn_rho_vxc_io_m
  implicit none

  integer, pointer :: atyp(:)
  real(DP),allocatable :: kw(:)
  real(DP),allocatable :: tkw(:)    !< total kweights
  real(DP),allocatable :: kw_n(:)    !< total kweights
  real(DP),allocatable :: kpt(:,:)
  real(DP),allocatable :: kpt_n(:,:)  !< total list of kpoints
  integer,allocatable  :: ifmin(:,:)
  integer, allocatable :: ifmin_n(:,:)  !< total ifmin
  integer,allocatable  :: ifmax(:,:)
  integer, allocatable :: ifmax_n(:,:)  !< total ifmax
  real(DP),allocatable :: occ(:,:,:)
  real(DP),allocatable :: occ_n(:,:,:)  !< total occupancy
  real(DP),allocatable :: en(:,:,:)
  real(DP),allocatable :: en_n(:,:,:)
  real(DP),allocatable :: ten(:,:,:)   !< total energies
  integer, allocatable :: ngk(:)
  integer, allocatable :: ngk_n(:)      !< total number of gvectors
  real(DP),pointer :: apos(:,:)
  integer, allocatable :: gvec(:,:)
  integer, allocatable :: kptfile(:)  !< kpts per file
  integer, allocatable :: tngkmax(:)  !< ngkmax for each file
  real(DP), allocatable :: rdata(:,:)
  complex(DPC), allocatable :: cdata(:,:)
  integer :: cell_symmetry,iflavor,ns,ng,nsym,nat,nk,nb,&
       ngkmax,kmax(3),kgrid(3),rot(3,3,48),expkpt,ik_init,ik_fin,nbnd_n,&
       is,ngkmax_n,ib_init,ifmax_set
  real(DP) :: ecutrho,ecutwfn,celvol,recvol,al,bl,a(3,3),b(3,3),adot(3,3),&
       bdot(3,3),kshift(3),tau(3,48),thesum

  integer :: total_kgrid(3),ifil,ik,ikindex,ib
  integer :: mngkmax  !< maximum number of gvectors for all kpoints
  real(DP) :: grid_shift(3)
  character*80 :: outfile, infile
  real(DP), allocatable :: relweight(:)
  character*3 :: sheader
  integer :: nfil,iunit,index, nspinor
  integer, allocatable :: nkfil(:)
  integer :: tnk, nk_n ! total number of kpoints

! Read wfnmerge.inp
! WK: description
  call open_file(55,file='wfn_modify.inp',form='formatted',status='old')
  read(55,'(a80)') infile     ! Input wfn file                                        
  read(55,'(a80)') outfile    ! Ouput wfn file
  read(55,*) ik_init          ! The initial index of k-points to be extracted
  read(55,*) ik_fin           ! The final index of k-points to be extracted
  read(55,*) ib_init          ! The lowest band index to be extracted (start from 1) 
  read(55,*) nbnd_n           ! The highest band index to be extracted
  read(55,*) ns               ! 'nspin' = 1 unless for a spin-polarized calculation 
  read(55,*) ifmax_set        ! 'ifmax' = # of occupied bands in the new wfn file   
  call close_file(55)

  write(6,*) 'Output -> ', TRUNC(outfile)
  write(6,*) "ib_init,nbnd_n,ik_init,ik_fin", ib_init, nbnd_n,ik_init,ik_fin
  nk_n = ik_fin - ik_init + 1
  SAFE_ALLOCATE(kw_n,(nk_n))
  SAFE_ALLOCATE(kpt_n,(3,nk_n))
  SAFE_ALLOCATE(ngk_n,(nk_n))
  SAFE_ALLOCATE(ifmin_n,(nk_n,ns))
  SAFE_ALLOCATE(ifmax_n,(nk_n,ns))
  SAFE_ALLOCATE(en_n,(nbnd_n - ib_init + 1,nk_n,ns))
  SAFE_ALLOCATE(occ_n,(nbnd_n - ib_init + 1,nk_n,ns))
  
! Open new file to be written
  call open_file(8,file=outfile,form='unformatted',status='replace')
! Open input file
  call open_file(18,file=infile,form='unformatted',status='old')
  sheader = 'WFN' 
  write(6,*) "Reading header"
  iflavor = 2
  call read_binary_header(18,sheader,iflavor,ns,ng,nsym,&
      cell_symmetry,nat,nk,nb,ngkmax,ecutrho,ecutwfn,kmax, &
      kgrid,kshift,celvol,al,a,adot,recvol,bl,b,bdot, &
      rot,tau,atyp,apos,ngk,kw,kpt,ifmin,ifmax,en,occ,nspinor,dont_warn_kgrid=.true.)
  write(6,*) "iflavor, nspinor:", iflavor, nspinor
  kw_n(:) = kw(ik_init:ik_fin)
  kpt_n(:,:) = kpt(:,ik_init:ik_fin)
  ngk_n(:) = ngk(ik_init:ik_fin)
  ifmin_n(:,:) = ifmin(ik_init:ik_fin,:)
  ifmax_n(:,:) = ifmax(ik_init:ik_fin,:)
   
  thesum=sum(kw_n(1:nk_n))
  do ik=1,nk_n
    kw_n(ik)=kw_n(ik)/thesum
  end do
  ifmax_n(:,:) = ifmax_set
  write(6,*) "New ifmin, ifmax:", ifmin_n, ifmax_n
  
  do is=1,ns 
    do ik=1,nk
      if (ik.ge.ik_init.and.ik.le.ik_fin) then
        en_n(:,ik-ik_init+1,is) = en(ib_init:nbnd_n,ik,is)
        occ_n(:,ik-ik_init+1,is) = occ(ib_init:nbnd_n,ik,is)
      endif
    enddo
  enddo
  write(6,*) "Eigenvalues in new WFN:",en_n
  ngkmax_n = maxval(ngk_n(:))
  call write_binary_header(8,sheader,iflavor,ns,ng,nsym,&
    cell_symmetry,nat,nk_n,nbnd_n - ib_init + 1,ngkmax_n,ecutrho,ecutwfn,kmax,&
    kgrid,kshift,celvol,al,a,adot,recvol,bl,b,bdot,&
    rot,tau,atyp,apos,ngk_n,kw_n,kpt_n,ifmin_n,ifmax_n,en_n,occ_n,nspinor,dont_warn_kgrid=.true.)  
  SAFE_ALLOCATE(gvec,(3,ng))
  call read_binary_gvectors(18,ng,ng,gvec)
  call write_binary_gvectors(8,ng,ng,gvec)
  SAFE_DEALLOCATE(gvec)

  SAFE_ALLOCATE(cdata,(ngkmax,ns*nspinor))
  SAFE_ALLOCATE(gvec,(3,ngkmax))
  
  do ik=1,nk
    if (ik.ge.ik_init.and.ik.le.ik_fin) then
      ! read g-vectors for current k-point
      call read_binary_gvectors(18,ngk(ik),ngkmax,gvec)
      ! write these g-vectors back to the outfile
      call write_binary_gvectors(8,ngk(ik),ngkmax_n,gvec)
      ! now read/write band data
      write(6,*) ngk(ik), "<", ngkmax_n
      do ib=1,ib_init-1
        write(6,*) "Reading ik, ib:", ik, ib
        cdata(:,:) = CMPLX(0d0, 0d0)
        call read_binary_data(18,ngk(ik),ngkmax,ns*nspinor,cdata(1:ngkmax,:))
      enddo
      do ib=ib_init,nbnd_n
        write(6,*) "Writing ik, ib:", ik, ib
        cdata(:,:) = CMPLX(0d0, 0d0)
        call read_binary_data(18,ngk(ik),ngkmax,ns*nspinor,cdata(1:ngkmax,:))
        call write_binary_data(8,ngk(ik),ngkmax_n,ns*nspinor,cdata(1:ngkmax_n,:))
      enddo
      do ib=nbnd_n+1,nb
        write(6,*) "Reading ik, ib:", ik, ib
        cdata(:,:) = CMPLX(0d0, 0d0)
        call read_binary_data(18,ngk(ik),ngkmax,ns*nspinor,cdata(1:ngkmax,:))
      enddo
    endif
  enddo
  write(6,*) "HERE"
  SAFE_DEALLOCATE(cdata)
  SAFE_DEALLOCATE(kw_n)
  SAFE_DEALLOCATE(ngk_n)
  SAFE_DEALLOCATE(kpt_n)
  SAFE_DEALLOCATE(ifmin_n)
  SAFE_DEALLOCATE(ifmax_n)
  SAFE_DEALLOCATE(en_n)
  SAFE_DEALLOCATE(occ_n)
  call close_file(8)
  call close_file(18)

end program wfn_modify_spinor
