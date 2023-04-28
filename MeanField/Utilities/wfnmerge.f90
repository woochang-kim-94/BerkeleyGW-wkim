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

program wfnmerge

  use global_m
  use wfn_rho_vxc_io_m
  implicit none

  integer, pointer :: atyp(:)
  real(DP),allocatable :: kw(:)
  real(DP),allocatable :: tkw(:)    !< total kweights
  real(DP),allocatable :: kpt(:,:)
  real(DP),allocatable :: tkpt(:,:)  !< total list of kpoints
  integer,allocatable  :: ifmin(:,:)
  integer, allocatable :: tifmin(:,:)  !< total ifmin
  integer,allocatable  :: ifmax(:,:)
  integer, allocatable :: tifmax(:,:)  !< total ifmax
  real(DP),allocatable :: occ(:,:,:)
  real(DP),allocatable :: tocc(:,:,:)  !< total occupancy
  real(DP),allocatable :: en(:,:,:)
  real(DP),allocatable :: ten(:,:,:)   !< total energies
  integer, allocatable :: ngk(:)
  integer, allocatable :: tngk(:)      !< total number of gvectors
  real(DP),pointer :: apos(:,:)
  integer, allocatable :: gvec(:,:)
  integer, allocatable :: kptfile(:)  !< kpts per file
  integer, allocatable :: tngkmax(:)  !< ngkmax for each file
  real(DP), allocatable :: rdata(:,:)
  complex(DPC), allocatable :: cdata(:,:)
  integer :: cell_symmetry,iflavor,ns,nspinor,ng,nsym,nat,nk,nb,&
       ngkmax,kmax(3),kgrid(3),rot(3,3,48),expkpt
  real(DP) :: ecutrho,ecutwfn,celvol,recvol,al,bl,a(3,3),b(3,3),adot(3,3),&
       bdot(3,3),kshift(3),tau(3,48),thesum

  integer :: total_kgrid(3),ifil,ik,ikindex,ib
  integer :: mngkmax  !< maximum number of gvectors for all kpoints
  real(DP) :: grid_shift(3)
  character*80 :: outfile
  character*80, allocatable :: infile(:)
  real(DP), allocatable :: relweight(:)
  character*3 :: sheader
  integer :: nfil,iunit,index
  integer, allocatable :: nkfil(:)
  integer :: tnk ! total number of kpoints

! Read wfnmerge.inp
  call open_file(55,file='wfnmerge.inp',form='formatted',status='old')
  read(55,'(a80)') outfile
  read(55,*) total_kgrid(1),total_kgrid(2),total_kgrid(3)
  read(55,*) grid_shift(1),grid_shift(2),grid_shift(3)
  read(55,*) nfil
  read(55,*) expkpt  ! the expected number of kpoints
  SAFE_ALLOCATE(infile,(nfil))
  SAFE_ALLOCATE(relweight,(nfil))
  SAFE_ALLOCATE(nkfil,(nfil))
  do ifil=1,nfil
    read(55,'(a80)') infile(ifil)
  enddo
  do ifil=1,nfil
    read(55,*) relweight(ifil)
  enddo
  call close_file(55)

  write(6,*) 'Output -> ', TRUNC(outfile)

  SAFE_ALLOCATE(tkw,(expkpt))
  SAFE_ALLOCATE(tkpt,(3,expkpt))
  SAFE_ALLOCATE(tngk,(expkpt))
  SAFE_ALLOCATE(kptfile,(nfil))
  SAFE_ALLOCATE(tngkmax,(nfil))

! Open all the WFN files to be read

  call open_file(8,file=outfile,form='unformatted',status='replace')
  do ifil=1,nfil
    iunit=128+ifil
    call open_file(iunit,file=infile(ifil),form='unformatted',status='old')
  end do

! Read headers
  tnk=0
  index=1
  sheader='WFN'
  iflavor=-1
  mngkmax=-1000
  do ifil=1,nfil
    iunit=128+ifil
    call read_binary_header(iunit,sheader,iflavor,ns,ng,nsym,&
      cell_symmetry,nat,nk,nb,ngkmax,ecutrho,ecutwfn,kmax, &
      kgrid,kshift,celvol,al,a,adot,recvol,bl,b,bdot, & 
      rot,tau,atyp,apos,ngk,kw,kpt,ifmin,ifmax,en,occ,nspinor,dont_warn_kgrid=.true.)
    if (ifil.eq.1) then
      SAFE_ALLOCATE(tifmin,(expkpt,ns))
      SAFE_ALLOCATE(tifmax,(expkpt,ns))
      SAFE_ALLOCATE(ten,(nb,expkpt,ns))
      SAFE_ALLOCATE(tocc,(nb,expkpt,ns))
    end if
    tnk=tnk+nk
    kptfile(ifil)=nk
    if (ngkmax>mngkmax) mngkmax=ngkmax
    tngkmax(ifil)=ngkmax
    tkw(index:index+nk-1)=kw(1:nk)*relweight(ifil)
    tkpt(:,index:index+nk-1)=kpt(:,1:nk)
    tifmin(index:index+nk-1,:)=ifmin(1:nk,:)
    tifmax(index:index+nk-1,:)=ifmax(1:nk,:)
    ten(:,index:index+nk-1,:)=en(:,1:nk,:)
    tocc(:,index:index+nk-1,:)=occ(:,1:nk,:)
    tngk(index:index+nk-1)=ngk(1:nk)
    index=index+nk
  end do

! Renormalize the weights.
! See note in the header for more about this

  thesum=sum(tkw(1:tnk))
  do ik=1,tnk
    tkw(ik)=tkw(ik)/thesum
  end do

  write(6,*) 'Total number of kpoints found=',tnk
  if (tnk.ne.expkpt) then
    call die('Unexpected number of kpoints found!')
  end if

! Headers read. Write output header

  call write_binary_header(8,sheader,iflavor,ns,ng,nsym,&
    cell_symmetry,nat,tnk,nb,mngkmax,ecutrho,ecutwfn,kmax,&
    total_kgrid,grid_shift,celvol,al,a,adot,recvol,bl,b,bdot,&
    rot,tau,atyp,apos,tngk,tkw,tkpt,tifmin,tifmax,ten,tocc,nspinor,dont_warn_kgrid=.true.)

! Read charge density gvectors

  SAFE_ALLOCATE(gvec,(3,ng))
  do ifil=1,nfil
    iunit=128+ifil
    call read_binary_gvectors(iunit,ng,ng,gvec)
    if (ifil.eq.1) then
      call write_binary_gvectors(8,ng,ng,gvec)
    end if
  end do
  SAFE_DEALLOCATE(gvec)

! Now we loop over kpoints, get their respective gvectors, and 
! read/write data band-by-band

  SAFE_ALLOCATE(gvec,(3,mngkmax))
  if (iflavor .eq. 1) then
     SAFE_ALLOCATE(rdata,(mngkmax,ns*nspinor))
  else
     SAFE_ALLOCATE(cdata,(mngkmax,ns*nspinor))
  end if
  ikindex=0
  do ifil=1,nfil
    iunit=128+ifil
    do ik=1,kptfile(ifil)
      ikindex=ikindex+1
      write(6,*) 'Working on kpoint #', ikindex
      ! read g-vectors for current k-point
      call read_binary_gvectors(iunit,tngk(ikindex),tngkmax(ifil),gvec)
      ! write these g-vectors back to the outfile
      call write_binary_gvectors(8,tngk(ikindex),mngkmax,gvec)
      ! now read/write band data
      do ib=1,nb
        if (iflavor .eq. 1) then
           rdata(:,:) = 0d0
           call read_binary_data(iunit,tngk(ikindex),tngkmax(ifil),ns*nspinor,rdata(1:tngkmax(ifil),:))
           call write_binary_data(8,tngk(ikindex),mngkmax,ns*nspinor,rdata)
        else
           cdata(:,:) = CMPLX(0d0, 0d0)
           call read_binary_data(iunit,tngk(ikindex),tngkmax(ifil),ns*nspinor,cdata(1:tngkmax(ifil),:))
           call write_binary_data(8,tngk(ikindex),mngkmax,ns*nspinor,cdata)
        endif
      end do
    end do ! ik   
  end do ! ifil

  write(6,*) 'JOB DONE. #', ikindex

! Close files

  do ifil=1,nfil
    iunit=128+ifil
    call close_file(iunit)
  end do
  call close_file(8)

! Deallocate arrays

  if (iflavor .eq. 1) then
     SAFE_DEALLOCATE(rdata)
  else
     SAFE_DEALLOCATE(cdata)
  endif
  SAFE_DEALLOCATE(gvec)
  SAFE_DEALLOCATE(tngkmax)
  SAFE_DEALLOCATE(kptfile)
  SAFE_DEALLOCATE(tkw)
  SAFE_DEALLOCATE(tngk)
  SAFE_DEALLOCATE(tkpt)
  SAFE_DEALLOCATE(tifmin)
  SAFE_DEALLOCATE(tifmax)
  SAFE_DEALLOCATE(ten)
  SAFE_DEALLOCATE(tocc)
  SAFE_DEALLOCATE(infile)
  SAFE_DEALLOCATE(relweight)
  SAFE_DEALLOCATE(nkfil)

end program wfnmerge
