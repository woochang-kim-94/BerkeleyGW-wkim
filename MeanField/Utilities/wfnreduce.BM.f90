!=============================================================================
!
! Utilities:
!
! (1) wfnreduce.BM     PD has work on this in wfnreduce.f90 that is under
!                      devel. This file is a separate attempt, but based 
!                      more along the lines of wfnmerge.f90. Sadly, neither
!                      version completely works at the moment. 
!     Reduces the cutoff of WFN files, Real/Complex flavor
!
!==============================================================================

! NOTES ON THIS UTILITY:
!  =====IMPORTANT======= 
!  1). This utility currently does NOT work. It *almost* works and gives
! answers that are about 30 meV off from what it should give, but that is 
! based on tests on a single system and there's no guarantee it will not be 
! larger in others. This needs to be be fixed before anyone uses this. I 
! do not know why it currently does not work. Please help fix or email me 
! with the problem. wfnreduce is a very important utility. For most systems 
! one is likely to waste a huge amount of time with the excessive wavefunction 
! cutoffs that are needed to produce enough empty states to converge the 
! calculations. 
! ======OTHER TODOS=========
!  2). Use kinetic energy routine in Common/ to calculate kinetic energy 
!     (development of this was done in 1.0.x where this routine has not yet 
!      been added. 
!  3). Be nice if FFT dimensions did not have to be manually specified. 
!  4). More human-readable output and options.
!  5). Improved estimate of energy savings (see bottom) and relation to 
!      feel-good quantities like trees and pandas. 
!  6). The ability to reduce the RHO file as well would be nice, but not 
!      critical since this file can always be regenerated pretty rapidly.
!
! =======NOTES===========
!  -If you want to run this just copy it over the other wfnreduce.f90 and 
!       make as normal.
!---------------------------------------------------------------------



#include "f_defs.h"

program wfnreduce_BM

  use global_m
  use wfn_rho_vxc_io_m

  implicit none


  integer, pointer :: atyp(:)
  real(DP),pointer :: kw(:)

  real(DP), allocatable :: rdatatemp(:,:),rdata(:,:)
  complex(DPC), allocatable :: cdatatemp(:,:),cdata(:,:)
  integer :: cell_symmetry,iflavor,ns,ng, newng,nsym,nat,nk,nb,&
       ngkmax,newngkmax,kmax(3),newfft(3),kgrid(3),rot(3,3,48),expkpt
  real(DP) :: ecutrho,ecutwfn,celvol,recvol,al,bl,a(3,3),b(3,3),adot(3,3),&
       bdot(3,3),kshift(3),tau(3,48),thesum, newecutrho
  real(DP), pointer :: kpt(:,:)
  integer, pointer :: ifmax(:,:)
  integer,pointer :: ifmin(:,:)
  real(DP),pointer :: occ(:,:,:)
  real(DP),pointer :: en(:,:,:)
  integer, pointer :: ngk(:)
  integer, pointer :: newngk(:) 
  real(DP),pointer :: apos(:,:)

  integer, allocatable :: gvecCD(:,:)
  integer, allocatable :: gvecK(:,:,:)
  integer, allocatable :: newgvecK(:,:,:)
  integer, allocatable :: gvecKtemp(:,:)
  integer, allocatable :: gmapping(:,:) 


  character*3 :: sheader, scutoff,sfft1,sfft2,sfft3
  character*7 :: sflavor, sbandcut
  character*256 :: infile, outfile, usage
  integer :: ik, ib, ig,newig
  integer :: nargs, bandcut
  real(DP) :: ecutw
  real(DP), allocatable :: ekin(:)  ! kinetic energies

  usage = 'Usage: wfnreduce.x newCutoff(Ry) lastBand fft1 fft2 fft3 infile outfile'


! Get file names from command-line arguments

  nargs = command_argument_count()

  if (nargs .ne. 7) then
    call die(usage)
  endif

  call get_command_argument(1, scutoff)
  call get_command_argument(2, sbandcut)
  call get_command_argument(3,sfft1)
  call get_command_argument(4,sfft2)
  call get_command_argument(5,sfft3)
  call get_command_argument(6, infile)
  call get_command_argument(7, outfile)

  read(scutoff, *) ecutw
  read(sbandcut, *) bandcut
  read(sfft1,*) newfft(1)
  read(sfft2,*) newfft(2)
  read(sfft3,*) newfft(3)
  newecutrho = ecutw * 4




  call open_file(unit=7,file=TRUNC(infile),form='unformatted',status='old')


  ! read headers
  iflavor=-1 
  sheader='WFN'
  call read_binary_header(7,sheader,iflavor,ns,ng,nsym,&
       cell_symmetry,nat,nk,nb,ngkmax,ecutrho,ecutwfn,kmax, &
       kgrid,kshift,celvol,al,a,adot,recvol,bl,b,bdot, & 
       rot,tau,atyp,apos,ngk,kw,kpt,ifmin,ifmax,en,occ,dont_warn_kgrid=.true.)


! read charge density gvectors
  SAFE_ALLOCATE(gvecCD,(3,ng))
  SAFE_ALLOCATE(ekin,(ng))
  call read_binary_gvectors(7,ng,ng,gvecCD)
  ! compute kinetic energies
  newng=0
  DO ig =1,ng
     ekin(ig)=DOT_PRODUCT(gvecCD(:,ig),MATMUL(bdot,gvecCD(:,ig)))
     if (ekin(ig)< newecutrho) then 
        newng=newng+1
     end if
  END DO

  if (iflavor .eq. 1) then
     SAFE_ALLOCATE(rdatatemp,(ngkmax,ns))
     SAFE_ALLOCATE(rdata,(ngkmax,ns))
  else
     SAFE_ALLOCATE(cdatatemp,(ngkmax,ns))
     SAFE_ALLOCATE(cdata,(ngkmax,ns))
  end if



  
  write(*,*) 'before loop over kpoints'
  ! now we loop over kpoints, and get their respective gvectors so that we 
  ! can trim those lists. Once we have all the info we will just read through
  ! the file again and write out the reduced wavefunctions
  SAFE_ALLOCATE(gvecK,(nk,3,ngkmax))  ! gvectors for each k, full list
  SAFE_ALLOCATE(gmapping,(nk,ngkmax))   ! the mapping between ig in the 
                                        ! full list and the ig in the
                                        ! reduced list
  SAFE_ALLOCATE(newgvecK,(nk,3,ngkmax)) ! gvectors for each k, reduced list
  SAFE_ALLOCATE(newngk,(nk))  ! new array for # of G vectors for each k
  newngk(:)=ZERO
  gmapping(:,:)=ZERO
  newgvecK(:,:,:)=ZERO
  ekin(:)=ZERO
  newngkmax=0

  do ik=1,nk
     newig=0
     ! read g-vectors for current k-point
     call read_binary_gvectors(7,ngk(ik),ngkmax,gvecK(ik,:,:))
     ! determine new arrays by seeing if the kinetic energy is below the 
     ! cutoff
     DO ig=1,ngkmax
        ekin(ig)=DOT_PRODUCT(gvecK(ik,:,ig)+kpt(:,ik),MATMUL(bdot,gvecK(ik,:,ig)+kpt(:,ik)))
        if (ekin(ig)< ecutw) then
           newig=newig+1  ! pronounced 'new ig'
           newgvecK(ik,:,newig)=gvecK(ik,:,ig)
           gmapping(ik,newig)=ig
           newngk(ik)=newngk(ik)+1
        end if
     end do
     if (newngkmax < newngk(ik)) then
        newngkmax = newngk(ik)
     end if
     ! now we "read" in the wavefunctions, but just to move correctly through
     ! the file. we don't care about no wavefunctions right na na na. 
     do ib=1,nb
        if (iflavor .eq. 1) then 
           rdata(:,:) =ZERO
           call read_binary_data(7,ngk(ik),ngkmax,ns,rdata(1:ngkmax,:))
        else
           cdata(:,:)=ZERO
           call read_binary_data(7,ngk(ik),ngkmax,ns,cdata(1:ngkmax,:))
        end if
     end do ! ib 
  end do

  write(*,*) 'newngk=',newngk
  ! now let us close file and reopen again, this time using our new 
  ! knowledge of the WFN world
  call close_file(7)
 
  call open_file(unit=7,file=TRUNC(infile),form='unformatted',status='old')
  call open_file(unit=8,file=TRUNC(outfile),form='unformatted',status='replace')

  ! read/write headers                                                                
  iflavor=-1
  sheader='WFN'
  call read_binary_header(7,sheader,iflavor,ns,ng,nsym,&
       cell_symmetry,nat,nk,nb,ngkmax,ecutrho,ecutwfn,kmax, &
       kgrid,kshift,celvol,al,a,adot,recvol,bl,b,bdot, &
       rot,tau,atyp,apos,ngk,kw,kpt,ifmin,ifmax,en,occ,dont_warn_kgrid=.true.)
! write new header 
    call write_binary_header(8,sheader,iflavor,ns,newng,nsym,&
    cell_symmetry,nat,nk,nb,newngkmax,newecutrho,ecutw,newfft,&
    kgrid,kshift,celvol,al,a,adot,recvol,bl,b,bdot,&
    rot,tau,atyp,apos,newngk,kw,kpt,ifmin,ifmax,en,occ,dont_warn_kgrid=.true.)

! read/write charge density gvectors  
    call read_binary_gvectors(7,ng,ng,gvecCD)
    call write_binary_gvectors(8,newng,newng,gvecCD(:,1:newng))
! now we loop over kpoints, read/write gvectors, and read/write data band-by
! band

  do ik=1,nk
     ! read g-vectors for current k-point                                       
     call read_binary_gvectors(7,ngk(ik),ngkmax,gvecK(ik,:,:))
!     write(*,*) ngk,ngkmax,size(gvecK(ik,:,:))
     call write_binary_gvectors(8,newngk(ik),newngkmax,newgvecK(ik,:,1:newngkmax))
     do ib=1,nb
        if (iflavor .eq. 1) then 
           rdata(:,:) =ZERO
           rdatatemp(:,:)=ZERO
           call read_binary_data(7,ngk(ik),ngkmax,ns,rdatatemp(1:ngkmax,:))
           do ig=1,newngkmax
              rdata(ig,:)=rdatatemp(gmapping(ik,ig),:)
           end do
           call write_binary_data(8,newngk(ik),newngkmax,ns,rdata(1:newngkmax,:))
        else
           cdata(:,:)=ZERO
           cdatatemp(:,:)=ZERO
           call read_binary_data(7,ngk(ik),ngkmax,ns,cdatatemp(1:ngkmax,:))
           ! now construct cdata from cdatatemp using gmapping
           do ig=1,newngkmax
              cdata(ig,:)=cdatatemp(gmapping(ik,ig),:)
           end do
           call write_binary_data(8,newngk(ik),newngkmax,ns,cdata(1:newngkmax,:))
        end if
     end do ! ib 

  end do    ! ik

  call close_file(8)

! TODO: write out estimate of time savings in terms of N^4 scaling of BGW. 
! Relating this quantity to some number of trees would be helpful. Relating
! it from trees to some number of pandas saved would be especially 
! motivating, but arguably includes some level of deception.
  write(*,*) 'Thanks for using wfnreduce. Your use of this script helps conserve energy.'

  if (iflavor .eq. 1) then
     SAFE_DEALLOCATE(rdata)
     SAFE_DEALLOCATE(rdatatemp)
  else
     SAFE_DEALLOCATE(cdata)
     SAFE_DEALLOCATE(cdatatemp)
  endif


  SAFE_DEALLOCATE(gmapping)
  SAFE_DEALLOCATE(newgvecK)
  SAFE_DEALLOCATE(newngk)
  SAFE_DEALLOCATE(gvecCD)    
  SAFE_DEALLOCATE(gvecK)
  SAFE_DEALLOCATE(ekin)
end program wfnreduce_BM
