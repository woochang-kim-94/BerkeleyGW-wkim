!===============================================================================
!
! Module:
!
! abi2bgw_m             Written by  Tonatiuh Rangel         Last Modified 07/30/2014    
!
! Main module used by ABI2BGW: reads ABINIT output files and 
! writes BGW input files.
!
! Beta version.
!
! Developers:
!   Tonatiuh Rangel, Berkeley CA, trangel@lbl.gov
!   Felipe Jornada, Berkeley CA
!
!===============================================================================


#include "f_defs.h"
module abi2bgw_m

 use global_m
 use wfn_rho_vxc_io_m
 use read_abi_m, only: abi_hdr_type,abi_geom_type,abi_basis_type,abi_kpt_type,abi_elec_type
 use reciprocal_lattice_m, only: reciprocal_lattice
 use print_error_m, only : print_error

 implicit none

 private
 public fill_kp_object
 public show_kp_object
 public fill_syms_object
 public show_syms_object
 public fill_crys_object
 public show_crys_object
 public fill_gvec_object
 public show_gvec_object
 public read_write_kpt_loop
 public update_symmetries_from_file

contains
!
subroutine update_symmetries_from_file(syms)
 implicit none
 type(symmetry),intent(inout)::syms
!
 integer,parameter:: iunit=1,ounit=2,tmpunit=3
 real(DP),parameter::zero=0.d0
 integer::iost,isym
 character(80)::aa,infile

 PUSH_SUB(update_symmetries_from_file)

 infile="symrel.in"
 write(*,*)'Reading symmetries from file ',trim(infile)
 call open_file(tmpunit, FILE=trim(infile), iostat=iost, STATUS="old", FORM="FORMATTED")

 read(tmpunit,*) aa, syms%ntran
 read(tmpunit,*) aa
 do isym=1,syms%ntran
   read(tmpunit,*)syms%mtrx(:,:,isym)
!  Assuming all translations are zero
   syms%tnp(:,isym)=zero
 end do

 call close_file(tmpunit)

 POP_SUB(update_symmetries_from_file)

end subroutine update_symmetries_from_file

!
!Reads WFN coefficients from ABINIT output and writes them into BGW format.
subroutine read_write_kpt_loop(abi_hdr,abi_basis,iunit,ounit,outformat,save_cg,pert_wfn)
 implicit none
 logical, intent(in)::outformat
 logical, intent(in)::save_cg  !if true, saves all cg into memory
 integer,intent(in)::iunit,ounit
 type( abi_hdr_type),intent(inout):: abi_hdr
 type( abi_basis_type),intent(inout):: abi_basis
 logical, intent(in) :: pert_wfn ! ZL adds
!
 integer:: ipw,jpw,ik,isppol,iband,ii,ij,ispin, icheck ! ZL adds icheck
 integer:: npw,mpw,nspinor,nkpt,nbandk,mband,nspin,nsppol
 real(DP):: occ_factor
 real(DP),allocatable::occk(:),cg(:,:),eigen(:),eigentmp(:)
 real(DP),allocatable::cg_save(:,:,:,:,:)
 complex(DPC),allocatable :: zwfn(:,:)
 integer, allocatable, target  :: red_coord_pw_k(:, :)
 integer, allocatable :: red_coord_pw_save(:,:,:,:)
 !
 PUSH_SUB(read_write_kpt_loop)
!
 !Safety check: save_cg only for nsppol>1
 if ( save_cg ) then
   if ( abi_hdr%number_of_spin_polarizations == 1 .or. abi_hdr%number_of_spinor_components > 1 ) &
       & call print_error('abi2bgw','save_cg only for nsppol>1',1)
 end if
 !Pass object variables to shorter variables:
 nsppol =abi_hdr%number_of_spin_polarizations
 nspinor=abi_hdr%number_of_spinor_components
 nkpt   =abi_hdr%number_of_kpoints
 nspin  = nsppol * nspinor
 mpw    = abi_hdr%max_number_of_coefficients
 mband  = abi_hdr%max_number_of_states
 SAFE_ALLOCATE(zwfn,(mpw, nspin))
 SAFE_ALLOCATE(red_coord_pw_k,(3, mpw))
 SAFE_ALLOCATE(cg,(2,nspinor*mpw))
 SAFE_ALLOCATE(eigen,(mband))
 SAFE_ALLOCATE(occk,(mband))
 if ( save_cg ) then
   SAFE_ALLOCATE(cg_save, (2,mpw,mband,nkpt,nsppol))
   SAFE_ALLOCATE(red_coord_pw_save, (3,mpw,nkpt,nsppol))
 end if
!
!Need to rescale spins to [0-1] for BGW
 if(nspin==1) then
   occ_factor=0.5d0
 else
   occ_factor=1.0d0
 end if
!
! big kpt loop to read G vectors, wave coefficients, 
! eigenvalues and occupations
!
 DO isppol=1, nsppol
   DO ik=1, nkpt
     write(*,*)'spin ',isppol,'kpt ',ik,'/',nkpt
!     kpttmp(:)=red_coord_kpt(:,ik)
     READ(iunit)npw,nspinor,nbandk      ! for each k point
     !check that the number of plane waves for this k-point
     !corresponds to the value already stored in number_of coefficients
     if(npw .ne. abi_basis%number_of_coefficients(ik))&
         & call print_error('abi2bgw','npw read is different from the stored value',1)
     if(nspinor .ne. abi_hdr%number_of_spinor_components)&
         & call print_error('abi2bgw','npw read is different from the stored value',1)
     if(nbandk .ne. abi_basis%number_of_states(ik,isppol))&
         & call print_error('abi2bgw','npw read is different from the stored value',1)

!    Read G vectors
     READ(iunit)red_coord_pw_k(1:3,1:npw)
     ! ZL check
     !do icheck=1, npw
     !  write(*,*) "check Gvec", red_coord_pw_k(1:3, icheck)
     !end do
!
     if ( save_cg ) then
       red_coord_pw_save(1:3,1:npw,ik,isppol)=red_coord_pw_k(1:3,1:npw)
     else
       call write_gvectors(ounit,outformat,npw,mpw,red_coord_pw_k)
     end if
 
!    Read eigenvalues and occupations
!    Here we just read 1 element, since they are not used here:
!    READ(iunit)eigen(1:nbandk),occk(1:nbandk)
     ! ZL modifies to be compatible with both ground-state and first-order perturbed wfns
     if (.not. pert_wfn) then
       READ(iunit) ! eigen(1),occk(1) ! ZL: they are not used, so not read at all
     end if
!
!    reads the wavefunction coefficients for each band
     DO iband=1,nbandk
       if (pert_wfn) then
         READ(iunit) ! ZL: for pert_wfn, eigen is written in this loop.
         ! No need to read these numbers, but need to skip
       end if
       READ(iunit)((cg(ii,ij),ii=1,2),&
         &ij=1,npw*nspinor)
       if ( save_cg ) then
         cg_save(1:2,1:npw,iband,ik,isppol)=cg(1:2,1:npw)
       else
         jpw=0
         do ispin=1,nspinor
           do ipw=1,npw
             jpw=jpw+1
             zwfn(ipw,ispin)=CMPLX(cg(1,jpw),cg(2,jpw))
             ! ZL: check
             !write(*,*) "check coeff:", zwfn(ipw,ispin)
           end do
         end do
         call write_complex_data(ounit,outformat,npw,mpw,nspinor,zwfn)
       end if
     END DO!iband
   END DO !ik
 END DO !dims%number_of_spin_polarizations
!
!For the option save_cg, write cg data:
 if ( save_cg ) then
   DO ik=1, nkpt
     npw=abi_basis%number_of_coefficients(ik)
     nbandk=abi_basis%number_of_states(ik,1)
     call write_gvectors(ounit,outformat,npw,mpw,red_coord_pw_save(:,:,ik,1))
     DO iband=1,nbandk
       do ispin=1,nspin
         do ipw=1,npw
           zwfn(ipw,ispin)=CMPLX(cg_save(1,ipw,iband,ik,ispin),cg_save(2,ipw,iband,ik,ispin))
         end do
       end do
       call write_complex_data(ounit,outformat,npw,mpw,nspin,zwfn)
     END DO !iband
   END DO !ikpt
 end if

 SAFE_DEALLOCATE(zwfn)
 SAFE_DEALLOCATE(cg)
 SAFE_DEALLOCATE(red_coord_pw_k)
 SAFE_DEALLOCATE(eigen)
 SAFE_DEALLOCATE(occk)
 if (save_cg) then
   SAFE_DEALLOCATE(cg_save)
   SAFE_DEALLOCATE(red_coord_pw_save)
 end if
!
 POP_SUB(read_write_kpt_loop)
end subroutine read_write_kpt_loop
!

!Fills BGW gvec object from ABINIT variables:
subroutine fill_gvec_object(gvec,kp,crys,abi_hdr,red_coord_pw_rho)
 implicit none
 type(gspace),intent(inout)::gvec
 type(kpoints),intent(in)::kp
 type(crystal),intent(in)::crys
 type(abi_hdr_type),intent(in)::abi_hdr
 integer, pointer  :: red_coord_pw_rho(:, :)
!
 real(DP), PARAMETER :: ha2ry=2.0d0
 real(DP), parameter :: four=4.0d0
 real(DP), parameter :: five=5.0d0
 real(DP), parameter :: six=6.0d0
 real(DP), parameter :: eight=8.0d0
 real(DP), parameter :: ten=10.0d0
!
 PUSH_SUB(fill_gvec_object)
!
 !real(DP) :: ecutrho !< charge-density cutoff, in Ry
 gvec%ecutrho = abi_hdr%energy_cutoff*ha2ry*four
 !integer :: nFFTgridpts !< number in FFT grid = product(FFTgrid(1:3))
 gvec%nFFTgridpts=product(abi_hdr%fft_grid(1:3))
 !integer :: FFTgrid(3)  !< gsm: FFTgrid is the size of the FFT grid, not the maximum G-vector   
 gvec%FFTgrid(1:3)=abi_hdr%fft_grid(1:3)
 !
 !integer :: ng       !< number of G-vectors
 call generate_gsphere(gvec,crys,red_coord_pw_rho)

 !integer, pointer :: components(:,:) !< the G-vectors, in units of 2 pi / a
 gvec%components => red_coord_pw_rho
 !integer, pointer :: index_vec(:) ! mapping to FFT grid
 !real(DP), pointer :: ekin(:) !< kinetic energy of G-vectors

 POP_SUB(fill_gvec_object)

end subroutine fill_gvec_object

subroutine generate_gsphere(gvec,crys,red_coord_pw_rho)
 implicit none
 type(gspace),intent(inout)::gvec
 type(crystal),intent(in)::crys
 integer, pointer  :: red_coord_pw_rho(:, :)
!
 integer,parameter:: iunit=1,ounit=2,tmpunit=3
 character*80 :: gsphere_file
 real(DP), parameter    :: two   = 2.0d0
 real(DP), parameter    :: four   = 4.0d0
 real(DP), parameter    :: twopi = 6.2831853071795864769252867665590D0
 integer::ig1,ig2,ig3,ng1,ng2,ng3,iost
 integer:: n1a,n1b,n2a,n2b,n3a,n3b
 real(DP)::ecut_,ecutoff,bvec(3,3),gmet(3,3)
 real(DP)::vecx,vecy,vecz
 integer ::ix,iy,iz,ig

 PUSH_SUB(generate_gsphere)

 bvec=crys%bvec*crys%blat
 gmet=crys%bdot !/twopi**2
 ecutoff=gvec%ecutrho/two  !in Hartree
 gsphere_file="gsphere_bgw.out"

 write(*,*)'#Generated list of g-vectors for density\/potential writen in file: ',trim(gsphere_file)
 ng1=gvec%FFTgrid(1);  ng2=gvec%FFTgrid(2);  ng3=gvec%FFTgrid(3); 
!
!Generating vectors in range [-n/2,n/2]
 n1a=ceiling(real(-ng1,DP)/two)
 n1b=floor(real(ng1,DP)/two)
 n2a=ceiling(real(-ng2,DP)/two)
 n2b=floor(real(ng2,DP)/two)
 n3a=ceiling(real(-ng3,DP)/two)
 n3b=floor(real(ng3,DP)/two)
!Check gvectors:
 if(2*n1b >= ng1) n1b=n1b-1
 if(2*n2b >= ng2) n2b=n2b-1
 if(2*n3b >= ng3) n3b=n3b-1
 if(-2*n1a >= ng1) n1a=n1a+1
 if(-2*n2a >= ng2) n2a=n2a+1
 if(-2*n3a >= ng3) n3a=n3a+1


 ig=0
 do ig3=n3a,n3b
   do ig2=n2a,n2b
     do ig1=n1a,n1b
       vecx=bvec(1,1)*real(ig1,DP) + bvec(1,2)*real(ig2,DP) + bvec(1,3)*real(ig3,DP) 
       vecy=bvec(2,1)*real(ig1,DP) + bvec(2,2)*real(ig2,DP) + bvec(2,3)*real(ig3,DP)
       vecz=bvec(3,1)*real(ig1,DP) + bvec(3,2)*real(ig2,DP) + bvec(3,3)*real(ig3,DP)

!      ecut= \hbar^2/2 |G|^2
       ecut_=(vecx**2+vecy**2+vecz**2)/two

       if( ecut_ < ecutoff)    ig=ig+1
     end do
   end do
 end do 
!
!Do this again but this time save the vectors in an array:
 call open_file(tmpunit, FILE=gsphere_file, iostat=iost, STATUS="unknown", FORM="FORMATTED")

 gvec%ng=ig
 write(tmpunit,'(1x,"#ng=",1x,i12)')gvec%ng
 write(tmpunit,'(1x,"g1",1x,"g2",1x,"g3",3x,"gx (bohrs)",1x,"gy (bohrs)",1x,"gz (bohrs)",3x,"ecut (Ha)")')
 SAFE_ALLOCATE(red_coord_pw_rho,(3,gvec%ng))
 ig=0
 do ig3=n3a,n3b
   do ig2=n2a,n2b
     do ig1=n1a,n1b
       vecx=bvec(1,1)*real(ig1,DP) + bvec(1,2)*real(ig2,DP) + bvec(1,3)*real(ig3,DP) 
       vecy=bvec(2,1)*real(ig1,DP) + bvec(2,2)*real(ig2,DP) + bvec(2,3)*real(ig3,DP)
       vecz=bvec(3,1)*real(ig1,DP) + bvec(3,2)*real(ig2,DP) + bvec(3,3)*real(ig3,DP)

!      ecut= \hbar^2/2 |G|^2
       ecut_=(vecx**2+vecy**2+vecz**2)/two

       if( ecut_ < ecutoff) then
         ig=ig+1
         write(tmpunit,'(3(1x,i4),3(1x,f12.6),2x,1(1x,f12.6))')ig1,ig2,ig3,vecx,vecy,vecz,ecut_
         red_coord_pw_rho(:,ig)=[ig1,ig2,ig3]
       end if
     end do
   end do
 end do 
 if (ig /= gvec%ng) call print_error("generate_gsphere","Inconsistency in number of g-vectors",1)

 call close_file(tmpunit)
 !
 POP_SUB(generate_gsphere)

end subroutine generate_gsphere
!
subroutine show_gvec_object(gvec)
!Fills BGW gvec object from ABINIT variables:
  implicit none
  type(gspace),intent(in)::gvec
  integer::ipw,jpw
  logical::lwrite_gvectors=.false.
!
  PUSH_SUB(show_gvec_object)
!
  write(*,*)
  write(*,'("Showing components of gvec object:")')
  write(*,'("Number of G-vectors",1x,i12)') gvec%ng
  write(*,'("Number of FFT points",1x,i6)')gvec%nFFTgridpts
  write(*,'("Charge-density cutoff in Ry",1x,f12.6)') gvec%ecutrho
  if( lwrite_gvectors) then
 !  integer, pointer :: components(:,:) !< the G-vectors, in units of 2 pi / a
    write(300,'("G-vectors(ik) in reduced units")')
    do ipw=1,gvec%ng
      write(300,'(3(1x,i6))',advance='yes')gvec%components(:,ipw)
    end do
  end if
 !integer :: FFTgrid(3)  !< gsm: FFTgrid is the size of the FFT grid, not the maximum G-vector   
 write(*,'("FFTgrid",1x,3(1x,i6))') gvec%FFTgrid(1:3)
 !integer, pointer :: index_vec(:) ! mapping to FFT grid
 !real(DP), pointer :: ekin(:) !< kinetic energy of G-vectors

  POP_SUB(show_gvec_object)
end subroutine show_gvec_object

subroutine fill_crys_object(crys,abi_hdr)
  implicit none
  type(crystal),intent(out)::crys
  type(abi_hdr_type),intent(in)::abi_hdr

  real(DP), parameter    :: twopi = 6.2831853071795864769252867665590D0
  real(DP), parameter :: one=1._DP
  integer, parameter:: iout=1 !if < 0 then do not write messages
  integer::iat
  real(DP):: blat, bdot(3,3), recvol,celvol
  real(DP):: xcart(3),xred(3)
  real(DP):: bvec(3,3) !reciprocal space vectors in blat units
  real(DP):: gprimd(3,3) !bvec*blat/2pi (a.u.)
  real(DP):: rprimd(3,3) !avec/alat (a.u.)
  real(DP):: gmet(3,3)   !bdot/(2pi)**2 (a.u)
  real(DP):: rmet (3,3)  !adot (a.u)

  PUSH_SUB(fill_crys_object)

!  call read_crystal_file(crys)
!
  crys%avec=abi_hdr%primitive_vectors
  crys%alat=one
  rprimd=crys%alat*crys%avec

  call reciprocal_lattice(gmet,gprimd,rmet,rprimd,recvol,celvol)

!Check consistency in unit cell vol:
!  if(abs(crys%celvol-celvol)>1.0d-3 ) then
!   write(*,*)'Unit cell vol from lattice parameters: ',celvol
!   write(*,*)'Not equal from celvol in crystal file: ',crys%celvol
!   call die('')
!  end if

! Convert reciprocal space quantities to blat units
! Twopi factors enter due to a distinct definition of 
! reciprocal space vectors:
!  e.g., b1 = 2pi ( a2 x a3)/(a1 . (a2 x a3)) or
!        b1 = ( a2 x a3)/(a1 . (a2 x a3))
  blat=twopi/crys%alat
  bvec=gprimd*crys%alat 
  recvol=recvol*twopi**3
  bdot=gmet*twopi**2

!Read ABINIT output and print a file with all of these:
!    real(DP) :: celvol !< cell volume in real space (a.u.)
  crys%celvol=celvol
!    real(DP) :: recvol !< cell volume in reciprocal space (a.u.)
  crys%recvol=recvol
!    real(DP) :: alat !< lattice constant in real space (a.u.)
!  crys%alat=crys_alat
!    real(DP) :: blat !< lattice constant in reciprocal space (a.u.)
  crys%blat=blat
!    real(DP) :: avec(3,3) !< lattice vectors in real space (alat)
!  crys%avec=crys%avec
!    real(DP) :: bvec(3,3) !< lattice vectors in reciprocal space (blat)
  crys%bvec=bvec
!    real(DP) :: adot(3,3) !< metric tensor in real space (a.u.)
  crys%adot=rmet
!    real(DP) :: bdot(3,3) !< metric tensor in reciprocal space (a.u.)
  crys%bdot=bdot
!    integer :: nat !< number of atoms
!  crys%nat=crys_nat
  crys%nat=abi_hdr%number_of_atoms
!    integer, pointer :: atyp(:) !< atomic species, atyp(1:nat)
!  crys%atyp=>crys_atyp
  SAFE_ALLOCATE(crys%atyp,(crys%nat))
  crys%atyp=abi_hdr%znucl
!    real(DP), pointer :: apos(:,:) !< atomic positions, apos(1:3,1:nat) (alat)
  SAFE_ALLOCATE(crys%apos,(3,crys%nat))
!
! convert xred to xcart
  do iat=1,crys%nat
    xred=abi_hdr%xred(:,iat)
    xcart=xred(:)*crys%avec(:,1)+&
&         xred(:)*crys%avec(:,2)+&
&         xred(:)*crys%avec(:,3)  
    crys%apos(:,iat)=xcart 
  end do
!  crys%apos=>crys_apos
  POP_SUB(fill_crys_object)
!
end subroutine fill_crys_object

subroutine show_crys_object(crys)
 implicit none
 type(crystal),intent(in)::crys
!
 integer::ii
!
 PUSH_SUB(show_crys_object)

 write(*,*)"*************************"
 write(*,'("Components on crys object")')
!    real(DP) :: celvol !< cell volume in real space (a.u.)
 write(*,'("celvol",1x,f20.12)')crys%celvol
!    real(DP) :: recvol !< cell volume in reciprocal space (a.u.)
 write(*,'("recvol",1x,f20.12)')crys%recvol
!    real(DP) :: alat !< lattice constant in real space (a.u.)
 write(*,'("alat",1x,f20.12)')crys%alat
!    real(DP) :: blat !< lattice constant in reciprocal space (a.u.)
 write(*,'("blat= 2 x pi / alat",1x,f20.12)')crys%blat
!    real(DP) :: avec(3,3) !< lattice vectors in real space (alat)
 write(*,'("avec \(alat\)")')
 do ii=1, 3
  write(*,'(3(1x,f20.12))')crys%avec(:,ii)
 end do
!    real(DP) :: bvec(3,3) !< lattice vectors in reciprocal space (blat)
 write(*,'("bvec \(blat\)")')
 do ii=1, 3
  write(*,'(3(1x,f20.12))')crys%bvec(:,ii)
 end do
!    real(DP) :: adot(3,3) !< metric tensor in real space (a.u.)
 write(*,'("adot \(a. u.\)")')
 do ii=1, 3
  write(*,'(3(1x,f20.12))')crys%adot(:,ii)
 end do
!    real(DP) :: bdot(3,3) !< metric tensor in reciprocal space (a.u.)
 write(*,'("bdot \(a. u.\)")')
 do ii=1, 3
  write(*,'(3(1x,f20.12))')crys%bdot(:,ii)
 end do
!    integer :: nat !< number of atoms
 write(*,'("nat",1x,i6)')crys%nat
!    integer, pointer :: atyp(:) !< atomic species, atyp(1:nat)
!    real(DP), pointer :: apos(:,:) !< atomic positions, apos(1:3,1:nat) (alat)
 do ii=1,crys%nat
   write(*,'(1x,i3,3(1x,f12.6))')crys%atyp(ii),crys%apos(:,ii)
 end do
 write(*,*)"*************************"

 POP_SUB(show_crys_object)

end subroutine show_crys_object


subroutine read_crystal_file(crys)
 implicit none
 integer,parameter::iunit=1
 type(crystal),intent(out)::crys
!
 integer::ii,iost
 character(80):: crystal_file,aa

 PUSH_SUB(read_crystal_file)

 crystal_file="crystal.dat"

 write(*,*)'Opening ',trim(crystal_file)
 call open_file(iunit, FILE=trim(crystal_file), iostat=iost, STATUS="old", FORM="FORMATTED")
 IF (iost .NE. 0) call print_error('abi2bgw','Error opening file '//trim(crystal_file),1)

 read(iunit,*)aa,crys%celvol
 if( trim(aa) .ne. "celvol") call print_error('read_crystal_file','celvol expected instead of'//trim(aa),1)
 read(iunit,*)aa,crys%alat
 if( trim(aa) .ne. "alat") call print_error('read_crystal_file','alat expected instead of'//trim(aa),1)
 read(iunit,*)aa
 if( trim(aa) .ne. "avec") call print_error('read_crystal_file','avec expected instead of'//trim(aa),1)
 do ii=1,3
   read(iunit,*)crys%avec(:,ii)
 end do
 read(iunit,*)aa,crys%nat
 if( trim(aa) .ne. "nat") call print_error('read_crystal_file','nat expected instead of'//trim(aa),1)
!Allocate crys object:
 SAFE_ALLOCATE(crys%atyp,(crys%nat))
 SAFE_ALLOCATE(crys%apos,(3,crys%nat))

 read(iunit,*)aa
 if( trim(aa) .ne. "atyp") call print_error('read_crystal_file','atyp expected instead of'//trim(aa),1)
 do ii=1,crys%nat
   read(iunit,*)crys%atyp(ii)
 end do
 read(iunit,*)aa
 if( trim(aa) .ne. "apos") call print_error('read_crystal_file','apos expected instead of'//trim(aa),1)
 do ii=1,crys%nat
   read(iunit,*)crys%apos(:,ii)
 end do
 
 call close_file(iunit)

 POP_SUB(read_crystal_file)
end subroutine read_crystal_file

subroutine show_syms_object(syms)
  implicit none
  type(symmetry),intent(in)::syms
  integer::isym
!
  PUSH_SUB(show_syms_object)
!
  write(*,*)"**********************************"
  write(*,'("Showing components of syms object:")')
!    integer :: ntran         !< number of operations in full group
  write(*,'("Number of sym. operations",1x,i2)') syms%ntran
!    integer :: mtrx(3,3,48)  !< symmetry matrix
  write(*,'("Symmetry matrix")')
  do isym=1,syms%ntran
    write(*,'(3(3x,3(1x,i4)))',advance='yes') syms%mtrx(1:3,1:3,isym)
  end do
!    real(DP) :: tnp(3,48)    !< fractional translations
  write(*,'("Symmetry translations")')
  do isym=1,syms%ntran
    write(*,'(3(1x,f12.6))',advance='yes') syms%tnp(1:3,isym)
  end do
!    integer :: indsub(48)    !< symmetry operations in subgroup of q
!    integer :: kgzero(3,48)  !< Umklapp vectors for subgroup symmetry operations
!    integer :: cell_symmetry !< 0 = cubic, 1 = hexagonal
!PENDING:
  write(*,'("Cell symmetry",1x,i2)')syms%cell_symmetry
  write(*,*)"**********************************"
  POP_SUB(show_syms_object)
end subroutine show_syms_object

subroutine fill_syms_object(syms,abi_hdr,abi_geom,cell_symmetry)

  use misc_m, only : invert_matrix
 
  implicit none
  real(DP), parameter    :: twopi = 6.2831853071795864769252867665590D0
  integer,intent(in)::cell_symmetry
  type(symmetry),intent(out)::syms
  type(abi_hdr_type),intent(in)::abi_hdr
  type(abi_geom_type),intent(in)::abi_geom
!
  integer::isym,irow,icol
  real(DP)::matrix_aux(3,3),matrix_aux2(3,3)
!
  PUSH_SUB(fill_syms_object)
!    integer :: ntran         !< number of operations in full group
  syms%ntran=abi_hdr%number_of_symmetry_operations
!    integer :: ntranq        !< number of operations in small group of q
!    real(DP) :: rq(3)        !< The q-point this ntranq belongs to
!    integer :: mtrx(3,3,48)  !< symmetry matrix
! Get inverse transpose: ABINIT rotation matrices are with respect to real 
!   space lattice vectors, whereas for BGW we need them with respect to 
!   reciprocal space:   inv(R)=trans(R`), R=rotation matrix in real space, 
!   R`= rotation matrix in reciprocal space
  do isym=1,abi_hdr%number_of_symmetry_operations !transpose matrices
    matrix_aux(:,:)=real(abi_geom%reduced_symmetry_matrices(:,:,isym),DP)
    matrix_aux2=TRANSPOSE(matrix_aux)
    call invert_matrix(matrix_aux2,matrix_aux)
    syms%mtrx(:,:,isym)=matrix_aux
  end do
!    real(DP) :: tnp(3,48)    !< fractional translations
  syms%tnp(1:3,1:abi_hdr%number_of_symmetry_operations)=&
&  abi_geom%reduced_symmetry_translations(1:3,1:abi_hdr%number_of_symmetry_operations)&
&  *twopi
!    integer :: indsub(48)    !< symmetry operations in subgroup of q
!    integer :: kgzero(3,48)  !< Umklapp vectors for subgroup symmetry operations
!    integer :: cell_symmetry !< 0 = cubic, 1 = hexagonal
!PENDING:
  syms%cell_symmetry=cell_symmetry
!
  POP_SUB(fill_syms_object)
end subroutine fill_syms_object

subroutine show_kp_object(kp,abi_hdr,abi_basis)
  implicit none
  type(kpoints),intent(in)::kp
  type(abi_hdr_type),intent(in)::abi_hdr
  type(abi_basis_type),intent(in)::abi_basis
  integer:: ik,ib,isppol,nbandk

  PUSH_SUB(show_kp_object)

  write(*,*)"**********************************"
  write(*,'("Showing components of kp object:")')
  write(*,'("Number of spinors",1x,i4)') kp%nspinor
  write(*,'("Number of spin polarizations",1x,i4)') kp%nspin
  write(*,'("Number of k-points",1x,i6)') kp%nrk
  write(*,'("Maximum number of bands",1x,i6)') kp%mnband
!!    integer :: nvband  !< number of valence bands
!!    integer :: ncband  !< number of conduction bands
  write(*,'("K-points grid",3(1x,i6))') kp%kgrid(:)
  write(*,'("K-points shift",3(1x,f12.6))') kp%shift(:)
  write(*,'("Wavefunction energy cutoff",1x,f12.6)')  kp%ecutwfc
  write(*,'("Number of planewaves for each k-point")')
  do ik=1, abi_hdr%number_of_kpoints
   write(*,'(1x,i6)',advance='no') kp%ngk(ik)
  end do
  write(*,*)
  write(*,'("Maximum number of planewaves",1x,i6)')  kp%ngkmax
  write(*,'("Lowest occupied band for each kpoint and spin")')
  do isppol=1,abi_hdr%number_of_spin_polarizations
    do ik=1, abi_hdr%number_of_kpoints
     write(*,'(1x,i6)',advance='no') kp%ifmin(ik,isppol)
    end do
    write(*,*)
  end do
  write(*,'("Highest occupied band for each kpoint and spin")')
  do isppol=1, abi_hdr%number_of_spin_polarizations
    do ik=1, abi_hdr%number_of_kpoints
     write(*,'(1x,i6)',advance='no') kp%ifmax(ik,isppol)
    end do
    write(*,*)
  end do
  write(*,'("Weights for each kpoint")')
  do ik=1,abi_hdr%number_of_kpoints
    write(*,'(1x,f12.6)',advance='no') kp%w(ik)
  end do
  write(*,*)
  write(*,'("K-vectors in crystal coordinates")')
  do ik=1,abi_hdr%number_of_kpoints
    write(*,'(i6,3(1x,f12.6))') ik, kp%rk(:,ik)
  end do
  write(*,'("Band energies for each band, kpoint and spin")')
  do isppol=1,abi_hdr%number_of_spin_polarizations
    write(*,'("spin=",1x,i6)')isppol
    do ik=1,abi_hdr%number_of_kpoints
      write(*,'(1x,i6)',advance='no')ik
      nbandk=abi_basis%number_of_states(ik,isppol)
      do ib=1,nbandk
        write(*,'(1x,f12.6)',advance='no') kp%el(ib,ik,isppol)
      end do
    write(*,*)
    end do
  end do
! elda is the same as above
!!    real(DP), pointer :: occ(:,:,:)  !< occupations (between 0 and 1)
  write(*,'("Occupatoins for each band, kpoint and spin")')
  do isppol=1,abi_hdr%number_of_spin_polarizations
    write(*,'("spin=",1x,i6)')isppol
    do ik=1,abi_hdr%number_of_kpoints
      write(*,'(1x,i6)',advance='no')ik
      nbandk=abi_basis%number_of_states(ik,isppol)
      do ib=1,nbandk
        write(*,'(1x,f6.4)',advance='no') kp%occ(ib,ik,isppol)
      end do
    write(*,*)
    end do
  end do
  write(*,*)"**********************************"

!     kp%occ => occupations
  POP_SUB(show_kp_object)
end subroutine show_kp_object


subroutine fill_kp_object(kp,abi_hdr,abi_basis,abi_elec,abi_kpt,wfng_nk,wfng_dk)
 implicit none
 type(kpoints),intent(inout)::kp
 type(abi_hdr_type), intent(in):: abi_hdr
 type(abi_basis_type), intent(in):: abi_basis
 type(abi_elec_type), intent(in):: abi_elec
 type(abi_kpt_type), intent(in):: abi_kpt
 integer,intent(in) :: wfng_nk(3)
 real(DP),intent(in) :: wfng_dk(3)
 real(DP), PARAMETER :: ha2ry=2.0d0

 PUSH_SUB(fill_kp_object)
!Fills BGW kp object from ABINIT variables:
  kp%nspinor = abi_hdr%number_of_spinor_components
  kp%nspin = abi_hdr%number_of_spin_polarizations
!    integer :: nrk     !< number of k-points
  kp%nrk=abi_hdr%number_of_kpoints
!    integer :: mnband  !< max number of bands
  kp%mnband=abi_hdr%max_number_of_states/abi_hdr%number_of_spin_polarizations
!    integer :: nvband  !< number of valence bands
!    integer :: ncband  !< number of conduction bands
!    integer  :: kgrid(3) !< Monkhorst-Pack number of k-points in each direction
  kp%kgrid=wfng_nk
!    real(DP) :: shift(3) !< Monkhorst-Pack shift of grid
  kp%shift=wfng_dk
!    real(DP) :: ecutwfc            !< wave-function cutoff, in Ry
  kp%ecutwfc=abi_hdr%energy_cutoff*ha2ry
!    integer, pointer :: ngk(:)     !< number of g-vectors for each k-point
  SAFE_ALLOCATE(kp%ngk, (kp%nrk))
  kp%ngk = abi_basis%number_of_coefficients
!    integer :: ngkmax              !< max(ngk(:))
  kp%ngkmax=abi_hdr%max_number_of_coefficients
!    integer, pointer :: ifmin(:,:) !< lowest occupied band (kpoint,spin)
!    integer, pointer :: ifmax(:,:) !< highest occupied band (kpoint,spin)
  call fill_ifminmax(kp,abi_hdr,abi_basis,abi_elec)
!    real(DP), pointer :: w(:)      !< weights (kpoint) (between 0 and 1)
  SAFE_ALLOCATE(kp%w, (kp%nrk))
  kp%w = abi_kpt%kpoint_weights
!    real(DP), pointer :: rk(:,:)   !< k-vector (3, kpoint) in crystal coords
  SAFE_ALLOCATE(kp%rk, (3, kp%nrk))
  kp%rk = abi_kpt%red_coord_kpt
!    real(DP), pointer :: el(:,:,:) !< band energies (band, kpoint, spin)
  SAFE_ALLOCATE( kp%el,(abi_hdr%max_number_of_states,abi_hdr%number_of_kpoints,abi_hdr%number_of_spin_polarizations))
  kp%el = abi_elec%eigenvalues*ha2ry
!    real(DP), pointer :: elda(:,:,:) !< band energies before eqp correction
!  kp%elda => eigenvalues
!    real(DP), pointer :: occ(:,:,:)  !< occupations (between 0 and 1)
  SAFE_ALLOCATE(kp%occ, (kp%mnband, kp%nrk, kp%nspin))
  kp%occ = abi_elec%occupations
!    integer, pointer :: degeneracy(:,:,:) !< size of deg. subspace for (band, kpoint, spin)
 POP_SUB(fill_kp_object)
end subroutine fill_kp_object

subroutine fill_ifminmax(kp,abi_hdr,abi_basis,abi_elec)
 implicit none
 type(abi_hdr_type),intent(in)::abi_hdr
 type(abi_basis_type),intent(in)::abi_basis
 type(abi_elec_type),intent(in)::abi_elec
 type(kpoints),intent(inout)::kp
 real(DP),parameter::tol4=0.0001_DP
!
 integer::ik,isppol,nbandk,ib
 integer::max_occ_band,min_occ_band
 real(DP)::occ_
 logical::found_ifmin,found_ifmax

!    integer, pointer :: ifmin(:,:) !< lowest occupied band (kpoint,spin)
!    integer, pointer :: ifmax(:,:) !< highest occupied band (kpoint,spin)

 PUSH_SUB(fill_ifminmax)

 SAFE_ALLOCATE(kp%ifmin,(abi_hdr%number_of_kpoints,abi_hdr%number_of_spin_polarizations))
 SAFE_ALLOCATE(kp%ifmax,(abi_hdr%number_of_kpoints,abi_hdr%number_of_spin_polarizations))


 do isppol=1,abi_hdr%number_of_spin_polarizations
   do ik=1,abi_hdr%number_of_kpoints
     found_ifmin=.false.
     found_ifmax=.false.
     nbandk=abi_basis%number_of_states(ik,isppol)
     min_occ_band=1
     max_occ_band=nbandk
     kp%ifmax(ik,isppol)=max_occ_band
     kp%ifmin(ik,isppol)=min_occ_band
     do ib=1,nbandk
       occ_=abi_elec%occupations(ib,ik,isppol)
       if ( .not. found_ifmin ) then
         if( occ_ > 0.999d0 ) then
           min_occ_band=ib
           found_ifmin=.true.
         end if
       end if
     end do
     do ib=nbandk,1,-1
       occ_=abi_elec%occupations(ib,ik,isppol)
       if ( .not. found_ifmax ) then
         if ( occ_ > 0.999d0 ) then
            max_occ_band=ib
            found_ifmax=.true. !last occupied band
         end if
       end if
     end do
     kp%ifmax(ik,isppol)=max_occ_band
     kp%ifmin(ik,isppol)=min_occ_band
   end do !kpoints
 end do !spin

 POP_SUB(fill_ifminmax)

end subroutine fill_ifminmax


end module abi2bgw_m
