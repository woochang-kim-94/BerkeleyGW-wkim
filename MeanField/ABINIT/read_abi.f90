!===============================================================================
!
! Module:
!
! read_abi_m             Written by  Tonatiuh Rangel         Last Modified 07/30/2014    
!
! Reads ABINIT objects: wavefunctions (WFK), density (DEN)
! and potential (VXC) files
!
!
! Developers:
!   Tonatiuh Rangel, Berkeley CA, trangel@lbl.gov
!   Felipe Jornada, Berkeley CA
!
!===============================================================================

#include "f_defs.h"

module read_abi_m

 use global_m
 use print_error_m, only : print_error
 implicit none

 private
 public :: abi_read_header
 public :: abi_allocate_elec
 public :: abi_read_kpt_elec_loop
 public :: abi_read_denpot
!
 private :: abi_read_header_first
 private :: abi_read_header_second
 private :: abi_show_basic_variables
 private :: abi_read_denpot_grid
 private :: abi_check_variables

!Object containing header informations (basically dimensions):
 type,public :: abi_hdr_type
!  variables in dims object:
   integer::number_of_atoms
   integer::number_of_kpoints,number_of_density_components,number_of_spinor_components
   integer::number_of_spin_polarizations,number_of_symmetry_operations,number_of_atom_species
   integer::max_number_of_coefficients,max_number_of_states,npsp
   integer:: headform,fft_grid(3)
   real(DP)::primitive_vectors(3,3)
   real(DP):: energy_cutoff
   real(DP):: smearing_width
   integer,allocatable:: znucl(:)
   real(DP),allocatable:: xred(:,:)
 end type abi_hdr_type

 type, public:: abi_kpt_type
!  Variables that will be used in the kpoints group.
   real(DP), pointer :: red_coord_kpt(:, :)
   real(DP), pointer :: kpoint_weights(:)
 end type abi_kpt_type

 type, public:: abi_basis_type
!   Variables that will be used in the basisdata group.
   integer, allocatable  :: number_of_states(:,:)
   integer, pointer  :: number_of_coefficients(:)
   integer, allocatable  :: red_coord_pw(:, :)
 end type abi_basis_type

 type, public:: abi_wf_type
   ! Variables that are declared in the main program in a real case
   real(DP), allocatable :: coef_pw_k(:, :)
 end type abi_wf_type

 type, public :: abi_geom_type
!   Variables for the geometry group
    integer, allocatable :: reduced_symmetry_matrices(:,:,:)
    real(DP), allocatable :: reduced_symmetry_translations(:,:)
    integer, allocatable :: atom_species(:)
    real(DP), allocatable :: reduced_atom_positions(:,:)
    real(DP), allocatable :: atomic_numbers(:)
    character(3), allocatable :: chemical_symbols(:)
 end type abi_geom_type

!type for density/potential
 type, public :: abi_denpot_type
   real(DP),allocatable::data_3d(:,:,:,:)
 end type abi_denpot_type

 type, public :: abi_elec_type
 ! Variables for the electrons group
  !character(len=etsf_charlen) :: smearing_scheme
  real(DP),allocatable:: eigenvalues(:,:,:)
  real(DP),pointer:: occupations(:,:,:)
  integer, allocatable:: number_of_states(:,:)
  real(DP) :: smearing_width, fermi_energy
 end type abi_elec_type


CONTAINS

subroutine abi_read_denpot(hdr,kpt,basis,geom,elec,denpot_file_abi,denpot,extract,denformat,vxc_flavor)
 implicit none
 integer, intent(in)::denformat
 logical, intent(in)::extract
 type(abi_hdr_type),intent(inout)::hdr
 type(abi_basis_type),intent(inout)::basis
 type(abi_geom_type),intent(inout)::geom
 type(abi_kpt_type),intent(inout)::kpt
 type(abi_elec_type),intent(inout)::elec
 type(abi_denpot_type),intent(out)::denpot
 character ( len = 256 ),intent(in) :: denpot_file_abi
 integer, intent(in), optional :: vxc_flavor ! ZL reads
!
 integer,parameter:: iunit=1,ounit=2,tmpunit=3
 logical::lcheck
 integer,parameter:: file_format= 1  !binary
 integer::nr1,nr2,nr3,nspden,iost
 integer :: cplex ! ZL adds
!
 PUSH_SUB(abi_read_denpot)
!
 if(present(vxc_flavor)) then
   cplex = vxc_flavor
 else
   cplex = 1  ! ZL: default is real
 endif
!
 lcheck=.false.
 if( .not. extract) lcheck=.true.
!
 nspden=hdr%number_of_density_components

 call open_file(iunit, FILE=denpot_file_abi, iostat=iost, STATUS="old", FORM="UNFORMATTED")

!Read header
 write(*,'("Reading header of density file")')
!Read header informations
 call abi_read_header(lcheck,extract,hdr,geom,elec,kpt,basis,iunit)

 nr1=hdr%fft_grid(1)
 nr2=hdr%fft_grid(2)
 nr3=hdr%fft_grid(3)

!Read the function on the 3D grid
 SAFE_ALLOCATE(denpot%data_3d,(cplex*nr1,nr2,nr3,nspden)) ! ZL adds cplex*
!
 write(*,'(" Reading density array",1x,i6," x ",i6," x ",i6," x ",i6)')nr1,nr2,nr3,nspden
 write(*,*) "Check vxc_flavor/cplex", cplex ! ZL adds
 call abi_read_denpot_grid(file_format,denpot%data_3d,hdr%fft_grid(1:3),nspden,iunit,denformat,cplex)
!
 call close_file(iunit)
!
 POP_SUB(abi_read_denpot)
end subroutine abi_read_denpot

!
subroutine abi_read_denpot_grid(file_format,grid,ngrid,nspden,iunit,den_format,cplex)
 implicit none
!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: den_format   ! 1 up, dn components, 2 total, up components 
 integer,intent(in) :: file_format ! 0 ACII, 1 binary
 integer,intent(in) :: ngrid(3)         ! grid dimensions
 integer,intent(in) :: nspden           ! Number of density components
 integer,intent(in) :: iunit            ! File unit 
 integer,intent(in) :: cplex  ! ZL adds
!arrays
 real(DP),intent(out) :: grid(cplex*ngrid(1),ngrid(2),ngrid(3),nspden)

!Local variables--------------------------------------------------------
!scalars
 character(500)::message
 integer :: ierr,i1,i2,i3,ispden
 real(DP),allocatable :: grid_up(:,:,:),grid_dn(:,:,:)

! *************************************************************************
 PUSH_SUB(abi_read_denpot_grid)

 if( file_format == 0 ) then
   if( nspden /= 1) then
     call die('read_abi: nspden /= 1 not allowed for formatted files')
   end if
   do i3=1,ngrid(3)
     do i2=1,ngrid(2)
       do i1=1,ngrid(1)*cplex
         read(unit=iunit,fmt=*) grid(i1,i2,i3,1)
       end do
     end do
   end do
 elseif ( file_format == 1)  then
   do ispden=1,nspden
     read(unit=iunit) grid(1:ngrid(1)*cplex,1:ngrid(2),1:ngrid(3),ispden)
   end do
!  Convert density to spin up, spin down
!  ABINIT format: if spin polarized, array contains total density in first half and spin-up density in second half)
   if ( nspden == 2 .and. den_format == 2 ) then
     SAFE_ALLOCATE( grid_up, (ngrid(1)*cplex,ngrid(2),ngrid(3)))
     SAFE_ALLOCATE( grid_dn, (ngrid(1)*cplex,ngrid(2),ngrid(3)))
     grid_up=grid(:,:,:,2)
     grid_dn=grid(:,:,:,1)-grid_up
     grid(:,:,:,1)=grid_up
     grid(:,:,:,2)=grid_dn
     SAFE_DEALLOCATE( grid_up )
     SAFE_DEALLOCATE( grid_dn )
   end if
 else
   write(message,'(a,i0)')"value for density/potential file format is invalid: ",file_format
   call print_error("abi_read_denpot_grid",message,1)
 end if

 POP_SUB(abi_read_denpot_grid)

end subroutine abi_read_denpot_grid
!


!Loop to read kpt info, and save only abi_elec_type informations
subroutine abi_read_kpt_elec_loop(hdr,basis,elec,iunit,pert_wfn)
 implicit none
 integer,intent(in)::iunit
 type(abi_elec_type),intent(inout)::elec
 type(abi_hdr_type),intent(in)::hdr
 type(abi_basis_type),intent(in)::basis
 logical, intent(in) :: pert_wfn

 integer:: ipw,isppol,nspinor,nbandk,ibandk,ik,npw
 real(DP):: occ_factor
 real(DP),allocatable::occk(:),cg(:,:),eigen(:),eigentmp(:)
 integer, allocatable  :: red_coord_pw_k(:, :)
!
 PUSH_SUB(abi_read_kpt_elec_loop)

 SAFE_ALLOCATE(cg,(2,hdr%number_of_spinor_components*hdr%max_number_of_coefficients))
 SAFE_ALLOCATE(red_coord_pw_k,(3, hdr%max_number_of_coefficients))
 SAFE_ALLOCATE(eigen,(hdr%max_number_of_states))
 SAFE_ALLOCATE(occk,(hdr%max_number_of_states))
!
!Need to rescale spins to [0-1] for BGW
 if(hdr%number_of_spin_polarizations==1 .and. hdr%number_of_spinor_components==1 ) then
   occ_factor=0.5d0
 else
   occ_factor=1.0d0
 end if
!
! big kpt loop to read G vectors, wave coefficients, 
! eigenvalues and occupations
!
 DO isppol=1, hdr%number_of_spin_polarizations
  ipw=0;
  DO ik=1, hdr%number_of_kpoints
   !write(*,*)'spin ',isppol,'kpt ',ik,'/',number_of_kpoints
   READ(iunit)npw,nspinor,nbandk      ! for each k point
   !check that the number of plane waves for this k-point
   !corresponds to the value already stored in number_of coefficients
   if(npw .ne. basis%number_of_coefficients(ik))&
       & call print_error('abi2bgw','npw read is different from the stored value [1]',1)
!
!Read G vectors (only one element, this is not used here:)
!   READ(iunit)red_coord_pw_k(1:3,1:npw)
   READ(iunit)red_coord_pw_k(1:3,1)
!
!Read eigenvalues and occupationsa
! ZL: for ground-state wfn
   if (.not. pert_wfn) then
     READ(iunit)eigen(1:nbandk),occk(1:nbandk)
   endif
!Fill the etsf variables containing the eigenvalues and the occupations   
   do ibandk=1,nbandk
    elec%eigenvalues( ibandk,ik,isppol)=eigen(ibandk)
    elec%occupations( ibandk,ik,isppol)=occk(ibandk)*occ_factor
   end do

!  reads the wavefunction coefficients for each band
   DO ibandk=1,nbandk
     if (pert_wfn) then
       READ(iunit)  ! ZL: for first-order pert wfn, no need to read eigen
     end if
!    READ(iunit)((cg(ii,ij),ii=1,2),&
!      &ij=1,npw*number_of_spinor_components)
!Probably we can just read 1 element here? we don`t use cg here:
    READ(iunit)cg(1:2,1)
   END DO!iband

  END DO !ik
 END DO !dims%number_of_spin_polarizations
!
 SAFE_DEALLOCATE(cg)
 SAFE_DEALLOCATE(red_coord_pw_k)
 SAFE_DEALLOCATE(eigen)
 SAFE_DEALLOCATE(occk)
!
 POP_SUB(abi_read_kpt_elec_loop)
end subroutine abi_read_kpt_elec_loop

subroutine abi_allocate_elec(elec,hdr)
  implicit none
  type(abi_elec_type),intent(inout)::elec
  type(abi_hdr_type),intent(inout)::hdr

  PUSH_SUB(abi_allocate_elec)

  SAFE_ALLOCATE( elec%eigenvalues,(hdr%max_number_of_states,hdr%number_of_kpoints,hdr%number_of_spin_polarizations))
  SAFE_ALLOCATE( elec%occupations,(hdr%max_number_of_states,hdr%number_of_kpoints,hdr%number_of_spin_polarizations))

  POP_SUB(abi_allocate_elec)
   
end subroutine abi_allocate_elec

subroutine abi_read_header(lcheck,extract,hdr,geom,elec,kpt,basis,iunit)
 implicit none
 integer, intent(in):: iunit
 logical, intent(in):: lcheck  !check consistency towards previously stored values
 logical, intent(in):: extract !whether to fill in bgw objects 
 type(abi_hdr_type),intent(inout)::hdr
 type(abi_geom_type),intent(inout)::geom
 type(abi_elec_type),intent(inout)::elec
 type(abi_basis_type),intent(inout)::basis
 type(abi_kpt_type),intent(inout)::kpt
!
 logical,parameter::debug=.true.

  PUSH_SUB(abi_read_header)

!Read header
!Start reading wavefunction
 call abi_read_header_first(iunit,hdr,lcheck)
!
!Read wavefunction basic quantities
!as nband,atoms information, symmetries, k-points and pseudopotentials.
 call abi_read_header_second(hdr,kpt,elec,basis,geom,extract,iunit)
!Check that variables are good:
 call abi_check_variables(hdr)
!
!  Show basic variables:
 if( debug ) call abi_show_basic_variables(hdr,geom,elec)
!
  POP_SUB(abi_read_header)
end subroutine abi_read_header


subroutine abi_show_basic_variables(hdr,geom,elec)
 implicit none
 type(abi_hdr_type),intent(in)::hdr
 type(abi_geom_type),intent(in)::geom
 type(abi_elec_type),intent(in)::elec
!
 character(len=1), parameter :: ch10 = char(10)
!
  PUSH_SUB(abi_show_basic_variables)
 write(*,'("hdr%fft_grid=",3i6)')&
  &hdr%fft_grid(1),&
  &hdr%fft_grid(2),&
  &hdr%fft_grid(3)
 write(*,'("number of atoms = ",i3)')hdr%number_of_atoms
 write(*,'("number of atom species = ",i3)')hdr%number_of_atom_species
!
 write(*,'("chemical_symbols ",10000a3)')geom%chemical_symbols(:)
!
 write(*,'("max number of coefficients (plane waves)=",&
  &i6,a,"max number of states (bands)=",i5,a,"number of kpoints=",&
  &i5,a,"number of spinor components=",i2,a,&
  &"number of spins=",i2)')  hdr%max_number_of_coefficients,ch10 &
  &,hdr%max_number_of_states,ch10,hdr%number_of_kpoints,ch10,&
  &hdr%number_of_spinor_components&
  &,ch10,hdr%number_of_spin_polarizations

 write(*,'("Fermi energy=",1x,f20.12)')elec%fermi_energy

  POP_SUB(abi_show_basic_variables)

end subroutine abi_show_basic_variables

subroutine abi_check_variables(hdr)
 implicit none
 type(abi_hdr_type), intent(in):: hdr

 if ( hdr%number_of_symmetry_operations > 48 ) then
   call print_error('abi_check_variables',&
& 'Number of symmetries >48, Please use symmorph=0 in your ABINIT input file',1)
 end if
 
end subroutine

subroutine abi_read_header_second(hdr,kpt,elec,basis,geom,extract,iunit)
 implicit none
 logical, intent(in) :: extract
 integer, intent(in) :: iunit
 type(abi_hdr_type), intent(inout):: hdr
 type(abi_kpt_type), intent(inout):: kpt
 type(abi_elec_type), intent(inout):: elec
 type(abi_basis_type), intent(inout):: basis
 type(abi_geom_type),intent(inout)::geom
!
 integer:: ii, ik, isppol, iatom,iband,ipsp,itypat,kptopt,lmn_size,icoulomb
 integer :: pawcpxocc,pspso,pspdat,pspcod,pspxc,lmax,lloc
 integer :: kptrlatt(3,3),kptrlatt_orig(3,3)
 character*132 :: title
 real(DP) :: znuclpsp,zionpsp,fermi_energy
 real(DP) :: charge,residm,etotal,nelect
 integer,allocatable::istwfk(:),nband(:),so_psp(:),rhoijselect(:,:)
 integer,allocatable:: symafm(:),nrhoijsel(:),npwarr(:),typat(:)
 integer,allocatable::syms(:,:,:)
 real(DP),allocatable ::amu(:),occ(:),tnons(:,:),occ3d(:,:,:)
 real(DP),allocatable::znucltypat(:),wk(:),kpts(:,:),xred(:,:)
 real(dp),allocatable :: shiftk_orig(:,:)
 real(dp),allocatable :: shiftk(:,:)
  ! shiftk(3,nshiftk), shiftks after inkpts   

  PUSH_SUB(abi_read_header_second)

 SAFE_ALLOCATE( amu, (hdr%number_of_atom_species))
 SAFE_ALLOCATE( istwfk,(hdr%number_of_kpoints))
 SAFE_ALLOCATE( so_psp,(hdr%npsp))
 SAFE_ALLOCATE( symafm,(hdr%number_of_symmetry_operations))
 SAFE_ALLOCATE(npwarr,(hdr%number_of_kpoints) )
 SAFE_ALLOCATE(nband,(hdr%number_of_kpoints*hdr%number_of_spin_polarizations))
 SAFE_ALLOCATE(syms,(3,3,hdr%number_of_symmetry_operations))
 SAFE_ALLOCATE(typat,(hdr%number_of_atoms) )
 SAFE_ALLOCATE(tnons,(3,hdr%number_of_symmetry_operations))
 SAFE_ALLOCATE(znucltypat,(hdr%number_of_atom_species))
 SAFE_ALLOCATE(xred,(3,hdr%number_of_atoms) )
 SAFE_ALLOCATE(kpts,(3, hdr%number_of_kpoints))
 SAFE_ALLOCATE(wk,(hdr%number_of_kpoints) )
 SAFE_ALLOCATE(occ,(hdr%max_number_of_states*hdr%number_of_kpoints))
 SAFE_ALLOCATE(occ3d,(hdr%max_number_of_states,hdr%number_of_kpoints,hdr%number_of_spin_polarizations))


  if ( hdr%headform >= 80 ) then
    read(iunit) &
    & istwfk(:), nband(:), npwarr(:), &
    & so_psp(:), symafm(:), syms(:,:,:), &
    & typat(:), kpts(:,:), occ3d, &
    & tnons(:,:), znucltypat(:), wk(:)
! Transfer occ3d to occ:
   ii = 0
   do ipsp=1,hdr%number_of_spin_polarizations
     do ik=1,hdr%number_of_kpoints
       do iband=1,nband(ik + (ipsp-1) * hdr%number_of_kpoints)
           ii = ii +1
           occ(ii) = occ3d(iband,ik,ipsp)
       end do
     end do
   end do

 else
   read(iunit)  istwfk(:),nband(:),npwarr(:),&
    &so_psp(:),symafm(:),syms(:,:,:),&
    &typat(:),kpts(:,:),occ(:),&
    &tnons(:,:),znucltypat(:),wk(:)
 end if
 ! FHJ&ZL: TODO - check that this is indeed the right thing to do
 wk(:) = wk(:) / sum(wk)

 if ( hdr%headform >= 80 ) then
!   Reading the final record of the header  ---------------------------------
   read(iunit) residm, xred(:,:), etotal, fermi_energy,amu
   
   read(iunit)&
      kptopt,pawcpxocc,nelect,charge,icoulomb,&
      kptrlatt,kptrlatt_orig !, hdr%shiftk_orig,hdr%shiftk

!   Reading the records with psp information ---------------------------------
   do ipsp=1,hdr%npsp
     read(iunit) &
&      title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc,lmn_size!,md5pseudos
   end do

! PAW is not supported, so we skip this
!   if (hdr%usepaw==1) then ! Reading the Rhoij tab if the PAW method was used.
!     call
!     pawrhoij_io(hdr%pawrhoij,unit,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,hdr%headform,"Read")
!   end if

 else

!  Read pseudopotentials info
   do ipsp=1,hdr%npsp
!   (npsp lines, 1 for each pseudopotential ; npsp=dims%number_of_atom_species, except if alchemical pseudo-atoms)
    read(iunit) title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc,lmn_size
   enddo
!  
!  (in case of usepaw==0, final record: residm, coordinates, total energy, Fermi energy)
   read(iunit) residm,xred(1:3,1:hdr%number_of_atoms),etotal,fermi_energy
 end if


 if( extract ) then
   SAFE_ALLOCATE(hdr%znucl,(hdr%number_of_atoms))
   SAFE_ALLOCATE(hdr%xred,(3,hdr%number_of_atoms))
   SAFE_ALLOCATE(basis%number_of_coefficients,(hdr%number_of_kpoints) )
   SAFE_ALLOCATE(basis%number_of_states,(hdr%number_of_kpoints,hdr%number_of_spin_polarizations) )
   SAFE_ALLOCATE(geom%reduced_symmetry_matrices,(3,3,hdr%number_of_symmetry_operations))
   SAFE_ALLOCATE(geom%atom_species,(hdr%number_of_atoms) )
   SAFE_ALLOCATE(geom%chemical_symbols,(hdr%number_of_atom_species) )
   SAFE_ALLOCATE(geom%reduced_symmetry_translations,(3,hdr%number_of_symmetry_operations))
   SAFE_ALLOCATE(geom%atomic_numbers,(hdr%number_of_atom_species))
   SAFE_ALLOCATE(geom%reduced_atom_positions,(3,hdr%number_of_atoms) )
   SAFE_ALLOCATE(kpt%red_coord_kpt,(3, hdr%number_of_kpoints))
   SAFE_ALLOCATE(kpt%kpoint_weights,(hdr%number_of_kpoints) )
!   SAFE_ALLOCATE(elec%occupations,(hdr%max_number_of_states,hdr%number_of_kpoints,hdr%number_of_spin_polarizations))
!
   basis%number_of_coefficients=npwarr
   geom%reduced_symmetry_matrices=syms
   geom%atom_species=typat
   geom%reduced_symmetry_translations=tnons
   geom%atomic_numbers=znucltypat
   geom%reduced_atom_positions=xred
   kpt%red_coord_kpt=kpts
   kpt%kpoint_weights=wk
   elec%fermi_energy=fermi_energy
!   elec%occupations=occ
   hdr%xred=xred
   do iatom=1,hdr%number_of_atoms
     itypat=typat(iatom)
     hdr%znucl(iatom)=znucltypat(itypat)
   end do
   !Get the chemical symbols from the atomic numbers
    do iatom=1,hdr%number_of_atom_species
      call get_symbol_by_z(nint(geom%atomic_numbers(iatom)), geom%chemical_symbols(iatom))
    end do
!   basis%number_of_states=reshape(nband,(/hdr%number_of_kpoints,hdr%number_of_spin_polarizations/))
   ii=0
   do isppol=1,hdr%number_of_spin_polarizations
    do ik=1,hdr%number_of_kpoints
      ii=ii+1
      basis%number_of_states(ik,isppol)=nband(ii)
    end do
   end do
   hdr%max_number_of_coefficients = maxval(basis%number_of_coefficients(:))
 end if
!
 SAFE_DEALLOCATE(amu)
 SAFE_DEALLOCATE(npwarr )
 SAFE_DEALLOCATE(nband)
 SAFE_DEALLOCATE(occ3d)
 SAFE_DEALLOCATE(syms)
 SAFE_DEALLOCATE(typat)
 SAFE_DEALLOCATE(tnons)
 SAFE_DEALLOCATE(znucltypat)
 SAFE_DEALLOCATE(xred)
 SAFE_DEALLOCATE(kpts)
 SAFE_DEALLOCATE(wk)
 SAFE_DEALLOCATE(occ)
  POP_SUB(abi_read_header_second)
end subroutine abi_read_header_second


!Check m_hdr in ABINIT
subroutine abi_read_header_first(iunit,hdr,lcheck)
 implicit none
 logical,intent(in) :: lcheck !check consistency towards previously stored values
 integer,intent(in) :: iunit
 type(abi_hdr_type),intent(inout):: hdr

 logical, parameter:: debug=.false.
 integer::ierr,iwarn
! Variables to read _WFK file
 character*6 :: codvsn
 integer:: ipsp,isppol,ii,ij,iband,iatom
 integer :: headform,fform
 integer :: bantot,date,intxc,ixc,mband,ngfft(3),npsp
 integer :: occopt,pertcase,usepaw
 integer :: natom,nkpt,nspden,nspinor,nsym,ntypat,nsppol
 integer :: npw,nshift_orig,nshiftk,nshiftk_orig,usewvl
 real(DP) :: acell(3),ecut,ecutdg,ecutsm,ecut_eff
 real(DP) :: kpttmp(3),qptn(3),stmbias,tphysel,tsmear
 real(DP) :: rprimd(3,3)

  PUSH_SUB(abi_read_header_first)

 READ(iunit) codvsn,headform,fform

 if ( headform >= 80 ) then
   READ(iunit) bantot, date, intxc, ixc, natom, ngfft,&
   &  nkpt, nspden, nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase,&
   &  usepaw, ecut, ecutdg, ecutsm, ecut_eff, qptn(1:3),rprimd,&
   &  stmbias, tphysel, tsmear, usewvl, nshiftk_orig, nshiftk, mband

 else
   READ(iunit) bantot,date,intxc,ixc,natom,ngfft,&
    &nkpt,nspden,nspinor,nsppol,&
    &nsym,npsp,ntypat,occopt,pertcase,usepaw,&
    &ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd&
    &,stmbias,tphysel,tsmear,usewvl
 end if

!Check that this agree with stored values
 if( lcheck ) then
   ierr=0; iwarn=0
   if(hdr%number_of_kpoints/=nkpt) iwarn=iwarn+1
   if(hdr%max_number_of_states/=bantot/nkpt) iwarn=iwarn+1
!
   if(hdr%number_of_spinor_components/=nspinor) ierr=ierr+1
   if(hdr%number_of_density_components/=nspden) ierr=ierr+1
   if(hdr%number_of_spin_polarizations/=nsppol) ierr=ierr+1
   if(hdr%number_of_symmetry_operations/=nsym) ierr=ierr+1
   if(hdr%number_of_atom_species/=ntypat) ierr=ierr+1
   if(sum(hdr%fft_grid(:))/=sum(ngfft(1:3))) ierr=ierr+1
   if(hdr%energy_cutoff/=ecut_eff) ierr=ierr+1
   if(hdr%npsp/=npsp) ierr=ierr+1
   if(ierr/=0) then
     write(*,*)'Inconsistency found in basic input parameters'
     write(*,*)'Please check your WFK/DEN file'
     write(*,*)'List of stored \& read values:'
     write(*,'(a,1x,i6,1x,"/=",1x,i6,1x)')'nspinor',hdr%number_of_spinor_components,nspinor
     write(*,'(a,1x,i6,1x,"/=",1x,i6,1x)')'nspden ',hdr%number_of_density_components,nspden
     write(*,'(a,1x,i6,1x,"/=",1x,i6,1x)')'nspin  ',hdr%number_of_spin_polarizations,nsppol
     write(*,'(a,1x,i6,1x,"/=",1x,i6,1x)')'nsym   ',hdr%number_of_symmetry_operations,nsym
     write(*,'(a,1x,i6,1x,"/=",1x,i6,1x)')'ntypat ',hdr%number_of_atom_species,ntypat
     write(*,'(a,3(1x,i6),1x,"/=",3(1x,i6),1x)')'ngfft  ',hdr%fft_grid(:),ngfft(1:3)
     write(*,'(a,1x,f12.6,1x,"/=",1x,f12.6,1x)')'ecut   ',hdr%energy_cutoff,ecut_eff
     write(*,'(a,1x,i6,1x,"/=",1x,i6,1x)')'npsp   ',hdr%npsp,npsp
     call die('Inconsistency found in basic input parameters')
   end if
   if(iwarn/=0) then
     write(*,*)'Inconsistency found in the following input parameters'
     write(*,*)'List of stored \& read values:'
     write(*,'(a,1x,i6,1x,"/=",1x,i6,1x)')'nkpt   ',hdr%number_of_kpoints,nkpt
     write(*,'(a,1x,i6,1x,"/=",1x,i6,1x)')'nband  ',hdr%max_number_of_states,bantot/nkpt
     write(*,*)'There may be a difference in kpoints/bands in the DEN/WFK calculations'
     write(*,*)'Assuming an experienced user, the calculation will continue'
   end if
 end if
!
!Transfer to hdr array
 hdr%number_of_atoms=natom
 hdr%number_of_kpoints=nkpt
 hdr%number_of_spinor_components=nspinor
 hdr%number_of_density_components=nspden
 hdr%number_of_spin_polarizations=nsppol
 hdr%number_of_symmetry_operations=nsym
 hdr%number_of_atom_species=ntypat
 hdr%fft_grid(:)=ngfft(1:3)
 hdr%energy_cutoff=ecut_eff
 hdr%primitive_vectors=rprimd
 hdr%smearing_width=tsmear
 hdr%npsp=npsp
 hdr%headform=headform
 if (headform >= 80 ) then
   hdr%max_number_of_states=mband
 else
   hdr%max_number_of_states=bantot/nkpt
 end if

 if(debug) then
   write(*,*)'**********************************'
   write(*,*)'Header components:'
   write(*,*)'Codvsn,headform,fform', codvsn,headform,fform
   write(*,*)'max_number_of_states ', bantot
   write(*,*)'date ',date
   write(*,*)'intxc,ixc ',intxc,ixc
   write(*,*)'number_of_atoms ',natom
   write(*,*)'ngfft ', ngfft(1:3)
   write(*,*)'number_of_kpoints', nkpt
   write(*,*)'number_of_components         ', nspden
   write(*,*)'number_of_spinor_components  ', nspinor
   write(*,*)'number_of_spin_polarizations ', nsppol
   write(*,*)'number_of_symmetry_operations', nsym
   write(*,*)'npsp                         ', npsp
   write(*,*)'number_of_atom_species       ', ntypat
   write(*,*)'occopt                       ', occopt
   write(*,*)'pertcase                     ', pertcase
   write(*,*)'usepaw                       ', usepaw
   write(*,*)'ecut                         ', ecut
   write(*,*)'ecutdg                       ', ecutdg
   write(*,*)'ecutsm                       ', ecutsm
   write(*,*)'ecut_eff                     ', ecut_eff
   write(*,*)'qptn(1:3)                    ', qptn(1:3)
   write(*,*)'primitive_vectors(1:3,1:3)   ', rprimd(1:3,1:3)
   write(*,*)'stmbias                      ', stmbias
   write(*,*)'tphysel                      ', tphysel
   write(*,*)'smearing_width               ', tsmear
   write(*,*)'**********************************'
 end if
  
  POP_SUB(abi_read_header_first)

end subroutine abi_read_header_first
!

end module read_abi_m
