!===============================================================================
!
! Program:
!
! ABI2BGW             Written by  Tonatiuh Rangel         Last Modified 07/30/2014    
!
! ABI2BGW reads ABINIT output files and 
! writes BGW input files.
!
! Beta version.
! ZL adds reading complex VXC for dVXC (or VXC1 in abinit language)
! ZL adds reading perturbed wfn for abinit.               10/04/2017
!
! Developers:
!   Tonatiuh Rangel, Berkeley CA, trangel@lbl.gov
!   Felipe Jornada, Berkeley CA
!   Zhenglu Li (ZL), Berkeley CA, lzlwell@berkeley.edu
!
!===============================================================================

#include "f_defs.h"

 program abi2bgw

 use global_m
 use wfn_rho_vxc_io_m
 use read_abi_m, only: abi_hdr_type,abi_geom_type,abi_basis_type,abi_kpt_type,abi_elec_type
 use abi2bgw_m

 implicit none
 character(len=1), parameter :: ch10 = char(10)
 integer,parameter:: iunit=1,ounit=2,tmpunit=3
 logical,parameter:: debug=.false.
 integer:: iost
!
!BGW objects:
 type(crystal) :: crys
 type(symmetry) :: syms
 type(kpoints) :: kp
 type(gspace) :: gvec
 character(len=3) :: sheader
 character(len=11) :: outform
 character(len=80) :: infile !, outfile, usage
 integer :: iflavor
 logical :: outformat
!
!Variables for ABINIT
 type(abi_hdr_type):: abi_hdr
 type(abi_geom_type):: abi_geom
 type(abi_basis_type):: abi_basis
 type(abi_kpt_type):: abi_kpt
 type(abi_elec_type):: abi_elec

!Input file variables
 character ( len = 256 ) :: wfng_file,rhog_file,vxcg_file
 character ( len = 256 ) :: wfng_file_abi,rhog_file_abi,vxcg_file_abi
 logical :: wfng_flag,rhog_flag,symrel_file_flag,vxcg_flag
 integer :: wfng_nk(3)
 real (DP) :: wfng_dk(3)
 integer :: cell_symmetry
 integer :: vxc_flavor ! ZL adds to deal with complex dVXC
 logical :: pert_wfn ! ZL adds to deal with first-order perturbative wavefunction
 real (DP) :: phonon_q(3) ! ZL adds
 integer :: ik

!Variables for density:
 integer, pointer :: red_coord_pw_rho(:, :)

!Check and read command line arguments
 CALL checkCommandLineInputs(infile)

!Read input file
 call read_input_file()

!Read ABINIT WFK file 
!(read only basic quantities, skip wfk coefficients)
 call read_abinit_basic_wfk(abi_hdr,abi_geom,abi_basis,abi_elec,abi_kpt)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Fill BGW objects:

!Read file with information on the crystal structure
 call fill_crys_object(crys,abi_hdr)

 if ( debug ) then
   call show_crys_object(crys)
 end if

!Fills kp object
 call fill_kp_object(kp,abi_hdr,abi_basis,abi_elec,abi_kpt,wfng_nk,wfng_dk)
 if ( debug ) then
   call show_kp_object(kp,abi_hdr,abi_basis)
 end if

!Fills syms object
 call fill_syms_object(syms,abi_hdr,abi_geom,cell_symmetry)
!Read symmetries from file
 if(symrel_file_flag) then
   call update_symmetries_from_file(syms)
   call show_syms_object(syms)
 end if
 if ( debug ) then
   call show_syms_object(syms)
 end if
! 
!Important to fill crys before this
!Creates sphere of gvectors in gvec object
 call fill_gvec_object(gvec,kp,crys,abi_hdr,red_coord_pw_rho)

 if ( debug ) then
   call show_gvec_object(gvec)
 end if

!Write BGW file:
! outformat = .true.
! outform = 'formatted'
!This is for binary:
 outformat = .false.
 outform = 'unformatted'
!
 sheader='WFN' !For wavefunctions
 iflavor=2 !1=real, 2=cplx ZL: in abinit, we always use cplx flavor

 ! ZL: shift back rk by phonon_q for BGW convention
 if (pert_wfn) then
   do ik=1, kp%nrk
     kp%rk(:,ik) = kp%rk(:,ik) - phonon_q(:)
   end do
 end if

 if(wfng_flag) then
!Write wfng file
   call open_file(ounit,file=trim(wfng_file),iostat=iost,form=outform,status='unknown')

   ! ZL: if doing pert_wfn, then shift kgrid by phonon_q,
   !     because abinit labels ( d_q psi_k ) by k+q for ngk issue,
   !     however, in BerkeleyGW, we want it to be on k (not k+q) because of the
   !     phase issue
   !     so we need to shift it back: (k_abi_pert - phonon_q) is really what we want 
   call write_header_type(ounit, outformat, sheader, iflavor, kp, gvec, syms, crys)
   call write_gvectors(ounit, outformat, gvec%ng, gvec%ng, gvec%components)
!  Now read again ABINIT WFK:
   if ( abi_hdr%number_of_spin_polarizations == 1 ) then
!    read and write the cg coefficients into the BGW file:
!    this saves memory, but only works for nspin=1
     call read_write_cg()
   else
!    Needs to save all WFN in memory
!    With the HDF5 format this can be improved
     call read_write_cg_spin()
   end if
!
   call close_file(ounit)
 end if
!
 if(rhog_flag) then
   Write(*,*)'Writing DEN file:'
   ! ZL: avoid potential inconsistency
   if ((vxc_flavor.eq.2) .or. pert_wfn) then
     call die("RHO cannot be calculated with vxc_flavor=2 or pert_wfn=.true.")
   end if
!Write rhog file
   call open_file(ounit,file=trim(rhog_file),iostat=iost,form=outform,status='unknown')
!
   sheader='RHO'
   call write_header_type(ounit, outformat, sheader, iflavor, kp, gvec, syms, crys)
   call write_gvectors(ounit, outformat, gvec%ng, gvec%ng, gvec%components)
   call read_write_den()
!
   call close_file(ounit)
 end if
!
 if(vxcg_flag) then
   Write(*,*)'Writing VXC file:'
!Write vxcg file
   call open_file(ounit,file=trim(vxcg_file),iostat=iost,form=outform,status='unknown')
!
   sheader='VXC'
   call write_header_type(ounit, outformat, sheader, iflavor, kp, gvec, syms, crys)
   call write_gvectors(ounit, outformat, gvec%ng, gvec%ng, gvec%components)
   call read_write_vxc()
!
   call close_file(ounit)
 end if

 write(*,*)'Calculation completed'


!
!Deallocate
! call deallocate_abi()
! call deallocate_crys()
!!  call dealloc_header_type(sheader, crys, kp)
!!  PENDING
!!  deallocate gvec% object
!
contains
!
subroutine read_write_den()
  use read_abi_m, only: abi_read_denpot,abi_denpot_type
  use fftw_m
  implicit none
  type(abi_denpot_type)::abi_den
!
  logical,parameter::extract=.false.
  integer, parameter :: denformat=2 !total and up components
  real(DP),parameter::zero=0.d0, one=1.d0
  complex(DPC),parameter::czero=(0d0,0d0)
  integer::ispden,ig1,ig2,ig3
  complex(DPC), dimension(:,:,:), allocatable :: fftbox
  complex(DPC), allocatable :: charge_density(:,:)
  integer, dimension(3) :: nfft
  real(DP) :: scale,norm,raux1,raux2

  PUSH_SUB(read_write_den)

!  Now read ABINIT DEN:
  call abi_read_denpot(abi_hdr,abi_kpt,abi_basis,abi_geom,&
&  abi_elec,rhog_file_abi,abi_den,extract,denformat)
! Warning: abi_* objects (except for abi_hdr) were not upgraded, these are not used anyway...
!   Also BGW objects were not upgraded

!  FFT density
! Check: Common/fftw_inc.f90 and Sigma/mtxel.f90
  call setup_FFT_sizes(gvec%FFTgrid,nfft,scale)

!  write(*,'(" FFT of density")')
!  write(*,'(" Sizes for FFT:",3(1x,i6))')nfft
!  write(*,'(" Scale:",1x,f12.6)')scale
!    
! Allocate FFT boxes
!  
  SAFE_ALLOCATE(fftbox, (Nfft(1),Nfft(2),Nfft(3)))
  SAFE_ALLOCATE(gvec%index_vec, (gvec%ng))
  SAFE_ALLOCATE(charge_density, (gvec%ng,abi_hdr%number_of_density_components))
  do ig1=1,gvec%ng
    gvec%index_vec(ig1)=ig1
  end do

  norm=zero
  fftbox=czero
  do ispden=1,abi_hdr%number_of_density_components
!   1) cp density into fftbox
    do ig3=1,nfft(3)
      do ig2=1,nfft(2)
        do ig1=1,nfft(1)
          fftbox(ig1,ig2,ig3)=abi_den%data_3d(ig1,ig2,ig3,ispden)
          norm=norm+abi_den%data_3d(ig1,ig2,ig3,ispden)
        end do
      end do
    end do
    raux1=crys%celvol; raux2=one/product(gvec%FFTgrid(:))
    norm=norm*raux1*raux2
    fftbox=fftbox*raux1
!    write(*,'( "Norm of rho(R,spin ",i1," = ",f10.2)')ispden,norm
!2) perform the FFT on fftbox: call do_FFT, using sign=-1
    call do_FFT(fftbox,nfft,-1)
   
!3) get all the FFT of the CD up to ng: 
    call get_from_fftbox(gvec%ng,charge_density(:,ispden),gvec%components,gvec%index_vec,fftbox,nfft,scale)
    do ig1=1,gvec%ng
      if( gvec%components(1,ig1)==0 .and. gvec%components(2,ig1)==0 .and. gvec%components(3,ig1)==0 ) ig2=ig1
    end do
!4) Check norm: rho(G=0)=Ne
    raux1=real(charge_density(ig2,ispden))
    raux2=aimag(charge_density(ig2,ispden))
    norm= sqrt(raux1**2+raux2**2)
!    write(*,'(" Norm of rho(G=0)=",1x,f12.6)')norm
    write(*,'( " Norm of rho( spin = ",i1," ) = ",f10.2)')ispden,norm

  end do

! Write density in BGW file
  call write_complex_data(ounit, outformat, gvec%ng, gvec%ng, &
&  abi_hdr%number_of_density_components, charge_density)
!
! DEALLOCATIONS
 SAFE_DEALLOCATE(charge_density)
 SAFE_DEALLOCATE(fftbox)
 SAFE_DEALLOCATE(abi_den%data_3d)
 SAFE_DEALLOCATE_P(gvec%index_vec)

 POP_SUB(read_write_den)

end subroutine read_write_den
!
subroutine read_write_vxc()
  use read_abi_m, only: abi_read_denpot,abi_denpot_type
  use fftw_m
!  use misc_m, only: set_jspinor
  implicit none
  type(abi_denpot_type)::abi_vxc
!
  logical, parameter :: extract=.false.
  integer, parameter :: denformat=1 !up and dn components
  real(DP),parameter :: zero=0.d0, one=1.d0
  real(DP),parameter :: pi=3.14159265359d0
  complex(DPC),parameter::czero=(0d0,0d0)
  integer::ispden,ig1,ig2,ig3,nvxc
  complex(DPC), dimension(:,:,:), allocatable :: fftbox
  complex(DPC), allocatable :: vxc(:,:)
  integer, dimension(3) :: nfft
  real(DP) :: scale,norm,raux1,raux2
!
  PUSH_SUB(read_write_vxc)

!  Now read ABINIT DEN:
  call abi_read_denpot(abi_hdr,abi_kpt,abi_basis,abi_geom,abi_elec,&
&  vxcg_file_abi,abi_vxc,extract,denformat,vxc_flavor=vxc_flavor)

!if ( abi_hdr%number_of_spinor_components == 2 ) then
!  abi_vxc%data_3d(:,:,:,1)=abi_vxc%data_3d(:,:,:,1)+abi_vxc%data_3d(:,:,:,2)
!  abi_vxc%data_3d=abi_vxc%data_3d/(2.d0,0.d0)
!  abi_vxc%data_3d(:,:,:,2:)=(0.d0,0.d0)
!end if

! Warning: abi_* objects (except for abi_hdr) were not upgraded, these are not used anyway...
!   Also BGW objects were not upgraded

!  FFT density
! Check: Common/fftw_inc.f90 and Sigma/mtxel.f90
  call setup_FFT_sizes(gvec%FFTgrid,nfft,scale)

!  write(*,'(" FFT of vxc")')
!  write(*,'(" Sizes for FFT:",3(1x,i6))')nfft
!  write(*,'(" Scale:",1x,f12.6)')scale
!    
! Allocate FFT boxes
! 
  nvxc=abi_hdr%number_of_density_components*abi_hdr%number_of_spinor_components
  SAFE_ALLOCATE(fftbox, (Nfft(1),Nfft(2),Nfft(3)))
  SAFE_ALLOCATE(gvec%index_vec, (gvec%ng))
  SAFE_ALLOCATE(vxc,(gvec%ng,nvxc))
  vxc=(0.d0,0.d0) !initialization
  do ig1=1,gvec%ng
    gvec%index_vec(ig1)=ig1
  end do

  norm=zero
  fftbox=czero
  do ispden=1,abi_hdr%number_of_spinor_components
!   1) cp density into fftbox
    do ig3=1,nfft(3)
      do ig2=1,nfft(2)
        do ig1=1,nfft(1)
          if (vxc_flavor .eq. 1) then
            fftbox(ig1,ig2,ig3)=abi_vxc%data_3d(ig1,ig2,ig3,ispden)
            norm=norm+abi_vxc%data_3d(ig1,ig2,ig3,ispden)
          else if (vxc_flavor .eq. 2) then
            fftbox(ig1,ig2,ig3)=CMPLX(abi_vxc%data_3d(ig1*2-1,ig2,ig3,ispden),abi_vxc%data_3d(ig1*2,ig2,ig3,ispden))
            norm=norm+abi_vxc%data_3d(ig1*2-1,ig2,ig3,ispden)  ! ZL: sum only real part
!            write(*,*) "real",real(fftbox(ig1,ig2,ig3))
!            write(*,*) "imag",aimag(fftbox(ig1,ig2,ig3))
!            write(*,*) "fft", real(fftbox(ig1,ig2,ig3)), aimag(fftbox(ig1,ig2,ig3))
          endif
        end do
      end do
    end do
   !raux1=crys%celvol/(4.0*pi); 
   raux1=2.d0; !ha2ry
   raux2=one/product(gvec%FFTgrid(:))
   norm=norm*raux1*raux2
   fftbox=fftbox*raux1
   write(*,*)'Norm of vxc(R)= ',norm, 'nvxc=',nvxc, 'raux1,raux2=',raux1,raux2
!2) perform the FFT on fftbox: call do_FFT, using sign=-1
    call do_FFT(fftbox,nfft,-1)
!    do ig3=1,nfft(3)
!      do ig2=1,nfft(2)
!        do ig1=1,nfft(1)
!          write(*,*)'ig3,ig2,ig1',ig3,ig2,ig1,real(fftbox(ig1,ig2,ig3)),aimag(fftbox(ig1,ig2,ig3))
!        end do
!      end do
!    end do
!3) get all the FFT of the CD up to ng: 
    call get_from_fftbox(gvec%ng,vxc(:,ispden),gvec%components,gvec%index_vec,fftbox,nfft,scale)
    do ig1=1,gvec%ng
      if( gvec%components(1,ig1)==0 .and. gvec%components(2,ig1)==0 .and. gvec%components(3,ig1)==0 ) ig2=ig1
    end do
!4) Check norm: vxc(G=0)
    raux1=real(vxc(ig2,ispden))
    raux2=aimag(vxc(ig2,ispden))
    norm= sqrt(raux1**2+raux2**2)
!    write(*,'(" Norm of vxc(G=0)=",1x,f12.6)')norm
    if (vxc_flavor .eq. 1) then
      write(*,'( " Norm of vxc( spin = ",i1," ) = ",f10.2)')ispden,norm
    else if (vxc_flavor .eq. 2) then
      write(*,'( " Norm of real part of vxc( spin = ",i1," ) = ",f10.2)')ispden,norm
!      write(*,*) "raux1,raux2",raux1,raux2
    endif
  end do

! Write density in BGW file
  call write_complex_data(ounit, outformat, gvec%ng, gvec%ng, nvxc, vxc)
!
! DEALLOCATIONS
 SAFE_DEALLOCATE(vxc)
 SAFE_DEALLOCATE(fftbox)
 SAFE_DEALLOCATE(abi_vxc%data_3d)
 SAFE_DEALLOCATE_P(gvec%index_vec)

 POP_SUB(read_write_vxc)
end subroutine read_write_vxc
!
subroutine read_write_cg()
 use read_abi_m, only: abi_read_header
 implicit none
 logical,parameter::extract=.false.
 logical,parameter::lcheck=.true.
 logical,parameter::save_cg=.false.
 integer::iost

 PUSH_SUB(read_write_cg)

 call open_file(iunit, FILE=wfng_file_abi, iostat=iost, STATUS="old", FORM="UNFORMATTED")

!Read header
 call abi_read_header(lcheck,extract,abi_hdr,abi_geom,abi_elec,abi_kpt,abi_basis,iunit)
!END DEBUG
!
!Read & write G vectors and wavefunction coefficients
!in a big kpt_loop
 call read_write_kpt_loop(abi_hdr,abi_basis,iunit,ounit,outformat,save_cg,pert_wfn)
!
!Close input file
 call close_file(iunit)

 POP_SUB(read_write_cg)

end subroutine read_write_cg

!Option to read and write cg, when spin > 1
!this needs to save all WFN in memory !!
subroutine read_write_cg_spin()
 use read_abi_m, only: abi_read_header
 implicit none
 logical,parameter::extract=.false.
 logical,parameter::lcheck=.true.
 logical,parameter::save_cg=.true.  !save all cg into memory and then write
 integer::iost

 PUSH_SUB(read_write_cg_spin)

 call open_file(iunit, FILE=wfng_file_abi, iostat=iost, STATUS="old", FORM="UNFORMATTED")

!Read header
 call abi_read_header(lcheck,extract,abi_hdr,abi_geom,abi_elec,abi_kpt,abi_basis,iunit)
!END DEBUG
!
!Read & write G vectors and wavefunction coefficients
!in a big kpt_loop
 call read_write_kpt_loop(abi_hdr,abi_basis,iunit,ounit,outformat,save_cg,pert_wfn)
!
!Close input file
 call close_file(iunit)

 POP_SUB(read_write_cg_spin)

end subroutine read_write_cg_spin

!
!
!
!
!
!
!
!
subroutine read_input_file()
 implicit none
 integer::iost
 character(80)::aa

 PUSH_SUB(read_input_file)

 ! assign defaults
wfng_file_abi='o_WFK'
  wfng_flag = .FALSE.
  wfng_file = 'WFN'
 wfng_nk(:) = 0
 wfng_dk(:) = 0.0D0
  rhog_flag = .FALSE.
  rhog_file = 'RHO'
 rhog_file_abi='o_DEN'
 symrel_file_flag=.FALSE.
  vxcg_flag = .FALSE.
  vxcg_file = 'VXC'
 rhog_file_abi='o_VXC'
 vxc_flavor = 1  ! ZL adds: 1 for real and 2 for complex. Abinit deals with cplx by two real(dp)
 pert_wfn = .FALSE.  ! ZL adds: if this is first-order perturbative wavefunction. Default is .false.  
 phonon_q(:) = 0.0D0  ! ZL adds: phonon_q coordinates, default is (0.0, 0.0, 0.0)
!  vxc_flag = .FALSE.
!  vxc_file = 'vxc.dat'

 write(*,*)'Opening ',trim(infile)
 call open_file(iunit, FILE=trim(infile), iostat=iost, STATUS="old", FORM="FORMATTED")
!
 read(iunit,*) aa, wfng_file_abi
 read(iunit,*) aa, wfng_flag
 read(iunit,*) aa, wfng_file
 read(iunit,*) aa, wfng_nk(1)
 read(iunit,*) aa, wfng_nk(2)
 read(iunit,*) aa, wfng_nk(3)
 read(iunit,*) aa, wfng_dk(1)
 read(iunit,*) aa, wfng_dk(2)
 read(iunit,*) aa, wfng_dk(3)
 read(iunit,*) aa, rhog_file_abi
 read(iunit,*) aa, rhog_flag
 read(iunit,*) aa, rhog_file
 read(iunit,*) aa, cell_symmetry
 read(iunit,*) aa, symrel_file_flag
 read(iunit,*) aa, vxcg_file_abi
 read(iunit,*) aa, vxcg_flag
 read(iunit,*) aa, vxcg_file
 ! This part is optional
 read(iunit,*, iostat=iost) aa, vxc_flavor  ! ZL adds
 read(iunit,*, iostat=iost) aa, pert_wfn  ! ZL adds
 read(iunit,*, iostat=iost) aa, phonon_q(1)
 read(iunit,*, iostat=iost) aa, phonon_q(2)
 read(iunit,*, iostat=iost) aa, phonon_q(3)
 call close_file(iunit)

 write(*,*)'Summary of input variables'
 write(*,*) "wfng_file_abi ", trim(wfng_file_abi)
 if(wfng_flag) then
   write(*,*) "wfng_flag     ", wfng_flag
   write(*,*) "wfng_file     ", trim(wfng_file)
   write(*,*) "wfng_nk1      ", wfng_nk(1)
   write(*,*) "wfng_nk2      ", wfng_nk(2)
   write(*,*) "wfng_nk3      ", wfng_nk(3)
   write(*,*) "wfng_dk1      ", wfng_dk(1)
   write(*,*) "wfng_dk2      ", wfng_dk(2)
   write(*,*) "wfng_dk3      ", wfng_dk(3)
 end if
 if(rhog_flag) then
   write(*,*) "rhog_file_abi ", trim(rhog_file_abi)
   write(*,*) "rhog_flag     ", rhog_flag
   write(*,*) "rhog_file     ", trim(rhog_file)
 end if
 write(*,*) "cell_symmetry   ",cell_symmetry
 write(*,*) "symrel_file_flag", symrel_file_flag
 if(vxcg_flag) then
   write(*,*) "vxcg_flag     ", vxcg_flag
   write(*,*) "vxcg_file     ", trim(vxcg_file)
 end if
 write(*,*) "For electron-phonon linear response"
 write(*,*) "      vxc_flavor      ", vxc_flavor
 write(*,*) "      pert_wfn        ", pert_wfn
 write(*,*) "      phonon_q1        ", phonon_q(1)
 write(*,*) "      phonon_q2        ", phonon_q(2)
 write(*,*) "      phonon_q3        ", phonon_q(3)
 if ((vxc_flavor .ne. 1) .and. (vxc_flavor .ne.2)) then
   call die('Invalid value for vxc_flavor')
 endif

 POP_SUB(read_input_file)

end subroutine read_input_file
!
! Reads basic geom data, occupations and eigenvalues.
! Fills the geom, basis, elec and kpt objects
subroutine read_abinit_basic_wfk(hdr,geom,basis,elec,kpt)

 use read_abi_m, only: abi_read_header, abi_allocate_elec,&
& abi_read_kpt_elec_loop
 implicit none
 type( abi_hdr_type),intent(out):: hdr
 type( abi_geom_type),intent(out):: geom
 type( abi_basis_type),intent(out):: basis
 type( abi_elec_type),intent(out):: elec
 type( abi_kpt_type),intent(out):: kpt
 logical,parameter::extract=.true.
 logical,parameter::lcheck=.false.
 integer::iost

 PUSH_SUB(read_abinit_basic_wfk)

 call open_file(iunit, FILE=wfng_file_abi, iostat=iost, STATUS="old", FORM="UNFORMATTED")

!Read header informations
 call abi_read_header(lcheck,extract,hdr,geom,elec,kpt,basis,iunit)

!Allocage elec object
 call abi_allocate_elec(abi_elec,abi_hdr)
!
!Read eigenvalues and occupations
 call abi_read_kpt_elec_loop(abi_hdr,abi_basis,abi_elec,iunit,pert_wfn)

!Close input file
 call close_file(iunit)

 POP_SUB(read_abinit_basic_wfk)

end subroutine read_abinit_basic_wfk
!!
!subroutine deallocate_basic_abi()
!
! SAFE_DEALLOCATE( istwfk)
! SAFE_DEALLOCATE ( nband)
! SAFE_DEALLOCATE(number_of_coefficients )
! SAFE_DEALLOCATE(number_of_states)
! SAFE_DEALLOCATE( so_psp )
! SAFE_DEALLOCATE( symafm)
! SAFE_DEALLOCATE( reduced_symmetry_matrices)
! SAFE_DEALLOCATE( atom_species)
! SAFE_DEALLOCATE( chemical_symbols)
! SAFE_DEALLOCATE( red_coord_kpt)
! SAFE_DEALLOCATE( occ)
! SAFE_DEALLOCATE( reduced_symmetry_translations)
! SAFE_DEALLOCATE( atomic_numbers)
! SAFE_DEALLOCATE( kpoint_weights)
! SAFE_DEALLOCATE( reduced_atom_positions )
!
!end subroutine deallocate_basic_abi
!!
!
!
!subroutine deallocate_abi()
! SAFE_DEALLOCATE( istwfk )
! SAFE_DEALLOCATE(nband)
! SAFE_DEALLOCATE( so_psp)
! SAFE_DEALLOCATE( symafm)
! SAFE_DEALLOCATE( occ)
! SAFE_DEALLOCATE(number_of_coefficients)
! SAFE_DEALLOCATE(red_coord_kpt)
! SAFE_DEALLOCATE(kpoint_weights)
! SAFE_DEALLOCATE(number_of_states)
! SAFE_DEALLOCATE(reduced_symmetry_matrices)
! SAFE_DEALLOCATE(reduced_symmetry_translations)
! SAFE_DEALLOCATE(atomic_numbers)
! SAFE_DEALLOCATE(atom_species)
! SAFE_DEALLOCATE(reduced_atom_positions)
! SAFE_DEALLOCATE(chemical_symbols)
! SAFE_DEALLOCATE(eigenvalues)
! SAFE_DEALLOCATE(occupations)
! SAFE_DEALLOCATE(red_coord_pw)
!end subroutine deallocate_abi
!
!subroutine deallocate_crys
! SAFE_DEALLOCATE(crys_atyp)
! SAFE_DEALLOCATE(crys_apos)
!end subroutine deallocate_crys

!!!#############################
SUBROUTINE checkCommandLineInputs(infile)
!!!#############################

  IMPLICIT NONE
  character(80)::infile
 
  PUSH_SUB(checkCommandLineInputs)

  if (command_argument_count()<1) then
    infile = "abi2bgw.inp"
  else
    call get_command_argument(1, infile)
    if (len_trim(infile)==0) then
      infile = "abi2bgw.inp"
    endif
  endif

  POP_SUB(checkCommandLineInputs)

END SUBROUTINE checkCommandLineInputs

end program abi2bgw

