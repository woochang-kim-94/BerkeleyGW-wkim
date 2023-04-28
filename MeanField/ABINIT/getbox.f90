!===============================================================================
!
! Program:
!
! GETBOX             Written by  Tonatiuh Rangel         Last Modified 07/30/2014    
!
! This finds a box size for a finite system that contains a given 
! percentage of the density (for GW purposes)
!
! Developers:
!   Tonatiuh Rangel, Berkeley CA, trangel@lbl.gov
!
!===============================================================================

#include "f_defs.h"

 program getbox

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
! type(symmetry) :: syms
! type(kpoints) :: kp
! type(gspace) :: gvec
! character(len=3) :: sheader
! character(len=11) :: outform
 character(len=80) :: infile !, outfile, usage
! integer :: iflavor
! logical :: outformat
!
!Variables for ABINIT
 type(abi_hdr_type):: abi_hdr
 type(abi_geom_type):: abi_geom
 type(abi_basis_type):: abi_basis
 type(abi_kpt_type):: abi_kpt
 type(abi_elec_type):: abi_elec

!Input file variables
 character ( len = 256 ) :: rho_file_abi
 real(DP):: den_percentage
 logical::slab


!Check and read command line arguments
 CALL checkCommandLineInputs(rho_file_abi,den_percentage,slab)

!Read input file
! call read_input_file()

 call read_den()
!
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
subroutine read_den()
  use read_abi_m, only: abi_read_denpot,abi_denpot_type
  use fftw_m
  implicit none
  type(abi_denpot_type)::abi_den
!
  logical,parameter::extract=.true.
  integer,parameter::den_format=1 !total, up components
  real(DP),parameter::zero=0.d0, one=1.d0, two=2.d0
  real(DP),parameter:: bohr2ang=0.52917721d0
  complex(DPC),parameter::czero=(0d0,0d0)
  integer::ispden,ig1,ig2,ig3
  integer:: ig1_1,ig1_2,ig1_3,ig2_1,ig2_2,ig2_3,ig3_1,ig3_2,ig3_3
  integer:: ig_left, ig_right, ig_half
  integer, dimension(3) :: nfft
  real(DP) :: scale,norm,norm1,norm2,norm3,norm_axis
  real(DP) :: raux1,raux2
  real(DP) :: den_goal,per_axis,norm_slide,slide_vol
  real(DP) :: slide_size1,slide_size2,slide_size3
  real(DP) :: rprimd_final(3,3)
  logical :: found

  PUSH_SUB(read_den)

  abi_hdr%number_of_density_components=1 !not coded otherwise
!  Now read ABINIT DEN:
  call abi_read_denpot(abi_hdr,abi_kpt,abi_basis,abi_geom,abi_elec,&
&  rho_file_abi,abi_den,extract,den_format)
! Warning: abi_* objects (except for abi_hdr) were not upgraded, these are not used anyway...
!   Also BGW objects were not upgraded

! Fill crystal object:
  call fill_crys_object(crys,abi_hdr)
  if ( debug ) then
    call show_crys_object(crys)
  end if

!  FFT density
! Check: Common/fftw_inc.f90 and Sigma/mtxel.f90
  call setup_FFT_sizes(abi_hdr%fft_grid,nfft,scale)

  write(*,'(" FFT of density")')
  write(*,'(" Sizes for FFT:",3(1x,i6))')nfft
  write(*,'(" Scale:",1x,f12.6)')scale

! Integrate total density:
  norm=zero
  do ispden=1,abi_hdr%number_of_density_components
    do ig3=1,nfft(3)
      do ig2=1,nfft(2)
        do ig1=1,nfft(1)
          norm=norm+abi_den%data_3d(ig1,ig2,ig3,ispden)
        end do
      end do
    end do
  end do
  norm=norm*scale*crys%celvol
  write(*,*)'Integral of rho(R)= ',norm

! Now try to find smaller cell with N% of density:

! The total density will be:
! (per_axis)^3 x den, we will remove a factor of per_axis at each axis
  if ( slab ) then
    per_axis=(den_percentage/100.0)
  else
    per_axis=(den_percentage/100.0)**(one/3.d0)
  end if
  norm_axis=per_axis*norm
!
! rprimd_final will contain the final box
  rprimd_final=abi_hdr%primitive_vectors
!
  if (.not. slab) then
!   1) a1 axis:
    found=.false.
    ig_half=nfft(1)/2
    slide_vol=zero; norm1=zero
    do ig1=1,ig_half
      if( found ) cycle
      norm_slide=zero
      ig_left=ig_half-(ig1-1)
      ig_right=ig_half+ig1
!  
      do ispden=1,abi_hdr%number_of_density_components
        do ig3=1,nfft(3)
          do ig2=1,nfft(2)
            norm_slide=norm_slide+abi_den%data_3d(ig_left,ig2,ig3,ispden)
            norm_slide=norm_slide+abi_den%data_3d(ig_right,ig2,ig3,ispden)
          end do
        end do
      end do
      norm1=norm1+norm_slide*scale*crys%celvol
      if ( norm1 >= norm_axis ) then
        slide_size1=real(ig_right-ig_left+1,DP)/real(nfft(1),DP)
        rprimd_final(:,1)=rprimd_final(:,1)*slide_size1
        found=.true.
      end if
    end do

!  2) a2 axis:
    found=.false.
    ig_half=nfft(2)/2
    slide_vol=zero; norm1=zero
    do ig2=1,ig_half
      if (found) cycle
      norm_slide=zero
      ig_left=ig_half-(ig2-1)
      ig_right=ig_half+ig2
!  
      do ispden=1,abi_hdr%number_of_density_components
        do ig3=1,nfft(3)
          do ig1=1,nfft(1)
            norm_slide=norm_slide+abi_den%data_3d(ig1,ig_left,ig3,ispden)
            norm_slide=norm_slide+abi_den%data_3d(ig1,ig_right,ig3,ispden)
          end do
        end do
      end do
      norm1=norm1+norm_slide*scale*crys%celvol
      if ( norm1 >= norm_axis ) then
        slide_size2=real(ig_right-ig_left+1,DP)/real(nfft(2),DP)
        rprimd_final(:,2)=rprimd_final(:,2)*slide_size2
        found=.true.
      end if
    end do
  else
    slide_size1=1.d0
    slide_size2=1.d0
  end if !not slab

!3) a3 axis:
  found=.false.
  ig_half=nfft(3)/2
  slide_vol=zero; norm1=zero
  do ig3=1,ig_half
    if (found) cycle
    norm_slide=zero
    ig_left=ig_half-(ig3-1)
    ig_right=ig_half+ig3
!
    do ispden=1,abi_hdr%number_of_density_components
      do ig2=1,nfft(2)
        do ig1=1,nfft(1)
          norm_slide=norm_slide+abi_den%data_3d(ig1,ig2,ig_left,ispden)
          norm_slide=norm_slide+abi_den%data_3d(ig1,ig2,ig_right,ispden)
        end do
      end do
    end do
    norm1=norm1+norm_slide*scale*crys%celvol
    if ( norm1 >= norm_axis ) then
      slide_size3=real(ig_right-ig_left+1,DP)/real(nfft(3),DP)
      rprimd_final(:,3)=rprimd_final(:,3)*slide_size3
      found=.true.
    end if
  end do

!Debug, check that the above is actually true:
!Integrate on the final box:
  ig1=real(nfft(1),DP)/two*slide_size1
  ig1_1=nfft(1)/2-(ig1-1)
  ig1_2=nfft(1)/2+ig1
!
  ig2=real(nfft(2),DP)/two*slide_size2
  ig2_1=nfft(2)/2-(ig2-1)
  ig2_2=nfft(2)/2+ig2
!
  ig3=real(nfft(3),DP)/two*slide_size3
  ig3_1=nfft(3)/2-(ig3-1)
  ig3_2=nfft(3)/2+ig3
!
  norm1=zero
  do ispden=1,abi_hdr%number_of_density_components
    do ig3=ig3_1,ig3_2
      do ig2=ig2_1,ig2_2
        do ig1=ig1_1,ig1_2
          norm1=norm1+abi_den%data_3d(ig1,ig2,ig3,ispden)
        end do
      end do
    end do
  end do
  norm1=norm1*scale*crys%celvol

  write(*,*)'Integral of rho(R) in reduced cell= ',norm1
  write(*,*)"New lattice vectors (bohr):"
  write(*,'(3(1x,f20.12))')rprimd_final(:,1)
  write(*,'(3(1x,f20.12))')rprimd_final(:,2)
  write(*,'(3(1x,f20.12))')rprimd_final(:,3)

  write(*,*)"New lattice vectors (Angstrom):"
  write(*,'(3(1x,f20.12))')rprimd_final(:,1)*bohr2ang
  write(*,'(3(1x,f20.12))')rprimd_final(:,2)*bohr2ang
  write(*,'(3(1x,f20.12))')rprimd_final(:,3)*bohr2ang

  write(*,*)"Box size for GW (Angstrom):"
  write(*,*)"(twice as big along non-periodic directions)"
  if ( slab ) then
    write(*,'(3(1x,f20.12))')rprimd_final(:,1)*bohr2ang
    write(*,'(3(1x,f20.12))')rprimd_final(:,2)*bohr2ang
    write(*,'(3(1x,f20.12))')rprimd_final(:,3)*bohr2ang*two
  else
    write(*,'(3(1x,f20.12))')rprimd_final(:,1)*bohr2ang*two
    write(*,'(3(1x,f20.12))')rprimd_final(:,2)*bohr2ang*two
    write(*,'(3(1x,f20.12))')rprimd_final(:,3)*bohr2ang*two
  end if
! DEALLOCATIONS
 SAFE_DEALLOCATE(abi_den%data_3d)

 POP_SUB(read_den)

end subroutine read_den
!
subroutine read_input_file()
 implicit none
 integer::iost
 character(80)::aa

 PUSH_SUB(read_input_file)

 ! assign defaults
 rho_file_abi='o_DEN'
 den_percentage=90

 write(*,*)'Opening ',trim(infile)
 call open_file(iunit, FILE=trim(infile), iostat=iost, STATUS="old", FORM="FORMATTED")
!
 read(iunit,*) aa, rho_file_abi
 read(iunit,*) aa, den_percentage
 call close_file(iunit)

 write(*,*)'Summary of input variables'
 write(*,*) "rho_file_abi ", trim(rho_file_abi)
 write(*,*) "den_percentage   ",den_percentage

 POP_SUB(read_input_file)

end subroutine read_input_file
!
subroutine calculate_volume(volume,aa)
  implicit none
  real(DP),intent(in)::aa(3,3)
  real(DP),intent(out)::volume
 
  PUSH_SUB(getbox.calculate_volume)


  volume=aa(1,1)*(aa(2,2)*aa(3,3)-aa(3,2)*aa(2,3))+&
&        aa(2,1)*(aa(3,2)*aa(1,3)-aa(1,2)*aa(3,3))+&
&        aa(3,1)*(aa(1,2)*aa(2,3)-aa(2,2)*aa(1,3))

  POP_SUB(getbox.calculate_volume)

end subroutine calculate_volume


!!!#############################
SUBROUTINE checkCommandLineInputs(rho_file_abi,percentage,slab)
!!!#############################

  IMPLICIT NONE
  real(DP),intent(out)::percentage
  logical,intent(out)::slab
  character(256),intent(out)::rho_file_abi
  character(80)::infile,arg2,arg3
  integer::int1
 
  PUSH_SUB(checkCommandLineInputs)

  call get_command_argument(1,rho_file_abi)
  call get_command_argument(2,arg2)
  call get_command_argument(3,arg3)
  write(*,*)'rho_file_abi',trim(rho_file_abi)
  write(*,*)'arg2: ',trim(arg2)
  write(*,*)'arg3: ',trim(arg3)
  if( trim(arg2).ne.'') then
    read(arg2,*)percentage
  end if
  slab=.FALSE.
  if( trim(arg3).ne.'') then
    read(arg3,*)int1
    if(int1==1)slab=.TRUE.
  end if
  if (slab) write(*,*)'Slab calculation'

  if ( trim(rho_file_abi).eq.'') then

     write(*,*) 'Usage'
     write(*,*) 'getbox.x o_DEN percentage [slab_keyword]'
     write(*,*) 'o_DEN: abinit density file'
     write(*,*) 'percentage range [0-100) '
     write(*,*) 'slab_keyword: set to 1 for slab'
     write(*,*) 'Important: along the non-periodic directions, atoms should be centered at the middle of the unit cell.'
     call die("Exiting")
  end if

  POP_SUB(checkCommandLineInputs)

END SUBROUTINE checkCommandLineInputs

end program getbox

