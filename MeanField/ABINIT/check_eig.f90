!===============================================================================
!
! Program:
!
! check_eig            Written by  Tonatiuh Rangel         Last Modified 07/30/2014    
!
! This program is used to debug ABI2BGW, it will compare the
! eigenvalues read in an ABINIT and a BGW files.
!
!===============================================================================

#include "f_defs.h"
 program check_eig

 use global_m
 use wfn_rho_vxc_io_m
 use read_abi_m, only: abi_hdr_type,abi_geom_type,abi_basis_type,abi_kpt_type,abi_elec_type
 use abi2bgw_m

 implicit none
 character(len=1), parameter :: ch10 = char(10)
 integer,parameter:: iunit=1,ounit=2,tmpunit=3
!
!BGW objects:
 type(crystal) :: crys,crys_bgw
 type(symmetry) :: syms,syms_bgw
 type(kpoints) :: kp,kp_bgw
 type(gspace) :: gvec,gvec_bgw
 character(len=3) :: sheader
 character(len=11) :: outform
 character(len=80) :: infile !, outfile, usage
 integer :: iflavor
 logical :: informat
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
!DEBUG
! call show_crys_object(crys)
!END DEBUG
 call fill_kp_object(kp,abi_hdr,abi_basis,abi_elec,abi_kpt,wfng_nk,wfng_dk)
!DEBUG
! call show_kp_object(kp,abi_hdr,abi_basis)
!END DEBUG
!
 call fill_syms_object(syms,abi_hdr,abi_geom,cell_symmetry)
!Read symmetries from file
 if(symrel_file_flag) then
   call update_symmetries_from_file(syms)
   call show_syms_object(syms)
 end if
!DEBUG
! call show_syms_object(syms)
!END DEBUG
! 
!Important to fill crys before this
 call fill_gvec_object(gvec,kp,crys,abi_hdr,red_coord_pw_rho)
!DEBUG
! call show_gvec_object(gvec)
!END DEBUG




!Write BGW file:
! outformat = .true.
! outform = 'formatted'
!This is for binary:
! outformat = .false.
 outform = 'unformatted'
!
! sheader='WFN' !For wavefunctions
! iflavor=2 !1=real, 2=cplx

!Read BGW file:
 call open_file(tmpunit,file=TRIM(wfng_file),form=outform,status='old')
  sheader = 'GET'
  iflavor = -1
  call read_header_type(tmpunit, informat, sheader, iflavor, kp_bgw, gvec_bgw, syms_bgw, crys_bgw, dont_warn_kgrid = .true.)

  close(tmpunit)

  call compare_eigenvalues()

!
contains
!
subroutine compare_eigenvalues()
 implicit none
 integer::is,ik,ib
 real(dp)::e1,e2,error,avg_error,max_error

 if(kp%nspin .ne. kp_bgw%nspin) then
   write(*,'("Error found, nspin is not the same:")')
   write(*,'("nspin (ABI)=",i2,"nspin (BGW)=",i2)')kp%nspin,kp_bgw%nspin
   call die('nspin mismatch')
 end if
 if(kp%nrk .ne. kp_bgw%nrk) then
   write(*,'("Error found, nkpt is not the same:")')
   write(*,'("nkpt (ABI)=",i2,"nkpt (BGW)=",i2)')kp%nrk,kp_bgw%nrk
   call die('nkpt mismatch')
 end if
 if(kp%mnband .ne. kp_bgw%mnband) then
   write(*,'("Error found, nband is not the same:")')
   write(*,'("nband (ABI)=",i2,"nband (BGW)=",i2)')kp%mnband,kp_bgw%mnband
   call die('nband mismatch')
 end if

!Now compare eigenvalues:
  max_error=0.d0;avg_error=0.d0
  do is=1,kp%nspin
    do ik=1,kp%nrk
      do ib=1,kp%mnband
        e1=kp%el(ib,ik,is)
        e2=kp_bgw%el(ib,ik,is)
        error=abs(e2-e1)
        avg_error=avg_error+error
        if(error> max_error) max_error=error
      end do
    end do
  end do
  avg_error=avg_error/(real(kp%nspin)*real(kp%nrk)*real(kp%mnband))
  
  write(*,'("Max. error=",f20.12)')max_error
  write(*,'("Avg. error=",f20.12)')avg_error

end subroutine compare_eigenvalues
!
subroutine read_input_file()
 implicit none
 integer::iost
 character(80)::aa

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
 close(iunit)

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

end subroutine read_input_file
!
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

 call open_file(iunit, FILE=wfng_file_abi, iostat=iost, STATUS="old", FORM="UNFORMATTED")


!Read header informations
 call abi_read_header(lcheck,extract,hdr,geom,elec,kpt,basis,iunit)

 call abi_allocate_elec(abi_elec,abi_hdr)
!
!Read eigenvalues and occupations
 call abi_read_kpt_elec_loop(abi_hdr,abi_basis,abi_elec,iunit,.false.)

!Close input file
 call close_file(iunit)

end subroutine read_abinit_basic_wfk

!!!#############################
SUBROUTINE checkCommandLineInputs(infile)
!!!#############################

  implicit none
  character(80) :: infile

  call get_command_argument(1,infile)

  if ( trim(infile).eq.'') then
     write(*,*) 'Usage:'
     write(*,*) '  abi2bgw input_file'
     write(*,*) 'where:'
     write(*,*) '  input_file: see abi2bgw.inp example'
     stop
  end if
END SUBROUTINE checkCommandLineInputs

end program check_eig

