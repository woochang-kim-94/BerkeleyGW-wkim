!===============================================================================
!
! Program:
!
! EQP_CONVERT      Written by  Tonatiuh Rangel         Last Modified 06/2016    
!
! Converts BGW eqp1.dat to ABINIT o_EQP file for post-processing purposes.
!
! Developers:
!   Tonatiuh Rangel, Berkeley CA, trangel@lbl.gov
!
!===============================================================================

#include "f_defs.h"

 program eqp_convert

 use global_m
 use wfn_rho_vxc_io_m
 use read_abi_m, only: abi_hdr_type,abi_geom_type,abi_basis_type,abi_kpt_type,abi_elec_type, abi_denpot_type
 use abi2bgw_m
 use read_eqp_m
 use print_error_m, only : print_error

 implicit none

 character(len=1), parameter :: ch10 = char(10)
 integer,parameter:: iunit=1,ounit=2,tmpunit=3
 logical,parameter:: debug=.false.
 integer:: iost
!
!BGW objects:
 type(crystal) :: crys
 character(len=80) :: infile !, outfile, usage
!
!Variables for ABINIT
 type(abi_hdr_type):: abi_hdr
 type(abi_kpt_type):: abi_kpt
 type(abi_elec_type):: abi_elec

!Input file variables
 character ( len = 256 ) :: rho_file_abi, wfng_file_abi,eqp1_file

!Objects containing the DEN and QPS data:
 type(qps_type)::qps
 type(abi_denpot_type)::abi_den

!Check and read command line arguments
 CALL checkCommandLineInputs(rho_file_abi,wfng_file_abi,eqp1_file)

!Reads ABINIT DEN file:
 write(*,*)'Read ABINIT density file:'
 call read_den(abi_den)

!Reads ABINIT WFK file to extract eigenvalues and header info:
 write(*,*)'Read basic ABINIT data for WFK file:'
 call read_abinit_basic_wfk(abi_hdr,abi_elec,abi_kpt,wfng_file_abi)


!Read BGW eqp1 file:
 call read_eqp1(qps)

!Check consistency in k-points:
 call check_kpoints()

!Write QPS file:
 call write_qps()

!Deallocate
 SAFE_DEALLOCATE(abi_den%data_3d)
! call deallocate_abi()
! call deallocate_crys()
!!  call dealloc_header_type(sheader, crys, kp)
!!  PENDING
!!  deallocate gvec% object
!
contains
!
subroutine check_kpoints()
  implicit none
  real(dp), parameter::tol=1.d-4
  integer:: ik,ii
  real(dp):: raux
  character(80)::message


  PUSH_SUB(check_kpoints)
  do ik=1,abi_hdr%number_of_kpoints
    do ii=1,3
      raux=abs(abi_kpt%red_coord_kpt(ii,ik)-qps%k_points(ii,ik))
      if( raux > tol ) then
        write(message,'(a,3f8.4,a,3f8.4)') &
&        "Inconsistency in k-points found in WFK and eqp1.dat files",&
&        abi_kpt%red_coord_kpt(:,ik)," /=", qps%k_points(:,ik)
        call print_error("check_kpoints: ",message,1)
      end if
    end do
  end do

  POP_SUB(check_kpoints)
end subroutine check_kpoints
!
subroutine read_abinit_basic_wfk(hdr,elec,kpt,wfng_file_abi)

 use read_abi_m, only: abi_read_header, abi_allocate_elec,&
& abi_read_kpt_elec_loop
 implicit none
 type( abi_hdr_type),intent(out):: hdr
 type( abi_elec_type),intent(out):: elec
 type( abi_kpt_type),intent(out):: kpt
 character ( len = 256 ),intent(in) :: wfng_file_abi
 logical,parameter::extract=.true.
 logical,parameter::lcheck=.false.
 integer,parameter::iunit=3
 type( abi_geom_type):: abi_geom
 type( abi_basis_type):: abi_basis
 integer::iost

 PUSH_SUB(read_abinit_basic_wfk)
 call open_file(iunit, FILE=wfng_file_abi, iostat=iost, STATUS="old", FORM="UNFORMATTED")

!Read header informations
 call abi_read_header(lcheck,extract,abi_hdr,abi_geom,abi_elec,abi_kpt,abi_basis,iunit)

 call abi_allocate_elec(abi_elec,abi_hdr)
!
!Read eigenvalues and occupations
 call abi_read_kpt_elec_loop(abi_hdr,abi_basis,abi_elec,iunit,.false.)

!Close input file
 call close_file(iunit)

 POP_SUB(read_abinit_basic_wfk)
end subroutine read_abinit_basic_wfk

!!!#############################
!
subroutine write_qps()
  integer:: iost,isppol,nfft
  integer:: ib,jb,ik,max_band_gw
  integer:: map_sigma(abi_hdr%max_number_of_states,abi_hdr%number_of_kpoints)
  complex(dp),allocatable :: mtmp(:,:)
  real(dp),allocatable :: rhor_tmp(:)
  real(dp)::edft,eqp,sigma
  PUSH_SUB(write_qps)

! First we need to extrapolate bands which we don`t compute GW:
  do ik=1, abi_hdr%number_of_kpoints
!   by default we will apply the corrections to the first band we calculated
    map_sigma(:,ik)=qps%index_eigenvalues(1,ik)
    max_band_gw=qps%index_eigenvalues(qps%number_of_bands,ik)
    map_sigma(max_band_gw:,ik)=qps%number_of_bands

    do ib=1, qps%number_of_bands
      jb=qps%index_eigenvalues(ib,ik)
      map_sigma(jb,ik)=ib
    end do
    
  end do


  call open_file(3, FILE='o_QPS', iostat=iost, STATUS="new", FORM="FORMATTED")
  write(3,'("1")')
  write(3,*)qps%number_of_k_points
  write(3,*)abi_hdr%max_number_of_states !qps%number_of_bands
!  SAFE_ALLOCATE(mtmp,(qps%number_of_bands,qps%number_of_bands))
  SAFE_ALLOCATE(mtmp,(abi_hdr%max_number_of_states,abi_hdr%max_number_of_states))
  mtmp=CMPLX(0.0d0,0.0d0)
  do ib=1,abi_hdr%max_number_of_states !qps%number_of_bands
    mtmp(ib,ib)=CMPLX(1.0d0,0.0d0)
  end do
  write(3,*)qps%number_of_spins
  do isppol=1,qps%number_of_spins
     do ik=1,qps%number_of_k_points
       write(3,*)qps%k_points(:,ik)
       do ib=1,abi_hdr%max_number_of_states !qps%number_of_bands
         jb=map_sigma(ib,ik)
         eqp=qps%eigenvalues_eqp1(jb,ik,isppol)
         edft=qps%eigenvalues_dft(jb,ik,isppol)
         sigma=eqp-edft
         edft=abi_elec%eigenvalues(ib,ik,isppol)
         eqp=edft+sigma
         write(3,*)eqp
         write(3,*)mtmp(:,ib)
       end do
     end do
  end do
  SAFE_DEALLOCATE(mtmp) 
  write(3,*)abi_hdr%fft_grid(:)

  
  nfft=PRODUCT(abi_hdr%fft_grid(:))
  SAFE_ALLOCATE(rhor_tmp,(nfft))
  rhor_tmp=reshape(abi_den%data_3d,(/nfft/))
  write(3,*)rhor_tmp
  SAFE_DEALLOCATE(rhor_tmp)
!
  call close_file(3)
  POP_SUB(write_qps)

end subroutine write_qps
!
subroutine read_eqp1(qps)

 IMPLICIT NONE
!
 real(dp),parameter::ev2ha=0.0367493088676916
 type(qps_type),intent(out)::qps
 integer::ii,iost,nband_,error
 integer::ib,jb,ik,nkpt
 real(dp)::edft,eqp,kpt_tmp(3)
 character(120)::s

  PUSH_SUB(read_eqp1)

  nkpt=abi_hdr%number_of_kpoints
  qps%number_of_spins= abi_hdr%number_of_spinor_components*abi_hdr%number_of_spin_polarizations
  if ( qps%number_of_spins /= 1)  call print_error("eqp_convert","not coded for nspin > 1",1)
  ii=0
!
  call open_file(3, FILE=eqp1_file, iostat=iost, STATUS="old", FORM="FORMATTED")
!Get number of rows:
  do while(.true.)
    ii=ii+1
    read (3,'(a)',iostat=error) s
    if ( error.ne.0) exit
  end do
  call close_file(3)
  !write(*,*)'nrows=',ii

  nband_=ii/nkpt-1
  write(*,*)'Number of k-points in wfn/eqp file:',nkpt
  write(*,*)'Number of bands in eqp file:',nband_
  qps%number_of_k_points=nkpt
  qps%number_of_bands=nband_
  SAFE_ALLOCATE(qps%k_points,(3,nkpt))
  SAFE_ALLOCATE(qps%eigenvalues_dft,(nband_,nkpt,1))
  SAFE_ALLOCATE(qps%eigenvalues_eqp1,(nband_,nkpt,1))
  SAFE_ALLOCATE(qps%index_eigenvalues,(nband_,nkpt))
!
!Read the file again, and save data:
  call open_file(3, FILE=eqp1_file, iostat=iost, STATUS="old", FORM="FORMATTED")
  do ik=1,nkpt
    read(3,*)kpt_tmp(:),nband_
    !write(*,*)kpt_tmp,nband_
    qps%k_points(:,ik)=kpt_tmp
    do ib=1,nband_
      read (3,*,iostat=error) ii,jb,edft,eqp
      !write(*,*)ii,jb,edft,eqp
      if ( error.ne.0) then
        call die('Error found when reading eqp file')
      end if
      qps%eigenvalues_eqp1(ib,ik,1)=eqp*ev2ha
      qps%eigenvalues_dft(ib,ik,1)=edft*ev2ha
      qps%index_eigenvalues(ib,ik)=jb
    end do
  end do
  call close_file(3)
  POP_SUB(read_eqp1)

end subroutine read_eqp1

!
subroutine read_den(abi_den)
  use read_abi_m, only: abi_read_denpot,abi_denpot_type
  use fftw_m
  implicit none
!
  type(abi_denpot_type),intent(out)::abi_den
  type(abi_basis_type):: abi_basis
  type(abi_geom_type):: abi_geom
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

  PUSH_SUB(read_den)

  abi_hdr%number_of_density_components=1 !not coded otherwise
!  Now read ABINIT DEN:
  call abi_read_denpot(abi_hdr,abi_kpt,abi_basis,abi_geom,abi_elec,&
&  rho_file_abi,abi_den,extract,den_format)

! Fill crystal object:
  call fill_crys_object(crys,abi_hdr)
  if ( debug ) then
    call show_crys_object(crys)
  end if

! DEALLOCATIONS
! SAFE_DEALLOCATE(abi_den%data_3d)

 POP_SUB(read_den)

end subroutine read_den
!
subroutine calculate_volume(volume,aa)
  implicit none
  real(DP),intent(in)::aa(3,3)
  real(DP),intent(out)::volume
 
  PUSH_SUB(eqp_convert.calculate_volume)


  volume=aa(1,1)*(aa(2,2)*aa(3,3)-aa(3,2)*aa(2,3))+&
&        aa(2,1)*(aa(3,2)*aa(1,3)-aa(1,2)*aa(3,3))+&
&        aa(3,1)*(aa(1,2)*aa(2,3)-aa(2,2)*aa(1,3))

  POP_SUB(eqp_convert.calculate_volume)

end subroutine calculate_volume


!!!#############################
SUBROUTINE checkCommandLineInputs(rho_file_abi,wfng_file_abi,eqp1_file)
!!!#############################

  IMPLICIT NONE
  character(256),intent(out)::rho_file_abi,eqp1_file, wfng_file_abi

  PUSH_SUB(checkCommandLineInputs)

  call get_command_argument(1,rho_file_abi)
  call get_command_argument(2,wfng_file_abi)
  call get_command_argument(3,eqp1_file)

  if ( trim(eqp1_file).eq.'') then

     write(*,*) 'Usage'
     write(*,*) 'eqp_convert.x o_DEN o_WFK eqp1.dat '
     write(*,*) 'o_DEN: abinit density file'
     write(*,*) 'o_WFK: abinit wavefunction file'
     write(*,*) 'eqp1.dat: BGW eqp1.dat file '
     call die("Exiting")
  end if

  write(*,*)'rho file:',trim(rho_file_abi)
  write(*,*)'WFK file:',trim(wfng_file_abi)
  write(*,*)'eqp1.dat ',trim(eqp1_file)

  POP_SUB(checkCommandLineInputs)

END SUBROUTINE checkCommandLineInputs

end program eqp_convert

