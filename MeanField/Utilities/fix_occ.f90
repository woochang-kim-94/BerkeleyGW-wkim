!=============================================================================
! Utilities:
! 
! fix_occ         Originally By FHJ      Last Modified 06/04/2014 (FHJ)
!
! Fix the occupations of a WFN file. This is useful if you have split a
! large nscf calculations into smaller blocks, and you want the occupations
! of the merged WFN to be consistent. It can also fix inconsistencies in the
! occupations and in the FE for metals.
!
! Additionally, the script can also manually add an extra charge to the WFN
! file (useful for internal testings).
!
! Note: this file is not compiled by default, you should type: make fix_occ.x
!
!=============================================================================

#include "f_defs.h"

program fix_occ

  use global_m
  use wfn_rho_vxc_io_m
  use sort_m

  implicit none

  character :: ierror_C
  integer :: iflavor

  type(kpoints) :: kp
  type(gspace) :: gvec
  type(symmetry) :: syms
  type(crystal) :: crys

  real(DP) :: n_el, charge
  integer :: n_occ, n_sort
  real(DP) :: n_el2
  integer, allocatable :: sort(:)
  real(DP), pointer :: en_ord(:), weight_orig(:,:,:), weight_ord(:)
  real(DP) :: E_F, E_M
  integer :: is,ik,ib
  integer :: ii, ii2

  character(len=80) :: outfile, infile
  character(len=30) :: sheader
  integer, parameter :: file_out = 11, file_in = 10

  ! Read fix_occ.inp
  call open_file(15, file='fix_occ.inp', form='formatted', status='old')
  read(15,*) infile
  read(15,*) outfile
  read(15,*) n_el    ! orig number of electrons (leave as -1 to auto detect)
  read(15,*) charge  ! additional number of electrons to put in/remove
  call close_file(15)

  write(6,*) 'Input  -> ' , TRUNC(infile)
  write(6,*) 'Output -> ' , TRUNC(outfile)

  ! Open all the WFN files to be read
  call open_file(file_in,  file=infile,  form='unformatted', status='old')
  call open_file(file_out, file=outfile, form='unformatted', status='unknown')

  ! Read original header
  sheader='WFN'
  iflavor=0
  call read_header_type(file_in, .false., sheader, iflavor, kp, gvec, syms, crys, .true.)

  ! calculate number of electrons, if the user requests
  if (n_el <= TOL_SMALL) then
    n_el = calc_nel(kp)
  endif

  write(6,*) 'Original number of electrons:', n_el
  n_el = n_el + charge
  write(6,*) 'Number of electrons with extra charge:', n_el
  n_occ = IDNINT(n_el*0.5d0 * kp%nspin * kp%nrk)
  write(6,*) 'Target # of occ. states:', n_occ

  n_sort = kp%mnband * kp%nrk * kp%nspin
  SAFE_ALLOCATE(sort, (n_sort))
  SAFE_ALLOCATE(en_ord, (n_sort))
  SAFE_ALLOCATE(weight_orig, (kp%mnband, kp%nrk, kp%nspin))
  SAFE_ALLOCATE(weight_ord, (n_sort))

  ! flattens the energy array and sort by energy
  en_ord = RESHAPE(kp%el, (/n_sort/))
  call sortrx(n_sort, en_ord, sort)
  en_ord = en_ord(sort)

  ! make sure sum(kw)==1
  if (abs(sum(kp%w)-1)>TOL_SMALL) then
    write(0,'(a,f0.12)') 'ERROR: sum of the k-weights = ', sum(kp%w)
    call die('k-weights don`t add up to 1.')
  endif

  do ik=1, kp%nrk
    weight_orig(:, ik, :) = kp%w(ik)
  enddo
  weight_ord(:) = RESHAPE(weight_orig, (/n_sort/))
  weight_ord(:) = weight_ord(sort)*(2.0d0/kp%nspin)

  !Determine Fermi energy by counting electrons
  n_el2 = 0.0d0
  !looping over states
  ii2 = 0
  do ii=1, n_sort
    if (weight_ord(ii) < (n_el-n_el2)) then
      n_el2 = n_el2 + weight_ord(ii)
    else
      ii2 = ii - 1
      exit
    endif
  enddo

  if (ii2==0) ii2 = n_sort

  E_F = en_ord(ii2)
  if (ii2==n_sort) then
    E_M = E_F
    write(0,*) 'WARNING: could not find Fermi Energy! Assuming all states are occupied.'
  else
    E_M = (E_F + en_ord(ii2+1))*0.5d0
  endif

  write(6,*) 'Fermi Energy  (eV):', E_F*ryd
  write(6,*) 'Middle Energy (eV):', E_M*ryd

  kp%ifmin(:,:)=1
  kp%ifmax(:,:)=0
  kp%occ(:,:,:)=0.0d0

  write(6,*) 'Resetting occupations' 
  !TODO: deal with case where Fermi level spans a degenerate subspace
  do is=1,kp%nspin
    do ik=1,kp%nrk
      do ib=1,kp%mnband
        if (kp%el(ib,ik,is)<E_M) then
          kp%ifmax(ik,is) = kp%ifmax(ik,is) + 1
          kp%occ(ib,ik,is) = 1.0d0
        endif
      enddo
    enddo
  enddo

  n_el = calc_nel(kp)

  write(6,*) 'Final number of electrons:', n_el

  write(6,*) 'Saving new header' 
  call write_header_type(file_out, .false., sheader, iflavor, kp, gvec, syms, crys, .true.)

  call close_file(file_in)
  call close_file(file_out)

  write(6,*) 'Copying rest of the file'
  ! this is an external C subroutine
  call copy_tail(ierror_C)

  write(6,*) 'Done' 

contains

  ! We can`t use calc_qtot here, b/c that subroutine assumes that there is a 
  ! partial occupation in some cases, which would not be true here.
  ! However, we can simply sum the occupations from the WFN (assuming it`s correct)
  real(DP) function calc_nel (kp_)
    type(kpoints), intent(in) :: kp_
    integer :: ik_

    calc_nel = 0
    do ik_ = 1, kp_%nrk
      calc_nel = calc_nel + kp_%w(ik_) * sum(kp_%occ(:,ik_,:))
    enddo
    ! note: occ is normalized to 1
    calc_nel = calc_nel * 2.0d0 / dble(kp_%nspin)
    return
  end function calc_nel

end program fix_occ

