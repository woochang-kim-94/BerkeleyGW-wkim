!=============================================================================
!
! Utilities:
!
! (1) createfullbz.x      Originally By MJ      Last Modified 10/02/2011 (MJ)
!
!     Creates fullbz maps from WFN file, Real/Complex flavor.
!
!==============================================================================

#include "f_defs.h"

program analyzebz

  use global_m
  use misc_m
  use wfn_rho_vxc_io_m
  use fullbz_m
  use irrbz_m
  use subgrp_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp
  type(gspace) :: gvec
  type(grid) :: gr
  character*3 :: sheader, sformat
  character*7 :: sflavor
  character*11 :: inform 
  character*256 :: infile, usage
  integer :: iflavor, if, ii, jj, itran, nq, nrq, iq
  integer :: nargs
  logical :: informat, do_irrbz, skip_checkbz
  integer, allocatable :: neq(:),indrq(:)
  real(DP), allocatable :: rq(:,:)

  usage = 'Usage: analyze.x ASC|BIN infile'

! Get file names from command-line arguments

  nargs = command_argument_count()

  if (nargs .ne. 2) then
    call die(usage)
  endif

  call get_command_argument(1, sformat)
  call get_command_argument(2, infile)

! Open units

  if (sformat .eq. 'ASC') then
    informat = .true.
    inform = 'formatted'
  elseif (sformat .eq. 'BIN') then
    informat = .false.
    inform = 'unformatted'
  else  
    call die(usage)
  endif

  call open_file(unit=7,file=TRUNC(infile),form=inform,status='old')

  write(6,'(/,3x,a,/)') 'Reading file ' // TRUNC(infile) // '  ' // TRUNC(inform) 

  sheader = 'GET'
  iflavor = -1
  call read_header_type(7, informat, sheader, iflavor, kp, gvec, syms, crys)

! Output info

  write(6,'(3x,"File header:",17x,a,/)') sheader // '-' // TRUNC(sflavor)
  write(6,'(3x,"Crystal volume:",f32.14)') crys%celvol
  write(6,'(3x,"Number of G-vectors:",i12)') gvec%ng
  write(6,'(3x,"Number of spins:",i16)') kp%nspin
  if (sheader .eq. 'WFN') then
    write(6,'(3x,"Number of bands:",i16)') kp%mnband
    write(6,'(3x,"Number of k-points:",i13)') kp%nrk
    write(6,'(3x,"Number of symmetry elements:",i4)') syms%ntran
    write(6,'(3x,"Grid reduced from:",i12,i4,i4)') kp%kgrid(:)
    write(6,'(3x,"            shift:",f14.3,f8.3,f8.3)') kp%shift(:)
  else
    write(6,'(3x,"k-point information is only in the WFN file")')
    call die('Use WFN file')
  endif
  write(6,*)
  call close_file(7)

  gr%nr = kp%nrk
  SAFE_ALLOCATE(gr%r, (3, gr%nr))
  gr%r = kp%rk
  call fullbz(crys,syms,gr,syms%ntran,skip_checkbz,wigner_seitz=.false.,paranoid=.false.)
  SAFE_DEALLOCATE_P(gr%r)

  call open_file(unit=8,file='fullbz.dat',form='formatted',status='replace')
  write(8,'(i6)') gr%nf
  do if = 1, gr%nf
    write(8,'(f13.9, f13.9, f13.9, i6, i6)') gr%f(1, if), gr%f(2,if),gr%f(3,if),gr%itran(if), gr%indr(if)
  enddo
  write(8,'(i6)') syms%ntran
  write(8,'(9i3)') (((syms%mtrx(ii,jj,itran),ii=1,3),jj=1,3),itran=1,syms%ntran)
  write(8,'(3f13.9)') ((syms%tnp(jj,itran),jj=1,3),itran=1,syms%ntran)
  write(8,'(i6)') kp%nrk
  do if = 1, kp%nrk
    write(8,'(3f13.9)') kp%rk(1:3, if)
  enddo
  call close_file(8)
  write(6,'(3x,"Grid unfolded:",i16)') gr%nf
  write(6,'(3x,"Grids difference:",i13)') gr%nf-(kp%kgrid(1)*kp%kgrid(2)*kp%kgrid(3))
  write(6,'(3x,"Done",/)')
  
  inquire(file='analyzebz.inp',exist=do_irrbz)
  if (do_irrbz) then
    write(6,'(3x,"Reading analyzebz.inp for calculating irrbz.dat at q pts")')
    call open_file(unit=9,file='analyzebz.inp',status='old')
    read(9,*) nq
    allocate(rq(3,nq))
    do iq = 1, nq
      read(9,*) rq(1:3,iq)
    enddo
    call close_file(9)

    call open_file(unit=8,file='irrbz_create.dat',form='formatted',status='replace')
    write(8,*) nq
    write(8,*) (rq(1:3,iq), iq = 1,nq)

    SAFE_ALLOCATE(indrq, (gr%nf))
    SAFE_ALLOCATE(neq, (gr%nf))

    do iq = 1, nq
      write(6,'(3x,"Grid reduced from:",i12,i4,i4)') kp%kgrid(:)
      write(6,'(3x,"Analyzing:",f16.4,f8.4,f8.4)') rq(1:3,iq)
      call subgrp(rq(1:3,iq),syms)
      call irrbz(syms,gr%nf,gr%f,nrq,neq,indrq)
      write(8,'(i6)') nrq
      write(8,*) (neq(ii), ii=1,nrq)
      write(8,*) (indrq(ii), ii=1,nrq)
    enddo

    SAFE_DEALLOCATE(indrq)
    SAFE_DEALLOCATE(neq)

    call close_file(8)
  endif

  call dealloc_grid(gr)

end program analyzebz
