!=============================================================================
!
! Utilities:
!
! (1) scissors2eqp    Originally By DAS      Last Modified 3/25/2012 (DAS)
!
!     Write an eqp.dat file based on a WFN file and scissors parameters.
!     For testing equivalence of eqp_corrections and scissors.
!
!==============================================================================

#include "f_defs.h"

program scissors2eqp

  use global_m
  use scissors_m
  use wfn_rho_vxc_io_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp
  type(gspace) :: gvec
  character*3 :: sheader
  character*256 :: infile, usage, tmpstr
  integer :: iflavor, nargs, ik, is, ib, ii
  type(scissors_t) :: scis
  type(spline_tck) :: spl
  logical :: has_splines

  usage = 'Usage: scissors2eqp.x WFN evs ev0 evdel ecs ec0 ecdel [spline_scissors.dat]'

! Get file name from command-line arguments

  nargs = command_argument_count()

  if (nargs<7 .or. nargs>8) then
    call die(usage)
  endif
  has_splines = nargs==8

  call get_command_argument(1, infile)
  call open_file(unit=7,file=TRUNC(infile),form='unformatted',status='old')

  sheader = 'WFN'
  iflavor = -1
  call read_binary_header_type(7, sheader, iflavor, kp, gvec, syms, crys, warn = .false., &
    dont_warn_kgrid = .true.)
  call close_file(7)

  call get_command_argument(2, tmpstr)
  read(tmpstr,*) scis%val%es
  call get_command_argument(3, tmpstr)
  read(tmpstr,*) scis%val%e0
  call get_command_argument(4, tmpstr)
  read(tmpstr,*) scis%val%edel
  call get_command_argument(5, tmpstr)
  read(tmpstr,*) scis%cond%es
  call get_command_argument(6, tmpstr)
  read(tmpstr,*) scis%cond%e0
  call get_command_argument(7, tmpstr)
  read(tmpstr,*) scis%cond%edel

  call scissors_write(6, scis)

  SAFE_ALLOCATE(kp%elda, (kp%mnband, kp%nrk, kp%nspin))
  kp%elda = kp%el

  if (has_splines) then
    call get_command_argument(8, tmpstr)
    call open_file(20,file=TRUNC(tmpstr),form='formatted',status='old')
    read(20,*) spl%n
    SAFE_ALLOCATE(spl%t, (spl%n))
    SAFE_ALLOCATE(spl%c, (spl%n))
    read(20,*) (spl%t(ii), ii=1,spl%n)
    read(20,*) (spl%c(ii), ii=1,spl%n)
    read(20,*) spl%k
    close(20)
    call scissors_shift(kp, scis, spl)
  else
    call scissors_shift(kp, scis)
  endif

  call open_file(unit=9,file='eqp.dat',form='formatted',status='replace')

  do ik = 1, kp%nrk
    write(9,'(3f13.9,i8)') kp%rk(1:3, ik),kp%mnband*kp%nspin 
    do is = 1, kp%nspin 
      do ib = 1, kp%mnband
        write(9,'(2i8,2f15.9)') is, ib, kp%elda(ib, ik, is)*RYD, kp%el(ib, ik, is)*RYD
      enddo
    enddo
  enddo

  call close_file(9)

  call dealloc_header_type(sheader, crys, kp)
  SAFE_DEALLOCATE(kp%elda)

  if (has_splines) then
    SAFE_DEALLOCATE_P(spl%t)
    SAFE_DEALLOCATE_P(spl%c)
  endif

end program scissors2eqp
