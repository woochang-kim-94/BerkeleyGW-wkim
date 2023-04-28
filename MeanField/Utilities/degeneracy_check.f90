!=============================================================================
!
! Utilities:
!
! (1) degeneracy_check    Originally By DAS      Last Modified 12/13/2010 (DAS)
!
!     Determines numbers of bands that can be used, compatible with the
!     degenerate subspaces. Any number of wavefunctions can be given,
!     in new format.
!
!==============================================================================

#include "f_defs.h"

program degeneracy_check

  use global_m
  use wfn_rho_vxc_io_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp
  type(gspace) :: gvec
  character*3 :: sheader
  character*256 :: infile, usage
  integer :: ib, iflavor, minband, minvband, mincband, ik, is
  integer :: nargs, ifile
  logical, allocatable :: ok(:), ok_vband(:), ok_cband(:)

  usage = 'Usage: degeneracy_check.x wfn [wfn2 [...]]'

! Get file names from command-line arguments

  nargs = command_argument_count()

  if (nargs .lt. 1) then
    call die(usage)
  endif

  do ifile = 1, nargs

    call get_command_argument(ifile, infile)
    call open_file(unit=7,file=TRUNC(infile),form='unformatted',status='old')

    sheader = 'WFN'
    iflavor = -1
    call read_binary_header_type(7, sheader, iflavor, kp, gvec, syms, crys, warn = .false., dont_warn_kgrid = .true.)

    ! Output info
    write(6,'(a)') 'Reading eigenvalues from file ' // TRUNC(infile)
    write(6,'("Number of spins:",i16)') kp%nspin
    write(6,'("Number of bands:",i16)') kp%mnband
    write(6,'("Number of k-points:",i13)') kp%nrk

    kp%nvband=minval(kp%ifmax(:,:)-kp%ifmin(:,:))+1
    kp%ncband=kp%mnband-maxval(kp%ifmax(:,:))

    if(ifile == 1) then
      SAFE_ALLOCATE(ok,(kp%mnband - 1))
      ok(1:kp%mnband - 1) = .true.
      minband = kp%mnband

      SAFE_ALLOCATE(ok_vband,(kp%nvband - 1))
      ok_vband(:) = .true.
      minvband = kp%nvband

      if(kp%ncband - 1 > 0) then
        ! avoid allocation <= 0
        SAFE_ALLOCATE(ok_cband,(kp%ncband - 1))
        ok_cband(:) = .true.
      endif
      mincband = kp%ncband
    else
      minband = min(minband, kp%mnband)
      minvband = min(minvband, kp%nvband)
      mincband = min(mincband, kp%ncband)
    endif

    do ib = 1, minband - 1
      ok(ib) = ok(ib) .and. &
        all(abs(kp%el(ib, 1:kp%nrk, 1:kp%nspin) - kp%el(ib + 1, 1:kp%nrk, 1:kp%nspin)) .gt. TOL_Degeneracy)
    enddo

    do ib = 1, minvband - 1
      do is = 1, kp%nspin
        do ik = 1, kp%nrk
          ok_vband(ib) = ok_vband(ib) .and. &
            (abs(kp%el(kp%ifmax(ik, is) - ib + 1, ik, is) &
            - kp%el(kp%ifmax(ik, is) - ib, ik, is)) .gt. TOL_Degeneracy)
        enddo
      enddo
    enddo

    do ib = 1, mincband - 1
      do is = 1, kp%nspin
        do ik = 1, kp%nrk
          ok_cband(ib) = ok_cband(ib) .and. &
            (abs(kp%el(kp%ifmax(ik, is) + ib, ik, is) &
            - kp%el(kp%ifmax(ik, is) + ib + 1, ik, is)) .gt. TOL_Degeneracy)
        enddo
      enddo
    enddo

    call dealloc_header_type(sheader, crys, kp)
    call close_file(7)
  enddo

  if(nargs > 1) write(6,'(a,i6)') 'Minimum number of bands in files: ', minband

  write(6,'(a)')
  write(6,'(a)') '== Degeneracy-allowed numbers of bands (for epsilon and sigma) =='

  do ib = 1, minband - 1
    if(ok(ib)) then
      write(6,*) ib
    endif
  enddo

  write(6,'(a,i6,a)') 'Note: cannot assess whether or not highest band ', minband, ' is degenerate.'

  write(6,'(a)')
  write(6,'(a)') '== Degeneracy-allowed numbers of valence bands (for inteqp, kernel, and absorption) =='

  do ib = 1, minvband - 1
    if(ok_vband(ib)) then
      write(6,*) ib
    endif
  enddo
  write(6,*) minvband ! using all bands is always allowed ... ?

  if(mincband - 1 > 0) then
    write(6,'(a)')
    write(6,'(a)') '== Degeneracy-allowed numbers of conduction bands (for inteqp, kernel, and absorption) =='
    
    do ib = 1, mincband - 1
      if(ok_cband(ib)) then
        write(6,*) ib
      endif
    enddo
    
    write(6,'(a,i6,a)') 'Note: cannot assess whether or not highest conduction band ', mincband, ' is degenerate.'
    SAFE_DEALLOCATE(ok_cband)
  endif

  SAFE_DEALLOCATE(ok)
  SAFE_DEALLOCATE(ok_vband)

end program degeneracy_check
