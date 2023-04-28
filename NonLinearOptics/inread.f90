!=======================================================================================
!
! Routines:
!
! (1) inread()          Originally By                   Last Modified 6/6/2008 (JRD)
!
!     Read input parameters from file nonlinearoptics.inp
!
!     input: none
!
!     output: xct%nvb_fi        number of valence bands (fine grid)
!             xct%ncb_fi        number of conduction bands (fine grid)
!             xct%eta           energy resolution (used to calculate the
!                 optical spectrum)
!             xct%read_kpoints          same as in Kernel
!             xct%nktotal       total number of points in fine grid
!                  (if read_kpoints, then xct%nktotal=xct%nkpt_fi, set later in input.f90)
!             xct%shift(1:3)    shift vector (this is the small shift, used
!                               to generate WFNq_fi)
!
!             neig : number of eigenvectors/eigenvalues to be computed
!
!========================================================================================

#include "f_defs.h"

module inread_m

  use global_m
  use input_utils_m
  use scissors_m
  implicit none

  private

  public :: &
    inread

contains

subroutine inread(eqp,xct,flag,nmax,neig,win,other)
  type (eqpinfo), intent(out) :: eqp
  type (xctinfo), intent(out) :: xct
  type (flags), intent(out) :: flag
  type (windowinfo), intent(out) :: win
  type (otherinfo), intent(out) :: other
  integer, intent(out) :: nmax,neig

  character*256 :: blockword,keyword,line,errmsg
  integer :: ii,iostat
  logical :: found

  integer, parameter :: nmaxr=10000
  integer, allocatable :: istates_read(:)
  real(DP), allocatable :: cstates_read(:), estates_read(:)
  
  PUSH_SUB(inread)

#ifdef MPI
  ! Non-root nodes should wait for root to read the whole file.
  ! That way, we can be sure root gets a chance to write errors before
  ! any die call is issued by another node. Root calls MPI_Barrier below.
  if(peinf%inode /= 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  SAFE_ALLOCATE(istates_read, (nmaxr))
  SAFE_ALLOCATE(cstates_read, (nmaxr))
  SAFE_ALLOCATE(estates_read, (nmaxr))

! JRD: Dumb debugging

!      if (peinf%inode .eq. 0) then
!        write(6,*) 'Inreading!'
!        write(6,*) 'Bout to start default values'
!      endif

!------------------------------
! Set default values

  other%ithreeD=0
  other%keta=0D0
  other%knx=0
  other%kny=0
  other%knz=0
  xct%nvb_fi=0
  xct%ncb_fi=0
  call scissors_zero(eqp%scis)
  flag%opr=999
  flag%job=1
  flag%vm=0
  flag%bz0=1
  flag%bzq=1
  flag%lor=0
  nmax=0
  xct%read_kpoints = .false.
  xct%shift(:)=0.d0
  xct%pol(:)=0.d0
  xct%nktotal=0
  neig=0
  xct%vmin=1
  xct%vmax=0
  xct%rgrid=0
  win%nstates=0
  xct%freplacebz=.false.
  xct%fwritebz=.false.
  xct%efermi_input=0.0d0      
  xct%rfermi=.true.
#ifdef HDF5
  xct%use_hdf5=.true.
#else
  xct%use_hdf5=.false.
#endif


! JRD: Dumb debugging

!      write(6,*) 'Starting read loop'

!-----------------------------------
! Never ending loop...

  do while(0.eq.0)

! Actually the loop ends when the end of the file is reached

    read(8,'(a256)',iostat=iostat) line
    if(iostat < 0) exit
    
! Skip comment lines

    if(len_trim(line).eq.0) cycle
    if(line(1:1).eq.'#') cycle

! Determine keyword:

    keyword=line(1:scan(line," ")-1)
    line=adjustl(line(scan(line," ")+1:256))
    
    if(trim(keyword).eq.'begin') then
      blockword=line(1:scan(line," ")-1)
      ii=0
      do while(trim(line).ne.'end')
        read(8,'(a256)',end=105) line
        if(trim(line).ne.'end') then
          ii=ii+1
          if(trim(blockword).eq.'initial_states') then
            read(line,*,iostat=iostat) istates_read(ii), estates_read(ii),cstates_read(ii)
            if(iostat /= 0) then
              write(errmsg,'(3a)') 'Unexpected characters while reading elements of the ', &
                trim(blockword),' block. '
              call die(errmsg, only_root_writes = .true.)
            endif
          else
            write(errmsg,'(3a)') 'Unexpected blockword ', trim(blockword), ' was found in nonlinearoptics.inp.'
            call die(errmsg, only_root_writes = .true.)
          end if
        end if
      end do
      if(trim(blockword).eq.'initial_states') then
        win%nstates=ii
      endif
    elseif(trim(keyword).eq.'verbosity') then
      read(line,*,err=110) peinf%verbosity
    elseif(trim(keyword).eq.'number_of_final_states') then
      read(line,*,err=110) win%nstates
    elseif(trim(keyword).eq.'number_val_bands_fine') then
      read(line,*,err=110) xct%nvb_fi
    elseif(trim(keyword).eq.'number_cond_bands_fine') then
      read(line,*,err=110) xct%ncb_fi
    elseif(trim(keyword).eq.'absorp3D') then
      other%ithreeD=1
    elseif(trim(keyword).eq.'keta') then
      read(line,*,err=110) other%keta
    elseif(trim(keyword).eq.'knx') then
      read(line,*,err=110) other%knx
    elseif(trim(keyword).eq.'kny') then
      read(line,*,err=110) other%kny
    elseif(trim(keyword).eq.'knz') then
      read(line,*,err=110) other%knz
!          write(6,*) 'Read knz'
    elseif(trim(keyword).eq.'energy_max') then
      read(line,*,err=110) win%emax
    elseif(trim(keyword).eq.'energy_resolution') then
      read(line,*,err=110) xct%eta
    elseif(trim(keyword).eq.'two_photon_job') then
      flag%job=1
    elseif(trim(keyword).eq.'ultrafast_job') then
      flag%job=0
    elseif(trim(keyword).eq.'use_velocity') then
      flag%opr=0
    elseif(trim(keyword).eq.'use_momentum') then
      flag%opr=1
!        elseif(trim(keyword).eq.'number_eigenvalues') then
!          read(line,*,err=110) neig
    elseif(trim(keyword).eq.'q_shift') then
      read(line,*,err=110) (xct%shift(ii),ii=1,3)
    elseif(trim(keyword).eq.'polarization') then
      read(line,*,err=110) (xct%pol(ii),ii=1,3)
    elseif(trim(keyword).eq.'read_vmtxel_nl') then
      flag%vm=1
    elseif(trim(keyword).eq.'no_symmetries_fine_grid') then
      flag%bz0 = 1
    elseif(trim(keyword).eq.'use_symmetries_fine_grid') then
      flag%bz0 = 0
    elseif(trim(keyword).eq.'no_symmetries_shifted_grid') then
      flag%bzq = 1
    elseif(trim(keyword).eq.'use_symmetries_shifted_grid') then
      flag%bzq = 0
    elseif(trim(keyword).eq.'lorentzian_broadening') then
      flag%lor=1
    elseif(trim(keyword).eq.'gaussian_broadening') then
      flag%lor=0
    elseif(trim(keyword).eq.'fullbz_replace') then
      xct%freplacebz=.true.
    elseif(trim(keyword).eq.'fullbz_write') then
      xct%fwritebz=.true.
    elseif(trim(keyword).eq.'dont_use_hdf5') then
      xct%use_hdf5 = .false.
    else
      call scissors_inread(keyword, line, eqp%scis, found)
      if(.not. found) then
        write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found in nonlinearoptics.inp.'
        call die(errmsg, only_root_writes = .true.)
      endif
    end if
  enddo

  if(xct%nvb_fi.eq.0) then
    call die('The number_val_bands_fine keyword could not be found.', only_root_writes = .true.)
  endif
  
  if(xct%ncb_fi.eq.0) then
    call die('The number_cond_bands_fine keyword could not be found.', only_root_writes = .true.)
  endif
  
  if(win%nstates.ne.0 .and. flag%job .eq. 0) then
    SAFE_ALLOCATE(win%istates, (win%nstates))
    SAFE_ALLOCATE(win%cstates, (win%nstates))
    SAFE_ALLOCATE(win%estates, (win%nstates))
    do ii = 1, win%nstates
      win%istates(ii)=istates_read(ii)
      win%cstates(ii)=cstates_read(ii)
      win%estates(ii)=estates_read(ii)
    enddo
  endif
  
  if(flag%opr.eq.0) then
    if(peinf%inode.eq.0) then
      write(0,*) 'Using the velocity operator '
      write(0,*)
      write(0,*) '************************************************'
      write(0,*) '************************************************'
      write(0,*) '**                                            **'
      write(0,*) '**                WARNING !!!!                **'
      write(0,*) '**                                            **'
      write(0,*) '**      VELOCITY OPERATOR DOESN"T WORK YET    **'
      write(0,*) '**             IN TWO-PHOTON CODES            **'
      write(0,*) '**                                            **'
      write(0,*) '**         YOUR ANSWER WILL BE GARBAGE        **'
      write(0,*) '**                                            **'
      write(0,*) '************************************************'
      write(0,*) '************************************************'
      write(0,*)
    endif
  elseif(flag%opr.eq.1) then
    if(peinf%inode.eq.0) then
      write(6,*) 'Using the momentum operator '
    endif
  else
    call die("No flag for velocity/momentum operator", only_root_writes = .true.)
  endif
  
  if(xct%eta.eq.0.d0) then
    call die('The energy_resolution keyword could not be found.', only_root_writes = .true.)
  endif
  
  if(xct%read_kpoints .and. xct%nktotal.eq.0) then
    call die("missing total number of unit cells")
  endif
  
  if (peinf%inode .eq. 0) then
    if (flag%job .eq. 0) then
      write(6,*) ' '
      write(6,*) 'Doing an Ultrafast Calculation.'
    else
      write(6,*) ' '
      write(6,*) 'Doing a Two-Photon Calculation.'
    endif
  end if
  
  if(peinf%inode.eq.0) write(6,*)

! JRD: More dumb debugging

!      write(6,*) 'Inread nkpt,ncband,nvband,nspin,neig',
!     > xct%nkpt_fi,xct%ncb_fi,xct%nvb_fi,xct%nspin,neig

  SAFE_DEALLOCATE(istates_read)
  SAFE_DEALLOCATE(cstates_read)
  SAFE_DEALLOCATE(estates_read)

#ifdef MPI
  ! root lets the others go after it is done reading (see beginning of function)
  if(peinf%inode == 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  POP_SUB(inread)
  
  return
  
105 write(errmsg,'(3a)') 'The end of the file was reached while reading elements of the ', &
      trim(blockword),' block. '
  call die(errmsg, only_root_writes = .true.)
  
110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', &
      trim(keyword), '. '
  call die(errmsg, only_root_writes = .true.)
  
end subroutine inread

end module inread_m
