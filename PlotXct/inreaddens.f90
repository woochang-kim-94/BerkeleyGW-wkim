!==============================================================================
!
! Routine:
!
! (1) inreaddens        Originally By Peter Doak 
!
!     Read information from input file plotxctdens.inp.
!
!==============================================================================

#include "f_defs.h"

module inreaddens_m

  use global_m
  use plotxct_common_m
  implicit none
  
  public :: inreaddens

contains

subroutine inreaddens(xct,pxct)
  type (xctinfo), intent(out) :: xct
  type (plotxct_t), intent(out) :: pxct

  character*256 :: blockword,keyword,line,errmsg
  integer :: ii,iostat
  
  PUSH_SUB(inreaddens)

#ifdef MPI
  ! Non-root nodes should wait for root to read the whole file.
  ! That way, we can be sure root gets a chance to write errors before
  ! any die call is issued by another node. Root calls MPI_Barrier below.
  if(peinf%inode /= 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  call open_file(8,file='plotxctdens.inp',form='formatted',status='old')

!     Set default values
  pxct%nsuper(:)=0
  pxct%pstateHigh=0
  pxct%pstateLow=0
  pxct%ispin=1
  xct%shift(:)=0.d0
  xct%freplacebz=.false.
  xct%fwritebz=.false.
  xct%nn=0
  pxct%int_dim(:)=.false.
  pxct%unfold=.true.
  pxct%unfoldq=.true.
  pxct%downsample(:)=1
  pxct%bz_paranoid=.true.


!     Never ending loop...
  do while(0.eq.0)

!       Actually the loop ends when the end of the file is reached
    read(8,'(a256)',iostat=iostat) line
    if(iostat < 0) exit

!       Skip comment lines
    if(len_trim(line).eq.0) cycle
    if(line(1:1).eq.'#') cycle

!       Determine keyword:
    keyword=line(1:scan(line," ")-1)
    line=adjustl(line(scan(line," ")+1:256))

    if(trim(keyword).eq.'begin') then
      blockword=line(1:scan(line," ")-1)
      ii=0
      do while(trim(line).ne.'end')
        read(8,'(a256)',iostat=iostat) line
        if(iostat /= 0) then
          write(errmsg,'(3a)') 'The end of the file was reached while reading elements of the ', &
            trim(blockword),' block.'
           call die(errmsg, only_root_writes = .true.)
        endif
        if(trim(line).ne.'end') then
          write(errmsg,'(3a)') 'Unexpected blockword ', trim(blockword), ' was found in plotxct.inp.'
           call die(errmsg, only_root_writes = .true.)
        end if
      end do
    elseif(trim(keyword).eq.'verbosity') then
      read(line,*,err=110) peinf%verbosity
    elseif(trim(keyword).eq.'restrict_kpoints') then
      read(line,*,err=110) xct%nn
    elseif(trim(keyword).eq.'no_symmetries_fine_grid') then
      pxct%unfold = .false.
    elseif(trim(keyword).eq.'use_symmetries_fine_grid') then
      pxct%unfold = .true.
    elseif(trim(keyword).eq.'no_symmetries_shifted_grid') then
      pxct%unfoldq = .false.
    elseif(trim(keyword).eq.'use_symmetries_shifted_grid') then
      pxct%unfoldq = .true.
!#BEGIN_INTERNAL_ONLY
    elseif(trim(keyword).eq.'no_bz_check') then
      pxct%bz_paranoid = .false.
!#END_INTERNAL_ONLY
    elseif(trim(keyword).eq.'plot_state_low') then
      read(line,*,err=110) pxct%pstateLow
    elseif(trim(keyword).eq.'plot_state_high') then
      read(line,*,err=110) pxct%pstateHigh
    elseif(trim(keyword).eq.'plot_spin') then
      read(line,*,err=110) pxct%ispin
    elseif(trim(keyword).eq.'q_shift') then
      read(line,*,err=110) (xct%shift(ii),ii=1,3)
    elseif(trim(keyword).eq.'supercell_size') then
      read(line,*,err=110) (pxct%nsuper(ii),ii=1,3)
    elseif(trim(keyword).eq.'fullbz_replace') then
      xct%freplacebz=.true.
    elseif(trim(keyword).eq.'fullbz_write') then
      xct%fwritebz=.true.
    else
      write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found in plotxctdens.inp.'
      call die(errmsg, only_root_writes = .true.)
    end if
  enddo
  
  call close_file(8)
  
  call peinfo_set_verbosity()
  if (peinf%inode==0) then
    if (pxct%unfold) then
      write(6,'(1x,a)') 'Unfolding fine BZ'
    else
      write(6,'(1x,a)') 'Not unfolding fine BZ'
    endif
    if (pxct%unfoldq) then
      write(6,'(1x,a)') 'Unfolding shifted BZ'
    else
      write(6,'(1x,a)') 'Not unfolding shifted BZ'
    endif
    if (.not.pxct%bz_paranoid) then
      write(6,'(1x,a)') 'We won`t be paranoid when unfolding the BZ'
    endif
    write(6,'(1x,a,3(i2,1x))') 'Real space grid will be downsampled by:', pxct%downsample
  endif

  if(pxct%ispin < 0 .or. pxct%ispin > 2) call die("plot_spin out of bounds", only_root_writes = .true.)
  ! we don`t know if we are spin-polarized yet so we cannot check that consistency here
  
  if(peinf%inode.eq.0) write(6,*)
  
#ifdef MPI
  ! root lets the others go after it is done reading (see beginning of function)
  if(peinf%inode == 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  POP_SUB(inreaddens)
  
  return
  
110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', &
      trim(keyword), '. '
  call die(errmsg, only_root_writes = .true.)
  
112 write(errmsg,'(3a)') 'Unexpected characters were found while reading elements of the ', &
      trim(blockword),' block.'
  call die(errmsg, only_root_writes = .true.)

end subroutine inreaddens

end module inreaddens_m

