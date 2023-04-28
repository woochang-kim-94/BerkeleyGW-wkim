!==============================================================================
!
! Routine:
!
! (1) inread        Originally By ?
!
!     Read information from input file plotxct.inp.
!
!==============================================================================

#include "f_defs.h"

module inread_m

  use global_m
  use plotxct_common_m
  implicit none
  
  public :: inread

contains

subroutine inread(xct,pxct)
  type (xctinfo), intent(out) :: xct
  type (plotxct_t), intent(out) :: pxct

  character*256 :: blockword,keyword,line,errmsg
  integer :: ii,iostat
  
  PUSH_SUB(inread)

#ifdef MPI
  ! Non-root nodes should wait for root to read the whole file.
  ! That way, we can be sure root gets a chance to write errors before
  ! any die call is issued by another node. Root calls MPI_Barrier below.
  if(peinf%inode /= 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  call open_file(8,file='plotxct.inp',form='formatted',status='old')

!     Set default values
  pxct%nsuper(:)=0
  pxct%pstate=0
  pxct%plot_hole=.false.
  pxct%plot_electron=.false.
  pxct%ispin=1
  pxct%nspinor=1
  pxct%rhole(:)=0.d0
  pxct%iwann=0
  pxct%seedname=''
  xct%shift(:)=0.d0
  xct%freplacebz=.false.
  xct%fwritebz=.false.
  xct%nn=0
  pxct%int_dim(:)=.false.
  pxct%unfold=.true.
  pxct%unfoldq=.true.
  pxct%downsample(:)=2
  pxct%bz_paranoid=.true.
  pxct%only_psi2=.false.
  pxct%hspinor=1
  pxct%espinor=1
  pxct%use_wfn_hdf5=.false.
  xct%wfn_hdf5_min_band_block = -1

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
    elseif(trim(keyword).eq.'integrate_x') then
      pxct%int_dim(1) = .true.
    elseif(trim(keyword).eq.'integrate_y') then
      pxct%int_dim(2) = .true.
    elseif(trim(keyword).eq.'integrate_z') then
      pxct%int_dim(3) = .true.
    elseif(trim(keyword).eq.'only_psi2') then
      pxct%only_psi2 = .true.
    elseif(trim(keyword).eq.'use_wfn_hdf5') then
      pxct%use_wfn_hdf5 = .true.
    elseif(trim(keyword).eq.'downsample') then
      read(line,*,err=110) pxct%downsample
!#BEGIN_INTERNAL_ONLY
    elseif(trim(keyword).eq.'no_bz_check') then
      pxct%bz_paranoid = .false.
    elseif(trim(keyword).eq.'plot_hole') then
      pxct%plot_hole = .true.
    elseif(trim(keyword).eq.'plot_electron') then
      pxct%plot_electron = .true.
    elseif(trim(keyword).eq.'wannier_index') then
      if(any(abs(pxct%rhole(1:3)) > TOL_Zero)) &
        call die("Keywords hole_position and wannier_index are incompatible.", only_root_writes = .true.)
      read(line,*,err=110) pxct%iwann
    elseif(trim(keyword).eq.'wannier_seedname') then
      if(any(abs(pxct%rhole(1:3)) > TOL_Zero)) &
        call die("Keywords hole_position and wannier_seedname are incompatible.", only_root_writes = .true.)
      read(line,*,err=110) pxct%seedname
!#END_INTERNAL_ONLY
    elseif(trim(keyword).eq.'plot_state') then
      read(line,*,err=110) pxct%pstate
    elseif(trim(keyword).eq.'plot_spin') then
      read(line,*,err=110) pxct%ispin
    elseif(trim(keyword).eq.'spinor') then
      pxct%nspinor = 2
    elseif(trim(keyword).eq.'hole_spin') then
      read(line,*,err=110) pxct%hspinor
    elseif(trim(keyword).eq.'electron_spin') then
      read(line,*,err=110) pxct%espinor
    elseif(trim(keyword).eq.'hole_position') then
      if(pxct%iwann /= 0) call die("Keywords hole_position and wannier_index are incompatible.", only_root_writes = .true.)
      read(line,*,err=110) (pxct%rhole(ii),ii=1,3)
    elseif(trim(keyword).eq.'q_shift') then
      read(line,*,err=110) (xct%shift(ii),ii=1,3)
    elseif(trim(keyword).eq.'supercell_size') then
      read(line,*,err=110) (pxct%nsuper(ii),ii=1,3)
    elseif(trim(keyword).eq.'fullbz_replace') then
      xct%freplacebz=.true.
    elseif(trim(keyword).eq.'fullbz_write') then
      xct%fwritebz=.true.
    else
      write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found in plotxct.inp.'
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
    if (pxct%int_dim(1)) then
      write(6,'(1x,a)') 'Integrating wavefunction modulus squared in x-direction'
    endif
    if (pxct%int_dim(2)) then
      write(6,'(1x,a)') 'Integrating wavefunction modulus squared in y-direction'
    endif
    if (pxct%int_dim(3)) then
      write(6,'(1x,a)') 'Integrating wavefunction modulus squared in z-direction'
    endif
    if (pxct%only_psi2) then
      write(6,'(1x,a)') 'We will only output the exciton wavefunction modulus squared'
    endif
    if (.not.pxct%bz_paranoid) then
      write(6,'(1x,a)') 'We won`t be paranoid when unfolding the BZ'
    endif
    write(6,'(1x,a,3(i2,1x))') 'Real-space grid will be downsampled by:', pxct%downsample
  endif

  if (any(pxct%downsample<1)) then
    call die('Invalid value for real-space grid subsampling.', only_root_writes = .true.)
  endif

  if (.not.pxct%only_psi2.and.any(pxct%int_dim)) then
    call die('integrate_* flags require only_psi2 flag.', only_root_writes = .true.)
  endif
  if (all(pxct%int_dim)) then
    call die('It makes no sense to integrate all dimensions!', only_root_writes = .true.)
  endif

  if(pxct%pstate.eq.0) call die("The plot_state keyword could not be found.", only_root_writes = .true.)
  if(pxct%ispin < 0 .or. pxct%ispin > 2) call die("plot_spin out of bounds", only_root_writes = .true.)
  ! we don`t know if we are spin-polarized yet so we cannot check that consistency here

  if(pxct%iwann /= 0 .and. trim(pxct%seedname) == '') &
    call die("If wannier_index is specified, wannier_seedname must be set.", only_root_writes = .true.)

  if(pxct%iwann == 0 .and. trim(pxct%seedname) /= '') &
    call die("If wannier_seedname is set, wannier_index must be specified.", only_root_writes = .true.)
  
  if(product(pxct%nsuper(:))==0) then
    if(peinf%inode == 0) write(0,'(a,3i6)') 'supercell = ', pxct%nsuper(1:3)
    call die("supercell cannot be zero", only_root_writes = .true.)
  endif
  
  if (pxct%ispin .ne. 1 .and. pxct%nspinor .eq. 2) then
    call die('Only one spin index when using spinor wavefunctions', only_root_writes = .true.)
  endif

  if (pxct%iwann .ne. 0 .and. pxct%nspinor .eq. 2) then
    call die('PlotXct with spinors and Wannier functions not implemented', only_root_writes = .true.)
  endif

  if (pxct%nspinor .ne. 2 .and. (pxct%espinor .eq. 2 .or. pxct%hspinor .eq. 2)) then
    call die('Use keyword spinor in input file and make sure wavefunctions have spinor components', only_root_writes = .true.)
  endif

  if(peinf%inode.eq.0) write(6,*)
  
#ifdef MPI
  ! root lets the others go after it is done reading (see beginning of function)
  if(peinf%inode == 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  POP_SUB(inread)
  
  return
  
110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', &
      trim(keyword), '. '
  call die(errmsg, only_root_writes = .true.)
  
end subroutine inread

end module inread_m
