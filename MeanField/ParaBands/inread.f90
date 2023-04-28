!==============================================================================
!
! Routine:
!
! (1) inread        Originally By ?
!
!     Read information from input file parabands.inp.
!
!==============================================================================

#include "f_defs.h"

module inread_m

  use global_m
  use kpoint_pool_m,   only: kpoint_pool_t
  implicit none

  private

    type pb_params_t
      integer :: cur_ik = 0
      character(len=256) :: fname_wfn_in = ''
      character(len=256) :: fname_wfn_out = ''
      character(len=256) :: fname_vsc = ''
      character(len=256) :: fname_vkb = ''
      character(len=256) :: fname_qkb = ''
      logical :: has_wfn_in = .false.
      logical :: has_wfn_out = .false.
      logical :: has_vsc = .false.
      logical :: has_vkb = .false.
      logical :: has_qkb = .false.
      integer :: nb = -1
      integer :: block_sz = 64
      integer :: solver_alg = -1
      integer :: wfn_io_driver = 0
      integer :: wfn_io_mpiio_mode = 0
      type(kpoint_pool_t) :: kpp
      real(DP) :: scalapack_mem_rescale = 1d0
      real(DP) :: primme_tol = 1d-8
      real(DP) :: primme_ritz_shift_min = 1d0/TOL_ZERO
      real(DP) :: primme_ritz_shift = 0.1d0
      integer :: primme_exact_diag_algo = -1
      integer :: primme_exact_diag_size = 1000
      integer :: primme_max_block_sz = 8
      integer :: primme_max_basis_sz = 256
      logical :: has_restart = .false.
    endtype pb_params_t

  public :: pb_params_t, inread

contains

subroutine inread(params)
  type (pb_params_t), intent(out) :: params

  character(len=256) :: keyword, line, errmsg
  integer :: iostat

  PUSH_SUB(inread)

#ifdef MPI
  ! Non-root nodes should wait for root to read the whole file.
  ! That way, we can be sure root gets a chance to write errors before
  ! any die call is issued by another node. Root calls MPI_Barrier below.
  if (peinf%inode/=0 .and. peinf%npes>1) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  call open_file(8, file='parabands.inp', form='formatted', status='old')

!     Set default values
  params%kpp%inode = peinf%inode
  params%kpp%npes = peinf%npes
  params%kpp%npes_diag = peinf%npes

!     Never ending loop...
  do while(.true.)

!       Actually the loop ends when the end of the file is reached
    read(8,'(a256)',iostat=iostat) line
    if(iostat < 0) exit

!       Skip comment lines
    if(len_trim(line).eq.0) cycle
    if(line(1:1).eq.'#') cycle

!       Determine keyword:
    keyword=line(1:scan(line," ")-1)
    line=adjustl(line(scan(line," ")+1:256))

    if(trim(keyword).eq.'input_wfn_file') then
      read(line,*,err=110) params%fname_wfn_in
      params%has_wfn_in = .true.
    elseif(trim(keyword).eq.'output_wfn_file') then
      read(line,*,err=110) params%fname_wfn_out
      params%has_wfn_out = .true.
    elseif(trim(keyword).eq.'vsc_file') then
      read(line,*,err=110) params%fname_vsc
      params%has_vsc = .true.
    elseif(trim(keyword).eq.'vkb_file') then
      read(line,*,err=110) params%fname_vkb
      params%has_vkb = .true.
    elseif(trim(keyword).eq.'qkb_file') then
      read(line,*,err=110) params%fname_qkb
      params%has_qkb = .true.
    elseif(trim(keyword).eq.'number_bands') then
      read(line,*,err=110) params%nb
    elseif(trim(keyword).eq.'number_pools') then
      read(line,*,err=110) params%kpp%npools
    elseif(trim(keyword).eq.'block_size') then
      read(line,*,err=110) params%block_sz
    elseif(trim(keyword).eq.'solver_algorithm') then
      read(line,*,err=110) params%solver_alg
    elseif(trim(keyword).eq.'scalapack_mem_rescale') then
      read(line,*,err=110) params%scalapack_mem_rescale
    !elseif(trim(keyword).eq.'wfn_io_driver') then
    !  read(line,*,err=110) params%wfn_io_driver
    elseif(trim(keyword).eq.'wfn_io_mpiio_mode') then
      read(line,*,err=110) params%wfn_io_mpiio_mode
    elseif(trim(keyword).eq.'tolerance') then
      read(line,*,err=110) params%primme_tol
    elseif(trim(keyword).eq.'ritz_shift') then
      read(line,*,err=110) params%primme_ritz_shift
    elseif(trim(keyword).eq.'exact_diag_algo') then
      read(line,*,err=110) params%primme_exact_diag_algo
    elseif(trim(keyword).eq.'exact_diag_size') then
      read(line,*,err=110) params%primme_exact_diag_size
    elseif(trim(keyword).eq.'max_block_size') then
      read(line,*,err=110) params%primme_max_block_sz
    elseif(trim(keyword).eq.'max_basis_size') then
      read(line,*,err=110) params%primme_max_basis_sz
    elseif(trim(keyword).eq.'restart') then
      params%has_restart = .true.
    elseif(trim(keyword).eq.'verbosity') then
      read(line,*,err=110) peinf%verbosity
    else
      write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found in parabands.inp.'
      call die(errmsg, only_root_writes = .true.)
    end if
  enddo
  call close_file(8)

  call peinfo_set_verbosity()
  if (peinf%inode==0) then
    write(6,*)
    if (.not.params%has_wfn_in) call die('Input WFN file was not specified')
    if (.not.params%has_vsc) call die('Input VSC file was not specified')
#ifndef HDF5
    write(0,'(/a)') 'WARNING: Code was not compiled with HDF5 support.'
    write(0,'(a/)') '         We will not be able to write the output WFN file!'
#endif
    write(6,'(1x,a)') 'Options for writing output WFN file:'
    select case (params%wfn_io_driver)
      case (0)
        write(6,'(1x,a)') '- Using regular HDF5 write operation'
      case (1)
        write(6,'(1x,a)') '- Using two-stage process'
      case default
        call die('Invalid option for wfn_io_driver', only_root_writes=.true.)
    endselect
    select case (params%wfn_io_mpiio_mode)
      case (0)
        write(6,'(1x,a)') '- Using collective MPI-IO'
      case (1)
        write(6,'(1x,a)') '- Using independent MPI-IO'
      case default
        call die('Invalid option for wfn_io_mpiio_mode', only_root_writes=.true.)
    endselect
    write(6,'()')
  endif

  ! FHJ: FIXME! Pools are not working!
  if (params%kpp%npools/=1) then
    params%kpp%npools = 1
    if (peinf%inode==0) write(0,'(/1x,a/)') &
      'WARNING: k-point pools are broken. Overwriting `number_pools` to 1.'
  endif

#ifndef HDF5
    params%has_wfn_out = .false.
#endif
#ifdef MPI
  ! root lets the others go after it is done reading (see beginning of function)
  if (peinf%inode==0 .and. peinf%npes>1) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  POP_SUB(inread)

  return

110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', &
      trim(keyword), '. '
  call die(errmsg, only_root_writes = .true.)

end subroutine inread

end module inread_m

