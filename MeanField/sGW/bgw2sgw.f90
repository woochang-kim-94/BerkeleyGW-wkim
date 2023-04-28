!===============================================================================
!
! bgw2sgw           Originally by FHJ             Last Modified Jul/19 (FHJ)
!
!     Converts wavefunctions from BerkeleyGW to StochasticGW format
!
!===============================================================================

#include "f_defs.h"

program bgw2sgw
  use fftw_m
  use global_m
#ifdef HDF5
  use hdf5
#endif
  use io_utils_m
  use wfn_io_hdf5_m
  use wfn_rho_vxc_io_m

  implicit none

  character(len=128) :: input_rho_file, output_rho_file, &
                        input_wfn_file, output_wfn_file, &
                        output_structure_file
  logical :: has_rho, has_wfn, write_wfn_bin, write_bin
  character(len=18) :: fmtstr
  character(len=18), parameter :: &
    fmtstr_low = '(*(e13.07e2,:,1x))', &
    fmtstr_high = '(*(e22.16e2,:,1x))'
  integer :: number_states, ierr

#ifdef HDF5
  call h5open_f(ierr)
#endif

  call inread()
  if (has_rho) call convert_rho()
  if (has_wfn) call convert_wfn()

#ifdef HDF5
  call h5close_f(ierr)
#endif


contains

  subroutine inread
    character(len=256) :: keyword, line, errmsg
    integer :: iostat, len_wfn

    has_rho = .false.
    has_wfn = .false.
    output_rho_file = 'dens.txt'
    output_wfn_file = 'wf.bin'
    write_wfn_bin = .true.
    output_structure_file = 'cnt.ini'
    fmtstr = fmtstr_low
    number_states = -1

    write(*,*) 'Opening input file bgw2sgw.inp'

    call open_file(8, file='bgw2sgw.inp', form='formatted', status='old')
    do while (.true.)

      read(8, '(a256)', iostat=iostat) line
      if (iostat < 0) exit

      ! Skip comment lines
      if (len_trim(line).eq.0) cycle
      if (line(1:1).eq.'#') cycle

      ! Determine keyword:
      keyword = line(1:scan(line," ")-1)
      line = adjustl(line(scan(line," ")+1:256))

      if(trim(keyword).eq.'input_rho_file') then
        read(line,*,err=110) input_rho_file
        has_rho = .true.
      elseif(trim(keyword).eq.'output_rho_file') then
        read(line,*,err=110) output_rho_file
      elseif(trim(keyword).eq.'input_wfn_file') then
        read(line,*,err=110) input_wfn_file
        has_wfn = .true.
      elseif(trim(keyword).eq.'output_wfn_file') then
        read(line,*,err=110) output_wfn_file
      elseif(trim(keyword).eq.'output_structure_file') then
        read(line,*,err=110) output_structure_file
      elseif(trim(keyword).eq.'number_states') then
        read(line,*,err=110) number_states
      elseif(trim(keyword).eq.'precision_high') then
        fmtstr = fmtstr_high
      elseif(trim(keyword).eq.'precision_low') then
        fmtstr = fmtstr_low
      else
        write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found in bgw2sgw.inp.'
        call die(errmsg)
      endif
    enddo
    call close_file(8)

    if (.not.has_rho .and. .not.has_wfn) then
      call die('Must specify either input_rho_file or input_wfn_file (or both).')
    endif

    if (has_wfn) then
      len_wfn = len_trim(output_wfn_file)
      select case (output_wfn_file(len_wfn-3:len_wfn))
        case ('.bin', '.Bin', '.BIN')
          write_wfn_bin = .true.
        case ('.txt', '.Txt', '.TXT')
          write_wfn_bin = .false.
        case default
          write(0,*)
          write(0,*) 'WARNING: Could not defect output file format for wfn file.'
          write(0,*) 'Assuming a binary output.'
          write(0,*)
          write_wfn_bin = .true.
      end select
    endif

    write(*,*)
    if (has_rho) then
      write(*,'(1x,a,a)') 'Input rho file: ', trim(input_rho_file)
      write(*,'(1x,a,a)') 'Output rho file: ', trim(output_rho_file)
    endif
    if (has_wfn) then
      write(*,'(1x,a,a)') 'Input wfn file: ', trim(input_wfn_file)
      write(*,'(1x,a,a)') 'Output wfn file: ', trim(output_wfn_file)
      write(*,'(1x,a,i0)') 'Number of states to write: ', number_states
      write(*,'(1x,a,l1)') 'Write output wfn file in binary format: ', write_wfn_bin
      write(*,'(1x,a,a)') 'Output structure file: ', trim(output_structure_file)
    endif
    write(*,*)

    return
110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', &
      trim(keyword), '. '
    call die(errmsg)

  end subroutine inread


  !============================================================================
  ! Convert RHO -> dens.txt
  !============================================================================
  subroutine convert_rho
    type(mf_header_t) :: mf
    integer :: nx, ny, nz, ii, loc(3)
    integer, allocatable :: isort_iden(:)
    complex(DPC) :: phase
    complex(DPC), allocatable :: fftbox(:,:,:)
    complex(DPC), allocatable :: zrho(:,:)
    real(DP), allocatable :: drho(:,:)
    real(DP), allocatable :: rho_r(:,:,:)

    write_bin = .false.
    write(*,*) 'Converting charge-density file ', trim(input_rho_file), ' to ', trim(output_rho_file)
    call open_file(unit=7, file=trim(input_rho_file), form='unformatted', status='old')

    call read_mf_header(7, mf, sheader='RHO', warn=.false., dont_warn_kgrid=.true.)
    call check_file(mf)
    call write_structure(mf%crys)

    SAFE_ALLOCATE(mf%gvec%components, (3, mf%gvec%ng))
    call read_binary_gvectors(7, mf%gvec%ng, mf%gvec%ng, mf%gvec%components)

    if (mf%iflavor == 1) then
      SAFE_ALLOCATE(drho, (mf%gvec%ng, 1))
      call read_binary_data(7, mf%gvec%ng, mf%gvec%ng, mf%kp%nspin, drho)
    else
      SAFE_ALLOCATE(zrho, (mf%gvec%ng, 1))
      call read_binary_data(7, mf%gvec%ng, mf%gvec%ng, mf%kp%nspin, zrho)
    endif

    call close_file(7)

    call open_file(unit=11, file=trim(output_rho_file), form='formatted', status='replace')
    call write_sgw_header(mf)

    nx = mf%gvec%FFTgrid(1)
    ny = mf%gvec%FFTgrid(2)
    nz = mf%gvec%FFTgrid(3)
    SAFE_ALLOCATE(fftbox, (nx, ny, nz))
    SAFE_ALLOCATE(rho_r, (nx, ny, nz))
    SAFE_ALLOCATE(isort_iden, (mf%gvec%ng))
    isort_iden(:) = [(ii, ii=1,mf%gvec%ng)]

    ! Perform FFT
    if (mf%iflavor == 1) then
      call put_into_fftbox(mf%gvec%ng, drho(:,1), mf%gvec%components, &
        isort_iden, fftbox, mf%gvec%FFTgrid)
    else
      call put_into_fftbox(mf%gvec%ng, zrho(:,1), mf%gvec%components, &
        isort_iden, fftbox, mf%gvec%FFTgrid)
    endif
    call do_FFT(fftbox, mf%gvec%FFTgrid, 1)

    ! Make data real
    loc = maxloc(abs(fftbox))
    phase = fftbox(loc(1), loc(2), loc(3))
    phase = conjg(phase) / abs(phase)
    ! Put in normalization factor of 1/V
    phase = phase / mf%crys%celvol
    rho_r = dble(phase * fftbox)
    call write_label('dens')
    call write_real_array(rho_r, size(rho_r))

    SAFE_DEALLOCATE(isort_iden)
    SAFE_DEALLOCATE(rho_r)
    SAFE_DEALLOCATE(fftbox)
    if (mf%iflavor == 1) then
      SAFE_DEALLOCATE(drho)
    else
      SAFE_DEALLOCATE(zrho)
    endif
    call close_file(11)

    write(*,*) 'Finished converting charge-density file'

  end subroutine convert_rho


  !============================================================================
  ! Convert WFN -> wf.bin/wf.txt
  !============================================================================
  subroutine convert_wfn
    logical :: use_hdf5
    type(mf_header_t) :: mf
    integer :: nx, ny, nz, ii, ib, loc(3), ngk
    integer, allocatable :: isort_iden(:), gvecs_wfn(:,:)
    complex(DPC) :: phase
    complex(DPC), allocatable :: fftbox(:,:,:)
    complex(DPC), allocatable :: zwfn(:,:)
    real(DP), allocatable :: dwfn(:,:)
    real(DP), allocatable :: wfn_r(:,:,:)
    type(progress_info) :: prog_info !< a user-friendly progress report

    write_bin = write_wfn_bin
    write(*,*) 'Converting wavefunction file ', trim(input_wfn_file), ' to ', trim(output_wfn_file)

    ii = len_trim(input_wfn_file)
    use_hdf5 = .false.
    if (input_wfn_file(ii-2:ii)=='.h5' .or. input_wfn_file(ii-2:ii)=='.H5') then
#ifdef HDF5
      use_hdf5 = .true.
#else
      call die('Not compiled with HDF5 support!')
#endif
    endif

    if (use_hdf5) then
#ifdef HDF5
      call read_hdf5_mf_header(trim(input_wfn_file), mf, sheader='WFN', &
        warn=.false., dont_warn_kgrid=.true.)
#endif
    else
      call open_file(unit=7, file=trim(input_wfn_file), form='unformatted', status='old')
      call read_mf_header(7, mf, sheader='WFN', warn=.false., dont_warn_kgrid=.true.)
    endif

    call check_file(mf)
    call write_structure(mf%crys)
    ngk = mf%kp%ngk(1)

    ! Master list of G-vectors
    SAFE_ALLOCATE(mf%gvec%components, (3, mf%gvec%ng))
    ! WFN G vectors at Gamma
    SAFE_ALLOCATE(gvecs_wfn, (3, ngk))

    ! The code below only works because nk==1 here.
    if (use_hdf5) then
#ifdef HDF5
      call read_hdf5_gvectors(trim(input_wfn_file), mf%gvec%ng, mf%gvec%components)
      call read_hdf5_wfn_gvectors(trim(input_wfn_file), gvecs_wfn, ngk)
#endif
    else
      call read_binary_gvectors(7, mf%gvec%ng, mf%gvec%ng, mf%gvec%components)
      call read_binary_gvectors(7, ngk, ngk, gvecs_wfn)
    endif

    if (mf%iflavor == 1) then
      SAFE_ALLOCATE(dwfn, (ngk, 1))
    else
      SAFE_ALLOCATE(zwfn, (ngk, 1))
    endif

    if (write_bin) then
      call open_file(unit=11, file=trim(output_wfn_file), form='unformatted', status='replace')
    else
      call open_file(unit=11, file=trim(output_wfn_file), form='formatted', status='replace')
    endif
    call write_sgw_header(mf)

    if (number_states<0) number_states = mf%kp%mnband
    call write_int('nstates', number_states)
    call write_label('evls')
    ! Divide by two to write in Hartree instead of ryd.
    call write_real_array(mf%kp%el(1:number_states,1,1)/2d0, number_states)

    nx = mf%gvec%FFTgrid(1)
    ny = mf%gvec%FFTgrid(2)
    nz = mf%gvec%FFTgrid(3)
    SAFE_ALLOCATE(fftbox, (nx, ny, nz))
    SAFE_ALLOCATE(wfn_r, (nx, ny, nz))
    SAFE_ALLOCATE(isort_iden, (ngk))
    isort_iden(:) = [(ii, ii=1,ngk)]

    call write_label('orbitals')
    call progress_init(prog_info, 'conversion of wavefunctions', 'bands', number_states)
    do ib = 1, number_states
      call progress_step(prog_info)
      if (use_hdf5) then
#ifdef HDF5
        if (mf%iflavor == 1) then
          call read_hdf5_band(trim(input_wfn_file), dwfn(:,:), ngk, 1, 0, ib-1)
        else
          call read_hdf5_band(trim(input_wfn_file), zwfn(:,:), ngk, 1, 0, ib-1)
        endif
#endif
      else
        if (mf%iflavor == 1) then
          call read_binary_data(7, ngk, ngk, 1, dwfn(:,:))
        else
          call read_binary_data(7, ngk, ngk, 1, zwfn(:,:))
        endif
      endif

      ! Perform FFT
      if (mf%iflavor == 1) then
        call put_into_fftbox(ngk, dwfn(:,1), gvecs_wfn, isort_iden, fftbox, mf%gvec%FFTgrid)
      else
        call put_into_fftbox(ngk, zwfn(:,1), gvecs_wfn, isort_iden, fftbox, mf%gvec%FFTgrid)
      endif
      call do_FFT(fftbox, mf%gvec%FFTgrid, 1)

      ! Make data real
      loc = maxloc(abs(fftbox))
      phase = fftbox(loc(1), loc(2), loc(3))
      phase = conjg(phase) / abs(phase)
      ! Put in normalization factor of 1/sqrt(V)
      phase = phase / sqrt(mf%crys%celvol)
      wfn_r = dble(phase * fftbox)
      call write_int_array([ib, 1], 2)
      call write_real_array(wfn_r, size(wfn_r))
    enddo
    call progress_free(prog_info)

    SAFE_DEALLOCATE(isort_iden)
    SAFE_DEALLOCATE(wfn_r)
    SAFE_DEALLOCATE(fftbox)
    if (mf%iflavor == 1) then
      SAFE_DEALLOCATE(dwfn)
    else
      SAFE_DEALLOCATE(zwfn)
    endif
    if (.not.use_hdf5) then
      call close_file(7)
    endif
    call close_file(11)

    write(*,*) 'Finished converting wavefunction file'

  end subroutine convert_wfn


  !============================================================================
  ! Convert RHO/WFN -> cnt.ini
  !============================================================================
  subroutine write_structure(crys)
    !void get_symbol_by_z(int *zz, char symbol[3])
    interface
      subroutine get_symbol_by_z(zz, symbol) bind(c)
        use, intrinsic :: iso_c_binding
        integer(c_int), intent(in) :: zz
        character(c_char), intent(out) :: symbol(3)
      end subroutine get_symbol_by_z
    end interface

    type(crystal), intent(in) :: crys

    logical, save :: written = .false.
    character(len=3) :: symbol
    integer :: iat
    real(DP) :: cell(3)

    if (written) then
      return
    endif

    call open_file(unit=12, file=trim(output_structure_file), form='formatted', status='replace')
    write(*,*)
    write(*,*) 'Writing structure file ', trim(output_structure_file)

    cell = crys%alat * [crys%avec(1,1), crys%avec(2,2), crys%avec(3,3)]

    do iat = 1, crys%nat
      call get_symbol_by_z(crys%atyp(iat), symbol)
      write(12,'(a,3(1x,f15.9))') symbol, mod(crys%alat * crys%apos(:,iat), cell)
    enddo
    call close_file(12)

    write(*,*) 'Ok'
    write(*,*)
    written = .true.

  end subroutine write_structure


  !============================================================================
  ! Auxiliary functions
  !============================================================================
  subroutine check_file(mf)
    type(mf_header_t), intent(in) :: mf

    integer :: ii, jj
    logical :: full_check
    real(DP) :: diag_max, offdiag_max

    select case (mf%sheader)
      case ('RHO')
        full_check = .false.
      case ('WFN')
        full_check = .true.
      case default
        call die('Unknown mf header: '//mf%sheader)
    end select

    if (mf%kp%nspin/=1) call die('Only support nspin=1')
    if (mf%kp%nspinor/=1) call die('Only support nspinor=1')
    if (full_check) then
      if (mf%kp%nrk/=1) call die('Only support a single k-point')
      if (mf%kp%ngk(1)/=mf%kp%ngkmax) call die('Inconsistencies in ngkmax')
    endif

    diag_max = maxval(abs([mf%crys%avec(1,1), mf%crys%avec(2,2), mf%crys%avec(3,3)]))
    offdiag_max = 0d0
    do jj = 1, 3
      do ii = 1,3
        if (ii/=jj) offdiag_max = max(offdiag_max, mf%crys%avec(ii,jj))
      enddo
    enddo
    if (offdiag_max/diag_max > TOL_SMALL) then
      call die('Only support orthogonal lattice vectors!')
    endif

  end subroutine check_file


  subroutine write_sgw_header(mf)
    type(mf_header_t), intent(in) :: mf

    call write_int('nx', mf%gvec%FFTgrid(1))
    call write_int('ny', mf%gvec%FFTgrid(2))
    call write_int('nz', mf%gvec%FFTgrid(3))

    call write_real('dx', mf%crys%alat * mf%crys%avec(1,1) / mf%gvec%FFTgrid(1))
    call write_real('dy', mf%crys%alat * mf%crys%avec(2,2) / mf%gvec%FFTgrid(2))
    call write_real('dz', mf%crys%alat * mf%crys%avec(3,3) / mf%gvec%FFTgrid(3))

    call write_int('nsp', mf%kp%nspin)

  end subroutine write_sgw_header

  subroutine write_label(label)
    character(len=*), intent(in) :: label

    character(len=9) :: ch
    ch = label
    if (write_bin) then
      write(11) ch
    else
      write(11,'(a)') ch
    endif

  end subroutine write_label

  subroutine write_int(label, val)
    character(len=*), intent(in) :: label
    integer, intent(in) :: val

    character(len=9) :: ch
    ch = label
    if (write_bin) then
      write(11) ch, val
    else
      write(11,'(a,i0)') ch, val
    endif

  end subroutine write_int

  subroutine write_real(label, val)
    character(len=*), intent(in) :: label
    real(DP), intent(in) :: val

    character(len=9) :: ch
    ch = label
    if (write_bin) then
      write(11) ch, val
    else
      write(11,'(a,g0)') ch, val
    endif

  end subroutine write_real

  subroutine write_int_array(vals, siz)
    integer, intent(in) :: siz
    integer, intent(in) :: vals(siz)

    if (write_bin) then
      write(11) vals
    else
      write(11,'(*(i0,:,1x))') vals
    endif

  end subroutine write_int_array

  subroutine write_real_array(vals, siz)
    integer, intent(in) :: siz
    real(DP), intent(in) :: vals(siz)

    if (write_bin) then
      write(11) vals
    else
      write(11,fmtstr) vals
    endif

  end subroutine write_real_array

end program bgw2sgw
