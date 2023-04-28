!=======================================================================================
!
! Module inread_m
!
! (1) inread()          Originally By MLT               Last Modified 7/8/2008 (JRD)
!
!     Read input parameters from file absorption.inp / inteqp.inp.
!
!     input: none
!
!     output: many components of xct
!             nmax :            number of iterations to be performed
!                               (Haydock only!)
!             neig :            number of eigenvectors/eigenvalues
!                               to be computed (diagonalization only!)
!
!=============================================================================================

#include "f_defs.h"

module inread_m

  use global_m
  use scissors_m
  use inread_common_m
  use references_m

  implicit none

  private

  public :: inread

contains

!> supply optionals for absorption, and do not for inteqp
subroutine inread(eqp,xct,flag,nmax,neig)
  type (eqpinfo), intent(out) :: eqp
  type (xctinfo), intent(out) :: xct
  type (flags), intent(out) :: flag
  integer, optional, intent(out) :: nmax, neig

  character*256 :: blockword,keyword,line,errmsg,filename
  integer :: ii,iostat
  logical :: unknown_keyword, found
  
  PUSH_SUB(inread)

  xct%is_absorption = present(nmax) .and. present(neig)
  if(xct%is_absorption .neqv. (present(nmax) .or. present(neig))) then
    call die("inread internal error: all or no optional arguments must be supplied.")
  endif

  if(xct%is_absorption) then
    filename = 'absorption.inp'
  else
    filename = 'inteqp.inp'
  endif

#ifdef MPI
  ! Non-root nodes should wait for root to read the whole file.
  ! That way, we can be sure root gets a chance to write errors before
  ! any die call is issued by another node. Root calls MPI_Barrier below.
  if(peinf%inode /= 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

!--------------------------
! Set default values

  xct%nvb_fi=0
  xct%nvb_co=0
  xct%ncb_fi=0
  xct%inteqp = .not. xct%is_absorption
  xct%ncb_co=0
  xct%nkpt_co = 0
  call scissors_zero(eqp%scis)
  eqp%spl_tck%n=0 
  if(xct%is_absorption) then
    flag%opr = -1 ! no default
  else
    flag%opr = 1  ! momentum for inteqp
  endif
  flag%spec=0
  if(xct%is_absorption) flag%lor=1
  flag%read_dtmat=.false.
  flag%read_dtmat_sub=.false.
  flag%vm=0
  flag%read_epsdiag=.false.
  flag%eig=0
  flag%krnl=1
  flag%bz0=1
  flag%bzq=1
  flag%bzc=1
  flag%lanczos_gauss_quad = .true.
  flag%debug_lanczos = .false.
  if(xct%is_absorption) nmax = 100
  xct%read_kpoints = .false.
  xct%shift(:)=0.d0
  xct%pol(:)=0.d0
  xct%npol=3
  xct%icutv=TRUNC_NONE
  xct%iscreen=SCREEN_SEMICOND
  xct%renorm_transf=.true.
  xct%iwritecoul=0
  xct%skipinterp=.false.
  xct%iabsorp0=0
  xct%truncval(:)=0.0d0
  xct%wplasmon=0.0d0
  if(xct%is_absorption) neig=0
  xct%vmin=1
  xct%vmax=0
  xct%rgrid=0
  xct%qflag=1
  xct%finiteq(:)=0.d0
  xct%no_mtxel=.false.
  xct%avgpot=0d0
  xct%efermi_input=0.0d0
  xct%rfermi=.true.
  xct%avgcut=-1 !FHJ: <0 means "auto", see code below
  xct%scaling=1.0d0
  xct%freplacebz=.false.
  xct%fwritebz=.false.
  xct%delaunay_interp=.true.
  xct%kpt_interp = KPT_INTERP_EXP_TRANSF
  xct%degeneracy_check_override = .false.
  if(xct%is_absorption) then
    xct%averagew=.true.
    xct%eta=0.0d0
    xct%sigma=0.0d0
    xct%gamma=0.0d0
  endif
  xct%eqp_corrections=.false.
  ! inteqp does nothing without eqp_co_corrections
  xct%eqp_co_corrections = .not. xct%is_absorption
  xct%eqp_co_q_corrections = .false.
  xct%algo = BSE_ALGO_DIAG
  xct%extended_kernel=.false.
  xct%unrestricted_transf=.false.
  xct%zero_unrestricted_contrib=.false.
  xct%patched_sampling=.false.
  xct%patched_sampling_co=.false.
  xct%zero_q0_element=0
  xct%theory=0
  xct%coulomb_mod%short_range_frac_fock=1.0d0
  xct%coulomb_mod%long_range_frac_fock=1.0d0
  xct%coulomb_mod%screening_length=0.0d0
  xct%coul_mod_flag=.false.
  xct%npts_intp_kernel = 1
  xct%subsample_cutoff = -10.d0
  xct%subsample_line = .false.
  xct%subsample_algo = 0
  xct%screen_exchange = .false.
  xct%exchange_fact= 1.d0
  xct%direct_fact= 1.d0
  xct%delta_frequency = 1.d-2
#ifdef HDF5
  xct%use_hdf5 = .true.
  xct%use_hdf5_output = .true.
#else
  xct%use_hdf5 = .false.
  xct%use_hdf5_output = .false.
#endif
  xct%use_wfn_hdf5 = .false.
  ! default is to read at least 128 Mb of data, using wfn_hdf5_min_band_block a different
  ! minimum value of band can be set
  xct%wfn_hdf5_min_band_block = -1
  xct%tda = .true.
  xct%zero_coupling_block=.false.
  xct%nspin = -1 ! Determined in input_fi.f90
  xct%nspinor = -1 ! Determined in input_fi.f90
  xct%use_elpa = .false.
#ifdef USEELPA
  xct%use_elpa = .true.
#endif
!#BEGIN_INTERNAL_ONLY
#ifdef USEPRIMME
  xct%primme_tol = 1d-6
  xct%primme_max_block_size = 8
  xct%primme_max_basis_size = 256
#endif
!#END_INTERNAL_ONLY

!-----------------------------
! Never ending loop...

  do while(0.eq.0)

! Actually the loop ends when the end of the file is reached

    read(8,'(a256)',iostat=iostat) line
    if(iostat < 0) exit

! Skip comment lines

    if(len_trim(line).eq.0 .or. line(1:1).eq.'#') cycle

! Determine keyword:
    
    keyword=line(1:scan(line," ")-1)
    line=adjustl(line(scan(line," ")+1:256))
    
    unknown_keyword = .false.
    if(trim(keyword).eq.'begin') then
      blockword=line(1:scan(line," ")-1)
      ii=0
      do while(trim(line).ne.'end')
        read(8,'(a256)',iostat=iostat) line
        if(iostat /= 0) then
          write(errmsg,'(3a)') 'The end of the file was reached while reading elements of the ', &
            trim(blockword),' block. '
          call die(errmsg, only_root_writes = .true.)
        endif
        if(trim(line).ne.'end') then
          ii=ii+1
          if(trim(blockword).ne.'XXXX') then
            write(errmsg,'(3a)') 'Unexpected blockword ', trim(blockword), ' was found in ' // trim(filename) // '.'
            call die(errmsg, only_root_writes = .true.)
          end if
        end if
      end do
    elseif(trim(keyword).eq.'verbosity') then
      read(line,*,err=110) peinf%verbosity
    elseif(trim(keyword).eq.'number_val_bands_fine') then
      read(line,*,err=110) xct%nvb_fi
    elseif(trim(keyword).eq.'number_val_bands_coarse') then
      read(line,*,err=110) xct%nvb_co
    elseif(trim(keyword).eq.'number_cond_bands_fine') then
      read(line,*,err=110) xct%ncb_fi
    elseif(trim(keyword).eq.'number_cond_bands_coarse') then
      read(line,*,err=110) xct%ncb_co

! Spline scissors
    elseif(trim(keyword).eq.'spline_scissors') then
      if (peinf%inode==0) then
        write(6,*) 'Reading spline coefficients from `spline_scissors.dat`'
        write(6,*)
        ! read number of pts, knots, coefficients, degree
        call open_file(20,file='spline_scissors.dat',form='formatted',status='old')
        read(20,*) eqp%spl_tck%n
        SAFE_ALLOCATE(eqp%spl_tck%t, (eqp%spl_tck%n))
        SAFE_ALLOCATE(eqp%spl_tck%c, (eqp%spl_tck%n))
        read(20,*) (eqp%spl_tck%t(ii), ii=1,eqp%spl_tck%n)
        read(20,*) (eqp%spl_tck%c(ii), ii=1,eqp%spl_tck%n)
        read(20,*) eqp%spl_tck%k
        close(20)
      endif

    elseif(trim(keyword).eq.'coarse_grid_points') then
      read(line,*,err=110) xct%nkpt_co
    elseif(trim(keyword).eq.'read_kpoints') then
      xct%read_kpoints = .true.
    elseif(trim(keyword).eq.'avgpot') then
      read(line,*,err=110) xct%avgpot
    elseif(trim(keyword).eq.'fermi_level') then
      read(line,*,err=110) xct%efermi_input
      if (peinf%inode==0 .and..not. xct%is_absorption) then
        write(0,*) 'WARNING: Interpolating eqp_co.dat using a shifted Fermi level.'
        write(0,*) '         The resulting eqp.dat may be unsuitable for absorption.x'
      endif
    elseif(trim(keyword).eq.'fermi_level_absolute') then
      xct%rfermi=.false.
    elseif(trim(keyword).eq.'fermi_level_relative') then
      xct%rfermi=.true.
    elseif(trim(keyword).eq.'no_symmetries_fine_grid') then
      flag%bz0 = 1
    elseif(trim(keyword).eq.'use_symmetries_fine_grid') then
      flag%bz0 = 0
    elseif(trim(keyword).eq.'no_symmetries_shifted_grid') then
      flag%bzq = 1
    elseif(trim(keyword).eq.'use_symmetries_shifted_grid') then
      flag%bzq = 0
    elseif(trim(keyword).eq.'no_symmetries_coarse_grid') then
      flag%bzc = 1
    elseif(trim(keyword).eq.'use_symmetries_coarse_grid') then
      flag%bzc = 0
    elseif(trim(keyword).eq.'lowest_occupied_band') then
      read(line,*,err=110) xct%vmin
    elseif(trim(keyword).eq.'highest_occupied_band') then
      read(line,*,err=110) xct%vmax
    elseif(trim(keyword).eq.'regular_grid') then
      read(line,*,err=110) xct%rgrid
    elseif(trim(keyword).eq.'fullbz_replace') then
      xct%freplacebz=.true.
    elseif(trim(keyword).eq.'fullbz_write') then
      xct%fwritebz=.true.
    elseif(trim(keyword).eq.'use_velocity') then
      xct%npol=1
      flag%opr=0
    elseif(trim(keyword).eq.'q_shift') then
      read(line,*,err=110) (xct%shift(ii),ii=1,3)
    elseif(trim(keyword).eq.'use_momentum') then
      flag%opr=1
    elseif(trim(keyword).eq.'degeneracy_check_override') then
      xct%degeneracy_check_override=.true.
    elseif(trim(keyword).eq.'greedy_interpolation') then
      xct%delaunay_interp=.false.
    elseif(trim(keyword).eq.'delaunay_interpolation') then
      xct%delaunay_interp=.true.
    elseif(trim(keyword).eq.'kpt_interpolation_linear') then
      xct%kpt_interp=KPT_INTERP_LINEAR
    elseif(trim(keyword).eq.'kpt_interpolation_exp_transform') then
      xct%kpt_interp=KPT_INTERP_EXP_TRANSF
    elseif(trim(keyword).eq.'dont_use_hdf5') then
      xct%use_hdf5 = .false.
    elseif(trim(keyword).eq.'dont_use_hdf5_output') then
      xct%use_hdf5_output = .false.
    elseif(trim(keyword).eq.'use_wfn_hdf5') then
#ifdef HDF5
      xct%use_wfn_hdf5 = .true.
#endif
    elseif(trim(keyword).eq.'wfn_hdf5_min_band_block') then
      read(line,*,err=110)  xct%wfn_hdf5_min_band_block
    elseif(trim(keyword).eq.'extended_kernel') then
      xct%extended_kernel=.true.
    elseif(trim(keyword).eq.'unrestricted_transformation') then
      xct%unrestricted_transf=.true.
    elseif(trim(keyword).eq.'zero_unrestricted_contribution') then  
      xct%zero_unrestricted_contrib=.true.
    elseif(trim(keyword).eq.'delta_frequency') then
      read(line,*,err=110) xct%delta_frequency
    elseif(trim(keyword).eq.'zero_q0_element') then
      read(line,*,err=110) xct%zero_q0_element
    elseif(trim(keyword).eq.'subsample_line') then
      xct%subsample_line  = .true.
      read(line,*,err=110) xct%subsample_cutoff
      if(peinf%inode.eq.0) write(6,*) "Will read subsampled BSE matrix elements from file"
    elseif(trim(keyword).eq.'subsample_algo') then
      read(line,*,err=110) xct%subsample_algo
!#BEGIN_INTERNAL_ONLY
    ! FHJ: in metallic systems, one should not use the TDA. In the meantime,
    ! we can turn off the renormalization of the interpolation transformation
    ! coefficients for debugging purposes.
    elseif(trim(keyword).eq.'no_renormalization') then
      xct%renorm_transf=.false.
    elseif(trim(keyword).eq.'patched_sampling') then  
      xct%patched_sampling=.true.
    elseif(trim(keyword).eq.'patched_sampling_co') then
      xct%patched_sampling_co=.true.
    elseif(trim(keyword).eq.'bse') then
      xct%theory=0
    elseif(trim(keyword).eq.'tddft') then
      xct%theory=1
    elseif(trim(keyword).eq.'short_range_frac_fock') then
      read(line,*,err=110) xct%coulomb_mod%short_range_frac_fock
    elseif(trim(keyword).eq.'long_range_frac_fock') then
      read(line,*,err=110) xct%coulomb_mod%long_range_frac_fock
    elseif(trim(keyword).eq.'screening_length') then
      read(line,*,err=110) xct%coulomb_mod%screening_length
    elseif(trim(keyword).eq.'q0vec') then
        read(line,*,err=110) (xct%q0vec(ii),ii=1,3)
    elseif(trim(keyword).eq.'dont_check_norms') then
      peinf%check_norms=.false.
!#END_INTERNAL_ONLY
    else
      call scissors_inread(keyword, line, eqp%scis, found)
      unknown_keyword = .not. found
    endif

    if(unknown_keyword .and. xct%is_absorption) then
      unknown_keyword = .false.
      if(trim(keyword).eq.'eqp_corrections') then
        xct%eqp_corrections=.true.
      elseif(trim(keyword).eq.'eqp_co_corrections') then
        xct%eqp_co_corrections=.true.
      elseif(trim(keyword).eq.'eqp_co_q_corrections') then
        xct%eqp_co_q_corrections=.true.
      elseif(trim(keyword).eq.'energy_resolution') then
        read(line,*,err=110) xct%eta
      elseif(trim(keyword).eq.'energy_resolution_sigma') then
        read(line,*,err=110) xct%sigma
      elseif(trim(keyword).eq.'energy_resolution_gamma') then
        read(line,*,err=110) xct%gamma
      elseif(trim(keyword).eq.'use_dos') then
        flag%opr=2
      elseif(trim(keyword).eq.'tda_bse') then
        xct%tda=.true.
      elseif(trim(keyword).eq.'full_bse') then
        xct%tda=.false.
      elseif(trim(keyword).eq.'exciton_Q_shift') then
        xct%qflag=2
        read(line,*,err=110) xct%qflag,xct%finiteq(1),xct%finiteq(2),xct%finiteq(3)
!#BEGIN_INTERNAL_ONLY
      elseif(trim(keyword).eq.'zero_coupling_block') then
        xct%zero_coupling_block=.true.
      elseif(trim(keyword).eq.'dont_use_matrix_elements') then
        xct%no_mtxel=.true.
      elseif(trim(keyword).eq.'exchange_factor') then
        read(line,*,err=110) xct%exchange_fact
      elseif(trim(keyword).eq.'direct_factor') then
        read(line,*,err=110) xct%direct_fact
!#END_INTERNAL_ONLY
      elseif(trim(keyword).eq.'read_eigenvalues') then
        flag%spec=1
      elseif(trim(keyword).eq.'read_dtmat') then
        flag%read_dtmat=.true.
      elseif(trim(keyword).eq.'number_iterations') then
        read(line,*,err=110) nmax
      elseif(trim(keyword).eq.'number_eigenvalues') then
        read(line,*,err=110) neig
      elseif(trim(keyword).eq.'cell_average_cutoff') then
        read(line,*,err=110) xct%avgcut
      elseif(trim(keyword).eq.'lorentzian_broadening') then
        flag%lor=0
      elseif(trim(keyword).eq.'gaussian_broadening') then
        flag%lor=1
      elseif(trim(keyword).eq.'voigt_broadening') then
        flag%lor=2
      elseif(trim(keyword).eq.'polarization') then
        xct%npol=1
        read(line,*,err=110) (xct%pol(ii),ii=1,3)
      elseif(trim(keyword).eq.'read_vmtxel') then
        flag%vm=1
      elseif(trim(keyword).eq.'read_eps2_moments') then
        flag%vm=2
      elseif(trim(keyword).eq.'read_epsdiag') then
        flag%read_epsdiag=.true.
      elseif(trim(keyword).eq.'write_eigenvectors') then
        flag%eig=-1
        read(line,*) flag%eig
      elseif(trim(keyword).eq.'diagonalization') then
        xct%algo = BSE_ALGO_DIAG
!#BEGIN_INTERNAL_ONLY
#ifdef USEPRIMME
      elseif(trim(keyword).eq.'diagonalization_primme') then
        xct%algo = BSE_ALGO_DIAG_PRIMME
#endif
!#END_INTERNAL_ONLY
      elseif(trim(keyword).eq.'lanczos') then
        xct%algo = BSE_ALGO_LANCZOS
      elseif(trim(keyword).eq.'haydock') then
        xct%algo = BSE_ALGO_HAYDOCK
!#BEGIN_INTERNAL_ONLY
#ifdef USEPRIMME
      elseif(trim(keyword).eq.'primme_tolerance') then
        read(line,*) xct%primme_tol
      elseif(trim(keyword).eq.'primme_max_basis_size') then
        read(line,*) xct%primme_max_basis_size
      elseif(trim(keyword).eq.'primme_max_block_size') then
        read(line,*) xct%primme_max_block_size
#endif
!#END_INTERNAL_ONLY
      elseif(trim(keyword).eq.'spin_triplet') then
        flag%krnl=0
      elseif(trim(keyword).eq.'spin_singlet') then
        flag%krnl=1
      elseif(trim(keyword).eq.'local_fields') then
        flag%krnl=2
      elseif(trim(keyword).eq.'spinor') then
        flag%krnl=3
      elseif(trim(keyword).eq.'dont_use_elpa') then
        xct%use_elpa = .false.
      elseif(trim(keyword).eq.'use_elpa') then
        xct%use_elpa = .true.
      elseif(trim(keyword).eq.'no_lanczos_gauss_quadrature') then
        flag%lanczos_gauss_quad = .false.
      elseif(trim(keyword).eq.'lanczos_gauss_quadrature') then
        flag%lanczos_gauss_quad = .true.
!#BEGIN_INTERNAL_ONLY
      elseif(trim(keyword).eq.'debug_lanczos') then
        flag%debug_lanczos = .true.
!#END_INTERNAL_ONLY
      elseif(trim(keyword).eq.'noeh_only') then
        xct%iabsorp0=1
      elseif(trim(keyword).eq.'kernel_k_interpolation') then
        xct%npts_intp_kernel = -1
      elseif(trim(keyword).eq.'write_vcoul') then
        xct%iwritecoul=1
      elseif(trim(keyword).eq.'skip_interpolation') then
        xct%skipinterp=.true.
      elseif(trim(keyword).eq.'kernel_scaling') then
        read(line,*,err=110) xct%scaling
      elseif(trim(keyword).eq.'average_w') then
        xct%averagew=.true.
      elseif(try_inread_truncation(trim(keyword), trim(line), xct%icutv, xct%truncval(1))) then
        ! subroutine already does the job
      elseif(try_inread_screening(trim(keyword), trim(line), xct%iscreen)) then
        ! subroutine already does the job
      else
        unknown_keyword = .true.
      endif
    endif ! xct%is_absorption

    if(unknown_keyword) then
      write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found in ' // trim(filename) // '.'
      call die(errmsg, only_root_writes = .true.)
    end if
  enddo
  
  xct%krnl = flag%krnl

  call peinfo_set_verbosity()
  if(xct%nvb_fi.eq.0) then
    call die('The number_val_bands_fine keyword could not be found.', only_root_writes = .true.)
  endif
  
  if(xct%ncb_fi.eq.0) then
    call die('The number_cond_bands_fine keyword could not be found.', only_root_writes = .true.)
  endif
  
  if(xct%nvb_co.eq.0) then
    call die('The number_val_bands_coarse keyword could not be found.', only_root_writes = .true.)
  endif
  
  if(xct%ncb_co.eq.0) then
    call die('The number_cond_bands_coarse keyword could not be found.', only_root_writes = .true.)
  endif

  if (xct%is_absorption) then
    if (xct%algo==BSE_ALGO_HAYDOCK) then
      ! FHJ: Haydock will only work with 1 polarization for now, even for momentum.
      xct%npol = 1
    endif
    if (flag%opr==1 .and. xct%npol==1 .and. all(xct%pol(:)==0d0)) then
      call die('Polarization direction not set.', only_root_writes=.true.)
    endif
    if (xct%algo/=BSE_ALGO_DIAG .and. xct%algo/=BSE_ALGO_DIAG_PRIMME .and. nmax<1) then
      call die('Invalid number of iterations.', only_root_writes=.true.)
    endif
  else
    ! FHJ: inteqp doesn`t care about npol
    xct%npol = 1
  endif

!#BEGIN_INTERNAL_ONLY
  if (xct%algo==BSE_ALGO_DIAG_PRIMME) then
#ifndef MPI
    call die('PRIMME iterative solver only implemented for MPI builds.')
#endif
    if (neig==0) call die('Number of eigenvalues must be set for iterative diagonalization.')
    if (.not.xct%tda) call die('Iterative diagonalization not implemented with full BSE solver.')
  endif
!#END_INTERNAL_ONLY

  if (xct%algo/=BSE_ALGO_DIAG .and. &
      xct%algo/=BSE_ALGO_DIAG_PRIMME .and. &
      xct%algo/=BSE_ALGO_LANCZOS .and. &
      xct%algo/=BSE_ALGO_HAYDOCK) then
      call die('Invalid solver algorithm.', only_root_writes=.true.)
  endif

#ifndef USEELPA
  if (xct%is_absorption) then
    xct%use_elpa = .true.
  endif
#endif

#ifndef USESCALAPACK
  if (xct%algo==BSE_ALGO_LANCZOS) then
    call die('Can only use lanczos algo. when code is compiled with ScaLAPACK support.',&
      only_root_writes=.true.)
  endif
#endif

  if (xct%subsample_line) then
    if(.not.xct%use_hdf5) then
      call die('Subsampling not implemented for non-hdf5 file format',only_root_writes=.true.)
    elseif ( xct%npts_intp_kernel.eq.-1 ) then
      call die('Subsampling not implemented for kernel_k_interpolation',only_root_writes=.true.)
    endif
  endif

  if (xct%subsample_algo/=0 .and. &
      xct%subsample_algo/=1) then
      call die('Invalid algorithm to find vertices for subsampling.', only_root_writes=.true.)
  endif


! FHJ: Default values for full BSE calculation
  if (.not.xct%tda) then
    xct%extended_kernel = .true.
    xct%unrestricted_transf = .true.
    !flag%lor = 0
    if (peinf%inode==0) then
      write(0,'(/,a,/)') 'WARNING: full BSE calculation is under development.'
    endif
  endif
! FHJ: How many kernel blocks did we compute? (inteqp doesn`t care about this)
  if (xct%is_absorption) then
    if (xct%extended_kernel) then
      xct%unrestricted_transf = .true.
      if (peinf%inode==0) then
        write(6,'(1x,a)') "Kernel was calculated for all possible (n1,n2)->(n1',n2') transitions."
        write(6,'(1x,a)') "This introduces interpolation artifacts for semiconducting systems."
        if (xct%zero_unrestricted_contrib) then
          write(6,'(1x,a)') "However, only the restricted portion of the BSE kernel will be used."
          !write(6,'(1x,a)') "Use this option for testing purposes only."
        endif
        write(6,*)
      endif
    else
      if (peinf%inode==0) then
        write(6,'(1x,a)') "Kernel was calculated only for (v,c)->(v',c') transitions."
        write(6,'(1x,a)') "This may not be a good idea for metals!"
        write(6,*)
      endif
    endif
  else
    ! FHJ: for inteqp, there`s no point for extended_kernel and unrestricted_transf
    ! to be different, that will only cause extra buffers to be allocated.
    xct%extended_kernel = xct%unrestricted_transf
  endif

  ! FHJ: Note that intwfn will redefine xct%n1b_co temporarily, but below
  ! are the final for values xct%n1b_co that we will use in intkernel.
  if (xct%extended_kernel) then
    xct%n1b_co = xct%nvb_co + xct%ncb_co
    xct%n2b_co = xct%n1b_co
  else
    xct%n1b_co = xct%nvb_co
    xct%n2b_co = xct%ncb_co
  endif

  if (peinf%inode==0) then
!#BEGIN_INTERNAL_ONLY
    if (xct%is_absorption) then
      if (xct%theory.eq.0) then
        write(6,'(1x,a)') 'This is a Bethe-Salpeter-equation calculation.'
      elseif (xct%theory.eq.1) then
        write(6,'(1x,a)') 'This is a time-dependent density-functional-theory calculation.'
      endif
    endif
    ! such a distinction does not apply to inteqp
!#END_INTERNAL_ONLY
    if (xct%unrestricted_transf) then
      ! FHJ: We actually do this unrestricted transformation even if there`s no
      ! eqp to be interpolated. This might sound silly if we are not using an
      ! extended kernel, because we`ll have to restrict the coefficients later
      ! on. However, by comparing the error between the unrestricted/restricted
      ! interpolations, we can determine the validity of the restricted kernel.
      write(6,'(1x,a)') "The unrestricted fi/co WFN transformation dvn',dcn' will be calculated."
      if (.not.xct%extended_kernel) then
        write(6,'(1x,a)') "However, the coefficients will be later projected onto the dvv',dcc' subspaces."
      endif
    else
      write(6,'(1x,a)') "The restricted fi/co WFN transformation dvv',dcc' will be calculated."
      !write(6,'(1x,a)') "This is not a good idea for metals and narrow-bandgap semiconductors."
      !write(6,'(1x,a)') "Consider using the 'unrestricted_transformation' flag."
    endif
    write(6,*)
    if (xct%tda) then
      write(6,'(1x,a,/)') 'Using the Tamm-Dancoff approximation.'
    else
      write(6,'(1x,a,/)') 'Solving the full BSE without the Tamm-Dancoff approximation.'
    endif
  endif

  if(flag%opr.eq.0) then
    if(peinf%inode.eq.0) then
      write(6,*) 'Using the velocity operator'
      write(6,*)
    endif
    if (xct%skipinterp) then
      call die('skip_interpolation is incompatible with velocity operator', only_root_writes = .true.)
    endif
  elseif(flag%opr.eq.1) then
    if(peinf%inode.eq.0) then
      write(6,*) 'Using the momentum operator'
      write(6,*)

      if (flag%bz0 /= flag%bzq) then
        write(0,*)
        write(0,*) 'ERROR: When using the momentum operator, you must use the same'
        write(0,*) 'symmetry operators for the fine ({use,no}_symmetries_fine_grid)'
        write(0,*) 'and q-shifted ({use,no}_symmetries_shifted_grid) grids'
        write(0,*)
        call die('Inconsistent symmetry treatment of the fine and shifted grids with the momentum operator')
      endif

      ! for inteqp, this just means we are only interpolating to WFN_fi, which is not
      ! an approximation, so we do not need to write a warning
      if(xct%is_absorption) then
        write(0,999)
        write(0,*)
999     format(1x,'**************************************************',/, &
          1x,'**                                              **',/, &
          1x,'**                   WARNING:                   **',/, &
          1x,'**      THE MOMENTUM OPERATOR IS NOT EXACT      **',/, &
          1x,'**  ONLY USE IT IF YOU KNOW WHAT YOU ARE DOING  **',/, &
          1x,'**  OTHERWISE YOU SHOULD USE VELOCITY OPERATOR  **',/, &
          1x,'**    for details see equations (29) and (30)   **',/, &
          1x,'**    of Rohlfing & Louie PRB 62, 4927 (2000)   **',/, &
          1x,'**                                              **',/, &
          1x,'**************************************************')
      endif
    endif
  elseif(flag%opr.eq.2) then
    if(peinf%inode.eq.0) then
      write(6,*) 'Using the Joint density of states operator'
      write(6,*)
    endif
  else
    call die("No flag for operator", only_root_writes = .true.)
  endif

  if (peinf%inode==0 .and. .not. xct%skipinterp) then
    if (xct%delaunay_interp) then
      write(6,'(1x,a,/)') 'Using the Delaunay interpolation method.'
    else
      write(6,'(1x,a)') 'Using the greedy interpolation method [DEPRECATED].'
      write(0,'(1x,a)') 'WARNING: This interpolation scheme is deprecated, and might lead to'
      write(0,'(1x,a)') 'non-continuous interpolants. Please make sure the new option'
      write(0,'(1x,a,/)') '`delaunay_interpolation` works in your system, and make the migration.'
    endif
  endif

  if (xct%nkpt_co>0 .and. peinf%inode==0) then
    write(0,*) 'WARNING: `coarse_grid_points` keyword is deprecated'
  endif

  ! FHJ: Let`s not pretend that QP corrections with shifted Fermi energy works,
  ! it doesn`t work with inteqp nor absorption. However, you can in principle run
  ! inteqp without any shift in the Fermi energy, which should give consistent results.
  ! See ticket #197.
  if ((abs(xct%efermi_input) > TOL_Zero .or. .not. xct%rfermi) .and. &
    (xct%eqp_corrections .or. xct%eqp_co_corrections)) then
    if (peinf%inode==0) then
      write(0,'(/a)') "WARNING: Use of shifted Fermi Energy and eqp*corrections will"
      write(0,'(a)') "likely give meaningless results. Either use scissors operators,"
      write(0,'(a)') "or manually fix the mean-field energies in the WFN files."
    endif
  endif

  if(xct%is_absorption) then ! rest of file does not apply to inteqp
    ! no fundamental reason, but the way the code is written, energy interpolation only occurs if kernel is interpolated
    if(xct%eqp_co_corrections .and. xct%skipinterp) then
      call die("Cannot do eqp_co_corrections with skip_interpolation. Use eqp_corrections instead.", &
        only_root_writes = .true.)
    endif

! FHJ: noeh_only doesn`t interpolate onto fine grid, so eqp_co won`t work
    if(xct%eqp_co_corrections .and. (.not.xct%eqp_corrections) .and. (xct%iabsorp0==1) .and. (flag%spec==0)) then
      if (peinf%inode==0) then
        write(0,*) 'Can`t use `eqp_co_corrections` with `noeh_only`. Either:'
        write(0,*) ' (1) Run inteqp, then run absorption with eqp_corrections, or'
        write(0,*) ' (2) Remove the `noeh_only` flag and run absorption directly.'
      endif 
      call die('Can`t use `eqp_co_corrections` and `noeh_only` simultaneously.',  only_root_writes = .true.)
    endif

    if ((abs(xct%efermi_input) > TOL_Zero .or. .not. xct%rfermi)) then ! Fermi level is being adjusted
      if((xct%eqp_corrections .neqv. xct%eqp_co_corrections) .and. .not. xct%skipinterp) then
        if(peinf%inode == 0) then
          write(0,'(a)') "If Fermi level is adjusted and interpolation is used, you must use both or neither"
          write(0,'(a)') "of eqp_corrections and eqp_co_corrections. (You can get eqp.dat by running inteqp.)"
        endif
        call die("Cannot use this combination of Fermi adjustment / eqp_corrections / interpolation.", &
          only_root_writes = .true.)
      endif
    endif
    ! eqp_corrections only cannot be used because we would need to 'anti-interpolate' to the coarse grid
    ! in order to correctly reset occupations (ifmax) on the coarse grid.
    ! eqp_co_corrections only cannot be used because then the Fermi level is set by the coarse grid, and
    ! it is cumbersome to reset the fine grid Fermi level since WFN_fi was already handled. --DAS
    ! But if we are not interpolating, then it is ok, since eqp.dat will be used on both.

    if(flag%read_dtmat .and. xct%skipinterp) then
      call die("Options read_dtmat and skip_interpolation are incompatible.", only_root_writes = .true.)
    endif

    if(flag%read_dtmat .and. xct%eqp_co_corrections) then
      if(peinf%inode == 0) then
        write(0,'(a)') "read_dtmat flag causes the eqp_co_corrections flag to be ignored, even if the"
        write(0,'(a)') "eqp_co_corrections flag has been set to .true. (leading to no eqp_co_corrections)."
        write(0,'(a)') "If you need eqp_co_corrections = .true., please set read_dtmat = .false."
      endif
      call die("Options read_dtmat and eqp_co_corrections cannot be both .true. .", only_root_writes = .true.)
    endif

    ! momentum and JDOS do not have q-shift. specifying one mistakenly will break things without this.
    if(flag%opr /= 0) then
      xct%shift(1:3) = ZERO
    endif
  
    if(abs(xct%eta).lt.TOL_Zero) then
      call die('The energy_resolution keyword could not be found.', only_root_writes = .true.)
    endif
  
    if(flag%lor.eq.2.and.abs(xct%sigma).lt.TOL_Small) then
      call die('Voigt requires energy_resolution_sigma > 0 and energy_resolution_gamma >= 0', &
        only_root_writes = .true.)
    endif

! JRD: Make a note if we have finite Q
    if (peinf%inode==0) then
      if (xct%qflag==0 .or. xct%qflag==2) then
        write(6,'(1x,a)') 'We are doing a finite-Q calculation'
        write(6,'(1x,a,3(1x,f12.9))') 'Q_shift =', xct%finiteq(:)
        write(6,'(1x,a)') 'where Q_shift is the NEGATIVE of the exciton center-of-mass (COM) momentum'
      else
        write(6,'(1x,a)') 'We are doing a calculation with zero exciton momentum'
      endif
    endif

! JRD: Check compatibility

    if ((xct%qflag==0 .and. flag%opr==1) .and..not. xct%no_mtxel) then
      call die('Momentum and Finite_q are incompatible',  only_root_writes = .true.)
    end if
    
    if (xct%qflag .eq. 0 .and. xct%skipinterp) then
      call die("Use of WFNq_co in finite momentum not compatible with skip_interpolation.",&
      only_root_writes=.true.)
    endif

    if (xct%qflag .eq. 0 .and. (.not. xct%delaunay_interp)) then
      call die("Use of WFNq_co in finite momentum only implemented with delaunay interpolation",&
      only_root_writes=.true.)
    endif

    if (xct%qflag .eq. 0 .and. flag%read_dtmat) then
      call die("read_dtmat not yet implemented for finite momentum with WFNq_co.",&
      only_root_writes=.true.)
    endif

    if (xct%qflag .eq. 0 .and. (.not.xct%eqp_co_q_corrections) .and. (xct%eqp_co_corrections)) then
      call die("Must use eqp_co_q_corrections with finite momentum with WFNq_co.",&
      only_root_writes=.true.)
    endif
    
    ! JRD: What screening is present?
    if (peinf%inode==0) then
      select case (xct%iscreen)
        case (SCREEN_SEMICOND)
          write(6,'(1x,a/)') 'Running with semiconductor screening'
        case (SCREEN_GRAPHENE)
          write(6,'(1x,a/)') 'Running with graphene screening'
        case (SCREEN_METAL)
          write(6,'(1x,a/)') 'Running with metal screening'
        case default
          call die('Unknown screening type', only_root_writes=.true.)
      endselect
    endif

    if(peinf%inode == 0) then
      if(peinf%npes > 1) then
        write(6,803)
      else
        write(6,805)
      endif
    endif
803 format(1x,'We are communicating via MPI',/)
805 format(1x,'We are not communicating',/)

    ! Make a note if we have finite Q
    ! WHERE (?)

    call print_truncation_summary(xct%icutv, xct%truncval(1))

    ! FHJ: set cell_average_cutoff, if not set by user
    if (xct%avgcut<0) then
      xct%avgcut = TOL_ZERO
      if (xct%icutv==TRUNC_NONE .and. xct%iscreen==SCREEN_SEMICOND) then
        xct%avgcut = 1d0/TOL_ZERO
      endif
    endif
    if (peinf%inode==0) then
      write(6,'(1x,a,es10.3e3,a/)') &
        'Cutoff for Monte-Carlo average of Coulomb potential: ', xct%avgcut, ' Ry'
    endif

  endif ! xct%is_absorption 

  call require_reference(REF_Deslippe2012)
  call require_reference(REF_Rohlfing2000)
  ! Truncation of the Coulomb potential: slab and write
  if (xct%icutv==TRUNC_SLAB .or. xct%icutv==TRUNC_WIRE) call require_reference(REF_IsmailBeigi2016)
  ! Finite-Q excitons
  if (xct%qflag/=1) call require_reference(REF_Qiu2015)
  ! Non-TDA
  if (xct%extended_kernel .or. .not.xct%tda) call require_reference(REF_Shao2016)
  ! New Lanczos algorithm
  if (xct%algo==BSE_ALGO_LANCZOS) call require_reference(REF_Shao2018)

#ifdef MPI
  ! root lets the others go after it is done reading (see beginning of function)
  if(peinf%inode == 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  ! FHJ: Broadcast spline scissors (this must be called after the barrier)
  call MPI_BCAST(eqp%spl_tck%n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  if ( eqp%spl_tck%n > 0 ) then
    if(peinf%inode/=0) then
      SAFE_ALLOCATE(eqp%spl_tck%t, (eqp%spl_tck%n))
      SAFE_ALLOCATE(eqp%spl_tck%c, (eqp%spl_tck%n))
    endif
    call MPI_BCAST(eqp%spl_tck%t, eqp%spl_tck%n, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(eqp%spl_tck%c, eqp%spl_tck%n, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(eqp%spl_tck%k, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  endif
#endif

  POP_SUB(inread)
  
  return
  
110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', &
      trim(keyword), '. '
  call die(errmsg, only_root_writes = .true.)
  
end subroutine inread

end module inread_m
