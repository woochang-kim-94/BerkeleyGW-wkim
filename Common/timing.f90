!==============================================================================
!
! Routines:
!
! (1) date_time()       Originally by ?         Last Modified: 5/12/2008 (JRD)
!
!     Gets current date and time.
!
! (2) timget()          Originally by gsm       Last Modified: 4/29/2010 (gsm)
!
!     Gets current cpu and wall time.
!     Note: it`s almost a private subroutine, if not for io_utils.f90
!
! (3) timacc(n,option,tsec,nslices)   Originally by ?
!                                               Last Modified: 6/17/2009 (PWD)
!     DEPRECATED
!
!     Timing subroutine.  Calls machine-dependent subroutine timget
!     which returns elapsed cpu and wall clock times in seconds
!     Also return the number of times the counter has been called
!
!     Depending on value of "option" routine will:
!       (0) zero all accumulators
!       (1) start with new incremental time slice for accumulator n
!           also increase by one the counter for this accumulator
!       (2) stop time slice; add time to accumlator n
!       (3) report accumulated time for accumulator n
!           and number of time that the routine has been called
!       (4) report time slice for accumulator n (not full time accumulated)
!
!       If, on first entry, subroutine is not being initialized, it
!       will automatically initialize as well as rezero accumulator n.
!       However, initialization SHOULD be done explicitly by the user
!       so that it can be done near the top of his/her main routine.
!
!       Input:
!         n=index of accumulator (distinguish what is being timed); not used if
!         option=0 option=see comment above
!       Output:
!         on option=3:
!         tottim(2,n)=accumulated time for accumulator n; otherwise
!         tottim is a dummy variable.
!         nslices is optional variable that give number of slices collected
!
! (4) logit()    Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!
!    Write out a debugging message with an inputed string and write time.
!
! (5) logitint() Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!
!    Same as logit but with an integer constant.
!
! (6) Timing class. See at the type definition for its (verbose) description.
!
!==============================================================================
!
! Todo: move logit and logitint to another place?

#include "f_defs.h"

module timing_m

  use os_m
  use message_m
  use nrtype_m
  use peinfo_m
  use push_pop_m

  implicit none

  private

  public :: date_time, timget, timacc, logit, logitint

  !> MTIM determines the maximum number of "timing slots" available
  integer, parameter, private :: MTIM=100
  real(DP), private, save :: acctim(2,MTIM), tzero(2,MTIM)
  integer, private, save :: ncount(MTIM)
  !
  !----------------------------------------------------------------------------
  !
  ! Timing class
  !
  ! Reference object to handle the timing of subroutines and of the program.
  ! The object is implemented as following.
  ! type(timing_class) is a baseclass that should not be used directly, as it
  ! only contains the definition of class methods.
  ! Below, the base class is extended for each code
  !
  ! First, the type timing_class is an abstract implementation of the
  ! timing class. It contains the common methods that are used to
  ! time the duration.
  ! The subclasses are the types that should be used in the code, that is:
  ! - timing_epsilon_class: timing of subroutines in ./Epsilon
  ! - timing_sigma_class: timing of subroutines in ./Sigma
  ! - timing_bse_class: timing of subroutines in ./BSE
  ! - timing_common_class: timing of subroutines in ./Common
  ! - timing_extra_class: the rest (PlotXct, ...)
  !
  ! The module defines 5 objects, that should be loaded when needed:
  ! - timing_epsilon, timing_sigma, timing_bse, timing_common, timing_extra
  !
  ! The schematic usage is as follow (using epsilon as example)
  !
  ! program epsilon
  !   use timing_m, only: timing => timing_epsilon

  !   call timing%init()
  !
  !   call timing%start(timing%sub1)
  !   ...
  !   call timing%stop(timing%sub1)
  !
  !   call timing%print()
  ! end program epsilon
  !
  ! where:
  ! - we renamed timing_epsilon as timing, so that it looks as in the code we
  !   always use the same object
  ! - we initialized the timing object (zeroing reference timings)
  ! - timing%sub1 is an attribute of the timing class
  ! - we kept track of the elapsed time between start() and stop().
  !   One can make several calls that refer to the same attribute (e.g. %sub1).
  !   The total time associated to %sub1 will be the sum of all the timings.
  ! - we printed the timing information to screen.
  !
  ! Note: print() should generally be at the end of the program
  !
  ! Note: if you call %start twice without calling %stop before, for example:
  !   call timing%start(timing%sub1)
  !   ...
  !   call timing%start(timing%sub1)
  !   ...
  !   call timing%stop(timing%sub1)
  !   Then, the timing will keep track of the time from the 2nd start to stop.
  !
  ! Note: if you call stop() without preceding it by a start(), timings will be
  !   almost random numbers. Note also that safe-proofing this case requires
  !   introducing a few `if` instructions that would slow the code.
  !   At the moment, the choice is to not correct this.
  !
  ! Modifications / New tags:
  !
  ! It`s recommended to not add new subclasses, rather, one should try to
  ! merge everything, if possible.
  ! This is an example to add a new tag (%sub2) to timing_sigma
  ! (or timing_epsilon, ...). The procedure is:
  !
  ! 1) Locate the definition of timing_sigma_class, e.g. the line with
  !    type, extends(timing_class) :: timing_sigma_class
  !
  !    Add a new variable 'sub2' and assign an integer value to it between
  !    1 and 100. MAKE SURE that the integer value hasn`t been taken already,
  !    or that the integer namespace hasn`t been taken already.
  !    As an example, suppose we can writeL
  !    integer :: sub2 = 9
  !
  ! 2) Locate the subroutine init_sigma_labels (or init_epsilon_labels, ...).
  !    Add a label to sub2:
  !    instance%labels(34) = 'Subroutine 2'
  !    This string will be used only for a nice formatting of the printing info
  !
  ! Devel note: one could have implemented a sort of dictionary in fortran.
  ! Then, one could have simply used a syntax "call timing%start('sub2')"
  ! However, it may lead to poor performance, since parts of the code call
  ! the timing information several thousand times
  !
  ! Devel note: due to a limitation of Fortran, it`s more convenient to
  ! initialize the chartacter array of labels at runtime, rather than
  ! writing them in the class description
  !
  type :: timing_class
     !
     ! Picky note: wall time is probably incorrect when there is a change in
     ! the OS time: for example, when daylight saving time changes or if the
     ! computer is moved across time zones.
     !
     integer :: num_times = MTIM ! max number of timing attributes
     character(len=100) :: labels(MTIM) ! array of timing labels, for printing
     !
     ! arrays used to make timing measurements.
     real(DP), dimension(MTIM) :: wall_times, cpu_times, &
          tmp_cpu_times, tmp_wall_times
     integer :: call_numbers(MTIM) ! the number of start() calls for each
     ! attribute

   contains
     ! Initialize arrays and timings
     procedure :: init => timing_class_init
     ! Print to screen a summary of timings
     procedure :: print => timing_class_print
     ! Start the chronometer for a tag
     procedure :: start => timing_class_start
     ! Stop the chronometer for a tag
     procedure :: stop => timing_class_stop
     ! Initialize the object
     procedure :: init_labels => bare_init_labels
  end type timing_class
  !
  type, extends(timing_class) :: timing_epsilon_class
     ! class definition for timing of subrourines in ./Epsilon
     integer :: io_total = 1
     integer :: input = 2
     integer :: input_q = 3
     integer :: fullbz = 4
     integer :: gvec = 5
     integer :: subgrp = 6
     integer :: irrbz = 8
     integer :: genwf = 9
     integer :: mtxel = 10
     integer :: rqstar = 11
     integer :: gmap = 12
     integer :: epsinv_total = 13
     integer :: chi_sum_comm = 14
     integer :: chi_sum_total = 15
     integer :: genwf_val = 16
     integer :: genwf_cond = 17
     integer :: epsinv_vcoul = 18
     integer :: job_setup = 19
     integer :: q_loop_setup = 20
     integer :: init_cutoff = 21
     integer :: init_scalapack = 22
     integer :: init_arrays = 23
     integer :: converge_tests = 24
     integer :: mtxel_denom = 25
     integer :: mtxel_fft = 26
     integer :: genwf_ekin = 28
     integer :: genwf_sort = 29
     integer :: chi_sum_gemm = 30
     integer :: chi_sum_prep = 31
     integer :: mtxel_exp_denom = 32
     integer :: mtxel_exp_fft = 33
     integer :: chi_sum_sub_vcoul = 34
     integer :: chi_sum_sub_diag = 35
     integer :: chi_sum_sub_omega_0 = 36
     integer :: chi_sum_sub_eigvet_comm = 37
     integer :: chi_sum_sub_transf = 38
     integer :: chi_sum_sub_omega_neq_0 = 39
     integer :: opt_fft = 40
     integer :: opt_fft_init = 41
     integer :: opt_fft_comm_fft = 42
     integer :: opt_fft_fft = 43
     integer :: chi_sum_array_alloc = 44
     integer :: epsinv_i_o = 45
     integer :: epsinv_invert = 46
     integer :: chi_read = 47
     integer :: chi_sum_bar = 49
     integer :: chi_sum_flt = 50
     integer :: chi_sum_row = 51
     integer :: chi_sum_column = 52
     integer :: chi_sum_ht_nb = 53
     integer :: subspace_pgemm = 60
     integer :: epsinv_omega_0 = 61
     integer :: epsinv_omega_neq_0 = 62
     ! Epsilon doesn`t use the put/mltply of Common.
     integer :: fft_put = 92
     integer :: fft_mltply = 95
     integer :: total = 100
   contains
     procedure :: init_labels => epsilon_init_labels
  end type timing_epsilon_class
  !
  !----------------------------------------------------------------------------
  !
  type, extends(timing_class) :: timing_sigma_class
     ! class definition for timing of subrourines in ./Sigma
     integer :: io_total = 1
     integer :: input = 2
     integer :: epscopy = 3
     integer :: fullbz = 4
     integer :: vxc = 5
     integer :: subgrp = 6
     integer :: irrbz = 7
     integer :: gmap = 8
     integer :: genwf = 9
     integer :: mtxel = 10
     integer :: mtxel_cor_tot = 11
     integer :: vcoul = 13
     integer :: epsread = 14
     integer :: input_outer = 15
     integer :: mtxel_ch = 6
     integer :: mtxel_comm =17
     integer :: bare_x = 18
     integer :: wf_comm = 19
     integer :: wf_ch_comm = 20
     integer :: input_read = 21
     integer :: input_write = 22
     integer :: sub_transf_tot = 31
     integer :: sub_transf_com = 32
     integer :: sub_transf_gemm = 33
     integer :: m_cor_init = 41
     integer :: m_cor_epsinit = 42
     integer :: m_cor_comm = 43
     integer :: m_cor_pp_prep = 44
     integer :: m_cor_sx_ch = 45
     integer :: m_cor_ra_sx = 46
     integer :: m_cor_ra_ch = 47
     integer :: m_cor_ra_ch2 = 48
     integer :: m_cor_ra_sum = 49
     integer :: m_cor_cd_res = 50
     integer :: m_cor_cd_int = 51
     integer :: m_cor_cd_sum = 52
     integer :: m_cor_cd_gemm = 53
     integer :: m_cor_remain = 55
     integer :: m_cor_sub_wings = 56
     integer :: read_neps = 59
     integer :: epscopy_io = 61
     integer :: epscopy_comm = 62
     integer :: epscopy_sub = 63
     integer :: epscopy_pgemm = 64
     integer :: epscopy_redstr = 65
     integer :: sub_io_vec = 66
     integer :: sub_prep_vec = 67
     integer :: sub_comm_vec = 68
     integer :: sub_io_eps = 69
     integer :: sub_prep_eps = 70
     integer :: sub_comm_eps = 71
     integer :: epscopy_vcoul = 72
     integer :: eqp_loop_total = 75
     integer :: total = 100
   contains
     procedure :: init_labels => sigma_init_labels
  end type timing_sigma_class
  !
  !----------------------------------------------------------------------------
  !
  type, extends(timing_class) :: timing_bse_class
     ! class definition for timing of subrourines in ./BSE
     integer :: input = 2
     integer :: input_q = 3
     integer :: intwfn = 4
     integer :: intkernel = 5
     integer :: epsdiag = 7
     integer :: eps_comm = 8
     integer :: absorp0 = 9
     integer :: vmtxel = 10
     integer :: trans_mtxel = 11
     integer :: absorp = 12
     integer :: write_eig = 13
     integer :: iw_input_co = 41
     integer :: iw_interp = 42
     integer :: iw_genwf = 43
     integer :: iw_genwf_co = 44
     integer :: iw_mtxel_t = 45
     integer :: iw_write = 46
     integer :: iw_reduce = 47
     integer :: ik_setup = 51
     integer :: ik_c_check = 52
     integer :: ik_input = 53
     integer :: ik_inteps = 54
     integer :: ik_vcoul = 55
     integer :: ik_cache = 56
     integer :: ik_interp = 57
     integer :: ik_sum = 58
     integer :: diagonalize = 61
     integer :: lanczos = 62
     integer :: iterate = 63
     integer :: peig_inter = 64
     integer :: total = 100
   contains
     procedure :: init_labels => bse_init_labels
  end type timing_bse_class
  !
  !----------------------------------------------------------------------------
  !
  type, extends(timing_class) :: timing_common_class
     ! class definition for timing of subrourines in ./Common
     integer :: eps_i_o_comm = 47
     integer :: eps_i_o_io = 48
     integer :: epscopy_comm = 62
     integer :: input_i_o = 81
     integer :: input_comm = 82
     integer :: fft_zero = 91
     integer :: fft_put = 92
     integer :: fft_plan = 93
     integer :: fft_exec = 94
     integer :: fft_mltply = 95
     integer :: fft_conjg = 96
     integer :: fft_get = 97
   contains
     procedure :: init_labels => common_init_labels
  end type timing_common_class
  !
  !----------------------------------------------------------------------------
  !
  type, extends(timing_class) :: timing_extra_class
    ! class definition for timing of subrourines in various parts of BGW
    ! that do not fall in folders ./Sigma, ./BSE, ./Common or ./Epsilon
    integer :: total = 100
    integer :: input = 2
    integer :: input_q = 3
    integer :: vmtxel = 4
    integer :: readasvck = 5
    integer :: os_comm = 6
    integer :: os_sums = 7
    integer :: genwf = 8
    integer :: genwf_q = 9
    integer :: summing = 9
    integer :: gather = 9
   contains
     procedure :: init_labels => extra_init_labels
  end type timing_extra_class
  !
  type, extends(timing_class) :: timing_tdgw_class
    ! class definition for timing of subroutines in ./TDGW
    integer :: input = 2
    integer :: init_gf_hamil = 3
    integer :: solve_bse = 4
    integer :: get_interaction = 5
    integer :: commutation = 6
    integer :: propagator = 7
    integer :: measurement = 8
    integer :: output = 9
    integer :: communication = 10
    integer :: test = 99  ! for testing purpose
    integer :: total = 100
   contains 
    procedure :: init_labels => tdgw_init_labels
  end type timing_tdgw_class
  !----------------------------------------------------------------------------
  !
  ! After the definition of the classes, these are the object instances used in
  ! the code.
  ! These should be used throughotu the BGW code by importing them as
  ! "use timing_m, only: timing => epsilon_timing"
  !
  type(timing_epsilon_class), save, public :: epsilon_timing
  type(timing_sigma_class), save, public :: sigma_timing
  type(timing_bse_class), save, public :: bse_timing
  type(timing_common_class), save, public :: common_timing
  type(timing_extra_class), save, public :: extra_timing
  type(timing_tdgw_class), save, public :: tdgw_timing
  !
contains
  !
  subroutine bare_init_labels(instance)
    ! abstract implementation
    implicit none
    class(timing_class), intent(inout) :: instance
    call die("Need a specific implementation of init_labels")
    return
  end subroutine bare_init_labels
  !
  !----------------------------------------------------------------------------
  !
  subroutine epsilon_init_labels(instance)
    ! Labels for attributes of timing_epsilon
    implicit none
    class(timing_epsilon_class), intent(inout) :: instance
    PUSH_SUB(epsilon_init_labels)

    instance%labels(1) = 'I/O TOTAL'
    instance%labels(2) = 'INPUT'
    instance%labels(3) = 'INPUT_Q'
    instance%labels(4) = 'FULLBZ'
    instance%labels(5) = 'GVEC'
    instance%labels(6) = 'SUBGRP'
    instance%labels(8) = 'IRRBZ'
    instance%labels(9) = 'GENWF'
    instance%labels(10) = 'MTXEL'
    instance%labels(11) = 'RQSTAR'
    instance%labels(12) = 'GMAP'
    instance%labels(13) = 'EPSINV (TOTAL)'
    instance%labels(14) = 'CHI SUM (COMM)'
    instance%labels(15) = 'CHI SUM (TOTAL)'
    instance%labels(16) = 'GENWF (VAL)'
    instance%labels(17) = 'GENWF (COND)'
    instance%labels(18) = 'EPSINV (VCOUL)'
    instance%labels(19) = 'JOB SETUP'
    instance%labels(20) = 'Q LOOP SETUP'
    instance%labels(21) = 'INIT CUTOFF'
    instance%labels(22) = 'INIT SCALAPACK'
    instance%labels(23) = 'INIT ARRAYS'
    instance%labels(24) = 'CONVERGE TESTS'
    instance%labels(25) = 'MTXEL (DENOM)'
    instance%labels(26) = 'MTXEL (FFT)'
    instance%labels(28) = 'GENWF (Ekin)'
    instance%labels(29) = 'GENWF (Sort)'
    instance%labels(30) = 'CHI SUM (' + TOSTRING(X(GEMM)) + ')'
    instance%labels(31) = 'CHI SUM (PREP)'
    instance%labels(32) = 'MTXEL EXP(DENOM)'
    instance%labels(33) = 'MTXEL EXP (FFT)'
    instance%labels(34) = 'CHI SUM SUB (VCOUL)'
    instance%labels(35) = 'CHI SUM SUB DIAG'
    instance%labels(36) = 'CHI SUM SUB OMEGA=0'
    instance%labels(37) = 'CHI SUM SUB EIGVET COMM'
    instance%labels(38) = 'CHI SUM SUB TRANSF'
    instance%labels(39) = 'CHI SUM SUB OMEGA neq 0'
    instance%labels(40) = 'OPT FFT'
    instance%labels(41) = 'OPT FFT (INIT)'
    instance%labels(42) = 'OPT FFT (COMM_FFT)'
    instance%labels(43) = 'OPT FFT (FFT)'
    instance%labels(44) = 'CHI SUM (ARRAY ALLOC)'
    instance%labels(45) = 'EPSINV (I/O)'
    instance%labels(46) = 'EPSINV (INVERT)'
    instance%labels(47) = 'READ CHI'
    instance%labels(49) = 'CHI SUM (BAR)'
    instance%labels(50) = 'CHI SUM (FLT)'
    instance%labels(51) = 'CHI SUM (ROW)'
    instance%labels(52) = 'CHI SUM (COLUMN)'
    instance%labels(53) = 'CHI SUM (HT/NB)'
    instance%labels(60) = 'SUBSPACE (P' + TOSTRING(X(GEMM)) + ')'
    instance%labels(61) = 'EPSINV OMEGA=0'
    instance%labels(62) = 'EPSINV OMEGA neq 0'
    instance%labels(92) = 'FFT PUT'
    instance%labels(95) = 'FFT MLTPLY'
    instance%labels(100) = 'TOTAL'
    !
    POP_SUB(epsilon_init_labels)
    return
  end subroutine epsilon_init_labels
  !
  !----------------------------------------------------------------------------
  !
  subroutine sigma_init_labels(instance)
    ! Labels for attributes of timing_sigma
    implicit none
    class(timing_sigma_class), intent(inout) :: instance
    PUSH_SUB(sigma_init_labels)
    !
    instance%labels(1) = 'I/O TOTAL'
    instance%labels(2) = 'INPUT'
    instance%labels(3) = 'EPSCOPY'
    instance%labels(4) = 'FULLBZ'
    instance%labels(5) = 'VXC'
    instance%labels(6) = 'SUBGRP'
    instance%labels(7) = 'IRRBZ'
    instance%labels(8) = 'GMAP'
    instance%labels(9) = 'GENWF'
    instance%labels(10) = 'MTXEL'
    instance%labels(11) = 'MTXEL_COR TOT'
    instance%labels(13) = 'VCOUL'
    instance%labels(14) = 'EPSREAD'
    instance%labels(15) = 'INPUT_OUTER'
    instance%labels(16) = 'MTXEL_CH'
    instance%labels(17) = 'MTXEL COMM'
    instance%labels(18) = 'BARE X'
    instance%labels(19) = 'WF COMM'
    instance%labels(20) = 'WF_CH COMM'
    instance%labels(21) = 'INPUT (READ)'
    instance%labels(22) = 'INPUT (WRITE)'
    instance%labels(31) = 'SUB-TRANSF TOT'
    instance%labels(32) = 'SUB-TRANSF COM'
    instance%labels(33) = 'SUB-TRANSF GEMM'
    instance%labels(41) = 'M.COR INIT'
    instance%labels(42) = 'M.COR EPSINIT'
    instance%labels(43) = 'M.COR COMM'
    instance%labels(44) = 'M.COR PP PREP'
    instance%labels(45) = 'M.COR SX+CH'
    instance%labels(46) = 'M.COR RA SX'
    instance%labels(47) = 'M.COR RA CH'
    instance%labels(48) = 'M.COR RA CH2'
    instance%labels(49) = 'M.COR RA SUM'
    instance%labels(50) = 'M.COR CD RES'
    instance%labels(51) = 'M.COR CD INT'
    instance%labels(52) = 'M.COR CD SUM'
    instance%labels(53) = 'M.COR CD GEMM'
    instance%labels(55) = 'M.COR REMAIN'
    instance%labels(56) = 'M.COR SUB WINGS'
    instance%labels(59) = 'READ NEPS'
    instance%labels(61) = 'EPSCOPY IO'
    ! Epscopy comm is a duplicate of common_timing:
    ! the common timing has the HDF5, this one has the binary
    instance%labels(62) = 'EPSCOPY COMM'
    instance%labels(63) = 'EPSCOPY SUB'
    instance%labels(64) = 'EPSCOPY PGEMM'
    instance%labels(65) = 'EPSCOPY REDSTR'
    instance%labels(66) = 'SUB IO Vec'
    instance%labels(67) = 'SUB Prep Vec'
    instance%labels(68) = 'SUB COMM Vec'
    instance%labels(69) = 'SUB IO Eps'
    instance%labels(70) = 'SUB Prep Eps'
    instance%labels(71) = 'SUB COMM Eps'
    instance%labels(72) = 'EPSCOPY VCOUL'
    instance%labels(75) = 'EQP LOOP TOTAL'
    instance%labels(100) = 'TOTAL'
    !
    POP_SUB(sigma_init_labels)
    return
  end subroutine sigma_init_labels
  !
  !----------------------------------------------------------------------------
  !
  subroutine bse_init_labels(instance)
    ! Labels for attributes of timing_bse
    implicit none
    class(timing_bse_class), intent(inout) :: instance
    PUSH_SUB(bse_init_labels)
    !
    instance%labels(2)='Input'
    instance%labels(3)='Input q'
    instance%labels(4)='Intwfn'
    instance%labels(5)='Intkernel'
    instance%labels(7)='Epsdiag'
    instance%labels(8)='Eps Comm'
    instance%labels(9)='Absorp0'
    instance%labels(10)='Vmtxel'
    instance%labels(11)='Trans Mtxel'
    instance%labels(12)='Absorp'
    instance%labels(13)='Write Eig'
    instance%labels(41)='Iw Input_co'
    instance%labels(42)='Iw Interp'
    instance%labels(43)='Iw Genwf'
    instance%labels(44)='Iw Gwnwf_Co'
    instance%labels(45)='Iw Mtxel_t'
    instance%labels(46)='Iw Write'
    instance%labels(47)='Iw Reduce'
    instance%labels(51)='Ik Setup'
    instance%labels(52)='Ik C-Check'
    instance%labels(53)='Ik Input'
    instance%labels(54)='Ik Inteps'
    instance%labels(55)='Ik Vcoul'
    instance%labels(56)='Ik Cache'
    instance%labels(57)='Ik Interp'
    instance%labels(58)='Ik Sum'
    instance%labels(61)='Diagonalize'
    instance%labels(62)='Lanczos'
    instance%labels(63)='Iterate'
    instance%labels(64)='Peig_Inter'
    instance%labels(100) = 'TOTAL'
    !
    POP_SUB(bse_init_labels)
    return
  end subroutine bse_init_labels
  !
  !----------------------------------------------------------------------------
  !
  subroutine common_init_labels(instance)
    ! Labels for attributes of timing_common
    implicit none
    class(timing_common_class), intent(inout) :: instance
    PUSH_SUB(common_init_labels)
    !
    instance%labels(47) = 'Eps (I/O) Comm'
    instance%labels(48) = 'Eps (I/O) IO'
    instance%labels(62) = 'Epscopy Comm'
    instance%labels(81) = 'Input I/O'
    instance%labels(82) = 'Input Comm'
    instance%labels(91) = 'FFT Zero'
    instance%labels(92) = 'FFT Put'
    instance%labels(93) = 'FFT Plan'
    instance%labels(94) = 'FFT Exec'
    instance%labels(95) = 'FFT Mltply'
    instance%labels(96) = 'FFT Conjg'
    instance%labels(97) = 'FFT Get'
    !
    POP_SUB(common_init_labels)
    return
  end subroutine common_init_labels
  !
  !----------------------------------------------------------------------------
  !
  subroutine extra_init_labels(instance)
    ! Labels for attributes of timing_extra
    implicit none
    class(timing_extra_class), intent(inout) :: instance
    PUSH_SUB(extra_init_labels)
    !
    instance%labels(2) = 'Input'
    instance%labels(3) = 'Input_q'
    instance%labels(4) = 'Vmtxel'
    instance%labels(5) = 'Readasvck'
    instance%labels(6) = 'OS - Comm'
    instance%labels(7) = 'OS - Sums'
    instance%labels(8) = 'Genwf'
    instance%labels(9) = 'Genwf_q'
    instance%labels(10) = 'Summing'
    instance%labels(11) = 'Gather'
    instance%labels(100) = 'TOTAL'
    !
    POP_SUB(extra_init_labels)
    return
  end subroutine extra_init_labels
  !
  !----------------------------------------------------------------------------
  !
  subroutine tdgw_init_labels(instance)
    ! Labels for attributes of timing_tdgw
    implicit none
    class(timing_tdgw_class), intent(inout) :: instance
    PUSH_SUB(tdgw_init_labels)
    !
    instance%labels(2) = 'INPUT'
    instance%labels(3) = 'INIT_GF_HAMIL'
    instance%labels(4) = 'SOLVE_BSE'
    instance%labels(5) = 'GET_INTERACTION'
    instance%labels(6) = 'COMMUTATION'
    instance%labels(7) = 'PROPAGATOR'
    instance%labels(8) = 'MEASUREMENT'
    instance%labels(9) = 'OUTPUT'
    instance%labels(10) = 'COMMUNICATION'
    instance%labels(99) = 'TEST'
    instance%labels(100) = 'TOTAL'
    !
    POP_SUB(tdgw_init_labels)
    return
  end subroutine tdgw_init_labels
  !
  !----------------------------------------------------------------------------
  !
  subroutine timing_class_init(instance)
    ! Initialize the timing methods
    ! Essentially, sets labels and sets times to zero
    implicit none
    class(timing_class), intent(inout) :: instance
    integer :: cm
    PUSH_SUB(timing_class_init)
    !
    instance%cpu_times = 0.0d0
    instance%wall_times = 0.0d0
    instance%tmp_wall_times = 0.0d0
    instance%tmp_cpu_times = 0.0d0
    instance%call_numbers = 0
    instance%labels = ''
    call instance%init_labels()
    !
    POP_SUB(timing_class_init)
    return
  end subroutine timing_class_init
  !
  !----------------------------------------------------------------------------
  !
  subroutine timing_class_print(instance, c_timing)
    ! Print to screen all the timing information.
    ! specifically, we will print the max, min (over the timing info of the
    ! MPI processes) and the root time associated to that tag
    ! Args:
    ! c_timing, optional: one could pass another timing object, e.g. the one
    !     for timing the calls in ./Common. The printing info will be merged.
    implicit none
    class(timing_class), intent(inout) :: instance
    type(timing_common_class), optional, intent(inout) :: c_timing
    !
    integer :: i, error, call_numbers(2*instance%num_times), N, N2
    integer, allocatable :: buffer_i(:)
    real(DP), allocatable :: buffer_r(:)
    real(DP) :: min_cpu_times(2*instance%num_times), &
         max_cpu_times(2*instance%num_times),             &
         root_cpu_times(2*instance%num_times),            &
         min_wall_times(2*instance%num_times),            &
         max_wall_times(2*instance%num_times),            &
         root_wall_times(2*instance%num_times)
    character(len=100) :: labels(2*instance%num_times)
    !
    PUSH_SUB(timing_class_print)
    !
    labels = ''
    N2 = 2*instance%num_times
    N = instance%num_times
    !
    ! store times and labels in a temporary array.
    !
    ! The simplest thing to do to merge with c_timing, given also that these
    ! arrays are small, is to create a bigger array that accomodates everything
    !
    min_cpu_times(N+1:) = instance%cpu_times
    max_cpu_times(N+1:) = instance%cpu_times
    root_cpu_times(N+1:) = instance%cpu_times
    min_wall_times(N+1:) = instance%wall_times
    max_wall_times(N+1:) = instance%wall_times
    root_wall_times(N+1:) = instance%wall_times
    call_numbers(N+1:) = instance%call_numbers
    labels(N+1:) = instance%labels
    !
    if ( present(c_timing) ) then
       min_cpu_times(:N) = c_timing%cpu_times
       max_cpu_times(:N) = c_timing%cpu_times
       root_cpu_times(:N) = c_timing%cpu_times
       min_wall_times(:N) = c_timing%wall_times
       max_wall_times(:N) = c_timing%wall_times
       root_wall_times(:N) = c_timing%wall_times
       call_numbers(:N) = c_timing%call_numbers
       labels(:N) = c_timing%labels
    end if
    !
#ifdef MPI
     SAFE_ALLOCATE(buffer_r,(N2))
     SAFE_ALLOCATE(buffer_i,(N2))
     call MPI_allreduce(max_cpu_times, buffer_r, N2, MPI_double_precision,  &
          MPI_max, mpi_comm_world, error)
     max_cpu_times = buffer_r
     call MPI_allreduce(min_cpu_times, buffer_r, N2, MPI_double_precision,  &
          MPI_min, mpi_comm_world, error)
     min_cpu_times = buffer_r
     call MPI_allreduce(max_wall_times, buffer_r, N2, MPI_double_precision, &
          MPI_max, mpi_comm_world, error)
     max_wall_times = buffer_r
     call MPI_allreduce(min_wall_times, buffer_r, N2, MPI_double_precision, &
          MPI_min, mpi_comm_world, error)
     min_wall_times = buffer_r
     call MPI_allreduce(call_numbers, buffer_i, N2, MPI_integer, MPI_max,   &
          mpi_comm_world, error)
     call_numbers = buffer_i
     SAFE_DEALLOCATE(buffer_r)
     SAFE_DEALLOCATE(buffer_i)
#endif

    ! print to screen
    if ( peinf%inode == 0 ) then
       print*, ''
       print*, 'Timing information'
       print*, ''

!                                CPU time (s)          WALL time (s)       Number
! Routine                       min        max        min        max     of calls
! -------------------------  ---------  ---------  ---------  ---------  --------
! - Input I/O                999999.99  999999.99  999999.99  999999.99  88888888


      write(6,'(a)') &
'                                CPU time (s)          WALL time (s)       Number'
      write(6,'(a)') &
' Routine                       min        max        min        max     of calls'
      write(6,'(a)') &
' -------------------------  ---------  ---------  ---------  ---------  --------'

       do i = 1, N2
          ! if the label is empty, or it`s never been called, we don`t print it
          if ( len_trim(labels(i)) == 0 ) cycle
          if ( call_numbers(i) == 0 ) cycle
           write(6,"(1x,'- ',a23,4(2x,f9.2),2x,i8)") labels(i), &
             min_cpu_times(i), max_cpu_times(i), &
             min_wall_times(i), max_wall_times(i), call_numbers(i)
       end do
       print*, ''
       print*, 'Job Done'
       print*, ''
    end if
    !
    POP_SUB(timing_class_print)
    return
  end subroutine timing_class_print
  !
  !----------------------------------------------------------------------------
  !
  subroutine timing_class_start(instance, tag)
    ! start measuring time for tag
    implicit none
    class(timing_class), intent(inout) :: instance
    integer, intent(in) :: tag
    !
    integer :: i, j, values(8), wt0
    real(DP) :: cpu, wall
    !
    ! increase the counter by 1
    instance%call_numbers(tag) = instance%call_numbers(tag) + 1
    ! save the initial time
    call timget(cpu,wall)
    instance%tmp_cpu_times(tag) = cpu
    instance%tmp_wall_times(tag) = wall
    !
    return
  end subroutine timing_class_start
  !
  !----------------------------------------------------------------------------
  !
  subroutine timing_class_stop(instance, tag)
    ! Stop measuring time for this tag
    !
    implicit none
    class(timing_class), intent(inout) :: instance
    integer, intent(in) :: tag
    real(DP) :: cpu, wall, t0_c, t0_w
    !
    call timget(cpu, wall)
    t0_c = instance%tmp_cpu_times(tag)
    t0_w = instance%tmp_wall_times(tag)
    instance%cpu_times(tag)  = instance%cpu_times(tag)  + (cpu  - t0_c)
    instance%wall_times(tag) = instance%wall_times(tag) + (wall - t0_w)
    !
    return
  end subroutine timing_class_stop
  !
  !----------------------------------------------------------------------------
  !
  subroutine date_time(bdate,btime)
    ! returns:
    ! - bdate: string with date
    ! - btime: string with time
    character(len=11), intent(out) :: bdate
    character(len=14), intent(out) :: btime
    !
    integer :: lmonth
    integer :: idate(8)
    character(len=10) :: atime
    character(len=8) :: adate
    character(len=5) :: azone
    character(len=4) :: year
    character(len=3) :: month(12)
    character(len=2) :: hour, min, sec, day
    !
    DATA month/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep', &
         'Oct','Nov','Dec'/
    !
    PUSH_SUB(date_time)

    call date_and_time(adate,atime,azone,idate)
    read(adate,"(a4,i2,a2)") year, lmonth, day
    write(bdate,"(a2,a1,a3,a1,a4)") day, '-', month(lmonth), '-', year
    read(atime,'(a2,a2,a2,a4)') hour, min, sec
    write(btime,"(a2,a1,a2,a1,a2,1x,a5)") hour, ':', min, ':', sec, azone
    !
    POP_SUB(date_time)
    return    
  end subroutine date_time

  !============================================================================

  subroutine timget(cpu, wall)
    real(DP), intent(out) :: cpu, wall
    integer :: values(8)
    ! no push_sub, called too frequently
    !
    call cpu_time(cpu)
    call date_and_time(VALUES=values)
    wall=((values(3)*24.0d0+values(5))*60.0d0 &
         +values(6))*60.0d0+values(7)+values(8)*1.0d-3
    !
    return
  end subroutine timget

  !============================================================================

  subroutine timacc(n, option, tottim, nslices)
    ! DEPRECATED
    ! old subroutine for measuring execution time
    integer, intent(in) :: n !< not used for option = 0
    integer, intent(in) :: option !< 0, 1, 2, 3, 4
    real(DP), intent(out), optional :: tottim(2) !present if option=3 or 4
    integer, intent(out), optional :: nslices !< optionally used when option=3
    !
    real(DP) :: cpu,wall
    character(len=100) :: tmpstr
    !
    ! no push_sub, called too frequently

    ! Check that n lies in sensible bounds

    if (n .lt. 0 .or. n .gt. MTIM) then
       write(tmpstr,'(a,i6,a,i8)')'timacc: dim MTIM = ',MTIM,' but input n =',n
       call die(tmpstr)
    end if

    if (option==0) then

       ! Zero out all accumulators of time and init timers

       acctim(:,:)=0.0d0
       tzero(:,:)=0.0d0
       ncount(:)=0

    else if (option==1) then

       ! Initialize timepw for n

       call timget(cpu,wall)
       tzero(1,n)=cpu
       tzero(2,n)=wall

    else if (option==2) then

       ! Accumulate time for n

       call timget(cpu,wall)
       acctim(1,n)=acctim(1,n)+cpu -tzero(1,n)
       acctim(2,n)=acctim(2,n)+wall-tzero(2,n)
       ncount(n)=ncount(n)+1

    else if (option==3) then

       ! Return accumulated time for n

       if(.not. present(tottim)) call die("timacc requires tottim for option 3.")

       tottim(1)=acctim(1,n)
       tottim(2)=acctim(2,n)
       if(present(nslices)) then
          nslices=ncount(n)
       end if

    else if (option==4) then

       ! Return elapsed time for n (do not accumulate)

       if(.not. present(tottim)) call die("timacc requires tottim for option 4.")

       call timget(cpu,wall)
       tottim(1)=cpu-tzero(1,n)
       tottim(2)=wall-tzero(2,n)

    else

       write(tmpstr,'(a,i10,a)') 'timacc: input option = ', option, 'not valid.'
       call die(tmpstr)

    end if

    return
  end subroutine timacc

  !============================================================================

  subroutine logit(str, should_print, iunit)
    character (len=*), intent(in) :: str
    logical, intent(in), optional :: should_print
    integer, intent(in), optional :: iunit

    character(len=15) :: mydate,mytime,tmpstr
    logical :: should_print_
    integer :: iunit_

    if ( .not. peinf%verb_log ) return

    iunit_ = 6
    if (present(iunit)) iunit_ = iunit
    should_print_ = peinf%inode==0
    if (present(should_print)) should_print_ = should_print

    if (should_print_) then
       call date_and_time(mydate,mytime)
       tmpstr = mytime(1:2)//':'//mytime(3:4)//':'//mytime(5:6)//'.'//mytime(8:10)
       mytime = tmpstr
       write(iunit_,'(4a)') '*** LOG: ', TRUNC(str),'  time = ', TRUNC(mytime)
    endif

  end subroutine logit
  !
  !============================================================================
  !
  subroutine logitint(str,i)
    character(len=*), intent(in) :: str
    integer, intent(in) :: i
    character(len=100) :: tmpstr
    !
    if (.not.peinf%verb_log) return
    write(tmpstr,'(a,i5)') str(1:len_trim(str)),i
    call logit(tmpstr)
    !
    return
  end subroutine logitint

  !============================================================================

end module timing_m
