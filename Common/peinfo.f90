!================================================================================
!
! Modules:
!
! (1) peinfo_m      Originally by DAS 8/20/2010
!
!     Defines type and global instance of object for "processor equivalent" info.
!     Use mpi module to define interfaces for MPI calls.
!     [For F77, MPI header 'mpif.h' was included.]
!
!================================================================================

#include "f_defs.h"

module peinfo_m

#ifdef MPI
  use mpi
! include "mpif.h" -- the old way, which will not provide interfaces
#endif
  USEOMPLIB
#ifdef ACC
  use openacc
#endif
  use nrtype_m
  implicit none

  ! default must not be private, or else the types defined in mpi module will not be available.

  public ::      &
    peinfo,      &
    peinfo_init, &
    create_mpi_group, &
    get_thread_id

!-------------------------------

  type peinfo
    !> default values for serial
    integer :: npes = 1
    integer :: ngpus = 0
    integer :: npes_freqgrp = 1
    integer :: nthreads = 1
    integer :: nthreads_sort = 1 !< Set with `export BGW_NUM_THREADS_SORT=<num>`
    integer :: inode = 0
    !> Verbosity level, not to be used directly. Use the verbosity flags instead.
    integer :: verbosity=1
    logical :: verb_medium=.false.
    logical :: verb_high=.false.
    logical :: verb_log=.false.
    logical :: verb_debug=.false.
    logical :: verb_max=.false.
    !> Use new algorithm throughout BGW when writing HDF5 FF matrices?
    logical :: hdf5_write_redist = .false.
    !> Use H5FD_MPIO_INDEPENDENT to read WFNs
    logical :: use_wfn_hdf5_independent = .false.
    !> Use to change max memory to read contigous chuncks of bands for HDF5 matrices (see Common/read_matrix.f90)
    real(DP) :: maxmem_read_h5_mat = 512d0
    !> Communication style when reading bands in the HDF5 format. Controlled by
    !! the BGW_HDF5_BANDS_COMM_STYLE environment variable. Options are:
    !! 0=native HDF5; 1=manual group comm; 2=manual send/recvs
    integer :: hdf5_bands_comm_style = 1
    !> initialize to zero, then keep track of memory
    real(DP) :: mymem = 0d0
    real(DP) :: mymaxmem = 0d0
    integer :: nckmem
    integer :: nkpe  !< number of k-points per processor, used in absorption only
    !> kernel: total number of block-transitions ( nk^2, (nk*nc)^2 or (nk*nc*nv)^2)
    !! Each block-transition has iholdperown 
    integer :: nck
    !> kernel: number of block-transitions that I own
    integer :: nckpe
    integer :: myown !< Kernel: number of unique (k,kp) pairs I own; BSE: number of blocks I own
    integer :: mypown !< in BSE, number of unprimed indices I own for all my blocks
    integer :: npown !< in BSE, max number of unprimed indices owned by any proc in my pool
    integer :: jobtypeeval  !< Coulomb potential generator
                            !! 0 = "dumb generator" ; 1 = "smart consumer"
                            !! See vcoul_generator_m for more details
    !> BSE: number of blocks I own in the one-dimentional block-cyclic distributed
    !! matrices hmtx_a/evecs_r.
    integer :: nblocks
    !> BSE: size of each block in the one-dimentional block-cyclic distributed
    !! matrices hmtx_a/evecs_r = ns*nc_block*nv_block, which varies according to ipar.
    integer :: block_sz
    !> kernel: (nv,nc,nk,nv,nc,nk) offset in the bse_matrix for the
    !! block-transition identified by (ivp,icp,ikp,iv,ic,ik)
    integer, pointer :: wown(:,:,:,:,:,:)
    integer, pointer :: ciown(:)
    integer, pointer :: ik(:,:) !< (inode,j) index of jth k owned by inode
    integer, pointer :: ic(:,:) !< (inode,j) index of jth cband owned by inode
    integer, pointer :: iv(:,:) !< (inode,j) index of jth vband owned by inode
    integer, pointer :: ikp(:,:) !< (inode,j) index of jth kp owned by inode
    integer, pointer :: icp(:,:) !< (inode,j) index of jth cpband owned by inode
    integer, pointer :: ivp(:,:) !< (inode,j) index of jth vpband owned by inode
    integer, pointer :: ib(:,:)
    integer, pointer :: ick(:,:)
    integer, pointer :: ipe(:)
    !> (inode,iv,ik) Maps the global index for valence band (iv) at kpt (ik) to
    !! the local list of valence band the proc owns. (ik) is defined in the 
    !! reducible wedge. ipec is 0 if the proc doesn`t own that band/kpt
    integer, pointer :: ipec(:,:,:)
    integer, pointer :: ipev(:,:,:) !< See ipec
    integer, pointer :: ipek(:,:)   !< Same as ipec, but w/o band index
    integer, pointer :: ipekq(:,:)  !< Local index of valence band k-points only used 
                                    !< for finite momemtnum calculations
    integer, pointer :: ipecb(:,:)
    integer, pointer :: ivckpe(:)
    !> (npes) Number of k-points in the full fine grid that each processors owns.
    !! This parallelization is only used for the WFN interpolation in BSE, and
    !! it has nothing to do with the ikt array used in the hbse_a matrix.
    integer, pointer :: ikt(:)
    !> (npes) Number of block-columns of the hbse_a matrix each processors owns.
    !! Used in BSE only. The size of each block is block_sz.
    integer, pointer :: ibt(:)
    !> (nblocks) ikb(ib) is the k-point associated to the ib-th block of the
    !! distributed BSE Hamiltonian that I own.
    integer, pointer :: ikb(:)
    !> (nblocks) icb(ib) is the cond band associated to the ib-th block of the
    !! distributed BSE Hamiltonian that I own. Used only if ipar==2 or ipar==3.
    integer, pointer :: icb(:)
    !> (nblocks) ivb(ib) is the val band associated to the ib-th block of the
    !! distributed BSE Hamiltonian that I own. Used only if ipar==3.
    integer, pointer :: ivb(:)
    !> Number of cond bands in each block of the distributed BSE Hamiltonian.
    !! This is xct%ncb_fi for ipar<2, and 1 for ipar>=2
    integer :: nc_block
    !> Number of val bands in each block of the distributed BSE Hamiltonian.
    !! This is xct%nvb_fi for ipar<3, and 1 for ipar>=3
    integer :: nv_block
    integer, allocatable :: neig(:)
    integer, allocatable :: peig(:,:)
    integer :: npools !< number of pools for the valence bands in Epsilon or outer bands in sigma
    integer :: npes_pool !< number of processors per pool
    integer :: pool_group !< mpi_group for pools
    integer :: pool_comm !< mpi_comm for pools
    integer :: pool_rank !< rank within pool
    integer :: my_pool !< what pool this processor is in 
    integer :: nvownmax  !< max. number of valence bands that I can own
    integer :: ncownmax  !< max. number of conduction bands that I can own
    integer :: nvownactual !< (total) number of valence bands that I *really* own
    integer :: ncownactual !< (total) number of conduction bands that I *really* own
    !> Who owns a particular pair of bands (v,c)?
    integer, pointer :: global_pairowner(:,:)
    !> (total) number of valence bands that a particular MPI process owns
    integer, pointer :: global_nvown(:)
    !> (total) number of conduction bands that a particular MPI process owns
    integer, pointer :: global_ncown(:)
    !> indexv(i) is the local index (in terms of bands that I own) of the ith 
    !! (global) valence band. It is zero if I don`t own valence band #i.
    integer, pointer :: indexv(:)
    integer, pointer :: global_indexv(:,:) !< local indices for all processes
    integer, pointer :: indexc(:) !< see indexv
    !> Given a local band #i that I own, invindexv(i) is the global index of
    !! that band. If i>nvownt, the result is zero.
    integer, pointer :: invindexv(:)
    integer, pointer :: invindexc(:) !< see invindexv
    logical, pointer :: doiownv(:) !< do I own a particular valence band?
    logical, pointer :: doiownc(:) !< do I own a particular conduction band?
    logical, pointer :: does_it_ownc(:,:) !< (band,node) does a particular node own a cond. band?
    logical, pointer :: does_it_ownv(:,:) !< (band,node) does a particular node own a val. band?
    integer, pointer :: iownwfv(:) !< number of val. WFNs each proc. owns
    integer, pointer :: iownwfc(:) !< number of cond WFNs each proc. owns
    integer, pointer :: iownwfk(:) !< number of distinct k-points each proc. (partially) owns
    integer, pointer :: iownwfkq(:) !< Same as iownwfk, but refers to k+Q point when using finite momentum Q
    integer, pointer :: nxqown(:)
    integer, pointer :: nxqi(:)
    integer :: ndiag_max
    integer :: noffdiag_max
    integer :: ntband_max
    integer :: ntband_node
    integer :: nvband_node
    integer, pointer :: indext(:)
    integer, pointer :: ntband_dist(:)
    integer, pointer :: indext_dist(:,:)
    integer, pointer :: index_diag(:)
    logical, pointer :: flag_diag(:)
    integer, pointer :: index_offdiag(:)
    logical, pointer :: flag_offdiag(:)
    !> Parallel frequencies mpi group variables
    !! igroup = your group number
    !! rank = your processor number in your group
    !! _f = frequency evaluation group
    !! _mtxel = matrix element communication group
    integer :: igroup_f
    integer :: rank_f
    integer :: igroup_mtxel
    integer :: rank_mtxel
    integer :: mtxel_comm         !< mtxel group communicator
    integer :: freq_comm          !< frequency group communicator
    integer :: npes_orig          !< original number of processors
                                  !! for when nfreq_group does not
                                  !! divide total number of procs
    integer :: mtxel_group        !< mtxel group handle
    integer :: freq_group         !< frequency group handle
    integer, pointer :: ranks(:)  !< ranks of processors to include in mpi group
    logical :: check_norms=.true. !< Whether to check norms, .true. unless doing pseudobands
    !> default values for accelerator run (GPU)
    integer :: accel_num_devices=0        !< number of accelerator (GPU)
    integer :: accel_my_active_device=0  !< my active accelerator MOD(para_env%mepos, accel_num_devices)
    integer :: total_number_devices=0
    integer :: mpi_tasks_per_device=-1
  end type peinfo
  
  type(peinfo), save, public :: peinf
#ifdef MPI
  integer, public :: mpistatus(MPI_STATUS_SIZE)
  integer, public :: mpierr
#endif

contains


  !> FHJ: Set verbosity flags, such as peinf%verb_medium, based on peinf%verbosity.
  !! Note that verbosity flags are cumulative.
  subroutine peinfo_set_verbosity()
    character(len=8) :: verb_str(6)
    ! cannot use push_pop because that module uses this one

    if (peinf%verbosity<1) peinf%verbosity = 1
    if (peinf%verbosity>6) peinf%verbosity = 6
    if (peinf%verbosity>=2) peinf%verb_medium = .true.
    if (peinf%verbosity>=3) peinf%verb_high = .true.
    if (peinf%verbosity>=4) peinf%verb_log = .true.
    if (peinf%verbosity>=5) peinf%verb_debug = .true.
    if (peinf%verbosity>=6) peinf%verb_max = .true.
#ifdef VERBOSE
    ! FHJ: -DVERBOSE flag overwrites everything. This is useful for buildbots.
    peinf%verb_medium = .true.
    peinf%verb_high = .true.
    peinf%verb_log = .true.
    peinf%verb_debug = .true.
    peinf%verb_max = .true.
#endif
    if (peinf%inode==0) then
      verb_str(1) = "default"
      verb_str(2) = "medium"
      verb_str(3) = "high"
      verb_str(4) = "log"
      verb_str(5) = "debug"
      verb_str(6) = "max"
      write(6,'(1x,a,i0,3a/)') 'Running with verbosity level ', &
        peinf%verbosity,' (', trim(verb_str(peinf%verbosity)), ').'
      if (peinf%verbosity>3) then
        write(0,'(/a)') 'WARNING: you are running the calculation with a high level of verbosity.'
        write(0,'(a/)') 'This will impact the performance of the code.'
      endif
    endif

  end subroutine peinfo_set_verbosity


  subroutine peinfo_init()
    character(len=256) :: nthreads_sort_str, tmp_str
    integer :: err
#ifdef OMP
    integer, external :: omp_get_num_threads
#endif
    ! cannot use push_pop because that module uses this one

#ifdef MPI
    mpierr = 0
    call MPI_Init(mpierr)
    if(mpierr .ne. MPI_SUCCESS) then
      write(0,'(a)') 'ERROR: MPI initialization failed!'
      stop 999
    endif
    call MPI_Comm_rank(MPI_COMM_WORLD, peinf%inode, mpierr)
    call MPI_Comm_size(MPI_COMM_WORLD, peinf%npes, mpierr)
#endif

! Why put OMP pragmas here?
! JRD: I want to make sure our code has a parallel region before that of any library. This affects
! performance when the libraries are using a different implementation of threads or OpenMP build.
! if serial, default values set in type peinfo above are left alone
#ifdef ACC
    peinf%ngpus = acc_get_num_devices(acc_device_not_host)
#endif

#ifdef OMP
!$OMP PARALLEL
    peinf%nthreads = omp_get_num_threads()
!$OMP END PARALLEL
#endif

#ifdef OMP
    ! FHJ: The user can specify the number of threads used for sorting functions.
    ! We allow for this granularity because some functions, such as LAPACK,
    ! gain some performance with hyperthreads, but not our threaded sorting.
    call get_environment_variable("BGW_NUM_THREADS_SORT", nthreads_sort_str)
    read(nthreads_sort_str,*,iostat=err) peinf%nthreads_sort
    ! Default value if BGW_NUM_THREADS_SORT is not set
    ! FHJ: FIXME - currently disabled by default due to segfaults
    !if (err/=0) peinf%nthreads_sort = peinf%nthreads
    if (err/=0) peinf%nthreads_sort = 1
    ! Default value if BGW_NUM_THREADS_SORT is invalid
    !if (peinf%nthreads_sort<1.or.peinf%nthreads_sort>peinf%nthreads) &
    !  peinf%nthreads_sort = peinf%nthreads
    if (peinf%nthreads_sort<1.or.peinf%nthreads_sort>peinf%nthreads) &
      peinf%nthreads_sort = 1
#endif

    ! FHJ: When writing parallel FF matrices with HDF5:
    ! use new algorithm based on data redistribution?
    call get_environment_variable("BGW_HDF5_WRITE_REDIST", tmp_str)
    if (len_trim(tmp_str)>0) then
      if ((tmp_str(1:1)=='t').or.(tmp_str(1:1)=='T').or.(tmp_str(1:1)=='y').or.&
          (tmp_str(1:1)=='Y').or.(tmp_str(1:1)=='1')) then
        peinf%hdf5_write_redist = .true.
      endif
      if ((tmp_str(1:1)=='f').or.(tmp_str(1:1)=='F').or.(tmp_str(1:1)=='n').or.&
          (tmp_str(1:1)=='N').or.(tmp_str(1:1)=='0')) then
        peinf%hdf5_write_redist = .false.
      endif
    endif

    call get_environment_variable("BGW_WFN_HDF5_INDEPENDENT ", tmp_str)
    if (len_trim(tmp_str)>0) then
      if ((tmp_str(1:1)=='t').or.(tmp_str(1:1)=='T').or.(tmp_str(1:1)=='y').or.&
          (tmp_str(1:1)=='Y').or.(tmp_str(1:1)=='1')) then
        peinf%use_wfn_hdf5_independent = .true.
      endif
      if ((tmp_str(1:1)=='f').or.(tmp_str(1:1)=='F').or.(tmp_str(1:1)=='n').or.&
          (tmp_str(1:1)=='N').or.(tmp_str(1:1)=='0')) then
        peinf%use_wfn_hdf5_independent = .false.
      endif
    endif

    ! FHJ: Which algorithm to use when reading WFN in HDF5 format?
    peinf%hdf5_bands_comm_style = 1
    call get_environment_variable("BGW_HDF5_BANDS_COMM_STYLE", tmp_str)
    if (len_trim(tmp_str)>0) then
      read(tmp_str,*,iostat=err) peinf%hdf5_bands_comm_style
      if (err/=0) peinf%hdf5_bands_comm_style = 1
      if (peinf%hdf5_bands_comm_style<0 .or. peinf%hdf5_bands_comm_style>2) then
        peinf%hdf5_bands_comm_style = 1
      endif
    endif

  end subroutine peinfo_init


  subroutine create_mpi_group(orig_group,group_size,ranks,group_handle,group_comm)
    integer, intent(in) :: orig_group    !< Handle for original MPI group, which you are breaking into smaller groups
    integer,intent(in)  :: group_size    !< number of processors in new mpi group
    integer,intent(in)  :: ranks(:)      !< (group_size) array specifying ranks of processors to include in MPI group 
    integer,intent(out) :: group_handle  !< handle for new MPI group 
    integer,intent(out) :: group_comm    !< communicator for new MPI group 

#ifdef MPI
! DVF : create new group from original group, using ranks specified in `ranks` array 
    call MPI_Group_incl(orig_group, group_size,ranks(1:group_size), group_handle, mpierr)
    if(mpierr .ne. MPI_SUCCESS) write(0,'(a)') "ERROR: mpi_group_incl failed!"
! DVF : create communicator for new group
    call MPI_Comm_create(MPI_COMM_WORLD,group_handle,group_comm,mpierr)
    if(mpierr .ne. MPI_SUCCESS) write(0,'(a)') "ERROR: mpi_comm_create failed!"
#else
    group_handle = -1
    group_comm = -1
#endif

    return
  end subroutine create_mpi_group


  !> FHJ: Gets the thread number of the caller. Should be called inside an OMP
  !! construct. Returns 0 if code was compiled without OMP support.
  integer function get_thread_id()
#ifdef OMP
    integer, external :: omp_get_thread_num
#endif

    ! cannot use push_pop because that module uses this one

#ifdef OMP
    get_thread_id = omp_get_thread_num()
#else
    get_thread_id = 0
#endif

    ! cannot use push_pop because that module uses this one

  end function get_thread_id

end module peinfo_m
