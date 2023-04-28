!===============================================================================
!
! Routines:
!   Contains all input routines for dynamics simulation      
!
! (1) read_input()
!     read dynamics.inp file
!
! (2) load_bse_kernel()   
!     load electron-electron interaction BSE kernel matrix elements, 
!     including finite q in excitons (and excitations beyond)
!
! (3) load_eph_matrix()
!     load electron-phonon matrix elements from GWPT or DFPT, 
!     phonon q overlays with exciton finite q 
!
!===============================================================================

#include "f_defs.h"

module input_dynamics_m

#ifdef CPLX

  use global_m
  use hdf5
  use timing_m, only: timing => tdgw_timing
  use equation_of_motion_m
  use utility_m, only: kindex

  implicit none

  private

  public :: &
    read_input, &
    load_bse_kernel, &
    load_eqp_vmtxel, &
    load_eph_matrix

contains

!=========================================================================
!> read input file dynamics.inp
subroutine read_input(dyn)

  type(dynaminfo), intent(out) :: dyn
  integer :: iostat
  character*256 :: line, keyword, blockword
  integer :: iq
  real(DP) :: qptsum(3)

  PUSH_SUB(read_input)

  call timing%start(timing%input)

!-----------------------------
! set default values

  dyn%calc_dyn = 1
  dyn%n_finite_q = 0
  dyn%tmax = 100.0d0
  dyn%dt = 0.1d0
  dyn%pulse_type = 1
  dyn%eamp = 0.1d0
  dyn%efreq = 0.0d0
  dyn%ewidth = 1.0d0
  dyn%tpulse = 10.0d0
  dyn%iq_efield = (/1, 1/)
  dyn%band_index_min = 0
  dyn%band_index_max = 0
  dyn%prop_algo = 2

!-----------------------------
! all nodes read dynamics.inp

  call open_file(unit=10, file='dynamics.inp', form='formatted', status='old')

  if(peinf%inode == 0) then
    write(*,'(1x,a)') 'Reading dynamics.inp file'
    write(*,*)
  endif

  do
    ! all processors read input
    read(10,'(a256)',iostat=iostat) line
    if (iostat < 0) exit ! exit do loop at end-of-file

    ! skip empty lines
    if(len_trim(line) .eq. 0) cycle

    ! skip comment lines
    if(line(1:1) .eq. '#') cycle

    ! determine keyword
    keyword=line(1:scan(line," ")-1)
    line=adjustl(line(scan(line," ")+1:256))

    ! match and find keyword
    if (trim(keyword).eq.'begin') then
      blockword=line(1:scan(line," ")-1)
      if (blockword.ne.'finite_q_points') then
        call die('Please list finite_q_points in dynamics.inp.', only_root_writes=.false.)
      endif

      if (dyn%n_finite_q .eq. 0) then
        call die('Please set n_finite_q value before the finite_q_points block.', only_root_writes=.false.)
      endif

      SAFE_ALLOCATE(dyn%qpts, (3, dyn%n_finite_q))

      iq = 0
      do while(trim(line).ne.'end')
        read(10,'(a256)',iostat=iostat) line

        if (trim(line).ne.'end') then
          iq = iq + 1
          read(line,*) dyn%qpts(:,iq)
        endif
      enddo

      if (iq .ne. dyn%n_finite_q) then
        call die('Number of finite q points mismatch between n_finite_q flag and finite_q_points block.', only_root_writes=.false.)
      endif
      ! done reading finite_q_points block
    elseif(trim(keyword).eq.'calc_dyn') then
      read(line,*) dyn%calc_dyn
    elseif(trim(keyword).eq.'n_finite_q') then
      read(line,*) dyn%n_finite_q
    elseif(trim(keyword).eq.'tmax') then
      read(line,*) dyn%tmax
    elseif(trim(keyword).eq.'dt') then
      read(line,*) dyn%dt
    elseif(trim(keyword).eq.'pulse_type') then
      read(line,*) dyn%pulse_type
    elseif(trim(keyword).eq.'eamp') then
      read(line,*) dyn%eamp
    elseif(trim(keyword).eq.'efreq') then
      read(line,*) dyn%efreq
    elseif(trim(keyword).eq.'ewidth') then
      read(line,*) dyn%ewidth
    elseif(trim(keyword).eq.'tpulse') then
      read(line,*) dyn%tpulse
    elseif(trim(keyword).eq.'iq_efield') then
      read(line,*) dyn%iq_efield(:)
    elseif(trim(keyword).eq.'band_index_min') then
      read(line,*) dyn%band_index_min
    elseif(trim(keyword).eq.'band_index_max') then
      read(line,*) dyn%band_index_max
    elseif(trim(keyword).eq.'prop_algo') then
      read(line,*) dyn%prop_algo
    else
      call die('Unrecogonized input flag '//keyword, only_root_writes = .false.)
    endif
  enddo

  dyn%nt = int(dyn%tmax / dyn%dt) + 1  ! contains the starting equilibrium status

  call close_file(10)

  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  ! here we print some input parameters
  if(peinf%inode == 0) then
    write(*,'(3x,a)') 'Input parameters:'
    write(*,'(3x,a,2x,i5)')    'dynamics calculation type', dyn%calc_dyn
    write(*,'(3x,a,2x,i5)')    'number of finite q', dyn%n_finite_q
    write(*,'(3x,a,2x,f10.3)') 'simulation time', dyn%tmax
    write(*,'(3x,a)') 'list of finite_q_points' 
    do iq = 1, dyn%n_finite_q
      write(*,'(5x, 3f15.10)'), dyn%qpts(:, iq)
    enddo
    write(*,*)
  endif

  ! check hermitian field
  qptsum(:) = dyn%qpts(:, dyn%iq_efield(1)) + dyn%qpts(:, dyn%iq_efield(2))
  if (sqrt(dot_product(qptsum, qptsum)) .ge. TOL_Small) then
    call die('Input two components iq_efield do not form Hermitian light field.')
  endif

  call timing%stop(timing%input)

  POP_SUB(read_input)

end subroutine read_input


!=========================================================================
!> load BSE kernel matrix elements, which are e-e scattering amplitudes
subroutine load_bse_kernel(dyn, iq, tdevol)

  type(dynaminfo), intent(in) :: dyn          ! input parameters
  integer, intent(in) :: iq  ! q-point label 
                             ! TODO add reading paths for different q
  type(td_equation), intent(inout) :: tdevol  ! time-dependent evolution object 

  character*256 :: kernel_fname   ! temporarily reading arranged file from Yang-Hao's preprocessing script
  character*256 :: bsemat_fname   ! reading header information
  character*256 :: iq_str         ! string format of iq
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dataspace_id  ! dataspace identifier 
  integer(HID_T) :: memspace_id   ! memspace identifier 
  integer(HSIZE_T), allocatable :: dims(:), cnt(:), offset(:)
  integer(HSIZE_T) :: hcntm(1)    ! count for memory dataspace

  real(DP), allocatable :: array_real(:,:,:,:,:,:), array_imag(:,:,:,:,:,:)
  complex(DPC), allocatable :: array_cplx(:,:,:,:,:,:)
  integer :: err
  integer :: ip
  integer :: ik, jk, ib, jb, mb, nb
  integer :: i, j
  real(DP) :: kqpt(3)
  integer :: idx_kq
  integer :: my_nk, my_nkp
  logical :: read_slice

  PUSH_SUB(load_bse_kernel)

  call timing%start(timing%input)

  write(iq_str,*) iq
  kernel_fname = 'bse_'//trim(adjustl(iq_str))//'/kernel.h5'
  bsemat_fname = 'bse_'//trim(adjustl(iq_str))//'/bsemat.h5'

  ! hdf5 I/O example:
  !   https://www.asc.ohio-state.edu/wilkins.5/computing/HDF/hdf5tutorial/examples/F90/rwdsetexample.f90

  tdevol%iq = iq   ! label this tdevol
  tdevol%qpt(:) = dyn%qpts(:,iq)
  ! TODO: modify file path w.r.t. iq

  !====== read bsemat.h5 for header information ======
  ! root node reads and then broadcast
  if (peinf%inode == 0) then
    write(*,*) 'Reading bsemat.h5 header'

    SAFE_ALLOCATE(dims, (1))
    dims(1) = 1

    ! initialize fortran interface
    call h5open_f(err)
    ! open existing file
    call h5fopen_f(bsemat_fname, H5F_ACC_RDONLY_F, file_id, err)

    ! read nk
    ! open existing dataset
    call h5dopen_f(file_id, 'bse_header/kpoints/nk', dset_id, err)
    ! read dataset
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, tdevol%nk, dims, err)
    ! close dataset
    call h5dclose_f(dset_id, err)

    ! read ncb
    call h5dopen_f(file_id, 'bse_header/bands/ncb', dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, tdevol%ncb, dims, err)
    call h5dclose_f(dset_id, err)

    ! read nvb
    call h5dopen_f(file_id, 'bse_header/bands/nvb', dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, tdevol%nvb, dims, err)
    call h5dclose_f(dset_id, err)

    ! read nspin
    call h5dopen_f(file_id, 'bse_header/bands/ns', dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, tdevol%nspin, dims, err)
    call h5dclose_f(dset_id, err)

    ! read nspinor
    call h5dopen_f(file_id, 'bse_header/bands/nspinor', dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, tdevol%nspinor, dims, err)
    call h5dclose_f(dset_id, err)

    ! read celvol
    call h5dopen_f(file_id, 'mf_header/crystal/celvol', dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tdevol%celvol, dims, err)
    call h5dclose_f(dset_id, err)

    ! read blat
    call h5dopen_f(file_id, 'mf_header/crystal/blat', dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tdevol%blat, dims, err)
    call h5dclose_f(dset_id, err)

    SAFE_DEALLOCATE(dims)
    SAFE_ALLOCATE(dims, (2))
    dims = (/3, 3/)

    ! read bvec
    call h5dopen_f(file_id, 'mf_header/crystal/bvec', dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tdevol%bvec, dims, err)
    call h5dclose_f(dset_id, err)

    dims = (/3, tdevol%nk/)
    SAFE_ALLOCATE(tdevol%kpts, (3, tdevol%nk))
    SAFE_ALLOCATE(tdevol%kqpts, (3, tdevol%nk))
    SAFE_ALLOCATE(tdevol%kqmap, (tdevol%nk))

    ! read kpts
    call h5dopen_f(file_id, 'bse_header/kpoints/kpts', dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tdevol%kpts, dims, err)
    call h5dclose_f(dset_id, err)

    ! close file
    call h5fclose_f(file_id, err)
    ! close interface
    call h5close_f(err)

    ! find k+q points and kqmap
    ! find k+q points and k+q maps to kpoints
    do ik = 1, tdevol%nk
      kqpt(:) = tdevol%kpts(:,ik) + tdevol%qpt(:)
      tdevol%kqpts(:,ik) = kqpt(:)
      ! index of k+q point
      idx_kq = kindex(kqpt, tdevol%kpts, tdevol%nk)
      tdevol%kqmap(ik) = idx_kq
    enddo

    SAFE_DEALLOCATE(dims)
  endif

  ! broadcast
  call MPI_Bcast(tdevol%nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(tdevol%ncb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(tdevol%nvb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(tdevol%nspin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(tdevol%nspinor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(tdevol%celvol, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(tdevol%blat, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(tdevol%bvec, 3*3, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)

  if (peinf%inode .ne. 0) then
    SAFE_ALLOCATE(tdevol%kpts, (3, tdevol%nk))
    SAFE_ALLOCATE(tdevol%kqpts, (3, tdevol%nk))
    SAFE_ALLOCATE(tdevol%kqmap, (tdevol%nk))
  endif
  call MPI_Bcast(tdevol%kpts, 3*tdevol%nk, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(tdevol%kqpts, 3*tdevol%nk, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(tdevol%kqmap, tdevol%nk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

  ! end reading bsemat.h5 header information

  !====== set up parallelization information ======
  call init_parallelization(dyn%n_finite_q, tdevol)
  ! after this call, the distributed k_start, k_end, kp_start, kp_end are fixed
  ! below we read kernel.h5 for slices only relevant to the core`s own part in interaction
  my_nk  = tdevol%mpi_k_end  - tdevol%mpi_k_start  + 1
  my_nkp = tdevol%mpi_kp_end - tdevol%mpi_kp_start + 1

  !====== read kernel.h5 for interaction matrix ======
  ! read slices, data distribution at read in
  tdevol%nb = tdevol%nvb + tdevol%ncb  ! TODO now this is for insulator only

  ! allocate array
  ! serial version: SAFE_ALLOCATE(tdevol%vwint, (tdevol%nb*tdevol%nb*tdevol%nk, tdevol%nb*tdevol%nb*tdevol%nk))
  if ((my_nk .gt. 0) .and. (my_nkp .gt. 0)) then
    read_slice = .true.
  else
    read_slice = .false.
  endif

  if (read_slice) then
    SAFE_ALLOCATE(tdevol%vwint, (tdevol%nb*tdevol%nb*my_nk, tdevol%nb*tdevol%nb*my_nkp))
  endif

  do ip=1, peinf%npes
    ! we ask every processor to read sequentially to avoid conflict in hdf5 I/O use MPI_Barrier
    !call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    if ((peinf%inode .eq. ip-1) .and. read_slice) then
      if (peinf%inode .eq. 0) then
        write(*,*) 'Reading kernel.h5 from ', trim(kernel_fname)
        write(*,*) 'Data of kernel matrix elements distributed at reading'
      endif

      ! initialize fortran interface
      call h5open_f(err)
      ! open existing file
      call h5fopen_f(kernel_fname, H5F_ACC_RDONLY_F, file_id, err)
      ! open existing dataset
      call h5dopen_f(file_id, 'kernel', dset_id, err)
      ! get a copy of dataspace
      call h5dget_space_f(dset_id, dataspace_id, err)

      ! all nodes read full data, and later can be modified into distributed
      !serial version: SAFE_ALLOCATE(array_real, (tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nk, tdevol%nk))
      !serial version: SAFE_ALLOCATE(array_imag, (tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nk, tdevol%nk))
      !serial version: SAFE_ALLOCATE(array_cplx, (tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nk, tdevol%nk))

      SAFE_ALLOCATE(array_real, (tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nb, my_nk, my_nkp))
      SAFE_ALLOCATE(array_imag, (tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nb, my_nk, my_nkp))
      SAFE_ALLOCATE(array_cplx, (tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nb, my_nk, my_nkp))

      ! serial version: cnt = (/1, tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nk, tdevol%nk/)
      cnt = (/1, tdevol%nb, tdevol%nb, tdevol%nb, tdevol%nb, my_nk, my_nkp/)
      hcntm = product(cnt)

      ! read real part
      ! serial version: offset = (/0, 0, 0, 0, 0, 0, 0/)   ! for real part, first index is 0
      offset = (/0, 0, 0, 0, 0, tdevol%mpi_k_start-1, tdevol%mpi_kp_start-1/)
      call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, cnt, err)
      call h5screate_simple_f(1, hcntm, memspace_id, err)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array_real, hcntm, err, memspace_id, dataspace_id)

      ! read imaginary part
      ! serial version: offset = (/1, 0, 0, 0, 0, 0, 0/)   ! for imaginary part, first index is 1
      offset = (/1, 0, 0, 0, 0, tdevol%mpi_k_start-1, tdevol%mpi_kp_start-1/)
      call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, cnt, err)
      call h5screate_simple_f(1, hcntm, memspace_id, err)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array_imag, hcntm, err, memspace_id, dataspace_id)

      array_cplx(:,:,:,:,:,:) = cmplx(array_real(:,:,:,:,:,:), array_imag(:,:,:,:,:,:))

      !write(*,*) 'chk array_cplx'
      !write(*,*) array_cplx(1,1,1,1,1,1)

      SAFE_DEALLOCATE(array_real)
      SAFE_DEALLOCATE(array_imag)

      ! assign input data to interaction matrix vwint
      tdevol%vwint(:,:) = (0.0d0, 0.0d0)

      ! serial version:
      !do ik = 1, tdevol%nk  ! k
      !  do jk = 1, tdevol%nk  ! k`
      !
      ! here we only need slices
      do ik = 1, my_nk  ! k
        do jk = 1, my_nkp  ! k`
          ! each (k, k`) reads full band matrix
          do ib = 1, tdevol%nb  ! c
            do jb = 1, tdevol%nb  ! v
              do mb = 1, tdevol%nb  ! c`
                do nb = 1, tdevol%nb  ! v`
                  i = (ik - 1) * tdevol%nb * tdevol%nb + (ib - 1) * tdevol%nb + jb  ! kcv
                  j = (jk - 1) * tdevol%nb * tdevol%nb + (mb - 1) * tdevol%nb + nb  ! k`c`v`

                  ! format of array_cplx from kernel.py, in Fortran index:
                  !                              (v`, c`, v,  c,  k,  k`)
                  tdevol%vwint(i, j) = array_cplx(nb, mb, jb, ib, ik, jk)
                enddo ! nb
              enddo ! mb
            enddo ! jb
          enddo ! ib
        enddo ! jk
      enddo ! ik

      SAFE_DEALLOCATE(array_cplx)

    endif ! inode == ip and read_slice

  enddo ! ip

  call timing%stop(timing%input)

  POP_SUB(load_bse_kernel)

end subroutine load_bse_kernel

!=========================================================================
!> load eqp.dat and vmtxel (binary) file
subroutine load_eqp_vmtxel(dyn, tdevol)

  type(dynaminfo), intent(in) :: dyn          ! input parameters
  type(td_equation), intent(inout) :: tdevol  ! time-dependent evolution object

  character*256 :: eqp_fname                  ! eqp file containing dft and gw eigenvalues 
  character*256 :: vmtxel_A_fname             ! vmtxel.A file, easily switchable to binary
  character*256 :: vmtxel_B_fname             ! vmtxel.B file
                                              ! here we need two file because for finite-Q
                                              ! the two blocks are NOT complex conjugate
                                              ! there is finite Q shift involved
                                              ! definition should be consistent with extended_kernel
                                              ! vmtxel.A is the vmtxel from Q calculation
                                              ! vmtxel.B is the vmtxel from -Q calculation
  character*256 :: iq_str

  integer :: ik, ib, jb, is, ib_read, nb_read
  integer :: nk_read, ncb_read, nvb_read, ns_read      ! for vmtxel.A
  integer :: nk_read_, ncb_read_, nvb_read_, ns_read_  ! for vmtxel.B
  integer :: icb, jvb
  integer :: opr  !< 0 = use velocity operator, 1 = use momentum operator
  integer :: opr_

  real(DP) :: kpt(3), kdiff(3)
  real(DP) :: emf, eqp                        ! mean-field energy and quasiparticle (GW) energy
  integer :: bse_index_A, bse_index_B         ! index conversion
  complex(DPC), allocatable :: array_A(:)
  complex(DPC), allocatable :: array_B(:)
  real(DP) :: kmqpt(3)   ! k-q point
  integer :: ikmq
  complex(DPC) :: imag_one = (0.0d0, 1.0d0)

  PUSH_SUB(load_eqp_vmtxel)

  call timing%start(timing%input)

  write(iq_str,*) tdevol%iq
  eqp_fname = 'bse_'//trim(adjustl(iq_str))//'/eqp.dat'
  vmtxel_A_fname = 'bse_'//trim(adjustl(iq_str))//'/vmtxel.A'
  vmtxel_B_fname = 'bse_'//trim(adjustl(iq_str))//'/vmtxel.B'

  !====== read eqp.dat for quasiparticle band energies ======
  call open_file(11, file=eqp_fname, form='formatted', status='old')

  if (peinf%inode .eq. 0) write(*,*) 'Reading eqp.dat from ', trim(eqp_fname)

  SAFE_ALLOCATE(tdevol%band_energy, (tdevol%nb, tdevol%nk))

  if (dyn%band_index_max - dyn%band_index_min + 1 .ne. tdevol%nb) then
    call die('Number of bands mismatch from dynamics.inp and bsemat.h5', only_root_writes = .true.)
  endif

  if (tdevol%nspin == 2) then
    call die('nspin = 2 case not supported for eqp.dat yet.', only_root_writes = .true.)  ! TODO
  endif

  ! eqp.dat file required to contain the same set of kpoints as in bsemat.h5
  ! this is controlled automatically by the code below:
  !  - if eqp.dat contains less kpoints, then the read command will reach end-of-file since ik loops over tdevol%nk
  !  - if eqp.dat contains a kpoint (in the first tdevol%nk points) not in the tdevol%kpts list, kindex will initiate error
  !
  do ik = 1, tdevol%nk
    read(11,*) kpt, nb_read
    
    do ib = 1, nb_read
      read(11,*) is, ib_read, emf, eqp

      if ((ib_read .ge. dyn%band_index_min) .and. (ib_read .le. dyn%band_index_max)) then
        ! BerkeleyGW/BSE orders bands differently
        ! the valence/conduction bands are counted from Fermi level
        ! therefore, for valence bands, we need to invert the band ordering
        !
        if (ib_read .lt. dyn%band_index_min + tdevol%nvb) then
          ! valence bands, invert ordering
          tdevol%band_energy(tdevol%nvb - (ib_read - dyn%band_index_min), kindex(kpt, tdevol%kpts, tdevol%nk)) = eqp / ryd
        else
          ! conduction bands, keep ordering
          tdevol%band_energy(ib_read - dyn%band_index_min + 1, kindex(kpt, tdevol%kpts, tdevol%nk)) = eqp / ryd
        endif ! invert valence band index

      endif ! dyn%band_index_min <= ib_read <= dyn%band_index_max
    enddo ! ib_read
  enddo ! ik

  call close_file(11)

  !====== read vmtxel for velocity matrix elements ======
  call open_file(15, file=vmtxel_A_fname, form='unformatted', status='old')
  call open_file(16, file=vmtxel_B_fname, form='unformatted', status='old')

  if (peinf%inode .eq. 0) write(*,*) 'Reading vmtxel from ', trim(vmtxel_A_fname), ' and ', trim(vmtxel_B_fname)

  read(15) nk_read, ncb_read, nvb_read, ns_read, opr
  read(16) nk_read_, ncb_read_, nvb_read_, ns_read_, opr_

  if(peinf%inode == 0) then
    if (opr == 0) then
      write(*,'(1x,a)') 'Dipole matrix elements using velocity operator'
    elseif (opr == 1) then
      write(*,'(1x,a)') 'Dipole matrix elements using momentum operator'
    endif
  endif

  ! sanity check
  if (nk_read .ne. tdevol%nk) call die('Number of kpoints mismatch in vmtxel and bsemat.h5.')
  if (ncb_read .ne. tdevol%ncb) call die('Number of conduction bands mismatch in vmtxel and bsemat.h5.')
  if (nvb_read .ne. tdevol%nvb) call die('Number of valence bands mismatch in vmtxel and bsemat.h5.')
  if (ns_read .ne. tdevol%nspin) call die('Number of spins mismatch in vmtxel and bsemat.h5.')
  if ((opr .ne. 0) .and. (opr .ne. 1)) call die('Illegal dipole operator type in vmtxel.')

  if (ns_read == 2) then
    call die('nspin = 2 case not supported for vmtxel yet.', only_root_writes = .true.)  ! TODO
  endif

  if (abs(nk_read-nk_read_) + abs(ncb_read-ncb_read_) + abs(nvb_read-nvb_read_) + abs(ns_read-ns_read_) + abs(opr-opr_) .ne. 0) then
    call die('Header information inconsistent between vmtxel.A and vmtxel.B', only_root_writes = .true.)
  endif

  SAFE_ALLOCATE(array_A, (ncb_read * nvb_read * nk_read))
  SAFE_ALLOCATE(array_B, (ncb_read * nvb_read * nk_read))
  SAFE_ALLOCATE(tdevol%dipole, (tdevol%nb, tdevol%nb, tdevol%nk))
  tdevol%dipole(:,:,:) = (0.0d0,0.0d0)

  read(15) array_A(:)
  read(16) array_B(:)

  do ik = 1, nk_read
    do icb = 1, ncb_read
      ib = tdevol%nvb + icb

      do jvb = 1, nvb_read
        jb = jvb

        ! convert to one dimensional index
        !   following ../Common/misc.f90
        bse_index_A = jvb + (icb - 1 + (ik - 1)*ncb_read)*nvb_read

        ! first complete the A part, (ck, vk+q)
        ! multiply by imaginary one
        tdevol%dipole(ib, jb, ik) = imag_one * array_A(bse_index_A)

        ! then complete the B part, (vk, ck+q) <=> (vk-q, ck)
        ! we take these elements from the -Q calculation
        ! now need to find the correct kpoint mapping of k-q point
        kmqpt(:) = tdevol%kpts(:,ik) - tdevol%qpt(:)
        ikmq = kindex(kmqpt, tdevol%kpts, tdevol%nk)

        ! the current ik in vmtxel_b corresponds to ikmq label in tdevol
        bse_index_B = jvb + (icb - 1 + (ik - 1)*ncb_read)*nvb_read

        ! then we do complex conjugate and transpose
        tdevol%dipole(jb, ib, ikmq) = MYCONJG(imag_one * array_B(bse_index_B))

      enddo ! jvb
    enddo ! icb
  enddo ! ik

  SAFE_DEALLOCATE(array_A)
  SAFE_DEALLOCATE(array_B)

  call close_file(15)
  call close_file(16)

  call timing%stop(timing%input)

  POP_SUB(load_eqp_vmtxel)

end subroutine load_eqp_vmtxel
!=========================================================================
!> load electron-phonon matrix elements 
subroutine load_eph_matrix()

  PUSH_SUB(load_eph_matrix)

  call timing%start(timing%input)

  call timing%stop(timing%input)

  POP_SUB(load_eph_matrix)

end subroutine load_eph_matrix

#endif

end module input_dynamics_m
