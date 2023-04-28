!==============================================================================
!
! Routines:
!
! (1) intkernel()       Originally By MLT        Last Modified: 6/16/2011 (FHJ)
!
!     input: crys, kg, kp, syms, xct, peinf types
!            hbse_a     effective Hamiltonian, with only the diagonal part
!            dcc, dvv  transformation matrices
!            kco       coordinates of k-points in the coarse grid
!            fi2co_wfn map of k-points fine grid -> coarse grid;
!                      defines which point in the coarse grid is close to
!                      one in the fine grid
!     output: hbse_a    effective Hamiltonian, with diagonal part and
!                      interaction part (direct + exchange)
!
!     Build the interaction Kernel using the interpolation scheme
!     (see Rohlfing & Louie). The kernel on the coarse grid is read
!     from unit, "bsedmat" and "bsexmat", and used immediately
!     (using temporary files may speed up the calculation...)
!
!     interaction kernel K =  K_d + 2*K_x   <---- spin singlet
!                        K =  K_d           <---- spin triplet
!                        K =  2*K_x         <---- local fields
!                        K =  K_d + K_x     <---- spinor
!
!     imatrix code: 1 dhead, 2 dwing, 3 dbody, 4 x (BSE)
!     imatrix code: 1 dhead, 2 dbody, 3 x, 4 fxc (TDDFT)
!
! (2) interpolate()     Originally By MLT        Last Modified: 7/1/2008 (JRD)
!
!     input: xct type
!            bse_co    kernel in the coarse grid
!            dcckp     dcc transf. matrix at point k`
!            dvvkp     dvv transf. matrix at point k`
!     output: bse_fi   (interpolated) kernel in the fine grid
!
!==============================================================================

#include "f_defs.h"

module intkernel_m

#ifdef HDF5
  use hdf5
#endif
  use hdf5_io_m
  use global_m
  use blas_m
  use checkgriduniformity_m
  use epsdiag_m
  use intpts_m
  use interp_m
  use misc_m
  use vcoul_generator_m
  use cells_m
  use io_utils_m
  use wfn_rho_vxc_io_m
  use kernel_io_m
  use timing_m, only: timing => bse_timing
  implicit none

  private

  public :: intkernel

contains

subroutine intkernel(crys,kg,kp,syms,xct,hbse_a,dcc,dvv,kco,fi2co_wfn,flag,gvec,intp_coefs,hbse_b,dcc_sub,dvv_sub,closepts_sub)
  type (crystal), intent(in) :: crys
  type (grid), intent(in) :: kg
  type (kpoints), intent(in) :: kp
  type (symmetry), intent(in) :: syms
  type (xctinfo), intent(inout) :: xct
  SCALAR, intent(inout) :: hbse_a(:,:) !< (xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin,peinf%nblocks*peinf%block_sz)
  !> (xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel) Band expansion coefficients
  SCALAR, intent(in) :: dcc(:,:,:,:,:), dvv(:,:,:,:,:)
  real(DP), intent(in) :: kco(:,:) !< (3,xct%nkpt_co)
  !> (xct%npts_intp_kernel, xct%nkpt_fi) For a given k-point in the fine grid
  !! returns the closest k-point in the coarse WFN k-grid.
  integer, intent(in) :: fi2co_wfn(:,:)
  type(flags), intent(in) :: flag
  type (gspace), intent(in) :: gvec
  !> (xct%npts_intp_kernel, xct%nkpt_fi) Delaunay/greedy interpolation coefficients
  real(DP), intent(in) :: intp_coefs(:,:)
  SCALAR, intent(inout), optional :: hbse_b(:,:)

  type (twork_scell) :: work_scell
  real(DP) :: vcoul, oneoverq
  real(DP) :: vcoul0(1), closeweights(4), q0temp(3)
  SCALAR   :: epshead_temp
  type (epsinfo) :: epsi

! FHJ: these arrays depend only on the distance w.r.t qq=0 (inside the BZ)
  real(DP), allocatable :: dist_array(:) !length of dq for a given index
  real(DP), allocatable :: vcoul_array(:), oneoverq_array(:)

  !> FHJ: cells structure that keeps all fine (k-kp) transitions
 type(cells_t) :: cells_fi
  ! FHJ: these variables store the cached value of |q|^2, v(q), epsinv, for
  ! each dq in the same order as it appears in the cell structure.
  real(DP) :: dq(3), dq_bz(3), abs_q2
  real(DP), allocatable :: dqs(:,:), dqs_bz(:,:), abs_qs2(:), eps_intp(:)
  integer, allocatable :: vq_map(:)
  integer :: ik_cells, jkp_offset

  integer :: isrtrq(1), closepts(4)
  type(interp_t) :: eps_interp

  integer :: iscreentemp,iparallel,ijk,icb,ivb,imatrix, &
    ik,ic,iv,ikp,icp,ivp,iblock,ikcvs,ikcvsd,jj,jc,jv,js, &
    jcp,jvp,jsp,jk,jkp,dimbse,nmatrices,ifile, &
    i1,i2,i3,iold,icout,ivout,inew,tda_sz,ii
  !> (xct%nkpt_co) maps a k-point in the WFN_co to the bse*mat files.
  integer, allocatable :: wfn2bse(:)
  !> (xct%npts_intp_kernel, xct%nkpt_fi) For a given k-point in the fine grid,
  !! returns the closest k-point in the (coarse) bse*mat files.
  integer, allocatable :: fi2co_bse(:,:)
  integer :: ivertp, ivert
  integer, allocatable :: jkp2offset(:)

  real(DP) :: bsemat_fac,fac_d,fac_x,w_eff,eps,tol,factor,fac_t
  type(progress_info) :: prog_info
  type(kernel_header_t) :: kern, kern2, kern3, kern_sub

#ifdef HDF5
  ! JIM: variables used for HDF5 interface
  integer(HID_T) :: file_id
  integer(HID_T) :: file_id_sub
  integer :: error
#endif

  ! DYQ: variables used for interpolation with subsample_line flag
  !SCALAR, allocatable :: hbse_a_temp(:,:)
  ! subsampled bsemat: bsedmt_sub(iv,ic,ivp,icp,ik_co,ik_sub,imatrix)
  SCALAR, allocatable :: bsedmt_sub(:,:,:,:,:,:,:),bsedmt_sub_temp(:,:,:,:,:,:)
  real(DP), allocatable :: kpoints_sub(:,:,:),kpoints_sub_len(:)  ! k-points and lengths read from subsampled bsemat files
  real(DP) :: k_temp(3),k_temp_bz(3),abs_ksub2,dist_min,qq(3)
  !integer, allocatable :: sub2bse(:),bse2sub(:) ! Maps kpt indices between subsampled bsemat files and coarse bsemat file
  integer :: nk_sub,nsub_files,isub,jk_sub,jkp_sub
  character(len=80) :: fname_sub, tmpstr
  SCALAR,intent(inout), optional ::dcc_sub(:,:,:,:,:,:),dvv_sub(:,:,:,:,:,:)
  integer,intent(in),optional :: closepts_sub(:,:)
  integer :: ivert_sub,ivertp_sub,jk_temp,jkp_temp,ivert_temp,ivertp_temp
  SCALAR,allocatable :: dcckp_sub(:,:,:),dvvkp_sub(:,:,:)
  logical :: found_vertex, old_csi_behavior


!------------------------------
! kernel and transformation matrices

  SCALAR, allocatable :: &
    bsedmatrix(:,:,:,:,:,:,:),bsedmt(:,:,:,:,:,:,:), &
    bsedmatrix_loc(:,:,:,:,:,:,:), &
    dcckp(:,:,:), dvvkp(:,:,:)

  if(abs(xct%scaling) < TOL_Zero) return ! everything will be zero

  PUSH_SUB(intkernel)
  !
  ! CSO: setting to .true. will reproduce old behavior (and e.g., absorption
  ! eigenvalues), which selects the shared vertex that is first identified by
  ! the algorithm. Setting to false selects the shared vertex that is closest to
  ! the fine k-point. The latter breaks the testsuite because it becomes
  ! possible to select different vertices of the same simplex which are just as
  ! close to the fine k-point, i.e. degenerate. When vertices are degenerate,
  ! different compilers may select different vertices. A testsuite (that insists
  ! on high precision) may break because using different (but degenerate)
  ! vertices leads to slightly different (but non-degenerate) interpolated
  ! kernels.
  !
  if (xct%subsample_algo==0) then
    old_csi_behavior = .true.
  else if (xct%subsample_algo==1) then
    old_csi_behavior = .false.
  else
    call die("Illegal input for subsample_algo.")
  endif

  if (.not.xct%tda.and..not.present(hbse_b)) &
    call die('Doing non-TDA but array hbse_b not present.', only_root_writes=.true.)

  if (xct%theory==0) then
!---------------------------------
! Read head of dielectric function from file (either 'epsdiag.dat' or
! 'eps0mat'/'epsmat')

    call logit('Calling epsdiag')

    call timing%start(timing%epsdiag)
    if(peinf%inode==0) then
      if(flag%read_epsdiag) then
        call read_epsdiag(epsi)
      else
        call epsdiag(crys,gvec,syms,epsi,xct)
      endif
    endif
    call timing%stop(timing%epsdiag)

    call timing%start(timing%eps_comm)
#ifdef MPI
    call MPI_BCAST(epsi%nq,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(epsi%q0vec,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(epsi%epshead,1,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(epsi%emax,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
    if(peinf%inode.ne.0) then
      SAFE_ALLOCATE(epsi%eps, (epsi%nq))
    endif
    call MPI_BCAST(epsi%eps,epsi%nq,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
    if(peinf%inode.ne.0) then
      SAFE_ALLOCATE(epsi%q, (3,epsi%nq))
    endif
    call MPI_BCAST(epsi%q,3*epsi%nq,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif
    call timing%stop(timing%eps_comm)
  endif

  SAFE_ALLOCATE(wfn2bse, (xct%nkpt_co))
  SAFE_ALLOCATE(fi2co_bse, (xct%npts_intp_kernel, xct%nkpt_fi))

  tol=TOL_Small
  factor = -8.d0*PI_D/(crys%celvol*xct%nktotal)*xct%scaling
  if (flag%krnl.eq.0) then
    if(peinf%inode.eq.0) write(6,'(/1x,a)') 'Assembling spin-triplet kernel.'
    fac_d = factor
    fac_x = 0d0
    fac_t = factor
  elseif (flag%krnl.eq.1) then
    if(peinf%inode.eq.0) write(6,'(/1x,a)') 'Assembling spin-singlet kernel.'
    fac_d = factor
    fac_x = factor
    fac_t = factor
  elseif (flag%krnl.eq.2) then
    if(peinf%inode.eq.0) write(6,'(/1x,a)') 'Assembling local-fields kernel.'
    fac_d = 0d0
    fac_x = factor
    fac_t = 0d0
  elseif (flag%krnl.eq.3) then
    if(peinf%inode.eq.0) write(6,'(/1x,a)') 'Assembling spinor kernel.'
    fac_d = factor
    fac_x = factor
    fac_t = factor
    if(xct%nspinor.ne.2) then
      call die("Spinor kernel requires nspinor=2.")
    endif
  else
    write(0,*) 'kernel type = ', flag%krnl
    call die("Illegal value of kernel type in intkernel.")
  endif

!-----------------------------------------------------------------------------!
!              EPSILON INTERPOLATION AND V(Q) PRE-CALCULATION
! FHJ: The interpolation of epsinv scales as O(nkpt) for regular grids.
! We pre-calculate all dqs = k(:,:) - k(:,1), move this quantity to the [0, 1)
! range, and put them into a cell structure. For each dq, we interpolate epsinv.
! In a separate memory copy, we move all dqs to the 1st BZ (dqs_bz), and
! calculate v(dq_bz) and related quantities.
! When we reach the actual kernel interpolation, we calculate dq=k-kp, move it
! to the same [0, 1) range, and get the mapping to the dqs, which is done in
! O(1) using the cell structure.
! This algorithm should also works for coarse-grained grids (a fine grid with
! holes), as long as the dqs are compute with all possible fine transitions.


! MJ: Only if we have an epsilon
  if (xct%theory .eq. 0) then
! FHJ: Prepare the cells for the interpolation of epsilon
    if (xct%delaunay_interp) then
      call interp_init(eps_interp, crys, epsi%q, dims=xct%idimensions, &
        periodic=.false., active_dims=xct%is_periodic)
    else
      call alloc_intpts(epsi%nq, epsi%q(:,:), periodic=.false.)
    endif

    SAFE_ALLOCATE(eps_intp, (xct%nktotal))
  endif

! MJ : All of these are needed for V(Q) as well
  SAFE_ALLOCATE(dist_array, (xct%nktotal))
  SAFE_ALLOCATE(vcoul_array, (xct%nktotal))
  SAFE_ALLOCATE(oneoverq_array, (xct%nktotal))
  SAFE_ALLOCATE(dqs, (3, xct%nktotal))
  SAFE_ALLOCATE(dqs_bz, (3, xct%nktotal))
  SAFE_ALLOCATE(abs_qs2, (xct%nktotal))
  SAFE_ALLOCATE(vq_map, (xct%nktotal))

  ! FHJ: This block is equivalent to dqs(:,ik) = UNIT_RANGE(kg%f(:,ik)-kg%f(:,1))
  ! if nkpt_fi==nktotal, but only up to TOL_SMALL. In order for the testsuite not to
  ! fail, we explicitly break this into standard and patched_sampling cases.
  if (xct%patched_sampling) then
    ik=0
    do i1=0,kp%kgrid(1)-1
      do i2=0,kp%kgrid(2)-1
        do i3=0,kp%kgrid(3)-1
          ik=ik+1
          ! FHJ: This is already in the [0,1) range
          dqs(1,ik)=(dble(i1))/dble(kp%kgrid(1))
          dqs(2,ik)=(dble(i2))/dble(kp%kgrid(2))
          dqs(3,ik)=(dble(i3))/dble(kp%kgrid(3))
        enddo
      enddo
    enddo
  else
    do ik = 1, xct%nkpt_fi
      dqs(:,ik) = UNIT_RANGE(kg%f(:,ik) - kg%f(:,1))
    enddo
  endif

  ! FHJ: some debugging stuff.
  !if (peinf%inode==0) then
  !  write(601,'(a)') '#iq  dq  dq_bz  abs_q2'
  !  write(602,'(a)') '#iq vcoul'
  !  write(603,'(a)') '#iq eps'
  !endif
  ! FHJ: TODO: parallelize this loop (trivial)
  call timing%start(timing%ik_inteps)
  call cells_init(cells_fi, dqs, periodic=.true.)
  do ik = 1, xct%nktotal
    dq(:) = dqs(:, ik)

    ! FHJ: Move pt to 1st BZ and save it
    call point_to_bz(crys, dq, dq_bz, abs_q2)
    dqs_bz(:, ik) = dq_bz(:)
    abs_qs2(ik) = abs_q2
    if (xct%theory .eq. 0) then
      if (abs_q2 > epsi%emax) write(0,*) ' WARNING: emax too small, ', abs_q2, epsi%emax

      ! Interpolate epsinv - Just the head
      if (xct%icutv/=TRUNC_NONE .and. xct%iscreen==SCREEN_SEMICOND .and. abs_q2<TOL_ZERO) then
        ! SIB:  if q=0 and we are truncating the coulomb interaction, then
        !       the dielectric is 1 at q=0 for semiconductors.
        eps = 1.0d0
      else
        if (xct%delaunay_interp) then
          call interp_eval(eps_interp, dq_bz, closepts(:), closeweights(:))
        else
          call intpts_local(crys, dq_bz, epsi%nq, epsi%q(:,:), xct, closepts(:), &
            closeweights(:), periodic=.false.)
        endif
        eps=0D0
        do ijk = 1, xct%idimensions + 1
          eps = eps + closeweights(ijk)*epsi%eps(closepts(ijk))
        enddo
      endif
      eps_intp(ik) = eps
      ! FHJ: some debugging stuff.
      !if (peinf%inode==0) then
      !  write(603,'(1x,i6,1x,es25.16)') ik, dble(eps)
      !endif
    endif ! theory=0
  enddo !ik
  if (xct%theory .eq. 0) then
    if (xct%delaunay_interp) then
      call interp_free(eps_interp)
    else
      call dealloc_intpts()
    endif
  endif ! theory=0
  call timing%stop(timing%ik_inteps)
! FHJ: Done with epsinv interpolation
!-----------------------------------

  if (xct%theory .eq. 0) then
    SAFE_DEALLOCATE_P(epsi%q)
    SAFE_DEALLOCATE_P(epsi%eps)
  endif

! FHJ: v(q) pre-calculation
!-------------------------------
! We should all contribute to calculating vcoul0 since it might involve minibz average

  call timing%start(timing%ik_vcoul)
  if (peinf%inode .eq. 0) then
    write(6,'(1x,a,i0,a)') 'Calculating Coulomb potential v(q) for ', xct%nktotal, ' q-points.'
  endif
  iscreentemp = 0
  iparallel = 1
  inew = 0

  if(peinf%inode == 0) call checkgriduniformity(kp%kgrid, crys, xct%icutv)

  if (peinf%verb_high .and. peinf%inode==0) then
    write(6,'(a)')
    write(6,1000) 'Point','kpt','dq  = kg%f(:,ik)-kg%f(:,1)','   |dq|^2   ','index'
    write(6,1000) '-----','---','--------------------------','------------','-----'
1000  format(2x, a5,1x, a5,1x, a26,1x, a12,1x, a5)
1001  format(2x, a5,1x, i5,1x, 3(f8.5,1x), f12.9,1x, i5)
  endif

  do ik = 1, xct%nktotal
    dq(:) = dqs(:, ik)
    dq_bz(:) = dqs_bz(:, ik)
    abs_q2 = abs_qs2(ik)

    ! FHJ: Take a look if v(q) for this |q| has already been calculated
    ! We assume that v(q) = v(|q|), which is fine since we only include q`s
    ! in the periodic directions, so truncation doesn`t affect.
    iold=0
    do jj=inew,1,-1
      if (dabs(dist_array(jj)-abs_q2) < TOL_Small) then
        iold=jj
        exit
      endif
    enddo

    if (iold==0) then
      inew=inew+1
      dist_array(inew)=abs_q2
      vcoul0(1)=0.0d0
      if (peinf%verb_high .and. peinf%inode==0) then
        write(6,1001) '~NEW~', ik, dq_bz, abs_q2, inew
      endif
      isrtrq(1) = 1
      if (xct%theory .eq. 0) then
        q0temp(:) = epsi%q0vec(:)
        epshead_temp = epsi%epshead
      else if (xct%theory .eq. 1) then
        q0temp(:) = xct%q0vec(:)
        epshead_temp = (1.0D0, 0.0d0)
        if ((xct%theory == 1) .and. ((1.0d0 - xct%coulomb_mod%long_range_frac_fock > TOL_SMALL) .or. &
           (1.0d0 - xct%coulomb_mod%short_range_frac_fock > TOL_SMALL))) then
            xct%coul_mod_flag=.true.
        endif
      endif

      if (.not. xct%coul_mod_flag) then
        call vcoul_generator(xct%icutv,xct%truncval,gvec, &
          crys%bdot,crys%celvol,xct%nktotal,1,isrtrq(:),iscreentemp,dq_bz,q0temp, &
          vcoul0(:),xct%iwritecoul,iparallel,xct%avgcut,oneoverq, &
          kp%kgrid,epshead_temp,work_scell,xct%averagew,xct%wcoul0)
      else
        call vcoul_generator(xct%icutv,xct%truncval,gvec, &
          crys%bdot,crys%celvol,xct%nktotal,1,isrtrq(:),iscreentemp,dq_bz,q0temp, &
          vcoul0(:),xct%iwritecoul,iparallel,xct%avgcut,oneoverq, &
          kp%kgrid,epshead_temp,work_scell,xct%averagew,xct%wcoul0, &
          coulomb_mod=xct%coulomb_mod)
      endif

      vcoul_array(inew) = vcoul0(1)/(8.0*Pi_D)
      oneoverq_array(inew) = oneoverq/(8.0*Pi_D)
      vq_map(ik) = inew
    else
      if (peinf%verb_high .and. peinf%inode==0) then
        write(6,1001) ' old ', ik, dq_bz, abs_q2, iold
      endif
      vq_map(ik) = iold
    endif !iold == 0
    ! FHJ: some debugging stuff.
    !if (peinf%inode==0) then
    !  write(601,'(1x,i6,1x,2(3(f12.9,1x),1x),es25.16)') ik, dq, dq_bz, abs_q2
    !  write(602,'(1x,i6,1x,es25.16)') ik, dble(vcoul_array(vq_map(ik)))
    !endif
  enddo !ik
  ! FHJ: some debugging stuff.
  !if (peinf%inode==0) then
  !  write(600,'("wcoul0=",es25.16)') dble(xct%wcoul0)
  !  write(600,'("factor=",es25.16)') factor
  !  write(600,'("celvol=",es25.16)') crys%celvol
  !endif

! FHJ: Done with v(q) pre-calculation
!-----------------------------------

  call destroy_qran()

  if (peinf%verb_high .and. peinf%inode==0) then
    write(6,'(a)')
    write(6,'(1x,a,i0,a)') 'Finished calculating Vcoul with ',inew,' unique points'
  endif

  if (xct%icutv==TRUNC_SUPERCELL) then
    SAFE_DEALLOCATE_P(work_scell%fftbox_1D)
  endif

  call timing%stop(timing%ik_vcoul)
  call timing%start(timing%ik_setup)


!-----------------------------------------------------------------------------
! DYQ: If using interpolation with subsample_line read the bsemat files containing the subsampled points
!

  if (xct%subsample_line) then
    if (peinf%inode==0) then
      write(6,'(/,1x,a)') 'Reading subsampled matrices'
      ! Read list of file names and nk_sub here
      call open_file(unit=77, file='subsample.inp',form='formatted',status='old')
      read(77,*) nk_sub
      write(6,'(1x,a,i6)') "Number of k-points per subsampled bsemat file = ",nk_sub
      read(77,*) nsub_files
      write(6,'(1x,a,i6)') "Number of bsemat files = ",nsub_files
      SAFE_ALLOCATE(bsedmt_sub, (xct%n1b_co,xct%n2b_co,xct%n1b_co,xct%n2b_co,xct%nkpt_co,nk_sub,4))
      bsedmt_sub=0.0
      SAFE_ALLOCATE(kpoints_sub,(3,nk_sub,xct%nkpt_co))
      SAFE_ALLOCATE(kpoints_sub_len,(nk_sub))
      ! skip reading coarse k-points
      do ii=1,nsub_files
        read(77,*)
      enddo
      do ii=1,nsub_files
        read(77,'(a80)') fname_sub
        write(6,'(1x,a,a)') "Opening ", fname_sub
#ifdef HDF5
        ! Read kpoints. Don`t read the mf header, because we don`t use it.
        call read_kernel_header_hdf5(fname_sub, kern_sub, rw_mf_header=.false.)
        kpoints_sub(:,:,ii) = kern_sub%kpts
        ! FHJ: TODO - deallocate kern_sub
        !
        if (ii.eq.1) then
          do jj=1,nk_sub
            k_temp = kpoints_sub(:,1,ii) - kpoints_sub(:,jj,ii)
            call point_to_bz(crys, k_temp, k_temp_bz, abs_ksub2)
            kpoints_sub_len(jj) = sqrt(abs_ksub2)
          enddo
        else
          ! CSO: self consistency check between different bsemats.h5
          do jj=1,nk_sub
            k_temp = kpoints_sub(:,1,ii) - kpoints_sub(:,jj,ii)
            call point_to_bz(crys, k_temp, k_temp_bz, abs_ksub2)
            if ( abs( sqrt(abs_ksub2) - kpoints_sub_len(jj)) .ge. TOL_ZERO ) then
              call die ( "Why are lengths of subsampled points from the coarse point are different for different coarse points?")
            endif
          enddo
        endif
        !
        ! FHJ: TODO - support collective read
        call hdf5_open_file(trim(fname_sub), 'r', file_id_sub)
        ! Read Matrices
        do imatrix=1,4
          call bse_hdf5_read(xct,nk_sub,file_id_sub,imatrix,ii,bsedmt_sub(:,:,:,:,:,:,:),subsample=.true.)
          !do jj=1,nk_sub
          !  do jv=1,xct%nvb_co
          !    do jc=1,xct%ncb_co
          !      do jvp=1,xct%nvb_co
          !        do jcp=1,xct%ncb_co
          !          write(6,*) ii,jj,imatrix,jv,jc,jvp,jcp,bsedmt_sub(imatrix,jv,jc,jvp,jcp,ii,jj)
          !        enddo
          !      enddo
          !    enddo
          !  enddo
          !enddo
        enddo
        call hdf5_close_file(file_id_sub)
#endif
      enddo ! ik_co
      write(6,'(1x,a)') 'Finished reading subsampled bsemat files'
    endif

#ifdef MPI
    call MPI_BCAST(nk_sub, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
    if (peinf%inode.ne.0) then
      SAFE_ALLOCATE(bsedmt_sub, (xct%n1b_co,xct%n2b_co,xct%n1b_co,xct%n2b_co,xct%nkpt_co,nk_sub,4))
      SAFE_ALLOCATE(kpoints_sub,(3,nk_sub,xct%nkpt_co))
      SAFE_ALLOCATE(kpoints_sub_len,(nk_sub))
    endif
#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(bsedmt_sub(1,1,1,1,1,1,1), 4*xct%n1b_co*xct%n2b_co*xct%n1b_co*xct%n2b_co*xct%nkpt_co*nk_sub, &
      MPI_SCALAR, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kpoints_sub(1,1,1), xct%nkpt_co*3*nk_sub, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kpoints_sub_len(1), nk_sub, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
    SAFE_ALLOCATE(bsedmt_sub_temp,(xct%n1b_co,xct%n2b_co,xct%n1b_co,xct%n2b_co,xct%nspin,xct%nspin))
  endif


!--------------------------------
! Allocate data

  if (peinf%inode .eq. 0) then
    if (.not. xct%skipinterp) then
      write(6,'(1x,a)') 'Performing kernel interpolation.'
    else
      write(6,'(1x,a)') 'Building interaction kernel.'
    endif
  endif

!  if(peinf%inode.eq.0) write(6,'(a)') 'Data allocation: intkernel'
  tda_sz = 1
  if (.not.xct%tda) tda_sz = 2
  SAFE_ALLOCATE(bsedmt, (xct%nvb_fi,xct%ncb_fi,peinf%nv_block,peinf%nc_block,xct%nspin,xct%nspin,tda_sz))
  SAFE_ALLOCATE(dcckp, (xct%n2b_co,xct%ncb_fi,xct%nspin))
  SAFE_ALLOCATE(dvvkp, (xct%n1b_co,xct%nvb_fi,xct%nspin))

  if(xct%subsample_line) then
    SAFE_ALLOCATE(dcckp_sub, (xct%n2b_co,xct%ncb_fi,xct%nspin))
    SAFE_ALLOCATE(dvvkp_sub, (xct%n1b_co,xct%nvb_fi,xct%nspin))
  endif

  dimbse=xct%nkpt_co*(xct%n2b_co*xct%n1b_co*xct%nspin)**2
  wfn2bse=0

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

  call timing%stop(timing%ik_setup)
  call timing%start(timing%ik_c_check)

  if(.not. xct%use_hdf5) then
    if (peinf%inode==0) then
      call open_file(unit=11,file='bsedmat',form='unformatted',status='old')
      call open_file(unit=12,file='bsexmat',form='unformatted',status='old')
      if (xct%theory .eq. 1) call open_file(unit=13,file='bsetmat',form='unformatted',status='old')
    endif
    call read_binary_kernel_header(11, kern)
    call require_version('bsedmat', kern%mf%version, VER_BSE_FORT)
    call check_xctinfo_kernel_header('bsedmat', xct, kern)
    call read_binary_kernel_header(12, kern2)
    call require_version('bsexmat', kern2%mf%version, VER_BSE_FORT)
    call check_xctinfo_kernel_header('bsexmat', xct, kern2)

    if (.not.all(abs(kern%kpts(:,:)-kern2%kpts(:,:))<tol)) then
      call die('k-grids from bsedmat and bsexmat don`t match.', &
        only_root_writes=.true.)
    endif

    if (xct%theory .eq. 1) then
      call read_binary_kernel_header(13, kern3)
      call require_version('bsetmat', kern3%mf%version, VER_BSE_FORT)
      call check_xctinfo_kernel_header('bsetmat', xct, kern3)
      if (.not.all(abs(kern%kpts(:,:)-kern3%kpts(:,:))<tol)) then
        call die('k-grids from bsedmat and bsexmat don`t match.', &
          only_root_writes=.true.)
      endif
    endif

#ifdef HDF5
  else
    call read_kernel_header_hdf5('bsemat.h5', kern)
    !call check_xctinfo_kernel_header('bsemat.h5', xct, kern)
    if (peinf%inode==0) then
      !FHJ: TODO - change to collective I/O
      call hdf5_open_file('bsemat.h5', 'r', file_id)
    endif
#endif
  endif

  if (peinf%inode==0) then
    i1 = 0
    ! Check consistency of files
    ! wfn2bse(ikp): index of coarse-grid point ikp in files bsedmat,bsexmat
    do jj=1,xct%nkpt_co
      ikp = 0
      do ik=xct%nkpt_co,1,-1
        ! FHJ: it`s all right if the two k-grids differ by a reciprocal lattice vector
        if (all(abs(MIN_RANGE(kern%kpts(:,jj) - kco(:,ik))) < tol)) then
          ikp=ik
          exit
        endif
      enddo
      if (ikp==0) then
        write(0,'(a)') 'ERROR: found a k-point in the bsemat file that is not present in WFN_co.'
        write(0,'(a)') 'Are you using the exact same WFN_co from your kernel calculation?'
        write(0,'(a,i0)') 'jj = ', jj
        write(0,'(a,3(f10.6,1x))') 'k = ', kern%kpts(:,jj)
        call die('found a k-point in the bsemat file that is not present in WFN_co', &
          only_root_writes=.true.)
      endif
      wfn2bse(ikp)=jj
      if (jj/=ikp) i1=i1+1
    enddo

    if (i1.ne.0) then
      write(6,'(1x,a)')
      write(6,'(1x,a)') 'The k-points in the kernel file has a different order than in the'
      write(6,'(1x,a)') 'WFN file. The mapping between the k-points was updated.'
      write(6,'(1x,a)')
      if (peinf%verb_high) then
        do jk=1,xct%nkpt_co
          write(6,'(1x,2i8)') jk, wfn2bse(jk)
        enddo
        write(6,'()')
      endif
    endif
  endif !peinf%inode = 0

#ifdef MPI
  call MPI_BCAST(wfn2bse, xct%nkpt_co, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
  do ik = 1, xct%nkpt_fi
    fi2co_bse(:,ik) = wfn2bse(fi2co_wfn(:,ik))
  enddo
  SAFE_DEALLOCATE(wfn2bse)

  call timing%stop(timing%ik_c_check)
  call logit('Done calculating Vcoul')

!------------- Read in coarse matrices: head, wings, body, exchange --------------------
! PE # 0 reads and broadcasts to everybody
  call logit('Read in coarse matrices: head, wings, body, exchange')
  ! FHJ: this is to generate nice output / time estimate
  ! MJ: This is true irrespective of BSE/TDDFT
  if (flag%krnl==0) then ! triplet (only W)
    nmatrices = 3
  elseif (flag%krnl==2) then ! local fields (only V)
    nmatrices = 1
  elseif (flag%krnl==3) then ! spinor (W + V)
    nmatrices = 4
  else ! singlet (W + V)
    nmatrices = 4
  endif
  call progress_init(prog_info, 'interpolating BSE kernel', 'block', &
    maxval(peinf%ibt(:))*nmatrices)

  ! FHJ: For each PE, figure out which coarse k-points jkp it needs in order to
  ! interpolate the BSE matrix on the fine k-points it owns.
  SAFE_ALLOCATE(jkp2offset, (xct%nkpt_co))
  jkp2offset = -1 ! offset=-1 means we don`t need this coarse k-point
  jkp_offset = 0
  do jkp=1,xct%nkpt_co
    do iblock = 1, peinf%ibt(peinf%inode+1)
      ikp = peinf%ikb(iblock)
      if (any(fi2co_bse(:,ikp)==jkp)) then
        jkp2offset(jkp) = jkp_offset
        jkp_offset = jkp_offset + xct%nkpt_co
        exit
      endif
    enddo
  enddo

  SAFE_ALLOCATE(bsedmatrix_loc, (xct%n1b_co,xct%n2b_co,xct%n1b_co,xct%n2b_co,xct%nspin,xct%nspin,jkp_offset))
  nmatrices = 4
  do imatrix = 1, nmatrices
    ! triplet kernel has no exchange term
    if(xct%theory == 0 .and. (imatrix == 4 .and. flag%krnl == 0)) cycle
    if(xct%theory == 1 .and. (imatrix == 3 .and. flag%krnl == 0)) cycle
    ! local-field kernel has no direct term
    if(xct%theory == 0 .and. ((imatrix >= 1 .and. imatrix <= 3) .and. flag%krnl == 2)) cycle
    if(xct%theory == 1 .and. ((imatrix == 1 .or. imatrix == 2 .or. imatrix == 4) .and. flag%krnl == 2)) cycle

    SAFE_ALLOCATE(bsedmatrix, (xct%n1b_co,xct%n2b_co,xct%n1b_co,xct%n2b_co,xct%nspin,xct%nspin,xct%nkpt_co))
    do jkp=1,xct%nkpt_co

      call timing%start(timing%ik_input)

      bsedmatrix = 0.0
      if (peinf%inode.eq.0) then
        if(.not. xct%use_hdf5) then
          if (xct%theory==0) then
            if (imatrix .lt. 4) then
              ifile=11
            else if (imatrix .eq. 4) then
              ifile=12
            endif
          else if (xct%theory == 1) then
            if (imatrix .lt. 3) then
              ifile=11
            else if (imatrix .eq. 3) then
              ifile=12
            else if (imatrix .eq. 4) then
              ifile=13
            endif
          endif

! JRD: READ EVERYTHING FOR A GIVEN K

          do jj=1, xct%n2b_co * xct%n1b_co
            read(ifile) ikp,jcp,jvp,(((((bsedmatrix(jv,jc,jvp,jcp,js,jsp,jk), &
              jsp=1,xct%nspin),js=1,xct%nspin),jv=1,xct%n1b_co), &
              jc=1,xct%n2b_co),jk=1,xct%nkpt_co)
          enddo
#ifdef HDF5
        else
! XXX: Change to parallel
! JRD: HDF5 read in serial for now
! FHJ: FIXME - THIS IS BAD!!! We are doing serial read + BCast!!
          call bse_hdf5_read(xct,xct%nkpt_co,file_id,imatrix,jkp,bsedmatrix,subsample=.false.)
#endif
        endif
      endif

#ifdef MPI
      call MPI_BCAST(bsedmatrix(1,1,1,1,1,1,1), dimbse, MPI_SCALAR, 0, MPI_COMM_WORLD, mpierr)
#endif
      call timing%stop(timing%ik_input)

      if (jkp2offset(jkp)/=-1) then
        bsedmatrix_loc(:,:,:,:,:,:,jkp2offset(jkp)+1:jkp2offset(jkp)+xct%nkpt_co) = &
            bsedmatrix(:,:,:,:,:,:,1:xct%nkpt_co)
        if (xct%zero_q0_element==2) then
          ! FHJ: This can improve the interpolation in some cases: zeroing out
          ! the kernel matrix elements for q_co==0
          bsedmatrix_loc(:,:,:,:,:,:,jkp2offset(jkp)+jkp) = ZERO
        endif
      endif
    enddo
    SAFE_DEALLOCATE(bsedmatrix)

!==============================================================================
! FHJ: Beginning of the big loop over k-points
!      For each fine k-points ik/ikp and expansion vertices ivert/iterp, find
!      the corresponding coarse-grid points jk/jkp and interpolate the kernel.
!==============================================================================
    do iblock = 1, peinf%ibt(peinf%inode+1)
      call progress_step(prog_info)
      ikp = peinf%ikb(iblock)

      if (peinf%verb_debug) then
        write(tmpstr,'(a,i0)') 'new ikp = ', ikp
        call logit(tmpstr)
      endif
      do ivertp = 1, xct%npts_intp_kernel
        jkp = fi2co_bse(ivertp, ikp)
        jkp_offset = jkp2offset(jkp)

          do js=1,xct%nspin
            dcckp(:,:,js) = transpose(dcc(:,:,js,ikp,ivertp))
            if (xct%qflag.eq.1) then
              dvvkp(:,:,js) = transpose(dvv(:,:,js,ikp,ivertp))
            else
              if (.not. xct%skipinterp) then
                dvvkp(:,:,js) = transpose(dvv(:,:,js,xct%indexq_fi(ikp),ivertp))
              else
                dvvkp(:,:,js) = transpose(dvv(:,:,js,ikp,ivertp))
              endif
            endif
          enddo

            ! FHJ: please don`t fix the indentation yet
            do ik = 1, xct%nkpt_fi

              do ivert = 1, xct%npts_intp_kernel
                jk = fi2co_bse(ivert, ik)
                bsedmt = ZERO

                if (ivert==1) then
                ! FHJ: No need to update this block for different jk`s
                dq(:) = kg%f(:,ik) - kg%f(:,ikp)
                ! FHJ: and now we can reuse all the information from the cell
                ! structure to get the mapping to the BZ, v(q) and the interpolated
                ! epsilon in O(1) operations.
                call timing%start(timing%ik_cache)
                call cells_find_exactly(cells_fi, dq, ik_cells)
                if (ik_cells==0) then
                  write(0,'(a)') 'Found a point that was not calculated before:'
                  write(0,'(3(f13.9,1x))') dq
                  write(0,'(a)') 'Are you using a non-uniform grid?'
                  call die('Found a point that was not calculated before.')
                endif
                abs_q2 = abs_qs2(ik_cells)
                if (xct%theory==0) then
                  eps = eps_intp(ik_cells)
                else
                  eps = 1.0d0
                endif
                iold = vq_map(ik_cells)

!------------------- Calculate Coulomb Interaction-------------------------------------
!
! SIB: Here we calculate the values of 1/q^2 and 1/q (no truncation case)
! for the fine q vector q=k-kp.  Things to note:
!
! (1) If q=0, then we want instead the averages of 1/q^2 and 1/q over
!     the "mini-BZ" (mini-BZ is the 1st BZ but scaled down by the
!     number of total q-points).
!
!JRD:
! (2) If we are using truncation, the wing is DEFINED as being
!     some smooth function multiplied by |q|*Vtrunc(G=0,q)*eps^-1(G=0,G`=0,q)
!     This is because the wings of the dielectric matrix are /propto |q| and pick up a factor
!     of eps^-1(G=0,G`=0,q) during the inversion process (see Epsilon/epsinv.f90).
!     We include this factor here and not above because eps^-1(G=0,G`=0,q) varies quicker in truncated case.
!     The Vtrunc(G=0,q) factor comes from the bare (truncated) coulomb interaction.

! We do it as if it were a semiconductor because we want q=0 here since the divergent
! parts were already removed in kernel

                ! The 1/(8*pi) is already included
                vcoul = vcoul_array(iold)
                oneoverq = oneoverq_array(iold)
                call timing%stop(timing%ik_cache)

                if (xct%theory == 0) then
                  w_eff = vcoul * eps
                  if (abs_q2<TOL_ZERO .and. xct%iscreen==SCREEN_SEMICOND) then
                    w_eff = xct%wcoul0/(8.0*PI_D)
                  endif
                else if (xct%theory == 1) then
                  w_eff = vcoul
                  if (abs_q2<TOL_ZERO .and. xct%iscreen==SCREEN_SEMICOND) then
                    w_eff = vcoul
                  endif
                endif

!---------------- Start interpolation for the present pair (k,kp)----------------------------

!-----------------------------
! Head (spin diagonal)

                if (xct%ipar==1) then
                  icout=1
                  ivout=1
                else if (xct%ipar==2) then
                  icout=peinf%icb(iblock)
                  ivout=1
                else if (xct%ipar==3) then
                  icout=peinf%icb(iblock)
                  ivout=peinf%ivb(iblock)
                endif
                endif !ivert==1

                if (xct%zero_q0_element==1 .and. abs_q2<TOL_Zero) then
                  ! FHJ: If we zero out W(q=0), we get the continuum right, but
                  ! the first exciton is typically wrong.
                  ! See Figure 7 of PRB 93, 235435 (2016).
                  cycle
                endif
                call timing%start(timing%ik_interp)

                if (.not. xct%skipinterp) then
                  if (xct%subsample_line .and. sqrt(abs_q2) < xct%subsample_cutoff .and. sqrt(abs_q2) > TOL_Small) then
                    ! Can only expand around a vertex shared by ik and ikp:
                    ! jkp_sub = closepts_sub(ivertp_sub,ikp) = closepts_sub(ivert_sub,ik)
                    ! Find shared vertex that minimizes the distance between k_fi and k_co.
                    ! The closest point in closepts_sub is the first element.
                    jk_temp = 0
                    ivert_sub = 0
                    jkp_temp = 0
                    ivertp_temp = 0
                    jk_sub = 0
                    jkp_sub = 0
                    found_vertex=.false.
                    do ivertp_temp=1,xct%idimensions+1
                      jkp_temp = closepts_sub(ivertp_temp,ikp)
                      do ivert_temp=1,xct%idimensions+1
                        jk_temp = closepts_sub(ivert_temp,ik)
                        if (jk_temp.eq.jkp_temp) then
                          ivert_sub = ivert_temp
                          ivertp_sub = ivertp_temp
                          jkp_sub = jkp_temp
                          if (old_csi_behavior) found_vertex=.true.
                          exit
                        endif
                        if (found_vertex) exit
                      enddo
                    enddo
#ifdef DEBUG
                    if (peinf%inode.eq.0) write(6,*) "found jkp_sub = ",jkp_sub
#endif
                    if (jkp_sub.eq.0) call die("Could not find shared vertex. Subsampling cutoff may be too large",&
                      only_root_writes=.true.)
                    ! Find jk_sub s.t. |k_co-k_sub| is as close as possible to |q|.
                    dist_min=INF
                    do isub=1,nk_sub
#ifdef DEBUG
                      if (peinf%inode.eq.0) write(6,*) isub, kpoints_sub_len(isub), sqrt(abs_q2)
#endif
                      if( abs(kpoints_sub_len(isub) - sqrt(abs_q2)) < dist_min) then
                        dist_min = abs(kpoints_sub_len(isub) - sqrt(abs_q2))
                        jk_sub = isub
                      endif
                    enddo
#ifdef DEBUG
                    if (peinf%inode.eq.0) then
                      write(6,*) "found jk_sub = ",jk_sub
                    endif
#endif
                    ! Find dcckp_sub and dvvkp_sub analogous to dcckp and dvvkp
                    do js=1,xct%nspin
                      dcckp_sub(:,:,js) = transpose(dcc_sub(:,:,js,ikp,ivertp_sub,1))
                      if (xct%qflag.eq.1) then
                        dvvkp_sub(:,:,js) = transpose(dvv_sub(:,:,js,ikp,ivertp_sub,1))
                      else
                        if (.not. xct%skipinterp) then
                          dvvkp_sub(:,:,js) = transpose(dvv_sub(:,:,js,xct%indexq_fi(ikp),ivertp_sub,1))
                        else
                          dvvkp_sub(:,:,js) = transpose(dvv_sub(:,:,js,ikp,ivertp_sub,1))
                        endif
                      endif
                    enddo
                    ! Interpolate
                    do js=1,xct%nspin
                      do jsp=1,xct%nspin
                        bsedmt_sub_temp(:,:,:,:,js,jsp) = bsedmt_sub(:,:,:,:,jkp_sub,jk_sub,imatrix)
                      enddo
                    enddo
                    if (xct%qflag.eq.1) then
                      call interpolate(xct, &
                        bsedmt_sub_temp(:,:,:,:,:,:), bsedmt(:,:,:,:,:,:,1), &
                        dcc_sub(:,:,:,ik,ivert_sub,jk_sub), dcckp_sub, dvv_sub(:,:,:,ik,ivert_sub,jk_sub), dvvkp_sub, &
                        ivout, icout, imatrix, flag%krnl, .true.)
                    else
                      call interpolate(xct, &
                        bsedmt_sub_temp(:,:,:,:,:,:), bsedmt(:,:,:,:,:,:,1), &
                        dcc_sub(:,:,:,ik,ivert_sub,jk_sub), dcckp_sub, dvv_sub(:,:,:,xct%indexq_fi(ik),ivert_sub,jk_sub), &
                        dvvkp_sub, ivout, icout, imatrix, flag%krnl, .true.)
                    endif
                  else
                    if (xct%qflag.eq.1) then
                      call interpolate(xct, &
                        bsedmatrix_loc(:,:,:,:,:,:,jkp_offset+jk), bsedmt(:,:,:,:,:,:,1), &
                        dcc(:,:,:,ik,ivert), dcckp, dvv(:,:,:,ik,ivert), dvvkp, &
                        ivout, icout, imatrix, flag%krnl, .true.)
                    else
                      call interpolate(xct, &
                      bsedmatrix_loc(:,:,:,:,:,:,jkp_offset+jk), bsedmt(:,:,:,:,:,:,1), &
                      dcc(:,:,:,ik,ivert), dcckp, dvv(:,:,:,xct%indexq_fi(ik),ivert), dvvkp, &
                      ivout, icout, imatrix, flag%krnl, .true.)
                    endif
                  endif ! if (xct%subsample_line)
                  if (.not.xct%tda) then
                    if (xct%subsample_line .and. sqrt(abs_q2) < xct%subsample_cutoff .and. sqrt(abs_q2) > TOL_Small) then
                      if (xct%qflag.eq.1) then
                        call interpolate(xct, &
                          bsedmt_sub_temp(:,:,:,:,:,:), bsedmt(:,:,:,:,:,:,2), &
                          dcc_sub(:,:,:,ik,ivert_sub,jk_sub), dcckp_sub, dvv_sub(:,:,:,ik,ivert_sub,jk_sub), dvvkp_sub, &
                          ivout, icout, imatrix, flag%krnl, .false.)
                      else
                        call interpolate(xct, &
                          bsedmt_sub_temp(:,:,:,:,:,:), bsedmt(:,:,:,:,:,:,2), &
                          dcc_sub(:,:,:,ik,ivert_sub,jk_sub), dcckp_sub, dvv_sub(:,:,:,xct%indexq_fi(ik),ivert_sub,jk_sub), &
                          dvvkp_sub, ivout, icout, imatrix, flag%krnl, .false.)
                      endif
                    else
                      if (xct%qflag.eq.1) then
                        call interpolate(xct, &
                          bsedmatrix_loc(:,:,:,:,:,:,jkp_offset+jk), bsedmt(:,:,:,:,:,:,2), &
                          dcc(:,:,:,ik,ivert), dcckp, dvv(:,:,:,ik,ivert), dvvkp, &
                          ivout, icout, imatrix, flag%krnl, .false.)
                      else
                        call interpolate(xct, &
                          bsedmatrix_loc(:,:,:,:,:,:,jkp_offset+jk), bsedmt(:,:,:,:,:,:,2), &
                          dcc(:,:,:,ik,ivert), dcckp, dvv(:,:,:,xct%indexq_fi(ik),ivert), dvvkp, &
                          ivout, icout, imatrix, flag%krnl, .false.)
                      endif
                    endif
                  endif
                else
! XXX Thread?
                  if (.not.xct%extended_kernel) then
                    do jsp=1,xct%nspin
                      do js=1,xct%nspin
                        do jcp=1,peinf%nc_block
                          do jvp=1,peinf%nv_block
                            do jc=1,xct%ncb_fi
                              do jv=1,xct%nvb_fi
                                bsedmt(jv,jc,jvp,jcp,js,jsp,1) = bsedmatrix_loc( &
                                  jv, jc, jvp+ivout-1, jcp+icout-1, &
                                  js, jsp, jkp_offset+jk)
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                  else
! FHJ: A note on the indices:
! For restricted kernel, we distinguish valence and conduction states in the bsemat file:
! valence -> (vbm, vbm-1, ..., vbm-nvb_co+1)
! conduction -> (cbm, cbm+1, ..., cbm+ncb_co-1)
! For extended kernel, we stack both band indices together in the bsemat file in this funny ordering:
! band -> (vbm, ...,  vbm-nvb_co+1, cbm, ..., cbm+ncb_co-1)
! This only applies to bands from the coarse grid. Indices from fine grids
! (such as those in the bsedmt/hbse_a arrays) are either valence or conduction,
! but never mixed.
                    do jsp=1,xct%nspin
                      do js=1,xct%nspin
                        do jcp=1,peinf%nc_block
                          do jvp=1,peinf%nv_block
                            do jc=1,xct%ncb_fi
                              do jv=1,xct%nvb_fi
                                ! FHJ: (vc) -> (v`c`)
                                bsedmt(jv,jc,jvp,jcp,js,jsp,1) = bsedmatrix_loc( &
                                  jv, xct%nvb_co+jc, jvp+ivout-1, xct%nvb_co+jcp+icout-1, &
                                  js, jsp, jkp_offset+jk)
                                if (.not.xct%tda) then
                                  ! FHJ: (vc) -> (c`v`)
                                  bsedmt(jv,jc,jvp,jcp,js,jsp,2) = bsedmatrix_loc( &
                                    jv, xct%nvb_co+jc, xct%nvb_co+jcp+icout-1, jvp+ivout-1, &
                                    js, jsp, jkp_offset+jk)
                                endif
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                  endif
                endif

                call timing%stop(timing%ik_interp)

!------------------------------
! Add interaction kernel to the Hamiltonian

! JRD: bsemat_fac is the prefactor for the various head/wing/body/exchange and screening

! XXX: Some of this should be moved outside of the jk,jkp,ik,ikp loops right?

                if (imatrix .eq. 1) then
                  if (xct%iscreen==SCREEN_SEMICOND) then
                    if (xct%theory==0) then
                      bsemat_fac = fac_d * w_eff
                    else
                      bsemat_fac = fac_d * w_eff !* xct%coulomb_mod%short_range_frac_fock
                    endif
                  elseif (xct%iscreen==SCREEN_GRAPHENE) then
                    if (xct%icutv==TRUNC_NONE .and. xct%theory .eq. 0) then
                      bsemat_fac = fac_d * oneoverq
                    else
                      bsemat_fac = fac_d * w_eff
                    endif
                  else
                    bsemat_fac = fac_d
                  endif
                else if (imatrix .eq. 2) then
                  if (xct%iscreen==SCREEN_SEMICOND .and. xct%icutv==TRUNC_NONE .and. xct%theory==0) then
                    bsemat_fac = fac_d * oneoverq
                  elseif (xct%theory .eq. 0) then
                    bsemat_fac = fac_d
                  else
                    bsemat_fac = fac_d !* xct%coulomb_mod%short_range_frac_fock
                  endif
                else if (imatrix .eq. 3) then
                  if (xct%theory .eq. 0) then
                    bsemat_fac = fac_d
                  else
                    ! This matrix is exchange
                    bsemat_fac = -fac_x
                    if ((xct%nspin .eq. 1) .and. (flag%krnl .ne. 3)) bsemat_fac = bsemat_fac * 2D0
                  endif
                else if (imatrix .eq. 4) then
                  if (xct%theory .eq. 0) then
                    bsemat_fac = -fac_x
                    if ((xct%nspin .eq. 1) .and. (flag%krnl .ne. 3)) bsemat_fac = bsemat_fac * 2D0
                  else
                    ! change back
                    ! This is the fxc matrix
                    bsemat_fac = -fac_t !* ( 1.0d0 - xct%coulomb_mod%short_range_frac_fock)
                    if ((xct%nspin .eq. 1) .and. (flag%krnl .ne. 3)) bsemat_fac = bsemat_fac * 2D0
                  endif
                endif
                bsemat_fac = bsemat_fac * intp_coefs(ivert, ik) * intp_coefs(ivertp, ikp)

!#BEGIN_INTERNAL_ONLY
                if (imatrix.eq.4) then
                  bsemat_fac = bsemat_fac*xct%exchange_fact
                endif
                if (imatrix.ne.4) then
                  bsemat_fac = bsemat_fac*xct%direct_fact
                endif
!#END_INTERNAL_ONLY

                if(abs(bsemat_fac) < TOL_Zero) cycle
                call timing%start(timing%ik_sum)


! JRD: When comparing these loops with the loops in diag.f90 notice that:
! iv -> ivp
! ic -> icp
! ...
!
! hbse_a array is accessed the same way since ikcvs here is based on ik,ic,iv and
! ikcvsd is based on ikp,icp,ivp. This is the opposite of diag.f90 which makes the
! loops identical. Someone should fix the variable names to match at some point...

! XXX: Thread?
! XXX: WHY IS BSE_INDEX FASTER ON THE SPIN INDEX??


                do icb = 1, peinf%nc_block
                    do ivb = 1, peinf%nv_block
                        do jsp = 1, xct%nspin
                          ikcvsd = bse_index(iblock, icb, ivb, jsp, xct, &
                            ncband=peinf%nc_block, nvband=peinf%nv_block)

                          do ic = 1, xct%ncb_fi
                            do iv = 1, xct%nvb_fi
                              do js = 1,xct%nspin
                                ikcvs = bse_index(ik, ic, iv, js, xct)
                                hbse_a(ikcvs,ikcvsd) = hbse_a(ikcvs,ikcvsd) + &
                                  bsemat_fac * bsedmt(iv,ic,ivb,icb,js,jsp,1)
                                ! DYQ: Print exciton continuum
                                !if (imatrix.eq.1) then
                                !  if (ik.eq.peinf%ikb(peinf%inode+1,ikt) .and. (iv.eq.ivp) .and. (ic.eq.icp)) then
                                !    write(6,*) "Head Diag = ", bsemat_fac * bsedmt(iv,ic,ivb,icb,js,jsp,1)
                                !  endif
                                !endif
                                if (.not.xct%tda) then
                                  hbse_b(ikcvs,ikcvsd) = hbse_b(ikcvs,ikcvsd) + &
                                    bsemat_fac * bsedmt(iv,ic,ivb,icb,js,jsp,2)
                                endif
                              enddo              ! js
                            enddo              ! iv
                          enddo              ! ic
                        enddo              ! jsp
                    enddo              ! ivb
                enddo              ! icb
                call timing%stop(timing%ik_sum)
              enddo ! ivert / jk
            enddo ! ik
      enddo ! ivertp / jkp
    enddo ! ikp
  enddo ! imatrix

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif



  call progress_free(prog_info)
!==============================================================================
! FHJ: End of the big loop over k-points
!==============================================================================

  if (peinf%inode.eq.0) then
    if(.not. xct%use_hdf5) then
      call close_file(11)
      call close_file(12)
#ifdef HDF5
    else
      ! FHJ: TODO - change to collective I/O
      call h5fclose_f(file_id, error)
#endif
    endif
  endif

  call cells_free(cells_fi)

  SAFE_DEALLOCATE(bsedmatrix_loc)
  SAFE_DEALLOCATE(bsedmt)
  SAFE_DEALLOCATE(dcckp)
  SAFE_DEALLOCATE(dvvkp)
  SAFE_DEALLOCATE(dist_array)
  SAFE_DEALLOCATE(vcoul_array)
  SAFE_DEALLOCATE(oneoverq_array)
  SAFE_DEALLOCATE(dqs)
  SAFE_DEALLOCATE(dqs_bz)
  SAFE_DEALLOCATE(abs_qs2)
  if (xct%theory==0) then
    SAFE_DEALLOCATE(eps_intp)
  endif
  SAFE_DEALLOCATE(vq_map)
  SAFE_DEALLOCATE(fi2co_bse)
  SAFE_DEALLOCATE(jkp2offset)
  if (xct%subsample_line) then
    SAFE_DEALLOCATE(dcckp_sub)
    SAFE_DEALLOCATE(dvvkp_sub)
    SAFE_DEALLOCATE(bsedmt_sub)
    SAFE_DEALLOCATE(bsedmt_sub_temp)
    SAFE_DEALLOCATE(kpoints_sub)
    SAFE_DEALLOCATE(kpoints_sub_len)
  endif

  POP_SUB(intkernel)

  return
end subroutine intkernel

!=======================================================================================

subroutine interpolate(xct,bse_co,bse_fi,dcck,dcckp,dvvk,dvvkp,ivin,icin,imatrix,flag,resonant)
  type (xctinfo), intent(inout) :: xct
     ! xct should not be changed but ifort 12.0.0 -O3 compiles
     ! incorrectly if it is intent(in)
  SCALAR, intent(in) :: bse_co(:,:,:,:,:,:) !< (xct%n1b_co,xct%n2b_co,xct%n1b_co,xct%n2b_co,xct%nspin,xct%nspin)
  SCALAR, intent(out) :: bse_fi(:,:,:,:,:,:) !< (xct%nvb_fi,xct%ncb_fi,peinf%nv_block,peinf%nc_block,xct%nspin,xct%nspin)
  SCALAR, intent(in) :: dcck(:,:,:)  !< (xct%ncb_fi,xct%n2b_co,xct%nspin)
  SCALAR, intent(in) :: dcckp(:,:,:) !< (xct%n2b_co,xct%ncb_fi,xct%nspin)
  SCALAR, intent(in) :: dvvk(:,:,:)  !< (xct%nvb_fi,xct%n1b_co,xct%nspin)
  SCALAR, intent(in) :: dvvkp(:,:,:) !< (xct%n1b_co,xct%nvb_fi,xct%nspin)
  integer, intent(in) :: ivin, icin, imatrix, flag
  !> Interpolate to resonant or coupling block? We calculate the coupling block
  !! <vck|K|c`v`k> by applying complex conjugation to dvvkp and dcckp and switching
  !! the jvp and jcp indices in bse_co, which is easy b/c we have an extended kernel.
  logical, intent(in) :: resonant

  integer :: js,jsp,iv,ivp,jc,jcp,icb,ivb,jv,jvp,ic,icp, &
    js_dvvk, js_dvvkp, js_dcck, js_dcckp, bse_co_js, bse_co_jsp

  SCALAR, allocatable :: mat_vcvc(:,:,:,:),mat_vfvc(:,:,:,:), &
    mat_vfvf(:,:,:,:),mat_cccc(:,:,:,:), &
    mat_cfcc(:,:,:,:),mat_cfcf(:,:,:,:), &
    dummy(:,:,:,:),dummyp(:,:,:,:), &
    dummy2(:,:,:,:),dummy3(:,:,:,:), dvvkn(:,:)

  PUSH_SUB(interpolate)

  bse_fi=0.0

  do js=1,xct%nspin
    do jsp=1,xct%nspin

      if (xct%theory .eq. 0) then
        if (xct%nspin .eq. 1) then
          js_dcck=1
          js_dcckp=1
          js_dvvk=1
          js_dvvkp=1
          bse_co_js=1
          bse_co_jsp=1
        elseif (flag .eq. 0 .and. imatrix .ne. 4 .and. js .eq. jsp ) then
          js_dcck=js
          js_dcckp=js
          bse_co_js=js
          if (js .eq. 1) then
            js_dvvk=js+1
            js_dvvkp=js+1
            bse_co_jsp=js+1
          else
            js_dvvk=js-1
            js_dvvkp=js-1
            bse_co_jsp=js-1
          end if
        else if((flag == 1 .and. imatrix /= 4 .and. js == jsp) .or. &
                (flag == 1 .and. imatrix == 4) .or. (flag == 2 .and. imatrix == 4)) then
          js_dcck=js
          js_dcckp=jsp
          js_dvvk=js
          js_dvvkp=jsp
          bse_co_js=js
          bse_co_jsp=jsp
        else
          cycle
        end if
      else ! theory = 2
        if (xct%nspin .eq. 1) then
          js_dcck=1
          js_dcckp=1
          js_dvvk=1
          js_dvvkp=1
          bse_co_js=1
          bse_co_jsp=1
        elseif (flag .eq. 0 .and. imatrix .ge. 3 .and. js .eq. jsp ) then
          js_dcck=js
          js_dcckp=js
          bse_co_js=js
          if (js .eq. 1) then
            js_dvvk=js+1
            js_dvvkp=js+1
            bse_co_jsp=js+1
          else
            js_dvvk=js-1
            js_dvvkp=js-1
            bse_co_jsp=js-1
          end if
        else if((flag == 1 .and. imatrix .le. 2 .and. js == jsp) .or. &
                (flag == 1 .and. imatrix .ge. 3) .or. (flag == 2 .and. imatrix .ge. 3)) then
          js_dcck=js
          js_dcckp=jsp
          js_dvvk=js
          js_dvvkp=jsp
          bse_co_js=js
          bse_co_jsp=jsp
        else
          cycle
        end if
      endif

      if (xct%ipar >= 2) then

! Faster and better on memory when peinf%nv_block is not 1

        SAFE_ALLOCATE(dummy, (xct%n1b_co,xct%ncb_fi,xct%n1b_co,peinf%nc_block))
        dummy = 0.0

        do icp=1,peinf%nc_block
          icb = icin
          do jcp=1,xct%n2b_co
            do jvp=1,xct%n1b_co
              if (resonant) then
                call X(gemm)('n','c',xct%n1b_co,xct%ncb_fi,xct%n2b_co, &
                  dcckp(jcp,icb,js_dcckp),bse_co(:,:,jvp,jcp,bse_co_js,bse_co_jsp), &
                  xct%n1b_co,dcck(:,:,js_dcck),xct%ncb_fi,ONE,dummy(:,:,jvp,icp),xct%n1b_co)
              else
                call X(gemm)('n','c',xct%n1b_co,xct%ncb_fi,xct%n2b_co, &
                  MYCONJG(dcckp(jcp,icb,js_dcckp)),bse_co(:,:,jcp,jvp,bse_co_js,bse_co_jsp), &
                  xct%n1b_co,dcck(:,:,js_dcck),xct%ncb_fi,ONE,dummy(:,:,jvp,icp),xct%n1b_co)
              endif
            enddo
          enddo
        enddo

        SAFE_ALLOCATE(dummyp, (xct%n1b_co,xct%n1b_co,xct%ncb_fi,peinf%nc_block))
        SAFE_ALLOCATE(dummy2, (xct%n1b_co,xct%ncb_fi,peinf%nv_block,peinf%nc_block))
        SAFE_ALLOCATE(dummy3, (xct%n1b_co,peinf%nv_block,xct%ncb_fi,peinf%nc_block))
        SAFE_ALLOCATE(dvvkn, (xct%n1b_co,peinf%nv_block))
        dummy3 = 0.0
        dummy2 = 0.0

        do icp = 1, peinf%nc_block
          do ic = 1, xct%ncb_fi
            do jvp = 1, xct%n1b_co
              dummyp(:,jvp,ic,icp) = dummy(:,ic,jvp,icp)
            enddo
          enddo
        enddo

        if (resonant) then
          if (xct%ipar .eq. 2) then
            do ivp = 1,peinf%nv_block
              dvvkn(:,ivp) = MYCONJG(dvvkp(:,ivp,js_dvvkp))
            enddo
          else
            dvvkn(:,1) = MYCONJG(dvvkp(:,ivin,js_dvvkp))
          endif
        else
          if (xct%ipar .eq. 2) then
            do ivp = 1,peinf%nv_block
              dvvkn(:,ivp) = dvvkp(:,ivp,js_dvvkp)
            enddo
          else
            dvvkn(:,1) = dvvkp(:,ivin,js_dvvkp)
          endif
        endif

        do icp=1,peinf%nc_block
          icb = icin
          do ic=1,xct%ncb_fi
            call X(gemm)('n','n',xct%n1b_co,peinf%nv_block,xct%n1b_co, &
              ONE,dummyp(1,1,ic,icp),xct%n1b_co,dvvkn(1,1),xct%n1b_co,ONE,dummy3(1,1,ic,icp),xct%n1b_co)
          enddo
        enddo

        do icp = 1, peinf%nc_block
          do ic = 1, xct%ncb_fi
            dummy2(:,ic,:,icp) = dummy3(:,:,ic,icp)
          enddo
        enddo

        SAFE_DEALLOCATE(dummy)
        SAFE_DEALLOCATE(dummyp)
        SAFE_DEALLOCATE(dummy3)
        SAFE_DEALLOCATE(dvvkn)

        do icp=1,peinf%nc_block
          icb = icin
          do ivp=1,peinf%nv_block
            call X(gemm)('n','n',xct%nvb_fi,xct%ncb_fi,xct%n1b_co,ONE,dvvk(:,:,js_dvvk),xct%nvb_fi, &
              dummy2(:,:,ivp,icp),xct%n1b_co,ONE,bse_fi(:,:,ivp,icp,js,jsp),xct%nvb_fi)
          enddo
        enddo

        SAFE_DEALLOCATE(dummy2)

      else if (xct%ipar .eq. 1) then

! Fastest but worst on memory

        SAFE_ALLOCATE(mat_vcvc, (xct%n1b_co,xct%n1b_co,xct%n2b_co,xct%n2b_co))

        if (resonant) then
          do jcp=1,xct%n2b_co
            do jc=1,xct%n2b_co
              mat_vcvc(:xct%n1b_co,1:xct%n1b_co,jc,jcp) = &
                bse_co(1:xct%n1b_co,jc,1:xct%n1b_co,jcp,bse_co_js,bse_co_jsp)
            enddo
          enddo
        else
          do jcp=1,xct%n2b_co
            do jc=1,xct%n2b_co
              mat_vcvc(:xct%n1b_co,1:xct%n1b_co,jc,jcp) = &
                bse_co(1:xct%n1b_co,jc,jcp,1:xct%n1b_co,bse_co_js,bse_co_jsp)
            enddo
          enddo
        endif

! Interpolate v,v`

        SAFE_ALLOCATE(mat_vfvc, (xct%nvb_fi,xct%n1b_co,xct%n2b_co,xct%n2b_co))

        do jc=1,xct%n2b_co
          do jcp=1,xct%n2b_co
            mat_vfvc(1:xct%nvb_fi,1:xct%n1b_co,jc,jcp) = &
              MATMUL((dvvk(1:xct%nvb_fi,1:xct%n1b_co,js_dvvk)), &
              (mat_vcvc(1:xct%n1b_co,1:xct%n1b_co,jc,jcp)))
          enddo
        enddo

        SAFE_DEALLOCATE(mat_vcvc)
        SAFE_ALLOCATE(mat_vfvf, (xct%nvb_fi,xct%nvb_fi,xct%n2b_co,xct%n2b_co))

        if (resonant) then
          do jc=1,xct%n2b_co
            do jcp=1,xct%n2b_co
              mat_vfvf(1:xct%nvb_fi,1:xct%nvb_fi,jc,jcp) = &
                MATMUL((mat_vfvc(1:xct%nvb_fi,1:xct%n1b_co,jc,jcp)), &
                MYCONJG(dvvkp(1:xct%n1b_co,1:xct%nvb_fi,js_dvvkp)))
            enddo
          enddo
        else
          do jc=1,xct%n2b_co
            do jcp=1,xct%n2b_co
              mat_vfvf(1:xct%nvb_fi,1:xct%nvb_fi,jc,jcp) = &
                MATMUL((mat_vfvc(1:xct%nvb_fi,1:xct%n1b_co,jc,jcp)), &
                dvvkp(1:xct%n1b_co,1:xct%nvb_fi,js_dvvkp))
            enddo
          enddo
        endif

! Reorder from v,v` to c,c`

        SAFE_DEALLOCATE(mat_vfvc)
        SAFE_ALLOCATE(mat_cccc, (xct%n2b_co,xct%n2b_co,xct%nvb_fi,xct%nvb_fi))

        do jcp=1,xct%n2b_co
          do jc=1,xct%n2b_co
            mat_cccc(jc,jcp,1:xct%nvb_fi,1:xct%nvb_fi) = &
              mat_vfvf(1:xct%nvb_fi,1:xct%nvb_fi,jc,jcp)
          enddo
        enddo

! Interpolate c,c`

        SAFE_DEALLOCATE(mat_vfvf)
        SAFE_ALLOCATE(mat_cfcc, (xct%ncb_fi,xct%n2b_co,xct%nvb_fi,xct%nvb_fi))

        do iv=1,xct%nvb_fi
          do ivp=1,xct%nvb_fi
            mat_cfcc(1:xct%ncb_fi,1:xct%n2b_co,iv,ivp) = &
              MATMUL(MYCONJG(dcck(1:xct%ncb_fi,1:xct%n2b_co,js_dcck)), &
              (mat_cccc(1:xct%n2b_co,1:xct%n2b_co,iv,ivp)))
          enddo
        enddo

        SAFE_DEALLOCATE(mat_cccc)
        SAFE_ALLOCATE(mat_cfcf, (xct%ncb_fi,xct%ncb_fi,xct%nvb_fi,xct%nvb_fi))

        if (resonant) then
          do iv=1,xct%nvb_fi
            do ivp=1,xct%nvb_fi
              mat_cfcf(1:xct%ncb_fi,1:xct%ncb_fi,iv,ivp) = &
                MATMUL((mat_cfcc(1:xct%ncb_fi,1:xct%n2b_co,iv,ivp)), &
                (dcckp(1:xct%n2b_co,1:xct%ncb_fi,js_dcckp)))
            enddo
          enddo
        else
          do iv=1,xct%nvb_fi
            do ivp=1,xct%nvb_fi
              mat_cfcf(1:xct%ncb_fi,1:xct%ncb_fi,iv,ivp) = &
                MATMUL((mat_cfcc(1:xct%ncb_fi,1:xct%n2b_co,iv,ivp)), &
                MYCONJG(dcckp(1:xct%n2b_co,1:xct%ncb_fi,js_dcckp)))
            enddo
          enddo
        endif

        SAFE_DEALLOCATE(mat_cfcc)

! Reorder matrix

        do ivp=1,xct%nvb_fi
          do iv=1,xct%nvb_fi
            bse_fi(iv,1:xct%ncb_fi,ivp,1:xct%ncb_fi,js,jsp) = &
              mat_cfcf(1:xct%ncb_fi,1:xct%ncb_fi,iv,ivp)
          enddo
        enddo

        SAFE_DEALLOCATE(mat_cfcf)

      endif
    enddo
  enddo

  POP_SUB(interpolate)

  return
end subroutine interpolate

#ifdef HDF5
subroutine bse_hdf5_read(xct,nk,file_id,imatrix,ikp,bsedmatrix,subsample)
  type(xctinfo), intent(in) :: xct
  integer(HID_T), intent(in) :: file_id
  integer, intent(in) :: imatrix
  integer, intent(in) :: ikp, nk
  logical, intent(in) :: subsample
  SCALAR, intent(out) :: bsedmatrix(:,:,:,:,:,:,:)

  integer :: ik, ic, icp, iv, ivp, is, isp,ik_sub,ikprime
!  integer(HID_T) :: plist_id      ! Property list identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in mem
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HSIZE_T) :: count(7), offset(7)
  integer :: rank, error
  real(DP), allocatable :: data(:,:,:,:,:,:,:)

  PUSH_SUB(bse_hdf5_read)

  !write(6,*) "About to read matrix ", imatrix

  if(subsample) then
    ! Extracting <k_sub|H|k_co=1>
    ik_sub = ikp ! This is actually the file number
    ikprime = 1  ! coarse point is always the first point in the file
  else
    ikprime = ikp
  endif

  if (xct%theory == 0) then
    if (imatrix .eq. 1) then
      call h5dopen_f(file_id, 'mats/head', dset_id, error)
    else if (imatrix .eq. 2) then
      call h5dopen_f(file_id, 'mats/wing', dset_id, error)
    else if (imatrix .eq. 3) then
      call h5dopen_f(file_id, 'mats/body', dset_id, error)
    else if (imatrix .eq. 4) then
      call h5dopen_f(file_id, 'mats/exchange', dset_id, error)
    endif
  else
    if (imatrix .eq. 1) then
      call h5dopen_f(file_id, 'mats/head', dset_id, error)
    else if (imatrix .eq. 2) then
      call h5dopen_f(file_id, 'mats/body', dset_id, error)
    else if (imatrix .eq. 3) then
      call h5dopen_f(file_id, 'mats/exchange', dset_id, error)
    else if (imatrix .eq. 4) then
      call h5dopen_f(file_id, 'mats/fxc', dset_id, error)
    endif
  endif

  rank = 7
  count(1) = SCALARSIZE
  count(2) = xct%n1b_co
  count(3) = xct%n1b_co
  count(4) = xct%n2b_co
  count(5) = xct%n2b_co
  count(6) = nk
  count(7) = 1

  SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6),count(7)))

  do isp = 1, xct%nspin
    do is = 1, xct%nspin

! Get filespace

      call h5screate_simple_f(rank, count, memspace, error)
      call h5dget_space_f(dset_id,filespace,error)

! Construct data and offset

      offset(1)=0
      offset(2)=0
      offset(3)=0
      offset(4)=0
      offset(5)=0
      offset(6)=(is-1)*xct%nkpt_co
      offset(7)=(ikprime-1)+(isp-1)*xct%nkpt_co     ! RJW

! Select hyperslab

      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

! Create property list for collective dataset read
! Read is serial for now

      !call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

! Collectively read the file
! Serial for now

      !call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, filespace, &
      !                  xfer_prp = plist_id)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, filespace)

! copy from data into bsedmatrix which holds everything for one kpoint
! XXX This isn`t optimal because bsedmatrix layout is funny.

      do ik=1,count(6)
        do icp=1,count(5)
          do ic=1,count(4)
            do ivp=1,count(3)
              do iv=1,count(2)
                if (subsample) then
                  bsedmatrix(iv,ic,ivp,icp,ik_sub,ik,imatrix)=SCALARIFY2(data(1,iv,ivp,ic,icp,ik,1),data(2,iv,ivp,ic,icp,ik,1))
                else
                  bsedmatrix(iv,ic,ivp,icp,is,isp,ik)=SCALARIFY2(data(1,iv,ivp,ic,icp,ik,1),data(2,iv,ivp,ic,icp,ik,1))
                endif
              enddo
            enddo
          enddo
        enddo
      enddo

! Close property list and filespace. Do we have to do this each time?

      call h5sclose_f(memspace, error)
      call h5sclose_f(filespace, error)

! JRD: Read is serial for now
      !call h5pclose_f(plist_id, error)

    enddo
  enddo

  call h5dclose_f(dset_id, error)

  SAFE_DEALLOCATE(data)
  POP_SUB(bse_hdf5_read)

end subroutine bse_hdf5_read
#endif

!> FHJ: Move a point (qq) to the 1st BZ (WS cell). Output vector is qq_bz,
!! with length abs_qq2.
subroutine point_to_bz(crys, qq, qq_bz, abs_qq2)
  type(crystal), intent(in) :: crys
  real(DP), intent(in) :: qq(3)
  real(DP), intent(out) :: qq_bz(3)
  real(DP), intent(out) :: abs_qq2

  real(DP) :: abs_tmp, qq_tmp(3)
  integer :: i1, i2, i3

  !no push/pop, called too often
  abs_qq2 = INF
  do i1=-ncell,ncell+1
    qq_tmp(1) = qq(1) - i1
    do i2=-ncell,ncell+1
      qq_tmp(2) = qq(2) - i2
      do i3=-ncell,ncell+1
        qq_tmp(3) = qq(3) - i3
        abs_tmp = DOT_PRODUCT(qq_tmp, MATMUL(crys%bdot, qq_tmp))
        if (abs_tmp < abs_qq2) then
          abs_qq2 = abs_tmp
          qq_bz(:) = qq_tmp(:)
        endif
      enddo
    enddo
  enddo

end subroutine point_to_bz

end module intkernel_m
