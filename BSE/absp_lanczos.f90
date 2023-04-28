!===============================================================================
!
! Routines:
!
! (1) absp_lanczos()     Originally By: FHJ     Last Modified   May/2016 (FHJ)
!
!     Calculate the absorption with electron-hole interaction
!     using the Lanczos algorithm developed by M. Shao and C. Yang.
!
!     See notes in absp.f90 and absp0.f90 about the units and prefactors.
!
!===============================================================================

#include "f_defs.h"

module absp_lanczos_m

#ifdef USESCALAPACK
  use global_m
  use diagonalize_m
  use misc_m
  use scalapack_m
  use timing_m, only: timing => bse_timing
  implicit none

  ! FHJ: FIXME - remove these module vars after Meiyue and I debug the solver
  real(DP) :: dt
  integer :: nmax_eff

  public :: absp_lanczos

  private

contains

subroutine absp_lanczos(xct, nmat, nmax, s1, vol, omega_plasma, &
                        flag, ipol, en_diff_max, hbse_a, hbse_b)
  type(xctinfo), intent(in) :: xct
  integer, intent(in) :: nmat
  integer, intent(in) :: nmax
  SCALAR, intent(in) :: s1(nmat)
  real(DP), intent(in) :: vol, omega_plasma
  type (flags), intent(in) :: flag
  integer, intent(in) :: ipol
  real(DP), intent(in) :: en_diff_max
  SCALAR, intent(inout) :: hbse_a(nmat,peinf%nblocks*peinf%block_sz)
  SCALAR, intent(inout), optional :: hbse_b(nmat,peinf%nblocks*peinf%block_sz)

  integer :: ii,iw,nwstep
  integer :: n_omega
  real(DP), allocatable :: omega(:), eps2(:)
  real(DP) :: emin,emax,iemax,domega
  real(DP) :: fac,fac1,fac2,sum0,sum1,sum2,pref
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
  ! FHJ: FIXME - flags to test the solver. Keep them for a while..
  integer :: nmax_iter, ialgo
  type (flags) :: flag2
  character(len=64) :: fname_base
  character(len=7) :: theories(3) = (/'gauss  ', 'nogauss', 'paired '/)
  
  PUSH_SUB(absp_lanczos)

  if (.not.xct%tda.and..not.present(hbse_b)) then
    call die('Internal parameter error in absp_lanczos.', only_root_writes=.true.)
  endif
  
  pref = 16.d0 * PI_D**2 / (vol * dble(xct%nspinor))
  ! the factor 1/sqrt(nspin) in absp0 is handled by the normalization of the exciton here

  emin = 0.d0
  emax = en_diff_max * ryd
  iemax = int(emax + 10.d0 * xct%eta) + 1
  emax = dble(iemax)
  domega = xct%delta_frequency
  nwstep = int(iemax / domega)

  ! FHJ: Note that omega is in ryd because the solver needs it like that.
  n_omega = nwstep + 1
  SAFE_ALLOCATE(omega, (n_omega))
  SAFE_ALLOCATE(eps2, (n_omega))
  omega(:) = (/ (emin + (emax - emin) * dble(iw) / (ryd*nwstep), iw=0,nwstep) /)

  if (peinf%inode==0) then

    write(6,'()')
    if (xct%npol==1) then
      if (flag%opr==0) then
        write(6,'(1x,a,3(f10.6,1x))') 'Polarization: ', xct%shift
      else
        write(6,'(1x,a,3(f10.6,1x))') 'Polarization: ', xct%pol
      endif
    else
      write(6,'(1x,a)') 'Polarization: '//suffix(ipol)
    endif
  endif

  if (flag%debug_lanczos) then
    call timing%start(timing%lanczos)
    flag2 = flag
    ! gauss, nogauss, paired
    do ialgo = 1, 3
      if (peinf%inode==0) write(6,'(1x,a,i0,a)') 'Algorthm ', ialgo, ' / 3'
      flag2%lanczos_gauss_quad = ialgo==1
      if (xct%npol==1) then
        fname_base = 'lanczos_'//trim(theories(ialgo))//'_'
      else
        fname_base = 'lanczos_'//suffix(ipol)//'_'//trim(theories(ialgo))//'_'
      endif

      do nmax_iter = 1, nmax
        if (peinf%inode==0) write(6,'(3x,a,i0,a,i0)') 'Iteration ', nmax_iter, ' / ', nmax
        if (ialgo<3) then
          call bsabsp(xct, nmat, nmax_iter, s1, n_omega, omega, eps2, flag2, hbse_a, hbse_b)
        else
          call bsabsp_paired(xct, nmat, nmax_iter, s1, n_omega, omega, eps2, flag2, hbse_a, hbse_b)
        endif
        write(fname,'(a,i0,a)') trim(fname_base), nmax_iter ,'.dat'
        if (peinf%inode==0) then
          call open_file(10,file=trim(fname),form='formatted',status='replace')
          write(10,'(a,i0)') "# Nmax_eff: ", nmax_eff
          write(10,'(a,es15.6e3)') "# Time: ", dt
          write(10,'(a)') "# Column 1: omega"
          write(10,'(a)') "# Column 2: eps2(omega)"
          do iw=1,n_omega
            write(10,'(2f16.9)') omega(iw)*ryd, eps2(iw)*pref
          enddo
          call close_file(10)
        endif
        if (nmax_iter>nmax_eff) exit
      enddo !nmax_iter
    enddo !ialgo
    call timing%stop(timing%lanczos)
    return
  else
    call timing%start(timing%lanczos)
    !call bsabsp(xct, nmat, nmax, s1, n_omega, omega, eps2, flag, hbse_a, hbse_b)
    call bsabsp_paired(xct, nmat, nmax, s1, n_omega, omega, eps2, flag, hbse_a, hbse_b)
    call timing%stop(timing%lanczos)
  endif

  if (peinf%inode==0) then

    if (xct%npol==1) then
      fname = 'absorption_eh.dat'
    else
      fname = 'absorption_'//suffix(ipol)//'_eh.dat'
    endif
    call open_file(10,file=trim(fname),form='formatted',status='replace')
    
    write(10,*) "# Column 1: omega"
    write(10,*) "# Column 2: eps2(omega)"

    do iw=1,n_omega
      write(10,'(2f16.9)') omega(iw)*ryd, eps2(iw)*pref
    enddo
    
    call close_file(10)

  endif

  POP_SUB(absp_lanczos)
  
end subroutine absp_lanczos

!> FHJ: this routines puts the Hamiltonian in the correct block-cyclic layout,
!! and calls PDBSABSP/PZBSABSP to calculate the optical absorption spectrum.
subroutine bsabsp(xct, nmat, nmax, s1, n_omega, omega, eps2, flag, hbse_a, hbse_b)
  type(xctinfo), intent(in) :: xct
  integer, intent(in) :: nmat
  integer, intent(in) :: nmax
  SCALAR, intent(in) :: s1(nmat)
  integer :: n_omega
  real(DP), intent(in) :: omega(n_omega)
  real(DP), intent(out) :: eps2(n_omega)
  type (flags), intent(in) :: flag
  SCALAR, intent(inout) :: hbse_a(nmat,peinf%nblocks*peinf%block_sz)
  SCALAR, intent(inout), optional :: hbse_b(nmat,peinf%nblocks*peinf%block_sz)

  ! Scalapack and blacs arrays
  type(scalapack) :: scal_h, scal_x, scal_d
  integer, dimension(9) :: desc_h_1d, desc_h_2d, desc_x_2d, desc_d_2d
  integer :: lwork, liwork
#ifdef CPLX
  integer :: lrwork
  real(DP), allocatable :: rwork(:)
#endif
  SCALAR, allocatable :: hbse_a_bl(:), hbse_b_bl(:), x_bl(:), s1_bl(:)
  complex(DPC), allocatable :: work(:)
  integer, allocatable :: iwork(:)
  real(DP) :: lambda(nmat), sigma
  integer :: cntxt_1d,  usermap(1,peinf%npes), lld_h, lld_x, lld_d, ii
  integer :: solver, irow_l, irow_g, info

  include 'full_solver/solver.f90'

  PUSH_SUB(bsabsp)

  if (.not.xct%tda.and..not.present(hbse_b)) then
    call die('Internal parameter error in absp_lanczos.', only_root_writes=.true.)
  endif

  eps2(:) = 0d0

  ! FHJ: Create BLACS context for 1d cyclic matrices. Processors: all aboard!
  call blacs_get(-1, 0, cntxt_1d)
  do ii=1,peinf%npes
    usermap(1,ii) = ii-1
  enddo
  call blacs_gridmap(cntxt_1d, usermap, 1, 1, peinf%npes) ! This will modify cntxt_1d
  call descinit(desc_h_1d, nmat, nmat, nmat, peinf%block_sz, 0, 0, cntxt_1d, nmat, info)
  if (info/=0) call die('got info/=0 in descinit')

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices hbse_{a,b}_bl.
  ! Some processors might get excluded for the sake of load balancing.
  call create_blacs_square_grid(scal_h)
  call create_blacs_rect_desc(nmat, nmat, scal_h, desc_h_2d, lld_h)

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices evecs_{r,l}_bl.
  ! The processor grid doesn`t change, so we reuse the context from scal_h%icntxt.  
  scal_x = scal_h
  scal_d = scal_h
  call create_blacs_rect_desc(2*nmat, nmax+1, scal_x, desc_x_2d, lld_x)
  call create_blacs_rect_desc(nmat, 1, scal_d, desc_d_2d, lld_d)

  SAFE_ALLOCATE(hbse_a_bl, (lld_h*max(1,scal_h%npc)))
  hbse_a_bl(:) = ZERO
  if (.not.xct%tda) then
    SAFE_ALLOCATE(hbse_b_bl, (lld_h*max(1,scal_h%npc)))
  else
    SAFE_ALLOCATE(hbse_b_bl, (1))
  endif
  hbse_b_bl(:) = ZERO
  SAFE_ALLOCATE(x_bl, (lld_x*max(1,scal_x%npc)))
  x_bl(:) = ZERO
  SAFE_ALLOCATE(s1_bl, (lld_d))
  s1_bl(:) = ZERO

  ! FHJ: Copy matrices hbse_{a,b} from 1d cyclic layout into hbse_{a,b}_bl (2d cyclic layout)
  call pX(gemr2d)(nmat, nmat, hbse_a, 1, 1, desc_h_1d, hbse_a_bl, 1, 1, desc_h_2d, cntxt_1d)
  if (.not.xct%tda) then
    call pX(gemr2d)(nmat, nmat, hbse_b, 1, 1, desc_h_1d, hbse_b_bl, 1, 1, desc_h_2d, cntxt_1d)
  endif
  ! FHJ: Copy global D into local D by looping over all of the local columns.
  ! Note that scal_d%npr==0 automatically if scal_d%icntxt<0.
  do irow_l = 1, scal_d%npr
    irow_g = INDXL2G(irow_l, scal_d%nbl, scal_d%myprow, 0, scal_d%nprow)
    s1_bl(irow_l) = s1(irow_g)
  enddo

  ! FHJ: zero coupling block if the user requests. Useful to test the solver.
  if (xct%zero_coupling_block) hbse_b_bl = ZERO
  lambda(:) = 0d0

  solver = BSE_LANCZOS
  if (xct%tda) then
    solver = solver + BSE_TDA
  else
    solver = solver + BSE_FULLBSE
  endif
  if (flag%lanczos_gauss_quad) solver = solver + BSE_QUADAVGGAUSS
  if (flag%lor==0) then
    solver = solver + BSE_LORENTZIAN
  else
    solver = solver + BSE_GAUSSIAN
  endif
  sigma = xct%eta/ryd
  if (peinf%inode==0 .and. peinf%verb_high) then
    write(6,'(/1x,a)') 'Omega grid:'
    write(6,'(1(5(1x,f15.9)))') omega
    write(6,'(1x,a,f0.9/)') 'Sigma: ', sigma
  endif

  if (scal_h%icntxt>=0) then
    if (peinf%inode==0) write(6,*)
    if (.not.flag%debug_lanczos) then
      if (.not.xct%tda) then
        call blacs_calc_broken_herm(hbse_a_bl, nmat, size(hbse_a_bl), desc_h_2d, 'A')
        call blacs_calc_broken_herm(hbse_b_bl, nmat, size(hbse_b_bl), desc_h_2d, 'B')
      else
        call blacs_calc_broken_herm(hbse_a_bl, nmat, size(hbse_a_bl), desc_h_2d, 'H')
      endif
    endif
    if (peinf%inode==0) write(6,*)

    lwork = -1
    SAFE_ALLOCATE(work, (10))
#ifdef CPLX
    lrwork = -1
    SAFE_ALLOCATE(rwork, (10))
#endif
    liwork = -1
    SAFE_ALLOCATE(iwork, (10))

    ! FHJ: Call p*bsabsp in query mode
    CALL pX(BSABSP)(solver, nmat, &
      n_omega, sigma, omega, eps2, &
      hbse_a_bl, 1, 1, desc_h_2d, &
      hbse_b_bl, 1, 1, desc_h_2d, lambda, &
      x_bl, 1, 1, desc_x_2d, &
      s1_bl, 1, 1, desc_d_2d, nmax, &
      work, lwork, &
#ifdef CPLX
      rwork, lrwork, &
#endif
      iwork, liwork, info)
    if (info/=0) then
      if(peinf%inode==0) then
        write(0,'(/,a,i0)') 'ERROR: Query mode for p*bsabsp failed with info=',info
        write(0,'(a,/)') 'Please, report this error to the BerkeleyGW developers!'
      endif
      call die("p*bseig internal error", only_root_writes=.true.)
    endif

    lwork = int(work(1))
    SAFE_DEALLOCATE(work)
    SAFE_ALLOCATE(work, (lwork))
#ifdef CPLX
    lrwork = int(rwork(1))
    SAFE_DEALLOCATE(rwork)
    SAFE_ALLOCATE(rwork, (lrwork))
#endif
    liwork = iwork(1)
    SAFE_DEALLOCATE(iwork)
    SAFE_ALLOCATE(iwork, (liwork))

    ! FHJ: Call p*bsabsp for realz
    if (peinf%inode==0) write(6,'(1x,a,i0,a,i0)') &
      'Beginning Lanczos iteration with ',nmax, ' steps. Matrix size: ', nmat
    dt = -MPI_WTIME()
    CALL pX(BSABSP)(solver, nmat, &
      n_omega, sigma, omega, eps2, &
      hbse_a_bl, 1, 1, desc_h_2d, &
      hbse_b_bl, 1, 1, desc_h_2d, lambda, &
      x_bl, 1, 1, desc_x_2d, &
      s1_bl, 1, 1, desc_d_2d, nmax, &
      work, lwork, &
#ifdef CPLX
      rwork, lrwork, &
#endif
      iwork, liwork, info)
    dt = dt + MPI_WTIME()
    if (peinf%inode==0) write(6,'(1x,a)') 'Done with Lanczos iteration.'

    SAFE_DEALLOCATE(work)
#ifdef CPLX
    SAFE_DEALLOCATE(rwork)
#endif
    SAFE_DEALLOCATE(iwork)

    if (info<0) then
      if (peinf%inode==0) then
        write(0,'(/,a,i0)') 'ERROR: Call to p*bsabsp failed with info=',info
        write(0,'(a,/)') 'Please, report this error to the BerkeleyGW developers!'
      endif
      call die("p*bseig internal error", only_root_writes=.true.)
    elseif (peinf%inode==0) then
      if (info==0) then
        write(6,'(1x,a,i0,a)') 'Lanczos algorithm performed all ', nmax, ' steps.'
      else
        write(6,'(1x,a,i0,a)') &
          'Lanczos algorithm converged to an invariant subspace after ', info, ' steps.'
      endif
    endif
    ! FHJ: TODO: save ritz vectors and values to reuse them it requested
  endif !scal_h%icntxt>=0
  nmax_eff = nmax
  if (info>0) nmax_eff = info

  ! FHJ: PEs outside the grid don`t have eps2!
  call MPI_Bcast(eps2, n_omega, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(nmax_eff, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

  ! FHJ: Cleanup
  if (scal_h%icntxt>=0) call blacs_gridexit(scal_h%icntxt)
  call blacs_gridexit(cntxt_1d)

  SAFE_DEALLOCATE(hbse_a_bl)
  SAFE_DEALLOCATE(hbse_b_bl)
  SAFE_DEALLOCATE(x_bl)
  SAFE_DEALLOCATE(s1_bl)

  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  if (peinf%inode==0) write(6,'(1x,a)') 'Lanczos step done.'

  POP_SUB(bsabsp)

end subroutine bsabsp

subroutine bsabsp_paired(xct, nmat, nmax, s1, n_omega, omega, eps2, flag, hbse_a, hbse_b)
  type(xctinfo), intent(in) :: xct
  integer, intent(in) :: nmat
  integer, intent(in) :: nmax
  SCALAR, intent(in) :: s1(nmat)
  integer :: n_omega
  real(DP), intent(in) :: omega(n_omega)
  real(DP), intent(out) :: eps2(n_omega)
  type (flags), intent(in) :: flag
  SCALAR, intent(inout) :: hbse_a(nmat,peinf%nblocks*peinf%block_sz)
  SCALAR, intent(inout), optional :: hbse_b(nmat,peinf%nblocks*peinf%block_sz)

  ! Scalapack and blacs arrays
  type(scalapack) :: scal_h, scal_x, scal_d
  integer, dimension(9) :: desc_h_1d, desc_h_2d, desc_x_2d, desc_d_2d
  integer :: liwork, lwork
  SCALAR, allocatable :: hbse_a_bl(:), hbse_b_bl(:), x_bl(:), s1_bl(:)
  complex(DPC), allocatable :: work(:)
  real(DP) :: m_mat(4*nmax*nmax), hx(4*nmax*nmax), ha(2*nmax)
  complex(DPC) :: hb(2*nmax)
  integer, allocatable :: iwork(:)
  real(DP) :: lambda(nmat), sigma, dtmp
  integer :: cntxt_1d,  usermap(1,peinf%npes), lld_h, lld_x, lld_d, ii
  integer :: irow_l, irow_g, info, nmax_real
  double precision, external :: pzlange

  PUSH_SUB(bsabsp_paired)

  if (.not.xct%tda.and..not.present(hbse_b)) then
    call die('Internal parameter error in absp_lanczos.', only_root_writes=.true.)
  endif

  eps2(:) = 0d0

  ! FHJ: Create BLACS context for 1d cyclic matrices. Processors: all aboard!
  call blacs_get(-1, 0, cntxt_1d)
  do ii=1,peinf%npes
    usermap(1,ii) = ii-1
  enddo
  call blacs_gridmap(cntxt_1d, usermap, 1, 1, peinf%npes) ! This will modify cntxt_1d
  call descinit(desc_h_1d, nmat, nmat, nmat, peinf%block_sz, 0, 0, cntxt_1d, nmat, info)
  if (info/=0) call die('got info/=0 in descinit')

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices hbse_{a,b}_bl.
  ! Some processors might get excluded for the sake of load balancing.
  call create_blacs_square_grid(scal_h)
  call create_blacs_rect_desc(nmat, nmat, scal_h, desc_h_2d, lld_h)

  ! FHJ: Initialize BLACS grid for 2d cyclic matrices evecs_{r,l}_bl.
  ! The processor grid doesn`t change, so we reuse the context from scal_h%icntxt.  
  scal_x = scal_h
  scal_d = scal_h
  call create_blacs_rect_desc(2*nmat, nmax+1, scal_x, desc_x_2d, lld_x)
  call create_blacs_rect_desc(nmat, 1, scal_d, desc_d_2d, lld_d)

  SAFE_ALLOCATE(hbse_a_bl, (lld_h*max(1,scal_h%npc)))
  hbse_a_bl(:) = ZERO
  if (.not.xct%tda) then
    SAFE_ALLOCATE(hbse_b_bl, (lld_h*max(1,scal_h%npc)))
  else
    SAFE_ALLOCATE(hbse_b_bl, (1))
  endif
  hbse_b_bl(:) = ZERO
  SAFE_ALLOCATE(x_bl, (lld_x*max(1,scal_x%npc)))
  x_bl(:) = ZERO
  SAFE_ALLOCATE(s1_bl, (lld_d))
  s1_bl(:) = ZERO

  hb = (0d0, 0d0)
  ha = 0d0
  hx = 0d0
  m_mat = 0d0
  ! FHJ: Copy matrices hbse_{a,b} from 1d cyclic layout into hbse_{a,b}_bl (2d cyclic layout)
  call pX(gemr2d)(nmat, nmat, hbse_a, 1, 1, desc_h_1d, hbse_a_bl, 1, 1, desc_h_2d, cntxt_1d)
  if (.not.xct%tda) then
    call pX(gemr2d)(nmat, nmat, hbse_b, 1, 1, desc_h_1d, hbse_b_bl, 1, 1, desc_h_2d, cntxt_1d)
  endif
  ! FHJ: Copy global D into local D by looping over all of the local columns.
  ! Note that scal_d%npr==0 automatically if scal_d%icntxt<0.
  do irow_l = 1, scal_d%npr
    irow_g = INDXL2G(irow_l, scal_d%nbl, scal_d%myprow, 0, scal_d%nprow)
    s1_bl(irow_l) = s1(irow_g)
  enddo

  ! FHJ: zero coupling block if the user requests. Useful to test the solver.
  if (xct%zero_coupling_block) hbse_b_bl = ZERO
  lambda(:) = 0d0

  sigma = xct%eta/ryd
  if (peinf%inode==0 .and. peinf%verb_high) then
    write(6,'(/1x,a)') 'Omega grid:'
    write(6,'(1(5(1x,f15.9)))') omega
    write(6,'(1x,a,f0.9/)') 'Sigma: ', sigma
  endif

  if (scal_h%icntxt>=0) then
    if (peinf%inode==0) write(6,*)
    if (.not.flag%debug_lanczos) then
      if (.not.xct%tda) then
        call blacs_calc_broken_herm(hbse_a_bl, nmat, size(hbse_a_bl), desc_h_2d, 'A')
        call blacs_calc_broken_herm(hbse_b_bl, nmat, size(hbse_b_bl), desc_h_2d, 'B')
      else
        call blacs_calc_broken_herm(hbse_a_bl, nmat, size(hbse_a_bl), desc_h_2d, 'H')
      endif
    endif
    if (peinf%inode==0) write(6,*)

    lwork = -1
    SAFE_ALLOCATE(work, (10))
    liwork = -1
    SAFE_ALLOCATE(iwork, (10))

    if (peinf%inode==0) write(6,'(1x,a)') 'Calling zbssolver in query mode'
    call zbssolver1(nmax, m_mat, 2*nmax, lambda, hx, 2*nmax, work, lwork, iwork, &
      liwork, info)

    if (info/=0) then
      if(peinf%inode==0) then
        write(0,'(/,a,i0)') 'ERROR: Query mode for zbssolver1 failed with info=',info
        write(0,'(a,/)') 'Please, report this error to the BerkeleyGW developers!'
      endif
      call die("zbssolver1 internal error", only_root_writes=.true.)
    endif

    lwork = max(int(work(1)), 2*nmat)
    SAFE_DEALLOCATE(work)
    SAFE_ALLOCATE(work, (lwork))
    liwork = iwork(1)
    SAFE_DEALLOCATE(iwork)
    SAFE_ALLOCATE(iwork, (liwork))

    dt = -MPI_WTIME()
    if (peinf%inode==0) write(6,'(1x,a)') 'Calling pzcopy'
    call pzcopy(nmat, &
      s1_bl, 1, 1, desc_d_2d, 1, &
      x_bl, 1, 1, desc_x_2d, 1)

    if (peinf%inode==0) write(6,'(1x,a)') 'Calling pzsplanczos'
    call pzsplanczos(nmat, nmax, &
      hbse_a_bl, 1, 1, desc_h_2d, &
      hbse_b_bl, 1, 1, desc_h_2d, &
      x_bl, 1, 1, desc_x_2d, &
      ha, hb, work, desc_x_2d, info)

    if (info<0) then
      if (peinf%inode==0) then
        write(0,'(/,a,i0)') 'ERROR: Call to pzsplanczos failed with info=',info
        write(0,'(a,/)') 'Please, report this error to the BerkeleyGW developers!'
      endif
      call die("pzsplanczos internal error", only_root_writes=.true.)
    elseif (peinf%inode==0) then
      if (info==0) then
        write(6,'(1x,a,i0,a)') 'Lanczos algorithm performed all ', nmax, ' steps.'
      else
        write(6,'(1x,a,i0,a)') &
          'Lanczos algorithm converged to an invariant subspace after ', info, ' steps.'
      endif
    endif

    nmax_real = nmax
    if (info>0) nmax_real = info
    nmax_eff = nmax_real

    if (peinf%inode==0) write(6,'(1x,a)') 'Calling zembed3'
    call zembed3(nmax_real, ha, hb, m_mat, 2*nmax)

    if (peinf%inode==0) write(6,'(1x,a)') 'Calling zbssolver1 (real call)'
    call zbssolver1(nmax_real, m_mat, 2*nmax, lambda, hx, 2*nmax, work, lwork, iwork, &
      liwork, info)

    if (info<0) then
      if (peinf%inode==0) then
        write(0,'(/,a,i0)') 'ERROR: Call to zbssolver1 failed with info=',info
        write(0,'(a,/)') 'Please, report this error to the BerkeleyGW developers!'
      endif
      call die("zbssolver1 internal error", only_root_writes=.true.)
    elseif (peinf%inode==0) then
      if (info==0) then
        write(6,'(1x,a,i0,a)') 'Lanczos algorithm performed all ', nmax, ' steps.'
      else
        write(6,'(1x,a,i0,a)') &
          'Lanczos algorithm converged to an invariant subspace after ', info, ' steps.'
      endif
    endif

    if (peinf%inode==0) write(6,'(1x,a)') 'Calling pzlange'
    dtmp = pzlange('F', nmat, 1, s1_bl, 1, 1, desc_d_2d, work)
    call zbsabsp2(nmax_real, n_omega, sigma, omega, eps2, dtmp**2, hx, 2*nmax, lambda)

    if (peinf%inode==0) write(6,'(1x,a)') 'Done with Lanczos iteration.'
    dt = dt + MPI_WTIME()

    SAFE_DEALLOCATE(work)
    SAFE_DEALLOCATE(iwork)

    ! FHJ: TODO: save ritz vectors and values to reuse them it requested
  endif !scal_h%icntxt>=0

  ! FHJ: PEs outside the grid don`t have eps2!
  call MPI_Bcast(eps2, n_omega, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
  call MPI_Bcast(nmax_eff, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

  ! FHJ: Cleanup
  if (scal_h%icntxt>=0) call blacs_gridexit(scal_h%icntxt)
  call blacs_gridexit(cntxt_1d)

  SAFE_DEALLOCATE(hbse_a_bl)
  SAFE_DEALLOCATE(hbse_b_bl)
  SAFE_DEALLOCATE(x_bl)
  SAFE_DEALLOCATE(s1_bl)

  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  if (peinf%inode==0) write(6,'(1x,a)') 'Lanczos step done.'

  POP_SUB(bsabsp_paired)

end subroutine bsabsp_paired

#endif

end module absp_lanczos_m
