!==============================================================================
!
! Routines:
!
! (1) diagonalize()      Originally AC           Last Modified 01/2015 (FHJ)
!
! Generalised parallel Hermitian eigenvalue solver,
! Returns eigenvectors and eigenvalues of the BSE Hamiltonian hbse_a.
!
!        neig = number of eigenvectors
!        nmat = lda of hbse_a and sis (=nv in code)
!
! interface routine for the pzheevx/pdsyevx, generalised parallel eigenvalue
! solver. Starts with the distributed S and H matrices. Redistributes
! them to blacs layout, then calls pzheevx/pdsyevx and
! then collects all the eigenvectors onto
! all pes. For more details of the scalapack routines and data layouts
! see http://www.netlib.org/scalapack/scalapack_home.html
!
! based on pssyevx_inter/pcheevx_inter, written by Andrew Canning
! NERSC/LBL 1998
!
!      distributed solver (block cyclic blacs layout):
!
!           nbl = blocksize
!           nprow = processor grid row
!           npcol = processor grid column
!
! double-precision eigenvalue solvers, pzheevx/pdsyevx
!
!===================================================================================

#include "f_defs.h"
#if defined MPI && !defined USESCALAPACK
  #error ScaLAPACK is required for MPI builds.
#endif

module diagonalize_m

  use global_m
  use scalapack_m
  use scalapack_hl_m
  use elpa_interface_m
  use primme_interface_m
  implicit none

  private

  public :: &
#ifdef USESCALAPACK
    blacs_calc_broken_herm, create_blacs_square_grid, create_blacs_rect_desc, &
#endif
    diagonalize

contains

subroutine diagonalize(xct, neig, nmat, hbse_a, evals, evecs_r, hbse_b, evecs_l)
  type(xctinfo), intent(in) :: xct
  integer, intent(in) :: neig
  integer, intent(in) :: nmat
  SCALAR, intent(inout) :: hbse_a(nmat,peinf%nblocks*peinf%block_sz)
  real(DP), intent(out) :: evals(neig)
  SCALAR, intent(out) :: evecs_r(:,:) !<({nmat,2*nmat},pblock)
  SCALAR, intent(inout), optional :: hbse_b(nmat,peinf%nblocks*peinf%block_sz)
  SCALAR, intent(out), optional :: evecs_l(:,:) !<({nmat,2*nmat},pblock)

  integer :: ii, info
  character :: range
  character*100 :: tmpstr

  PUSH_SUB(diagonalize)

  if (.not.xct%tda.and.(.not.present(hbse_b).or..not.present(evecs_l))) then
    call die('Internal parameter error in diagonalize.', only_root_writes=.true.)
  endif

! If neig > nmat, only the first neig eigenvectors/eigenvalues will
! be computed. Otherwise, calculate all of them.
! Set range
  range = 'A'
  if (xct%tda) then
    if (neig==nmat) then
      range='A'
    else if (neig<nmat) then
      range='I'
    else
      write(tmpstr,'(a,i0,a,i0)') 'diagonalize: ', neig, ' eigenvalues requested > matrix size ', nmat
      call die(tmpstr)
    endif
  endif


  ! FHJ: Full BSE
  if (.not.xct%tda) then
    ! FHJ: Even if we have only 1 processor, we use the ScaLAPACK routines from
    ! SSEIG because it`s faster than our pure LAPACK algorithm.
#ifdef USESCALAPACK
    call blacs_full_sseig()
    !call blacs_full_gvx()
#else
    call serial_full()
    !call serial_full_gvx()
#endif
    POP_SUB(diagonalize)
    return
  endif

  ! FHJ: Serial TDA
  if (peinf%npes==1) then
    call serial_tda()
    POP_SUB(diagonalize)
    return
  endif


  ! FHJ: Parallel TDA below
#ifdef USESCALAPACK

  ! FHJ: PRIMME parallel diag.
#ifdef USEPRIMME
  if (xct%algo==BSE_ALGO_DIAG_PRIMME) then
    call primme_tda()
    POP_SUB(diagonalize)
    return
  endif
#endif

  ! FHJ: ELPA parallel diag.
#ifdef USEELPA
  if (xct%use_elpa) then
    ! FHJ: cannot use ELPA for small block sizes, i.e., when
    ! nblk*(max(np_rows,np_cols)-1) >= na
    ! Block size is currently hard-coded to 64.
    !if (64*int(sqrt(dble(peinf%npes))) >= nmat) then
    select case (int(nmat/sqrt(dble(peinf%npes))))
      case (32:) ! Use ELPA with regular block size of 32
        call blacs_elpa_tda(32)
      case (16:31) ! Use ELPA with smaller block size, useful for tests
        call blacs_elpa_tda(16)
      case (8:15) ! Use ELPA with smaller block size, useful for tests
        call blacs_elpa_tda(8)
      case default ! Don`t bother using ELPA
        if (peinf%inode==0) then
          write(0,'(/a)') 'WARNING: block size too large for this matrix size and process grid!'
          write(0,'(a/)') '         Switching to ScaLAPACK solver.'
        endif
        call blacs_tda()
    endselect
    POP_SUB(diagonalize)
    return
  endif
#endif

  ! FHJ: Generic ScaLAPACK parallel diag.
  !call mpi_tda()
  call blacs_tda()
  POP_SUB(diagonalize)
  return
#endif
! End of USESCALAPACK region

  ! FHJ: If we got here, we had a weird combination of npes>0 but no ScaLAPACK
  call die('Running with more than one processor without ScaLAPACK!')

  POP_SUB(diagonalize)

contains


  !============================================================================
  !> Serial code to diagonalize the TDA Hamiltonian.
  !============================================================================
  subroutine serial_tda()
    real(DP) :: evals_t(nmat)
    real(DP) :: abstol
    integer, allocatable :: iwork(:), ifail(:)
#ifdef CPLX
    real(DP), allocatable :: rwork(:)
#endif
    SCALAR, allocatable :: work(:)
    integer :: iloop, jloop
#ifdef USEESSL
    SCALAR, allocatable :: hupper(:)
    integer :: hup_counter
#endif
    integer :: nfound, ilow, iup, info
    real(DP) :: ellow, elup

    PUSH_SUB(diagonalize.serial_tda)

    if (nmat.ne.peinf%nblocks*peinf%block_sz) then
      call die('Hamiltonian matrix does not seem to be square!')
    endif

    if (peinf%inode==0) write(6,*)
    call calc_broken_herm(hbse_a, nmat, 'H')
    if (peinf%inode==0) write(6,*)

    abstol=0.0
    SAFE_ALLOCATE(ifail, (nmat))
    SAFE_ALLOCATE(work, (10*nmat))
    SAFE_ALLOCATE(iwork, (5*nmat))
    ilow = 1
    iup = neig
#ifdef CPLX
    SAFE_ALLOCATE(rwork, (7*nmat))

#ifdef USEESSL
    SAFE_ALLOCATE(hupper, (nmat*(nmat+1)/2))
    hup_counter = 1
    do jloop=1,nmat
      do iloop=1,jloop
        hupper(hup_counter) =  hbse_a(iloop,jloop)
        hup_counter = hup_counter + 1
      enddo
    enddo
    call zhpev(21,hupper,evals_t,evecs_r,nmat,nmat,iwork,0)
    SAFE_DEALLOCATE(hupper)
#else
    call zheevx('V',range,'U',nmat,hbse_a,nmat,ellow,elup,ilow,iup, &
      abstol,nfound,evals_t,evecs_r,nmat,work,8*nmat,rwork,iwork,ifail,info)
    if(nfound .lt. neig) then
      write(tmpstr,'(a, i10, a, i10, a)') 'Diagonalization with zheevx failed: only ', &
        nfound, ' of ', neig, ' eigenvalues found.'
      call die(tmpstr)
    endif
#endif
    SAFE_DEALLOCATE(rwork)

#else

#ifdef USEESSL
    SAFE_ALLOCATE(hupper, (nmat*(nmat+1)/2))
    hup_counter = 1
    do jloop=1,nmat
      do iloop=1,jloop
        hupper(hup_counter) =  hbse_a(iloop,jloop)
        hup_counter = hup_counter + 1
      enddo
    enddo
    call dspev(21,hupper,evals_t,evecs_r,nmat,nmat,iwork,0)
    SAFE_DEALLOCATE(hupper)
#else
    call dsyevx('V',range,'U',nmat,hbse_a,nmat,ellow,elup,ilow,iup, &
      abstol,nfound,evals_t,evecs_r,nmat,work,10*nmat,iwork,ifail,info)
    if(nfound .lt. neig) then
      write(tmpstr,'(a, i10, a, i10, a)') 'Diagonalization with dsyevx failed: only ', &
        nfound, ' of ', neig, ' eigenvectors found.'
      call die(tmpstr)
    endif
#endif
#endif

    if(info.lt.0) then
      write(tmpstr,*) "Problem in input parameters for zheevx/dsyevx: info = ",info
      call die(tmpstr)
    endif

    if(info.gt.0) then
      write(0,*) "Convergence problems in zheevx/dsyevx: info = ",info
      write(tmpstr,*) 'The following eigenvector failed to converge: ifail = ',ifail
      call die(tmpstr)
    endif

    evals(1:neig) = evals_t(1:neig)
    SAFE_DEALLOCATE(iwork)
    SAFE_DEALLOCATE(work)
    SAFE_DEALLOCATE(ifail)

    POP_SUB(diagonalize.serial_tda)

  end subroutine serial_tda


  !============================================================================
  !> Serial code to diagonalize the full BSE Hamiltonian (Originally by FHJ).
  !============================================================================
  subroutine serial_full()
    SCALAR, allocatable :: work(:)
#ifdef CPLX
    real(DP) :: rwork(4*nmat)
    SCALAR :: evals_t(2*nmat)
#else
    real(DP) :: evals_r(2*nmat), evals_i(2*nmat)
#endif
    SCALAR :: hbse(2*nmat,2*nmat), tmp_norm
    integer :: lwork, info, jj

    PUSH_SUB(diagonalize.serial_full)

    if (nmat.ne.peinf%nblocks*peinf%block_sz) then
      call die('Hamiltonian matrix does not seem to be square!')
    endif

    if (peinf%inode==0) write(6,*)
    call calc_broken_herm(hbse_a, nmat, 'A')
    call calc_broken_herm(hbse_b, nmat, 'B')
    if (peinf%inode==0) write(6,*)

    ! FHJ: symmetrize hbse_a and hbse_b using the lower triangular portions.
    ! We need to to this to obtain the same result as SSEIG solver.
    do ii=1,nmat
      hbse_a(ii, ii+1:nmat) = MYCONJG(hbse_a(ii+1:nmat, ii))
      hbse_b(ii, ii+1:nmat) = hbse_b(ii+1:nmat, ii)
    enddo
    hbse(1:nmat, 1:nmat) = hbse_a(1:nmat, 1:nmat)
    hbse(nmat+1:2*nmat, nmat+1:2*nmat) = -MYCONJG(hbse_a(1:nmat, 1:nmat))
    hbse(1:nmat, nmat+1:2*nmat) = hbse_b(1:nmat, 1:nmat)
    hbse(nmat+1:2*nmat, 1:nmat) = -MYCONJG(hbse_b(1:nmat, 1:nmat))

    SAFE_ALLOCATE(work, (10))

    write(6,'(1x,a,i0)') 'Beginning LAPACK diagonalization. Size: ', 2*nmat
    call X(geev)('V', 'V', 2*nmat, hbse, 2*nmat, &
#ifdef CPLX
      evals_t, &
#else
      evals_r, evals_i, &
#endif
      evecs_l, 2*nmat, evecs_r, 2*nmat, work, -1, &
#ifdef CPLX
      rwork, &
#endif
      info)

    if(info/=0) then
      write(tmpstr,*) "problem in xgeev/query mode: got info = ", info
      call die(tmpstr)
    endif

    lwork = max(1, int(work(1)))
    SAFE_DEALLOCATE(work)
    SAFE_ALLOCATE(work, (lwork))

    call X(geev)('V', 'V', 2*nmat, hbse, 2*nmat, &
#ifdef CPLX
      evals_t, &
#else
      evals_r, evals_i, &
#endif
      evecs_l, 2*nmat, evecs_r, 2*nmat, work, lwork, &
#ifdef CPLX
      rwork, &
#endif
      info)

    if(info/=0) then
      write(tmpstr,*) "problem in xgeev: got info = ", info
      call die(tmpstr)
    endif

    write(6,*) 'Done with LAPACK diagonalization.'

#ifdef CPLX
    evals(1:neig) = dble(evals_t(1:neig))
#else
    evals(1:neig) = evals_r(1:neig)
#endif

    ! FHJ: No, the left/right eigenvectors that come out of LAPACK are *not*
    ! correctly normalized against each other!
    do jj=1,2*nmat
      tmp_norm = sum(evecs_l(1:2*nmat,jj)*evecs_r(1:2*nmat,jj))
      evecs_l(1:2*nmat,jj) = evecs_l(1:2*nmat,jj) / tmp_norm
    enddo

    PUSH_SUB(diagonalize.serial_full)

  end subroutine serial_full


  !============================================================================
  !> Serial code to diagonalize the full BSE Hamiltonian (Orig. by FHJ).
  !! This version calls a generalized eigensolver from LAPACK.
  !============================================================================
  subroutine serial_full_gvx()
    SCALAR :: hbse(2*nmat,2*nmat)
    integer :: lwork, liwork, lrwork, nfound, ii
    integer :: ifail(2*nmat), iwork(10*nmat)
    SCALAR :: b_mat(2*nmat,2*nmat), s_mat(2*nmat,2*nmat)
    SCALAR, allocatable :: work(:)
#ifdef CPLX
    real(DP) :: rwork(14*nmat)
#endif
    real(DP) :: abstol, evals_t(2*nmat)

    PUSH_SUB(diagonalize.serial_full_gvx)

    ! FHJ: zero coupling block if the user requests. Useful to test the solver.
    if (xct%zero_coupling_block) hbse_b = ZERO

    if (nmat.ne.peinf%nblocks*peinf%block_sz) then
      call die('Hamiltonian matrix does not seem to be square!')
    endif


    ! FHJ: symmetrize hbse_a and hbse_b using the lower triangular portions.
    ! We need to to this to obtain the same result as SSEIG solver.
    do ii=1,nmat
      hbse_a(ii, ii+1:nmat) = MYCONJG(hbse_a(ii+1:nmat, ii))
      hbse_b(ii, ii+1:nmat) = hbse_b(ii+1:nmat, ii)
    enddo
    hbse(1:nmat, 1:nmat) = hbse_a(1:nmat, 1:nmat)
    hbse(nmat+1:2*nmat, nmat+1:2*nmat) = MYCONJG(hbse_a(1:nmat, 1:nmat))
    hbse(1:nmat, nmat+1:2*nmat) = hbse_b(1:nmat, 1:nmat)
    hbse(nmat+1:2*nmat, 1:nmat) = MYCONJG(hbse_b(1:nmat, 1:nmat))

    ! Create matrix B
    ! b_mat = [ 1  0 ]
    !         [ 0 -1 ]
    call X(laset)('A', nmat, nmat, ZERO, ONE, b_mat(1,1), 2*nmat)
    call X(laset)('A', nmat, nmat, ZERO, -ONE, b_mat(nmat+1,nmat+1), 2*nmat)

    ! FHJ: Call *evx in query mode
    abstol = 0.0d0
    SAFE_ALLOCATE(work, (10))

#ifdef CPLX
    call zhegvx &
#else
    call dsygvx &
#endif
      (2, 'V', range, 'U', 2*nmat, b_mat, 2*nmat, hbse, 2*nmat, &
      0d0, 0d0, 1, neig, abstol, nfound, &
      evals_t, evecs_r, 2*nmat, work, -1, &
#ifdef CPLX
      rwork, &
#endif
      iwork, ifail, info)
    if (info/=0) then
      !if(peinf%inode==0) then
        write(0,'(/,a,i0,/)') 'ERROR: Query mode for ???gvx failed with info=', info
      !endif
      call die("???gvx query failed", only_root_writes=.true.)
    endif

    lwork = max(1,int(work(1))) + 10*nmat
    SAFE_DEALLOCATE(work)
    SAFE_ALLOCATE(work, (lwork))

    ! FHJ: Call p*evx for realz
    if (peinf%inode==0) write(6,'(1x,a,i0)') 'Beginning LAPACK diagonalization. Size: ', 2*nmat
#ifdef CPLX
    call zhegvx &
#else
    call dsygvx &
#endif
      (2, 'V', range, 'U', 2*nmat, b_mat, 2*nmat, hbse, 2*nmat, &
      0d0, 0d0, 1, neig, abstol, nfound, &
      evals_t, evecs_r, 2*nmat, work, lwork, &
#ifdef CPLX
      rwork, &
#endif
      iwork, ifail, info)

    if (peinf%inode==0) write(6,*) 'Done LAPACK diagonalization'
    SAFE_DEALLOCATE(work)

    if(nfound<neig) then
      if (peinf%inode==0) then
        write(0,'(/,a)') 'ERROR: Diagonalization with zhegvx/dsygvx failed:'
        write(0,'(3(a,i0),/)') 'only ', nfound, ' out of ', neig, ' eigenvalues found.'
      endif
        write(0,'(/,a,i0,/)') 'ERROR: Diagonalization with ???evx failed with info=',info
        write(0,'(3(a,i0),/)') 'only ', nfound, ' out of ', neig, ' eigenvalues found.'
      call die("LAPACK found wrong number of eigenvalues", only_root_writes=.true.)
    endif
    if (info/=0.and.info/=2) then
      !if(peinf%inode==0) then
        write(0,'(/,a,i0,/)') 'ERROR: Diagonalization with ???evx failed with info=',info
      !endif
      call die("p???evx diagonalization failed", only_root_writes=.true.)
    endif

    ! FHJ: Copy eigenvalues/eigenvectors
    evals(1:neig) = evals_t(1:neig)

    if (peinf%inode==0) write(6,*) 'Calculating overlap matrix'
    ! FHJ: Calculate overlap S = evecs_r evecs_r^H.
    !      Note: other codes/works define S by evecs_r^H evecs_r
#ifdef CPLX
    call zherk&
#else
    call dsyrk&
#endif
      ('L', 'N', 2*nmat, 2*nmat, ONE, evecs_r, 2*nmat, ZERO, s_mat, 2*nmat)

    if (peinf%inode==0) write(6,*) 'Performing Cholesky decomposition'
    ! FHJ: Invert overlap matrix. First, do Cholesky decomposition.
    call X(potrf)('L', 2*nmat, s_mat, 2*nmat, info)
    if (info/=0) then
      if (peinf%inode==0) write(0,*) 'ERROR: got info=', info
      call die('Cholesky decomposition failed', only_root_writes=.true.)
    endif
    if (peinf%inode==0) write(6,*) 'Inverting overlap matrix'
    ! FHJ: Now, invert the matrix
    call X(potri)('L', 2*nmat, s_mat, 2*nmat, info)
    if (info/=0) then
      if (peinf%inode==0) write(0,*) 'ERROR: got info=', info
      call die('matrix inversion failed', only_root_writes=.true.)
    endif
    if (peinf%inode==0) write(6,*) 'Calculating left evecs'
    ! FHJ: Multiply S^{-1} by evecs_r to get evecs_l^H
#ifdef CPLX
    call zhemm&
#else
    call dsymm&
#endif
      ('L', 'L', 2*nmat, 2*nmat, ONE, s_mat, 2*nmat, &
      evecs_r, 2*nmat, ZERO, evecs_l, 2*nmat)

    ! FHJ: Apply complex conjugation to get evecs_l^T
    evecs_l = MYCONJG(evecs_l)

    if (peinf%inode==0) write(6,*) 'Diagonalization done'

    POP_SUB(diagonalize.serial_full_gvx)

  end subroutine serial_full_gvx


#ifdef USEELPA
  !============================================================================
  !> Parallel code to diagonalize the TDA Hamiltonian (Originally by FHJ).
  !! It uses BLACS (P?GEMR2D) to change the matrix distribution and ELPA
  !============================================================================
  subroutine blacs_elpa_tda(block_sz)
    integer, intent(in) :: block_sz

    ! Scalapack and blacs arrays
    type(scalapack) :: scal
    integer :: desc_1d(9), desc_2d(9)
    SCALAR, allocatable :: hbse_a_bl(:), evecs_r_bl(:)
    real(dp), allocatable :: evals_all(:)
    integer :: cntxt_1d, usermap(1,peinf%npes), lld
    integer :: comm, inode, npes, color

    PUSH_SUB(diagonalize.blacs_elpa_tda)

    ! FHJ: Create BLACS context for 1d cyclic matrix (hbse_a). Processors: all aboard!
    call blacs_get(-1, 0, cntxt_1d)
    do ii=1,peinf%npes
      usermap(1,ii) = ii-1
    enddo
    call blacs_gridmap(cntxt_1d, usermap, 1, 1, peinf%npes) ! This will modify cntxt_1d
    call descinit(desc_1d, nmat, nmat, nmat, peinf%block_sz, 0, 0, cntxt_1d, nmat, info)
    if (info/=0) call die('got info/=0 in descinit')

    ! FHJ: Initialize BLACS grid for 2d cyclic matrices hbse_a_bl, evecs_r_bl.
    ! Some processors might get excluded for the sake of load balancing.
    call create_blacs_square_grid(scal, block_sz)
    call create_blacs_rect_desc(nmat, nmat, scal, desc_2d, lld)

    SAFE_ALLOCATE(hbse_a_bl, (lld*max(1,scal%npc)))
    hbse_a_bl(:) = ZERO
    SAFE_ALLOCATE(evecs_r_bl, (lld*max(1,scal%npc)))
    evecs_r_bl(:) = ZERO

    ! FHJ: symmetrize the matrix. This is important for consistency: ScaLAPACK
    ! only reference the upper tringular part, but not ELPA.
    call blacs_symmetrize_u2l(hbse_a, nmat, size(hbse_a), desc_1d, .true.)

    ! FHJ: Copy matrix hbse_r from 1d cyclic layout into hbse_a_bl (2d cyclic layout)
    call pX(gemr2d)(nmat, nmat, hbse_a, 1, 1, desc_1d, hbse_a_bl, 1, 1, desc_2d, cntxt_1d)

    ! Split the MPI communicator and only include the processors with data.
    ! If inode==-1, the processor is not in the group.
    color = 1
    if (scal%icntxt<0) color = 2
    inode = peinf%inode
    call MPI_Comm_split(MPI_COMM_WORLD, color, inode, comm, mpierr)
    if (peinf%inode==0) call MPI_Comm_size(comm, npes, mpierr)
    call MPI_Bcast(npes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    if (color==2) inode = -1

    ! MDB need to allocate all eigenvalues to avoid invalid memory reference in ELPA if
    ! number_eigenvalues is given in input smaller than nmat
    SAFE_ALLOCATE( evals_all, (nmat))
    evals_all = 0.0D+00

    ! FHJ: Call ELPA
    if (scal%icntxt<0) then
      call diagonalize_elpa(scal, nmat, neig, hbse_a_bl, evals_all, evecs_r_bl, comm=MPI_COMM_NULL)
    else
      call diagonalize_elpa(scal, nmat, neig, hbse_a_bl, evals_all, evecs_r_bl, comm=comm)
    endif

    ! MDB copy the eigenvalues back
    evals(1:neig) = evals_all(1:neig)
    SAFE_DEALLOCATE( evals_all )

    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    ! FHJ: PEs outside the grid don`t have the evals!
    call MPI_Bcast(evals, neig, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
    if (peinf%inode==0) write(6,*) 'Sharing eigenvectors'
    call pX(gemr2d)(nmat, neig, evecs_r_bl, 1, 1, desc_2d, evecs_r, 1, 1, desc_1d, cntxt_1d)

    ! FHJ: Cleanup
    if (scal%icntxt>=0) call blacs_gridexit(scal%icntxt)
    call blacs_gridexit(cntxt_1d)
    SAFE_DEALLOCATE(hbse_a_bl)
    SAFE_DEALLOCATE(evecs_r_bl)
    call MPI_Comm_free(comm, mpierr)
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    if (peinf%inode==0) write(6,*) 'Diagonalization done'

    POP_SUB(diagonalize.blacs_elpa_tda)

  end subroutine blacs_elpa_tda
#endif

!#BEGIN_INTERNAL_ONLY
#ifdef USEPRIMME
  subroutine primme_tda()
    SCALAR, allocatable :: evecs_primme(:,:), ham_primme(:,:)
    SCALAR :: evec_column(nmat)
    integer :: Nb, Nown_bse, Nown_primme, ieval, ipe, ieval_loc, ii_l, ii_g
    !integer :: recvcounts(peinf%npes), displs(peinf%npes)
    integer :: desc_bse(9), desc_primme(9), desc_evecs_primme(9)
    integer :: cntxt_1d, cntxt_1d_r, usermap(1,peinf%npes)

    PUSH_SUB(diagonalize.primme_tda)

    ! WARNING: peinf%nblocks is the max. number of blocks each PE owns
    ! peinf%ibt(peinf%inode+1) == Nown_bse is the number of blocks the current
    ! PE actually owns.
    Nb = DIVUP(nmat, peinf%npes)
    Nown_primme = NUMROC(nmat, Nb, peinf%inode, 0, peinf%npes)
    Nown_bse = NUMROC(nmat, peinf%block_sz, peinf%inode, 0, peinf%npes)
    !print '(a,5(1x,i0))', '!!!: ', peinf%inode, Nown_bse, Nown_primme, &
    !  peinf%block_sz*peinf%nblocks, peinf%block_sz*peinf%ibt(peinf%inode+1)
    if (Nown_bse/=peinf%block_sz*peinf%ibt(peinf%inode+1)) then
      print '(a,5(1x,i0))', 'ERROR: ', peinf%inode, Nown_bse, Nown_primme, &
        peinf%block_sz*peinf%nblocks, peinf%block_sz*peinf%ibt(peinf%inode+1)
      call die('Inconsistency betwen Nown_bse and peinf%ibt')
    endif

    ! Distribute original Hamiltonian from block size peinf%block_sz to Nb.
    SAFE_ALLOCATE(ham_primme, (nmat,max(1,Nown_primme)))
    ham_primme(:,:) = ZERO

    ! FHJ: Create BLACS context for 1d column-cyclic distributed matrices
    call blacs_get(-1, 0, cntxt_1d)
    do ii=1,peinf%npes
      usermap(1,ii) = ii-1
    enddo
    call blacs_gridmap(cntxt_1d, usermap, 1, 1, peinf%npes)

    ! FHJ: Create descriptor for original distribution of hbse_a with block
    ! size peinf%block_sz
    call descinit(desc_bse, nmat, nmat, nmat, peinf%block_sz, 0, 0, cntxt_1d, nmat, info)
    if (info/=0) call die('got info/=0 in descinit')

    ! FHJ: symmetrize matrix. This is very important for PRIMME!
    if (peinf%inode==0) write(6,*)
    call blacs_calc_broken_herm(hbse_a, nmat, size(hbse_a), desc_bse, 'H')
    if (peinf%inode==0) write(6,*)
    call blacs_symmetrize_u2l(hbse_a, nmat, size(hbse_a), desc_bse, .true.)

    ! FHJ: Create descriptor for disribution of ham_primme with block size Nb
    call descinit(desc_primme, nmat, nmat, nmat, Nb, 0, 0, cntxt_1d, nmat, info)
    if (info/=0) call die('got info/=0 in descinit')

    ! FHJ: Copy matrix hbse_a with block size peinf%block_sz to ham_primme, with
    ! block size Nb
    if (peinf%inode==0) write(6,*) 'Redistributing Hamiltonian'
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    call pX(gemr2d)(nmat, nmat, hbse_a, 1, 1, desc_bse, ham_primme, 1, 1, desc_primme, cntxt_1d)

    ! Note that evecs_primme is distributed over the first index (vck) instead
    ! of over the eigenvalue index.
    SAFE_ALLOCATE(evecs_primme, (max(1,Nown_primme),neig))

    ! Call PRIMME
    if (peinf%inode==0) write(6,*) 'Calling PRIMME solver interface'
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    call solve_primme_explicit(nmat, Nown_primme, neig, evals, evecs_primme, &
      ham_primme, max_block_size=xct%primme_max_block_size, &
      max_basis_size=xct%primme_max_basis_size, tol=xct%primme_tol)

    ! Redistribute eigenvectors
    ! First, create BLACS context for block row distributed matrices
    call blacs_get(-1, 0, cntxt_1d_r)
    call blacs_gridmap(cntxt_1d_r, transpose(usermap), peinf%npes, peinf%npes, 1)

    ! Then, create descriptor for eigenvectors computed by PRIMME
    call descinit(desc_evecs_primme, nmat, nmat, Nb, nmat, 0, 0, cntxt_1d_r, max(1,Nown_primme), info)
    if (info/=0) call die('got info/=0 in descinit')

    ! FHJ: Copy matrix evecs_primme to evecs_a
    ! The final distribution of evecs_a is the same as hbse_a
    if (peinf%inode==0) write(6,*) 'Redistributing eigenvectors'
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    call pX(gemr2d)(nmat, neig, evecs_primme, 1, 1, desc_evecs_primme, &
      evecs_r, 1, 1, desc_bse, cntxt_1d)

    if (peinf%inode==0) write(6,*) 'Freeing BLACS context structures'
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    if (cntxt_1d>=0) call blacs_gridexit(cntxt_1d)
    if (cntxt_1d_r>=0) call blacs_gridexit(cntxt_1d_r)

    SAFE_DEALLOCATE(evecs_primme)
    SAFE_DEALLOCATE(ham_primme)

    POP_SUB(diagonalize.primme_tda)

  end subroutine primme_tda
#endif
!#END_INTERNAL_ONLY

#ifdef USESCALAPACK
  !============================================================================
  !> Parallel code to diagonalize the TDA Hamiltonian (Originally by FHJ and MS).
  !! It uses BLACS (P?GEMR2D) to change the matrix distribution.
  !============================================================================
  subroutine blacs_tda()
    ! Scalapack and blacs arrays
    type(scalapack) :: scal
    integer :: desc_1d(9), desc_2d(9)
    SCALAR, allocatable :: hbse_a_bl(:), evecs_r_bl(:)
    integer :: cntxt_1d, usermap(1,peinf%npes), lld

    PUSH_SUB(diagonalize.blacs_tda)

    ! FHJ: Create BLACS context for 1d cyclic matrix (hbse_a). Processors: all aboard!
    call blacs_get(-1, 0, cntxt_1d)
    do ii=1,peinf%npes
      usermap(1,ii) = ii-1
    enddo
    call blacs_gridmap(cntxt_1d, usermap, 1, 1, peinf%npes) ! This will modify cntxt_1d
    call descinit(desc_1d, nmat, nmat, nmat, peinf%block_sz, 0, 0, cntxt_1d, nmat, info)
    if (info/=0) call die('got info/=0 in descinit')

    ! FHJ: Initialize BLACS grid for 2d cyclic matrices hbse_a_bl, evecs_r_bl.
    ! Some processors might get excluded for the sake of load balancing.
    call create_blacs_square_grid(scal)
    call create_blacs_rect_desc(nmat, nmat, scal, desc_2d, lld)

    SAFE_ALLOCATE(hbse_a_bl, (lld*max(1,scal%npc)))
    hbse_a_bl(:) = ZERO
    SAFE_ALLOCATE(evecs_r_bl, (lld*max(1,scal%npc)))
    evecs_r_bl(:) = ZERO

    ! FHJ: Copy matrix hbse_r from 1d cyclic layout into hbse_a_bl (2d cyclic layout)
    call pX(gemr2d)(nmat, nmat, hbse_a, 1, 1, desc_1d, hbse_a_bl, 1, 1, desc_2d, cntxt_1d)

    ! FHJ: Call ScaLAPACK
    call diagonalize_scalapack(desc_2d, neig, hbse_a_bl, evals, evecs_r_bl)

    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    ! FHJ: PEs outside the grid don`t have the evals!
    call MPI_Bcast(evals, neig, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
    if (peinf%inode==0) write(6,*) 'Sharing eigenvectors'
    call pX(gemr2d)(nmat, neig, evecs_r_bl, 1, 1, desc_2d, evecs_r, 1, 1, desc_1d, cntxt_1d)

    ! FHJ: Cleanup
    if (scal%icntxt>=0) call blacs_gridexit(scal%icntxt)
    call blacs_gridexit(cntxt_1d)
    SAFE_DEALLOCATE(hbse_a_bl)
    SAFE_DEALLOCATE(evecs_r_bl)
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    if (peinf%inode==0) write(6,*) 'Diagonalization done'

    POP_SUB(diagonalize.blacs_tda)

  end subroutine blacs_tda


  !============================================================================
  !> Parallel code to diagonalize the full BSE Hamiltonian (Orig. by FHJ).
  !! It uses BLACS (P?GEMR2D) to change the matrix distribution.
  !! This version calls a generalized eigensolver from ScaLAPACK and does
  !! some manual work to get the left eigenvectors. It also doesn`t respect
  !! the symmetries of the BSE Hamiltonian, but it doesn`t require SSEIG.
  !! NOTE: This routine gives slightly different results than blacs_full_sseig
  !! because we don`t symmetrize the matrix.
  !============================================================================
  subroutine blacs_full_gvx()
    ! Scalapack and blacs arrays
    type(scalapack) :: scal
    integer, dimension(9) :: desc_h_1d, desc_e_1d, desc_2d
    integer :: lwork, liwork, lrwork, nfound, nzfound
    integer, allocatable :: iclustr(:), iwork(:), ifail(:)
    SCALAR, allocatable :: hbse_bl(:), evecs_r_bl(:), evecs_l_bl(:)
    SCALAR, allocatable :: b_bl(:), s_bl(:)
    SCALAR, allocatable :: work(:)
#ifdef CPLX
    real(DP), allocatable :: rwork(:)
#endif
    integer :: cntxt_1d,  usermap(1,peinf%npes), lld
    real(DP), allocatable :: gap(:)
    real(DP) :: abstol, orfac, evals_t(2*nmat)

    PUSH_SUB(diagonalize.blacs_full_gvx)

    ! FHJ: Create BLACS context for 1d cyclic matrices. Processors: all aboard!
    call blacs_get(-1, 0, cntxt_1d)
    do ii=1,peinf%npes
      usermap(1,ii) = ii-1
    enddo
    call blacs_gridmap(cntxt_1d, usermap, 1, 1, peinf%npes) ! This will modify cntxt_1d
    call descinit(desc_h_1d, nmat, nmat, nmat, peinf%block_sz, 0, 0, cntxt_1d, nmat, info)
    if (info/=0) call die('got info/=0 in descinit')
    call descinit(desc_e_1d, 2*nmat, 2*nmat, 2*nmat, peinf%block_sz, 0, 0, cntxt_1d, 2*nmat, info)
    if (info/=0) call die('got info/=0 in descinit')

    ! FHJ: Initialize BLACS grid for 2d cyclic matrices hbse_bl, evecs_{r,l}_bl.
    ! Some processors might get excluded for the sake of load balancing.
    call create_blacs_square_grid(scal)
    call create_blacs_rect_desc(2*nmat, 2*nmat, scal, desc_2d, lld)

    SAFE_ALLOCATE(hbse_bl, (lld*max(1,scal%npc)))
    hbse_bl(:) = ZERO
    SAFE_ALLOCATE(b_bl, (lld*max(1,scal%npc)))
    b_bl(:) = ZERO
    SAFE_ALLOCATE(s_bl, (lld*max(1,scal%npc)))
    s_bl(:) = ZERO
    SAFE_ALLOCATE(evecs_r_bl, (lld*max(1,scal%npc)))
    evecs_r_bl(:) = ZERO
    SAFE_ALLOCATE(evecs_l_bl, (lld*max(1,scal%npc)))
    evecs_l_bl(:) = ZERO

    if (peinf%inode==0) write(6,*)
    call blacs_calc_broken_herm(hbse_a, nmat, size(hbse_a), desc_h_1d, 'A')
    call blacs_calc_broken_herm(hbse_b, nmat, size(hbse_b), desc_h_1d, 'B')
    if (peinf%inode==0) write(6,*)

    ! FHJ: zero coupling block if the user requests. Useful to test the solver.
    if (xct%zero_coupling_block) hbse_b = ZERO

    ! FHJ: Copy matrices hbse_{a,b} from 1d cyclic layout into hbse_bl (2d cyclic layout)
    ! hbse_bl = [ A  B  ]
    !           [ B* A* ]
    if (peinf%inode==0) write(6,*) 'Preparing matrices'
    call blacs_symmetrize_l2u(hbse_a, nmat, size(hbse_a), desc_h_1d, .true.)
    call blacs_symmetrize_l2u(hbse_b, nmat, size(hbse_b), desc_h_1d, .false.)
    call pX(gemr2d)(nmat, nmat, hbse_a, 1, 1, desc_h_1d, hbse_bl, 1, 1, desc_2d, cntxt_1d)
    call pX(gemr2d)(nmat, nmat, hbse_b, 1, 1, desc_h_1d, hbse_bl, 1, nmat+1, desc_2d, cntxt_1d)
    hbse_a = MYCONJG(hbse_a)
    hbse_b = MYCONJG(hbse_b)
    call pX(gemr2d)(nmat, nmat, hbse_a, 1, 1, desc_h_1d, hbse_bl, nmat+1, nmat+1, desc_2d, cntxt_1d)
    call pX(gemr2d)(nmat, nmat, hbse_b, 1, 1, desc_h_1d, hbse_bl, nmat+1, 1, desc_2d, cntxt_1d)

    if (scal%icntxt>=0) then
      ! Create distributed matrix B
      ! b_bl = [ 1  0 ]
      !        [ 0 -1 ]
      call pX(laset)('A', nmat, nmat, ZERO, ONE, b_bl, 1, 1, desc_2d)
      call pX(laset)('A', nmat, nmat, ZERO, -ONE, b_bl, nmat+1, nmat+1, desc_2d)

      ! FHJ: Call p*evx in query mode
      abstol = 0.0d0
      orfac = 5d-7
      SAFE_ALLOCATE(iclustr, (2*scal%nprow*scal%npcol))
      SAFE_ALLOCATE(gap, (scal%nprow*scal%npcol))
      SAFE_ALLOCATE(ifail, (2*nmat))
      SAFE_ALLOCATE(work, (10))
#ifdef CPLX
      SAFE_ALLOCATE(rwork, (10))
#endif
      SAFE_ALLOCATE(iwork, (10))

#ifdef CPLX
      call pzhegvx &
#else
      call pdsygvx &
#endif
        (2, 'V', range, 'U', 2*nmat, b_bl, 1, 1, desc_2d, &
        hbse_bl, 1, 1, desc_2d, &
        0d0, 0d0, 1, neig, abstol, nfound, &
        nzfound, evals_t, orfac, evecs_r_bl, 1, 1, desc_2d, work, -1, &
#ifdef CPLX
        rwork, -1, &
#endif
        iwork, -1, ifail, iclustr, gap, info)
      if (info/=0) then
        !if(peinf%inode==0) then
          write(0,'(/,a,i0,/)') 'ERROR: Query mode for p???evx failed with info=', info
        !endif
        call die("p???evx query failed", only_root_writes=.true.)
      endif

#ifdef CPLX
      lrwork = max(3,int(rwork(1))) + (2*nmat)*10
      SAFE_DEALLOCATE(rwork)
      SAFE_ALLOCATE(rwork, (lrwork))
      lwork = max(3,int(work(1)))
      SAFE_DEALLOCATE(work)
      SAFE_ALLOCATE(work, (lwork))
#else
      lwork = max(3,int(work(1))) + (2*nmat)*10
      SAFE_DEALLOCATE(work)
      SAFE_ALLOCATE(work, (lwork))
#endif
      liwork = iwork(1)
      SAFE_DEALLOCATE(iwork)
      SAFE_ALLOCATE(iwork, (liwork))

      ! FHJ: Call p*evx for realz
      if (peinf%inode==0) write(6,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', 2*nmat
#ifdef CPLX
      call pzhegvx &
#else
      call pdsygvx &
#endif
        (2, 'V', range, 'U', 2*nmat, b_bl, 1, 1, desc_2d, &
        hbse_bl, 1, 1, desc_2d, &
        0d0, 0d0, 1, neig, abstol, nfound, &
        nzfound, evals_t, orfac, evecs_r_bl, 1, 1, desc_2d, work, lwork, &
#ifdef CPLX
        rwork, lrwork, &
#endif
        iwork, liwork, ifail, iclustr, gap, info)
      if (peinf%inode==0) write(6,*) 'Done ScaLAPACK diagonalization'

      SAFE_DEALLOCATE(work)
#ifdef CPLX
      SAFE_DEALLOCATE(rwork)
#endif
      SAFE_DEALLOCATE(iwork)

      if(nfound<neig) then
        if (peinf%inode==0) then
          write(0,'(/,a)') 'ERROR: Diagonalization with pzheevx/pdsyevx failed:'
          write(0,'(3(a,i0),/)') 'only ', nfound, ' out of ', neig, ' eigenvalues found.'
        endif
          write(0,'(/,a,i0,/)') 'ERROR: Diagonalization with p???evx failed with info=',info
          write(0,'(3(a,i0),/)') 'only ', nfound, ' out of ', neig, ' eigenvalues found.'
        call die("ScaLAPACK found wrong number of eigenvalues", only_root_writes=.true.)
      endif
      if(nzfound<neig) then
        if (peinf%inode==0) then
          write(0,'(/,a)') 'ERROR: Diagonalization with pzheevx/pdsyevx failed:'
          write(0,'(3(a,i0),/)') 'only ', nzfound, ' out of ', neig, ' eigenvectors found.'
        endif
        call die("ScaLAPACK found wrong number of eigenvectors", only_root_writes=.true.)
      endif
      if (info/=0.and.info/=2) then
        !if(peinf%inode==0) then
          write(0,'(/,a,i0,/)') 'ERROR: Diagonalization with p???evx failed with info=',info
        !endif
        call die("p???evx diagonalization failed", only_root_writes=.true.)
      endif

      ! FHJ: Copy eigenvalues/eigenvectors
      evals(1:neig) = evals_t(1:neig)
      SAFE_DEALLOCATE(iclustr)
      SAFE_DEALLOCATE(gap)
      SAFE_DEALLOCATE(ifail)

      if (peinf%inode==0) write(6,*) 'Calculating overlap matrix'
      ! FHJ: Calculate overlap S = evecs_r evecs_r^H.
      !      Note: other codes/works define S by evecs_r^H evecs_r
#ifdef CPLX
      call pzherk&
#else
      call pdsyrk&
#endif
        ('L', 'N', 2*nmat, 2*nmat, ONE, &
          evecs_r_bl, 1, 1, desc_2d, ZERO, &
          s_bl, 1, 1, desc_2d)

      if (peinf%inode==0) write(6,*) 'Performing Cholesky decomposition'
      ! FHJ: Invert overlap matrix. First, do Cholesky decomposition.
      call pX(potrf)('L', 2*nmat, s_bl, 1, 1, desc_2d, info)
      if (info/=0) then
        if (peinf%inode==0) write(0,*) 'ERROR: got info=', info
        call die('Cholesky decomposition failed', only_root_writes=.true.)
      endif
      if (peinf%inode==0) write(6,*) 'Inverting overlap matrix'
      ! FHJ: Now, invert the matrix
      call pX(potri)('L', 2*nmat, s_bl, 1, 1, desc_2d, info)
      if (info/=0) then
        if (peinf%inode==0) write(0,*) 'ERROR: got info=', info
        call die('matrix inversion failed', only_root_writes=.true.)
      endif
      if (peinf%inode==0) write(6,*) 'Calculating left evecs'
      ! FHJ: Multiply S^{-1} by evecs_r to get evecs_l^H
#ifdef CPLX
      call pzhemm&
#else
      call pdsymm&
#endif
        ('L', 'L', 2*nmat, 2*nmat, ONE, &
        s_bl, 1, 1, desc_2d, &
        evecs_r_bl, 1, 1, desc_2d, ZERO, &
        evecs_l_bl, 1, 1, desc_2d)

      ! FHJ: Apply complex conjugation to get evecs_l^T
      evecs_l_bl = MYCONJG(evecs_l_bl)

    endif !scal%icntxt>=0

    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    ! FHJ: PEs outside the grid don`t have the evals!
    call MPI_Bcast(evals, neig, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
    if (peinf%inode==0) write(6,*) 'Sharing eigenvectors'
    call pX(gemr2d)(2*nmat, 2*nmat, evecs_r_bl, 1, 1, desc_2d, evecs_r, 1, 1, desc_e_1d, cntxt_1d)
    call pX(gemr2d)(2*nmat, 2*nmat, evecs_l_bl, 1, 1, desc_2d, evecs_l, 1, 1, desc_e_1d, cntxt_1d)

    ! FHJ: Cleanup
    if (scal%icntxt>=0) call blacs_gridexit(scal%icntxt)
    call blacs_gridexit(cntxt_1d)
    SAFE_DEALLOCATE(hbse_bl)
    SAFE_DEALLOCATE(evecs_r_bl)
    SAFE_DEALLOCATE(evecs_l_bl)
    SAFE_DEALLOCATE(s_bl)
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    if (peinf%inode==0) write(6,*) 'Diagonalization done'

    POP_SUB(diagonalize.blacs_full_gvx)

  end subroutine blacs_full_gvx


  !============================================================================
  !> Parallel code to diagonalize the full BSE Hamiltonian (Orig. by FHJ and MS).
  !! It uses BLACS (P?GEMR2D) to change the matrix distribution.
  !============================================================================
  subroutine blacs_full_sseig()
    ! Scalapack and blacs arrays
    type(scalapack) :: scal_h, scal_e
    integer, dimension(9) :: desc_h_1d, desc_e_1d, desc_h_2d, desc_e_2d
    integer :: lwork, liwork
#ifdef CPLX
    integer :: lrwork
    real(DP), allocatable :: rwork(:)
#endif
    SCALAR, allocatable :: hbse_a_bl(:), hbse_b_bl(:), evecs_r_bl(:)
    complex(DPC), allocatable :: work(:)
    integer, allocatable :: iwork(:)
    real(DP) :: evals_t(nmat)
    integer :: cntxt_1d,  usermap(1,peinf%npes), lld_h, lld_e

    include 'full_solver/solver.f90'

    PUSH_SUB(diagonalize.blacs_full_sseig)

    ! FHJ: Create BLACS context for 1d cyclic matrices. Processors: all aboard!
    call blacs_get(-1, 0, cntxt_1d)
    do ii=1,peinf%npes
      usermap(1,ii) = ii-1
    enddo
    call blacs_gridmap(cntxt_1d, usermap, 1, 1, peinf%npes) ! This will modify cntxt_1d
    call descinit(desc_h_1d, nmat, nmat, nmat, peinf%block_sz, 0, 0, cntxt_1d, nmat, info)
    if (info/=0) call die('got info/=0 in descinit')
    ! FHJ: Note that the eigenvectors are 2*nmat by nmat!
    call descinit(desc_e_1d, 2*nmat, nmat, 2*nmat, peinf%block_sz, 0, 0, cntxt_1d, 2*nmat, info)
    if (info/=0) call die('got info/=0 in descinit')

    ! FHJ: Initialize BLACS grid for 2d cyclic matrices hbse_{a,b}_bl.
    ! Some processors might get excluded for the sake of load balancing.
    call create_blacs_square_grid(scal_h)
    call create_blacs_rect_desc(nmat, nmat, scal_h, desc_h_2d, lld_h)

    ! FHJ: Initialize BLACS grid for 2d cyclic matrices evecs_{r,l}_bl.
    ! The processor grid doesn`t change, so we reuse the context from scal_h%icntxt.
    scal_e = scal_h
    ! FHJ: There`s a small bug in PZBSEIG, as it actually require the descriptor
    ! of X to be (2N,2N).
    call create_blacs_rect_desc(2*nmat, nmat, scal_e, desc_e_2d, lld_e)
    !call create_blacs_rect_desc(2*nmat, 2*nmat, scal_e, desc_e_2d, lld_e)

    SAFE_ALLOCATE(hbse_a_bl, (lld_h*max(1,scal_h%npc)))
    hbse_a_bl(:) = ZERO
    SAFE_ALLOCATE(hbse_b_bl, (lld_h*max(1,scal_h%npc)))
    hbse_b_bl(:) = ZERO
    SAFE_ALLOCATE(evecs_r_bl, (lld_e*max(1,scal_e%npc)))
    evecs_r_bl(:) = ZERO

    ! FHJ: Copy matrices hbse_{a,b} from 1d cyclic layout into hbse_{a,b}_bl (2d cyclic layout)
    call pX(gemr2d)(nmat, nmat, hbse_a, 1, 1, desc_h_1d, hbse_a_bl, 1, 1, desc_h_2d, cntxt_1d)
    call pX(gemr2d)(nmat, nmat, hbse_b, 1, 1, desc_h_1d, hbse_b_bl, 1, 1, desc_h_2d, cntxt_1d)
    ! FHJ: zero coupling block if the user requests. Useful to test the solver.
    if (xct%zero_coupling_block) hbse_b_bl = ZERO
    evals_t(:) = 0d0

    if (scal_h%icntxt>=0) then
      if (peinf%inode==0) write(6,*)
      call blacs_calc_broken_herm(hbse_a_bl, nmat, size(hbse_a_bl), desc_h_2d, 'A')
      call blacs_calc_broken_herm(hbse_b_bl, nmat, size(hbse_b_bl), desc_h_2d, 'B')
      if (peinf%inode==0) write(6,*)

      lwork = -1
      SAFE_ALLOCATE(work, (10))
#ifdef CPLX
      lrwork = -1
      SAFE_ALLOCATE(rwork, (10))
#endif
      liwork = -1
      SAFE_ALLOCATE(iwork, (10))
      ! FHJ: Call p*bseig in query mode
      CALL pX(BSEIG)(BSE_FULLBSE+BSE_DIRECT, nmat, &
        hbse_a_bl, 1, 1, desc_h_2d, &
        hbse_b_bl, 1, 1, desc_h_2d, evals_t, &
        evecs_r_bl, 1, 1, desc_e_2d, &
        work, lwork, &
#ifdef CPLX
        rwork, lrwork, &
#endif
        iwork, liwork, info)
      if (info/=0) then
        if(peinf%inode==0) then
          write(0,'(/,a,i0)') 'ERROR: Query mode for p*bseig failed with info=',info
          write(0,'(a,/)') 'Please, report this error to the BerkeleyGW developers!'
        endif
        call die("p*bseig internal error", only_root_writes=.true.)
      endif

      lwork = max(int(work(1)), 6*size(evecs_r_bl))
      SAFE_DEALLOCATE(work)
      SAFE_ALLOCATE(work, (lwork))
#ifdef CPLX
      lrwork = max(int(rwork(1)), 6*size(evecs_r_bl))
      SAFE_DEALLOCATE(rwork)
      SAFE_ALLOCATE(rwork, (lrwork))
#endif
      liwork = iwork(1)
      SAFE_DEALLOCATE(iwork)
      SAFE_ALLOCATE(iwork, (liwork))

      ! FHJ: Call p*bseig for realz
      if (peinf%inode==0) write(6,'(1x,a,i0)') 'Beginning ScaLAPACK diagonalization. Size: ', nmat
      CALL pX(BSEIG)(BSE_FULLBSE+BSE_DIRECT, nmat, &
        hbse_a_bl, 1, 1, desc_h_2d, &
        hbse_b_bl, 1, 1, desc_h_2d, evals_t, &
        evecs_r_bl, 1, 1, desc_e_2d, &
        work, lwork, &
#ifdef CPLX
        rwork, lrwork, &
#endif
        iwork, liwork, info)
      if (peinf%inode==0) write(6,*) 'Done ScaLAPACK diagonalization'

      SAFE_DEALLOCATE(work)
#ifdef CPLX
      SAFE_DEALLOCATE(rwork)
#endif
      SAFE_DEALLOCATE(iwork)

      if (info/=0) then
        if(peinf%inode==0) then
          write(0,'(/,a,i0)') 'ERROR: Diagonalization with p*bseig failed with info=',info
          write(0,'(a,/)') 'Please, report this error to the BerkeleyGW developers!'
        endif
        call die("p*bseig internal error", only_root_writes=.true.)
      endif
      ! FHJ: Copy eigenvalues/eigenvectors
      evals(1:neig) = evals_t(1:neig)
    endif !scal_h%icntxt>=0

    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    ! FHJ: PEs outside the grid don`t have the evals!
    call MPI_Bcast(evals, neig, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    ! FHJ: Copy eigenvectors from 2d block to 1d block-column format
    if (peinf%inode==0) write(6,*) 'Sharing eigenvectors'
    call pX(gemr2d)(2*nmat, nmat, evecs_r_bl, 1, 1, desc_e_2d, evecs_r, 1, 1, desc_e_1d, cntxt_1d)

    ! FHJ: Cleanup
    if (scal_h%icntxt>=0) call blacs_gridexit(scal_h%icntxt)
    call blacs_gridexit(cntxt_1d)
    SAFE_DEALLOCATE(hbse_a_bl)
    SAFE_DEALLOCATE(hbse_b_bl)
    SAFE_DEALLOCATE(evecs_r_bl)

    ! FHJ: We can calculate left eigenvectors (Y) from the block structure of
    ! the right eigenvectors (X):
    !        X = [ X_1, conj(X_2);     Y = [ X_1, -conj(X_2);
    !              X_2, conj(X_1) ],        -X_2,  conj(X_1) ].
    ! Since we only compute positive eigenvalues, i.e., the first block column
    ! of X, we only need to compute the first block column of Y:
    !        X_pos = [ X_1;     Y_pos = [ X_1;
    !                  X_2 ],            -X_2 ].
    ! The second block column of X and Y are associated with negative
    ! eigenvalues and are not computed explicitly here.
    if (size(evecs_l,2)>0) then
      evecs_l(1:nmat,:) = evecs_r(1:nmat,:)
      evecs_l(nmat+1:,:) = -evecs_r(nmat+1:,:)
    endif

    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    if (peinf%inode==0) write(6,*) 'Diagonalization done'

    POP_SUB(diagonalize.blacs_full_sseig)

  end subroutine blacs_full_sseig
#endif

end subroutine diagonalize


!> FHJ: Internal subroutine that actually prints ||M - M^H||_F / ||M + M^H||_F.
subroutine print_broken_herm(norm, block)
  real(DP), intent(in) :: norm
  character, intent(in) :: block

  logical, save :: warned=.false.

  PUSH_SUB(print_broken_herm)

  if (peinf%inode==0) then
    write(6,'(1x,a)',advance='no') 'Degree of broken Hermiticity of the '
    if (block=='A'.or.block=='B') then
      write(6,'(a)',advance='no') 'subblock '//block//' of the BSE Hamiltonian: '
    else
      write(6,'(a)',advance='no') 'BSE Hamiltonian: '
    endif
    if (norm>=1.and..not.warned) then
      write(6,'(f0.4,a)') norm, ' %'
      write(0,*)
      write(0,*) 'WARNING: non-Hermiticity of the BSE Hamiltonian is large!'
      write(0,*) 'Some possible reasons include:'
      write(0,*) '- Coarse grid that is too coarse'
      write(0,*) '- Large q0 shift (>0.001)'
      write(0,*) '- Small dielectric cutoff (<8 Ry)'
      write(0,*)
      warned = .true.
    else
      write(6,'(f6.4,a)') norm, ' %'
    endif
  endif

  POP_SUB(print_broken_herm)

end subroutine print_broken_herm


!> FHJ: Calculate the degree of broken Hermiticity: ||M - M^H||_F / ||M + M^H||_F.
subroutine calc_broken_herm(mat, nmat, block)
  integer, intent(in) :: nmat
  SCALAR, intent(in) :: mat(nmat,nmat)
  character, intent(in) :: block

  SCALAR :: mat2(nmat,nmat)
  real(DP) :: norm, norm2
  real(DP), external :: X(lange)

  PUSH_SUB(calc_broken_herm)

  ! Calc ||A - A^H||_F
  mat2(:,:) = transpose(mat(:,:))
  if (block=='B') then
    mat2(:,:) = mat(:,:) - mat2(:,:)
  else
    mat2(:,:) = mat(:,:) - MYCONJG(mat2(:,:))
  endif
  norm = X(lange)('F', nmat, nmat, mat2, nmat, mat2)
  ! Calc ||A + A^H||_F
  mat2(:,:) = transpose(mat(:,:))
  if (block=='B') then
    mat2(:,:) = mat(:,:) + mat2(:,:)
  else
    mat2(:,:) = mat(:,:) + MYCONJG(mat2(:,:))
  endif
  norm2 = X(lange)('F', nmat, nmat, mat2, nmat, mat2)

  norm = norm / norm2 * 1d2
  call print_broken_herm(norm, block)

  POP_SUB(calc_broken_herm)

end subroutine calc_broken_herm


#ifdef USESCALAPACK
!> FHJ: Calculate the degree of broken Hermiticity: ||M - M^H||_F / ||M + M^H||_F.
subroutine blacs_calc_broken_herm(mat, ldm, buf_sz, descm, block)
  integer, intent(in) :: buf_sz
  SCALAR, intent(in) :: mat(buf_sz)
  integer, intent(in) :: ldm
  integer, intent(in) :: descm(9)
  character, intent(in) :: block

  SCALAR :: mat2(buf_sz)
  real(DP) :: norm, norm2
  SCALAR, parameter :: MINUSONE=-ONE
  real(DP), external :: pX(lange)
  character :: op

  PUSH_SUB(blacs_calc_broken_herm)

  op = 'C'
  if (block=='B') op = 'T'

  ! Calc ||A - A^H||_F
  mat2(:) = mat(:)
  call pX(geadd)(op, ldm, ldm, MINUSONE, mat, 1, 1, descm, ONE, mat2, 1, 1, descm)
  norm = pX(lange)('F', ldm, ldm, mat2, 1, 1, descm, mat2)
  ! Calc ||A + A^H||_F
  mat2(:) = mat(:)
  call pX(geadd)(op, ldm, ldm, ONE, mat, 1, 1, descm, ONE, mat2, 1, 1, descm)
  norm2 = pX(lange)('F', ldm, ldm, mat2, 1, 1, descm, mat2)

  norm = norm / norm2 * 1d2
  call print_broken_herm(norm, block)

  POP_SUB(blacs_calc_broken_herm)

end subroutine blacs_calc_broken_herm


!> FHJ: Make diagonals of a distributed matrix real
subroutine blacs_make_diagonal_real(mat, ldm, buf_sz, descm)
  integer, intent(in) :: buf_sz
  SCALAR, intent(inout) :: mat(buf_sz)
  integer, intent(in) :: ldm
  integer, intent(in) :: descm(9)

  integer :: j_l, j_g, i_l, prow
  integer :: ctxt, M, N, MB, NB
  integer :: npr, npc, nprow, npcol, myprow, mypcol

  PUSH_SUB(blacs_make_diagonal_real)

  ctxt = descm(2)
  M = descm(3)
  N = descm(4)
  MB = descm(5)
  NB = descm(6)
  call blacs_gridinfo(ctxt, nprow, npcol, myprow, mypcol)
  npr = numroc(M, MB, myprow, 0, nprow)
  npc = numroc(N, NB, mypcol, 0, npcol)
  do j_l = 1, npc
    j_g = indxl2g(j_l, NB, mypcol, 0, npcol)
    prow = indxg2p(j_g, MB, myprow, 0, nprow) ! Since i_g = j_g
    if (prow==myprow) then
      i_l = indxg2l(j_g, MB, myprow, 0, nprow)
      if (i_l>=1 .and. i_l<=npr) then
        mat(i_l + ldm*(j_l-1)) = dble(mat(i_l + ldm*(j_l-1)))
      endif
    endif
  enddo

  POP_SUB(blacs_make_diagonal_real)

end subroutine blacs_make_diagonal_real


!> FHJ: Symmetrize a distributed matrix by copying the lower triangular
!! portion to the upper triangular. Applies complex conjugation if necessary.
subroutine blacs_symmetrize_l2u(mat, ldm, buf_sz, descm, do_conjg)
  integer, intent(in) :: buf_sz
  SCALAR, intent(inout) :: mat(buf_sz)
  integer, intent(in) :: ldm
  integer, intent(in) :: descm(9)
  logical :: do_conjg

  SCALAR :: mat_T(buf_sz)
#ifdef CPLX
  integer :: j_l, j_g, i_l, i_g, prow
  integer :: ctxt, M, N, MB, NB
  integer :: npr, npc, nprow, npcol, myprow, mypcol
#endif

  PUSH_SUB(blacs_symmetrize_l2u)

  ! Zero strictly upper triangular part of mat
  call pX(laset)('U', ldm-1, ldm-1, ZERO, ZERO, mat, 1, 2, descm)
  ! Transpose mat
#ifdef CPLX
  if (do_conjg) then
    call pztranc(ldm, ldm, ONE, mat, 1, 1, descm, ZERO, mat_T, 1, 1, descm)
  else
    call pztranu(ldm, ldm, ONE, mat, 1, 1, descm, ZERO, mat_T, 1, 1, descm)
  endif
#else
  call pdtran(ldm, ldm, ONE, mat, 1, 1, descm, ZERO, mat_T, 1, 1, descm)
#endif
  ! Zero lower triangular part of mat_T
  call pX(laset)('L', ldm, ldm, ZERO, ZERO, mat_T, 1, 1, descm)
  ! Add mat_T to mat
  call pX(geadd)('N', ldm, ldm, ONE, mat_T, 1, 1, descm, ONE, mat, 1, 1, descm)

#ifdef CPLX
  if (do_conjg) call blacs_make_diagonal_real(mat, ldm, buf_sz, descm)
#endif

  POP_SUB(blacs_symmetrize_l2u)

end subroutine blacs_symmetrize_l2u


!> FHJ: Symmetrize a distributed matrix by copying the upper triangular
!! portion to the lower triangular. Applies complex conjugation if necessary.
subroutine blacs_symmetrize_u2l(mat, ldm, buf_sz, descm, do_conjg)
  integer, intent(in) :: buf_sz
  SCALAR, intent(inout) :: mat(buf_sz)
  integer, intent(in) :: ldm
  integer, intent(in) :: descm(9)
  logical :: do_conjg

  SCALAR :: mat_T(buf_sz)

  PUSH_SUB(blacs_symmetrize_u2l)

  ! Zero strictly lower triangular part of mat
  call pX(laset)('L', ldm-1, ldm-1, ZERO, ZERO, mat, 2, 1, descm)
  ! Transpose mat
#ifdef CPLX
  if (do_conjg) then
    call pztranc(ldm, ldm, ONE, mat, 1, 1, descm, ZERO, mat_T, 1, 1, descm)
  else
    call pztranu(ldm, ldm, ONE, mat, 1, 1, descm, ZERO, mat_T, 1, 1, descm)
  endif
#else
  call pdtran(ldm, ldm, ONE, mat, 1, 1, descm, ZERO, mat_T, 1, 1, descm)
#endif
  ! Zero upper triangular part of mat_T
  call pX(laset)('U', ldm, ldm, ZERO, ZERO, mat_T, 1, 1, descm)
  ! Add mat_T to mat
  call pX(geadd)('N', ldm, ldm, ONE, mat_T, 1, 1, descm, ONE, mat, 1, 1, descm)

#ifdef CPLX
  if (do_conjg) call blacs_make_diagonal_real(mat, ldm, buf_sz, descm)
#endif

  POP_SUB(blacs_symmetrize_u2l)

end subroutine blacs_symmetrize_u2l


subroutine create_blacs_square_grid(scal, block_sz)
  type(scalapack), intent(out) :: scal
  integer, intent(in), optional :: block_sz

  PUSH_SUB(create_blacs_square_grid)

  call blacs_get(-1, 0, scal%icntxt)
  scal%nprow = int(sqrt(dble(peinf%npes)))
  scal%npcol = peinf%npes/scal%nprow
  call blacs_gridinit(scal%icntxt, 'R', scal%nprow, scal%npcol)
  call blacs_gridinfo(scal%icntxt, scal%nprow, scal%npcol, scal%myprow, scal%mypcol)
  scal%nbl = 64
  if (present(block_sz)) scal%nbl = block_sz
  if (peinf%inode==0) then
    write(6,'(/,1x,3(a,i0))') 'BLACS processor grid: ', &
      scal%nprow, ' x ', scal%npcol, '; BLOCKSIZE = ', scal%nbl
    write(6,'(1x,a,i0,/)') 'Number of idle processors: ', peinf%npes - scal%nprow*scal%npcol
  endif

  POP_SUB(create_blacs_square_grid)

end subroutine create_blacs_square_grid


subroutine create_blacs_rect_desc(MM, NN, scal, desc, lld)
  integer, intent(in) :: MM !< Number of rows in the global matrix
  integer, intent(in) :: NN !< Number of columns in the global matrix
  type(scalapack), intent(inout) :: scal
  integer, intent(out) :: desc(9)
  integer, intent(out) :: lld

  PUSH_SUB(create_blacs_rect_desc)

  if (scal%icntxt>=0) then
    scal%npr = numroc(MM, scal%nbl, scal%myprow, 0, scal%nprow)
    scal%npc = numroc(NN, scal%nbl, scal%mypcol, 0, scal%npcol)
  else
    scal%npr = 0
    scal%npc = 0
  endif
  lld = max(1, scal%npr)
  call descset(desc, MM, NN, scal%nbl, scal%nbl, 0, 0, scal%icntxt, lld)

  POP_SUB(create_blacs_rect_desc)

end subroutine create_blacs_rect_desc
#endif


end module diagonalize_m
