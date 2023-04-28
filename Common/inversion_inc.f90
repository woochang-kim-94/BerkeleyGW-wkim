!===========================================================================
!
! Included from inversion.F90
!
!============================================================================

!---------------------- Use scaLAPACK For Inversion -----------------------------------


#ifdef USESCALAPACK

subroutine X(invert_with_scalapack)(nmtx, scal, matrix,integ_rpa_val)
  integer, intent(in) :: nmtx
  type (scalapack), intent(in) :: scal
  SCALAR, intent(inout) :: matrix(scal%npr,scal%npc)
  real(DP), optional, intent(out) :: integ_rpa_val

  integer :: info, desca(9), lwork, liwork, ipiv(scal%npr+scal%nbl)
  integer, allocatable :: iwork(:)
  SCALAR, allocatable :: work(:)

  SCALAR, allocatable :: diag_pol(:), diag_LU(:)
  integer :: i, j, irow, icol, icurr, irowm, icolm
  logical :: integrate_rpa

  PUSH_SUB(X(invert_with_scalapack))

  integrate_rpa = .false.
  if(present(integ_rpa_val)) integrate_rpa = .true.

  if(integrate_rpa) then
    allocate(diag_pol(nmtx))
    diag_pol = ZERO

    icurr=0
    do i=1, nmtx
      irow=mod(int(((i-1)/scal%nbl)+TOL_SMALL),scal%nprow)
      if(irow.ne.scal%myprow) cycle
      do j=1, nmtx
        icol=mod(int(((j-1)/scal%nbl)+TOL_SMALL),scal%npcol)
        if(icol .eq. scal%mypcol) then
          icurr=icurr+1
          irowm=int((icurr-1)/scal%npc+TOL_SMALL)+1
          icolm=mod((icurr-1),scal%npc)+1

          if(i == j) diag_pol(i) = matrix(irowm,icolm) - ONE

        endif
      end do
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE,diag_pol(1),nmtx,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)

  endif

  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*)
    write(6,*) 'We are doing inversion with ScaLAPACK.'
    write(6,*)
    write(6,*) 'Peinf: ', peinf%inode,' Scalapack Grid:'
    write(6,*) scal%nprow,scal%npcol,scal%myprow,scal%mypcol
    write(6,*) nmtx,scal%nbl,scal%npr,scal%npc
    write(6,*)
  endif

  call descinit(desca, nmtx, nmtx, scal%nbl, scal%nbl, 0, 0, &
    scal%icntxt, max(1,scal%npr), info)
  if(info < 0) then
    write(0,'(a,i3,a)') 'Argument number ', -info, ' had an illegal value on entry.'
    call die("descinit error for descaA in inversion")
  else if(info > 0) then
    write(0,*) 'info = ', info
    call die("descinit error for descaA in inversion")
  endif

  ! FHJ: LU factorization of the matrix
  call pX(getrf)(nmtx, nmtx, matrix, 1, 1, desca, ipiv, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in p?getrf'
    call die('p?getrf failed')
  endif

  if(integrate_rpa) then
    allocate(diag_LU(nmtx))
    diag_LU = ZERO

    icurr=0
    do i=1, nmtx
      irow=mod(int(((i-1)/scal%nbl)+TOL_SMALL),scal%nprow)
      if(irow.ne.scal%myprow) cycle
      do j=1, nmtx
        icol=mod(int(((j-1)/scal%nbl)+TOL_SMALL),scal%npcol)
        if(icol .eq. scal%mypcol) then
          icurr=icurr+1
          irowm=int((icurr-1)/scal%npc+TOL_SMALL)+1
          icolm=mod((icurr-1),scal%npc)+1

          IF(i == j) diag_LU(i) = matrix(irowm,icolm) 

        endif
      end do
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE,diag_LU(1),nmtx,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)

    ! calculate integrand RPA value
    ! calculate determinant
    integ_rpa_val = 1.0D+00
    do i = 1, nmtx
      integ_rpa_val = diag_LU(i)*integ_rpa_val
    enddo
    ! take the log
    integ_rpa_val = log(integ_rpa_val)
    ! remove trace of polarizability
    do i = 1, nmtx
      integ_rpa_val = integ_rpa_val - diag_pol(i)
    enddo

  endif

  ! FHJ: tringular inversion of LU decomposition
  SAFE_ALLOCATE(work, (10))
  SAFE_ALLOCATE(iwork, (10))
  call pX(getri)(nmtx, matrix, 1, 1, desca, ipiv, work, -1, iwork, -1, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in p?getri'
    call die('p?getri failed for query mode')
  endif

  lwork = max(1,int(work(1)))
  liwork = max(1,iwork(1))
  SAFE_DEALLOCATE(work)
  SAFE_DEALLOCATE(iwork)
  SAFE_ALLOCATE(work, (lwork))
  SAFE_ALLOCATE(iwork, (liwork))

  call pX(getri)(nmtx, matrix, 1, 1, desca, ipiv, work, lwork, iwork, liwork, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in p?getri'
    call die('p?getri failed')
  endif

  SAFE_DEALLOCATE(iwork)
  SAFE_DEALLOCATE(work)
  
  if(integrate_rpa) then
    deallocate(diag_pol)
    deallocate(diag_LU)
  endif

  POP_SUB(X(invert_with_scalapack))
  return
end subroutine X(invert_with_scalapack)

#endif

!------------------------------------------------------------

subroutine X(invert_serial)(nmtx, matrix)
  integer, intent(in) :: nmtx
  SCALAR, intent(inout) :: matrix(nmtx,nmtx)

  integer :: ii, info, lwork, ipiv(nmtx)
  SCALAR, allocatable :: work(:)

  PUSH_SUB(X(invert_serial))

  ! FHJ: LU factorization of the matrix
  call X(getrf)(nmtx, nmtx, matrix, nmtx, ipiv, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getrf'
    call die('?getrf failed')
  endif

  ! FHJ: tringular inversion of LU decomposition
  SAFE_ALLOCATE(work, (10))
  call X(getri)(nmtx, matrix, nmtx, ipiv, work, -1, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
    call die('?getri failed for query mode')
  endif

  lwork = max(1,int(work(1)))
  SAFE_DEALLOCATE(work)
  SAFE_ALLOCATE(work, (lwork))

  call X(getri)(nmtx, matrix, nmtx, ipiv, work, lwork, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
    call die('?getri failed')
  endif

  SAFE_DEALLOCATE(work)
  
  POP_SUB(X(invert_serial))
  return
end subroutine X(invert_serial)
