! This file is based on the work of Michel Olagnon.
! The original code for the MRGRNK subroutine is available at:
! http://fortran-2000.com/rank/mrgrnk.f90

! MRGRNK - Copyright (c) Michel Olagnon
! Copying and distribution of this file, with or without modification,
! are permitted in any medium without royalty provided the copyright
! notice and this notice are preserved.  This file is offered as-is,
! without any warranty.


! FHJ: WARNING - make sure you don`t change the following lines too much,
! otherwise they will be longer than 120 characters after the preprocessors kicks in.
! Note that, if there the extra "gvec" argument, we use a tolerance to figure
! out if the two items AA(ii) and AA(jj) are degenerate.
#ifdef HAS_GVEC
 #define is_greater(ii,jj) \
  (AA(ii)-AA(jj)>TOL.or.(AA(ii)-AA(jj)>-TOL.and.GK(ii)>GK(jj)))
#else
 #define is_greater(ii,jj) \
  (AA(ii)>AA(jj))
#endif

! FHJ: Note - we need to nest the "JOIN" macro otherwise the symbol LABEL
! doesn`t get expanded by the C preprocessor.
#define JOIN2(x,y) x ## _ ## y
#define JOIN(x,y) JOIN2(x,y)
! LABEL_GK is the sorting function without automatically generating the GK
! array out of the gvec array. This is useful for threaded sorting, when we
! first generate the GK array, but then call the sorting routine many times.
! LABEL_GK is the same as LABEL if there is no gvec input.
#ifdef HAS_GVEC
#define LABEL_GK JOIN(LABEL,GK)
#else
#define LABEL_GK LABEL
#endif
! This is the kernel of the sourting routine
#define LABEL_INSERTSORT JOIN(LABEL,insertsort)

! These functions are threaded wrappers for the sorting routines
#ifdef HAS_GVEC
#define LABEL_THREADED_GK JOIN(LABEL,threaded_GK)
#else
#define LABEL_THREADED_GK JOIN(LABEL,threaded)
#endif
#define LABEL_THREADED_GK_ACTUAL JOIN(LABEL_THREADED_GK,actual)
#define LABEL_THREADED JOIN(LABEL,threaded)
#define LABEL_THREADED_MERGE JOIN(LABEL,threaded_merge)

!> Sorts (actually, ranks) a small array AA using the insert sort algorithm.
subroutine LABEL_INSERTSORT(NVAL, AA, ord&
#ifdef HAS_GVEC
, GK&
#endif
)
  integer, intent(in) :: NVAL
  DTYPE, intent(in) :: AA(NVAL)
  integer, intent(inout) :: ord(NVAL)
#ifdef HAS_GVEC
  integer, intent(in) :: GK(NVAL)
  DTYPE, parameter :: TOL=TOL_ZERO
#endif

  integer :: ii, jj, tord

  PUSH_SUB(LABEL_INSERTSORT)

  do ii = 2, NVAL
    tord = ord(ii)
    jj = ii - 1
    do while (jj>0)
      if (.not.is_greater(ord(jj),tord)) exit
      ord(jj+1) = ord(jj)
      jj = jj - 1
    enddo
    ord(jj+1) = tord
  enddo

  POP_SUB(LABEL_INSERTSORT)

end subroutine LABEL_INSERTSORT

#ifdef HAS_GVEC
subroutine LABEL(NVAL, AA, ord, gvec)
  integer, intent(in) :: NVAL
  DTYPE, intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer, intent(in) :: gvec(3,NVAL) !< (3, N) G-vectors, used to break tie for equal data

  integer :: GK(NVAL)

  PUSH_SUB(LABEL)

  call get_GK_array_from_gvecs(NVAL, gvec, GK)
  call LABEL_GK(NVAL, AA, ord, GK)

  POP_SUB(LABEL)

end subroutine LABEL
#endif

!> Sorts (actually, ranks) a real/integer array AA.
!! The rank is written to the output array ord.
!! This subroutine is based on the routine MRGRNK by Michel Olagnon, which
!! uses the merge sort algorithm.
subroutine LABEL_GK(NVAL, AA, ord&
#ifdef HAS_GVEC
, GK&
#endif
)
  integer, intent(in) :: NVAL
  DTYPE, intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
#ifdef HAS_GVEC
  integer, intent(in) :: GK(NVAL)
#endif

#ifdef HAS_GVEC
  DTYPE, parameter :: TOL=TOL_ZERO
#endif
  integer :: JT(NVAL)
  integer :: LMTNA, LMTNC, IRNG1, IRNG2
  integer :: IIND, ID, IWRK, IWRKF, JINDA, IA, IB

  PUSH_SUB(LABEL_GK)
!
!  Fill-in the index array, creating ordered couples
!
  Do IIND = 2, NVAL, 2
    If (&
is_greater(IIND,IIND-1)&
    ) Then
      ord (IIND-1) = IIND - 1
      ord (IIND) = IIND
    Else
      ord (IIND-1) = IIND
      ord (IIND) = IIND - 1
    End If
  End Do
  If (Modulo(NVAL, 2) /= 0) Then
    ord (NVAL) = NVAL
  End If

  ! FHJ - shortcut if the array is small enough
  if (NVAL<16) then
    call LABEL_INSERTSORT(NVAL, AA, ord&
#ifdef HAS_GVEC
, GK&
#endif
)
    POP_SUB(LABEL_GK)
    return
  endif

!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into    C  -  C  - ...
!
  LMTNA = 2
  LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
  Do
    If (NVAL <= 2) Exit
!
!  Loop on merges of A and B into C
!
    Do ID = 0, NVAL - 1, 4
      If ((ID+4) > NVAL) Then
        If ((ID+2) >= NVAL) Exit
!
!  1 2 3
!
        If (&
is_greater(ord(ID+3),ord(ID+2))&
        ) Exit
!
!  1 3 2
!
        If (&
is_greater(ord(ID+3),ord(ID+1))&
        ) Then
          IRNG2 = ord (ID+2)
          ord (ID+2) = ord (ID+3)
          ord (ID+3) = IRNG2
!
!  3 1 2
!
        Else
          IRNG1 = ord (ID+1)
          ord (ID+1) = ord (ID+3)
          ord (ID+3) = ord (ID+2)
          ord (ID+2) = IRNG1
        End If
        Exit
      End If
!
!  1 2 3 4
!
      If (&
is_greater(ord(ID+3),ord(ID+2))&
      ) Cycle
!
!  1 3 x x
!
      If (&
is_greater(ord(ID+3),ord(ID+1))&
      ) Then
        IRNG2 = ord (ID+2)
        ord (ID+2) = ord (ID+3)
        If (&
is_greater(ord(ID+4),IRNG2)&
        ) Then
!  1 3 2 4
          ord (ID+3) = IRNG2
        Else
!  1 3 4 2
          ord (ID+3) = ord (ID+4)
          ord (ID+4) = IRNG2
        End If
!
!  3 x x x
!
      Else
        IRNG1 = ord (ID+1)
        IRNG2 = ord (ID+2)
        ord (ID+1) = ord (ID+3)
        If (&
is_greater(ord(ID+4),IRNG1)&
        ) Then
          ord (ID+2) = IRNG1
          If (&
is_greater(ord(ID+4),IRNG2)&
          ) Then
!  3 1 2 4
            ord (ID+3) = IRNG2
          Else
!  3 1 4 2
            ord (ID+3) = ord (ID+4)
            ord (ID+4) = IRNG2
          End If
        Else
!  3 4 1 2
          ord (ID+2) = ord (ID+4)
          ord (ID+3) = IRNG1
          ord (ID+4) = IRNG2
        End If
      End If
    End Do
!
!  The Cs become As and Bs
!
    LMTNA = 4
    Exit
  End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
  Do
    If (LMTNA >= NVAL) Exit
    IWRKF = 0
    LMTNC = 2 * LMTNC
!
!  Loop on merges of A and B into C
!
    Do
      IWRK = IWRKF
      ID = IWRKF + 1
      JINDA = IWRKF + LMTNA
      IWRKF = IWRKF + LMTNC
      If (IWRKF >= NVAL) Then
        If (JINDA >= NVAL) Exit
        IWRKF = NVAL
      End If
      IA = 1
      IB = JINDA + 1
!
!  Shortcut for the case when the max of A is smaller
!  than the min of B. This line may be activated when the
!  initial set is already close to sorted.
!
      IF (&
is_greater(ord(IB),ord(JINDA))&
      ) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
      JT (1:LMTNA) = ord (ID:JINDA)
!
      Do
        IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
        If (&
is_greater(JT(IA),ord(IB))&
        ) Then
          ord (IWRK) = ord (IB)
          IB = IB + 1
          If (IB > IWRKF) Then
!  Only A still with unprocessed values
            ord (IWRK+1:IWRKF) = JT (IA:LMTNA)
            Exit
          End If
        Else
          ord (IWRK) = JT (IA)
          IA = IA + 1
          If (IA > LMTNA) Exit! Only B still with unprocessed values
        End If
!
      End Do
    End Do
!
!  The Cs become As and Bs
!
    LMTNA = 2 * LMTNA
  End Do
!
  POP_SUB(LABEL_GK)
!
End Subroutine LABEL_GK


#ifdef HAS_GVEC
subroutine LABEL_THREADED(NVAL, AA, ord, gvec)
  integer, intent(in) :: NVAL
  DTYPE, intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer, intent(in) :: gvec(3,NVAL) !< (3, N) G-vectors, used to break tie for equal data

  integer :: GK(NVAL)

  PUSH_SUB(LABEL_THREADED)

  call get_GK_array_from_gvecs(NVAL, gvec, GK)
  call LABEL_THREADED_GK(NVAL, AA, ord, GK)

  POP_SUB(LABEL_THREADED)

end subroutine LABEL_THREADED
#endif


!> Threaded sort (actually, ranking) of a real/integer array AA.
!! Public code adapted from:
!! https://github.com/cphyc/Fortran-parallel-sort/blob/master/mod_sort.f90
subroutine LABEL_THREADED_GK(NVAL, AA, ord&
#ifdef HAS_GVEC
, GK&
#endif
)
  integer, intent(in) :: NVAL
  DTYPE, intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
#ifdef HAS_GVEC
  integer, intent(in) :: GK(NVAL)
#endif

  integer :: nthreads
#ifdef OMP
  integer, allocatable :: ord_loc(:)
  integer :: from, middle, to, thread, chunk, chunk2, ii, lim
  logical :: old_is_dynamic
  integer :: old_num_threads
#endif

  PUSH_SUB(LABEL_THREADED_GK)

  ! FHJ: TODO - nthreads should be different that peinf%nthreads_sort if we
  ! have nested OMPs. However, starting a $OMP PARALLEL region to figure out
  ! th number of threads is just too much overhead.
  nthreads = peinf%nthreads_sort
  if (nthreads==1.or.NVAL<nthreads*100) then
    call LABEL_GK(NVAL, AA, ord&
#ifdef HAS_GVEC
, GK&
#endif
)
    POP_SUB(LABEL_THREADED_GK)
    return
  endif

#ifdef OMP

  ! Size of each chunk
  chunk = NVAL / nthreads

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(from, middle, to, thread, lim, ii, chunk2, ord_loc) &
  !$OMP NUM_THREADS(nthreads)

  ! Sort each chunk
  !$OMP DO SCHEDULE(STATIC,1)
  do thread = 0, nthreads
    from = thread * chunk + 1
    to = min((thread+1) * chunk, NVAL)
    if (from<=to) then
      ! Initialize order
      SAFE_ALLOCATE(ord_loc, (to-from+1))
      !ord(from:to) = ord(from:to) - from + 1
      call LABEL_GK(to-from+1, AA(from:to), ord_loc&
#ifdef HAS_GVEC
, GK(from:to)&
#endif
)
      ord(from:to) = ord_loc + from - 1
      SAFE_DEALLOCATE(ord_loc)
    endif
  enddo
  !$OMP END DO

  ! Merge pieces
  ii = 1
  chunk2 = chunk
  do while (chunk2 < NVAL)
    lim = DIVUP(NVAL, 2*chunk2)
    !$OMP DO SCHEDULE(STATIC,1)
    do thread = 0, lim
      from = thread*2 * chunk2 + 1
      middle = min((thread*2+1) * chunk2, NVAL)
      to = min((thread*2+2) * chunk2, NVAL)
      if (from<to) then
        call LABEL_THREADED_MERGE(NVAL, AA, ord, &
#ifdef HAS_GVEC
GK, &
#endif
        from, middle, to)
      endif
    enddo
    !$OMP END DO
    chunk2 = chunk2 * 2
    ii = ii + 1
  enddo
  !$OMP END PARALLEL

#endif

  POP_SUB(LABEL_THREADED_GK)

end subroutine LABEL_THREADED_GK


! Merge two parts of A, ordered by order from left to right around middle.
subroutine LABEL_THREADED_MERGE(NVAL, AA, ord, &
#ifdef HAS_GVEC
GK, &
#endif
left, middle, right)
  integer, intent(in) :: NVAL
  DTYPE, intent(in) :: AA(NVAL)
  integer, intent(inout) :: ord(NVAL)
#ifdef HAS_GVEC
  integer, intent(in) :: GK(NVAL)
#endif
  integer, intent(in) :: left, middle, right

  integer :: leftA, rightA, leftB, rightB
  integer :: iA, iB, i
  integer :: lenA, lenB
  integer :: idxA, idxB
  integer :: orderA(left:middle)
  integer :: orderB(middle+1:right)
#ifdef HAS_GVEC
  DTYPE, parameter :: TOL=TOL_ZERO
#endif

  PUSH_SUB(LABEL_THREADED_MERGE)

  ! copy order
  orderA = ord(left:middle)
  orderB = ord(middle+1:right)

  ! more explicit variables
  leftA = left
  rightA = middle
  leftB = middle+1
  rightB = right

  ! initialize iA, iB to their leftmost position
  iA = leftA
  iB = leftB
  i = leftA

  do while ((iA <= rightA) .and. (iB <= rightB))
    idxA = orderA(iA)
    idxB = orderB(iB)
    if (is_greater(idxA,idxB)) then
       ord(i) = idxB
       iB = iB + 1
    else
       ord(i) = idxA
       iA = iA + 1
    end if
    i = i + 1
  end do

  ! either A or B still have elements, append them to the new order
  do while (iA <= rightA)
    ord(i) = orderA(iA)
    iA = iA + 1
    i  = i + 1
  end do
  do while (iB <= rightB)
    ord(i) = orderB(iB)
    iB = iB + 1
    i = i + 1
  end do

  POP_SUB(LABEL_THREADED_MERGE)

end subroutine LABEL_THREADED_MERGE

#undef LABEL_THREADED_MERGE
#undef LABEL_THREADED
#undef LABEL_THREADED_GK
#undef LABEL_INSERTSORT
#undef LABEL_GK
#undef JOIN
#undef JOIN2
#undef is_greater
