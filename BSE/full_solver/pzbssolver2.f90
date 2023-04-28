SUBROUTINE PZBSSOLVER2( N, H, IH, JH, DESCH, LAMBDA, X, IX, JX, &
                        DESCX, Y, IY, JY, DESCY, WORK, LWORK, &
                        INFO )
!
IMPLICIT NONE
!
!     .. Scalar Arguments ..
INTEGER            N, IH, JH, IX, JX, IY, JY, LWORK, INFO
!     ..
!     .. Array Arguments ..
INTEGER            DESCH( * ), DESCX( * ), DESCY( * )
COMPLEX*16         H( * ), X( * ), Y( * ), LAMBDA( * )
DOUBLE PRECISION   WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  PZBSSOLVER2() computes all eigenvalues and (both right and
!  left) eigenvectors of 2n-by-2n complex matrix
!
!     H = [       A,        B;
!          -conj(B), -conj(A) ],
!
!  where A is n-by-n Hermitian, and B is n-by-n symmetric.
!
!  On exit, the outputs X, Y, and lambda satisfy that
!
!     Y**H * X = I,
!     H * X = X * diag(lambda),
!     Y**H * H = diag(lambda) * Y**H.
!
!  Notes
!  =====
!
!  Each global data object is described by an associated description
!  vector.  This vector stores the information required to establish
!  the mapping between an object element and its corresponding process
!  and memory location.
!
!  Let A be a generic term for any 2D block cyclicly distributed array.
!  Such a global array has an associated description vector DESCA.
!  In the following comments, the character _ should be read as
!  "of the global array".
!
!  NOTATION        STORED IN      EXPLANATION
!  --------------- -------------- --------------------------------------
!  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
!                                 DTYPE_A = 1.
!  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
!                                 the BLACS process grid A is distribu-
!                                 ted over. The context itself is glo-
!                                 bal, but the handle (the integer
!                                 value) may vary.
!  M_A    (global) DESCA( M_ )    The number of rows in the global
!                                 array A.
!  N_A    (global) DESCA( N_ )    The number of columns in the global
!                                 array A.
!  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
!                                 the rows of the array.
!  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
!                                 the columns of the array.
!  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
!                                 row of the array A is distributed.
!  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
!                                 first column of the array A is
!                                 distributed.
!  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
!                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
!
!  Let K be the number of rows or columns of a distributed matrix,
!  and assume that its process grid has dimension p x q.
!  LOCr( K ) denotes the number of elements of K that a process
!  would receive if K were distributed over the p processes of its
!  process column.
!  Similarly, LOCc( K ) denotes the number of elements of K that a
!  process would receive if K were distributed over the q processes of
!  its process row.
!  The values of LOCr() and LOCc() may be determined via a call to the
!  ScaLAPACK tool function, NUMROC:
!          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
!          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
!  An upper bound for these quantities may be computed by:
!          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
!          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
!
!  Arguments
!  =========
!
!  N       (global input) INTEGER
!          The number of rows and columns of A and B.
!          The order of the distributed submatrices sub( M ), sub( X ),
!          and sub( Y ) are 2*N.
!          N >= 0.
!
!  H       (local input) COMPLEX*16 pointer into the local
!          memory to an array of dimension (LLD_H, LOCc(JH+2*N-1)).
!          The distributed matrix to be diagonalized.
!
!  IH      (global input) INTEGER
!          The row index in the global array H indicating the first
!          row of sub( H ).
!
!  JH      (global input) INTEGER
!          The column index in the global array H indicating the
!          first column of sub( H ).
!
!  DESCH   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix H.
!
!  LAMBDA  (global output) COMPLEX*16 array, dimension (2*N)
!
!  X       (local output) COMPLEX*16 array,
!          global dimension (2N, 2N),
!          local dimension ( LLD_X, LOCc(JX+2*N-1) )
!          On normal exit X contains the right eigenvectors of H.
!
!  IX      (global input) INTEGER
!          X`s global row index, which points to the beginning of the
!          submatrix which is to be operated on.
!
!  JX      (global input) INTEGER
!          X`s global column index, which points to the beginning of
!          the submatrix which is to be operated on.
!
!  DESCX   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix X.
!          DESCX( CTXT_ ) must equal DESCM( CTXT_ )
!
!  Y       (local output) COMPLEX*16 array,
!          global dimension (2N, 2N),
!          local dimension ( LLD_Y, LOCc(JY+2*N-1) )
!          On normal exit Y contains the left eigenvectors of H.
!
!  IY      (global input) INTEGER
!          Y`s global row index, which points to the beginning of the
!          submatrix which is to be operated on.
!
!  JY      (global input) INTEGER
!          Y`s global column index, which points to the beginning of
!          the submatrix which is to be operated on.
!
!  DESCY   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix Y.
!          DESCY( CTXT_ ) must equal DESCM( CTXT_ )
!
!  WORK    (local workspace/output) DOUBLE PRECISION array,
!          dimension (LWORK)
!          On output, WORK( 1 ) returns the minimal amount of workspace
!          needed to guarantee completion.
!          If the input parameters are incorrect, WORK( 1 ) may also be
!          incorrect.
!
!  LWORK   (local input) INTEGER
!          The length of the workspace array WORK.
!          If LWORK = -1, the LWORK is global input and a workspace
!          query is assumed; the routine only calculates the minimum
!          size for the WORK array. The required workspace is returned
!          as the first element of WORK and no error message is issued
!          by PXERBLA.
!
!  INFO    (global output) INTEGER
!          = 0:  successful exit
!          < 0:  If the i-th argument is an array and the j-entry had
!                an illegal value, then INFO = -(i*100+j), if the i-th
!                argument is a scalar and had an illegal value, then
!                INFO = -i.
!          > 0:  The eigensolver did not converge.
!
!  =====================================================================
!
!     .. Parameters ..
INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, &
                   LLD_, MB_, M_, NB_, N_, RSRC_
PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
COMPLEX*16         ZERO, ONE
PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
                     ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
LOGICAL            LQUERY
INTEGER            ICTXT, NPROCS, NPROW, NPCOL, MYROW, MYCOL, &
                   TWON, LWKOPT, LLWORK, LOCALMAT, LROWS, LCOLS, &
                   LLD, INDH, INDX, INDY, INDRWORK, INDWORK, I
DOUBLE PRECISION   DTMP
COMPLEX*16         ZTMP
!     ..
!     .. Local Arrays ..
INTEGER            DESCHL( DLEN_ ), DESCXL( DLEN_ ), &
                   DESCYL( DLEN_ )
!     ..
!     .. Intrinsic Functions ..
INTRINSIC          DBLE
!     ..
!     .. External Functions ..
EXTERNAL           ZDOTC
COMPLEX*16         ZDOTC
!     ..
!     .. External Subroutines ..
EXTERNAL           PZGEMR2D, PXERBLA, BLACS_GRIDINFO, CHK1MAT, &
                   DESCINIT, ZGSUM2D, ZGEEV, ZSCAL
!     ..
!     .. External Functions ..
EXTERNAL           NUMROC
INTEGER            NUMROC
!     ..
!     .. Executable Statements ..
!
INFO = 0
TWON = 2*N
ICTXT = DESCH( CTXT_ )
CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
NPROCS = NPROW*NPCOL
IF ( NPROW .EQ. -1 ) THEN
   INFO = -( 200+CTXT_ )
END IF
!
!     Test the input arguments.
!
IF ( INFO .EQ. 0 .AND. N .LT. 0 ) &
   INFO = -1
IF ( INFO .EQ. 0 ) &
   CALL CHK1MAT( TWON, 3, TWON, 3, IH, JH, DESCH, 5, INFO )
IF ( INFO .EQ. 0 ) &
   CALL CHK1MAT( TWON, 3, TWON, 3, IX, JX, DESCX, 10, INFO )
IF ( INFO .EQ. 0 .AND. DESCX( CTXT_ ) .NE. ICTXT ) &
   INFO = -( 1000+CTXT_ )
IF ( INFO .EQ. 0 ) &
   CALL CHK1MAT( TWON, 3, TWON, 3, IY, JY, DESCY, 14, INFO )
IF ( INFO .EQ. 0 .AND. DESCY( CTXT_ ) .NE. ICTXT ) &
   INFO = -( 1400+CTXT_ )
!
!     Compute required workspace.
!
IF ( INFO .EQ. 0 ) THEN
!
!        Set up local indices for the workspace.
!
   LQUERY = LWORK .EQ. -1
   LROWS = NUMROC( TWON, TWON, MYROW, 0, NPROW )
   LCOLS = NUMROC( TWON, TWON, MYCOL, 0, NPCOL )
   LLD = MAX( LROWS, 1 )
   LOCALMAT = LLD*LCOLS
   CALL DESCINIT( DESCHL, TWON, TWON, TWON, TWON, 0, 0, ICTXT, &
        LLD, INFO )
   DESCXL( 1:DLEN_ ) = DESCHL( 1:DLEN_ )
   DESCYL( 1:DLEN_ ) = DESCHL( 1:DLEN_ )
   INDH = 1
   INDX = INDH + 2*MAX( LOCALMAT, 1 )
   INDY = INDX + 2*MAX( LOCALMAT, 1 )
   INDRWORK = INDY + 2*MAX( LOCALMAT, 1 )
   INDWORK = INDRWORK + 2*TWON
   LLWORK = ( LWORK - INDWORK + 1 ) / 2
!
!        Estimate the workspace required by external subroutines.
!
   IF ( MYROW+MYCOL .EQ. 0 ) THEN
      CALL ZGEEV( 'V', 'V', TWON, WORK, LLD, LAMBDA, WORK, LLD, &
           WORK, LLD, ZTMP, -1, DTMP, INFO )
   ELSE
      ZTMP = ONE
   END IF
   LWKOPT = INDWORK - 1 + 2*INT( ZTMP )
!
   IF ( .NOT. LQUERY .AND. LWORK .LT. LWKOPT ) &
      INFO = -16
END IF
!
IF ( INFO .NE. 0 ) THEN
   CALL PXERBLA( ICTXT, 'PZBSSOLVER2', -INFO )
   RETURN
END IF
WORK( 1 ) = DBLE( LWKOPT )
IF ( LQUERY ) &
   RETURN
!
!     Quick return if possible.
!
IF ( N .EQ. 0 ) &
   RETURN
!
CALL PZGEMR2D( TWON, TWON, H, IH, JH, DESCH, WORK( INDH ), 1, 1, &
     DESCHL, ICTXT )
IF ( MYROW+MYCOL .EQ. 0 ) THEN
   CALL ZGEEV( 'V', 'V', TWON, WORK( INDH ), LLD, LAMBDA, &
        WORK( INDY ), LLD, WORK( INDX ), LLD, WORK( INDWORK ), &
        LLWORK, WORK( INDRWORK ), INFO )
   DO I = 1, TWON
      ZTMP = ZDOTC( TWON, WORK( INDX + 2*(I-1)*LLD ), 1, &
           WORK( INDY + 2*(I-1)*LLD ), 1 )
      IF ( ZTMP .NE. ZERO ) &
         CALL ZSCAL( TWON, ONE/ZTMP, WORK( INDY + 2*(I-1)*LLD ), &
              1 )
   END DO
ELSE
   LAMBDA( 1:TWON ) = ZERO
END IF
IF ( NPROW*NPCOL .GT. 1 ) &
   CALL ZGSUM2D( ICTXT, 'A', ' ', TWON, 1, LAMBDA, TWON, -1, -1 )
CALL PZGEMR2D( TWON, TWON, WORK( INDX ), 1, 1, DESCXL, X, IX, JX, &
     DESCX, ICTXT )
CALL PZGEMR2D( TWON, TWON, WORK( INDY ), 1, 1, DESCYL, Y, IY, JY, &
     DESCY, ICTXT )
!
WORK( 1 ) = DBLE( LWKOPT )
!
RETURN
!
!     End of PZBSSOLVER2().
!
END
