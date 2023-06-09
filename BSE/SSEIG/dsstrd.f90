SUBROUTINE DSSTRD( UPLO, N, A, LDA, E, TAU, WORK, LWORK, INFO )
!
IMPLICIT NONE
!
!     .. Scalar Arguments ..
CHARACTER          UPLO
INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSSTRD reduces a real skew-symmetric matrix A to skew-symmetric
!  tridiagonal form T by an orthogonal similarity transformation:
!     Q**T * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          skew-symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the skew-symmetric matrix A.  If UPLO = 'U', the
!          leading n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the first superdiagonal of A is
!          overwritten by the corresponding elements of the tridiagonal
!          matrix T, and the elements above the first superdiagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors; if UPLO = 'L', the first
!          subdiagonal of A is overwritten by the corresponding elements
!          of the tridiagonal matrix T, and the elements below the first
!          subdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The superdiagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U',
!          E(i) = -A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) DOUBLE PRECISION array,
!          dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= 1.
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v**T,
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v**T,
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  0   e   v2  v3  v4 )              (  0                  )
!    (      0   e   v3  v4 )              (  -e  0              )
!    (          0   e   v4 )              (  v1  -e  0          )
!    (              0   e  )              (  v1  v2  -e  0      )
!    (                  0  )              (  v1  v2  v3  -e  0  )
!
!  where e denotes super-diagonal elements of T, and vi denotes an
!  element of the vector defining H(i).
!
!  =====================================================================
!
!  Level 3 LAPACK-like routine.
!
!  DSSTRD is modified from DSYTRD in LAPACK version 3.5.0.
!
!  Written by Meiyue Shao, Lawrence Berkeley National Laboratory.
!  Last change: October 2014
!
!  =====================================================================
!
!     .. Parameters ..
DOUBLE PRECISION   ONE
PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
LOGICAL            LQUERY, UPPER
INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB, &
                   NBMIN, NX
!     ..
!     .. External Subroutines ..
EXTERNAL           DLASTD, DSSR2K, DSSTD2, XERBLA
!     ..
!     .. External Functions ..
INTEGER            ILAENV
LOGICAL            LSAME
EXTERNAL           ILAENV, LSAME
!     ..
!     .. Intrinsic Functions ..
INTRINSIC          MAX, DBLE
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
INFO = 0
UPPER = LSAME( UPLO, 'U' )
LQUERY = ( LWORK .EQ. -1 )
IF ( .NOT. UPPER .AND. .NOT. LSAME( UPLO, 'L' ) ) THEN
   INFO = -1
ELSE IF ( N .LT. 0 ) THEN
   INFO = -2
ELSE IF ( LDA .LT. MAX( 1, N ) ) THEN
   INFO = -4
ELSE IF ( LWORK .LT. 1 .AND. .NOT. LQUERY ) THEN
   INFO = -8
END IF
!
IF ( INFO .EQ. 0 ) THEN
!
!        Determine the block size.
!
   NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
   LWKOPT = N*NB
   WORK( 1 ) = DBLE( LWKOPT )
END IF
!
IF ( INFO .NE. 0 ) THEN
   CALL XERBLA( 'DSSTRD', -INFO )
   RETURN
ELSE IF ( LQUERY ) THEN
   RETURN
END IF
!
!     Quick return if possible.
!
IF ( N .LE. 1 ) THEN
   WORK( 1 ) = ONE
   RETURN
END IF
!
NX = N
IWS = 1
IF ( NB .GT. 1 .AND. NB .LT. N ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).
!
   NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
   IF ( NX .LT. N ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
      LDWORK = N
      IWS = LDWORK*NB
      IF ( LWORK .LT. IWS ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code by setting NX = N.
!
         NB = MAX( LWORK / LDWORK, 1 )
         NBMIN = ILAENV( 2, 'DSYTRD', UPLO, N, -1, -1, -1 )
         IF ( NB .LT. NBMIN ) &
            NX = N
      END IF
   ELSE
      NX = N
   END IF
ELSE
   NB = 1
END IF
!
IF ( UPPER ) THEN
!
!        Reduce the upper triangle of A.
!        Columns 1:kk are handled by the unblocked method.
!
   KK = N - ( ( N-NX+NB-1 ) / NB )*NB
   DO I = N-NB+1, KK+1, -NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix.
!
      CALL DLASTD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK, &
           LDWORK )
!
!           Update the unreduced submatrix A(1:i-1,1:i-1), using an
!           update of the form:  A := A + V*W**T - W*V**T.
!
      CALL DSSR2K( UPLO, 'N', I-1, NB, ONE, A( 1, I ), LDA, WORK, &
           LDWORK, ONE, A, LDA )
!
!           Copy superdiagonal elements back into A.
!
      DO J = I, I+NB-1
         A( J-1, J ) = E( J-1 )
      END DO
   END DO
!
!        Use unblocked code to reduce the last or only block
!
   CALL DSSTD2( UPLO, KK, A, LDA, E, TAU, IINFO )
ELSE
!
!        Reduce the lower triangle of A.
!
   DO I = 1, N-NX, NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix.
!
      CALL DLASTD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ), &
                   TAU( I ), WORK, LDWORK )
!
!           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
!           an update of the form:  A := A + V*W**T - W*V**T.
!
      CALL DSSR2K( UPLO, 'N', N-I-NB+1, NB, ONE, A( I+NB, I ), &
           LDA, WORK( NB+1 ), LDWORK, ONE, A( I+NB, I+NB ), LDA )
!
!           Copy subdiagonal elements back into A.
!
      DO J = I, I+NB-1
         A( J+1, J ) = -E( J )
      END DO
   END DO
!
!        Use unblocked code to reduce the last or only block.
!
   CALL DSSTD2( UPLO, N-I+1, A( I, I ), LDA, E( I ), TAU( I ), &
        IINFO )
END IF
!
WORK( 1 ) = DBLE( LWKOPT )
RETURN
!
!     End of DSSTRD.
!
END
