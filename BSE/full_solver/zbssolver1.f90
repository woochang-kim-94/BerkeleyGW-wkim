SUBROUTINE ZBSSOLVER1( N, M, LDM, LAMBDA, X, LDX, WORK, LWORK, &
                       IWORK, LIWORK, INFO )
!
IMPLICIT NONE
!
!     .. Scalar Arguments ..
INTEGER            N, LDM, LDX, LWORK, LIWORK, INFO
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION   M( LDM, * ), LAMBDA( * ), WORK( * )
COMPLEX*16         X( LDX, * )
INTEGER            IWORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZBSSOLVER1() computes all eigenvalues and (both right and
!  left) eigenvectors of 2n-by-2n complex matrix
!
!     H = [       A,        B;
!          -conj(B), -conj(A) ],
!
!  where A is n-by-n Hermitian, and B is n-by-n symmetric.
!
!  On entry, the information of H is provided in the lower triangular
!  part of
!
!     M = [ real(A+B), imag(A-B);
!          -imag(A+B), real(A-B) ].
!
!  The matrix M is required to be positive definite.
!
!  The structure of H leads to the following properties.
!
!     1) H is diagonalizable with n pairs of real eigenvalues
!        (lambda_i, -lambda_i).
!
!     2) The eigenvectors of H has the block structure
!
!           X = [ X_1, conj(X_2);     Y = [ X_1, -conj(X_2);
!                 X_2, conj(X_1) ],        -X_2,  conj(X_1) ],
!
!        and satisfy that
!
!           X_1**T * X_2 = X_2**T * X_1,
!           X_1**H * X_1 - X_2**H * X_2 = I,
!           Y**H * X = I,
!           H * X = X * diag(lambda, -lambda),
!           Y**H * H = diag(lambda, -lambda) * Y**H.
!
!  On exit, only the positive eigenvalues and the corresponding right
!  eigenvectors are returned.  The eigenvalues are sorted in ascending
!  order.  The eigenvectors are normalized (i.e., X = [ X_1; X_2 ] with
!  X_1**H * X_1 - X_2**H * X_2 = I).
!
!  M is destroyed on exit.
!
!  Arguments
!  =========
!
!  N       (input) INGEGER
!          2*N is the number of rows and columns of M.
!
!  M       (input/output) DOUBLE PRECISION array, dimension (LDM,2*N)
!          On entry, the symmetric positive definite matrix M. Only the
!          lower triangular part of M is used to define the elements of
!          the symmetric matrix.
!          On exit, all entries of M are destroyed.
!
!  LDM     (input) INGEGER
!          The leading dimension of the array M. LDM >= max(1,2*N).
!
!  LAMBDA  (output) DOUBLE PRECISION array, dimension (N)
!          On normal exit LAMBDA contains the positive eigenvalues of H
!          in ascending order.
!
!  X       (output) COMPLEX*16 array, dimension (LDX, N)
!          On normal exit X contains the normalized right eigenvectors
!          of H corresponding to the positive eigenvalues.
!
!  LDX     (input) INGEGER
!          The leading dimension of the array X. LDX >= max(1,2*N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array,
!          dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK( 1 ) returns the optimal LWORK.
!
!  LWORK   (input) INGEGER
!          The dimension of the array WORK.
!          If LWORK = -1 or LIWORK = -1, then a workspace query is
!          assumed; the routine only calculates the optimal sizes of
!          WORK and IWORK, returns these values as the first entries of
!          WORK and IWORK, respectively, and no error message related to
!          LWORK or LIWORK is issued by XERBLA.
!
!  IWORK   (workspace/output) INTEGER array,
!          dimension (MAX(1,LIWORK))
!          On exit, if INFO = 0, IWORK( 1 ) returns the optimal LIWORK.
!
!  LWORK   (input) INGEGER
!          The dimension of the array WORK.
!          If LWORK = -1 or LIWORK = -1, then a workspace query is
!          assumed; the routine only calculates the optimal sizes of
!          WORK and IWORK, returns these values as the first entries of
!          WORK and IWORK, respectively, and no error message related to
!          LWORK or LIWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  the eigensolver did not converge.
!
!  =====================================================================
!
!     .. Parameters ..
DOUBLE PRECISION   ZERO, ONE, TWO, HALF
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0, &
                     HALF = 0.5D+0 )
!     ..
!     .. Local Scalars ..
LOGICAL            LQUERY
INTEGER            I, J, ITMP, LWKOPT, LIWKOPT, LLWORK, LLIWORK, &
                   TWON, DIMV, NSPLIT, INDD, INDE, INDTAU, INDU, &
                   INDV, INDW, INDWORK, INDIFAIL, INDIBLOCK, &
                   INDISPLIT, INDIWORK
DOUBLE PRECISION   ABSTOL, DTMP
!     ..
!     .. Intrinsic Functions ..
INTRINSIC          DBLE, DSQRT, INT, MAX
!     ..
!     .. External Subroutines ..
EXTERNAL           ZBSAUX1, DSORTEIG, DSCAL, DTRMM, DLACPY, &
                   DLASET, DORMTR, DPOTRF, DSSTRD, DSTEBZ, DSTEIN, &
                   XERBLA, ZLASCL
!     ..
!     .. External Functions ..
EXTERNAL           DLAMCH
DOUBLE PRECISION   DLAMCH
!     ..
!     .. Executable Statements ..
!
INFO = 0
LQUERY = LWORK .EQ. -1 .OR. LIWORK .EQ. -1
TWON = 2*N
!
!     Test the input arguments.
!
IF ( N .LT. 0 ) THEN
   INFO = -1
ELSE IF ( LDM .LT. MAX( 1, TWON ) ) THEN
   INFO = -3
ELSE IF ( LDX .LT. MAX( 1, TWON ) ) THEN
   INFO = -6
END IF
!
!     Compute required workspace.
!
IF ( INFO .EQ. 0 ) THEN
!
!        Set the indices of temporary matrices in the workspace.
!
   ABSTOL = TWO*DLAMCH( 'S' )
   INDE = 1
   INDW = INDE + TWON
   INDU = INDW + TWON*TWON
   INDV = INDU + TWON*TWON
   INDWORK = INDV + TWON*TWON
   INDD = INDW
   INDTAU = INDV
   LLWORK = LWORK - INDWORK + 1
   INDIBLOCK = 1
   INDISPLIT = INDIBLOCK + TWON
   INDIFAIL = INDISPLIT + TWON
   INDIWORK = INDIFAIL + N
   LLIWORK = LIWORK - INDIWORK + 1
!
!        Estimate the workspace required by external subroutines.
!
   CALL DSSTRD( 'U', TWON, DTMP, TWON, DTMP, DTMP, WORK, -1, &
        ITMP )
   LWKOPT = INT( WORK( 1 ) )
   CALL DORMTR( 'R', 'U', 'N', TWON, TWON, DTMP, TWON, DTMP, DTMP, &
        TWON, WORK, -1, ITMP )
   LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
   LWKOPT = MAX( LWKOPT, MAX( 5, TWON )*TWON )
   LIWKOPT = 3*TWON
!
   LWKOPT = INDWORK - 1 + LWKOPT
   LIWKOPT = INDIWORK - 1 + LIWKOPT
   IF ( .NOT. LQUERY ) THEN
      IF ( LWORK .LT. LWKOPT ) THEN
         INFO = -8
      ELSE IF ( LIWORK .LT. LIWKOPT ) THEN
         INFO = -10
      END IF
   END IF
END IF
!
IF ( INFO .NE. 0 ) THEN
   CALL XERBLA( 'ZBSSOLVER1', -INFO )
   RETURN
END IF
!
WORK( 1 ) = DBLE( LWKOPT )
IWORK( 1 ) = LIWKOPT
IF ( LQUERY ) RETURN
!
!     Quick return if possible.
!
IF ( N .EQ. 0 ) RETURN
!
!     Compute the Cholesky factorization M = L * L**T.
!     In case of failure, a general solver is needed.
!
CALL DPOTRF( 'L', TWON, M, LDM, INFO )
IF ( INFO .NE. 0 ) THEN
   INFO = -2
   RETURN
END IF
!
!     Explicitly formulate W = L**T * J * L,
!     where
!
!        J = [   0, I_n;
!             -I_n,   0 ]
!
!     is the standard symplectic matrix.
!     Only the upper triangular part of W is filled by the correct
!     values.
!
CALL DLASET( 'U', TWON, TWON, ZERO, ZERO, WORK( INDW ), TWON )
CALL DLACPY( 'A', N, N, M( N+1, 1 ), LDM, WORK( INDW ), TWON )
CALL DLACPY( 'L', N, N, M( N+1, N+1 ), LDM, &
     WORK( INDW + TWON*N ), TWON )
CALL DTRMM( 'L', 'L', 'T', 'N', N, TWON, ONE, M, LDM, &
     WORK( INDW ), TWON )
DO J = 0, N-1
   DO I = 0, J
      WORK( INDW + I + J*TWON ) = WORK( INDW + I + J*TWON ) &
           -WORK( INDW + J + I*TWON )
   END DO
END DO
!
!     Tridiagonalization: U**T * W * U = T.
!     W is overwritten by the orthogonal matrix U.
!
CALL DSSTRD( 'U', TWON, WORK( INDW ), TWON, WORK( INDE ), &
     WORK( INDTAU ), WORK( INDWORK ), LLWORK, ITMP )
!
!     Formulate L * U.
!
CALL DLASET( 'A', TWON, TWON, ZERO, ZERO, WORK( INDU ), TWON )
CALL DLACPY( 'L', TWON, TWON, M, LDM, WORK( INDU ), TWON )
CALL DORMTR( 'R', 'U', 'N', TWON, TWON, WORK( INDW ), TWON, &
     WORK( INDTAU ), WORK( INDU ), TWON, WORK( INDWORK ), LLWORK, &
     ITMP )
!
!     Diagonalize the tridiagonal matrix
!
!        D**H * ( -i*T ) * D = V * diag( lambda ) * V**T,
!
!     where D = diag{1,i,i^2,i^3,...}.
!     Only the positive half of the eigenpairs are computed.
!
DO I = 0, TWON-1
   WORK( INDD + I ) = ZERO
END DO
!
CALL DSTEBZ( 'I', 'B', TWON, ZERO, ZERO, N+1, TWON, ABSTOL, &
     WORK( INDD ), WORK( INDE ), DIMV, NSPLIT, LAMBDA, &
     IWORK( INDIBLOCK ), IWORK( INDISPLIT ), WORK( INDWORK ), &
     IWORK( INDIWORK ), ITMP )
IF ( ITMP .NE. 0 ) THEN
   INFO = ITMP
   WRITE( *,* ) 'DSTEBZ fails with INFO =', INFO
   RETURN
END IF
CALL DSTEIN( TWON, WORK( INDD ), WORK( INDE ), DIMV, LAMBDA, &
     IWORK( INDIBLOCK ), IWORK( INDISPLIT ), WORK( INDV ), &
     TWON, WORK( INDWORK ), IWORK( INDIWORK ), IWORK( INDIFAIL ), &
     ITMP )
IF ( ITMP .NE. 0 ) THEN
   INFO = ITMP
   WRITE( *,* ) 'DSTEIN fails with INFO =', INFO
   RETURN
END IF
CALL DSORTEIG( TWON, DIMV, LAMBDA, WORK( INDV ), TWON, 0 )
!
!     Restore the eigenvectors of H:
!
!        X = Q * L**{-T} * U * D * V * Lambda**{1/2},
!        Y = Q * L * U * D * V * Lambda**{-1/2},
!
!     where
!
!        Q = sqrt(1/2) * [ I_n, -i*I_n;
!                          I_n,  i*I_n ].
!
!     Scale V by (2*Lambda)**{-1/2}.
!
DO I = 0, N-1
   CALL DSCAL( TWON, DSQRT( HALF/LAMBDA( I+1 ) ), &
        WORK( INDV + I*TWON ), 1 )
END DO
!
!     Formulate
!
!        Y = Q * L * U * D * V * Lambda**{-1/2},
!        X = Q * L**{-T} * U * D * V * Lambda**{1/2},
!
!     Once Y is calculated, X is constructed from Y through
!
!        Y = [ Y1; Y2 ], X = [ Y1; -Y2 ].
!
!     Use WORK and M as workspace for real and imaginary parts,
!     respectively.
!
CALL ZBSAUX1( N, WORK( INDV ), TWON, X, LDX, WORK( INDU ), TWON, &
     WORK( INDWORK ), TWON, M, LDM )
CALL ZLASCL( 'G', N, N, -ONE, ONE, N, N, X( N+1, 1 ), LDX, ITMP )
!
WORK( 1 ) = DBLE( LWKOPT )
IWORK( 1 ) = LIWKOPT
!
RETURN
!
!     End of ZBSSOLVER1().
!
END
