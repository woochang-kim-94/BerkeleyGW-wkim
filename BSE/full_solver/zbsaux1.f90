SUBROUTINE ZBSAUX1( N, V, LDV, X, LDX, Z, LDZ, ZR, LDZR, ZI, &
                    LDZI )
!
IMPLICIT NONE
!
!     .. Scalar Arguments ..
INTEGER            N, LDV, LDX, LDZ, LDZR, LDZI
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION   V( LDV, * ), Z( LDZ, * ), ZR( LDZR, * ), &
                   ZI( LDZI, * )
COMPLEX*16         X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  ZBSAUX1() is an auxiliary routine called by ZBSSOLVER1().
!  It applies three 2n-by-2n or 2n-by-n unitary transformations
!
!     X = Q * Z * D * V,
!
!  where D = diag{1,i,i^2,i^3,...},
!
!        Q = 1/sqrt(2) * [ I_n, -i*I_n;
!                          I_n,  i*I_n ],
!
!  and V is a given 2n-by-2n real orthogonal matrix.
!  On entry, V has already been prescaled by the factor 1/sqrt(2) in Q.
!
!  ZR and ZI are provided as workspace to store the real and imaginary
!  parts of Q * Z * D, respectively, i.e.,
!
!     Q * Z * D = ZR + i*ZI.
!
!  Z is destroyed on exit.
!
!  No argument check is performed, i.e.,
!  all arguments are assumed to be valid.
!
!  Arguments
!  =========
!
!  N       (input) INGEGER
!          2*N is the number of rows and columns of Z.
!
!  V       (input) DOUBLE PRECISION array, dimension (LDV, N)
!          The eigenvectors before similarity transformations.
!
!  LDV     (input) INGEGER
!          The leading dimension of the array V. LDV >= max(1,2*N).
!
!  X       (output) COMPLEX*16 array, dimension (LDX, N)
!          The eigenvectors after similarity transformations.
!
!  LDX     (input) INGEGER
!          The leading dimension of the array X. LDX >= max(1,2*N).
!
!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, 2*N)
!          The transformation matrix.
!
!  LDZ     (input) INGEGER
!          The leading dimension of the array Z. LDZ >= max(1,2*N).
!
!  ZR      (workspacek/output) DOUBLE PRECISION array,
!          dimension (LDZR, 2*N)
!          Work space.
!
!  LDZR    (input) INGEGER
!          The leading dimension of the array ZR. LDZR >= max(1,2*N).
!
!  ZI      (workspacek/output) DOUBLE PRECISION array,
!          dimension (LDZI, 2*N)
!          Work space.
!
!  LDZI    (input) INGEGER
!          The leading dimension of the array ZR. LDZI >= max(1,2*N).
!
!  =====================================================================
!
!     .. Parameters ..
DOUBLE PRECISION   ZERO, ONE
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
INTEGER            I, J
!     ..
!     .. Intrinsic Functions ..
INTRINSIC          DBLE, DCMPLX, MOD
!     ..
!     .. External Subroutines ..
EXTERNAL           DCOPY, DGEMM, DSCAL
!     ..
!     .. Executable Statements ..
!
!     Apply Q from the left and D from the right.
!
DO J = 1, 2*N
   IF ( MOD( J-1, 4 ) .EQ. 0 ) THEN
      CALL DCOPY( N, Z( 1, J ), 1, ZR( 1, J ), 1 )
      CALL DCOPY( N, Z( 1, J ), 1, ZR( N+1, J ), 1 )
      CALL DCOPY( N, Z( N+1, J ), 1, ZI( 1, J ), 1 )
      CALL DSCAL( N, -ONE, ZI( 1, J ), 1 )
      CALL DCOPY( N, Z( N+1, J ), 1, ZI( N+1, J ), 1 )
   ELSE IF ( MOD( J-1, 4 ) .EQ. 1 ) THEN
      CALL DCOPY( N, Z( 1, J ), 1, ZI( 1, J ), 1 )
      CALL DCOPY( N, Z( 1, J ), 1, ZI( N+1, J ), 1 )
      CALL DCOPY( N, Z( N+1, J ), 1, ZR( 1, J ), 1 )
      CALL DCOPY( N, Z( N+1, J ), 1, ZR( N+1, J ), 1 )
      CALL DSCAL( N, -ONE, ZR( N+1, J ), 1 )
   ELSE IF ( MOD( J-1, 4 ) .EQ. 2 ) THEN
      CALL DCOPY( N, Z( 1, J ), 1, ZR( 1, J ), 1 )
      CALL DSCAL( N, -ONE, ZR( 1, J ), 1 )
      CALL DCOPY( N, Z( 1, J ), 1, ZR( N+1, J ), 1 )
      CALL DSCAL( N, -ONE, ZR( N+1, J ), 1 )
      CALL DCOPY( N, Z( N+1, J ), 1, ZI( 1, J ), 1 )
      CALL DCOPY( N, Z( N+1, J ), 1, ZI( N+1, J ), 1 )
      CALL DSCAL( N, -ONE, ZI( N+1, J ), 1 )
   ELSE !IF ( MOD( J-1, 4 ) .EQ. 3 ) THEN
      CALL DCOPY( N, Z( 1, J ), 1, ZI( 1, J ), 1 )
      CALL DSCAL( N, -ONE, ZI( 1, J ), 1 )
      CALL DCOPY( N, Z( 1, J ), 1, ZI( N+1, J ), 1 )
      CALL DSCAL( N, -ONE, ZI( N+1, J ), 1 )
      CALL DCOPY( N, Z( N+1, J ), 1, ZR( 1, J ), 1 )
      CALL DSCAL( N, -ONE, ZR( 1, J ), 1 )
      CALL DCOPY( N, Z( N+1, J ), 1, ZR( N+1, J ), 1 )
   END IF
END DO
!
!     Apply V from the right.
!
CALL DGEMM( 'N', 'N', 2*N, N, 2*N, ONE, ZR, LDZR, V, LDV, ZERO, &
     Z, LDZ )
DO J = 1, N
   DO I = 1, 2*N
      X( I, J ) = DCMPLX( Z( I, J ) )
   END DO
END DO
CALL DGEMM( 'N', 'N', 2*N, N, 2*N, ONE, ZI, LDZI, V, LDV, ZERO, &
     Z, LDZ )
DO J = 1, N
   DO I = 1, 2*N
      X( I, J ) = DCMPLX( DBLE( X( I, J ) ), Z( I, J ) )
   END DO
END DO
!
RETURN
!
!     End of ZBSAUX1().
!
END
