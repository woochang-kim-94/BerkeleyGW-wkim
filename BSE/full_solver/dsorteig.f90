SUBROUTINE DSORTEIG( N, M, LAMBDA, V, LDV, MODE )
!
IMPLICIT NONE
!
!     .. Scalar Arguments ..
INTEGER            N, M, LDV, MODE
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION   LAMBDA( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  DSORTEIG() is an auxiliary routine.
!  It sorts the eigenvalues stored in lambda and the corresponding
!  eigenvectors stored in V.
!
!  On entry, lambda is of length m, and V is of size n-by-m.
!
!  On exit,
!  if MODE = 0, then
!     lambda( 1 ) <= lambda( 2 ) <= ... <= lambda( m );
!  if MODE = 1, then
!     lambda( 1 ) >= lambda( 2 ) >= ... >= lambda( m ).
!
!  V( :, j ) always contains the eigenvector corresponds to lambda( j ).
!
!  No argument check is performed, i.e.,
!  all arguments are assumed to be valid.
!
!  Arguments
!  =========
!
!  N       (input) INGEGER
!          The length of each eigenvector.
!
!  M       (input) INGEGER
!          The number of eigenpairs.
!
!  LAMBDA  (input/output) DOUBLE PRECISION array, dimension (M)
!          The eigenvalues to be sorted.
!
!  V       (input/output) DOUBLE PRECISION array, dimension (LDV, M)
!          On entry/exit, V( :, j ) always contains the eigenvector
!          corresponds to lambda( j ).
!
!  LDV     (input) INGEGER
!          The leading dimension of the array V. LDV >= max(1,N).
!
!  MODE    (input) INGEGER
!          = 0:  lambda is sorted in ascending order.
!          = 1:  lambda is sorted in descending order.
!
!  =====================================================================
!
!     .. Local Scalars ..
INTEGER            I, J, JJ
DOUBLE PRECISION   TMP
!     ..
!     .. External Subroutines ..
EXTERNAL           DSWAP
!
IF ( MODE .EQ. 0 ) THEN
   DO J = 1, M
      I = J
      TMP = LAMBDA( J )
      DO JJ = J+1, M
         IF ( LAMBDA( JJ ) .LT. TMP ) I = JJ
      END DO
      IF ( I .NE. J ) THEN
         LAMBDA( J ) = LAMBDA( I )
         LAMBDA( I ) = TMP
         CALL DSWAP( N, V( 1, J ), 1, V( 1, I ), 1 )
      END IF
   END DO
ELSE
   DO J = 1, M
      I = J
      TMP = LAMBDA( J )
      DO JJ = J+1, M
         IF ( LAMBDA( JJ ) .GT. TMP ) I = JJ
      END DO
      IF ( I .NE. J ) THEN
         LAMBDA( J ) = LAMBDA( I )
         LAMBDA( I ) = TMP
         CALL DSWAP( N, V( 1, J ), 1, V( 1, I ), 1 )
      END IF
   END DO
END IF
!
RETURN
!
!     End of DSORTEIG().
!
END
