PROGRAM TEST_DSSMM
!
IMPLICIT NONE
!
INTEGER            NMAX, OFFSET, MATLEN
PARAMETER          ( NMAX = 32, OFFSET = 3, &
                     MATLEN = (NMAX+OFFSET)*NMAX )
DOUBLE PRECISION   A( MATLEN ), B( MATLEN ), &
                   C( MATLEN ), CC( MATLEN )
INTEGER            I, J, K, L, M, N, LDA, LDB, LDC, PASS, TOTAL
CHARACTER          SIDE, UPLO
INTEGER            ISEED( 4 )
DOUBLE PRECISION   ALPHA, BETA, EPS, ERR, TOL, AMAX, BMAX, CMAX
INTRINSIC          DBLE
EXTERNAL           DCOPY, DGEMM, DLACPY, DLAMCH, &
                   DGENMAT, DGENSSMAT, DSSMM, DERRMAT
DOUBLE PRECISION   DLAMCH, DERRMAT
!
ISEED( 1 ) = 0
ISEED( 2 ) = 1
ISEED( 3 ) = 2
ISEED( 4 ) = 3
!
TOTAL = 0
PASS = 0
EPS = DLAMCH( 'P' )
!
DO M = 1, 28, 9
DO N = 1, 28, 9
DO LDA = MAX( M, N ), MAX( M, N )+OFFSET, OFFSET
DO LDB = M, M+OFFSET, OFFSET
DO LDC = M, M+OFFSET, OFFSET
DO I = -1, 2
DO J = -1, 2
DO K = 0, 1
DO L = 0, 1
   IF ( K .LT. 1 ) THEN
      UPLO = 'L'
   ELSE
      UPLO = 'U'
   END IF
   IF ( L .LT. 1 ) THEN
      SIDE = 'L'
      CALL DGENSSMAT( M, A, LDA, ISEED, AMAX )
   ELSE
      SIDE = 'R'
      CALL DGENSSMAT( N, A, LDA, ISEED, AMAX )
   END IF
   ALPHA = DBLE( I )
   BETA = DBLE( J )
!
   CALL DGENMAT( M, N, B, LDB, ISEED, BMAX )
   CALL DGENMAT( M, N, C, LDC, ISEED, CMAX )
   CALL DLACPY( 'A', M, N, C, LDC, CC, LDC )
!
   CALL DSSMM( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, &
        C, LDC )
   IF ( L .LT. 1 ) THEN
      CALL DGEMM( 'N', 'N', M, N, M, ALPHA, A, LDA, B, LDB, &
           BETA, CC, LDC )
   ELSE
      CALL DGEMM( 'N', 'N', M, N, N, ALPHA, B, LDB, A, LDA, &
           BETA, CC, LDC )
   END IF
!
   ERR = DERRMAT( 'A', M, N, C, LDC, CC, LDC )
   TOL = 20 * EPS * N**1.5 * ( N * AMAX * BMAX + CMAX )
   IF ( ERR .LE. TOL ) PASS = PASS + 1
   TOTAL = TOTAL + 1
END DO
END DO
END DO
END DO
END DO
END DO
END DO
END DO
END DO
!
WRITE(*,*), '%', PASS, 'out of', TOTAL, 'tests passed!'
!
END
