PROGRAM TEST_SSSR2K
!
IMPLICIT NONE
!
REAL               ZERO, ONE
PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
INTEGER            NMAX, OFFSET, MATLEN
PARAMETER          ( NMAX = 32, OFFSET = 3, &
                     MATLEN = (NMAX+OFFSET)*NMAX )
REAL               A( MATLEN ), B( MATLEN ), &
                   C( MATLEN ), CC( MATLEN )
INTEGER            I, J, K, L, M, N, LDA, LDB, LDC, PASS, TOTAL
CHARACTER          UPLO, TRANS, OPP
INTEGER            ISEED( 4 )
REAL               ALPHA, BETA, EPS, ERR, TOL, AMAX, BMAX, CMAX
INTRINSIC          REAL
EXTERNAL           SGEMM, SLACPY, SLAMCH, &
                   SGENMAT, SGENSSMAT, SSSR2K, SERRMAT
REAL               SLAMCH, SERRMAT
!
ISEED( 1 ) = 0
ISEED( 2 ) = 1
ISEED( 3 ) = 2
ISEED( 4 ) = 3
!
TOTAL = 0
PASS = 0
EPS = SLAMCH( 'P' )
!
DO K = 1, 28, 9
DO N = 1, 28, 9
DO LDA = MAX( K, N ), MAX( K, N )+OFFSET, OFFSET
DO LDB = MAX( K, N ), MAX( K, N )+OFFSET, OFFSET
DO LDC = N, N+OFFSET, OFFSET
DO I = -1, 2
DO J = -1, 2
DO M = 0, 1
DO L = 0, 1
   IF ( M .LT. 1 ) THEN
      UPLO = 'L'
      OPP = 'U'
   ELSE
      UPLO = 'U'
      OPP = 'L'
   END IF
   IF ( L .LT. 1 ) THEN
      TRANS = 'N'
      CALL SGENMAT( N, K, A, LDA, ISEED, AMAX )
      CALL SGENMAT( N, K, B, LDB, ISEED, BMAX )
   ELSE
      TRANS = 'T'
      CALL SGENMAT( K, N, A, LDA, ISEED, AMAX )
      CALL SGENMAT( K, N, B, LDB, ISEED, BMAX )
   END IF
   ALPHA = REAL( I )
   BETA = REAL( J )
   CALL SGENSSMAT( N, C, LDC, ISEED, CMAX )
   CALL SLACPY( 'A', N, N, C, LDC, CC, LDC )
!
   CALL SSSR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA, &
        C, LDC )
!
   ERR = SERRMAT( OPP, N, N, C, LDC, CC, LDC )
   IF ( ERR .NE. ZERO ) ERR = SLAMCH( 'O' )
!
   IF ( L .LT. 1 ) THEN
      CALL SGEMM( 'N', 'T', N, N, K, ALPHA, A, LDA, B, LDB, &
           BETA, CC, LDC )
      CALL SGEMM( 'N', 'T', N, N, K, -ALPHA, B, LDB, A, LDA, &
           ONE, CC, LDC )
   ELSE
      CALL SGEMM( 'T', 'N', N, N, K, ALPHA, A, LDA, B, LDB, &
           BETA, CC, LDC )
      CALL SGEMM( 'T', 'N', N, N, K, -ALPHA, B, LDB, A, LDA, &
           ONE, CC, LDC )
   END IF
!
   ERR = SERRMAT( UPLO, N, N, C, LDC, CC, LDC )
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
