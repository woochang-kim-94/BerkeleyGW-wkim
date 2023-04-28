PROGRAM TEST_DSSR2K
!
IMPLICIT NONE
!
DOUBLE PRECISION   ZERO, ONE
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
INTEGER            NMAX, OFFSET, MATLEN
PARAMETER          ( NMAX = 32, OFFSET = 3, &
                     MATLEN = (NMAX+OFFSET)*NMAX )
DOUBLE PRECISION   A( MATLEN ), B( MATLEN ), &
                   C( MATLEN ), CC( MATLEN )
INTEGER            I, J, K, L, M, N, LDA, LDB, LDC, PASS, TOTAL
CHARACTER          UPLO, TRANS, OPP
INTEGER            ISEED( 4 )
DOUBLE PRECISION   ALPHA, BETA, EPS, ERR, TOL, AMAX, BMAX, CMAX
INTRINSIC          DBLE
EXTERNAL           DGEMM, DLACPY, DLAMCH, &
                   DGENMAT, DGENSSMAT, DSSR2K, DERRMAT
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
      CALL DGENMAT( N, K, A, LDA, ISEED, AMAX )
      CALL DGENMAT( N, K, B, LDB, ISEED, BMAX )
   ELSE
      TRANS = 'T'
      CALL DGENMAT( K, N, A, LDA, ISEED, AMAX )
      CALL DGENMAT( K, N, B, LDB, ISEED, BMAX )
   END IF
   ALPHA = DBLE( I )
   BETA = DBLE( J )
   CALL DGENSSMAT( N, C, LDC, ISEED, CMAX )
   CALL DLACPY( 'A', N, N, C, LDC, CC, LDC )
!
   CALL DSSR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA, &
        C, LDC )
!
   ERR = DERRMAT( OPP, N, N, C, LDC, CC, LDC )
   IF ( ERR .NE. ZERO ) ERR = DLAMCH( 'O' )
!
   IF ( L .LT. 1 ) THEN
      CALL DGEMM( 'N', 'T', N, N, K, ALPHA, A, LDA, B, LDB, &
           BETA, CC, LDC )
      CALL DGEMM( 'N', 'T', N, N, K, -ALPHA, B, LDB, A, LDA, &
           ONE, CC, LDC )
   ELSE
      CALL DGEMM( 'T', 'N', N, N, K, ALPHA, A, LDA, B, LDB, &
           BETA, CC, LDC )
      CALL DGEMM( 'T', 'N', N, N, K, -ALPHA, B, LDB, A, LDA, &
           ONE, CC, LDC )
   END IF
!
   ERR = DERRMAT( UPLO, N, N, C, LDC, CC, LDC )
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
