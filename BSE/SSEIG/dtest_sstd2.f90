PROGRAM TEST_DSSTD2
!
IMPLICIT NONE
!
DOUBLE PRECISION   ZERO, ONE
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
INTEGER            NMAX, OFFSET, MATLEN, VECLEN
PARAMETER          ( NMAX = 128, OFFSET = 3, &
                     MATLEN = (NMAX+OFFSET)*NMAX, &
                     VECLEN = NMAX-1 )
DOUBLE PRECISION   A( MATLEN ), AC( MATLEN ), T( MATLEN ), &
                   W( MATLEN ), &
                   E( VECLEN ), TAU( VECLEN )
INTEGER            K, N, LDA, LDT, LDW, INFO, OFFDIAG, PASS, TOTAL
CHARACTER          UPLO, OPP
INTEGER            ISEED( 4 )
DOUBLE PRECISION   ALPHA, EPS, ERR, TOL, AMAX
INTRINSIC          MIN
EXTERNAL           DAXPY, DGEMM, DLACPY, DLASET, &
                   DLAMCH, DLANGE, DORGTR, &
                   DGENSSMAT, DSSTD2, DERRMAT, DERRVEC
DOUBLE PRECISION   DLAMCH, DLANGE, DERRMAT, DERRVEC
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
DO N = 1, 78, 11
DO LDA = N, N+OFFSET
DO K = 0, 1
   IF ( K .LT. 1 ) THEN
      UPLO = 'L'
      OPP = 'U'
      OFFDIAG = MIN( 2, N )
      ALPHA = -ONE
   ELSE
      UPLO = 'U'
      OPP = 'L'
      OFFDIAG = MIN( LDA+1, LDA*N )
      ALPHA = ONE
   END IF
!
   CALL DGENSSMAT( N, A, LDA, ISEED, AMAX )
   CALL DLACPY( 'A', N, N, A, LDA, AC, LDA )
!
   CALL DSSTD2( UPLO, N, A, LDA, E, TAU, INFO )
!
   ERR = DERRMAT( OPP, N, N, A, LDA, AC, LDA )
   ERR = ERR + DERRVEC( N-1, A( OFFDIAG ), LDA+1, ALPHA, E, 1 )
   IF ( INFO .NE. 0 ) ERR = DLAMCH( 'O' )
   IF ( ERR .NE. ZERO ) ERR = DLAMCH( 'O' )
!
   CALL DORGTR( UPLO, N, A, LDA, TAU, W, MATLEN, INFO )
   IF ( INFO .NE. 0 ) ERR = DLAMCH( 'O' )
   LDW = N
   CALL DLASET( 'A', N, N, ZERO, -ONE, W, LDW )
   CALL DGEMM( 'T', 'N', N, N, N, ONE, A, LDA, A, LDA, ONE, W, &
        LDW )
   ERR = ERR + DLANGE( 'F', N, N, W, LDW, T )
!
   LDT = N
   CALL DLASET( 'A', N, N, ZERO, ZERO, T, LDT )
   CALL DAXPY( N-1, ONE, E, 1, T( LDT+1 ), LDT+1 )
   CALL DAXPY( N-1, -ONE, E, 1, T( 2 ), LDT+1 )
   CALL DGEMM( 'N', 'N', N, N, N, ONE, AC, LDA, A, LDA, ZERO, W, &
        LDW )
   CALL DGEMM( 'T', 'N', N, N, N, -ONE, A, LDA, W, LDW, ONE, T, &
        LDT )
   ERR = ERR + DLANGE( 'F', N, N, T, LDT, W )
!
   TOL = 10 * EPS * N * AMAX
   IF ( ERR .LE. TOL ) PASS = PASS + 1
   TOTAL = TOTAL + 1
END DO
END DO
END DO
!
WRITE(*,*), '%', PASS, 'out of', TOTAL, 'tests passed!'
!
END
