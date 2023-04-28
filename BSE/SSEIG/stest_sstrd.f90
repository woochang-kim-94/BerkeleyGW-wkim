PROGRAM TEST_SSSTRD
!
IMPLICIT NONE
!
REAL               ZERO, ONE
PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
INTEGER            NMAX, OFFSET, MATLEN, VECLEN
PARAMETER          ( NMAX = 128, OFFSET = 3, &
                     MATLEN = (NMAX+OFFSET)*NMAX, &
                     VECLEN = NMAX-1 )
REAL               A( MATLEN ), AC( MATLEN ), T( MATLEN ), &
                   W( MATLEN ), &
                   E( VECLEN ), TAU( VECLEN )
INTEGER            K, N, LDA, LDT, LDW, LWORK, INFO, OFFDIAG, &
                   PASS, TOTAL
CHARACTER          UPLO, OPP
INTEGER            ISEED( 4 )
REAL               ALPHA, EPS, ERR, TOL, AMAX
INTRINSIC          INT, MIN
EXTERNAL           SAXPY, SGEMM, SLACPY, SLASET, &
                   SLAMCH, SLANGE, SORGTR, &
                   SGENSSMAT, SSSTRD, SERRMAT, SERRVEC
REAL               SLAMCH, SLANGE, SERRMAT, SERRVEC
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
   CALL SGENSSMAT( N, A, LDA, ISEED, AMAX )
   CALL SLACPY( 'A', N, N, A, LDA, AC, LDA )
!
   CALL SSSTRD( UPLO, N, A, LDA, E, TAU, W, -1, INFO )
   LWORK = INT( W( 1 ) )
   CALL SSSTRD( UPLO, N, A, LDA, E, TAU, W, LWORK, INFO )
!
   ERR = SERRMAT( OPP, N, N, A, LDA, AC, LDA )
   ERR = ERR + SERRVEC( N-1, A( OFFDIAG ), LDA+1, ALPHA, E, 1 )
   IF ( INFO .NE. 0 ) ERR = SLAMCH( 'O' )
   IF ( ERR .NE. ZERO ) ERR = SLAMCH( 'O' )
!
   CALL SORGTR( UPLO, N, A, LDA, TAU, W, MATLEN, INFO )
   IF ( INFO .NE. 0 ) ERR = SLAMCH( 'O' )
   LDW = N
   CALL SLASET( 'A', N, N, ZERO, -ONE, W, LDW )
   CALL SGEMM( 'T', 'N', N, N, N, ONE, A, LDA, A, LDA, ONE, W, &
        LDW )
   ERR = ERR + SLANGE( 'F', N, N, W, LDW, T )
!
   LDT = N
   CALL SLASET( 'A', N, N, ZERO, ZERO, T, LDT )
   CALL SAXPY( N-1, ONE, E, 1, T( LDT+1 ), LDT+1 )
   CALL SAXPY( N-1, -ONE, E, 1, T( 2 ), LDT+1 )
   CALL SGEMM( 'N', 'N', N, N, N, ONE, AC, LDA, A, LDA, ZERO, W, &
        LDW )
   CALL SGEMM( 'T', 'N', N, N, N, -ONE, A, LDA, W, LDW, ONE, T, &
        LDT )
   ERR = ERR + SLANGE( 'F', N, N, T, LDT, W )
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
