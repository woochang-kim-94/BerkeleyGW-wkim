PROGRAM TEST_SSSR2
!
IMPLICIT NONE
!
REAL               ZERO
PARAMETER          ( ZERO = 0.0E+0 )
INTEGER            NMAX, OFFSET, MATLEN, VECLEN
PARAMETER          ( NMAX = 32, OFFSET = 3, &
                     MATLEN = (NMAX+OFFSET)*NMAX, &
                     VECLEN = NMAX*3+1 )
REAL               A( MATLEN ), AC( MATLEN ), &
                   X( VECLEN ), Y( VECLEN )
INTEGER            I, K, N, LDA, INCX, INCY, PASS, TOTAL
CHARACTER          UPLO, OPP
INTEGER            ISEED( 4 )
REAL               ALPHA, EPS, ERR, TOL, AMAX, XMAX, YMAX
INTRINSIC          REAL
EXTERNAL           SCOPY, SGER, SLACPY, SLAMCH, &
                   SGENSSMAT, SGENVEC, SSSR2, SERRMAT
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
DO N = 1, 28, 9
DO LDA = N, N+OFFSET, OFFSET
DO INCX = -3, 3, 2
DO INCY = -3, 3, 2
DO I = -1, 2
DO K = 0, 1
   IF ( K .LT. 1 ) THEN
      UPLO = 'L'
      OPP = 'U'
   ELSE
      UPLO = 'U'
      OPP = 'L'
   END IF
   ALPHA = REAL( I )
!
   CALL SGENSSMAT( N, A, LDA, ISEED, AMAX )
   CALL SGENVEC( N, X, INCX, ISEED, XMAX )
   CALL SGENVEC( N, Y, INCY, ISEED, YMAX )
   CALL SLACPY( 'A', N, N, A, LDA, AC, LDA )
!
   CALL SSSR2( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!
   ERR = SERRMAT( OPP, N, N, A, LDA, AC, LDA )
   IF ( ERR .NE. ZERO ) ERR = SLAMCH( 'O' )
!
   CALL SGER( N, N, ALPHA, X, INCX, Y, INCY, AC, LDA )
   CALL SGER( N, N, -ALPHA, Y, INCY, X, INCX, AC, LDA )
!
   ERR = ERR + SERRMAT( UPLO, N, N, A, LDA, AC, LDA )
   TOL = 10 * EPS * N * ( AMAX + XMAX * YMAX )
   IF ( ERR .LE. TOL ) PASS = PASS + 1
   TOTAL = TOTAL + 1
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
