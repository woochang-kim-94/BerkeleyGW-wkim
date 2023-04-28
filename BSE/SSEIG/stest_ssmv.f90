PROGRAM TEST_SSSMV
!
IMPLICIT NONE
!
REAL               ONE
PARAMETER          ( ONE = 1.0E+0 )
INTEGER            NMAX, OFFSET, MATLEN, VECLEN
PARAMETER          ( NMAX = 32, OFFSET = 3, &
                     MATLEN = (NMAX+OFFSET)*NMAX, &
                     VECLEN = NMAX*3+1 )
REAL               A( MATLEN ), X( VECLEN ), &
                   Y( VECLEN ), YC( VECLEN )
INTEGER            I, J, K, N, LDA, INCX, INCY, PASS, TOTAL
CHARACTER          UPLO
INTEGER            ISEED( 4 )
REAL               ALPHA, BETA, EPS, ERR, TOL, AMAX, XMAX, YMAX
INTRINSIC          REAL
EXTERNAL           SCOPY, SGEMV, SLAMCH, &
                   SGENSSMAT, SGENVEC, SSSMV, SERRVEC
REAL               SLAMCH, SERRVEC
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
DO J = -1, 2
DO K = 0, 1
   IF ( K .LT. 1 ) THEN
      UPLO = 'L'
   ELSE
      UPLO = 'U'
   END IF
   ALPHA = REAL( I )
   BETA = REAL( J )
!
   CALL SGENSSMAT( N, A, LDA, ISEED, AMAX )
   CALL SGENVEC( N, X, INCX, ISEED, XMAX )
   CALL SGENVEC( N, Y, INCY, ISEED, YMAX )
   CALL SCOPY( N, Y, INCY, YC, INCY )
!
   CALL SSSMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
   CALL SGEMV( 'N', N, N, ALPHA, A, LDA, X, INCX, BETA, &
        YC, INCY )
!
   ERR = SERRVEC( N, Y, INCY, ONE, YC, INCY )
   TOL = 10 * EPS * N * ( N * AMAX * XMAX + YMAX )
   IF ( ERR .LE. TOL ) PASS = PASS + 1
   TOTAL = TOTAL + 1
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
