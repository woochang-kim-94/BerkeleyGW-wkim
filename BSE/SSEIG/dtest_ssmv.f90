PROGRAM TEST_DSSMV
!
IMPLICIT NONE
!
DOUBLE PRECISION   ONE
PARAMETER          ( ONE = 1.0D+0 )
INTEGER            NMAX, OFFSET, MATLEN, VECLEN
PARAMETER          ( NMAX = 32, OFFSET = 3, &
                     MATLEN = (NMAX+OFFSET)*NMAX, &
                     VECLEN = NMAX*3+1 )
DOUBLE PRECISION   A( MATLEN ), X( VECLEN ), &
                   Y( VECLEN ), YC( VECLEN )
INTEGER            I, J, K, N, LDA, INCX, INCY, PASS, TOTAL
CHARACTER          UPLO
INTEGER            ISEED( 4 )
DOUBLE PRECISION   ALPHA, BETA, EPS, ERR, TOL, AMAX, XMAX, YMAX
INTRINSIC          DBLE
EXTERNAL           DCOPY, DGEMV, DLAMCH, &
                   DGENSSMAT, DGENVEC, DSSMV, DERRVEC
DOUBLE PRECISION   DLAMCH, DERRVEC
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
   ALPHA = DBLE( I )
   BETA = DBLE( J )
!
   CALL DGENSSMAT( N, A, LDA, ISEED, AMAX )
   CALL DGENVEC( N, X, INCX, ISEED, XMAX )
   CALL DGENVEC( N, Y, INCY, ISEED, YMAX )
   CALL DCOPY( N, Y, INCY, YC, INCY )
!
   CALL DSSMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
   CALL DGEMV( 'N', N, N, ALPHA, A, LDA, X, INCX, BETA, &
        YC, INCY )
!
   ERR = DERRVEC( N, Y, INCY, ONE, YC, INCY )
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
