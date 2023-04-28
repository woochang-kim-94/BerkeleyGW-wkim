PROGRAM TEST_PDSSMV
!
IMPLICIT NONE
!
DOUBLE PRECISION   ONE
PARAMETER          ( ONE = 1.0D+0 )
INTEGER            NMAX, MATLEN
PARAMETER          ( NMAX = 32, MATLEN = NMAX*NMAX )
DOUBLE PRECISION   A( MATLEN ), X( MATLEN ), &
                   Y( MATLEN ), YC( MATLEN ), &
                   WORK( MATLEN )
INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM, &
                   SYS_NPROCS
INTEGER            I, J, K, N, NB, LLD, INCX, INCY, PASS, TOTAL
CHARACTER          UPLO
INTEGER            DESCA( 9 ), DESCX( 9 ), DESCY( 9 ), ISEED( 4 )
DOUBLE PRECISION   ALPHA, BETA, EPS, ERR, TOL, AMAX, XMAX, YMAX
INTRINSIC          DBLE
EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_EXIT, &
                   BLACS_GRIDINIT, BLACS_GRIDINFO, BLACS_GRIDEXIT, &
                   NOMROC, DLAMCH, PDCOPY, PDGEMV, &
                   PDGENMAT, PDGENSSMAT, PDSSMV, PDERRVEC
INTEGER            NUMROC
DOUBLE PRECISION   DLAMCH, PDERRVEC
!
CALL BLACS_PINFO( IAM, SYS_NPROCS )
NPROW = INT( SQRT( DBLE( SYS_NPROCS ) ) )
NPCOL = SYS_NPROCS / NPROW
CALL BLACS_GET( 0, 0, ICTXT )
CALL BLACS_GRIDINIT( ICTXT, '2D', NPROW, NPCOL )
CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
ISEED( 1 ) = 0
ISEED( 2 ) = 1
ISEED( 3 ) = 2
ISEED( 4 ) = 3
!
IF ( ICTXT .LT. 0 ) GOTO 9999
!
TOTAL = 0
PASS = 0
EPS = DLAMCH( 'P' )
!
DO N = 1, 28, 9
DO NB = 1, 31, 6
DO INCX = 1, N
IF ( INCX .NE. 1 .AND. INCX .NE. N ) CYCLE
DO INCY = 1, N
IF ( INCY .NE. 1 .AND. INCY .NE. N ) CYCLE
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
   LLD = MAX( 1, NUMROC( N, NB, MYROW, 0, NPROW ) )
   CALL DESCSET( DESCA, N, N, NB, NB, 0, 0, ICTXT, LLD )
   CALL DESCSET( DESCX, N, N, NB, NB, 0, 0, ICTXT, LLD )
   CALL DESCSET( DESCY, N, N, NB, NB, 0, 0, ICTXT, LLD )
   CALL PDGENSSMAT( N, A, 1, 1, DESCA, ISEED, AMAX )
   CALL PDGENMAT( N, N, X, 1, 1, DESCX, ISEED, XMAX )
   CALL PDGENMAT( N, N, Y, 1, 1, DESCY, ISEED, YMAX )
   CALL PDCOPY( N, Y, 1, 1, DESCY, INCY, YC, 1, 1, DESCY, &
        INCY )
!
   CALL PDSSMV( UPLO, N, ALPHA, A, 1, 1, DESCA, X, 1, 1, DESCX, &
        INCX, BETA, Y, 1, 1, DESCY, INCY, WORK )
   CALL PDGEMV( 'N', N, N, ALPHA, A, 1, 1, DESCA, X, 1, 1, &
        DESCX, INCX, BETA, YC, 1, 1, DESCY, INCY )
!
   ERR = PDERRVEC( N, Y, 1, 1, DESCY, INCY, ONE, YC, 1, 1, &
         DESCY, INCY )
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
IF ( MYROW + MYCOL .EQ. 0 ) THEN
   WRITE(*,*), '%', PASS, 'out of', TOTAL, 'tests passed!'
END IF
!
CALL BLACS_GRIDEXIT( ICTXT )
!
9999 CONTINUE
CALL BLACS_EXIT( 0 )
!
END