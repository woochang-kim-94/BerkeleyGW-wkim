PROGRAM TEST_PSSSTRD
!
IMPLICIT NONE
!
REAL               ZERO, ONE
PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
INTEGER            NMAX, MATLEN
PARAMETER          ( NMAX = 128, MATLEN = NMAX*NMAX )
REAL               A( MATLEN ), AC( MATLEN ), T( MATLEN ), &
                   U( MATLEN ), W( MATLEN ), &
                   E( NMAX ), E2( NMAX ), TAU( NMAX ), &
                   WORK( 2*MATLEN )
INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM, &
                   SYS_NPROCS
INTEGER            I, K, N, NB, LLDA, LLDT, LLDU, LLDW, &
                   LWORK, INFO, PASS, TOTAL
CHARACTER          UPLO, TRANS, OPP
INTEGER            DESCA( 9 ), DESCT( 9 ), DESCU( 9 ), DESCW( 9 ), &
                   ISEED( 4 )
REAL               ALPHA, EPS, ERR, TOL, AMAX
INTRINSIC          REAL, INT, MAX
EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_EXIT, &
                   BLACS_GRIDINIT, BLACS_GRIDINFO, BLACS_GRIDEXIT, &
                   NOMROC, SLAMCH, PSGEMM, PSLACPY, &
                   PSORMTR, PSELGET, PSELSET, PSSSTRD, &
                   PSLANGE, PSGENSSMAT, PSERRMAT
INTEGER            NUMROC
REAL               SLAMCH, PSLANGE, PSERRMAT
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
EPS = SLAMCH( 'P' )
!
DO N = 1, 78, 11
DO NB = 1, 64, 9
DO K = 0, 1
   IF ( K .LT. 1 ) THEN
      UPLO = 'L'
      OPP = 'U'
      ALPHA = -ONE
   ELSE
      UPLO = 'U'
      OPP = 'L'
      ALPHA = ONE
   END IF
!
   LLDA = MAX( 1, NUMROC( N, NB, MYROW, 0, NPROW ) )
   CALL DESCSET( DESCA, N, N, NB, NB, 0, 0, ICTXT, LLDA )
   CALL PSGENSSMAT( N, A, 1, 1, DESCA, ISEED, AMAX )
   CALL PSLACPY( 'A', N, N, A, 1, 1, DESCA, AC, 1, 1, DESCA )
!
   CALL PSSSTRD( UPLO, N, A, 1, 1, DESCA, E, TAU, WORK, -1, INFO )
   LWORK = INT( WORK( 1 ) )
   CALL PSSSTRD( UPLO, N, A, 1, 1, DESCA, E, TAU, WORK, LWORK, &
        INFO )
   IF ( K .LT. 1 ) THEN
      DO I = 1, N-1
         CALL PSELGET( 'A', ' ', E2( I ), A, I+1, I, DESCA )
         E2( I ) = -E2( I )
      END DO
   ELSE
      DO I = 1, N-1
         CALL PSELGET( 'A', ' ', E2( I ), A, I, I+1, DESCA )
      END DO
   END IF
!
   ERR = PSERRMAT( OPP, N, N, A, 1, 1, DESCA, AC, 1, 1, DESCA )
   IF ( INFO .NE. 0 ) ERR = SLAMCH( 'O' )
   IF ( ERR .NE. ZERO ) ERR = SLAMCH( 'O' )
!
   LLDU = LLDA
   CALL DESCSET( DESCU, N, N, NB, NB, 0, 0, ICTXT, LLDU )
   CALL PSLASET( 'A', N, N, ZERO, ONE, U, 1, 1, DESCU )
   LWORK = 2*MATLEN
   CALL PSORMTR( 'R', UPLO, 'N', N, N, A, 1, 1, DESCA, TAU, U, 1, &
        1, DESCU, WORK, LWORK, INFO )
   IF ( INFO .NE. 0 ) ERR = SLAMCH( 'O' )
   LLDW = LLDA
   CALL DESCSET( DESCW, N, N, NB, NB, 0, 0, ICTXT, LLDW )
   CALL PSLASET( 'A', N, N, ZERO, -ONE, W, 1, 1, DESCW )
   CALL PSGEMM( 'T', 'N', N, N, N, ONE, U, 1, 1, DESCU, U, 1, 1, &
        DESCU, ONE, W, 1, 1, DESCW )
   ERR = ERR + PSLANGE( 'F', N, N, W, 1, 1, DESCW, WORK )
!
   LLDT = LLDA
   CALL DESCSET( DESCT, N, N, NB, NB, 0, 0, ICTXT, LLDT )
   CALL PSLASET( 'A', N, N, ZERO, ZERO, T, 1, 1, DESCT )
   DO I = 1, N-1
      CALL PSELSET( T, I, I+1, DESCT, E2( I ) )
      CALL PSELSET( T, I+1, I, DESCT, -E2( I ) )
   END DO
   CALL PSGEMM( 'N', 'N', N, N, N, ONE, AC, 1, 1, DESCA, U, 1, 1, &
        DESCU, ZERO, W, 1, 1, DESCW )
   CALL PSGEMM( 'T', 'N', N, N, N, -ONE, U, 1, 1, DESCU, W, 1, 1, &
        DESCW, ONE, T, 1, 1, DESCT )
   ERR = ERR + PSLANGE( 'F', N, N, T, 1, 1, DESCT, WORK )
!
   TOL = 10 * EPS * N * AMAX
   IF ( ERR .LE. TOL ) PASS = PASS + 1
   TOTAL = TOTAL + 1
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
