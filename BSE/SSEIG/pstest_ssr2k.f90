PROGRAM TEST_PSSSR2K
!
IMPLICIT NONE
!
REAL               ZERO, ONE
PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
INTEGER            NMAX, MATLEN
PARAMETER          ( NMAX = 32, MATLEN = NMAX*NMAX )
REAL               A( MATLEN ), B( MATLEN ), &
                   C( MATLEN ), CC( MATLEN ), &
                   WORK( MATLEN )
INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM, &
                   SYS_NPROCS
INTEGER            I, J, K, L, M, N, NB, LLDA, LLDB, LLDC, &
                   PASS, TOTAL
CHARACTER          UPLO, TRANS, OPP
INTEGER            DESCA( 9 ), DESCB( 9 ), DESCC( 9 ), ISEED( 4 )
REAL               ALPHA, BETA, EPS, ERR, TOL, AMAX, BMAX, CMAX
INTRINSIC          REAL
EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_EXIT, &
                   BLACS_GRIDINIT, BLACS_GRIDINFO, BLACS_GRIDEXIT, &
                   NOMROC, SLAMCH, PSGEMM, PSLACPY, &
                   PSGENMAT, PSGENSSMAT, PSSSR2K, PSERRMAT
INTEGER            NUMROC
REAL               SLAMCH, PSERRMAT
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
DO K = 1, 28, 9
DO N = 1, 28, 9
DO NB = 1, 31, 6
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
      LLDA = MAX( 1, NUMROC( N, NB, MYROW, 0, NPROW ) )
      LLDB = LLDA
      CALL DESCSET( DESCA, N, K, NB, NB, 0, 0, ICTXT, LLDA )
      CALL DESCSET( DESCB, N, K, NB, NB, 0, 0, ICTXT, LLDB )
      CALL PSGENMAT( N, K, A, 1, 1, DESCA, ISEED, AMAX )
      CALL PSGENMAT( N, K, B, 1, 1, DESCB, ISEED, BMAX )
   ELSE
      TRANS = 'T'
      LLDA = MAX( 1, NUMROC( K, NB, MYROW, 0, NPROW ) )
      LLDB = LLDA
      CALL DESCSET( DESCA, K, N, NB, NB, 0, 0, ICTXT, LLDA )
      CALL DESCSET( DESCB, K, N, NB, NB, 0, 0, ICTXT, LLDB )
      CALL PSGENMAT( K, N, A, 1, 1, DESCA, ISEED, AMAX )
      CALL PSGENMAT( K, N, B, 1, 1, DESCB, ISEED, BMAX )
   END IF
   ALPHA = REAL( I )
   BETA = REAL( J )
!
   LLDC = MAX( 1, NUMROC( N, NB, MYROW, 0, NPROW ) )
   CALL DESCSET( DESCC, N, N, NB, NB, 0, 0, ICTXT, LLDC )
   CALL PSGENSSMAT( N, C, 1, 1, DESCC, ISEED, CMAX )
   CALL PSLACPY( 'A', N, N, C, 1, 1, DESCC, CC, 1, 1, DESCC )
!
   CALL PSSSR2K( UPLO, TRANS, N, K, ALPHA, A, 1, 1, DESCA, B, 1, &
        1, DESCB, BETA, C, 1, 1, DESCC, WORK )
!
   ERR = PSERRMAT( OPP, N, N, C, 1, 1, DESCC, CC, 1, 1, &
        DESCC )
   IF ( ERR .NE. ZERO ) ERR = SLAMCH( 'O' )
!
   IF ( L .LT. 1 ) THEN
      CALL PSGEMM( 'N', 'T', N, N, K, ALPHA, A, 1, 1, DESCA, &
           B, 1, 1, DESCB, BETA, CC, 1, 1, DESCC )
      CALL PSGEMM( 'N', 'T', N, N, K, -ALPHA, B, 1, 1, DESCB, &
           A, 1, 1, DESCA, ONE, CC, 1, 1, DESCC )
   ELSE
      CALL PSGEMM( 'T', 'N', N, N, K, ALPHA, A, 1, 1, DESCA, &
           B, 1, 1, DESCB, BETA, CC, 1, 1, DESCC )
      CALL PSGEMM( 'T', 'N', N, N, K, -ALPHA, B, 1, 1, DESCB, &
           A, 1, 1, DESCA, ONE, CC, 1, 1, DESCC )
   END IF
!
   ERR = PSERRMAT( UPLO, N, N, C, 1, 1, DESCC, CC, 1, 1, &
        DESCC )
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
