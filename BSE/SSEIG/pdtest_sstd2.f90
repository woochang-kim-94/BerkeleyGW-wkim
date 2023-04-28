PROGRAM TEST_PDSSTD2
!
IMPLICIT NONE
!
DOUBLE PRECISION   ZERO, ONE
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
INTEGER            NMAX, MATLEN
PARAMETER          ( NMAX = 128, MATLEN = NMAX*NMAX )
DOUBLE PRECISION   A( MATLEN ), AC( MATLEN ), T( MATLEN ), &
                   U( MATLEN ), W( MATLEN ), &
                   E( NMAX ), E2( NMAX ), TAU( NMAX ), &
                   WORK( 2*MATLEN )
INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM, &
                   SYS_NPROCS
INTEGER            I, J, K, KK, N, NB, LLDA, LLDT, LLDU, LLDW, &
                   LWORK, INFO, PASS, TOTAL
CHARACTER          UPLO, TRANS, OPP
INTEGER            DESCA( 9 ), DESCT( 9 ), DESCU( 9 ), DESCW( 9 ), &
                   ISEED( 4 )
DOUBLE PRECISION   ALPHA, EPS, ERR, TOL, AMAX
INTRINSIC          DBLE, MAX
EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_EXIT, &
                   BLACS_GRIDINIT, BLACS_GRIDINFO, BLACS_GRIDEXIT, &
                   NOMROC, DLAMCH, PDGEMM, PDLACPY, &
                   PDORMTR, PDELGET, PDELSET, PDSSTD2, &
                   PDLANGE, PDGENSSMAT, PDERRMAT
INTEGER            NUMROC
DOUBLE PRECISION   DLAMCH, PDLANGE, PDERRMAT
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
   CALL PDGENSSMAT( N, A, 1, 1, DESCA, ISEED, AMAX )
   IF ( K .LT. 1 ) THEN
      KK = MOD( N, NB )
      IF ( KK .EQ. 0 ) KK = NB
      DO J = 1, N-KK
         DO I = J+2, N
            CALL PDELSET( A, I, J, DESCA, ZERO )
            CALL PDELSET( A, J, I, DESCA, ZERO )
         END DO
      END DO
   ELSE
      DO J = NB+1, N
         DO I = 1, J-2
            CALL PDELSET( A, I, J, DESCA, ZERO )
            CALL PDELSET( A, J, I, DESCA, ZERO )
         END DO
      END DO
   END IF
   CALL PDLACPY( 'A', N, N, A, 1, 1, DESCA, AC, 1, 1, DESCA )
!
   TAU( : ) = ZERO
   LWORK = 2*MATLEN
   IF ( K .LT. 1 ) THEN
      CALL PDSSTD2( UPLO, KK, A, N-KK+1, N-KK+1, DESCA, E, TAU, &
           WORK, LWORK, INFO )
   ELSE
      CALL PDSSTD2( UPLO, MIN( N, NB ), A, 1, 1, DESCA, E, TAU, &
           WORK, LWORK, INFO )
   END IF
   IF ( K .LT. 1 ) THEN
      DO I = 1, N-1
         CALL PDELGET( 'A', ' ', E2( I ), A, I+1, I, DESCA )
         E2( I ) = -E2( I )
      END DO
   ELSE
      DO I = 1, N-1
         CALL PDELGET( 'A', ' ', E2( I ), A, I, I+1, DESCA )
      END DO
   END IF
!
   ERR = PDERRMAT( OPP, N, N, A, 1, 1, DESCA, AC, 1, 1, DESCA )
   IF ( INFO .NE. 0 ) ERR = DLAMCH( 'O' )
   IF ( ERR .NE. ZERO ) ERR = DLAMCH( 'O' )
!
   LLDU = LLDA
   CALL DESCSET( DESCU, N, N, NB, NB, 0, 0, ICTXT, LLDU )
   CALL PDLASET( 'A', N, N, ZERO, ONE, U, 1, 1, DESCU )
   CALL PDORMTR( 'R', UPLO, 'N', N, N, A, 1, 1, DESCA, TAU, U, 1, &
        1, DESCU, WORK, LWORK, INFO )
   IF ( INFO .NE. 0 ) ERR = DLAMCH( 'O' )
   LLDW = LLDA
   CALL DESCSET( DESCW, N, N, NB, NB, 0, 0, ICTXT, LLDW )
   CALL PDLASET( 'A', N, N, ZERO, -ONE, W, 1, 1, DESCW )
   CALL PDGEMM( 'T', 'N', N, N, N, ONE, U, 1, 1, DESCU, U, 1, 1, &
        DESCU, ONE, W, 1, 1, DESCW )
   ERR = ERR + PDLANGE( 'F', N, N, W, 1, 1, DESCW, WORK )
!
   LLDT = LLDA
   CALL DESCSET( DESCT, N, N, NB, NB, 0, 0, ICTXT, LLDT )
   CALL PDLASET( 'A', N, N, ZERO, ZERO, T, 1, 1, DESCT )
   DO I = 1, N-1
      CALL PDELSET( T, I, I+1, DESCT, E2( I ) )
      CALL PDELSET( T, I+1, I, DESCT, -E2( I ) )
   END DO
   CALL PDGEMM( 'N', 'N', N, N, N, ONE, AC, 1, 1, DESCA, U, 1, 1, &
        DESCU, ZERO, W, 1, 1, DESCW )
   CALL PDGEMM( 'T', 'N', N, N, N, -ONE, U, 1, 1, DESCU, W, 1, 1, &
        DESCW, ONE, T, 1, 1, DESCT )
   ERR = ERR + PDLANGE( 'F', N, N, T, 1, 1, DESCT, WORK )
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
