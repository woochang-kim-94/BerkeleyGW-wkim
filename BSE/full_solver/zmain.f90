PROGRAM ZMAIN
!
IMPLICIT NONE
INCLUDE 'solver.f90'
!
!     .. Parameters ..
INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, &
                   LLD_, MB_, M_, NB_, N_, RSRC_
PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
COMPLEX*16         ZERO, ONE
PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
                     ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
INTEGER            ICTXT, IAM, SYS_NPROCS, &
                   NPROW, NPCOL, MYROW, MYCOL, SOLVER, INFO
INTEGER            N, NB, LLDA, LLDB, LLDX, LLDY, AROWS, ACOLS, &
                   XROWS, XCOLS, LWORK, LRWORK, LIWORK, &
                   INDW1, INDW2, INDW3
INTEGER            I, J, ITMP1(1), ITMP2(1), ITMP3(1), J0
DOUBLE PRECISION   DTMP1(1), DTMP2(1), ORTH, RES, T, T_IO, T_SOLVER, &
                   T_CHECK, ANORM, BNORM, HNORM, XNORM, &
                   ORTH1, ORTH2, RES1, RES2
COMPLEX*16         ZTMP1(1), ZTMP2(1)
!     ..
!     .. Local Arrays ..
INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), DESCX( DLEN_ ), &
                   DESCY( DLEN_ )
!
INTEGER         , ALLOCATABLE :: IWORK( : )
DOUBLE PRECISION, ALLOCATABLE :: LAMBDA( : ), RWORK( : )
COMPLEX*16,       ALLOCATABLE :: A( : ), B( : ), X( : ), Y( : ), &
                   WORK( : ), WORK0( : ), AA( : ), BB( : )
!     ..
!     .. Intrinsic Functions ..
INTRINSIC          DBLE, DCMPLX, DCONJG, IAND, INT, MAX, DSQRT
!     ..
!     .. External Functions ..
EXTERNAL           INDXL2G, NUMROC, PZLANGE, MPI_WTIME
INTEGER            INDXL2G, NUMROC
DOUBLE PRECISION   PZLANGE, MPI_WTIME
!     ..
!     .. External Subroutines ..
EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_EXIT, &
                   BLACS_GRIDINIT, BLACS_GRIDINFO, BLACS_GRIDEXIT, &
                   IGEBS2D, IGEBR2D, ZGEBS2D, ZGEBR2D, PZTRADD, &
                   PZELSET, PZBSEIG, PZEMBED2, PZGEMM, PZLACPY, &
                   PZLASET, PZGEADD, PZSCAL
EXTERNAL           PZLAPRNT
!     ..
!     .. Executable Statements ..
!
!     BLACS initialization.
!
CALL BLACS_PINFO( IAM, SYS_NPROCS )
NPROW = INT( DSQRT( DBLE( SYS_NPROCS ) ) )
NPCOL = SYS_NPROCS / NPROW
CALL BLACS_GET( 0, 0, ICTXT )
CALL BLACS_GRIDINIT( ICTXT, 'R', NPROW, NPCOL )
CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
NB = 64
SOLVER = BSE_FULLBSE + BSE_DIRECT
!      SOLVER = BSE_TDA + BSE_DIRECT + BSE_LAPACK_HEEVR
!
IF ( ICTXT .GE. 0 ) THEN
   INFO = 0
   T = MPI_WTIME()
!
!        Read the dimension of A and B.
!
   IF ( MYROW+MYCOL .EQ. 0 ) THEN
      OPEN( UNIT = 10, FILE = 'input' )
      READ( UNIT = 10, FMT = * ) N, ITMP1
      IF ( NPROW*NPCOL .GT. 1 ) &
         CALL IGEBS2D( ICTXT, 'A', ' ', 1, 1, N, 1 )
   ELSE
      CALL IGEBR2D( ICTXT, 'A', ' ', 1, 1, N, 1, 0, 0 )
   END IF
!
   IF ( N .GT. 0 ) THEN
!
!           Set up description vectors and local indices.
!
      AROWS = NUMROC( N, NB, MYROW, 0, NPROW )
      ACOLS = NUMROC( N, NB, MYCOL, 0, NPCOL )
      LLDA = MAX( AROWS, 1 )
      LLDB = LLDA
      CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, ICTXT, LLDA, &
           INFO )
      CALL DESCINIT( DESCB, N, N, NB, NB, 0, 0, ICTXT, LLDB, &
           INFO )
      XROWS = NUMROC( 2*N, NB, MYROW, 0, NPROW )
      XCOLS = NUMROC( 2*N, NB, MYCOL, 0, NPCOL )
      LLDX = MAX( XROWS, 1 )
      LLDY = LLDX
      CALL DESCINIT( DESCX, 2*N, 2*N, NB, NB, 0, 0, ICTXT, LLDX, &
           INFO )
      CALL DESCINIT( DESCY, 2*N, 2*N, NB, NB, 0, 0, ICTXT, LLDY, &
           INFO )
      CALL PZBSEIG( SOLVER, N, ZTMP1, 1, 1, DESCA, ZTMP1, 1, 1, &
           DESCB, DTMP1, ZTMP1, 1, 1, DESCX, ZTMP2, -1, DTMP2, -1, &
           ITMP3, -1, INFO )
      LWORK = INT( ZTMP2(1) )
      LRWORK = INT( DTMP2(1) )
      LIWORK = ITMP3(1)
      INDW1 = 1
      INDW2 = INDW1 + LLDX*XCOLS
      INDW3 = INDW2 + LLDX*XCOLS
      LWORK = MAX( LWORK, 3*LLDX*XCOLS )
!
!           Allocate memory.
!
      ALLOCATE( A( LLDA*ACOLS ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate A!'
      ALLOCATE( B( LLDA*ACOLS ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate B!'
      ALLOCATE( LAMBDA( 2*N ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate lambda!'
      ALLOCATE( X( LLDX*XCOLS ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate X!'
      ALLOCATE( Y( LLDX*XCOLS ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate Y!'
      ALLOCATE( WORK0( 4*N ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate WORK0!'
      ALLOCATE( WORK( LWORK ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate WORK!'
      ALLOCATE( RWORK( LRWORK ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate RWORK!'
      ALLOCATE( IWORK( LIWORK ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate IWORK!'
      ALLOCATE( AA( LLDA*ACOLS ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate AA!'
      ALLOCATE( BB( LLDA*ACOLS ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate BB!'
!
!           Read the entries of A and B.
!
      DO J = 1, N
         IF ( MYROW+MYCOL .EQ. 0 ) THEN
            DO I = 1, N
               READ( UNIT = 10, FMT = * ) DTMP1(1), DTMP2(1)
               WORK0( I ) = DCMPLX( DTMP1(1), DTMP2(1) )
            END DO
            IF ( NPROW*NPCOL .GT. 1 ) &
               CALL ZGEBS2D( ICTXT, 'A', ' ', N, 1, WORK0, N )
         ELSE
            CALL ZGEBR2D( ICTXT, 'A', ' ', N, 1, WORK0, N, 0, 0 )
         END IF
         DO I = 1, N
            CALL PZELSET( A, I, J, DESCA, WORK0( I ) )
         END DO
         CALL PZELSET( A, J, J, DESCA, &
              DCMPLX( DBLE( WORK0( J ) ) ) )
      END DO
      IF ( MYROW+MYCOL .EQ. 0 ) THEN
         READ( UNIT = 10, FMT = * ) ITMP1, ITMP2
      END IF
      DO J = 1, N
         IF ( MYROW+MYCOL .EQ. 0 ) THEN
            DO I = 1, N
               READ( UNIT = 10, FMT = * ) DTMP1(1), DTMP2(1)
               WORK0( I ) = DCMPLX( DTMP1(1), DTMP2(1) )
            END DO
            IF ( NPROW*NPCOL .GT. 1 ) &
               CALL ZGEBS2D( ICTXT, 'A', ' ', N, 1, WORK0, N )
         ELSE
            CALL ZGEBR2D( ICTXT, 'A', ' ', N, 1, WORK0, N, 0, 0 )
         END IF
         DO I = 1, N
            CALL PZELSET( B, I, J, DESCB, WORK0( I ) )
         END DO
      END DO
      IF ( MYROW+MYCOL .EQ. 0 ) THEN
         CLOSE( UNIT = 10 )
      END IF
      CALL PZTRADD( 'U', 'C', N-1, N-1, ONE, A, 2, 1, DESCA, ZERO, &
           A, 1, 2, DESCA )
      CALL PZTRADD( 'U', 'T', N-1, N-1, ONE, B, 2, 1, DESCB, ZERO, &
           B, 1, 2, DESCB )
      CALL PZLACPY( 'A', N, N, A, 1, 1, DESCA, AA, 1, 1, DESCA )
      CALL PZLACPY( 'A', N, N, B, 1, 1, DESCB, BB, 1, 1, DESCB )
!
!            IF ( MYROW+MYCOL .EQ. 0 ) WRITE( *, * ) 'clear all;'
!            CALL PZLAPRNT( N, N, A, 1, 1, DESCA, 0, 0, 'A', 6, WORK0 )
!            CALL PZLAPRNT( N, N, B, 1, 1, DESCB, 0, 0, 'B', 6, WORK0 )
      T_IO = MPI_WTIME() - T
!
      CALL PZLASET( 'A', 2*N, N, ZERO, ZERO, X, 1, 1, DESCX )
      CALL PZBSEIG( SOLVER, N, A, 1, 1, DESCA, B, 1, 1, DESCB, &
           LAMBDA, X, 1, 1, DESCX, WORK, -1, RWORK, -1, IWORK, -1, &
           INFO )
      CALL IGAMX2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, ITMP1, ITMP2, &
           -1, -1, -1 )
      IF ( INFO .EQ. 0 ) THEN
         T = MPI_WTIME()
         CALL PZBSEIG( SOLVER, N, A, 1, 1, DESCA, B, 1, 1, DESCB, &
              LAMBDA, X, 1, 1, DESCX, WORK, LWORK, RWORK, LRWORK, &
              IWORK, LIWORK, INFO )
         T_SOLVER = MPI_WTIME() - T
         T = MPI_WTIME()
         IF ( MYROW+MYCOL .EQ. 0 ) THEN
            DO I = 1, N
               WRITE( *, * ) 'lambda(', I, ', 1) =', &
                    LAMBDA( I ), ';'
            END DO
         END IF
!               CALL PZLAPRNT( 2*N, N, X, 1, 1, DESCX, 0, 0, 'X', 6,
!     $              WORK0 )
         T_IO = MPI_WTIME() - T + T_IO
!
!              Check the accuracy.
!
         T = MPI_WTIME()
         CALL PZLACPY( 'A', N, N, AA, 1, 1, DESCA, A, 1, 1, &
              DESCA )
         IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_TDA ) THEN
            CALL PZLASET( 'A', N, N, ZERO, ZERO, B, 1, 1, DESCB )
         ELSE
            CALL PZLACPY( 'A', N, N, BB, 1, 1, DESCB, B, 1, 1, &
                 DESCB )
         END IF
         ANORM = PZLANGE( 'F', N, N, A, 1, 1, DESCA, WORK0 )
         BNORM = PZLANGE( 'F', N, N, B, 1, 1, DESCB, WORK0 )
         XNORM = DSQRT( 2.0D+0 ) &
              *PZLANGE( 'F', 2*N, N, X, 1, 1, DESCX, WORK0 )
         HNORM = DSQRT( ANORM**2 + BNORM**2 )
!
         CALL PZGEADD( 'N', N, N, ONE, X, 1, 1, DESCX, ZERO, &
              WORK( INDW1 ), 1, 1, DESCX )
         DO I = 1, N
            CALL PZSCAL( N, DCMPLX( -LAMBDA( I ) ), WORK( INDW1 ), &
                 1, I, DESCX, 1 )
         END DO
         CALL PZGEMM( 'N', 'N', N, N, N, ONE, A, 1, 1, DESCA, &
              X, 1, 1, DESCX, ONE, WORK( INDW1 ), 1, 1, DESCX )
         CALL PZGEMM( 'N', 'N', N, N, N, ONE, B, 1, 1, DESCB, &
              X, N+1, 1, DESCX, ONE, WORK( INDW1 ), 1, 1, DESCX )
         RES1 = PZLANGE( 'F', N, N, WORK( INDW1 ), 1, 1, DESCX, &
              WORK0 )
         CALL PZGEADD( 'N', N, N, ONE, X, N+1, 1, DESCX, ZERO, &
              WORK( INDW1 ), 1, 1, DESCX )
         DO I = 1, N
            CALL PZSCAL( N, DCMPLX( LAMBDA( I ) ), WORK( INDW1 ), &
                 1, I, DESCX, 1 )
         END DO
         CALL PZGEMM( 'T', 'N', N, N, N, ONE, A, 1, 1, DESCA, &
              X, N+1, 1, DESCX, ONE, WORK( INDW1 ), 1, 1, DESCX )
         CALL PZGEMM( 'C', 'N', N, N, N, ONE, B, 1, 1, DESCB, &
              X, 1, 1, DESCX, ONE, WORK( INDW1 ), 1, 1, DESCX )
         RES2 = PZLANGE( 'F', N, N, WORK( INDW1 ), 1, 1, DESCX, &
              WORK0 )
         RES = DSQRT( 2.0D+0*( RES1**2 + RES2**2 ) )/HNORM/XNORM
         IF ( MYROW+MYCOL .EQ. 0 ) THEN
            WRITE( *, * ) &
                 '% ||A*X_1+BX_2-X_1*diag(lambda)||_F =', RES1
            WRITE( *, * ) &
                 '% ||A**T*X_2+B**H*X_1+X_2*diag(lambda)||_F =', &
                 RES2
            WRITE( *, * ) &
                 '% ||H*X-X*diag(lambda)||_F/(||H||_F||X||_F) =', &
                 RES
         END IF
!
         CALL PZLASET( 'A', N, N, ZERO, -ONE, WORK( INDW1 ), 1, 1, &
              DESCX )
         CALL PZGEMM( 'C', 'N', N, N, N, ONE, X, 1, 1, DESCX, &
              X, 1, 1, DESCX, ONE, WORK( INDW1 ), 1, 1, DESCX )
         CALL PZGEMM( 'C', 'N', N, N, N, -ONE, X, N+1, 1, DESCX, &
              X, N+1, 1, DESCX, ONE, WORK( INDW1 ), 1, 1, DESCX )
         ORTH1 = PZLANGE( 'F', N, N, WORK( INDW1 ), 1, 1, DESCX, &
              WORK0 )
         CALL PZLASET( 'A', N, N, ZERO, ZERO, WORK( INDW1 ), 1, 1, &
              DESCX )
         CALL PZGEMM( 'T', 'N', N, N, N, ONE, X, 1, 1, DESCX, &
              X, N+1, 1, DESCX, ONE, WORK( INDW1 ), 1, 1, DESCX )
         CALL PZGEMM( 'T', 'N', N, N, N, ONE, X, N+1, 1, DESCX, &
              X, 1, 1, DESCX, -ONE, WORK( INDW1 ), 1, 1, DESCX )
         ORTH2 = PZLANGE( 'F', N, N, WORK( INDW1 ), 1, 1, DESCX, &
              WORK0 )
         ORTH = DSQRT( ( ORTH1**2 + ORTH2**2 )/2.0D+0 )/N
         IF ( MYROW+MYCOL .EQ. 0 ) THEN
            WRITE( *, * ) &
                 '% ||X_1**T*X_1-X_2**T*X_2-I||_F =', ORTH1
            WRITE( *, * ) &
                 '% ||X_2**T*X_1-X_1**T*X_2||_F =', ORTH2
            WRITE( *, * ) &
                 '% ||Y**T*X-I||_F/(2*N) =', ORTH
         END IF
!
         DO I = 1, N
            LAMBDA( N+I ) = -LAMBDA( I )
         END DO
         CALL PZGEADD( 'N', N, N, ONE, X, 1, 1, DESCX, &
              ZERO, X, 1+N, 1+N, DESCX )
         CALL PZGEADD( 'N', N, N, ONE, X, 1+N, 1, DESCX, &
              ZERO, X, 1, 1+N, DESCX )
         DO J = 1, XCOLS
            J0 = INDXL2G( J, NB, MYCOL, DESCX( CSRC_ ), NPCOL )
            IF ( J0 .GT. N ) THEN
               DO I = 1, XROWS
                  X( I + (J-1)*DESCX( LLD_ ) ) = &
                       DCONJG( X( I + (J-1)*DESCX( LLD_ ) ) )
               END DO
            END IF
         END DO
         CALL PZGEADD( 'N', N, N, ONE, X, 1, 1, DESCX, &
              ZERO, Y, 1, 1, DESCY )
         CALL PZGEADD( 'N', N, N, -ONE, X, 1+N, 1, DESCX, &
              ZERO, Y, 1+N, 1, DESCY )
         CALL PZGEADD( 'N', N, N, -ONE, X, 1, 1+N, DESCX, &
              ZERO, Y, 1, 1+N, DESCY )
         CALL PZGEADD( 'N', N, N, ONE, X, 1+N, 1+N, DESCX, &
              ZERO, Y, 1+N, 1+N, DESCY )
!
         CALL PZLASET( 'A', 2*N, 2*N, ZERO, -ONE, WORK( INDW1 ), &
              1, 1, DESCX )
         CALL PZGEMM( 'C', 'N', 2*N, 2*N, 2*N, ONE, X, 1, 1, &
              DESCX, Y, 1, 1, DESCY, ONE, WORK( INDW1 ), 1, 1, &
              DESCX )
         ORTH = PZLANGE( 'F', 2*N, 2*N, WORK( INDW1 ), 1, 1, &
              DESCX, WORK0 )
         IF ( MYROW+MYCOL .EQ. 0 ) &
            WRITE( *, * ) '% || Y**H * X - I ||_F =', ORTH
!
         CALL PZEMBED2( N, A, 1, 1, DESCA, B, 1, 1, DESCB, &
              WORK( INDW1 ), 1, 1, DESCX, WORK( INDW2 ) )
         CALL PZLASET( 'A', 2*N, 2*N, ZERO, ZERO, WORK( INDW3 ), &
              1, 1, DESCX )
         DO I = 1, 2*N
            CALL PZELSET( WORK( INDW3 ), I, I, DESCX, &
                 DCMPLX( -LAMBDA( I ) ) )
         END DO
         CALL PZGEMM( 'C', 'N', 2*N, 2*N, 2*N, ONE, Y, 1, 1, &
              DESCY, WORK( INDW1 ), 1, 1, DESCX, ZERO, &
              WORK( INDW2 ), 1, 1, DESCX )
         CALL PZGEMM( 'N', 'N', 2*N, 2*N, 2*N, ONE, WORK( INDW2 ), &
              1, 1, DESCX, X, 1, 1, DESCX, ONE, WORK( INDW3 ), 1, &
              1, DESCX )
         RES = PZLANGE( 'F', 2*N, 2*N, WORK( INDW3 ), 1, 1, &
              DESCX, WORK0 )
         IF ( MYROW+MYCOL .EQ. 0 ) &
            WRITE( *, * ) &
                 '% || Y**H * H * X - diag(lambda) ||_F =', RES
         T_CHECK = MPI_WTIME() - T
      END IF
      IF ( INFO.NE. 0 ) THEN
         WRITE( *, * ) '%', MYROW, MYCOL, 'INFO =', INFO
      END IF
      IF ( MYROW+MYCOL .EQ. 0 ) THEN
         WRITE( *, * ), 't_io =', T_IO, ';'
         WRITE( *, * ), 't_solver =', T_SOLVER, ';'
         WRITE( *, * ), 't_check =', T_CHECK, ';'
      END IF
!
!           Release allocated memory.
!
      DEALLOCATE( A, B, LAMBDA, X, Y, WORK, RWORK, IWORK, WORK0, &
           AA, BB )
!
   END IF
!
   CALL BLACS_GRIDEXIT( ICTXT )
!
END IF
!
CALL BLACS_EXIT( 0 )
!
END
