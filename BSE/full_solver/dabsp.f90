PROGRAM DABSP
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
INTEGER            NPTS
PARAMETER          ( NPTS = 1024 )
INTEGER            ITMAX
PARAMETER          ( ITMAX = 100 )
DOUBLE PRECISION   ZERO, ONE
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
INTEGER            ICTXT, IAM, SYS_NPROCS, &
                   NPROW, NPCOL, MYROW, MYCOL, SOLVER, INFO
INTEGER            N, NB, LLDA, LLDX, AROWS, ACOLS, XROWS, XCOLS, &
                   K, LWORK, LIWORK
INTEGER            I, J, ITMP1(1), ITMP2
DOUBLE PRECISION   SIGMA, STEP, ANORM
DOUBLE PRECISION   DTMP1(1), DTMP2(1), T, T_IO, T_SOLVER
!     ..
!     .. Local Arrays ..
INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), DESCX( DLEN_ ), &
                   DESCD( DLEN_ )
DOUBLE PRECISION   OMEGA( NPTS ), EPS( NPTS )
!
INTEGER,          ALLOCATABLE :: IWORK( : )
DOUBLE PRECISION, ALLOCATABLE :: A( : ), B( : ), X( : ), D( : ), &
                   LAMBDA( : ), WORK( : ), WORK0( : )
!     ..
!     .. Intrinsic Functions ..
INTRINSIC          DBLE, IAND, INT, MAX, MIN, DSQRT
!     ..
!     .. External Functions ..
EXTERNAL           NUMROC, MPI_WTIME, PDLANGE
INTEGER            NUMROC
DOUBLE PRECISION   MPI_WTIME, PDLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_EXIT, &
                   BLACS_GRIDINIT, BLACS_GRIDINFO, BLACS_GRIDEXIT, &
                   DESCINIT, IGEBS2D, IGEBR2D, DGEBS2D, DGEBR2D, &
                   PDELSET, PDTRADD, PDBSABSP
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
SOLVER = BSE_FULLBSE + BSE_LANCZOS + BSE_GAUSSIAN &
     + BSE_QUADAVGGAUSS
!      SOLVER = BSE_TDA + BSE_DIRECT + BSE_GAUSSIAN + BSE_LAPACK_HEEVR
!
IF ( ICTXT .GE. 0 ) THEN
   INFO = 0
   T = MPI_WTIME()
!
!        Read the dimension of A and B.
!
   IF ( MYROW+MYCOL .EQ. 0 ) THEN
      OPEN( UNIT = 10, FILE = 'input_real' )
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
      K = MIN( N, ITMAX )
      AROWS = NUMROC( N, NB, MYROW, 0, NPROW )
      ACOLS = NUMROC( N, NB, MYCOL, 0, NPCOL )
      LLDA = MAX( AROWS, 1 )
      CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, ICTXT, LLDA, &
           INFO )
      CALL DESCINIT( DESCB, N, N, NB, NB, 0, 0, ICTXT, LLDA, &
           INFO )
      XROWS = NUMROC( 2*N, NB, MYROW, 0, NPROW )
      LLDX = MAX( XROWS, 1 )
      IF ( IAND( SOLVER, BSE_ALGORITHM ) .EQ. BSE_DIRECT ) THEN
         XCOLS = NUMROC( 2*N, NB, MYCOL, 0, NPCOL )
         CALL DESCINIT( DESCX, 2*N, 2*N, NB, NB, 0, 0, ICTXT, &
              LLDX, INFO )
      ELSE
         XCOLS = NUMROC( K + 1, NB, MYCOL, 0, NPCOL )
         CALL DESCINIT( DESCX, 2*N, K + 1, NB, NB, 0, 0, ICTXT, &
              LLDX, INFO )
      END IF
      CALL DESCINIT( DESCD, N, 1, NB, NB, 0, 0, ICTXT, LLDA, &
           INFO )
      CALL PDBSABSP( SOLVER, N, NPTS, ONE, OMEGA, EPS, &
           DTMP1, 1, 1, DESCA, DTMP1, 1, 1, DESCB, DTMP1, &
           DTMP1, 1, 1, DESCX, DTMP1, 1, 1, DESCD, ITMAX, &
           DTMP2, -1, ITMP1, -1, INFO )
      LWORK = INT( DTMP2(1) )
      LIWORK = ITMP1(1)
!
!           Allocate memory.
!
      ALLOCATE( A( LLDA*ACOLS ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate A!'
      ALLOCATE( B( LLDA*ACOLS ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate B!'
      ALLOCATE( LAMBDA( N ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate lambda!'
      ALLOCATE( X( LLDX*XCOLS ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate X!'
      ALLOCATE( D( LLDA ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate D!'
      ALLOCATE( WORK( LWORK ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate WORK!'
      ALLOCATE( IWORK( LIWORK ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate IWORK!'
      ALLOCATE( WORK0( N ), STAT = INFO )
      IF ( INFO .NE. 0 ) WRITE( *, * ) 'Fail to allocate WORK0!'
!
!           Read the entries of A and B.
!
      DO J = 1, N
         IF ( MYROW+MYCOL .EQ. 0 ) THEN
            DO I = 1, N
               READ( UNIT = 10, FMT = * ) WORK0( I )
            END DO
            IF ( NPROW*NPCOL .GT. 1 ) &
               CALL DGEBS2D( ICTXT, 'A', ' ', N, 1, WORK0, N )
         ELSE
            CALL DGEBR2D( ICTXT, 'A', ' ', N, 1, WORK0, N, 0, 0 )
         END IF
         DO I = 1, N
            CALL PDELSET( A, I, J, DESCA, WORK0( I ) )
         END DO
      END DO
      IF ( MYROW+MYCOL .EQ. 0 ) THEN
         READ( UNIT = 10, FMT = * ) ITMP1, ITMP2
      END IF
      DO J = 1, N
         IF ( MYROW+MYCOL .EQ. 0 ) THEN
            DO I = 1, N
               READ( UNIT = 10, FMT = * ) WORK0( I )
            END DO
            IF ( NPROW*NPCOL .GT. 1 ) &
               CALL DGEBS2D( ICTXT, 'A', ' ', N, 1, WORK0, N )
         ELSE
            CALL DGEBR2D( ICTXT, 'A', ' ', N, 1, WORK0, N, 0, 0 )
         END IF
         DO I = 1, N
            CALL PDELSET( B, I, J, DESCB, WORK0( I ) )
         END DO
      END DO
      CALL PDTRADD( 'U', 'T', N-1, N-1, ONE, A, 2, 1, DESCA, ZERO, &
           A, 1, 2, DESCA )
      CALL PDTRADD( 'U', 'T', N-1, N-1, ONE, B, 2, 1, DESCB, ZERO, &
           B, 1, 2, DESCB )
!
!           Read the first N entries of D.
!
      IF ( MYROW+MYCOL .EQ. 0 ) THEN
         READ( UNIT = 10, FMT = * ) ITMP1, ITMP2
         DO I = 1, N
            READ( UNIT = 10, FMT = * ) WORK0( I )
         END DO
         IF ( NPROW*NPCOL .GT. 1 ) &
            CALL DGEBS2D( ICTXT, 'A', ' ', N, 1, WORK0, N )
      ELSE
         CALL DGEBR2D( ICTXT, 'A', ' ', N, 1, WORK0, N, 0, 0 )
      END IF
      DO I = 1, N
         CALL PDELSET( D, I, 1, DESCD, WORK0( I ) )
      END DO
!
!           Read sigma.
!
      IF ( MYROW+MYCOL .EQ. 0 ) THEN
         READ( UNIT = 10, FMT = * ) ITMP1, ITMP2
         READ( UNIT = 10, FMT = * ) SIGMA
         CLOSE( UNIT = 10 )
         IF ( NPROW*NPCOL .GT. 1 ) &
            CALL DGEBS2D( ICTXT, 'A', ' ', 1, 1, SIGMA, 1 )
      ELSE
         CALL DGEBR2D( ICTXT, 'A', ' ', 1, 1, SIGMA, 1, 0, 0 )
      END IF
      T_IO = MPI_WTIME() - T
!
!           Compute the optical absorption spectra.
!
      ANORM = PDLANGE( 'F', N, N, A, 1, 1, DESCA, WORK0 )
      STEP = ( 2*ANORM/DSQRT( DBLE( N ) ) )/( NPTS - 1 )
      DO J = 1, NPTS
         OMEGA( J ) = ( J - 1 )*STEP
      END DO
      T = MPI_WTIME()
      CALL PDBSABSP( SOLVER, N, NPTS, SIGMA, OMEGA, EPS, &
           A, 1, 1, DESCA, B, 1, 1, DESCB, LAMBDA, X, 1, 1, DESCX, &
           D, 1, 1, DESCD, K, WORK, LWORK, IWORK, LIWORK, INFO )
      T_SOLVER = MPI_WTIME() - T
!
      IF ( MYROW+MYCOL .EQ. 0 ) THEN
         WRITE( *, * ), 'clear all; close all;'
         DO J = 1, NPTS
            WRITE( *, * ), 'omega(', J, ', 1 ) =', OMEGA( J ), ';'
         END DO
         DO J = 1, NPTS
            WRITE( *, * ), 'eps(', J, ', 1 ) =', EPS( J ), ';'
         END DO
         IF ( IAND( SOLVER, BSE_ALGORITHM ) .EQ. BSE_DIRECT ) THEN
            K = N
         ELSE IF ( INFO .NE. 0 ) THEN
            K = INFO
         END IF
         DO J = 1, K
            WRITE( *, * ) 'lambda(', J, ', 1) =', LAMBDA( J ), ';'
         END DO
         WRITE( *, * ), 't_io =', T_IO, ';'
         WRITE( *, * ), 't_solver =', T_SOLVER, ';'
         WRITE( *, * ), 'info =', INFO, ';'
         WRITE( *, * ), 'sigma =', SIGMA, ';'
         WRITE( *, * ), 'plot(omega, eps);'
      END IF
!
!           Release allocated memory.
!
      DEALLOCATE( A, B, LAMBDA, X, D, WORK, IWORK, WORK0 )
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
