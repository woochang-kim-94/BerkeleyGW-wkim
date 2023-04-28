SUBROUTINE PZSPLANCZOS( N, K, A, IA, JA, DESCA, B, IB, JB, DESCB, &
                        X, IX, JX, DESCX, HA, HB, W, DESCW, INFO )
!
IMPLICIT NONE
!
!     .. Scalar Arguments ..
INTEGER            N, K, IA, JA, IB, JB, IX, JX, INFO
!     ..
!     .. Array Arguments ..
COMPLEX*16         A( * ), B( * ), X( * ), HB( 2, * ), W( * )
DOUBLE PRECISION   HA( 2, * )
INTEGER            DESCA( * ), DESCB( * ), DESCX( * ), DESCW( * )
!
!     .. Parameters ..
INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, &
                   LLD_, MB_, M_, NB_, N_, RSRC_
PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
DOUBLE PRECISION   ZERO, ONE
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
COMPLEX*16         CPLX_ZERO, CPLX_ONE
PARAMETER          ( CPLX_ZERO = ( 0.0D+0, 0.0D+0 ), &
                     CPLX_ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
INTEGER            ICTXT, J
DOUBLE PRECISION   DTMP, TMAX, TOL
!     ..
!     .. Local Arrays ..
COMPLEX*16         ZTMP( 4 )
complex*16 work0(2*N)
!     ..
!     .. Intrinsic Functions ..
INTRINSIC          DABS, DBLE, DCMPLX, DCONJG, DIMAG, DSQRT
!     ..
!     .. External Functions ..
EXTERNAL           PDLAMCH
DOUBLE PRECISION   PDLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL           PZAXPY, PZDSCAL, PZDOTC_ALL, PZDOTU_ALL, &
                   PZGEMM, PZHEMM, PZSYMM, PZLASET, PZLACONJ
!
INFO = 0
ICTXT = DESCX( CTXT_ )
TOL = 1.0D+4*N*PDLAMCH( ICTXT, 'E' )
TMAX = ZERO
!
CALL PZLASET( 'A', N, 1, ZERO, ZERO, X, IX+N, JX, DESCX )
CALL PZDOTC_ALL( N, ZTMP(1), X, IX, JX, DESCX, 1, &
     X, IX, JX, DESCX, 1 )
DTMP = DSQRT( ONE/DBLE( ZTMP( 1 ) ) )
CALL PZDSCAL( 2*N, DTMP, X, IX, JX, DESCX, 1 )
!      call pzlaprnt(n, 1, x, ix, jx, descx, 0,0,'x', 6, work0)
!      call pzlaprnt(n, 1, x, ix+n, jx, descx, 0,0,'y', 6, work0)
!
!     Lanczos process
!
DO J = 1, K
!
!        Matrix vector multiplications
!
!      call pzlaprnt(n, 1, x, ix, jx+j-1, descx, 0,0,'x1', 6, work0)
!      call pzlaprnt(n, 1, x, ix+n, jx+j-1, descx, 0,0,'y1', 6, work0)
   CALL PZHEMM( 'L', 'L', N, 1, CPLX_ONE, A, IA, JA, DESCA, &
        X, IX, JX+J-1, DESCX, CPLX_ZERO, W, 1, 1, DESCW )
!      call pzlaprnt(n, 1, w, 1, 1, descw, 0,0,'wx1', 6, work0)
   CALL PZSYMM( 'L', 'L', N, 1, CPLX_ONE, B, IB, JB, DESCB, &
        X, IX+N, JX+J-1, DESCX, CPLX_ONE, W, 1, 1, DESCW )
!      call pzlaprnt(n, 1, w, 1, 1, descw, 0,0,'wx2', 6, work0)
   CALL PZLACONJ( 2*N, 1, X, IX, JX+J-1, DESCX )
   CALL PZHEMM( 'L', 'L', N, 1, -CPLX_ONE, A, IA, JA, DESCA, &
        X, IX+N, JX+J-1, DESCX, CPLX_ZERO, W, N+1, 1, DESCW )
!      call pzlaprnt(n, 1, w, n+1, 1, descw, 0,0,'wy1', 6, work0)
   CALL PZSYMM( 'L', 'L', N, 1, -CPLX_ONE, B, IB, JB, DESCB, &
        X, IX, JX+J-1, DESCX, CPLX_ONE, W, N+1, 1, DESCW )
!      call pzlaprnt(n, 1, w, n+1, 1, descw, 0,0,'wy2', 6, work0)
   CALL PZLACONJ( 2*N, 1, X, IX, JX+J-1, DESCX )
   CALL PZLACONJ( N, 1, W, N+1, 1, DESCW )
!      call pzlaprnt(n, 1, w, 1, 1, descw, 0,0,'wx', 6, work0)
!      call pzlaprnt(n, 1, w, n+1, 1, descw, 0,0,'wy', 6, work0)
!
!        Gram--Schmidt.
!
   IF ( J .GT. 1 ) THEN
!            CALL PZGEADD( 'N', N, 1, DCMPLX( -HA( 2, J-1 ), ZERO ),
!     $           X, IX, JX+J-2, DESCX, CPLX_ONE, W, 1, 1, DESCW )
!            CALL PZGEADD( 'N', N, 1, DCMPLX( -HA( 2, J-1 ), ZERO ),
!     $           X, IX+N, JX+J-2, DESCX, CPLX_ONE, W, N+1, 1, DESCW )
      CALL PZAXPY( N, DCMPLX( -HA( 2, J-1 ), ZERO ), &
           X, IX, JX+J-2, DESCX, 1, W, 1, 1, DESCW, 1 )
      CALL PZAXPY( N, DCMPLX( -HA( 2, J-1 ), ZERO ), &
           X, IX+N, JX+J-2, DESCX, 1, W, N+1, 1, DESCW, 1 )
      CALL PZLACONJ( 2*N, 1, W, 1, 1, DESCW )
!            CALL PZGEADD( 'N', N, 1, HB( 2, J-1 ),
!     $           X, IX+N, JX+J-2, DESCX, CPLX_ONE, W, 1, 1, DESCW )
!            CALL PZGEADD( 'N', N, 1, HB( 2, J-1 ),
!     $           X, IX, JX+J-2, DESCX, CPLX_ONE, W, N+1, 1, DESCW )
      CALL PZAXPY( N, HB( 2, J-1 ), &
           X, IX+N, JX+J-2, DESCX, 1, W, 1, 1, DESCW, 1 )
      CALL PZAXPY( N, HB( 2, J-1 ), &
           X, IX, JX+J-2, DESCX, 1, W, N+1, 1, DESCW, 1 )
      CALL PZLACONJ( 2*N, 1, W, 1, 1, DESCW )
   END IF
!      call pzlaprnt(n, 1, w, 1, 1, descw, 0,0,'vx1', 6, work0)
!      call pzlaprnt(n, 1, w, n+1, 1, descw, 0,0,'vy1', 6, work0)
   CALL PZDOTC_ALL( N, ZTMP( 1 ), X, IX, JX+J-1, DESCX, 1, &
        W, 1, 1, DESCW, 1 )
   CALL PZDOTC_ALL( N, ZTMP( 2 ), X, IX+N, JX+J-1, DESCX, 1, &
        W, N+1, 1, DESCW, 1 )
   CALL PZDOTU_ALL( N, ZTMP( 3 ), X, IX, JX+J-1, DESCX, 1, &
        W, N+1, 1, DESCW, 1 )
   CALL PZDOTU_ALL( N, ZTMP( 4 ), X, IX+N, JX+J-1, DESCX, 1, &
        W, 1, 1, DESCW, 1 )
   HA( 1, J ) = DBLE( ZTMP( 1 ) ) - DBLE( ZTMP( 2 ) )
   HB( 1, J ) = DCONJG( ZTMP( 4 ) - ZTMP( 3 ) )
   IF ( J .LT. K )  THEN
!            CALL PZGEADD( 'N', N, 1, DCMPLX( -HA( 1, J ), ZERO ),
!     $           X, IX, JX+J-1, DESCX, CPLX_ONE, W, 1, 1, DESCW )
!            CALL PZGEADD( 'N', N, 1, DCMPLX( -HA( 1, J ), ZERO ),
!     $           X, IX+N, JX+J-1, DESCX, CPLX_ONE, W, N+1, 1, DESCW )
      CALL PZAXPY( N, DCMPLX( -HA( 1, J ), ZERO ), &
           X, IX, JX+J-1, DESCX, 1, W, 1, 1, DESCW, 1 )
      CALL PZAXPY( N, DCMPLX( -HA( 1, J ), ZERO ), &
           X, IX+N, JX+J-1, DESCX, 1, W, N+1, 1, DESCW, 1 )
      CALL PZLACONJ( 2*N, 1, W, 1, 1, DESCW )
!            CALL PZGEADD( 'N', N, 1, HB( 1, J ),
!     $           X, IX+N, JX+J-1, DESCX, CPLX_ONE, W, 1, 1, DESCW )
!            CALL PZGEADD( 'N', N, 1, HB( 1, J ),
!     $           X, IX, JX+J-1, DESCX, CPLX_ONE, W, N+1, 1, DESCW )
      CALL PZAXPY( N, HB( 1, J ), &
           X, IX+N, JX+J-1, DESCX, 1, W, 1, 1, DESCW, 1 )
      CALL PZAXPY( N, HB( 1, J ), &
           X, IX, JX+J-1, DESCX, 1, W, N+1, 1, DESCW, 1 )
      CALL PZLACONJ( 2*N, 1, W, 1, 1, DESCW )
!      call pzlaprnt(n, 1, w, 1, 1, descw, 0,0,'vx2', 6, work0)
!      call pzlaprnt(n, 1, w, n+1, 1, descw, 0,0,'vy2', 6, work0)
!
!           Normalization.
!
      CALL PZDOTC_ALL( N, ZTMP( 1 ), W, 1, 1, DESCW, 1, &
           W, 1, 1, DESCW, 1 )
      CALL PZDOTC_ALL( N, ZTMP( 2 ), W, N+1, 1, DESCW, 1, &
           W, N+1, 1, DESCW, 1 )
      DTMP = DSQRT( DBLE( ZTMP( 1 ) )**2 + DBLE( ZTMP( 2 ) )**2 )
      IF ( DTMP .LE. TOL*TMAX ) THEN
         WRITE( *, * ) '% Lucky breakdown at step', J
         INFO = J
         RETURN
      END IF
      TMAX = MAX( TMAX, DTMP )
      DTMP = DBLE( ZTMP( 1 ) ) - DBLE( ZTMP( 2 ) )
      IF ( DSQRT( DABS( DTMP ) ) .LE. TOL*TMAX ) THEN
         WRITE( *, * ) '% Neutral vectors encountered at step', J
         INFO = J
         RETURN
      ELSE IF ( DTMP .GT. ZERO ) THEN
         DTMP = DSQRT( DTMP )
         CALL PZDSCAL( 2*N, ONE/DTMP, W, 1, 1, DESCW, 1 )
         CALL PZCOPY( 2*N, W, 1, 1, DESCW, 1, X, IX, JX+J, DESCX, &
              1 )
         HA( 2, J ) = DTMP
         HB( 2, J ) = DCMPLX( ZERO, ZERO )
      ELSE IF ( DTMP .LT. ZERO ) THEN
         DTMP = DSQRT( -DTMP )
         CALL PZDSCAL( 2*N, ONE/DTMP, W, 1, 1, DESCW, 1 )
         CALL PZLACONJ( 2*N, 1, W, 1, 1, DESCW )
         CALL PZCOPY( N, W, N+1, 1, DESCW, 1, X, IX, JX+J, DESCX, &
              1 )
         CALL PZCOPY( N, W, 1, 1, DESCW, 1, X, IX+N, JX+J, DESCX, &
              1 )
         HA( 2, J ) = ZERO
         HB( 2, J ) = DCMPLX( -DTMP, ZERO )
      END IF
!      call pzlaprnt(n, 1, x, ix, jx+j, descx, 0,0,'x', 6, work0)
!      call pzlaprnt(n, 1, x, ix+n, jx+j, descx, 0,0,'y', 6, work0)
   END IF
END DO
!
RETURN
!
!     End of PZSPLANCZOS().
!
!
END
