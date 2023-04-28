SUBROUTINE ZBSABSP2( K, NPTS, SIGMA, OMEGA, EPS2, NORMD2, X, LDX, &
                     MU )
!
IMPLICIT NONE
INCLUDE 'solver.f90'
!
!     .. Scalar Arguments ..
INTEGER            K, NPTS, LDX
DOUBLE PRECISION   SIGMA, NORMD2
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION   OMEGA( * ), EPS2( * ), MU( * )
COMPLEX*16         X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  ZBSABSP2() estimates the optical absorption spectra
!
!     epsilon_2( omega ) =
!        ||d||**2 * sum_i |x_i(1) - x_i(k+1)|^2 * delta( omega - mu_i )
!
!  using the eigenpairs (mu_j,x_j) of the 2k-by-2k projected
!  Bethe--Salpeter Hamiltonian
!
!     H_proj = [       HA,        HB;
!                -conj(HB), -conj(HA) ].
!
!  The summation is taken over all positive eigenvalues/eigenvectors.
!
!  The delta function is replaced by the Gaussian function
!
!     g( x ) = exp( -x**2/( 2*sigma**2 ) )/( sqrt( 2*pi )*sigma ).
!
!  The eigenvectors x_i = [ x_i1**H, x_i2**H ]**H are normalized such
!  that x_i1**H * x_i1 - x_i2**H * x_i2 = 1.
!
!  Arguments
!  =========
!
!  K       (input) INTEGER
!          2*N is the number of rows and columns of the projected
!          Hamiltonian.
!          K >= 0.
!
!  NPTS    (input) INTEGER
!          NPTS is the number of sampling points in OMEGA.
!          NPTS >= 0.
!
!  SIGMA   (input) DOUBLE PRECISION
!          Standard deviation of the Gaussian function.
!          SIGMA > 0.
!
!  OMEGA   (input) DOUBLE PRECISION array, dimension (NPTS)
!          Sampling points of omega.
!          All entries of OMEGA are nonnegative.
!
!  EPS2    (output) DOUBLE PRECISION array, dimension (NPTS)
!          Sampling points of epsilon_2, i.e.,
!          EPS2( I ) = epsilon_2( OMEGA( I ) ).
!
!  NORMD2  (input) DOUBLE PRECISION
!          The squared norm of the first block, d_1, of the dipole
!          vector d, i.e., ||d_1||**2.
!
!  X       (input) COMPLEX*16 pointer into the local memory to an
!          array of dimension (LLD_X, LOCc(JX+N-1)).
!          The normalized right eigenvectors of the projected
!          Hamiltonian corresponding to the positive eigenvalues.
!
!  LDX     (input) INGEGER
!          The leading dimension of the array X. LDX >= max(1,2*K).
!
!  MU      (input) DOUBLE PRECISION array, dimension (K)
!          The positive eigenvalues of the projected Hamiltonian.
!
!  =====================================================================
!
!     .. Parameters ..
DOUBLE PRECISION   ZERO
PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
INTEGER            I, J
DOUBLE PRECISION   DTMP
COMPLEX*16         ZTMP
!     ..
!     .. Intrinsic Functions ..
INTRINSIC          DBLE, DCONJG
!     ..
!     .. External Functions ..
EXTERNAL           APPROX_DELTA
DOUBLE PRECISION   APPROX_DELTA
!     ..
!     .. Executable Statements ..
!
DO J = 1, NPTS
   EPS2( J ) = ZERO
END DO
DO I = 1, K
   ZTMP = X( 1, I ) - X( K+1, I )
   DTMP = NORMD2*DBLE( ZTMP*DCONJG( ZTMP ) )
   DO J = 1, NPTS
      EPS2( J ) = EPS2( J ) + DTMP &
           *APPROX_DELTA( BSE_GAUSSIAN, OMEGA( J ) - MU( I ), SIGMA )
   END DO
END DO
!
RETURN
!
!     End of ZBSABSP2().
!
END
