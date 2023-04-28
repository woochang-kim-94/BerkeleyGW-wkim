SUBROUTINE ZEMBED3( N, A, B, M, LDM )
!
IMPLICIT NONE
!
!     .. Scalar Arguments ..
INTEGER            N, LDA, LDB, LDM
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION   A( 2, * ), M( LDM, * )
COMPLEX*16         B( 2, * )
!     ..
!
!  Purpose
!  =======
!
!  ZEMBED3() is an auxiliary routine post-processing the output of the
!  Lanczos process.
!  It generates a 2n-by-2n real symmetric matrix
!
!     M = [ real(A+B), imag(A-B);
!          -imag(A+B), real(A-B) ],
!
!  where A is n-by-n real symmetric tridiagonal, and B is n-by-n complex
!  symmetric tridiagonal with real subdiagonal entries.
!  The lower triangular parts of A and B are provided as band matrices.
!  Only the lower triangular part of M is generated.
!
!  No argument check is performed, i.e.,
!  all arguments are assumed to be valid.
!
!  Arguments
!  =========
!
!  N       (input) INGEGER
!          2*N is the number of rows and columns of M.
!
!  A       (input) DOUBLE PRECISION array, dimension (2,N)
!          The n-by-n real symmetric tridiagonal matrix A.
!          The lower triangular part of A is stored in band format,
!          i.e., the first row contains the diagonals and the second row
!          contains the subdiagonals.
!
!  B       (input) COMPLEX*16 array, dimension (2,N)
!          The n-by-n complex symmetric tridiagonal matrix B.
!          The lower triangular part of B is stored in band format,
!          i.e., the first row contains the diagonals and the second row
!          contains the subdiagonals.
!
!  M       (output) DOUBLE PRECISION array, dimension (LDM,2*N)
!          The 2n-by-2n matrix real symmetric matrix M.
!          Only the lower triangular part of M is referenced.
!
!  LDM     (input) INGEGER
!          The leading dimension of the array M. LDM >= max(1,2*N).
!
!  =====================================================================
!
!     .. Parameters ..
DOUBLE PRECISION   ZERO
PARAMETER          ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
INTEGER            I
!     ..
!     .. Intrinsic Functions ..
INTRINSIC          DBLE, DIMAG
!     ..
!     .. External Subroutines ..
EXTERNAL           DLASET
!     ..
!     .. Executable Statements ..
!
CALL DLASET( 'L', 2*N, 2*N, ZERO, ZERO, M, LDM )
DO I = 1, N
   M( I, I ) = A( 1, I ) + DBLE( B( 1, I ) )
   M( N+I, I ) = -DIMAG( B( 1, I ) )
   M( N+I, N+I ) = A( 1, I ) - DBLE( B( 1, I ) )
END DO
DO I = 1, N-1
   M( I+1, I ) = A( 2, I ) + DBLE( B( 2, I ) )
   M( N+I+1, N+I ) = A( 2, I ) - DBLE( B( 2, I ) )
END DO
!
RETURN
!
!     End of ZEMBED3().
!
END
