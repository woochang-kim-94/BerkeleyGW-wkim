
!===============================================================================
!
! Routines:
!
! 1. real_wfng()        Originally By gsm        Last Modified 9/1/2010 (gsm)
!
!    Constructs real wavefunctions in G-space for systems
!    with inversion symmetry by applying Gram-Schmidt process.
!    Based on paratecSGL/src/para/gwreal.f90
!
!===============================================================================

#include "f_defs.h"

module real_m

  use global_m

  implicit none

  private

  public :: real_wfng

contains

!-------------------------------------------------------------------------------

SUBROUTINE real_wfng ( ngk_l, nb, ns, nk, energy, wfn_d )

  integer, intent (in) :: ngk_l, nb, ns, nk
  real(DP), intent (in) :: energy ( :, :, : ) !< ( nb, nk, ns )
  complex(DPC), intent (inout) :: wfn_d ( :, :, :, : ) !< ( ngk_l, nb, ns, nk )
  
  real(DP), PARAMETER :: eps2 = 1.0D-2
  real(DP), PARAMETER :: eps5 = 1.0D-5
  real(DP), PARAMETER :: eps6 = 1.0D-6
  
  integer :: i, j, k, ik, is, ib, jb, ig, inum, deg, mdeg, inc
  integer :: dimension_span, reduced_span
  real(DP) :: x, y
  integer, allocatable :: imap ( :, : )
  integer, allocatable :: inums ( : )
  integer, allocatable :: inull ( : )
  integer, allocatable :: null_map ( :, : )
  real(DP), allocatable :: psi ( :, : )
  real(DP), allocatable :: phi ( :, : )
  real(DP), allocatable :: vec ( : )
  complex(DPC), allocatable :: wfc ( : )
  
  PUSH_SUB(real_wfng)
  
  DO ik = 1, nk
    
    mdeg = 1
    DO is = 1, ns
      DO ib = 1, nb - 1
        deg = 1
        DO jb = ib + 1, nb
          IF ( abs ( energy ( ib, ik, is ) - &
            energy ( jb, ik, is ) ) .LT. eps5 * &
            dble ( jb - ib + 1 ) ) deg = deg + 1
        ENDDO
        IF ( deg .GT. mdeg ) mdeg = deg
      ENDDO
    ENDDO
    mdeg = mdeg * 2
    
    SAFE_ALLOCATE ( imap, ( nb, ns ) )
    SAFE_ALLOCATE ( inums, ( ns ) )
    SAFE_ALLOCATE ( inull, ( nb ) )
    SAFE_ALLOCATE ( null_map, ( mdeg, nb ) )
    
    DO is = 1, ns
      inum = 1
      DO ib = 1, nb
        IF ( ib .EQ. nb ) THEN
          imap ( inum, is ) = ib
          inum = inum + 1
        ELSEIF ( abs ( energy ( ib, ik, is ) - &
          energy ( ib + 1, ik, is ) ) .GT. eps5 ) THEN
          imap ( inum, is ) = ib
          inum = inum + 1
        ENDIF
      ENDDO
      inum = inum - 1
      inums ( is ) = inum
    ENDDO
    
    SAFE_ALLOCATE ( wfc, ( ngk_l ) )
    SAFE_ALLOCATE ( psi, ( ngk_l, mdeg ) )
    SAFE_ALLOCATE ( phi, ( ngk_l, mdeg ) )
    SAFE_ALLOCATE ( vec, ( ngk_l ) )
    
    DO is = 1, ns
      inc = 1
      inum = inums ( is )
      DO i = 1, inum
        inull ( i ) = 1
        DO ib = inc, imap ( i, is )
          DO ig = 1, ngk_l
            wfc ( ig ) = wfn_d ( ig, ib, is, ik )
          ENDDO
          x = 0.0D0
          DO ig = 1, ngk_l
            x = x + dble ( wfc ( ig ) ) **2
          ENDDO
#ifdef MPI
          y = x
          CALL MPI_Allreduce ( y, x, 1, MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
          IF ( x .LT. eps2 ) null_map ( inull ( i ), i ) = 0
          IF ( x .GT. eps2 ) null_map ( inull ( i ), i ) = 1
          inull ( i ) = inull ( i ) + 1
          x = 0.0D0
          DO ig = 1, ngk_l
            x = x + aimag ( wfc ( ig ) ) **2
          ENDDO
#ifdef MPI
          y = x
          CALL MPI_Allreduce ( y, x, 1, MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
          IF ( x .LT. eps2 ) null_map ( inull ( i ), i ) = 0
          IF ( x .GT. eps2 ) null_map ( inull ( i ), i ) = 1
          inull ( i ) = inull ( i ) + 1
        ENDDO
        inull ( i ) = inull ( i ) - 1
        inc = imap ( i, is ) + 1
      ENDDO
      inc = 1
      ib = 1
      DO i = 1, inum
        k = 1
        DO j = 1, 2 * ( imap ( i, is ) - inc ) + 1, 2
          IF ( null_map ( j, i ) .EQ. 1 .OR. &
            null_map ( j + 1, i ) .EQ. 1 ) THEN
            DO ig = 1, ngk_l
              wfc ( ig ) = wfn_d ( ig, ib, is, ik )
            ENDDO
            IF ( null_map ( j, i ) .EQ. 1 ) THEN
              DO ig = 1, ngk_l
                phi ( ig, k ) = dble ( wfc ( ig ) )
              ENDDO
              k = k + 1
            ENDIF
            IF ( null_map ( j + 1, i ) .EQ. 1 ) THEN
              DO ig = 1, ngk_l
                phi ( ig, k ) = aimag ( wfc ( ig ) )
              ENDDO
              k = k + 1
            ENDIF
            ib = ib + 1
          ENDIF
        ENDDO
        dimension_span = k - 1
        IF ( dimension_span .EQ. 0 .AND. peinf%inode .EQ. 0 ) &
          WRITE ( 0, 201 ) ik, is, inc
        DO j = 1, dimension_span
          x = 0.0D0
          DO ig = 1, ngk_l
            x = x + phi ( ig, j ) **2
          ENDDO
#ifdef MPI
          y = x
          CALL MPI_Allreduce ( y, x, 1, MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
          x = sqrt ( x )
          DO ig = 1, ngk_l
            phi ( ig, j ) = phi ( ig, j ) / x
          ENDDO
        ENDDO
!
! the Gram-Schmidt process begins
!
        reduced_span = 1
        DO ig = 1, ngk_l
          psi ( ig, 1 ) = phi ( ig, 1 )
        ENDDO
        DO j = 1, dimension_span - 1
          DO ig = 1, ngk_l
            vec ( ig ) = phi ( ig, j + 1 )
          ENDDO
          DO k = 1, reduced_span
            x = 0.0D0
            DO ig = 1, ngk_l
              x = x + phi ( ig, j + 1 ) * psi ( ig, k )
            ENDDO
#ifdef MPI
            y = x
            CALL MPI_Allreduce ( y, x, 1, MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
            DO ig = 1, ngk_l
              vec ( ig ) = vec ( ig ) - psi ( ig, k ) * x
            ENDDO
          ENDDO
          x = 0.0D0
          DO ig = 1, ngk_l
            x = x + vec ( ig ) **2
          ENDDO
#ifdef MPI
          y = x
          CALL MPI_Allreduce ( y, x, 1, MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
          x = sqrt ( x )
          IF ( x .GT. eps6 ) THEN
            reduced_span = reduced_span + 1
            DO ig = 1, ngk_l
              psi ( ig, reduced_span ) = vec ( ig ) / x
            ENDDO
          ENDIF
        ENDDO
!
! the Gram-Schmidt process ends
!
        IF ( reduced_span .LT. imap ( i, is ) - inc + 1 .AND. peinf%inode .EQ. 0 ) &
          WRITE ( 0, 202 ) ik, is, inc
        DO ib = inc, imap ( i, is )
          DO ig = 1, ngk_l
            wfn_d ( ig, ib, is, ik ) = &
              CMPLX ( psi ( ig, ib - inc + 1 ), 0.0D0 )
          ENDDO
        ENDDO
        inc = imap ( i, is ) + 1
      ENDDO
    ENDDO
    
    SAFE_DEALLOCATE ( vec )
    SAFE_DEALLOCATE ( phi )
    SAFE_DEALLOCATE ( psi )
    SAFE_DEALLOCATE ( wfc )
    SAFE_DEALLOCATE ( null_map )
    SAFE_DEALLOCATE ( inull )
    SAFE_DEALLOCATE ( inums )
    SAFE_DEALLOCATE ( imap )
    
  ENDDO
  
  POP_SUB(real_wfng)
  
  RETURN
  
201 FORMAT(1x,"WARNING: failed Gram-Schmidt dimension span",/,14x, &
      "k-point =",1x,i4.4,1x,"spin =",1x,i1.1,1x,"band =",1x,i4.4,/)
202 FORMAT(1x,"WARNING: failed Gram-Schmidt reduced span",/,14x, &
      "k-point =",1x,i4.4,1x,"spin =",1x,i1.1,1x,"band =",1x,i4.4,/)
  
END SUBROUTINE real_wfng

!-------------------------------------------------------------------------------

end module real_m

