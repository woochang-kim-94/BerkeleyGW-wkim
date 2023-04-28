#include "f_defs.h"

!*************************************************************************
Module NormalizeMod
!*************************************************************************

  use push_pop_m
  implicit none

  contains

!*************************************************************************
  Subroutine EigenStatesOrth( eigenStates, iBand, kPoint )
!*************************************************************************

!   This Subroutine orthonormalizes the starting vector for iBand against
!   all the lower bands.

    use EigenStatesMod, only : eigenStatesT

    integer, intent(in) :: iBand, kPoint
    type(eigenStatesT), pointer :: eigenStates

    PUSH_SUB(EigenStatesOrth)

    ! Orthogonalize to all the LOWER bands
    Call LowerBandOrthogonalize( eigenStates, &
        & eigenStates%eigenVectors( :, iBand, kPoint ), iBand, kPoint )

    ! Normalize the state-vector
    Call VecNormalize( eigenStates%eigenvectors( :, iBand, kPoint ))

!    print *, 'OrthStartVec:',eigenStates%eigenvectors(:,iBand)

    POP_SUB(EigenStatesOrth)
  end Subroutine EigenStatesOrth

!*************************************************************************
  Subroutine LowerBandOrthogonalize( eigenStates, vector, iBand, kPoint)
!*************************************************************************

!   This Subroutine Orthogonalizes the vector corresponding to iBand against 
!   all the bands (in eigenStates) below iBand.

    use SysParams, only : double
    use EigenStatesMod, only : eigenStatesT

    type( eigenStatesT ), pointer :: eigenStates
    complex( double ), intent( inout ) :: vector(:)
    integer, intent( in ) :: iBand, kPoint
    
    integer :: jBand
    complex( double ) :: orthfac

    PUSH_SUB(LowerBandOrthogonalize)

    do jBand = 1, iBand-1

       !    Calculate the orthogonalising factor for the j th band wrt 
       ! the i th band.
       orthfac = Dot_Product( &
               & eigenStates%eigenvectors( :, jBand, kPoint), vector)

       !   Calculate the Orthogonalising vector between i and j bands
       vector = vector - orthfac * eigenStates%eigenvectors( :, jband, kPoint)
          
    end do
    
    vector(eigenStates%hSize+1:eigenStates%hDim) = cmplx(0.d0,0.d0,kind=double)

    POP_SUB(LowerBandOrthogonalize)
  end Subroutine LowerBandOrthogonalize

!*************************************************************************
  Subroutine VecNormalize( vector )
!*************************************************************************

!   This Subroutine normalizes vector

    use SysParams, only : double

    complex( double ), intent( inout ) :: vector(:)
    real( double ) :: norm

    PUSH_SUB(VecNormalize)

    norm = Dot_Product( vector, vector )
    norm = sqrt( norm ) + 1.d-32
    vector = vector / norm

    POP_SUB(VecNormalize)
  end Subroutine VecNormalize

!*************************************************************************
  Subroutine SelfBandOrthogonalize( eigenStates, vector, iBand, kPoint )
!*************************************************************************

!   This Subroutine orthogonalizes the vector corresponding to iBand, against
!   the state-vector iBand contained in eigenStates.

    use EigenStatesMod, only : eigenStatesT
    use SysParams, only : double

    integer, intent( in ) :: iBand, kPoint
    type( eigenStatesT ), pointer :: eigenStates
    complex( double ), dimension(:) :: vector
    complex( double ) :: orthfac

    PUSH_SUB(SelfBandOrthogonalize)

    orthfac = Dot_Product( eigenStates%eigenvectors( :, iBand, kPoint ), vector )

    vector = vector - orthfac * eigenStates%eigenvectors( :, iBand, kPoint )
       
    vector(eigenStates%hSize+1:eigenStates%hDim) = cmplx(0.d0,0.d0,kind=double)

    POP_SUB(SelfBandOrthogonalize)
  end Subroutine SelfBandOrthogonalize

!*************************************************************************
  Subroutine BandOrthogonalize( eigenStates, vector, iBand, kPoint )
!*************************************************************************

!   This Subroutine orthogonalizes the vector corresponding to iBand against
!   all the bands below and including iBand.

    use SysParams, only : double
    use EigenStatesMod, only : eigenStatesT

    integer, intent(in) :: iBand, kPoint
    type(eigenStatest), pointer :: eigenStates
    complex( double ), dimension(:) :: vector

    PUSH_SUB(BandOrthogonalize)

    ! First Orthogonalize against iBand
    Call SelfBandOrthogonalize( eigenStates, vector, iBand, kPoint )

    ! Next Orthogonalize that against all the lower bands
    Call LowerBandOrthogonalize( eigenStates, vector, iBand, kPoint )

    POP_SUB(BandOrthogonalize)
  end Subroutine BandOrthogonalize

end Module NormalizeMod
