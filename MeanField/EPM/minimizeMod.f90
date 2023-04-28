#include "f_defs.h"

!*************************************************************************
Module MinimizeMod
!*************************************************************************

  use NormalizeMod,  only : SelfBandOrthogonalize, VecNormalize
  use push_pop_m
  implicit none

  contains

!*************************************************************************
  Subroutine StDescentVec( eigenStates, minVec, iBand, kPoint )
!*************************************************************************

    use TypeMod, only : minVecT
    use EigenStatesMod, only : eigenStatesT

    integer, intent(in) :: iBand, kPoint
    type(eigenStatesT), pointer :: eigenStates
    type(minVecT), pointer :: minVec

    PUSH_SUB(StDescentVec)

    ! Evaluate the Steepest Descent direction for iBand
    minVec%StDesc = - eigenStates%ePsi( : , kPoint ) + &
         & eigenStates%eigenValues( iBand, kPoint)* &
         & eigenStates%eigenVectors( :, iBand, kPoint )

    POP_SUB(StDescentVec)
  end Subroutine StDescentVec

!*************************************************************************
  Subroutine ConjDirection( eigenStates, minVec, oldVec, iBand, kPoint )
!*************************************************************************

    use SysParams, only : double
    use TypeMod, only : minVecT
    use EigenStatesMod, only : eigenStatesT

    integer, intent(in) :: iBand, kPoint
    type(eigenStatesT), pointer :: eigenStates
    type(minVecT), pointer :: minVec,oldVec

    real( double ) :: Gamma,numer,denom

    PUSH_SUB(ConjDirection)

    !   Calculate Gamma, the phi is adjusted such that in the first 
    ! iteration Gamma becomes inconsequential
    numer = Dot_Product(( minVec%StDesc - oldVec%StDesc ), minVec%StDesc )
    denom = Dot_Product(( minVec%StDesc - oldVec%StDesc ), oldVec%Phi )

!   Gamma = 0.d0     ! This makes it the steepest descent method
    Gamma = -numer / ( denom + 10.d-40 )

    minVec%Phi = minVec%StDesc + Gamma * oldVec%Phi

    ! Orthogonalize Phi against the state-vector of the same band 'iBand'
    call SelfBandOrthogonalize( eigenStates, minVec%Phi, iBand, kPoint )

    ! Normalize Phi
    call VecNormalize( minVec%Phi )

!    print *, 'Conjgdir:',minVec%Phi

    ! Roll the current values into oldVec
    oldVec%sine = minVec%sine
    oldVec%cosine = minVec%cosine
    oldVec%Phi = minVec%Phi
    oldVec%EPhi = minVec%EPhi
    oldVec%StDesc = minVec%StDesc

    POP_SUB(ConjDirection)
  end Subroutine ConjDirection

end Module MinimizeMod
