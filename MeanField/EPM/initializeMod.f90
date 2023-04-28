#include "f_defs.h"

!*************************************************************************
Module InitializeMod
!*************************************************************************

  use push_pop_m
  implicit none

  contains

!*************************************************************************
  Subroutine MinVecAllocate( minVec, nSize )
!*************************************************************************

    use TypeMod, only : minVecT

    integer, intent( in ) :: nSize
    type( minVecT ), pointer :: minVec

    integer :: error

    PUSH_SUB(MinVecAllocate)
    if( associated( minVec )) stop 'Error minVec already allocated'

    allocate( minVec, stat = error )

    if( error /= 0 ) stop 'Error, Vector allocation failed.'

    allocate( minVec%Phi( NSize ), minVec%StDesc( nSize), minVec%EPhi(nSize), &
         & stat = error )

    if( error /= 0 ) stop 'Error in allocating minVec space'

    POP_SUB(MinVecAllocate)
  end Subroutine MinVecAllocate

!*************************************************************************
  Subroutine MinVecDestroy(minVec)
!*************************************************************************

    use TypeMod, only : minVecT

    type(minVecT), pointer :: minVec

    PUSH_SUB(MinVecDestroy)

    deallocate(minVec%StDesc,minVec%Phi,minVec%EPhi)
    deallocate(minVec)

    POP_SUB(MinVecDestroy)
  end Subroutine MinVecDestroy

!*************************************************************************
  Subroutine EigenStatesInit( eigenStates, iband, kPoint )
!*************************************************************************

    use SysParams, only : double
    use EigenStatesMod, only : eigenStatesT

    integer, intent(in) :: iband, kPoint
    type(EigenStatesT),  pointer :: eigenStates

    integer :: nSize,r
    real( double ) :: R1,rnd1,rnd2,del

    PUSH_SUB(EigenStatesInit)

    nSize = eigenStates%hSize
    del = 0.1

    eigenStates%eigenVectors( :, iband, kPoint) = cmplx(0.d0,0.d0,kind=double)

    do r = 1, nSize
       call Random_Number(R1)
       rnd1 = -1.d0 + 2.d0*R1
       call Random_Number(R1)
       rnd2 = -1.d0 + 2.d0*R1
       eigenStates%eigenVectors( r, iband, kPoint) = cmplx(rnd1,rnd2,kind=double)
    end do

    eigenStates%ePsi( :, kPoint ) = cmplx(0.d0,0.d0,kind=double)

    POP_SUB(EigenStatesInit)
  end Subroutine EigenStatesInit

!*************************************************************************
  Subroutine MinVecInit( minVec )
!*************************************************************************

    use SysParams, only : double
    use TypeMod, only : minVecT

    type(minVecT), pointer :: minVec

    PUSH_SUB(MinVecInit)

    ! Initialize all these vectors to 0.d0
    minVec%cosine = 1.d0
    minVec%sine = 0.d0
    minVec%Phi = cmplx(0.d0,0.d0,kind=double)
    minVec%EPhi = cmplx(0.d0,0.d0,kind=double)
    minVec%StDesc = cmplx(0.d0,0.d0,kind=double)

    POP_SUB(MinVecInit)
  end Subroutine MinVecInit

end Module InitializeMod

