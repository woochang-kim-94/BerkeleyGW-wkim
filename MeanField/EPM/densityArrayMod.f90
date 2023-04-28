#include "f_defs.h"

!*************************************************************************
Module densityArrayMod
!*************************************************************************
! Module for density defined on a grid
! densityArrayT is a type containing the wavefunction and density information
!*************************************************************************
! NEED TO CHANGE TO MAKE MORE UNIVERSAL
!*************************************************************************

  use SysParams,    only : double
  use push_pop_m

  type densityArrayT
    complex(double), pointer :: densityArray(:,:,:), waveFunction(:,:,:), rhog(:)
    integer                  :: gridDimensions(3), gridSize
  end type densityArrayT

  contains

!*************************************************************************
  Subroutine densityArrayInit( densityArray, ndim, FFTgrid, numGVectors )
!*************************************************************************

    implicit none
    
    type( densityArrayT ), pointer :: densityArray
    integer, intent(in) :: ndim
    integer, intent(in) :: FFTgrid(3)
    integer, intent(in) :: numGVectors
    integer :: i, error

    PUSH_SUB(densityArrayInit)

    if( associated( densityArray )) stop 'Error, densityArray already allocated.'
    allocate( densityArray, stat=error )
    if( error /= 0 ) stop 'Error, densityArray allocation failed.'
    nullify( densityArray%densityArray )
    nullify( densityArray%waveFunction )
    nullify( densityArray%rhog )

    densityArray%gridSize = 1

    densityArray%gridDimensions = 1
    do i = 1, ndim
      densityArray%gridDimensions(i) = FFTgrid(i)
      densityArray%gridSize = densityArray%gridSize * densityArray%gridDimensions(i)
    end do

    allocate( densityArray%densityArray(densityArray%gridDimensions(1), &
                                        densityArray%gridDimensions(2), densityArray%gridDimensions(3)), &
              densityArray%waveFunction(densityArray%gridDimensions(1), &
                                        densityArray%gridDimensions(2), densityArray%gridDimensions(3)), &
              densityArray%rhog(numGVectors), stat = error )
    if(error /= 0) stop 'Error allocating density arrays in densityArrayInit.'

    POP_SUB(densityArrayInit)
  end Subroutine densityArrayInit

!*************************************************************************
  Subroutine densityArrayDestroy( densityArray )
!*************************************************************************

    implicit none
    
    type( densityArrayT ), pointer :: densityArray
    integer :: error

    PUSH_SUB(densityArrayDestroy)

    if( .not. associated( densityArray )) stop 'Error, densityArray not allocated.'

    deallocate( densityArray%densityArray, &
              & densityArray%waveFunction, densityArray%rhog, stat = error )
    if(error /= 0) stop 'Error deallocating density arrays.'

    deallocate( densityArray, stat = error )
    if(error /= 0) stop 'Error, densityArray deallocation failed.'

    nullify( densityArray )

    POP_SUB(densityArrayDestroy)
  end Subroutine densityArrayDestroy

end Module densityArrayMod

