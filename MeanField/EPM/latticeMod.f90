#include "f_defs.h"

!*************************************************************************
Module LatticeMod
!*************************************************************************
! This module contains information on the lattice.
! The Lattice portion contains information on the direct and reciprocal 
!   Bravais lattices.
!*************************************************************************
! The module contains the public functions:
!    LatticeInit           allocates arrays; reads lattice constant and vectors;
!                          calculates reciprocal lattice, metric tensors
!    LatticeDestroy        deallocates arrays for the dynamic variables 
!*************************************************************************
! The module contains the public variables:
!    latticeConst          scale factor for lengths
!    aLatVec(i,j)          i-th component of the (scaled) j-th lattice vector
!    aLatVecLen(j)         length of the j-th lattice vector (a.u.)
!    aLatMetricTensor(i,j) = aLatMetricTensor(j,i) metric tensor in
!                          real space, used to scale the integer multiples
!                          of the reciprocal lattice vectors into real values
!    bLatVec(i,j)          i-th component of the j-th reciprocal lattice vector
!    bLatVecLen(j)         length of the j-th reciprocal lattice vector
!    bLatMetricTensor(i,j) = bLatMetricTensor(j,i): metric tensor in
!                          reciprocal space, used to scale the integer multiples
!                          of the reciprocal lattice vectors into real values
!*************************************************************************
! The module contains the private functions:
!    InvertMatrix(a)    returns the inverse of a.
!*************************************************************************
  
  use SysParams, only : double
  use push_pop_m
  implicit none

  type LatticeT
    real(double)             :: latticeConst
    real(double),    pointer :: aLatVec(:,:)
    real(double),    pointer :: aLatVecLen(:)
    real(double),    pointer :: aLatMetricTensor(:,:)
    real(double),    pointer :: bLatVec(:,:)
    real(double),    pointer :: bLatVecLen(:)
    real(double),    pointer :: bLatMetricTensor(:,:)
  end type LatticeT

  contains

!*************************************************************************
  Subroutine LatticeInit( lattice, ndim, tagHandler )
!*************************************************************************

    use TagHandlerMod,    only : FindTag, StringCompare, tagHandlerT
    use SysParams,     only : one, twopi, bohr2ang, double

    implicit none

    type( LatticeT ),    pointer :: lattice
    type( tagHandlerT ), pointer :: tagHandler
    integer, intent(in)  :: ndim

    integer :: i, error
    character :: c
    character(80) :: buff
    logical :: default
    real( double ) :: temp

    PUSH_SUB(LatticeInit)
    if( associated( lattice )) stop 'Error: lattice already associated.'

    allocate( lattice, stat = error )

    if( error /= 0 ) stop 'Error allocation lattice.'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate space for the lattice vectors, the reciprocal lattice
    ! vectors, and the metric tensor.
    allocate( lattice%aLatVec(ndim,ndim), &
            & lattice%bLatVec(ndim,ndim), &
            & lattice%aLatMetricTensor(ndim,ndim), &
            & lattice%bLatMetricTensor(ndim,ndim), &
            & lattice%aLatVecLen(ndim), lattice%bLatVecLen(ndim), &
            & stat = error)

    if(error /= 0) stop 'Error allocating vector space in LatticeInit.'
    ! Allocation succeeded.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read in the lattice constant.  If there is no lattice constant in
    ! the input file or if there is an error reading it then the default
    ! value of 1.0d0 will be used.

    print *, 'Lattice Constant = scale factor for lengths'
    default = .True.

    temp = 1.0D0
    call FindTag( tagHandler, 'LatticeConstantFormat', error)
    if(error .eq. 0) then
      read(unit = tagHandler%fileno, fmt = '(A)', iostat = error) buff
      if (error .eq. 0) then
        if (StringCompare(buff, 'Angstrom')) temp = bohr2ang
      endif
    endif

    call FindTag( tagHandler, 'LatticeConstant', error)
    if(error .eq. 0) then
      read(unit = tagHandler%fileno, fmt = *, iostat = error) lattice%latticeConst

      if(error .eq. 0) then
        lattice%latticeConst = lattice%latticeConst / temp
        print *, 'LatticeConstant ', lattice%latticeConst
        default = .False.
      end if

    end if

    if(default) then
      lattice%latticeConst = one
      print *, 'Using default value of LatticeConstant 1.0.'
    end if
    ! latticeConst has been initialized
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read in the lattice vectors and scale them by the lattice constant.
    ! Calculate the lengths of the lattice vectors.
    ! If the lattice vectors are missing or there is an error reading them
    ! the program will halt.

    call FindTag( tagHandler, 'LatticeVectors', error)
    print *, 'Lattice Vectors as read in (unscaled) LatticeVectors'

    if(error .eq. 0) then ! read dummy variable c to start 

      ! folowing read staement on next line

      read(unit = tagHandler%fileno, fmt ='(A1)') c  

      do i = 1, ndim
        read(unit = tagHandler%fileno, fmt = *, iostat = error) lattice%aLatVec(:,i)
        if(error /= 0) stop 'Error, couldn''t read LatticeVector.'
        print *, lattice%aLatVec(:,i)
      end do

    else
      stop 'Error, could not find lattice%aLatVec.'
    endif

    lattice%aLatVec = lattice%latticeConst * lattice%aLatVec ! rescale lattice vectors

    lattice%aLatMetricTensor = matmul(TRANSPOSE(lattice%aLatVec),lattice%aLatVec)

    do i = 1, ndim         ! determine lattice vector lengths
      lattice%aLatVecLen(i) = SQRT(dot_product(lattice%aLatVec(:,i),lattice%aLatVec(:,i)))
    end do

    ! lattice%aLatVec and lattice%aLatVecLen has been initialized and scaled.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate the reciprocal lattice vectors using Invertmatrix.  Then
    ! calculate the metric tensor for the reciprocal lattice vectors.

    ! calculate reciprocal lattice
    lattice%bLatVec=twopi*TRANSPOSE(InvertMatrix(lattice%aLatVec))

    ! reciprocal space metric tensor
    lattice%bLatMetricTensor = matmul(TRANSPOSE(lattice%bLatVec),lattice%bLatVec)

    do i = 1, ndim          ! determine lattice vector lengths
      lattice%bLatVecLen(i) = SQRT(dot_product(lattice%bLatVec(:,i),lattice%bLatVec(:,i)))
    end do

    !  lattice%bLatVec and lattice%aLatMetricTensor have been initialized.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    POP_SUB(LatticeInit)
  end Subroutine LatticeInit

!*************************************************************************
  Subroutine LatticeDestroy( lattice )
!*************************************************************************

! Deallocate space for the lattice information 

    implicit none

    type( LatticeT ), pointer :: lattice
    integer :: error

    PUSH_SUB(LatticeDestroy)

    deallocate(lattice%aLatVec, lattice%bLatVec, lattice%aLatMetricTensor, lattice%bLatmetricTensor, &
             & lattice%aLatVecLen, lattice%bLatVecLen, stat = error)
    deallocate(lattice, stat = error)
    if(error /= 0) stop 'Error deallocating in LatticeDestroy.'

    POP_SUB(LatticeDestroy)
  end Subroutine LatticeDestroy

!*************************************************************************
  Subroutine LatticePrint( lattice )
!*************************************************************************

    use SysParams,    only : twopi

    implicit none

    type( LatticeT ), pointer :: lattice

    PUSH_SUB(LatticePrint)

    print *, ' '
    print *, 'The scaled lattice vectors are:'
    print *, lattice%aLatVec
    print *, ' '
    print *, 'The reciprocal lattice vectors (in units of 2pi/a) are:'
    print *, lattice%bLatVec/twopi

    POP_SUB(LatticePrint)
  end Subroutine LatticePrint

!*************************************************************************
  Function InvertMatrix(m)
!*************************************************************************

    implicit none

    real(double), intent(in) :: m(3, 3)

    real(double) :: InvertMatrix(3, 3)

    real(double) :: a(3, 3)
    real(double) :: b

    PUSH_SUB(InvertMatrix)

    a(1, 1) =   m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2)
    a(2, 1) = - m(2, 1) * m(3, 3) + m(2, 3) * m(3, 1)
    a(3, 1) =   m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1)
    a(1, 2) = - m(1, 2) * m(3, 3) + m(1, 3) * m(3, 2)
    a(2, 2) =   m(1, 1) * m(3, 3) - m(1, 3) * m(3, 1)
    a(3, 2) = - m(1, 1) * m(3, 2) + m(1, 2) * m(3, 1)
    a(1, 3) =   m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2)
    a(2, 3) = - m(1, 1) * m(2, 3) + m(1, 3) * m(2, 1)
    a(3, 3) =   m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)

    b = m(1, 1) * a(1, 1) + m(1, 2) * a(2, 1) + m(1, 3) * a(3, 1)

    if (abs(b) .lt. 1.0d-6) stop 'Singular Matrix in InvertMatrix.'

    InvertMatrix(1, 1) = a(1, 1) / b
    InvertMatrix(2, 1) = a(2, 1) / b
    InvertMatrix(3, 1) = a(3, 1) / b
    InvertMatrix(1, 2) = a(1, 2) / b
    InvertMatrix(2, 2) = a(2, 2) / b
    InvertMatrix(3, 2) = a(3, 2) / b
    InvertMatrix(1, 3) = a(1, 3) / b
    InvertMatrix(2, 3) = a(2, 3) / b
    InvertMatrix(3, 3) = a(3, 3) / b

    POP_SUB(InvertMatrix)
  end function InvertMatrix

end module LatticeMod
