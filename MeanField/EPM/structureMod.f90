#include "f_defs.h"

!*************************************************************************
Module StructureMod
!*************************************************************************
! This module contains information on the crystal structure.
! The Dimensions portion reads in the space dimension ndim used throughout 
! The Lattice portion contains information on the direct and reciprocal 
! Bravais lattices.
! The Atoms portion contains information on the atoms in the basis, 
! mainly their type and position.
!*************************************************************************
! The module contains the public functions:
!    StructureInit         calls DimensionsRead, LatticeInit, AtomsInit 
!    StructureDestory      deallocates the StructureMod module memory
!*************************************************************************
! The module contains the public variables:
!*************************************************************************
  
  use SysParams,    only : double
  use AtomBasisMod, only : AtomBasisT, AtomBasisInit, AtomBasisDestroy
  use LatticeMod,   only : LatticeT, LatticeInit, LatticeDestroy
  use push_pop_m
  implicit none

  type StructureT
    integer ndim         ! dimension of space used for lattce and atom positions
    type( AtomBasisT ),  pointer :: atomBasis
    type( LatticeT ),    pointer :: lattice
  end type StructureT

  contains

!*************************************************************************
  Subroutine StructureInit( structure1, tagHandler )
!*************************************************************************

    use  TagHandlerMod, only : tagHandlerT
    implicit none

    type( StructureT ), pointer :: structure1
    type( tagHandlerT ), pointer :: tagHandler

    integer :: error

    PUSH_SUB(StructureInit)

    if( associated( structure1 )) stop 'Error: structure already allocated.'

    allocate( structure1, stat = error )

    if( error /= 0 ) stop 'Error allocating structure.'

    nullify( structure1%atomBasis )
    nullify( structure1%lattice )

    call DimensionsRead( structure1%ndim, tagHandler )

    call LatticeInit( structure1%lattice, structure1%ndim, tagHandler )

    call AtomBasisInit( structure1%atomBasis, structure1%ndim, structure1%lattice, tagHandler )

    POP_SUB(StructureInit)
  end Subroutine StructureInit

!*************************************************************************
  Subroutine StructureDestroy( structure1 )
!*************************************************************************

    implicit none

    type( StructureT ), pointer :: structure1
    integer :: error

    PUSH_SUB(StructureDestroy)

    if( .not. associated( structure1 )) stop 'Error structure not allocated.'

    call LatticeDestroy( structure1%lattice )

    call AtomBasisDestroy( structure1%atomBasis )

    deallocate( structure1, stat = error )

    if( error /= 0 ) stop 'Error deallocating structure.'

    nullify( structure1 )

    POP_SUB(StructureDestroy)
  end Subroutine StructureDestroy

!*************************************************************************
  Subroutine DimensionsRead(ndim, tagHandler )
!*************************************************************************

    use TagHandlerMod, only : FindTag, tagHandlerT
    implicit none

    type( tagHandlerT ), pointer :: taghandler
    integer, intent(out) :: ndim
    integer :: error1, error
    logical :: default

    PUSH_SUB(DimensionsRead)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read in the number of dimensions.  If none is present or their is
    ! an error reading the number then the default value of 3 is used.  The
    ! number of dimensions must be from 1 to 3. 

    default = .True.

    call FindTag( tagHandler, 'NumberOfDimensions', error1)
    if(error1 .eq. 0) then
      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) ndim

      if(error.eq. 0 .and.ndim .gt. 0 .and.ndim.lt. 4) then
        print *, 'NumberOfDimensions ', ndim
        default = .False.
      end if

    end if

    if(default) then
      ndim = 3
      print *, 'Using default value of NumberOfDimensions ', ndim
    end if

    ! ndim = space dimension has been initialized.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    POP_SUB(DimensionsRead)
  end Subroutine DimensionsRead

end Module StructureMod

