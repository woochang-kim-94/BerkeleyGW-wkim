#include "f_defs.h"

!*************************************************************************
Module AtomBasisMod
!*************************************************************************
! This module contains information on the atom basis.
! The Atoms portion contains information on the atoms in the basis, 
! mainly their type and position.
!*************************************************************************
! The module contains the public functions:
!    AtomBasisInit            reads information on atom positions and types
!                         allocates the space for the dynamic variables
!    AtomBasisDestroy         deallocates space for the dynamic variables
!*************************************************************************
! The module contains the type:
! AtomBasisT which contains the variables:
!    totalNumAtoms        total number of atoms in the unit cell
!    atomCart(k,n)        k-th component (in cartesian coordinates) of
!                         the position of the n-th atom, n = 1, .., totalNumAtoms
!    atomLat(k,n)         k-th component (in primitive lattice 
!                         coordinates) of the position of the n-th atom
!    speciesOfAtom(n)     The species number ( = 1, .., numSpecies) of the n-th atom
!    numSpecies           number of species of atoms
!    numAtomsOfSpecies(i) number of atoms of species i = 1, .., numSpecies in the unit cell
!    elementLabel(i)      chemical symbol for the atom species i
!    atomicNumber(i)      atomic number for the atom species i -- added by DAS
!*************************************************************************

  use SysParams,    only : double
  use push_pop_m
  implicit none

  type AtomBasisT  
    integer                      :: totalNumAtoms
    real(double),    pointer     :: atomCart(:,:)
    real(double),    pointer     :: atomLat(:,:)
    integer,         pointer     :: speciesOfAtom(:) 

    integer                      :: numSpecies   
    integer,         pointer     :: numAtomsOfSpecies(:) 
    character*4,     pointer     :: elementLabel(:)
    integer,         pointer     :: atomicNumber(:) ! added by DAS
  end type AtomBasisT


  integer, parameter           :: scaledByLatVec = 1
  integer, parameter           :: scaledCart = 2

  private scaledByLatVec
  private scaledCart

  contains

!*************************************************************************
  Subroutine AtomBasisInit( atomBasis, ndim, lattice, tagHandler )
!*************************************************************************

    use TagHandlerMod, only : FindTag, StringCompare, tagHandlerT
    use Sysparams,  only : one, twopi
    use LatticeMod,    only : LatticeT

    implicit none

    type( AtomBasisT ),  pointer :: atomBasis
    type( LatticeT ),    pointer :: lattice
    type( TagHandlerT ), pointer :: tagHandler
    integer, intent(in)  :: ndim

    integer       :: atomFormat ! Format for positions in input file  
    integer       :: i, j, error
    character(80) :: buff
    logical       :: default
    real(double), allocatable :: atomCoorTemp(:,:)

    PUSH_SUB(AtomBasisInit)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read in the number of atomic species.  If there is no entry in the 
    ! input file or there is an error reading the number a default of 1 will
    ! be used.

    allocate(atomBasis, stat = error)

    if(error /= 0 ) stop 'Error, atomBasis allocation failed.'

    default = .True.
    call FindTag(tagHandler, 'NumberOfSpecies', error)

    if(error .eq. 0) then

    read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) atomBasis%numSpecies

      if(error .eq. 0) then

        print *, 'NumberOfSpecies ', atomBasis%numSpecies
        default = .False.

      end if
    end if

    if(default) then

      atomBasis%numSpecies = 1
      print *, 'Using default value of NumberOfSpecies  1.'

    end if
    ! numSpecies has been initialized.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate the space for the number of atoms for each element and the
    ! chemical label for each element.
    allocate(atomBasis%numAtomsOfSpecies(atomBasis%numSpecies), &
           & atomBasis%elementLabel(atomBasis%numSpecies), &
           & atomBasis%atomicNumber(atomBasis%numSpecies), stat = error)
    if(error /= 0) stop 'Error allocating atom space in AtomBasisInit.'
    ! Allocation succeeded.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read in the number of atoms for each element.  If number is not in
    ! the file or there is an error reading the program will halt.
    call FindTag( tagHandler, 'NumberOfAtoms', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt ='(I10)', iostat = error) atomBasis%totalNumAtoms
      if(error /= 0) stop 'Error, couldn''t read NumberOfAtoms.'
      print *, 'NumberOfAtoms ', atomBasis%totalNumAtoms

    else

      stop 'Error, could not find the NumberOfAtoms.'

    end if
    ! totalNumAtoms initialized.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read in the Labels for each element.  If there is an error, or
    ! the labels are not in the input file then the program will exit.
    call FindTag( tagHandler, 'ChemicalSpeciesLabel', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)') i

      print *, 'ChemicalSpeciesLabel '
      do i = 1, atomBasis%numSpecies

        read(unit = tagHandler%fileno, fmt = *, iostat = error) j, atomBasis%atomicNumber(i), atomBasis%elementLabel(i)
        if(error /= 0) stop 'Error reading ChemicalSpeciesLabel.'
        print *, j, ' ', atomBasis%atomicNumber(i), ' ', atomBasis%elementLabel(i)
        if(j /= i) stop 'Wrong order of ChemicalSpeciesLabel.'

      end do
    else

      stop 'Error, could not find ChemicalSpeciesLabel.'

    end if
    ! elementLabel has been initialized.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 1.  Allocate temporary space to store the atom coordinates before they
    ! are sorted into atomCart.  
    ! 2.  Then read in the atomic coordinates into the array atoms. 
    ! The format of the atom 
    ! coordinates in the input file is a set of ndim = numDimensions coordinates
    ! followed by an integer that refers to the place of the element label 
    ! in the list of labels in the input file.  
    ! 3.  Read 'AtomicCoordinatesFormat' and assign the integer atomFormat
    ! to treat the coordinates input file 
    ! 4.  If the coordinates are not in the input file or there is an
    ! error reading them the program will exit.  
    ! 5. atomCart atomLat are then allocated and the
    ! atomic coordinates are sorted and stored in them.

    allocate(atomCoorTemp(ndim, atomBasis%totalNumAtoms), &
           & atomBasis%speciesOfAtom(atomBasis%totalNumAtoms), stat = error)
    if(error /= 0) stop 'Error allocating atoms in AtomBasisInit.'
    
    call FindTag( tagHandler, 'AtomicCoordinatesFormat', error)
    
    default = .True.

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(A80)', iostat = error) buff

      if(error .eq. 0) then

        if(StringCompare(buff, 'ScaledCartesian')) then

          atomFormat = scaledCart
          print *, 'AtomicCoordinatesFormat ScaledCartesian'
          default = .False.

        else if(StringCompare(buff, 'ScaledByLatticeVectors')) then

          atomFormat = scaledByLatVec
          print *, 'AtomicCoordinatesFormat ScaledByLatticeVectors'
          default = .False.

        end if
      end if
    end if

    if(default) then
      atomFormat = scaledByLatVec
      print *, 'Using default AtomicCoordinatesFormat ScaledByLatticeVectors'
    end if

    call FindTag( tagHandler, 'AtomicCoordinatesAndAtomicSpecies', error)
    atomBasis%numAtomsOfSpecies = 0

    if(error .eq. 0 .and. (atomFormat .ne. scaledByLatVec &
           & .or. atomFormat .ne. scaledCart)) then

      read(unit = tagHandler%fileno, fmt = '(I10)') i

      do i = 1, atomBasis%totalNumAtoms

        read( unit=tagHandler%fileno, fmt=*, iostat=error) atomCoorTemp(:,i), atomBasis%speciesOfAtom(i)

        if(error /= 0) stop 'Error reading AtomicCoord(Cart,Lat).'

        atomBasis%numAtomsOfSpecies(atomBasis%speciesOfAtom(i)) = &
                 & atomBasis%numAtomsOfSpecies(atomBasis%speciesOfAtom(i)) + 1

      end do

    else

      stop 'Error, could not find AtomicCoord(Cart,Lat).'

    end if

    allocate(atomBasis%atomCart(ndim, atomBasis%totalNumAtoms), &
           & atomBasis%atomLat(ndim, atomBasis%totalNumAtoms), &
           & stat = error)
    if(error /= 0) stop 'Error allocating atomCart in AtomBasisInit.'

    if(atomFormat .eq. scaledCart) then

      call AtomSort(atomBasis, atomBasis%atomCart, atomCoorTemp)
      atomBasis%atomCart = atomBasis%atomCart * lattice%latticeConst
      atomBasis%atomLat = matmul(transpose(lattice%bLatVec),atomBasis%atomCart) / twopi

      ! do i = 1, ndim
      !   do j = 1, totalNumAtoms
      !     atomLat(i,j) = sum(lattice%bLatVec(:,i) * atomCart(:,j)) / twopi
      !   end do
      ! end do

    else if(atomFormat .eq. scaledByLatVec) then

      call AtomSort(atomBasis, atomBasis%atomLat, atomCoorTemp)
      atomBasis%atomCart = matmul(lattice%aLatVec,atomBasis%atomLat)

      ! do i = 1, ndim
      !   do j = 1, totalNumAtoms
      !     atomCart(i,j) = sum(lattice%aLatVec(i,:) * atomLat(:,j))
      !   end do
      ! end do

    end if

    print *, 'AtomicCoordinatesAndAtomicSpecies'
    do i = 1, atomBasis%totalNumAtoms
      print *, atomCoorTemp(:,i), atomBasis%speciesOfAtom(i)
    end do
    ! atomCart and atomLat have been initialized.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate temp space for reading in atomic coordinates.
    deallocate(atomCoorTemp, stat = error)
    if(error /= 0) stop 'Error deallocating in AtomBasisInit.'
    ! All temp space has been freed.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    POP_SUB(AtomBasisInit)

  end Subroutine AtomBasisInit

!*************************************************************************
  Subroutine AtomBasisDestroy( atomBasis )
!*************************************************************************

!  Deallocate space for atomic positions and information

    implicit none

    type(AtomBasisT), pointer :: atomBasis
    integer :: error

    PUSH_SUB(AtomBasisDestroy)

    if( .not. associated( atomBasis )) stop 'Error: atomBasis not allocated.'
    deallocate(atomBasis%numAtomsOfSpecies, atomBasis%elementLabel, &
             & atomBasis%speciesOfAtom, atomBasis%atomicNumber, stat = error)
    if(error /= 0) stop 'Error deallocating numAtoms in AtomBasisDestroy.'

    deallocate(atomBasis%atomCart, atomBasis%atomLat, stat = error)
    if(error /= 0) stop 'Error deallocating atomCart in AtomBasisDestroy.'

    deallocate(atomBasis, stat = error)
    if(error /= 0) stop 'Error deallocating atomCart in AtomBasisDestroy.'

    nullify( atomBasis )

    POP_SUB(AtomBasisDestroy)

  end Subroutine AtomBasisDestroy

!*************************************************************************
  Subroutine AtomSort(atomBasis, atomCoor, atomCoorTemp)
!*************************************************************************

!  Sort atoms so that they are grouped by element

  implicit none
  type( atomBasisT ), pointer :: atomBasis
  real( double ),     pointer :: atomCoor(:,:)
  real( double ),     intent( in )  :: atomCoorTemp(:,:)
  integer :: i, j, k

    PUSH_SUB(AtomSort)
    
    k = 0

    do i = 1, atomBasis%numSpecies
      do j = 1, atomBasis%totalNumAtoms
        if(atomBasis%speciesOfAtom(j) .eq. i) then
          k = k + 1
          atomCoor(:, k) = atomCoorTemp(:, j)
        end if
      end do
    end do

    k = 0

    do i = 1, atomBasis%numSpecies
      do j = 1, atomBasis%numAtomsOfSpecies(i)
        k = k + 1
        atomBasis%speciesOfAtom(k) = i
      end do
    end do

    POP_SUB(AtomSort)

  end Subroutine AtomSort

!*************************************************************************
  Subroutine AtomsPrint( atomBasis )
!*************************************************************************

    implicit none
    type(AtomBasisT), pointer :: atomBasis
    integer :: i

    PUSH_SUB(AtomsPrint)

    print *, 'Atomic cartesian coordinates:'
      do i = 1, atomBasis%totalNumAtoms
        print *, atomBasis%atomCart(:,i),  &
          & atomBasis%elementLabel(atomBasis%speciesOfAtom(i))
      end do

    print *, 'Atomic lattice vector coordinates:'
      do i = 1, atomBasis%totalNumAtoms
        print *, atomBasis%atomLat(:,i), &
          & atomBasis%elementLabel(atomBasis%speciesOfAtom(i))
        end do

    POP_SUB(AtomsPrint)

  end Subroutine AtomsPrint

!*************************************************************************
  Subroutine XYZCoordPrint( atomBasis, ndim, lattice )
!*************************************************************************
  
!RMM    use DimensionsMod, only : DimensionsT
        use LatticeMod,    only : LatticeT

    use message_m
    implicit none

    type( AtomBasisT ), pointer :: atomBasis
!RMM    type( DimensionsT ), pointer :: dimensions
        type( LatticeT ),    pointer :: lattice
        integer, intent(in)  :: ndim
    integer :: i, k

    PUSH_SUB(XYZCoordPrint)

    call open_file(7, file = 'coord.xyz',status='replace')
    write(7, '(I4)') 10 * atomBasis%totalNumAtoms

    do k = 1, atomBasis%totalNumAtoms

      write(7, '(A2, I1)', advance = "no") &
         & atomBasis%elementLabel(atomBasis%speciesOfAtom(k)), k
      write(7, '(3F13.6)') atomBasis%atomCart(:,k)

    end do

    do i = 1, ndim

      do k = 1, atomBasis%totalNumAtoms

        write(7, '(A2, I1)', advance = "no") &
            & atomBasis%elementLabel(atomBasis%speciesOfAtom(k)), k
        write(7, '(3F13.6)') atomBasis%atomCart(:,k) + lattice%aLatVec(:,i)

      end do

      do k = 1, atomBasis%totalNumAtoms

        write(7, '(A2, I1)', advance = "no") &
            & atomBasis%elementLabel(atomBasis%speciesOfAtom(k)), k
        write(7, '(3F13.6)') atomBasis%atomCart(:,k) - lattice%aLatVec(:,i)

      end do
    end do

    call close_file(7)

    POP_SUB(XYZCoordPrint)

  end Subroutine XYZCoordPrint

end Module AtomBasisMod

