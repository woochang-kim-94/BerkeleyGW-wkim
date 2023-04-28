#include "f_defs.h"

!*************************************************************************
Module KPointsMod
!*************************************************************************
! This module contains information about the k-points.
!*************************************************************************
! The module contains the public functions:
!    KPointsInit         reads information on k points and 
!                        allocates the space for the dynamic variables
!    KPointsDestroy      deallocates space for the dynamic variables
!*************************************************************************
! The module contains the type:
! kPointsT which contains the variables:
!    numKPoints          total number of k-points used
!    kPointsLat(i,j)     i-th component (reciprocal lattice coordinates)
!                        of the j-th k-point
!    kPointsCart(i,j)    i-th component (cartesian coordinates * 2pi / a)
!                        of the j-th k-point
!    numBands            number of bands to be calculated
!*************************************************************************
! Information for graphs
!*************************************************************************

  use SysParams,     only : double
  use GraphMod,      only : graphT
  use push_pop_m
  implicit none

  save

  type kPointsT
    integer               :: numKPoints
    integer               :: kPointNumbers(3)
    real(double)          :: kPointOffsets(3)
    real(double), pointer :: kPointsCart(:,:)
    real(double), pointer :: kPointsLat(:,:)
    real(double), pointer :: kPointWeights(:)
    integer               :: numBands
    integer               :: numOccupiedBands
  end type kPointsT

  integer, parameter      :: scaledByRecipLatVec = 1
  integer, parameter      :: scaledByTwoPi = 2

  private scaledByRecipLatVec
  private scaledByTwoPi

  contains

!*************************************************************************
  Subroutine KPointsInit( kPoints, graphInfo, ndim, lattice, tagHandler )
!*************************************************************************

!    generate k-points around a BZ for plotting band structure plots such that
!    each line segment has a constant density of points to make the graphics
!    look good

    use SysParams,     only : one, zero, twopi
    use TagHandlerMod, only : FindTag, StringCompare, tagHandlerT
    use LatticeMod,    only : latticeT
    use GraphMod,      only : graphT

    implicit none

    type( kPointsT ),    pointer :: kPoints
    type( LatticeT ),    pointer :: lattice
    type( graphT ),      pointer :: graphInfo
    type( tagHandlerT ), pointer :: tagHandler
    integer, intent(in)          :: ndim
    real(double), allocatable    :: ktemp(:,:) ! temporary list of k-points

    integer       :: j,error,kPointsScale
    real(double)  :: normKPointWeights
    character(80) :: buff
    logical       :: default

    PUSH_SUB(KPointsInit)

    if( associated( kPoints )) stop 'Error, kPoints already allocated.'

    allocate( kPoints, stat = error )
    if( error /= 0 ) stop 'Error, kPoints allocation failed.'
    nullify( kPoints%kPointsCart )
    nullify( kPoints%kPointsLat )
    nullify( kPoints%kPointWeights )

    nullify( graphInfo )  ! used only if associated by another subroutine

    call FindTag( tagHandler, 'NumberOfBands', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) kPoints%numBands
      if(error /= 0) stop 'Error, couldn''t read NumberOfBands.'

      print *, 'NumberOfBands ', kPoints%numBands

    else

      stop 'Error, could not find NumberOfBands.'

    endif

!   Input parameter NumberOfOccupiedBands for BerkeleyGW wavefunction file

    call FindTag( tagHandler, 'NumberOfOccupiedBands', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) kPoints%numOccupiedBands
      if(error /= 0) stop 'Error, couldn''t read NumberOfOccupiedBands.'

      print *, 'NumberOfOccupiedBands ', kPoints%numOccupiedBands

    else

      kPoints%numOccupiedBands = kPoints%numBands

    endif

! Determine whether kpoints are read from input file are in units of ReciprocalLatticeVectors
! or in Cartesian units.  Default is scaledByRecipLatVec
 
    default = .True.

    call FindTag(tagHandler, 'KPointsScale', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(A80)', iostat = error) buff

      if(error .eq. 0) then

        if(StringCompare(buff, 'ReciprocalLatticeVectors')) then

          kPointsScale = scaledByRecipLatVec
          print *, 'KPointsScale ReciprocalLatticeVectors'
          default = .False.

        else if(StringCompare(buff, 'TwoPi/a')) then

          kPointsScale = scaledByTwoPi
          print *, 'KPointsScale TwoPi/a'
          default = .False.

        end if
      end if
    end if

    if(default) then
      kPointsScale = scaledByRecipLatVec
      print *, 'Using default KPointsScale ReciprocalLatticeVectors'
    end if

!   Start Various Possible types of generation of k ponts

!   Input parameter KPointNumbers for BerkeleyGW wavefunction file
    call FindTag( tagHandler, 'KPointNumbers', error)
    if(error .eq. 0) then
      read(unit = tagHandler%fileno, fmt = *, iostat = error) (kPoints%KPointNumbers(j),j=1,3)
    else
      do j=1,3
        kPoints%KPointNumbers(j)=0
      enddo
    endif

!   Input parameter KPointOffsets for BerkeleyGW wavefunction file
    call FindTag( tagHandler, 'KPointOffsets', error)
    if(error .eq. 0) then
      read(unit = tagHandler%fileno, fmt = *, iostat = error) (kPoints%KPointOffsets(j),j=1,3)
    else
      do j=1,3
        kPoints%KPointOffsets(j)=0.0d0
      enddo
    endif

    call FindTag( tagHandler, 'KPointsList', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) kPoints%numKPoints
      if(error /= 0) stop 'Error: Could not read number of KPoints after KPointsList.'

      print *, 'Use k point list read from input file '
      print *, 'Calculation will be done for ' , kPoints%numKPoints, ' k points' 

      allocate(kPoints%kPointsCart(ndim, kPoints%numKPoints), &
        & kPoints%kPointsLat(ndim, kPoints%numKPoints), &
        & kPoints%kPointWeights(kPoints%numKPoints), &
        & ktemp(ndim, kPoints%numKPoints), &
        & stat=error) ! allocate memory for k-points, and temporary ktemp 
      if(error /= 0) stop 'could not allocate memory for k-points'

!     Needed to read data from following line
!     read(unit = tagHandler%fileno,fmt = '(A1)', iostat = error) c

!     Input parameter kPointWeights for BerkeleyGW wavefunction file
      do j = 1, kPoints%numKPoints
        read(unit = tagHandler%fileno, fmt = *, iostat = error) ktemp(:,j), kPoints%kPointWeights(j)
        if(error /= 0) stop 'Error reading KPoints.'
        print *, j, ' ', ktemp(:,j)
      enddo

      normKPointWeights = 0.0d0
      do j = 1, kPoints%numKPoints
        normKPointWeights = normKPointWeights + kPoints%kPointWeights(j)
      enddo
      if (abs(normKPointWeights).gt.1.0d-6) then
        do j = 1, kPoints%numKPoints
          kPoints%kPointWeights(j) = kPoints%kPointWeights(j) / normKPointWeights
        enddo
      else
        stop "k-weights sum is close to zero."
      endif

      print *, ' HERE WE ASSUME INPUT K ARE IN LATTICE VECTOR UNITS - NEED TO CHANGE'

      do j = 1, kPoints%numKPoints
        kPoints%kPointsLat(:,j) = ktemp(:,j)
        print *, j, ' ', ktemp(:,j)
      enddo

      deallocate(ktemp, stat = error)
      if(error /= 0) stop 'Error deallocating ktemp in KPointsInit.'

    else

      !Checking for next option - read lines and divisions for plotting

      call FindTag( tagHandler, 'NumberOfLines', error)

      if(error .eq. 0) then

      call KPointsLinesandGraphInit( kPoints, graphInfo, ndim, lattice, tagHandler )

      else
        if(error /= 0) stop 'Error, could not find any tag for kpoints'
      endif

    endif

    POP_SUB(KPointsInit)
   end Subroutine KPointsInit

!*************************************************************************
  Subroutine KPointsDestroy( kPoints )
!*************************************************************************

    use GraphMod, only : GraphDestroy

    implicit none

    type( kPointsT ), pointer :: kPoints
    integer :: error

    PUSH_SUB(KPointsDestroy)

    if( .not. associated( kPoints )) stop 'Error, kPoints not allocated.'

!    if( associated( kPoints%graph ))  call GraphDestroy( kPoints%graph )

!    deallocate(kPoints%kPointsLat, kPoints%kPointsCart,kPoints%xAxis, stat=error)
    deallocate(kPoints%kPointsLat, kPoints%kPointsCart, kPoints%kPointWeights, stat=error)
    if(error /= 0) stop 'Error deallocating kPointsLat in KPointsDestroy.'

    deallocate( kPoints, stat = error )
    if(error /= 0) stop 'Error deallocating kPoints in KPointsDestroy.'

    nullify( kPoints )

    POP_SUB(KPointsDestroy)
  end Subroutine KPointsDestroy

!*************************************************************************
  Subroutine KPointsLinesandGraphInit( kPoints, graphInfo, ndim, lattice, tagHandler )
!*************************************************************************

!    generate k-points around a BZ for plotting band structure plots such that
!    each line segment has a constant density of points to make the graphics
!    look good

    use SysParams,        only : one, zero, twopi
    use TagHandlerMod,    only : FindTag, StringCompare, tagHandlerT
    use LatticeMod,       only : LatticeT
    use GraphMod,         only : AddLabel, GraphInit, graphT

    implicit none

    type( LatticeT ),    pointer :: lattice
    type( kPointsT ),    pointer :: kPoints
    type( tagHandlerT ), pointer :: tagHandler
    type( graphT ), pointer      :: graphInfo
    integer, intent(in)  :: ndim

    integer                   :: nlines       ! number of lines in k-space to plot
    integer                   :: ndiv         ! number of divisions per line
    integer,      allocatable :: num(:)       ! number of divisions per line
    real(double), allocatable :: length(:)    ! length of line
    real(double), allocatable :: kc(:,:)      ! list of circuit k-points
    real(double), allocatable :: tick(:)      ! list of positions of lables
    character(2), allocatable :: label(:)     ! list of labels for plot
    real(double), allocatable :: dk(:,:)      ! line vector
    character(2)  :: c
    character(80) :: buff
    integer       :: i,j,npt,error, kPointsScale
    real(double)  :: delta_k
    logical       :: default

    PUSH_SUB(KPointsLinesandGraphInit)

    call FindTag( tagHandler, 'NumberOfLines', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) nlines
      if(error /= 0) stop 'Error, couldn''t read NumberOfLines.'

      print *, 'NumberOfLines ', nlines

    else

      stop 'Error, could not find NumberOfLines.'

    endif

    call FindTag( tagHandler, 'NumberOfDivisions', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) ndiv
      if(error /= 0) stop 'Error, couldn''t read NumberOfDivisions.'

      print *, 'NumberOfDivisions ', ndiv

    else

      stop 'Error, could not find NumberOfDivisions.'

    endif

    allocate(kc(ndim, nlines + 1), dk(ndim, nlines), &
      & label(nlines + 1), tick(nlines + 1), &
      & num(nlines), length(nlines), stat = error) 
    if(error /= 0) stop 'Error allocating k-points space in KPointsInit.'

! Determine whether kpoints are read from input file are in units of ReciprocalLatticeVectors
! or in Cartesian units.  Default is scaledByRecipLatVec
 
    default = .True.

    call FindTag(tagHandler, 'KPointsScale', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(A80)', iostat = error) buff

      if(error .eq. 0) then

        if(StringCompare(buff, 'ReciprocalLatticeVectors')) then

          kPointsScale = scaledByRecipLatVec
          print *, 'KPointsScale ReciprocalLatticeVectors'
          default = .False.

        else if(StringCompare(buff, 'TwoPi/a')) then

          kPointsScale = scaledByTwoPi
          print *, 'KPointsScale TwoPi/a'
          default = .False.

        end if
      end if
    end if

    if(default) then
      kPointsScale = scaledByRecipLatVec
      print *, 'Using default KPointsScale ReciprocalLatticeVectors'
    end if

    call FindTag( tagHandler,'KPointsAndLabels', error)
    print *, 'KPointsAndLabels'
    if(error .eq. 0) then

      read(unit = tagHandler%fileno,fmt = '(A1)', iostat = error) c

      read(unit = tagHandler%fileno,fmt =*, iostat = error) kc(:,1), label(1)

      print *, kc(:,1),' ', label(1)
!      call AddLabel( graphInfo, label, zero )
      tick(1) = 0.0d0

      do i = 1, nlines

        ! input circuit information
        read(unit = tagHandler%fileno, fmt = *, iostat = error) kc(:,i+1),label(i+1)
        dk(:,i)    = kc(:,i+1) - kc(:,i)
        length(i)  = SQRT(dot_product(matmul(dk(:,i),lattice%bLatMetricTensor),dk(:,i)))
        tick(i+1) = tick(i) + length(i)
!        call AddLabel( graphInfo, label(i+1), tick(i+1) )
        print *, kc(:,i+1),' ', label(i+1)

      end do

    else

      stop 'Error, could not find KPointsAndLabels.'

    end if

    num = (length/minval(length))*ndiv          ! number of divisions per line

    kPoints%numKPoints  = SUM(num) + 1          ! total number of k-points

    allocate(kPoints%kPointsCart(ndim, kPoints%numKPoints), &
        & kPoints%kPointsLat(ndim, kPoints%numKPoints), &
        & kPoints%kPointWeights(kPoints%numKPoints), &
        stat=error) ! allocate memory for k-points
    if(error /= 0) stop 'Error, could not allocate memory for k-points.'

! Generate a continuous set of k-points around the circuit 

    npt = 0
    do i = 1, nlines                             
      delta_k     = one/num(i)                
      npt         = npt + 1
      kPoints%kPointsLat(:,npt) = kc(:,i)
      do j = 1, num(i) - 1
        npt         = npt + 1
        kPoints%kPointsLat(:,npt) = kc(:,i) + j * delta_k * dk(:,i)
      end do
    end do

    kPoints%kPointsLat(:, npt+1) = kc(:, nlines+1) ! get the last k-point

! Generate information for plotting bands

   call GraphInit(graphInfo, kpoints%numKpoints, nlines)

      do i = 1, nlines + 1
        call AddLabel( graphInfo, label(i), tick(i) )
      end do

    graphInfo%xAxis = 0

    do i = 2, kPoints%numKPoints

      dk(:,1)        = kPoints%kPointsLat(:,i) - kPoints%kPointsLat(:,i-1)
      graphInfo%xAxis(i) = graphInfo%xAxis(i-1) + &
          & sqrt(dot_product(dk(:,1),matmul(lattice%bLatMetricTensor,dk(:,1))))

    end do

! Calculate k points in Cartesian coordinates for possible use

    if(kPointsScale .eq. scaledByRecipLatVec) then
      kPoints%kPointsCart = matmul(lattice%bLatVec, kPoints%kPointsLat)  &
                          & * lattice%latticeConst / twopi
    else if(kPointsScale .eq. scaledByTwoPi) then
      kPoints%kPointsCart = kPoints%kPointsLat
      kPoints%kPointsLat = matmul(transpose(lattice%aLatVec), kPoints%kPointsCart) / & 
                 & lattice%latticeConst * twopi
    end if

    deallocate(kc, dk, label, tick, num, length, stat = error)
    if(error /= 0) stop 'Error deallocating kc, dk, label, tick, num, length in KPointsInit.'

    POP_SUB(KPointsLinesandGraphInit)
  end Subroutine KPointsLinesandGraphInit

!*************************************************************************
  Subroutine PrintKPoints( kPoints )
!*************************************************************************

    implicit none

    type( kPointsT ), pointer :: kPoints
    integer :: i

    PUSH_SUB(PrintKPoints)

    do i = 1, kPoints%numKPoints
      print *, kPoints%kPointsLat(:,i)
    end do

    POP_SUB(PrintKPoints)
  end Subroutine PrintKPoints

end Module KPointsMod

