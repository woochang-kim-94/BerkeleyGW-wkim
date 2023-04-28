#include "f_defs.h"

!*************************************************************************
Module pwHamMod
!*************************************************************************
! The generic name of the module "SpecificHamMod" and the subroutines
! "HGenerate", "HamInfoInit", "HamInfoDestroy",  must be used so they
! are called properly by the general hamiltonian routines
!*************************************************************************
! This module contains all the routines needs for seeting up the
! Plane wave hamiltonian h (a 1d array which is upper triangular part
! of the H matrix of size = hSize)  Size is fixed by energy cutoff
! read from input file in the GVectors module
!*************************************************************************
! Function atomPotential(q2,elementLabel) is called here - It is
! in module AtomPotentialMod - it is the only subroutine
! needs to be changed for different potentials 
!*************************************************************************
! This module uses the following modules:
!    SysParams
!    GVectorsMod
!    Hamiltonian
!    kPointsMod
!    LatticeMod (only for bLatMetricTensor used to find lengths of k + G)
!*************************************************************************
! The module contains the public functions:
!    HGenerate       Generates the plane wave hamiltonian for
!                    the given k point K_Points$kPointsLat(:,k).
!    HamInfoInit     Sets up relevant information for the hamiltonian
!    HamInfoDestroy
!    GridInit
!    GridDestroy
!*************************************************************************
! The module contains the private functions:
!   PotentialInit    Allocates space for initializes potentialV.
!   PotentialDestroy Deallocates the space for potentialV
!*************************************************************************
! The module contains the private variables:
!   potentialV(i)    Fourier component of the crystal (pseudo)potential
!                    for G-vector gVectors(i,:)
!*************************************************************************

  use SysParams,        only : double
  use GVectorsMod,      only : gVectorObjectsT
  Use EigenStatesMod,   only : eigenStatesT
  use typeMod,          only : GridT
  use fftw_m
  use push_pop_m

  implicit none

  save

  type hamInfoT
    integer                    :: hMaxSize
    complex(double),   pointer :: potentialV(:)
    type( gVectorObjectsT ), pointer :: gVectors
  end type hamInfoT
  
  private PotentialInit
  private PotentialDestroy

  contains

!*************************************************************************
  Subroutine HamInfoInit( hamInfo, struct, tagHandler, v_of_q, ff )
!*************************************************************************

    use SysParams,     only : double
    use GVectorsMod,   only : GVectorsInit
    use StructureMod,  only : structureT
    use TagHandlerMod, only : tagHandlerT, FindTag

    implicit none

    type( hamInfoT ),    pointer :: hamInfo
    type( structureT ),  pointer :: struct
    type( tagHandlerT ), pointer :: tagHandler
    logical :: v_of_q
    real ( double ) :: ff(6)

    integer :: error

    PUSH_SUB(HamInfoInit)

    if( associated( hamInfo )) stop 'Error: SpecificHamiltonian already allocated.'
    allocate( hamInfo, stat = error )
    if( error /= 0 ) stop 'Error: SpecificHamiltonian allocation failed.'
    nullify( hamInfo%potentialV )
    nullify( hamInfo%gVectors )

    call GVectorsInit( hamInfo%gVectors, struct, tagHandler, v_of_q )
    call PotentialInit( hamInfo, struct, v_of_q, ff )

    POP_SUB(HamInfoInit)
  end Subroutine HamInfoInit

!*************************************************************************
  Subroutine HamInfoDestroy( hamInfo )
!*************************************************************************

    use GVectorsMod, only : GVectorsDestroy
    implicit none

    type( hamInfoT ), pointer :: hamInfo

    integer :: error

    PUSH_SUB(HamInfoDestroy)
    if( .not. associated( hamInfo )) stop 'Error: SpecificHamiltonian not allocated.'

    call GVectorsDestroy( hamInfo%gVectors )
    call PotentialDestroy( hamInfo )

    deallocate( hamInfo, stat = error )
    if( error /= 0 ) stop 'Error: specificHamiltonian deallocation failed.'

    POP_SUB(HamInfoDestroy)
  end Subroutine HamInfoDestroy

!*************************************************************************
  Subroutine PotentialInit( hamInfo, struct, v_of_q, ff )
!*************************************************************************

!  For all local potentials the Hamiltonian matrix depends on k only through
!  the diagonal terms. The Fourier components of the potential terms 
!  are independent of k and can be set up only once for all k.

    use SysParams,    only : two, twopi, half, double
    use StructureMod, only : StructureT
    use AtomPotentialMod

    implicit none

    real( double ), parameter :: eps6 = 1.0d-6

    type( hamInfoT ),   pointer :: hamInfo
    type( StructureT ), pointer :: struct
    logical :: v_of_q
    real( double ) :: ff(6)

    integer :: error, i, j, i2
    real( double ) :: q2, ffr, ffi, sfr, sfi

    PUSH_SUB(PotentialInit)

    allocate( hamInfo%potentialV(hamInfo%gVectors%numGVectors), stat=error )
    if( error /= 0 ) stop 'Error allocating potential space.'
                
    hamInfo%potentialV(:) = 0

    if (v_of_q) then

      do j = 1, struct%atomBasis%numSpecies
        do i = 1, hamInfo%gVectors%numGVectors

          ! Compute potential at G-vectors
          hamInfo%potentialV(i) = hamInfo%potentialV(i) + &
                  & atomPotential(two * hamInfo%gVectors%kineticE(i), struct%atomBasis%elementLabel(j)) * &
                  & hamInfo%gVectors%structureFactor(i,j)

        end do
      end do

    else

      do i = 1, hamInfo%gVectors%numGVectors

        q2 = two * hamInfo%gVectors%kineticE(i) &
               & * (struct%lattice%latticeConst / twopi)**2
        i2 = int(q2 + eps6)
        if (abs(q2 - dble(i2)) .gt. eps6) &
        stop 'Error, q2 is not integer.'

        select case (i2)
        case(3)
          ffr = ff(1)
          ffi = ff(4)
        case(4)
          ffr = 0.0d0
          ffi = ff(5)
        case(8)
          ffr = ff(2)
          ffi = 0.0d0
        case(11)
          ffr = ff(3)
          ffi = ff(6)
        case default
          ffr = 0.0d0
          ffi = 0.0d0
        end select
        ffr = ffr * half
        ffi = ffi * half

        sfr = dble(hamInfo%gVectors%structureFactor(i,1))
        sfi = -aimag(hamInfo%gVectors%structureFactor(i,1))

        hamInfo%potentialV(i) = hamInfo%potentialV(i) + &
                & cmplx(ffr * sfr, ffi * sfi, kind = double)

      end do

    endif

    ! Explicitly put the potential at the origin to be zero.
    hamInfo%PotentialV(hamInfo%gVectors%gVecIndex(0,0,0)) = cmplx(0.d0,0.d0,kind=double)

    POP_SUB(PotentialInit)
  end Subroutine PotentialInit

!*************************************************************************
  Subroutine PotentialDestroy( hamInfo )
!*************************************************************************

    implicit none

    type( hamInfoT ), pointer :: hamInfo
    integer :: error

    PUSH_SUB(PotentialDestroy)

    deallocate( hamInfo%potentialV, stat=error )
    if(error /= 0) stop 'Error deallocating potential space.'

    POP_SUB(PotentialDestroy)
  end Subroutine PotentialDestroy


!*************************************************************************
  Subroutine FindhMaxSize( hamInfo, struct, kPoints )
!*************************************************************************

!   This subroutine goes through the same loops as for generating the Hamiltonian
!   in order to find the maximum memory needed 
!
!   Output is hamInfo%hMaxSize which is defined in module variables and is public

    use SysParams,     only : double, two
    use StructureMod,  only : StructureT
    use kPointsMod,    only : kPointsT

    implicit none

    type( StructureT ),   pointer :: struct
    type( kPointsT ),     pointer :: kPoints
    type( hamInfoT ),     pointer :: hamInfo
    integer                        :: i
    integer                        :: k, hSize, ktemp 
    real(double)                   :: ek
    real(double), dimension(3)     :: kPlusG

    PUSH_SUB(FindhMaxSize)

    hamInfo%hMaxSize = 0

    do k = 1, kPoints%numKPoints

    hSize = 0
    kPlusG = 0

! This section is
! inefficient since their are a lot of G vectors in gVectors that
! will never be within the energy cut off radius, because all of the
! G vectors within two times the cut off radius are stored in gVectors
! for the structure factor calculation.

       do i = 1, hamInfo%gVectors%numGVectors

         kPlusG(1:struct%ndim) = kPoints%kPointsLat(:,k) + hamInfo%gVectors%gVectors(:,i)
         ek = dot_product(kPlusG(1:struct%ndim), &
             & matmul(struct%lattice%bLatMetricTensor, kPlusG(1:struct%ndim)))/two

         if(ek <= hamInfo%gVectors%energyCutOff) then
           hSize        = hSize + 1
         end if

       end do


       if(hSize < kPoints%numBands) then
         print *, 'Error number of bands = ', kPoints%numBands, &
             & ' > size of the hamiltonian = ', hSize, &
             & ' for kPoint = ', k
         stop
       end if


       if(hSize > hamInfo%hMaxSize) then
         hamInfo%hMaxSize = hSize ! hamInfo%hMaxSize set to largest hSize
         ktemp = k
       end if 

    end do

!  print *, 'Maximum size of the hamiltonian = ', hamInfo%hMaxSize, &
!                & ' for k point = ', ktemp, numKPoints
!  print *, 'In recip. lattice basis vectors, this k point is:'
!  print *, kPointsLat(:,ktemp)

    POP_SUB(FindhMaxSize)
  end Subroutine FindhMaxSize

!*************************************************************************
  Subroutine KineticGenerate(hamInfo,struct,kPoints,eigenStates,k,Kinetic, & 
       & hSize)
!*************************************************************************

!   This subroutine will fill up the kinetic energy T_{G} array of the 
!   hamiltonian. T_{G} = |k+G|^2/2, and the range of G is such that the 
!   energy cutoff is satisfied.

    use EigenStatesMod,   only : eigenStatesT
    use SysParams,        only : double, two
    use StructureMod,     only : StructureT
    use kPointsMod,       only : kPointsT

    Implicit None

    type( StructureT ),   pointer :: struct
    type( kPointsT ),     pointer :: kPoints
    type( eigenStatesT ), pointer :: eigenStates
    type( hamInfoT ),     pointer :: hamInfo

    integer, intent(in)           :: k
    integer                       :: hSize,i
    real(double)                  :: ek
    real(double), dimension(3)    :: kPlusG
    real(double), dimension(:)    :: Kinetic

    PUSH_SUB(KineticGenerate)

    hSize = 0
    kPlusG = 0

    do i = 1, hamInfo%gVectors%numGVectors
       kPlusG(1:struct%ndim) = kPoints%kPointsLat(:,k) + & 
            & hamInfo%gVectors%gVectors(:,i) 
       ek = dot_product(kPlusG(1:struct%ndim), &
            & matmul(struct%lattice%bLatMetricTensor, & 
            & kPlusG(1:struct%ndim)))/two

       if(ek <= hamInfo%gVectors%energyCutOff) then
          hSize        = hSize + 1
          Kinetic(hSize) = ek
          eigenStates%basisIndex(hSize,k) = i

       end if
    end do

    POP_SUB(KineticGenerate)
  end Subroutine KineticGenerate

!*************************************************************************
  Subroutine PotentialRGenerate(hamInfo,Grid)
!*************************************************************************

!   This subroutine will put the Potential_{G} array on the grid, the backward
!   fourier transform into V_{r}

    Type(hamInfoT),           pointer :: hamInfo
    Type(GridT),              pointer :: Grid
    integer                           :: i,j,k,s,ii

    PUSH_SUB(PotentialRGenerate)

    ! size of the VGrid has to be up to maxGLen
    Grid%V = cmplx(0.d0,0.d0,kind=double)

    ! Translation
    do i = 1, 2*Grid%Size(1)
       do j = 1, 2*Grid%Size(2)
          do k = 1, 2*Grid%Size(3)
             ii = Grid%GInvIndex(i-1-Grid%Size(1), &
                  & j-1-Grid%Size(2),k-1-Grid%Size(3))
             if (ii.gt.0) Grid%V(i,j,k) = hamInfo%PotentialV(ii)
          end do
       end do
    end do

    ! This is so that the final grid VPsi_G is centered at the origin.
    do i = 1, 2*Grid%Size(1)
       do j = 1, 2*Grid%Size(2)
          do k = 1, 2*Grid%Size(3)
             s = i+j+k-3
             Grid%V(i,j,k) = Grid%V(i,j,k)*(-1)**s
          end do
       end do
    end do

    POP_SUB(PotentialRGenerate)
  end Subroutine PotentialRGenerate

!*************************************************************************
  Subroutine HGenerate( hamInfo, struct, kPoints, eigenStates, &
                      & k, h, sqH, hSize )
!*************************************************************************

    use EigenStatesMod,   only : eigenStatesT
    use SysParams,        only : double, two
    use StructureMod,     only : StructureT
    use kPointsMod,       only : kPointsT

    implicit none

    type( StructureT ),   pointer :: struct
    type( kPointsT ),     pointer :: kPoints
    type( eigenStatesT ), pointer :: eigenStates
    type( hamInfoT ),     pointer :: hamInfo
    integer, intent(in)            :: k
    complex(double), intent(inout) :: h(:)     
    complex(double), intent(inout) :: SqH(:,:)     
    integer, intent(inout)         :: hSize

    integer, dimension(3)          :: gg       
    integer                        :: i,j,indx
    real(double)                   :: ek
    real(double), dimension(3)     :: kPlusG

    !!!!!!!!!!!!!!!! CHANGED ddas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: N

    PUSH_SUB(HGenerate)

    N = 180
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    hSize = 0
    kPlusG = 0

    ! For each G vector calculate the kinetic energy if it is less than the 
    ! cut off energy add a row for it on the hamiltonian insert the kinetic
    ! energy along the diagonal.  BasisIndex is used to store the indices 
    ! of the G vectors that are in the hamiltonian.  This section is
    ! inefficient since their are a lot of G vectors in gVectors that
    ! will never be within the energy cut off radius, because all of the
    ! G vectors within two times the cut off radius are stored in gVectors
    ! for the structure factor calculation.

    do i = 1, hamInfo%gVectors%numGVectors

      kPlusG(1:struct%ndim) = kPoints%kPointsLat(:,k) + hamInfo%gVectors%gVectors(:,i) 
      ek = dot_product(kPlusG(1:struct%ndim), &
        & matmul(struct%lattice%bLatMetricTensor, kPlusG(1:struct%ndim)))/two

      if(ek <= hamInfo%gVectors%energyCutOff) then

        hSize        = hSize + 1
        h(hSize * (hSize + 1) / 2) = ek

        SqH(hSize,hSize) = ek                       ! For the square matrix

        eigenStates%basisIndex(hSize,k) = i

      end if
    end do

!!$    !!!!!!!!!!!!!!! CHANGED ddas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    do i = 1, N
!!$       SqH(i,i) = cmplx(1.d0,0.d0,kind=double)    
!!$    end do
!!$
!!$    do i = 1, N
!!$       do j = i+1, N
!!$          SqH(i,j) = cmplx(0.d0,0.d0,kind=double)
!!$          SqH(j,i) = conjg(SqH(i,j))
!!$       end do
!!$    end do
!!$    !!!!!!!!!!!!!!! CHANGED ddas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    eigenStates%numBasisVectors(k) = hSize

    if(hSize < kPoints%numBands) then
      print *, 'Error number of bands = ', kPoints%numBands, &
             & ' > size of the hamiltonian = ', hSize
      stop
    end if

    ! For each row in the hamiltonian store the kinetic energy along the
    ! diagonal and calculate the off diagonal terms.
    do i = 1, hSize

      do j = i+1, hSize

        gg(1:struct%ndim) = &
      & hamInfo%gVectors%gVectors(1:struct%ndim, eigenStates%basisIndex(i, k))&
    & - hamInfo%gVectors%gVectors(1:struct%ndim, eigenStates%basisIndex(j, k))
        if(struct%ndim < 3) gg(struct%ndim+1:3) = 0
        indx = hamInfo%gVectors%gVecIndex(gg(1),gg(2),gg(3))    ! index for G-G` vector
        h(i+j*(j-1)/2) = hamInfo%potentialV(indx)

        SqH(i,j) = hamInfo%potentialV(indx)
        SqH(j,i) = conjg(SqH(i,j))

      end do
    end do

    POP_SUB(HGenerate)
  end Subroutine HGenerate

!*************************************************************************
  Subroutine CalculateDensity( hamInfo, kPoints, eigenStates, densityArray )
!*************************************************************************

    use EigenStatesMod,  only : eigenStatesT
    use densityArrayMod, only : densityArrayT
    use StructureMod,    only : StructureT
    use kPointsMod,      only : kPointsT

    implicit none

    type( kPointsT ),       pointer :: kPoints
    type( densityArrayT ),  pointer :: densityArray
    type( eigenStatesT ),   pointer :: eigenStates
    type( hamInfoT ),       pointer :: hamInfo
    integer :: n, k, ig

    integer, allocatable :: identity(:)

    PUSH_SUB(CalculateDensity)

    densityArray%densityArray(:,:,:) = 0

    do k = 1, kPoints%numKPoints
      do n = 1, kPoints%numOccupiedBands
        call put_into_fftbox(eigenStates%numBasisVectors(k), eigenStates%eigenVectors(:, n, k), hamInfo%gVectors%gVectors, &
          eigenStates%basisIndex(:, k), densityArray%waveFunction, densityArray%gridDimensions)
        call do_FFT(densityArray%waveFunction, densityArray%gridDimensions, -1)
        ! normalization
        densityArray%waveFunction = densityArray%waveFunction / sqrt(dble(densityArray%gridSize))
        densityArray%densityArray = densityArray%densityArray + &
          & densityArray%waveFunction*conjg(densityArray%waveFunction) &
          & *kPoints%kPointWeights(k)
      end do
    end do

    call do_FFT(densityArray%densityArray, densityArray%gridDimensions, +1)
    allocate(identity(hamInfo%gVectors%numGVectors))
    do ig = 1, hamInfo%gVectors%numGVectors
      identity(ig) = ig
    enddo
    call get_from_fftbox(hamInfo%gVectors%numGVectors, densityArray%rhog, hamInfo%gVectors%gVectors, identity, &
      densityArray%densityArray, densityArray%gridDimensions, 1d0)
    deallocate(identity)
!    write(0,*) 'RHO G=0 : ', densityArray%rhog(hamInfo%gVectors%gVecIndex(0, 0, 0)), &
!      hamInfo%gVectors%gVectors(:, hamInfo%gVectors%gVecIndex(0, 0, 0))

    call destroy_fftw_plans()

    POP_SUB(CalculateDensity)
  end Subroutine CalculateDensity

!*************************************************************************
  Subroutine GridInit(Grid,hamInfo,eigenStates,ndim)
!*************************************************************************

    type(GridT),                               pointer :: Grid
    type(hamInfoT),                            pointer :: hamInfo
    Type(eigenStatesT),                        pointer :: eigenStates    
    integer,                                   intent(in) :: ndim      

    PUSH_SUB(GridInit)

    !  Allocate space for the Grid Type
    allocate(Grid)
    nullify(Grid%V)
    nullify(Grid%Psi)
    nullify(Grid%VPsi)
    nullify(Grid%VPsiG)
    nullify(Grid%Work)
    nullify(Grid%GVecIndex)
    nullify(Grid%GInvIndex)
    nullify(Grid%BasisIndex)

    !  Length of the Grid in each direction
    Grid%Size = hamInfo%gVectors%maxGLen

    Grid%fac = dble(8*Grid%Size(1)*Grid%Size(2)*Grid%Size(3))

    !  Number of GVectors
    Grid%NumGVec = hamInfo%gVectors%numGVectors

    !  Allocate space for the VGrid, and initialise it to zero
    allocate(Grid%V(2*Grid%Size(1),2*Grid%Size(2),2*Grid%Size(3)))
    allocate(Grid%Psi(2*Grid%Size(1),2*Grid%Size(2),2*Grid%Size(3)))
    allocate(Grid%VPsi(2*Grid%Size(1),2*Grid%Size(2),2*Grid%Size(3)))

    allocate(Grid%VPsiG(-Grid%Size(1):Grid%Size(1),-Grid%Size(2):Grid%Size(2),&
         &  -Grid%Size(3):Grid%Size(3)))


    Grid%V   = cmplx(0.d0,0.d0,kind=double)

    !  Allocate and Initialise the Grid Index
    allocate(Grid%GVecIndex(ndim,Grid%NumGVec))
    Grid%GVecIndex = hamInfo%gVectors%gVectors

    allocate(Grid%GInvIndex( &
         & -hamInfo%gVectors%maxGLen(1):hamInfo%gVectors%maxGLen(1), &
         & -hamInfo%gVectors%maxGLen(2):hamInfo%gVectors%maxGLen(2), &
         & -hamInfo%gVectors%maxGLen(3):hamInfo%gVectors%maxGLen(3)))
    Grid%GInvIndex = hamInfo%gVectors%gVecIndex

    allocate(Grid%BasisIndex(EigenStates%maxBasisVectors, & 
         & EigenStates%numKPoints))

    POP_SUB(GridInit)
  end Subroutine GridInit

!*************************************************************************
  Subroutine GridDestroy(Grid)
!*************************************************************************
! by DAS. original developers did not see fit to include such a thing.

    type(GridT), pointer :: Grid !< intent not permitted for pointer

    integer :: error

    PUSH_SUB(GridDestroy)

    deallocate(Grid%V, Grid%Psi, Grid%VPsi, Grid%VPsiG, Grid%GVecIndex, &
         & Grid%GInvIndex, Grid%BasisIndex, stat = error)
    if( error /= 0 ) stop 'Error deallocating Grid space.'
    deallocate(Grid, stat = error)
    if( error /= 0 ) stop 'Error deallocating Grid.'

    POP_SUB(GridDestroy)
  end Subroutine GridDestroy

end Module pwHamMod

