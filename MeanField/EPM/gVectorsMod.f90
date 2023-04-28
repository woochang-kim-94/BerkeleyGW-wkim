#include "f_defs.h"

!*************************************************************************
Module GVectorsMod
!*************************************************************************
! Information related to the reciprocal lattice vectors (G-vectors).
!*************************************************************************
! The module contains the public functions:
!    GVectorsInit         reads in EnergyCutOff and allocates space for gVectors,
!                         gVecIndex, and kineticE. It initializes gVectors,
!                         and gVecIndex. Finally it calls StructureFactorInit.
!    GVectorsDestroy      deallocates space for gVectors, gVecIndex, and kineticE.
!*************************************************************************
! The module contains the public variables:
!    energyCutOff         plane-wave energy cutoff (a.u.) for wavefunctions
!    numGVectors          total number of G-vectors  
!    gVectors(i,j)        i-th component (reciprocal lattice coordinates)
!                         of the j-th G-vector
!    gVecIndex(i,j,k)     index for the G-vectors: 
!                         gVectors(gVecIndex(i,j,k),:) is the G-vector
!                         with reciprocal lattice coordinates (i,j,k)
!    kineticE(i)          kineticE(gVecIndex(i,j,k)) is the kinetic energy
!                         of the G-vector with reciprocal lattice
!                         coordinates (i,j,k)
!    structureFactor(i,j) complex structure factor for atom type j and
!                         G-vector i
!*************************************************************************
! The module contains the private functions:
!    StructureFactorInit  allocates space for structureFactor and then
!                         initializes structureFactor.
!*************************************************************************

  use SysParams, only : double
  use push_pop_m
  implicit none

  save

  type gVectorObjectsT
    integer                  :: numGVectors     
    integer, pointer         :: gVectors(:,:)   
    integer, pointer         :: gVecIndex(:,:,:) 
    real(double)             :: energyCutOff    
    integer, dimension(3)    :: maxGLen      
    real(double), pointer    :: kineticE(:)     
    complex(double), pointer :: structureFactor(:,:)
  end type gVectorObjectsT

  private StructureFactorInit
    
  contains

!*************************************************************************
  Subroutine GVectorsInit( gVectors, structure1, tagHandler, v_of_q )
!*************************************************************************

! GVectorInit reads in the "EnergyCutOff" and creates the G-vectors lying
! in a sphere of radius sqrt(two * maxEnergy) where maxEnergy is the
! cutoff for the potential which is 4 * EnergyCutOff for the wavefunctions.
! It creates the corresponding index gVecIndex, and 
! calls StructureFactorinit.

! Defined Tags:
!     EnergyCutOff

    use SysParams,       only : four, two, twopi, electronMass, &
                              & hBar2, hartree2eV
    use TagHandlerMod,   only : FindTag, tagHandlerT
    use StructureMod,    only : StructureT
    
    implicit none

    type( gVectorObjectsT ), pointer :: gVectors
    type( StructureT ), pointer :: structure1
    type( tagHandlerT ), pointer :: tagHandler
    logical :: v_of_q

    real(double)              :: kinetic, maxEnergy, maxGradius
    integer                   :: i, j, k, maxNumGVec, error
    integer,      allocatable :: gVectors_temp(:,:)
    integer,      dimension(3):: g_trial 
    real(double), allocatable :: kineticE_temp(:)

    PUSH_SUB(GVectorsInit)

    allocate( gVectors, stat = error )
    if( error /= 0 ) stop 'Error: gVectors allocation failed.'

    !**********************************************************************
    ! Read in EnergyCutOff
    call FindTag( tagHandler,'EnergyCutOff', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = *, iostat = error) gVectors%energyCutOff
      if(error /= 0) stop 'Error, couldn''t read EnergyCutOff.'
    else
      stop 'Error, could not find EnergyCutOff.'
    endif

    ! Energies read in Rydberg must be converted to Hartree
    gVectors%energyCutoff = gVectors%energyCutoff * 0.5d0

    ! Energies read in eV must be converted to Hartree
    call FindTag( tagHandler, 'InputEnergiesInEV', error)
    if(error .eq. 0) then
       print *, 'InputEnergiesInEV'
       gVectors%energyCutoff = gVectors%energyCutoff / hartree2eV
    else
       print *, 'EnergyCutOff in Rydberg ', gVectors%energyCutOff * 2.0d0
    endif

    ! energyCutOff has been initialized.
    !**********************************************************************

    !**********************************************************************
    ! Initialize gVectors, kineticE and gVecIndex.
    gVectors%numGVectors = 0
    gVectors%maxGLen     = 0
    
    ! The maximum energy for a G vector is four times the cut off energy.
    maxEnergy   = four * gVectors%energyCutOff 
    maxGradius =  sqrt(two * maxEnergy)

    ! The maximum length maxGLen(i) of a G vector along the ith reciprocal
    ! lattice vector bLatVec(i).  Think about it.
    gVectors%maxGLen(1:structure1%ndim) = int(maxGradius * structure1%lattice%aLatVecLen / twopi) + 1

    ! The number of G vectors in the parallelepiped constructed with
    ! int(2*maxGLen(i)) + 1 cells in each direction.
    maxNumGVec              = product(2 * gVectors%maxGLen(1:structure1%ndim) + 1)

    allocate(gVectors%gVecIndex(-gVectors%maxGLen(1):gVectors%maxGLen(1), &
           & -gVectors%maxGLen(2):gVectors%maxGLen(2), &
           & -gVectors%maxGLen(3):gVectors%maxGLen(3)),&
           & stat=error)
    if(error /= 0) stop 'Error allocating gVectors%gVecIndex in GVectorsInit.'

    gVectors%gVecIndex   = 0

    allocate(gVectors_temp(structure1%ndim, maxNumGVec), &
           & kineticE_temp(maxNumGVec), &
           & stat=error)
    if(error /= 0) stop 'Error allocating temporary space in GVectorInit.'

    ! Loop over all cells in the parallelepiped containing the sphere of 
    ! G vectors with energy less than maxEnergy.
    do i = -gVectors%maxGLen(1), gVectors%maxGLen(1)
      do j = -gVectors%maxGLen(2), gVectors%maxGLen(2)
        do k = -gVectors%maxGLen(3), gVectors%maxGLen(3)

          g_trial = (/ i, j, k /)
          kinetic = dot_product(g_trial(1:structure1%ndim), &
                   & matmul(structure1%lattice%bLatMetricTensor,g_trial(1:structure1%ndim)))/two 

          ! If the energy is less than maxEnergy store the G vector for
          ! later StructureMod.
          if(kinetic < maxEnergy) then

            gVectors%numGVectors = gVectors%numGVectors + 1
            gVectors_temp(:,gVectors%numGVectors) = g_trial(1:structure1%ndim)
            kineticE_temp(gVectors%numGVectors) = kinetic
            gVectors%gVecIndex(i,j,k) = gVectors%numGVectors  ! determine index array  

          end if
        end do
      end do
    end do

    ! allocate final gVectors array and copy over from temporary storage
    allocate(gVectors%gVectors(structure1%ndim, gVectors%numGVectors), &
           & gVectors%kineticE(gVectors%numGVectors), stat=error)
    if(error /= 0) stop 'Error allocating space for G vectors in GVectorsInit.'

    ! transfer vectors to permanent storage
    gVectors%gVectors(:,1:gVectors%numGVectors) = gVectors_temp(:,1:gVectors%numGVectors) 
    gVectors%kineticE(1:gVectors%numGVectors)   = kineticE_temp(1:gVectors%numGVectors)

    deallocate(gVectors_temp,kineticE_temp) ! deallocate associated work space
    ! gVectors, kineticE, and gVecIndex initialized.
    !**********************************************************************

    ! Initialize the structure factors.
    call StructureFactorInit( gVectors, structure1, v_of_q )

    POP_SUB(GVectorsInit)
  end Subroutine GVectorsInit

!*************************************************************************
  Subroutine GVectorsDestroy( gVectors )
!*************************************************************************

! Deallocates the space allocated for the GVectors.

    implicit none

    type( gVectorObjectsT ), pointer :: gVectors
    integer :: error

    PUSH_SUB(GVectorsDestroy)

    if( .not. associated( gVectors )) stop 'Error: gVectors not allocated.'

    deallocate(gVectors%gVectors, gVectors%kineticE, gVectors%structureFactor, stat = error)
    if(error /= 0) stop 'Error gVectors deallocating in GVectorsDestroy.'

    deallocate(gVectors%gVecIndex, stat = error)
    if(error /= 0) stop 'Error gVectors%gVecIndex deallocating in GVectorsDestroy.'

    deallocate(gVectors, stat = error)
    if(error /= 0) stop 'Error gVectors deallocating in GVectorsDestroy.'

    POP_SUB(GVectorsDestroy)
  end Subroutine GVectorsDestroy

!*************************************************************************
  Subroutine StructureFactorInit( gVectors, structure1, v_of_q )
!*************************************************************************

!   calculates the structure factors for each atom type at each G-vector

    use SysParams,        only : one, twopi, double
    use StructureMod, only : StructureT

    implicit none

    type( gVectorObjectsT ), pointer :: gVectors
    type( StructureT ), pointer :: structure1
    logical :: v_of_q

    real(double) :: arg
    integer      :: iGVec, iElem, iAtom, iAtomMin, error, atomI

    PUSH_SUB(StructureFactorInit)

    allocate (gVectors%structureFactor(gVectors%numGVectors, structure1%atomBasis%numSpecies),stat=error)
    if(error /= 0) stop 'structureFactor allocation error'

    atomI = 0
    gVectors%structureFactor = 0

    ! Structure Factor exist for each G vector and each element.  Each
    ! structure factor has a contribution from each atom of that element =
    ! exp(i* 2 * pi (G * R))
    do iElem = 1, structure1%atomBasis%numSpecies

      if (v_of_q) then
        iAtomMin = 1
      else
        iAtomMin = structure1%atomBasis%numAtomsOfSpecies(iElem)
      endif

      do iAtom = iAtomMin, structure1%atomBasis%numAtomsOfSpecies(iElem)

        atomI = atomI + 1

        do iGVec = 1, gVectors%numGVectors                             

          arg = twopi * dot_product(gVectors%gVectors(:,iGVec), structure1%atomBasis%atomLat(:,atomI))
          gVectors%structureFactor(iGVec,iElem) = gVectors%structureFactor(iGVec,iElem) + &
                                   cmplx(cos(arg),-sin(arg),kind=double)

        end do
      end do
    end do

    ! Normalize the structure factors.
    if (v_of_q) then
      gVectors%structureFactor = gVectors%structureFactor / structure1%atomBasis%totalNumAtoms
    endif

    POP_SUB(StructureFactorInit)
  end Subroutine StructureFactorInit

! These legacy routines are not currently used and are wrong for anisotropic systems.
!!*************************************************************************
!  Subroutine GVectorScatter( gVectors, structure1, eigenStates, densityArray, k, n )
!!*************************************************************************
!
!    use StructureMod,    only : StructureT
!    use EigenStatesMod,  only : eigenStatesT
!    use densityArrayMod, only : densityArrayT
!
!    implicit none
!
!    type( gVectorObjectsT ),      pointer :: gVectors
!    type( StructureT ),     pointer :: structure1
!    type( eigenStatesT ),   pointer :: eigenStates
!    type( densityArrayT ), pointer :: densityArray
!
!    integer, intent(in) :: k, n
!    integer             :: i, j, index
!
!    PUSH_SUB(GVectorScatter)
!
!    densityArray%waveFunction = 0
!
!    do i = 1, eigenStates%numBasisVectors(k)  
!
!      if(gVectors%gVectors(structure1%ndim, eigenStates%basisIndex(i,k)) .ge. 0) then
!        index = gVectors%gVectors(structure1%ndim, eigenStates%basisIndex(i,k))
!      else
!        index = gVectors%gVectors(structure1%ndim, eigenStates%basisIndex(i,k)) + &
!              & densityArray%gridDimensions(structure1%ndim)
!      end if
!
!      do j = structure1%ndim-1, 1, -1
!        if(gVectors%gVectors(j, eigenStates%basisIndex(i,k)) .ge. 0) then
!          index = index * densityArray%gridDimensions(j+1) + gVectors%gVectors(j, eigenStates%basisIndex(i,k))
!        else
!          index = index * densityArray%gridDimensions(j+1) + gVectors%gVectors(j, eigenStates%basisIndex(i,k)) &
!                & + densityArray%gridDimensions(j)
!        end if
!      end do  
!      index = index + 1
!      densityArray%waveFunction(index) = eigenStates%eigenVectors(i, n, k)
!
!    end do
!
!    POP_SUB(GVectorScatter)
!  end Subroutine GVectorScatter
!
!!*************************************************************************
!  Subroutine GVectorCollect( gVectors, structure1, densityArray)
!!*************************************************************************
!
!    use StructureMod,    only : StructureT
!    use densityArrayMod, only : densityArrayT
!
!    implicit none
!
!    type( gVectorObjectsT ),      pointer :: gVectors
!    type( StructureT ),     pointer :: structure1
!    type( densityArrayT ), pointer :: densityArray
!
!    integer             :: i, j, index
!
!    PUSH_SUB(GVectorCollect)
!
!    densityArray%rhog = 0
!
!    do i = 1, gVectors%numGVectors
!
!      if(gVectors%gVectors(structure1%ndim, i) .ge. 0) then
!        index = gVectors%gVectors(structure1%ndim, i)
!      else
!        index = gVectors%gVectors(structure1%ndim, i) + &
!              & densityArray%gridDimensions(structure1%ndim)
!      end if
!
!      do j = structure1%ndim-1, 1, -1
!        if(gVectors%gVectors(j, i) .ge. 0) then
!          index = index * densityArray%gridDimensions(j+1) + gVectors%gVectors(j, i)
!        else
!          index = index * densityArray%gridDimensions(j+1) + gVectors%gVectors(j, i) &
!                & + densityArray%gridDimensions(j)
!        end if
!      end do  
!      index = index + 1
!      densityArray%rhog(i) = densityArray%densityArray(index)
!
!    end do
!
!    POP_SUB(GVectorCollect)
!  end Subroutine GVectorCollect

end Module GVectorsMod

