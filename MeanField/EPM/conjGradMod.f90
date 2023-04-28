#include "f_defs.h"

!*************************************************************************
Module ConjGradMod
!*************************************************************************

  use push_pop_m
  implicit none

  contains

!*************************************************************************
  Subroutine CGDiagonalize( hMatrix, eigenStates, kPoint,tol,period,ReportNo)
!*************************************************************************

    use SysParams,         only : double
    use TypeMod,           only : minVecT
    use EigenStatesMod,    only : EigenStatesT
    use HamiltonianMod,    only : hamiltonianArrayT
    use InitializeMod,     only : MinVecAllocate, MinVecInit, &
                                & MinVecDestroy, EigenStatesInit
    use NormalizeMod,      only : EigenStatesOrth
    use tiseOpsMod,        only : CalculateEPsi


    integer, intent( in ) :: kPoint,ReportNo,period
    real(double), intent(in) :: tol

    type( eigenStatesT ), pointer :: eigenStates
    type( hamiltonianArrayT ), pointer :: hMatrix
    type( MinVecT ), pointer :: minVec, oldVec
    complex( double ), dimension(:,:), allocatable :: h

    integer :: iBand, cntr, oldband
    real(double) :: tolerance

    PUSH_SUB(CGDiagonalize)
    
    allocate( h( hMatrix%hSize, hMatrix%hSize ))

    eigenStates%hSize = hMatrix%hSize
    eigenStates%hDim = hMatrix%hDim

    ! The size of HMatrix might be bigger than hMatrix%hSize, with only the 
    ! hMatrix%hSize * hMatrix%hSize part filled
    Call ConvertUTToFullMatrix( hMatrix%h, h, hMatrix%hSize )

    Call MinVecAllocate( minVec, hMatrix%hDim )
    Call MinVecAllocate( oldVec, hMatrix%hDim )

    cntr = 1
    iBand = 1
    tolerance = tol !eigenStates%tolerance

    do while( iBand <= eigenStates%numStates )

       oldband = iBand

       Call MinVecInit( minVec )
       Call MinVecInit( oldVec )

       Call EigenStatesInit( eigenStates, iBand, kPoint )

       Call EigenStatesOrth( eigenStates, iBand, kPoint )

       Call CalculateEPsi( eigenStates, h, iBand, kPoint )

       ! Calculate the Ritz values
       eigenStates%eigenValues( iBand, kPoint ) = &
           & Dot_Product( eigenStates%eigenVectors( :, iBand, kPoint ), &
           & eigenStates%ePsi( :, kPoint ))



       Call Diagonalize( hMatrix%hSize, h, eigenStates,& 
       & iBand, kPoint, tolerance, minVec, oldVec ,ReportNo)

       iBand = iBand + 1

       !   If the calculation for the band has converged then iBand has
       ! increased by one; if so then reset counter to 1 and tolerance to
       ! tol; else increase counter, if the counter crosses 10 (or some 
       ! user-defined number), then reset counter to 1 and decrease 
       ! tolerance. 
       if( oldband /= iBand) then
          cntr = 1
          tolerance = tol !eigenStates%tolerance
       else
          cntr = cntr + 1
          if( cntr > period )then
             cntr = 1
             tolerance = 100 * tolerance
             write( ReportNo, * ) 'Failed to converge in', period, &
                  & ' steps: Decreasing tolerance to ', tolerance
          end if
       end if

    end do

    Call ConvertFullToUTMatrix( hMatrix%h, h, hMatrix%hSize)

    deallocate( h )
    Call MinVecDestroy( minVec )
    Call MinVecDestroy( oldVec )

    POP_SUB(CGDiagonalize)

  end Subroutine CGDiagonalize

!*************************************************************************
  Subroutine ConvertUTToFullMatrix( hMatrix, h, hSize )
!*************************************************************************

    use SysParams, only : double

    integer, intent(in) :: hSize
    complex(double), dimension(:,:) :: h
    complex(double), pointer :: hMatrix(:)

    integer :: i, j

    PUSH_SUB(ConvertUTToFullMatrix)

    do i = 1, hSize
       h( i, i ) = hMatrix( i + i * ( i - 1 ) / 2 )
       do j = i + 1, hSize
          h(i,j) = hMatrix( i + j * ( j - 1 ) / 2 )
          h(j,i) = conjg(h(i,j))
       end do
    end do

    POP_SUB(ConvertUTToFullMatrix)
  end Subroutine ConvertUTToFullMatrix

!*************************************************************************
  Subroutine ConvertFullToUTMatrix( hMatrix, h, hSize )
!*************************************************************************

    use SysParams, only : double

    integer, intent(in) :: hSize
    complex(double), dimension(:,:) :: h
    complex(double), pointer :: hMatrix(:)

    integer :: i, j

    PUSH_SUB(ConvertFullToUTMatrix)

    do i = 1, hSize
       hMatrix( i + i * ( i - 1 ) / 2 ) = h( i, i )
       do j = i + 1, hSize
          hMatrix( i + j * ( j - 1 ) / 2 ) = h(i,j) 
       end do
    end do

    POP_SUB(ConvertFullToUTMatrix)

  end Subroutine ConvertFullToUTMatrix

!*************************************************************************
  Subroutine Diagonalize( hSize, hMatrix, eigenStates, & 
                        & iBand, kPoint, tol, minVec, oldVec, ReportNo)
!*************************************************************************

    use SysParams,         only : double
    use EigenStatesMod,    only : eigenStatesT
    use TypeMod,           only : minVecT
    use MinimizeMod,       only : StDescentVec, ConjDirection
    use NormalizeMod,      only : LowerBandOrthogonalize, BandOrthogonalize
    use tiseOpsMod,        only : PreCondition, ExactEnergyMinimize

    ! Define the types of the subroutine arguments
    integer, intent( in ) :: hSize, kPoint, ReportNo
    integer, intent( inout ) :: iBand
    real( double ), intent( in ) :: tol
    complex( double ), dimension(:,:) :: hMatrix

    type(eigenStatesT), pointer :: eigenStates
    type(MinVecT), pointer :: MinVec,OldVec

    integer :: iter
    real(double) :: residue, tolsq

    PUSH_SUB(Diagonalize)

    tolsq = 10.d0*tol**2

    !   Iterate the Eigen Vector until convergence before moving over to the
    ! next higher band.
    EigenVectorLoop: do iter = 1, hSize

       !   Evaluate the Steepest Descent diection as the negative gradient
       ! In this case it is StDesc = -(H-lambda)|Psi>        
       Call StDescentVec( eigenStates, MinVec, iBand, kPoint ) ! In "Minimize.f90"


       !   The residue is the Steepst Descent Vector, so its square is 
       ! calculated, and the iterations are carried out until userdefined
       ! convergence is reached, OR 'HamSize' iterations have been carried out.
       residue = Dot_Product( minVec%StDesc, minVec%StDesc )

       if( residue > tolsq ) then
          
          !   Orhtogonalize the steepst descent vector against all the LOWER
          ! Eigen Vectors. (In "NormalizeMod.f90")
          Call LowerBandOrthogonalize( eigenStates, minVec%StDesc, &
                                     & iBand, kPoint )

          !   PreCondition the Steepest Descent Vector In "HamiltonMod.f90"
          Call PreCondition( eigenStates, minVec, hMatrix, iBand, kPoint )
          
          !   Orthogonalize the PreConditioned Vector against all bands, lower 
          ! than and including the current band.  In "NormalizeMod.f90"
          Call BandOrthogonalize( eigenStates, minVec%StDesc, iBand, kPoint )
          
          !   Evaluate the conjugate direction orthogonal to the current band 
          ! and normalized. (In "MinimizeMod.f90")
          Call ConjDirection( eigenStates, minVec, oldVec, iBand, kPoint )

          !   Do the line minimization by searching for the energy minimum, and
          ! update the band.
          Call ExactEnergyMinimize( eigenStates, minVec, &
                                  & hMatrix, iBand, kPoint )
       else
          
          exit EigenVectorLoop

       end if

    end do EigenVectorLoop
    
    ! Print the eigenvalue and residue
    write( ReportNo, 11 ) iBand, eigenStates%eigenValues( iBand, kPoint ),&
         & sqrt( residue ), iter

11     format(1X,'Band:',I3,1X,'EigenValues:',E15.5,1X,'Norm Resed:',E15.5, &
            & 1X,'Iterations:',I3)

    if( residue > tolsq ) then
       iBand = iBand - 1
    end if

    POP_SUB(Diagonalize)
  end Subroutine Diagonalize

!*************************************************************************
  Subroutine cgEigenSystem( hamiltonian,eigenStates,k,tol,period,ReportNo )
!*************************************************************************

! Hamiltonian matrix is diagonalized by a call to conjugate gradient

    use SysParams,      only : double
    use EigenStatesMod, only : eigenStatesT
    use HamiltonianMod, only : hamiltonianArrayT

    implicit none
    
    integer, intent(in)                :: period,ReportNo,k
    real(double), intent(in)           :: tol
    type( hamiltonianArrayT ), pointer :: hamiltonian
    type( eigenStatesT ),      pointer :: eigenStates

    PUSH_SUB(cgEigenSystem)

    Call CGDiagonalize(hamiltonian,eigenStates,k,tol,period,ReportNo)

    POP_SUB(cgEigenSystem)
  end Subroutine cgEigenSystem

!*************************************************************************
  Subroutine cgFFTDiagonalise(hamArray,eigenStates,Grid,kPoint,tol,period, &
       & ReportNo)
!*************************************************************************

    use SysParams,           only : double
    use TypeMod,             only : minVecT,GridT
    use EigenStatesMod,      only : EigenStatesT
    use HamiltonianMod,      only : hamiltonianArrayT
    use InitializeMod,       only : MinVecAllocate, MinVecInit, &
                                  & MinVecDestroy, EigenStatesInit
    use NormalizeMod,        only : EigenStatesOrth
    use tiseOpsMod,          only : CalculateEPsi,CalcHPsi


    integer, intent( in )              :: kPoint,ReportNo,period
    real(double), intent(in)           :: tol

    type(GridT),               pointer :: Grid 
    type(eigenStatesT),        pointer :: eigenStates
    type(hamiltonianArrayT),   pointer :: hamArray
    type(MinVecT),             pointer :: minVec, oldVec
    
    integer :: iBand, cntr, oldband
    real(double) :: tolerance

    PUSH_SUB(cgFFTDiagonalise)

    eigenStates%hSize = hamArray%hSize
    eigenStates%hDim = hamArray%hDim

    Call MinVecAllocate( minVec, hamArray%hDim )
    Call MinVecAllocate( oldVec, hamArray%hDim )

    cntr = 1
    iBand = 1
    tolerance = tol 

    do while( iBand <= eigenStates%numStates )

       oldband = iBand

       Call MinVecInit( minVec )
       Call MinVecInit( oldVec )

       Call EigenStatesInit(eigenStates,iBand,kPoint)

       Call EigenStatesOrth(eigenStates,iBand,kPoint)

       Call CalcHPsi(eigenStates%hSize,hamArray%KineticG, &
            & eigenStates%eigenVectors(:,iBand,kPoint),eigenStates%epsi(:,kPoint), &
            & Grid)

       ! Calculate the Ritz values
       eigenStates%eigenValues( iBand, kPoint ) = &
           & Dot_Product( eigenStates%eigenVectors( :, iBand, kPoint ), &
           & eigenStates%ePsi( :, kPoint ))

       Call FFTDiagonalise(hamArray%KineticG,eigenStates,Grid,iBand,kPoint, &
       & tolerance,minVec,oldVec,ReportNo)

       iBand = iBand + 1

       !   If the calculation for the band has converged then iBand has
       ! increased by one; if so then reset counter to 1 and tolerance to
       ! tol; else increase counter, if the counter crosses 10(or some 
       ! user defined number), then reset counter to 1 and decrease 
       ! tolerance. 
       if( oldband /= iBand) then
          cntr = 1
          tolerance = tol !eigenStates%tolerance
       else
          cntr = cntr + 1
          if( cntr > period )then
             cntr = 1
             tolerance = 100 * tolerance
             write( ReportNo, * ) 'Failed to converge in', period, &
                  & ' steps: Decreasing tolerance to ', tolerance
          end if
       end if

    end do

    Call MinVecDestroy( minVec )
    Call MinVecDestroy( oldVec )

    POP_SUB(cgFFTDiagonalise)
  end Subroutine cgFFTDiagonalise

!*************************************************************************
  Subroutine FFTDiagonalise(T_G,eigenStates,Grid,iBand,kPoint, &
       & tol,minVec,oldVec,ReportNo)
!*************************************************************************

    use SysParams,           only : double
    use TypeMod,             only : minVecT,GridT
    use EigenStatesMod,      only : EigenStatesT
    use MinimizeMod,         only : StDescentVec, ConjDirection
    use NormalizeMod,        only : LowerBandOrthogonalize, BandOrthogonalize
    use tiseOpsMod,          only : FFTPreCondition, FFTEnergyMinimize

    type(GridT),           pointer :: Grid
    type(eigenStatesT),    pointer :: eigenStates
    type(minVecT),         pointer :: minVec,oldVec
    integer, intent(in)            :: kPoint,ReportNo
    real(double)                   :: tol
    real(double), dimension(:)     :: T_G

    integer :: iter,iBand
    real(double) :: residue, tolsq

    PUSH_SUB(FFTDiagonalise)

    tolsq = 10.d0*tol**2

    !   Iterate the Eigen Vector until convergence before moving over to the
    ! next higher band.
    EigenVectorLoop: do iter = 1, eigenStates%hSize

       !   Evaluate the Steepest Descent diection as the negative gradient
       ! In this case it is StDesc = -(H-lambda)|Psi>        
       Call StDescentVec( eigenStates, MinVec, iBand, kPoint )

       !   The residue is the Steepst Descent Vector, so its square is 
       ! calculated, and the iterations are carried out until userdefined
       ! convergence is reached, OR 'HamSize' iterations have been carried out.
       residue = Dot_Product( minVec%StDesc, minVec%StDesc )

       if( residue > tolsq ) then
          
          !   Orhtogonalize the steepst descent vector against all the LOWER
          ! Eigen Vectors. (In "NormalizeMod.f90")
          Call LowerBandOrthogonalize( eigenStates, minVec%StDesc, &
                                     & iBand, kPoint )

          !   PreCondition the Steepest Descent Vector In "HamiltonMod.f90"
          Call FFTPreCondition(eigenStates,minVec,T_g,iBand,kPoint)

          !   Orthogonalize the PreConditioned Vector against all bands, lower 
          ! than and including the current band.  In "NormalizeMod.f90"
          Call BandOrthogonalize( eigenStates, minVec%StDesc, iBand, kPoint )
          
          !   Evaluate the conjugate direction orthogonal to the current band 
          ! and normalized. (In "MinimizeMod.f90")
          Call ConjDirection( eigenStates, minVec, oldVec, iBand, kPoint )

          !   Do the line minimization by searching for the energy minimum, and
          ! update the band.
          Call FFTEnergyMinimize(eigenStates,MinVec,T_g,Grid,iBand,kPoint)

       else
          
          exit EigenVectorLoop

       end if

    end do EigenVectorLoop
    
    ! Print the eigenvalue and residue
    write( ReportNo, 11 ) iBand, eigenStates%eigenValues( iBand, kPoint ),&
         & sqrt( residue ), iter

11     format(1X,'Band:',I3,1X,'EigenValues:',E15.5,1X,'Norm Resed:',E15.5, &
            & 1X,'Iterations:',I3)


    if( residue > tolsq ) then
       iBand = iBand - 1
    end if

    POP_SUB(FFTDiagonalise)
  end Subroutine FFTDiagonalise

end Module ConjGradMod

