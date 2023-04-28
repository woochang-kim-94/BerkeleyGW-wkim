#include "f_defs.h"

!*************************************************************************
Module HamiltonianMod
!*************************************************************************
! gsm [2010-10-31] added zheevx and workspace query for ACML,
! kept zhpevx (implemented in the original EPM code) for MKL,
! replaced abstol = -1.0 with AbsoluteTolerance (default = -1.0)
!-------------------------------------------------------------------------
! on hbar, eigenvectors are reproduced in the following configurations:
! ACML/zheevx with abstol = -1.0, MKL/zhpevx with abstol = -1.0,
! ACML/zhpevx with abstol = 1.0d-12, MKL/zheevx with abstol = 1.0d-14
!*************************************************************************
! HDiagonalize calls LAPACK routines
!*************************************************************************
! The module contains the public functions:
!    HInit
!    HDestroy
!    HDiagonalize (calls LAPACK routines)
!    HPrint
!    HCheck
!*************************************************************************
! The module contains the public variables:
!*************************************************************************
! The module contains the private functions:
!*************************************************************************
! The module contains the private variables:
!*************************************************************************
! The module contains the type:
! hamiltonianArrayT which contains the variables:
!    hDim        Maximum dimension of the Hamiltonian matrix that gets used.
!    hSize       Actual dimension of the Hamiltonian matrix that gets used.
!                The variable name is misleading, but we`re stuck with it for now.
!    h(i)        I-th element of the Hamiltonian matrix in packed storage
!    SqH(i,j)    The complete Hamiltonian matrix 
!    work        work arrays for LAPACK/Conjugate gradient
!*************************************************************************

  use SysParams, only : double
  use push_pop_m
  implicit none

  save

  type hamiltonianArrayT
    integer                   :: hSize, hDim, lwork
    real(double)              :: abstol
    complex(double), pointer  :: h(:)
    complex(double), pointer  :: SqH(:,:)
    complex(double), pointer  :: PotentialR(:)
    real(double),    pointer  :: KineticG(:)
    complex(double), pointer  :: work(:)
    real(double),    pointer  :: rwork(:)
    integer,         pointer  :: iwork(:)
    integer,         pointer  :: fail(:)
  end type hamiltonianArrayT

  contains

!*************************************************************************
  Subroutine HInit(hamiltonian, maxHDim)
!*************************************************************************

!   Allocate space for the Hamiltonian and the work matrices.

    use SysParams,      only : one, zero
    use lapack_m
    implicit none

    type( hamiltonianArrayT ), pointer :: hamiltonian
    integer, intent(in) :: maxHDim

    integer :: error
    ! Workspace variables
#ifdef UNPACKED
    integer             :: dummym, dummyi(1), dummyf(1)
    real(double)        :: dummyw(1), dummyr(1)
    complex(double)     :: dummyh(1,1), dummyz(1,1), dummywork(1)
#endif

    PUSH_SUB(HInit)
    
    if( associated( hamiltonian )) stop 'Error, hamiltonian already allocated.'
    allocate( hamiltonian, stat = error )
    if( error /= 0 ) stop 'Error, hamiltonian allocation failed.'
    nullify( hamiltonian%h )
    nullify( hamiltonian%SqH )
    nullify( hamiltonian%PotentialR )
    nullify( hamiltonian%KineticG )
    nullify( hamiltonian%work )
    nullify( hamiltonian%rwork )
    nullify( hamiltonian%iwork )
    nullify( hamiltonian%fail )

    hamiltonian%hDim  = maxHDim
    
    ! Assume we used the entire Hamiltonian matrix.
    hamiltonian%hSize = maxHDim

    ! workspace query
#ifdef UNPACKED
    ! for ACML
    hamiltonian%lwork = -1
    call zheevx( 'V','I','U', hamiltonian%hDim, dummyh, &
               & hamiltonian%hDim, zero, zero, 1, hamiltonian%hDim, &
               & hamiltonian%abstol, dummym, dummyw, dummyz, hamiltonian%hDim, &
               & dummywork, hamiltonian%lwork, dummyr, dummyi, dummyf, error )
    hamiltonian%lwork = nint(dble(dummywork(1)))
#else
    ! for MKL and Netlib
    hamiltonian%lwork = 2 * hamiltonian%hDim
#endif

    ! The Hamiltonian matrix is stored in packed form. Only the upper triangle
    ! plus half the main diagonal is stored. This is why the allocate statement
    ! looks funny.
    allocate( hamiltonian%h( hamiltonian%hDim * ( hamiltonian%hDim + 1 ) / 2 ), &
            & hamiltonian%sqH( hamiltonian%hDim, hamiltonian%hDim ), &
            & hamiltonian%KineticG(hamiltonian%hDim), &
            & hamiltonian%PotentialR(hamiltonian%hDim), &
            & hamiltonian%fail( hamiltonian%hDim ), stat=error )
    if(error /= 0) stop 'Error in allocating Hamiltonian space, part 1.'
    allocate( hamiltonian%work( hamiltonian%lwork ), &
            & hamiltonian%rwork( 7 * hamiltonian%hDim ), &
            & hamiltonian%iwork( 5 * hamiltonian%hDim ), stat=error )
    if(error /= 0) stop 'Error in allocating Hamiltonian space, part 2.'

    POP_SUB(HInit)
  end Subroutine HInit

!*************************************************************************
  Subroutine HDestroy( hamiltonian )
!*************************************************************************

    implicit none

    type( hamiltonianArrayT ), pointer :: hamiltonian
    integer :: error

    PUSH_SUB(HDestroy)

    if( .not. associated( hamiltonian )) stop 'Error, hamiltonian not allocated.'

    deallocate( hamiltonian%h, hamiltonian%sqH, &
              & hamiltonian%work, hamiltonian%rwork, &
              & hamiltonian%iwork, hamiltonian%fail, &
              & hamiltonian%KineticG, hamiltonian%PotentialR, stat=error )
    if(error /= 0) stop 'Error in deallocating Hamiltonian space.'

    deallocate( hamiltonian, stat = error )
    if(error /= 0) stop 'Error, hamiltonian deallocation failed.'

    nullify( hamiltonian )

    POP_SUB(HDestroy)
  end Subroutine Hdestroy

!*************************************************************************
  Subroutine HDiagonalize( hamiltonian, eigenStates, k )
!*************************************************************************

! Hamiltonian matrix is diagonalized by a call to LAPACK routines.

    use SysParams,      only : one, zero
    use EigenStatesMod, only : eigenStatesT
    use lapack_m

    implicit none

    type( hamiltonianArrayT ), pointer :: hamiltonian
    type( eigenStatesT ),      pointer :: eigenStates
    integer, intent( in )         :: k

    integer :: neig, error, eigenvecdim
    real(double), allocatable :: eigenvalues(:)
    complex(double), allocatable :: eigenvectors(:,:)

#ifdef DEBUG
!#ifdef UNPACKED
    complex(double), allocatable :: sqH(:,:)
    complex(double), allocatable :: packed(:)
    integer :: packed_size
    PUSH_SUB(HDiagonalize)
    allocate(sqH(hamiltonian%hDim, hamiltonian%hDim))
    sqH(1:hamiltonian%hDim, 1:hamiltonian%hDim) = hamiltonian%sqH(1:hamiltonian%hDim, 1:hamiltonian%hDim)
#ifndef UNPACKED
    packed_size = hamiltonian%hSize * (hamiltonian%hSize + 1)/2
    allocate(packed(1:packed_size))
    packed(1:packed_size) = hamiltonian%h(1:packed_size)
#endif
#else
    PUSH_SUB(HDiagonalize)
#endif

    allocate(eigenvalues(1:hamiltonian%hSize))
    eigenvecdim = max(hamiltonian%hSize, eigenStates%maxBasisVectors)
    allocate(eigenvectors(1:eigenvecdim, eigenStates%numStates))

#ifdef UNPACKED
    ! for ACML
    call zheevx( 'V','I','U', hamiltonian%hSize, hamiltonian%SqH, hamiltonian%hDim, &
#else
    ! for MKL and Netlib
    call zhpevx( 'V','I','U', hamiltonian%hSize, hamiltonian%h, &
#endif
               & zero, zero, 1, eigenStates%numStates, &
               & hamiltonian%abstol, neig, eigenvalues(:), &
               & eigenvectors(:,:), eigenvecdim, hamiltonian%work, &
#ifdef UNPACKED
               & hamiltonian%lwork, &
#endif
               & hamiltonian%rwork, hamiltonian%iwork, hamiltonian%fail, error )

    if(neig .lt. eigenStates%numStates) then
       write(0,'(a,i6,a,i6)') 'Error: failed in zheevx: neig = ', neig, ' < numStates = ', eigenStates%numStates
       stop
    endif
    if(error /= 0) then
       write(0,'(a,i6)') 'Error: failed in zheevx: code = ', error
       stop
    endif

    eigenStates%eigenValues(1:eigenStates%numStates, k) = eigenvalues(1:eigenStates%numStates)
    eigenStates%eigenVectors(1:eigenStates%maxBasisVectors, 1:eigenStates%numStates, k) = &
      & eigenvectors(1:eigenStates%maxBasisVectors, 1:eigenStates%numStates)
    deallocate(eigenvalues, eigenvectors)

#ifdef DEBUG
!#ifdef UNPACKED
    hamiltonian%sqH(1:hamiltonian%hDim, 1:hamiltonian%hDim) = sqH(1:hamiltonian%hDim, 1:hamiltonian%hDim)
    deallocate(sqH)
#ifndef UNPACKED
    hamiltonian%h(1:packed_size) = packed(1:packed_size)
    deallocate(packed)
#endif
#endif

    POP_SUB(HDiagonalize)
  end Subroutine HDiagonalize

!*************************************************************************
  Subroutine HPrint( hamiltonian )
!*************************************************************************

    use message_m
    implicit none

    type( hamiltonianArrayT ), pointer :: hamiltonian
    integer :: i, j

    PUSH_SUB(HPrint)

    call open_file(7, file = 'h.dat', status='replace')

    do i = 1, hamiltonian%hSize
      do j = i, hamiltonian%hSize
        write( unit = 7, fmt = '(F8.4)', advance = 'NO') &
             & dble( hamiltonian%h( i + j * ( j - 1 ) / 2 ))
      end do
      write(unit = 7, fmt = '()')
    end do

    POP_SUB(HPrint)
  end Subroutine HPrint

!*************************************************************************
  Subroutine HCheck( hamiltonian, eigenStates, ik )
! Added by DAS, 9 Nov 2010
!*************************************************************************

    use SysParams,      only : one, zero
    use EigenStatesMod, only : eigenStatesT
    use blas_m
    implicit none

    type( hamiltonianArrayT ), pointer :: hamiltonian
    type( eigenStatesT ),      pointer :: eigenStates
    integer, intent(in) :: ik

    integer :: ist, ist2, ng, ig, ig2
    real(double) :: diff, norm
    complex(double) :: energy, overlap
    complex(double), allocatable :: Hpsi(:), Hpsimat(:,:), Hdiag(:,:)
    logical, parameter :: noblas = .true. ! enable naive ways to do it without BLAS

    PUSH_SUB(HCheck)

    ng = hamiltonian%hSize
    allocate(Hpsi(ng))

    do ist = 1, eigenStates%numStates
      if(noblas) then

        energy = 0d0
        do ig = 1, ng
          do ig2 = 1, ng
            energy = energy + conjg(eigenStates%eigenVectors(ig2, ist, ik)) * &
              hamiltonian%SqH(ig2, ig) * eigenStates%eigenVectors(ig, ist, ik)
          enddo
        enddo

      else

#ifdef UNPACKED
        call zhemv('U', ng, (one, zero), hamiltonian%SqH(:,:), hamiltonian%hDim, &
          eigenStates%eigenVectors(:, ist, ik), 1, (zero, zero), Hpsi(:), 1)
#else
        call zhpmv('U', ng, (one, zero), hamiltonian%h(:), &
          eigenStates%eigenVectors(:, ist, ik), 1, (zero, zero), Hpsi(:), 1)
#endif
        energy = zdotc(ng, eigenStates%eigenVectors(:, ist, ik), 1, Hpsi(:), 1)
      endif

! can`t multiply these ways since the leading dimension of SqH is not ng
!      energy = dot_product(conjg(eigenStates%eigenVectors(:, ist, ik)), &
!        matmul(hamiltonian%SqH(:,:), eigenStates%eigenVectors(:, ist, ik)))
!      energy = dot_product(conjg(eigenStates%eigenVectors(1:ng, ist, ik)), &
!        matmul(hamiltonian%SqH(1:ng, 1:ng), eigenStates%eigenVectors(1:ng, ist, ik)))


      if(abs(aimag(energy)) > 1d-10) then
        write(0,'(a)') 'ERROR: EPM diagonalization failure -- expectation value is complex.'
        write(0,'(a,i4,a,i4,a,e13.6)') 'ik = ', ik, '   ist = ', ist, '   aimag(energy) = ', aimag(energy)
        write(0,*) 'energy = ', energy
        write(0,*) 'abs(energy) = ', abs(energy)
        write(0,*) 'eigenvalue = ', eigenStates%eigenValues(ist, ik)
        stop
      endif

      diff = dble(energy) - eigenStates%eigenValues(ist, ik)
      if(abs(diff) > 1d-12) then
        write(0,'(a)') 'ERROR: EPM diagonalization failure -- expectation value differs from eigenvalue.'
        write(0,'(a,i4,a,i4,a,e13.6)') 'ik = ', ik, '   ist = ', ist, '   diff = ', diff
        stop
      endif

      do ist2 = 1, eigenStates%numStates
        if(noblas) then
          overlap = 0d0
          do ig = 1, ng
            overlap = overlap + conjg(eigenStates%eigenVectors(ig, ist, ik)) * eigenStates%eigenVectors(ig, ist2, ik)
          enddo
        else
          overlap = zdotc(ng, eigenStates%eigenVectors(:, ist, ik), 1, eigenStates%eigenVectors(:, ist2, ik), 1)
        endif

        if(ist == ist2) then
! This should be more efficient, but is seg-faulting for some reason
!            norm = dznrm2(ng, eigenStates%eigenVectors(:, ist, ik), 1)
!          norm = sqrt(zdotc(ng, eigenStates%eigenVectors(:, ist, ik), 1, eigenStates%eigenVectors(:, ist2, ik), 1))
          norm = overlap
          if(abs(norm - 1d0) > 4d-14) then
            write(0,'(a)') 'ERROR: EPM diagonalization failure -- eigenvector is not normalized.'
            write(0,'(a,i4,a,i4,a,e13.6,3x,e13.6)') 'ik = ', ik, '   ist = ', ist, '   norm - 1 = ', norm - 1d0
            stop
          endif
        else
          if(abs(overlap) > 4d-14) then
            write(0,'(a)') 'ERROR: EPM diagonalization failure -- eigenvectors are not orthogonal.'
            write(0,'(a,i4,a,i4,a,i4,a,e13.6,3x,e13.6)') 'ik = ', ik, '  ist1 = ', ist, '  ist2 = ', ist2,'  overlap = ', overlap
            stop
          endif
        endif
      enddo
    enddo

    deallocate(Hpsi)

    allocate(Hdiag(eigenStates%numStates, eigenStates%numStates))
    if(.not. noblas) then
      allocate(Hpsimat(ng, eigenStates%numStates))
! DAS -- not clear if there is a way to do this with the packed representation
      call zhemm('L', 'U', ng, eigenStates%numStates, (one, zero), hamiltonian%SqH(:,:), hamiltonian%hDim, &
        eigenStates%eigenVectors(:, :, ik), eigenStates%maxBasisVectors, (zero, zero), Hpsimat(:, :), ng)
    endif
    
    do ist = 1, eigenStates%numStates
      do ist2 = 1, eigenStates%numStates
        if(noblas) then
          Hdiag(ist2, ist) = 0d0
          do ig = 1, ng
            do ig2 = 1, ng
              Hdiag(ist2, ist) = Hdiag(ist2, ist) + conjg(eigenStates%eigenVectors(ig2, ist2, ik)) * &
                hamiltonian%SqH(ig2, ig) * eigenStates%eigenVectors(ig, ist, ik)
            enddo
! DAS -- way to do it with ZHEMM but not ZDOTC
!          Hdiag(ist2, ist) = Hdiag(ist2, ist) + conjg(eigenStates%eigenVectors(ig, ist2, ik)) * Hpsimat(ig, ist)
          enddo
        else
          Hdiag(ist2, ist) = zdotc(ng, eigenStates%eigenVectors(:, ist2, ik), 1, Hpsimat(:, ist), 1)
        endif
      enddo
    enddo

    if(.not. noblas) deallocate(Hpsimat)

    do ist = 1, eigenStates%numStates
      do ist2 = 1, eigenStates%numStates
        if(ist == ist2) then
          diff = Hdiag(ist, ist) - eigenStates%eigenValues(ist, ik)
          if(abs(diff) > 2d-12) then
            write(0,'(a)') 'ERROR: EPM diagonalization failure -- Hdiag(ist, ist) != eigenvalue(ist).'
            write(0,'(a,i4,a,i4,a,e13.6,3x,e13.6)') 'ik = ', ik, '  ist1 = ', ist, ' diff = ', diff
            stop
          endif
        else
          if(abs(Hdiag(ist2, ist)) > 2d-10) then
            write(0,'(a)') 'ERROR: EPM diagonalization failure -- Hamiltonian not diagonalized.'
            write(0,'(a,i4,a,i4,a,i4,a,e13.6,3x,e13.6)') 'ik = ', ik, '  ist1 = ', ist, &
              '  ist2 = ', ist2,'  element = ', Hdiag(ist2, ist)
            stop
          endif
        endif
      enddo
    enddo
    deallocate(Hdiag)

    POP_SUB(HCheck)
  end Subroutine HCheck

end Module HamiltonianMod

