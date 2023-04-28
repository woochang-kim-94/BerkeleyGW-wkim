!==============================================================================
!
! Modules:
!
! chi_convergence_m   (Originally by DVF)            Last Modified: 05/08/2015 (DVF)
!
!   This module has functions to create the structures for the convergence tests 
!   with respect to to bands in epsilon, to do the convergence tests, and to 
!   print the results.
!
!==============================================================================

#include "f_defs.h"

module chi_convergence_m

  use global_m

  implicit none

  private

  public :: create_chi_converger, free_chi_converger,&
    chi_convergence_test,chi_convergence_print, &
    chi_converger_t

  type chi_converger_t
    SCALAR, allocatable :: head(:)
    SCALAR, allocatable :: maxG(:)
    SCALAR, allocatable :: head2(:)
    SCALAR, allocatable :: maxG2(:)
    SCALAR, allocatable :: headtotal(:)
    SCALAR, allocatable :: maxGtotal(:)
  end type chi_converger_t

contains

  subroutine create_chi_converger(conv,nvalbands,ntotbands)
    type(chi_converger_t), intent(inout) :: conv !< the chi_converger_t object
    integer, intent(in) :: nvalbands  !< number of valence bands
    integer, intent(in) :: ntotbands  !< total number of bands

    PUSH_SUB(create_chi_converger)

    SAFE_ALLOCATE(conv%head,(ntotbands-nvalbands))
    SAFE_ALLOCATE(conv%maxG,(ntotbands-nvalbands))
    SAFE_ALLOCATE(conv%head2,(ntotbands-nvalbands))
    SAFE_ALLOCATE(conv%maxG2,(ntotbands-nvalbands))
    SAFE_ALLOCATE(conv%headtotal,(ntotbands-nvalbands))
    SAFE_ALLOCATE(conv%maxGtotal,(ntotbands-nvalbands))      

    conv%head = ZERO
    conv%head2 = ZERO
    conv%headtotal = ZERO
    conv%maxG = ZERO
    conv%maxG2 = ZERO
    conv%maxGtotal = ZERO

    POP_SUB(create_chi_converger)
    return

  end subroutine create_chi_converger

  subroutine free_chi_converger(conv)
    type(chi_converger_t), intent(inout) :: conv !<the chi_converger_t object

    PUSH_SUB(free_chi_converger)

    SAFE_DEALLOCATE(conv%head)
    SAFE_DEALLOCATE(conv%maxG)
    SAFE_DEALLOCATE(conv%head2)
    SAFE_DEALLOCATE(conv%maxG2)
    SAFE_DEALLOCATE(conv%headtotal)
    SAFE_DEALLOCATE(conv%maxGtotal)      

    POP_SUB(free_chi_converger)
    return

  end subroutine free_chi_converger


  subroutine chi_convergence_test(pol,pht,indt,kp,nrk,nst,nvalbands,ntotbands,fact,conv)
    type(polarizability), intent(in) :: pol !< polarizability
    SCALAR, intent(in) :: pht(:, :, :) !< (pol%nmtx,neqmax,nrk) phases
      !! neqmax = max number of kpoints equivalent by symmetry (across all irr. k-points)
    integer, intent(in) :: indt(:, :, :) !< (pol%nmtx,neqmax,nrk) indexing array 
    type(kpoints),intent(in) :: kp !<kpoint array
    integer, intent(in) :: nrk !<number of kpoints
    integer, intent(in) :: nst(:) !< (nrk) !<number of elements in the star
    integer, intent(in) :: nvalbands !< number of valence bands
    integer, intent(in) :: ntotbands !< number of total bands
    real(DP), intent(in) :: fact !< volume factor used throughout the code
    type(chi_converger_t), intent(inout) :: conv

    SCALAR :: gmehead,gmemaxG  !<matrix elements for head and tail (maxG)
    SCALAR :: mod_square_phthead,mod_square_phtmaxG !<mod squared of phase
    integer :: ispin,irk,iv,ic,istar,isend,iown

    ! JRD/DVF: Test convergence with bands for G=G`=0 and G=G`=pol%nmtx
    ! The result is chi(0,0), chi(Gmax,Gmax) as a function of conduction
    ! bands, i.e. for these gvectors the chi summation is done for all
    ! valence bands and k-points for certain numbers of conduction bands
    ! and the convergence is tracked as the number of conduction bands
    ! is increased. These will be printed in chi_convergence_print.

    PUSH_SUB(chi_convergence_test)

    do ispin =1, kp%nspin
      do irk = 1, nrk
        do istar = 1, nst(irk)

          mod_square_phthead=pht(1,istar,irk)*MYCONJG(pht(1,istar,irk))
          mod_square_phtmaxG=pht(pol%nmtx,istar,irk)*MYCONJG(pht(pol%nmtx,istar,irk))

          do iv=1,(nvalbands+pol%ncrit)
            iown =1 
            do ic=1,ntotbands-nvalbands
              isend=peinf%global_pairowner(iv,ic)-1
              if (isend .eq. peinf%inode) then
                if (iown .gt. peinf%ncownactual) write(6,*) 'iown bigger than ncownactual'
                gmehead = pol%gme(indt(1,istar,irk),iown,peinf%indexv(iv), &
                  ispin,irk,pol%nfreq_group)
                gmemaxG = pol%gme(indt(pol%nmtx,istar,irk),iown,peinf%indexv(iv), &
                  ispin,irk,pol%nfreq_group)
                conv%head2(ic) = conv%head2(ic) + gmehead*MYCONJG(gmehead)*mod_square_phthead*fact*(-1D0)
                conv%maxG2(ic) = conv%maxG2(ic) + gmemaxG*MYCONJG(gmemaxG)*mod_square_phtmaxG*fact*(-1D0)
                iown=iown+1
              endif
            enddo ! ic
          enddo ! iv

        enddo ! istar
      enddo ! irk
    enddo ! ispin

#ifdef MPI
    call MPI_reduce(conv%head2,conv%head,ntotbands-nvalbands,MPI_SCALAR,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
    call MPI_reduce(conv%maxG2,conv%maxG,ntotbands-nvalbands,MPI_SCALAR,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
#else
    conv%head = conv%head2
    conv%maxG = conv%maxG2
#endif

    POP_SUB(chi_convergence_test)

  end subroutine chi_convergence_test

  subroutine chi_convergence_print(pol,iq,nvalbands,ntotbands,conv)

    type(polarizability), intent(in) :: pol !<polarizability
    integer, intent(in) :: iq !<current qpoint
    integer, intent(in) :: nvalbands !<number of valence bands
    integer, intent(in) :: ntotbands !<total number of bands
    type(chi_converger_t), intent(inout) :: conv

    SCALAR :: extrap_head,extrap_maxG
    integer :: ic,idis,ncondbands
    character*24 :: strhead,strtail

    ! JRD: Print out convergence

    PUSH_SUB(chi_convergence_print)

    conv%headtotal(1)=conv%head(1)
    conv%maxGtotal(1)=conv%maxG(1)
    ncondbands = ntotbands-nvalbands
    do ic = 2, ncondbands
      conv%headtotal(ic)=conv%headtotal(ic-1)+conv%head(ic)
      conv%maxGtotal(ic)=conv%maxGtotal(ic-1)+conv%maxG(ic)
    enddo

    ! DVF: Get the value of the matrix elements extrapolated with respect to bands
    ! idis is some sort of extrapolation "distance"

    idis = (ncondbands)/10
    if (idis.gt.0) then
      extrap_head = (ncondbands*conv%headtotal(ncondbands) - (ncondbands-idis)*conv%headtotal(ncondbands-idis)) / idis
      extrap_maxG = (ncondbands*conv%maxGtotal(ncondbands) - (ncondbands-idis)*conv%maxGtotal(ncondbands-idis)) / idis
    else
      extrap_head = ZERO
      extrap_maxG = ZERO
    endif

    write(strhead,701)1
    write(strtail,701)pol%nmtx

    if (pol%fullConvLog .eq. 0) then
      write(17,'(a,3e16.8,a,i0,a)') '# q=', pol%qpt(:,iq), '  (qpt ', iq, ')'
      write(17,'(a1,a7,4a20)') '#', 'ncbands', 'Re chi(0,0)', 'extrap', 'Re chi(Gmax,Gmax)', 'extrap'

      do ic = 1, ncondbands
        write(17,'(i8,4e20.8)') ic, dble(conv%headtotal(ic)), &
          dble(extrap_head), dble(conv%maxGtotal(ic)), dble(extrap_maxG)
      enddo
      write(17,*)
    elseif (pol%fullConvLog .eq. 1) then
      write(17,801) pol%qpt(:,iq), iq
      write(17,802)

      write(17,803)TRUNC(strhead)
      do ic = 1, ncondbands
        write(17,805) ic, dble(conv%headtotal(ic))
      enddo

      write(17,804)TRUNC(strtail)
      do ic = 1, ncondbands
        write(17,805) ic, dble(conv%maxGtotal(ic))
      enddo
    elseif (pol%fullConvLog .eq. 2) then
      write(17,801) pol%qpt(:,iq), iq
      write(17,802)

      write(17,803)TRUNC(strhead)
      do ic = 1, ncondbands
        write(17,805) ic, conv%headtotal(ic)
      enddo

      write(17,804)TRUNC(strtail)
      do ic = 1, ncondbands
        write(17,805) ic, conv%maxGtotal(ic)
      enddo
    endif

701 format(i16)
801 format('#',1x,'q =',3f10.6,1x,'iq =',i4)
802 format(2x,'nbands',1x,'chi0')
803 format(2x,'head ig = igp =',1x,a)
804 format(2x,'tail ig = igp =',1x,a)
805 format(i8,2e16.8)

    SAFE_DEALLOCATE(conv%headtotal)
    SAFE_DEALLOCATE(conv%maxGtotal)

    POP_SUB(chi_convergence_print)

  end subroutine chi_convergence_print

end module chi_convergence_m
