#include "f_defs.h"

!*************************************************************************
Module cgParamsMod
!*************************************************************************

  use sysparams, only : double
  use push_pop_m

  use TagHandlerMod

  Type cgParamsT
     integer          :: period,switch
     real(double)     :: tol
  end Type cgParamsT

  contains

!*************************************************************************
  Subroutine cgParamsInit(Params,TagHandler)
!*************************************************************************

    type(TagHandlerT), pointer :: TagHandler
    type(cgParamsT) :: Params    
    
    integer :: error,error1
    logical :: default

    PUSH_SUB(cgParamsInit)
    
    ! Fill in the tolerance in Params%tol
    default = .True.
    Call FindTag(TagHandler,'CGTolerance',error1)
    if(error1.eq.0) then
       read(unit=tagHandler%fileno, fmt = *, iostat=error ) & 
            & Params%tol
       if(error.eq.0 .and. Params%tol.ne.0 ) default = .False.
    end if
    if(default) then
       Params%tol = 1.d-5
       print *, 'Using the default value of CGTolerance',Params%tol
    end if

    ! Fill in the tolerance in Params%switch
    default = .True.
    Call FindTag(TagHandler,'DiagonalizationSwitch',error1)
    if(error1.eq.0) then
       read(unit=tagHandler%fileno, fmt = '(I10)', iostat=error ) & 
            & Params%Switch
       if(error.eq.0 .and. Params%Switch.le.2) &
            & default = .False.
    end if
    if(default) then
       Params%Switch = 0
       print *, 'Using the Lapack routine for Diagonalisation'
    end if

    if(Params%Switch.eq.1) then

       ! Fill in the tolerance in Params%period
       default = .True.
       Call FindTag(TagHandler,'CGIterationPeriod',error1)
       if(error1.eq.0) then
          read(unit=tagHandler%fileno, fmt = '(I10)', iostat=error ) & 
               & Params%period
          if(error.eq.0 .and. Params%period.gt.0 ) default = .False.
       end if
       if(default) then
          Params%period = 3
          print *, 'Using the default value of CGIterationPeriod',Params%period
       end if

       ! These parameters are never used. Why were they here?? --DAS
       ! Fill in the preconditioning in Params%precond
!       default = .True.
!       Call FindTag(TagHandler,'CGPreConditioner',error1)
!       if(error1.eq.0) then
!          read(unit=tagHandler%fileno, fmt = '(I10)', iostat=error ) & 
!               & Params%precond
!          if(error.eq.0 .and. (Params%precond.eq.0 .or. Params%precond.eq.1)) default = .False.
!       end if
!       if(default) then
!          Params%precond = 0
!          print *, 'Using the default Preconditioner',Params%precond
!       end if

       ! Fill in the convergence switch in Params%conv
!       default = .True.
!       Call FindTag(TagHandler,'CGConvergenceSwitch',error1)
!       if(error1.eq.0) then
!          read(unit=tagHandler%fileno, fmt = '(I10)', iostat=error ) & 
!               & Params%conv
!          if(error.eq.0 .and. (Params%conv.eq.0 .or. Params%conv.eq.1)) default = .False.
!       end if
!       if(default) then
!          Params%conv = 0
!          print *, 'Not testing for convergence'
!       end if

    end if

    POP_SUB(cgParamsInit)

  end Subroutine cgParamsInit

end Module cgParamsMod

