!inspired by ABINIT/metric.f90

#include "f_defs.h"

module reciprocal_lattice_m

  use global_m
  use misc_m, only : invert_matrix
  use print_error_m, only : print_error
  implicit none

  private

  public :: &
    reciprocal_lattice

contains

subroutine reciprocal_lattice(gmet,gprimd,rmet,rprimd,gcvol,rcvol)
 logical, parameter :: debug=.false.
 real(DP), parameter :: tol12=0.000000000001_DP
 real(DP), parameter :: rad2deg=57.29577951308232d0

 real(DP),intent(out) :: rcvol,gcvol
 real(DP),intent(in) :: rprimd(3,3)
 real(DP),intent(out) :: gmet(3,3),gprimd(3,3),rmet(3,3)

 real(DP) :: angle(3),matrix_aux(3,3)

 PUSH_SUB(reciprocal_lattice)
! *************************************************************************

!Calculate real-space volume:
 call calculate_volume(rcvol,rprimd)

 if (abs(rcvol)<tol12) then 
   call print_error("reciprocal_lattice","the unit-cell volume is zero",1)
 end if

!Generates gprimd
 call invert_matrix(rprimd,matrix_aux)
 gprimd=TRANSPOSE(matrix_aux)

 call calculate_volume(gcvol,gprimd)

!Calculate metric tensor
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)   

!Calculate reciprocal-space metric tensor
 gmet = MATMUL(TRANSPOSE(gprimd),gprimd)

!Write out the angles
 if (debug) then
   angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))*rad2deg 
   angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))*rad2deg
   angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))*rad2deg
   write(*, '(a,3f16.8,a)' )' Angles (23,13,12)=',angle(1:3),' degrees'
 end if

 POP_SUB(reciprocal_lattice)
contains

subroutine calculate_volume(volume,aa)
  real(DP),intent(in)::aa(3,3)
  real(DP),intent(out)::volume
 
  PUSH_SUB(reciprocal_lattice.calculate_volume)


  volume=aa(1,1)*(aa(2,2)*aa(3,3)-aa(3,2)*aa(2,3))+&
&        aa(2,1)*(aa(3,2)*aa(1,3)-aa(1,2)*aa(3,3))+&
&        aa(3,1)*(aa(1,2)*aa(2,3)-aa(2,2)*aa(1,3))

  POP_SUB(reciprocal_lattice.calculate_volume)

end subroutine calculate_volume


end subroutine reciprocal_lattice
!!***


end module reciprocal_lattice_m
