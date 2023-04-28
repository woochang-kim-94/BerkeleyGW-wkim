!===============================================================================
!
! Module:
!
! read_eqp_m             Written by  Tonatiuh Rangel         Last Modified 07/30/2014    
!
! Only contains objects to read eqp1.dat file
!
! Developers:
!   Tonatiuh Rangel, Berkeley CA, trangel@lbl.gov
!
!===============================================================================

#include "f_defs.h"

module read_eqp_m

 use global_m
 implicit none

!Object containing QPS corrections
 type, public:: qps_type
   integer :: number_of_spins
   integer :: number_of_k_points
   integer :: number_of_bands
   integer , allocatable :: index_eigenvalues(:,:)
   real(dp)   , allocatable :: eigenvalues_dft(:,:,:)
   real(dp)   , allocatable :: eigenvalues_eqp1(:,:,:)
   real(dp)   , allocatable :: k_points(:,:)
   complex(dp), allocatable :: rotation_matrix(:,:,:,:)
 end type qps_type
end module read_eqp_m
