!==============================================================================
!
!  test_delaunay_f        Originally by FHJ     Last Modified 12/04/2012 (FHJ)
!
!    A simple Fortran program that tests the BerkeleyGW bindings for Qhull.
!
!    This file is part of the BerkeleyGW package.
!
!==============================================================================

#include "f_defs.h"

program test_delaunay_f

  use global_m
  use tile_m
  implicit none

  integer ::  nnei, nnei_max, dim_, npts, npts_orig, nrep
  integer :: ip, idim_, div, rem
  integer, allocatable ::  indices(:), indices2(:)
  real(DP), allocatable :: points(:,:), pt(:), coefs(:), coefs2(:)
  logical :: show_coords = .true.;
  integer :: ierr

  PUSH_SUB(test_delaunay_f)

  nrep = 1000
  write(6,'(a)') "Testing libtile_qhull - Fortran bindings."
  call open_file(unit=10, file='points.dat', form='formatted', status='old')
  read(10,*) dim_, npts, npts_orig
  SAFE_ALLOCATE(points, (dim_, npts))
  read(10,*) points
  call close_file(10)

  write(6,'(a,i0)') "Dimensions: ", dim_
  write(6,'(a,i0)') "Number of points: ", npts
  write(6,'(a,i0)') "Number of points w/o ghost points: ", npts_orig

  if (show_coords) then
    write(6,'(/,a)') "Input coordinates:";
    do ip = 1, npts
      call print_vect(6, ip, points(:,ip), dim_)
    enddo
  endif

  SAFE_ALLOCATE(indices, (dim_+1))
  SAFE_ALLOCATE(coefs, (dim_+1))
  SAFE_ALLOCATE(pt, (dim_))
  pt(:) = 0
  pt(1) = 0.125d0
  pt(2) = 0.75d0
  write(6,'(/,a)') "Test coordinate:"
  call print_vect(6, 0, pt, dim_)
  write(6,'()')

  write(6,'(a)') "Calling init_delaunay."
  ierr = init_delaunay(points, npts, dim_)

  write(6,'(/a)') "Finding max number of second neighbors:"
  FLUSH(0)
  FLUSH(6)
  ierr = init_second_neighbors()
  ierr = get_max_num_second_neighbors(nnei_max)
  write(6,'(i0/)') nnei_max

  SAFE_ALLOCATE(indices2, (dim_ + 1 + nnei_max))
  SAFE_ALLOCATE(coefs2, (dim_ + 1 + nnei_max))
  write(6,'(a)') 'Calling find_delaunay_simplex_with_second_neighbors.'
  ierr = find_delaunay_simplex_with_second_neighbors(pt, indices2, coefs2, nnei);
  write(6,'(/,a)') 'Found vertices:'
  do ip = 1, nnei
    call print_vect(6, indices2(ip), points(:,indices2(ip)), dim_)
  enddo
  write(6,'(/,a)',advance='no') "Coefficients: "
  call print_vect(6, -1, coefs2, nnei)

  write(6,'(/a)') 'Writing output data to file ''points_out_f.dat''.'
  ! Write header
  call open_file(unit=10, file='points_out_f.dat', form='formatted', status='replace')
  write(10,'(4(i0,1x))') dim_, npts, npts_orig, nnei
  ! Input coordinates
  do ip = 1, npts
    call print_vect(10, ip, points(:,ip), dim_)
  enddo
  write(10,*)
  ! First + second neighbors
  do ip = 1, nnei
    call print_vect(10, indices2(ip), points(:,indices2(ip)), dim_)
  enddo
  write(10,*)
  ! Original point
  call print_vect(10, 0, pt, dim_)
  call close_file(10)

  write(6,'()')
  write(6,'(a,i8,a)') "Calling find_delaunay_simplex ", nrep, " times."
  do ip = 1, nrep
    ierr = find_delaunay_simplex(pt, indices, coefs)
  enddo
  write(6,'(a)') "Calling free_delaunay."
  ierr = free_delaunay()

  write(6,'(/,a)') "Found vertices:"
  do ip = 1, dim_+1
    call print_vect(6, indices(ip), points(:,indices(ip)), dim_)
  enddo
  write(6,'(/,a)',advance='no') "Coefficients: "
  call print_vect(6, -1, coefs, dim_+1)

  SAFE_DEALLOCATE(indices2)
  SAFE_DEALLOCATE(coefs2)
  SAFE_DEALLOCATE(indices)
  SAFE_DEALLOCATE(coefs)
  SAFE_DEALLOCATE(pt)
  SAFE_DEALLOCATE(points)
  write(6,'(/,a)') "All Done!"

  POP_SUB(test_delaunay_f)

  contains

    subroutine print_vect(iunit, id, vec, dims)
      integer, intent(in) :: iunit, id
      real(DP), intent(in) :: vec(:)
      integer, intent(in) :: dims

      integer :: jdim_

      PUSH_SUB(print_vect)

      if (id>=0) write(iunit,'(1x,i4,1x)',advance='no') id
      do jdim_ = 1, dims
        write(iunit,'(f6.3,1x)',advance='no') vec(jdim_)
      enddo
      write(iunit,*)

      POP_SUB(print_vect)

    end subroutine print_vect

end program test_delaunay_f
