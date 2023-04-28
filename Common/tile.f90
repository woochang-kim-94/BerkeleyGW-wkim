!===============================================================================
!
! MODULE:
!
! tile_m                                          Last Modified: Dec/2012 (FHJ)
!
!>  Tiling and tessellation routines for the R^2/R^3. Currently, only Delaunay 
!!  triangulation is available, but Voronoi tessellation should be added soon.
!
! DESCRIPTION:
!
!   This module is an abstraction layer that provides tessellation methods in
!   arbitrary dimensions. The "API" declared here (in the `contains` section)
!   should be generic, and the methods may be implemented by in external 
!   libraries. In this case, please write appropriate wrappers so that the
!   routines declared here do not depend on the technical details of the
!   external back end.
!
!   We currently use Qhull as the back end for the Delaunay tessellation.
!
!===============================================================================

#include "f_defs.h"

module tile_m

  use global_m
  implicit none

  public :: &
#ifdef USEVORO
    pts2bz, get_kpts_volumes, &
#endif
    init_delaunay, &
    find_delaunay_simplex, &
    get_num_simplices, &
    get_simplices, &
    get_neighbors, &
    init_second_neighbors, &
    get_max_num_second_neighbors, &
    find_delaunay_simplex_with_second_neighbors, &
    free_delaunay

  private
  
  interface

    !---------------------------------------------------------------------------
    ! FHJ: these tessellation functions are implemented in qhull/libtile_qhull.a

    integer function qhull_init_delaunay(points, num_points, dimensions, inode)
      use global_m
      implicit none
      real(DP), intent(in) :: points
      integer, intent(in) :: num_points
      integer, intent(in) :: dimensions
      integer, intent(in) :: inode
    end function qhull_init_delaunay

    !> Finds the Delaunay triangle/tetrahedron that encloses `point`.
    !!
    !! @param point [in] array (ndim).
    !! @param indices [out] array (ndim+1) of the (Fortran) indices of the vertices.
    !! @param coefs [out] coefficients of point in barycentric coordinates.
    integer function qhull_find_delaunay_simplex(point, indices, coefs)
      use global_m
      implicit none
      real(DP), intent(in) :: point
      integer, intent(out) :: indices
      real(DP), intent(out) :: coefs
    end function qhull_find_delaunay_simplex

    integer function qhull_get_num_simplices(num_simplices)
      implicit none
      integer, intent(out) :: num_simplices
    end function qhull_get_num_simplices

    integer function qhull_get_simplices(indices)
      implicit none
      integer, intent(out) :: indices
    end function qhull_get_simplices

    integer function qhull_get_neighbors(neighbors)
      implicit none
      integer, intent(out) :: neighbors
    end function qhull_get_neighbors

    integer function qhull_init_second_neighbors()
      implicit none
    end function qhull_init_second_neighbors

    integer function qhull_get_max_num_second_neighbors(n_nei)
      implicit none
      integer, intent(out) :: n_nei
    end function qhull_get_max_num_second_neighbors

    !> Computes second-nearest neighbors and barycentric coefficients for point `point`.
    !!
    !! The arrrays `indices` and `coefs` must be ndim+1 + `max_num_second_neighbors`
    !! big. The first ndim+1 points will correspond to the first neighbors.
    !! `nnei` is set to the total number of neighbors found, with
    !! `nnei` <= ndim + 1 + max_num_second_neighbors.
    integer function qhull_find_delaunay_simplex_with_second_neighbors(point, indices, coefs, nnei)
      use global_m
      implicit none
      real(DP), intent(in) :: point
      integer, intent(out) :: indices
      real(DP), intent(out) :: coefs
      integer, intent(out) :: nnei
    end function qhull_find_delaunay_simplex_with_second_neighbors

    integer function qhull_free_delaunay()
      implicit none
    end function qhull_free_delaunay

#ifdef USEVORO
    !---------------------------------------------------------------------------
    ! FHJ: these functions are implemented in voro/libtile_voro.a

    integer function voro_pts2bz(M, pts_cart, npts, umks)
      use global_m
      implicit none
      real(DP), intent(in) :: M
      real(DP), intent(in) :: pts_cart
      integer, intent(in) :: npts
      integer, intent(out) :: umks
    end function voro_pts2bz

    integer function voro_get_kpts_volumes(M, pts_cart, npts, vols)
      use global_m
      implicit none
      real(DP), intent(in) :: M
      real(DP), intent(in) :: pts_cart
      integer, intent(in) :: npts
      real(DP), intent(out) :: vols
    end function voro_get_kpts_volumes

    !---------------------------------------------------------------------------
#endif

  end interface

contains

  !> Calculates the Delaunay triangulation/tetrahedralization of a set of points.
  !!
  !! For performance issues, only one Delaunay triangulation can be active at
  !! a given time.
  !!
  !! \param points [inout] coarse points (dim, npts)
  !! \param num_points [in] number of points.
  !! \param dim [in] number of dimensions
  !! \return 0 for success
  integer function init_delaunay(points, num_points, dimensions)
    real(DP), intent(in) :: points(:,:)
    integer, intent(in) :: num_points
    integer, intent(in) :: dimensions
    
    PUSH_SUB(init_delaunay)
    init_delaunay = qhull_init_delaunay(points(1,1), num_points, dimensions, peinf%inode)
    POP_SUB(init_delaunay)
    return

  end function init_delaunay

  !> Finds the Delaunay triangle/tetrahedron that encloses point.
  !!
  !! \param point [in] real array (dim) containing the coordinates.
  !! \param indices [out] integer array (dim+1) with the indices of the vertices.
  !! \param coefs [out] coefficients of point in barycentric coordinates.
  integer function find_delaunay_simplex(point, indices, coefs)
    real(DP), intent(in) :: point(:)
    integer, intent(out) :: indices(:)
    real(DP), intent(out) :: coefs(:)

    PUSH_SUB(find_delaunay_simplex)
    find_delaunay_simplex = qhull_find_delaunay_simplex(point(1), indices(1), coefs(1))
    POP_SUB(find_delaunay_simplex)
    return

  end function find_delaunay_simplex

  !> Returns the total number of simplices obtained from the Delaunay
  !! triangulation.
  !!
  !! \param num_simplices [out] number of simplices.
  integer function get_num_simplices(num_simplices)
    integer, intent(out) :: num_simplices

    PUSH_SUB(get_num_simplices)
    get_num_simplices = qhull_get_num_simplices(num_simplices)
    POP_SUB(get_num_simplices)
    return

  end function get_num_simplices

  !> Returns the indices of the points that define each simplex.
  !!
  !! \param simplices [out] (ndims+1, num_simplices) Vertex ivert of
  !! simplex isimp corresponds to the point simplices(ivert,isimp) of
  !! the original list of points.
  integer function get_simplices(simplices)
    integer, intent(out) :: simplices(:,:)

    PUSH_SUB(get_simplices)
    get_simplices = qhull_get_simplices(simplices(1,1))
    POP_SUB(get_simplices)
    return

  end function get_simplices

  !> Returns the indices of the neighbors for each simplex.
  !!
  !! \param neighbors [out] (ndims+1, num_simplices) The neighbor ivert
  !! of simplex isimp is neighbors(ivert,isimp), and corresponds to the
  !! simplex reached by following the ridge opposite to vertex ivert.
  integer function get_neighbors(neighbors)
    integer, intent(out) :: neighbors(:,:)

    PUSH_SUB(get_neighbors)
    get_neighbors = qhull_get_neighbors(neighbors(1,1))
    POP_SUB(get_neighbors)
    return

  end function get_neighbors

  !> Initializes extra buffers required to compute second nearest neighbors
  !! You don`t need to do any extra procedure to free this functin, just
  !! call `free_delaunay` at the end as usual.
  integer function init_second_neighbors()

    PUSH_SUB(init_second_neighbors)
    init_second_neighbors = qhull_init_second_neighbors()
    POP_SUB(init_second_neighbors)
    return

  end function init_second_neighbors

  !> Get maximum number of second nearest neighbors
  integer function get_max_num_second_neighbors(n_nei)
    integer, intent(out) :: n_nei

    PUSH_SUB(get_max_num_second_neighbors)
    get_max_num_second_neighbors = qhull_get_max_num_second_neighbors(n_nei)
    POP_SUB(get_max_num_second_neighbors)
    return

  end function get_max_num_second_neighbors

  !> Finds the Delaunay simplex that encloses a point, plus second neighbors.
  !!
  !! \param point [in] real array (dim) containing the coordinates.
  !! \param indices [out] integer array (dim+1) with the indices of the vertices.
  !!        The first ndim+1 points are the first neighbors, and remaining are
  !!        second neighbors. The size of the array indices MUST BE that given
  !!        by get_max_num_second_neighbors().
  !! \param coefs [out] coefficients of point in barycentric coordinates.
  !!        The first ndim+1 points are the first neighbors and will have the
  !!        regular coefficients for barycentric interpolation. The remaining
  !!        coefficients are zero.
  !! \param nnei [out] total number of neighbors (first+second) around this
  !!        particular input point. The value of the array indices(:) will
  !!        be zero for points outside nnei.
  integer function find_delaunay_simplex_with_second_neighbors(point, indices, coefs, nnei)
    real(DP), intent(in) :: point(:)
    integer, intent(out) :: indices(:)
    real(DP), intent(out) :: coefs(:)
    integer, intent(out) :: nnei

    PUSH_SUB(find_delaunay_simplex_with_second_neighbors)
    find_delaunay_simplex_with_second_neighbors = &
      qhull_find_delaunay_simplex_with_second_neighbors(point(1), indices(1), coefs(1), nnei)
    POP_SUB(find_delaunay_simplex_with_second_neighbors)
    return

  end function find_delaunay_simplex_with_second_neighbors

  !> Frees buffers associated to the Delaunay triangulation.
  integer function free_delaunay()

    PUSH_SUB(free_delaunay)
    free_delaunay = qhull_free_delaunay()
    POP_SUB(free_delaunay)
    return

  end function free_delaunay

#ifdef USEVORO

  !> Transform points from crystal coordinates to "Cholesky Cartesian", i.e.,
  !! using lattice vectors derived from the Cholesky decomposition of the
  !! metric bdot. We can view this as a form of rotated Cartesian points with
  !! some nice properties, and it`s the kind of metric that voro++ uses.
  !!
  !! \param bdot [in] (3,3) reciprocal metric
  !! \params pts_crys [in] (3,npts) k-points in crystal coordinates
  !! \params pts_cart [out] (3,npts) k-points in "Cholesky Cartesian" coordinates
  !! \param M [out] (6) non-zero elements of Cholesky decomposition of the
  !!   metric
  subroutine pts2cholesky_cartesian(bdot, pts_crys, pts_cart, M)
    real(DP), intent(in) :: bdot(3,3)
    real(DP), intent(in) :: pts_crys(:,:)
    real(DP), intent(out) :: pts_cart(:,:)
    real(DP), intent(out) :: M(6)

    real(DP), parameter :: FX=1d0+1d-12, FY=1d0+1d-13, FZ=1d0+1d-14
    real(DP) :: xx, xy, xz, yy, yz, zz
    integer :: npts

    PUSH_SUB(pts2cholesky_cartesian)

    npts = size(pts_crys, 2)
    if (npts/=size(pts_cart,2)) then
      write(0,*) 'pts_crys: ', size(pts_crys, 2)
      write(0,*) 'pts_cart: ', size(pts_cart, 2)
      call die('Incompatible arrays in function pts2cholesky_cartesian', &
        only_root_writes=.true.)
    endif

    ! FHJ: Cholesky decomposition of the metric, required by voro++
    xx = sqrt(bdot(1,1))
    xy = bdot(2,1)/xx
    xz = bdot(3,1)/xx
    yy = sqrt(bdot(2,2)-xy*xy)
    yz = (bdot(3,2)-xy*xz)/yy
    zz = sqrt(bdot(3,3)-(xz*xz+yz*yz))
    M = (/xx, xy, yy, xz, yz, zz/)
#ifdef DEBUG
    if (peinf%verb_max) then
      write(6,*) 'BDOT'
      write(6,'(3(f12.9,1x))') bdot
      write(6,'(a,f12.9)') 'xx=',xx
      write(6,'(a,f12.9)') 'xy=',xy
      write(6,'(a,f12.9)') 'xz=',xz
      write(6,'(a,f12.9)') 'yy=',yy
      write(6,'(a,f12.9)') 'yz=',yz
      write(6,'(a,f12.9)') 'zz=',zz
      write(6,*) 'M'
      write(6,'(6(f12.9,1x))') M
    endif
#endif

    ! FHJ: We slightly shift the crystal coords to make the umklapp vectors unique
    pts_cart(1,:) = xx*(pts_crys(1,:)*FX) + xy*(pts_crys(2,:)*FY) + xz*(pts_crys(3,:)*FZ)
    pts_cart(2,:) = yy*(pts_crys(2,:)*FY) + yz*(pts_crys(3,:)*FZ)
    pts_cart(3,:) = zz*(pts_crys(3,:)*FZ)

    POP_SUB(pts2cholesky_cartesian)

  end subroutine pts2cholesky_cartesian


  !> Puts points into Brillouin Zone
  !!
  !! \param bdot [in] (3,3) reciprocal metric
  !! \param pts_crys [in] (3,npts) input points, in crystal coords
  !! \param umks [out] (3,npts) umklapp vectors, defined such that:
  !!   pts_crys_bz = pts_crys + umks
  subroutine pts2bz(bdot, pts_crys, umks)
    real(DP), intent(in) :: bdot(3,3)
    real(DP), intent(in) :: pts_crys(:,:)
    integer, intent(out) :: umks(:,:)

    real(DP) :: M(6)
    real(DP), allocatable :: pts_cart(:,:)
    integer :: npts, ierr

    PUSH_SUB(pts2bz)

    npts = size(pts_crys, 2)
    SAFE_ALLOCATE(pts_cart, (3,npts))
    call pts2cholesky_cartesian(bdot, pts_crys, pts_cart, M)
    ierr = voro_pts2bz(M(1), pts_cart(1,1), npts, umks(1,1))
    SAFE_DEALLOCATE(pts_cart)

    POP_SUB(pts2bz)

  end subroutine pts2bz

  !> Get volume associated to each k-point by via Voronoi diagram.
  !!
  !! \param bdot [in] (3,3) reciprocal metric
  !! \param pts_crys [in] (3,npts) input points, in crystal coords
  !! \param vols [out] (npts) volume for each k-point
  subroutine get_kpts_volumes(bdot, pts_crys, vols)
    real(DP), intent(in) :: bdot(3,3)
    real(DP), intent(in) :: pts_crys(:,:)
    real(DP), intent(out) :: vols(:)

    real(DP) :: M(6)
    real(DP), allocatable :: pts_cart(:,:)
    integer :: npts, ierr

    PUSH_SUB(get_kpts_volumes)

    npts = size(pts_crys, 2)
    SAFE_ALLOCATE(pts_cart, (3,npts))
    call pts2cholesky_cartesian(bdot, pts_crys, pts_cart, M)
    ierr = voro_get_kpts_volumes(M(1), pts_cart(1,1), npts, vols(1))
    SAFE_DEALLOCATE(pts_cart)

    POP_SUB(get_kpts_volumes)

  end subroutine get_kpts_volumes
#endif

end module tile_m
