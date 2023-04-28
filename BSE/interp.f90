!===============================================================================
!
! MODULE:
!
! interp_m          Originally Dec/2012 by FHJ    Last Modified: Apr/2013 (FHJ)
!
!>  Bourne-again interpolation routines.
!
! DESCRIPTION:
!
!   This module contains routines to interpolate points. When initialized,
!   the code performs a Delaunay triangulation of the coarse points. When 
!   evaluating each interpolated point, the code finds the Delaunay simplex 
!   that encloses the point, and returns the appropriate coefficients.
!
!===============================================================================

#include "f_defs.h"

module interp_m

  use global_m
  use lapack_m
  use sort_m
  use tile_m

  implicit none

  !> The interpolator `object`
  type :: interp_t
    integer :: dims !< Number of dimensions
    integer :: npts_orig !< Number of points w/o periodicity
    integer :: npts !< Number of points including periodic images, if periodic
    !> Coarse points in "effective" crystal coords (this%dims, this%npts)
    !! "Effective" means that only the active dimensions are kept here, and
    !! they are selected/reordered according to dims_map(:).
    real(DP), pointer :: pts_crys(:,:)
    !> Coarse points in Cartesian coords (this%dims, this%npts)
    !! pts_cart = dot(bvec_, pts_crys). Only calculated for 2D and 3D.
    real(DP), pointer :: pts_cart(:,:)
    !> This is the effective transformation matrix from crystal to Cartesian,
    !! with dimensions this%dims x this%dims. Only calculated for 2D and 3D.
    real(DP), pointer :: bvec_(:,:)
    logical :: periodic !< If .true., add copy of pts in all directions
    !> Map each internal interpolated dimension to the crystal axes.
    !! This is needed b/c Qhull only accepts rank-1/-2 arrays for
    !! 1-/2-D interpolation.
    integer, pointer :: dims_map(:)
    !> For 1D systems only, this is a list of indices that sorts pts_crys.
    integer, pointer :: pts_crys_ind(:)
    !> For 1D systems, this is the sorted version of pts_crys, sorted. Same as
    !! pts_crys(1,pts_crys_ind(:)), but better for the cache.
    real(DP), pointer :: pts_crys_ord(:)
  end type interp_t

  public :: interp_init, interp_free, interp_eval

contains


  !> Initialize the interpolation object
  !! This routine:
  !! 1) Adds ghost copies of points in the periodic directions.
  !!    Reason: Qhull doesn`t deal with periodicity.
  !! 2) Finds an effective transformation from crystal coords to Cartesian,
  !!    so that the number of Cartesian coordinates = number of dimensions
  !!    This is trivial in 3D, and a bit tedious in 2D.
  !!    Reason: Qhull doesn`t accept "metrics", and requires #coordinates = 
  !!    dimensionality of the Delaunay tessellation.
  !! 3) Performs the Delaunay triangulation/tetrahedralization.
  subroutine interp_init(this, crys, pts, dims, periodic, active_dims)
    type(interp_t), intent(out) :: this
    type(crystal), intent(in) :: crys
    real(DP), intent(in) :: pts(:,:)
    integer, intent(in) :: dims
    logical, intent(in) :: periodic
    logical, intent(in) :: active_dims(3)

    integer :: ipt, ighost, n_ghost, jdim
    integer :: offset, ghost_G(3), ierr
    real(DP), pointer :: sub_pts(:,:)

    PUSH_SUB(interp_init)

    if (peinf%verb_debug .and. peinf%inode==0) then
      write(6,'(/1x,a)') 'Starting inter_init:'
      write(6,'(1x,a,i0)') '- dims = ', dims
      write(6,'(1x,a,l1)') '- periodic = ', periodic
      write(6,'(1x,a,3(1x,l1))') '- active_dims =', active_dims
      write(6,'(1x,a,i0," x ",i0)') '- shape(pts) = ', size(pts, dim=1), size(pts, dim=2)
    endif
    ! FHJ: Initialize number of dimensions, figure out dims_map, etc.
    if (dims<0 .or. dims>3) &
      call die("Wrong number of dimensions", only_root_writes=.true.)
    if (dims/=count(active_dims)) &
      call die("Inconsistent parameters dims and active_dims", &
      only_root_writes=.true.)
    this%dims = dims
    this%periodic = periodic
    this%npts_orig = size(pts, dim=2)
    if (this%npts_orig<this%dims+1) then
      ! If we have less points than dims+1, we can`t really interpolate...
      this%dims = 0
      if (peinf%inode==0) then
        write(0,'(/,a)') 'WARNING: Insufficient points to perform linear interpolation.'
        write(0,'(a)') 'Please, make sure your coarse k-grid is correct.'
        write(0,'(a,/)') 'Falling back to constant extrapolation using 1st coarse point.'
      endif
    endif
    if (this%dims==0) then
      ! No initialization needed for 0D systems or nearest neighbor interpolation
      this%npts = 0
      POP_SUB(interp_init)
      return
    endif
    SAFE_ALLOCATE(this%dims_map, (this%dims))
    call get_dims_map()

    ! FHJ: Deal with periodicity by adding ghost copies of the universe.
    ! Note: n_ghost=1 means only the original set of points
    n_ghost = 1
    if (this%periodic) n_ghost = 3**this%dims
    this%npts = this%npts_orig * n_ghost
    SAFE_ALLOCATE(this%pts_crys, (this%dims, this%npts))
    ghost_G(:) = 0
    do ighost = 0, n_ghost-1
      offset = this%npts_orig*ighost + 1
      ! FHJ: the outermost mod()`s map (0, 1, 2) to (0, 1, -1)
      ghost_G(1) = mod(mod(ighost, 3) + 1, 3) - 1
      if (this%dims>1) &
        ghost_G(2) = mod(mod(ighost/3, 3) + 1, 3) - 1
      if (this%dims>2) &
        ghost_G(3) = mod(mod(ighost/9, 3) + 1, 3) - 1

      sub_pts => this%pts_crys(1:this%dims, offset:offset+this%npts_orig-1)
      do ipt = 1, this%npts_orig
        do jdim = 1, this%dims
          sub_pts(jdim, ipt) = pts(this%dims_map(jdim), ipt) + ghost_G(jdim)
        enddo
      enddo
    enddo

    ! FHJ: For 1D, just sort the points in crystal coords. For 2D/3D, we
    ! need the effective metric and the Delaunay tessellation.
    if (dims==1) then ! 1D
      ! FHJ: Just sort the points.
      SAFE_ALLOCATE(this%pts_crys_ind, (this%npts))
      call sortrx(this%npts, this%pts_crys(1,:), this%pts_crys_ind)
      SAFE_ALLOCATE(this%pts_crys_ord, (this%npts))
      do ipt = 1, this%npts
        this%pts_crys_ord(ipt) = this%pts_crys(1, this%pts_crys_ind(ipt))
      enddo
    else ! 2D and 3D
      ! FHJ: Find the appropriate metric that transforms 3D/2D.
      ! crystal coordinates into Cartesian coordinates in the subspace.
      SAFE_ALLOCATE(this%bvec_, (this%dims, this%dims))
      if (dims==3) then
        call init_metric_3D()
      else if (dims==2) then
        call init_metric_2D()
      endif

      ! FHJ: Apply the transformation from the crystal coordinates to get the
      ! subspace Cartesian coordinates
      SAFE_ALLOCATE(this%pts_cart, (this%dims, this%npts))
      call dgemm('n', 'n', this%dims, this%npts, this%dims, 1D0, this%bvec_,&
        this%dims, this%pts_crys, this%dims, 0D0, this%pts_cart, this%dims)
    
      ! FHJ: Delaunay tessellation
      ierr = init_delaunay(this%pts_cart(:,:), this%npts, this%dims)
      if (ierr/=0) &
        call die('Could not perform Delaunay triangulation!', only_root_writes=.true.)
    endif

    POP_SUB(interp_init)

  contains

    !> Initialize the dims_map array with the periodic directions.
    !! E.g. #1: for a 2D system where z is truncated and x and y are 
    !! periodic, this%dims_map = (/1, 2/).
    !! E.g. #1: for a 1D system where z is periodic and x and y are 
    !! truncated, this%dims_map = (/3/).
    subroutine get_dims_map()
      integer :: jdim

      PUSH_SUB(interp_init.get_dims_map)

      if (this%dims==3) then
        this%dims_map(:) = (/1, 2, 3/)
      else if (this%dims==2) then
        if (.not.active_dims(3)) then
          this%dims_map = (/1, 2/)
        else if (.not.active_dims(2)) then
          this%dims_map = (/1, 3/)
        else
          this%dims_map = (/2, 3/)
        endif
      else
        do jdim = 1,3
          if (active_dims(jdim)) this%dims_map(1) = jdim
        enddo
      endif
      if (peinf%verb_debug .and. peinf%inode==0) then
        write(6,'(1x,a,3(1x,i1)/)') '- dims_map =', this%dims_map
      endif

      POP_SUB(interp_init.get_dims_map)
    end subroutine get_dims_map

    ! Initialize the metric for the 3D and 2D case.
    ! Note: we also slightly distort the metric so that the Delaunay
    ! triangulation is always the same for degenerate grids.
    !---------------------------------------------------------------------
    ! 3D: Effective vectors = crystal vectors
    !---------------------------------------------------------------------
    subroutine init_metric_3D()
      real(DP) :: delta

      PUSH_SUB(interp_init.init_metric_3D)

      this%bvec_(:,:) = crys%blat * crys%bvec(:,:)
      ! FHJ: Make a small perturbation to the off-diagonal terms of the
      ! transformation to make the Delaunay triangulation unique.
      delta = maxval(this%bvec_)*1d-10
      this%bvec_(1,2) = this%bvec_(1,2) + delta
      this%bvec_(1,3) = this%bvec_(1,3) + delta*sqrt(2d0)
      this%bvec_(2,3) = this%bvec_(2,3) + delta*sqrt(3d0)

      POP_SUB(interp_init.init_metric_3D)
    end subroutine init_metric_3D
    !---------------------------------------------------------------------
    ! 2D: Eff. vectors from Cholesky decomposition of leading 2x2 metric
    !---------------------------------------------------------------------
    subroutine init_metric_2D()
      real(DP) :: LT(2, 2)

      PUSH_SUB(interp_init.init_metric_2D)

      ! FHJ: Cholesky decomposition of the metric bdot = L * L^T
      LT(:,:) = crys%bdot(this%dims_map(1:2), this%dims_map(1:2))
      LT(1,1) = dsqrt(LT(1,1))
      LT(2,1) = 0d0
      LT(1,2) = LT(1,2)/LT(1,1)
      LT(2,2) = dsqrt(LT(2,2) - LT(1,2)**2)
      ! FHJ: Effective 2D transformation: V_eff = L^T
      this%bvec_(:,:) = LT(:,:)
      ! FHJ: Make a small perturbation to the off-diagonal term of the
      ! transformation to make the Delaunay triangulation unique.
      this%bvec_(1,2) = this%bvec_(1,2) + maxval(this%bvec_)*1d-10

      POP_SUB(interp_init.init_metric_2D)
    end subroutine init_metric_2D

  end subroutine interp_init


  !=========================================================================
  !=========================================================================


  !> Given a fine point `qq`, return the coarse points `ind_co` that enclose
  !! it and the coefficients `coefs_co` of the point in barycentric coordinates.
  subroutine interp_eval(this, qq, ind_co, coefs_co)
    type(interp_t), intent(in) :: this
    !> The fine point to be interpolated, in "regular" crystal coordinates.
    real(DP), intent(in) :: qq(:)
    !> Indices of the coarse point that form the simplex around `qq_fi`.
    integer, intent(out) :: ind_co(:)
    !> Coefficients of `qq_fi` in barycentric coordinates.
    real(DP), intent(out) :: coefs_co(:)

    real(DP) :: qq_crys(3), qq_cart(3), coefs_orig(4)
    real(DP) :: qq_co(3), dists(4)
    integer :: ind_orig(4), ierr, ipt, isort(4)

    PUSH_SUB(interp_eval)

    ind_co(:) = 1
    coefs_co(:) = 0d0
    if (this%dims==0) then ! 0D or nearest neighbor here
      coefs_co(1) = 1d0
      POP_SUB(interp_eval)
      return
    endif

    qq_crys(:) = 0d0
    qq_crys(1:this%dims) = qq(this%dims_map(:))
    if (this%dims==1) then ! 1D here
      call interp_1D()
      coefs_co = abs(coefs_co)
      POP_SUB(interp_eval)
      return
    endif

    ! FHJ: only 2D/3D from now on.
    qq_cart(:) = 0d0
    qq_cart(1:this%dims) = matmul(this%bvec_, qq_crys(1:this%dims))

    ! Find the simplex that encloses the fine point `qq_fi`.
    ! Notes: `coefs_co` is just the sorted version of `coefs_orig`, and 
    ! `ind_co` is `ind_orig` sorted + with points mapped back to 1st BZ.
    ierr = find_delaunay_simplex(qq_cart(1:this%dims), ind_orig, coefs_orig)
    if (ierr<0) then
      write(0,*) 'ERROR: could not find simplex:'
      write(0,'(a,3(f10.6,1x))') 'q(orig)=', qq(:)
      write(0,'(a,3(f10.6,1x))') 'q(crys)=', qq_crys(:)
      write(0,'(a,3(f10.6,1x))') 'q(cart)=', qq_cart(:)
      call die('Could not find Delaunay simplex.', only_root_writes=.true.)
    endif

    ! Sort: closest point should the first element of ind_co/coefs_co
    qq_co(:) = 0d0
    dists(:) = 0d0
    do ipt = 1, this%dims + 1
      qq_co(1:this%dims) = this%pts_cart(1:this%dims, ind_orig(ipt))
      dists(ipt) = sum( (qq_co(:)-qq_cart(:))**2 )
    enddo
    call sortrx(this%dims+1, dists, isort)
    do ipt = 1, this%dims + 1
      ! FHJ: This is how we map the ghost points to the original set
      ind_co(ipt) = mod(ind_orig(isort(ipt))-1, this%npts_orig) + 1
      coefs_co(ipt) = coefs_orig(isort(ipt))
    enddo

    coefs_co = abs(coefs_co)
    if (peinf%verb_max .and. peinf%inode==0) then
      write(6,'(1x,a,3(f10.6,1x))') 'q(orig)=', qq(:)
      write(6,'(1x,a,3(f10.6,1x))') 'q(crys)=', qq_crys(:)
      write(6,'(1x,a,3(f10.6,1x))') 'q(cart)=', qq_cart(:)
      write(6,'(1x,a,4(f7.3, 1x))') 'coefs=', coefs_co(:)
    endif

    POP_SUB(interp_eval)

  contains

    !> FHJ: For the 1D interpolation, just find the two neighboring points
    !! using a binary search algorithm, so we can claim we scale as O(log(N)).
    subroutine interp_1D()
      real(DP) :: dr
      integer :: imin, imid, imax
#ifdef DEBUG
      integer :: extrap_dir
#endif

      PUSH_SUB(interp_eval.interp_1D)

      imin = 1
      imax = this%npts
#ifdef DEBUG
      extrap_dir = 0
#endif

      if (this%pts_crys_ord(imax) < qq_crys(1) + TOL_SMALL) then
        ! FHJ: qq_crys is past the last point
        ind_co(1) = this%pts_crys_ind(this%npts)
        ind_co(2) = this%pts_crys_ind(this%npts-1)
        imax = this%npts
        imin = this%npts-1
        dr = this%pts_crys_ord(imax) - qq_crys(1)
#ifdef DEBUG
        extrap_dir = 1
#endif
      else if (this%pts_crys_ord(imin) > qq_crys(1) - TOL_SMALL) then
        ! FHJ: qq_crys is before the first point
        ind_co(1) = this%pts_crys_ind(1)
        ind_co(2) = this%pts_crys_ind(2)
        imax = 2
        imin = 1
        dr = qq_crys(1) - this%pts_crys_ord(imin)
#ifdef DEBUG
        extrap_dir = -1
#endif
      else
        ! FHJ: general case: use a binary search algorithm to locate one point
        ! larger and one smaller than the interpolant.
        do while (imin<imax-1)
          imid = (imin+imax)/2
          if (this%pts_crys_ord(imid) < qq_crys(1)) then
            imin = imid
          else
            imax = imid
          endif
        enddo
        if ( dabs(this%pts_crys_ord(imin)-qq_crys(1)) < &
             dabs(this%pts_crys_ord(imax)-qq_crys(1))) then
          dr = qq_crys(1) - this%pts_crys_ord(imin)
          ind_co(1) = this%pts_crys_ind(imin)
          ind_co(2) = this%pts_crys_ind(imax)
        else
          dr = this%pts_crys_ord(imax) - qq_crys(1)
          ind_co(1) = this%pts_crys_ind(imax)
          ind_co(2) = this%pts_crys_ind(imin)
        endif
      endif

      ind_co(1) = mod(ind_co(1)-1, this%npts_orig) + 1
      ind_co(2) = mod(ind_co(2)-1, this%npts_orig) + 1
      coefs_co(2) = dr / (this%pts_crys_ord(imax) - this%pts_crys_ord(imin))
      coefs_co(1) = 1d0 - coefs_co(2)
#ifdef DEBUG
      ! FHJ: Warn if we had to extrapolate the point, but don`t warn if the
      ! fine point is an almost "exact match".
      if (extrap_dir==1) then
        if (dabs(coefs_co(2))>TOL_SMALL) then
          write(0,'(a,f10.6,a)') 'WARNING: point ',qq_crys(1),' was extrapolated to positive values.'
        endif
      else if (extrap_dir==-1) then
        if (dabs(coefs_co(2))>TOL_SMALL) then
          write(0,'(a,f10.6,a)') 'WARNING: point ',qq_crys(1),' was extrapolated to negative values.'
        endif
      endif
#endif

      POP_SUB(interp_eval.interp_1D)
    end subroutine interp_1D

  end subroutine interp_eval

  !=========================================================================
  !=========================================================================


  !> Deallocate buffers used by interpolator object.
  subroutine interp_free(this)
    type(interp_t), intent(inout) :: this

    PUSH_SUB(interp_free)

    if (this%dims>0) then
      SAFE_DEALLOCATE_P(this%pts_crys)
      SAFE_DEALLOCATE_P(this%dims_map)
      if (this%dims==1) then
        SAFE_DEALLOCATE_P(this%pts_crys_ind)
        SAFE_DEALLOCATE_P(this%pts_crys_ord)
      else
        SAFE_DEALLOCATE_P(this%pts_cart)
        SAFE_DEALLOCATE_P(this%bvec_)
        if (free_delaunay()/=0) then
          call die("Error deallocating Delaunay triangulation.", only_root_writes=.true.)
        endif
      endif
    endif

    POP_SUB(interp_free)

  end subroutine interp_free

end module interp_m
