!=============================================================================
!
! Routines:
!
! (1) intpts_local()    Originally By JRD        Last Modified: 6/16/2011 (FHJ)
!  
!     Calculates the interpolating coefficients for a fine point qq given a
!     coarse grid. The function is optimized for a large number of points
!     using a local search algorithm. The routine only looks for coarse points
!     that are in neighborhood of the fine point, which is accomplished by
!     dividing the space into cells. This subroutine might is especially useful
!     if the coarse grid contains more than ~1000 points.
!
!     - Note 1: before using the function, you must prepare the cell structure
!     via alloc_intpts(), and later decommission it via dealloc_intpts().
!     - Note 2: results are mostly the same as those obtained by the original
!     function intpts_full(). But there might be some reproducible differences
!     in unshifted grids due to the way the points are sorted.
!
!     Changelog:
!       2010/28/06 [FHJ] Mostly rewritten to make control flow more logical.
!
!     intpts_full()     Originally By JRD        Last Modified: 6/16/2011 (FHJ)
!
!     Original "non-local" version of intpts written by JRD.
!
!     Changelog:
!       2010/14/06 [FHJ] Fixed the interpolation near the border of the BZ
!                  when umklapp = false (amat should be the smallest distance)
!     This routine was decommissioned at r3253.
!     We used to have ESSL support here. It could be put back if desired.
!
! (2) alloc_intpts()    Originally by FHJ        Last Modified: 6/16/2011 (FHJ)
!      
!     Allocates memory and initializes cell structure for intpts_local. You
!     should decommission and reinitialize the cells for each kind of
!     interpolation that you perform.
!
! (3) dealloc_intpts()  Originally by FHJ        Last Modified: 6/16/2011 (FHJ)
!     Free memory allocated by alloc_intpts().
!
! (4) get_ndims()       Originally by FHJ        Last Modified: 2/08/2012 (FHJ)
!
!     
!==============================================================================

#include "f_defs.h"

module intpts_m

  use global_m
  use lapack_m
  use sort_m

  implicit none

  public :: alloc_intpts, dealloc_intpts, intpts_local, get_ndims

  private

  real(DP), allocatable :: cell_head(:,:,:), cell_list(:)
! FHJ: cell_factor = 1/length of the individual cell
!      cell_shift = length of the individual cell / 2
  real(DP), dimension(3) :: cell_dmin, cell_dmax, cell_factor, cell_shift
  integer , dimension(3) :: cell_N
  !> size of the list that keep the nearest neighbors.
  integer, parameter :: NEI_BUFFER = 1000
  !> These variables represent what alloc_intpts *thinks* that the space
  !! (geometry/#dims) looks like
  integer :: cell_ndims
  logical :: cell_active_dim(3) !< .true. if system is non-deg. in a particular dim

contains

!> FHJ: interpolate a point qq using a coarse grid coarsepts. Call alloc_intpts
!!      before calling this function, and dealloc_intpts after.
  subroutine intpts_local(crys,qq,ncoarse,coarsepts,xct,closepts,closeweights,periodic)
    type (crystal), intent(in) :: crys
    real(DP), intent(in) :: qq(3)
    integer, intent(in) :: ncoarse
    real(DP), intent(in) :: coarsepts(3, ncoarse)
    type (xctinfo), intent(in) :: xct
    integer, intent(out) :: closepts(4)
    real(DP), intent(out) :: closeweights(4)
    logical, intent(in) :: periodic
    
    real(DP) :: qq_(3)
    real(DP) :: norm,den,closediffs(4)
    integer ii,jj,info,indx(4),iq,iii,jjj,kkk,kkk_max,iqarray(2)
    real(DP) :: diff(3),vol_sum,delta_r
    integer :: ipiv(4)
    real(DP) :: amat(4,4),vol(3,3),amatt(3,3),bvec(4,4)
    real(DP) :: fct1,fct2,fct3
! FHJ: The size of arrays like iqnear are now fixed. This way, we make use of
!      the stack memory, which is faster than the heap. This is important b/c
!      this function is called many times.
    real(DP) :: diffl(NEI_BUFFER), tempr(NEI_BUFFER)
    integer :: iqnear(NEI_BUFFER), idx_found(NEI_BUFFER), tempi(NEI_BUFFER)
! FHJ: We should find at least this number of neighbors for a given dimension
    integer :: min_pts

    integer :: nn_found, cell_idx(3), cell_min(3), cell_max(3)
    integer :: cell_dist, cell_dist_max
    integer :: j1,j2,j3, i1,i2,i3

    PUSH_SUB(intpts_local)

    closeweights(:) = 0d0

! FHJ: number of neighbors to find
    min_pts = xct%idimensions + 1

! FHJ: The qq_ is a vector inside the cubic reciprocal cell.
    qq_(:) = qq(:)
    if (periodic) then
      ! FHJ: this is guaranteed to be in [0, 1)
      qq_(:) = UNIT_RANGE(qq_(:))
    endif

! FHJ: Find the central cell for this qq_
    call get_cell_idx(qq_, cell_idx, periodic)

! FHJ: Initially, we only consider the +/- 1 neighboring cells, and if the cell
!      partition is good, the routine will return when cell_dist=1.
!      If we can`t find the minimum number of neighbors, include cells farther
!      away, up to cell_dist_max apart.
    cell_dist_max = idint(maxval(cell_N)*0.5 + TOL_SMALL)
    if (cell_dist_max<1) cell_dist_max=1
cell_dist_do:&
    do cell_dist = 1,cell_dist_max
      iqnear=0
      indx=0

! FHJ: Cell_min and cell_max are the starting/ending local neighboring cells
!      But start/end cells should not overlap!
      do ii=1,3
        if ((2*cell_dist+1).ge.cell_N(ii)) then
          cell_min(ii) = 1
          cell_max(ii) = cell_N(ii)
        else
          cell_min(ii) = cell_idx(ii) - cell_dist
          cell_max(ii) = cell_idx(ii) + cell_dist
        endif
      enddo
      if (.not.periodic) then
        cell_min(:) = max(cell_min(:), 1)
        cell_max(:) = min(cell_max(:), cell_N)
      endif

! FHJ: Search for points in neighboring cells.
      nn_found=0
do_:  do i1=cell_min(1),cell_max(1)
        j1 = fix_index(i1,1)
        do i2=cell_min(2),cell_max(2)
          j2 = fix_index(i2,2)
          do i3=cell_min(3),cell_max(3)
            j3 = fix_index(i3,3)
            call get_cell_pts()
            if (nn_found >= nei_buffer) exit do_
          enddo
        enddo
      enddo do_

! FHJ: Require a minimum number of pts, which depend on the # of dimensions
      if (nn_found < min_pts) then
        if (peinf%inode==0) then
          write(0,'(a,i1,a)') 'WARNING: Could not find ',min_pts,' points.'
          write(0,*) 'Increasing the number of neighboring cells'
        endif
        cycle cell_dist_do !continue to the next value of cell_dist
      endif

! FHJ: Sort the distances and break ties by the index on the coarse grid.
!      This should mimic the behavior of the original code, so that results are
!      almost compatible with the non-local routine. Deviations are due
!      to rounding errors, and b/c of the stability of the sorting algorithm.
!      Instead of doing one search with the two conditions, we break it into
!      2 searches (ideally, we should make the second search stable!)
!   1) Sort by indices and rearrange vectors.
      call sortix(nn_found, idx_found, iqnear)
      ! DAS: one shot without the temp arrays may be optimized incorrectly
      tempi(1:nn_found) = idx_found(iqnear(1:nn_found))
      tempr(1:nn_found) = diffl(iqnear(1:nn_found))
      idx_found(1:nn_found) = tempi(1:nn_found)
      diffl(1:nn_found) = tempr(1:nn_found)
!   2) Sort by distances.
      call sortrx(nn_found, diffl, iqnear)

! FHJ: iqnear stores the indices of the (nn_found) *nearest* q-points. Also,
!      idx_co stores the *coarse* point corresponding to a neighbor. Example:
!      to get the index of the *closest coarse* pt, use idx_co(iqnear(1))

! FHJ: For the interpolation scheme, we have to consider now triangles or
!      tetrahedrons. However, the points may be degenerate. So, we have
!      to test many combinations of 4 points. Obviously, we can always fix the
!      first coordinate idx_co(iqnear(1)) and loop only through the other ones.

! FHJ: Now, instead of looping through coordinates, we can loop through all the
!      possible combinations of indices. We define:
!      - indx(:): one particular combination of indices for the iqnear list of
!        near points. For 3D, some possible values are
!          {indx(:)} = {(1,2,3,4), (1,2,3,5), ..., (1,N-2,N-1,N)}
!        where N=nn_found. For 2D:
!          {indx(:)} = {(1,2,3), (1,2,4), ..., (1,N-1,N)}
!      - closepts(i) = idx_found(indx(i)): index of the coarse pt associated to
!        the ith-vertex of the tetrahedron.
!      - closediffs(i) = diffl(indx(i)): distance from the ith-vertex to the
!        point qq. 
!     Note: indx(1), closepts(1) and closediffs(1) are always fixed
      indx(1) = iqnear(1)
      closepts(1) = idx_found(indx(1))
      closediffs(1) = diffl(indx(1))

! FHJ: Before doing anything too fancy, are we within the numerical error of 
!     a point in the coarse grid?
!------------------------
      if (closediffs(1)<TOL_Small) then
        call return_first_pt()
        POP_SUB(intpts_local)
        return
      endif    

! FHJ: Special case: 1 dimension (just do a simple linear interpolation)
!------------------------
      if (xct%idimensions==1) then
        closepts(2) = idx_found(iqnear(2))
        closediffs(2) = diffl(iqnear(2))

        closeweights(1) = closediffs(2)/(closediffs(1)+closediffs(2))
        closeweights(2) = closediffs(1)/(closediffs(1)+closediffs(2))
! GKA: Avoid small negative numbers
        closeweights = abs(closeweights)
        POP_SUB(intpts_local)
        return
      endif

! FHJ: For 2D, if nn_found==3, we have to make up a third point to avoid
!       memory problems. Either this or writing lots of if`s later...
      if ((xct%idimensions==2).and.(nn_found==3)) iqnear(4) = iqnear(3)

! FHJ: 2 and 3 dimensions
!------------------------
      do iii = 2,nn_found-1
        do jjj = iii+1,nn_found

! FHJ: If there are only 2 dimensions, we don`t have to loop through kkk.
!      The simplest solution is to make kkk constant. Note that, for 3D,
!      kkk_max automatically sets the maximum values of iii and jjj to
!      nn_found-2 and nn_found-1, respectively.
          if (xct%idimensions==3) then
            kkk_max = nn_found
          else
            kkk_max = jjj+1
          endif
          do kkk = jjj+1,kkk_max
            indx(2) = iqnear(iii)
            indx(3) = iqnear(jjj)
            indx(4) = iqnear(kkk) !not used for 2D
            closepts(2:4) = idx_found(indx(2:4))
            closediffs(2:4) = diffl(indx(2:4))

! FHJ: amat(ii,:) is the coordinate of the ii-th vertex of the tetrah./triang.
            do ii=1,4
              do jj=1,3
                amat(ii,jj)=coarsepts(jj,closepts(ii))
! FHJ: If periodic, then center amat around the qq point
                if (periodic) then
                  delta_r = amat(ii,jj) - qq(jj)
                  amat(ii,jj) = MIN_RANGE(delta_r) + qq(jj)
                endif
              enddo
              amat(ii,4) = 1.d0
            enddo

! FHJ: Check if the four/three vertices form a tetrahedron/triangle with non-zero volume
            if (xct%idimensions==3) then
              do ii=1,3 !note: we are indeed using all the 4 vertices
                norm=0.d0
                do jj=1,3
                  vol(ii,jj) = amat(ii+1,jj) - amat(1,jj)
                  norm = norm + vol(ii,jj)**2
                enddo
                norm=sqrt(norm) !normalize the vectors 
                if(abs(norm) < TOL_Zero) call die("intpts polyhedron has zero norm")
                do jj=1,3
                  vol(ii,jj) = vol(ii,jj)/norm
                enddo
              enddo
              vol_sum = vol(1,1)*vol(2,2)*vol(3,3) + vol(1,2)*vol(2,3)*vol(3,1) + &
                vol(1,3)*vol(2,1)*vol(3,2) - vol(1,1)*vol(2,3)*vol(3,2) - &
                vol(1,2)*vol(2,1)*vol(3,3) - vol(1,3)*vol(2,2)*vol(3,1)
              vol_sum = dabs(vol_sum)
            else !2D case
              do ii=1,2
                norm=0.d0
                do jj=1,3
                vol(ii,jj) = amat(ii+1,jj) - amat(1,jj)
                norm = norm + vol(ii,jj)**2
                enddo

                norm=sqrt(norm) ! FHJ: Normalize the vectors 
                do jj=1,3
                  vol(ii,jj) = vol(ii,jj)/norm
                enddo
              enddo
! JRD: This does a cross product
              fct1 = vol(1,2)*vol(2,3) - vol(2,2)*vol(1,3)
              fct2 = vol(1,3)*vol(2,1) - vol(1,1)*vol(2,3)
              fct3 = vol(1,1)*vol(2,2) - vol(1,2)*vol(2,1)
              vol_sum = sqrt(fct1**2+fct2**2+fct3**2)
            endif !3D/2D cases

            if (vol_sum > TOL_Small) then! this is a valid tetrahedron/triangle
              if (xct%idimensions==3) then
                bvec(:,:)=0D0
                bvec(1,1)=1D0
                bvec(2,2)=1D0
                bvec(3,3)=1D0
                bvec(4,4)=1D0
                call DGESV(4,4,amat,4,ipiv,bvec,4,info)

                if (info.ne.0) then
                  write(0,*) 'WARNING: LAPACK failed at intpts_local, which is unexpected. Cycling indices.'
                  cycle
                else
                  do i1 = 1,4
                    closeweights(i1) = bvec(1,i1)*qq(1) + &
                      bvec(2,i1)*qq(2)+bvec(3,i1)*qq(3) + bvec(4,i1)
                  enddo
                endif
              else !2D
                do ii=1,3
                  iq=0
                  do jj=1,3
                    if (xct%is_periodic(jj)) then
                      iq=iq+1
                      iqarray(iq)=jj
                      amatt(ii,iq)=amat(ii,jj)
                    endif
                  enddo
                  amatt(ii,3)=1.d0
                enddo

                den=-1.d0*amatt(2,1)*amatt(1,2)+ &
                  amatt(3,1)*amatt(1,2)+ &
                  amatt(1,1)*amatt(2,2)- &
                  amatt(3,1)*amatt(2,2)- &
                  amatt(1,1)*amatt(3,2)+ &
                  amatt(2,1)*amatt(3,2)

                closeweights(1) = qq(iqarray(1))*(amatt(2,2)-amatt(3,2))/den + &
                            qq(iqarray(2))*(amatt(3,1)-amatt(2,1))/den + &
                            (amatt(2,1)*amatt(3,2)-amatt(3,1)*amatt(2,2))/den
                closeweights(2) = qq(iqarray(1))*(amatt(3,2)-amatt(1,2))/den + &
                            qq(iqarray(2))*(amatt(1,1)-amatt(3,1))/den + &
                            (amatt(3,1)*amatt(1,2)-amatt(1,1)*amatt(3,2))/den
                closeweights(3) = qq(iqarray(1))*(amatt(1,2)-amatt(2,2))/den + &
                            qq(iqarray(2))*(amatt(2,1)-amatt(1,1))/den + &
                            (amatt(1,1)*amatt(2,2)-amatt(2,1)*amatt(1,2))/den
              endif !3D/2D
              
! GKA: Avoid small negative numbers
              closeweights = abs(closeweights)

              POP_SUB(intpts_local)
              return !since we got a valid set of vertices
            endif !sum_vol>TOL_Small (valid vertices)
            
          enddo !kkk
        enddo !jjj
      enddo !iii

      if (peinf%inode==0) then
        write(0,*)
        write(0,'(/,a,i3,a)') ' WARNING: could not perform the interpolation using ',&
          cell_dist,'-order cell neighbors.'
        if (cell_dist<cell_dist_max) then
          write(0,'(/,a,i3,a)') ' Increasing to ',cell_dist+1,'-order neighbors.'
        else
          write(0,*) ' Interpolation failed! Using nearest point available.'
          write(0,*) ' Check if the coarse grid is well-behaved and if the cell'
          write(0,*) ' population is evenly distributed.'
        endif
      endif
    enddo cell_dist_do !the loop over cell_dist

! FHJ: Epic Fail! Just return the nearest-neighbor.
    call return_first_pt()
    POP_SUB(intpts_local)
    return

!------------------------------------------------------------------------------

  contains
  
    subroutine return_first_pt()
      PUSH_SUB(intpts_local.return_first_pt)
! FHJ: closepts(1) was already set, and memory for 2:4 has to be accessible
      closepts(2:4) = 1
      closeweights(1) = 1.0D0
      closeweights(2:4) = 0D0
      POP_SUB(intpts_local.return_first_pt)
      return
    end subroutine return_first_pt
  
    subroutine get_cell_pts()

      integer :: idx_co, ikb     

      ! no push/pop since called too frequently
      idx_co=cell_head(j1,j2,j3)
      do while (idx_co>0)
        nn_found = nn_found+1

        ! FHJ: we don`t need to punish b/c the coarse grid, by construction,
        ! doesn`t have points in non-periodic directions.
        ! TODO: die if WFN and epsilon grids are incompatible.
        diff(:) = qq_(:) - coarsepts(:, idx_co)
        if (periodic) then
          diff(:) = MIN_RANGE(diff(:))
        endif
        diffl(nn_found)=DOT_PRODUCT(diff,MATMUL(crys%bdot,diff))
        if (xct%idimensions==1) then
          diffl(nn_found)=sqrt(diffl(nn_found))
        endif
        idx_found(nn_found) = idx_co
        ! FHJ - Avoid buffer overflow
        if (nn_found >= NEI_BUFFER) return ! no push/pop since called too frequently
        ! FHJ - Move to next point in this cell
        idx_co=cell_list(idx_co)
      enddo

      ! no push/pop since called too frequently
    end subroutine get_cell_pts
          
    !Returns n such that 1 <= n <= N
    integer function fix_index(idx, dim_)
      integer, intent(in) :: idx, dim_

      ! no push/pop since called too frequently
      fix_index = idx
      if (periodic) then
        fix_index = modulo(fix_index-1, cell_N(dim_)) + 1
      endif
    end function fix_index
    
  end subroutine intpts_local

!------------------------------------------------------------------------------   
!> Calculate the kgrid geometry by looking at the space spawned by the kpoints
!! This is useful when the value passed as kgrid is (0,0,0), when doing
!! inteqp, for instance. This function actually works by allocating and 
!! deallocating a cell structure, so make sure there is no active cell
!! before calling this function!
  subroutine get_ndims(kg, xct)
    type (grid), intent(in) :: kg
    type (xctinfo), intent(inout) :: xct
    integer :: ii

    PUSH_SUB(get_ndims)

    call alloc_intpts(kg%nf, kg%f, periodic=.true.)
    xct%idimensions = cell_ndims
    do ii=1,3
      xct%is_periodic(ii) = .false.
      if (cell_active_dim(ii)) xct%is_periodic(ii) = .true.
    enddo
    call dealloc_intpts()

    POP_SUB(get_ndims)

  end subroutine get_ndims

!------------------------------------------------------------------------------   
  subroutine alloc_intpts(ncoarse, coarsepts_, periodic)
    integer, intent(in) :: ncoarse
    real(DP), intent(in) :: coarsepts_(3,ncoarse)
    logical, intent(in) :: periodic

! FHJ: These are the coarsepts_ mapped to the cubic cell (if periodic)
    real(DP) :: coarsepts(3,ncoarse)
    integer :: j1,j2,j3, ii, idx_co,cell_idx(3)
    real(DP) :: dim_vol, prop_const
    logical :: should_write

    PUSH_SUB(alloc_intpts)

    should_write = peinf%verb_debug

    if (peinf%inode==0.and.should_write) write(6,'(/,a)') 'Starting alloc_intpts'

    if (periodic) then
      ! FHJ: this is guaranteed to be in [0, 1)
      coarsepts(:,:) = UNIT_RANGE(coarsepts_(:,:))
    else
      coarsepts(:,:) = coarsepts_(:,:)
    endif
   
! FHJ: Avoid problems that a very small dimension may cause
    cell_active_dim(:) = .false.
    cell_ndims = 0
    cell_dmin(:) = minval(coarsepts(:,:), dim=2)
    cell_dmax(:) = maxval(coarsepts(:,:), dim=2)
    do ii=1,3
      if ((cell_dmax(ii)-cell_dmin(ii)) > TOL_Small) then
        cell_active_dim(ii) = .true.
        cell_ndims = cell_ndims + 1
        if (periodic) then
          cell_dmin(ii) = 0d0
          cell_dmax(ii) = 1d0
        endif
      endif
    enddo

    ! This is only to avoid division by zero later on
    cell_dmin(:) = cell_dmin(:) - TOL_Small
    cell_dmax(:) = cell_dmax(:) + TOL_Small

! FHJ: number of cells is proportional to the length of each active dimension
!      this should work for any number of dimensions. Proof:
!        d1*d2*... = V
!        n1 = c*d1*(N)^(1/D)
!        n1*n2*... = c^(D)*V*N = N => c = (1/V)^(1/D)
!      so, const = (N/V)^(1/D)
    cell_N(:) = 1
    if (cell_ndims>0) then
      dim_vol=1.0d0
      do ii=1,3
        if (cell_active_dim(ii)) dim_vol = dim_vol * (cell_dmax(ii)-cell_dmin(ii))
      enddo
      prop_const = (ncoarse/dim_vol)**(1.0d0/cell_ndims)
      do ii=1,3
        if (cell_active_dim(ii)) then
          cell_N(ii) = dnint(prop_const * (cell_dmax(ii)-cell_dmin(ii)))
          cell_N(ii) = min(max(cell_N(ii), 1), ncoarse)
        endif
      enddo
    endif
    
    cell_factor(:) = cell_N(:) / (cell_dmax(:) - cell_dmin(:))
    cell_shift(:) = (cell_dmax(:) - cell_dmin(:)) / (2.0d0*cell_N(:))

    if (peinf%inode==0.and.should_write) then
      write(6,*)
      write(6,'(a,i6,a)') ' Automatically creating a cell structure for ',ncoarse,' points'
      write(6,'(a,i1,a)') '  Found ',cell_ndims,' dimension(s)'
      write(6,'(a,i5,i5,i5)') '  Number of cells:', cell_N(1),cell_N(2),cell_N(3)
      write(6,'(a,i6)') '  Total number of cells:', cell_N(1)*cell_N(2)*cell_N(3)
    endif
   
    SAFE_ALLOCATE(cell_head, (cell_N(1),cell_N(2),cell_N(3)))
    SAFE_ALLOCATE(cell_list, (ncoarse))

! FHJ: Initialize and populate cells
    cell_head(:,:,:) = 0

    if (peinf%inode==0.and.should_write) then
      write(6,*)
      do ii=1,3
        write(6,801) ii,cell_dmin(ii),cell_dmax(ii),cell_shift(ii)*2.0d0
801     format(' Cells [',i1,'], dmin= ',f8.5,' dmax= ',f8.5,' length= ',f12.5)
      enddo
    endif

    do idx_co=1,ncoarse
      call get_cell_idx(coarsepts(:,idx_co), cell_idx, periodic)
      cell_list(idx_co) = cell_head(cell_idx(1),cell_idx(2),cell_idx(3))
      cell_head(cell_idx(1),cell_idx(2),cell_idx(3)) = idx_co
    enddo

    if (peinf%inode==0.and.should_write) then
      write(6,*)
      write(6,*) 'Cell Population Analysis'
      write(6,*)
      write(6,900) ' x ',' y ',' z ',' members '
      write(6,900) '---','---','---','---------'
  900 format(2x, 3(a3,1x), a9)
      do j1=1,cell_N(1)
        do j2=1,cell_N(2)
          do j3=1,cell_N(3)
            write(6,'(2x,3(i3,1x))',advance='no') j1,j2,j3      
            idx_co=cell_head(j1,j2,j3)
            do while (idx_co.gt.0)
              write(6,'(1x,i5)',advance='no') idx_co
              idx_co=cell_list(idx_co)
            enddo
            write(6,*)
          enddo
        enddo
      enddo
      write(6,'(/,a)') 'Finished alloc_intpts'
    endif

    POP_SUB(alloc_intpts)

  end subroutine alloc_intpts

  subroutine dealloc_intpts()
    PUSH_SUB(dealloc_intpts)

    SAFE_DEALLOCATE(cell_head)
    SAFE_DEALLOCATE(cell_list)

    POP_SUB(dealloc_intpts)
  end subroutine dealloc_intpts

  !> Returns the cell index for point `pt`.
  subroutine get_cell_idx(pt, cell_idx, periodic)
    real(DP), intent(in) :: pt(3)
    integer, intent(out) :: cell_idx(3)
    logical, intent(in) :: periodic

    ! no push/pop, called too frequently
    cell_idx(:) = idint((pt(:)-cell_dmin(:)+cell_shift(:))*cell_factor(:))
    if (periodic) then
      cell_idx(:) = modulo(cell_idx(:), cell_N(:)) + 1
    else
      cell_idx(:) = min(max(cell_idx(:)+1, 1), cell_N(:))
    endif

  end subroutine get_cell_idx

end module intpts_m
