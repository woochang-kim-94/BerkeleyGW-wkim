
#include "f_defs.h"

module intp_utils_m

  use global_m

  implicit none

  private

  public :: calc_max_dist, bleed_bz, alloc_cells, get_cell_idx, fix_cell_idx, &
    bled_grid, pp, cell_struc

  !type for the `bled grid`, i.e., grid with more kpts than the 1st BZ
  !note: `bleeding` is a term used in the context of printing documents:
  ! see http://en.wikipedia.org/wiki/Bleed_(printing)
  type bled_grid
    integer :: nq                !number of pts
    integer, pointer :: idx(:)   !index of the reduced-BZ point
    integer, pointer :: kg0(:,:) !umklapp
    real(DP), pointer :: q(:,:)  !all points (includes BZ)
  end type bled_grid

  !defines a pointer
  type pp
    real(DP), pointer :: p(:,:)
  end type pp

  !defines a cell structure
  type cell_struc
    integer, allocatable :: head(:,:,:), list(:)
    real(DP), dimension(3) :: dmin, dmax, factor, shift
    integer, dimension(3) :: N
    ! FHJ: cell_factor = 1/length of the individual cell
    !      cell_shift = length of the individual cell / 2
  end type cell_struc

  real(DP), parameter :: fact=1.15

contains

  subroutine get_cell_idx(cells, qq, cell_idx)
    type(cell_struc), intent(in) :: cells
    real(DP), intent(in) :: qq(3)
    integer, intent(out) :: cell_idx(3)

    integer :: ii

    ! no push/pop since called too frequently

    do ii=1,3
      cell_idx(ii) = idint((qq(ii)-cells%dmin(ii)-cells%shift(ii))*cells%factor(ii)+TOL_Small)+1
      if (cell_idx(ii).gt.cells%N(ii)) cell_idx(ii) = cell_idx(ii)-cells%N(ii)
      if (cell_idx(ii).gt.cells%N(ii) .or. cell_idx(ii).lt.1) then
        write(0,'(1x,a,3(f10.7,1x))') ' qq=',qq
        write(0,'(1x,a,3(i5,1x))') ' index cell= ',cell_idx
        call die('Invalid index for cell', only_root_writes = .true.)
      endif
    enddo

    ! no push/pop since called too frequently

  end subroutine get_cell_idx

  !put cell_idx inside [1, cells%N]
  subroutine fix_cell_idx(cells, cell_idx)
    type(cell_struc), intent(in) :: cells
    integer, intent(inout) :: cell_idx(3)
    integer :: ii

    ! no push/pop since called too frequently

    do ii=1,3
      !MOD vs. MODULO: http://www.nsc.liu.se/~boein/f77to90/a5.html
      cell_idx(ii) = modulo(cell_idx(ii), cells%N(ii))
      if(cell_idx(ii)==0) cell_idx(ii) = cell_idx(ii) + cells%N(ii)
    enddo

    ! no push/pop since called too frequently

  end subroutine fix_cell_idx

  subroutine alloc_cells(bg, cells, dx)
    type(bled_grid), intent(in) :: bg
    type(cell_struc), intent(out) :: cells
    real(DP), intent(in) :: dx

    integer :: ii, idx
    integer :: j1,j2,j3
    integer :: cell_idx(3)

    PUSH_SUB(alloc_cells)

    cells%shift = dx * 0.5d0
    cells%factor = 1.0d0/dx
    cells%dmin = minval(bg%q)-cells%shift
    cells%dmax = maxval(bg%q)+cells%shift
    do ii=1,3
      cells%N(ii) = idint((cells%dmax(ii)-cells%dmin(ii))/dx + TOL_Small)
    enddo

    write(6,'(1x,A,3(I5,1x))') 'Number of Cells:',cells%N
    SAFE_ALLOCATE(cells%head, (cells%N(1), cells%N(2), cells%N(3)) )
    SAFE_ALLOCATE(cells%list, (bg%nq) )
    cells%head(:,:,:) = 0

    write(6,*)
    do ii=1,3
      write(6,801) ii,cells%dmin(ii),cells%dmax(ii),cells%shift(ii)*2.0d0
801   format(' Cells [',i1,'], dmin= ',f8.5,' dmax= ',f8.5,' length= ',f12.5)
    enddo

    do idx=1, bg%nq
      call get_cell_idx(cells, bg%q(:,idx), cell_idx)
      cells%list(idx) = cells%head(cell_idx(1),cell_idx(2),cell_idx(3))
      cells%head(cell_idx(1),cell_idx(2),cell_idx(3)) = idx
    enddo

    if (all(cells%list==0)) then
      write(6,*) 'All cells were properly populated'
    else
      write(6,*) 'Some cells have occupancy>1. Details:'
      write(6,*)
      write(6,*) 'Cell Population Analysis'
      write(6,*)
      write(6,900) ' x ',' y ',' z ',' members '
      write(6,900) '---','---','---','---------'
900   format(2x, 3(a3,1x), a9)
      do j1=1,cells%N(1)
        do j2=1,cells%N(2)
          do j3=1,cells%N(3)
            idx=cells%head(j1,j2,j3)
            !if(idx>0) then
            write(6,'(2x,3(i3,1x))',advance='no') j1,j2,j3      
            do while (idx.gt.0)
              write(6,'(1x,i5)',advance='no') idx
              idx=cells%list(idx)
            enddo
            write(6,*)
            !endif
          enddo
        enddo
      enddo
    endif
    write(6,*)

    POP_SUB(alloc_cells)

  end subroutine alloc_cells

  !add extra kpts outside the BZ
  subroutine bleed_bz(crys, gr, bg, nq, q_in)
    type(crystal), intent(in) :: crys
    type(grid), intent(in) :: gr
    type(bled_grid), intent(out) :: bg
    integer, intent(in) :: nq
    real(DP), intent(in) :: q_in(3,nq)

    !using these vars while I don`t know the size of the array
    integer, allocatable ::  tmp_idx(:)
    real(DP), allocatable :: tmp_q(:,:)
    integer, allocatable ::  tmp_kg0(:,:)   

    integer n_bled, n_gvects

    integer :: idx_f
    real(DP) :: max_dist, l_q, qq(3)
    integer :: gx,gy,gz

    PUSH_SUB(bleed_bz)

    !calculate maximum distance |q|
    call calc_max_dist(crys, nq, q_in(:,:), max_dist)

    n_gvects = 3*3 -1 !because we are in 2D
    SAFE_ALLOCATE(tmp_idx, (n_gvects*gr%nf))
    SAFE_ALLOCATE(tmp_q, (3, n_gvects*gr%nf))
    SAFE_ALLOCATE(tmp_kg0, (3, n_gvects*gr%nf))

    n_bled = 0
    do gx = -1,1
      do gy = -1,1
        gz = 0
        if ((gx==0).and.(gy==0).and.(gz==0)) cycle

        do idx_f = 1,gr%nf
          qq = gr%f(:,idx_f) + ((/gx,gy,gz/)*1.0d0)
          l_q = DOT_PRODUCT(qq,MATMUL(crys%bdot,qq))
          if (l_q > max_dist*fact) cycle
          !ok, we accepted this distance...

          n_bled = n_bled + 1
          tmp_idx(   n_bled) = gr%indr(idx_f)
          tmp_q  (:, n_bled) = qq
          tmp_kg0(:, n_bled) = (/gx,gy,gz/)
        enddo
      enddo
    enddo

    !now we can properly allocate the arrays
    bg%nq = gr%nf + n_bled
    SAFE_ALLOCATE(bg%idx, (bg%nq))
    SAFE_ALLOCATE(bg%q,   (3, bg%nq))
    SAFE_ALLOCATE(bg%kg0, (3, bg%nq))

    !copy the older points
    bg%idx(   1:gr%nf) = gr%indr(1:gr%nf)
    bg%q  (:, 1:gr%nf) = gr%f(:, 1:gr%nf)
    bg%kg0(:, 1:gr%nf) = 0

    !and now copy the newer ones
    bg%idx(   gr%nf+1:bg%nq) = tmp_idx(1:n_bled)
    bg%q  (:, gr%nf+1:bg%nq) = tmp_q  (:, 1:n_bled)
    bg%kg0(:, gr%nf+1:bg%nq) = tmp_kg0(:, 1:n_bled)
    SAFE_DEALLOCATE(tmp_idx)
    SAFE_DEALLOCATE(tmp_q)
    SAFE_DEALLOCATE(tmp_kg0)

    POP_SUB(bleed_bz)

  end subroutine bleed_bz

  subroutine calc_max_dist(crys, nq, q_in, max_dist)
    type(crystal), intent(in) :: crys
    integer, intent(in) :: nq
    real(DP), intent(in) :: q_in(3,nq)
    real(DP), intent(out) :: max_dist

    integer :: inq
    real(DP) :: length

    ! no push/pop since called too frequently

    max_dist=0.0d0
    do inq=1,nq
      length = DOT_PRODUCT(q_in(:,inq),MATMUL(crys%bdot,q_in(:,inq)))
      if (length>max_dist) max_dist=length
    enddo

    ! no push/pop since called too frequently

  end subroutine calc_max_dist

end module intp_utils_m
