!=============================================================================
!
! Utilities:
!
! (1) epsmat_intp()  Originally by FHJ       Last Modified 11/04/2011 (FHJ)
!
!<    Doubles the number of points in an epsmat file by interpolating it.
!
!     This utility *doubles* the number of q-points in an epsmat file by doing
!     several interpolation/extrapolation techniques. Note that you *can`t*
!     interpolate onto *arbitrary grids*. This was a design choice to
!     increase the accuracy of the interpolation.
!
!     You must provide a file epsmat_intp.inp (see example one).
!     
!     Current limitations:
!      - Only REAL version
!      - Only works for 1D and 2D systems
!      - Only GPP
!
!=============================================================================

#include "f_defs.h"

!Z symmetry speeds up the algorithm to find vectors. Use for graphene only
! (the whole program only works for graphene anyways...)
#define Z_SYMMETRY

!The following macros should be used inside interpolate_1D and interpolate_2D,
! and the parameter is always the *index* of the eps_pt array.
!The macros assume that the spacing b/w two coarse points is twice that of the
! interpolated point and the closest coarse point. The diagrams use the
! following convention:
! number=given coarse point, x=intp point, o=unused coarse point

!Linear extrapolation [x,1,2]
! o x 1   2   o 
#define EXTP_2(i1,i2) \
  (1.5d0*eps_pt(i1)%eps - 0.5d0*eps_pt(i2)%eps)

!Quadratic extrapolation [x,1,2,3]
! o x 1   2   3
#define EXTP_3(i1,i2,i3) \
  0.125d0*(15d0*eps_pt(i1)%eps - 10d0*eps_pt(i2)%eps + 3d0*eps_pt(i3)%eps)

!Cubic symmetric interpolation [1,2,x,3,4]
! 1   2 x 3   4
#define INTP_4(i1,i2,i3,i4) \
  0.0625d0*(9d0*(eps_pt(i2)%eps+eps_pt(i3)%eps) - (eps_pt(i1)%eps+eps_pt(i4)%eps))

!Quadratic asymmetric interpolation [1,x,2,3]
! o   1 x 2   3
#define INTP_3(i1,i2,i3) \
  0.125d0*(3d0*eps_pt(i1)%eps + 6d0*eps_pt(i2)%eps - eps_pt(i3)%eps)


!This macro should be used inside interpolate_1D and interpolate_2D
#define AVG_2(i1,i2) ((eps_pt(i1)%eps + eps_pt(i2)%eps)*0.5d0)

program epsmat_intp

  use global_m
  use wfn_rho_vxc_io_m
  use sort_m
  use intp_utils_m
  use fullbz_m

  implicit none

  type eps_point
    integer :: idx_bled, idx_fold
    integer :: ieps, jeps
    real(DP) :: eps
    logical :: avail
  end type eps_point

  character(len=6) :: ajname
  character(len=11) :: adate
  real(DP) :: div,ecuts
  integer :: i_1,ii,iq,j_1,ngq,nq, qgrid(3), &
             nq_intp, icol, col_bar, last_col_bar

  integer :: nq_intp_save, iq_
  integer, pointer :: nc_map(:) !maps the non-coarse points into indices

  character(len=80) :: wfn, file_in, file_out
  real(DP), allocatable :: q_eps(:,:), q_intp(:,:)

  !eps is an array of pointers
  type(pp), pointer :: eps(:)

  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvec_wfn
  type(kpoints) :: kp
  type(crystal) :: crys
  type(symmetry) :: syms
  type(grid) :: gr
  type(bled_grid) :: bg
  type(cell_struc) :: cells
  real(DP) :: dx
  real(DP), allocatable :: ekin(:,:), ekin_intp(:) ! to get # of G vectors for interpolated pts
  integer, allocatable :: isort_intp(:), isort_i_intp(:), gvec(:,:)
  integer :: ng
  integer, pointer :: nmtx_intp(:)
  integer :: cell_idx(3)
  real(DP) :: qq(3), tmp_qq(3)
  logical :: match_dim(3)
  integer :: line_dir, ig,jg, ieps,jeps, intp_kind
  real(DP) :: line_dx(3) ! direction of dx (times dx)
  real(DP), pointer :: eps_intp(:,:)
  integer :: ax(3)
  integer :: ierr
  logical :: skip_check

  integer :: gvec_i(3), gvec_j(3) !these are the actual g vectors associated to ig & jg

  integer, allocatable :: isort(:,:), isort_i(:,:), nmtx(:)
  integer, parameter :: F_WFN=10, F_INP=11, F_EPSIN=12, F_EPSOUT=13
  integer, parameter :: MAX_NG=100000
  integer, parameter :: MAX_NQ=1000

  PUSH_SUB(epsmat_intp)

  !read input file
  call open_file(F_INP, file='epsmat_intp.inp', form='formatted', status='old')
  read(F_INP,*) wfn
  read(F_INP,*) file_in
  read(F_INP,*) file_out
  read(F_INP,*) nq_intp
  SAFE_ALLOCATE(q_intp, (3, nq_intp))
  do iq = 1,nq_intp
    read(F_INP,*) q_intp(1,iq),q_intp(2,iq),q_intp(3,iq),div
    q_intp(:,iq) = q_intp(:,iq)/div
  enddo
  call close_file(F_INP)
 
  SAFE_ALLOCATE(q_eps, (3, MAX_NQ))
  SAFE_ALLOCATE(gvec, (3, MAX_NG))

  !read WFN b/c we need the symmetries
  write(6,*) 'Reading WFN:  ', TRUNC(wfn)
  sheader = 'WFN'
  iflavor = 0
  call open_file(F_WFN, file=wfn, form='unformatted', status='old')
  call read_binary_header_type(F_WFN, sheader, iflavor, kp, gvec_wfn, syms, crys)
  call close_file(F_WFN)

  write(6,*) 'Reading epsmat file from: ', TRUNC(file_in)

  !read the epsmat file
  call open_file(F_EPSIN, file=file_in, form='unformatted', status='old')
  read(F_EPSIN) ajname,adate
  read(F_EPSIN) !freq_dep,nFreq
  read(F_EPSIN) (qgrid(j_1),j_1=1,3)
  read(F_EPSIN) !assuming not freq dep.
  read(F_EPSIN)
  read(F_EPSIN)
  read(F_EPSIN) ecuts
  read(F_EPSIN) nq,((q_eps(j_1,iq),j_1=1,3),iq=1,nq)
  read(F_EPSIN) ng,((gvec(j_1,ig),j_1=1,3),ig=1,ng)

  SAFE_ALLOCATE(isort  , (ng,nq))
  SAFE_ALLOCATE(isort_i, (ng,nq))
  SAFE_ALLOCATE(ekin   , (ng,nq))
  SAFE_ALLOCATE(nmtx   , (nq))
  !some compilers (cray) don`t know how to SIZEOF an array of types
#ifndef NOSIZEOF_TYPE
  SAFE_ALLOCATE(eps    , (nq))
#else
  allocate(eps(nq))
#endif
 
  write(6,'(1x,a,i5,a)') 'Reading ',nq,' q_eps-points'
  
  write(6,'(1x,a)',advance='yes') 'Reading: |                                                  |' 
  write(6,'(1x,a)',advance='no' ) '         +'
  FLUSH(6)
  last_col_bar = 0
  do iq=1, nq 
  !do iq=1, 1
    read(F_EPSIN) ngq, nmtx(iq), (isort(icol, iq),isort_i(icol, iq), icol=1,ng)
    read(F_EPSIN) (ekin(icol, iq), icol=1,ng)
    read(F_EPSIN) !(q2(j,iq),j=1,3)
    SAFE_ALLOCATE(eps(iq)%p, (nmtx(iq),nmtx(iq)))
    
    do j_1=1, nmtx(iq)
      read(F_EPSIN) (eps(iq)%p(icol,j_1), icol=1,nmtx(iq))
    enddo

    col_bar = int((iq*1.0d0/nq)*50.0d0)
    do icol=1,(col_bar-last_col_bar)
      write(6,'(a)',advance='no') '-'
      FLUSH(6)
    enddo
    last_col_bar = col_bar
  enddo
  write(6,'(a)') '+'
  write(6,*)

  call close_file(F_EPSIN)

  !unfold the q-grid from the epsmat file 
  write (6,*) 'Unfolding BZ'
  SAFE_ALLOCATE(gr%r, (3,nq))
  gr%nr = nq
  gr%r(:,:) = q_eps(:,:)
  skip_check=.false.
  call fullbz(crys, syms, gr, syms%ntran, skip_check, wigner_seitz=.true., paranoid=.false.) 

  !add extra points at the border
  write(6,*) 'Bleeding BZ'
  call bleed_bz(crys, gr, bg, nq, q_eps)
  write(6,*) 'BZ has',bg%nq,'points'

  call dealloc_grid(gr)

  !create cell structure
  dx = sum(q_eps(:,2)-q_eps(:,1))
  write(6,'(1x,A,F8.5)') 'dx=',dx
  call alloc_cells(bg, cells, dx)

  SAFE_ALLOCATE(nmtx_intp,    (nq_intp))
  SAFE_ALLOCATE(nc_map,       (nq_intp))
  SAFE_ALLOCATE(ekin_intp,    (ng))
  SAFE_ALLOCATE(isort_intp,   (ng))
  SAFE_ALLOCATE(isort_i_intp, (ng))

  write(6,*) 'Calculating rank of eps matrix per q_eps-pt'
  !nc_map: map a non-coarse pt into a interpolated point
  !nq_intp_save:  number of interpolated points that will be saved (they are
  !               either fine points *or* coarse points that exists in the grid)
  !nq_intp: total number of interpolated points (including those that will be ignored)
  nq_intp_save = 0
  nc_map(:) = 0
  nmtx_intp(:) = 0
  !this loop is only needed to construct the header of the epsmat file
  !TODO: compute these things only once!
  do iq=1, nq_intp
    !Find if this is a coarse point
    do ii=1,3
      tmp_qq(ii) = (q_intp(ii,iq)-cells%dmin(ii)-cells%shift(ii))*cells%factor(ii)
      match_dim(ii) = dabs(tmp_qq(ii)-dnint(tmp_qq(ii)))<TOL_Small
    enddo

    if (all(match_dim)) then 
      !see if the coarse point exists
      call get_cell_idx(cells, q_intp(:,iq), cell_idx(:))
      iq_ = cells%head(cell_idx(1), cell_idx(2), cell_idx(3))
      if (iq_==0) cycle      
      !ok, reuse the values for kinetic energy and nmtx
      iq_ = bg%idx(iq_)
      nmtx_intp(iq) = nmtx(iq_)
    else   
      !calculate kinetic energy
      do i_1=1,ng
        qq = q_intp(:,iq) + gvec(:,i_1)
        ekin_intp(i_1) = DOT_PRODUCT(qq,MATMUL(crys%bdot,qq))
      enddo
      call sortrx(ng, ekin_intp, isort_intp, gvec)
      do i_1=1,ng
        isort_i_intp(isort_intp(i_1)) = i_1
      enddo
      nmtx_intp(iq) = count(ekin_intp<ecuts)
    endif
    !we accepte this point, save it
    nq_intp_save = nq_intp_save + 1
    nc_map(nq_intp_save) = iq
  enddo

  write(6,*) 'Writing epsmat file header to: ', TRUNC(file_out)
  call open_file(F_EPSOUT, file=file_out, form='unformatted', status='unknown')
  write(F_EPSOUT) ajname,adate
  write(F_EPSOUT) 0,1
  write(F_EPSOUT) (qgrid(j_1),j_1=1,3)
  write(F_EPSOUT) !assuming not freq dep.
  write(F_EPSOUT)
  write(F_EPSOUT)
  write(F_EPSOUT) ecuts
  write(F_EPSOUT) nq_intp_save,((q_intp(j_1,nc_map(iq)),j_1=1,3),iq=1,nq_intp_save)
  write(F_EPSOUT) ng,((gvec(j_1,ig),j_1=1,3),ig=1,ng)

  SAFE_ALLOCATE(eps_intp, (maxval(nmtx_intp),maxval(nmtx_intp)))

  !loop through the q-points to be interpolated and do the fun part
  do iq=1, nq_intp
    write(6,*)
    write(6,'(1x,A,I4,": ",3(F10.7,1x),A)') '*** Dealing with point ', iq, q_intp(:,iq),' ***'
    write(6,*)

    !Find if the interpolated point is lying on a line or in the center of a square
    do ii=1,3
      tmp_qq(ii) = (q_intp(ii,iq)-cells%dmin(ii)-cells%shift(ii))*cells%factor(ii)
      match_dim(ii) = dabs(tmp_qq(ii)-dnint(tmp_qq(ii)))<TOL_Small
    enddo
    
    select case (count(match_dim))
      case (3)
        intp_kind = 0
      case (2)
        line_dx = 0
        do ii=1,3
          if(.not.match_dim(ii)) then
            line_dir = ii
            line_dx(ii) = dx
            exit
          endif
        enddo
        write(6,'(A,I1)') ' Doing 1D interpolation, dir=',line_dir
        intp_kind = 1
      case (1)
        write(6,*) 'Doing 2D interpolation'        
        intp_kind = 2
    end select

    !store the axis where the interpolated plane lives
    j_1 = 1
    do ii=1,3
      if(.not.match_dim(ii)) then
        ax(j_1) = ii
        j_1 = j_1 + 1
      endif
    enddo


    if (intp_kind==0) then 
      !see if the coarse point exists
      call get_cell_idx(cells, q_intp(:,iq), cell_idx(:))
      iq_ = cells%head(cell_idx(1), cell_idx(2), cell_idx(3))
      if (iq_==0) then
        write(6,*) 'Ignoring q0'
        cycle
      endif
      write(6,*) 'Copying point'
      !ok, reuse the values for kinetic energy and nmtx
      iq_ = bg%idx(iq_)
      isort_intp(:) = isort(:,iq_)
      isort_i_intp(:) = isort_i(:,iq_)
      ekin_intp(:) = ekin(:,iq_)
    else   
      !calculate kinetic energy
      do i_1=1,ng
        qq = q_intp(:,iq) + gvec(:,i_1)
        ekin_intp(i_1) = DOT_PRODUCT(qq,MATMUL(crys%bdot,qq))
      enddo
      call sortrx(ng, ekin_intp, isort_intp, gvec)
      do i_1=1,ng
        isort_i_intp(isort_intp(i_1)) = i_1
      enddo      
    endif

    !Write headers
    write(F_EPSOUT) ng, nmtx_intp(iq), (isort_intp(icol),isort_i_intp(icol), icol=1,ng)
    write(F_EPSOUT) (ekin_intp(icol), icol=1,ng)
    write(F_EPSOUT) (q_intp(j_1,iq),j_1=1,3)

    !ekin, gvec: sorted by index of G
    !eps: sorted by column of eps (increasing value of ekin)
    !isort -> get (index of) G  given (column of) eps
    !isort_i -> get (column of) eps given (index of) G
    !nmtx: size of eps matrix = number values of ekin less than ecuts

    write(6,'(1x,a,i5,a)') 'Writing ',nmtx_intp(iq),' matrix elements into epsmat file'
    if (intp_kind==0) then 
      !no need to interpolate, just copy from coarse grid
      do j_1=1, nmtx_intp(iq)
        write(F_EPSOUT) (eps(iq_)%p(icol,j_1), icol=1,nmtx_intp(iq))
      enddo
    else
      !calculate
      do ieps = 1,nmtx_intp(iq)
        do jeps = 1,nmtx_intp(iq)
          ig = isort_intp(ieps)
          jg = isort_intp(jeps)
          gvec_i = gvec(:,ig)
          gvec_j = gvec(:,jg)

          if (intp_kind==1) then
            call interpolate_1D(ierr)
          else
            call interpolate_2D(ierr,ax)
          endif

          !if the interpolation failed, try extrapolating by the G vectors
          if (ierr/=0) then
            call extrapolate_by_gvec(ierr)
          endif
        enddo
      enddo

      do j_1=1, nmtx_intp(iq)
        write(F_EPSOUT) (eps_intp(icol,j_1), icol=1,nmtx_intp(iq))
      enddo
    endif

  enddo

  write(6,*) 'Done'
  write(6,*)

  call close_file(F_EPSOUT)
  PUSH_SUB(epsmat_intp)

contains

  !return the eps_pt associate to the point at cell # cell_idx_
  !this structure contains the eps matrix element for the corresponding
  !gvecs ig,jg (for the columns ieps,jeps). Note that we point may not be
  !available (check for eps_pt%avail)
  subroutine get_eps_point(cell_idx_, eps_pt)
    integer, intent(in) :: cell_idx_(3)
    type(eps_point), intent(out) :: eps_pt

    integer :: ibled, ifold
    integer ig_, jg_
    integer :: gvec_i_(3), gvec_j_(3)

    ! no push/pop since called too frequently

    eps_pt%idx_bled = cells%head(cell_idx_(1), cell_idx_(2), cell_idx_(3)) ! idx of bled point
    if (eps_pt%idx_bled==0) then
      !point doesn`t exists, probably its the q=0
      eps_pt%avail=.false.
      eps_pt%idx_bled=0
      eps_pt%idx_fold=0
      return
    endif
    eps_pt%idx_fold = bg%idx(eps_pt%idx_bled) ! and the corresponding index in the folded BZ, used for eps
    ibled = eps_pt%idx_bled ! just an alias
    ifold = eps_pt%idx_fold ! just an alias
    if (all(bg%kg0(:, ibled)==0)) then
      !write(6,*) 'Kein umklapp!'
      ig_ = ig
      jg_ = jg
    else
      !write(6,'(1x,a,3(I3,x),"!")') 'Mit Umklapp =',bg%kg0(:,idx_(ipt))
      gvec_i_(:) = gvec_i + bg%kg0(:, ibled) !add umklapp to the gvecs
      gvec_j_(:) = gvec_j + bg%kg0(:, ibled)
      call find_gvec(ig_, gvec_i_(:)) !and find these new gvecs
      call find_gvec(jg_, gvec_j_(:))
    endif

    if (ig_==0.or.jg_==0) then
      !the gvecs are out of our space!
      !there is probably something really wrong going on...
      eps_pt%ieps = 0
      eps_pt%jeps = 0
      eps_pt%avail = .false.
    else
      !gvecs are valid. test if the corresponding eps columns were calculated
      !now, test if ieps_ and jeps_ are in the bounds for these points
      eps_pt%ieps = isort_i(ig_, ifold)
      eps_pt%jeps = isort_i(jg_, ifold)
      eps_pt%avail = (eps_pt%ieps<=nmtx(ifold)).and.(eps_pt%jeps<=nmtx(ifold))
      if (eps_pt%avail) then
        !great success! store the eps matrix element
        eps_pt%eps = eps(ifold)%p(eps_pt%ieps, eps_pt%jeps)
      endif
    endif

    ! no push/pop since called too frequently

  end subroutine get_eps_point

  ! use this if nothing else works...
  subroutine extrapolate_by_gvec(ierr)
    integer, intent(out) :: ierr

    integer :: tmp_gvec(3)
    integer :: sig, i_1, ii
    integer :: g_dir, gmax
    logical :: use_i
    integer :: tmp_epscol, tmp_gcol
    integer :: eps_var
    real(DP) :: tmp_eps(3)

    PUSH_SUB(epsmat_intp.extrapolate_by_gvec)

    ierr = 0
    !find the largest component of the g vector
    gmax = -1
    do ii=1,3
      if(abs(gvec_j(ii))>=gmax) then
        tmp_gvec(:) = gvec_j(:)
        gmax = abs(gvec_j(ii))
        g_dir = ii
        use_i = .false.
      endif
      if(abs(gvec_i(ii))>=gmax) then
        tmp_gvec(:) = gvec_i(:)
        gmax = abs(gvec_i(ii))
        g_dir = ii
        use_i = .true.
      endif
    enddo

    sig = 1
    if (tmp_gvec(g_dir)<1) sig=-1
    if (use_i) then
#ifdef VERBOSE
      write(6,*) 'Extrapolating i_1, dim=',g_dir
#endif
      eps_var = ieps
    else
#ifdef VERBOSE
      write(6,*) 'Extrapolating j_1, dim=',g_dir
#endif
      eps_var = jeps
    endif

    call find_gvec(tmp_gcol, tmp_gvec(:))
#ifdef VERBOSE
    write(6,'(1x,a,3(i3,1x))') 'G=', gvec(:, tmp_gcol)
    write(6,*)
#endif

    do i_1 = 1,3
      !we always start with the largest possible G vector, and go down
      tmp_gvec(g_dir) = tmp_gvec(g_dir) - sig
      call find_gvec(tmp_gcol, tmp_gvec(:))
      tmp_epscol = isort_i_intp(tmp_gcol)
#ifdef VERBOSE
      write(6,'(1x,a,(i5,1x),a,3(i3,1x))') 'epscol=', tmp_epscol,'  G=', gvec(:, tmp_gcol)
#endif
      if (tmp_epscol<eps_var) then
        !if the last G vector was already calculated, save the value
        if (use_i) then
          tmp_eps(i_1) = eps_intp(tmp_epscol, jeps)
        else
          tmp_eps(i_1) = eps_intp(ieps, tmp_epscol)
        endif
#ifdef VERBOSE
        write(6,*) ' eps(G,Gp)=',tmp_eps(i_1)
#endif
      else
        !if not, we don`t have enough points to do the full quadratic
        !extrapolation, so let`s do the best we can with what we`ve got
        select case (i_1)
        case (1)
#ifdef VERBOSE
          call die('no G vector available for extrapolation!')
#endif
        case (2)
          !use NN
#ifdef VERBOSE
          write(6,*) 'using NN'
#endif
          eps_intp(ieps,jeps) = tmp_eps(1)
          POP_SUB(epsmat_intp.extrapolate_by_gvec)
          return
        case (3)
          !linear extrapolation
#ifdef VERBOSE
          write(6,*) 'using linear extrapolation'
#endif
          eps_intp(ieps,jeps) = 2.0*tmp_eps(1) - tmp_eps(2)
          POP_SUB(epsmat_intp.extrapolate_by_gvec)
          return
        end select
      endif
#ifdef VERBOSE
      write(6,*)
#endif
    enddo
    !got all points, do quadratic extrapolation
    eps_intp(ieps,jeps) = 3d0*(tmp_eps(1) - tmp_eps(2)) + tmp_eps(3)
#ifdef VERBOSE
    write(6,*) 'extrapolation:',eps_intp(ieps,jeps)
    write(6,*)
#endif

    POP_SUB(epsmat_intp.extrapolate_by_gvec)

  end subroutine extrapolate_by_gvec
    

  ! Interpolate in a line. Used in situations like:
  ! 1   2 x 3   4
  !              
  ! o   o   o   o
  !------------------------------------------------------------------------+
  subroutine interpolate_1D(ierr)!                                         |
  !------------------------------------------------------------------------+
    integer, intent(out) :: ierr

    real(DP) :: qq1(3)
    integer :: cell_idx_(3,5) !pt #5 is extra
    logical :: pt_avail(4)
    type(eps_point) :: eps_pt(5) !pt #5 is extra
    integer :: ipt
    integer :: cel_num(3), sig !used to order extra point

    PUSH_SUB(interpolate_1D)

    pt_avail(:) = .false.

    !qq1 is the coarse point 1 in the diagram above
    qq1(:) = q_intp(:,iq)
    qq1(line_dir) = qq1(line_dir) - 1.5d0*dx
    call get_cell_idx(cells, qq1(:), cell_idx_(:,1))

    do ipt=1, 4
      !for each step, move 1 cell right/up
      if (ipt>1) then
        cell_idx_(:,ipt) = cell_idx_(:, ipt-1)
        cell_idx_(line_dir, ipt) = cell_idx_(line_dir, ipt) + 1
      endif
      call get_eps_point(cell_idx_(:,ipt), eps_pt(ipt))
      pt_avail(ipt) = eps_pt(ipt)%avail
    enddo

    select case (count(pt_avail))

      !---------------------------+
    case (4) !      4 VALID POINTS       |
      !---------------------------+

      !cubic symmetric interpolation

      eps_intp(ieps,jeps) = INTP_4(1,2,3,4)

      !---------------------------+
    case (3) !      3 VALID POINTS       |
      !---------------------------+

      !do we have the nearest points?
      if (pt_avail(2).and.pt_avail(3)) then
        !yes, do asymmetric quadratic interpolation

        if (pt_avail(1)) then
          ! 1   2 x 3   o
          eps_intp(ieps,jeps) = INTP_3(3,2,1)
        else
          ! o   2 x 3   4
          eps_intp(ieps,jeps) = INTP_3(2,3,4)
        endif
        POP_SUB(epsmat_intp.interpolate_1D)
        return
      endif

      !lets still try quadratic extrapolation. need a new point...
      if (pt_avail(2)) then
        ! 5   1   2 x o
        cel_num(:) = (/2, 1, 5/)
        sig = -1
      else
        ! o x 3   4   5
        cel_num(:) = (/3, 4, 5/)
        sig = +1
      endif

      cell_idx_(:,5) = cell_idx_(:,cel_num(2))
      cell_idx_(line_dir,5) = cell_idx_(line_dir,5) + sig
      call get_eps_point(cell_idx_(:,5), eps_pt(5))
      if (eps_pt(5)%avail) then
        !quadratic extrapolation is possible
        eps_intp(ieps,jeps) = EXTP_3(cel_num(1),cel_num(2),cel_num(3))
        POP_SUB(epsmat_intp.interpolate_1D)
        return
      endif

      !no extra point found, do linear extrapolation
      eps_intp(ieps,jeps) = EXTP_2(cel_num(1),cel_num(2))

      !---------------------------+
    case (2) !      2 VALID POINTS       |
      !---------------------------+

      !if we have no nearest neighbor, fail
      if (.not.(pt_avail(2).or.pt_avail(3))) then
        ierr=1
        POP_SUB(epsmat_intp.interpolate_1D)
        return
      endif

      !do we have both nearest pts?
      if (pt_avail(2).and.pt_avail(3)) then
        !just take the average
        eps_intp(ieps,jeps) = AVG_2(2,3)
#ifdef VERBOSE
        write(6,'(1x,a,3(F10.6,1x))') '[2] eps_intp(ieps,jeps), eps(1), eps(2) =', &
          eps_intp(ieps,jeps),eps_pt(1)%eps,eps_pt(2)%eps
#endif
        POP_SUB(epsmat_intp.interpolate_1D)
        return
      endif

      !take the value of the valid nearest neighbor
      if (pt_avail(2)) then
        eps_intp(ieps,jeps) = eps_pt(2)%eps
      else
        eps_intp(ieps,jeps) = eps_pt(3)%eps
      endif

      !--------------------------+
    case (1) !      1 VALID POINT       |
      !--------------------------+

      !if we have no nearest neighbor, fail
      if (.not.(pt_avail(2).or.pt_avail(3))) then
        ierr=1
        POP_SUB(epsmat_intp.interpolate_1D)
        return
      endif

      !take the value of the valid nearest neighbor
      if (pt_avail(2)) then
        eps_intp(ieps,jeps) = eps_pt(2)%eps
      else
        eps_intp(ieps,jeps) = eps_pt(3)%eps
      endif

      !--------------------------+
    case (0) !      NO VALID POINT      |
      !--------------------------+

#ifdef VERBOSE
      write(0,*) 'WARNING: interpolate_1D: could not find any suitable point'
      write(0,'(1x,a,i4,1x,a,3(F9.6,1x))') 'iq=',iq,'q_eps=',q_intp(:,iq)
      write(0,*)
      write(0,*) '>Fine Point'
      write(0,'(1x,a,3(i5,1x))') 'ieps, jeps, nmtx=',ieps,jeps,nmtx_intp(iq)
      write(0,'(1x,a,2(3(i5,1x),";"))') 'g,gp=',gvec(:,isort_intp(ieps)), gvec(:,isort_intp(jeps))
      write(0,*)
#endif
      ierr = 1
    end select
    POP_SUB(epsmat_intp.interpolate_1D)

  end subroutine interpolate_1D


  ! Interpolate in a plane. Used in situations like:
  ! 42  o   o   41
  ! 
  ! o   32  31  o
  !       x     
  ! o   21  22  o
  ! 
  ! 11  o   o   12
  !------------------------------------------------------------------------+
  subroutine interpolate_2D(ierr,ax)!                                      |
  !------------------------------------------------------------------------+
    integer, intent(out) :: ierr
    integer, intent(in) :: ax(3) !axes that span the interpolation plane

    real(DP) :: qq1(3,2)
    integer :: cell_idx_(3,4,2)
    logical, target  :: pt_avail_(4,2)
    logical, pointer :: pt_avail(:)

    type(eps_point), target  :: eps_pt_(4,2)
    type(eps_point), pointer :: eps_pt(:)
    integer :: ipt, idir
    integer :: sign_x(2) = (/1,-1/)

    real(DP) :: eps_diag(2) !value of eps for a given diagonal
    integer :: eps_cnt(2) !number of points actually used in the diagonal

    PUSH_SUB(epsmat_intp.interpolate_2D)

    pt_avail_(:,:) = .false.

    do idir=1,2
      !qq1 is the coarse point 1 in the diagram above
      qq1(:,idir) = q_intp(:,iq)
      qq1(ax(1),idir) = qq1(ax(1),idir) - 1.5d0*sign_x(idir)*dx !`x` anchors goes up or down
      qq1(ax(2),idir) = qq1(ax(2),idir) - 1.5d0*dx !`y` anchor always goes down
      call get_cell_idx(cells, qq1(:,idir), cell_idx_(:,1,idir))

      do ipt=1, 4
        !for each step, move 1 cell right/up
        if (ipt>1) then
          cell_idx_(:,ipt,idir) = cell_idx_(:, ipt-1,idir)
          cell_idx_(ax(1),ipt,idir) = cell_idx_(ax(1),ipt,idir) + sign_x(idir)
          cell_idx_(ax(2),ipt,idir) = cell_idx_(ax(2),ipt,idir) + 1
        endif

        call get_eps_point(cell_idx_(:,ipt,idir), eps_pt_(ipt,idir))
        pt_avail_(ipt,idir) = eps_pt_(ipt,idir)%avail
      enddo
    enddo

    !loop through diagonals
    do idir=1,2
      !use this so that the macros for eps_pt continue working 
      eps_pt => eps_pt_(:,idir)
      pt_avail => pt_avail_(:,idir)
      eps_cnt(idir) = 0
      eps_diag(idir) = 0d0

      select case (count(pt_avail))

        !--------------------------+
      case (4)  !      4 VALID POINTS      |
        !--------------------------+

        !cubic symmetric interpolation
        eps_diag(idir) = INTP_4(1,2,3,4)
        eps_cnt(idir) = 4

        !--------------------------+
      case (3)  !     3 VALID POINTS       |
        !--------------------------+

        !do we have the nearest points?
        if (pt_avail(2).and.pt_avail(3)) then
          !yes, do asymmetric quadratic interpolation

          if (pt_avail(1)) then
            ! 1   2 x 3   o
            eps_diag(idir) = INTP_3(3,2,1)
          else
            ! o   2 x 3   4
            eps_diag(idir) = INTP_3(2,3,4)
          endif
          eps_cnt(idir) = 3
          cycle
        endif

        !do linear extrapolation
        if(pt_avail(2)) then
          eps_diag(idir) = EXTP_2(2,1)
        else
          eps_diag(idir) = EXTP_2(3,4)
        endif
        eps_cnt(idir) = 2

        !--------------------------+
      case (2)  !     2 VALID POINTS       |
        !--------------------------+

        !if we have no nearest neighbor, fail
        if (.not.(pt_avail(2).or.pt_avail(3))) then
          cycle
        endif

        !do we have both nearest pts?
        if (pt_avail(2).and.pt_avail(3)) then
          !just take the average
          eps_diag(idir) = AVG_2(2,3)
          eps_cnt(idir) = 2
          cycle
        endif

        !take the value of the valid nearest neighbor
        if (pt_avail(2)) then
          eps_diag(idir) = eps_pt(2)%eps
        else
          eps_diag(idir) = eps_pt(3)%eps
        endif
        eps_cnt(idir) = 1

        !--------------------------+
      case (1)  !      1 VALID POINT       |
        !--------------------------+

        !if we have no nearest neighbor, fail
        if (.not.(pt_avail(2).or.pt_avail(3))) cycle

        !take the value of the valid nearest neighbor
        if (pt_avail(2)) then
          eps_diag(idir) = eps_pt(2)%eps
        else
          eps_diag(idir) = eps_pt(3)%eps
        endif
        eps_cnt(idir) = 1

        !--------------------------+
      case (0)  !      NO VALID POINT      |
        !--------------------------+

        ! do nothing, the sum will be zero anyways...
      end select

    enddo

    if (sum(eps_cnt(:))==0) then
      ierr = 1
      POP_SUB(epsmat_intp.interpolate_2D)
      return
    endif

    eps_intp(ieps,jeps) = (eps_diag(1)*eps_cnt(1) + eps_diag(2)*eps_cnt(2))/dble(sum(eps_cnt(:)))

    POP_SUB(epsmat_intp.interpolate_2D)

  end subroutine interpolate_2D


  subroutine find_gvec(iout, vec)
    integer, intent(out) :: iout
    integer, intent(in) :: vec(3)

    integer :: i_2

    PUSH_SUB(epsmat_intp.find_gvec)

    ! I should probably think of a smarter way to check for gvectors.
    ! Some sort of cell scheme could be used...
#ifdef Z_SYMMETRY
    iout = 0
    i_2 = 1
    do while (gvec(1,i_2)/=vec(1).or.gvec(2,i_2)/=vec(2))
      i_2 = i_2 + iabs(gvec(3,i_2)*2) + 1
      if (i_2.gt.ng) return
    enddo
    iout = i_2 + vec(3)-gvec(3,i_2)
#else
    iout = 0
    do i_2=1,ng
      if (all(gvec(:,i_2)==vec)) then
        iout=i_2
        exit
      endif
    enddo
#endif
    !couldn`t find gvector?
    if (iout==0) then
      write(0,*) 'WARNING: find_gvec'
      write(0,'(1x,a,3(I4,1x))') 'wanted=',vec
      write(0,*)
    endif

    POP_SUB(epsmat_intp.find_gvec)

  end subroutine find_gvec

end program epsmat_intp
