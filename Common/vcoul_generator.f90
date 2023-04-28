!==============================================================================
!
! Module:
!
! (1) vcoul_generator()     Originally by JRD   Last Modified: 06/21/2013 (FHJ)
!
!     Generates the (Truncated) Coulomb Interaction for all G at a particular
!     q.  Outputs what would be 8Pi/q^2 if not for truncation.
!
!==============================================================================

#include "f_defs.h"

module vcoul_generator_m

  use global_m
  use minibzaverage_m
  use random_m
  use trunc_cell_box_m
  use trunc_cell_box_d_m
  use trunc_cell_wire_m
  use trunc_scell_box_d_m

  implicit none

  private

  public ::          &
    vcoul_generator, &
    destroy_qran

  real(DP), allocatable, private :: qran(:,:)
  real(DP), private :: q0sph2 = 0d0
  integer, private :: ifirst = 1

contains

  !> FHJ: Calculates the (average) bare Coulomb interaction v(q+G) in Rydberg.
  !! This function uses different forms of v(q+G) (see below) depending on the
  !! dimensionality of the system and whether the potential is truncated. The
  !! Coulomb interaction is optionally averaged on the mini-BZ around the q
  !! vector, and the average screened interaction wcoul0 is also calculated for
  !! q->0 and G=0.
  !!
  !! Note: we use a special flag, peinf%jobtypeeval, to determine whether
  !! vcoul_generator should be a "dumb generator" (if jobtypeeval=0) or 
  !! "smart consumer" (if 1). The main reasons for this division are the following:
  !! 1) For sigma and absorption, the idea is that we are evaluating the
  !! "best estimate" for W(q=0) and v(q=0) for a given k-grid, using the known
  !! analytical limits for epsinv and performing a Monte Carlo integration to
  !! evaluate quantities such as <1/q^2>. Epsilon is always calculated with a
  !! finite but non-averaged v(q), and each program that calculates W determines
  !! if it`s necessary to performing such an averaging.
  !! 2) In sigma, we divide SX into (SX-X) + X, so we have to be consistent with
  !! the way we average the quantities v=<1/q^2> and w=<epsinv/q^2>. In particular,
  !! we "hack" v to make it consistent with W, otherwise the partition (SX-X) + X
  !! would not be correct for metals. We never perform these "hacks" on v(q) 
  !! if peinf%jobtypeeval=0.
  subroutine vcoul_generator(itruncflag, truncval, gvec, bdot, celvol, nfq, ncoul, &
    isrtq, iscreen, qvec, q0vec, vcoul, iwritecoul, iparallel, avgcut, &
    oneoverq, qgrid, epshead, work_scell, averagew, wcoul0, coulomb_mod,nfreq_group)
    !> Truncation flag. The current supported options are:
    !! 0: No truncation (for 3D systems)
    !! 4: Wire truncation (for 1D systems)
    !! 5: Box truncation (for 0D systems)
    !! 6: Slab truncation (for 2D systems)
    integer, intent(in) :: itruncflag
    !> Used for spherical truncation only. 
    real(DP), intent(in) :: truncval(3) 
    !> G space containing the G-vectors used to calculate v(q+G)
    type (gspace), intent(in) :: gvec
    !> Metric of the reciprocal lattice
    real(DP), intent(in) :: bdot(3,3)
    !> Volume of the unit cell
    real(DP), intent(in) :: celvol
    !> FHJ: This is used for miniBZ averages, and only if peinf%jobtypeeval==1
    !! This should always match product(qgrid)
    integer, intent(in) :: nfq
    !> Number of G vectors to calculate v(q+G)
    integer, intent(in) :: ncoul
    !> (ncoul) Indices of the G vectors from gvec which we are using to build v(q+G)
    integer, intent(in) :: isrtq(:)
    !> Screening type. 0 for semiconductor, 1 for graphene, 2 for metal
    integer, intent(in) :: iscreen
    !> The q vector that we use to calculate v(q+G).
    !! The null vector q=(/0,0,0/) is perfectly valid here.
    real(DP), intent(in) :: qvec(3)
    !> The q0 vector that was used in eps0mat, which can never be the null vector.
    real(DP), intent(in) :: q0vec(3)
    !> (ncoul) The output, 8*pi*<1/(q+G)^2> if there`s no truncation.
    real(DP), intent(out) :: vcoul(:)
    !> Write vcoul to a file?
    integer, intent(in) :: iwritecoul
    !> FIXME Not sure what 0 really means if you use MPI!
    integer, intent(in) :: iparallel
    !> If |q+G|^2 < avgcut, calculate <1/(q+G)^2>. Otherwise, calculate 1/(q+G)^2.
    !! The default value for avgcut is TOL_SMALL, i.e., average only done if G=0.
    real(DP), intent(in) :: avgcut
    !> Returns <1/q> or 1/q0
    real(DP), intent(out) :: oneoverq
    !> The q/k-grid being used, which determined the mini-BZ used to construct
    !! averages such as <1/(q+G)^2>. This is only used if peinf%jobtypeeval==1,
    !! except for supercell truncation (untested). In this case, qgrid actualy
    !! means the number of replicas of the system to include.
    integer, intent(in) :: qgrid(3)
    !> Head of epsilon^{-1}. Used to construct wcoul0.
    SCALAR, intent(in) :: epshead 
    !> FIXME A mysterious working buffer for "supercell box truncation"
    type (twork_scell), intent(inout) :: work_scell
    !! If .true., calculate W as <eps*v(q)>
    logical, intent(in) :: averagew
    !> If averagew==.true., wcoul0 is <epsinv*v(q)>. We use the analytical limit of
    !! epsinv(q->0) to perform the proper average. Examples:
    !! - For 3D semiconductors, epsinv(q->0)~const => wcoul0 = epsinv*vcoul(0)
    !! - For 3D metals, epsinv(q->0)~q^2 => wcoul0 = epsinv*8*pi/q0^2, where q0 is
    !! the same q used to evaluate epsinv(q->0).
    SCALAR, intent(inout) :: wcoul0
    !> Optional coulomb modifier parameters. Used in TDDFT and hybrid functional calculations
    !! they allow one to change the bare coulomb interaction. Can only 
    !! be used with no truncation 
    type(coulomb_modifier_t), optional, intent(in) :: coulomb_mod
    !> DVF: Optional number of parallel frequency groups. For when this computation of the 
    !! bare coulomb interaction is for a calculation of epsilon using parallel frequencies
    integer, intent(in), optional :: nfreq_group

    integer :: i1, i2, i3, ii, jj, nn, ig, nint, iCommFlag, iCommFlagGlobal
    integer :: gk(3), seed
    real(DP), allocatable :: vcoul2(:)
    real(DP) :: dvalue, temp_exp, screeninv
    real(DP) :: qk(3), qkxy(3), qlen, q0len, ekinx, kxy
    real(DP) :: kz, zc, qkz(3), dummy, qvec_mod(3)
    real(DP) :: dd(3),fq(3),vlength, lmin, qpp(3) 
    real(DP), allocatable :: ran(:,:)

    real(DP), save :: oneoverq0
    real(DP), save :: vcoul0
    !> Depending on avg_cut, the code needs to average the potential in the
    !! miniBZ. The following flag is set to true after we calculated <v(q)> for
    !! q->0, so we can reuse that integral later for q=0.
    logical, save :: minibz_done = .false.
    !> It there`s no truncation, we also need to calculate the average <1/q>
    !! for the wings. The following flag is set to true after we calculated
    !! <1/q> for q->0.
    logical, save :: oneoverq0_done = .false.
    !> Depending on avg_cut, the code needs to average the modified potential in the
    !! miniBZ. The following flag is set to true after we calculated <v(q)> for
    !! q->0 for the modified potential, so we can reuse that integral later for q=0.
    !! This is separate from minibz_done because you may need to recalculate vcoul 
    !! without the modifier for the exchange part of the kernel etc.
    logical, save :: minibz_mod_done = .false.
    !> DVF: This is needed to accommodate the use of parallel frequencies group in FF 
    !! epsilon calculations. It is just equal to peinf%npes unless you are using
    !! parallel frequency groups in FF epsilon calculations. See the epsilon code
    !! for an explanation of why you need a different value in this case. 
    integer :: npes_local

    ! FHJ: Variables used for minibz calculation
    real(DP) :: UU(3,3), vmid(3)
    integer :: qran_nb, qran_own, qran_first, qran_last, info
    integer, allocatable :: recv_cnts(:), displs(:)

    logical :: verbose

    PUSH_SUB(vcoul_generator)

    verbose = peinf%verb_debug
    if(size(isrtq) < ncoul) call die("vcoul_generator: isrtq not allocated to size ncoul")
    if(size(vcoul) < ncoul) call die("vcoul_generator: vcoul not allocated to size ncoul")

    if(ifirst == 2) call die("you destroyed qran already!")
    ! FHJ: TODO check that nfq == product(qgrid). This should *always*
    ! be the case, since qgrids are Gamma centered. So, we don`t even have
    ! to ask for nfq.
    if(nfq<1 .and. peinf%jobtypeeval==1) call die("Illegal value for nfq in vcoul_generator")

    if(.not. present(nfreq_group)) then
      npes_local=peinf%npes
    else
      npes_local=peinf%npes_orig
    endif


    iCommFlag = 0
    iCommFlagGlobal = 0

    vcoul = 0d0

    nint=100

    qlen  = sqrt(DOT_PRODUCT(qvec,MATMUL(bdot,qvec)))
    q0len = sqrt(DOT_PRODUCT(q0vec,MATMUL(bdot,q0vec)))
    
! JRD: If we have a metal then the divergent components of W cancel
!      and you get a constant.  Epsilon was calculated at q ~ q0vec
! MJ:  Not quite. We while W goes to a constant, V does not. See the
!      comment below about metals. Here we set the q to be qvec (i.e.
!      0 0 0) but later make sure we change W to be evaluated at this
!      q. We change epsiloninv such that epsinv*V(q0vec) = epsnewinv*V(q)
!      where V(q) is the averaged quantity.

    if (abs(qlen)<TOL_ZERO .and. iscreen==SCREEN_METAL .and. &
      peinf%jobtypeeval==1 .and. itruncflag/=TRUNC_NONE) then
      qvec_mod = q0vec
    else
      qvec_mod = qvec
    endif
    
!------------------------------------------------------------------
! Generate random numbers for minibz averaging

    nn = nmc
    
    if (peinf%jobtypeeval==1 .and. ifirst==1 .and. &
     (itruncflag==TRUNC_NONE .or. itruncflag==TRUNC_SLAB)) then

      if (any(qgrid(1:3) .eq. 0)) then
        if(peinf%inode == 0) then
          write(0,*) 'Error qgrid', qgrid
          write(0,*) 'You must specify qgrid in .inp file'
        endif
        call die ('Zero qgrid. Cannot determine minibz', only_root_writes = .true.)
      endif

      SAFE_ALLOCATE(ran, (3,nn))

      if (peinf%inode .eq. 0) then
!     We require a fixed seed for reproducibility of our 'random' numbers.
!     call date_and_time(VALUES=values)
!     seed=((values(3)*24+values(5))*60+values(6))*60+values(7)
        seed = 5000
        call genrand_init(put=seed)
        do jj = 1, 3*nn
          call genrand_real4(dummy)
        enddo
        do jj = 1, nn
          do ii = 1, 3
            call genrand_real4(ran(ii, jj))
          enddo
        enddo
      endif
      
#ifdef MPI
!   Our 'random' numbers must be the same on all processors.
      call MPI_Bcast(ran, 3*nn, MPI_REAL_DP, 0, MPI_COMM_WORLD,mpierr)
#endif
      
      dd(1:3) = 1D0 / dble(qgrid(1:3))

      UU(1:3, 1:3) = bdot(1:3, 1:3)
      ! FHJ: Cholesky decomposition of the metric: bdot = UU^T UU
      call dpotrf('U', 3, UU, 3, info)
      if (info/=0) call die('Could not compute Cholesky decomposition of the metric', &
          only_root_writes=.true.)

      SAFE_ALLOCATE(qran, (3,nn))
      qran = 0d0
      qran_nb = DIVUP(nn, npes_local) ! Block size for qran distrib.
      qran_own = NUMROC(nn, qran_nb, peinf%inode, 0, npes_local) ! How many qran I own
      qran_first = INDXL2G(1, qran_nb, peinf%inode, 0, npes_local) ! First qran I own
      qran_last = qran_first + qran_own - 1 ! Last qran I own
      ! FHJ: OpenMP threading actually hurts performance (at least on Intel i5)
      do jj = qran_first, qran_last
        lmin = INF
        qpp(:) = ran(:,jj)
        do i1 = -ncell+1, ncell
          fq(1) = qpp(1) - dble(i1)
          do i2 = -ncell+1, ncell
            fq(2) = qpp(2) - dble(i2)
            do i3 = -ncell+1, ncell
              fq(3) = qpp(3) - dble(i3)
              !vlength = DOT_PRODUCT(fq,MATMUL(bdot,fq))
              !FHJ: The following lines are mathematically equiv. to the one
              !before, but it`s much faster.
              vmid(1) = UU(1,1)*fq(1) + UU(1,2)*fq(2) + UU(1,3)*fq(3)
              vmid(2) =                 UU(2,2)*fq(2) + UU(2,3)*fq(3)
              vmid(3) =                                 UU(3,3)*fq(3)
              vlength = vmid(1)**2 + vmid(2)**2 + vmid(3)**2
              if (vlength < lmin) then
                lmin = vlength
                qran(:,jj) = fq(:)
              endif
            enddo
          enddo
        enddo
        qran(1:3,jj) = dd(1:3) * qran(1:3,jj)
      enddo ! jj

      if (itruncflag/=TRUNC_SLAB) then
        !
        ! FB: find the (square of the) radius of the sphere contained in the miniBZ
        ! 1/q**2 is going to be calculated analytically in the sphere
        ! and numerically outside
        q0sph2 = INF

        do i1 = -ncell+1, ncell
          fq(1) = dble(i1) * dd(1) * 0.5D0
          do i2 = -ncell+1, ncell
            fq(2) = dble(i2) * dd(2) * 0.5D0
            do i3 = -ncell+1, ncell
              fq(3) = dble(i3) * dd(3) * 0.5D0
              if( i1==0 .AND. i2==0 .AND. i3==0 ) cycle

              vlength = DOT_PRODUCT(fq,MATMUL(bdot,fq))
              if (vlength < q0sph2) then
                q0sph2 = vlength
              endif
            enddo
          enddo
        enddo
      else
        !
        ! FB: find the (square of the) radius of the sphere contained in the miniBZ
        ! 1/q**2 is going to be calculated analytically in the sphere
        ! and numerically outside
        q0sph2 = INF

        do i1 = -ncell+1, ncell
          fq(1) = dble(i1) * dd(1) * 0.5D0
          do i2 = -ncell+1, ncell
            fq(2) = dble(i2) * dd(2) * 0.5D0
            if( i1==0 .AND. i2==0 ) cycle

            vlength = DOT_PRODUCT(fq,MATMUL(bdot,fq))
            if (vlength < q0sph2) then
              q0sph2 = vlength
            endif
          enddo
        enddo

      endif

      SAFE_DEALLOCATE(ran)

#ifdef MPI
      if (iparallel .eq. 1) then
        SAFE_ALLOCATE(recv_cnts, (npes_local))
        SAFE_ALLOCATE(displs, (npes_local))
        do jj = 1, npes_local
          recv_cnts(jj) = NUMROC(3*nn, 3*qran_nb, jj-1, 0, npes_local)
          displs(jj) = INDXL2G(1, 3*qran_nb, jj-1, 0, npes_local) - 1
        enddo
        call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, qran, &
          recv_cnts, displs, MPI_REAL_DP, MPI_COMM_WORLD, mpierr)
        SAFE_DEALLOCATE(recv_cnts)
        SAFE_DEALLOCATE(displs)
      endif
#endif

    endif

!-------------------------------------------------------------------
! No Truncation

    if (itruncflag==TRUNC_NONE .and. (.not. present(coulomb_mod))) then

! Calculate Wing Correction Factor - this is not done for Epsilon and Kernel
! since avgcut is zero and qvec_mod is not

      if (peinf%inode==0 .and. qlen**2<avgcut .and. peinf%jobtypeeval==1) then
        if (qlen**2>TOL_ZERO .or. .not.oneoverq0_done) then
          call minibzaverage_3d_oneoverq(nn,bdot,dvalue,qran,qvec_mod)
          oneoverq=dvalue
          if (qlen**2<=TOL_ZERO) then
            oneoverq0=dvalue
            oneoverq0_done = .true.
          endif
        else
          oneoverq=oneoverq0
        endif
      endif
      ! otherwise we set oneoverq a little later

      do ig=1,ncoul
        if (iparallel .eq. 1) then
          if(mod(ig-1,npes_local).ne.peinf%inode) cycle
        endif
        
        gk(:)=gvec%components(:,isrtq(ig))
        qk(:)=gk+qvec_mod(:)
        ekinx=DOT_PRODUCT(qk,MATMUL(bdot,qk))

! We Do 3D Mini Brillouin Zone Average if q is exactly 0
! and G = 0. This should be the case when constructing W,
! but don`t want this in Epsilon code for example, where you
! have to use a finite Q.

        if (ekinx<avgcut .and. peinf%jobtypeeval==1) then
          ! FHJ: this branch means that we |G+k|^2 is small enough that
          ! we want to do miniBZ average, and that we are running a kind of
          ! job (jobtypeeval==1) where we want to do miniBZ averagesa.
          select case (iscreen)
            !-------------------------------------------------------------------
            case (SCREEN_SEMICOND)
              if (ekinx>TOL_ZERO .or. .not.minibz_done) then 
                if (ekinx<=TOL_ZERO) then
                  call minibzaverage_3d_oneoverq2(nn,bdot,dvalue,qran,qk,averagew,epshead,wcoul0,q0sph2,celvol,nfq)
                  vcoul(ig)=dvalue
                  vcoul0=dvalue
                  iCommFlag = peinf%inode+1
                  minibz_done = .true.
                else
                  call minibzaverage_3d_oneoverq2(nn,bdot,dvalue,qran,qk,averagew,epshead,wcoul0,q0sph2,celvol,nfq)
                  vcoul(ig)=dvalue
                endif
              else
                vcoul(ig)=vcoul0
              endif
            !-------------------------------------------------------------------
            case (SCREEN_GRAPHENE)
              if (ekinx>TOL_ZERO .or. .not.minibz_done) then
                call minibzaverage_3d_oneoverq2(nn,bdot,dvalue,qran,qk,averagew,epshead,wcoul0,q0sph2,celvol,nfq)
                if (ekinx<TOL_ZERO) then
                  vcoul0=dvalue
                  call minibzaverage_3d_oneoverq(nn,bdot,dvalue,qran,qk)
                  if (q0len .lt. TOL_ZERO) then
                    write(0,'(a)') 'You have q0vec=0 but a graphene-type system!!'
                    call die('Bad q0vec', only_root_writes = .true.)
                  endif
                  wcoul0=epshead*dvalue/q0len
                  iCommFlag = peinf%inode+1
                  minibz_done = .true.
                  vcoul(ig)=vcoul0
                else
! JRD: This seems wrong. We now use correct vcoul average.
!                vcoul(ig)=dvalue/sqrt(ekinx)
                endif
              else
                vcoul(ig)=vcoul0
              endif
            !-------------------------------------------------------------------
            case (SCREEN_METAL)
! MJ : W(q->0) goes to a constant. But the COH term (as well as the
!      Fock term) also require V(q->0). It is simpler to see this in the 
!      context of the COH term in COHSEX which is 0.5*(W-V). Setting V(q->0)
!      based on the q0vec is dangerous because it just adds a constant to all
!      the quasiparticle levels. So we keep both of these quantities
!      at hand. This is only needed in the Sigma code.
               ! write(0,'(a)') 'You want cell averaging on a metal!!'
               ! write(0,'(a)') 'Specify q0vec and cell averaging will automatically'
               ! write(0,'(a)') 'be done for vcoul. Wcoul will use the q0vec that you'
               ! write(0,'(a)') 'specify!!'
               ! call die('Bad Screening Options', only_root_writes = .true.)
              if (ekinx>TOL_ZERO .or. .not.minibz_done) then
                call minibzaverage_3d_oneoverq2(nn,bdot,dvalue,qran,qk,averagew,epshead,wcoul0,q0sph2,celvol,nfq)

!      Keep in mind that while the vcoul average out of minibzaverage_3d_oneoverq2
!      is correct, wcoul0 is not. wcoul0 is just modified so that when multiplied
!      by vcoul it gives the correct result.
                !if (abs(ekinx-q0len) .lt. TOL_ZERO) then
                if (ekinx<TOL_ZERO) then
                  vcoul0=dvalue
                  !wcoul0=epshead*(q0len*q0len)/(dvalue/8.0/PI_D)
                  wcoul0=epshead*8.0*PI_D/(q0len*q0len)
                  iCommFlag = peinf%inode+1
                  minibz_done = .true.
                  vcoul(ig)=vcoul0
                endif
              else
                vcoul(ig)=vcoul0
              endif
          endselect

        else if (ekinx<TOL_ZERO) then
          ! FHJ: |G+k|^2 is zero, but we don`t want to do miniBZ averages,
          ! because this is not an "evaluation job".
          ! So, just set v(q+G=0)=0. The outside code should do something about this value.
          vcoul(ig) = 0D0
        else
          ! FHJ: |G+q|^2 is large. Just use the regular expression for the
          ! Coulomb potential, v(q+G) = 4*Pi*e^2 / |q+G|^2
          ! MJ : Note that, for metals, we already changed qvec_mod to q0vec above
          ! FHJ: (note sure what that means...)
          vcoul(ig) = 8.0d0*PI_D/ekinx
        endif
      enddo

      if (qlen**2>=avgcut) then
        ! FHJ: For |q|^2 larger than the miniBZ average, 1/|q| is computed
        ! without miniBZ averages, as |q| * v(q).
        oneoverq=vcoul(1)*qlen
      endif

#ifdef MPI
      if (iparallel .eq. 1) then
        call MPI_Bcast(oneoverq, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)

        SAFE_ALLOCATE(vcoul2, (ncoul))
        vcoul2(:)=0.0d0
        call MPI_ALLREDUCE(vcoul,vcoul2,ncoul, &
          MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD,mpierr)
        call MPI_ALLREDUCE(iCommFlag,iCommFlagGlobal,1, &
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)
        if (iCommFlagGlobal .ne. 0) then
          minibz_done = .true.
          call MPI_BCAST(vcoul0,1,MPI_DOUBLE_PRECISION,iCommFlagGlobal-1, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(wcoul0,1,MPI_SCALAR,iCommFlagGlobal-1, MPI_COMM_WORLD, mpierr)
        endif
        vcoul(:)=vcoul2(:)
        SAFE_DEALLOCATE(vcoul2)
      endif
#endif

!#BEGIN_INTERNAL_ONLY
! This is what we do if we have TDDFT where there are some modifiers for the
! Coulomb interaction
    elseif (itruncflag==TRUNC_NONE .and. present(coulomb_mod)) then

      screeninv = 1.0D0/(4.0D0 * coulomb_mod%screening_length *coulomb_mod%screening_length)
      ! convert screening_length from A^{-1} to Bohr
      screeninv = screeninv/(BOHR*BOHR)
      do ig=1,ncoul
        if (iparallel .eq. 1) then
          if(mod(ig-1,npes_local).ne.peinf%inode) cycle
        endif
        
        gk(:)=gvec%components(:,isrtq(ig))
        qk(:)=gk+qvec_mod(:)
        ekinx=DOT_PRODUCT(qk,MATMUL(bdot,qk))

! We Do 3D Mini Brillouin Zone Average if q is exactly 0
! and G = 0. This should be the case when constructing W,
! but don`t want this in Epsilon code for example, where you
! have to use a finite Q.

        if (ekinx<avgcut .and. peinf%jobtypeeval==1) then
          if (iscreen==SCREEN_SEMICOND) then ! semiconductor
            if (ekinx>TOL_ZERO .or. .not.minibz_mod_done) then
              call minibzaverage_3d_oneoverq2_mod(nn,bdot,dvalue,qran,qk,coulomb_mod)
              vcoul(ig)=dvalue
              if (ekinx<=TOL_ZERO) then
                vcoul0=dvalue
                iCommFlag = peinf%inode+1
                minibz_mod_done = .true.
              endif
            else
              vcoul(ig)=vcoul0
            endif
         else ! iscreen /= SCREEN_SEMICOND
            call die('No trunc TDDFT with metal/graphene not supported', only_root_writes = .true.)
         endif

        else if (ekinx<TOL_ZERO) then

! We are not in an evaluation job and were sent q=0, we don`t average. We just zero. The outside code should
! do something about this value.

          vcoul(ig) = 0D0

        else

          temp_exp = exp(-screeninv*ekinx)
          vcoul(ig) = 8.0d0*PI_D/ekinx* (coulomb_mod%long_range_frac_fock*temp_exp + &
                      coulomb_mod%short_range_frac_fock*(1.0d0-temp_exp))
        endif
      enddo

      if (qlen**2 .ge. avgcut) then
        oneoverq=vcoul(1)*qlen
      endif

#ifdef MPI
      if (iparallel .eq. 1) then
        call MPI_Bcast(oneoverq, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)

        SAFE_ALLOCATE(vcoul2, (ncoul))
        vcoul2(:)=0.0d0
        call MPI_ALLREDUCE(vcoul,vcoul2,ncoul, &
          MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD,mpierr)
        call MPI_ALLREDUCE(iCommFlag,iCommFlagGlobal,1, &
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)
        if (iCommFlagGlobal .ne. 0) then
          minibz_done = .true.
          call MPI_BCAST(vcoul0,1,MPI_DOUBLE_PRECISION,iCommFlagGlobal-1, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(wcoul0,1,MPI_SCALAR,iCommFlagGlobal-1, MPI_COMM_WORLD, mpierr)
        endif
        vcoul(:)=vcoul2(:)
        SAFE_DEALLOCATE(vcoul2)
      endif
#endif
!#END_INTERNAL_ONLY
!--------------------------------------------------------------------
! Rectangular Box Truncation

    elseif (itruncflag .eq. 1) then
    
      call die('Rectangular Truncation is No Longer Supported', only_root_writes = .true.)
    
!------------------------------------------------------------------
! Spherical Truncation

    elseif (itruncflag==TRUNC_SPHERICAL .and. (.not. present(coulomb_mod))) then
      
      do ig=1,ncoul
        qk(:)=gvec%components(:,isrtq(ig))+qvec_mod(:)
        ekinx=DOT_PRODUCT(qk,MATMUL(bdot,qk))
        if ( ekinx.lt.TOL_ZERO .and. peinf%jobtypeeval .eq. 1 ) then
          vcoul(ig)=8.0D0*PI_D*((truncval(1))**2)/2D0
        else if ( ekinx.lt.TOL_ZERO ) then
          vcoul(ig) = 0D0
        else
          vcoul(ig)=8.0d0*PI_D/ekinx* &
            (1.0d0-cos(truncval(1)*sqrt(ekinx)))
        endif
      enddo
      
      oneoverq=vcoul(1)*qlen
      wcoul0 = vcoul(1)

!#BEGIN_INTERNAL_ONLY
    elseif (itruncflag==TRUNC_SPHERICAL .and. present(coulomb_mod)) then

      ! This is WRONG. Basically we need to do a convolution in G-space
      ! not a product

      screeninv = 1.0D0/(4.0D0 * coulomb_mod%screening_length *coulomb_mod%screening_length)
      ! convert screening_length from A^{-1} to Bohr
      screeninv = screeninv/(BOHR*BOHR)

      do ig=1,ncoul
        qk(:)=gvec%components(:,isrtq(ig))+qvec_mod(:)
        ekinx=DOT_PRODUCT(qk,MATMUL(bdot,qk))
        if ( ekinx.lt.TOL_ZERO .and. peinf%jobtypeeval .eq. 1 ) then
          vcoul(ig)=8.0D0*PI_D*coulomb_mod%long_range_frac_fock*((truncval(1))**2)/2D0
        else if ( ekinx.lt.TOL_ZERO ) then
          vcoul(ig) = 0D0
        else
          temp_exp = dexp(-ekinx*screeninv)

          vcoul(ig) = 8.0d0*PI_D/ekinx* (coulomb_mod%long_range_frac_fock*temp_exp + &
                      coulomb_mod%short_range_frac_fock*(1.0d0-temp_exp)) * (1.0d0-cos(truncval(1)*sqrt(ekinx)))
        endif
      enddo
      
      oneoverq=vcoul(1)*qlen
      wcoul0 = vcoul(1)
!#END_INTERNAL_ONLY

!------------------------------------------------------------------
! Cylindrical Box Truncation

    elseif (itruncflag .eq. 3) then

      call die('Cylindrical Truncation is No Longer Supported', only_root_writes = .true.)
      
!-----------------------------------------------------------------
! Cell Wire Truncation

    elseif (itruncflag==TRUNC_WIRE) then
      
!    if (peinf%inode .eq. 0) then
!      write(6,*) 'Generating Vcoul_cell_wire with FFT'
!    endif

! JRD: This is hopefully temporary.  q=0 gives NaN
! because the log diverges in Eq. 5 of Ismail-Beigi paper.  
! For all G =/ 0, using the q0vec (small shifted vector) is 
! probably good enough.  
!
! JRD: This is now fixed in trunc_cell_wire where we use the scheme
! of Ismail-Beigi

      if (iparallel .eq. 1) then
        call trunc_cell_wire(gvec,verbose,peinf%inode, &
          npes_local,bdot,qvec_mod(3),ncoul,isrtq,vcoul)
      else
        call trunc_cell_wire(gvec,verbose,0, &
          1,bdot,qvec_mod(3),ncoul,isrtq,vcoul)
      endif

! We Do 1D Mini Brillouin Zone Average if q=0, G=0 or less than avgcut
! May not be implemented correctly for graphene-type system... I`m not even sure there is a system with linear DOS in 1D...

      if (iscreen==SCREEN_SEMICOND) then ! semiconductor
        do ig=1,ncoul
          qk(:)=gvec%components(:,isrtq(ig))+qvec_mod(:)
          ekinx=DOT_PRODUCT(qk,MATMUL(bdot,qk))
          if ( ekinx .lt. avgcut .and. peinf%jobtypeeval .eq. 1) then 

            if (ekinx .gt. TOL_ZERO .or. .not.minibz_done) then
              call minibzaverage_1d(gvec,nfq,bdot,dvalue,iparallel,qk,epshead,q0len,averagew,wcoul0)
              vcoul(ig)=dvalue
              if (ekinx .lt. TOL_ZERO) then
                vcoul0=dvalue
                minibz_done = .true.
              endif
            else
              vcoul(ig)=vcoul0
            endif
          else if (ekinx<TOL_ZERO) then
            vcoul(ig) = 0
          endif
        enddo
      endif

      oneoverq=vcoul(1)*qlen

!----------------------------------------------------------------
! Cell Box Truncation

    elseif (itruncflag==TRUNC_BOX) then
      
      if (qlen>TOL_ZERO) then
        write(0,'(a)') 'You asked for cell box truncation but have more q-points than q=0!!'
        call die('Bad Truncation', only_root_writes = .true.)
      endif
      
      if (iparallel .eq. 1) then
        if (.not. present(coulomb_mod)) then
          call trunc_cell_box_d(gvec,verbose,peinf%inode, &
            npes_local,bdot,ncoul,isrtq,vcoul)
!#BEGIN_INTERNAL_ONLY
        else
          call trunc_cell_box_d(gvec,verbose,peinf%inode, &
            npes_local,bdot,ncoul,isrtq,vcoul,coulomb_mod)
!#END_INTERNAL_ONLY
        endif
      else
        if (.not. present(coulomb_mod)) then
          call trunc_cell_box(gvec,verbose,bdot,ncoul,isrtq,vcoul)
!#BEGIN_INTERNAL_ONLY
        else
          call trunc_cell_box(gvec,verbose,bdot,ncoul,isrtq,vcoul,coulomb_mod)
!#END_INTERNAL_ONLY
        endif
      endif
      
      oneoverq=vcoul(1)*qlen
      wcoul0 = vcoul(1)

!----------------------------------------------------------------
! Cell Slab Truncation

! JRD: This is easy because an analytic expression exists.  See Sohrab.

    elseif (itruncflag==TRUNC_SLAB) then
      
      if (abs(qvec_mod(3))>TOL_ZERO) then
        write(0,'(a)') 'You asked for cell slab truncation but have more q-points in z direction than qz=0!!'
        call die('Bad Truncation', only_root_writes = .true.)
      endif
      
      do ig=1,ncoul
        
        if (iparallel .eq. 1) then
          if(mod(ig,npes_local).ne.peinf%inode) cycle
        endif
        
        qk(:)=gvec%components(:,isrtq(ig))+qvec_mod(:)
        ekinx=DOT_PRODUCT(qk,MATMUL(bdot,qk))
        qkxy(1:2)=qk(1:2)
        qkxy(3)=0D0
        kxy=sqrt(DOT_PRODUCT(qkxy,MATMUL(bdot,qkxy)))
        qkz(1:2)=0D0
        qkz(3)=qk(3)
        kz=sqrt(DOT_PRODUCT(qkz,MATMUL(bdot,qkz)))
        zc=2D0*PI_D/(sqrt(bdot(3,3))*2D0)
        
!          write(6,*) "zc", zc

        if (ekinx<avgcut .and. peinf%jobtypeeval==1) then
          if (iscreen/=SCREEN_METAL) then ! semiconductor or graphene
            if (ekinx>=TOL_ZERO .or. .not.minibz_done) then
              call minibzaverage_2d_oneoverq2(nn,bdot, &
                dvalue,qran,qk,kz,zc,epshead,q0len,averagew,wcoul0)
              vcoul(ig)=dvalue
              if (ekinx .lt. TOL_ZERO) then
                vcoul0=dvalue
                if (iscreen==SCREEN_GRAPHENE) then
                  wcoul0=vcoul0*epshead
                endif
                minibz_done = .true.
                iCommFlag=peinf%inode+1
              endif
            else
              vcoul(ig)=vcoul0
            endif
          else
            write(0,'(a)') 'You have q0vec=0 but a metal!!'
            call die('Bad Screening Options', only_root_writes = .true.)
          endif
        else if (ekinx<TOL_ZERO) then
          vcoul(ig) = 0D0
        else
          vcoul(ig) = 8.0d0*PI_D/ekinx
          vcoul(ig) = vcoul(ig)*(1.0d0-exp(-kxy*zc)*cos(kz*zc))
        endif
      enddo
      
#ifdef MPI
      if (iparallel .eq. 1) then
        SAFE_ALLOCATE(vcoul2, (ncoul))
        vcoul2(:)=0.d0
        call MPI_ALLREDUCE(vcoul,vcoul2,ncoul, &
          MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD,mpierr)
        vcoul(:)=vcoul2(:)
        SAFE_DEALLOCATE(vcoul2)
        call MPI_ALLREDUCE(iCommFlag,iCommFlagGlobal,1, &
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)
        if (iCommFlagGlobal .ne. 0) then
          minibz_done = .true.
          call MPI_BCAST(vcoul0,1,MPI_DOUBLE_PRECISION,iCommFlagGlobal-1, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(wcoul0,1,MPI_SCALAR,iCommFlagGlobal-1, MPI_COMM_WORLD, mpierr)
        endif
      endif
#endif
      
! This is wrong too?

      oneoverq=vcoul(1)*qlen

!----------------------------------------------------------------
! Supercell Box Truncation

    elseif (itruncflag==TRUNC_SUPERCELL) then
      
      if (iparallel .eq. 1) then
        call trunc_scell_box_d(gvec,.true.,peinf%inode,npes_local, &
          bdot,qvec_mod,qgrid,ncoul,isrtq,vcoul,work_scell)
      else
        call trunc_scell_box_d(gvec,.true.,0,1, &
          bdot,qvec_mod,qgrid,ncoul,isrtq,vcoul,work_scell)
      endif
      
      oneoverq=vcoul(1)*qlen
      wcoul0 = vcoul(1)
      
    endif

! Saving qran between calls
!      deallocate(qran)

!-----------------------------------------------------------------
! Print vcoul to file

    if (iwritecoul .eq. 1) then
      if (peinf%inode.eq.0) then
        do ig=1,ncoul
          write(19,'(3f12.8,1x,3i7,1x,e20.8)') &
            qvec_mod(:),gvec%components(:,isrtq(ig)),vcoul(ig)
        enddo
      endif
    endif

    ifirst = 0
    
    POP_SUB(vcoul_generator)
    
    return
  end subroutine vcoul_generator

!-----------------------------------------------------------------
  subroutine destroy_qran()

    PUSH_SUB(destroy_qran)

    ifirst = 2
    SAFE_DEALLOCATE(qran)

    call logit('Deallocated random numbers.')

    POP_SUB(destroy_qran)
    return
  end subroutine destroy_qran

end module vcoul_generator_m
