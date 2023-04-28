!=================================================================================
!
! Routines:
!
! (1) mtxel_kernel()         Originally By MLT     Last Modified  02/02/2013 (FHJ)
!
!     input: crys, gvec, syms, qg, wfnc, wfncp, wfnv, wfnvp,
!            xct types
!            ii  = label of ck block
!            ipe = label of the PE that contains the ckp block
!     output: bsedbody,bsedhead,bsedwing,bsex = kernel matrix elements
!              between ck (all v) and ckp (all vp) blocks
!
!     Calculate the head, wings, body, exchange of the kernel (see
!     eq. 34, 35, and 41-46, Rohlfing & Louie). The exchange has just
!     the proper part, no divergent contribution.
!
!=================================================================================

#include "f_defs.h"

module mtxel_kernel_m

  use global_m
  use fftw_m
  use gmap_m
  use misc_m
  use w_sum_m
  use g_sum_m
  use gx_sum_m
  use accel_memory_m
  use accel_fft_m
  use algos_kernel_m

  implicit none

  public :: mtxel_kernel

  private

  !> FHJ: pointer to a wavefunction. We use this construction b/c Fortran, in
  !! all its infinite might and wisdom, doesn`t allow arrays of pointers.
  type wavefunction_ptr
    type(wavefunction), pointer :: p
  end type wavefunction_ptr

  !> FHJ: contains the indices of G, G-G0, -G and G0-G. Used by charge_matrix
  !! depending on whether there`s umklapp involved.
  type gvec_indices
    !> Index of G, G-G0, -G, and G0-G
    integer, pointer :: gmg0(:), g0mg(:), g(:), mg(:)
  end type gvec_indices

  type matrix_4
    SCALAR, pointer :: p(:,:,:,:) => Null()
  end type matrix_4

contains

  subroutine mtxel_kernel(crys,gvec,syms,qg,wfnc,wfncp,wfnvp, &
    wfnv,xct,leading_dim,bsedbody,bsedhead,bsedwing,bsex,ii,ik,ikp, &
    ic_in,icp_in,iv_in,ivp_in,vcoularray,fq,qq,g0,ifq,irq,q0len)

    type (crystal), intent(in) :: crys
    type (gspace), intent(in) :: gvec
    type (symmetry), intent(in) :: syms
    type (grid), intent(in) :: qg
    type (wavefunction), target, intent(in) :: wfnc,wfncp,wfnvp,wfnv
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: leading_dim
    SCALAR, intent(inout) :: bsedbody(:,:,:), bsedhead(:,:,:), bsedwing(:,:,:), &
      bsex(:,:,:) !< (leading_dim,xct%nspin,xct%nspin)
    integer, intent(in) :: ii,ik,ikp,ic_in,icp_in,iv_in,ivp_in
    real(DP), intent(in) :: vcoularray(:,:) !< (xct%ng,qg%nf)
    real(DP), intent(in) :: fq(3),qq
    integer, intent(in) :: g0(3)
    integer, intent(in) :: ifq !< q-pt index in full BZ
    integer, intent(in) :: irq !< k-pt index in reduced BZ
    real(DP), intent(in) :: q0len

    integer :: ipe
    character :: filename*20
    integer :: jj,sinv
    integer :: ig,igp,igp_loc
    integer :: neps,nepsmin,ngpown_ipe
    !> Number of g-vectors in the gvec G-space needed to describe the dielectric matrix.
    integer :: ng_eps
    !> Number of conduction/valence bands we are dealing with. nbands is an
    !! array with (/ invband, incband /)
    integer :: incband, invband, nbands(2)
    !> The "left" and "right" set of WFNs that we use to get charge density matrices.
    !! The first index represents a v=1 or c=2 type of WFN, and the second
    !! index is whether it is at k(=1) or kp(=2).
    integer :: wfn_l(2), wfn_r(2)
    real(DP) :: qlen,vq(3)

    integer, save :: ik_old = 0
    integer, save :: ikp_old = 0
    integer, save :: iv_old = 0
    integer, save :: ivp_old = 0
    integer, save :: ic_old = 0
    integer, save :: icp_old = 0
    logical :: ivsaved, icsaved, ixsaved, ixpsaved

    integer, allocatable :: isrtqi(:), isrtqt(:), ind(:), indinv(:)
    real(DP), allocatable :: vcoul(:)
    SCALAR :: wval,epshead
    SCALAR, allocatable :: &
      epscol(:),wptcol(:),ph(:),phinv(:),epscolt(:,:),epscolt2(:,:), &
      tempw(:,:,:,:),tempb(:,:,:,:),temph(:,:,:)
    SCALAR, pointer :: mccp(:,:,:,:),mvvp(:,:,:,:),mvc(:,:,:,:),mvpcp(:,:,:,:)
    SCALAR, save, allocatable :: &
      tempw_old(:,:,:,:),tempb_old(:,:,:,:),temph_old(:,:,:)
    SCALAR, save, allocatable :: &
      mccp_old(:,:,:,:),mvc_old(:,:,:,:),mvpcp_old(:,:,:,:)

    SCALAR, allocatable :: wptcol_mat(:,:)

    integer :: irqt
    real(DP) :: fact

    complex(DPC), dimension(:,:,:), allocatable :: fftbox1,fftbox2
    real(DP) :: scale
    type(gvec_indices) :: g_idx
    integer :: Nfft(3)
    integer :: itotj
    logical :: save_ffts !< reuse WFNs ffts within the subroutine?
    !> If we are saving the ffts, precalculate the fftboxes for all WFNs.
    !! Indices are (x, y, z, spin, band, ik or ikp?)
    complex(DPC), dimension(:,:,:,:,:,:), allocatable :: fftboxes
    ! (c/p, k/kp) points to the appropriate WFW (wfnv, wfnc, wfnvp, wfncp)
    type(wavefunction_ptr) :: wfns(2,2)
    type(matrix_4), target :: mats(2,2,2,2)

    PUSH_SUB(mtxel_kernel)

    !==============================================================================
    ! Initialization: determine if we can reuse previous matrix elements, etc.
    !==============================================================================

    save_ffts = .false.
    wfns(1,1)%p => wfnv
    wfns(1,2)%p => wfnvp
    wfns(2,1)%p => wfnc
    wfns(2,2)%p => wfncp
    ivsaved = .false.
    ixsaved = .false.
    ixpsaved = .false.
    icsaved = .false.
    if ( ik /= -1) then
      if (xct%ivpar .eq. 0) invband=xct%nvb_co
      if (xct%ivpar .eq. 1) invband=1
      if (xct%icpar .eq. 0) incband=xct%ncb_co
      if (xct%icpar .eq. 1) incband=1
      nbands = (/ invband, incband /)

      if (ik_old .eq. ik .and. ikp_old .eq. ikp) then
        if (iv_old .eq. iv_in .and. ic_old .eq. ic_in) then
          ixsaved = .true.
          !            if (peinf%inode .eq. 0) write(6,781) ii,iv_old,iv_in,ic_old,ic_in
        endif
        if (ivp_old .eq. ivp_in .and. icp_old .eq. icp_in) then
          ixpsaved = .true.
          !            if (peinf%inode .eq. 0) write(6,782) ii,ivp_old,ivp_in,icp_old,icp_in
        endif
        if (iv_old .eq. iv_in .and. ivp_old .eq. ivp_in) then
          ivsaved = .true.
          !            if (peinf%inode .eq. 0) write(6,783) ii,iv_old,iv_in,ivp_old,ivp_in
        else
          iv_old=iv_in
          ivp_old=ivp_in
        endif
        if (ic_old .eq. ic_in .and. icp_old .eq. icp_in) then
          icsaved = .true.
          !            if (peinf%inode .eq. 0) write(6,784) ii,ic_old,ic_in,icp_old,icp_in
        else
          ic_old=ic_in
          icp_old=icp_in
        endif
      else
        ic_old=ic_in
        icp_old=icp_in
        iv_old=iv_in
        ivp_old=ivp_in
        ik_old=ik
        ikp_old=ikp
      endif
      !781     format('Reusing Exchange Matrix Elements',5i6)
      !782     format('Reusing Exchange Matrix P Elements',5i6)
      !783     format('Reusing Valence Matrix Elements',5i6)
      !784     format('Reusing Conduction Matrix Elements',5i6)

      sinv = 1
      if (ik.gt.ikp) sinv = -1
      if (.not. ivsaved) then
        SAFE_ALLOCATE(isrtqi, (gvec%ng))
        isrtqi=0
      endif
      call init_g_idx(xct, gvec, g0, sinv, g_idx)
      call timacc(61,1)
    endif

    !==============================================================================
    ! Initialize epsilon and v(q)
    !==============================================================================
    ! Read dielectric matrix at q. The umklapp vector g0 is k - kp = q + g0
    ! JRD: All procs need to participate in the communication
    if (xct%bLowComm) then
      if ( ik /= -1 .and. (.not. ivsaved)) then
        isrtqi(:)=xct%isrtqi(:,irq)
      endif
      ! JRD NEED GLOBAL IVSAVED HERE
    else if (.true.) then
      SAFE_ALLOCATE(isrtqt, (gvec%ng))
      do irqt = 1, qg%nr
        if (peinf%inode .eq. 0) then
          isrtqt(:)=xct%isrtqi(:,irqt)
        endif
#ifdef MPI
        call MPI_BCAST(isrtqt,gvec%ng,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
        if (irq .eq. irqt .and. ik /= -1 .and. (.not. ivsaved)) isrtqi(:) = isrtqt(:)
      enddo
      SAFE_DEALLOCATE(isrtqt)
    endif

    if ( ik /= -1 .and. (.not. ivsaved)) then
      neps=xct%nmtxa(irq)
      nepsmin = max(xct%neps, neps)
      SAFE_ALLOCATE(epscol, (neps))
      call timacc(61,2)

      ! Compute Coulomb interaction at this q vector: vcoul(q+G) for q+G=0 we set
      ! the interaction to V(q0) where q0 is the small vector used for epsilon
      SAFE_ALLOCATE(vcoul, (xct%ng))
      call get_vcoul(.false., .false.,fq)

      ! FHJ: There are two G-spaces here, which historically caused several bugs:
      ! - The gvec space: vectors ig_gvec sorted wrt |G|^2; and
      ! - The epsilon space: vectors ig_eps sorted wrt |R(G+q)|^2
      ! Different vectors/matrices will be ordered wrt a different space:
      ! - gvec space: M12 matrices, temp{w/b}, wptcol, vcoul
      ! - eps space: epscol
      ! To get from one space to the other use the following mapping arrays:
      !  ig_eps = ind(ig_gvec); ig_gvec = indinv(ig_eps)
      ! The code below loops over all xct%ng vectors ig_gvec from gvec G-space to
      ! find the mapping to the eps G-space.
      ! IMPORTANT: ind and indinv will always return a valid index. If the mapping
      ! is actually missing, we zero out the phases ph(ig_gvec)/phinv(ig_eps).
      ! So, always use phinv(ig_eps), and not ph(indinv(ig_eps))!
      call timacc(69,1)
      igp=0
      SAFE_ALLOCATE(ind, (xct%ng))
      SAFE_ALLOCATE(indinv, (nepsmin))
      SAFE_ALLOCATE(ph, (xct%ng))
      SAFE_ALLOCATE(phinv, (nepsmin))
      ind(:)=0
      indinv(:)=0
      ph(:)=ZERO
      phinv(:)=ZERO

      call gmap(gvec, syms, xct%ng, qg%itran(ifq), qg%kg0(:,ifq), g_idx%g, &
        isrtqi, ind, ph, xct%die_outside_sphere)
      if(ind(1)==0 .or. ind(1)>neps) then
        call die('Could not map the head of the dielectric matrix!', &
          only_root_writes=.true.)
      endif

      ! FHJ: find the inverse mapping and the maximum number of G-vectors
      ! that we need in order to describe wptcoul(ig_gvec).
      ng_eps = 0
      do ig=1,xct%ng
        if (ind(ig)>0 .and. ind(ig)<=neps) then
          indinv(ind(ig)) = ig
          phinv(ind(ig)) = ph(ig)
          ng_eps = ig
        else
          ind(ig) = neps
          ph(ig) = ZERO
        endif
      enddo
      do ig=1,neps
        if (indinv(ig)>0 .and. indinv(ig)<=xct%ng) then
          if (ind(indinv(ig))/=ig) &
            call die('Error in mtxel_kernel: mapping from eps to g-space not consistent!')
        else
          indinv(ig) = xct%ng
          phinv(ig) = ZERO
          call warn_about_mapping()
        endif
      enddo

      SAFE_ALLOCATE(wptcol, (ng_eps))
      wptcol=ZERO

      if (w_sum_algo /= CPU_ALGO) then
        SAFE_ALLOCATE(wptcol_mat, (ng_eps, xct%ngpown_max))
        wptcol_mat = ZERO
      end if

      call timacc(69,2)
      call timacc(70,1) ! why are we timing retrieving a single value from an array...?

      ! We need to now find the head of epsilon
      epshead=xct%epsdiag(ind(1),irq)

      call timacc(70,2)
    endif

    !==============================================================================
    ! Prepare fftboxes and precalculate WFN FFTs, if there`s enough memory.
    !==============================================================================
    if ( ik /= -1 ) then
      ! Compute size of FFT box we need
      call timacc(34,1)
      call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)
      ! Allocate FFT boxes
      if (xct%ilowmem==-1) then
        save_ffts = .true.
        SAFE_ALLOCATE(fftboxes, (Nfft(1),Nfft(2),Nfft(3),xct%nspin*xct%nspinor,invband+incband,2))
        call precompute_ffts()
      endif
      SAFE_ALLOCATE(fftbox1, (Nfft(1),Nfft(2),Nfft(3)))
      SAFE_ALLOCATE(fftbox2, (Nfft(1),Nfft(2),Nfft(3)))
      call timacc(34,2)
    endif

    !==============================================================================
    ! Compute direct term matrix elements: <ck|e^{i(G-G0).r}|cpkp>, etc.
    !==============================================================================
    if ( ik /= -1 ) then

      ! If g0(:) is non-zero, the size of mccp,mvvp must increase: umklapp vector
      call timacc(65,1)
      call logit('         mtxel_kernel: direct term mats')

      call timacc(33,1)
      if (.not. ivsaved .and. g_sum_algo == CPU_ALGO) then
        if (ii .ne. 1) then
          SAFE_DEALLOCATE(tempb_old)
          SAFE_DEALLOCATE(tempw_old)
          SAFE_DEALLOCATE(temph_old)
        end if
      endif

      if (.not. icsaved .and. xct%icpar .eq. 1) then
        if (ii .ne. 1) then
          SAFE_DEALLOCATE(mccp_old)
        end if
        SAFE_ALLOCATE(mccp_old, (xct%ng,1,1,xct%nspin))
      endif
      call timacc(33,2)

      ! Compute matrix elements: <ck|exp(i(k-kp-G0+G).r)|ckp> -> mccp
      call timacc(30,1)
      wfn_l = (/2,1/) ! <c,k|
      wfn_r = (/2,2/) ! |c,kp>
      call get_charge_matrix(mats, wfn_l, wfn_r, mccp, .not.icsaved, &
        m12_old=mccp_old, should_save=xct%icpar==1)
      call timacc(30,2)
      ! Compute matrix element: <vk|exp(i(k-kp-G0+G).r)|vkp> -> mvvp
      call timacc(31,1)
      wfn_l = (/1,1/) ! <v,k|
      wfn_r = (/1,2/) ! |v,kp>
      call get_charge_matrix(mats, wfn_l, wfn_r, mvvp, .not.ivsaved)
      call timacc(31,2)

      call timacc(65,2)
    endif ! ik /= -1

    !==============================================================================
    ! Compute screened interaction and perform sum over G (head, wings and body)
    !==============================================================================
    ! FHJ: The direct term is proportional to:
    ! bsed(iv,ivp,ic,icp) = sum(ig,igp) { [Mvvp(igp)]^* * W(ig,igp) * Mccp(ig) }
    !
    ! The strategy is:
    ! (1) Loop over all (distributed) columns of epsilon (igp):
    ! (1.1) Construct W(:, igp) for all G.
    ! (1.2) Call subroutine w_sum to multiply [Mvvp]^* by W. Eventually, we`ll get:
    !     temp?(ig,iv,ivp,is) = \sum{igp} [Mvvp(igp)]^* * W(ig,igp)
    ! (2) After the loop over igp, call g_sum to perform the sum over G and get:
    !     bsed(it,s,sp) = \sum_{ig} temp?(ig,iv,ivp,is) * Mccp(ig)
    !
    ! Things look a little bit messy because:
    ! (a) Epsilon is distributed across all processors
    ! (b) We have to worry about head/wing/body parts separately
    if ( ik /= -1 ) then

      !--------------------------------------------------------------------------
      ! Allocate buffers
      if ( .not. ivsaved ) then
        call timacc(61,1)
        call logit('         mtxel_kernel: head-wings-body')
        fact = 16.0d0 * PI_D / crys%celvol
        call logit('         mtxel_kernel: computing W(g,gp)')
        call timacc(61,2)

        if (w_sum_algo /= CPU_ALGO) then
          if (allocated(temph_accel)) then
            call accel_exit_data_map_delete(temph_accel, w_sum_algo)
            deallocate(temph_accel)
          end if
          if (allocated(tempw_accel)) then
            call accel_exit_data_map_delete(tempw_accel, w_sum_algo)
            deallocate(tempw_accel)
          end if
          if (allocated(tempb_accel)) then
            call accel_exit_data_map_delete(tempb_accel, w_sum_algo)
            deallocate(tempb_accel)
          end if
          if (xct%ivpar .eq. 1) then
            allocate(temph_accel(1, 1, xct%nspin))
            allocate(tempw_accel(xct%ng, 1, 1, xct%nspin))
            allocate(tempb_accel(xct%ng, 1, 1, xct%nspin))
          else
            allocate(temph_accel(xct%n1b_co, xct%n1b_co, xct%nspin))
            allocate(tempw_accel(xct%ng, xct%n1b_co, xct%n1b_co, xct%nspin))
            allocate(tempb_accel(xct%ng, xct%n1b_co, xct%n1b_co, xct%nspin))
          end if
          ! Transferring up zeroes is *much* faster than manually setting zeroes
          temph_accel = 0.0d0
          tempw_accel = 0.0d0
          tempb_accel = 0.0d0
          call accel_enter_data_map_to(temph_accel, w_sum_algo)
          call accel_enter_data_map_to(tempw_accel, w_sum_algo)
          call accel_enter_data_map_to(tempb_accel, w_sum_algo)
          ! Dummy arrays
          SAFE_ALLOCATE(temph, (1,1,1))
          SAFE_ALLOCATE(tempw, (1,1,1,1))
          SAFE_ALLOCATE(tempb, (1,1,1,1))
          if(.not. allocated(tempw_old)) allocate(tempw_old (1,1,1,1))
        else
          if (xct%ivpar .eq. 1) then
            SAFE_ALLOCATE(temph, (1,1,xct%nspin))
            SAFE_ALLOCATE(tempw, (xct%ng,1,1,xct%nspin))
            SAFE_ALLOCATE(tempb, (xct%ng,1,1,xct%nspin))
          else
            SAFE_ALLOCATE(temph, (xct%n1b_co,xct%n1b_co,xct%nspin))
            SAFE_ALLOCATE(tempw, (xct%ng,xct%n1b_co,xct%n1b_co,xct%nspin))
            SAFE_ALLOCATE(tempb, (xct%ng,xct%n1b_co,xct%n1b_co,xct%nspin))
          endif
          temph(:,:,:) = 0.0d0
          tempw(:,:,:,:) = 0.0d0
          tempb(:,:,:,:) = 0.0d0
        end if
      endif ! ivsaved
    endif ! ik /= -1

    ! FHJ: we want to calculate \sum_igp Mvvp^*(igp) eps_ig,igp for all ig, but eps
    ! is distributed along all processors! So instead of looping over all global
    ! igp, we loop over all processors, broadcast its chunk of eps (epscolt(:,:))
    ! and loop over all its local eps columns (igp_loc=1..ngpown_ipe)
    SAFE_ALLOCATE(epscolt, (xct%nmtxmax,xct%ngpown_max))
    SAFE_ALLOCATE(epscolt2, (xct%nmtxmax,xct%ngpown_max))
    epscolt=0
    epscolt2=0
    itotj=0
    do ipe = 0, peinf%npes-1
      ! FHJ: Number of columns of epsinv that ipe owns
      ngpown_ipe = NUMROC(xct%nmtxmax, xct%nb, ipe, 0, peinf%npes)

      !----------------------------------------------------------------------
      ! Read and distribute epsilon. All processors must work here.
      call timacc(74,1)
      if ( xct%bLowComm ) then
        ! For LowComm, each PE has all epsilon, we we have to present it`s
        ! distributed so that we can still use the loop over ipe.
        epscolt = 0D0
        irqt = irq
        if (irqt .ne. 0) then
          do igp_loc = 1, ngpown_ipe
            igp = INDXL2G(igp_loc, xct%nb, ipe, 0, peinf%npes)
            epscolt(:, igp_loc) = xct%epscol(:, igp, irqt)
          enddo
        endif
      else ! .not. bLowComm?
        do irqt = 1, qg%nr
          if (ipe==peinf%inode) then
            epscolt2(:,:)=xct%epscol(:,:,irqt)
          endif
#ifdef MPI
          call MPI_BCAST(epscolt2,xct%nmtxmax*xct%ngpown_max,MPI_SCALAR,ipe,MPI_COMM_WORLD,mpierr)
#endif
          if (irq .eq. irqt) then
            epscolt=epscolt2
          endif
        enddo
      endif
      call timacc(74,2)

      ! FHJ: Epsilon was distributed and we have no work to do. Just cycle.
      if (ik .eq. -1 .or. ivsaved) cycle

      do igp_loc = 1, ngpown_ipe

        call timacc(71,1)
        itotj = itotj + 1
        igp = INDXL2G(igp_loc, xct%nb, ipe, 0, peinf%npes)

        epscol=0d0
        if (igp .le. nepsmin) then
          epscol(1:neps)=epscolt(1:neps,igp_loc)
        else
          call timacc(71,2)
          cycle
        endif
        call timacc(71,2)
        if (phinv(igp)==0) cycle

        !------------------------------------------------------------------------
        ! Compute W(g,gp) for q + g0 = kp - k

        ! JRD: What we actually interpolate is only the head matrix elements (i.e. excluding
        ! the 1/q^2 factor), the wing part that does not diverge (i.e. multiplied by |q|), and the
        ! complete body element.  If we are truncating, instead of multiplying the wings by |q|
        ! we divide them noting that in this case the wing is DEFINED as being
        ! some smooth function multiplied by |q|*Vtrunc(G=0,q)*eps^-1(G=0,G`=0,q)
        ! This is because the wings of the dielectric matrix are /propto |q| and pick up a factor
        ! of eps^-1(G=0,G`=0,q) during the inversion process (see Epsilon/epsinv.f90).
        ! We include this factor here and not above because eps^-1(G=0,G`=0,q) varies quicker in truncated case.
        ! The Vtrunc(G=0,q) factor comes from the bare (truncated) coulomb interaction.

        call timacc(63,1)
        if (indinv(igp) .eq. 1) then ! Gp = 0. Treat head and right wing here.
          ! HEAD for q=0 and q<>0. Interpolate different things depending on
          ! what screening/truncation we use. See CPC paper Table 2.
          select case (xct%iscreen)
            case (SCREEN_SEMICOND)
              wptcol(1) = 1.0d0
            case (SCREEN_GRAPHENE)
              ! Graphene Screening
              if (xct%icutv==TRUNC_NONE) then
                if ( qq .gt. Tol_Zero) then
                  wptcol(1)=vcoul(1)*epscol(ind(1))*qq
                else
                  wptcol(1)=vcoul(1)*epscol(ind(1))*q0len
                endif
              else
                wptcol(1)=1.0d0
              endif
            case (SCREEN_METAL)
              ! Metal Screening
              wptcol(1)=vcoul(1)*epscol(ind(1))
          endselect
          ! FIRST WING (g<>0, gp=0).
          call calc_wings()
        else ! Gp /= 0 follows.
          ! SECOND WING (g=0, gp<>0).
          call calc_wings()
          ! BODY (g<>0, gp<>0). Does not depend on truncation or screening.
          do ig=2,ng_eps
            wptcol(ig)=ph(ig)*MYCONJG(phinv(igp)) * &
              epscol(ind(ig))*vcoul(indinv(igp))
          enddo
        endif ! indinv(gp)

        ! Take complex conjugate of epsinv if ik>ikp
        if (sinv.eq.-1) wptcol = MYCONJG(wptcol)
        call timacc(63,2)
        ! Perform partial summation over ig
        if (w_sum_algo == CPU_ALGO) then
          call calc_direct_w_sum()
        else
          wptcol_mat(:,igp_loc) = wptcol(:)
        end if
      enddo !igp_loc (gp loop)

      if (w_sum_algo /= CPU_ALGO) then
        call calc_direct_w_sum()
      end if
    enddo !ipe
    SAFE_DEALLOCATE(epscolt)
    SAFE_DEALLOCATE(epscolt2)
    SAFE_DEALLOCATE(isrtqi)

    if ( .not. ivsaved) then

      if (ik .eq. -1) then
        ! We don`t have to work, and no one else requires our epsilon -> return
        POP_SUB(mtxel_kernel)
        return
      endif
      SAFE_DEALLOCATE(wptcol)
      SAFE_DEALLOCATE(epscol)
      if (w_sum_algo /= CPU_ALGO) then
        SAFE_DEALLOCATE(wptcol_mat)
      end if

      ! Add in contribution from unscreened coulomb at high G vectors
      if (w_sum_algo /= CPU_ALGO) then
        call accel_update_from(tempb_accel, w_sum_algo)
        call calc_direct_unscreened_contrib(tempb_accel)
        call accel_update_to(tempb_accel, w_sum_algo)
      else
        call calc_direct_unscreened_contrib(tempb)
      end if

      SAFE_DEALLOCATE(vcoul)
      SAFE_DEALLOCATE(ind)
      SAFE_DEALLOCATE(indinv)
      SAFE_DEALLOCATE(ph)
      call free_charge_matrix(mats, xct, mvvp)

      ! When we`re using an accelerated version, the "old" values are always
      ! stored on the GPU as module variables.  So we don`t need these
      ! variables.
      if (w_sum_algo == CPU_ALGO) then
        if (xct%ivpar .eq. 1) then
          SAFE_ALLOCATE(temph_old, (1,1,xct%nspin))
          SAFE_ALLOCATE(tempw_old, (xct%ng,1,1,xct%nspin))
          SAFE_ALLOCATE(tempb_old, (xct%ng,1,1,xct%nspin))
        else
          SAFE_ALLOCATE(temph_old, (xct%n1b_co,xct%n1b_co,xct%nspin))
          SAFE_ALLOCATE(tempw_old, (xct%ng,xct%n1b_co,xct%n1b_co,xct%nspin))
          SAFE_ALLOCATE(tempb_old, (xct%ng,xct%n1b_co,xct%n1b_co,xct%nspin))
        endif

        temph_old = temph
        tempw_old = tempw
        tempb_old = tempb
      end if
    endif ! ivsaved

    !==============================================================================
    ! SUM G vectors to get DIRECT kernel term.
    !==============================================================================
    call calc_direct_g_sum()
    call timacc(67,1)

    !==============================================================================
    ! Compute exchange matrix elements: <vk|e^{iG.r}|ck>, etc.
    !==============================================================================
    call logit('         mtxel_kernel: X term matrices')
    if (.not. ixsaved .and. xct%ivpar .eq. 1 .and. xct%icpar .eq. 1) then
      if (ii .ne. 1) then
        SAFE_DEALLOCATE(mvc_old)
      endif
      SAFE_ALLOCATE(mvc_old, (xct%ng,1,1,xct%nspin))
    endif
    if (.not. ixpsaved .and. xct%ivpar .eq. 1 .and. xct%icpar .eq. 1) then
      if (ii .ne. 1) then
        SAFE_DEALLOCATE(mvpcp_old)
      endif
      SAFE_ALLOCATE(mvpcp_old, (xct%ng,1,1,xct%nspin))
    endif

    call timacc(35,1)
    ! Compute matrix elements: <vk|e^{i*G.r}|ck> -> mvc
    wfn_l = (/1,1/) ! <v,k|
    wfn_r = (/2,1/) ! |c,k>
    call get_charge_matrix(mats, wfn_l, wfn_r, mvc, &
      .not.ixsaved.or.xct%icpar==0.or.xct%ivpar==0, m12_old=mvc_old, &
      should_save=(xct%icpar==1.and.xct%ivpar==1))

    ! Compute matrix elements: <vkp|e^{i*G.r}|ckp> -> mvpcp
    wfn_l = (/1,2/) ! <v,kp|
    wfn_r = (/2,2/) ! |c,kp>
    call get_charge_matrix(mats, wfn_l, wfn_r, mvpcp, &
      .not.ixpsaved.or.xct%icpar==0.or.xct%ivpar==0, m12_old=mvpcp_old, &
      should_save=(xct%icpar==1.and.xct%ivpar==1))  
    call timacc(35,2)
    call timacc(67,2)

    !==============================================================================
    ! SUM G vectors to get EXCHANGE kernel term.
    !==============================================================================
    call logit('         mtxel_kernel: computing bsex')
    call timacc(68,1)
    ! Calc. modified Coulomb potential vbar for q=0, where vbar(G=0)=0
    SAFE_ALLOCATE(vcoul, (xct%ng))
    if (xct%qflag .ne. 1) then
      ! DYQ: If there`s a finite center-of-mass momentum, Q, the integral of the exchange 
      ! term in real-space gives delta(k-Q-k+q), so q=Q.
      ! Do not zero v(G=0) if calculating energy loss
      ! FHJ: remember that:
      ! (1) Q = -xct%finiteq; and
      ! (2) While we should compute
      ! \sum_G M_cv(G) v(G+Q) [M_cpvp(G)]^*,
      ! for historical reasons we end up evaluating this via (see gx_sum.f90)
      ! \sum_G [M_vc(G)]^* v(G-Q) M_cpvp(G)
      ! So we need to compute v(G-Q), which is at q = -Q = xct%finiteq
      ! JBH : When calling get_vcoul, also include index of exhange element,
      ! we wish to zero out.
      ! DYQ: Note, we should never zero head of exchange when using finite Q.
      call get_vcoul(.not.xct%energy_loss, .true., xct%finiteq, xct%qpg0_ind)
    elseif (xct%qflag .eq. 1) then
      call get_vcoul(.true.,.false.)
    endif
    ! Sum over G-vectors to get exchange term.
    if (.not.xct%extended_kernel) then
      call gx_sum_TDA(xct,invband,incband,vcoul,mvc,mvpcp, &
        bsex,ivp_in,icp_in,ikp,iv_in,ic_in,ik)
    else
      call calc_exchange_extended()
    endif
    call timacc(68,2)
    call logit('         mtxel_kernel: done bsex')

    call free_charge_matrices(mats)
    SAFE_DEALLOCATE(fftbox1)
    SAFE_DEALLOCATE(fftbox2)
    SAFE_DEALLOCATE_P(g_idx%g)
    SAFE_DEALLOCATE_P(g_idx%gmg0)
    if (xct%extended_kernel) then
      SAFE_DEALLOCATE_P(g_idx%mg)
      SAFE_DEALLOCATE_P(g_idx%g0mg)
    endif
    SAFE_DEALLOCATE(vcoul)
    if (save_ffts) then
      SAFE_DEALLOCATE(fftboxes)
    endif

    POP_SUB(mtxel_kernel)
    return

  contains

    !> Calculates the partial sum over ig for the direct term.
    !! Works for TDA and extended kernels.
    subroutine calc_direct_w_sum()
      PUSH_SUB(mtxel_kernel.calc_direct_w_sum)

      select case (w_sum_algo)
      case (CPU_ALGO)
         call calc_direct_w_sum_cpu()
      case (OPENACC_ALGO)
         call calc_direct_w_sum_accel()
      case (OMP_TARGET_ALGO)
         call calc_direct_w_sum_accel()
      case default
        call die("Invald algorithm for calc_direct_w_sum", &
                 only_root_writes = .true.)
      end select

      POP_SUB(mtxel_kernel.calc_direct_w_sum)

    end subroutine calc_direct_w_sum

    subroutine calc_direct_w_sum_cpu()
      integer :: t1, t1p, t1_max
      integer :: ofs1, ofs1p
      integer :: n1, n1p
      SCALAR, pointer :: m11p(:,:,:,:)
      integer :: wfn_l(2), wfn_r(2)

      PUSH_SUB(mtxel_kernel.calc_direct_w_sum_cpu)

      call timacc(64,1)
      t1_max = 1
      if (xct%extended_kernel) t1_max = 2
      ! Loop over "generalized valence WFNs"
      do t1=1,t1_max
        ofs1 = invband*(t1-1)
        n1 = nbands(t1)
        do t1p=1,t1_max
          ofs1p = invband*(t1p-1)
          n1p = nbands(t1p)
          ! This would be Mvvp within TDA
          wfn_l = (/t1,1/)  ! <t1,k|
          wfn_r = (/t1p,2/) ! |t1p,kp>
          call get_charge_matrix(mats, wfn_l, wfn_r, m11p, .true.)
          call w_sum_cpu(xct,wptcol,ofs1,ofs1p,n1,n1p,temph,tempw,tempb,m11p,indinv(igp),ng_eps)
        enddo !t1p
      enddo !t1
      call timacc(64,2)

      POP_SUB(mtxel_kernel.calc_direct_w_sum_cpu)

    end subroutine calc_direct_w_sum_cpu

    subroutine calc_direct_w_sum_accel()
      integer :: t1, t1p, t1_max
      integer :: ofs1, ofs1p
      integer :: n1, n1p
      SCALAR, pointer :: m11p(:,:,:,:)
      integer :: wfn_l(2), wfn_r(2)
      integer :: my_igp_loc, my_igp, t1_all
      SCALAR, allocatable, target :: m11p_all(:,:,:,:,:)
      integer :: s1,s2,s3,s4

      PUSH_SUB(mtxel_kernel.calc_direct_w_sum_accel)

      call timacc(64,1)

      t1_all = 1
      if (xct%extended_kernel) t1_all = 2

      SAFE_ALLOCATE( m11p_all, (xct%ng,MAXVAL(nbands),MAXVAL(nbands),xct%nspin,t1_all*t1_all) )

      t1_max = 1
      if (xct%extended_kernel) t1_max = 2
      ! Loop over "generalized valence WFNs"
      t1_all = 0
      do t1=1,t1_max
        ofs1 = invband*(t1-1)
        n1 = nbands(t1)
        do t1p=1,t1_max
          t1_all = t1_all + 1
          ofs1p = invband*(t1p-1)
          n1p = nbands(t1p)
          ! This would be Mvvp within TDA
          wfn_l = (/t1,1/)  ! <t1,k|
          wfn_r = (/t1p,2/) ! |t1p,kp>
          call get_charge_matrix(mats, wfn_l, wfn_r, m11p, .true.)
          s1 = SIZE(m11p,1)
          s2 = SIZE(m11p,2)
          s3 = SIZE(m11p,3)
          s4 = SIZE(m11p,4)
          m11p_all(1:s1,1:s2,1:s3,1:s4,t1_all) = m11p(:,:,:,:)
        enddo !t1p
      enddo !t1

      select case (w_sum_algo)
      case (OPENACC_ALGO)
         call w_sum_openacc(xct,wptcol_mat,m11p_all,temph,tempw,tempb,indinv,&
                     ng_eps, phinv, nepsmin, ngpown_ipe, invband, nbands, ipe)
      case (OMP_TARGET_ALGO)
         call w_sum_omp_target(xct,wptcol_mat,m11p_all,temph,tempw,tempb,indinv,&
                     ng_eps, phinv, nepsmin, ngpown_ipe, invband, nbands, ipe)
      case default
        call die("Invald algorithm for calc_direct_w_sum_accel", &
                 only_root_writes = .true.)
      end select

      SAFE_DEALLOCATE( m11p_all )

      call timacc(64,2)

      POP_SUB(mtxel_kernel.calc_direct_w_sum_accel)
    end subroutine calc_direct_w_sum_accel

    !> Add in contribution from unscreened coulomb at high G vectors.
    !! We can write W in a blocked structure as:
    !! [ W1     Wd ]
    !! [ Wd^H   W2 ]
    !! W1 is that part that we know (from eps*mat * vcoul), Wd ~ 0,
    !! and W2 can be approximated by a diagonal matrix:
    !! W2 ~ delta(G,Gp) * vcoul(q+Gp).
    !! This routine adds in the contribution from W2 to the body term.
    subroutine calc_direct_unscreened_contrib(tempb)
      SCALAR, target, intent(in) :: tempb(:,:,:,:)

      integer :: t1, t1p, t1_max
      integer :: ofs1, ofs1p
      integer :: n1, n1p, isv
      SCALAR, pointer :: m11p(:,:,:,:)
      SCALAR, pointer :: tempb2(:,:,:,:)
      integer :: wfn_l(2), wfn_r(2)

      PUSH_SUB(mtxel_kernel.calc_direct_unscreened_contrib)

      call timacc(66,1)
      t1_max = 1
      if (xct%extended_kernel) t1_max = 2
      ! Loop over "generalized valence WFNs"
      do t1=1,t1_max
        ofs1 = invband*(t1-1)
        n1 = nbands(t1)
        do t1p=1,t1_max
          ofs1p = invband*(t1p-1)
          n1p = nbands(t1p)
          ! This would be Mvvp within TDA
          wfn_l = (/t1,1/)  ! <t1,k|
          wfn_r = (/t1p,2/) ! |t1p,kp>
          call get_charge_matrix(mats, wfn_l, wfn_r, m11p, .true.)
          tempb2 => tempb(:, ofs1+1:ofs1+n1, ofs1p+1:ofs1p+n1p, :)
          do isv=1,xct%nspin
            do ig=1,xct%ng
              ! FHJ: integers are represented *exactly* by floating points!
              if (ph(ig)==ZERO) then
                wval = vcoul(ig)
                tempb2(ig, :, :, isv) = tempb2(ig, :, :, isv) + &
                  wval*MYCONJG(m11p(ig, :, :, isv))
              endif
            enddo
          enddo
        enddo !t1p
      enddo !t1
      call timacc(66,2)

      POP_SUB(mtxel_kernel.calc_direct_unscreened_contrib)

    end subroutine calc_direct_unscreened_contrib

    !> Calculates the partial sum over ig for the direct term.
    !! Works for TDA and extended kernels.
    subroutine calc_direct_g_sum()
      integer :: t2, t2p, t2_min
      integer :: ofs2, ofs2p
      integer :: n2, n2p
      SCALAR, pointer :: m22p(:,:,:,:)
      integer :: wfn_l(2), wfn_r(2)

      PUSH_SUB(mtxel_kernel.calc_direct_g_sum)

      call timacc(64,1)
      t2_min = 2
      if (xct%extended_kernel) t2_min = 1
      ! Loop over "generalized conduction WFNs"
      do t2=t2_min,2
        ofs2 = invband*(t2-1)
        n2 = nbands(t2)
        do t2p=t2_min,2
          ofs2p = invband*(t2p-1)
          n2p = nbands(t2p)
          ! This would be Mccp within TDA
          wfn_l = (/t2,1/)  ! <t2,k|
          wfn_r = (/t2p,2/) ! |t2p,kp>
          call get_charge_matrix(mats, wfn_l, wfn_r, m22p, .true.)
          call timacc(73,1)
          if (.not.xct%extended_kernel) then
            if (ivsaved) then
              call g_sum_TDA(xct,invband,incband,temph_old,tempb_old,m22p, &
                bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik, &
                tempw=tempw_old, bsedwing=bsedwing)
            else 
              call g_sum_TDA(xct,invband,incband,temph,tempb,m22p, &
                bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik, &
                tempw=tempw, bsedwing=bsedwing)
            endif
            call free_charge_matrix(mats, xct, mccp)
          else
            if (ivsaved) then
              call g_sum_extended(xct,ofs2,ofs2p,n2,n2p,temph_old,tempb_old,m22p, &
                bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik, &
                tempw=tempw_old, bsedwing=bsedwing)
            else 
              call g_sum_extended(xct,ofs2,ofs2p,n2,n2p,temph,tempb,m22p, &
                bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik, &
                tempw=tempw, bsedwing=bsedwing)
            endif
          endif
          call timacc(73,2)
        enddo !t2p
      enddo !t2
      if (.not.ivsaved) then
        SAFE_DEALLOCATE(temph)
        SAFE_DEALLOCATE(tempw)
        SAFE_DEALLOCATE(tempb)
      endif
      call timacc(64,2)

      POP_SUB(mtxel_kernel.calc_direct_g_sum)

    end subroutine calc_direct_g_sum

    !> Calculates all possible kernel exchange blocks:
    !! \sum_G M_21(G) v(G+Q) [M_2p1p(G)]^*
    !! = \sum_G [M_12(G)]^* v(G-Q) M_1p2p(G)
    !! Note: within TDA, we would have: 1=1p=v, 2=2p=c
    subroutine calc_exchange_extended()
      integer :: t1, t1p, t2, t2p
      integer :: ofs1, ofs1p, ofs2, ofs2p
      integer :: n1, n1p, n2, n2p
      !> This is vcoul(:) * m1p2p(:,:,:,:)
      SCALAR, pointer :: m1p2p(:,:,:,:), m12(:,:,:,:)
      SCALAR, allocatable :: v_m1p2p(:,:,:,:)
      integer :: i1p,i2p,is
      integer :: wfn_l(2), wfn_r(2)

      PUSH_SUB(mtxel_kernel.calc_exchange_extended)

      ! Loop over WFNs at kp
      do t2p=1,2
        ofs2p = invband*(t2p-1)
        n2p = nbands(t2p)
        do t1p=1,2
          ofs1p = invband*(t1p-1)
          n1p = nbands(t1p)
          ! This would be Mvpcp within TDA
          wfn_l = (/t1p,2/) ! <t1p,kp|
          wfn_r = (/t2p,2/) ! |t2p,kp>
          call get_charge_matrix(mats, wfn_l, wfn_r, m1p2p, .true.)
          SAFE_ALLOCATE(v_m1p2p, (xct%ng,n1p,n2p,xct%nspin))
          ! Construct v * m1p2p. Note that we are using the modified
          ! Coulomb potential here, i.e., vbar(G=0)=0.
          ! TODO: BLASSIFY
          do is=1,xct%nspin
            do i2p=1,n2p
              do i1p=1,n1p
                v_m1p2p(:,i1p,i2p,is) = vcoul(:) * m1p2p(:,i1p,i2p,is)
              enddo
            enddo
          enddo

          ! Loop over WFNs at k
          do t2=1,2
            ofs2 = invband*(t2-1)
            n2 = nbands(t2)
            do t1=1,2
              ofs1 = invband*(t1-1)
              n1 = nbands(t1)
              ! This would be Mvc within TDA
              wfn_l = (/t1,1/) ! <t1,k|
              wfn_r = (/t2,1/) ! |t2,k>
              call get_charge_matrix(mats, wfn_l, wfn_r, m12, .true.)
              call gx_sum_extended(xct, ofs1, ofs2, ofs1p, ofs2p, n1, n2, n1p, n2p, &
                m12, v_m1p2p, bsex, ivp_in, icp_in, ikp, iv_in, ic_in, ik)
            enddo
          enddo
          SAFE_DEALLOCATE(v_m1p2p)
        enddo
      enddo

      POP_SUB(mtxel_kernel.calc_exchange_extended)

    end subroutine calc_exchange_extended

    !> Calculate either W(0,igp) or W(ig,0).
    !! Note: we interpolate different quantities depending on the
    !! truncation and screening, and we change W to take that into account.
    subroutine calc_wings()
      integer :: ig_start, ig_end, igp_gvec, ig_gvec, ig_eps
      SCALAR :: ph2_conjg, ph_tot

      PUSH_SUB(mtxel_kernel.calc_wings)

      igp_gvec = indinv(igp)
      ph2_conjg = MYCONJG(phinv(igp))
      if (igp_gvec==1) then ! First/right wing: gp=1, g>1
        ig_start = 2
        ig_end = ng_eps
      else ! Second/left wing: gp>1, g=1
        ig_start = 1
        ig_end = 1
      endif

      do ig_gvec = ig_start, ig_end
        ph_tot = ph(ig_gvec) * ph2_conjg
        ig_eps = ind(ig_gvec)
        if (xct%iscreen==SCREEN_SEMICOND) then
          ! Semiconductor Screening
          ! --No Truncation
          wptcol(ig_gvec) = ph_tot*epscol(ig_eps)*vcoul(igp_gvec)*qq
          if (xct%icutv/=TRUNC_NONE) then
            ! --With Truncation
            wptcol(ig_gvec) = ph_tot*epscol(ig_eps)*vcoul(igp_gvec)
          endif
          if (qq<Tol_ZERO) then
            ! Semiconductor => zero q=0 wings
            wptcol(ig_gvec) = ZERO
          endif
        else
          ! Metal/Graphene Screening
          wptcol(ig_gvec) = ph_tot*epscol(ig_eps)*vcoul(igp_gvec)
          if (xct%iscreen==SCREEN_GRAPHENE .and. xct%icutv==TRUNC_WIRE .and. qq<TOL_ZERO) then
            ! Graphene Wire => zero q=0 wings
            wptcol(ig_gvec) = ZERO
          endif
        endif
      enddo

      POP_SUB(mtxel_kernel.calc_wings)

    end subroutine calc_wings

    !> Populates the array vcoul with v_G(q), where q=qq, if qq is given,
    !! or q=q0, if qq is omitted. If vbar is set to true, we zero out the
    !! G=0 component to get the modified Coulomb potential vbar.
    subroutine get_vcoul(vbar, finiteQ, qq, iqpg0)
      logical, intent(in) :: vbar
      logical, intent(in) :: finiteQ
      real(DP), intent(in), optional :: qq(3)
      integer, intent(in), optional :: iqpg0

      integer :: ikpt, ik

      PUSH_SUB(mtxel_kernel.get_vcoul)

      call timacc(62,1)
      vcoul(:) = 0.0d0
      ikpt = 0
      if (present(qq)) then
        ! if qq is given, find the index of the k-point identical to qq
        ! DYQ:  qq and k-point might differ by an Umklapp
        do ik=1,qg%nf
          vq(:) = qg%f(:,ik) - qq(:)
          qlen = DOT_PRODUCT(vq,MATMUL(crys%bdot,vq))
          if (qlen < TOL_Zero) then
            ikpt = ik
            exit
          endif
        enddo
      else
        do ik = 1,qg%nf
          ! if no qq is given, find the index of the q0 point
          vq(:) = qg%f(:,ik)
          qlen = DOT_PRODUCT(vq,MATMUL(crys%bdot,vq))
          if (qlen < TOL_Zero) then
            ikpt = ik
            exit
          endif
        enddo
      endif
      if (ikpt == 0) then
        if (finiteQ) then
          ikpt = qg%nf+1
        else
          call die("Couldn't find q-point")
        endif
      endif
      vcoul(:) = vcoularray(:,ikpt)
      if (vbar) then
        ! JBH : if finiteQ, we zero out the index of the G vector with the samllest
        ! |-Q+G|^2
        if (present(iqpg0)) then
          vcoul(iqpg0) = 0d0
        else
          vcoul(1) = 0d0
        endif
      endif
      call timacc(62,2)

      POP_SUB(mtxel_kernel.get_vcoul)

    end subroutine get_vcoul

    !> Compute all real-space WFNs to speed up the calculation of the charge
    !! density matrix later on. Call this function for each type of wave function
    !! (valence/conduction) at each k-point (k or kp).
    subroutine precompute_ffts()
      type(wavefunction), pointer :: wfn
      integer :: wtype, wprime
      integer :: is, isp, ib, offset

      PUSH_SUB(mtxel_kernel.precompute_ffts)

      do wtype = 1,2 ! loop valence/conduction wfns
        offset = (wtype-1)*invband ! 0 for valence, invband for conduction
        do wprime = 1,2 ! at k/lp
          wfn => wfns(wtype, wprime)%p
          do is = 1, xct%nspin
            do isp = 1, xct%nspinor
              do ib = 1, nbands(wtype) ! offset + ib = local band index
                call put_into_fftbox(wfn%ng, wfn%cg(1:,ib,is*isp), gvec%components, &
                  wfn%isort, fftboxes(:,:,:,is*isp,offset+ib,wprime), Nfft)
                call do_FFT(fftboxes(:,:,:,is*isp,offset+ib,wprime), Nfft, 1)
              enddo
            enddo
          enddo
        enddo
      enddo

      POP_SUB(mtxel_kernel.precompute_ffts)   

    end subroutine precompute_ffts

    !> Calculates the matrix elements M_12(G) = <1|e^(i*G.r)|2>
    !! This the complex conjugate of the charge density matrix <2|e^(-i*G.r)|1>,
    !! and we calculate it for all G vectors and for nb1 wavefunctions of type
    !! <1| and nb2 of type |2>. Result is stored in array m12.
    subroutine calc_charge_matrix(w1, w2, m12, &
      should_calculate, m12_old, should_save)
      !> each WFN is a pair (v/c, k/kp)
      integer, intent(in) :: w1(2), w2(2)
      !> The resulting matrix elements (xct%ng, nb1, nb2, xct%nspin).
      SCALAR, intent(out) :: m12(:,:,:,:)
      !> If true, calculate m12, otherwise use m12_old, if present
      logical, intent(in) :: should_calculate
      !> Save/reuse matrix elements, if appropriate/possible
      SCALAR, intent(inout), allocatable, optional :: m12_old(:,:,:,:)
      !> If true, set m12_old from m12 (calculated in the function)
      logical, intent(in), optional :: should_save

      type (wavefunction), pointer :: wfn1, wfn2
      !> Order to get G-vectors. Should be (1,2,...) unless there`s umklapp involved.
      integer, pointer :: gvecs_order(:)

      PUSH_SUB(mtxel_kernel.calc_charge_matrix)

      wfn1 => wfns(w1(1),w1(2))%p
      wfn2 => wfns(w2(1),w2(2))%p
      if (w1(2)==w2(2)) then ! <k|..|k> or <kp|..|kp> => e^{iG.r}
        gvecs_order => g_idx%g(:)
      elseif(w1(2)<w2(2)) then ! <k|..|kp> => e^{i(G-G0).r}
        gvecs_order => g_idx%gmg0(:)
      else ! <kp|..|k> => e^{i(G0-G).r}
        gvecs_order => g_idx%g0mg(:)
      endif

      if (should_calculate) then
        select case (mtxel_algo)
        case (CPU_ALGO)
          call calc_charge_matrix_core_cpu(w1, w2, wfn1, wfn2, gvecs_order, m12)
        case (OPENACC_ALGO)
          if (accel_mtxel_eps%n_ffts_per_batch .gt. 1) then
            call calc_charge_matrix_core_accel_batched(w1, w2, wfn1, wfn2, gvecs_order, m12)
          else
            call calc_charge_matrix_core_accel(w1, w2, wfn1, wfn2, gvecs_order, m12)
          end if
        case (OMP_TARGET_ALGO)
          if (accel_mtxel_eps%n_ffts_per_batch .gt. 1) then
            call calc_charge_matrix_core_accel_batched(w1, w2, wfn1, wfn2, gvecs_order, m12)
          else
            call calc_charge_matrix_core_accel(w1, w2, wfn1, wfn2, gvecs_order, m12)
          end if
        case default
          call die("Invald algorithm for mtxel", only_root_writes = .true.)
        end select

        if (present(m12_old).and.present(should_save)) then
          if (should_save) m12_old(:,:,:,:) = m12(:,:,:,:)
        endif
      else
        if (present(m12_old)) m12(:,:,:,:) = m12_old(:,:,:,:)
      endif

      POP_SUB(mtxel_kernel.calc_charge_matrix)

    end subroutine calc_charge_matrix

    subroutine calc_charge_matrix_core_cpu(w1, w2, wfn1, wfn2, gvecs_order, m12)
      !> each WFN is a pair (v/c, k/kp)
      integer, intent(in) :: w1(2), w2(2)
      type (wavefunction), intent(in) :: wfn1, wfn2
      !> Order to get G-vectors. Should be (1,2,...) unless there`s umklapp involved.
      integer, intent(in) :: gvecs_order(:)
      !> The resulting matrix elements (xct%ng, nb1, nb2, xct%nspin).
      SCALAR, intent(out) :: m12(:,:,:,:)

      integer :: ib1, ib2, is, isp
      SCALAR, dimension(:), allocatable :: tmparray

      PUSH_SUB(mtxel_kernel.calc_charge_matrix_core_cpu)

      m12(:,:,:,:) = ZERO
      SAFE_ALLOCATE(tmparray, (xct%ng))

      do is = 1, xct%nspin
        do isp = 1, xct%nspinor
          do ib1 = 1, nbands(w1(1))
            if (save_ffts) then
              fftbox1(:,:,:) = CONJG(fftboxes(:,:,:,is*isp,invband*(w1(1)-1)+ib1,w1(2)))
            else
              call put_into_fftbox(wfn1%ng,wfn1%cg(1:,ib1,is*isp),gvec%components,wfn1%isort,fftbox1,Nfft)
              call do_FFT(fftbox1,Nfft,1)
              call conjg_fftbox(fftbox1,Nfft)
            endif
            do ib2 = 1, nbands(w2(1))
              if (save_ffts) then
                ! FHJ: If we use multiply_fftboxes directly we endup moving
                ! more data than necessary here.
                fftbox2 = fftbox1*fftboxes(:,:,:,is*isp,invband*(w2(1)-1)+ib2,w2(2))
              else
                call put_into_fftbox(wfn2%ng,wfn2%cg(1:,ib2,is*isp),gvec%components,wfn2%isort,fftbox2,Nfft)
                call do_FFT(fftbox2,Nfft,1)
                call multiply_fftboxes(fftbox1,fftbox2,Nfft)
              endif
              call do_FFT(fftbox2,Nfft,1)
              call get_from_fftbox(xct%ng,tmparray,gvec%components,gvecs_order,fftbox2,Nfft,scale)
              m12(1:,ib1,ib2,is) = m12(1:,ib1,ib2,is) + tmparray
            enddo !ib2
          enddo !ib1
        enddo !isp
      enddo !is

      SAFE_DEALLOCATE(tmparray)

      POP_SUB(mtxel_kernel.calc_charge_matrix_core_cpu)

    end subroutine calc_charge_matrix_core_cpu

    subroutine calc_charge_matrix_core_accel_batched(w1, w2, wfn1, wfn2, gvecs_order, m12)
      !> each WFN is a pair (v/c, k/kp)
      integer, intent(in) :: w1(2), w2(2)
      type (wavefunction), intent(in) :: wfn1, wfn2
      !> Order to get G-vectors. Should be (1,2,...) unless there`s umklapp involved.
      integer, intent(in) :: gvecs_order(:)
      !> The resulting matrix elements (xct%ng, nb1, nb2, xct%nspin).
      SCALAR, intent(out) :: m12(:,:,:,:)

      integer :: ib1, ib2, is, isp

      integer :: ierr

      integer :: i_batch, i_in_batch, n_ffts_in_batch

      PUSH_SUB(mtxel_kernel.calc_charge_matrix_core_accel)

      m12(:,:,:,:) = ZERO

      if (save_ffts) then
        call die("save_ffts not implemented for OpenACC MTXEL", only_root_writes = .true.)
      end if

      accel_mtxel_eps%outer_sort = 0
      accel_mtxel_eps%outer_sort(1:wfn1%ng) = wfn1%isort(1:wfn1%ng)
      call accel_update_to(accel_mtxel_eps%outer_sort, mtxel_algo)
      accel_mtxel_eps%inner_sort = 0
      accel_mtxel_eps%inner_sort(1:wfn2%ng) = wfn2%isort(1:wfn2%ng)
      call accel_update_to(accel_mtxel_eps%inner_sort, mtxel_algo)
      accel_mtxel_eps%result_sort = 0
      accel_mtxel_eps%result_sort(1:xct%ng) = gvecs_order(1:xct%ng)
      call accel_update_to(accel_mtxel_eps%result_sort, mtxel_algo)
      accel_mtxel_eps%Nfft(1:3) = Nfft(1:3)
      call accel_update_to(accel_mtxel_eps%Nfft, mtxel_algo)

      do is = 1, xct%nspin
        do isp = 1, xct%nspinor
          accel_mtxel_eps%outer_bands_2d = ZERO
          accel_mtxel_eps%outer_bands_2d(1:wfn1%ng, 1:nbands(w1(1))) = &
              wfn1%cg(1:wfn1%ng, 1:nbands(w1(1)), is * isp)
          call accel_update_to(accel_mtxel_eps%outer_bands_2d, mtxel_algo)

          accel_mtxel_eps%inner_bands_2d = ZERO
          accel_mtxel_eps%inner_bands_2d(1:wfn2%ng, 1:nbands(w2(1))) = &
              wfn2%cg(1:wfn2%ng, 1:nbands(w2(1)), is * isp)
          call accel_update_to(accel_mtxel_eps%inner_bands_2d, mtxel_algo)
          do ib1 = 1, nbands(w1(1))
            call accel_zero_box(accel_mtxel_eps%outer_fftbox, &
                              accel_mtxel_eps%Nfft, mtxel_algo)
            call accel_put_into_fftbox(wfn1%ng, &
                     accel_mtxel_eps%outer_bands_2d(:, ib1), &
                     accel_mtxel_eps%gvec_comp, accel_mtxel_eps%outer_sort, &
                     accel_mtxel_eps%outer_fftbox, accel_mtxel_eps%Nfft, &
                     1.0D+00, mtxel_algo)
            call accel_run_fft(accel_mtxel_eps%outer_fftbox, &
                     accel_mtxel_eps%outer_plan, 1, mtxel_algo)

            ! Do a round-robin distribution of the conduction band into batches
            do i_batch = 1, (nbands(w2(1)) - 1) / accel_mtxel_eps%n_ffts_per_batch + 1
              n_ffts_in_batch = 0

              ! Find size of batch and create batches of initial wavefunctions
              call accel_zero_box_batched(accel_mtxel_eps%inner_fftbox_batched, &
                                          accel_mtxel_eps%Nfft, accel_mtxel_eps%n_ffts_per_batch, mtxel_algo)
              do i_in_batch = 1, accel_mtxel_eps%n_ffts_per_batch
                ib2 = (i_batch - 1) * accel_mtxel_eps%n_ffts_per_batch + i_in_batch
                if (ib2 > nbands(w2(1))) cycle

                n_ffts_in_batch = n_ffts_in_batch + 1
              end do
              ! Offset for first band in current batch
              ib2 = (i_batch - 1) * accel_mtxel_eps%n_ffts_per_batch
              call accel_put_into_fftbox_batched(wfn2%ng, &
                   accel_mtxel_eps%inner_bands_2d, accel_mtxel_eps%gvec_comp, &
                   accel_mtxel_eps%inner_sort, &
                   accel_mtxel_eps%inner_fftbox_batched, &
                   accel_mtxel_eps%Nfft, 1.0D+00, ib2, n_ffts_in_batch, &
                   mtxel_algo)

              ! Do the hard work of FFTs across batch
              call accel_run_fft_batched(accel_mtxel_eps%inner_fftbox_batched, &
                                         accel_mtxel_eps%inner_plan_batched, 1, mtxel_algo)
              call accel_box_multiply_batched(accel_mtxel_eps%inner_fftbox_batched, &
                                              accel_mtxel_eps%outer_fftbox, accel_mtxel_eps%Nfft, 1, &
                                              n_ffts_in_batch, mtxel_algo)
              call accel_run_fft_batched(accel_mtxel_eps%inner_fftbox_batched, &
                                         accel_mtxel_eps%inner_plan_batched, 1, mtxel_algo)

              ! Extract finished FFTs and update results
              ! TODO:  In Epsilon, had to use tmparray to give better performance, needed here too?
              call accel_get_from_fftbox_batched(xct%ng, accel_mtxel_eps%result_array_batched, &
                                                 accel_mtxel_eps%gvec_comp, accel_mtxel_eps%result_sort, &
                                                 accel_mtxel_eps%inner_fftbox_batched, &
                                                 accel_mtxel_eps%Nfft, scale, n_ffts_in_batch, mtxel_algo)
              call accel_update_from(accel_mtxel_eps%result_array_batched, mtxel_algo)

              ! Accumulate results on CPU, as m12 too big to fit on GPU
              do i_in_batch = 1, n_ffts_in_batch
                ib2 = (i_batch - 1) * accel_mtxel_eps%n_ffts_per_batch + i_in_batch
                m12(:, ib1, ib2, is) = m12(:, ib1, ib2, is) + accel_mtxel_eps%result_array_batched(:, i_in_batch)
              end do
            end do ! i_batch over ib2
          enddo !ib1
        enddo !isp
      enddo !is

      POP_SUB(mtxel_kernel.calc_charge_matrix_core_accel_batched)

    end subroutine calc_charge_matrix_core_accel_batched

    subroutine calc_charge_matrix_core_accel(w1, w2, wfn1, wfn2, gvecs_order, m12)
      !> each WFN is a pair (v/c, k/kp)
      integer, intent(in) :: w1(2), w2(2)
      type (wavefunction), intent(in) :: wfn1, wfn2
      !> Order to get G-vectors. Should be (1,2,...) unless there`s umklapp involved.
      integer, intent(in) :: gvecs_order(:)
      !> The resulting matrix elements (xct%ng, nb1, nb2, xct%nspin).
      SCALAR, intent(out) :: m12(:,:,:,:)

      integer :: ib1, ib2, is, isp

      integer :: ierr

      PUSH_SUB(mtxel_kernel.calc_charge_matrix_core_accel)

      m12(:,:,:,:) = ZERO

      if (save_ffts) then
        call die("save_ffts not implemented for OpenACC MTXEL", only_root_writes = .true.)
      end if

      accel_mtxel_eps%outer_sort = 0
      accel_mtxel_eps%outer_sort(1:wfn1%ng) = wfn1%isort(1:wfn1%ng)
      call accel_update_to(accel_mtxel_eps%outer_sort, mtxel_algo)
      accel_mtxel_eps%inner_sort = 0
      accel_mtxel_eps%inner_sort(1:wfn2%ng) = wfn2%isort(1:wfn2%ng)
      call accel_update_to(accel_mtxel_eps%inner_sort, mtxel_algo)
      accel_mtxel_eps%result_sort = 0
      accel_mtxel_eps%result_sort(1:xct%ng) = gvecs_order(1:xct%ng)
      call accel_update_to(accel_mtxel_eps%result_sort, mtxel_algo)
      accel_mtxel_eps%Nfft(1:3) = Nfft(1:3)
      call accel_update_to(accel_mtxel_eps%Nfft, mtxel_algo)

      do is = 1, xct%nspin
        do isp = 1, xct%nspinor
          accel_mtxel_eps%outer_bands_2d = ZERO
          accel_mtxel_eps%outer_bands_2d(1:wfn1%ng, 1:nbands(w1(1))) = &
              wfn1%cg(1:wfn1%ng, 1:nbands(w1(1)), is * isp)
          call accel_update_to(accel_mtxel_eps%outer_bands_2d, mtxel_algo)
          accel_mtxel_eps%inner_bands_2d = ZERO
          accel_mtxel_eps%inner_bands_2d(1:wfn2%ng, 1:nbands(w2(1))) = &
              wfn2%cg(1:wfn2%ng, 1:nbands(w2(1)), is * isp)
          call accel_update_to(accel_mtxel_eps%inner_bands_2d, mtxel_algo)
          do ib1 = 1, nbands(w1(1))
            call accel_zero_box(accel_mtxel_eps%outer_fftbox, &
                              accel_mtxel_eps%Nfft, mtxel_algo)
            call accel_put_into_fftbox(wfn1%ng, &
                     accel_mtxel_eps%outer_bands_2d(:, ib1), &
                     accel_mtxel_eps%gvec_comp, accel_mtxel_eps%outer_sort, &
                     accel_mtxel_eps%outer_fftbox, accel_mtxel_eps%Nfft, &
                     1.0D+00, mtxel_algo)
            call accel_run_fft(accel_mtxel_eps%outer_fftbox, &
                     accel_mtxel_eps%outer_plan, 1, mtxel_algo)
            do ib2 = 1, nbands(w2(1))
              call accel_zero_box(accel_mtxel_eps%inner_fftbox, &
                                accel_mtxel_eps%Nfft, mtxel_algo)
              call accel_put_into_fftbox(wfn2%ng, &
                     accel_mtxel_eps%inner_bands_2d(:, ib2), &
                     accel_mtxel_eps%gvec_comp, accel_mtxel_eps%inner_sort, &
                     accel_mtxel_eps%inner_fftbox, accel_mtxel_eps%Nfft, &
                     1.0D+00, mtxel_algo)
              call accel_run_fft(accel_mtxel_eps%inner_fftbox, &
                       accel_mtxel_eps%outer_plan, 1, mtxel_algo)
              call accel_box_multiply(accel_mtxel_eps%inner_fftbox, &
                       accel_mtxel_eps%outer_fftbox, accel_mtxel_eps%Nfft, 1, &
                       mtxel_algo)
              call accel_run_fft(accel_mtxel_eps%inner_fftbox, &
                       accel_mtxel_eps%outer_plan, 1, mtxel_algo)
              call accel_get_from_fftbox(xct%ng, &
                       accel_mtxel_eps%result_array_2d(:, ib2), &
                       accel_mtxel_eps%gvec_comp, accel_mtxel_eps%result_sort, &
                       accel_mtxel_eps%inner_fftbox, accel_mtxel_eps%Nfft, &
                       scale, mtxel_algo)
              call accel_update_from(accel_mtxel_eps%result_array_2d(:, ib2), &
                       mtxel_algo)
            end do
            do ib2 = 1, nbands(w2(1))
              m12(:, ib1, ib2, is) = m12(:, ib1, ib2, is) + accel_mtxel_eps%result_array_2d(:, ib2)
            enddo !ib2
          enddo !ib1
        enddo !isp
      enddo !is

      POP_SUB(mtxel_kernel.calc_charge_matrix_core_accel)

    end subroutine calc_charge_matrix_core_accel

    !> Return the matrix element M_12(G) = <1|e^(i*G.r)|2>
    !! The matrix element will be either calculated directly, or obtained
    !! from complex conjugation of a known matrix element.
    !! Parameter m12 will point to the resulting matrix element, so do
    !! not free m12 directly, use free_charge_matri{x|ces} instead.
    subroutine get_charge_matrix(mats, w1, w2, m12, should_calculate, m12_old, should_save)
      type(matrix_4), intent(inout), target :: mats(2,2,2,2)
      integer, intent(in) :: w1(2), w2(2)
      SCALAR, pointer :: m12(:,:,:,:) !< intent not permitted for pointer
      !> If true, calculate m12, otherwise use m12_old, if present
      logical, intent(in) :: should_calculate
      !> Save/reuse matrix elements, if appropriate/possible
      SCALAR, intent(inout), allocatable, optional :: m12_old(:,:,:,:)
      !> If true, set m12_old from m12 (calculated in the function)
      logical, intent(in), optional :: should_save

      integer :: ig, nb1, nb2, i1, i2
      logical :: is_assoc

      PUSH_SUB(mtxel_kernel.get_charge_matrix)

      if (associated(mats(w1(1),w1(2),w2(1),w2(2))%p)) then
        m12 => mats(w1(1),w1(2),w2(1),w2(2))%p
      else
        nb1 = nbands(w1(1))
        nb2 = nbands(w2(1))

        ! store status in variable, in case some indices are equal and the allocation changes the status
        is_assoc = associated(mats(w2(1),w2(2),w1(1),w1(2))%p)
        SAFE_ALLOCATE(mats(w1(1),w1(2),w2(1),w2(2))%p, (xct%ng,nb1,nb2,xct%nspin))
        m12 => mats(w1(1),w1(2),w2(1),w2(2))%p

        if (is_assoc) then
          m12 = ZERO
          do i2 = 1, nb2
            do i1 = 1, nb1
              do ig=1,xct%ng
                m12(ig,i1,i2,:) = MYCONJG(mats(w2(1),w2(2),w1(1),w1(2))%p(g_idx%mg(ig),i2,i1,:))
              enddo
            enddo
          enddo
        else
          call calc_charge_matrix(w1, w2, m12, should_calculate, m12_old, should_save)
        endif
      endif

      POP_SUB(mtxel_kernel.get_charge_matrix)

    end subroutine get_charge_matrix

  end subroutine mtxel_kernel

  !> If there`s umklapp involved, e.g. <nk|e^{i(G-G0).r}|mkp>, we
  !! perform the FFTs as usual, but map each G vector to G-G0, etc.
  !! This subroutine initializes all these mappings in g_idx.
  subroutine init_g_idx(xct, gvec, g0, sinv, g_idx)
    type (xctinfo), intent(in) :: xct
    type(gspace), intent(in) :: gvec
    integer, intent(in) :: g0(3) !< umklapp vector
    integer, intent(in) :: sinv !< -1 if ik > ikp, 1 o.w.
    type(gvec_indices), intent(out) :: g_idx

    integer :: ig, gg(3)

    PUSH_SUB(init_g_idx)

    call timacc(32,1)
    ! Index of G => identity
    SAFE_ALLOCATE(g_idx%g, (gvec%ng))
    do ig=1,gvec%ng
      g_idx%g(ig)=ig
    enddo
    ! Index of G-G0
    SAFE_ALLOCATE(g_idx%gmg0, (gvec%ng))
    do ig=1,xct%ng
      gg(:) = sinv * (gvec%components(:,ig) - g0(:))
      call findvector(g_idx%gmg0(ig),gg,gvec)
      if (g_idx%gmg0(ig) == 0) call die('cannot find G-G0')
    enddo
    if (xct%extended_kernel) then
      ! Index of -G
      SAFE_ALLOCATE(g_idx%mg, (gvec%ng))
      do ig=1,xct%ng
        gg(:) = - gvec%components(:,ig)
        call findvector(g_idx%mg(ig),gg,gvec)
        if (g_idx%mg(ig) == 0) call die('cannot find -G')
        if (g_idx%mg(ig) > xct%ng) g_idx%mg(ig) = xct%ng
      enddo
      ! Index of G0-G
      ! TODO: remove meh, we never have kp,k
      SAFE_ALLOCATE(g_idx%g0mg, (gvec%ng))
      do ig=1,xct%ng
        gg(:) = sinv * (g0(:) - gvec%components(:,ig))
        call findvector(g_idx%g0mg(ig),gg,gvec)
        if (g_idx%g0mg(ig) == 0) call die('cannot find G0-G')
      enddo
    endif
    call timacc(32,2)

    POP_SUB(init_g_idx)

  end subroutine init_g_idx

  !> Frees and nullifies charge density pointed by m12.
  subroutine free_charge_matrix(mats, xct, m12)
    type(matrix_4), intent(inout) :: mats(2,2,2,2)
    type(xctinfo), intent(in) :: xct
    SCALAR, pointer, intent(inout) :: m12(:,:,:,:)

    integer :: t1, p1, t2, p2

    PUSH_SUB(free_charge_matrix)

    if (xct%extended_kernel) then
      POP_SUB(free_charge_matrix)
      return
    endif

    outer_loop: do t1=1,2
      do p1=1,2
        do t2=1,2
          do p2=1,2
            if (associated(m12, mats(t1,p1,t2,p2)%p)) then
              SAFE_DEALLOCATE_P(mats(t1,p1,t2,p2)%p)
              nullify(mats(t1,p1,t2,p2)%p)
              nullify(m12)
              exit outer_loop
            endif
            SAFE_DEALLOCATE_P(mats(t1,p1,t2,p2)%p)
            nullify(mats(t1,p1,t2,p2)%p)
          enddo
        enddo
      enddo
    enddo outer_loop
    
    POP_SUB(free_charge_matrix)

  end subroutine free_charge_matrix

  !> Deallocates all charge density matrices.
  subroutine free_charge_matrices(mats)
    type(matrix_4), intent(inout) :: mats(2,2,2,2)
    
    integer :: t1, p1, t2, p2
    
    PUSH_SUB(free_charge_matrices)
    
    do t1=1,2
      do p1=1,2
        do t2=1,2
          do p2=1,2
            SAFE_DEALLOCATE_P(mats(t1,p1,t2,p2)%p)
            nullify(mats(t1,p1,t2,p2)%p)
          enddo
        enddo
      enddo
    enddo
    
    POP_SUB(free_charge_matrices)
    
  end subroutine free_charge_matrices

  !> Warn the user that the g-space is too small and that mapping is no good.
  subroutine warn_about_mapping()
    logical, save :: warned=.false.
    
    PUSH_SUB(warn_about_mapping)
    
    if (peinf%inode==0 .and. .not.warned) then
      write(0,*)
      write(0,'(a)') 'WARNING: at least one vector from the epsilon G-space could not be mapped to '
      write(0,'(a)') ' the WFN G-space. This means that either:'
      write(0,'(a)') ' (1) The WFN cutoff is too small (most likely and dangerous); or'
      write(0,'(a)') ' (2) The cutoff of the dielectric matrix is simply huge.'
      write(0,'(a)') ' Consider using the gsphere.py utility to figure out the cause of this warning.'
      write(0,*)
      warned=.true.
    endif
    
    POP_SUB(warn_about_mapping)
    
  end subroutine warn_about_mapping
  
end module mtxel_kernel_m
