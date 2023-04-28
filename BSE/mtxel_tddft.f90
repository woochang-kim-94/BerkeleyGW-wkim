!=================================================================================
!
! Routines:
!
! (1) mtxel_tddft()         Originally By MJ     Last Modified  10/14/2013 (MJ)
!
!     input: crys, gvec, syms, qg, wfnc, wfncp, wfnv, wfnvp,
!            xct types
!            ii  = label of ck block
!            ipe = label of the PE that contains the ckp block
!     output: bsedbody,bsedhead,bsex,bset = kernel matrix elements
!              between ck (all v) and ckp (all vp) blocks
!
! FIXME: Please refactor against mtxel_kernel.f90! They are ~50% identical!
!
!     Calculate the head, body, exchange of the kernel (see
!     eq. 34, 35, and 41-46, Rohlfing & Louie). Set W(G, Gp) to
!     V(G) \delta(G,Gp) in those expressions. The exchange has just
!     the proper part, no divergent contribution.
!     Also calculate the fxc matrix elements
!
!=================================================================================

#include "f_defs.h"

module mtxel_tddft_m

  use global_m
  use fftw_m
  use misc_m
  use g_sum_m
  use gx_sum_m
  implicit none

  public :: mtxel_tddft

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

  subroutine mtxel_tddft(crys,gvec,syms,qg,wfnc,wfncp,wfnvp, &
    wfnv,xct,leading_dim,bsedbody,bsedhead,bsex,bset,ii,ik,ikp, &
    ic_in,icp_in,iv_in,ivp_in,vcoularray_mod,vcoularray,fq,qq,g0, &
    ifq,irq,q0len)

    type (crystal), intent(in) :: crys
    type (gspace), intent(in) :: gvec
    type (symmetry), intent(in) :: syms
    type (grid), intent(in) :: qg
    type (wavefunction), target, intent(in) :: wfnc,wfncp,wfnvp,wfnv
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: leading_dim
    SCALAR, intent(inout) :: bsedbody(:,:,:), bsedhead(:,:,:),  &
      bsex(:,:,:), bset(:,:,:) !< (leading_dim,xct%nspin,xct%nspin)
    integer, intent(in) :: ii,ik,ikp,ic_in,icp_in,iv_in,ivp_in
    real(DP), intent(in) :: vcoularray_mod(:,:) !< (xct%ng,qg%nf)
    real(DP), intent(in) :: vcoularray(:,:) !< (xct%ng,qg%nf)
    real(DP), intent(in) :: fq(3),qq
    integer, intent(in) :: g0(3)
    integer, intent(in) :: ifq !< q-pt index in full BZ
    integer, intent(in) :: irq !< k-pt index in reduced BZ
    real(DP), intent(in) :: q0len

    integer :: ipe
    character :: filename*20
    integer :: jj,sinv
    integer :: ig,igp,igpe
    !> Number of conduction/valence bands we are dealing with. nbands is an
    !! array with (/ invband, incband /)
    integer :: incband, invband, nbands(2)
    !> The "left" and "right" set of WFNs that we use to get charge density matrices.
    !! The first index represents a v=1 or c=2 type of WFN, and the second
    !! index is whether it is at k(=1) or kp(=2).
    integer :: wfn_l(2), wfn_r(2)
    real(DP) :: qlen,vq(3)

    integer :: ipeDumb
    integer, save :: ik_old = 0
    integer, save :: ikp_old = 0
    integer, save :: iv_old = 0
    integer, save :: ivp_old = 0
    integer, save :: ic_old = 0
    integer, save :: icp_old = 0
    logical :: ivsaved, icsaved, ixsaved, ixpsaved

    real(DP), allocatable :: vcoul(:)
    SCALAR, allocatable :: &
      ph(:), &
      temph(:,:,:), tempb(:,:,:,:)
    SCALAR, pointer :: mccp(:,:,:,:),mvvp(:,:,:,:),mvc(:,:,:,:),mvpcp(:,:,:,:)
    SCALAR, save, allocatable :: &
      temph_old(:,:,:), tempb_old(:,:,:,:)
    SCALAR, save, allocatable :: &
      mccp_old(:,:,:,:),mvc_old(:,:,:,:),mvpcp_old(:,:,:,:)

    integer :: irqt
    real(DP) :: fact

    complex(DPC), dimension(:,:,:), allocatable :: fftbox1,fftbox2,fxc
    complex(DPC), dimension(:,:,:), allocatable :: fftboxv,fftboxc,fftboxvp,fftboxcp
    complex(DPC) :: sumt
    integer :: iscp, isc, iv, ivp, ic, icp, iit, ix, iy, iz
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

    PUSH_SUB(mtxel_tddft)

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
      call init_g_idx(xct, gvec, g0, sinv, g_idx)
      call timacc(61,1)
    endif

    !==============================================================================
    ! Initialize v(q)
    !==============================================================================

    if ( ik /= -1 .and. (.not. ivsaved)) then

      ! Compute Coulomb interaction at this q vector: vcoul(q+G) for q+G=0 we set
      ! the interaction to V(q0) where q0 is the small vector used for epsilon
      SAFE_ALLOCATE(vcoul, (xct%ng))
      call get_vcoul(.false., fq)

      ! MJ: There are only one G-space here
      ! - The gvec space: vectors ig_gvec sorted wrt |G|^2;
      ! Different vectors/matrices will be ordered wrt this space:
      ! - gvec space: M12 matrices, vcoul,
      !               temph, tempb
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
        SAFE_ALLOCATE(fftboxes, (Nfft(1),Nfft(2),Nfft(3),xct%nspin,invband+incband,2))
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
      call logit('         mtxel_tddft: direct term mats')

      call timacc(33,1)
      if (.not. ivsaved) then
        if (ii .ne. 1) then
          SAFE_DEALLOCATE(temph_old)
          SAFE_DEALLOCATE(tempb_old)
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
    ! Compute bare interaction and perform sum over G (body)
    !==============================================================================
    ! MJ: The direct term (in this case coming from Fock) is proportional to:
    ! bsed(iv,ivp,ic,icp) = sum(ig,igp) { [Mvvp(igp)]^* * V(ig,igp) * Mccp(ig) }
    !                     = sum(ig) {[Mvvp(ig)]^* * V(ig) * Mccp(ig) }
    !
    ! (1.2) Call subroutine w_sum to multiply [Mvvp]^* by V. Eventually, we`ll get:
    !     temp?(ig,iv,ivp,is) = [Mvvp(ig)]^* * V(ig)
    ! (2) After the loop over igp, call g_sum to perform the sum over G and get:
    !     bsed(it,s,sp) = \sum_{ig} temp?(ig,iv,ivp,is) * Mccp(ig)
    !
    if ( ik /= -1 ) then

      !--------------------------------------------------------------------------
      ! Allocate buffers
      if ( .not. ivsaved ) then
        call timacc(61,1)
        call logit('         mtxel_tddft: head-body')
        call logit('         mtxel_tddft: computing V(g,gp)')
        call timacc(61,2)
        if (xct%ivpar .eq. 1) then
          SAFE_ALLOCATE(temph, (1,1,xct%nspin))
          SAFE_ALLOCATE(tempb, (xct%ng,1,1,xct%nspin))
        else
          SAFE_ALLOCATE(temph, (xct%n1b_co,xct%n1b_co,xct%nspin))
          SAFE_ALLOCATE(tempb, (xct%ng,xct%n1b_co,xct%n1b_co,xct%nspin))
        endif
        temph(:,:,:) = 0.0d0
        tempb(:,:,:,:) = 0.0d0

      endif ! ivsaved
    endif ! ik /= -1

    if ( .not. ivsaved) then

      if (ik .eq. -1) then
        ! We don`t have to work, and no one else requires our epsilon -> return
        POP_SUB(mtxel_tddft)
        return
      endif
      ! The first part of the sum -- this is not really a sum here... 
      call calc_direct_unscreened_contrib()
      SAFE_DEALLOCATE(vcoul)
      call free_charge_matrix(mats, xct, mvvp)

      if (xct%ivpar .eq. 1) then
        SAFE_ALLOCATE(tempb_old, (xct%ng,1,1,xct%nspin))
        SAFE_ALLOCATE(temph_old, (1,1,xct%nspin))
      else
        SAFE_ALLOCATE(tempb_old, (xct%ng,xct%n1b_co,xct%n1b_co,xct%nspin))
        SAFE_ALLOCATE(temph_old, (xct%n1b_co,xct%n1b_co,xct%nspin))
      endif
      temph_old = temph
      tempb_old = tempb
    endif ! ivsaved

    !==============================================================================
    ! SUM G vectors to get DIRECT kernel term.
    !==============================================================================
    call calc_direct_g_sum()
    call timacc(67,1)

    !==============================================================================
    ! Compute exchange matrix elements: <vk|e^{iG.r}|ck>, etc.
    !==============================================================================
    call logit('         mtxel_tddft: X term matrices')
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
    call logit('         mtxel_tddft: computing bsex')
    call timacc(68,1)
    ! Calc. modified Coulomb potential vbar for q=0, where vbar(G=0)=0
    SAFE_ALLOCATE(vcoul, (xct%ng))
    call get_vcoul(.not.xct%energy_loss)
    ! Sum over G-vectors to get exchange term.
    if (.not.xct%extended_kernel) then
      call gx_sum_TDA(xct,invband,incband,vcoul,mvc,mvpcp, &
        bsex,ivp_in,icp_in,ikp,iv_in,ic_in,ik)
    else
      call calc_exchange_extended()
    endif
    call timacc(68,2)
    call logit('         mtxel_tddft: done bsex')

    call free_charge_matrices(mats)

    !==============================================================================
    ! Compute fxc matrix elements: <vck|fxc|v`c`k`>, etc.
    !==============================================================================
    call logit('         mtxel_tddft: computing bset')
    if (xct%coul_mod_flag) then

      SAFE_ALLOCATE(fxc, (Nfft(1), Nfft(2), Nfft(3)))
      fxc(:,:,:) = ZERO
      if (xct%coulomb_mod%screening_length < TOL_SMALL) then
        call tdlda(xct,crys,gvec,Nfft,xct%wpg%rho,fxc)
      else
        call tdrsh(xct,crys,gvec,Nfft,xct%wpg%rho,fxc)
      endif

      ! For now we will do more work by recalculating the realspace wfns..
      ! unless they have already been computed before
      ! Note: in principle we could use mcv and mcpvp and do inverse 
      ! Fourier transform as well...  this would save on two FFT`s per
      ! matrix element

      if (.not.xct%extended_kernel) then

        SAFE_ALLOCATE(fftboxv, (Nfft(1),Nfft(2),Nfft(3)))
        SAFE_ALLOCATE(fftboxc, (Nfft(1),Nfft(2),Nfft(3)))
        SAFE_ALLOCATE(fftboxvp, (Nfft(1),Nfft(2),Nfft(3)))
        SAFE_ALLOCATE(fftboxcp, (Nfft(1),Nfft(2),Nfft(3)))

        fact = scale/(4.0D0*PI_D)
        do iscp=1,xct%nspin
          do isc=1,xct%nspin
             do icp=1,incband
               call put_into_fftbox(wfns(2,2)%p%ng, wfns(2,2)%p%cg(1:,icp,iscp), gvec%components, &
                       wfns(2,2)%p%isort, fftboxcp, Nfft)
               call do_FFT(fftboxcp, Nfft, 1)
               do ivp=1,invband
                 call put_into_fftbox(wfns(1,2)%p%ng, wfns(1,2)%p%cg(1:,ivp,iscp), gvec%components, &
                       wfns(1,2)%p%isort, fftboxvp, Nfft)
                 call do_FFT(fftboxvp, Nfft, 1)
                 call conjg_fftbox(fftboxvp,Nfft)
                 do ic=1,incband
                   call put_into_fftbox(wfns(2,1)%p%ng, wfns(2,1)%p%cg(1:,ic,isc), gvec%components, &
                       wfns(2,1)%p%isort, fftboxc, Nfft)
                   call do_FFT(fftboxc, Nfft, 1)
                   call conjg_fftbox(fftboxc,Nfft)
                   do iv=1,invband
                     call put_into_fftbox(wfns(1,1)%p%ng, wfns(1,1)%p%cg(1:,iv,isc), gvec%components, &
                       wfns(1,1)%p%isort, fftboxv, Nfft)
                     call do_FFT(fftboxv, Nfft, 1)

                     ! Multiply fftboxes
                     sumt = ZERO
                     !$OMP PARALLEL PRIVATE (ix,iy,iz)
                     !$OMP DO
                     do iz = 1, Nfft(3)
                       do iy = 1, Nfft(2)
                         do ix = 1, Nfft(1)
                           sumt = sumt + fftboxv(ix,iy,iz) * fftboxc(ix,iy,iz) * fftboxvp(ix, iy, iz)  &
                                         * fftboxcp(ix,iy,iz) * fxc(ix,iy,iz)
                         enddo
                       enddo
                     enddo
                     !$OMP END DO
                     !$OMP END PARALLEL

                     !write(6,*) 'iv, ic, ivp, icp, sumt = ',iv, ic, ivp, icp, sumt*(1/scale)/(4.0d0*PI_D)*(scale)**2.
                     !write(6,*) 'iv, ic, ivp, icp, sumt = ',iv, ic, ivp, icp, sumt*fact

                     if (xct%icpar .eq. 0) then
                       iit = peinf%wown(1,1,ikp,1,1,ik) + xct%nvb_co*xct%nvb_co*xct%ncb_co*(icp-1) &
                          + xct%nvb_co*xct%nvb_co*(ic-1) + xct%nvb_co*(ivp-1) + iv -1
                     else if (xct%ivpar .eq. 0) then
                       iit = peinf%wown(1,icp_in,ikp,1,ic_in,ik) + xct%nvb_co*(ivp-1) + iv -1
                     else
                       iit = peinf%wown(ivp_in,icp_in,ikp,iv_in,ic_in,ik)
                     endif
                     bset(iit,isc,iscp) = bset(iit,isc,iscp) + sumt*fact
                   enddo !iv
                 enddo !ic
               enddo !ivp
             enddo !icp
           enddo !isc
         enddo !iscp
        SAFE_DEALLOCATE(fftboxv)
        SAFE_DEALLOCATE(fftboxc)
        SAFE_DEALLOCATE(fftboxvp)
        SAFE_DEALLOCATE(fftboxcp)

!      else ! extended kernel
!
!         do iscp=1,xct%nspin
!           do isc=1,xct%nspin
!             do i2p=1,n2p
!               gi2p = i2p + ofs2p
!               do i1p=1,n1p
!                 gi1p = i1p + ofs1p
!                 do i2=1,n2
!                   gi2 = i2 + ofs2
!                   do i1=1,n1
!                     gi1 = i1 + ofs1
!                     if (xct%icpar .eq. 0) then
!                       iit = peinf%wown(1,1,ikp,1,1,ik) + xct%n1b_co*xct%n1b_co*xct%n2b_co*(gi2p-1) &
!                        + xct%n1b_co*xct%n1b_co*(gi2-1) + xct%n1b_co*(gi1p-1) + gi1 -1
!                     else if (xct%ivpar .eq. 0) then
!                       iit = peinf%wown(1,i2p_in,ikp,1,i2_in,ik) + xct%n1b_co*(gi1p-1) + gi1 -1
!                     else
!                       iit = peinf%wown(i1p_in,i2p_in,ikp,i1_in,i2_in,ik)
!                     endif
!                     bset(iit,isc,iscp) = bset(iit,isc,iscp) + outtemp(i1,i2,i1p,i2p)
!                   enddo !i1
!                 enddo !i2
!               enddo !i1p
!             enddo !i2p
!           enddo !isc
!         enddo !iscp
      endif

      SAFE_DEALLOCATE(fxc)

    endif
    call logit('         mtxel_tddft: done bset')

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

    POP_SUB(mtxel_tddft)
    return

  contains


    !> Calculate the contribution from unscreened coulomb at all G vectors.
    !! V ~ delta(G,Gp) * vcoul(q+Gp).
    !! This routine adds in the contribution from V to the body term.
    subroutine calc_direct_unscreened_contrib()

      integer :: t1, t1p, t1_max
      integer :: ofs1, ofs1p
      integer :: n1, n1p, isv
      integer :: k1, k1p, i1p, i1
     
      SCALAR, pointer :: m11p(:,:,:,:)

      integer :: wfn_l(2), wfn_r(2)

      PUSH_SUB(mtxel_tddft.calc_direct_unscreened_contrib)

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
          ! BLASSIFY IF YOU CAN
          do isv=1,xct%nspin
            do i1p=ofs1p+1,ofs1p+n1p
              k1p = i1p - ofs1p
              do i1=ofs1+1,ofs1+n1
                ! We have to exclude the head element here ie vcoul(1)...
                ! as that is handled separately -- vcoul(1) will be mutliplied
                ! later when we interpolate
                k1 = i1 - ofs1
                temph(i1, i1p, isv) = MYCONJG(m11p(1 ,k1, k1p, isv))
                tempb(2:xct%ng, i1, i1p, isv) = tempb(2:xct%ng, i1, i1p, isv) + &
                  vcoul(2:xct%ng)*MYCONJG(m11p(2:xct%ng, k1, k1p, isv))
              enddo
            enddo
          enddo
        enddo !t1p
      enddo !t1
      call timacc(66,2)

      POP_SUB(mtxel_tddft.calc_direct_unscreened_contrib)

    end subroutine calc_direct_unscreened_contrib

    !> Calculates the partial sum over ig for the direct term.
    !! Works for TDA and extended kernels.
    subroutine calc_direct_g_sum()
      integer :: t2, t2p, t2_min
      integer :: ofs2, ofs2p
      integer :: n2, n2p
      SCALAR, pointer :: m22p(:,:,:,:)
      integer :: wfn_l(2), wfn_r(2)

      PUSH_SUB(mtxel_tddft.calc_direct_g_sum)

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
                bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik)
            else 
              call g_sum_TDA(xct,invband,incband,temph,tempb,m22p, &
                bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik)
            endif
            call free_charge_matrix(mats, xct, mccp)
          else
            if (ivsaved) then
              call g_sum_extended(xct,ofs2,ofs2p,n2,n2p,temph_old,tempb_old,m22p, &
                bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik)
            else 
              call g_sum_extended(xct,ofs2,ofs2p,n2,n2p,temph,tempb,m22p, &
                bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik)
            endif
          endif
          call timacc(73,2)
        enddo !t2p
      enddo !t2
      if (.not.ivsaved) then
        SAFE_DEALLOCATE(temph)
        SAFE_DEALLOCATE(tempb)
      endif
      call timacc(64,2)

      POP_SUB(mtxel_tddft.calc_direct_g_sum)

    end subroutine calc_direct_g_sum

    !> Calculates all possible kernel exchange blocks:
    !! \sum_G [M_12(G)]^* v(G+q) M_1p2p(G)
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

      PUSH_SUB(mtxel_tddft.calc_exchange_extended)

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

      POP_SUB(mtxel_tddft.calc_exchange_extended)

    end subroutine calc_exchange_extended


    !> Populates the array vcoul with v_G(q), where q=qq, if qq is given,
    !! or q=q0, if qq is omitted. If vbar is set to true, we zero out the
    !! G=0 component to get the modified Coulomb potential vbar.
    subroutine get_vcoul(vbar, qq)
      logical, intent(in) :: vbar
      real(DP), intent(in), optional :: qq(3)

      integer :: ikpt, ik

      PUSH_SUB(mtxel_tddft.get_vcoul)

      call timacc(62,1)
      vcoul(:) = 0.0d0
      ikpt = 0
      if (present(qq)) then
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
          vq(:) = qg%f(:,ik)
          qlen = DOT_PRODUCT(vq,MATMUL(crys%bdot,vq))
          if (qlen < TOL_Zero) then
            ikpt = ik
            exit
          endif
        enddo
      endif
      if (ikpt == 0) then
        call die("Couldn't find q-point")
      endif
      if (vbar) then
        vcoul(:) = vcoularray(:,ikpt)
        vcoul(1) = 0D0
      else
        vcoul(:) = vcoularray_mod(:,ikpt)
      endif

      call timacc(62,2)

      POP_SUB(mtxel_tddft.get_vcoul)

    end subroutine get_vcoul

    !> Compute all real-space WFNs to speed up the calculation of the charge
    !! density matrix later on. Call this function for each type of wave function
    !! (valence/conduction) at each k-point (k or kp).
    subroutine precompute_ffts()
      type(wavefunction), pointer :: wfn
      integer :: wtype, wprime
      integer :: is, ib, offset

      PUSH_SUB(mtxel_tddft.precompute_ffts)

      do wtype = 1,2 ! loop valence/conduction wfns
        offset = (wtype-1)*invband ! 0 for valence, invband for conduction
        do wprime = 1,2 ! at k/lp
          wfn => wfns(wtype, wprime)%p
          do is = 1, xct%nspin
            do ib = 1, nbands(wtype) ! offset + ib = local band index
              call put_into_fftbox(wfn%ng, wfn%cg(1:,ib,is), gvec%components, &
                wfn%isort, fftboxes(:,:,:,is,offset+ib,wprime), Nfft)
              call do_FFT(fftboxes(:,:,:,is,offset+ib,wprime), Nfft, 1)
            enddo
          enddo
        enddo
      enddo

      POP_SUB(mtxel_tddft.precompute_ffts)   

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
      integer :: ib1, ib2, is
      !> Order to get G-vectors. Should be (1,2,...) unless there`s umklapp involved.
      integer, pointer :: gvecs_order(:)

      PUSH_SUB(mtxel_tddft.calc_charge_matrix)

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
        do is = 1,xct%nspin
          do ib1 = 1, nbands(w1(1))
            if (save_ffts) then
              fftbox1(:,:,:) = CONJG(fftboxes(:,:,:,is,invband*(w1(1)-1)+ib1,w1(2)))
            else
              call put_into_fftbox(wfn1%ng,wfn1%cg(1:,ib1,is),gvec%components,wfn1%isort,fftbox1,Nfft)
              call do_FFT(fftbox1,Nfft,1)
              call conjg_fftbox(fftbox1,Nfft)
            endif
            do ib2 = 1, nbands(w2(1))
              if (save_ffts) then
                ! FHJ: If we use multiply_fftboxes directly we endup moving
                ! more data than necessary here.
                fftbox2 = fftbox1*fftboxes(:,:,:,is,invband*(w2(1)-1)+ib2,w2(2))
              else
                call put_into_fftbox(wfn2%ng,wfn2%cg(1:,ib2,is),gvec%components,wfn2%isort,fftbox2,Nfft)
                call do_FFT(fftbox2,Nfft,1)
                call multiply_fftboxes(fftbox1,fftbox2,Nfft)
              endif
              call do_FFT(fftbox2,Nfft,1)
              call get_from_fftbox(xct%ng,m12(1:,ib1,ib2,is),gvec%components,gvecs_order,fftbox2,Nfft,scale)
            enddo !ib2
          enddo !ib1
        enddo !is
        if (present(m12_old).and.present(should_save)) then
          if (should_save) m12_old(:,:,:,:) = m12(:,:,:,:)
        endif
      else
        if (present(m12_old)) m12(:,:,:,:) = m12_old(:,:,:,:)
      endif

      POP_SUB(mtxel_tddft.calc_charge_matrix)

    end subroutine calc_charge_matrix

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

      PUSH_SUB(mtxel_tddft.get_charge_matrix)

      if (associated(mats(w1(1),w1(2),w2(1),w2(2))%p)) then
        m12 => mats(w1(1),w1(2),w2(1),w2(2))%p
      else
        nb1 = nbands(w1(1))
        nb2 = nbands(w2(1))
        ! FHJ: don`t put the SAFE_ALLOCATE before the "associated" statement!
        if (associated(mats(w2(1),w2(2),w1(1),w1(2))%p)) then
          SAFE_ALLOCATE(mats(w1(1),w1(2),w2(1),w2(2))%p, (xct%ng,nb1,nb2,xct%nspin))
          m12 => mats(w1(1),w1(2),w2(1),w2(2))%p
          m12 = ZERO
          do i2 = 1, nb2
            do i1 = 1, nb1
              do ig=1,xct%ng
                m12(ig,i1,i2,:) = MYCONJG(mats(w2(1),w2(2),w1(1),w1(2))%p(g_idx%mg(ig),i2,i1,:))
              enddo
            enddo
          enddo
        else
          SAFE_ALLOCATE(mats(w1(1),w1(2),w2(1),w2(2))%p, (xct%ng,nb1,nb2,xct%nspin))
          m12 => mats(w1(1),w1(2),w2(1),w2(2))%p
          call calc_charge_matrix(w1, w2, m12, should_calculate, m12_old, should_save)
        endif
      endif

      POP_SUB(mtxel_tddft.get_charge_matrix)

    end subroutine get_charge_matrix

  end subroutine mtxel_tddft

    subroutine tdlda(xct,crys,gvec,Nfft,rhog,fxc)
      implicit none

      type (xctinfo), intent(in) :: xct
      type (crystal), intent(in) :: crys
      type (gspace), intent(in) :: gvec 
      integer, dimension(3), intent(in) :: Nfft
      SCALAR, intent(in) :: rhog(gvec%ng)
      complex(DPC), intent(out) :: fxc(Nfft(1),Nfft(2),Nfft(3))

      integer :: ngrid,nx,ny,nz,i,j,k,l
      real(DP) :: r, rs, A, B, C, D, X, damma, Beta1, Beta2
      ! Rho(G), read from file 
      complex(DPC), allocatable :: rho(:,:,:)
      integer, allocatable :: isort(:)
      
      
      PUSH_SUB(mtxel_tddft.tdlda)

      ! Grid size
      nx = Nfft(1)
      ny = Nfft(2)
      nz = Nfft(3)
      

      !  Allocate
      SAFE_ALLOCATE(rho,(nx,ny,nz))
      SAFE_ALLOCATE(isort,(gvec%ng))

      ! Sort G-vectors, trivial in this case
      do j=1,gvec%ng
        isort(j)=j;
      enddo

      ! Put data in rhog(:) into rho(:,:,:)
      call put_into_fftbox(gvec%ng,rhog, gvec%components,isort,rho,Nfft)

      ! Convert to real space rho(r)
      call do_FFT(rho,Nfft,1)
      rho(:,:,:)=rho(:,:,:)/crys%celvol

      ! Now, calculate fxc from rho

      ! LDA parameter
      ! Taken from Vasiliev, et al., PRB 65, 115416 (2002)
      A =  0.0311
      B = -0.048
      C = -0.0015
      D = -0.0116
      X =  0.0036
      damma = -0.1423
      Beta1 =  1.0529
      Beta2 =  0.3334
    
      do j=1,nx
        do k=1,ny
          do l=1,nz
    !           r = rho(j,k,l)/(nx*ny*nz)
            r = real(rho(j,k,l), dp)
            rs = 3.0d0/(4.0d0*PI_D*r)
            rs = rs**(1./3.)
            !write(6,*) ' rho, rs', r, rs
! Add a minimum density check --- fxc is zero when density is less than 10^-6
            if ((r.le.TOL_ZERO)) then
              fxc(j,k,l)=0.0d0
            else
              if(rs.lt.1.0d0) then
                 fxc(j,k,l) = -1/(9.*r)*(3.*A+(2.*D+C)*rs + &
                     2.*C*rs*log(rs)+ X*(rs**2.)*(2.*log(rs)-1)) - (1/(9.*PI_D*r**2))**(1./3.) 
                 !fxc(j,k,l) = 2.*fxc(j,k,l)
              else
                 fxc(j,k,l) = 5.*Beta1*sqrt(rs)+(7.*(Beta1**2)+8.*Beta2)*rs+ &
                     21.*Beta1*Beta2*(rs**(3./2.))+16.*(Beta2**2)*(rs**2)
                 fxc(j,k,l) = damma*fxc(j,k,l)/(36.*r*(1+Beta1*sqrt(rs)+Beta2*rs)**3.)  &
                     - (1/(9.*PI_D*r**2))**(1./3.)
                 !fxc(j,k,l) = 2.*fxc(j,k,l)
              endif
            endif
          end do !l
        end do !k
      end do !j

      SAFE_DEALLOCATE(rho)
      SAFE_DEALLOCATE(isort)
      POP_SUB(mtxel_tddft.tdlda)
  
      return
  end subroutine tdlda

!SRA: adding lda-based tdrsh, namely short-range Slater fxc
  subroutine tdrsh(xct,crys,gvec,Nfft,rhog,fxc)
      implicit none

      type (xctinfo), intent(in) :: xct
      type (crystal), intent(in) :: crys
      type (gspace), intent(in) :: gvec 
      integer, dimension(3), intent(in) :: Nfft
      SCALAR, intent(in) :: rhog(gvec%ng)
      complex(DPC), intent(out) :: fxc(Nfft(1),Nfft(2),Nfft(3))
      complex(DPC) :: fc,fxsr,fxlr,fxlda

      integer :: ngrid,nx,ny,nz,i,j,k,l
      real(DP) :: r, rs, A, B, C, D, X, damma, Beta1, Beta2, gamma, thrpis
      ! Rho(G), read from file 
      complex(DPC), allocatable :: rho(:,:,:)
      integer, allocatable :: isort(:)
      
      
      PUSH_SUB(mtxel_tddft.tdrsh)

      ! Grid size
      nx = Nfft(1)
      ny = Nfft(2)
      nz = Nfft(3)

      !  Allocate
      SAFE_ALLOCATE(rho,(nx,ny,nz))
      SAFE_ALLOCATE(isort,(gvec%ng))

      ! Sort G-vectors, trivial in this case
      do j=1,gvec%ng
        isort(j)=j;
      enddo

      ! Put data in rhog(:) into rho(:,:,:)
      call put_into_fftbox(gvec%ng,rhog, gvec%components,isort,rho,Nfft)

      ! Convert to real space rho(r)
      call do_FFT(rho,Nfft,1)
      rho(:,:,:)=rho(:,:,:)/crys%celvol

      ! Now, calculate fxc from rho

      ! LDA parameter
      ! Taken from Vasiliev, et al., PRB 65, 115416 (2002)
      A =  0.0311
      B = -0.048
      C = -0.0015
      D = -0.0116
      X =  0.0036
      damma = -0.1423
      Beta1 =  1.0529
      Beta2 =  0.3334

      gamma = xct%coulomb_mod%screening_length*BOHR

      do j=1,nx
        do k=1,ny
          do l=1,nz
    !           r = rho(j,k,l)/(nx*ny*nz)
            r = real(rho(j,k,l), dp)
            rs = 3.0d0/(4.0d0*PI_D*r)
            rs = rs**(1./3.)
            thrpis = 3.0d0*PI_D**2
! Add a minimum density check --- fxc is zero when density is less than 10^-6
            if ((r.le.TOL_ZERO)) then
              fxc(j,k,l)=0.0d0
            else
              if(rs.lt.1.0d0) then
                 !fxc(j,k,l) = -1/(9.*r)*(3.*A+(2.*D+C)*rs + &
                 !    2.*C*rs*log(rs)+ X*(rs**2.)*(2.*log(rs)-1)) + (1/(9.*PI_D*r**2))**(1./3.)&
                 !    *(gamma**2/((thrpis*r)**2/3)*(1-exp(-((thrpis*r)**2/3)/(gamma**2)))-1)
                 fc = -1/(9.*r)*(3.*A+(2.*D+C)*rs + &
                     2.*C*rs*log(rs)+ X*(rs**2.)*(2.*log(rs)-1)) 
                 fxlda = - (1/(9.*PI_D*r**2))**(1./3.)
                 fxsr = (1/(9.*PI_D*r**2))**(1./3.)&
                     *(gamma**2/((thrpis*r)**2/3)*(1-exp(-((thrpis*r)**2/3)/(gamma**2)))-1)
                 fxlr=fxlda-fxsr
                 fxc(j,k,l)=fc+(1-xct%coulomb_mod%short_range_frac_fock)*fxsr+(1-xct%coulomb_mod%long_range_frac_fock)*fxlr
                 !fxc(j,k,l) = 2.*fxc(j,k,l)
              else
                 !fxc(j,k,l) = 5.*Beta1*sqrt(rs)+(7.*(Beta1**2)+8.*Beta2)*rs+ &
                 !    21.*Beta1*Beta2*(rs**(3./2.))+16.*(Beta2**2)*(rs**2)
                 !fxc(j,k,l) = damma*fxc(j,k,l)/(36.*r*(1+Beta1*sqrt(rs)+Beta2*rs)**3.)  &
                 !    + (1/(9.*PI_D*r**2))**(1./3.)&
                 !    *(gamma**2/((thrpis*r)**2/3)*(1-exp(-((thrpis*r)**2/3)/(gamma**2)))-1)
                 fc = 5.*Beta1*sqrt(rs)+(7.*(Beta1**2)+8.*Beta2)*rs+ &
                     21.*Beta1*Beta2*(rs**(3./2.))+16.*(Beta2**2)*(rs**2)
                 fc = damma*fc/(36.*r*(1+Beta1*sqrt(rs)+Beta2*rs)**3.)
                 fxlda = - (1/(9.*PI_D*r**2))**(1./3.)
                 fxsr = (1/(9.*PI_D*r**2))**(1./3.)&
                     *(gamma**2/((thrpis*r)**2/3)*(1-exp(-((thrpis*r)**2/3)/(gamma**2)))-1)
                 fxlr=fxlda-fxsr
                 fxc(j,k,l)=fc+(1-xct%coulomb_mod%short_range_frac_fock)*fxsr+(1-xct%coulomb_mod%long_range_frac_fock)*fxlr

                 !fxc(j,k,l) = 2.*fxc(j,k,l)
              endif
            endif
          end do !l
        end do !k
      end do !j

      SAFE_DEALLOCATE(rho)
      SAFE_DEALLOCATE(isort)
      POP_SUB(mtxel_tddft.tdrsh)
  
      return
  end subroutine tdrsh



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
      POP_SUB(free_charge_matrices)
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
  
end module mtxel_tddft_m
