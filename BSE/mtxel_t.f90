!==============================================================================
!
! Routines:
!
! (1) mtxel_t()         Originally By MLT       Last Modified 02/13/2013 (FHJ)
!
!     Calculated matrix elements needed for interpolation of kernel (dcn,dvn)
!
!     input: gvec, wfnc_co, wfnc_fi, wfnv_co, wfnvq_fi types
!            igumk        index of umklapp, translation, vector
!            ikpt         index of the current k-point
!
!     output: dcc, dvv     transformation matrices between a k-point in the
!                          fine grid and one in the coarse grid
!
!==============================================================================

#include "f_defs.h"

module mtxel_t_m

  use global_m

  implicit none

  public :: mtxel_t

  private

contains

  !> Calculates the transformation matrix between the coarse and fine
  !! wavefunction. Can use full or restricted transformation, where restricted
  !! means cond states are expanded in terms of other cond. states, etc.
  !! TODO: this whole function can be BLAS-ified
  subroutine mtxel_t(gvec,wfnc_co,wfnc_fi,wfnv_co,wfnvq_fi,dcn,dvn,igumk,igumkq,restricted,qflag)
    type (gspace), intent(in) :: gvec
    type (wavefunction), intent(in) :: wfnc_co,wfnc_fi,wfnv_co,wfnvq_fi
    !> Restricted: (wfnc_fi%nband, wfnc_co%nband, wfnc_fi%nspin)
    !! Unrestric.: (wfnc_fi%nband, wfnc_co%nband+wfnv_co%nband, wfnc_fi%nspin)
    SCALAR, intent(inout), target :: dcn(:,:,:)
    !> Restricted: (wfnvq_fi%nband, wfnv_co%nband, wfnc_fi%nspin)
    !! Unrestric.: (wfnvq_fi%nband, wfnv_co%nband+wfnc_co%nband, wfnc_fi%nspin)
    SCALAR, intent(inout), target :: dvn(:,:,:)
    integer, intent(in) :: igumk,igumkq, qflag
    logical, intent(in) :: restricted

    integer :: ig,ig_co
    integer, allocatable :: isorti(:), isortiq(:)
    SCALAR, pointer :: dmm(:,:,:)
    
    PUSH_SUB(mtxel_t)
    
    dcn(:,:,:)=ZERO
    dvn(:,:,:)=ZERO
    
    ! Compute inverse array to wfnc_co%isort
    
    SAFE_ALLOCATE(isorti, (gvec%ng))
    isorti=0
    do ig=1,gvec%ng
      isorti(wfnc_co%isort(ig))=ig
    enddo
    
    SAFE_ALLOCATE(isortiq, (gvec%ng))
    isortiq=0
    do ig=1,gvec%ng
      isortiq(wfnv_co%isort(ig))=ig
    enddo
    
    ! Loop over G-vectors
    ! Compute matrix element <c_co,k_co|exp(i(k_co-k_fi).r)|c_fi,k_fi>
    
    ! fine grid = valence
    ! v_fi <- v_co
    dmm => dvn(:,:wfnv_co%nband,:)
    do ig=1,wfnvq_fi%ng
      call map_ig(gvec, ig, ig_co, igumkq, isortiq, wfnvq_fi, wfnv_co)
      if (ig_co==0) cycle
      call mtxel_dmm(ig, ig_co, wfnvq_fi%nspinor, wfnvq_fi, wfnv_co, dmm)
    enddo
    if (qflag==1 .and. .not. restricted) then
      ! v_fi <- c_co
      dmm => dvn(:,wfnv_co%nband+1:,:)
      do ig=1,wfnvq_fi%ng
        call map_ig(gvec, ig, ig_co, igumkq, isortiq, wfnvq_fi, wfnc_co)
        if (ig_co==0) cycle
        call mtxel_dmm(ig, ig_co, wfnvq_fi%nspinor, wfnvq_fi, wfnc_co, dmm)
      enddo
    endif
    
    ! fine grid = conduction
    if (restricted) then
      ! c_fi <- c_co
      dmm => dcn(:,:,:)
      do ig=1,wfnc_fi%ng
        call map_ig(gvec, ig, ig_co, igumk, isorti, wfnc_fi, wfnc_co)
        if (ig_co==0) cycle
        call mtxel_dmm(ig, ig_co, wfnc_fi%nspinor, wfnc_fi, wfnc_co, dmm)
      enddo
    else
      ! c_fi <- c_co
      dmm => dcn(:,wfnv_co%nband+1:,:)
      do ig=1,wfnc_fi%ng
        call map_ig(gvec, ig, ig_co, igumk, isorti, wfnc_fi, wfnc_co)
        if (ig_co==0) cycle
        call mtxel_dmm(ig, ig_co, wfnc_fi%nspinor, wfnc_fi, wfnc_co, dmm)
      enddo
      if (qflag==1) then
        ! c_fi <- v_co
        dmm => dcn(:,:wfnv_co%nband,:)
        do ig=1,wfnc_fi%ng
          call map_ig(gvec, ig, ig_co, igumk, isorti, wfnc_fi, wfnv_co)
          if (ig_co==0) cycle
          call mtxel_dmm(ig, ig_co, wfnc_fi%nspinor, wfnc_fi, wfnv_co, dmm)
        enddo
      endif
    endif
    
    SAFE_DEALLOCATE(isorti)
    SAFE_DEALLOCATE(isortiq)
    
    POP_SUB(mtxel_t)
    
    return
    
  end subroutine mtxel_t


  !> Maps the G vector (ig) from a fine WFN to the coarse one, 
  !! including umklapp. Resulting index is saved to ig_co (0 if invalid)
  subroutine map_ig(gvec, ig, ig_co, ig_umklapp, isorti_co, wfn_fi, wfn_co)
    type(gspace), intent(in) :: gvec
    integer, intent(in) :: ig
    integer, intent(out) :: ig_co
    integer, intent(in) :: ig_umklapp
    integer, intent(in) :: isorti_co(:)
    type(wavefunction), intent(in) :: wfn_fi, wfn_co
    
    integer :: igadd, igaddn, kn(3), kaddn
    
    ! no push/pop, called too frequently
    ig_co = 0
    if (ig_umklapp==1) then
      ig_co = isorti_co(wfn_fi%isort(ig))
    else
      igadd = wfn_fi%isort(ig)
      
      ! Compute address of gn=g+gumk:
      ! If the result is beyond planewaves included in wavefunction, skip
      kn(1:3) = gvec%components(1:3,igadd) + &
        gvec%components(1:3,ig_umklapp) + gvec%FFTgrid(1:3)/2 + 1
      if (any(kn(1:3) < 1) .or. any(kn(1:3) > gvec%FFTgrid(1:3))) then
        ig_co = 0
        return
      endif
      kaddn = ((kn(1) - 1)*gvec%FFTgrid(2) + kn(2) - 1)*gvec%FFTgrid(3) + kn(3)
      igaddn = gvec%index_vec(kaddn) ! relate the cube to the sphere
      ig_co = isorti_co(igaddn)
    endif
    if (ig_co > wfn_co%ng) ig_co = 0
    
  end subroutine map_ig
  
  
  !> Calculates the overlap between the coarse and fine wavefunction, store
  !! info the dmm transformation matrix. Assumes you know the mapping between
  !! the the fine g-vector ig and the coarse ig_co
  subroutine mtxel_dmm(ig, ig_co, nspinor, wfn_fi, wfn_co, dmm)
    integer, intent(in) :: ig, ig_co
    integer, intent(in) :: nspinor
    type(wavefunction), intent(in) :: wfn_fi, wfn_co
    SCALAR, intent(inout) :: dmm(:,:,:) !< (wfn_fi%nband, wfn_co%nband, nspin)
    
    integer :: im_fi, im_co, is, ispinor
    
    ! no push/pop, called too frequently
    
    do is=1,wfn_co%nspin
      do im_co=1,wfn_co%nband
        do im_fi=1,wfn_fi%nband
          do ispinor=1,nspinor
            dmm(im_fi, im_co, is) = dmm(im_fi, im_co, is) + &
              wfn_fi%cg(ig, im_fi, is*ispinor) * MYCONJG(wfn_co%cg(ig_co, im_co, is*ispinor))
          enddo
        enddo
      enddo
    enddo
    
  end subroutine mtxel_dmm

end module mtxel_t_m

