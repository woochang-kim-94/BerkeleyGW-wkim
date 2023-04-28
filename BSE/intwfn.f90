!=============================================================================
!
! Module:
!
! (1) intwfn()  Originally By MLT               Last Modified 7/8/2008 (JRD)
!
!     Evaluate the transformation matrices and store them
!     in the arrays dcc,dvv.
!     The matrices are defined according to eq. (39) of Rohlfing & Louie:
!
!     dcc(ic,jc,ik,is,ivert) == d_{ic,ik}^{\tilde jc,\tilde jk)
!        ik: index of k-point in fine grid
!        ic: index of cond. band in fine grid
!        jk: index of k-point in coarse grid, fi2co_wfn(ivert,ik)=jk
!        jc: index of cond. band in coarse grid
!        is: spin index (the same in both grids)
!
!     dvv(iv,jv,ik,is,ivert) == d_{iv,ik}^{\tilde jv,\tilde jk)
!        ik: index of k-point in fine grid
!        iv: index of valence band in fine grid
!        jk: index of k-point in coarse grid, fi2co_wfn(ivert,ik)=jk
!        jv: index of valence band in coarse grid
!        is: spin index (the same in both grids)
!
!     Mixing between valence/conduction states is taken into account
!     by allowing jc to be also an occupied band, and similarly
!     for jv.
!
!     input: crys,syms,xct,gvec,kg_fi,kgq_fi types
!            indexq_fi (mapping between shifted and unshifted k-points)
!
!     output: dcc,dvv arrays
!             kco array
!             fi2co_wfn array
!
! (2) intwfn_sub()  DYQ
!
!     Based on intwfn, this module evaluates the transformation matrices 
!     between fine and subsampled wavefunctions and stores them
!     in the arrays dcc_sub,dvv_sub.
!     The matrices are defined according to eq. (39) of Rohlfing & Louie
!     with k_co` replaced by the subsampled point in every
!     instance.
!
!===============================================================================

#include "f_defs.h"

module intwfn_m
  use checkbz_m
  use fullbz_m
  use genwf_co_m
  use genwf_m
  use global_m
  use gmap_m
  use input_co_m
  use input_co_q_m
  use input_utils_m
  use interp_m
  use intpts_m
  use io_utils_m
  use misc_m
  use mtxel_t_m
  use timing_m, only: timing => bse_timing
  use wfn_rho_vxc_io_m
  implicit none

  private

  public :: intwfn, intwfn_sub

contains

subroutine intwfn(kp,crys,syms,xct,flag,gvec,kg_fi,kgq_fi,kg_co,kgq_co, &
  dcc,dvv,kco,fi2co_wfn,indexq_fi,eqp,intwfnv,intwfnc,intp_coefs)
  type (kpoints), intent(in) :: kp
  type (crystal), intent(in) :: crys
  type (symmetry), intent(in) :: syms
  type (xctinfo), intent(inout) :: xct
  type (flags), intent(in) :: flag
  type (gspace), intent(in) :: gvec
  type (grid), intent(in) :: kg_fi
  type (grid), intent(in) :: kgq_fi
  !> (xct%ncb_fi,xct%n2b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel) Interp coefs.
  SCALAR, intent(out), target :: dcc(:,:,:,:,:)
  !> (xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel) Interp coefs.
  SCALAR, intent(out), target :: dvv(:,:,:,:,:)
  real(DP), intent(out) :: kco(:,:) !< 3,xct%nkpt_co
  integer, intent(out) :: fi2co_wfn(:,:) !< (xct%npts_intp_kernel, xct%nkpt_fi)
  integer, intent(in) :: indexq_fi(:) !< xct%nkpt_fi
  type (eqpinfo), intent(inout) :: eqp
  type (int_wavefunction), intent(inout) :: intwfnc
  type (int_wavefunction), intent(inout) :: intwfnv
  !> (xct%npts_intp_kernel, xct%nkpt_fi) Delaunay/greedy interpolation coefficients
  real(DP), intent(out) :: intp_coefs(:,:)
  type(grid),intent(out) :: kg_co,kgq_co

  integer :: closepts(4),umk
  real(DP) :: closeweights(4),qq_temp(3),delta
  type(work_genwf) :: work, workq, workco, workcoq
  type (wavefunction) :: wfnc_co,wfnc_fi
  type (wavefunction) :: wfnv_co,wfnvq_fi
  type (tdistgwf) :: distgwfco,distgwfcoq

  integer :: ik,ikt,ik_co,ik_fi,ikq_fi,vunit,eqp_dat_size,ikq
  integer :: ik_co_q, igumkq, ijk, ivert
  integer :: ic,ic_co,ic_fi,iv,iv_co,iv_fi,is
  integer :: ii,jj,kk,gumk(3),igumk
  real(DP) :: dist,qg(3),dnorm,tempshift,tempsum,dweight, worst_normc, worst_normv
  real(DP), allocatable :: normc(:,:,:,:),normv(:,:,:,:)
  integer, allocatable :: idum(:,:)
  SCALAR, allocatable :: dummy(:,:,:,:,:)
  real(DP), allocatable :: edummy(:,:,:),  dummy_coefs(:,:)
  type(interp_t) :: wfn_interp,wfn_interp_q
  !> (xct%nvb_fi, xct%n1b, xct%nspin, xct%nkpt_fi, xct%npts_intp_kernel)
  !> dvn points to dvv if we are restricting the transf right away,
  !! or if we never restrict it (i.e., extended_kernel==.true.). Otherwise,
  !! we allocate a buffer, calculate the unrestricted transf, and restrict+copy
  ! dvn to dvv later on.
  SCALAR, pointer :: dvn(:,:,:,:,:) 
  SCALAR, pointer :: dcn(:,:,:,:,:) !< See: dvn
  !> .true. if we calculate the unrestricted transformation, but the kernel is
  !! restricted (i.e., we`ll have to truncate the dvn to dvv later on).
  logical :: truncate_coefs
  type(progress_info) :: prog_info
  type(kpoints) ::kp_co
 ! variables used in finite Q calculcation only
  integer :: closepts_q(4)
  real(DP) :: closeweights_q(4),dweight_q

  PUSH_SUB(intwfn)

  ! FHJ: Normally, xct%n1b_co depends only on xct%extended_kernel. However,
  ! we`ll temporarily overwrite this value depending on xct%unrestricted_transf
  ! till the end of the subroutine, because if xct%unrestricted_transf==.true,
  ! intwfn "thinks" it`s dealing with an extended kernel when calculating dvn/dcn.
  if (xct%unrestricted_transf) then
    xct%n1b_co = xct%nvb_co + xct%ncb_co
    xct%n2b_co = xct%n1b_co
  else
    xct%n1b_co = xct%nvb_co
    xct%n2b_co = xct%ncb_co
  endif
  truncate_coefs = xct%unrestricted_transf.and.(.not.xct%extended_kernel)
  if (truncate_coefs) then
    ! FHJ: we`ll restrict the coefficients later on, but we want to do the
    ! right expansion here => need a temp buffer.
    SAFE_ALLOCATE(dvn, (xct%nvb_fi, xct%n1b_co, xct%nspin, xct%nkpt_fi, xct%npts_intp_kernel))
    SAFE_ALLOCATE(dcn, (xct%ncb_fi, xct%n2b_co, xct%nspin, xct%nkpt_fi, xct%npts_intp_kernel))
    dvn = 0.d0
    dcn = 0.d0
  else
    ! FHJ: we either do not restrict, or restrict right from the beginning.
    dvn => dvv
    dcn => dcc
  endif
  dcc = 0.d0
  dvv = 0.d0

  SAFE_ALLOCATE(normc, (xct%nkpt_fi,xct%ncb_fi,xct%nspin,xct%npts_intp_kernel))
  SAFE_ALLOCATE(normv, (xct%nkpt_fi,xct%nvb_fi,xct%nspin,xct%npts_intp_kernel))
  intp_coefs = 0d0
  fi2co_wfn = 0
  if (flag%read_dtmat) then
    if (peinf%inode.eq.0) then

!-------------------------------
! Read transformation matrices, if they exist...
! FHJ: TODO - test this part for xct%unrestricted_transf=.true.

      if (xct%unrestricted_transf) then
        write(6,*) 'Reading dcn,dvn matrices from file'
      else
        write(6,*) 'Reading dcc,dvv matrices from file'
      endif
      
! Check if dtmat has the right parameters

      call open_file(unit=13,file='dtmat',form='unformatted',status='old')
      read(13)
      read(13) ik_co,ic_co,iv_co,ik,ic_fi,iv_fi,is

      if(ik_co /= xct%nkpt_co) then
        write(0,*) 'File has ', ik_co, '; we need ', xct%nkpt_co
        call die("dtmat does not have the correct number of k-points in coarse grid")
      endif
      if(ic_co /= xct%n2b_co) then
        write(0,*) 'File has ', ic_co, '; we need ', xct%n2b_co
        call die("dtmat does not have the correct number of conduction bands in coarse grid")
      endif
      if(iv_co /= xct%n1b_co) then
        write(0,*) 'File has ', iv_co, '; we need ', xct%n1b_co
        call die("dtmat does not have the correct number of valence bands in coarse grid")
      endif
      if(ik /= xct%nkpt_fi) then
        write(0,*) 'File has ', ik, '; we need ', xct%nkpt_fi
        call die("dtmat does not have the correct number of k-points in fine grid")
      endif
      if(ic_fi /= xct%ncb_fi) then
        write(0,*) 'File has ', ic_fi, '; we need ', xct%ncb_fi
        call die("dtmat does not have the correct number of conduction bands in fine grid")
      endif
      if(iv_fi /= xct%nvb_fi) then
        write(0,*) 'File has ', iv_fi, '; we need ', xct%nvb_fi
        call die("dtmat does not have the correct number of valence bands in fine grid")
      endif
      if(is /= xct%nspin) then
        write(0,*) 'File has ', is, '; we need ', xct%nspin
        call die("dtmat does not have the correct number of spins")
      endif
          
! kco(1:3,1:xct%nkpt_co) has the coordinates of k-points in
! the coarse grid used to generate "dtmat" (in general, this
! is equal to kg_co%fk(1:3,1:xct%nkpt_co) of a previous run)
!
! fi2co_wfn(1:xct%nkpt_fi) is a mapping array between a point in
! the fine grid and the closest ones in the coarse grid

      do jj = 1, xct%nkpt_co
        read(13) kco(1,jj),kco(2,jj),kco(3,jj)
      enddo
      do ivert = 1, xct%npts_intp_kernel
        do jj = 1, xct%nkpt_fi * xct%ncb_fi * xct%n2b_co * xct%nspin
          read(13) ik_fi,ic_fi,ik_co,ic_co,is,dcn(ic_fi,ic_co,is,ik_fi,ivert)
        enddo
      enddo
      do ivert = 1, xct%npts_intp_kernel
        do jj = 1, xct%nkpt_fi * xct%nvb_fi * xct%n1b_co * xct%nspin
          read(13) ik_fi,iv_fi,ik_co,iv_co,is,dvn(iv_fi,iv_co,is,ik_fi,ivert)
          fi2co_wfn(ivert,ik_fi) = ik_co
        enddo
      enddo

      intp_coefs(:,:) = 1d0
      if (xct%npts_intp_kernel>1) then
        read(13) intp_coefs(:,:)
      endif
      call close_file(13)
    endif

! PE # 0 distribute data

#ifdef MPI
    call MPI_BCAST(fi2co_wfn(1,1), xct%npts_intp_kernel*xct%nkpt_fi, MPI_INTEGER, 0, MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(kco, 3*xct%nkpt_co, MPI_REAL_DP, 0, MPI_COMM_WORLD,mpierr)
    call MPI_BCAST(dcn(1,1,1,1,1), xct%nkpt_fi*xct%ncb_fi*xct%n2b_co*xct%nspin*xct%npts_intp_kernel, &
      MPI_SCALAR, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(dvn(1,1,1,1,1), xct%nkpt_fi*xct%nvb_fi*xct%n1b_co*xct%nspin*xct%npts_intp_kernel, &
      MPI_SCALAR, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(intp_coefs(1,1), xct%npts_intp_kernel*xct%nkpt_fi, MPI_REAL_DP, 0, MPI_COMM_WORLD,mpierr)
#endif

  else

!-------------------------------
! ...Otherwise, start the computation of them

! Read the wavefunctions on the coarse grid and initialize dtmat

    call timing%start(timing%iw_input_co)
    call input_co(kp,kp_co,crys,gvec,kg_co,kgq_co,syms,xct,flag%bzc,distgwfco,eqp)

    if (xct%qflag .eq. 0) then
      if (peinf%inode .eq. 0) then
        write(6,999)
      endif
    endif
999 format(/,1x,'We are doing a Finite Q Calculation',/)

    if (xct%qflag.eq.0) then
      call input_co_q(kp,kp_co,crys,gvec,kg_co,kgq_co,syms,xct,flag%bzc,distgwfcoq,eqp)
    endif
    call timing%stop(timing%iw_input_co)

!---------------------------------------------------- 
! If skip_interpolation keyword is present don`t interpolate

    if (xct%skipinterp) then

      intp_coefs = 1d0
      if (peinf%inode.eq.0) then
        write(6,'(a)') ' Skipping interpolation'
        write(6,'(a)')
      endif

      do ii = 1, xct%nkpt_fi
        do jj = 1, xct%nkpt_co
          if (all(abs(kg_fi%f(1:3,ii)-kg_co%f(1:3,jj)) .lt. TOL_Small)) then
            fi2co_wfn(1,ii)=jj
            kco(1:3,jj)=kg_co%f(1:3,jj)
          endif
        enddo
      enddo

      dcn(:,:,:,:,:) = 0d0
      do ik = 1, xct%nkpt_fi
        do is = 1, xct%nspin
          do ic = 1, xct%ncb_fi
            dcn(ic,ic,is,ik,1) = 1d0
          enddo
        enddo
      enddo
      dvn(:,:,:,:,:) = 0d0
      do ik = 1, xct%nkpt_fi
        do is = 1, xct%nspin
          do iv = 1, xct%nvb_fi
            dvn(iv,iv,is,ik,1) = 1d0
          enddo
        enddo
      enddo

    else ! xct%skipinterp
      
      call timing%start(timing%iw_interp)
      if (xct%delaunay_interp) then
        call interp_init(wfn_interp, crys, kg_co%f, dims=xct%idimensions, &
          periodic=.true., active_dims=xct%is_periodic)
        !DYQ: For now, since we can only use one delaunay grid, interpolate WFNq using 
        ! old scheme
        if (xct%qflag.eq.0) then
          call alloc_intpts(kgq_co%nf, kgq_co%f(:,:), periodic=.true.)
          !call interp_init(wfn_interp_q, crys, kgq_co%f, dims=xct%idimensions, &
          !  periodic=.true., active_dims=xct%is_periodic)          
        endif
      else
        call alloc_intpts(kg_co%nf, kg_co%f(:,:), periodic=.true.)
      endif
      call timing%stop(timing%iw_interp)
    
      do kk=1,xct%nkpt_co
        kco(:,kk) = kg_co%f(:,kk)
      enddo
      
      if (peinf%inode.eq.0) then
        call open_file(unit=13,file='dtmat',form='unformatted',status='replace')
        write(13) xct%idimensions, xct%is_periodic(1:3), xct%npts_intp_kernel, xct%nkpt_co
        write(13) xct%nkpt_co,xct%n2b_co,xct%n1b_co, &
          xct%nkpt_fi,xct%ncb_fi,xct%nvb_fi,xct%nspin,xct%npts_intp_kernel
        do kk=1,xct%nkpt_co
          write(13) (kco(ii,kk),ii=1,3)
        enddo
        call close_file(13)
      endif

!-------------------------
! Loop over k-points in fine grid
  
      if (xct%inteqp) then
        SAFE_ALLOCATE(eqp%evshift, (xct%nvb_fi,kg_fi%nf,xct%nspin))
        SAFE_ALLOCATE(eqp%ecshift, (xct%ncb_fi,kg_fi%nf,xct%nspin))
        eqp%evshift=0D0
        eqp%ecshift=0D0
      endif

     call progress_init(prog_info, 'wavefunction interpolation', 'k-point', &
      peinf%nkpe)
     do ii=1,peinf%nkpe
       call progress_step(prog_info)
       ik_fi=peinf%ik(peinf%inode+1,ii)
       
! JRD: Find the points in coarse grid closest to each point in the fine grid
! We need 4 points here for 3D

        if (ik_fi .ne. 0) then
          if (xct%qflag.eq.1) then
            ikq_fi = indexq_fi(ik_fi)
          else
            ikq_fi = xct%indexq_fi(ik_fi)
          endif
        else
          ! this process has no work to do
          ikq_fi = 0
        endif
          
! JRD: In principle, we should be calling the below function twice. Once for ik_fi and once for ikq_fi. If using 
! momentum operator, it makes no difference. And for velocity operator, it makes insignificant difference for typical small q`s.

       if (kg_co%nf .gt. 1 .and. ik_fi > 0) then
          call timing%start(timing%iw_interp)
          if (xct%delaunay_interp) then
            ! Returns closest coarse points to ik_fi 
            call interp_eval(wfn_interp, kg_fi%f(:,ik_fi), closepts(:), closeweights(:))
          else
            call intpts_local(crys, kg_fi%f(:,ik_fi), kg_co%nf, kg_co%f(:,:), &
              xct, closepts(:), closeweights(:), periodic=.true.)
          endif
          ! DYQ: If using finite momentum, we need to call interpolation routines again for the shifted grid
          if (xct%qflag.eq.2) then
            if (xct%delaunay_interp) then
              ! Returns closest coarse points to ikq_fi 
              call interp_eval(wfn_interp, kgq_fi%f(:,ikq_fi), closepts_q(:), closeweights_q(:))
            else
              call intpts_local(crys, kgq_fi%f(:,ikq_fi), kg_co%nf, kg_co%f(:,:), &
                xct, closepts_q(:), closeweights_q(:), periodic=.true.)
            endif
          else if (xct%qflag .eq. 0) then
            if (xct%delaunay_interp) then
              ! TODO: Fix This!!! Using old interpolation scheme instead of delaunay
              !call interp_eval(wfn_interp_q, kgq_fi%f(:,ikq_fi), closepts_q(:), closeweights_q(:))
              call intpts_local(crys, kgq_fi%f(:,ikq_fi), kgq_co%nf, kgq_co%f(:,:), &
                xct, closepts_q(:), closeweights_q(:), periodic=.true.)
            else
              call intpts_local(crys, kgq_fi%f(:,ikq_fi), kgq_co%nf, kgq_co%f(:,:), &
                xct, closepts_q(:), closeweights_q(:), periodic=.true.)
            endif
          endif
          if (xct%npts_intp_kernel>1) then
            intp_coefs(:,ik_fi) = closeweights(1:xct%npts_intp_kernel)
            !if (xct%qflag.ne.1) then
            !  intp_coefs_q(:,ikq_fi) = closeweights_q(1:xct%npts_intp_kernel)
            !endif
          else
            intp_coefs(1,ik_fi) = 1d0
            !if (xct%qflag.ne.1) then
            !  intp_coefs_q(1,ikq_fi) = 1d0
            !endif
          endif
          call timing%stop(timing%iw_interp)
        else
          closepts(:)=1
          closeweights(:)=0D0
          closeweights(1)=1D0
          if (ik_fi>0) intp_coefs(1,ik_fi) = 1d0
          if (xct%qflag.ne.1) then
            closepts_q(:)=1
            closeweights_q(:)=0D0
            closeweights_q(1)=1D0
            !if (ikq_fi>0) intp_coefs_q(1,ikq_fi) = 1d0
          endif
        endif

! FHJ: we need to be careful when interpolating QP corrections near band extrema.
! If the closeweights go to 1 as ~ 1-q as we go away from a coarse point, a direct
! interpolation will introduce a linear term in the QP corrections and break
! the effective masses. A better approach is to use a transformed variable
! as the interpolated weights. The transformation below is simple, yet guarantees
! that the weights approach 1 much faster than 1-q or (1-q)^2 for small q and
! has a nice linear region where weights from different k-points blend in.
! FHJ: (a few years later) the band structure interpolation simply doesn`t
! look nice with this, even for MoS2.
! TODO: we need to fix this...
        if (xct%kpt_interp==KPT_INTERP_EXP_TRANSF) then
          where (abs(closeweights)>TOL_ZERO)
            closeweights = exp(1d0 - 1d0/abs(closeweights))
          elsewhere
            closeweights = 0d0
          endwhere
        endif
        closeweights = closeweights / sum(closeweights)

        do ijk = xct%idimensions + 1, 1, -1
          
          ik_co=closepts(ijk)
          dweight=closeweights(ijk)
          
          if (xct%qflag .eq. 1) then
            ik_co_q=ik_co
          elseif (xct%qflag .eq. 0) then
            !ik_co_q=xct%indexq(ik_co)
            ik_co_q = closepts_q(ijk)
            dweight_q=closeweights_q(ijk)
          elseif (xct%qflag .eq. 2) then
            if (xct%patched_sampling_co .and. xct%indexq(ik_co).eq.0) then
              if (peinf%inode.eq.0) write(6,*) "Skipping coarse point, ik_co",ik_co
              ik_co_q =0
              cycle
            else
              dweight_q=closeweights_q(ijk)
              ik_co_q = xct%indexq(ik_co)
            endif
          endif

          ! Store the umklapp G-vector, if not zero
          if (ik_fi .ne. 0) then
            qg(:)=kg_fi%f(:,ik_fi)-kg_co%f(:,ik_co)
            gumk(:)=nint(qg(:))
            call findvector(igumk, gumk, gvec)
            
            ! Store the umklapp G-vector, if not zero - Finite Q
            if (xct%qflag .eq. 1) then
            !JWJ add, we always do this mapping
               qg(:)=kgq_fi%f(:,ik_fi)-kg_co%f(:,ik_co)
               gumk(:)=nint(qg(:))
               call findvector(igumkq, gumk, gvec)
            !  igumkq=igumk
            elseif (xct%qflag .eq. 0) then
              qg(:)=kgq_fi%f(:,ikq_fi)-kgq_co%f(:,ik_co_q)
              gumk(:)=nint(qg(:))
              call findvector(igumkq, gumk, gvec)
            elseif (xct%qflag .eq.2) then
              qg(:)=kgq_fi%f(:,ikq_fi) - kg_co%f(:,ik_co_q) 
              gumk(:)=nint(qg(:))
              call findvector(igumkq, gumk, gvec)
            endif
          endif
      
!------------------------------
! Read needed wavefunctions: conduction bands from unshifted grid
! and valence bands from shifted grid

          if (ik_fi .ne. 0) then
            call timing%start(timing%iw_genwf)
            call genwf(crys,gvec,kg_fi,syms,wfnc_fi,ik_fi,ik_fi,xct%nspin,xct%ncb_fi,&
                       work,intwfnc,xct%iwriteint,is_cond=.true.)
            call genwf(crys,gvec,kgq_fi,syms,wfnvq_fi,ik_fi,ikq_fi,xct%nspin,xct%nvb_fi,&
                       workq,intwfnv,xct%iwriteint,is_cond=.false.)
            call timing%start(timing%iw_genwf)
          endif
          
          if (xct%qflag .eq. 0) then
            call timing%start(timing%iw_genwf_co)
            call genwf_co(crys,gvec,kg_co,kgq_co,syms,wfnc_co, &
              wfnv_co,xct,ik_co,ik_co_q,distgwfco,distgwfcoq,workco,workcoq)
            call timing%stop(timing%iw_genwf_co)
          else
            call timing%start(timing%iw_genwf_co)
            ! Unchanged for qflag=2
            call genwf_co(crys,gvec,kg_co,kgq_co,syms,wfnc_co, &
              wfnv_co,xct,ik_co,ik_co_q,distgwfco,distgwfco,workco,workcoq)
            call timing%stop(timing%iw_genwf_co)
          endif

          if (ik_fi .ne. 0) then
            call timing%start(timing%iw_mtxel_t)
            ivert = 1
            if (xct%npts_intp_kernel>1) ivert = ijk
            if (xct%qflag.eq.1) then
              call mtxel_t(gvec,wfnc_co, wfnc_fi, wfnv_co, wfnvq_fi, &
                dcn(:,:,:,ik_fi,ivert), dvn(:,:,:,ik_fi,ivert),  &
                igumk, igumkq, .not.xct%unrestricted_transf, xct%qflag)
            else
              call mtxel_t(gvec,wfnc_co, wfnc_fi, wfnv_co, wfnvq_fi, &
                dcn(:,:,:,ik_fi,ivert), dvn(:,:,:,ikq_fi,ivert),  &
                igumk, igumkq, .not.xct%unrestricted_transf, xct%qflag)
            endif
            call timing%start(timing%iw_mtxel_t)
          endif
          
          SAFE_DEALLOCATE_P(wfnc_co%cg)
          SAFE_DEALLOCATE_P(wfnc_co%isort)
          SAFE_DEALLOCATE_P(wfnv_co%cg)
          SAFE_DEALLOCATE_P(wfnv_co%isort)
          
          if (ik_fi .ne. 0) then
            SAFE_DEALLOCATE_P(wfnc_fi%cg)
            SAFE_DEALLOCATE_P(wfnc_fi%isort)
            SAFE_DEALLOCATE_P(wfnvq_fi%cg)
            SAFE_DEALLOCATE_P(wfnvq_fi%isort)
            fi2co_wfn(ivert, ik_fi) = ik_co
          endif

! Construct Interpolated eqp shifts
          
          if (xct%inteqp .and. ik_fi .ne. 0) then
            
            if (xct%qflag.eq.1) then
              do is=1,xct%nspin
                do iv_fi=1,xct%nvb_fi
                  tempshift=0d0
                  tempsum=0d0
                  do iv_co=1,xct%nvb_co
                    tempshift = tempshift + (abs(dvn(iv_fi,iv_co,is,ik_fi,ivert))**2)*&
                      eqp%evshift_co(iv_co,kg_co%indr(ik_co),is)
                    tempsum = tempsum + (abs(dvn(iv_fi,iv_co,is,ik_fi,ivert))**2)
                  enddo
                  do iv_co=xct%nvb_co+1,xct%n1b_co
                    tempshift = tempshift + (abs(dvn(iv_fi,iv_co,is,ik_fi,ivert))**2)*&
                      eqp%ecshift_co(iv_co-xct%nvb_co,kg_co%indr(ik_co),is)
                    tempsum = tempsum + (abs(dvn(iv_fi,iv_co,is,ik_fi,ivert))**2)
                  enddo
                  if (.not.xct%renorm_transf) tempsum = 1d0
                  eqp%evshift(iv_fi,ik_fi,is)= eqp%evshift(iv_fi,ik_fi,is)+tempshift*dweight/(tempsum)
                enddo
              enddo
            else if (xct%qflag.eq.0) then
              do is=1,xct%nspin
                do iv_fi=1,xct%nvb_fi
                  tempshift=0d0
                  tempsum=0d0
                  do iv_co=1,xct%nvb_co
                    tempshift = tempshift + (abs(dvn(iv_fi,iv_co,is,ikq_fi,ivert))**2)*&
                      eqp%evshift_co_q(iv_co,kgq_co%indr(ik_co_q),is)
                    tempsum = tempsum + (abs(dvn(iv_fi,iv_co,is,ikq_fi,ivert))**2)
                  enddo
                  do iv_co=xct%nvb_co+1,xct%n1b_co
                    tempshift = tempshift + (abs(dvn(iv_fi,iv_co,is,ikq_fi,ivert))**2)*&
                      eqp%ecshift_co_q(iv_co-xct%nvb_co,kgq_co%indr(ik_co_q),is)
                    tempsum = tempsum + (abs(dvn(iv_fi,iv_co,is,ikq_fi,ivert))**2)
                  enddo
                  if (.not.xct%renorm_transf) tempsum = 1d0
                  eqp%evshift(iv_fi,ik_fi,is)= eqp%evshift(iv_fi,ik_fi,is)+tempshift*dweight_q/(tempsum)
                enddo
              enddo
            else if (xct%qflag.eq.2) then
              do is=1,xct%nspin
                do iv_fi=1,xct%nvb_fi
                  tempshift=0d0
                  tempsum=0d0
                  do iv_co=1,xct%nvb_co
                    tempshift = tempshift + (abs(dvn(iv_fi,iv_co,is,ikq_fi,ivert))**2)*&
                      eqp%evshift_co(iv_co,kg_co%indr(ik_co_q),is)
                    tempsum = tempsum + (abs(dvn(iv_fi,iv_co,is,ikq_fi,ivert))**2)
                  enddo
                  do iv_co=xct%nvb_co+1,xct%n1b_co
                    tempshift = tempshift + (abs(dvn(iv_fi,iv_co,is,ikq_fi,ivert))**2)*&
                      eqp%ecshift_co(iv_co-xct%nvb_co,kg_co%indr(ik_co_q),is)
                    tempsum = tempsum + (abs(dvn(iv_fi,iv_co,is,ikq_fi,ivert))**2)
                  enddo
                  if (.not.xct%renorm_transf) tempsum = 1d0
                  eqp%evshift(iv_fi,ikq_fi,is)= eqp%evshift(iv_fi,ikq_fi,is)+tempshift*dweight_q/(tempsum)
                enddo
              enddo
            endif !qflag

            do is=1,xct%nspin
              do ic_fi=1,xct%ncb_fi
                tempshift=0d0
                tempsum=0d0
                if (xct%unrestricted_transf) then
                  do ic_co=1,xct%nvb_co
                    tempshift = tempshift + (abs(dcn(ic_fi,ic_co,is,ik_fi,ivert))**2)*&
                      eqp%evshift_co(ic_co,kg_co%indr(ik_co),is)
                    tempsum = tempsum + (abs(dcn(ic_fi,ic_co,is,ik_fi,ivert))**2)
                  enddo
                  do ic_co=xct%nvb_co+1,xct%n2b_co
                    tempshift = tempshift + (abs(dcn(ic_fi,ic_co,is,ik_fi,ivert))**2)*&
                      eqp%ecshift_co(ic_co-xct%nvb_co,kg_co%indr(ik_co),is)
                    tempsum = tempsum + (abs(dcn(ic_fi,ic_co,is,ik_fi,ivert))**2)
                  enddo
                else
                  do ic_co=1,xct%ncb_co
                    tempshift = tempshift + (abs(dcn(ic_fi,ic_co,is,ik_fi,ivert))**2)*&
                      eqp%ecshift_co(ic_co,kg_co%indr(ik_co),is)
                    tempsum = tempsum + (abs(dcn(ic_fi,ic_co,is,ik_fi,ivert))**2)
                  enddo
                endif
                if (.not.xct%renorm_transf) tempsum = 1d0
                eqp%ecshift(ic_fi,ik_fi,is)= eqp%ecshift(ic_fi,ik_fi,is)+ tempshift*dweight/(tempsum)
              enddo
            enddo
          endif
      
        enddo ! ijk
      enddo ! ii
      call progress_free(prog_info)
  
      if (xct%inteqp) then
        
        if (xct%nkpt_fi .ne. kg_fi%nf) then
          write(0,*) "WARNING: nkpt != nf", xct%nkpt_fi, kg_fi%nf
        endif
        
        call timing%start(timing%iw_write)
#ifdef MPI 
        SAFE_ALLOCATE(edummy, (xct%nvb_fi,xct%nkpt_fi,xct%nspin))
        edummy = eqp%evshift
        call MPI_Allreduce(edummy(1,1,1),eqp%evshift(1,1,1),xct%nkpt_fi*xct%nvb_fi* &
          xct%nspin,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr)
        SAFE_DEALLOCATE(edummy)
        
        SAFE_ALLOCATE(edummy, (xct%ncb_fi,xct%nkpt_fi,xct%nspin))
        edummy = eqp%ecshift
        call MPI_Allreduce(edummy(1,1,1),eqp%ecshift(1,1,1),xct%nkpt_fi*xct%ncb_fi* &
          xct%nspin,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr)
        SAFE_DEALLOCATE(edummy)
#endif
        
        if (peinf%inode .eq. 0) then
          
          ! Only write eqp_q.dat for velocity. For momentum, it is identical to valence bands of eqp.dat.
          if(flag%opr == 0 .or. xct%qflag.eq.0) then
            call open_file(unit=21,file='eqp_q.dat',form='formatted',status='replace')
            vunit = 21
            eqp_dat_size = xct%ncb_fi*xct%nspin
          else
            vunit = 22
            eqp_dat_size = (xct%nvb_fi+xct%ncb_fi)*xct%nspin
          endif
          call open_file(unit=22,file='eqp.dat',  form='formatted',status='replace')
          do ik_fi = 1, kg_fi%nf
            if(xct%qflag.eq.0 ) then
              write(21,'(3f13.9,i8)') kgq_fi%f(1,xct%indexq_fi(ik_fi)),kgq_fi%f(2,xct%indexq_fi(ik_fi)), &
                kgq_fi%f(3,xct%indexq_fi(ik_fi)),xct%nvb_fi*xct%nspin
            else if (flag%opr.eq.0) then
              write(21,'(3f13.9,i8)') kgq_fi%f(1,ik_fi),kgq_fi%f(2,ik_fi),kgq_fi%f(3,ik_fi),xct%nvb_fi*xct%nspin
            endif
            write(22,'(3f13.9,i8)') kg_fi%f(1,ik_fi),kg_fi%f(2,ik_fi),kg_fi%f(3,ik_fi),eqp_dat_size
            do is = 1, xct%nspin
              do iv_fi = xct%nvb_fi, 1, -1
                write(vunit,'(2i8,2f15.9)') is, xct%ifmaxq(ik_fi ,is) + 1 - iv_fi, eqp%evlda(iv_fi,ik_fi,is)*RYD, &
                  (eqp%evlda(iv_fi,ik_fi,is)+eqp%evshift(iv_fi,ik_fi,is))*RYD
              enddo
              do ic_fi = 1, xct%ncb_fi
                write(22,'(2i8,2f15.9)') is, xct%ifmax(ik_fi ,is) + ic_fi, eqp%eclda(ic_fi,ik_fi,is)*RYD, &
                  (eqp%eclda(ic_fi,ik_fi,is)+eqp%ecshift(ic_fi,ik_fi,is))*RYD
              enddo
            enddo
          enddo
          if(flag%opr == 0 .or. xct%qflag.eq.0) call close_file(21)
          call close_file(22)
          
          call open_file(unit=23,file='bandstructure.dat',form='formatted',status='replace')
          ! Write k-points in cart. coord. instead of cryst. coord. for the 
          ! band structure plot (as it is done in espresso and wannier90) --GSM
          write(23,'(a)') '# spin    band          kx          ky          kz          E(MF)          E(QP)        Delta E'
          write(23,'(a)') '#                        (Cartesian coordinates)             (eV)           (eV)           (eV)'

          do is = 1, xct%nspin
            do iv_fi = xct%nvb_fi, 1, -1
              ! For the valence bands, the "q" structures are the ones that contained the calculated info. NOTE: kgq%fi = kg%fi
              do ik_fi = 1,kgq_fi%nf
                write(23,'(i6,i8,3f12.5,3f15.9)') is,xct%ifmaxq(ik_fi,is) + 1 - iv_fi, &
                  matmul(crys%bvec(1:3,1:3),kgq_fi%f(1:3,ik_fi))*crys%blat, &
                  eqp%evlda(iv_fi,ik_fi,is)*RYD,(eqp%evlda(iv_fi,ik_fi,is)+eqp%evshift(iv_fi,ik_fi,is))*RYD, &
                  eqp%evshift(iv_fi,ik_fi,is)*RYD
              enddo
            enddo
            do ic_fi = 1,xct%ncb_fi
              do ik_fi = 1,kg_fi%nf
                write(23,'(i6,i8,3f12.5,3f15.9)') is,xct%ifmax(ik_fi,is) + ic_fi, &
                  matmul(crys%bvec(1:3,1:3),kg_fi%f(1:3,ik_fi))*crys%blat, &
                  eqp%eclda(ic_fi,ik_fi,is)*RYD,(eqp%eclda(ic_fi,ik_fi,is)+eqp%ecshift(ic_fi,ik_fi,is))*RYD, &
                  eqp%ecshift(ic_fi,ik_fi,is)*RYD
              enddo
            enddo
          enddo
          call close_file(23)

! JRD: We now update our fine-grid eigenvalues with the interpolated ones.
  
        endif ! peinf%inode .eq. 0

        call timing%stop(timing%iw_write)
        eqp%ecqp(:,:,:)=eqp%eclda(:,:,:)+eqp%ecshift(:,:,:)
        eqp%evqp(:,:,:)=eqp%evlda(:,:,:)+eqp%evshift(:,:,:)    
        
      endif ! xct%inteqp
      
      call dealloc_grid(kg_co)
      
      !if (xct%qflag .eq. 0) call dealloc_grid(kg_co_q)
      
      ! typedefs initializes all of these ikolds to 0
      if (work%ikold .ne. 0) then
        SAFE_DEALLOCATE_P(work%cg)
        SAFE_DEALLOCATE_P(work%ph)
        SAFE_DEALLOCATE_P(work%ind)
        SAFE_DEALLOCATE_P(work%isort)
      endif
      if (workq%ikold .ne. 0) then
        SAFE_DEALLOCATE_P(workq%cg)
        SAFE_DEALLOCATE_P(workq%ph)
        SAFE_DEALLOCATE_P(workq%ind)
        SAFE_DEALLOCATE_P(workq%isort)
      endif
      if (workco%ikold .ne. 0) then
        SAFE_DEALLOCATE_P(workco%cg)
        SAFE_DEALLOCATE_P(workco%ph)
        SAFE_DEALLOCATE_P(workco%ind)
        SAFE_DEALLOCATE_P(workco%isort)
      endif
      if (workcoq%ikold .ne. 0) then
        SAFE_DEALLOCATE_P(workcoq%cg)
        SAFE_DEALLOCATE_P(workcoq%ph)
        SAFE_DEALLOCATE_P(workcoq%ind)
        SAFE_DEALLOCATE_P(workcoq%isort)
      endif

!--------------------------------
! Write out dcn,dvn. All PEs access dtmat
! FHJ: TODO - only one PE should do the job!
      
      do ivert = 1, xct%npts_intp_kernel
      do ii=1,peinf%npes
        if(peinf%inode+1.eq.ii) then
          call open_file(unit=13,file='dtmat',position='append',form='unformatted',status='old')
          do ikt=1,peinf%ikt(peinf%inode+1)
            ik_fi=peinf%ik(peinf%inode+1,ikt)
            do is=1,xct%nspin
              do ic_co=1,xct%n2b_co
                do ic_fi=1,xct%ncb_fi
                  write(13) ik_fi,ic_fi,fi2co_wfn(ivert,ik_fi),ic_co,is, &
                    dcn(ic_fi,ic_co,is,ik_fi,ivert)
                enddo
              enddo
            enddo
          enddo
          call close_file(13)
        endif
#ifdef MPI
        call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

      enddo
      enddo
  
      do ivert = 1, xct%npts_intp_kernel
      do ii=1,peinf%npes
        if(peinf%inode+1.eq.ii) then
          call open_file(unit=13,file='dtmat',position='append',form='unformatted',status='old')
          do ikt=1,peinf%ikt(peinf%inode+1)
            ik_fi=peinf%ik(peinf%inode+1,ikt)
            do is=1,xct%nspin
              do iv_fi=1,xct%nvb_fi
                do iv_co=1,xct%n1b_co
                  write(13) ik_fi,iv_fi,fi2co_wfn(ivert,ik_fi),iv_co,is, &
                    dvn(iv_fi,iv_co,is,ik_fi,ivert)
                enddo
              enddo
            enddo
          enddo
          call close_file(13)
        endif
#ifdef MPI
        call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

      enddo
      enddo
      
      SAFE_DEALLOCATE_P(distgwfco%isort)
      if (xct%qflag.ne.0) then
        SAFE_DEALLOCATE_P(distgwfco%zv)
      endif
      SAFE_DEALLOCATE_P(distgwfco%zc)
      if (xct%qflag .eq. 0) then
        SAFE_DEALLOCATE_P(distgwfcoq%isort)
        SAFE_DEALLOCATE_P(distgwfcoq%zv)
      endif
      
!--------------------------------
! Finished the computation of transformation matrices, dcn/dvn
! Now, PEs share the information

      call timing%start(timing%iw_reduce)
#ifdef MPI
      SAFE_ALLOCATE(dummy, (xct%ncb_fi,xct%n2b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
      dummy = dcn
      call MPI_Allreduce(dummy(1,1,1,1,1),dcn(1,1,1,1,1),xct%nkpt_fi*xct%ncb_fi* &
        xct%n2b_co*xct%nspin*xct%npts_intp_kernel,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
      SAFE_DEALLOCATE(dummy)
      SAFE_ALLOCATE(dummy, (xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
      dummy = dvn
      call MPI_Allreduce(dummy(1,1,1,1,1),dvn(1,1,1,1,1),xct%nkpt_fi*xct%nvb_fi* &
        xct%n1b_co*xct%nspin*xct%npts_intp_kernel,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
      SAFE_DEALLOCATE(dummy)
      SAFE_ALLOCATE(idum, (xct%npts_intp_kernel, xct%nkpt_fi))
      idum = fi2co_wfn
      call MPI_Allreduce(idum(1,1),fi2co_wfn(1,1),xct%nkpt_fi*xct%npts_intp_kernel,MPI_INTEGER, &
        MPI_SUM,MPI_COMM_WORLD,mpierr)
      SAFE_DEALLOCATE(idum)
      SAFE_ALLOCATE(dummy_coefs, (xct%npts_intp_kernel, xct%nkpt_fi))
      dummy_coefs = intp_coefs
      call MPI_Allreduce(dummy_coefs(1,1), intp_coefs(1,1), xct%nkpt_fi*xct%npts_intp_kernel, MPI_REAL_DP, &
        MPI_SUM, MPI_COMM_WORLD, mpierr)
      SAFE_DEALLOCATE(dummy_coefs)
      !if (xct%qflag.ne.1) then
      !  SAFE_ALLOCATE(dummy_coefs, (xct%npts_intp_kernel, xct%nkpt_fi))
      !  dummy_coefs = intp_coefs_q
      !  call MPI_Allreduce(dummy_coefs(1,1), intp_coefs_q(1,1), xct%nkpt_fi*xct%npts_intp_kernel, MPI_REAL_DP, &
      !    MPI_SUM, MPI_COMM_WORLD, mpierr)
      !  SAFE_DEALLOCATE(dummy_coefs)
      !endif
#endif
      call timing%stop(timing%iw_reduce)

      ! I am not sure why but this barrier is necessary to avoid hanging. --DAS
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

      if (xct%npts_intp_kernel>1 .and. peinf%inode==0) then
        call open_file(unit=13,file='dtmat',position='append',form='unformatted',status='old')
        write(13) intp_coefs(:,:)
        call close_file(13)
        !if (xct%qflag.ne.1) then
        !  call open_file(unit=13,file='dtmat_q',position='append',form='unformatted',status='old')
        !  write(13) intp_coefs_q(:,:)
        !  call close_file(13)
        !endif
      endif

#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

      
      if (xct%delaunay_interp) then
        call interp_free(wfn_interp)
        if (xct%qflag.eq.0) then
          !call interp_free(wfn_interp_q)
          call dealloc_intpts()
        endif
      else
        call dealloc_intpts()
      endif
    
    endif ! xct%skipinterp

  endif ! flag%read_dtmat

!--------------------------------
! Check if transformation matrices are unitary and print out their
! norm (i.e., the norm of wavefunctions in the fine grid)
! The norm should be as close as possible to 1.d0. If it is not,
! important overlaps are being neglected.
  do ivert = 1, xct%npts_intp_kernel
  do is = 1, xct%nspin
    do ik_fi=1,xct%nkpt_fi
      do iv=1,xct%nvb_fi
        dnorm = 0.d0
        do iv_co=1,xct%n1b_co
          dnorm = dnorm + abs(dvn(iv,iv_co,is,ik_fi,ivert))**2
        enddo
        normv(ik_fi,iv,is,ivert) = dnorm
      enddo
      
      do ic=1,xct%ncb_fi
        dnorm = 0.d0
        do ic_co=1,xct%n2b_co
          dnorm = dnorm + abs(dcn(ic,ic_co,is,ik_fi,ivert))**2
        enddo
        normc(ik_fi,ic,is,ivert) = dnorm
      enddo
    enddo
  enddo
  enddo
      
  if (.not.flag%read_dtmat .and. peinf%inode==0) then
    call open_file(20,file='dcmat_norm.dat',form='formatted',status='replace')
    write(20,'(a,i12,a)') ' -------  Norm of dcc matrices : Spins = ', xct%nspin,'  -------'
    if (xct%nspin.eq.1) then
      write(20,'(a)') '           k-point           ik_co    c   dist    |dcc|^2      '
    else
      write(20,'(a)') '           k-point           ik_co    c   dist    |dcc|^2 (spin up)  |dcc|^2 (spin down) '
    endif
    write(20,'(21a3)') ('---',kk=1,21)
    do ik_fi=1,xct%nkpt_fi
      ik_co= fi2co_wfn(1,ik_fi)
      do kk=1,3
        qg(kk)=kg_fi%f(kk,ik_fi)-kco(kk,ik_co)
        qg(kk) = qg(kk) - anint( qg(kk) )
      enddo
      dist= DOT_PRODUCT(qg,MATMUL(crys%bdot,qg))
      do jj=1,xct%ncb_fi
        write(20,250) kg_fi%f(:,ik_fi),ik_co,jj,sqrt(dist),(normc(ik_fi,jj,is,1), is = 1, xct%nspin)
      enddo
    enddo
    call close_file(20)
    
    call open_file(20,file='dvmat_norm.dat',form='formatted',status='replace')
    write(20,'(a,i12,a)') ' -------  Norm of dvv matrices : Spins = ', xct%nspin,'  -------'
    if (xct%nspin.eq.1) then
      write(20,'(a)') '           k-point           ik_co    v   dist    |dvv|^2      '
    else
      write(20,'(a)') '           k-point           ik_co    v   dist    |dvv|^2 (spin up)  |dvv|^2 (spin down) '
    endif
    write(20,'(21a3)') ('---',kk=1,21)
    do ik_fi=1,xct%nkpt_fi
      ik_co= fi2co_wfn(1,ik_fi)
      do kk=1,3
        qg(kk)=kg_fi%f(kk,ik_fi)-kco(kk,ik_co)
        qg(kk) = qg(kk) - anint( qg(kk) )
      enddo
      dist= DOT_PRODUCT(qg,MATMUL(crys%bdot,qg))
      do jj=1,xct%nvb_fi
        write(20,250) kg_fi%f(:,ik_fi),ik_co,jj,sqrt(dist),(normv(ik_fi,jj,is,1), is = 1, xct%nspin)
      enddo
    enddo
    call close_file(20)
250   format(' ( ',f6.3,' , ',f6.3,' , ',f6.3, ' ) ',i4,1x,i4,1x,f6.3,1x,f10.6,9x,f10.6)
  endif    

!---------------------------------
! Renormalize transformation matrices: matrices dcn,dvn are
! renormalized to 1.d0 (ideally, this step should not be
! performed but the interpolation can be really bad if dcn,dvn
! do not have unitary norm).

  ! FHJ: Print worst norm before restriction/renormalization
  worst_normv = minval(normv)
  worst_normc = minval(normc)
  if (peinf%inode==0) then    
    write(6,'(1x,a)') 'Max. error in norm of transformation coefficients (1 - \sum_co |d_fi,co|^2):'
    write(6,'(1x,a,es9.3e2)') ' - For valence states: ', 1d0 - worst_normv
    write(6,'(1x,a,es9.3e2)') ' - For conduction states: ', 1d0 - worst_normc
    call warn_norm(xct, worst_normv, worst_normc)
  endif
  ! FHJ: truncate or zero out coefficients. User doesn`t care about the details...
  if (truncate_coefs) then
    if (peinf%inode==0) write(6,'(1x,a)') "Restricting transformation to dvv',dcc' subspaces."
    dvv(:,:,:,:,:) = dvn(:,:xct%nvb_co,:,:,:)
    dcc(:,:,:,:,:) = dcn(:,xct%nvb_co+1:,:,:,:)
  elseif (xct%zero_unrestricted_contrib) then
    if (peinf%inode==0) write(6,'(1x,a)') "Restricting transformation to dvv',dcc' subspaces."
    dvn(:,xct%nvb_co+1:,:,:,:) = ZERO
    dcn(:,:xct%nvb_co,:,:,:) = ZERO
  endif

  ! FHJ: And from now on, xct%n1b_co will depend only on xct%extended_kernel.
  if (.not.xct%extended_kernel) then
    xct%n1b_co = xct%nvb_co
    xct%n2b_co = xct%ncb_co
  else
    xct%n1b_co = xct%nvb_co + xct%ncb_co
    xct%n2b_co = xct%n1b_co
  endif

  ! FHJ: inteqp also renormalizes the coefs for consistency (and someone might
  ! need this in the future). And note that we go here if truncate_coefs==.true.
  ! only to recalculate the worst norm.
  if (truncate_coefs.or.xct%renorm_transf) then
    worst_normc = INF
    worst_normv = INF
    do ivert = 1, xct%npts_intp_kernel
    do is=1, xct%nspin
      do ik=1,xct%nkpt_fi
        do ic_fi=1,xct%ncb_fi
          dnorm = sum(abs(dcc(ic_fi,:,is,ik,ivert))**2)
          worst_normc = min(worst_normc, dnorm)
          if (xct%renorm_transf) then
            dnorm = sqrt(dnorm)
            dcc(ic_fi,:,is,ik,ivert) = dcc(ic_fi,:,is,ik,ivert)/dnorm
          endif
        enddo
        do iv_fi=1,xct%nvb_fi
          dnorm = sum(abs(dvv(iv_fi,:,is,ik,ivert))**2)
          worst_normv = min(worst_normv, dnorm)
          if (xct%renorm_transf) then
            dnorm = sqrt(dnorm)
            dvv(iv_fi,:,is,ik,ivert) = dvv(iv_fi,:,is,ik,ivert)/dnorm
          endif
        enddo
      enddo
    enddo
    enddo
  endif
  if (peinf%inode==0.and.(truncate_coefs.or.xct%zero_unrestricted_contrib)) then
    write(6,'(1x,a)') 'Max. error after restricting WFN transformation:'
    write(6,'(1x,a,es9.3e2)') ' - For valence states: ', 1d0 - worst_normv
    write(6,'(1x,a,es9.3e2)') ' - For conduction states: ', 1d0 - worst_normc
    call warn_norm(xct, worst_normv, worst_normc)
  endif

!-------------------------------------
! Print out information about the number of k-points in the fine grid
! interpolated from each point in the coarse grid
! ii : number of fine-grid points close to point ik_co

  if (peinf%verb_high .and. peinf%inode==0) then
    write(6,'(/1x,a)') 'K-points in coarse grid:'
    write(6,'(3x,a6,3(1x,a10),1x,a6)') 'ik', 'kx', 'ky', 'kz', 'Nfi'
    do ik_co=1,xct%nkpt_co
      ii = 0
      do ik=1,xct%nkpt_fi
        if (any(fi2co_wfn(:,ik)==ik_co)) ii = ii + 1
      enddo
      write(6,'(3x,i6,3(1x,f10.6),1x,i6)') ik_co, kco(:,ik_co), ii
    enddo
    write(6,'()')
  endif

  if (.not.associated(dvn,dvv)) then
    SAFE_DEALLOCATE_P(dvn)
  endif
  if (.not.associated(dcn,dcc)) then
    SAFE_DEALLOCATE_P(dcn)
  endif
  SAFE_DEALLOCATE(normc)
  SAFE_DEALLOCATE(normv)

  POP_SUB(intwfn)
  
  return
end subroutine intwfn

subroutine intwfn_sub(kp,crys,syms,xct,flag,gvec,kg_fi,kgq_fi, &
  dcc_sub,dvv_sub,kg_sub_co,kg_co,indexq_fi,nk_sub,intwfnv,intwfnc,closepts)
  type (kpoints), intent(inout) :: kp
  type (crystal), intent(in) :: crys
  type (symmetry), intent(in) :: syms
  type (xctinfo), intent(inout) :: xct
  type (flags), intent(in) :: flag
  type (gspace), intent(in) :: gvec
  type (grid), intent(in) :: kg_fi,kgq_fi
  !> (xct%ncb_fi,xct%n2b_co,xct%nspin,xct%nkpt_fi,4) Interp coefs.
  SCALAR, intent(out), target :: dcc_sub(:,:,:,:,:,:)
  !> (xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,4) Interp coefs.
  SCALAR, intent(out), target :: dvv_sub(:,:,:,:,:,:)
  type (grid), intent(inout) :: kg_sub_co,kg_co
  integer, intent(in) :: indexq_fi(:) !< xct%nkpt_fi
  integer, intent(in) :: nk_sub
  type (int_wavefunction), intent(inout) :: intwfnc
  type (int_wavefunction), intent(inout):: intwfnv
  integer, intent(out) :: closepts(:,:)

  type (grid) :: kg_sub,kgq_sub
  type (kpoints) :: kp_sub,kpq_sub
  type (crystal) :: crys_sub
  type (symmetry) :: syms_sub
  type (gspace) gvec_sub,gvec_kpt,gvec_kpt_q
  type (wavefunction) :: wfnc,wfnc_fi,wfnv,wfnvq_fi
  type(interp_t) :: wfn_interp

  character*80, allocatable :: fname(:),fnameq(:)
  character(len=3) :: sheader
  integer :: ik,ik_sub,ik_fi,ikq_fi,ig,ib,ik_co,ik_loc,iunit,ik_co_q,vunit
  integer :: ic,ic_co,ic_fi,iv,iv_co,iv_fi,is,isp,kk,jj,hh,ii
  integer :: ivert,iflavor,ngkmax
  integer :: gumk(3),igumk,igumkq,save_npes,save_inode
  integer, allocatable :: dummy_closepts(:,:),isorti(:)
  real(DP) :: closeweights(4),qg(3),qk(3)
  real(DP), allocatable :: normc(:,:,:,:,:),normv(:,:,:,:,:),ekin(:)
  real(DP) :: worst_normc, worst_normv,dnorm
  SCALAR, allocatable :: dummy(:,:,:,:,:,:)
  SCALAR, allocatable :: cg(:,:), cgarray(:),cgq(:,:)
  type(work_genwf) :: work, workq,work_co,workq_co
  !> (xct%nvb_fi, xct%n1b, xct%nspin, xct%nkpt_fi, 4, nk_sub)
  !> dvn points to dvv if we are restricting the transf right away,
  !! or if we never restrict it (i.e., extended_kernel==.true.). Otherwise,
  !! we allocate a buffer, calculate the unrestricted transf, and restrict+copy
  ! dvn to dvv later on.
  SCALAR, pointer :: dvn(:,:,:,:,:,:) 
  SCALAR, pointer :: dcn(:,:,:,:,:,:) !< See: dvn
  !> .true. if we calculate the unrestricted transformation, but the kernel is
  !! restricted (i.e., we`ll have to truncate the dvn to dvv later on).
  logical :: truncate_coefs,skip_checkbz
  type(progress_info) :: prog_info

  PUSH_SUB(intwfn_sub)

  ! FHJ: Normally, xct%n1b_co depends only on xct%extended_kernel. However,
  ! we`ll temporarily overwrite this value depending on xct%unrestricted_transf
  ! till the end of the subroutine, because if xct%unrestricted_transf==.true,
  ! intwfn "thinks" it`s dealing with an extended kernel when calculating dvn/dcn.
  if (xct%unrestricted_transf) then
    xct%n1b_co = xct%nvb_co + xct%ncb_co
    xct%n2b_co = xct%n1b_co
  else
    xct%n1b_co = xct%nvb_co
    xct%n2b_co = xct%ncb_co
  endif
  truncate_coefs = xct%unrestricted_transf.and.(.not.xct%extended_kernel)
  if (truncate_coefs) then
    ! FHJ: we`ll restrict the coefficients later on, but we want to do the
    ! right expansion here => need a temp buffer.
    SAFE_ALLOCATE(dvn, (xct%nvb_fi, xct%n1b_co, xct%nspin, xct%nkpt_fi,xct%idimensions+1, nk_sub))
    SAFE_ALLOCATE(dcn, (xct%ncb_fi, xct%n2b_co, xct%nspin, xct%nkpt_fi,xct%idimensions+1, nk_sub))
    dvn = 0.d0
    dcn = 0.d0
  else
    ! FHJ: we either do not restrict, or restrict right from the beginning.
    dvn => dvv_sub
    dcn => dcc_sub
  endif
  dcc_sub = 0.d0
  dvv_sub = 0.d0

  SAFE_ALLOCATE(normc, (xct%nkpt_fi,xct%ncb_fi,xct%nspin,xct%idimensions+1,nk_sub))
  SAFE_ALLOCATE(normv, (xct%nkpt_fi,xct%nvb_fi,xct%nspin,xct%idimensions+1,nk_sub))

  ! Move points in kg_sub_co to Wigner Seitz cell
  call fullbz(crys,syms,kg_sub_co,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  ! Use delaunay interpolation to find which points on the coarse grid are closest to each fine point
  call interp_init(wfn_interp, crys, kg_sub_co%f, dims=xct%idimensions, &
    periodic=.true., active_dims=xct%is_periodic)

  ! Read names of subsampled wavefunctions
  SAFE_ALLOCATE(fname,(kg_sub_co%nr))
  if (xct%qflag.eq.0) then
    SAFE_ALLOCATE(fnameq,(kg_sub_co%nr))
  endif
  call open_file(unit=77, file='subsample.inp',form='formatted',status='old')
  read(77,*) !nk_sub
  read(77,*) !nsub_files
  do ik_sub=1,kg_sub_co%nr*2 ! read and ignore list of k-points and bsemat files
    read(77,*)
  enddo
  do ik_sub=1,kg_sub_co%nr
    read(77,'(a80)') fname(ik_sub)
    if (peinf%inode==0.and.peinf%verb_medium) write(6,*) "Read fname ",trim(fname(ik_sub))
  enddo
  if (xct%qflag.eq.0) then
    do ik_sub=1,kg_sub_co%nr
      read(77,'(a80)') fnameq(ik_sub)
    enddo
  endif
  call close_file(77)

  ! Save the number of processes, so we can restore after setting to 1 later on
  save_npes=peinf%npes
  save_inode=peinf%inode
  closepts = 0

  ! allocating some arrays
  SAFE_ALLOCATE(gvec_sub%components, (3, gvec%ng))
  SAFE_ALLOCATE(wfnv%isort, (gvec%ng))
  SAFE_ALLOCATE(wfnc%isort, (gvec%ng))
  SAFE_ALLOCATE(work_co%isort, (gvec%ng))
  SAFE_ALLOCATE(work_co%ind,(gvec%ng))
  SAFE_ALLOCATE(work_co%ph,(gvec%ng))
  if (xct%qflag.eq.0)then
    SAFE_ALLOCATE(workq_co%isort, (gvec%ng))
    SAFE_ALLOCATE(workq_co%ind,(gvec%ng))
    SAFE_ALLOCATE(workq_co%ph,(gvec%ng))
  endif
  SAFE_ALLOCATE(isorti, (gvec%ng))
  SAFE_ALLOCATE(ekin, (gvec%ng))
          

  ! DYQ - TODO: add ability to read precalculated transformation matrices
  if (flag%read_dtmat_sub) then
    call die('Reading of precomputed dtmat matrices not implemented.')
  else

!-------------------------------
! Read subsampled wavefunctions

   call progress_init(prog_info, 'subsampled wavefunction interpolation', 'k-point', &
    peinf%nkpe)
    do ik_loc=1,peinf%nkpe ! Outermost loop over fine k-points
       call progress_step(prog_info)
      if (peinf%inode==0.and.peinf%verb_high) write(6,'(a,i6,a,i6)') "Proc 0: working on kpt",ik_loc,"/",peinf%nkpe
      ! Each proc calculates dvn, dcn, closepts for the fine k-points it owns
      ik_fi = peinf%ik(peinf%inode+1,ik_loc)
      if (ik_fi.eq.0) cycle
      if (xct%qflag.eq.1) then
        ikq_fi = indexq_fi(ik_fi)
      else
        ikq_fi = xct%indexq_fi(ik_fi)
      endif

      ! Find which points on the coarse grid are closest to each fine point
      call interp_eval(wfn_interp, kg_fi%f(:,ik_fi), closepts(:,ik_fi), closeweights(:))

      ! Read the current fine wavefunctions
      call timing%start(timing%iw_genwf)
      intwfnc%nspinor=kp%nspinor
      intwfnv%nspinor=kp%nspinor
      call genwf(crys,gvec,kg_fi,syms,wfnc_fi,ik_fi,ik_fi,xct%nspin,xct%ncb_fi,&
                 work,intwfnc,xct%iwriteint,is_cond=.true.)
      call genwf(crys,gvec,kgq_fi,syms,wfnvq_fi,ik_fi,ikq_fi,xct%nspin,xct%nvb_fi,&
                 workq,intwfnv,xct%iwriteint,is_cond=.false.)
      call timing%stop(timing%iw_genwf)

      do ivert=1,xct%idimensions+1 ! loop over closepts (2 points for 1D, 3 points for 2D)
        ik_co = closepts(ivert,ik_fi)
        ! Each proc reads the subsampled wavefunction if the coarse point is close to a fine point owned by the proc
        if (peinf%inode==0.and.peinf%verb_high) write(6,'(a,a)') "Reading CSI wavefunction ",trim(fname(ik_co))
        iunit=211+ik_fi ! each process has a different unit number
        ! Open and read WFN_sub header
        call open_file(unit=iunit,file=trim(fname(ik_co)),form='unformatted',status='old')
        sheader = 'WFN'
        iflavor = 0
        ! To avoid broadcasting anything while reading the subsampled wavefunctions, temporarily set number of processes to 1
        ! DYQ - TODO: Do something less stupid here
        peinf%npes=1
        peinf%inode=0
        call read_binary_header_type(iunit, sheader, iflavor, kp_sub, gvec_sub, syms_sub, crys_sub, warn = .false., &
          dont_warn_kgrid=.True.)
        peinf%inode=save_inode
        peinf%npes=save_npes
        call check_trunc_kpts(xct%icutv, kp_sub) ! this is only done by proc. 0
        call check_header('WFN_fi', kp, gvec, syms, crys, trim(fname(ik_co)), kp_sub, gvec_sub, syms_sub,&
          crys_sub, is_wfn = .true.)
        ! generate kgrid and move to 1st bz
        kg_sub%nr=kp_sub%nrk
        SAFE_ALLOCATE(kg_sub%r,(3,kg_sub%nr))
        kg_sub%r(1:3,1:kg_sub%nr)=kp_sub%rk(1:3,1:kp_sub%nrk)
        call fullbz(crys_sub,syms_sub,kg_sub,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
        
        ! Read global gvecs
        peinf%inode=0
        peinf%npes=1
        SAFE_ALLOCATE(gvec_sub%components, (3, gvec_sub%ng))
        call read_binary_gvectors(iunit, gvec%ng, gvec%ng, gvec_sub%components,bcast=.false.)
        peinf%inode=save_inode
        peinf%npes=save_npes

        if(any(kp_sub%ifmax(:,:) == 0)) &
          call die("BSE codes cannot handle a system where some k-points have no occupied bands.")
        kp_sub%nvband=minval(kp_sub%ifmax(:,:)-kp_sub%ifmin(:,:))+1
        kp_sub%ncband=kp_sub%mnband-maxval(kp_sub%ifmax(:,:))
        
        ! Check whether requested number of bands is available
        if(xct%nvb_co .gt. kp_sub%nvband) then
          call die("The requested number of valence bands is not available in WFN_sub.")
        endif
        if(xct%ncb_co .gt. kp_sub%ncband) then
          call die("The requested number of conduction bands is not available in WFN_sub.")
        endif

        ! For finite q, read valence wavefunctions from WFNq_sub
        if (xct%qflag.eq.0) then
          vunit=311+ik_fi
          if (peinf%inode==0.and.peinf%verb_medium) write(6,'(a,a)') "Reading CSI_q wavefunction ", trim(fnameq(ik_co))
          call open_file(unit=vunit,file=trim(fnameq(ik_co)),form='unformatted',status='old')
          sheader = 'WFN'
          iflavor = 0
          ! To avoid broadcasting anything while reading the subsampled wavefunctions, temporarily set number of processes to 1
          ! DYQ - TODO: Do something less stupid here
          peinf%npes=1
          peinf%inode=0
          call read_binary_header_type(vunit, sheader, iflavor, kpq_sub, gvec_sub, syms_sub, crys_sub, warn = .false., &
            dont_warn_kgrid=.True.)
          call read_binary_gvectors(vunit, gvec%ng, gvec%ng, gvec_sub%components,bcast=.false.)
          peinf%inode=save_inode
          peinf%npes=save_npes
          if(any(kpq_sub%ifmax(:,:) == 0)) &
            call die("BSE codes cannot handle a system where some k-points have no occupied bands.")
          kpq_sub%nvband=minval(kp_sub%ifmax(:,:)-kp_sub%ifmin(:,:))+1
          kpq_sub%ncband=kp_sub%mnband-maxval(kp_sub%ifmax(:,:))
        
          ! Check whether requested number of bands is available
          if(xct%nvb_co .gt. kpq_sub%nvband) then
            call die("The requested number of valence bands is not available in WFNq_sub.")
          endif
          ! generate kgrid and move to 1st bz
          kgq_sub%nr=kpq_sub%nrk
          SAFE_ALLOCATE(kgq_sub%r,(3,kgq_sub%nr))
          kgq_sub%r(1:3,1:kgq_sub%nr)=kpq_sub%rk(1:3,1:kpq_sub%nrk)
          call fullbz(crys_sub,syms_sub,kgq_sub,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
          if (peinf%inode==0.and.peinf%verb_high) write(6,'(a,a)') "Read header of ",trim(fnameq(ik_co))
        endif ! qflag.eq.0

        wfnv%nband=xct%nvb_co
        wfnv%nspin=kp_sub%nspin
        wfnv%nspinor=kp_sub%nspinor
        wfnc%nband=xct%ncb_co
        wfnc%nspin=kp_sub%nspin
        wfnc%nspinor=kp_sub%nspinor

!------------------------------------------------------------------
! Loop over the subsampled k-points
!
        do ik_sub=1,nk_sub

          ! Read and sort gvecs at each k-point
          SAFE_ALLOCATE(gvec_kpt%components, (3, kp_sub%ngk(ik_sub)))
          peinf%inode=0
          peinf%npes=1
          call read_binary_gvectors(iunit, kp_sub%ngk(ik_sub), kp_sub%ngk(ik_sub), gvec_kpt%components,bcast=.false.)
          if (xct%qflag.eq.0) then
            SAFE_ALLOCATE(gvec_kpt_q%components, (3, kpq_sub%ngk(ik_sub)))
            call read_binary_gvectors(vunit, kpq_sub%ngk(ik_sub), kpq_sub%ngk(ik_sub), gvec_kpt_q%components,bcast=.false.)
          endif
          peinf%inode=save_inode
          peinf%npes=save_npes
          SAFE_ALLOCATE(cg, (kp_sub%ngk(ik_sub),kp_sub%nspin*kp_sub%nspinor))
          if (xct%qflag.eq.0)then
            SAFE_ALLOCATE(cgq, (kpq_sub%ngk(ik_sub),kp_sub%nspin*kp_sub%nspinor))
          endif

          ! allocate arrays for wfnv and wfnc (This is different on each proc.)
          if (xct%qflag.ne.0) then
            wfnv%ng=kp_sub%ngk(ik_sub)
            SAFE_ALLOCATE(cgarray, (kp_sub%ngk(ik_sub)))
          else
            wfnv%ng=kpq_sub%ngk(ik_sub)
            ngkmax=max(kp_sub%ngk(ik_sub),kpq_sub%ngk(ik_sub))
            SAFE_ALLOCATE(cgarray, (ngkmax))
            workq_co%ng=kpq_sub%ngk(ik_sub)
            SAFE_ALLOCATE(workq_co%cg, (wfnv%ng,wfnv%nband,wfnv%nspin*wfnv%nspinor))
          endif
          wfnc%ng=kp_sub%ngk(ik_sub)
          work_co%ng=kp_sub%ngk(ik_sub)
          SAFE_ALLOCATE(wfnv%cg, (wfnv%ng,wfnv%nband,wfnv%nspin*wfnv%nspinor))
          SAFE_ALLOCATE(wfnc%cg, (wfnc%ng,wfnc%nband,wfnc%nspin*wfnc%nspinor))
          SAFE_ALLOCATE(work_co%cg, (wfnc%ng,wfnc%nband,wfnc%nspin*wfnc%nspinor))

          ! loop over bands
          do ib=1,kp_sub%mnband
            peinf%inode=0
            peinf%npes=1
            call read_binary_data(iunit, kp_sub%ngk(ik_sub), kp_sub%ngk(ik_sub), kp_sub%nspin*kp_sub%nspinor, cg, bcast=.false.)
            if ((ib.le.kpq_sub%mnband) .and. (xct%qflag.eq.0)) then
              call read_binary_data(vunit,kpq_sub%ngk(ik_sub),kpq_sub%ngk(ik_sub),kp_sub%nspin*kp_sub%nspinor,cgq,bcast=.false.)
            endif
            peinf%inode=save_inode
            peinf%npes=save_npes
            do is=1, kp_sub%nspin
              if (ib .gt. (kp_sub%ifmax(ik_sub,is)-xct%nvb_co) .and. ib .le. (kp_sub%ifmax(ik_sub,is)+xct%ncb_co)) then
                do isp=1, kp_sub%nspinor
                  ! valence bands
                  if ((ib.le.kp_sub%ifmax(ik_sub,is)).and. &
                    (ib.gt.kp_sub%ifmax(ik_sub,is)-xct%nvb_co)) then
                    do ig=1, wfnv%ng
                      if (xct%qflag.ne.0) then
                        cgarray(ig)=cg(ig, is*isp)
                      else
                        cgarray(ig)=cgq(ig, is*isp)
                      endif
                    enddo
                    wfnv%cg(1:wfnv%ng,kp_sub%ifmax(ik_sub,is)-ib+1,is*isp)=cgarray(1:wfnv%ng)
                  endif
                  ! conduction bands
                  if ((ib.gt.kp_sub%ifmax(ik_sub,is)).and. &
                    (ib.le.kp_sub%ifmax(ik_sub,is)+xct%ncb_co)) then
                    do ig=1, wfnc%ng
                      cgarray(ig)=cg(ig, is*isp)
                    enddo
                    wfnc%cg(1:wfnc%ng,ib-kp_sub%ifmax(ik_sub,is),is*isp)=cgarray(1:wfnc%ng)
                  endif
                enddo ! isp
                if (xct%qflag.ne.0) then
                  call checknorm('WFN_sub',ib,ik_sub,kp_sub%ngk(ik_sub),is,kp_sub%nspinor,cg(:,:))
                endif
              endif
            enddo ! is
          enddo !ib

          SAFE_DEALLOCATE(cgarray)
          SAFE_DEALLOCATE(cg)
          if (xct%qflag.eq.0) then
            SAFE_DEALLOCATE(cgq)
          endif

          ! Fixing phases
          work_co%ind = 0
          work_co%ph = ZERO
          if (xct%qflag.eq.0)then
            workq_co%ind = 0
            workq_co%ph = ZERO
          endif
          isorti(:)=0
          ekin = 0.0d0
          ! isort stored mapping between local and global gvecs
          do ig=1,work_co%ng
            call findvector(work_co%isort(ig), gvec_kpt%components(:, ig), gvec)
            if (work_co%isort(ig) == 0) call die('Could not find g-vector.')
            isorti(work_co%isort(ig))=ig
          enddo
          ! mtxel_t wants isort to sort the total g-space by kinetic energies? (same as genwf_co)
          if (xct%qflag .eq. 0) then
            qk(1:3)=kg_sub%f(1:3,ik_sub)
          else
            qk(1:3)=kg_sub%f(1:3,ik_sub)
          endif
          call kinetic_energies(gvec, crys%bdot, ekin, qvec = qk)
          call sortrx(gvec%ng, ekin, work_co%isort, gvec = gvec%components)
          call gmap(gvec,syms,work_co%ng,kg_sub%itran(ik_sub), &
            kg_sub%kg0(:,ik_sub),work_co%isort,isorti,work_co%ind,work_co%ph,.true.)
          !       Compute and renormalize wfnc wavefunctions 
          do kk=1,wfnc%nspin
            do jj=1,wfnc%nband
              do hh=1,wfnc%nspinor
                do ii=1,work_co%ng
                  if (work_co%ind(ii) .gt. 0) then
                    work_co%cg(ii,jj,kk*hh)=work_co%ph(ii)*wfnc%cg(work_co%ind(ii),jj,kk*hh)
                  else
                    work_co%cg(ii,jj,kk*hh)=ZERO
                  endif
                enddo
              enddo
              call checknorm('wfnc%cg',jj,ik,work_co%ng,kk,work_co%nspinor,work_co%cg(1:work_co%ng,jj,:))
            enddo
          enddo
          wfnc%cg=work_co%cg(:,1:wfnc%nband,:)
          wfnc%isort=work_co%isort
          !       Compute and renormalize wfnv wavefunctions 
          SAFE_DEALLOCATE_P(work_co%cg)
          SAFE_ALLOCATE(work_co%cg, (wfnv%ng,wfnv%nband,wfnv%nspin*wfnv%nspinor))
          if (xct%qflag.ne.0) then
            do kk=1,wfnv%nspin
              do jj=1,wfnv%nband
                do hh=1,wfnv%nspinor
                  do ii=1,work_co%ng
                    if (work_co%ind(ii) .gt. 0) then
                      work_co%cg(ii,jj,kk*hh)=work_co%ph(ii)*wfnv%cg(work_co%ind(ii),jj,kk*hh)
                    else
                      work_co%cg(ii,jj,kk*hh)=ZERO
                    endif
                  enddo
                enddo
                call checknorm('wfnv%cg',jj,ik,work_co%ng,kk,work_co%nspinor,work_co%cg(1:work_co%ng,jj,:))
              enddo
            enddo
            wfnv%cg=work_co%cg(:,1:wfnv%nband,:)
            wfnv%isort=work_co%isort
          else ! If reading valence wavefunction from WFNq_sub
            isorti(:)=0
            ekin = 0.0d0
            ! isort stored mapping between local and global gvecs
            do ig=1,workq_co%ng
              call findvector(workq_co%isort(ig), gvec_kpt_q%components(:, ig), gvec)
              if (workq_co%isort(ig) == 0) call die('Could not find g-vector.')
              isorti(workq_co%isort(ig))=ig
            enddo
            ! mtxel_t wants isort to sort the total g-space by kinetic energies? (same as genwf_co)
            qk(1:3)=kgq_sub%f(1:3,ik_sub)
            call kinetic_energies(gvec, crys%bdot, ekin, qvec = qk)
            call sortrx(gvec%ng, ekin, workq_co%isort, gvec = gvec%components)
            call gmap(gvec,syms,workq_co%ng,kgq_sub%itran(ik_sub), &
              kgq_sub%kg0(:,ik_sub),workq_co%isort,isorti,workq_co%ind,workq_co%ph,.true.)
            
            do kk=1,wfnv%nspin
              do jj=1,wfnv%nband
                do hh=1,wfnv%nspinor
                  do ii=1,workq_co%ng
                    if (workq_co%ind(ii) .gt. 0) then
                      workq_co%cg(ii,jj,kk*hh)=workq_co%ph(ii)*wfnv%cg(workq_co%ind(ii),jj,kk*hh)
                    else
                      workq_co%cg(ii,jj,kk*hh)=ZERO
                    endif
                  enddo
                enddo
                call checknorm('wfnv%cg',jj,ik,workq_co%ng,kk,workq_co%nspinor,workq_co%cg(1:workq_co%ng,jj,:))
              enddo
            enddo
            wfnv%cg=workq_co%cg(:,1:wfnv%nband,:)
            wfnv%isort=workq_co%isort          
          endif

!-----------------------------------------------
! Umklapps between fine and coase grids
!
          if (xct%qflag.eq.1) then
            ikq_fi = indexq_fi(ik_fi)
          else
            ikq_fi = xct%indexq_fi(ik_fi)
          endif

          ! Store the umklapp G-vector, if not zero 
          if (ik_fi .ne. 0) then
            qg(:)=kg_fi%f(:,ik_fi)-kg_sub%f(:,ik_sub)
            gumk(:)=nint(qg(:))
            call findvector(igumk, gumk, gvec)
            if (xct%qflag .eq. 1) then
              igumkq=igumk
            elseif (xct%qflag.eq.0) then
              qg(:)=kgq_fi%f(:,ikq_fi) - kgq_sub%f(:,ik_sub)
              gumk(:)=nint(qg(:))
              call findvector(igumkq, gumk, gvec)
            endif
          endif

!------------------------------ 
! Read needed fine wavefunctions: conduction bands from unshifted grid 
! and valence bands from shifted grid  
          ! Compute transformation matrices
          if (xct%qflag.eq.1) then
            call mtxel_t(gvec,wfnc, wfnc_fi, wfnv, wfnvq_fi, &
              dcn(:,:,:,ik_fi,ivert,ik_sub), dvn(:,:,:,ik_fi,ivert,ik_sub),  &
              igumk, igumkq, .not.xct%unrestricted_transf, xct%qflag)
          else
            call mtxel_t(gvec,wfnc, wfnc_fi, wfnv, wfnvq_fi, &
              dcn(:,:,:,ik_fi,ivert,ik_sub), dvn(:,:,:,ikq_fi,ivert,ik_sub),  &
              igumk, igumkq, .not.xct%unrestricted_transf, xct%qflag)
          endif
          SAFE_DEALLOCATE_P(gvec_kpt%components)
          SAFE_DEALLOCATE_P(wfnv%cg)
          SAFE_DEALLOCATE_P(wfnc%cg)
          SAFE_DEALLOCATE_P(work_co%cg)
          if (Xct%qflag.eq.0)then
            SAFE_DEALLOCATE_P(workq_co%cg)
          endif
          
        enddo !ik_sub
        
        SAFE_DEALLOCATE_P(gvec_sub%components)
        call close_file(iunit)
        if (xct%qflag.eq.0) then
          call close_file(vunit)
        endif

        SAFE_DEALLOCATE_P(kg_sub%r)
      enddo ! ivert
      
    enddo !ik_loc
    call progress_free(prog_info)

    SAFE_DEALLOCATE_P(gvec_sub%components)
    SAFE_DEALLOCATE_P(wfnv%isort)
    SAFE_DEALLOCATE_P(wfnc%isort)
    SAFE_DEALLOCATE_P(work%isort)
    SAFE_DEALLOCATE_P(work_co%ind)
    SAFE_DEALLOCATE_P(work_co%ph)
    if (xct%qflag.eq.0)then
      SAFE_DEALLOCATE_P(workq_co%ind)
      SAFE_DEALLOCATE_P(workq_co%ph)
    endif
    SAFE_DEALLOCATE(isorti)
    SAFE_DEALLOCATE(ekin)

    if (work%ikold .ne. 0) then
      SAFE_DEALLOCATE_P(work%cg)
      SAFE_DEALLOCATE_P(work%ph)
      SAFE_DEALLOCATE_P(work%ind)
      SAFE_DEALLOCATE_P(work%isort)
    endif
    if (workq%ikold .ne. 0) then
      SAFE_DEALLOCATE_P(workq%cg)
      SAFE_DEALLOCATE_P(workq%ph)
      SAFE_DEALLOCATE_P(workq%ind)
      SAFE_DEALLOCATE_P(workq%isort)
    endif

  endif ! if read_dtmat_sub
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
        
!-------------------------------------------------------------------------------
! Finished computation of transformation matrices, now have all processes share dvn, dcn
!
#ifdef MPI
  SAFE_ALLOCATE(dummy, (xct%ncb_fi,xct%n2b_co,xct%nspin,xct%nkpt_fi,xct%idimensions+1,nk_sub))
  dummy = dcn
  call MPI_Allreduce(dummy(1,1,1,1,1,1),dcn(1,1,1,1,1,1),xct%nkpt_fi*xct%ncb_fi* &
    xct%n2b_co*xct%nspin*(xct%idimensions+1)*nk_sub,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
  SAFE_DEALLOCATE(dummy)

  SAFE_ALLOCATE(dummy, (xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,xct%idimensions+1,nk_sub))
  dummy = dvn
  call MPI_Allreduce(dummy(1,1,1,1,1,1),dvn(1,1,1,1,1,1),xct%nkpt_fi*xct%nvb_fi* &
    xct%n1b_co*xct%nspin*(xct%idimensions+1)*nk_sub,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
  SAFE_DEALLOCATE(dummy)

  SAFE_ALLOCATE(dummy_closepts,(4,xct%nkpt_fi))
  dummy_closepts = closepts ! closepts(ivert,ik_fi) = 0 if ik_fi is not owned by proc
  call MPI_Allreduce(dummy_closepts(1,1),closepts(1,1),4*xct%nkpt_fi,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
  SAFE_DEALLOCATE(dummy_closepts)

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

!-------------------------------------------------------------------------------
! TODO: Save dvn and dcn to a file

  call interp_free(wfn_interp)

!--------------------------------
! Check if transformation matrices are unitary and print out their
! norm (i.e., the norm of wavefunctions in the fine grid)
! The norm should be as close as possible to 1.d0. If it is not,
! important overlaps are being neglected.
  do ivert = 1,xct%idimensions + 1
    do ik_sub = 1,nk_sub
      do is = 1, xct%nspin
        do ik_fi=1,xct%nkpt_fi
          do iv=1,xct%nvb_fi
            dnorm = 0.d0
            do iv_co=1,xct%n1b_co
              dnorm = dnorm + abs(dvn(iv,iv_co,is,ik_fi,ivert,ik_sub))**2
            enddo
            normv(ik_fi,iv,is,ivert,ik_sub) = dnorm
          enddo
          
          do ic=1,xct%ncb_fi
            dnorm = 0.d0
            do ic_co=1,xct%n2b_co
              dnorm = dnorm + abs(dcn(ic,ic_co,is,ik_fi,ivert,ik_sub))**2
            enddo
            normc(ik_fi,ic,is,ivert,ik_sub) = dnorm
          enddo
        enddo
      enddo
    enddo
  enddo

! TODO: Write out norms to file

  if (.not.flag%read_dtmat_sub .and. peinf%inode==0) then
    call open_file(20,file='dcmat_sub_norm.dat',form='formatted',status='replace')
    write(20,'(a,i12,a)') ' -------  Norm of dcc matrices : Spins = ', xct%nspin,'  -------'
    if (xct%nspin.eq.1) then
      write(20,'(a)') '           k-point           ik_co    ivert    ik_sub   c   |dcc|^2      '
    else
      write(20,'(a)') '           k-point           ik_co    ivert    ik_sub   c   |dcc|^2 (spin up)  |dcc|^2 (spin down) '
    endif
    write(20,'(21a3)') ('---',kk=1,21)
    do ik_fi=1,xct%nkpt_fi
      do ivert=1,xct%idimensions+1
        ik_co= closepts(ivert,ik_fi)
        do ik_sub=1,nk_sub
          do jj=1,xct%ncb_fi
            write(20,250) kg_fi%f(:,ik_fi),ik_co,ivert,ik_sub,jj,(normc(ik_fi,jj,is,ivert,ik_sub), is = 1, xct%nspin)
          enddo
        enddo
      enddo
    enddo
    call close_file(20)
    
    call open_file(20,file='dvmat_sub_norm.dat',form='formatted',status='replace')
    write(20,'(a,i12,a)') ' -------  Norm of dvv matrices : Spins = ', xct%nspin,'  -------'
    if (xct%nspin.eq.1) then
      write(20,'(a)') '           k-point           ik_co    ivert    ik_sub   c   dist    |dcc|^2      '
    else
      write(20,'(a)') '           k-point           ik_co    ivert    ik_sub   c   dist    |dcc|^2 (spin up)  |dcc|^2 (spin down) '
    endif
    write(20,'(21a3)') ('---',kk=1,21)
    do ik_fi=1,xct%nkpt_fi
      do ivert=1,xct%idimensions+1
        ik_co= closepts(ivert,ik_fi)
        do ik_sub=1,nk_sub
          do jj=1,xct%nvb_fi
            write(20,250) kg_fi%f(:,ik_fi),ik_co,ivert,ik_sub,jj,(normv(ik_fi,jj,is,ivert,ik_sub), is = 1, xct%nspin)
          enddo
        enddo
      enddo
    enddo
    call close_file(20)
250   format(' ( ',f6.3,' , ',f6.3,' , ',f6.3, ' ) ',i4,1x,i4,4x,i4,4x,i4,1x,f10.6,9x,f10.6)
  endif
!---------------------------------
! Renormalize transformation matrices: matrices dcn,dvn are
! renormalized to 1.d0 (ideally, this step should not be
! performed but the interpolation can be really bad if dcn,dvn
! do not have unitary norm).

  ! FHJ: Print worst norm before restriction/renormalization
  worst_normv = minval(normv)
  worst_normc = minval(normc)
  if (peinf%inode==0) then    
    write(6,'(1x,a)') 'Max. error in norm of transformation coefficients (1 - \sum |d_fi,sub|^2):'
    write(6,'(1x,a,es9.3e2)') ' - For valence states: ', 1d0 - worst_normv
    write(6,'(1x,a,es9.3e2)') ' - For conduction states: ', 1d0 - worst_normc
    call warn_norm(xct, worst_normv, worst_normc)
  endif
  ! FHJ: truncate or zero out coefficients. User doesn`t care about the details...
  if (truncate_coefs) then
    if (peinf%inode==0) write(6,'(1x,a)') "Restricting transformation to dvv',dcc' subspaces."
    dvv_sub(:,:,:,:,:,:) = dvn(:,:xct%nvb_co,:,:,:,:)
    dcc_sub(:,:,:,:,:,:) = dcn(:,xct%nvb_co+1:,:,:,:,:)
  elseif (xct%zero_unrestricted_contrib) then
    if (peinf%inode==0) write(6,'(1x,a)') "Restricting transformation to dvv',dcc' subspaces."
    dvn(:,xct%nvb_co+1:,:,:,:,:) = ZERO
    dcn(:,:xct%nvb_co,:,:,:,:) = ZERO
  endif

  ! FHJ: And from now on, xct%n1b_co will depend only on xct%extended_kernel.
  if (.not.xct%extended_kernel) then
    xct%n1b_co = xct%nvb_co
    xct%n2b_co = xct%ncb_co
  else
    xct%n1b_co = xct%nvb_co + xct%ncb_co
    xct%n2b_co = xct%n1b_co
  endif

  ! FHJ: inteqp also renormalizes the coefs for consistency (and someone might
  ! need this in the future). And note that we go here if truncate_coefs==.true.
  ! only to recalculate the worst norm.
  if (truncate_coefs.or.xct%renorm_transf) then
    worst_normc = INF
    worst_normv = INF
    do ivert = 1, xct%idimensions+1
      do ik_sub=1,nk_sub
        do is=1, xct%nspin
          do ik=1,xct%nkpt_fi
            do ic_fi=1,xct%ncb_fi
              dnorm = sum(abs(dcc_sub(ic_fi,:,is,ik,ivert,ik_sub))**2)
              worst_normc = min(worst_normc, dnorm)
              if (xct%renorm_transf) then
                dnorm = sqrt(dnorm)
                dcc_sub(ic_fi,:,is,ik,ivert,ik_sub) = dcc_sub(ic_fi,:,is,ik,ivert,ik_sub)/dnorm
              endif
            enddo
            do iv_fi=1,xct%nvb_fi
              dnorm = sum(abs(dvv_sub(iv_fi,:,is,ik,ivert,ik_sub))**2)
              worst_normv = min(worst_normv, dnorm)
              if (xct%renorm_transf) then
                dnorm = sqrt(dnorm)
                dvv_sub(iv_fi,:,is,ik,ivert,ik_sub) = dvv_sub(iv_fi,:,is,ik,ivert,ik_sub)/dnorm
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  endif
  if (peinf%inode==0.and.(truncate_coefs.or.xct%zero_unrestricted_contrib)) then
    write(6,'(1x,a)') 'Max. error after restricting WFN transformation:'
    write(6,'(1x,a,es9.3e2)') ' - For valence states: ', 1d0 - worst_normv
    write(6,'(1x,a,es9.3e2)') ' - For conduction states: ', 1d0 - worst_normc
    call warn_norm(xct, worst_normv, worst_normc)
  endif

  if (.not.associated(dvn,dvv_sub)) then
    SAFE_DEALLOCATE_P(dvn)
  endif
  if (.not.associated(dcn,dcc_sub)) then
    SAFE_DEALLOCATE_P(dcn)
  endif
  SAFE_DEALLOCATE(normc)
  SAFE_DEALLOCATE(normv)
  
  POP_SUB(intwfn_sub)

  return
end subroutine intwfn_sub

    subroutine warn_norm(xct, worst_normv, worst_normc)
      type(xctinfo), intent(in) :: xct
      real(DP), intent(in) :: worst_normv, worst_normc

      logical, save :: warned = .false.

      PUSH_SUB(warn_norm)

      if (.not.warned) then
        ! FHJ: just want to make sure there`s no band crossing that the user is missing
        if (worst_normv<0.5 .or. worst_normc<0.5) then
          write(0,'(/,a)') 'WARNING: there are fine/coarse transformation coefficients with 2-norm < 50%.'
          write(0,'(a)') 'Most likely, you are including the same number of bands in the coarse grid'
          write(0,'(a)') 'as in the fine grid, or your coarse k-grid is too coarse. You should always'
          write(0,'(a)') 'include more bands in the coarse grid to capture, for instance, band crossings.'
          write(0,'(a)') 'To improve convergence, you might want to consider:'
          write(0,'(a)') '- including more bands from the coarse WFN_co file (number_*_bands_coarse)'
          write(0,'(a)') '- decreasing the number of bands from the fine WFN_fi file (number_*_bands_fine)'
          write(0,'(a)') '- using a coarse WFN_co file with a denser k-mesh'
          if (xct%is_absorption) then
            write(0,'(a)') '- recalculating the kernel using the "extended_kernel" flag (if metal)'
          endif
          if (.not. xct%unrestricted_transf) then
            write(0,'(a)') '- using the "unrestricted_transformation" flag in the input file (if metal)'
          endif
          write(0,*)
          warned = .true.
        endif
      endif

      POP_SUB(warn_norm)

    end subroutine warn_norm

end module intwfn_m
