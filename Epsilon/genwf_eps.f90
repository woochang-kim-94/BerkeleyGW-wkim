!==========================================================================
!
! Module:
!
! genwf_eps      Originally By FHJ             Last Modified 04/2012 (FHJ)
!
!
!==========================================================================

#include "f_defs.h"

module genwf_eps_m

  use global_m
  use input_utils_m
  use fftw_m
  use genwf_mpi_m
  use gmap_m
  use susymmetries_m
  use sort_m
  use timing_m, only: timing => epsilon_timing
  use genwf_mpi_m

  implicit none

  private

  !> communicator object for the WFN FFTs
  type wfn_FFT_comm_t
    logical :: done !< Have we received all the buffers?
    integer, pointer :: req_recvv(:), req_recvc(:) !< Array of MPI_REQUEST
    integer, pointer :: req_sendv(:), req_sendc(:) !< Array of MPI_REQUEST
    integer :: recv_cntv, recv_cntc !< Number of requests
    integer :: send_cntv, send_cntc !< Number of requests
    integer :: nv !< Number of valence bands
  end type wfn_FFT_comm_t

  public :: genwf_gen, genwf_FFT, free_wfns, genwf_lvl2, &
    genwf_FFT_Isend, genwf_FFT_Wait, wfn_FFT_comm_t, &
    get_wfn_fftgrid, get_eps_fftgrid

contains  

  !> FHJ: Figure out what is the min/max gvec components give an isort array
  !! Note: you should manually initialize box_min and box_max to zero!
  subroutine get_gvecs_bounds(gvec, ng, isort, box_min, box_max)
    type(gspace), intent(in) :: gvec
    integer, intent(in) :: ng
    integer, intent(in) :: isort(:)
    integer, intent(inout) :: box_min(3), box_max(3)

    integer :: ig

    PUSH_SUB(get_gvecs_bounds)

    do ig=1,ng
      box_min(1:3) = min(box_min(1:3),  gvec%components(1:3, isort(ig)))
      box_max(1:3) = max(box_max(1:3),  gvec%components(1:3, isort(ig)))
    enddo

    POP_SUB(get_gvecs_bounds)
    return

  end subroutine get_gvecs_bounds

  !> FHJ: Figure out the minimum fftbox that holds all the WFNs
  subroutine get_wfn_fftgrid(pol, gvec, kp, intwfn)
    type(polarizability), intent(inout) :: pol
    type(gspace), intent(in) :: gvec
    type(kpoints), intent(in) :: kp
    type(int_wavefunction), intent(in) :: intwfn

    integer :: ik, wfn_box_min(3), wfn_box_max(3)

    PUSH_SUB(get_wfn_fftgrid)

    wfn_box_min(:) = 0; wfn_box_max(:) = 0
    do ik=1,kp%nrk
      call get_gvecs_bounds(gvec, intwfn%ng(ik), intwfn%isort(:,ik), wfn_box_min, wfn_box_max)
    enddo
    pol%WFN_FFTgrid(1:3) = wfn_box_max(1:3) - wfn_box_min(1:3) + 1
    if (peinf%verb_debug .and. peinf%inode==0) then
      write(6,*) 'WFN min. FFT grid:',pol%WFN_FFTgrid
    endif

    POP_SUB(get_wfn_fftgrid)
    return

  end subroutine get_wfn_fftgrid

  !> FHJ: Figure out the minimum fftbox that allows us to convolve the WFNs within
  !! an energy window of nmtx G vectors.
  subroutine get_eps_fftgrid(pol, gvec)
    type(polarizability), intent(inout) :: pol
    type(gspace), intent(in) :: gvec
    integer :: eps_box_min(3), eps_box_max(3)

    PUSH_SUB(get_eps_fftgrid)

    eps_box_min(:) = 0; eps_box_max(:) = 0
    call get_gvecs_bounds(gvec, pol%nmtx, pol%isrtx, eps_box_min, eps_box_max)
    ! FHJ: Note: the amount of padding is actually N_sig + N_window - 1
    pol%FFTgrid(1:3) = pol%WFN_FFTgrid(1:3) + (eps_box_max(1:3) - eps_box_min(1:3))
    pol%FFTgrid(1:3) = min(pol%FFTgrid(1:3), gvec%FFTgrid(1:3))
    if (peinf%verb_debug .and. peinf%inode==0) then
      write(6,'(1x,a,3(1x,i0))') 'Original FFT grid:', gvec%FFTgrid
      write(6,'(1x,a,3(1x,i0))') 'Minimal  FFT grid:', pol%FFTgrid
    endif

    POP_SUB(get_eps_fftgrid)
    return

  end subroutine get_eps_fftgrid

  !> FHJ: to be used internally with genwf_FFT_Isend
  subroutine do_my_FFTs(this,gvec,Nfft,wfn_fft,intwfn,my_bands,my_cnt,ng,tmp_wfn,isort,ind,ph,is_val)
    type (wfn_FFT_comm_t), intent(inout), target :: this !< communicator object for the WFN FFTs
    type(gspace), intent(in) :: gvec
    complex(DPC), intent(inout) :: wfn_fft(:,:,:,:)
    integer, intent(in) :: Nfft(3)
    type(int_wavefunction), intent(in) :: intwfn
    integer, intent(in) :: my_bands(:)  !< my bands
    integer, intent(in) :: my_cnt       !< number of FFTs to do = sizeof(my_bands)
    integer, intent(in) :: ng           !< number of gvectors for the k-pt in question
    SCALAR,  intent(inout) :: tmp_wfn(:)!< buffer to reorder WFN using ind and ph
    integer, intent(in) :: isort(:)     !< wfn isort
    integer, intent(in) :: ind(:)
    SCALAR,  intent(in) :: ph(:)    
    logical, intent(in) :: is_val       !< .true. to take the conjg_fftbox

    integer :: is, ig, fft_size
    integer :: ib_list, ib_local, ib, iproc, offset
    integer, pointer :: invindex(:), send_cnt, req_send(:)
    logical, pointer :: does_it_own(:,:)

    PUSH_SUB(do_my_FFTs)

    is = 1
    if ( is_val ) then
      invindex => peinf%invindexv
      does_it_own => peinf%does_it_ownv
      req_send => this%req_sendv
      send_cnt => this%send_cntv
      offset = 0
    else
      invindex => peinf%invindexc
      does_it_own => peinf%does_it_ownc
      req_send => this%req_sendc
      send_cnt => this%send_cntc
      offset = this%nv
    endif
    fft_size = product(Nfft(1:3))

    do ib_list = 1, my_cnt
      ib_local = my_bands(ib_list)
      do ig=1,ng
        tmp_wfn(ig) = intwfn%cg(ind(ig), ib_local, is)*ph(ig)
      enddo
      call timing%start(timing%opt_fft_fft)
      call put_into_fftbox(ng, tmp_wfn, gvec%components, &
        isort, wfn_fft(:,:,:,ib_local), Nfft)
      call do_FFT( wfn_fft(:,:,:,ib_local), Nfft, 1)
      if ( is_val ) call conjg_fftbox( wfn_fft(:,:,:,ib_local), Nfft )
      call timing%stop(timing%opt_fft_fft)

#ifdef MPI
      ! FHJ: distribute all real-space WFNs
      call timing%start(timing%opt_fft_comm_fft)
      ib = invindex(ib_local)
      do iproc=0, peinf%npes-1
        if (does_it_own(ib, iproc+1).and.iproc/=peinf%inode) then
          send_cnt = send_cnt + 1
          call MPI_Isend( wfn_fft(1,1,1,ib_local), fft_size, MPI_COMPLEX_DPC, iproc, &
            offset + ib, MPI_COMM_WORLD, req_send(send_cnt), mpierr )
        endif
      enddo
      call timing%stop(timing%opt_fft_comm_fft)
#endif        
    enddo

    POP_SUB(do_my_FFTs)
    return

  end subroutine do_my_FFTs

  !> Generates all the real-space wavefunctions. Used only if pol%os_opt_fft==2
  !! This version avoids communication, and it`s under development!
  !!TODO`s:
  !! (1) we are just supporting one spin and one kpt/qpt.
  !! (2) support serial code
  subroutine genwf_FFT_Isend(this,crys,gvec,syms,kp,kpq,vwfn,pol,cwfn,intwfnv,intwfnvq,intwfnc)
    type (wfn_FFT_comm_t), intent(inout) :: this !< communicator object for the WFN FFTs
    type (crystal), intent(in) :: crys
    type (gspace), intent(in) :: gvec
    type (symmetry), intent(in) :: syms
    type (kpoints), target, intent(in) :: kp
    type (kpoints), target, intent(in) :: kpq
    type (valence_wfns), intent(inout) :: vwfn
    type (polarizability), intent(inout) :: pol
    type (conduction_wfns), intent(inout) :: cwfn
    type (int_wavefunction), intent(inout) :: intwfnv
    type (int_wavefunction), intent(inout) :: intwfnvq
    type (int_wavefunction), intent(inout) :: intwfnc

    integer :: npes_per_pool, nc_groups, nproc_max, nc
    integer, allocatable :: nv_bands(:), v_owners(:)
    integer :: Nfft(3), fft_size
    real(DP) :: scale
    integer, allocatable :: my_vbands(:)
    integer :: my_vcnt, iv, iv_local
    integer :: ipool, isubrank, inode, ioffset
    integer, allocatable :: grp_global_ranks(:), ntot_bands(:), &
      grp_local_owners(:) 
    integer, allocatable :: my_cbands(:)
    integer :: my_ccnt, ic, ic_local
    integer :: my_grp_rank, grp_nprocs
    integer :: min_bands, iproc, iworker
    integer :: ik, is, ib

    !sort stuff
    SCALAR, allocatable :: tmp_wfn(:)
    SCALAR, allocatable :: ph(:)
    real(DP), allocatable :: ekin(:)
    integer, allocatable :: ind(:), isorti(:)
    integer, allocatable, target :: wfn_isort(:)
    integer :: ng0, ig

    PUSH_SUB(genwf_FFT_Isend)

    call logit('generating all real-space wavefunctions')

    call timing%start(timing%opt_fft)

    if(pol%nq>1.or.kp%nrk>1.or.kpq%nrk>0.or.pol%need_WFNq) &
      call die('FFT opt. level 2 only works for 1 qpt and 1 kpt, and without WFNq',&
      only_root_writes=.true.)

#ifdef MPI

    ik = 1 ! Fixing one kpt for now
    is = 1 ! And the spin
    vwfn%idx_kp = ik
    cwfn%idx_kp = ik    
    nc = cwfn%nband - vwfn%nband  ! `real` number of conduction bands
    npes_per_pool = (peinf%npes/peinf%npools) ! number of processors per pool
    nc_groups = (nc + peinf%ncownmax - 1)/(peinf%ncownmax) ! number of conduction groups     
    nproc_max = (peinf%npes + nc_groups - 1)/(nc_groups) ! max num. of proc. per conduction group

    ! Basic initialization of the communicator object
    this%recv_cntv = 0; this%recv_cntc = 0
    this%send_cntv = 0; this%send_cntc = 0
    this%done = .false.
    this%nv = vwfn%nband
    SAFE_ALLOCATE(this%req_recvv, (vwfn%nband*npes_per_pool))
    SAFE_ALLOCATE(this%req_recvc, (nc*nproc_max))
    SAFE_ALLOCATE(this%req_sendv, (vwfn%nband*npes_per_pool))
    SAFE_ALLOCATE(this%req_sendc, (nc*nproc_max))

    call timing%start(timing%opt_fft_init)

    ! FHJ: prepare the sorting buffers
    ng0 = intwfnv%ng(ik) ! We all have at least one valence band
    SAFE_ALLOCATE(ind,(ng0))
    SAFE_ALLOCATE(ph, (ng0))
    SAFE_ALLOCATE(wfn_isort, (gvec%ng))
    SAFE_ALLOCATE(isorti, (gvec%ng))
    SAFE_ALLOCATE(ekin, (gvec%ng))
    call kinetic_energies(gvec, crys%bdot, ekin)
    call sortrx(gvec%ng, ekin, wfn_isort, gvec=gvec%components)

    do ig=1,gvec%ng
      isorti(wfn_isort(ig)) = ig
    enddo
    do ig=1,ng0
      isorti(intwfnv%isort(ig, ik)) = ig
    enddo
    call gmap(gvec,syms,ng0,1,(/0,0,0/),wfn_isort,isorti,ind,ph,.true.)
    vwfn%ngv=ng0
    cwfn%ngc=ng0

    ! FHJ: this will be more complicated to implement for more than 1 kpt
    if (pol%min_fftgrid) then
      pol%isrtx => wfn_isort
      pol%nmtx = gcutoff(gvec%ng, ekin, pol%isrtx, pol%ecuts)
      call get_eps_fftgrid(pol, gvec)
      nullify(pol%isrtx)
    endif
    call setup_FFT_sizes(pol%FFTgrid,Nfft,scale)
    fft_size = Nfft(1)*Nfft(2)*Nfft(3)

    ! FHJ: Distribute VALENCE bands !
    !--------------------------------
    ! NOTE: we index the processors wrt their real MPI ranks (peinf%inode)
    !       we refer all the val bands wrt the local index in my_vbands

    ! Number of valence bands that a particular processor owns      
    SAFE_ALLOCATE(nv_bands, (peinf%npes))
    ! Who owns a particular valence band?
    SAFE_ALLOCATE(v_owners, (vwfn%nband + pol%ncrit))
    ! Which valence bands do I own?
    SAFE_ALLOCATE(my_vbands, (peinf%nvownmax))
    nv_bands(:) = 0; v_owners(:) = 0
    my_vbands(:) = 0
    my_vcnt = 0

    do iv=1, vwfn%nband + pol%ncrit
      ipool = (iv-1)/peinf%nvownmax ! Get pool for the band
      ! Distribute bands to the processors in a round-robin way. We add an
      ! offset so that we get a good load balance for the cond. bands later on.
      isubrank = mod(iv-1, peinf%nvownmax)
      ioffset = ipool*peinf%nvownmax
      inode = mod(isubrank+ioffset, npes_per_pool) + ipool*npes_per_pool
      v_owners(iv) = inode
      ! Keep track of the # of times that each processors got a band:
      nv_bands(inode+1) = nv_bands(inode+1) + 1
      if (peinf%inode == inode) then
        my_vcnt = my_vcnt + 1
        my_vbands(my_vcnt) = isubrank + 1
      endif
    enddo

    ! FHJ: Distribute CONDUCTION bands !
    !----------------------------------!
    ! NOTE: we use *local ranks* to organize the processors in a particular cond. group
    !       we refer all the cond. bands wrt the local index in my_cbands

    ! Ranks of the processors that participate in this conduction group
    SAFE_ALLOCATE(grp_global_ranks, (nproc_max))
    grp_global_ranks(:) = -1
    SAFE_ALLOCATE(ntot_bands, (nproc_max)) ! Total number of bands per processor
    ntot_bands(:) = 0
    ! Who owns a particular conduction band? (in terms of local workers)
    SAFE_ALLOCATE(grp_local_owners, (peinf%ncownactual))
    grp_local_owners(:) = 0
    ! Which conduction bands do I own?
    SAFE_ALLOCATE(my_cbands, (peinf%ncownactual))
    my_cbands(:) = 0
    my_ccnt = 0

    ! This will be set to my rank within the group of conduction bands, starting at 0.
    my_grp_rank = -1

    ! Create list of all processors in that same group
    ! want: grp_global_ranks(:) = [icgroup, icgroup + npes_per_group, ...]
    grp_global_ranks(:) = 0
    grp_nprocs = 0
    ib = peinf%invindexc(1) ! index of first band that I own.
    if (ib>0) then
      do iproc=0, peinf%npes-1
        if (peinf%does_it_ownc(ib, iproc+1)) then
          grp_nprocs = grp_nprocs + 1
          grp_global_ranks(grp_nprocs) = iproc
          if (iproc==peinf%inode) then
            my_grp_rank = grp_nprocs - 1
          endif
        endif
      enddo
    endif

    ! Initialize list of the number of bands that each processor has
    do iproc = 1, grp_nprocs
      ntot_bands(iproc) = nv_bands(grp_global_ranks(iproc)+1)
    enddo

    do ic_local = 1, peinf%ncownactual
      ! Get smallest index of ntot_bands array -> iworker
      ! note: rank = grp_local_ranks(iworker)
      iworker = 0
      min_bands = cwfn%nband + 1
      do iproc=1, grp_nprocs
        if (ntot_bands(iproc)<min_bands) then
          min_bands = ntot_bands(iproc)
          iworker = iproc - 1
        endif
      enddo

      ! add ic_local to list of bands that iworker owns
      grp_local_owners(ic_local) = iworker
      ntot_bands(iworker+1) = ntot_bands(iworker+1) + 1
      if (my_grp_rank == iworker) then
        my_ccnt = my_ccnt + 1
        my_cbands(my_ccnt) = ic_local
      endif
    enddo

    call timing%stop(timing%opt_fft_init)
    call timing%start(timing%opt_fft_comm_fft)

    ! FHJ: Non-blocking Recvs !
    !-------------------------!

    ! Recv valence bands
    SAFE_ALLOCATE( vwfn%wfn_fft, (Nfft(1), Nfft(2), Nfft(3), peinf%nvownactual) )
    do iv_local=1, peinf%nvownactual
      iv = peinf%invindexv(iv_local)
      if (v_owners(iv)/=peinf%inode) then 
        this%recv_cntv = this%recv_cntv + 1
        call MPI_Irecv( vwfn%wfn_fft(1,1,1,iv_local), fft_size, MPI_COMPLEX_DPC, v_owners(iv), &
          iv, MPI_COMM_WORLD, this%req_recvv(this%recv_cntv), mpierr )
      endif
    enddo
    ! Recv conduction bands
    SAFE_ALLOCATE( cwfn%wfn_fft, (Nfft(1), Nfft(2), Nfft(3), peinf%ncownactual) )
    do ic_local=1, peinf%ncownactual
      ic = peinf%invindexc(ic_local)
      if (grp_local_owners(ic_local)/=my_grp_rank) then
        this%recv_cntc = this%recv_cntc + 1
        iproc = grp_global_ranks(grp_local_owners(ic_local)+1)
        call MPI_Irecv( cwfn%wfn_fft(1,1,1,ic_local), fft_size, MPI_COMPLEX_DPC, iproc, &
          vwfn%nband + ic, MPI_COMM_WORLD, this%req_recvc(this%recv_cntc), mpierr )
      endif
    enddo

    call timing%stop(timing%opt_fft_comm_fft)

    ! FHJ: Do the FFTs + Isend`s !
    !----------------------------!

    call logit('doing FFTs')

    SAFE_ALLOCATE(tmp_wfn, (ng0))
    call do_my_FFTs(this,gvec,Nfft,vwfn%wfn_fft,intwfnv,my_vbands,my_vcnt,&
      ng0,tmp_wfn,wfn_isort,ind,ph,.true.)
    call do_my_FFTs(this,gvec,Nfft,cwfn%wfn_fft,intwfnc,my_cbands,my_ccnt,&
      ng0,tmp_wfn,wfn_isort,ind,ph,.false.)
    SAFE_DEALLOCATE(tmp_wfn)

    ! FHJ: Deallocate original wavefunctions (but only the wfn!)
    call free_wfns(pol, intwfnv, intwfnvq, intwfnc, .false.)

    call logit('done with FFTs')

    ! Free comm buffers
    SAFE_DEALLOCATE(nv_bands)
    SAFE_DEALLOCATE(v_owners)
    SAFE_DEALLOCATE(my_vbands)
    SAFE_DEALLOCATE(grp_global_ranks)
    SAFE_DEALLOCATE(ntot_bands)
    SAFE_DEALLOCATE(grp_local_owners)
    SAFE_DEALLOCATE(my_cbands)

    ! and free sorting stuff
    SAFE_DEALLOCATE(ind)
    SAFE_DEALLOCATE(ph)
    SAFE_DEALLOCATE(wfn_isort)
    SAFE_DEALLOCATE(isorti)
    SAFE_DEALLOCATE(ekin)

#endif

    call logit('done generating real-space wavefunctions')

    call timing%stop(timing%opt_fft)

    POP_SUB(genwf_FFT_Isend)
    return

  end subroutine genwf_FFT_Isend


  ! FHJ: call me after genwf_FFT_Isend, but just before you actually need the data
  subroutine genwf_FFT_Wait(this)
    type (wfn_FFT_comm_t), intent(inout) :: this !< communicator object for the WFN FFTs

    PUSH_SUB(genwf_FFT_Wait)

#ifdef MPI      
    call timing%start(timing%opt_fft)
    call timing%start(timing%opt_fft_comm_fft)
    if (this%recv_cntv>0) then
      call MPI_Waitall(this%recv_cntv, this%req_recvv,  MPI_STATUSES_IGNORE, mpierr)
    endif
    if (this%send_cntv>0) then
      call MPI_Waitall(this%send_cntv, this%req_sendv,  MPI_STATUSES_IGNORE, mpierr)
    endif
    if (this%recv_cntc>0) then
      call MPI_Waitall(this%recv_cntc, this%req_recvc,  MPI_STATUSES_IGNORE, mpierr)
    endif
    if (this%send_cntc>0) then
      call MPI_Waitall(this%send_cntc, this%req_sendc,  MPI_STATUSES_IGNORE, mpierr)
    endif
    call timing%stop(timing%opt_fft_comm_fft)

    SAFE_DEALLOCATE_P(this%req_recvv)
    SAFE_DEALLOCATE_P(this%req_sendv)
    SAFE_DEALLOCATE_P(this%req_recvc)
    SAFE_DEALLOCATE_P(this%req_sendc)
    this%done = .true.
    call timing%stop(timing%opt_fft)
#endif

    POP_SUB(genwf_FFT_Wait)
    return

  end subroutine genwf_FFT_Wait


  subroutine genwf_lvl2(kp,kpq,vwfn,pol,cwfn)
    type (kpoints), target, intent(in) :: kp
    type (kpoints), target, intent(in) :: kpq
    type (valence_wfns), intent(inout) :: vwfn
    type (polarizability), intent(in) :: pol
    type (conduction_wfns), intent(inout) :: cwfn

    type(kpoints), pointer :: kp_point

    PUSH_SUB(genwf_lvl2)

    if(pol%need_WFNq) then ! FIXME I think this is wrong if pol%nq1>0
      kp_point => kpq
    else
      kp_point => kp
    endif

    SAFE_ALLOCATE(vwfn%ev, (vwfn%nband+pol%ncrit, kp%nspin))
    SAFE_ALLOCATE(cwfn%ec, (cwfn%nband,kp%nspin))
    vwfn%ev(1:vwfn%nband+pol%ncrit,1:kp%nspin) = &
      kp_point%el(1:vwfn%nband+pol%ncrit, vwfn%idx_kp, 1:kp%nspin)
    cwfn%ec(1:cwfn%nband,1:kp%nspin) = &
      kp%el(1:cwfn%nband, cwfn%idx_kp, 1:kp%nspin)

    POP_SUB(genwf_lvl2)
    return
  end subroutine genwf_lvl2


  !> Generates all the real-space wavefunctions. Used only if pol%os_opt_fft==2
  !!TODO`s:
  !! (1) we are just supporting one spin and one kpt/qpt.
  !! (2) communication can be reduced if we distribute the WFNs in a smarter way
  !! (3) support serial code
  subroutine genwf_FFT(crys,gvec,syms,kp,kpq,vwfn,pol,cwfn,intwfnv,intwfnvq,intwfnc)
    type (crystal), intent(in) :: crys
    type (gspace), intent(in) :: gvec
    type (symmetry), intent(in) :: syms
    type (kpoints), target, intent(in) :: kp
    type (kpoints), target, intent(in) :: kpq
    type (valence_wfns), intent(inout) :: vwfn
    type (polarizability), intent(inout) :: pol
    type (conduction_wfns), intent(inout) :: cwfn
    type (int_wavefunction), intent(inout) :: intwfnv
    type (int_wavefunction), intent(inout) :: intwfnvq
    type (int_wavefunction), intent(inout) :: intwfnc

    integer :: ik, is, ipe
    integer :: own_max
    integer, allocatable :: band_owners(:)
    integer :: ib_loc, ib, ng0
    integer :: receiver, recv_cnt, send_cnt
    integer :: local_band_idx
    integer, allocatable :: req_send(:), req_recv(:)
    integer, allocatable :: my_bands(:), isort0(:)
    SCALAR, allocatable :: bufs_wfn(:,:)
    complex(DPC), allocatable :: work_ffts(:,:,:,:)
    integer :: Nfft(3), fft_size
    real(DP) :: scale
    SCALAR, allocatable :: tmp_wfn(:)
    SCALAR, allocatable :: ph(:)
    real(DP), allocatable :: ekin(:)
    integer, allocatable :: ind(:), isorti(:)
    integer, allocatable, target :: wfn_isort(:)
    integer :: ig

    PUSH_SUB(genwf_FFT)

    call logit('generating all real-space wavefunctions')

    call timing%start(timing%opt_fft)

    if(pol%nq>1.or.kp%nrk>1.or.kpq%nrk>0.or.pol%need_WFNq) &
      call die('FFT opt. level 2 only works for 1 qpt and 1 kpt, and without WFNq',&
      only_root_writes=.true.)

#ifdef MPI

    own_max = (cwfn%nband + peinf%npes - 1)/(peinf%npes)

    ik = 1 ! Fixing one kpt for now
    is = 1 ! And the spin

    vwfn%idx_kp = ik
    cwfn%idx_kp = ik

    call timing%start(timing%opt_fft_init)

    ! FHJ: Only root is guaranteed to have at least one band.
    if (peinf%inode==0) ng0 = intwfnv%ng(ik)
    call MPI_Bcast(ng0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    SAFE_ALLOCATE(isort0, (ng0))
    if (peinf%inode==0) isort0(1:ng0) = intwfnv%isort(1:ng0, ik)
    call MPI_Bcast(isort0(1), ng0, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

    ! FHJ: prepare the sorting stuff
    SAFE_ALLOCATE(ind,(ng0))
    SAFE_ALLOCATE(ph, (ng0))
    SAFE_ALLOCATE(wfn_isort, (gvec%ng))
    SAFE_ALLOCATE(isorti, (gvec%ng))
    SAFE_ALLOCATE(ekin, (gvec%ng))
    SAFE_ALLOCATE(tmp_wfn, (ng0))
    call kinetic_energies(gvec, crys%bdot, ekin)
    call sortrx(gvec%ng, ekin, wfn_isort, gvec=gvec%components)
    do ig=1,gvec%ng
      isorti(wfn_isort(ig)) = ig
    enddo
    do ig=1,ng0
      isorti(isort0(ig)) = ig
    enddo
    call gmap(gvec,syms,ng0,1,(/0,0,0/),wfn_isort,isorti,ind,ph,.true.)
    vwfn%ngv=ng0
    cwfn%ngc=ng0

    ! FHJ: Work buffers
    SAFE_ALLOCATE(band_owners, (cwfn%nband))
    SAFE_ALLOCATE(my_bands, (cwfn%nband))
    SAFE_ALLOCATE(req_send, (cwfn%nband))
    SAFE_ALLOCATE(req_recv, (cwfn%nband))
    SAFE_ALLOCATE(bufs_wfn, (ng0, own_max))

    ! FHJ: Re-distribute WFNs
    recv_cnt = 0
    send_cnt = 0
    do ib = 1, cwfn%nband
      receiver = mod(ib-1, peinf%npes)
      band_owners(ib) = receiver

      ! FHJ: Receiving part
      if (receiver==peinf%inode) then
        recv_cnt = recv_cnt + 1
        my_bands(recv_cnt) = ib
        call MPI_Irecv( bufs_wfn(1, recv_cnt), ng0, MPI_SCALAR, MPI_ANY_SOURCE, &
          ib, MPI_COMM_WORLD, req_recv(recv_cnt), mpierr )
      endif

      ! FHJ: Sending part
      if (ib<=vwfn%nband) then
        local_band_idx = peinf%indexv(ib)+(ik-1)*peinf%nvownactual
        if (should_send(peinf%does_it_ownv)) &
          call MPI_Isend( intwfnv%cg(1, local_band_idx, is), ng0, MPI_SCALAR, receiver, ib, &
          MPI_COMM_WORLD, req_send(send_cnt), mpierr )
      else
        local_band_idx = peinf%indexc(ib-vwfn%nband)+(ik-1)*peinf%ncownactual
        if (should_send(peinf%does_it_ownc)) &
          call MPI_Isend( intwfnc%cg(1, local_band_idx, is), ng0, MPI_SCALAR, receiver, ib, &
          MPI_COMM_WORLD, req_send(send_cnt), mpierr )
      endif
    enddo

    if (recv_cnt>0) then
      call MPI_Waitall(recv_cnt, req_recv,  MPI_STATUSES_IGNORE, mpierr)
    endif
    if (send_cnt>0) then
      call MPI_Waitall(send_cnt, req_send,  MPI_STATUSES_IGNORE, mpierr)
    endif

    call timing%stop(timing%opt_fft_init)
    call timing%start(timing%opt_fft_fft)

    ! FHJ: TODO - change genwf_mpi.f90 so that it doesn`t use cwfn%el
    ! (see genwf_mpi.f90)

    ! FHJ: Deallocate original wavefunctions (but only the wfn!)
    call free_wfns(pol, intwfnv, intwfnvq, intwfnc, .false.)

    call logit('doing FFTs')

    ! FHJ: this will be more complicated to implement for more than 1 kpt
    if (pol%min_fftgrid) then
      pol%isrtx => wfn_isort
      pol%nmtx = gcutoff(gvec%ng, ekin, pol%isrtx, pol%ecuts)
      call get_eps_fftgrid(pol, gvec)
      nullify(pol%isrtx)
    endif
    ! FHJ: Do all the FFTs
    call setup_FFT_sizes(pol%FFTgrid,Nfft,scale)
    fft_size = Nfft(1)*Nfft(2)*Nfft(3)

    SAFE_ALLOCATE( vwfn%wfn_fft, (Nfft(1), Nfft(2), Nfft(3), peinf%nvownactual) )
    SAFE_ALLOCATE( cwfn%wfn_fft, (Nfft(1), Nfft(2), Nfft(3), peinf%ncownactual) )
    SAFE_ALLOCATE( work_ffts, (Nfft(1), Nfft(2), Nfft(3), recv_cnt) )

    if (peinf%inode==0) write(0,*)
    do ib = 1, recv_cnt
      do ig=1,ng0
        tmp_wfn(ig) = bufs_wfn(ind(ig), ib)*ph(ig)
      enddo
      call put_into_fftbox(ng0, tmp_wfn, gvec%components, wfn_isort, work_ffts(:,:,:,ib), Nfft)
      call do_FFT(work_ffts(:,:,:,ib), Nfft, 1)
      if (my_bands(ib)<=vwfn%nband) then
        call conjg_fftbox(work_ffts(:,:,:,ib), Nfft)
      endif
    enddo

    call timing%stop(timing%opt_fft_fft)
    call timing%start(timing%opt_fft_comm_fft)

    call logit('done with FFTs')

    do ib_loc = 1,peinf%nvownactual
      ib = peinf%invindexv(ib_loc)
      call MPI_Irecv( vwfn%wfn_fft(1,1,1,ib_loc), fft_size, MPI_COMPLEX_DPC, band_owners(ib), &
        ib, MPI_COMM_WORLD, req_recv(ib_loc), mpierr )
    enddo
    do ib_loc = 1,peinf%ncownactual
      ib = peinf%invindexc(ib_loc)+vwfn%nband
      call MPI_Irecv( cwfn%wfn_fft(1,1,1,ib_loc), fft_size, MPI_COMPLEX_DPC, band_owners(ib), &
        ib, MPI_COMM_WORLD, req_recv(peinf%nvownactual+ib_loc), mpierr )
    enddo

    ! FHJ: And send them back using point-to-point communication.
    send_cnt = 0
    do ib_loc = 1, recv_cnt
      ib = my_bands(ib_loc)

      if (ib<=vwfn%nband) then
        do ipe = 0, peinf%npes-1
          if ( peinf%does_it_ownv(ib,ipe+1) ) then
            send_cnt = send_cnt + 1
            call MPI_Isend( work_ffts(1,1,1,ib_loc), fft_size, MPI_COMPLEX_DPC, ipe, &
              ib, MPI_COMM_WORLD, req_send(send_cnt), mpierr )
          endif
        enddo
      else
        do ipe = 0, peinf%npes-1
          if ( peinf%does_it_ownc(ib-vwfn%nband,ipe+1) ) then
            send_cnt = send_cnt + 1
            call MPI_Isend( work_ffts(1,1,1,ib_loc), fft_size, MPI_COMPLEX_DPC, ipe, &
              ib, MPI_COMM_WORLD, req_send(send_cnt), mpierr )
          endif
        enddo
      endif
    enddo

    recv_cnt = peinf%nvownactual + peinf%ncownactual
    if (recv_cnt>0) then
      call MPI_Waitall(recv_cnt, req_recv,  MPI_STATUSES_IGNORE, mpierr)
    endif
    if (send_cnt>0) then
      call MPI_Waitall(send_cnt, req_send,  MPI_STATUSES_IGNORE, mpierr)
    endif

    call timing%stop(timing%opt_fft_comm_fft)

    SAFE_DEALLOCATE(work_ffts)
    SAFE_DEALLOCATE(bufs_wfn)
    SAFE_DEALLOCATE(req_send)
    SAFE_DEALLOCATE(req_recv)
    SAFE_DEALLOCATE(my_bands)
    SAFE_DEALLOCATE(band_owners)
    SAFE_DEALLOCATE(isort0)

    ! and free sorting stuff
    SAFE_DEALLOCATE(ind)
    SAFE_DEALLOCATE(ph)
    SAFE_DEALLOCATE(wfn_isort)
    SAFE_DEALLOCATE(isorti)
    SAFE_DEALLOCATE(ekin)
    SAFE_DEALLOCATE(tmp_wfn)

    call logit('done generating real-space wavefunctions')

    call timing%stop(timing%opt_fft)

#endif

    POP_SUB(genwf_FFT)
    return

  contains

    logical function should_send(own_arr)
      logical, intent(in) :: own_arr(:,:)

      integer :: sender, ib_

      PUSH_SUB(genwf_FFT.should_send)

      ib_ = ib
      if (ib_ > vwfn%nband) ib_ = ib_ - vwfn%nband

      should_send = .false.
#ifdef MPI
      sender = -1
      do ipe = 0, peinf%npes-1
        if (own_arr(ib_, ipe+1)) then
          sender = ipe
          exit
        endif
      enddo

      if (sender==-1) call die('No sender found!')

      if (sender==peinf%inode) then
        send_cnt = send_cnt + 1
        should_send = .true.
      endif
#endif

      POP_SUB(genwf_FFT.should_send)
      return

    end function should_send

  end subroutine genwf_FFT

  !> Generic routine that generates wavefunctions 
  subroutine genwf_gen(syms,gvec,crys,kp,kpq,irk,rk,qq,vwfn,pol,cwfn,use_wfnq,intwfnv,intwfnvq,intwfnc,iv)
    type (symmetry), intent(in) :: syms
    type (gspace), intent(in) :: gvec
    type (crystal), intent(in) :: crys
    type (kpoints), target, intent(in) :: kp
    type (kpoints), target, intent(in) :: kpq
    integer, intent(in) :: irk
    real(DP), intent(in) :: rk(3)
    real(DP), intent(in) :: qq(3)
    type (valence_wfns), intent(inout) :: vwfn
    type (polarizability), intent(in) :: pol
    type (conduction_wfns), intent(inout) :: cwfn
    logical, intent(in) :: use_wfnq
    type (int_wavefunction), intent(in) :: intwfnv
    type (int_wavefunction), intent(in) :: intwfnvq
    type (int_wavefunction), intent(in) :: intwfnc
    integer, intent(in) :: iv

    integer :: ivr

    PUSH_SUB(genwf_gen)

    call logit('calling genwf')
    call timing%start(timing%genwf)
    if (iv .le. peinf%nvownactual) then
      call genwf_mpi(syms,gvec,crys,kp,kpq,irk,rk,qq,vwfn,pol,cwfn,use_wfnq,intwfnv,intwfnvq,intwfnc,iv)
    endif
    call timing%stop(timing%genwf)
    call logit('done genwf')

    POP_SUB(genwf_gen)
    return

  end subroutine genwf_gen

  !> Deallocate all "intermediate" wavefunctions
  subroutine free_wfns(pol, intwfnv, intwfnvq, intwfnc, free_all)
    type(polarizability), intent(in) :: pol
    type(int_wavefunction), intent(inout) :: intwfnv, intwfnvq, intwfnc
    logical, intent(in) :: free_all !< if .false., then only %ng and %isort will be preserved

    PUSH_SUB(free_wfns)

    if (free_all) then
      SAFE_DEALLOCATE_P(intwfnv%ng)
      SAFE_DEALLOCATE_P(intwfnv%isort)
    endif
    SAFE_DEALLOCATE_P(intwfnv%cg)
    SAFE_DEALLOCATE_P(intwfnv%qk)
    if (pol%need_WFNq) then
      if (free_all) then
        SAFE_DEALLOCATE_P(intwfnvq%ng)
        SAFE_DEALLOCATE_P(intwfnvq%isort)
      endif
      SAFE_DEALLOCATE_P(intwfnvq%cg)
      SAFE_DEALLOCATE_P(intwfnvq%qk)
    endif
    if (free_all) then
      SAFE_DEALLOCATE_P(intwfnc%ng)
      SAFE_DEALLOCATE_P(intwfnc%isort)
    endif
    SAFE_DEALLOCATE_P(intwfnc%cg)
    SAFE_DEALLOCATE_P(intwfnc%cbi)
    SAFE_DEALLOCATE_P(intwfnc%qk)

    POP_SUB(free_wfns)
    return

  end subroutine free_wfns

end module genwf_eps_m
