! ========================================================================
!
! Routines:
!
! (1) plotxctdens(main) By PW Doak              Last Modified 5/1/2014
!
!     Reads files WFN_fi/WFNq_fi (from mean field) and eigenvectors (output
!     from BSE/absorption code) and plots the change excited state charge
!     density in real space from the expression
!
!     delta rho(r) = Sum_cvk (A_svck * Conj(A_svck) * (|phi_sck(r) * exp(ik.r)|^2  - 
!         |phi_svk(r) * exp(i[k+q].r)|^2 - Sum_vp.neq.v A_svck*Conj(Asvpck)*
!           phi_sv(r) * Conj(phi_svp(r)) + Sum_cp.neq.c A_svck*Conj(Asvcpk)*
!           phi_sc(r) * Conj(phi_sc(r)) )
!
!     where r runs over the supercell. Observables such as the change in
!     multipole moments on excitation can be calculated from this signed 
!     charge density.
!
!     Output is written in format suitable for Data Explorer (DX).
!     Mesh is defined by the supercell size and the FFT parameters
!     used to generate WFN_fi/WFNq_fi.
!
!     input option 'restrict_kpoints' reduces the sum over k-points
!     above to a sum over the ones that give most of the contribution
!     to the norm of eigenvectors. This is handy if there are many k-points
!     but only a few of them give sizable contribution. 
!     This feature only works for the first state plotted.
!
!     input: WFN_fi WFNq_fi       wavefunction files
!            eigenvectors         output from BSE/absorption code
!            plotxctdens.inp          input parameters
!            (seedname.chk         checkpoint from Wannier90 code)
!
!     output: xct.<ie>_s<spin>_tfull.a3Dr      real-space data file full delta rho(r)
!     output: xct.<ie>_s<spin>_tdiag.a3Dr      real-space data file diag delta rho(r)
!
!     compilation: make plotxctdens
!
!     Dependencies: uses FFTW library
!
!     Developers:
!
!        Peter Doak, Berkeley CA, pdoak@lbl.gov
!        Sahar Sharifzadeh, Berkeley CA, ssharifzadeh@lbl.gov
!           
!=================================================================================

#include "f_defs.h"

program plotxctdens

  use global_m
  use fftw_m
  use genwf_m
  use inreaddens_m
  use misc_m
  use sort_m
  use write_xct_dens_m
  use plotxct_common_m
  use fullbz_m
  use timing_m, only: timing => extra_timing
  use input_m
  use input_q_m
  use write_program_header_m
  implicit none

  type (crystal) :: crys
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (xctinfo) :: xct
  type (grid) :: kg,kgq
  type (wavefunction) :: wfnc,wfnv
  type (work_genwf) :: work, workq
  type (int_wavefunction) :: intwfn
  type (plotxct_t) :: pxct
!PD: ixct is state not spin
  integer :: ii,jj,ikq,ik,ikt,ic,iv,ixct,ikcvs,i1,i2,i3,ir1,ir2,ir3,ip1,ip2,ip3
!PD
  integer :: ivp, icp, ikcvsp
  integer :: ncount,ntim,nmat,iunit_c,iunit_v
  integer :: nstate
  integer :: nw(3),nfft(3),step !no more subsampling
  real(DP) :: sum,sumf,scale
  real(DP) :: phel,phhole,kk(3),kkq(3),rel(3)
  real(DP) :: tsec(2)
!PD
  real(DP), allocatable  :: normHole(:),normElec(:) !< (nv), (nc)

  complex(DPC) :: wfnvint
  character :: filename*20 
  character*16, allocatable :: routnam(:)
  character :: AsvckString*256
  integer, allocatable :: indexq(:),index_k(:)
  real(DP), allocatable :: wwk(:),kgr(:,:)
  complex(DPC) :: Abs2Asvck
  complex(DPC) :: phfel, phfhole !< phase factors so we don`t recalculate these many times 
  complex(DPC) :: diagPoint, crossPoint
  complex(DPC), allocatable :: fel(:,:,:,:) !< (nfft1, nfft2, nfft3, ncb)
!SS,PD:
  complex(DPC), allocatable :: fhole(:,:,:,:) !< (nfft1, nfft2, nfft3,nvb)
!  complex(DPC), allocatable :: totalHole(:), totalElec(:) ! Do we get 1 electron and 1 hole total?

  complex(DPC), allocatable :: coefficients(:,:) !< (nvb, pxct%nsuper)
!  complex(DPC), allocatable :: ucfftv(:,:,:) !< (nfft1, nfft2, nfft3)
  complex(DPC), allocatable :: scfft(:,:,:,:) !< (exciton_state,nw1, nw2, nw3)
! PD:
  complex(DPC), allocatable :: scfftDiag(:,:,:,:) !< (excition_state,nw1, nw2, nw3)
!  SCALAR, allocatable :: ucfftv_sum(:,:,:), ucfftc_sum(:,:,:) !< (nfft1, nfft2, nfft3)
  real(DP), allocatable :: hole_prob(:,:,:), elec_prob(:,:,:) !< (nfft1, nfft2, nfft3)

  integer  :: nmpinode
  real(DP) :: mem,rmem,scalarsize

  SCALAR, allocatable :: Asvck(:,:)

  call peinfo_init()

  call timing%init()
  call timing%start(timing%total)

!  call set_debug_level(1)

  call write_program_header('PlotXctDens', .false.)

  !---------------------------------
  ! Read plotxctdens.inp

  call inreaddens(xct,pxct)

  !---------------------------------
  ! Read eigenvectors

  call logit('plotxct: reading file eigenvectors')
  if(peinf%inode.eq.0) then
    call open_file(unit=10,file='eigenvectors',form='unformatted',status='old')
    read(10) xct%nspin
    read(10) xct%nvb_fi
    read(10) xct%ncb_fi
    read(10) xct%nkpt_fi
    SAFE_ALLOCATE(kgr, (3,xct%nkpt_fi))
    read(10) ((kgr(ii,ik),ii=1,3),ik=1,xct%nkpt_fi)
  endif

  if(pxct%ispin == 2 .and. xct%nspin == 1) then
    call die("plot_spin = 2 but calculation is not spin_polarized", only_root_writes = .true.)
  endif

#ifdef MPI
  call MPI_BCAST(xct%nspin, 1, MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
  call MPI_BCAST(xct%nvb_fi, 1, MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
  call MPI_BCAST(xct%ncb_fi, 1, MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
  call MPI_BCAST(xct%nkpt_fi, 1, MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
  if (peinf%inode.ne.0) then
    SAFE_ALLOCATE(kgr, (3,xct%nkpt_fi))
  endif
  call MPI_BCAST(kgr, 3*xct%nkpt_fi,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif

  !-------------------------------
  ! Read selected state
  if(pxct%pstateLow.eq.0) then
    pxct%pstateLow=1
  endif
  if(pxct%pstateHigh.eq.0) then
    pxct%pstateHigh=xct%nvb_fi*xct%ncb_fi
  endif
  pxct%offset = pxct%pstateLow-1
  nstate = pxct%pstateHigh-pxct%offset
  if(peinf%inode.eq.0) then
    write(6,*) "processing ", nstate, " excitons"
  endif
  nmat= xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin
  SAFE_ALLOCATE(Asvck, (nstate,nmat))
  SAFE_ALLOCATE(pxct%e_Ses, (nstate))
  if(peinf%inode.eq.0) then
    if(pxct%pstateLow.gt.1) then
      do ii=1,pxct%offset
        read(10)
        read(10)
      enddo
    endif
    do ii=pxct%pstateLow,pxct%pstateHigh
      read(10) pxct%e_Ses(ii-pxct%offset)
      read(10) (Asvck(ii-pxct%offset,jj),jj=1,nmat)
      write(6,*) ' Reading state # ',ii
      write(6,*) ' Excitation energy = ',pxct%e_Ses(ii-pxct%offset),' eV'
    enddo
  endif

#ifdef MPI
  call MPI_BCAST(Asvck,(nstate*nmat),MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#endif
  call logit('plotxct: done reading file eigenvectors')

  !-------------------------------
  ! Restrict the # of k-points to xct%nn
  !   not sure if this will work well for CD
  !   we only check and sort based on the the first state
  !   and this should be refactored into an _m shared function
  !   for plotxct and plotxctdens

  if (xct%nn.eq.0) xct%nn=xct%nkpt_fi
  SAFE_ALLOCATE(wwk, (xct%nkpt_fi))
  wwk=0.d0
  do ik=1,xct%nkpt_fi
    do iv=1,xct%nvb_fi
      do ic=1,xct%ncb_fi
        ikcvs = bse_index(ik, ic, iv, pxct%ispin, xct)
        wwk(ik) = wwk(ik) + abs(Asvck(1,ikcvs))**2 !checking on the first exciton
      enddo
    enddo
  enddo

  SAFE_ALLOCATE(index_k, (xct%nkpt_fi))
  index_k=0

  if(xct%nkpt_fi.gt.1 .and. nstate.gt.1 .and. xct%nn/=xct%nkpt_fi) then
    call die("ERROR: kpoint reduction will produce incorrect for states past first plotted", only_root_writes = .true.)
  endif
    
  call sortrx_D(xct%nkpt_fi, wwk, index_k)
  ! cannot use gvec here, not read in yet!
  do ii=1,xct%nkpt_fi/2
    ik=index_k(ii)
    index_k(ii)=index_k(xct%nkpt_fi+1-ii)
    index_k(xct%nkpt_fi+1-ii)=ik
  enddo

  ! Index_k reversed (check if it works!)

  if(peinf%inode.eq.0) then
    write(6,*) ' Reducing k-points to ',xct%nn
    sumf=0.d0
    do ik=1,xct%nkpt_fi
      sumf=sumf + wwk(ik)
    enddo
    sum=0.d0
    do ii=1,xct%nn
      sum=sum + wwk(index_k(ii))
    enddo
    write(6,*) ' Wfn norm (restricted sum): ',sum
    write(6,*) ' Wfn norm (full): ',sumf
  endif

  SAFE_DEALLOCATE(wwk)

  !------------------------------
  ! Read WFN_fi

  call timing%start(timing%input)
  call input(crys,gvec,kg,syms,xct,kgr,index_k,pxct)

  if(peinf%inode.eq.0) then
    write(6,*) ' # valence bands = ',xct%nvb_fi
    write(6,*) ' # cond. bands   = ',xct%ncb_fi
    write(6,*) ' Plotting spin component ',pxct%ispin
    write(6,*) ' lattice box = ',pxct%nsuper(:)
  endif
  call timing%stop(timing%input)
  SAFE_DEALLOCATE(kgr)

  !----------------------------------
  ! Read WFNq_fi

  SAFE_ALLOCATE(indexq, (kg%nf))
  call timing%start(timing%input_q)
  call input_q(crys,gvec,kg,kgq,syms,xct,indexq,pxct)
  call timing%stop(timing%input_q)

  !----------------------------------
  ! Compute size of FFT box we need
  call setup_FFT_sizes(gvec%FFTgrid,nfft,scale)

  nw(:) = nfft(:)*pxct%nsuper(:)

  !  Allocate FFT boxes

!  SAFE_ALLOCATE(ucfftv_sum, (nfft(1),nfft(2),nfft(3)))
!  SAFE_ALLOCATE(hole_prob, (nfft(1),nfft(2),nfft(3)))
!  hole_prob(:,:,:) = 0d0
!  SAFE_ALLOCATE(ucfftc_sum, (nfft(1),nfft(2),nfft(3)))
!  SAFE_ALLOCATE(elec_prob, (nfft(1),nfft(2),nfft(3)))
!  elec_prob(:,:,:) = 0d0
  SAFE_ALLOCATE(fel, (nfft(1),nfft(2),nfft(3),xct%ncb_fi))
  SAFE_ALLOCATE(fhole, (nfft(1),nfft(2),nfft(3),xct%nvb_fi))

  !FHJ : effective supercell taking into consideration the subsampling

  call procmem(mem,nmpinode)

  if(peinf%inode.eq.0) then
    write(6,998) mem / 1024**2
    write(7,998) mem / 1024**2
  endif
998 format(1x,'Memory available:',f10.1,1x,'MB per PE')

    ! Determine the available memory
    ! from epsilon (candidate for refactor?)

  ! required memory
  ! think this is in bytes
  scalarsize = sizeof_scalar()
  rmem=0.0d0
  ! the grid full and diag
  !rmem = rmem + dble(nw(1)) * dble(nw(2)) * dble(nw(3)) * 8 * 2 * nstate
  ! Asvck for all the states
  !rmem = rmem + nstate * nmat * 8
  ! pxct
  !rmem = rmem + nstate * 8
  ! kpoints
  !rmem = rmem + 5 * xct%nkpt_fi
  ! fhole and felec
  rmem = rmem + dble(nw(1)) * dble(nw(2)) * dble(nw(3)) * scalarsize * nstate * 2
  if(peinf%inode.eq.0) then
    write(6,989) rmem/1024.0d0**2
    write(7,989) rmem/1024.0d0**2
  endif
989 format(1x,'Memory required for execution:',f7.1,1x,'MB per PE')
  if(peinf%inode.eq.0 .and. rmem > mem) then
    write(6,*) "mem needed is too big"
  endif
  SAFE_ALLOCATE(scfft, (nstate,nw(1),nw(2),nw(3)))
  scfft=0.d0
  SAFE_ALLOCATE(scfftDiag, (nstate,nw(1),nw(2),nw(3)))
  scfftDiag=0.d0
  ! SAFE_ALLOCATE(totalHole, (nstate))
  ! totalHole=0.d0
  ! SAFE_ALLOCATE(totalElec, (nstate))
  ! totalElec=0.d0
  
  if (peinf%inode.eq.0) then
    write(6,'(a,3i6)') ' FFT box size = ',nfft
    write(6,'(a,3i6)') ' Supercell grid = ',nw
!    write(6,'(a,3i6)') ' Effective supercell grid = ',nw_eff
  endif

  step = max(1,(peinf%nkpe/20))

    !---------------------------------
  do ikt=1,peinf%ikt(peinf%inode+1)
    ik=peinf%ik(peinf%inode+1,ikt)
    ikq=indexq(ik)
    call timing%start(timing%genwf)
    call genwf(crys,gvec,kg,syms,wfnc,ik,ik,xct%nspin,xct%ncb_fi,&
               work,intwfn,0,is_cond=.true.) ! do not use intwfn
    call timing%stop(timing%genwf)


    call timing%start(timing%genwf_q)
    call genwf(crys,gvec,kgq,syms,wfnv,ik,ikq,xct%nspin,xct%nvb_fi,&
               workq,intwfn,0,is_cond=.false.) ! do not use intwfn
    call timing%stop(timing%genwf_q)


    kk(:) = kg%f(:,ik)
    kkq(:) = kgq%f(:,ikq) 

    fhole=0.d0

    do iv=1,xct%nvb_fi
        call put_into_fftbox(wfnv%ng,wfnv%cg(1:,iv,pxct%ispin),gvec%components,wfnv%isort,fhole(:,:,:,iv),nfft)
        call do_FFT(fhole(:,:,:,iv),nfft,1)
    enddo

    fel=0.d0

    !---------------------------------
    ! Calculate all needed conduction wavefunctions in the unit cell

    do ic=1,xct%ncb_fi
      call put_into_fftbox(wfnc%ng,wfnc%cg(1:,ic,pxct%ispin),gvec%components,wfnc%isort,fel(:,:,:,ic),nfft)
      call do_FFT(fel(:,:,:,ic),nfft,1)
    enddo !ic

    !---------------------------------
    ! Calculate contribution from each kcvs quadruplet in the excited state

    SAFE_DEALLOCATE_P(wfnc%cg)
    SAFE_DEALLOCATE_P(wfnc%isort)
    SAFE_DEALLOCATE_P(wfnv%cg)
    SAFE_DEALLOCATE_P(wfnv%isort)

    ! rel(:) = 0.1d0
    ! phel = 2.0d0*PI_D*DOT_PRODUCT( kk, rel )
    ! phhole = 2.0d0*PI_D*DOT_PRODUCT( kkq, rel )
    ! phfel = CMPLX(cos(phel),sin(phel))
    ! phfhole = CMPLX(cos(phhole),sin(phhole))
    ! write(6,*) "phases and factors",phel,phhole,phfel,phfhole

    call procmem(mem,nmpinode)

    if(peinf%inode.eq.0) then
      write(6,998) mem / 1024**2
      write(7,998) mem / 1024**2
    endif

    call timing%start(timing%summing)
    do ixct=pxct%pstateLow-pxct%offset,pxct%pstateHigh-pxct%offset      

      !$omp parallel default(none) shared(scfft,scfftDiag,nw,pxct,ixct,Asvck,nfft,fel,fhole,xct,index_k,ik,kk,kkq) &
      !$omp private(ikcvs,iv,ic,ir3,ir2,ir1,Abs2Asvck,rel,diagPoint,ikcvsp,phel,phhole,phfel,phfhole)
      !$omp do        
      do i3=1,nw(3)
        rel(3) = dble(i3)/dble(nfft(3))
        ir3 = mod(i3-1, nfft(3)) + 1
        do i2=1,nw(2)
          rel(2) = dble(i2)/dble(nfft(2))
          ir2 = mod(i2-1, nfft(2)) + 1
          do i1=1,nw(1)
            rel(1) = dble(i1)/dble(nfft(1))
            ir1 = mod(i1-1, nfft(1)) + 1
            phel = 2.0d0*PI_D*DOT_PRODUCT( kk, rel )
            phhole = 2.0d0*PI_D*DOT_PRODUCT( kkq, rel )
            phfel = CMPLX(cos(phel),sin(phel))
            phfhole = CMPLX(cos(phhole),sin(phhole))
            do ic=1,xct%ncb_fi
               do iv=1,xct%nvb_fi
                ikcvs = bse_index(index_k(ik), ic, iv, pxct%ispin, xct)
                Abs2Asvck = Asvck(ixct,ikcvs) * MYCONJG(Asvck(ixct,ikcvs))
                diagPoint = Abs2Asvck * &
                 (CONJG(fel(ir1,ir2,ir3,ic)*phfel)*fel(ir1,ir2,ir3,ic)*phfel - &
                 CONJG(fhole(ir1,ir2,ir3,iv)*phfhole)*fhole(ir1,ir2,ir3,iv)*phfhole)
                scfft(ixct,i1,i2,i3) = scfft(ixct,i1,i2,i3) + diagPoint
                scfftDiag(ixct,i1,i2,i3) = scfftDiag(ixct,i1,i2,i3) + diagPoint
                do ivp=1,xct%nvb_fi
                  if (ivp /= iv) then
                    ikcvsp = bse_index(index_k(ik), ic, ivp, pxct%ispin, xct)
                    scfft(ixct,i1,i2,i3) = scfft(ixct,i1,i2,i3) - &
                     Asvck(ixct,ikcvs) * MYCONJG(Asvck(ixct,ikcvsp)) * &
                     fhole(ir1,ir2,ir3,iv)*phfhole * CONJG(fhole(ir1,ir2,ir3,ivp)*phfhole)
                  endif
                enddo ! ivp
                do icp=1,xct%ncb_fi
                  if (icp /= ic) then
                    ikcvsp = bse_index(index_k(ik), icp, iv, pxct%ispin, xct)
                    scfft(ixct,i1,i2,i3) = scfft(ixct,i1,i2,i3) + &
                     Asvck(ixct,ikcvs) * MYCONJG(Asvck(ixct,ikcvsp)) * &
                     CONJG(fel(ir1,ir2,ir3,ic)*phfel) * fel(ir1,ir2,ir3,icp)*phfel
                  endif
                enddo ! icp
              enddo ! iv
            enddo ! ic

          enddo ! i1
        enddo ! i2
      enddo ! i3
      !$omp end do
      !$omp end parallel
    enddo ! ixct
    call timing%stop(timing%summing)

    if (peinf%inode.eq.0.and.mod(ikt,step).eq.0) &
      write(6,*)'PE #0 working at k-point # ',ikt,' out of ',peinf%ikt(1)

  enddo ! ik

  call dealloc_grid(kg)
  call dealloc_grid(kgq)

  call destroy_fftw_plans()

  ! typedefs initializes all of these ikolds to 0
  if(work%ikold.ne.0) then
    SAFE_DEALLOCATE_P(work%cg)
    SAFE_DEALLOCATE_P(work%ph)
    SAFE_DEALLOCATE_P(work%ind)
    SAFE_DEALLOCATE_P(work%isort)
  endif
  if(workq%ikold.ne.0) then
    SAFE_DEALLOCATE_P(workq%cg)
    SAFE_DEALLOCATE_P(workq%ph)
    SAFE_DEALLOCATE_P(workq%ind)
    SAFE_DEALLOCATE_P(workq%isort)
  endif

  SAFE_DEALLOCATE(fel)
  SAFE_DEALLOCATE(fhole)
  SAFE_DEALLOCATE(indexq)
  SAFE_DEALLOCATE(index_k)

  call timing%start(timing%gather)
  !-----------------------------
  ! Add up information on supercell 'scfft' and 'scfftdiag' (dest = PE 0)
  call gather_data()
  call timing%stop(timing%gather)


  !-----------------------------------
  ! Write output (only PE 0)

  if (peinf%inode.eq.0) call write_xct_dens(pxct, crys, nfft, 'full',  scfft)
  if (peinf%inode.eq.0) call write_xct_dens(pxct, crys, nfft, 'diag',  scfftDiag)

  !------------------------------------
  ! Time accounting
  call timing%stop(timing%total)
  call timing%print()

  call write_memory_usage()

  SAFE_DEALLOCATE(scfft)
  SAFE_DEALLOCATE(scfftDiag)

  ! SAFE_DEALLOCATE(normHole)
  ! SAFE_DEALLOCATE(normElec)

  ! SAFE_DEALLOCATE(totalHole)
  ! SAFE_DEALLOCATE(totalElec)

  iunit_v = 128+2*(peinf%inode-1)+2
  write(filename,'(a,i4.4)') 'INT_VWFNQ_', peinf%inode
  
  call open_file(iunit_v, filename, status = 'old')
  call close_file(iunit_v, delete = .true.) ! files INT_VWFNQ*

  iunit_c = 128+2*(peinf%inode-1)+1
  write(filename,'(a,i4.4)') 'INT_CWFN_', peinf%inode
  call open_file(iunit_c, filename, status = 'old')
  call close_file(iunit_c, delete = .true.) ! files INT_CWFN_*

#ifdef MPI
  call MPI_FINALIZE(mpierr)
#endif

contains

  !> Gather the scfft array to node 0 by slicing it in the larger axis and
  !! performing multiple MPI_Gather`s. We could do only one large MPI_Gather
  !! without temporary buffers if we could use in-place MPI routines...
  subroutine gather_data()

    complex(DPC), allocatable :: buf_send(:,:), buf_recv(:,:)
! PD second buffer for Diag scfft
    complex(DPC), allocatable :: buf_sendDiag(:,:), buf_recvDiag(:,:)
    integer :: slice_idx  !< the index that we will use to break the scfft array
    integer :: aux_idx(2) !< the other two indices
    integer :: islice

    PUSH_SUB(gather_data)

    if (nw(3)>=nw(1).and.nw(3)>=nw(2)) then
      slice_idx = 3
      aux_idx = [1,2]
    else if (nw(2)>=nw(1).and.nw(2)>=nw(3)) then
      slice_idx = 2
      aux_idx = [1,3]
    else
      slice_idx = 1
      aux_idx = [2,3]
    endif

    SAFE_ALLOCATE(buf_send, (nw(aux_idx(1)),nw(aux_idx(2))))
    SAFE_ALLOCATE(buf_recv, (nw(aux_idx(1)),nw(aux_idx(2))))
! PD
    SAFE_ALLOCATE(buf_sendDiag, (nw(aux_idx(1)),nw(aux_idx(2))))
    SAFE_ALLOCATE(buf_recvDiag, (nw(aux_idx(1)),nw(aux_idx(2))))

    do ixct=pxct%pstateLow-pxct%offset,pxct%pstateHigh-pxct%offset
      if (peinf%inode==0) then
        write(6,'(/,1x,a,i4,a,i5,a)') 'Gathering data for slice', ixct+pxct%offset, &
         ': breaking array in ',nw(slice_idx),' arrays containing'
        write(6,'(1x,i8,a,/)') SIZE(buf_recv), ' elements each.'
      endif

      do islice=1,nw(slice_idx)

        if (slice_idx==3) then
          buf_send(:,:) = scfft(ixct,:,:,islice)
          buf_sendDiag(:,:) = scfftDiag(ixct,:,:,islice)
        else if (slice_idx==2) then
          buf_send(:,:) = scfft(ixct,:,islice,:)
          buf_sendDiag(:,:) = scfftDiag(ixct,:,islice,:)
        else
          buf_send(:,:) = scfft(ixct,islice,:,:)
          buf_sendDiag(:,:) = scfftDiag(ixct,islice,:,:)
        endif

        buf_recv = ZERO

#ifdef MPI
        !FHJ: only node 0 wants this, so don`t Allreduce!
        call MPI_Reduce(buf_send(1,1), buf_recv(1,1), SIZE(buf_recv), &
         MPI_COMPLEX_DPC, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
        call MPI_Reduce(buf_sendDiag(1,1), buf_recvDiag(1,1), SIZE(buf_recv), &
         MPI_COMPLEX_DPC, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
#else
        buf_recv = buf_send
        buf_recvDiag = buf_sendDiag
#endif

        if (peinf%inode==0) then
          if (slice_idx==3) then
            scfft(ixct,:,:,islice) = buf_recv(:,:)
            scfftDiag(ixct,:,:,islice) = buf_recvDiag(:,:)
          else if (slice_idx==2) then
            scfft(ixct,:,islice,:) = buf_recv(:,:)
            scfftDiag(ixct,:,islice,:) = buf_recvDiag(:,:)
          else
            scfft(ixct,islice,:,:) = buf_recv(:,:)
            scfftDiag(ixct,islice,:,:) = buf_recvDiag(:,:)
          endif
        endif

      enddo
    enddo

    SAFE_DEALLOCATE(buf_send)
    SAFE_DEALLOCATE(buf_recv)
    SAFE_DEALLOCATE(buf_sendDiag)
    SAFE_DEALLOCATE(buf_recvDiag)

! PD: apply volume factor
    scfft = scfft/crys%celvol
    scfftDiag = scfftDiag/crys%celvol

    POP_SUB(gather_data)

  end subroutine gather_data

end program plotxctdens
