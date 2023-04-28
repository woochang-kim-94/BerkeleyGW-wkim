!==========================================================================
!
! Routines:
!
! (1) plotxct(main)     Originally By MLT       Last Modified 6/2008 FJR
!
!     Original developers:
!
!        Murilo Tiago, Berkeley CA, mtiago@civet.berkeley.edu
!           original creator
!
!        Sohrab Ismail-Beigi, Berkeley CA, sohrab@civet.berkeley.edu
!           FFTW interface routines
!
!        Filipe Ribeiro, Berkeley CA, fribeiro@alum.berkeley.edu
!           single output xct.<ie>.r.asc
!
!=================================================================================

#include "f_defs.h"

program plotxct

#ifdef HDF5
  use hdf5
#endif

  use global_m
  use fftw_m
  use genwf_m
  use susymmetries_m
  use random_m
  use inread_m
  use misc_m
  use read_wannier_m
  use sort_m
  use write_xct_m
  use plotxct_common_m
  use references_m
  use fullbz_m
  use timing_m, only: timing => extra_timing
  use input_m
  use input_q_m
  use interpol_m
  use write_program_header_m
  use evecs_m
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
  type (evecs_t) :: evecs

  integer :: ii,ikq,ik,ikt,ic,iv,ikcvs,i1,i2,i3,ir1,ir2,ir3,ip1,ip2,ip3,jsp,jspp
  integer :: iunit_c,iunit_v
  integer :: error
  integer :: nw(3),nw_eff(3),nfft(3),step
  real(DP) :: sum,sumf,scale
  real(DP) :: phhole,phel,kk(3),rel(3)
  complex(DPC) :: wfnvint
  character :: filename*20 
  integer, allocatable :: indexq(:),index_k(:)
  real(DP), allocatable :: wwk(:),kgr(:,:)

  complex(DPC), allocatable :: fel(:,:,:,:) !< (nfft1, nfft2, nfft3, ncb)
  complex(DPC), allocatable :: fhole(:) !< (nvb)
  complex(DPC), allocatable :: coefficients(:,:) !< (nvb, pxct%nsuper)
  complex(DPC), allocatable :: ucfftv(:,:,:) !< (nfft1, nfft2, nfft3)
  complex(DPC), allocatable :: scfft(:,:,:) !< (nw_eff1, nw_eff2, nw_eff3)
  SCALAR, allocatable :: ucfftv_sum(:,:,:), ucfftc_sum(:,:,:) !< (nfft1, nfft2, nfft3)
  real(DP), allocatable :: hole_prob(:,:,:), elec_prob(:,:,:) !< (nfft1, nfft2, nfft3)

#ifdef HDF5
  call h5open_f(error)
#endif

  call peinfo_init()

  call timing%init()
  call timing%start(timing%total)

  call write_program_header('PlotXct', .false.)

  !---------------------------------
  ! Read plotxct.inp

  call inread(xct,pxct)

  !---------------------------------
  ! Read eigenvectors

  call logit('plotxct: reading eigenvectors')

  ! Read the file, allocate, and broadcast to all nodes
  if(peinf%inode.eq.0) then
    write(6,'(1x,a,i0)') 'Reading state # ', pxct%pstate
  endif
  call evecs%read_file(neig_read=1, ieig_offset=pxct%pstate-1)
  call evecs%print_out_dimensions(6)

  xct%nspin = evecs%ns
  xct%nvb_fi = evecs%nv
  xct%ncb_fi = evecs%nc
  xct%nkpt_fi = evecs%nk
  SAFE_ALLOCATE(kgr, (3,xct%nkpt_fi))
  kgr(:,:) = evecs%kpts(:,:)

  pxct%e_S = evecs%evals(1)

  if(pxct%ispin == 2 .and. xct%nspin == 1) then
    call die("plot_spin = 2 but calculation is not spin_polarized", only_root_writes = .true.)
  endif

  if(peinf%inode.eq.0) then
    write(6,'(1x,a,f0.9,a)') 'Excitation energy: ', pxct%e_S, ' eV'
  endif

  call logit('plotxct: done reading file eigenvectors')

  !-------------------------------
  ! Restrict the # of k-points to xct%nn

  if (xct%nn.eq.0) xct%nn=xct%nkpt_fi
  SAFE_ALLOCATE(wwk, (xct%nkpt_fi))
  wwk=0.d0
  do ik=1,xct%nkpt_fi
    do iv=1,xct%nvb_fi
      do ic=1,xct%ncb_fi
        wwk(ik) = wwk(ik) + abs(evecs%Avc(pxct%ispin,iv,ic,ik,1))**2
      enddo
    enddo
  enddo

  SAFE_ALLOCATE(index_k, (xct%nkpt_fi))
  index_k=0
  call sortrx(xct%nkpt_fi, wwk, index_k)
  ! cannot use gvec here, not read in yet!
  do ii=1,xct%nkpt_fi/2
    ik=index_k(ii)
    index_k(ii)=index_k(xct%nkpt_fi+1-ii)
    index_k(xct%nkpt_fi+1-ii)=ik
  enddo

  ! index_k reversed (check if it works!)

  if(peinf%inode.eq.0) then
    write(6,'(1x,a,i0)') 'Reducing k-points to ', xct%nn
    sumf=0.d0
    do ik=1,xct%nkpt_fi
      sumf=sumf + wwk(ik)
    enddo
    sum=0.d0
    do ii=1,xct%nn
      sum=sum + wwk(index_k(ii))
    enddo
    write(6,'(1x,a)') 'Norm squared of the exciton expansion coefficient:'
    write(6,'(1x,a,f0.15)') '- Restricted sum: ', sum
    write(6,'(1x,a,f0.15)') '- Full sum: ', sumf
  endif

  SAFE_DEALLOCATE(wwk)

  !------------------------------
  ! Read WFN_fi

  call timing%start(timing%input)
  call input(crys,gvec,kg,syms,xct,kgr,index_k,pxct)

  if(peinf%inode.eq.0) then
    write(6,'(/1x,a)') 'Plotting parameters:'
    write(6,'(1x,a,i0)') '- Number of valence bands: ', xct%nvb_fi
    write(6,'(1x,a,i0)') '- Number of cond. bands: ', xct%ncb_fi
    if(pxct%plot_hole) write(6,'(1x,a)') '- Plotting hole probability.'
    if(pxct%plot_electron) write(6,'(1x,a)') '- Plotting electron probability.'
    if(pxct%iwann == 0) then
      ! FHJ: note: position not in [0,1) range
      write(6,'(1x,a,3(1x,f0.6))') '- Coordinates of hole (supercell units):', pxct%rhole(:)
    else
      write(6,'(1x,a,i0)') '- Projecting onto Wannier function # ', pxct%iwann
    endif
    write(6,'(1x,a,i0)') '- Plotting spin component: ', pxct%ispin
    write(6,'(1x,a,3(1x,i0))') '- Unit cell repetitions / supercell size:', pxct%nsuper(:)
  endif
  call timing%stop(timing%input)
  SAFE_DEALLOCATE(kgr)

  !----------------------------------
  ! Read WFNq_fi

  SAFE_ALLOCATE(indexq, (kg%nf))
  call timing%start(timing%input_q)
  call input_q(crys,gvec,kg,kgq,syms,xct,indexq,pxct)
  call timing%stop(timing%input_q)

!#BEGIN_INTERNAL_ONLY
  if(pxct%iwann /= 0) then
    SAFE_ALLOCATE(coefficients, (xct%nvb_fi, kgq%nf))
    call read_wannier(pxct%seedname, crys, kgq, pxct%iwann, xct%nvb_fi, coefficients)
  endif
!#END_INTERNAL_ONLY

  !----------------------------------
  ! Compute size of FFT box we need
  call setup_FFT_sizes(gvec%FFTgrid,nfft,scale)

  ! nw defines the mesh in supercell. Should be defined so that
  ! nfft(ii)*pxct%nsuper(ii)/nw(ii) = integer for ii=1,3

  nw(:) = nfft(:)*pxct%nsuper(:)

  !  Allocate FFT boxes

  if(pxct%iwann == 0 .or. pxct%plot_hole) then
    SAFE_ALLOCATE(ucfftv, (nfft(1),nfft(2),nfft(3)))
  endif
  if(pxct%plot_hole) then
    SAFE_ALLOCATE(ucfftv_sum, (nfft(1),nfft(2),nfft(3)))
    SAFE_ALLOCATE(hole_prob, (nfft(1),nfft(2),nfft(3)))
    hole_prob(:,:,:) = 0d0
  endif
  if(pxct%plot_electron) then
    SAFE_ALLOCATE(ucfftc_sum, (nfft(1),nfft(2),nfft(3)))
    SAFE_ALLOCATE(elec_prob, (nfft(1),nfft(2),nfft(3)))
    elec_prob(:,:,:) = 0d0
  endif
  SAFE_ALLOCATE(fel, (nfft(1),nfft(2),nfft(3),xct%ncb_fi))
  SAFE_ALLOCATE(fhole, (xct%nvb_fi))

  !FHJ : effective supercell taking into consideration the subsampling
  nw_eff = (nw + pxct%downsample - 1) / pxct%downsample

  SAFE_ALLOCATE(scfft, (nw_eff(1),nw_eff(2),nw_eff(3)))
  scfft=0.d0
  if (peinf%inode.eq.0) then
    write(6,'(/1x,a)') 'Grids information:'
    write(6,'(1x,a,3(1x,i0))') '- FFT box size:', nfft
    write(6,'(1x,a,3(1x,i0))') '- Supercell grid:', nw
    write(6,'(1x,a,3(1x,i0))') '- Effective supercell grid:', nw_eff
  endif

  if (peinf%inode==0) write(6,'()')
  step = max(1,(peinf%nkpe/20))
  do ikt=1,peinf%ikt(peinf%inode+1)
    ik=peinf%ik(peinf%inode+1,ikt)
    ikq=indexq(ik)
    call timing%start(timing%genwf)
    call genwf(crys,gvec,kg,syms,wfnc,ik,ik,xct%nspin,xct%ncb_fi,&
               work,intwfn,0,is_cond=.true.) ! do not use intwfn
    call timing%stop(timing%genwf)

    if(pxct%iwann == 0 .or. pxct%plot_hole) then
      call timing%start(timing%genwf_q)
      call genwf(crys,gvec,kgq,syms,wfnv,ik,ikq,xct%nspin,xct%nvb_fi,&
                 workq,intwfn,0,is_cond=.false.) ! do not use intwfn
      call timing%stop(timing%genwf_q)
    endif

    phhole = 2.0d0 * PI_D * DOT_PRODUCT(kgq%f(:,ikq), pxct%rhole)
    kk(:) = kg%f(:,ik)

    !-------------------------------
    ! Calculate all needed valence wavefunctions in real space and
    ! store their value at hole position in fhole

    do iv=1,xct%nvb_fi
      if(pxct%iwann == 0) then
        call put_into_fftbox(wfnv%ng,wfnv%cg(1:,iv,pxct%ispin*pxct%hspinor),gvec%components,wfnv%isort,ucfftv,nfft)
        call do_FFT(ucfftv,nfft,1)
        call interpol(pxct%rhole,nfft,wfnvint,ucfftv)
        fhole(iv) = CONJG(wfnvint)*CMPLX(cos(phhole),-sin(phhole))
      else
        fhole(iv) = CONJG(coefficients(xct%nvb_fi - iv + 1,ikq))
      endif
    enddo ! iv

    if(pxct%plot_hole) then
      do ic=1,xct%ncb_fi
        ucfftv_sum(:,:,:) = ZERO
        do iv=1,xct%nvb_fi
          call put_into_fftbox(wfnv%ng,wfnv%cg(1:,iv,pxct%ispin*pxct%hspinor),gvec%components,wfnv%isort,ucfftv,nfft)
          call do_FFT(ucfftv,nfft,1)
          ucfftv_sum(:,:,:) = ucfftv_sum(:,:,:) + MYCONJG(evecs%Avc(pxct%ispin,iv,ic,ik,1)) * ucfftv(:,:,:)
        enddo
        ! this is accumulating over k-points too since we are inside that loop
        hole_prob(:,:,:) = hole_prob(:,:,:) + abs(ucfftv_sum(:,:,:))**2
      enddo
    endif

    fel=0.d0

    !---------------------------------
    ! Calculate all needed conduction wavefunctions in the unit cell

    do ic=1,xct%ncb_fi
      call put_into_fftbox(wfnc%ng,wfnc%cg(1:,ic,pxct%ispin*pxct%espinor),gvec%components,wfnc%isort,fel(:,:,:,ic),nfft)
      call do_FFT(fel(:,:,:,ic),nfft,1)
    enddo !ic

    !---------------------------------
    ! Calculate contribution from each kcvs quadruplet in the excited state

    do ic=1,xct%ncb_fi
      do iv=1,xct%nvb_fi

        !$omp parallel default(shared) private(i1,i2,i3,ip1,ip2,ip3,ir1,ir2,ir3,rel,phel)
        !$omp do schedule(dynamic)
        do i3=1,nw(3),pxct%downsample(3)
          ip3 = ((i3-1) / pxct%downsample(3)) + 1
          rel(3) = dble(i3-1)/dble(nfft(3))
          ir3 = mod(i3-1, nfft(3)) + 1
          ip2 = 0
          do i2=1,nw(2),pxct%downsample(2)
            ip2 = ip2 + 1
            rel(2) = dble(i2-1)/dble(nfft(2))
            ir2 = mod(i2-1, nfft(2)) + 1
            ip1 = 0
            do i1=1,nw(1),pxct%downsample(1)
              ip1 = ip1 + 1
              rel(1) = dble(i1-1)/dble(nfft(1))
              ir1 = mod(i1-1, nfft(1)) + 1

              phel = 2.0d0 * PI_D * DOT_PRODUCT(kk, rel)

              scfft(ip1,ip2,ip3) = scfft(ip1,ip2,ip3) + &
                evecs%Avc(pxct%ispin,iv,ic,index_k(ik),1) * fhole(iv) * fel(ir1,ir2,ir3,ic) * &
                CMPLX(cos(phel),sin(phel))

            enddo ! i1
          enddo ! i2
        enddo ! i3
        !$omp end do

        !$omp end parallel

      enddo ! iv
    enddo ! ic

    SAFE_DEALLOCATE_P(wfnc%cg)
    SAFE_DEALLOCATE_P(wfnc%isort)
    SAFE_DEALLOCATE_P(wfnv%cg)
    SAFE_DEALLOCATE_P(wfnv%isort)

    if (peinf%inode.eq.0.and.mod(ikt,step).eq.0) &
      write(6,'(1x,a,i0,a,i0)') 'PE #0 working at k-point # ', ikt, ' out of ', peinf%ikt(1)

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

  if(pxct%iwann == 0) then
    SAFE_DEALLOCATE(ucfftv)
  else
    SAFE_DEALLOCATE(coefficients)
  endif
  SAFE_DEALLOCATE(fel)
  SAFE_DEALLOCATE(indexq)
  SAFE_DEALLOCATE(index_k)

  !-----------------------------
  ! Add up information on supercell 'scfft' (dest = PE 0)

  call gather_data()

  !-----------------------------------
  ! Write output (only PE 0)

  if (peinf%inode.eq.0) call write_xct(pxct, crys, nfft, scfft)

  call require_reference(REF_DESLIPPE2012)
  call require_reference(REF_ROHLFING2000)
  call show_references()

  !------------------------------------
  ! Time accounting
  call timing%print()

  call evecs%free()

  SAFE_DEALLOCATE(scfft)

  iunit_v = 128+2*(peinf%inode-1)+2
  write(filename,'(a,i4.4)') 'INT_VWFNQ_', peinf%inode
  call open_file(iunit_v, filename, status = 'old')
  call close_file(iunit_v, delete = .true.) ! files INT_VWFNQ*

  iunit_c = 128+2*(peinf%inode-1)+1
  write(filename,'(a,i4.4)') 'INT_CWFN_', peinf%inode
  call open_file(iunit_c, filename, status = 'old')
  call close_file(iunit_c, delete = .true.) ! files INT_CWFN_*

  call write_memory_usage()

#ifdef HDF5
  call h5close_f(error)
#endif
  
#ifdef MPI
  call MPI_FINALIZE(mpierr)
#endif

contains

  !> Gather the scfft array to node 0 by slicing it in the larger axis and
  !! performing multiple MPI_Gather`s. We could do only one large MPI_Gather
  !! without temporary buffers if we could use in-place MPI routines...
  subroutine gather_data()

    complex(DPC), allocatable :: buf_send(:,:), buf_recv(:,:)
    integer :: slice_idx  !< the index that we will use to break the scfft array
    integer :: aux_idx(2) !< the other two indices
    integer :: islice

    PUSH_SUB(gather_data)

    if (nw_eff(3)>=nw_eff(1).and.nw_eff(3)>=nw_eff(2)) then
      slice_idx = 3
      aux_idx = [1,2]
    else if (nw_eff(2)>=nw_eff(1).and.nw_eff(2)>=nw_eff(3)) then
      slice_idx = 2
      aux_idx = [1,3]
    else
      slice_idx = 1
      aux_idx = [2,3]
    endif

    SAFE_ALLOCATE(buf_send, (nw_eff(aux_idx(1)),nw_eff(aux_idx(2))))
    SAFE_ALLOCATE(buf_recv, (nw_eff(aux_idx(1)),nw_eff(aux_idx(2))))

    if (peinf%inode==0) then
      write(6,'(/,1x,a,2(i0,a))') 'Gathering data: breaking array in ', &
        nw_eff(slice_idx), ' arrays containing ', SIZE(buf_recv), ' elements each.'
    endif

    do islice=1,nw_eff(slice_idx)

      if (slice_idx==3) then
        buf_send(:,:) = scfft(:,:,islice)
      else if (slice_idx==2) then
        buf_send(:,:) = scfft(:,islice,:)
      else
        buf_send(:,:) = scfft(islice,:,:)
      endif

      buf_recv = ZERO

#ifdef MPI
      !FHJ: only node 0 wants this, so don`t Allreduce!
      call MPI_Reduce(buf_send(1,1), buf_recv(1,1), SIZE(buf_recv), &
        MPI_COMPLEX_DPC, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
#else
      buf_recv = buf_send
#endif

      if (peinf%inode==0) then
        if (slice_idx==3) then
          scfft(:,:,islice) = buf_recv(:,:)
        else if (slice_idx==2) then
          scfft(:,islice,:) = buf_recv(:,:)
        else
          scfft(islice,:,:) = buf_recv(:,:)
        endif
      endif

    enddo

    SAFE_DEALLOCATE(buf_send)
    SAFE_DEALLOCATE(buf_recv)

    POP_SUB(gather_data)

  end subroutine gather_data

end program plotxct
