!==============================================================================
!
! Modules:
!
! chi_summation_m                               Last Modified: 04/19/2012 (FHJ)
!
!   This module creates the (distributed) polarizability matrix chi by summing
!   the (distributed) pol%gme matrices. There are routines that communicate the
!   gme matrix using either matrix or element partitioning scheme.
!
!==============================================================================

#include "f_defs.h"

module chi_summation_m

  use accel_linalg_m
  use accel_memory_m
  use algos_epsilon_m
  use global_m
  use blas_m
  use mtxelmultiply_m
  use scalapack_m
  use lin_denominator_m
  use io_utils_m
  use vcoul_generator_m
  use elpa_interface_m
  use misc_m, only : procmem

  use timing_m, only: timing => epsilon_timing

  implicit none

  private

  public :: create_chi_summator, free_chi_summator,&
    chi_summation_comm_matrix, chi_summation_comm_elements
#if defined MPI && defined USESCALAPACK
  public :: chi_summation_sub_trunc
#endif

  !> FHJ: the chi_summator "object"
  type chi_summator_t
    real(DP) :: fact
    !> DVF: below are some temporary buffers needed for the chi summation. They are
    !! described in detail in this comment.
    !! gme = g-vector matrix element
    !! gmetempX where X = n,r,c are temporary arrays for static calculations
    !! n = normalized by the proper factor used in BGW
    !! r = row, meaning the matrix elements for all bands (nv*nc*nk) that the proc owns
    !! for the G-vectors in the rows of chi currently being computed
    !! c = column, the matrix elements for all bands (nv*nc*nk) that the proc owns
    !! for the G`-vectors in the rows of chi currently being computed
    !! the RDyn arrays are needed for full-frequency (FF) calculations, real and complex
    !! r2 is used in FF with matrix communication because you have a different denominator for
    !! each frequency. Only the r2 array (not the r) array is used for element communication
    !! the denominators are built into the gme`s for static calculations
    !! eden arrays hold the energy denominators for FF
    !! chilocal holds the contribution of a given processor to the GG` chunk of epsilon
    !! being computed
    !! chilocal2 holds the MPI reduced GG` chunk of epsilon being computed
    SCALAR, allocatable :: chilocal(:,:)
    SCALAR, allocatable :: chilocal2(:,:,:)
    complex(DPC), allocatable :: chilocalRDyn(:,:,:)
    complex(DPC), allocatable :: chilocal2RDyn(:,:,:,:)
    complex(DPC), allocatable :: chiRDyntmp(:)
    complex(DPC), allocatable :: chiRDynOld(:,:,:,:)
    SCALAR, allocatable :: gmetempr(:,:),gmetempc(:,:)
    SCALAR, allocatable :: gmetempn(:)
    complex(DPC), allocatable :: gmeRDyntempn(:)
    complex(DPC), allocatable :: gmeRDyntempr2(:,:,:)
    complex(DPC), allocatable :: gmeRDyntempr3(:,:)
    complex(DPC), allocatable :: gmeRDyntempc(:,:)
    complex(DPC), allocatable :: gmeRDyntempcs(:,:)
    integer, allocatable :: deltaCount(:,:)
    integer, allocatable :: deltaCountReduce(:,:)
  end type chi_summator_t

  public :: chi_summator_t

contains

  subroutine create_chi_summator(this, pol, scal, fact, nspin)
    type(chi_summator_t), intent(INOUT) :: this !<the chi_summator_t object
    type(polarizability), intent(IN) :: pol
    type(scalapack), intent(IN) :: scal
    real(DP), intent(IN) :: fact
    integer, intent(IN) :: nspin

    PUSH_SUB(create_chi_summator)

    this%fact = fact

    if (pol%freq_dep .eq. 0) then
      SAFE_ALLOCATE(this%gmetempn, (pol%nmtx))
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
      SAFE_ALLOCATE(this%gmeRDyntempn, (pol%nmtx))
      SAFE_ALLOCATE(this%chiRDyntmp, (pol%nfreq_in_group))
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      SAFE_ALLOCATE(this%gmeRDyntempn, (pol%nmtx))
      SAFE_ALLOCATE(this%chiRDyntmp, (pol%nfreq_in_group))
      SAFE_ALLOCATE(this%chiRDynOld, (pol%nfreq_in_group,scal%npr,scal%npc,nspin))
    endif
    if (pol%freq_dep .eq. 3) then
      SAFE_ALLOCATE(this%gmeRDyntempn, (pol%nmtx))
      SAFE_ALLOCATE(this%chiRDyntmp, (pol%nfreq_in_group))
    endif

    POP_SUB(create_chi_summator)
    return

  end subroutine create_chi_summator

  subroutine free_chi_summator(this, pol)
    type(chi_summator_t), intent(INOUT) :: this !<the chi_summator_t object
    type(polarizability), intent(IN) :: pol

    PUSH_SUB(free_chi_summator)

    if (pol%freq_dep .eq. 0) then
      SAFE_DEALLOCATE(this%gmetempn)
    endif
    if (pol%freq_dep .eq. 2) then
      SAFE_DEALLOCATE(this%gmeRDyntempn)
      SAFE_DEALLOCATE(this%chiRDyntmp)
      if ((pol%freq_dep_method .eq. 1)) then
        SAFE_DEALLOCATE(this%chiRDynOld)
      endif
    endif
    if (pol%freq_dep .eq. 3) then
      SAFE_DEALLOCATE(this%gmeRDyntempn)
      SAFE_DEALLOCATE(this%chiRDyntmp)
    endif

    POP_SUB(free_chi_summator)
    return

  end subroutine free_chi_summator

  !-----------------------------------------------------------------------------
  !                              GCOMM_MATRIX
  !-----------------------------------------------------------------------------

  !> Create the pol%chi matrix using gcomm_matrix sceheme
  subroutine chi_summation_comm_matrix(this,pol,scal,kp,kpq,vwfn,&
    nst,nrk,indt,pht)
    type(chi_summator_t), intent(INOUT) :: this
    type(polarizability), intent(INOUT) :: pol
    type(scalapack), intent(in) :: scal
    type(kpoints), intent(IN) :: kp,kpq
    type(valence_wfns), intent(IN) :: vwfn

    integer, intent(IN) :: nst(:)
    integer, intent(IN) :: nrk
    integer, intent(INOUT) :: indt(:,:,:)
    SCALAR,  intent(INOUT) :: pht(:,:,:)

    integer :: ntot_members(pol%nfreq_group)
    integer :: icurr,ntot,ntot2,itot,ntotmax,ifreq_para
    integer :: ipe, ilimit, ii, jj, kk, ll, iv, irk, it, ispin,im,mytot
    integer :: i_myband,tag,irank,tmp_iv,im_proc
    complex(DPC) :: negfact
    real(DP) :: zvalue, occ_diff,cv_energy ! DVF partial occupations
    type(cvpair_info) :: cvpair_temp
    integer, allocatable :: tmprowindex(:),tmpcolindex(:)
    SCALAR, allocatable :: tmprowph(:),tmpcolph(:), sendbuf(:,:)
    complex(DPC), allocatable :: gmeRDyntempr(:,:)
    complex(DPC), allocatable :: edenDRtemp(:,:)

    ! frequency points for the spectral functions of the polarizability
    integer :: isfreql, isfreqr, nwarn
    real(DP) :: sfreql, sfreqr, wr, wl
    ! Hilbert transform coefficients
    complex(DPC) :: htwR(pol%nfreq,pol%nsfreq), htwA(pol%nfreq,pol%nsfreq)
    complex(DPC) :: c11,c12,c13,c14,c21,c22,c23,c24
    complex(DPC) :: cA11,cA12,cA13,cA14,cA21,cA22,cA23,cA24

    integer :: isf,nf,nsf, request
    real(DP) :: sf1,sf2
    real(DP) :: step1,step2,fqt,eta
    complex(DPC) :: c1,c2,j_dpc,cA1,cA2
    logical :: first_reduce

    ! variables for non-blocking cyclic scheme
    integer :: isend_static, irec_static
    integer :: actual_send, actual_rec
    integer :: nsend_row, nsend_col, nrec_row, nrec_col
    integer :: req_send, tag_send, req_rec, tag_rec
    integer :: nrow_max, ncol_max
#ifdef MPI
    integer :: stat(MPI_STATUS_SIZE)
#endif
    SCALAR, allocatable :: buf_send(:,:,:), buf_rec(:,:,:), buf_temp(:,:,:)

    type(progress_info) :: prog_info

    PUSH_SUB(chi_summation_comm_matrix)

    !call alloc_summation_buffers(pol, this%fact)

    if (pol%freq_dep.eq.0) then
      SAFE_ALLOCATE(this%chilocal2, (scal%npr,scal%npc,kp%nspin))
!$OMP PARALLEL DO collapse(3)
      do ii = 1, kp%nspin
        do jj = 1, scal%npc
          do kk = 1, scal%npr
            this%chilocal2(kk,jj,ii)=0D0
          enddo
        enddo
      enddo
    endif

    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
      SAFE_ALLOCATE(this%chilocal2RDyn, (scal%npr,scal%npc,pol%nfreq_in_group,kp%nspin))
!$OMP PARALLEL DO collapse(4)
      do ii = 1, kp%nspin
        do jj = 1, pol%nfreq_in_group
          do kk = 1, scal%npc
            do ll = 1, scal%npr
              this%chilocal2RDyn(ll,kk,jj,ii)=0D0
            enddo
          enddo
        enddo
      enddo
    endif

! At this moment the spectral method only works for gcomm_elements.
! call die in inread.f90
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      SAFE_ALLOCATE(this%chilocal2RDyn, (scal%npr,scal%npc,pol%nfreq_in_group,kp%nspin))
!$OMP PARALLEL DO collapse(4)
      do ii = 1, kp%nspin
        do jj = 1, pol%nfreq_in_group
          do kk = 1, scal%npc
            do ll = 1, scal%npr
              this%chilocal2RDyn(ll,kk,jj,ii)=0D0
            enddo
          enddo
        enddo
      enddo

    endif

    if (pol%freq_dep .eq. 3) then
      SAFE_ALLOCATE(this%chilocal2RDyn, (scal%npr,scal%npc,pol%nfreq_in_group,kp%nspin))
!$OMP PARALLEL DO collapse(4)
      do ii = 1, kp%nspin
        do jj = 1, pol%nfreq_in_group
          do kk = 1, scal%npc
            do ll = 1, scal%npr
              this%chilocal2RDyn(ll,kk,jj,ii)=0D0
            enddo
          enddo
        enddo
      enddo
    endif

    ntot=0
    ntot2=0

    ntot = peinf%nvownactual*peinf%ncownactual
    do irk = 1, nrk
      ntot2=ntot2 + nst(irk)
    enddo
    ntot=ntot*ntot2

    !-------------------------------------------------------------------
    ! Static Be Here


    if (pol%freq_dep .eq. 0) then

#ifdef MPI
      call timing%start(timing%chi_sum_bar)
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call timing%stop(timing%chi_sum_bar)
#endif

      !XXXXXXXXXXXXXXXX
      ! Temporary fix!!!!
      if (chi_summation_algo == OMP_TARGET_ALGO) then
        pol%accel_chisum_fulloffload = .false.
      end if
      !XXXXXXXXXXXXXXXX

      ! setup the blas handle if existing
      call create_common_blas_handle( chi_summation_algo )
      call set_queue_to_common_blas_handle( chi_summation_algo, 1 )

      first_reduce = .true.

      if ( chi_summation_algo /= CPU_ALGO .and. peinf%inode == 0 ) then
        if( pol%accel_chisum_fulloffload ) then
          write(6,'(1x,A)') 'CHI summation GPU algorithm offloading all matrix elements to device.'
        else
          write(6,'(1x,A)') 'CHI summation GPU algorithm offloading only submatrices to device.'
        end if
      end if

      call progress_init(prog_info, 'building polarizability matrix', 'processor', &
        peinf%npes)

      if(pol%nonblocking_cyclic) then
        ! initialize buffer for non-blocking cyclic scheme
        ! static process for the communication
        isend_static = MOD(peinf%inode + 1 + peinf%npes, peinf%npes)
        irec_static = MOD(peinf%inode - 1 + peinf%npes, peinf%npes)
        ! allocate my size for the first send
        nsend_row = scal%nprd(peinf%inode+1)
        nsend_col = scal%npcd(peinf%inode+1)
        SAFE_ALLOCATE(buf_send, (nsend_row, nsend_col, kp%nspin))
        do ispin = 1 , kp%nspin
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, nsend_col
            do jj = 1, nsend_row
              buf_send(jj,ii,ispin) = ZERO
            enddo
          enddo
        end do
      end if

      nrow_max = MAXVAL(scal%nprd)
      ncol_max = MAXVAL(scal%npcd)
      SAFE_ALLOCATE(this%chilocal, (nrow_max, ncol_max))
      SAFE_ALLOCATE(this%gmetempr, (nrow_max, ntot))
      SAFE_ALLOCATE(this%gmetempc, (ncol_max, ntot))

      if ( chi_summation_algo /= CPU_ALGO ) then
        select case(chi_summation_algo)
         case( OPENACC_ALGO )
#ifdef OPENACC
           if ( pol%accel_chisum_fulloffload ) then
             !$acc enter data copyin(pol) async(1)
             !$acc enter data copyin(pol%gme) async(1)
             !$acc enter data copyin(scal) async(1)
             !$acc enter data copyin(scal%nprd, scal%npcd, scal%imyrowd, scal%imycold) async(1)
             !$acc enter data copyin(indt, pht) async(1)
             !$acc enter data copyin(peinf) async(1)
           end if
           !$acc enter data create(this) async(1)
           !$acc enter data create(this%chilocal, this%gmetempr, this%gmetempc) async(1)
#else
           call die_algos("OpenACC")
#endif
         case( OMP_TARGET_ALGO )
#ifdef OMP_TARGET
           !$omp target enter data map (to: this)
           !$omp target enter data map (to: this%chilocal, this%gmetempr, this%gmetempc)
#else
           call die_algos("OpenMP Target")
#endif
         case default
           call die("Invald algorithm for chi_summation.", only_root_writes = .true.)
        end select
      end if

      do ipe = 1, peinf%npes
        call progress_step(prog_info)

        if(pol%nonblocking_cyclic) then
          ! calculate the actual process we have to send and we are receiving
          actual_send = MOD(peinf%inode + ipe + peinf%npes, peinf%npes)
          actual_rec  = MOD(peinf%inode - ipe + peinf%npes, peinf%npes)
          nrec_row = scal%nprd(actual_rec+1)
          nrec_col = scal%npcd(actual_rec+1)
          ! allocate reciving buffer
          SAFE_ALLOCATE(buf_rec, (nrec_row,nrec_col,kp%nspin))
          do ispin = 1 , kp%nspin
            !$OMP PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_rec(jj,ii,ispin) = ZERO
              enddo
            enddo
          end do
#ifdef MPI
          ! post message
          call timing%start(timing%chi_sum_comm)
          tag_rec = 1
          tag_send = 1
          req_rec = MPI_REQUEST_NULL
          req_send = MPI_REQUEST_NULL
          CALL MPI_Irecv(buf_rec, nrec_row*nrec_col*kp%nspin,MPI_SCALAR,irec_static,&
                         tag_rec, MPI_COMM_WORLD, req_rec, mpierr)
          CALL MPI_Isend(buf_send, nsend_row*nsend_col*kp%nspin,MPI_SCALAR,isend_static,&
                         tag_send, MPI_COMM_WORLD, req_send, mpierr)
#endif
          call timing%stop(timing%chi_sum_comm)

          ! allocate working array
          call timing%start(timing%chi_sum_array_alloc)
!$OMP PARALLEL DO collapse(2)
          do ii = 1, nrec_col
            do jj = 1, nrec_row
              this%chilocal(jj,ii)=0D0
            enddo
          enddo
          SAFE_ALLOCATE(buf_temp, (nrec_row, nrec_col, kp%nspin))
          do ispin = 1 , kp%nspin
!$OMP PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_temp(jj,ii,ispin) = ZERO
              enddo
            enddo
          end do
          call timing%stop(timing%chi_sum_array_alloc)
        else

          call timing%start(timing%chi_sum_array_alloc)
!$OMP PARALLEL DO collapse(2)
          do ii = 1, scal%npcd(ipe)
            do jj = 1, scal%nprd(ipe)
              this%chilocal(jj,ii)=0D0
            enddo
          enddo
! JRD XXX Should change order here like do in FF case. I think I did...
          call timing%stop(timing%chi_sum_array_alloc)
        end if  ! pol%nonblocking_cyclic

        do ispin = 1 , kp%nspin
          if(pol%nonblocking_cyclic) then
            call mtxelmultiply(scal,ntot,nrk,nst,this%fact,vwfn, &
              this%gmetempr,this%gmetempc,this%chilocal,pol%gme,pol,indt,pht,actual_rec+1,ispin,&
               nrow_max, ncol_max)
            if ( chi_summation_algo /= CPU_ALGO ) then
#ifdef MPI
              call timing%start(timing%chi_sum_comm)
              ! make sure the buffer is received
              CALL MPI_Wait(req_rec,stat,mpierr)
              CALL MPI_Wait(req_send,stat,mpierr)
              call timing%stop(timing%chi_sum_comm)
#endif

              if ( chi_summation_algo == OPENACC_ALGO ) then
#ifdef OPENACC
                !$acc wait(1)
#else
                call die_algos("OpenACC")
#endif
              end if
            end if
            !$OMP PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_temp(jj,ii,ispin) = this%chilocal(jj,ii)
              enddo
            enddo
          else

            call mtxelmultiply(scal,ntot,nrk,nst,this%fact,vwfn, &
              this%gmetempr,this%gmetempc,this%chilocal,pol%gme,pol,indt,pht,ipe,ispin,&
               nrow_max, ncol_max)

            if ( chi_summation_algo == OPENACC_ALGO ) then
#ifdef OPENACC
              !$acc wait(1)
#else
              call die_algos("OpenACC")
#endif
            end if
#ifdef NONBLOCKING
          call timing%start(timing%chi_sum_comm)
          if (.not. first_reduce) then
            call MPI_WAIT(request,mpistatus,mpierr)
            SAFE_DEALLOCATE(sendbuf)
          endif
          call timing%stop(timing%chi_sum_comm)

          call timing%start(timing%chi_sum_flt)
          SAFE_ALLOCATE(sendbuf, (scal%nprd(ipe),scal%npcd(ipe)))
!$OMP PARALLEL DO collapse(2)
          do ii = 1, scal%npcd(ipe)
            do jj = 1, scal%nprd(ipe)
              sendbuf(jj,ii) = this%chilocal(jj,ii)
            enddo
          enddo
          call timing%stop(timing%chi_sum_flt)

          call timing%start(timing%chi_sum_ht_nb)
#ifdef MPI
! JRD XXX Barrier probably implicit in IReduce
          call MPI_IReduce(sendbuf(1,1),this%chilocal2(1,1,ispin),scal%nprd(ipe)*scal%npcd(ipe),MPI_SCALAR, &
            MPI_SUM,ipe-1,MPI_COMM_WORLD,request,mpierr)
#else
          this%chilocal2(:,:,ispin)=this%chilocal(1:scal%nprd(ipe),1:scal%npcd(ipe))
#endif
          first_reduce = .false.
          call timing%stop(timing%chi_sum_ht_nb)

#else
          call timing%start(timing%chi_sum_comm)
#ifdef MPI
          !call MPI_barrier(MPI_COMM_WORLD,mpierr)
          if ( peinf%inode == (ipe-1) ) then
            call MPI_Reduce(MPI_IN_PLACE,this%chilocal(1,1),ncol_max*nrow_max,MPI_SCALAR, &
              MPI_SUM,ipe-1,MPI_COMM_WORLD,mpierr)
            ! update local buffer
            this%chilocal2(:,:,ispin)=this%chilocal(1:scal%nprd(ipe),1:scal%npcd(ipe))
          else
            call MPI_Reduce(this%chilocal(1,1),this%chilocal(1,1),ncol_max*nrow_max,MPI_SCALAR, &
              MPI_SUM,ipe-1,MPI_COMM_WORLD,mpierr)
          end if
#else
          this%chilocal2(:,:,ispin)=this%chilocal(1:scal%nprd(ipe),1:scal%npcd(ipe))
#endif
          call timing%stop(timing%chi_sum_comm)
#endif
          ! in case nonblocking_cyclic communication will be finalize outside the spin-loop
          end if ! pol%nonblocking_cyclic

        enddo ! ispin

        if(pol%nonblocking_cyclic) then
#ifdef MPI
          call timing%start(timing%chi_sum_comm)
          ! make sure the buffer is received
          CALL MPI_Wait(req_rec,stat,mpierr)
#endif
          ! accumulate contribution into receiving buffer
          ! buf_rec(:,:,:) = buf_rec(:,:,:) + buf_temp
          do ispin = 1 , kp%nspin
!$OMP PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_rec(jj,ii,ispin) = buf_rec(jj,ii,ispin) + buf_temp(jj,ii,ispin)
              enddo
            enddo
          end do
          SAFE_DEALLOCATE(buf_temp)
#ifdef MPI
          ! wait for the massage to be sent
          CALL MPI_Wait(req_send,stat,mpierr)
#endif
          ! copy the messega to the sending buffer for the next cycle
          SAFE_DEALLOCATE(buf_send)
          SAFE_ALLOCATE(buf_send, (nrec_row, nrec_col, kp%nspin))
          buf_send = buf_rec
          nsend_row = nrec_row
          nsend_col = nrec_col
          ! deallocate receiving buffer
          SAFE_DEALLOCATE(buf_rec)
#ifdef MPI
          call timing%stop(timing%chi_sum_comm)
#endif
        end if

      enddo ! ipe

      if ( chi_summation_algo /= CPU_ALGO ) then
        select case(chi_summation_algo)
         case( OPENACC_ALGO )
#ifdef OPENACC
           if ( pol%accel_chisum_fulloffload ) then
             !$acc exit data delete(pol%gme) async(1)
             !$acc exit data delete(pol) async(1)
             !$acc exit data delete(scal%nprd, scal%npcd, scal%imyrowd, scal%imycold) async(1)
             !$acc exit data delete(scal) async(1)
             !$acc exit data delete(indt, pht) async(1)
             !$acc exit data delete(peinf) async(1)
           end if
           !$acc exit data delete(this%chilocal, this%gmetempr, this%gmetempc) async(1)
           !$acc exit data delete(this) async(1)
           !$acc wait(1)
#else
           call die_algos("OpenACC")
#endif
         case( OMP_TARGET_ALGO )
#ifdef OMP_TARGET
           !$omp target exit data map(delete: this%chilocal, this%gmetempr, this%gmetempc)
           !$omp target exit data map(delete: this)
#else
           call die_algos("OpenMP Target")
#endif
         case default
           call die("Invald algorithm for chi_summation.", only_root_writes = .true.)
        end select
      end if

      SAFE_DEALLOCATE(this%chilocal)
      SAFE_DEALLOCATE(this%gmetempr)
      SAFE_DEALLOCATE(this%gmetempc)

      if(pol%nonblocking_cyclic) then
        ! done
        this%chilocal2(:,:,:) = buf_send(:,:,:)
        SAFE_DEALLOCATE(buf_send)
      else
#ifdef NONBLOCKING
        call timing%start(timing%chi_sum_comm)

        call MPI_WAIT(request,mpistatus,mpierr)
        SAFE_DEALLOCATE(sendbuf)
        first_reduce = .true.

        call timing%stop(timing%chi_sum_comm)
#endif
      end if

      call progress_free(prog_info)

      do ispin =1, kp%nspin
        pol%chi(:,:,ispin) = this%chilocal2(:,:,ispin)
      enddo
      SAFE_DEALLOCATE(this%chilocal2)

      call destroy_common_blas_handle( chi_summation_algo )

      !XXXXXXXXXXx
      ! call diagonalize_scalapack(pol%nmtx, scal, pol%chi(:,:,1), 1.0D-3, icurr, this%chilocal)
      !XXXXXXXXXXX

    endif ! pol%freq_dep .eq. 0

    !-------------------------------------
    ! Full Frequency Be Here

    negfact = -1D0*this%fact

    if ( ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) .or. pol%freq_dep .eq. 3 ) then
      if(pol%nfreq_group .gt. 1 .and. peinf%inode .lt. peinf%npes) then
#ifdef MPI
        tag = 10 ! MPI tag
        if(peinf%rank_mtxel .eq. 0) then
          ntot_members(1) = ntot
        endif
        do irank = 1, pol%nfreq_group-1
          if(peinf%rank_mtxel .eq. irank) then
            call MPI_Send(ntot,1,MPI_INTEGER,0,tag,peinf%mtxel_comm,mpierr)
          endif
          if(peinf%rank_mtxel .eq. 0) then
            call MPI_Recv(ntot_members(irank+1),1,MPI_INTEGER,irank,tag,peinf%mtxel_comm,mpistatus,mpierr)
          endif
        enddo
        if(peinf%rank_mtxel .eq. 0) then
          ntot = sum(ntot_members(1:pol%nfreq_group))
        endif
        call MPI_Bcast(ntot_members,pol%nfreq_group,MPI_INTEGER,0,peinf%mtxel_comm,mpierr)
        call MPI_Bcast(ntot,1,MPI_INTEGER,0,peinf%mtxel_comm,mpierr)
        call MPI_Bcast(ntotmax,1,MPI_INTEGER,0,peinf%mtxel_comm,mpierr)
        ntotmax = peinf%nvownmax*peinf%ncownmax*nrk
        !Communicate gme`s and energy denominators among processors in matrix element communication groups
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_comm)
        endif
        call MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,pol%gme(1,1,1,1,1,1),&
                          ntotmax*kp%nspin*pol%nmtx,MPI_SCALAR,peinf%mtxel_comm,mpierr)
        call MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,pol%edenDyn(1,1,1,1,1),&
                          ntotmax*kp%nspin,MPI_REAL_DP,peinf%mtxel_comm,mpierr)
        call MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,pol%occ_diff(1,1,1,1,1),&
                          ntotmax*kp%nspin,MPI_REAL_DP,peinf%mtxel_comm,mpierr)
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_comm)
        endif
#endif
      endif

      call progress_init(prog_info, 'building polarizability matrix', 'processor', &
        peinf%npes_freqgrp)
      do ipe = 1, peinf%npes_freqgrp
        call progress_step(prog_info)
        if (peinf%verb_debug .and. peinf%inode==0) then
          write(6,'(A,i8,6x,A,i8,A)') '### ipe=',ipe,'(npes=',peinf%npes,')'
        endif
#ifdef MPI
        call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif
        SAFE_ALLOCATE(this%gmeRDyntempr3, (scal%nprd(ipe),ntot))
        SAFE_ALLOCATE(this%gmeRDyntempc, (scal%npcd(ipe),ntot))
        SAFE_ALLOCATE(this%chilocalRDyn, (scal%nprd(ipe),scal%npcd(ipe),pol%nfreq_in_group))

        SAFE_ALLOCATE(gmeRDyntempr, (scal%nprd(ipe),ntot))
        SAFE_ALLOCATE(edenDRtemp, (ntot, pol%nfreq_in_group))


!$OMP PARALLEL DO collapse(2)
        do ii = 1, pol%nfreq_in_group
          do jj = 1, scal%npcd(ipe)
            do kk = 1, scal%nprd(ipe)
              this%chilocalRDyn(kk,jj,ii)=0D0
            enddo
          enddo
        enddo

        do ispin = 1 , kp%nspin

          itot = 0

          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_prep)
          endif

          SAFE_ALLOCATE(tmprowindex,(scal%nprd(ipe)))
          SAFE_ALLOCATE(tmpcolindex,(scal%npcd(ipe)))
          SAFE_ALLOCATE(tmprowph,(scal%nprd(ipe)))
          SAFE_ALLOCATE(tmpcolph,(scal%npcd(ipe)))

          do im=1,pol%nfreq_group                          ! im labels which member of the mtxel comm group are you
            im_proc=peinf%rank_f+1+(im-1)*peinf%npes_freqgrp ! im_proc gives this mtxel comm group member`s global
            do irk = 1, nrk                                  ! proc number
              do it = 1, nst(irk)

                do icurr=1,scal%nprd(ipe)
                  tmprowindex(icurr) = indt(scal%imyrowd(icurr,ipe),it,irk)
                  tmprowph(icurr) = pht(scal%imyrowd(icurr,ipe),it,irk)
                enddo
                do icurr=1,scal%npcd(ipe)
                  tmpcolindex(icurr) = indt(scal%imycold(icurr,ipe),it,irk)
                  tmpcolph(icurr) = pht(scal%imycold(icurr,ipe),it,irk)
                enddo

! JRD XXX - Would rather put OMP here with collapse 2 - but unbalanced... need to loop over tmp_iv = 1, nvown instead
                do iv = 1,vwfn%nband+pol%ncrit

                  tmp_iv = peinf%global_indexv(iv,im_proc)

                  if (peinf%does_it_ownv(iv,im_proc)) then
                    ilimit = peinf%global_ncown(im_proc)
                  else
                    ilimit = 0
                  endif

                  !$OMP PARALLEL DO private (mytot,zvalue,jj,icurr,ifreq_para,cvpair_temp,cv_energy)
                  do i_myband = 1, ilimit
                    mytot = itot + i_myband
                    zvalue = pol%edenDyn(peinf%global_indexv(iv,im_proc),i_myband,ispin,irk,im)
! DVF partial occupations; might need to change the OMP directive above for
! threading to work
                    occ_diff = pol%occ_diff(peinf%global_indexv(iv,im_proc),i_myband,ispin,irk,im)
                    if(pol%lin_denominator<TOL_Zero) then
                      ! this is when the lin_denominator mode is not active.

                      ! DVF: These factors of (peinf%inode)/peinf%npes*100000 are so that you don`t access the
                      ! array edenDRtemp and other similar arrays for the processors that are not used
                      ! when pol%nfreq_group does not divide the original number of processors in the
                      ! calculation (peinf%npes_orig). This factor makes the loop start at some huge
                      ! number, so the loop will in fact not execute at all

                      do jj=1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
                        ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
                        if (abs(zvalue) .gt. Tol_Zero) then
                          ! DYQ: To get a polarizability consistent with the Tamm-Dancoff Approximation,
                          ! we only include positive frequency terms in the full frequency calculation
                          if(pol%tda) then
                            edenDRtemp(mytot,ifreq_para)= -0.5d0*( &
                              1d0/(zvalue+(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd))
                          else
                            ! JRD XXX - INDEX ORDER here worries me
                            edenDRtemp(mytot,ifreq_para)= -0.5d0*occ_diff*( &
                              1d0/(zvalue-(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd)+ &
                              1d0/(zvalue+(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd))
                          endif
                        else
                          edenDRtemp(mytot,ifreq_para)= 0D0
                        endif
                      enddo

!#BEGIN_INTERNAL_ONLY
                    else
                      !this is when lin_denominator mode is active

                      cvpair_temp = pol%lin_edenDyn(peinf%global_indexv(iv,im_proc),i_myband,ispin,irk,im)
                      cv_energy = cvpair_temp%ev - cvpair_temp%ec
                      do jj=1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
                      !check if denominator is not small. If so, do the usual
                      !thing
                        ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
                        if(abs(cv_energy+pol%dFreqGrid(jj)/ryd)-pol%lin_denominator>-TOL_Zero .and. &
                           abs(cv_energy-pol%dFreqGrid(jj)/ryd)-pol%lin_denominator>-TOL_Zero) then

                          if (abs(zvalue) .gt. Tol_Zero) then
                            if(pol%tda) then
                              edenDRtemp(mytot,ifreq_para)= -0.5d0*( &
                                1d0/(zvalue+(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd))                              
                            else
                              edenDRtemp(mytot,ifreq_para)= -0.5d0*( &
                                1d0/(zvalue-(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd)+ &
                                1d0/(zvalue+(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd))
                            endif
                          else
                            edenDRtemp(mytot,ifreq_para)= 0D0
                          endif

                        else ! if denominator is small, do linearized energy denominator
                          if(cvpair_temp%vltc) then
                            if(pol%tda) then
                              edenDRtemp(mytot,ifreq_para) = -0.5d0*(&
                                integrate_mbz(kp,kpq%rk(1:3,cvpair_temp%idx_kp),&
                                cvpair_temp%ev,cvpair_temp%ec,&
                                cvpair_temp%vv,cvpair_temp%vc,&
                                pol%dFreqGrid(jj)/ryd,pol%efermi/ryd))
                            else
                              edenDRtemp(mytot,ifreq_para) = -0.5d0*(&
                                conjg(integrate_mbz(kp,kpq%rk(1:3,cvpair_temp%idx_kp),&
                                cvpair_temp%ev,cvpair_temp%ec,&
                                cvpair_temp%vv,cvpair_temp%vc,&
                                -pol%dFreqGrid(jj)/ryd,pol%efermi/ryd))&
                                
                                +integrate_mbz(kp,kpq%rk(1:3,cvpair_temp%idx_kp),&
                                cvpair_temp%ev,cvpair_temp%ec,&
                                cvpair_temp%vv,cvpair_temp%vc,&
                                pol%dFreqGrid(jj)/ryd,pol%efermi/ryd))
                            endif
                            else
                            edenDRtemp(mytot,ifreq_para)= 0D0
                          endif
                        endif

                      enddo
!#END_INTERNAL_ONLY
                    endif
                    do icurr=1,scal%nprd(ipe)
                      gmeRDyntempr(icurr,mytot)=pol%gme(tmprowindex(icurr), &
                        i_myband,tmp_iv,ispin,irk,im) * tmprowph(icurr)
                    enddo
                    !do jj = 1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
                    !  ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
                    !  this%gmeRDyntempr2(:,mytot,ifreq_para)=gmeRDyntempr(:)*edenDRtemp(ifreq_para)
                    !enddo
                    do icurr=1,scal%npcd(ipe)
                      this%gmeRDyntempc(icurr,mytot) = &
                        MYCONJG( pol%gme(tmpcolindex(icurr),i_myband,tmp_iv,ispin,irk,im) * tmpcolph(icurr) )
                    enddo
                  enddo ! i_myband
                  !$OMP END PARALLEL DO
                  itot = itot+ilimit
                enddo ! iv
              enddo ! it
            enddo ! irk
          enddo ! im

          SAFE_DEALLOCATE(tmprowindex)
          SAFE_DEALLOCATE(tmpcolindex)
          SAFE_DEALLOCATE(tmprowph)
          SAFE_DEALLOCATE(tmpcolph)

          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_prep)
          endif


          !Do the zgemm`s
          if(ntot > 0) then
            do jj =1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
              ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_gemm)
              endif

              !$omp parallel do private(icurr)
              do mytot = 1, ntot
                do icurr=1,scal%nprd(ipe)
                  this%gmeRDyntempr3(icurr,mytot)=gmeRDyntempr(icurr,mytot)*edenDRtemp(mytot,ifreq_para)
                enddo
              enddo

              call zgemm('n','t',scal%nprd(ipe),scal%npcd(ipe),ntot, &
                negfact,this%gmeRDyntempr3(:,:),scal%nprd(ipe),this%gmeRDyntempc(:,:),scal%npcd(ipe),&
                (0D0,0D0),this%chilocalRDyn(:,:,ifreq_para),scal%nprd(ipe))

              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_gemm)
              endif
            enddo
          endif

          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_comm)
          endif

          if(pol%nfreq_group .eq. 1) then
#ifdef MPI
            call MPI_Reduce(this%chilocalRDyn(1,1,1),this%chilocal2RDyn(1,1,1,ispin), &
              pol%nfreq_in_group*scal%npcd(ipe)*scal%nprd(ipe),MPI_COMPLEX_DPC,&
              MPI_SUM,ipe-1,MPI_COMM_WORLD,mpierr)
#endif
          elseif(pol%nfreq_group .gt. 1 .and. peinf%inode .lt. peinf%npes) then
#ifdef MPI
            call MPI_Reduce(this%chilocalRDyn(1,1,1),this%chilocal2RDyn(1,1,1,ispin), &
              pol%nfreq_in_group*scal%npcd(ipe)*scal%nprd(ipe),MPI_COMPLEX_DPC,&
              MPI_SUM,ipe-1,peinf%freq_comm,mpierr)
#endif
          endif
#ifndef MPI
          this%chilocal2RDyn(:,:,:,ispin)=this%chilocalRDyn(:,:,:)
#endif
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_comm)
          endif

        enddo ! ispin
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_array_alloc)
        endif
        SAFE_DEALLOCATE(this%chilocalRDyn)
        SAFE_DEALLOCATE(edenDRtemp)
        SAFE_DEALLOCATE(gmeRDyntempr)
        SAFE_DEALLOCATE(this%gmeRDyntempr3)
        SAFE_DEALLOCATE(this%gmeRDyntempc)
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_array_alloc)
        endif
      enddo ! ipe
      call progress_free(prog_info)

      do ispin =1, kp%nspin
        do jj=1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
          ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
! JRD XXX This assignment is now a waste of time and memory. Should just set pol%chiR(A)Dyn
! directly above
          pol%chiRDyn(:,:,ifreq_para,ispin) = this%chilocal2RDyn(:,:,ifreq_para,ispin)
        enddo ! jj
      enddo ! ispin
      SAFE_DEALLOCATE(this%chilocal2RDyn)
      !XXXXXXXXXXX
      ! DO jj = 1, pol%nfreq
      !   IF(peinf%inode .eq. 0) WRITE(2000,*) jj
      !   call diagonalize_scalapack(pol%nmtx, scal, pol%chiRDyn(:,:,jj,1), 1.0D-3, icurr, this%chilocal, eigenval)
      !   DEALLOCATE(eigenval)
      !   DEALLOCATE(this%chilocal)
      ! END DO
      !XXXXXXXXXXX
    endif ! (pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0)

    !-------------------------------------
    ! Full Frequency Be Here.
    ! M. Shishkin and G. Kresse, Implementation and performance of the
    ! frequency-dependent GW method within the PAW framework,
    ! PHYSICAL REVIEW B 74, 035101, 2006.

    negfact = -1D0*this%fact

    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      if(pol%nfreq_group .gt. 1 .and. peinf%inode .lt. peinf%npes) then
#ifdef MPI
        tag = 10 ! MPI tag
        if(peinf%rank_mtxel .eq. 0) then
          ntot_members(1) = ntot
        endif
        do irank = 1, pol%nfreq_group-1
          if(peinf%rank_mtxel .eq. irank) then
            call MPI_Send(ntot,1,MPI_INTEGER,0,tag,peinf%mtxel_comm,mpierr)
          endif
          if(peinf%rank_mtxel .eq. 0) then
            call MPI_Recv(ntot_members(irank+1),1,MPI_INTEGER,irank,tag,peinf%mtxel_comm,mpistatus,mpierr)
          endif
        enddo
        if(peinf%rank_mtxel .eq. 0) then
          ntot = sum(ntot_members(1:pol%nfreq_group))
        endif
        call MPI_Bcast(ntot_members,pol%nfreq_group,MPI_INTEGER,0,peinf%mtxel_comm,mpierr)
        call MPI_Bcast(ntot,1,MPI_INTEGER,0,peinf%mtxel_comm,mpierr)
        call MPI_Bcast(ntotmax,1,MPI_INTEGER,0,peinf%mtxel_comm,mpierr)
        ntotmax = peinf%nvownmax*peinf%ncownmax*nrk
        !Communicate gme`s among processors in matrix element communication groups
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_comm)
        endif
                  call MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,pol%gme(1,1,1,1,1,1),&
                                   ntotmax*kp%nspin*pol%nmtx,MPI_SCALAR,peinf%mtxel_comm,mpierr)
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_comm)
        endif
#endif
      endif

      ! -------------------------------------------
      ! compute the Hilbert transform coefficients
      ! -------------------------------------------
      j_dpc=(0.0,1.0)
      htwR(:,:)=(0.0,0.0)

      nf=pol%nfreq
      nsf=pol%nsfreq
      do jj=1,nf
        eta=pol%dBrdning/ryd
        fqt=pol%dFreqGrid(jj)/ryd
        do isf=1,nsf
          if (isf==1) then
            c1=(0.0,0.0)
            step1=1.0
          else
            sf1=pol%dSFreqGrid(isf-1)/ryd
            sf2=pol%dSFreqGrid(isf)/ryd
            step1=sf2-sf1

            c11=((sf1-fqt)*(-1.0*j_dpc)-eta)*(atan((fqt-sf2)/(eta))&
              &          -atan((fqt-sf1)/(eta)))
            c12=0.50*(sf1-fqt+(-1.0*j_dpc)*eta)*log(((fqt-sf2)**2.0+eta*eta)&
              &          /(((fqt-sf1)**2.0+eta*eta)))

            c13=-((sf1+fqt)*j_dpc-eta)*(atan((fqt+sf2)/(eta))&
              &          -atan((fqt+sf1)/(eta)))
            c14=0.50*(sf1+fqt+j_dpc*eta)*log(((fqt+sf2)**2.0+eta*eta)&
              &          /(((fqt+sf1)**2.0+eta*eta)))

            c1=c11+c12+c13+c14
            c1=c1/step1
          endif

          if (isf==nsf) then
            c2=(0.0,0.0)
            step2=1.0
          else
            sf1=pol%dSFreqGrid(isf)/ryd
            sf2=pol%dSFreqGrid(isf+1)/ryd

            step2=sf2-sf1

            c21=((sf2-fqt)*(-1.0*j_dpc)-eta)*(atan((fqt-sf1)/(eta))&
              &          -atan((fqt-sf2)/(eta)))
            c22=0.50*(sf2-fqt+(-1.0*j_dpc)*eta)*log(((fqt-sf1)**2.0+eta*eta)&
              &          /(((fqt-sf2)**2.0+eta*eta)))

            c23=-((sf2+fqt)*j_dpc-eta)*(atan((fqt+sf1)/(eta))&
              &          -atan((fqt+sf2)/(eta)))
            c24=0.50*(sf2+fqt+j_dpc*eta)*log(((fqt+sf1)**2.0+eta*eta)&
              &          /(((fqt+sf2)**2.0+eta*eta)))

            c2=c21+c22+c23+c24

            c2=c2/step2
          endif

          if (isf==1.or.isf==nsf) then
            htwR(jj,isf)=0.5d0*(c1/step1+c2/step2)
          else
            htwR(jj,isf)=1.0d0*(c1+c2)/(step1+step2)
          endif

        enddo
      enddo

      ! ----------------------------------------------------
      ! compute the spectral functions of the polarizability
      ! ----------------------------------------------------

      call progress_init(prog_info, 'building polarizability matrix', 'processor', &
        peinf%npes_freqgrp)

      nwarn=0
      do ipe = 1, peinf%npes_freqgrp
        call progress_step(prog_info)
        if (peinf%verb_debug .and. peinf%inode==0) then
          write(6,'(A,i8,6x,A,i8,A)') '### ipe=',ipe,'(npes=',peinf%npes,')'
        endif
#ifdef MPI
        call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif

        SAFE_ALLOCATE(gmeRDyntempr, (scal%nprd(ipe),1))
        SAFE_ALLOCATE(this%gmeRDyntempr2, (scal%nprd(ipe),ntot,pol%os_nsfreq_para))
        SAFE_ALLOCATE(this%gmeRDyntempc, (ntot,scal%npcd(ipe)))
        SAFE_ALLOCATE(this%chilocalRDyn, (scal%nprd(ipe),scal%npcd(ipe),pol%os_nsfreq_para))

! JRD XXX should thread if keep COMM elements
        this%chilocalRDyn=0

! JRD XXX if we keep COMM elements we should thread this
        gmeRDyntempr=(0.0,0.0)
        this%gmeRDyntempr2=(0.0,0.0)
        this%gmeRDyntempc=(0.0,0.0)

        do ispin = 1 , kp%nspin

          itot = 0

          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_prep)
          endif

          SAFE_ALLOCATE(tmprowindex,(scal%nprd(ipe)))
          SAFE_ALLOCATE(tmpcolindex,(scal%npcd(ipe)))
          SAFE_ALLOCATE(tmprowph,(scal%nprd(ipe)))
          SAFE_ALLOCATE(tmpcolph,(scal%npcd(ipe)))

          do im=1,pol%nfreq_group                          ! im labels which member of the mtxel comm group are you
            im_proc=peinf%rank_f+1+(im-1)*peinf%npes_freqgrp ! im_proc gives this mtxel comm group member`s global
            do irk = 1, nrk                                  ! proc number
              do it = 1, nst(irk)

                do icurr=1,scal%nprd(ipe)
                  tmprowindex(icurr) = indt(scal%imyrowd(icurr,ipe),it,irk)
                  tmprowph(icurr) = pht(scal%imyrowd(icurr,ipe),it,irk)
                enddo
                do icurr=1,scal%npcd(ipe)
                  tmpcolindex(icurr) = indt(scal%imycold(icurr,ipe),it,irk)
                  tmpcolph(icurr) = pht(scal%imycold(icurr,ipe),it,irk)
                enddo

                do iv = 1,vwfn%nband+pol%ncrit

                  tmp_iv = peinf%global_indexv(iv,im_proc)

                  if (peinf%does_it_ownv(iv,im_proc)) then
                    ilimit = peinf%global_ncown(im_proc)
                  else
                    ilimit = 0
                  endif

                  !$OMP PARALLEL private (mytot,zvalue,gmeRDyntempr,icurr, &
                  !$OMP                   isfreql,isfreqr,sfreql,sfreqr, &
                  !$OMP                   wl,wr,i_myband,jj)
                  !$OMP DO
                  do i_myband = 1, ilimit

                    zvalue = -pol%edenDyn(peinf%global_indexv(iv,im_proc),i_myband,ispin,irk,im)
                    if (abs(zvalue) .gt. Tol_Zero) then

                      mytot = itot + i_myband
                      isfreql=-1

                      do jj=pol%nsfreq,1,-1
                        if ((pol%dSFreqGrid(jj)/ryd)<zvalue) then
                          isfreql=jj
                          EXIT
                        endif
                      enddo

                      if (isfreql.eq.pol%nsfreq) then
                        nwarn=nwarn+1
                        if (nwarn==1.and.peinf%inode.eq.0) then
                          write(0,*) 'WARNING: for accuracy, sfrequency_high_cutoff should be '
                          write(0,*) 'larger than energy difference between highest unoccupied '
                          write(0,*) 'state and lowest occupied state.'
                        endif
                        cycle
                      endif

                      sfreql=pol%dSFreqGrid(isfreql)/ryd
                      isfreqr=isfreql+1
                      sfreqr=pol%dSFreqGrid(isfreqr)/ryd

                      wr=  (zvalue-sfreql)/(sfreqr-sfreql)
                      wl= -(zvalue-sfreqr)/(sfreqr-sfreql)

                      do icurr=1,scal%nprd(ipe)
                        gmeRDyntempr(icurr,1)=pol%gme(tmprowindex(icurr), &
                          i_myband,tmp_iv,ispin,irk,im) * tmprowph(icurr)
                      enddo

                      this%gmeRDyntempr2(:,mytot,isfreql)=gmeRDyntempr(:,1)*wl
                      this%gmeRDyntempr2(:,mytot,isfreqr)=gmeRDyntempr(:,1)*wr

                      do icurr=1,scal%npcd(ipe)
                        this%gmeRDyntempc(mytot,icurr) = &
                          MYCONJG( pol%gme(tmpcolindex(icurr),i_myband,tmp_iv,ispin,irk,im) * tmpcolph(icurr) )
                      enddo
                    endif
                  enddo ! i_myband
                  !$OMP END DO
                  !$OMP END PARALLEL
                  itot = itot+ilimit

                enddo ! iv
              enddo ! it
            enddo ! irk
          enddo ! im

          SAFE_DEALLOCATE(tmprowindex)
          SAFE_DEALLOCATE(tmpcolindex)
          SAFE_DEALLOCATE(tmprowph)
          SAFE_DEALLOCATE(tmpcolph)

          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_prep)
          endif

          !Do the zgemm`s
          if(ntot > 0) then
            do jj =1+peinf%rank_mtxel, pol%nsfreq,pol%nfreq_group

              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_gemm)
              endif

              call zgemm('n','n',scal%nprd(ipe),scal%npcd(ipe),ntot, &
                negfact,this%gmeRDyntempr2(:,:,jj),scal%nprd(ipe),this%gmeRDyntempc(:,:),ntot,&
                (0D0,0D0),this%chilocalRDyn(:,:,jj),scal%nprd(ipe))

              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_gemm)
              endif

            enddo
          endif

          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_comm)
          endif

          if(pol%nfreq_group .eq. 1) then
#ifdef MPI
            call MPI_Reduce(this%chilocalRDyn(1,1,1),this%chilocal2RDyn(1,1,1,ispin), &
              pol%os_nsfreq_para*scal%npcd(ipe)*scal%nprd(ipe),MPI_COMPLEX_DPC,&
              MPI_SUM,ipe-1,MPI_COMM_WORLD,mpierr)
#endif
          elseif(pol%nfreq_group .gt. 1 .and. peinf%inode .lt. peinf%npes) then
#ifdef MPI
            call MPI_Reduce(this%chilocalRDyn(1,1,1),this%chilocal2RDyn(1,1,1,ispin), &
              pol%os_nsfreq_para*scal%npcd(ipe)*scal%nprd(ipe),MPI_COMPLEX_DPC,&
              MPI_SUM,ipe-1,peinf%freq_comm,mpierr)
#endif
          endif
#ifndef MPI
          this%chilocal2RDyn(:,:,:,ispin)=this%chilocalRDyn(:,:,:)
#endif
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_comm)
          endif


        enddo ! ispin
          
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_array_alloc)
        endif

        SAFE_DEALLOCATE(this%chilocalRDyn)
        SAFE_DEALLOCATE(gmeRDyntempr)
        SAFE_DEALLOCATE(this%gmeRDyntempr2)
        SAFE_DEALLOCATE(this%gmeRDyntempc)

        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_array_alloc)
        endif

      enddo ! ipe
      call progress_free(prog_info)

      do ispin =1, kp%nspin
        do jj=1+peinf%rank_mtxel,pol%nsfreq,pol%nfreq_group
          pol%chiTDyn(jj,:,:,ispin) = this%chilocal2RDyn(:,:,jj,ispin)
        enddo ! jj
      enddo ! ispin
      SAFE_DEALLOCATE(this%chilocal2RDyn)

      ! -----------------------------
      ! Hilbert transform
      ! ------------------------------


! JRD XXX - pol%chi... arrays need to be reordered here.
! I think we specifying LDC wrong here - should pol%nfreq_in_group not pol%nfreq right?
      call zgemm('n','n',pol%nfreq,scal%npr*scal%npc*kp%nspin,pol%os_nsfreq_para, &
        (-1D0,0D0),htwR(:,:),pol%nfreq,pol%chiTDyn(:,:,:,:),pol%os_nsfreq_para, &
        (0D0,0D0),this%chiRDynOld(:,:,:,:),pol%nfreq)

      do ispin =1, kp%nspin
        do jj=1,pol%nfreq_in_group
          pol%chiRDyn(:,:,jj,ispin) = this%chiRDynOld(jj,:,:,ispin)
        enddo ! jj
      enddo ! ispin


    endif ! pol%freq_dep.eq.2.and.pol%freq_dep_method.eq.1

    !call free_summation_buffers(pol)

    POP_SUB(chi_summation_comm_matrix)
    return

  end subroutine chi_summation_comm_matrix

  !-----------------------------------------------------------------------------
  !                              GCOMM_ELEMENTS
  !-----------------------------------------------------------------------------

  !> Create the pol%chi matrix using gcomm_elements sceheme
  subroutine chi_summation_comm_elements(this,pol,scal,kp,vwfn,cwfn,&
    nst,nrk,indt,pht)
    type(chi_summator_t), intent(INOUT) :: this
    type(polarizability), intent(INOUT) :: pol
    type(scalapack), intent(in) :: scal
    type(kpoints), intent(IN) :: kp
    type(valence_wfns), intent(IN) :: vwfn
    type(conduction_wfns), intent(IN) :: cwfn

    integer, intent(IN) :: nst(:)
    integer, intent(IN) :: nrk
    integer, intent(INOUT) :: indt(:,:,:)
    SCALAR,  intent(INOUT) :: pht(:,:,:)

    integer :: icurr,ntot,itot
    integer :: iv, ic, irk, it, ispin

    SCALAR :: temp_gme
    integer, allocatable :: iowna(:)
    integer :: isend

    complex(DPC), allocatable :: edenDRtemp(:)

    real(DP) :: zvalue
    ! frequency points for the spectral functions of the polarizability
    integer :: isfreql, isfreqr
    real(DP) :: sfreql, sfreqr, wr, wl
    ! Hilbert tranform coefficients
    complex(DPC) :: htwR(pol%nfreq,pol%nsfreq), htwA(pol%nfreq,pol%nsfreq)
    complex(DPC) :: c11,c12,c13,c14,c21,c22,c23,c24
    complex(DPC) :: cA11,cA12,cA13,cA14,cA21,cA22,cA23,cA24

    integer :: isf,nf,nsf,ifreq_para
    real(DP) :: sf1,sf2
    real(DP) :: step1,step2,fqt,eta
    complex(DPC) :: c1,c2,j_dpc,cA1,cA2

    integer :: ii, jj
    type(progress_info) :: prog_info

    integer :: nsftot, il, ir, n1, n2, n3, max_nv, nwarn
    integer, allocatable :: count_v(:), ind_v(:,:), ind_sf(:)

    PUSH_SUB(chi_summation_comm_elements)

    !call alloc_summation_buffers(pol, this%fact)

    if(peinf%inode .eq. 0) then
      call timing%start(timing%chi_sum_array_alloc)
    endif

    SAFE_ALLOCATE(iowna, (vwfn%nband+pol%ncrit))

    if (pol%freq_dep .eq. 0) then
      SAFE_ALLOCATE(this%chilocal, (scal%npr,scal%npc))
    endif

    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
      SAFE_ALLOCATE(this%chilocalRDyn, (scal%npr,scal%npc,pol%nfreq_in_group))
    endif

    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      SAFE_ALLOCATE(this%chilocalRDyn, (scal%npr,scal%npc,pol%os_nsfreq_para))

      SAFE_ALLOCATE(count_v, (pol%os_nsfreq_para))
      SAFE_ALLOCATE(ind_sf, (pol%os_nsfreq_para))
    endif

    if(peinf%inode .eq. 0) then
      call timing%stop(timing%chi_sum_array_alloc)
    endif

    call logit("Starting chi Sum")
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      ! -------------------------------------------
      ! compute the Hilbert transform coefficients
      ! -------------------------------------------

      if(peinf%inode .eq. 0) then
        call timing%start(timing%chi_sum_ht_nb)
      endif
  
      j_dpc=(0.0,1.0)
      htwR(:,:)=(0.0,0.0)

      nf=pol%nfreq
      nsf=pol%nsfreq
      do jj=1,nf
        eta=pol%dBrdning/ryd
        fqt=pol%dFreqGrid(jj)/ryd
        do isf=1,nsf
          if (isf==1) then
            c1=(0.0,0.0)
            step1=1.0
          else
            sf1=pol%dSFreqGrid(isf-1)/ryd
            sf2=pol%dSFreqGrid(isf)/ryd
            step1=sf2-sf1

            c11=((sf1-fqt)*(-1.0*j_dpc)-eta)*(atan((fqt-sf2)/(eta))&
              &          -atan((fqt-sf1)/(eta)))
            c12=0.50*(sf1-fqt+(-1.0*j_dpc)*eta)*log(((fqt-sf2)**2.0+eta*eta)&
              &          /(((fqt-sf1)**2.0+eta*eta)))

            c13=-((sf1+fqt)*j_dpc-eta)*(atan((fqt+sf2)/(eta))&
              &          -atan((fqt+sf1)/(eta)))
            c14=0.50*(sf1+fqt+j_dpc*eta)*log(((fqt+sf2)**2.0+eta*eta)&
              &          /(((fqt+sf1)**2.0+eta*eta)))

            c1=c11+c12+c13+c14
            c1=c1/step1
          endif

          if (isf==nsf) then
            c2=(0.0,0.0)
            step2=1.0
          else
            sf1=pol%dSFreqGrid(isf)/ryd
            sf2=pol%dSFreqGrid(isf+1)/ryd

            step2=sf2-sf1

            c21=((sf2-fqt)*(-1.0*j_dpc)-eta)*(atan((fqt-sf1)/(eta))&
              &          -atan((fqt-sf2)/(eta)))
            c22=0.50*(sf2-fqt+(-1.0*j_dpc)*eta)*log(((fqt-sf1)**2.0+eta*eta)&
              &          /(((fqt-sf2)**2.0+eta*eta)))

            c23=-((sf2+fqt)*j_dpc-eta)*(atan((fqt+sf1)/(eta))&
              &          -atan((fqt+sf2)/(eta)))
            c24=0.50*(sf2+fqt+j_dpc*eta)*log(((fqt+sf1)**2.0+eta*eta)&
              &          /(((fqt+sf2)**2.0+eta*eta)))

            c2=c21+c22+c23+c24

            c2=c2/step2
          endif

          if (isf==1.or.isf==nsf) then
            htwR(jj,isf)=0.5d0*(c1/step1+c2/step2)
          else
            htwR(jj,isf)=1.0d0*(c1+c2)/(step1+step2)
          endif

        enddo
      enddo

      if(peinf%inode .eq. 0) then
        call timing%stop(timing%chi_sum_ht_nb)
      endif
  
    endif!(pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)

    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      pol%chiTDyn=(0.0,0.0)
    endif

    call progress_init(prog_info, 'building polarizability matrix', 'blocks', &
      nrk*kp%nspin*(cwfn%nband-vwfn%nband))
    nwarn=0
    do irk=1,nrk
      if (peinf%verb_debug .and. peinf%inode==0) then
        write(6,'(A,i8,6x,A,i8,A)') '### irk=',irk,'(nrk=',nrk,')'
      endif
#ifdef MPI
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif
      do ispin=1,kp%nspin
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_array_alloc)
        endif
        iowna(:)=1
        ntot=(vwfn%nband+pol%ncrit)*nst(irk)
        if (pol%freq_dep .eq. 0) then
!JRD XXX Should thread
          this%chilocal=0
          SAFE_ALLOCATE(this%gmetempr, (scal%npr,ntot))
          SAFE_ALLOCATE(this%gmetempc, (ntot,scal%npc))
        endif
        if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
!JRD XXX Should thread
          this%chilocalRDyn=0
          SAFE_ALLOCATE(this%gmeRDyntempr2, (scal%npr,ntot,pol%nfreq_in_group))
          SAFE_ALLOCATE(this%gmeRDyntempc, (ntot,scal%npc))
        endif

        if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
!JRD XXX Should thread
          this%chilocalRDyn=(0.0,0.0)

          max_nv=0
          do ic=1,cwfn%nband-vwfn%nband
            count_v=0
            do iv=1,(vwfn%nband+pol%ncrit)

              isend=peinf%global_pairowner(iv,ic)-1
              if (isend .lt. 0) then
                write(0,*) 'Illegal value for mpi proc, isend: ',iv,ic
                call die("internal error in chi_summation")
              endif
              if (isend .eq. peinf%inode) then
                zvalue=-pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1)
                iowna(iv)=iowna(iv)+1
              endif

#ifdef MPI
              call MPI_Bcast(zvalue,1,MPI_REAL_DP,isend,MPI_COMM_WORLD,mpierr)
              call MPI_Bcast(pol%dSFreqGrid,pol%os_nsfreq_para,MPI_REAL_DP,isend,MPI_COMM_WORLD,mpierr)
#endif

              if (scal%npr*scal%npc .ne. 0) then
                do it=1, nst(irk)

                  if (abs(zvalue) .gt. Tol_Zero) then
                    isfreql=-1
                    do jj=pol%nsfreq,1,-1
                      if ((pol%dSFreqGrid(jj)/ryd)<zvalue) then
                        isfreql=jj
                        EXIT
                      endif
                    enddo

                    if (isfreql.eq.pol%nsfreq) then
                      nwarn=nwarn+1
                      if (nwarn==1.and.peinf%inode.eq.0) then
                        write(0,*) 'WARNING: for accuracy, sfrequency_high_cutoff should be '
                        write(0,*) 'larger than energy difference between highest unoccupied '
                        write(0,*) 'state and lowest occupied state.'
                      endif
                      cycle
                    endif

                    isfreqr=isfreql+1

                    count_v(isfreql)=count_v(isfreql)+1
                    count_v(isfreqr)=count_v(isfreqr)+1
                  endif
                enddo !it
              endif
            enddo !iv

            if (max_nv<maxval(count_v(:))) then
              max_nv=maxval(count_v(:))
            endif
          enddo !ic

          SAFE_ALLOCATE(this%gmeRDyntempr2, (scal%npr,max_nv,pol%os_nsfreq_para))
          SAFE_ALLOCATE(this%gmeRDyntempc, (ntot,scal%npc))
          SAFE_ALLOCATE(this%gmeRDyntempcs, (max_nv,scal%npc))

          this%gmeRDyntempr2=(0.0,0.0)
          this%gmeRDyntempc=(0.0,0.0)
          this%gmeRDyntempcs=(0.0,0.0)

          SAFE_ALLOCATE(ind_v, (pol%os_nsfreq_para, max_nv))

          this%gmeRDyntempr2=(0.0,0.0)
          iowna(:)=1
        endif!(pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)

        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_array_alloc)
        endif

        do ic=1,cwfn%nband-vwfn%nband
          call progress_step(prog_info)

          ! We do two giant loops here for freq_dep cases
          if (pol%freq_dep .eq. 0) then
            itot=0
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_comm)
              endif
              do iv=1,(vwfn%nband+pol%ncrit)
              isend=peinf%global_pairowner(iv,ic)-1
              if (isend .lt. 0) then
                write(0,*) 'Illegal value for mpi proc, isend:',&
                  peinf%inode,iv,ic,peinf%global_pairowner(iv,ic)
                call die("internal error in chi_summation")
              endif
              if (isend .eq. peinf%inode) then
                if (iowna(iv) .gt. peinf%ncownactual) call die('iowna(iv) bigger than ncownactual')
                this%gmetempn(:) = pol%gme(:,iowna(iv),peinf%indexv(iv), &
                  ispin,irk,1) * sqrt(this%fact)
                iowna(iv)=iowna(iv)+1
              endif
#ifdef MPI
              call MPI_Bcast(this%gmetempn,pol%nmtx,MPI_SCALAR,isend,MPI_COMM_WORLD,mpierr)
#endif
              if (scal%npr*scal%npc .ne. 0) then

                do it =1, nst(irk)
                  itot = itot + 1

                  do icurr=1,scal%npr
                    this%gmetempr(icurr,itot)=this%gmetempn(indt(scal%imyrow(icurr),it,irk)) * &
                      pht(scal%imyrow(icurr),it,irk)
                  enddo

                  do icurr=1,scal%npc
                    temp_gme = this%gmetempn(indt(scal%imycol(icurr),it,irk))
                    this%gmetempc(itot,icurr)=MYCONJG(temp_gme * pht(scal%imycol(icurr),it,irk) )
                  enddo
                enddo ! it
              endif
            enddo ! iv
            if(peinf%inode .eq. 0) then
              call timing%stop(timing%chi_sum_comm)
            endif


            ! JRD: Using Level3 BLAS here for better performance

            if (scal%npr*scal%npc .ne. 0 .and. ntot > 0) then
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_gemm)
              endif

              call X(gemm)('n','n',scal%npr,scal%npc,ntot, &
                -ONE,this%gmetempr,scal%npr,this%gmetempc,ntot,ONE,this%chilocal,scal%npr)

              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_gemm)
              endif
  
            endif

          endif ! pol%freq_dep .eq. 0

          !---------------------
          ! JRD: Full Frequency Be Here

          if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then

            if(peinf%inode .eq. 0) then
              call timing%start(timing%chi_sum_array_alloc)
            endif
            SAFE_ALLOCATE(edenDRtemp, (pol%nfreq))
            if(peinf%inode .eq. 0) then
              call timing%stop(timing%chi_sum_array_alloc)
            endif

            itot=0
!JRD XXX Should thread
            do iv=1,(vwfn%nband+pol%ncrit)
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_comm)
              endif
  
              isend=peinf%global_pairowner(iv,ic)-1
              if (isend .lt. 0) then
                write(0,*) 'Illegal value for mpi proc, isend: ',iv,ic
                call die("internal error in chi_summation")
              endif
              if (isend .eq. peinf%inode) then
                this%gmeRDyntempn(:) = pol%gme(:,iowna(iv),peinf%indexv(iv), &
                  ispin,irk,1) * sqrt(this%fact)
                if (abs(pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1)) .gt. Tol_Zero) then
                  do jj=1,pol%nfreq
                    if (pol%tda) then
                      edenDRtemp(jj)= -0.50d0 * ( 1d0/(pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1)+&
                        (pol%dFreqGrid(jj)+pol%dFreqBrd(jj))/ryd))
                    else
                      ! DVF partial occupations: occ_diff factor added here
                      edenDRtemp(jj)= -0.50d0 * pol%occ_diff(peinf%indexv(iv),iowna(iv),ispin,irk,1) * &
                        ( 1d0/(pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1) &
                        -(pol%dFreqGrid(jj)+pol%dFreqBrd(jj))/ryd)+&
                        1d0/(pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1)+&
                        (pol%dFreqGrid(jj)+pol%dFreqBrd(jj))/ryd))
                    endif
                  enddo
                else
                  edenDRtemp(:)=0D0
                endif
                iowna(iv)=iowna(iv)+1
              endif

#ifdef MPI
              call MPI_Bcast(this%gmeRDyntempn,pol%nmtx,MPI_COMPLEX_DPC,isend,MPI_COMM_WORLD,mpierr)
              call MPI_Bcast(edenDRtemp,pol%nfreq,MPI_COMPLEX_DPC,isend,MPI_COMM_WORLD,mpierr)
#endif
              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_comm)
              endif

              if (scal%npr*scal%npc .ne. 0) then
                do it =1, nst(irk)
                  if(peinf%inode .eq. 0) then
                    call timing%start(timing%chi_sum_row)
                  endif
    
                  itot = itot + 1
                  do icurr=1,scal%npr
                    do jj =1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
                      ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
                      this%gmeRDyntempr2(icurr,itot,ifreq_para)= &
                        (this%gmeRDyntempn(indt(scal%imyrow(icurr),it,irk))*pht(scal%imyrow(icurr),it,irk))*edenDRtemp(jj)
                    enddo
                  enddo

                  if(peinf%inode .eq. 0) then
                    call timing%stop(timing%chi_sum_row)
                  endif
  
                  if(peinf%inode .eq. 0) then
                    call timing%start(timing%chi_sum_column)
                  endif
  
                  do icurr=1,scal%npc
                    this%gmeRDyntempc(itot,icurr) = &
                      CONJG(this%gmeRDyntempn(indt(scal%imycol(icurr),it,irk))*pht(scal%imycol(icurr),it,irk))
                  enddo

                  if(peinf%inode .eq. 0) then
                    call timing%stop(timing%chi_sum_column)
                  endif
  
                enddo ! it
              endif

            enddo ! iv


            ! JRD: Using Level3 BLAS here for better performance

            if (scal%npr*scal%npc .ne. 0 .and. ntot > 0) then
              do jj =1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
                if(peinf%inode .eq. 0) then
                  call timing%start(timing%chi_sum_gemm)
                endif

                ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
                call zgemm('n','n',scal%npr,scal%npc,ntot,(-1D0,0D0),this%gmeRDyntempr2(:,:,ifreq_para),scal%npr, &
                  this%gmeRDyntempc(:,:),ntot,(1D0,0D0),this%chilocalRDyn(:,:,ifreq_para),scal%npr)
                if(peinf%inode .eq. 0) then
                  call timing%stop(timing%chi_sum_gemm)
                endif
  
              enddo
            endif

            if(peinf%inode .eq. 0) then
              call timing%start(timing%chi_sum_array_alloc)
            endif

            SAFE_DEALLOCATE(edenDRtemp)

            if(peinf%inode .eq. 0) then
              call timing%stop(timing%chi_sum_array_alloc)
            endif

          endif ! (pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0)

          !---------------------
          ! Full Frequency Be Here(shishkin and Kresse 2006)

          if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
            count_v=0
            ind_v=0
            ind_sf=0
            nsftot=0
            itot=0

            do iv=1,(vwfn%nband+pol%ncrit)

              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_comm)
              endif
  
              isend=peinf%global_pairowner(iv,ic)-1
              if (isend .lt. 0) then
                write(0,*) 'Illegal value for mpi proc, isend: ',iv,ic
                call die("internal error in chi_summation")
              endif
              if (isend .eq. peinf%inode) then
                this%gmeRDyntempn(:) = pol%gme(:,iowna(iv),peinf%indexv(iv), &
                  ispin,irk,1) * sqrt(this%fact)
                zvalue=-pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1)
                iowna(iv)=iowna(iv)+1
              endif

#ifdef MPI
              call MPI_Bcast(this%gmeRDyntempn,pol%nmtx,MPI_COMPLEX_DPC,isend,MPI_COMM_WORLD,mpierr)
              call MPI_Bcast(zvalue,1,MPI_REAL_DP,isend,MPI_COMM_WORLD,mpierr)
              call MPI_Bcast(pol%dSFreqGrid,pol%os_nsfreq_para,MPI_REAL_DP,isend,MPI_COMM_WORLD,mpierr)
#endif

              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_comm)
              endif
              ! compute spectral functions of the polarizability

              if (scal%npr*scal%npc .ne. 0) then
                do it=1, nst(irk)

                  if (abs(zvalue) .gt. Tol_Zero) then
                    
                    if(peinf%inode .eq. 0) then
                      call timing%start(timing%chi_sum_row)
                    endif

                    itot=itot+1
                    isfreql=-1
                    do jj=pol%nsfreq,1,-1
                      if ((pol%dSFreqGrid(jj)/ryd)<zvalue) then
                        isfreql=jj
                        EXIT
                      endif
                    enddo

                    if (isfreql.eq.pol%nsfreq) then
                      cycle
                    endif

                    isfreqr=isfreql+1

                    count_v(isfreql)=count_v(isfreql)+1
                    count_v(isfreqr)=count_v(isfreqr)+1

                    il=count_v(isfreql)
                    ir=count_v(isfreqr)

                    ind_v(isfreql,il)=itot
                    ind_v(isfreqr,ir)=itot

                    sfreql=pol%dSFreqGrid(isfreql)/ryd
                    sfreqr=pol%dSFreqGrid(isfreqr)/ryd

                    wl=-(zvalue-sfreqr)/(sfreqr-sfreql)
                    wr=(zvalue-sfreql)/(sfreqr-sfreql)

                    do icurr=1,scal%npr
                      this%gmeRDyntempr2(icurr,il,isfreql)=this%gmeRDyntempn( &
                        indt(scal%imyrow(icurr),it,irk))*pht(scal%imyrow(icurr),it,irk)*wl
                      this%gmeRDyntempr2(icurr,ir,isfreqr)=this%gmeRDyntempn( &
                        indt(scal%imyrow(icurr),it,irk))*pht(scal%imyrow(icurr),it,irk)*wr
                    enddo

                    if(peinf%inode .eq. 0) then
                      call timing%stop(timing%chi_sum_row)
                    endif
  
                    if(peinf%inode .eq. 0) then
                      call timing%start(timing%chi_sum_column)
                    endif

                    do icurr=1,scal%npc
                      this%gmeRDyntempc(itot,icurr) = &
                        CONJG(this%gmeRDyntempn(indt(scal%imycol(icurr),it,irk))*pht(scal%imycol(icurr),it,irk))
                    enddo

                    if(peinf%inode .eq. 0) then
                      call timing%stop(timing%chi_sum_column)
                    endif
  
                  endif

                enddo ! it
              endif

            enddo ! iv

            if(peinf%inode .eq. 0) then
              call timing%start(timing%chi_sum_array_alloc)
            endif

            jj=0
            do ii=1+peinf%rank_mtxel,pol%nsfreq,pol%nfreq_group
              if (count_v(ii)>0) then
                jj=jj+1
                ind_sf(jj)=ii
              endif
            enddo
            nsftot=jj
            if(peinf%inode .eq. 0) then
              call timing%stop(timing%chi_sum_array_alloc)
            endif


            if (scal%npr*scal%npc .ne. 0 .and. ntot > 0) then
              do ii=1, nsftot
                if(peinf%inode .eq. 0) then
                  call timing%start(timing%chi_sum_gemm)
                endif
    
                n1=ind_sf(ii)
                n2=count_v(n1)

                do jj=1,n2
                  n3=ind_v(n1,jj)
                  this%gmeRDyntempcs(jj,:)=this%gmeRDyntempc(n3,:)
                enddo

                call zgemm('n','n',scal%npr,scal%npc,n2,(-1D0,0D0),this%gmeRDyntempr2(:,:,n1),scal%npr, &
                  this%gmeRDyntempcs(:,:),max_nv,(1D0,0D0),this%chilocalRDyn(:,:,n1),scal%npr)

                if(peinf%inode .eq. 0) then
                  call timing%stop(timing%chi_sum_gemm)
                endif
  
              enddo
            endif
          endif ! (pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)

        enddo ! ic (loop over conduction bands)


        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_array_alloc)
        endif
        if (pol%freq_dep .eq. 0) then
          if (scal%npr*scal%npc .ne. 0) then
            pol%chi(:,:,ispin) = pol%chi(:,:,ispin) + this%chilocal(:,:)
          endif
          SAFE_DEALLOCATE(this%gmetempr)
          SAFE_DEALLOCATE(this%gmetempc)
        endif

        if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
          if (scal%npr*scal%npc .ne. 0) then
            do jj = 1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
              ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
! JRD XXX This copy is now a waste of time. Should set pol%chi directly above in the zgemm
              pol%chiRDyn(:,:,ifreq_para,ispin) = pol%chiRDyn(:,:,ifreq_para,ispin) + this%chilocalRDyn(:,:,ifreq_para)
            enddo
          endif
          SAFE_DEALLOCATE(this%gmeRDyntempr2)
          SAFE_DEALLOCATE(this%gmeRDyntempc)
        endif
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_array_alloc)
        endif


        if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then

          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_array_alloc)
          endif
  
          if (scal%npr*scal%npc .ne. 0) then
            do jj = 1+peinf%rank_mtxel, pol%nsfreq,pol%nfreq_group
              pol%chiTDyn(jj,:,:,ispin) = pol%chiTDyn(jj,:,:,ispin) + this%chilocalRDyn(:,:,jj)
            enddo
          endif
          SAFE_DEALLOCATE(this%gmeRDyntempr2)
          SAFE_DEALLOCATE(this%gmeRDyntempc)
          SAFE_DEALLOCATE(this%gmeRDyntempcs)
          SAFE_DEALLOCATE(ind_v)
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_array_alloc)
          endif

        endif

      enddo ! ispin (loop over spins)
    enddo ! irk (loop over k-points in set rk)
    call progress_free(prog_info)

    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then

      ! -------------------------
      ! Hilbert transform
      ! -------------------------

      if(peinf%inode .eq. 0) then
        call timing%start(timing%chi_sum_ht_nb)
      endif


! JRD XXX pol%chi arrays out of order below
      call zgemm('n','n',pol%nfreq,scal%npr*scal%npc*kp%nspin,pol%os_nsfreq_para, &
        (-1D0,0D0),htwR(:,:),pol%nfreq,pol%chiTDyn(:,:,:,:),pol%os_nsfreq_para, &
        (0D0,0D0),this%chiRDynOld(:,:,:,:),pol%nfreq)

      do ispin =1, kp%nspin
        do jj=1,pol%nfreq_in_group
          pol%chiRDyn(:,:,jj,ispin) = this%chiRDynOld(jj,:,:,ispin)
        enddo ! jj
      enddo ! ispin

      if(peinf%inode .eq. 0) then
        call timing%stop(timing%chi_sum_ht_nb)
      endif
  

    endif !(pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)

    if(peinf%inode .eq. 0) then
      call timing%start(timing%chi_sum_array_alloc)
    endif

    if (pol%freq_dep .eq. 0) then
      SAFE_DEALLOCATE(this%chilocal)
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
      SAFE_DEALLOCATE(this%chilocalRDyn)
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      SAFE_DEALLOCATE(this%chilocalRDyn)
      SAFE_DEALLOCATE(count_v)
      SAFE_DEALLOCATE(ind_sf)
    endif

    SAFE_DEALLOCATE(iowna)
    !call free_summation_buffers(pol)

    if(peinf%inode .eq. 0) then
      call timing%stop(timing%chi_sum_array_alloc)
    endif


    POP_SUB(chi_summation_comm_elements)
    return

  end subroutine chi_summation_comm_elements

#if defined MPI && defined USESCALAPACK
  subroutine chi_summation_sub_trunc(this,pol,scal,kp,kpq,vwfn,cwfn,&
                                     nst,nrk,indt,pht,gvec,crys,q0,iq)
    use inversion_m
    type(chi_summator_t), intent(INOUT) :: this
    type(polarizability), intent(INOUT) :: pol
    type(scalapack), intent(in) :: scal
    type(kpoints), intent(IN) :: kp, kpq
    type(valence_wfns), intent(IN) :: vwfn
    type(conduction_wfns), intent(IN) :: cwfn

    integer, intent(IN) :: nst(:)
    integer, intent(IN) :: nrk
    integer, intent(INOUT) :: indt(:,:,:)
    SCALAR,  intent(INOUT) :: pht(:,:,:)
    type (gspace), intent(in) :: gvec
    type (crystal), intent(in) :: crys
    real(DP), intent(in) :: q0(3)
    integer, intent(in)  :: iq
    real(DP) :: zvalue, cv_energy, occ_diff

    SCALAR, ALLOCATABLE :: chilocal(:,:), chiRDyn_local(:,:)
    SCALAR, allocatable :: chi_omega0_save(:,:,:)
    INTEGER :: ntot, ntot2
    INTEGER :: neig_sub
    INTEGER :: irk, ipe, ispin, itot, it, iv, i, j, ii, jj
    INTEGER :: irow, icol, irowm, icolm, icurr, ilimit, freq_idx, &
               i_myband, mytot
    integer :: ig_l, ig_g, igp_l, igp_g
    INTEGER :: ipe_real

    integer, allocatable :: tmprowindex(:),tmpcolindex(:)
    SCALAR, allocatable :: tmprowph(:),tmpcolph(:)
    type(cvpair_info) :: cvpair_temp
    SCALAR :: negfact
    SCALAR :: edenDRtemp

    SCALAR :: epsheaddummy, wcoul0
    type (twork_scell) :: work_scell
    real(DP), allocatable :: vcoul(:)
    integer, allocatable :: isorti(:)
    real(DP) :: vc, oneoverq, avgcut
    integer :: iscreen, nfk, iparallel
    integer :: qgrid(3)
    real(DP) :: q0vec(3)

    ! variables for non-blocking cyclic scheme
    integer :: isend_static, irec_static
    integer :: actual_send, actual_rec
    integer :: nsend_row, nsend_col, nrec_row, nrec_col
    integer :: req_send, tag_send, req_rec, tag_rec
    integer :: stat(MPI_STATUS_SIZE)
    SCALAR, allocatable :: buf_send(:,:,:), buf_rec(:,:,:), buf_temp(:,:,:)

    type(progress_info) :: prog_info
    ! variables for fixed buffer size
    logical :: keep_fix_buf
    integer :: size_N, size_M, size_K
    ! variables for matrix dumping
    logical :: do_read_dump_matrix, do_dump_matrix

    PUSH_SUB(chi_summation_sub_trunc)

    ! Keep these variables hardcoded, can be used by advance
    ! user to dump the block cyclic distributed polarizability
    ! as check point, for specific format see the 2 routines
    ! at the bottom of the module (dump_matrix and read_dump_matrix)
    do_dump_matrix = .false.
    do_read_dump_matrix = .false.
    ! if we read dumped matrix make sure we don't write, so we don't
    ! override the results
    if ( do_read_dump_matrix ) do_dump_matrix = .false.

    ! this is the new communication scheme which avoid continuos allocation
    ! deallocation of the communication / computation buffers. it will be
    ! set default to .true. and turn off by using the correspondinf input key
    keep_fix_buf = .false.
    ! keep_fix_buf only work for pol%nonblocking_cyclic
    if( pol%nonblocking_cyclic ) then
      keep_fix_buf = .true.
      if( pol%dont_keep_fix_buffers ) then
         keep_fix_buf = .false.
      end if
    end if

    !XXXXXXXXX
    ! WRITE(2000,*)  pol%gme
    ! WRITE(2000,*)  'XXXX'
    ! WRITE(2000,*)  pht
    ! WRITE(*,*) SIZE(pol%gme,1)
    ! WRITE(*,*) SIZE(pol%gme,2)
    ! WRITE(*,*) SIZE(pol%gme,3)
    ! WRITE(*,*) SIZE(pol%gme,4)
    ! WRITE(*,*) SIZE(pol%gme,5)
    ! WRITE(*,*) SIZE(pol%gme,6)
    ! DO irk = 1, nrk
    !   DO j = 1, SIZE(pol%gme,3)
    !     DO i = 1, SIZE(pol%gme,2)
    !       WRITE(2000,*) pht(i,j,irk)
    !       WRITE(2000,*) pol%gme(:,i,j,1,irk,1)
    !     END DO
    !   END DO
    ! END DO
    !XXXXXXXXX

    ! first calculate polarizability at omega=0 using all processes (similar to static case)
    ntot=0
    ntot2=0
    ntot = peinf%nvownactual*peinf%ncownactual
    do irk = 1, nrk
      ntot2=ntot2 + nst(irk)
    enddo
    ntot=ntot*ntot2

    SAFE_ALLOCATE(this%chilocal2, (scal%npr,scal%npc,kp%nspin))
    this%chilocal2=0

    freq_idx = pol%nfreq_group
    IF(pol%nfreq_group .gt. 1) THEN
      freq_idx = peinf%rank_mtxel+1
    END IF

    negfact = -1D0*this%fact

    ! trace the time for chi omega=0
    if(peinf%inode .eq. 0) then
      call timing%start(timing%chi_sum_sub_omega_0)
    endif


    ! here we want to use all processes for computing omega=0, we
    ! do not take care of the parallelization over frequencies
    call progress_init(prog_info, 'building polarizability matrix omega=0', 'processor', &
      peinf%npes)

    if(pol%nonblocking_cyclic) then
      ! initialize buffer for non-blocking cyclic scheme
      ! static process for the communication
      isend_static = MOD(peinf%inode + 1 + peinf%npes, peinf%npes)
      irec_static = MOD(peinf%inode - 1 + peinf%npes, peinf%npes)
      ! allocate my size for the first send
      nsend_row = scal%nprd(peinf%inode+1)
      nsend_col = scal%npcd(peinf%inode+1)

      if(keep_fix_buf) then
        ! precompute max sizes
        size_N = MAXVAL(scal%nprd(:))
        size_M = MAXVAL(scal%npcd(:))
        size_K = ntot
        call mpi_allreduce(MPI_IN_PLACE, size_K, 1, &
                           MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
        ! allocate only once
        SAFE_ALLOCATE(buf_send, (size_N, size_M, kp%nspin))
        SAFE_ALLOCATE(buf_rec,  (size_N, size_M, kp%nspin))
        SAFE_ALLOCATE(buf_temp, (size_N, size_M, kp%nspin))
        do ispin = 1 , kp%nspin
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, size_M
            do jj = 1, size_N
              buf_send(jj,ii,ispin) = ZERO
            enddo
          enddo
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, size_M
            do jj = 1, size_N
              buf_rec(jj,ii,ispin)  = ZERO
            end do
          end do
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, size_M
            do jj = 1, size_N
              buf_temp(jj,ii,ispin) = ZERO
            end do
          end do
        end do
        ! allocate for index mapping
        SAFE_ALLOCATE(tmprowindex, (size_N)) ! (scal%nprd(ipe_real)))
        SAFE_ALLOCATE(tmpcolindex, (size_M)) ! (scal%npcd(ipe_real)))
        SAFE_ALLOCATE(tmprowph,  (size_N)) ! (scal%nprd(ipe_real)))
        SAFE_ALLOCATE(tmpcolph,  (size_M)) ! (scal%npcd(ipe_real)))
        ! allocate for matrix multiplication
        SAFE_ALLOCATE(this%chilocal, (size_N, size_M)) ! (scal%nprd(ipe),scal%npcd(ipe))
        SAFE_ALLOCATE(this%gmetempr, (size_N, size_K)) ! (scal%nprd(ipe),ntot)
        SAFE_ALLOCATE(this%gmetempc, (size_K, size_M)) ! (ntot,scal%npcd(ipe))
        !$OMP PARALLEL DO collapse(2)
        do ii = 1, size_M
          do jj = 1, size_N
            this%chilocal(jj,ii) = ZERO
          end do
        end do
        !$OMP PARALLEL DO collapse(2)
        do ii = 1, size_K
          do jj = 1, size_N
            this%gmetempr(jj,ii) = ZERO
          end do
        end do
        !$OMP PARALLEL DO collapse(2)
        do ii = 1, size_M
          do jj = 1, size_K
            this%gmetempc(jj,ii) = ZERO
          end do
        end do
        !
      else  ! keep_fix_buf
        !
        ! standard case
        SAFE_ALLOCATE(buf_send, (nsend_row, nsend_col, kp%nspin))
        do ispin = 1 , kp%nspin
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, nsend_col
            do jj = 1, nsend_row
              buf_send(jj,ii,ispin) = ZERO
            enddo
          enddo
        end do
        !
      end if  ! keep_fix_buf
    end if

    do ipe = 1, peinf%npes
      call progress_step(prog_info)

      if ( do_read_dump_matrix ) cycle

      if(pol%nonblocking_cyclic) then
        ! calculate the actual process we have to send and we are receiving
        actual_send = MOD(peinf%inode + ipe + peinf%npes, peinf%npes)
        actual_rec  = MOD(peinf%inode - ipe + peinf%npes, peinf%npes)
        nrec_row = scal%nprd(actual_rec+1)
        nrec_col = scal%npcd(actual_rec+1)
        !
        if( keep_fix_buf ) then
          ! just make sure sizes will always be the same
          nsend_row = size_N
          nrec_row  = size_N
          nsend_col = size_M
          nrec_col  = size_M
        else
          ! allocate reciving buffer
          SAFE_ALLOCATE(buf_rec, (nrec_row,nrec_col,kp%nspin))
        end if
        !
        do ispin = 1 , kp%nspin
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, nrec_col
            do jj = 1, nrec_row
              buf_rec(jj,ii,ispin) = ZERO
            enddo
          enddo
        end do
        ! post message
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_comm)
        endif
        tag_rec = 1
        tag_send = 1
        req_rec = MPI_REQUEST_NULL
        req_send = MPI_REQUEST_NULL
        CALL MPI_Irecv(buf_rec, nrec_row*nrec_col*kp%nspin,MPI_SCALAR,irec_static,&
                       tag_rec, MPI_COMM_WORLD, req_rec, mpierr)
        CALL MPI_Isend(buf_send, nsend_row*nsend_col*kp%nspin,MPI_SCALAR,isend_static,&
                       tag_send, MPI_COMM_WORLD, req_send, mpierr)
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_comm)
        endif
              
        ! allocate working array
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_array_alloc)
        endif
        !
        if( .not. keep_fix_buf ) then
          SAFE_ALLOCATE(this%chilocal, (nrec_row, nrec_col))
        end if
        !
!$OMP PARALLEL DO collapse(2)
        do ii = 1, nrec_col
          do jj = 1, nrec_row
            this%chilocal(jj,ii)=0D0
          enddo
        enddo
        !
        if( .not. keep_fix_buf ) then
          SAFE_ALLOCATE(buf_temp, (nrec_row, nrec_col, kp%nspin))
        end if
        !
        do ispin = 1 , kp%nspin
!$OMP PARALLEL DO collapse(2)
          do ii = 1, nrec_col
            do jj = 1, nrec_row
              buf_temp(jj,ii,ispin) = ZERO
            enddo
          enddo
        end do
        !
        if( keep_fix_buf ) then
          ! make sure to put buffer to zero
          !
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, size_K
            do jj = 1, size_N
              this%gmetempr(jj,ii) = ZERO
            end do
          end do
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, size_M
            do jj = 1, size_K
              this%gmetempc(jj,ii) = ZERO
            end do
          end do
          !
        else
          SAFE_ALLOCATE(this%gmetempr, (nrec_row, ntot))
          SAFE_ALLOCATE(this%gmetempc, (ntot, nrec_col))
        end if
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_array_alloc)
        endif
        
        ! in this case the ipe is the real process from which we should receive
        ipe_real = actual_rec+1
      else
        ! standard calculation
        SAFE_ALLOCATE(this%chilocal, (scal%nprd(ipe),scal%npcd(ipe)))
        this%chilocal=0D0
        SAFE_ALLOCATE(this%gmetempr, (scal%nprd(ipe),ntot))
        SAFE_ALLOCATE(this%gmetempc, (ntot,scal%npcd(ipe)))
        ! this is the standard case
        ipe_real = ipe
      end if

      do ispin = 1 , kp%nspin

        itot = 0

        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_prep)
        endif

        if ( .not. keep_fix_buf ) then
          SAFE_ALLOCATE(tmprowindex,(scal%nprd(ipe_real)))
          SAFE_ALLOCATE(tmpcolindex,(scal%npcd(ipe_real)))
          SAFE_ALLOCATE(tmprowph,(scal%nprd(ipe_real)))
          SAFE_ALLOCATE(tmpcolph,(scal%npcd(ipe_real)))
        end if

        do irk = 1, nrk
          do it = 1, nst(irk)

            do icurr=1,scal%nprd(ipe_real)
              tmprowindex(icurr) = indt(scal%imyrowd(icurr,ipe_real),it,irk)
              tmprowph(icurr) = pht(scal%imyrowd(icurr,ipe_real),it,irk)
            enddo
            do icurr=1,scal%npcd(ipe_real)
              tmpcolindex(icurr) = indt(scal%imycold(icurr,ipe_real),it,irk)
              tmpcolph(icurr) = pht(scal%imycold(icurr,ipe_real),it,irk)
            enddo

            do iv = 1, vwfn%nband+pol%ncrit
               ! check if I own this partuicular band
               if (peinf%doiownv(iv)) then
                 ilimit = peinf%ncownactual
               else
                 ilimit = 0
               endif

               !$OMP PARALLEL private (mytot, zvalue, edenDRtemp, icurr, &
               !$OMP                   cvpair_temp, cv_energy )
               !$OMP DO
               do i_myband = 1, ilimit
                 mytot = itot + i_myband
                 zvalue = pol%edenDyn(peinf%indexv(iv),i_myband,ispin,irk,freq_idx)
! DVF partial occupations: Also, occ_diff factor added here. Same comment
! applies about OMP directive
                 occ_diff = pol%occ_diff(peinf%indexv(iv),i_myband,ispin,irk,freq_idx)
                 if(pol%lin_denominator<TOL_Zero) then
                   ! this is when the lin_denominator mode is not active.
                   if(abs(zvalue) .gt. Tol_Zero) then
                     if(pol%tda) then
                       edenDRtemp = -0.5D+00 * ( &
                         1d0/(zvalue + pol%dFreqBrd(1)/ryd))
                     else
                       edenDRtemp = -0.5D+00 * occ_diff * ( &
                         1d0/(zvalue - pol%dFreqBrd(1)/ryd)+ &
                         1d0/(zvalue + pol%dFreqBrd(1)/ryd))
                     endif
                   else
                     edenDRtemp = 0.0D+00
                   end if
                 else
                   !this is when lin_denominator mode is active
                   cvpair_temp = pol%lin_edenDyn(peinf%indexv(iv),i_myband,ispin,irk,freq_idx)
                   cv_energy = cvpair_temp%ev - cvpair_temp%ec

                   if(abs(cv_energy) - pol%lin_denominator > (-TOL_Zero)) then
                     !
                     if (abs(zvalue) .gt. Tol_Zero) then
                       if(pol%tda) then
                         edenDRtemp = -0.5d0*( &
                           1d0/(zvalue+pol%dFreqBrd(1)/ryd))
                       else
                         edenDRtemp = -0.5d0*( &
                           1d0/(zvalue-pol%dFreqBrd(1)/ryd)+ &
                           1d0/(zvalue+pol%dFreqBrd(1)/ryd))
                       endif
                     else
                       edenDRtemp = 0D0
                     endif
                     !
                   else ! if denominator is small, do linearized energy denominator
                     !
                     if(cvpair_temp%vltc) then
                       if (pol%tda) then
                         edenDRtemp = -0.5d0*(&
                           
                           integrate_mbz(kp,kpq%rk(1:3,cvpair_temp%idx_kp),&
                           cvpair_temp%ev,cvpair_temp%ec,&
                           cvpair_temp%vv,cvpair_temp%vc,&
                           pol%dFreqGrid(1)/ryd,pol%efermi/ryd))
                       else
                         edenDRtemp = -0.5d0*(&
                           
                           conjg(integrate_mbz(kp,kpq%rk(1:3,cvpair_temp%idx_kp),&
                           cvpair_temp%ev,cvpair_temp%ec,&
                           
                           cvpair_temp%vv,cvpair_temp%vc,&
                           -pol%dFreqGrid(1)/ryd,pol%efermi/ryd))&
                           
                           +integrate_mbz(kp,kpq%rk(1:3,cvpair_temp%idx_kp),&
                           cvpair_temp%ev,cvpair_temp%ec,&
                           cvpair_temp%vv,cvpair_temp%vc,&
                           pol%dFreqGrid(1)/ryd,pol%efermi/ryd))
                       endif
                     else
                       edenDRtemp = 0D0
                     endif
                     !
                   end if
                   !
                 end if ! lin denominator

                 ! divede the matrix elements by the corresponding energy denominator
                 ! right hand side matrix
                 do icurr=1,scal%nprd(ipe_real)
                   this%gmetempr(icurr,mytot)=pol%gme(tmprowindex(icurr), &
                     i_myband,peinf%indexv(iv),ispin,irk,freq_idx) * tmprowph(icurr) * edenDRtemp
                 enddo
                 ! left hand side matrix
                 do icurr=1,scal%npcd(ipe_real)
                   this%gmetempc(mytot,icurr) = &
                     MYCONJG( pol%gme(tmpcolindex(icurr),i_myband,peinf%indexv(iv),ispin,irk,freq_idx) * tmpcolph(icurr) )
                 enddo

               end do ! i_myband
               !$OMP END DO
               !$OMP END PARALLEL
               itot = itot+ilimit

            end do  ! iv

          end do   ! it
        end do   ! irk

        if ( .not. keep_fix_buf ) then
          SAFE_DEALLOCATE(tmprowindex)
          SAFE_DEALLOCATE(tmpcolindex)
          SAFE_DEALLOCATE(tmprowph)
          SAFE_DEALLOCATE(tmpcolph)
        end if

        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_prep)
        endif


        ! Go with matrix multiplication
        ! if we keep fixed buffer size we always multiply
        if( keep_fix_buf ) then
          !
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_gemm)
          endif

          call accel_enter_data_map_to(this%gmetempr, chi_summation_algo)
          call accel_enter_data_map_to(this%gmetempc, chi_summation_algo)
          call accel_enter_data_map_alloc(this%chilocal, chi_summation_algo)
          call accel_xgemm('n', 'n', &
                         size_N, size_M, size_K, &
                         negfact, &
                         this%gmetempr, size_N, &
                         this%gmetempc, size_K,&
                         ZERO, &
                         this%chilocal, size_N, &
                         chi_summation_algo)
          call accel_exit_data_map_from(this%chilocal, chi_summation_algo)
          call accel_exit_data_map_delete(this%gmetempc, chi_summation_algo)
          call accel_exit_data_map_delete(this%gmetempr, chi_summation_algo)

          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_gemm)
          endif

          !
        else
        !
        if(ntot > 0) then
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_gemm)
          endif

          ! in the case pol%nonblocking_cyclic, ipe is the shift wrt inode
          ! position, ipe_real is the absolute position of the ipe from which
          ! we are receiving. in the standard case ipe == ipe_real
#ifdef CPLX
          call zgemm('n','n',scal%nprd(ipe_real),scal%npcd(ipe_real),ntot, &
            negfact,this%gmetempr(:,:),scal%nprd(ipe_real),this%gmetempc(:,:),ntot,&
            (0D0,0D0),this%chilocal(:,:),scal%nprd(ipe_real))
#else
          call dgemm('n','n',scal%nprd(ipe_real),scal%npcd(ipe_real),ntot, &
                     negfact,this%gmetempr(:,:),scal%nprd(ipe_real),this%gmetempc(:,:),ntot,&
                     0D0,this%chilocal(:,:),scal%nprd(ipe_real))
#endif
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_gemm)
          endif

        endif
        !
        end if !keep_fix_buf

        if(pol%nonblocking_cyclic) then
          ! just copy the result
!$OMP PARALLEL DO collapse(2)
          do ii = 1, nrec_col
            do jj = 1, nrec_row
              buf_temp(jj,ii,ispin) = this%chilocal(jj,ii)
            enddo
          enddo
        else
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_comm)
          endif

#ifdef MPI
          ! here ipe == ipe_real
          call MPI_Reduce(this%chilocal(1,1),this%chilocal2(1,1,ispin),scal%npcd(ipe)*scal%nprd(ipe),MPI_SCALAR, &
                          MPI_SUM,ipe-1,MPI_COMM_WORLD,mpierr)
#else
          this%chilocal2(:,:,ispin)=this%chilocal(:,:)
#endif
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_comm)
          endif

        end if

      enddo ! ispin
      if(peinf%inode .eq. 0) then
        call timing%start(timing%chi_sum_array_alloc)
      endif

      !
      if( .not. keep_fix_buf ) then
        SAFE_DEALLOCATE(this%chilocal)
        SAFE_DEALLOCATE(this%gmetempr)
        SAFE_DEALLOCATE(this%gmetempc)
      end if
      !
      if(peinf%inode .eq. 0) then
        call timing%stop(timing%chi_sum_array_alloc)
      endif


      ! finalize non-blocking cyclic communication
      if(pol%nonblocking_cyclic) then
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_comm)
          endif
  
          ! make sure the buffer is received
          CALL MPI_Wait(req_rec,stat,mpierr)
          ! accumulate contribution into receiving buffer
          ! buf_rec(:,:,:) = buf_rec(:,:,:) + buf_temp
          do ispin = 1 , kp%nspin
!$OMP PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_rec(jj,ii,ispin) = buf_rec(jj,ii,ispin) + buf_temp(jj,ii,ispin)
              enddo
            enddo
          end do
          !
          if( .not. keep_fix_buf ) then
            SAFE_DEALLOCATE(buf_temp)
          end if
          !
          ! wait for the massage to be sent
          CALL MPI_Wait(req_send,stat,mpierr)
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_comm)
          endif

          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_array_alloc)
          endif

          ! copy the messega to the sending buffer for the next cycle
          !
          if( .not. keep_fix_buf ) then
            SAFE_DEALLOCATE(buf_send)
            SAFE_ALLOCATE(buf_send, (nrec_row, nrec_col, kp%nspin))
          end if
          !
          do ispin = 1 , kp%nspin
!$OMP PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_send(jj,ii,ispin) = buf_rec(jj,ii,ispin)
              enddo
            enddo
          end do
          !
          !XXX buf_send = buf_rec
          !
          nsend_row = nrec_row
          nsend_col = nrec_col
          ! deallocate receiving buffer
          if( .not. keep_fix_buf ) then
            SAFE_DEALLOCATE(buf_rec)
          end if
          !
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_array_alloc)
          endif

        end if

    enddo ! ipe
    call progress_free(prog_info)

    ! time for chi omega = 0
    if(peinf%inode .eq. 0) then
      call timing%stop(timing%chi_sum_sub_omega_0)
    endif

    if(pol%nonblocking_cyclic) then
      ! done
      do ispin = 1 , kp%nspin
        !$OMP PARALLEL DO collapse(2)
        do ii = 1, scal%npc
          do jj = 1, scal%npr
            this%chilocal2(jj,ii,ispin) = buf_send(jj,ii,ispin)
          enddo
        enddo
      end do
      !XXXXX this%chilocal2(:,:,:) = buf_send(:,:,:)
      !
      if(keep_fix_buf) then
        ! clean-up
        SAFE_DEALLOCATE(buf_send)
        SAFE_DEALLOCATE(buf_rec)
        SAFE_DEALLOCATE(buf_temp)
        !
        SAFE_DEALLOCATE(tmprowindex)
        SAFE_DEALLOCATE(tmpcolindex)
        SAFE_DEALLOCATE(tmprowph)
        SAFE_DEALLOCATE(tmpcolph)
        !
        SAFE_DEALLOCATE(this%chilocal)
        SAFE_DEALLOCATE(this%gmetempr)
        SAFE_DEALLOCATE(this%gmetempc)
      else
        ! standard case
        SAFE_DEALLOCATE(buf_send)
      end if
    end if

    SAFE_ALLOCATE(chiRDyn_local, (scal%npr,scal%npc))
    chiRDyn_local(:,:) = this%chilocal2(:,:,1)
    if(kp%nspin.eq.2)  chiRDyn_local(:,:) = chiRDyn_local(:,:) + this%chilocal2(:,:,2)

    ! make a copy of chi at omega zero (only in case full chi is needed)
    if(pol%need_full_chi) then
      SAFE_ALLOCATE(chi_omega0_save, (scal%npr,scal%npc,kp%nspin))
      chi_omega0_save(:,:,:) = this%chilocal2(:,:,:)
    end if

    SAFE_DEALLOCATE(this%chilocal2)

    ! Generator Coulomb Interaction Array Vcoul
    SAFE_ALLOCATE(isorti, (gvec%ng))
    SAFE_ALLOCATE(vcoul, (pol%nmtx))

    vcoul(:)=0.0d0
    do i=1,gvec%ng
      isorti(pol%isrtx(i)) = i
    end do

    avgcut=TOL_ZERO
    iscreen=0
    nfk=product(kp%kgrid(1:3))
    q0vec=0d0
    iparallel=1

    if(peinf%inode .eq. 0) then
      call timing%start(timing%chi_sum_sub_vcoul)
    endif


    qgrid(:)=1

    epsheaddummy=0.0d0
    wcoul0=0.0d0

    if(pol%nfreq_group .eq. 1) then
      call vcoul_generator(pol%icutv,pol%truncval,gvec,crys%bdot,crys%celvol, &
        nfk,pol%nmtx,pol%isrtx,iscreen,q0,q0vec,vcoul, &
        pol%iwritecoul,iparallel,avgcut,oneoverq,qgrid,epsheaddummy, &
        work_scell,.false.,wcoul0)
    else
      call vcoul_generator(pol%icutv,pol%truncval,gvec,crys%bdot,crys%celvol, &
        nfk,pol%nmtx,pol%isrtx,iscreen,q0,q0vec,vcoul, &
        pol%iwritecoul,iparallel,avgcut,oneoverq,qgrid,epsheaddummy, &
        work_scell,.false.,wcoul0,nfreq_group=pol%nfreq_group)
    endif

    if(peinf%inode .eq. 0) then
      call timing%stop(timing%chi_sum_sub_vcoul)
    endif


    ! calculate symmetrized polarizability
    do igp_l = 1, scal%npc
      igp_g = indxl2g(igp_l, scal%nbl, scal%mypcol, 0, scal%npcol)
      vc = sqrt(vcoul(igp_g))
      do ig_l = 1, scal%npr
        ig_g = indxl2g(ig_l, scal%nbl, scal%myprow, 0, scal%nprow)
        chiRDyn_local(ig_l,igp_l) = chiRDyn_local(ig_l,igp_l) * sqrt(vcoul(ig_g)) * vc
      enddo
    enddo

    ! Dump write/read matrices if needed
    if( do_dump_matrix ) then
      call dump_matrix(scal, chiRDyn_local, peinf%inode, 100000*iq)
    else
      if( do_read_dump_matrix ) then
        chiRDyn_local = ZERO
        call read_dump_matrix(scal, chiRDyn_local, peinf%inode, 100000*iq)
      end if
    end if

    ! copy the symmetrized chi for omega = 0 only if we work within the subspace
    IF(.NOT.pol%need_full_chi) THEN
      SAFE_ALLOCATE(pol%chiRDyn_sym_omega0, (scal%npr,scal%npc))
      pol%chiRDyn_sym_omega0(:,:) = chiRDyn_local(:,:)
    END IF

    !XXXXXXXXXXXXXXXXx
    ! ALLOCATE(chilocal(pol%nmtx,pol%nmtx))
    ! chilocal = ZERO
    ! icurr=0
    ! do i=1,pol%nmtx
    !   irow=MOD(INT(((i-1)/scal%nbl)+TOL_SMALL),scal%nprow)
    !   if(irow.ne.scal%myprow) cycle
    !   vc = SQRT(vcoul(i))
    !   do j=1,pol%nmtx
    !     icol=MOD(INT(((j-1)/scal%nbl)+TOL_SMALL),scal%npcol)
    !     if(icol .eq. scal%mypcol) then
    !       icurr=icurr+1
    !       irowm=INT((icurr-1)/scal%npc+TOL_SMALL)+1
    !       icolm=MOD((icurr-1),scal%npc)+1

    !       chilocal(i,j) = chiRDyn_local(irowm,icolm)
    !     endif
    !   end do
    ! end do
    ! CALL MPI_ALLREDUCE(MPI_IN_PLACE,chilocal(1,1),pol%nmtx*pol%nmtx,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
    ! IF(peinf%inode.eq.0) THEN
    !   WRITE(3000,*) pol%nmtx
    !   DO j = 1, pol%nmtx
    !     DO i = 1, pol%nmtx
    !       WRITE(3000,*) i, j, chilocal(i,j)
    !     END DO
    !   END DO
    ! END IF
    ! DEALLOCATE(chilocal)
    !XXXXXXXXXXXXXXXXx

    ! diagonalize
    if(peinf%inode .eq. 0) then
      call timing%start(timing%chi_sum_sub_diag)
    endif

    SAFE_ALLOCATE(pol%eigenvect_omega0, (scal%npr,scal%npc))
    SAFE_ALLOCATE(pol%eigenval_omega0, (pol%nmtx))

#ifdef USEELPA
#ifdef CPLX
    IF(pol%use_elpa) THEN
      CALL diagonalize_elpa_chi(pol%nmtx, scal, chiRDyn_local, &
                            pol%chi_eigenvalue_cutoff, pol%neig_sub_input, neig_sub, &
                            pol%eigenvect_omega0, pol%eigenval_omega0)
    ELSE
#endif
#endif
      CALL diagonalize_scalapack(pol%nmtx, scal, chiRDyn_local, &
                                 pol%chi_eigenvalue_cutoff, pol%neig_sub_input, neig_sub, &
                                 pol%eigenvect_omega0, pol%eigenval_omega0)
#ifdef USEELPA
#ifdef CPLX
    END IF
#endif
#endif

    SAFE_DEALLOCATE(chiRDyn_local)
    if(peinf%inode .eq. 0) then
      call timing%stop(timing%chi_sum_sub_diag)
    endif


#ifdef MPI
    call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif

    ! here we have the eigenvector at omega=0, redistribute gme into scalapack format
    ! easy communication scheme
    CALL calc_pol_sub_trunc_easy(this,pol,scal,kp,kpq,vwfn,cwfn,&
                                 nst,nrk,indt,pht,gvec,crys,q0,&
                                 neig_sub, pol%eigenvect_omega0, pol%eigenval_omega0, vcoul)

    ! replace the polarizability at omega zero with the non-truncated chi
    ! (only in case full chi is needed)
    if(pol%need_full_chi) then
      DO ispin = 1, kp%nspin
        pol%chiRDyn(:,:,1,ispin) = chi_omega0_save(:,:,ispin)
      END DO
      SAFE_DEALLOCATE(chi_omega0_save)
      ! in this case we don`t need eigenvector and eigenvalues
      SAFE_DEALLOCATE(pol%eigenvect_omega0)
      SAFE_DEALLOCATE(pol%eigenval_omega0)
    else
      ! this is the case in which we work only within the subspace (except for omega = 0)
      ! leave the possibility to do something
      ! save coulomb potential to be reused in epsinv
      SAFE_ALLOCATE(pol%vcoul_sub, (pol%nmtx))
      pol%vcoul_sub(:) = vcoul(:)
    end if
#ifdef MPI
    call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif

    POP_SUB(chi_summation_sub_trunc)

  end subroutine chi_summation_sub_trunc

#ifdef USEELPA
#ifdef CPLX
  subroutine diagonalize_elpa_chi(nmtx, scal, matrix, chi_eigenvalue_cutoff, neig_sub_input,&
                              neig_sub_out, eigen_vect, eigenval)
  integer, intent(in) :: nmtx
  type (scalapack), intent(in) :: scal
  SCALAR, intent(inout) :: matrix(scal%npr*scal%npc)
  SCALAR, intent(out)   :: eigen_vect(scal%npr*scal%npc)
  real(DP), intent(out) :: eigenval(nmtx)
  real(DP)              :: chi_eigenvalue_cutoff, max_abs_eigen
  INTEGER               :: neig_sub_input, neig_sub_out

  real(DP)              :: scaling
  integer               :: neig_sub, neig_sub_temp, i

  PUSH_SUB(diagonalize_elpa_chi)

  ! if not given in input use 25% of eigenvectors (fixed number, move in input?)
  scaling = 0.25D+00
  neig_sub = MIN(nmtx,INT(nmtx/scaling)+1)
  IF(neig_sub_input > 0) THEN
    neig_sub = MIN(neig_sub_input,nmtx)
  END IF
  IF(peinf%inode .eq. 0) then
    write(6,'(1x,A)')      'Diagonalization with ELPA'
    write(6,'(1x,A,e8.2)') 'EPS Truncation = ', chi_eigenvalue_cutoff
    write(6,'(1x,A,F8.2)') 'Scaling =        ', scaling
    write(6,'(1x,A,I8)')   'Scaled basis =   ', neig_sub
  END IF

  eigenval(:) = 0d0
  CALL diagonalize_elpa(scal, nmtx, neig_sub, matrix, eigenval, eigen_vect)

  ! truncate
  neig_sub_temp = 0
  max_abs_eigen = MAXVAL(ABS(eigenval))
  DO i = 1, neig_sub
    IF(ABS(eigenval(i))/max_abs_eigen >= chi_eigenvalue_cutoff) neig_sub_temp = neig_sub_temp + 1
    ! IF(ABS(eigen(i)) >= chi_eigenvalue_cutoff) neig_sub_temp = neig_sub_temp + 1
    ! IF(peinf%inode .eq. 0) WRITE(2000,*) eigen(i)
  END DO
  neig_sub = neig_sub_temp
  IF(peinf%inode .eq. 0) then
     write(6,*) 'N eigen with relative MAX ABS value >= chi_eigenvalue_cutoff', neig_sub
  END IF
  neig_sub_out = neig_sub

  POP_SUB(diagonalize_elpa_chi)

  end subroutine diagonalize_elpa_chi
#endif
#endif

  subroutine diagonalize_scalapack(nmtx, scal, matrix, chi_eigenvalue_cutoff, neig_sub_input,&
                                   neig_sub_out, eigen_vect, eigenval)
  integer, intent(in) :: nmtx
  type (scalapack), intent(in) :: scal
  SCALAR, intent(inout) :: matrix(scal%npr,scal%npc)
  complex(DPC), intent(out)  :: eigen_vect(scal%npr,scal%npc)
  real(DP), intent(out) :: eigenval(nmtx)
  real(DP)              :: chi_eigenvalue_cutoff
  INTEGER               :: neig_sub_input, neig_sub_out

  integer :: info, desca(9), lwork, lrwork, liwork, ipiv(scal%npr+scal%nbl)
  integer :: descb(9)
  integer, allocatable :: iwork(:)
  SCALAR, allocatable :: work(:)

  SCALAR, allocatable :: work_2(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: eigen(:), rwork(:)
  DOUBLE PRECISION :: vl, vu
  INTEGER :: iii, jjj, jj, i, j, icurr, irow, icol, irowm, icolm
  INTEGER :: neig_sub, neig_sub_temp
  DOUBLE PRECISION :: scaling, max_abs_eigen
  DOUBLE PRECISION :: ABSTOL, ORFAC
  DOUBLE PRECISION, ALLOCATABLE :: gap(:)
  INTEGER, ALLOCATABLE :: iclustr(:), ifail(:)
  INTEGER :: Neigen_found, Neigenvect_found, clustersize, locsize, nbc, nbr
  INTEGER :: nn, nnp, np0, mq0, nq0
  character(len=100) :: tmpstr
  logical :: find_eigen_in_range

  PUSH_SUB(diagonalize_scalapack)

  call descinit(desca, nmtx, nmtx, scal%nbl, scal%nbl, 0, 0, &
                scal%icntxt, max(1,scal%npr), info)
  if(info < 0) then
    write(0,'(a,i3,a)') 'Argument number ', -info, ' had an illegal value on entry.'
    call die("descinit error for descaA in inversion")
  else if(info > 0) then
    write(0,*) 'info = ', info
    call die("descinit error for descaA in inversion")
  endif

  SAFE_ALLOCATE(work_2, (scal%npr,scal%npc))
  work_2 = matrix

  find_eigen_in_range = .TRUE.

  scaling = 1.0D+00
  neig_sub = MIN(nmtx,INT(nmtx/scaling)+1)
  IF(neig_sub_input > 0) THEN
    neig_sub = MIN(neig_sub_input,nmtx)
    find_eigen_in_range = .FALSE.
  END IF
  IF(peinf%inode .eq. 0) then
    write(6,'(1x,A)')      'Diagonalization with Scalapack'
    write(6,'(1x,A,e8.2)') 'EPS Truncation = ', chi_eigenvalue_cutoff
    write(6,'(1x,A,F8.2)') 'Scaling =        ', scaling
    write(6,'(1x,A,I8)')   'Scaled basis =   ', neig_sub
  END IF

  vl = -100.0D+00
  vu =   -ABS(chi_eigenvalue_cutoff)
  ABSTOL = 0.0D+00
  ORFAC  = 5.0D-7
  ! ABSTOL = 1.0D-7
  ! ORFAC  = 5.0D-7
  Neigen_found = 0
  Neigenvect_found = 0

  clustersize = MIN(50,neig_sub)
  ! clustersize = nmtx/MAX(INT(SQRT(DBLE(scal%npcol*scal%nprow))),1)
  nbc = nmtx/(scal%nbl * scal%npcol) + 1
  nbr = nmtx/(scal%nbl * scal%nprow) + 1
  locsize = scal%nbl * scal%nbl * nbc * nbr
  nn  = MAX(nmtx, scal%nbl, 2)
  nnp = MAX(nmtx, scal%nprow*scal%npcol+1, 4)
  np0 = numroc(nn, scal%nbl, 0, 0, scal%nprow)
  mq0 = numroc(MAX(neig_sub,scal%nbl,2), scal%nbl, 0, 0, scal%npcol)

#ifdef CPLX
  nq0    = numroc(nn, scal%nbl, 0, 0, scal%npcol)
  lwork  = nmtx + (np0+nq0+scal%nbl)*scal%nbl
  lrwork = 4*nmtx + MAX(5*nn,np0*mq0) + iceil(neig_sub,scal%nprow*scal%npcol)*nn + &
           (clustersize-1)*nmtx
#else
  lwork  = 5*nmtx + MAX(5*nn,np0*mq0+2*scal%nbl*scal%nbl) + &
           iceil(neig_sub,scal%nprow*scal%npcol)*nn + (clustersize-1)*nmtx
  lrwork = 0
#endif
  liwork = 6*nnp

  SAFE_ALLOCATE(work, (lwork))
#ifdef CPLX
  SAFE_ALLOCATE(rwork, (lrwork))
#endif
  SAFE_ALLOCATE(iwork, (liwork))
  SAFE_ALLOCATE(iclustr, (2*scal%nprow*scal%npcol))
  SAFE_ALLOCATE(gap, (scal%nprow*scal%npcol))
  SAFE_ALLOCATE(ifail, (nmtx))
  SAFE_ALLOCATE(eigen, (nmtx))
  eigen = 0.0D+00

  info = 0
#ifdef CPLX
  IF(find_eigen_in_range) THEN
    CALL pzheevx('V', 'V', 'U', nmtx, work_2, 1, 1, desca, vl, vu, 1, neig_sub, ABSTOL, Neigen_found, Neigenvect_found, &
                 eigen, ORFAC,  matrix, 1, 1, desca, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, &
                 gap, info)
  ELSE
    CALL pzheevx('V', 'I', 'U', nmtx, work_2, 1, 1, desca, vl, vu, 1, neig_sub, ABSTOL, Neigen_found, Neigenvect_found, &
                 eigen, ORFAC,  matrix, 1, 1, desca, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, &
                 gap, info)
  END IF
#else
  CALL pdsyevx('V', 'I', 'U', nmtx, work_2, 1, 1, desca, vl, vu, 1, neig_sub, ABSTOL, Neigen_found, Neigenvect_found, &
               eigen, ORFAC,  matrix, 1, 1, desca, work, lwork, iwork, liwork, ifail, iclustr, &
               gap, info)
#endif

  IF(info/=0) THEN
    IF(info < 0) THEN
      call die(" error in parameters for pzheevx/pdsyevx")
    ELSE
      ! Source: http://www.netlib.org/lapack-dev/Patch/SRC/pzheevx.f
      ! if (MOD(INFO,2).NE.0), then one or more eigenvectors failed to converge.  Their indices are stored
      ! in IFAIL.  Ensure ABSTOL=2.0*PDLAMCH( 'U' ). Send e-mail to scalapack@cs.utk.edu
      if(mod(info, 2) /= 0) then
        if(peinf%inode == 0) then
          write(0,*) " info = ", info
          write(0,'(a)',advance='no') "Unconverged eigenvalues from ifail: "
          do jj = 1, size(ifail)
            if(ifail(jj) == 0) exit
            write(0,'(i8)',advance='no') ifail(jj)
          enddo
        endif
        call die("Convergence problems in pzheevx/pdsyevx.")
      endif

      ! if (MOD(INFO/2,2).NE.0),then eigenvectors corresponding to one or more clusters of eigenvalues could not be
      ! reorthogonalized because of insufficient workspace. The indices of the clusters are stored in the array ICLUSTR.
      if(mod(info / 2, 2) /= 0) then
        if(peinf%inode == 0) then
          write(0,*) " info = ", info
          do jj = 1, size(iclustr)
            if(iclustr(jj) == 0) exit
            write(0,'(i8)',advance='no') iclustr(jj)
          enddo
        endif
        call die("Could not reorthogonalize due to insufficient workspace in pzheevx/pdsyevx.")
      endif

      ! if (MOD(INFO/4,2).NE.0), then space limit prevented PZHEEVX from computing all of the eigenvectors
      ! between VL and VU.  The number of eigenvectors computed is returned in NZ.
      if(mod(info / 4, 2) /= 0) then
        if(peinf%inode == 0) then
          write(0,*) " info = ", info
          write(0,*) " number eigenvalues found = ", Neigenvect_found
        endif
        call die("Space limit prevented computing all eigenvectors in pzheevx/pdsyevx.")
      endif

      ! if (MOD(INFO/8,2).NE.0), then PZSTEBZ failed to compute eigenvalues.  Ensure ABSTOL=2.0*PDLAMCH( 'U' )
      ! Send e-mail to scalapack@cs.utk.edu
      if(mod(info / 8, 2) /= 0) then
        if(peinf%inode == 0) then
          write(0,*) " info = ", info
        endif
        call die("PZSTEBZ failed to compute eigenvalues in pzheevx/pdsyevx.")
      endif
    END IF
  END IF

  SAFE_DEALLOCATE(work)
#ifdef CPLX
  SAFE_DEALLOCATE(rwork)
#endif
  SAFE_DEALLOCATE(iwork)
  SAFE_DEALLOCATE(iclustr)
  SAFE_DEALLOCATE(gap)
  SAFE_DEALLOCATE(ifail)

  SAFE_DEALLOCATE(work_2)

  IF(Neigen_found/=neig_sub) THEN
    neig_sub = Neigen_found
    ! print a warning
    IF(peinf%inode .eq. 0) then
      write(6,*) 'Eigenvalue found = ', neig_sub
    END IF
  END IF
  IF(Neigenvect_found < neig_sub) THEN
    write(tmpstr,'(a, i10, a, i10, a)') 'Diagonalization with pzheevx/pdsyevx failed: only ', &
       Neigenvect_found, ' of ', neig_sub, ' eigenvectors found.'
    call die(tmpstr)
  END IF

  ! truncate
  neig_sub_temp = 0
  max_abs_eigen = MAXVAL(ABS(eigen))
  DO i = 1, neig_sub
    IF(ABS(eigen(i))/max_abs_eigen >= chi_eigenvalue_cutoff) neig_sub_temp = neig_sub_temp + 1
    ! IF(ABS(eigen(i)) >= chi_eigenvalue_cutoff) neig_sub_temp = neig_sub_temp + 1
    ! IF(peinf%inode .eq. 0) WRITE(2000,*) eigen(i)
  END DO
  neig_sub = neig_sub_temp

  IF(peinf%inode .eq. 0) then
     write(6,*) 'N eigen with relative MAX ABS value >= chi_eigenvalue_cutoff', neig_sub
  END IF
  neig_sub_out = neig_sub

  ! matrix at this point contains the eigenvectors contain here the e
  ! make a copy into eigenvect (eigenvectors are defined complex)
  eigen_vect = (0D0,0D0)
  eigen_vect = matrix
  eigenval = eigen

  SAFE_DEALLOCATE(eigen)

  POP_SUB(diagonalize_scalapack)

  end subroutine diagonalize_scalapack

  ! This is going to be the same communication scheme as comm_matrix
  ! for which the eigenvector of the symmetrized polarizability are
  ! replicated on each MPI task
  subroutine calc_pol_sub_trunc_easy(this,pol,scal,kp,kpq,vwfn,cwfn,&
                                     nst,nrk,indt,pht,gvec,crys,q0,&
                                     neig_sub, eigen_vect, eigenval, vcoul)
    type(chi_summator_t), intent(INOUT) :: this
    type(polarizability), intent(INOUT) :: pol
    type(scalapack), intent(in) :: scal
    type(kpoints), intent(IN) :: kp, kpq
    type(valence_wfns), intent(IN) :: vwfn
    type(conduction_wfns), intent(IN) :: cwfn
    INTEGER :: neig_sub
    complex(DPC), allocatable :: eigen_vect(:,:)
    real(DP), allocatable :: eigenval(:)
    real(DP), allocatable :: vcoul(:)

    integer, intent(IN) :: nst(:)
    integer, intent(IN) :: nrk
    integer, intent(INOUT) :: indt(:,:,:)
    SCALAR,  intent(INOUT) :: pht(:,:,:)
    type (gspace), intent(in) :: gvec
    type (crystal), intent(in) :: crys
    real(DP), intent(in) :: q0(3)
    real(DP) :: zvalue, cv_energy, occ_diff
    complex(DPC), allocatable :: edenDRtemp(:)

    INTEGER :: ntot, ntot2
    INTEGER :: irk, ipe, ipe_row, ipe_col, ipe_real, ispin, itot, it, iv, ic, i, j, iv_local
    INTEGER :: irow, icol, irowm, icolm, icurr, ilimit, freq_idx, &
               i_myband, mytot

    integer, allocatable :: tmprowindex(:),tmpcolindex(:)
    SCALAR, allocatable :: tmprowph(:),tmpcolph(:)
    type(cvpair_info) :: cvpair_temp
    complex(DPC) :: negfact

    SCALAR :: epsheaddummy, wcoul0
    type (twork_scell) :: work_scell
    integer, allocatable :: isorti(:)
    real(DP) :: vc, oneoverq, avgcut
    integer :: iscreen, nfk, iparallel
    integer :: qgrid(3)
    real(DP) :: q0vec(3)

    type(progress_info) :: prog_info

    INTEGER :: n_vc_tot, n_v_tot, n_c_tot, nmtx, kk, jj, ii, info
    INTEGER :: my_sub_nlocal_col, my_sub_nlocal_row, itot_rk2, &
      irow_global, icol_global

    complex(DPC), allocatable :: C_aux(:,:), gme_sub(:,:,:,:,:), gme_temp(:,:,:)
    complex(DPC), allocatable :: C_Pgemm(:,:)
    integer, allocatable :: grid2D_2inode(:,:), row_col_loc(:,:), inode_2grid2D(:,:)
    INTEGER :: pe_col_size, pe_row_size, desc_sub(9), desca(9)
    INTEGER, ALLOCATABLE :: row_indx(:), col_indx(:)
    ! variables for block transformation
    integer :: iproc_dum
    integer :: nmpinode
    integer :: iblock, block_size_max, ib_start, ib_end, ib_size, Nblocks
    real(DP) :: mem, avail_mem, mem_offset, mem_single_block
    ! variables for non-blocking cyclic scheme
    integer :: isend_static, irec_static
    integer :: actual_send, actual_rec
    integer :: nsend_row, nsend_col, nrec_row, nrec_col
    integer :: req_send, tag_send, req_rec, tag_rec
    integer :: stat(MPI_STATUS_SIZE)
    complex(DPC), allocatable :: buf_send(:,:,:,:), buf_rec(:,:,:,:), buf_temp(:,:,:,:)
    integer :: size_N, size_M, size_K
    ! variables for fast eigenvectors communication
    logical :: fast_eig_comm
    complex(DPC), allocatable :: C_elem(:), C_elem_send(:), C_elem_rec(:)
    real(DP), allocatable :: sqrt_vcoul_local(:)
    integer, allocatable :: elem_per_block(:)
    integer, allocatable :: glob_C_indec(:,:)
    integer, allocatable :: glob_C_indec_send(:,:), glob_C_indec_rec(:,:)
    integer :: num_loc_col_within_Neig, max_num_elem, iii, jjj
    integer :: req_send_info, tag_send_info, req_rec_info, tag_rec_info

    PUSH_SUB(calc_pol_sub_trunc_easy)

    fast_eig_comm = .true.
    if( pol%sub_collective_eigen_redistr ) then
      ! slow redistribution scheme, for debug only
      fast_eig_comm = .false.
    end if

    ! precompute the dimensions of the different objects
    ntot=0
    ntot2=0
    ntot = peinf%nvownactual*peinf%ncownactual
    do irk = 1, nrk
      ntot2=ntot2 + nst(irk)
    enddo
    ntot=ntot*ntot2

    ! WRITE(*,*) ntot2, nrk, maxval(nst), minval(nst)

    n_c_tot = cwfn%nband
    n_v_tot = vwfn%nband+pol%ncrit
    n_vc_tot = n_c_tot * n_v_tot
    nmtx = pol%nmtx

    ! first replicate the eigenvector of the symmetrized polarizability at omega=0
    ! on all processes (do it by blocks)
    ! calculate block size according to available memory
    call procmem(mem,nmpinode)
    avail_mem = mem
    mem_offset = dble(neig_sub) * dble(peinf%ncownmax) * dble(peinf%nvownmax) * &
                 dble(ntot2) * dble(kp%nspin) * 16D+00
    mem_single_block = MAX(dble(neig_sub)*2.0D+00, dble(neig_sub + peinf%ncownmax*peinf%nvownmax)) * 16D+00
    !XXXXX  factor 2 (probably the allreduce inplace allocate internaly an additional buffer)
    mem_single_block = mem_single_block * 2.0D+00
    !XXXXX
    block_size_max = int((avail_mem - mem_offset) / mem_single_block)
    block_size_max = MAX(block_size_max,1)
    block_size_max = MIN(block_size_max,nmtx)
    if(peinf%inode .eq. 0) then
      write(6,*)
      write(6,'(1x,"Basis Transformation info:")')
      write(6,'(1x,"Memory available: ",f0.1," MB per PE")') avail_mem / 1024**2
      write(6,'(1x,"Memory offset: ",f0.1," MB")') mem_offset / 1024**2
      write(6,'(1x,"Memory for a sigle block: ",f0.1," kB")') mem_single_block / 1024
      write(6,'(1x,"Max Block Size: ",i0)') block_size_max
    end if
    ! broadcast the result to make sure all processes have same value
#ifdef MPI
    call MPI_Bcast(block_size_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
    Nblocks = (nmtx-1) / block_size_max + 1
    SAFE_ALLOCATE(gme_sub, (neig_sub,peinf%ncownmax,peinf%nvownmax,ntot2,kp%nspin))
    gme_sub = ZERO

    call accel_enter_data_map_to(gme_sub, chi_summation_algo)

    ! initialize some stuff for fast_eig_comm
    if ( fast_eig_comm ) then
      ! calculate how many local elements have a global index smaller than Neig
      num_loc_col_within_Neig = 0
      do icolm = 1, scal%npc
        j = INDXL2G (icolm, scal%nbl, scal%mypcol, 0, scal%npcol)
        if ( j <= neig_sub ) then
          num_loc_col_within_Neig = num_loc_col_within_Neig + 1
        end if
      end do
      ! figure the biggest number of elements will be exchanged and the actual
      ! number of elements for each block
      SAFE_ALLOCATE(elem_per_block, (Nblocks))
      elem_per_block = 0
      iii = 0
      do iblock = 1, nmtx, block_size_max
        ib_start = iblock
        ib_end = ib_start + block_size_max - 1
        ib_end = MIN(ib_end, nmtx)
        ib_size = ib_end - ib_start + 1
        ! increase block index
        iii = iii + 1
        ! compute how many rows are in range and [ib_start,ib_end] and
        ! accumulate relative number of cols
        do i = ib_start, ib_end
          irow = INDXG2P( i, scal%nbl, iproc_dum, 0, scal%nprow)
          if ( irow .ne. scal%myprow) cycle
          elem_per_block(iii) = elem_per_block(iii) + num_loc_col_within_Neig
        end do
      end do
      max_num_elem = MAXVAL(elem_per_block(:))
      call mpi_allreduce(MPI_IN_PLACE, max_num_elem, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
      ! allocate communication buffers
      SAFE_ALLOCATE(C_elem, (max_num_elem))
      SAFE_ALLOCATE(C_elem_send, (max_num_elem))
      SAFE_ALLOCATE(C_elem_rec, (max_num_elem))
      SAFE_ALLOCATE(glob_C_indec, (2,max_num_elem+1))
      SAFE_ALLOCATE(glob_C_indec_send, (2,max_num_elem+1))
      SAFE_ALLOCATE(glob_C_indec_rec, (2,max_num_elem+1))
      ! precoumpute SQRT of coulomb potential
      SAFE_ALLOCATE(sqrt_vcoul_local, (scal%npr))
      do irowm = 1, scal%npr
        i = INDXL2G (irowm, scal%nbl, scal%myprow, 0, scal%nprow)
        sqrt_vcoul_local(irowm) = SQRT(vcoul(i))
      end do
      ! get process for communication
      isend_static = MOD(peinf%inode + 1 + peinf%npes, peinf%npes)
      irec_static = MOD(peinf%inode - 1 + peinf%npes, peinf%npes)
    end if

    call progress_init(prog_info, 'transforming basis', 'block', Nblocks)
    icurr=0
    iii = 0
    do iblock = 1, nmtx, block_size_max
      call progress_step(prog_info)
      ! actual block indeces
      ib_start = iblock
      ib_end = ib_start + block_size_max - 1
      ib_end = MIN(ib_end, nmtx)
      ib_size = ib_end - ib_start + 1
      ! allocate C_aux
      SAFE_ALLOCATE(C_aux, (ib_size,neig_sub))
      C_aux = ZERO

      if ( fast_eig_comm ) then
        ! increase block index
        iii = iii + 1
        ! copy information into communication buffers
        jjj = 0
        do i = ib_start, ib_end
          irow = INDXG2P( i, scal%nbl, iproc_dum, 0, scal%nprow)
          if ( irow .ne. scal%myprow) cycle
          irowm = INDXG2L(i, scal%nbl, scal%myprow, 0, scal%nprow)
          vc = sqrt_vcoul_local(irowm)
          do icolm = 1, scal%npc
            j = INDXL2G (icolm, scal%nbl, scal%mypcol, 0, scal%npcol)
            if ( j <= neig_sub ) then
              jjj = jjj + 1
              C_elem(jjj) = eigen_vect(irowm,icolm) * vc
              glob_C_indec(1,jjj) = i
              glob_C_indec(2,jjj) = j
            end if
          end do
        end do
        glob_C_indec(1,max_num_elem+1) = elem_per_block(iii)

        !XXX if ( jjj .ne. elem_per_block(iii) ) STOP

        ! loop over process and start ring communication
        glob_C_indec_send(:,:) = glob_C_indec(:,:)
        C_elem_send(:) = C_elem(:)

        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_sub_eigvet_comm)
        endif
    
        do ipe = 1, peinf%npes
          !
          tag_rec_info = 1
          tag_send_info = 1
          req_rec_info = MPI_REQUEST_NULL
          req_send_info = MPI_REQUEST_NULL
          CALL MPI_Irecv(glob_C_indec_rec, 2*(max_num_elem+1), &
                         MPI_INTEGER, irec_static, &
                         tag_rec_info, MPI_COMM_WORLD, req_rec_info, mpierr)
          CALL MPI_Isend(glob_C_indec_send, 2*(max_num_elem+1), &
                         MPI_INTEGER, isend_static, &
                         tag_send_info, MPI_COMM_WORLD, req_send_info, mpierr)
          !
          tag_rec = 1
          tag_send = 1
          req_rec = MPI_REQUEST_NULL
          req_send = MPI_REQUEST_NULL
          CALL MPI_Irecv(C_elem_rec, max_num_elem, &
                         MPI_COMPLEX_DPC, irec_static, &
                         tag_rec, MPI_COMM_WORLD, req_rec, mpierr)
          CALL MPI_Isend(C_elem_send, max_num_elem, &
                         MPI_COMPLEX_DPC, isend_static, &
                         tag_send, MPI_COMM_WORLD, req_send, mpierr)


          do jjj = 1, glob_C_indec(1,max_num_elem+1)
            i = glob_C_indec(1,jjj)
            j = glob_C_indec(2,jjj)
            C_aux(i-ib_start+1,j) = C_elem(jjj)
          end do

          ! make sure the buffer is received/send
          CALL MPI_Wait(req_rec_info,stat,mpierr)
          CALL MPI_Wait(req_send_info,stat,mpierr)
          !
          glob_C_indec_send(:,:) = glob_C_indec_rec(:,:)
          glob_C_indec(:,:)      = glob_C_indec_rec(:,:)
          !
          CALL MPI_Wait(req_rec,stat,mpierr)
          CALL MPI_Wait(req_send,stat,mpierr)
          !
          C_elem_send(:) = C_elem_rec(:)
          C_elem(:)      = C_elem_rec(:)
          !
        end do
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_sub_eigvet_comm)
        endif
  
      else ! fast_eig_comm
        !
        if ( .true. ) then
          !
          do i = ib_start, ib_end
            irow = INDXG2P( i, scal%nbl, iproc_dum, 0, scal%nprow)
            if ( irow .ne. scal%myprow) cycle
            vc = SQRT(vcoul(i))
            irowm = INDXG2L(i, scal%nbl, scal%myprow, 0, scal%nprow)
            !
            do icolm = 1, scal%npc
              j = INDXL2G (icolm, scal%nbl, scal%mypcol, 0, scal%npcol)
              if ( j <= neig_sub ) then
                C_aux(i-ib_start+1,j) = eigen_vect(irowm,icolm) * vc
              end if
            end do
            !
          end do
          !
        else
          !
          ! old and slow
          do i = ib_start, ib_end
            irow=MOD(INT(((i-1)/scal%nbl)+TOL_SMALL),scal%nprow)
            if(irow.ne.scal%myprow) cycle
            vc = SQRT(vcoul(i))
            do j=1, nmtx
              icol=MOD(INT(((j-1)/scal%nbl)+TOL_SMALL),scal%npcol)
              if(icol .eq. scal%mypcol) then
                icurr=icurr+1
                irowm=INT((icurr-1)/scal%npc+TOL_SMALL)+1
                icolm=MOD((icurr-1),scal%npc)+1

                IF(j<=neig_sub) THEN
                  C_aux(i-ib_start+1,j) = eigen_vect(irowm,icolm) * vc
                END IF

              endif
            end do
          end do
          !
        end if
        ! WRITE(*,*) nmtx,neig_sub
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_sub_eigvet_comm)
        endif

#ifdef MPI
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,C_aux(1,1),ib_size*neig_sub,MPI_COMPLEX_DPC,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_sub_eigvet_comm)
        endif

        !
      end if ! fast_eig_comm

      ! transform basis
      if(peinf%inode .eq. 0) then
        call timing%start(timing%chi_sum_sub_transf)
      endif

      SAFE_ALLOCATE(gme_temp, (ib_size,peinf%ncownmax,peinf%nvownmax))
      gme_temp = ZERO
      call accel_enter_data_map_to(C_aux, chi_summation_algo)
      call accel_enter_data_map_alloc(gme_temp, chi_summation_algo)
      do ispin = 1 , kp%nspin
        itot_rk2 = 0
        do irk = 1, nrk
          do it = 1, nst(irk)
            itot_rk2 = itot_rk2 + 1

            gme_temp = ZERO
            !$OMP PARALLEL DO
            do i = ib_start, ib_end
              gme_temp(i-ib_start+1,1:peinf%ncownactual,1:peinf%nvownactual) = &
                      pol%gme(indt(i,it,irk),1:peinf%ncownactual,1:peinf%nvownactual,ispin,irk,1) * &
                      pht(i,it,irk)
            end do
            !$OMP END PARALLEL DO
            call accel_update_to(gme_temp, chi_summation_algo)

            do iv = 1, vwfn%nband+pol%ncrit
              ! check if I own this partuicular band
              if (peinf%doiownv(iv)) then
                ilimit = peinf%ncownactual
              else
                cycle
              endif
              iv_local = peinf%indexv(iv)

              call accel_zgemm('c', 'n', &
                             neig_sub, peinf%ncownactual, ib_size, &
                             (1D0,0D0), &
                             C_aux, ib_size, &
                             gme_temp(:, :, peinf%indexv(iv)), ib_size, &
                             (1D0,0D0), &
                             gme_sub(:, :, peinf%indexv(iv), itot_rk2, ispin), neig_sub, &
                             chi_summation_algo)
            end do
          end do
        end do
      end do
      call accel_exit_data_map_delete(gme_temp, chi_summation_algo)
      call accel_exit_data_map_delete(C_aux, chi_summation_algo)
      SAFE_DEALLOCATE(gme_temp)
      SAFE_DEALLOCATE(C_aux)
      if(peinf%inode .eq. 0) then
        call timing%stop(timing%chi_sum_sub_transf)
      endif

    end do !iblock
    call progress_free(prog_info)
    SAFE_DEALLOCATE_P(pol%gme)
    ! this is going to be deallocated outside
    SAFE_ALLOCATE(pol%gme, (1,1,1,1,1,1))

    call accel_exit_data_map_from(gme_sub, chi_summation_algo)

    if ( fast_eig_comm ) then
      SAFE_DEALLOCATE(elem_per_block)
      SAFE_DEALLOCATE(C_elem)
      SAFE_DEALLOCATE(C_elem_send)
      SAFE_DEALLOCATE(C_elem_rec)
      SAFE_DEALLOCATE(glob_C_indec)
      SAFE_DEALLOCATE(glob_C_indec_send)
      SAFE_DEALLOCATE(glob_C_indec_rec)
      SAFE_DEALLOCATE(sqrt_vcoul_local)
    end if


    !XXXXX
    ! initialize the descriptor for the subspace matrix
    my_sub_nlocal_row = numroc(neig_sub, scal%nbl, scal%myprow, 0, scal%nprow)
    my_sub_nlocal_col = numroc(neig_sub, scal%nbl, scal%mypcol, 0, scal%npcol)
    call descinit(desc_sub, neig_sub, neig_sub, scal%nbl, scal%nbl, 0, 0, &
                  scal%icntxt, max(1,my_sub_nlocal_row), info)
    if(info < 0) then
      write(0,'(a,i3,a)') 'Argument number ', -info, ' had an illegal value on entry.'
      call die("descinit error for desca sub in sub_chi_sum")
    end if
    ! exchange some info
    SAFE_ALLOCATE(row_col_loc, (2,peinf%npes))
    SAFE_ALLOCATE(grid2D_2inode, (scal%nprow,scal%npcol))
    SAFE_ALLOCATE(inode_2grid2D, (2,peinf%npes))
    row_col_loc   = 0
    grid2D_2inode = 0
    inode_2grid2D = 0
    row_col_loc(1,peinf%inode+1) = my_sub_nlocal_row
    row_col_loc(2,peinf%inode+1) = my_sub_nlocal_col
    grid2D_2inode(scal%myprow+1,scal%mypcol+1) = peinf%inode
    inode_2grid2D(1,peinf%inode+1) = scal%myprow
    inode_2grid2D(2,peinf%inode+1) = scal%mypcol
#ifdef MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,row_col_loc(1,1),2*peinf%npes,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,grid2D_2inode(1,1),scal%nprow*scal%npcol,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,inode_2grid2D(1,1),2*peinf%npes,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
    !XXXXX

    freq_idx = pol%nfreq_group
    IF(pol%nfreq_group .gt. 1) THEN
      freq_idx = peinf%rank_mtxel+1
    END IF

    ! avoid allocate zero size array (for mpi reduce?)
    SAFE_ALLOCATE(this%chilocal2RDyn, (MAX(1,my_sub_nlocal_row),MAX(1,my_sub_nlocal_col),pol%nfreq_in_group,kp%nspin))

    negfact = -1D0*this%fact

    ! setup non-blocking cyclic communication (directly keep fix allocated buffers)
    if(pol%nonblocking_cyclic) then
      ! initialize buffer for non-blocking cyclic scheme
      ! static process for the communication
      isend_static = MOD(peinf%inode + 1 + peinf%npes, peinf%npes)
      irec_static = MOD(peinf%inode - 1 + peinf%npes, peinf%npes)
      ! get sizes
      size_N = my_sub_nlocal_row
      size_M = my_sub_nlocal_col
      size_K = ntot
      call mpi_allreduce(MPI_IN_PLACE, size_N, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
      call mpi_allreduce(MPI_IN_PLACE, size_M, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
      call mpi_allreduce(MPI_IN_PLACE, size_K, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
      ! allocate only once
      SAFE_ALLOCATE(buf_send, (size_N, size_M, pol%nfreq_in_group, kp%nspin))
      SAFE_ALLOCATE(buf_rec,  (size_N, size_M, pol%nfreq_in_group, kp%nspin))
      SAFE_ALLOCATE(buf_temp, (size_N, size_M, pol%nfreq_in_group, kp%nspin))
      !
      do ispin = 1 , kp%nspin
        do kk = 1, pol%nfreq_in_group
          !
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, size_M
            do jj = 1, size_N
              buf_send(jj,ii,kk,ispin) = (0.0D+00,0.0D+00)
            enddo
          enddo
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, size_M
            do jj = 1, size_N
              buf_rec(jj,ii,kk,ispin)  = (0.0D+00,0.0D+00)
            end do
          end do
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, size_M
            do jj = 1, size_N
              buf_temp(jj,ii,kk,ispin) = (0.0D+00,0.0D+00)
            end do
          end do
          !
        end do
      end do
      !
      ! allocate for matrix multiplication
      SAFE_ALLOCATE(this%chilocalRDyn,  (size_N, size_M, pol%nfreq_in_group))  ! (pe_row_size,pe_col_size,pol%nfreq_in_group))
      SAFE_ALLOCATE(this%gmeRDyntempr2, (size_N, size_K, pol%nfreq_in_group))  ! (pe_row_size,ntot,pol%nfreq_in_group))
      SAFE_ALLOCATE(this%gmeRDyntempc,  (size_K, size_M))  ! (ntot,pe_col_size))
      !
      do kk = 1, pol%nfreq_in_group
        !$OMP PARALLEL DO collapse(2)
        do ii = 1, size_M
          do jj = 1, size_N
            this%chilocalRDyn(jj,ii,kk) = (0.0D+00,0.0D+00)
          end do
        end do
        !$OMP PARALLEL DO collapse(2)
        do ii = 1, size_K
          do jj = 1, size_N
            this%gmeRDyntempr2(jj,ii,kk) = (0.0D+00,0.0D+00)
          end do
        end do
      end do
      !$OMP PARALLEL DO collapse(2)
      do ii = 1, size_M
        do jj = 1, size_K
          this%gmeRDyntempc(jj,ii) = (0.0D+00,0.0D+00)
        end do
      end do
      !
    end if


    ! total time for omega neq zero
    if(peinf%inode .eq. 0) then
      call timing%start(timing%chi_sum_sub_omega_neq_0)
    endif


    call progress_init(prog_info, 'building polarizability matrix within the subspace of omega=0', 'processor', &
                       peinf%npes)
    ! go with the big loop over processes calculations
    DO ipe = 1, peinf%npes
      call progress_step(prog_info)
      !
      ! in the standard case there is no difference between ipe_real and ipe
      ! this only matters for pol%nonblocking_cyclic
      ipe_real = ipe
      !
      if(pol%nonblocking_cyclic) then
        ! calculate the actual process we have to send and we are receiving
        actual_send = MOD(peinf%inode + ipe + peinf%npes, peinf%npes)
        actual_rec  = MOD(peinf%inode - ipe + peinf%npes, peinf%npes)
        ! fixed size
        nsend_row = size_N
        nrec_row  = size_N
        nsend_col = size_M
        nrec_col  = size_M
        !
        do ispin = 1 , kp%nspin
          do kk = 1, pol%nfreq_in_group
            !$OMP PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_rec(jj,ii,kk,ispin) = ZERO
              enddo
            enddo
          end do
        end do
        !
        ! post message
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_comm)
        endif

        tag_rec = 1
        tag_send = 1
        req_rec = MPI_REQUEST_NULL
        req_send = MPI_REQUEST_NULL
        CALL MPI_Irecv(buf_rec, nrec_row*nrec_col*pol%nfreq_in_group*kp%nspin, &
                       MPI_COMPLEX_DPC, irec_static, &
                       tag_rec, MPI_COMM_WORLD, req_rec, mpierr)
        CALL MPI_Isend(buf_send, nsend_row*nsend_col*pol%nfreq_in_group*kp%nspin, &
                       MPI_COMPLEX_DPC, isend_static, &
                       tag_send, MPI_COMM_WORLD, req_send, mpierr)
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_comm)
        endif

        ! make sure to put buffer to zero
        do kk = 1, pol%nfreq_in_group
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, size_M
            do jj = 1, size_N
              this%chilocalRDyn(jj,ii,kk) = (0.0D+00,0.0D+00)
            end do
          end do
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, size_K
            do jj = 1, size_N
              this%gmeRDyntempr2(jj,ii,kk) = (0.0D+00,0.0D+00)
            end do
          end do
        end do
        !$OMP PARALLEL DO collapse(2)
        do ii = 1, size_M
          do jj = 1, size_K
            this%gmeRDyntempc(jj,ii) = (0.0D+00,0.0D+00)
          end do
        end do
        !
        ! in this case the ipe is the real process from which we should receive
        ipe_real = actual_rec+1
      end if  ! pol%nonblocking_cyclic
      !

      !
      pe_row_size = row_col_loc(1,ipe_real)
      pe_col_size = row_col_loc(2,ipe_real)
      IF(pe_row_size==0 .OR. pe_col_size==0) CYCLE
      ! precompute indices
      ! row
      SAFE_ALLOCATE(row_indx, (pe_row_size))
      ipe_row = inode_2grid2D(1,ipe_real)
      row_indx = 0
      DO irow = 1, pe_row_size
        row_indx(irow) = INDXL2G( irow, scal%nbl, ipe_row, 0, scal%nprow)
      END DO
      ! col
      SAFE_ALLOCATE(col_indx, (pe_col_size))
      ipe_col = inode_2grid2D(2,ipe_real)
      col_indx = 0
      DO icol = 1, pe_col_size
        col_indx(icol) = INDXL2G( icol, scal%nbl, ipe_col, 0, scal%npcol)
      END DO

      ! IF(peinf%inode==0) WRITE(*,*) ipe, ipe_real-1
      ! IF(peinf%inode==0) WRITE(*,*) row_indx
      ! IF(peinf%inode==0) WRITE(*,*) col_indx
      ! IF(peinf%inode==0) WRITE(*,*)

      do ispin = 1, kp%nspin
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_prep)
        endif

        if ( .not. pol%nonblocking_cyclic ) then
          ! standard calculation
          SAFE_ALLOCATE(this%chilocalRDyn, (pe_row_size,pe_col_size,pol%nfreq_in_group))
          SAFE_ALLOCATE(this%gmeRDyntempr2, (pe_row_size,ntot,pol%nfreq_in_group))
          SAFE_ALLOCATE(this%gmeRDyntempc, (ntot,pe_col_size))
          this%chilocalRDyn  = ZERO
          this%gmeRDyntempr2 = ZERO
          this%gmeRDyntempc  = ZERO
          !
        end if

        itot = 0
        itot_rk2 = 0
        do irk = 1, nrk
          do it = 1, nst(irk)
            itot_rk2 = itot_rk2 + 1
            do iv = 1,vwfn%nband+pol%ncrit
              if (peinf%doiownv(iv)) then
                ilimit = peinf%ncownactual
              else
                ilimit = 0
              endif

              !$OMP PARALLEL private (mytot, zvalue, edenDRtemp, jj, icurr, &
              !$OMP                   cvpair_temp, cv_energy, irow, irow_global, icol, icol_global)
              SAFE_ALLOCATE(edenDRtemp, (pol%nfreq_in_group))
              !$OMP DO
              do i_myband = 1, ilimit
                 mytot = itot + i_myband
! DVF partial occupations: Also, occ_diff factor added here. Same comment about
! OMP directive.
                 zvalue = pol%edenDyn(peinf%indexv(iv),i_myband,ispin,irk,freq_idx)
                 occ_diff = pol%occ_diff(peinf%indexv(iv),i_myband,ispin,irk,freq_idx)
                 if(pol%lin_denominator<TOL_Zero) then
                   ! this is when the lin_denominator mode is not active.

                   do jj=1, pol%nfreq_in_group
                     if (abs(zvalue) .gt. Tol_Zero) then
                       if (pol%tda) then
                         edenDRtemp(jj)= -0.5d0*( &
                           1d0/(zvalue+(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd))
                       else
                         edenDRtemp(jj)= -0.5d0*occ_diff*( &
                           1d0/(zvalue-(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd)+ &
                           1d0/(zvalue+(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd))
                       endif
                     else
                       edenDRtemp(jj)= 0D0
                     endif
                   enddo

                 else
                   !this is when lin_denominator mode is active

                   cvpair_temp = pol%lin_edenDyn(peinf%indexv(iv),i_myband,ispin,irk,freq_idx)
                   cv_energy = cvpair_temp%ev - cvpair_temp%ec
                   do jj=1, pol%nfreq_in_group
                     !check if denominator is not small. If so, do the usual thing

                     if(abs(cv_energy+pol%dFreqGrid(jj)/ryd)-pol%lin_denominator>-TOL_Zero .and. &
                        abs(cv_energy-pol%dFreqGrid(jj)/ryd)-pol%lin_denominator>-TOL_Zero) then

                       if (abs(zvalue) .gt. Tol_Zero) then
                         if (pol%tda) then
                           edenDRtemp(jj)= -0.5d0*( &
                             1d0/(zvalue+(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd))
                         else
                           edenDRtemp(jj)= -0.5d0*( &
                             1d0/(zvalue-(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd)+ &
                             1d0/(zvalue+(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd))
                         endif
                         else
                         edenDRtemp(jj)= 0D0
                       endif

                     else ! if denominator is small, do linearized energy denominator
                       if(cvpair_temp%vltc) then
                         if (pol%tda) then
                           edenDRtemp(jj) = -0.5d0*(&
                             
                             integrate_mbz(kp,kpq%rk(1:3,cvpair_temp%idx_kp),&
                             cvpair_temp%ev,cvpair_temp%ec,&
                             cvpair_temp%vv,cvpair_temp%vc,&
                             pol%dFreqGrid(jj)/ryd,pol%efermi/ryd))
                         else
                           edenDRtemp(jj) = -0.5d0*(&
                             
                             conjg(integrate_mbz(kp,kpq%rk(1:3,cvpair_temp%idx_kp),&
                             cvpair_temp%ev,cvpair_temp%ec,&
                        
                             cvpair_temp%vv,cvpair_temp%vc,&
                             -pol%dFreqGrid(jj)/ryd,pol%efermi/ryd))&
                             
                             +integrate_mbz(kp,kpq%rk(1:3,cvpair_temp%idx_kp),&
                             cvpair_temp%ev,cvpair_temp%ec,&
                             cvpair_temp%vv,cvpair_temp%vc,&
                             pol%dFreqGrid(jj)/ryd,pol%efermi/ryd))
                         endif
                       else
                         edenDRtemp(jj)= 0D0
                       endif
                     endif
                   end do

                 end if

                 ! copy matrix elements
                 DO irow = 1, pe_row_size
                   irow_global = row_indx(irow)
                   DO jj = 1, pol%nfreq_in_group
                     this%gmeRDyntempr2(irow,mytot,jj) = gme_sub(irow_global,i_myband,peinf%indexv(iv),itot_rk2,ispin)*&
                                                         edenDRtemp(jj)
                   END DO
                 END DO
                 DO icol = 1, pe_col_size
                   icol_global = col_indx(icol)
                   this%gmeRDyntempc(mytot,icol) = MYCONJG(gme_sub(icol_global,i_myband,peinf%indexv(iv),itot_rk2,ispin))
                 END DO

              end do ! i_band
              !$OMP END DO
              SAFE_DEALLOCATE(edenDRtemp)
              !$OMP END PARALLEL
              itot = itot+ilimit
            end do ! iv
          end do ! it
        end do ! irk

        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_prep)
        endif
  

        if( pol%nonblocking_cyclic ) then
          !
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_gemm)
          endif
          do jj = 1, pol%nfreq_in_group
            call zgemm('n','n', size_N, size_M, size_K, &
                negfact,this%gmeRDyntempr2(:,:,jj), size_N, &
                        this%gmeRDyntempc(:,:), size_K,&
                (0D0,0D0),this%chilocalRDyn(:,:,jj), size_N)
          end do
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_gemm)
          endif
          !
        else ! pol%nonblocking_cyclic
          !
          ! go with zgemm (standard calculation)
          if(ntot > 0) then
            call accel_enter_data_map_to(this%gmeRDyntempr2, chi_summation_algo)
            call accel_enter_data_map_to(this%gmeRDyntempc, chi_summation_algo)
            call accel_enter_data_map_alloc(this%chilocalRDyn, chi_summation_algo)
            do jj =1, pol%nfreq_in_group
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_gemm)
              endif
              call accel_zgemm('n', 'n', &
                             pe_row_size, pe_col_size, ntot, &
                             negfact, &
                             this%gmeRDyntempr2(:,:,jj), pe_row_size, &
                             this%gmeRDyntempc(:,:), ntot, &
                             (0D0,0D0), &
                             this%chilocalRDyn(:,:,jj), pe_row_size, &
                             chi_summation_algo)
              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_gemm)
              endif
            enddo
            call accel_exit_data_map_from(this%chilocalRDyn, chi_summation_algo)
            call accel_exit_data_map_delete(this%gmeRDyntempc, chi_summation_algo)
            call accel_exit_data_map_delete(this%gmeRDyntempr2, chi_summation_algo)
          endif
          !
        end if ! pol%nonblocking_cyclic

        if( pol%nonblocking_cyclic ) then
          !
          ! just copy the result
          do kk = 1, pol%nfreq_in_group
!$OMP PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_temp(jj,ii,kk,ispin) = this%chilocalRDyn(jj,ii,kk)
              enddo
            enddo
          end do
          !
        else ! pol%nonblocking_cyclic
        !
        ! standard calculation
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_comm)
        endif

#ifdef MPI
        ! all reduce
        call MPI_Reduce(this%chilocalRDyn(1,1,1),this%chilocal2RDyn(1,1,1,ispin), &
                        pol%nfreq_in_group*pe_row_size*pe_col_size,MPI_COMPLEX_DPC,&
                        MPI_SUM,ipe_real-1,MPI_COMM_WORLD,mpierr)
#endif
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_comm)
        endif

#ifndef MPI
          this%chilocal2RDyn(:,:,:,ispin)=this%chilocalRDyn(:,:,:)
#endif
          !
        end if ! pol%nonblocking_cyclic


        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_array_alloc)
        endif
        if ( .not. pol%nonblocking_cyclic ) then
          ! deallocate
          SAFE_DEALLOCATE(this%chilocalRDyn)
          SAFE_DEALLOCATE(this%gmeRDyntempr2)
          SAFE_DEALLOCATE(this%gmeRDyntempc)
        end if
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_array_alloc)
        endif

      end do ! ispin

      SAFE_DEALLOCATE(row_indx)
      SAFE_DEALLOCATE(col_indx)

      ! finalize non-blocking cyclic communication
      if(pol%nonblocking_cyclic) then
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_comm)
        endif

          ! make sure the buffer is received
          CALL MPI_Wait(req_rec,stat,mpierr)
          ! accumulate contribution into receiving buffer
          !XXXX  buf_rec(:,:,:) = buf_rec(:,:,:) + buf_temp
          do ispin = 1 , kp%nspin
            do kk = 1, pol%nfreq_in_group
!$OMP PARALLEL DO collapse(2)
              do ii = 1, nrec_col
                do jj = 1, nrec_row
                  buf_rec(jj,ii,kk,ispin) = buf_rec(jj,ii,kk,ispin) + buf_temp(jj,ii,kk,ispin)
                enddo
              enddo
            end do
          end do
          ! wait for the massage to be sent
          CALL MPI_Wait(req_send,stat,mpierr)
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_comm)
          endif
  
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_array_alloc)
          endif
  
          ! copy the messega to the sending buffer for the next cycle
          !
          do ispin = 1 , kp%nspin
            do kk = 1, pol%nfreq_in_group
!$OMP PARALLEL DO collapse(2)
              do ii = 1, nrec_col
                do jj = 1, nrec_row
                  buf_send(jj,ii,kk,ispin) = buf_rec(jj,ii,kk,ispin)
                enddo
              enddo
            end do
          end do

          !XXXX buf_send = buf_rec
          !
          nsend_row = nrec_row
          nsend_col = nrec_col
          !
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_array_alloc)
          endif

      end if ! pol%nonblocking_cyclic

    END DO ! ipe
    call progress_free(prog_info)
    !
#ifdef MPI
    call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif
    !
    if(peinf%inode .eq. 0) then
      call timing%stop(timing%chi_sum_sub_omega_neq_0)
    endif


    if(pol%nonblocking_cyclic) then
      ! done
      do ispin = 1 , kp%nspin
        do kk = 1, pol%nfreq_in_group
          !$OMP PARALLEL DO collapse(2)
          do ii = 1, my_sub_nlocal_col
            do jj = 1, my_sub_nlocal_row
              this%chilocal2RDyn(jj,ii,kk,ispin) = buf_send(jj,ii,kk,ispin)
            enddo
          enddo
        end do
      end do
      !XXXX this%chilocal2RDyn(:,:,:,:) = buf_send(:,:,:,:)
      !
      ! clean-up
      SAFE_DEALLOCATE(buf_send)
      SAFE_DEALLOCATE(buf_rec)
      SAFE_DEALLOCATE(buf_temp)
      SAFE_DEALLOCATE(this%chilocalRDyn)
      SAFE_DEALLOCATE(this%gmeRDyntempr2)
      SAFE_DEALLOCATE(this%gmeRDyntempc)
      !
    end if

    ! deallocate stuff
    SAFE_DEALLOCATE(gme_sub)
    SAFE_DEALLOCATE(row_col_loc)
    SAFE_DEALLOCATE(grid2D_2inode)

    IF(pol%need_full_chi) THEN
      ! here allocate chiRdyson (case for which we need full chi)
      SAFE_ALLOCATE(pol%chiRDyn, (scal%npr,scal%npc,pol%nfreq_in_group,kp%nspin))
      pol%chiRDyn = ZERO

      ! initialize descriptor for the eigenvector matrix
      call descinit(desca, nmtx, nmtx, scal%nbl, scal%nbl, 0, 0, &
                    scal%icntxt, max(1,scal%npr), info)
      if(info < 0) then
        write(0,'(a,i3,a)') 'Argument number ', -info, ' had an illegal value on entry.'
        call die("descinit error for descaA in sub_chi_sum")
      end if

      ! transform to original basis use C_Pgemm as temporary array
      SAFE_ALLOCATE(C_Pgemm, (scal%npr,scal%npc))
      C_Pgemm = ZERO
      do ispin = 1, kp%nspin
        do jj = 1, pol%nfreq_in_group
          if(peinf%inode .eq. 0) then
            call timing%start(timing%subspace_pgemm)
          endif
      
          CALL pzgemm('N','N', nmtx, neig_sub, neig_sub, (1.0d0,0.0d0), eigen_vect, 1, 1, desca, &
                      this%chilocal2RDyn(:,:,jj,ispin), 1, 1, desc_sub, (0.0d0,0.0d0), &
                      C_Pgemm, 1, 1, desca)
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%subspace_pgemm)
          endif
          if(peinf%inode .eq. 0) then
            call timing%start(timing%subspace_pgemm)
          endif
          CALL pzgemm('N','C', nmtx, nmtx, neig_sub, (1.0d0,0.0d0), C_Pgemm, 1, 1, desca, &
                      eigen_vect, 1, 1, desca, (0.0d0,0.0d0), &
                      pol%chiRDyn(:,:,jj,ispin), 1, 1, desca)
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%subspace_pgemm)
          endif
        end do
      end do
      SAFE_DEALLOCATE(C_Pgemm)
      SAFE_DEALLOCATE(this%chilocal2RDyn)

      ! remove coulomb potential
      do ispin = 1, kp%nspin
        ! go back to polarizability (remove v^{1/2})
        icurr=0
        do i=1, nmtx
          irow=MOD(INT(((i-1)/scal%nbl)+TOL_SMALL),scal%nprow)
          if(irow.ne.scal%myprow) cycle

          vc = SQRT(vcoul(i))

          do j=1, nmtx
            icol=MOD(INT(((j-1)/scal%nbl)+TOL_SMALL),scal%npcol)
            if(icol .eq. scal%mypcol) then
              icurr=icurr+1
              irowm=INT((icurr-1)/scal%npc+TOL_SMALL)+1
              icolm=MOD((icurr-1),scal%npc)+1

              pol%chiRDyn(irowm,icolm,1:pol%nfreq,ispin) = pol%chiRDyn(irowm,icolm,1:pol%nfreq,ispin) / vc / SQRT(vcoul(j))

            endif
          end do
        end do
      end do
    ELSE
      ! save information for epsinv_sub
      pol%nrow_local_sub = my_sub_nlocal_row
      pol%ncol_local_sub = my_sub_nlocal_col
      pol%neig_sub = neig_sub

      SAFE_ALLOCATE(pol%chiRDyn, (MAX(1,my_sub_nlocal_row),MAX(1,my_sub_nlocal_col),pol%nfreq_in_group,kp%nspin))
      pol%chiRDyn = ZERO
      pol%chiRDyn(:,:,:,:) = this%chilocal2RDyn(:,:,:,:)
      SAFE_DEALLOCATE(this%chilocal2RDyn)
    END IF

    POP_SUB(calc_pol_sub_trunc_easy)

  end subroutine calc_pol_sub_trunc_easy

  subroutine dump_matrix(scal, matrix, inode, unit)
    type (scalapack), intent(in) :: scal
    SCALAR, intent(in)           :: matrix(scal%npr,scal%npc)
    integer :: inode, unit
    character(len=7) :: nam
    integer :: jjj

    PUSH_SUB(dump_matrix)

    ! write(nam,'(I5)') unit+inode
    ! open(unit+inode,file='chi.'//TRIM(nam),form='unformatted',status='replace')
    do jjj = 1, scal%npc
      write(unit+inode) matrix(1:scal%npr,jjj)
    end do
    flush(unit+inode)
    ! close(unit+inode)

    POP_SUB(dump_matrix)

  end subroutine dump_matrix

  subroutine read_dump_matrix(scal, matrix, inode, unit)
    type (scalapack), intent(in) :: scal
    SCALAR, intent(inout)        :: matrix(scal%npr,scal%npc)
    integer :: inode, unit
    character(len=7) :: nam
    integer :: jjj

    PUSH_SUB(read_dump_matrix)

    ! write(nam,'(I5)') unit+inode
    ! open(unit+inode,file='chi.'//TRIM(nam),form='unformatted',status='old')
    do jjj = 1, scal%npc
      read(unit+inode) matrix(1:scal%npr,jjj)
    end do
    ! close(unit+inode)

    POP_SUB(read_dump_matrix)

  end subroutine

#endif
! close the condition ( if defined MPI && defined USESCALAPACK) for subspace method

end module chi_summation_m
