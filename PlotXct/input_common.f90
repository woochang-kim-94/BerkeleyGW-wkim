#include "f_defs.h"

module input_common_m

  use global_m
  use wfn_rho_vxc_io_m
  use misc_m
  implicit none

  private

  !> the input reader `object`
  type input_reader_t
    character(len=7) :: fname !< WFN_fi or WFNq_fi
    logical :: is_cond   !< .true. if using WFN_fi, .false. if WFNq_fi
    integer :: iunit     !< unit corresponding to the WFN file
    integer :: iunit_xct !< unit corresponding to the INT* file
    integer :: nk        !< number of k-points in the WFN file
    integer :: ng        !< number of g vectors in the current k-point
    integer :: nb_wfn    !< number of bands in the WFN file
    integer :: nv_wfn    !< number of valence bands in the WFN file
    integer :: nv_xct    !< number of valence bands in the XCT file
    integer :: nc_xct    !< number of conduction bands in the XCT file
    integer :: ns        !< number of spin components
    integer :: nspinor   !< number of spinor components
    integer :: gvec_ng
    logical, pointer :: send_to(:)
    integer, pointer :: isort(:)
  end type input_reader_t

  public :: input_reader_t, init_reader, free_reader, skip_kpt, read_gvecs, read_bands

contains

  !> Initialize the input reader `object`
  subroutine init_reader(inp, iunit, iunit_xct, fname, kpts, xct, gvec)
    type (input_reader_t), intent(out) :: inp !< the input reader `object`
    integer, intent(in) :: iunit              !< unit of the WFN file
    integer, intent(in) :: iunit_xct          !< unit of the WFN file
    character(len=*), intent(in) :: fname     !< WFN file name
    type (kpoints), intent(in) :: kpts        !< kp/kpq
    type (xctinfo), intent(in) :: xct
    type (gspace), intent(in) :: gvec

    PUSH_SUB(init_reader)

    inp%iunit = iunit
    inp%fname = fname
    inp%is_cond = fname=='WFN_fi'
    inp%iunit_xct = iunit_xct
    inp%nk = kpts%nrk
    inp%nb_wfn = kpts%mnband
    inp%nv_wfn = kpts%nvband
    inp%nv_xct = xct%nvb_fi
    inp%nc_xct = xct%ncb_fi
    inp%ns = kpts%nspin
    inp%nspinor = kpts%nspinor
    inp%gvec_ng = gvec%ng
    SAFE_ALLOCATE(inp%isort, (gvec%ng))
    SAFE_ALLOCATE(inp%send_to, (peinf%npes))
    inp%isort(:) = 0

    POP_SUB(init_reader)

  end subroutine init_reader

  !> Deallocates buffers from input reader `object`
  subroutine free_reader(inp)
    type (input_reader_t), intent(inout) :: inp !< the input reader `object`

    PUSH_SUB(free_reader)

    SAFE_DEALLOCATE_P(inp%isort)
    SAFE_DEALLOCATE_P(inp%send_to)

    POP_SUB(free_reader)

  end subroutine free_reader


  !> Skip the WFNs of the current k-point.
  subroutine skip_kpt(inp)
    type (input_reader_t), intent(in) :: inp !< the input reader `object`

    integer, allocatable :: gvecs_comp(:,:)
    SCALAR, allocatable :: cg(:,:)
    integer :: ib

    PUSH_SUB(skip_kpt)

    SAFE_ALLOCATE(gvecs_comp, (3, inp%ng))
    call read_binary_gvectors(inp%iunit, inp%ng, inp%ng, gvecs_comp, dont_read=.true.)
    SAFE_DEALLOCATE(gvecs_comp)
    SAFE_ALLOCATE(cg, (inp%ng, inp%ns*inp%nspinor))
    do ib=1,inp%nb_wfn
      call read_binary_data(inp%iunit, inp%ng, inp%ng, inp%ns*inp%nspinor, cg, dont_read=.true.)
    enddo
    SAFE_DEALLOCATE(cg)

    POP_SUB(skip_kpt)

  end subroutine skip_kpt

  !> Read all gvectors for the current k-point.
  subroutine read_gvecs(inp, ik, gvec, use_wfn_hdf5, gvec_kpt_components_all,kp)
    type (input_reader_t), intent(inout) :: inp !< the input reader `object`
    integer, intent(in) :: ik    !< kpt in question

    type (gspace), intent(in) :: gvec
    logical, intent(in) :: use_wfn_hdf5
    integer, intent(in), optional :: gvec_kpt_components_all(:,:) ! For HDF5 read
    type (kpoints), intent(in), optional :: kp
    integer, allocatable :: gvecs_comp(:,:)
    integer :: ipe, ig, istart
    logical :: should_write

    PUSH_SUB(read_gvecs)

    should_write = inp%send_to(peinf%inode+1)
    ! FHJ: Idle nodes can leave the function
    if (.not. (peinf%inode==0 .or. should_write)) then 
      POP_SUB(read_gvecs)
      return
    endif

    ! FHJ: If we need the kpt, read gvector indices.
    ! Note: we don`t bcast, and only the receiving nodes call findvector
    SAFE_ALLOCATE(gvecs_comp, (3, inp%ng))
    if (use_wfn_hdf5) then
      if (ik>1) then
        istart = sum(kp%ngk(1:ik-1)) + 1
      else
        istart = 1
      endif
      gvecs_comp(:,:) = gvec_kpt_components_all(1:3,istart:istart+kp%ngk(ik)-1)
    else
      call read_binary_gvectors(inp%iunit, inp%ng, inp%ng, gvecs_comp, bcast=.false.)
#ifdef MPI
      if (peinf%inode==0) then
        do ipe=2,peinf%npes !FHJ: root already has the data!
          if (inp%send_to(ipe)) then
            call MPI_Send(gvecs_comp, 3*inp%ng, MPI_INTEGER, ipe-1, ik, MPI_COMM_WORLD, mpierr)
          endif
        enddo
      else 
        call MPI_Recv(gvecs_comp, 3*inp%ng, MPI_INTEGER, &
          0, ik, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
      endif
#endif
    endif
    if (should_write) then
      do ig = 1, inp%ng
        call findvector(inp%isort(ig), gvecs_comp(:, ig), gvec)
        if(inp%isort(ig) == 0) call die('could not find gvec')
      enddo
    endif
    SAFE_DEALLOCATE(gvecs_comp)

    POP_SUB(read_gvecs)

  end subroutine read_gvecs

  !> Read all WFNs the current k-point.
  subroutine read_bands(inp, ik, use_wfn_hdf5, npes_hdf5, kp,ib_size_array, ib_size_max, iband_min, iband_max, band_block, my_cg_all)
    type (input_reader_t), intent(inout) :: inp !< the input reader `object`
    integer, intent(in) :: ik    !< kpt in question
    logical, intent(in) :: use_wfn_hdf5 !< read hdf5
    integer, intent(in) :: npes_hdf5 !< number of processes for hdf5 read

    !MHN: Optional inputs for hdf5 read
    type (kpoints), intent(in), optional :: kp
    integer, intent(in), optional :: ib_size_array(:) !< number of bands per process
    integer, intent(in), optional :: ib_size_max !< max number of bands per process
    integer, intent(in), optional :: iband_min !< min band of the block
    integer, intent(in), optional :: iband_max !< max band of the block
    integer, intent(in), optional :: band_block !< number of bands in block
    SCALAR, intent(in), optional :: my_cg_all(:,:,:) !< Distributed cgs (band wise)

    SCALAR, allocatable :: cg_xct(:,:,:)
    SCALAR, allocatable :: ipe_cg_all(:,:,:)
    SCALAR :: cg_dummy(1,1)
    integer :: ipe, ig, ib_wfn, ib_xct, ib_xct_max, is, tag
    logical :: should_write
    logical :: cond
    integer :: ipe_hdf5, ii_loc, ii_end, ii_start, istart, ngktot
    integer :: ib_size, iend, ngkmax

    PUSH_SUB(read_bands)

    should_write = inp%send_to(peinf%inode+1)
    ! FHJ: Idle nodes can leave the function
    if (.not. use_wfn_hdf5) then
      if (.not. (peinf%inode==0 .or. should_write)) then 
        POP_SUB(read_bands)
        return
      endif
    endif

    if (inp%is_cond) then !input
      ib_xct_max = inp%nc_xct       
    else !input_q
      ib_xct_max = inp%nv_xct
    endif

    ! FHJ: note that we use a non-default order for the indices of cg_all.
    ! This makes MPI_Sends/Recvs easier
    SAFE_ALLOCATE(cg_xct, (inp%ng,inp%ns*inp%nspinor,ib_xct_max))
    if (use_wfn_hdf5)  then
      ngkmax = MAXVAL(kp%ngk)
      SAFE_ALLOCATE( ipe_cg_all, (ngkmax, kp%nspin*kp%nspinor, ib_size_max) )
    endif

    ! MHN: For HDF5 read, loop over the HDF5 processes and 
    ! broadcast the owned band block for each ipe_hdf5
  do ipe_hdf5 = 0, npes_hdf5
    if ( use_wfn_hdf5 ) then
      istart = sum(kp%ngk(1:ik-1)) + 1
      iend = istart+kp%ngk(ik)-1
      if ( ib_size_array(ipe_hdf5+1) > 0) then
        ipe_cg_all = ZERO
        if ( ipe_hdf5 == peinf%inode ) then
          ib_size = ib_size_array(ipe_hdf5+1)
          ! ipe_cg_all stores the block of bands on ipe_hdf5 process 
          ! for the current k-point, ik
          ipe_cg_all(1:kp%ngk(ik), :, 1:ib_size) = my_cg_all(istart:iend, :, 1:ib_size)
        end if
#ifdef MPI

        call MPI_Bcast(ipe_cg_all, ngkmax * kp%nspin*kp%nspinor * ib_size_max, &
                     MPI_SCALAR, ipe_hdf5, MPI_COMM_WORLD, mpierr)
#endif
      else
         ! no elements on this MPI task, cycle
         cycle
      end if
    end if

    ii_start = 1
    ii_end   = inp%nb_wfn
    if ( use_wfn_hdf5 ) then
      ! MHN: ii_start to ii_end are the band indices in the current 
      ! band block (for hdf5 read)
      ii_start = iband_min + band_block * ipe_hdf5
      ii_end   = min(ii_start+ band_block - 1, iband_max)
    end if
    ii_loc = 0
    do ib_wfn=ii_start,ii_end
      ii_loc = ii_loc + 1
      if (inp%is_cond) then !input
        ib_xct = ib_wfn - inp%nv_wfn
      else !input_q
        ib_xct = inp%nv_wfn - ib_wfn + 1
      endif

      cond = (ib_xct > 0) .and. (ib_xct <= ib_xct_max)

      ! If ib_xct is one of the selected bands...
      if (cond) then
        if ( use_wfn_hdf5 ) then
          if (ik>1) then
            istart = sum(kp%ngk(1:ik-1)) + 1
          else
            istart = 1
          endif
          cg_xct(1:kp%ngk(ik),:,ib_xct) = ipe_cg_all(1:kp%ngk(ik), 1:kp%nspin*kp%nspinor, ii_loc)
        else
          call read_binary_data(inp%iunit, inp%ng, inp%ng, inp%ns*inp%nspinor, &
            cg_xct(:,:,ib_xct), bcast=.false.)
#ifdef MPI
          tag = (ib_wfn-1)*inp%nk + ik
          if (peinf%inode==0) then
            do ipe=2,peinf%npes
              if (inp%send_to(ipe)) then
                call MPI_Send(cg_xct(1,1,ib_xct), inp%ng*inp%ns*inp%nspinor, MPI_SCALAR, &
                  ipe-1, tag, MPI_COMM_WORLD, mpierr)
              endif
            enddo
          else
            call MPI_Recv(cg_xct(1,1,ib_xct), inp%ng*inp%ns*inp%nspinor, MPI_SCALAR, &
              0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
          endif
#endif
        endif

        if (should_write) &
          call checknorm(inp%fname,ib_wfn,ik,inp%ng,inp%ns,inp%nspinor,cg_xct(:,:,ib_xct))

      else

        ! FHJ: not a selected band => ignore. cg_xct should not be touched!!
        if (.not. use_wfn_hdf5) then
          call read_binary_data(inp%iunit, 1, 1, 1, cg_dummy(:,:), dont_read=.true.)
        endif
      endif !ib_xct is one of the selected valence band

    enddo
  enddo ! ipe_hdf5 loop ends

    if (should_write) then
      write(inp%iunit_xct) ik,inp%ng,ib_xct_max,inp%ns,inp%nspinor
      write(inp%iunit_xct) (inp%isort(ig),ig=1,inp%gvec_ng), &
        (((cg_xct(ig,is,ib_xct),ig=1,inp%ng),ib_xct=1,ib_xct_max),is=1,inp%ns*inp%nspinor)
    endif

    SAFE_DEALLOCATE(cg_xct)
    if (use_wfn_hdf5) then
      SAFE_DEALLOCATE(ipe_cg_all)
    endif
    POP_SUB(read_bands)


  end subroutine read_bands

end module input_common_m
