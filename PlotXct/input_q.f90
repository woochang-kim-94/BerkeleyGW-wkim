#include "f_defs.h"

module input_q_m

  use checkbz_m
  use fullbz_m
  use global_m
  use input_common_m
  use misc_m
  use plotxct_common_m
#ifdef HDF5
  use hdf5
  use wfn_io_hdf5_m
#endif
  use wfn_rho_vxc_io_m
  implicit none

  private

  public :: &
    input_q

contains

!-----------------------------------------------------------------------
subroutine input_q(crys,gvec,kg,kgq,syms,xct,indexq,pxct)
!-----------------------------------------------------------------------
!
!     input: crys, gvec, kg,  syms, xct types
!
!     output: kgq type
!             indexq(1:xct%nkpt_fi_fi) : mapping of points in shifted grid
!             INT_VWFNQ_* files
!
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (grid), intent(in) :: kg
  type (grid), intent(out) :: kgq
  type (symmetry), intent(in) :: syms
  type (xctinfo), intent(inout) :: xct
  integer, intent(out) :: indexq(xct%nkpt_fi)
  type (plotxct_t), intent(in) :: pxct

  type (input_reader_t) :: inp !< the input reader `object`
  type (crystal) :: crysq
  type (kpoints) :: kpq
  type (symmetry) :: symsq
  character :: filenamev*20
  character :: tmpfn*16, errmsg*100
  integer :: iunit_v,dummygvec(1,1)
  logical :: send_any !< = any(send_to)
  integer :: ii,jj,ik,ikq,irkq
  logical :: irkq_match
  real(DP) :: delta,qq(3),tol

  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvecq

  logical :: skip_checkbz

  integer, allocatable :: gvec_kpt_components_all(:,:)
  integer, allocatable :: ib_size_array(:)
  integer :: iband_min, iband_max, tot_nbands, band_block, min_band_block
  integer :: ib_first, ib_last, ib_size, ib_size_max
  integer :: ipe, npes_hdf5, ii_start, ii_end, ii_loc, error, ngktot
#ifdef HDF5
  integer(HID_T) :: file_id
  integer(HID_T) :: plist_id
#endif
  SCALAR, allocatable :: my_cg_all(:,:,:)

  PUSH_SUB(input_q)

!-----------------------------------------------------------------------
!     Read data in WFNq_fi


  sheader = 'WFN'
  iflavor = 0
  if ( pxct%use_wfn_hdf5 ) then
    if (peinf%inode==0) write(6,'(a)') ' Reading header of WFNq_fi.h5'
#ifdef HDF5
    call read_hdf5_header_type('WFNq_fi.h5', sheader, iflavor, kpq, gvecq, symsq, crysq)
#endif
  else
    if(peinf%inode == 0) call open_file(26,file='WFNq_fi',form='unformatted',status='old')
    call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, symsq, crysq, warn = .false.)
  endif

  if (pxct%nspinor .ne. kpq%nspinor) then
    call die('Mismatch between expected number of spinors from input and WFNq_fi', only_root_writes = .true.)
  endif

  if (pxct%use_wfn_hdf5) then
    call check_header('WFNq_fi.h5', kpq, gvec, syms, crys, 'WFNq_fi.h5', kpq, gvecq, symsq, crysq, is_wfn = .true.)
  else
    call check_header('WFN_fi', kpq, gvec, syms, crys, 'WFNq_fi', kpq, gvecq, symsq, crysq, is_wfn = .true.)
  endif
  ! there is no kp object in this code??

  kpq%nvband=minval(kpq%ifmax(:,:)-kpq%ifmin(:,:))+1
  if(xct%nvb_fi.gt.kpq%nvband) then
    write(errmsg,'(a,i6,a,i6,a)') 'You requested ', xct%nvb_fi, ' valence bands but WFNq_fi contains only ', kpq%nvband, '.'
    call die(errmsg, only_root_writes = .true.)
  endif
  
  
  if( .not. pxct%use_wfn_hdf5 ) then
    call read_binary_gvectors(26, gvec%ng, gvec%ng, dummygvec, dont_read = .true.)
  endif

!
!     Define shifted irreducible BZ, kgq%r, and define the shift vector
!     (if not done before). Make sure it is right!
!
  kgq%nr=kpq%nrk
  SAFE_ALLOCATE(kgq%r, (3,kgq%nr))
  kgq%r(1:3,1:kgq%nr)=kpq%rk(1:3,1:kpq%nrk)
  xct%qshift= sqrt( DOT_PRODUCT(xct%shift(:), &
    MATMUL(crys%bdot,xct%shift(:) )) )
  if (xct%qshift.eq.0.d0) then
    xct%shift(:)=kgq%r(:,1)-kg%r(:,1)
    xct%qshift= sqrt( DOT_PRODUCT(xct%shift(:), &
      MATMUL(crys%bdot,xct%shift(:) )) )
  endif
  if(peinf%inode.eq.0) write(6,90) xct%shift(:),xct%qshift
90 format(/,1x,'Shift vector:',3(1x,f9.6),' ; Length = ',f0.6)

!-----------------------------------------------------------------------
!     Generate full brillouin zone from irreducible wedge, rk -> fk

  if (.not.pxct%unfoldq) then
    call fullbz(crys,symsq,kgq,1,skip_checkbz,wigner_seitz=.true.,paranoid=pxct%bz_paranoid)
  else
    call fullbz(crys,symsq,kgq,symsq%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=pxct%bz_paranoid)
  endif
  tmpfn='WFNq_fi'
  if (.not. skip_checkbz) then
    call checkbz(kgq%nf,kgq%f,kpq%kgrid,kpq%shift,crys%bdot, &
      tmpfn,'k',.true.,xct%freplacebz,xct%fwritebz)
  endif
  call logit('input_q:  done unfolding/checking BZ')

!-----------------------------------------------------------------------
!     Find correspondence with fk from WFNq_fi
!
!     indexq : correspondence between a k-point in the full BZ, kg%f, and
!       its shifted vector, in kgq%f
!     tol : tolerance

  tol = 1.d-6
  do ik=1,kg%nf
    ikq=0
    delta=0.1d0
    do while((delta.gt.tol).and.(ikq.lt.kgq%nf))
      ikq=ikq+1
      qq(:) = kg%f(:,ik)-(kgq%f(:,ikq)-xct%shift(:))
      qq(:) = qq(:) - anint( qq(:) )
      delta=sqrt(sum(qq(:)**2))
    enddo
    if(delta.gt.tol) then
      if(peinf%inode.eq.0) then
        write(0,'(a,3f12.6)') 'Could not find point equivalent to ', (kg%f(ii,ik),ii=1,3)
      endif
      call die('k-point mismatch between WFN_fi and WFNq_fi.', only_root_writes = .true.)
    else
!
!     make sure that kgq%f(:,ikq)-kg%f(:,ik) = shift vector
!     Near the zone edge, they may differ by a lattice vector
!
      do jj=1,3
        ii = nint( kgq%f(jj,ikq)-kg%f(jj,ik) )
        kgq%f(jj,ikq) = kgq%f(jj,ikq) - dble(ii)
        kgq%kg0(jj,ikq) = kgq%kg0(jj,ikq) - ii
      enddo
      qq(:) = kg%f(:,ik)-(kgq%f(:,ikq)-xct%shift(:))
      delta=sqrt(sum(qq(:)**2))
      if (delta.gt.tol) then
        call die('k-point mismatch between WFN_fi and WFNq_fi. Wrong shift', only_root_writes = .true.)
      endif
      indexq(ik)=ikq
    endif
  enddo
  call logit('input_q:  done mapping kpts (1)')

!-----------------------------------------------------------------------
!     Read the wavefunctions and create INT_VWFNQ_*

  call logit('input_q:  reading WFNq_fi')

  if(peinf%inode.lt.10000) then
    write(filenamev,'(a,i4.4)') 'INT_VWFNQ_', peinf%inode
  else
    call die('input_q: cannot use more than 10000 nodes')
  endif
  iunit_v=128+(2*peinf%inode)+2
  call open_file(iunit_v,file=filenamev,form='unformatted',status='replace')

  !FHJ : set up input reader `object`
  call init_reader(inp, 26, iunit_v, 'WFNq_fi', kpq, xct, gvec)

  if( pxct%use_wfn_hdf5 ) then
    ngktot = SUM(kpq%ngk)
    SAFE_ALLOCATE( gvec_kpt_components_all, (3,ngktot) )
    ! MHN : Read gvecs for all k-points from WFN_fi.h5
#ifdef HDF5
    call read_hdf5_wfn_gvectors('WFNq_fi.h5', gvec_kpt_components_all, ngktot)
#endif
    iband_min      = MINVAL(kpq%ifmax(:,:)) - xct%nvb_fi + 1
    iband_min      = MAX(iband_min, 1)
    ! MHN: We only need to read valence states from WFNq_fi
    iband_max      = MAXVAL(kpq%ifmax(:,:)) !+ xct%ncb_fi
    iband_max      = MIN(iband_max, kpq%mnband)
    tot_nbands     = iband_max - iband_min + 1
    band_block     = (tot_nbands + peinf%npes - 1) / peinf%npes
    ! read at least 128 Mb per MPI task
    min_band_block = 128.0D+00 / ( dble(ngktot) * dble(kpq%nspin*kpq%nspinor) * 16.0D+00 /1024.0D+00/1024D+00 )
    min_band_block = MAX(min_band_block, 1)
    if ( xct%wfn_hdf5_min_band_block > 0 )  min_band_block = xct%wfn_hdf5_min_band_block
    if ( min_band_block > band_block ) then
      band_block = min_band_block
    end if

    ib_first = iband_min + band_block * peinf%inode
    ib_last  = min(ib_first + band_block - 1, iband_max)
    ib_size  = ib_last - ib_first + 1
    if ( ib_size < 1 ) then
      ! don`t read
      ib_first = -1
      ib_last  = -1
      ib_size  =  0
    end if
    SAFE_ALLOCATE( ib_size_array, (peinf%npes) )
    ib_size_array = 0
    ib_size_array(peinf%inode+1) = ib_size

#ifdef HDF5
#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, mpierr)
    call h5fopen_f('WFNq_fi.h5', H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id,error)
    call MPI_AllReduce(MPI_IN_PLACE, ib_size_array, peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
    call h5fopen_f('WFNq_fi.h5', H5F_ACC_RDONLY_F, file_id, error)
#endif
#endif
    
    ib_size_max = MAXVAL( ib_size_array )
    SAFE_ALLOCATE( my_cg_all, (ngktot, kpq%nspin*kpq%nspinor, MAX(ib_size,1)) )
    my_cg_all = ZERO
#ifdef HDF5
    ! dummy allocation
    SAFE_ALLOCATE( peinf%does_it_ownv, (1, 1) )
    peinf%does_it_ownv = .true.
    ! MHN : Read bands for all k-points from WFN_fi.h5
    ! Each process reads a band block of ib_size
    call read_hdf5_bands_block(file_id, kpq, ib_size_max, ib_size, &
                               peinf%does_it_ownv, ib_first, my_cg_all, is_bse = .true.)

    SAFE_DEALLOCATE_P(peinf%does_it_ownv)

    call h5fclose_f(file_id, error)
#endif
  endif

  npes_hdf5 = 0
  if ( pxct%use_wfn_hdf5 ) npes_hdf5 = peinf%npes-1

  do irkq = 1, kpq%nrk

    inp%ng = kpq%ngk(irkq)

    ! FHJ: this scales as nk^2, but could scale as nk*log(nk)
    irkq_match = .false.
    do ii=1,kg%nf
      if (irkq == kgq%indr(indexq(ii))) then
        irkq_match = .true.
        exit
      endif
    enddo

    ! FHJ: Determine to which PE the wavefunctions for this k-point need to be sent.
    send_any = .false.
    inp%send_to(:) = .false.
    do jj=1,peinf%npes
      do ii=1, peinf%ikt(jj)
        if(kgq%indr(indexq(peinf%ik(jj,ii))).eq.irkq) then
          inp%send_to(jj) = .true.
          send_any = .true.
          exit
        endif
      enddo
    enddo

    ! FHJ: Don`t bother reading data if we don`t need the kpt
    if((.not. send_any) .or. (.not. irkq_match)) then
        ! MHN: Call skip_kpt only for binary read
        if (.not. pxct%use_wfn_hdf5) call skip_kpt(inp)
        cycle
    endif

    ! FHJ: Read gvectors` indices, store in inp%isort
    if( pxct%use_wfn_hdf5 ) then
      call read_gvecs(inp, irkq, gvec, pxct%use_wfn_hdf5, gvec_kpt_components_all, kpq)
    else
      call read_gvecs(inp, irkq, gvec, pxct%use_wfn_hdf5)
    endif

    ! FHJ: Read all bands and distribute
    if( pxct%use_wfn_hdf5 ) then
      call read_bands(inp, irkq, pxct%use_wfn_hdf5, npes_hdf5, kpq, ib_size_array, ib_size_max, iband_min,& 
                      iband_max, band_block, my_cg_all)
    else
      call read_bands(inp, irkq, pxct%use_wfn_hdf5, npes_hdf5)
    endif

  enddo !end loop over k-points

  call free_reader(inp)
  if( pxct%use_wfn_hdf5 ) then
    SAFE_DEALLOCATE(gvec_kpt_components_all)
  endif
  
  if (peinf%inode.eq.0) then
    call kpq%free()
  endif
  
!-------------------------------
! Write out info about xtal

  if(peinf%inode.eq.0) then
    write(6,4004)
4004 format(/1x,'Crystal wavefunctions read from WFNq_fi')
    write(6,3007) kgq%nr
3007 format(1x,'- nrk = ',i0)
    write(6,'(12x,3f10.4)') ((kg%r(ii,jj),ii=1,3),jj=1,kg%nr)
    write(6,3070) kgq%nf
3070 format(1x,'- nfk = ',i0)
    if (.not. pxct%use_wfn_hdf5) then
      call close_file(26)
    endif
  endif !end if(inode.eq.0)
  
  call close_file(iunit_v)

  ! only needed for comm_disk
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif
  call logit('input_q:  done reading WFNq_fi')

  POP_SUB(input_q)

  return
end subroutine input_q

end module input_q_m
