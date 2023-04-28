#include "f_defs.h"

module input_m

  use checkbz_m
  use fullbz_m
  use global_m
  use input_common_m
  use input_utils_m
  use misc_m
  use plotxct_common_m
  use sort_m
#ifdef HDF5
  use hdf5
  use wfn_io_hdf5_m
#endif
  use wfn_rho_vxc_io_m
  use distrib_m
  implicit none

  private

  public :: &
    input

contains

!-----------------------------------------------------------------------
subroutine input(crys,gvec,kg,syms,xct,kgr,index_k,pxct)
!-----------------------------------------------------------------------
!
!       Read data from file WFN_fi and initialize variables
!
!     input: xct types
!
!     output: crys,gvec,kg,syms types
!             peinf type (from distrib.f90)
!             INT_VWFN_* and INT_CWFN_* files
!
!     Copied from BSE/input.f90, without eqp arrays.
!
  type (crystal), intent(out) :: crys
  type (gspace), intent(out) :: gvec
  type (grid), intent(out) :: kg
  type (symmetry), intent(out) :: syms
  type (xctinfo), intent(inout) :: xct
  real(DP), intent(in) :: kgr(3,xct%nkpt_fi)
  integer, intent(in) :: index_k(xct%nkpt_fi)
  type (plotxct_t), intent(in) :: pxct

  type (input_reader_t) :: inp !< the input reader `object`
  type (kpoints) :: kp
  type (grid) :: kgt

  character :: filenamec*20
  character :: tmpfn*16, errmsg*100
  integer :: iunit_c
  logical :: send_any !< = any(send_to)
  integer :: ii,jj,irk
  integer :: irks
  real(DP) :: diffvol,vcell,kt(3)
  real(DP) :: tol
  real(DP), allocatable :: ek_tmp(:)
  integer, allocatable :: indxk(:),k_tmp(:,:)
  integer, allocatable :: index(:)
  integer, allocatable :: gvec_kpt_components_all(:,:)
  integer, allocatable :: ib_size_array(:)
  integer :: iband_min, iband_max, tot_nbands, band_block, min_band_block
  integer :: ib_first, ib_last, ib_size, ib_size_max
  integer :: ipe, npes_hdf5, ii_start, ii_end, ii_loc, error
  ! WFN HDF5 stuff
#ifdef HDF5
  integer(HID_T) :: file_id
  integer(HID_T) :: plist_id
#endif
  SCALAR, allocatable :: my_cg_all(:,:,:)
  
  character(len=3) :: sheader
  integer :: iflavor
  integer :: ngktot

  logical :: skip_checkbz

  PUSH_SUB(input)

  call logit('entering input')


  sheader = 'WFN'
  iflavor = 0
  if ( pxct%use_wfn_hdf5 ) then
    if (peinf%inode==0) write(6,'(a)') ' Reading header of WFN_fi.h5'
#ifdef HDF5
    call read_hdf5_header_type('WFN_fi.h5', sheader, iflavor, kp, gvec, syms, crys)
#endif
  else
    if(peinf%inode == 0) call open_file(25,file='WFN_fi',form='unformatted',status='old')
    call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys)
  endif

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  if( pxct%use_wfn_hdf5 ) then
#ifdef HDF5
    call read_hdf5_gvectors('WFN_fi.h5', gvec%ng, gvec%components)
#endif
  else
    call read_binary_gvectors(25, gvec%ng, gvec%ng, gvec%components)
  endif
  
  call get_volume(vcell,crys%bdot)
  diffvol=abs(crys%celvol-vcell)
  if (diffvol.gt.1.0d-6) then
    call die('volume mismatch.', only_root_writes = .true.)
  endif
  
  kp%nvband=minval(kp%ifmax(:,:)-kp%ifmin(:,:))+1
  kp%ncband=kp%mnband-maxval(kp%ifmax(:,:))
  
  if(xct%nvb_fi.gt.kp%nvband) then
    write(errmsg,'(a,i6,a,i6,a)') 'You requested ', xct%nvb_fi, ' valence bands but WFN_fi contains only ', kp%nvband, '.'
    call die(errmsg, only_root_writes = .true.)
  endif
  if(xct%ncb_fi.gt.kp%ncband) then
    write(errmsg,'(a,i6,a,i6,a)') 'You requested ', xct%ncb_fi, ' conduction bands but WFN_fi contains only ', kp%ncband, '.'
    call die(errmsg, only_root_writes = .true.)
  endif
  if(xct%nspin.ne.kp%nspin) then
    write(errmsg,'(a,2i6)') 'Number of spins mismatch: ', xct%nspin, kp%nspin
    call die(errmsg, only_root_writes = .true.)
  endif
  if(pxct%nspinor.ne.kp%nspinor) then
    write(errmsg,'(a,2i6)') 'Number of spinors mismatch: ', pxct%nspinor, kp%nspinor
    call die(errmsg, only_root_writes = .true.)
  endif
!-----------------------------------------------------------------------
!     Check if all k-points are available and define grid
!
  tol = 1.d-4
  kgt%nr=kp%nrk
  SAFE_ALLOCATE(kgt%r, (3,kgt%nr))
  kgt%r(1:3,1:kgt%nr)=kp%rk(1:3,1:kp%nrk)
  if (.not.pxct%unfold) then
    call fullbz(crys,syms,kgt,1,skip_checkbz,wigner_seitz=.true.,paranoid=pxct%bz_paranoid)
  else
    call fullbz(crys,syms,kgt,syms%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=pxct%bz_paranoid)
  endif
  tmpfn='WFN_fi'
  if (.not. skip_checkbz) then
    call checkbz(kgt%nf,kgt%f,kp%kgrid,kp%shift,crys%bdot, &
      tmpfn,'k',.true.,xct%freplacebz,xct%fwritebz)
  endif
  call logit('input:  done unfolding/checking BZ')
!
  ! FHJ: map the xct grid kgr(:,:) to the WFN grid kgt%f(:,:)
  SAFE_ALLOCATE(indxk, (xct%nkpt_fi))
  indxk=0
  do jj=1,xct%nkpt_fi
    do ii=1,kgt%nf
      kt(:) = mod(kgr(:,jj) - kgt%f(:,ii)+10.0,1.0d0)
      if (all(abs(kt(1:3)).lt.tol)) then
        if (indxk(jj).ne.0) write(0,*) 'WARNING: multiple definition of k-point',jj,indxk(jj),kgr(:,jj)
        indxk(jj)=ii
      endif
    enddo
!
!     If some k-point listed in kgr is not found in WFN_fi, indxk
!     will store zero.
!
    if (indxk(jj).eq.0) then
      write(errmsg,'(a,3f12.6,a)') 'Could not find vector ', kgr(:,jj), ' in WFN_fi'
      call die(errmsg, only_root_writes = .true.)
    endif
  enddo
  call logit('input:  done mapping kpts (1)')
!
!   update kgt -> kg
!
  kg%nr = kgt%nr
  kg%nf = xct%nn
  SAFE_ALLOCATE(kg%r, (3,kg%nr))
  kg%r = kgt%r
  kg%sz = kgt%sz
  SAFE_ALLOCATE(kg%itran, (kg%nf))
  SAFE_ALLOCATE(kg%indr, (kg%nf))
  SAFE_ALLOCATE(kg%f, (3,kg%nf))
  SAFE_ALLOCATE(kg%kg0, (3,kg%nf))
  do jj=1,xct%nn
    kg%itran(jj) = kgt%itran(indxk(index_k(jj)))
    kg%indr(jj) = kgt%indr(indxk(index_k(jj)))
    kg%kg0(:,jj) = kgt%kg0(:,indxk(index_k(jj)))
    kg%f(:,jj) = kgt%f(:,indxk(index_k(jj)))
  enddo
  call dealloc_grid(kgt)
  SAFE_DEALLOCATE(indxk)
!
!     indxk : stores the correspondence between k-points kg%r and kp%rk
!     (it is used to select the set of wavefunctions to be stored)
!     tol : tolerance in the coordinates of k-points
!
  SAFE_ALLOCATE(indxk, (kg%nr))
  indxk=0
  do jj=1,kg%nr
    do ii=1,kp%nrk
      kt(:) = kg%r(:,jj) - kp%rk(:,ii)
      if (all(abs(kt(1:3)).lt.tol)) then
        if (indxk(jj).ne.0) write(0,*) 'WARNING: multiple definition of k-point',jj,indxk(jj),kg%r(:,jj)
        indxk(jj)=ii
      endif
    enddo
!
!     If some k-point listed in kg%r is not found in WFN_fi, indxk
!     will store zero. Later, the job will stop in genwf.
!
    if (indxk(jj).eq.0) write(0,'(a,3f12.6,a)') 'WARNING: could not find vector ',kg%r(:,jj),' in WFN_fi'
  enddo
  call logit('input:  done mapping kpts (2)')

!-----------------------------------------------------------------------
!       Distribute kpoints among the PEs
!
  call logit('input:  calling distrib')
  call distrib(xct)
!

!-----------------------------------------------------------------------
!     Order g-vectors with respect to their kinetic energy
!
  call logit('input:  reordering gvecs')
  SAFE_ALLOCATE(index, (gvec%ng))
  SAFE_ALLOCATE(gvec%ekin, (gvec%ng))
  call kinetic_energies(gvec, crys%bdot, gvec%ekin)
  call sortrx(gvec%ng, gvec%ekin, index, gvec = gvec%components)
  
  SAFE_ALLOCATE(ek_tmp, (gvec%ng))
  ek_tmp = gvec%ekin
  SAFE_ALLOCATE(k_tmp, (3,gvec%ng))
  k_tmp = gvec%components
  do ii=1,gvec%ng
    gvec%ekin(ii) = ek_tmp(index(ii))
    gvec%components(:,ii) = k_tmp(:,index(ii))
  enddo
  SAFE_DEALLOCATE(ek_tmp)
  SAFE_DEALLOCATE(k_tmp)
  SAFE_DEALLOCATE(index)

  call gvec_index(gvec)

!-----------------------------------------------------------------------
!     Read the wavefunctions and create INT_CWFN_*
!

  call logit('input:  reading WFN_fi')

  if(peinf%inode.lt.10000) then
    write(filenamec,'(a,i4.4)') 'INT_CWFN_', peinf%inode
  else
    call die('input: cannot use more than 10000 nodes')
  endif
  iunit_c=128+(2*peinf%inode)+1
  call open_file(iunit_c,file=filenamec,form='unformatted',status='replace')

  !FHJ : set up input reader `object`
  call init_reader(inp, 25, iunit_c, 'WFN_fi', kp, xct, gvec)

  if( pxct%use_wfn_hdf5 ) then
    ! MHN : Read gvecs for all k-points from WFN_fi.h5 
    ngktot = SUM(kp%ngk)
    SAFE_ALLOCATE( gvec_kpt_components_all, (3,ngktot) )
#ifdef HDF5
    call read_hdf5_wfn_gvectors('WFN_fi.h5', gvec_kpt_components_all, ngktot)
#endif
    !iband_min      = MINVAL(kp%ifmax(:,:)) - xct%nvb_fi + 1
    ! MHN: We only read conduction states
    iband_min      = MINVAL(kp%ifmax(:,:))  + 1
    iband_min      = MAX(iband_min, 1)
    iband_max      = MAXVAL(kp%ifmax(:,:)) + xct%ncb_fi
    iband_max      = MIN(iband_max, kp%mnband)
    tot_nbands     = iband_max - iband_min + 1
    band_block     = (tot_nbands + peinf%npes - 1) / peinf%npes
    ! read at least 128 Mb per MPI task
    min_band_block = 128.0D+00 / ( dble(ngktot) * dble(kp%nspin*kp%nspinor) * 16.0D+00 /1024.0D+00/1024D+00 )
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
    call h5fopen_f('WFN_fi.h5', H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id,error)
    call MPI_AllReduce(MPI_IN_PLACE, ib_size_array, peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
    call h5fopen_f('WFN_fi.h5', H5F_ACC_RDONLY_F, file_id, error)
#endif
#endif

    ib_size_max = MAXVAL( ib_size_array )
    SAFE_ALLOCATE( my_cg_all, (ngktot, kp%nspin*kp%nspinor, MAX(ib_size,1)) )
    my_cg_all = ZERO
#ifdef HDF5
    ! dummy allocation
    SAFE_ALLOCATE( peinf%does_it_ownv, (1, 1) )
    ! MHN : Read bands for all k-points from WFN_fi.h5 
    ! Each process reads a band block of ib_size
    peinf%does_it_ownv = .true.
    call read_hdf5_bands_block(file_id, kp, ib_size_max, ib_size, &
                               peinf%does_it_ownv, ib_first, my_cg_all, is_bse = .true.)

    SAFE_DEALLOCATE_P(peinf%does_it_ownv)

    call h5fclose_f(file_id, error)
#endif
  endif

  npes_hdf5 = 0
  if ( pxct%use_wfn_hdf5 ) npes_hdf5 = peinf%npes-1

  do irk=1,kp%nrk


    inp%ng = kp%ngk(irk)

    ! FHJ: this scales as nk^2, but could scale as nk*log(nk)
    irks = 0
    do ii=1,kg%nr
      if (irk.eq.indxk(ii)) then
        irks=ii
        exit
      endif
    enddo
    if (irks==0) then
      call die('Cannot find k-point!', only_root_writes=.true.)
    endif

    ! FHJ: Determine to which PE the wavefunctions for this k-point need to be sent.
    send_any = .false.
    inp%send_to(:) = .false.
    do jj=1,peinf%npes
      do ii=1, peinf%ikt(jj)
        if(kg%indr(peinf%ik(jj,ii))==irks) then
          inp%send_to(jj) = .true.
          send_any = .true.
          exit
        endif
      enddo
    enddo

    ! FHJ: Don`t bother reading data if we don`t need the kpt
    if ( .not. send_any) then
      ! MHN: Call skip_kpt only for binary read
      if (.not. pxct%use_wfn_hdf5) call skip_kpt(inp)
      cycle
    endif

    ! FHJ: Read gvectors` indices, store in wfnv%isort
    if( pxct%use_wfn_hdf5 ) then
      call read_gvecs(inp, irk, gvec, pxct%use_wfn_hdf5,gvec_kpt_components_all, kp)
    else  
      call read_gvecs(inp, irk, gvec, pxct%use_wfn_hdf5)
    endif

    ! FHJ: Read all bands and distribute
    if( pxct%use_wfn_hdf5 ) then
      call read_bands(inp, irk, pxct%use_wfn_hdf5, npes_hdf5, kp, &
        ib_size_array, ib_size_max, iband_min, iband_max, band_block, my_cg_all)
    else
      call read_bands(inp, irk, pxct%use_wfn_hdf5, npes_hdf5)
    endif

  enddo !end loop over k-points

  call free_reader(inp)
  if( pxct%use_wfn_hdf5 ) then
    SAFE_DEALLOCATE(gvec_kpt_components_all)
  endif

  SAFE_DEALLOCATE(indxk)
  call close_file(iunit_c)

  if(peinf%inode.eq.0) then
    write(6,3004)
3004 format(/1x,'Crystal wavefunctions read from WFN_fi')
    write(6,3007) kg%nr
3007 format(1x,'- nrk = ',i0)
    write(6,'(12x,3f10.4)') ((kg%r(ii,jj),ii=1,3),jj=1,kg%nr)
    write(6,3070) kg%nf
3070 format(1x,'- nfk = ',i0)
    if (.not. pxct%use_wfn_hdf5) then
      call close_file(25)
    endif
  endif !end if(inode.eq.0)
  
  call kp%free()

  ! only needed for comm_disk
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif
  call logit('input:  done reading WFN_fi')
  call logit('leaving input')

  POP_SUB(input)

  return
end subroutine input

end module input_m
