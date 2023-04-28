!=============================================================================
!
! Utilities:
!
! (1) epsmat_old2hdf5     Originally by DVF      Last Modified 08/2016 (FHJ)
!
!     This utility converts epsmat files from the old epsilon format to
!     the new HDF5 format, to facilitate people switching to HDF5.
!     The script currently writes the epsmat.h5 file with version 3.
!
!=============================================================================

#include "f_defs.h"

program epsmat_old2hdf5

  use global_m
  use hdf5
  use write_matrix_m
  use epswrite_hdf5_m
  use wfn_rho_vxc_io_m
  use wfn_io_hdf5_m
  use input_utils_m
  use vcoul_generator_m

  implicit none
  
  character :: ajname*6,adate*11,outfile*80,infile*80,fname_wfn*80
  real(DP)  :: qk(3)
  real(DP), allocatable :: epsdiag(:,:,:),ekin(:),qpts(:,:)
  complex(DPC), allocatable :: epsRDyn(:,:,:)
  integer   :: ig,ng,ngq, ii,iq,qgrid(3), &
    nq,nmtxmax,jj,ifreq,error
  SCALAR, allocatable :: eps(:,:)
  integer, allocatable :: isort(:),isorti(:),gcomponents(:,:),nmtx0_of_q(:),nmtx0_of_q_file(:)
  logical, allocatable :: is_imag(:)
  type(mf_header_t) :: mf
  type(polarizability) :: pol
 
  write(6,*) MYFLAVOR // ' version is used to convert old epsmat format to hdf5 format'
  peinf%jobtypeeval = 0

  call h5open_f(error)
  
!------------------
! DVF : get gvector & qpoint information from input file, and the input/output file names
! Want number of gvectors so that we can read in g-vector components before the call
! to write the hdf5 epsmat header - call eps_hdf5_setup
! Want number of qpoints so that we can read in nmtx0_of_q. We want to read in this
! array because otherwise we`d have to read get it from the file by looping through
! the whole file, which can take a long time for big files

  call open_file(unit=10,file='epsmat_old2hdf5.inp',form='formatted',status='old')
  pol%truncval(:) = 0
  read(10,*) infile
  read(10,*) outfile
  read(10,*) fname_wfn
  read(10,*) pol%icutv
  if (pol%icutv==TRUNC_SPHERICAL) then
    read(10,*) pol%truncval(1)
  endif
  close(10)

  ! FHJ: Figure out nq
  call open_file(unit=11,file=TRUNC(infile),form='unformatted',status='old')
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11) nq
  read(11) ng
  call close_file(11)

  SAFE_ALLOCATE(qpts,(3,nq))
  SAFE_ALLOCATE(nmtx0_of_q,(nq))
  SAFE_ALLOCATE(nmtx0_of_q_file,(nq))
  SAFE_ALLOCATE(gcomponents, (3,ng))

  nmtxmax=0

  ! FHJ: Read WFN file, get the header
!#BEGIN_INTERNAL_ONLY
  if (index(fname_wfn,'.h5')/=0) then
    call read_hdf5_mf_header(TRUNC(fname_wfn), mf, sheader='WFN', iflavor=SCALARSIZE, &
      warn=.false., dont_warn_kgrid=.true.)
    SAFE_ALLOCATE(mf%gvec%components, (3, mf%gvec%ng))
    call read_hdf5_gvectors(TRUNC(fname_wfn), mf%gvec%ng, mf%gvec%components)
  else
!#END_INTERNAL_ONLY
    call open_file(unit=25,file=TRUNC(fname_wfn),form='unformatted',status='old')
    call read_mf_header(25, mf, sheader='WFN', iflavor=SCALARSIZE, &
      warn=.false., dont_warn_kgrid=.true.)
    SAFE_ALLOCATE(mf%gvec%components, (3, mf%gvec%ng))
    call read_binary_gvectors(25, mf%gvec%ng, mf%gvec%ng, mf%gvec%components)
    call close_file(25)
!#BEGIN_INTERNAL_ONLY
  endif
!#END_INTERNAL_ONLY

!--------------------------
!  Read in the epsmat file
!
  call open_file(unit=11,file=TRUNC(infile),form='unformatted',status='old')
  read(11) ajname,adate
  read(11) pol%freq_dep, pol%nFreq
  read(11) (qgrid(ii),ii=1,3)
  if (pol%freq_dep==2.or.pol%freq_dep==3) then
    SAFE_ALLOCATE(pol%dFreqGrid,(pol%nFreq))
    SAFE_ALLOCATE(pol%dFreqBrd,(pol%nFreq))
    read(11) (pol%dFreqGrid(ii),ii=1,pol%nFreq),(pol%dFreqBrd(ii),ii=1,pol%nFreq)
    SAFE_ALLOCATE(is_imag, (pol%nFreq))
    is_imag = (dabs(pol%dFreqGrid)<TOL_ZERO)
    is_imag(1) = .false.
    pol%nfreq_imag = count(is_imag)
    SAFE_DEALLOCATE(is_imag)
  else
    pol%nFreq = 1
    pol%nfreq_imag = 0
    SAFE_ALLOCATE(pol%dFreqGrid,(1))
    SAFE_ALLOCATE(pol%dFreqBrd,(1))
    read(11)
  endif
  read(11)
  read(11)
  read(11) pol%ecuts
  read(11) nq,((qpts(jj,iq),jj=1,3),iq=1,nq)
  read(11) ng
  do iq=1,nq
    read(11) ngq,nmtx0_of_q(iq)
    if(nmtx0_of_q(iq) .gt. nmtxmax) then
      nmtxmax=nmtx0_of_q(iq)
    endif
    read(11)
    read(11)
    if (pol%freq_dep/=0) then    
      do jj = 1, nmtx0_of_q(iq)
        do ii = 1, nmtx0_of_q(iq)
          read(11)
        enddo
#ifdef CPLX
        do ii = 1, nmtx0_of_q(iq)
          read(11)
        enddo
#endif
      enddo
    else
      do jj = 1, nmtx0_of_q(iq)
        read(11)
      enddo
    endif
  enddo

! JRD - LETS READ the nmtx here... Putting it in .inp is too much to ask user

  call close_file(11)

  SAFE_ALLOCATE(isort, (ng))  
  SAFE_ALLOCATE(isorti, (ng)) 
  SAFE_ALLOCATE(ekin, (ng))   
                              
! Write header for hdf5 file

! Read body of epsmat from old format and write to hdf5 file

  call open_file(unit=11,file=TRUNC(infile),form='unformatted',status='old')
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11) nq
  read(11) ng,((gcomponents(jj,ig),jj=1,3),ig=1,ng)

  if (any(gcomponents/=mf%gvec%components)) then
    write(0,*) 'ERROR: gspaces between WFN and epsmat files don`t match!'
    call die('gspace mismatch', only_root_writes=.true.)
  endif
  
  ! FHJ: some defaults..
  if (index(outfile,'epsmat')/=0 .or. index(outfile,'eps0mat')/=0) then
    pol%matrix_type = 0
  elseif (index(outfile,'chimat')/=0 .or. index(outfile,'chi0mat')/=0) then
    pol%matrix_type = 2
    call die('chimat files cannot be converted at this time. Sorry!')
  else
    pol%matrix_type = 0
    write(0,*)
    write(0,*) 'WARNING: could not determine if we have a chimat or epsmat file.'
    write(0,*) 'Assuming we are dealing with an epsmat file.'
    write(0,*)
  endif
  call eps_setup_sizes(pol, SCALARSIZE, mf%kp%nspin)
  pol%efermi = 0d0
  pol%intraband_flag = 0
  pol%intraband_overlap_min = 0.9d0
  pol%subsample = .false.
  call eps_hdf5_setup(mf%kp, mf%gvec, mf%syms, mf%crys, pol, qgrid, nq, qpts, &
    nmtx0_of_q, nmtxmax, -1, trim(outfile))
  SAFE_DEALLOCATE_P(pol%dFreqBrd)
  SAFE_DEALLOCATE_P(pol%dFreqGrid)

  do iq=1,nq
    write(6,'(1x,a,i0)') "iq = ", iq
    if(pol%freq_dep/=0) then
      ! FHJ: Note: we no longer store epsADyn, we compute it on-the-fly
      SAFE_ALLOCATE(epsRDyn, (nmtx0_of_q(iq),nmtx0_of_q(iq),pol%nFreq))
    else
      SAFE_ALLOCATE(eps, (nmtx0_of_q(iq),nmtx0_of_q(iq))) 
    endif
    read(11) ngq,nmtx0_of_q_file(iq),(isort(ig),isorti(ig),ig=1,ngq)
    if(nmtx0_of_q_file(iq) .ne. nmtx0_of_q(iq)) then
      call die('nmtx0_of_q in input file not the same as what is in epsmat file')
    endif
    read(11) (ekin(ig),ig=1,ngq)
    read(11) (qk(jj),jj=1,3)
    if (pol%freq_dep/=0) then    
      do jj = 1, nmtx0_of_q(iq)
        do ii = 1, nmtx0_of_q(iq)
          read(11) (epsRDyn(ii, jj, ifreq), ifreq= 1, pol%nFreq)
        enddo
#ifdef CPLX
        ! FHJ: Note: we no longer store epsADyn, we compute it on-the-fly
        do ii = 1, nmtx0_of_q(iq)
          read(11)
        enddo
#endif
      enddo
    else
      do jj = 1, nmtx0_of_q(iq)
         read(11) (eps(ii,jj),ii=1,nmtx0_of_q(iq))
      enddo
    endif
    call write_gvec_indices_hdf(ng, isort, isorti, ekin, iq, trim(outfile))
    call write_vcoul()

    if(pol%freq_dep/=0) then
      call write_matrix_f_ser_hdf(pol%nFreq, epsRDyn, nmtx0_of_q(iq), iq, 1, trim(outfile))
    else
      call write_matrix_ser_hdf(eps, nmtx0_of_q(iq), iq, 1, trim(outfile))
    endif

    SAFE_ALLOCATE(epsdiag,(pol%matrix_flavor,nmtx0_of_q(iq),1))

    epsdiag(:,:,:) = 0d0
    if (pol%freq_dep==0) then
      do jj = 1, nmtx0_of_q(iq)
        epsdiag(1,jj,1) = dble(eps(jj,jj))
#ifdef CPLX
        epsdiag(2,jj,1) = IMAG(eps(jj,jj))
#endif
      enddo
    else
      do jj = 1, nmtx0_of_q(iq)
        epsdiag(1,jj,1) = dble(epsRDyn(jj,jj,1))
        epsdiag(2,jj,1) = IMAG(epsRDyn(jj,jj,1))
      enddo
    endif
    call write_matrix_diagonal_hdf(epsdiag, nmtx0_of_q(iq), iq, pol%matrix_flavor, trim(outfile))

    SAFE_DEALLOCATE(epsdiag)
    if(pol%freq_dep/=0) then
      SAFE_DEALLOCATE(epsRDyn)
    else
      SAFE_DEALLOCATE(eps)
    endif
    call set_qpt_done(trim(outfile), iq)
  enddo
  write(6,*) 'All done!'

  call close_file(11)

contains

  ! FHJ: Compute vcoul and write to epsmat file
  subroutine write_vcoul
    integer :: nmtx
    type (twork_scell) :: work_scell
    real(DP) :: oneoverq, q0(3), qq(3)
    real(DP), allocatable :: vcoul(:)
    SCALAR :: epshead, wcoul0

    PUSH_SUB(write_vcoul)

    epshead = ZERO
    wcoul0 = ZERO
    nmtx = nmtx0_of_q(iq)
    q0(:) = 0d0
    qq(:) = qpts(:,iq)
    SAFE_ALLOCATE(vcoul, (nmtx))
    call vcoul_generator(pol%icutv, pol%truncval, mf%gvec, mf%crys%bdot, mf%crys%celvol, &
      0, nmtx, isort, 0, qq, q0, vcoul, &
      0, 0, 0d0, oneoverq, qgrid, epshead, work_scell, .false., wcoul0)
    call write_vcoul_hdf(vcoul, iq, trim(outfile))
    SAFE_DEALLOCATE(vcoul)

    POP_SUB(write_vcoul)

  end subroutine write_vcoul


end program epsmat_old2hdf5
