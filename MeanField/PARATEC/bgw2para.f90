
!-------------------------------------------------------------------------------
!
!   bgw2para.f90
!   converts BerkeleyGW file WFN to PARATEC files WFN$n.$s and BAND
!   written by M. Jain (May 2010)
!
!-------------------------------------------------------------------------------

#include "f_defs.h"

program bgw2para

  use global_m
  use band_m
  use input_utils_m
  use iolib_m
  use wfn_rho_vxc_io_m

  implicit none

  type(band) :: bands
  type(crystal) :: crys
  type(kpoints) :: kp
  type(symmetry) :: syms
  type(gspace) :: gvec

  character(len=11) :: file_wfn
  character(len=24) :: bdate
  character(len=8) :: btime
  character(len=11) :: cdate
  character(len=14) :: stime
  integer :: iproc

  integer :: i, ib
  integer :: ik, ig, ngk_local, is
  integer ::  nstart_local
  integer, allocatable :: gvec_local(:,:)
  integer, allocatable :: ngk_dist(:), ngk_start(:), ngk_end(:), ngk_t(:)
  complex(DPC), allocatable :: zc(:,:), wfn(:,:,:)
  SCALAR, allocatable :: temp(:,:)
  real(DP), allocatable :: ekin(:)
  integer :: nargs, iflavor
  character(len=80) :: infile, usage
  character(len=3) :: sheader

  usage = 'Usage: bgw2para.x infile'

  call peinfo_init()

  if (peinf%inode == 0) then
    nargs = command_argument_count()

    if (nargs .ne. 1) then
      call die(usage)
    endif

    call get_command_argument(1, infile)

    call open_file(unit=10, file=infile, form='unformatted', status='old')
  endif

  sheader = 'WFN'
  iflavor = 0
  call read_binary_header_type(10, sheader, iflavor, kp, gvec, syms, crys, warn = .false.)

! Assign the band values to the corresponding ones...
  bands%nspin = kp%nspin
  bands%nrk = kp%nrk
  bands%min(1) = kp%mnband
  bands%min(2) = kp%mnband
  bands%max = kp%mnband

  SAFE_ALLOCATE(bands%nband, (kp%nrk, kp%nspin)) 
  bands%nband = 0 
  SAFE_ALLOCATE(bands%occup, (kp%mnband, kp%nrk, kp%nspin)) 
  bands%occup = 0d0 
  SAFE_ALLOCATE(bands%energy, (kp%mnband, kp%nrk, kp%nspin)) 
  bands%energy = 0d0 
  SAFE_ALLOCATE(bands%ekn, (kp%mnband, kp%nrk, kp%nspin)) 
  bands%ekn = 0d0 
  SAFE_ALLOCATE(bands%ifmax, (kp%nrk, kp%nspin)) 
  bands%ifmax = 0 
  do is = 1, kp%nspin
    do ik = 1, kp%nrk
      do i = 1, kp%mnband
         bands%energy(i, ik, is) = kp%el(i, ik, is)
         bands%occup(i, ik, is)  = kp%occ(i, ik, is)*kp%w(ik)*2.0d0
      enddo
      bands%nband(ik, is) = kp%mnband
      bands%ifmax(ik,is) = kp%ifmax(ik,is)
    enddo
  enddo

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  call read_binary_gvectors(10, gvec%ng, gvec%ng, gvec%components)

  SAFE_ALLOCATE(zc, (kp%ngkmax, kp%nspin))
  SAFE_ALLOCATE(ngk_dist, (peinf%npes))
  SAFE_ALLOCATE(ngk_start, (peinf%npes))
  SAFE_ALLOCATE(ngk_end, (peinf%npes))
  SAFE_ALLOCATE(ngk_t, (peinf%npes))

  do ik = 1, kp%nrk
      call read_binary_gvectors(10, kp%ngk(ik), gvec%ng, gvec%components)

    ! Calculate the local dimensions, start and end (wrt to global index)
    if (peinf%inode .lt. mod(kp%ngk(ik), peinf%npes)) then
      ngk_local = kp%ngk(ik)/peinf%npes + 1
      nstart_local = peinf%inode*ngk_local
    else
      ngk_local = kp%ngk(ik)/peinf%npes 
      nstart_local = peinf%inode*ngk_local + mod(kp%ngk(ik),peinf%npes)
    endif
    ngk_dist = 0
    ngk_start = 0
    ngk_end = 0
    ngk_dist(peinf%inode+1) = ngk_local
    ngk_start(peinf%inode+1) = nstart_local
#ifdef MPI
    ngk_t = 0
    call mpi_allreduce(ngk_dist,ngk_t,peinf%npes,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
    ngk_dist = ngk_t
    ngk_t = 0
    call mpi_allreduce(ngk_start,ngk_t,peinf%npes,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
    ngk_start = ngk_t
#endif

    ! Now everyone knows where everyone begins and how many G-vectors they hold

    SAFE_ALLOCATE(gvec_local, (3,ngk_local))
    if (peinf%inode == 0) then
      SAFE_ALLOCATE(ekin, (kp%ngk(ik)))
      call kinetic_energies(gvec, crys%bdot, ekin, qvec = kp%rk(:, ik))
    endif
#ifdef MPI
    call mpi_scatterv(gvec%components(1,:),ngk_dist,ngk_start,MPI_INTEGER,gvec_local(1,:), &
      ngk_local,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_barrier(MPI_COMM_WORLD,mpierr)
    call mpi_scatterv(gvec%components(2,:),ngk_dist,ngk_start,MPI_INTEGER,gvec_local(2,:), &
      ngk_local,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_barrier(MPI_COMM_WORLD,mpierr)
    call mpi_scatterv(gvec%components(3,:),ngk_dist,ngk_start,MPI_INTEGER,gvec_local(3,:), &
      ngk_local,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_barrier(MPI_COMM_WORLD,mpierr)
#else
    gvec_local(:,:) = gvec%components(:,:)
#endif
    
    SAFE_ALLOCATE(wfn, (ngk_local,kp%mnband,kp%nspin))
    SAFE_ALLOCATE(temp, (kp%ngk(ik),kp%nspin))

    do ib = 1, kp%mnband
      call read_binary_data(10, kp%ngk(ik), kp%ngk(ik), kp%nspin, temp, bcast = .false.)
      if (peinf%inode == 0) then
        zc(1:kp%ngk(ik), 1:kp%nspin) = COMPLEXIFY(temp(1:kp%ngk(ik), 1:kp%nspin))
        do is = 1, kp%nspin
          do ig = 1, kp%ngk(ik)
            bands%ekn(ib,ik,is) = bands%ekn(ib,ik,is) + &
              ekin(ig)*real(zc(ig,is)*conjg(zc(ig,is)),dp)
          enddo
        enddo
      endif
#ifdef MPI
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif
      do is = 1, kp%nspin
#ifdef MPI
        call mpi_scatterv(zc(:,is),ngk_dist,ngk_start,MPI_DOUBLE_COMPLEX,wfn(:,ib,is), &
          ngk_local,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
#else
        wfn(:,ib,is)=zc(:,is)
#endif
      enddo
#ifdef MPI
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif
    enddo ! ib

    ! For this k point all the data is now read and stored parallel fashion
    ! Write the corresponding WFN file
    
    do is = 1, kp%nspin
      write(file_wfn, 120) ik, is
      if (peinf%inode == 0) then
        call date_time(cdate, stime)
        write(bdate, '(a11,13x)') cdate
        write(btime, '(a8)') stime(1:8)
        write(6,*) "Opening file : ", file_wfn, " for writing!"
        call open_file(unit=21, file=file_wfn, form = 'unformatted', status = 'unknown')
        rewind(21)
        write(21) bdate, btime
        write(21) kp%mnband, 1
        call close_file(21)
      endif
#ifdef MPI
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif
      do iproc = 0, peinf%npes-1
        if (peinf%inode == iproc) then
          call open_file(unit = 21, file = file_wfn, position = 'append', status = 'unknown', form = 'unformatted')
          
          do ig = 1, ngk_local
            write(21) gvec_local(1,ig),gvec_local(2,ig),gvec_local(3,ig)
            write(21) (wfn(ig,ib,is), ib = 1, kp%mnband)
          enddo
          call close_file(21)
        endif
#ifdef MPI
        call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif
      enddo
#ifdef MPI
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif
      if(peinf%inode == 0) then
        call open_file(unit = 21, file = file_wfn, position = 'append', status = 'unknown', form = 'unformatted')
        write(21) -1234567, 0, 0
        call close_file(21)
      endif
    enddo
    
    SAFE_DEALLOCATE(wfn)
    SAFE_DEALLOCATE(gvec_local)
    if (peinf%inode .eq. 0) then
      SAFE_DEALLOCATE(ekin)
    end if
    SAFE_DEALLOCATE(temp)

  enddo

  call dealloc_header_type(sheader, crys, kp)
  SAFE_DEALLOCATE_P(gvec%components)

  SAFE_DEALLOCATE(zc)
  SAFE_DEALLOCATE(ngk_dist)
  SAFE_DEALLOCATE(ngk_start)
  SAFE_DEALLOCATE(ngk_end)
  SAFE_DEALLOCATE(ngk_t)
  
  if (peinf%inode == 0) then
    call write_band('BAND',bands)
    call close_file(10)
  endif
  
120 format('WFN',i5.5,'.',i1)
  
#ifdef MPI
  call MPI_Finalize(mpierr)
#endif
  
end program bgw2para
