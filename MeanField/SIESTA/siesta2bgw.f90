!===============================================================================
!
! Program:
!
! (1) siesta2bgw         Originally By gsm        Last Modified 3/03/2010 (gsm)
!
! The SIESTA wrapper reads lattice vectors, atomic species and
! positions from SystemLabel.XV file, k-points from SystemLabel.KP
! file, eigenvalues from SystemLabel.EIG file, wavefunctions in R-space
! from SystemLabel{.K$k}.WF$n{.UP|.DOWN}{.REAL|.IMAG}.cube file, symmetry
! group, kinetic energy cutoffs for wavefunctions and charge density,
! k-grid dimensions and offsets from the input file, generates reciprocal
! lattice vectors, rotation matrices, fractional translations, FFT grid and
! G-vectors, performs FFT from R-space to G-space, and writes the selected
! range of bands to the output WFN file. Make sure to set up the real
! space grid in SIESTA/DENCHAR utility same as FFT grid. If the two grids
! differ the SIESTA wrapper will return an error. The output WFN file
! can be used as an auxiliary wavefunction for the SAPO code.
!
! A full GW/BSE calculation with SIESTA wavefunctions also requires RHO and
! VXC files. These files can be generated from SystemLabel.RHO{.UP|.DN}.cube,
! SystemLabel.VT{.UP|.DN}.cube and SystemLabel.VH.cube files obtained with
! SIESTA/GRID2CUBE utility on SIESTA real space grid. Make sure SIESTA real
! space grid is equivalent to FFT grid by setting parameter MeshCutoff
! in SIESTA input file same as kinetic energy cutoff for charge density.
!
! Input is read from file siesta2bgw.inp.
!
!===============================================================================

#include "f_defs.h"

program siesta2bgw

  use global_m
  use check_inversion_m
  use fftw_m
  use fft_parallel_m
  use read_cube_m
  use symmetries_m
  use wfn_rho_vxc_io_m
  use name_m
  use real_m
  use norm_m
  use write_program_header_m
  implicit none

!------------------------
! Allocatable arrays

  integer, pointer :: as(:)
  real(DP), pointer :: ap(:,:), ap_lat(:,:)
  real(DP), allocatable :: kw(:)
  real(DP), allocatable :: kpt(:,:)
  real(DP), allocatable :: enlo(:,:)
  real(DP), allocatable :: en(:,:,:)
  real(DP), allocatable :: occ(:,:,:)
  integer, allocatable :: ifmin(:,:)
  integer, allocatable :: ifmax(:,:)
  integer, allocatable :: gvec(:,:)
  integer, allocatable :: ngk(:)
  integer, allocatable :: idummy(:)
  complex(DPC), allocatable :: wfn(:,:)
  complex(DPC), allocatable :: wfn_d(:,:,:,:)
  integer, allocatable :: isort(:)
  integer, allocatable :: isort_d(:,:)
  complex(DPC), pointer :: fftbox_3D(:,:,:)
  complex(DPC), pointer :: fftbox_2D(:,:,:)
  complex(DPC), allocatable :: fftbox_1D(:,:)
  complex(DPC), allocatable :: buffer_2D(:,:,:)
  complex(DPC), allocatable :: buffer_1D(:,:,:)
  complex(DPC), allocatable :: buffer_rod(:)
  complex(DPC), allocatable :: buffer_sph(:)
  integer, allocatable :: inv_indx(:,:,:)
  integer, allocatable :: crod(:)
  integer, allocatable :: csph(:)
  integer, allocatable :: drod(:)
  integer, allocatable :: dsph(:)
  integer, allocatable :: irod(:,:)
  integer, allocatable :: isph(:)
  complex(DPC), allocatable :: rho(:,:)
  complex(DPC), allocatable :: rho_buf(:)
  complex(DPC), allocatable :: rho_d(:)
  complex(DPC), allocatable :: vxc(:,:)
  complex(DPC), allocatable :: vxc_buf(:)
  complex(DPC), allocatable :: vxc_d(:)
  SCALAR, allocatable :: data_buf(:,:)

  logical :: wfng_ref_flag,wfng_band_flag,wfng_energy_flag, &
    wfng_gamma_real,rhog_flag,vxcg_flag,fparafft
  integer :: cell_symmetry,wfng_nk1,wfng_nk2,wfng_nk3, &
    wfng_ref_kpoint,wfng_ref_spin,wfng_ref_band,wfng_band_min, &
    wfng_band_max,i,j,k,ierr,is,ik,ib,ig,igk,ig0, &
    nsymc,ns,nk,nb,ng,ngkmax,ng_l,ng_g,ngk_l,ngk_g, &
    nbmin,nbmax,nbnum,nbocc,nbstart,nbend,nFFTgridpts,ncount,na, &
    k0,k1,k2,i1,j1,j2,j3,tmpnum,ip,np,Nplane,Nrod,iflavor, &
    nelec,mrod,msph,kgrid(3),FFTgrid(3), &
    rotc(3,3,48),real_or_complex,spacegroup
  real(DP) :: ecutwfn,ecutrho,wfng_dk1,wfng_dk2,wfng_dk3, &
    wfng_ref_energy,wfng_energy_min,wfng_energy_max,efermi,de, &
    x,gcutm,ekin,scale,celvol,recvol,al,bl,relec,rdummy, &
    tsec(2),kvec(3),a(3,3),b(3,3),adot(3,3), &
    bdot(3,3),kshift(3),tauc(3,48)
  character(len=256) :: systemlabel,wfng_output_file, &
    rhog_output_file,vxcg_output_file,siesta_input_file,tmpstr
  character(len=16) :: s1,s2,s3,s4,s5,s6
  character*16 :: routnam(7)
  character :: sheader*3
  character*21 :: symbol

!------------------------
! Parameters

  integer, parameter :: Nfac = 3
  character(len=256), parameter :: input_file = 'siesta2bgw.inp'

!------------------------
! Input namelist

  namelist /input_siesta2bgw/ systemlabel, wfng_output_file, &
    rhog_output_file, vxcg_output_file, ecutwfn, ecutrho, &
    wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3, &
    wfng_ref_flag, wfng_ref_kpoint, wfng_ref_spin, wfng_ref_band, &
    wfng_ref_energy, wfng_band_flag, wfng_band_min, wfng_band_max, &
    wfng_energy_flag, wfng_energy_min, wfng_energy_max, wfng_gamma_real
  
!------------------------
! Initialize MPI

  ierr = 0
  call peinfo_init()
  
#ifdef MPI
  fparafft=.true.
#else
  fparafft=.false.
#endif

#ifdef CPLX
  real_or_complex = 2
#else
  real_or_complex = 1
#endif

!------------------------
! Initialize timer

  call timacc(0,0)
  
  call timacc(1,1)
  
  call timacc(2,1)

!------------------------
! Write header

  call write_program_header('SIESTA2BGW', .false.)

!------------------------
! Read input file

  if (peinf%inode.eq.0) then
    write(6,801) TRUNC(input_file)
    systemlabel = 'prefix'
    wfng_output_file = 'wfng.lo'
    rhog_output_file = 'rhog.lo'
    vxcg_output_file = 'vxcg.lo'
    ecutwfn = 0.0d0
    ecutrho = 0.0d0
    wfng_nk1 = 0
    wfng_nk2 = 0
    wfng_nk3 = 0
    wfng_dk1 = 0.0d0
    wfng_dk2 = 0.0d0
    wfng_dk3 = 0.0d0
    wfng_ref_flag = .false.
    wfng_ref_kpoint = 0
    wfng_ref_spin = 0
    wfng_ref_band = 0
    wfng_ref_energy = 0.0
    wfng_band_flag = .false.
    wfng_band_min = 0
    wfng_band_max = 0
    wfng_energy_flag = .false.
    wfng_energy_min = 0.0d0
    wfng_energy_max = 0.0d0
    wfng_gamma_real = .false.
    call open_file(55,file=input_file,status='old',form='formatted')
    read(55,input_siesta2bgw,iostat=ierr)
    if (ierr.eq.0) call close_file(55)
    if (ecutwfn.lt.TOL_Zero.or. &
      ecutrho.lt.TOL_Zero.or. &
      wfng_nk1.lt.0.or. &
      wfng_nk2.lt.0.or. &
      wfng_nk3.lt.0.or. &
      wfng_ref_kpoint.lt.0.or. &
      wfng_ref_spin.lt.0.or. &
      wfng_ref_band.lt.0.or. &
      wfng_band_min.lt.0.or. &
      wfng_band_max.lt.wfng_band_min.or. &
      wfng_energy_max.lt.(wfng_energy_min-TOL_Zero)) &
      ierr=1
    if (ierr.eq.0) then
      if (wfng_nk1.eq.0.or.wfng_nk2.eq.0.or.wfng_nk3.eq.0) then
        write(0,701)
      endif
      if (wfng_band_flag.and.wfng_energy_flag) then
        write(0,702)
        wfng_energy_flag=.false.
      endif
      kgrid(1)=wfng_nk1
      kgrid(2)=wfng_nk2
      kgrid(3)=wfng_nk3
      kshift(1)=wfng_dk1
      kshift(2)=wfng_dk2
      kshift(3)=wfng_dk3
    endif
  endif
#ifdef MPI
  call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  if (ierr.ne.0) then
    call die('failed to read input file ' // TRUNC(input_file), only_root_writes = .true.)
  endif
#ifdef MPI
  call MPI_Bcast(systemlabel,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_output_file,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(rhog_output_file,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(vxcg_output_file,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(ecutwfn,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(ecutrho,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(kgrid,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(kshift,3,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_ref_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_ref_kpoint,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_ref_spin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_ref_band,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_ref_energy,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_band_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_band_min,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_band_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_energy_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_energy_min,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_energy_max,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_gamma_real,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif

!------------------------
! Initialize flags for output files

  if (len(TRUNC(rhog_output_file)).gt.0) then
    rhog_flag=.true.
  else
    rhog_flag=.false.
  endif
  if (len(TRUNC(vxcg_output_file)).gt.0) then
    vxcg_flag=.true.
  else
    vxcg_flag=.false.
  endif

!------------------------
! Read parameters from SIESTA lattice file

  call xvfilename(systemlabel,siesta_input_file)
  
  if (peinf%inode.eq.0) then
    write(6,802) TRUNC(siesta_input_file)
    a(:,:)=0.0d0
    na=0
    call open_file(7,file=siesta_input_file,status='old',form='formatted')
    do i=1,3
      if (ierr.eq.0) read(7,*,iostat=ierr)(a(j,i),j=1,3)
    enddo
    if (ierr.eq.0) read(7,*,iostat=ierr)na
    if (ierr.eq.0) call close_file(7)
    do i=1,3
      al=sqrt(a(1,i)**2+a(2,i)**2+a(3,i)**2)
      if (al.lt.TOL_Small) ierr=1
    enddo
    if (na.lt.1) ierr=1
  endif
#ifdef MPI
  call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  if (ierr.ne.0) then
    call die('failed to read SIESTA file ' // TRUNC(siesta_input_file), only_root_writes = .true.)
  endif
#ifdef MPI
  call MPI_Bcast(a,9,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(na,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif

  al = sqrt(a(1,1)**2 + a(2,1)**2 + a(3,1)**2)
  a(:,:) = a(:,:) / al
  celvol = a(1,1) * (a(2,2) * a(3,3) - a(2,3) * a(3,2)) - &
    a(2,1) * (a(1,2) * a(3,3) - a(1,3) * a(3,2)) + &
    a(3,1) * (a(1,2) * a(2,3) - a(1,3) * a(2,2))
  bl = 2.0d0 * PI_D / al
  b(1,1) = (a(2,2) * a(3,3) - a(3,2) * a(2,3)) / celvol
  b(2,1) = (a(3,2) * a(1,3) - a(1,2) * a(3,3)) / celvol
  b(3,1) = (a(1,2) * a(2,3) - a(2,2) * a(1,3)) / celvol
  b(1,2) = (a(2,3) * a(3,1) - a(3,3) * a(2,1)) / celvol
  b(2,2) = (a(3,3) * a(1,1) - a(1,3) * a(3,1)) / celvol
  b(3,2) = (a(1,3) * a(2,1) - a(2,3) * a(1,1)) / celvol
  b(1,3) = (a(2,1) * a(3,2) - a(3,1) * a(2,2)) / celvol
  b(2,3) = (a(3,1) * a(1,2) - a(1,1) * a(3,2)) / celvol
  b(3,3) = (a(1,1) * a(2,2) - a(2,1) * a(1,2)) / celvol
  celvol = abs(celvol) * al**3
  recvol = (2.0d0 * PI_D)**3 / celvol
  adot(:,:) = 0.0d0
  do i=1,3
    do j=1,3
      do k=1,3
        adot(j,i) = adot(j,i) + a(k,j) * a(k,i)
      enddo
    enddo
  enddo
  adot(:,:) = adot(:,:) * al**2
  bdot(:,:) = 0.0d0
  do i=1,3
    do j=1,3
      do k=1,3
        bdot(j,i) = bdot(j,i) + b(k,j) * b(k,i)
      enddo
    enddo
  enddo
  bdot(:,:) = bdot(:,:) * bl**2
  
  gcutm=ecutrho/bl**2
  do i=1,3
    FFTgrid(i)=int(2.0d0*sqrt(gcutm)*sqrt(a(1,i)**2+a(2,i)**2+a(3,i)**2))+1
    do while (.not.check_FFT_size(FFTgrid(i),Nfac))
      FFTgrid(i)=FFTgrid(i)+1
    enddo
  enddo
  
!------------------------
! Read SIESTA lattice file

  SAFE_ALLOCATE(as, (na))
  SAFE_ALLOCATE(ap, (3,na))
  SAFE_ALLOCATE(ap_lat, (3,na))
  
  if (peinf%inode.eq.0) then
    call open_file(7,file=siesta_input_file,status='old',form='formatted')
    do i=1,4
      read(7,*)
    enddo
    do i=1,na
      read(7,*)k,as(i),(ap(j,i),j=1,3)
    enddo
    call close_file(7)
  endif
#ifdef MPI
  call MPI_Bcast(as,na,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(ap,3*na,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif

  ap_lat = matmul(transpose(b),ap) / al ! convert to crystal coords 
  call get_symmetries(na, as, ap_lat, a, FFTgrid, cell_symmetry, nsymc, rotc, tauc, spacegroup, symbol)
  ap(:,:) = ap(:,:) / al

!------------------------
! Read parameters from SIESTA kpoint file

  call kpfilename(systemlabel,siesta_input_file)
  
  if (peinf%inode.eq.0) then
    write(6,802) TRUNC(siesta_input_file)
    nk=0
    call open_file(7,file=siesta_input_file,status='old',form='formatted')
    read(7,*,iostat=ierr)nk
    if (ierr.eq.0) call close_file(7)
    if (nk.lt.1) ierr=1
  endif
#ifdef MPI
  call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  if (ierr.ne.0) then
    call die('failed to read SIESTA file ' // TRUNC(siesta_input_file), only_root_writes = .true.)
  endif
#ifdef MPI
  call MPI_Bcast(nk,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  
  SAFE_ALLOCATE(kw, (nk))
  SAFE_ALLOCATE(kpt, (3,nk))

!------------------------
! Read SIESTA kpoint file

  if (peinf%inode.eq.0) then
    call open_file(7,file=siesta_input_file,status='old',form='formatted')
    read(7,*)
    do ik=1,nk
      read(7,*)k,(kpt(j,ik),j=1,3),kw(ik)
    enddo
    call close_file(7)
  endif
#ifdef MPI
  call MPI_Bcast(kpt,3*nk,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(kw,nk,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif

  ! convert k-points from cartesian to crystal
  do ik=1,nk
    kvec(:)=kpt(:,ik)
    do j=1,3
      kpt(j,ik)=(kvec(1)*a(1,j)+kvec(2)*a(2,j)+kvec(3)*a(3,j)) &
        *al/(2.0d0*PI_D)
    enddo
  enddo

  ! normalize weights of k-points
  x=0.0d0
  do ik=1,nk
    x=x+kw(ik)
  enddo
  if (abs(x).gt.TOL_Small) then
    do ik=1,nk
      kw(ik)=kw(ik)/x
    enddo
  endif

!------------------------
! Read parameters from SIESTA eigenvalue file

  call eigfilename(systemlabel,siesta_input_file)
  
  if (peinf%inode.eq.0) then
    write(6,802) TRUNC(siesta_input_file)
    efermi=0.0d0
    nb=0
    ns=0
    call open_file(7,file=siesta_input_file,status='old',form='formatted')
    if (ierr.eq.0) read(7,*,iostat=ierr)efermi
    if (ierr.eq.0) read(7,*,iostat=ierr)nb,ns,i
    if (ierr.eq.0) call close_file(7)
    if (i.ne.nk) ierr=1
  endif
#ifdef MPI
  call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  if (ierr.ne.0) then
    call die('failed to read SIESTA file ' // TRUNC(siesta_input_file), only_root_writes = .true.)
  endif
#ifdef MPI
  call MPI_Bcast(efermi,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(nb,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(ns,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif

!------------------------
! Read SIESTA eigenvalue file

  SAFE_ALLOCATE(enlo, (nb*ns,nk))
  
  if (mod(nb*ns,10).eq.0) then
    tmpnum=(nb*ns)/10
  else
    tmpnum=(nb*ns)/10+1
  endif
  
  if (peinf%inode.eq.0) then
    call open_file(7,file=siesta_input_file,status='old',form='formatted')
    read(7,*)
    read(7,*)
    do i=1,nk
      do k=1,tmpnum
        read(7,103)tmpstr
        k1=10*(k-1)+1
        k2=10*(k-1)+10
        if (k2.gt.nb*ns) k2=nb*ns
        if (k.eq.1) then
          read(tmpstr,*)ik,(enlo(k0,i),k0=k1,k2)
        else
          read(tmpstr,*)(enlo(k0,i),k0=k1,k2)
        endif
      enddo
    enddo
    call close_file(7)
  endif
#ifdef MPI
  call MPI_Bcast(enlo,nb*ns*nk,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif

!------------------------
! Check if real/complex version is appropriate

  call check_inversion(real_or_complex, nsymc, rotc, ns, .true., .true., tnp = tauc)
  
  call timacc(2,2)

  call timacc(3,1)

!------------------------
! Count numbers of G-vectors for rho and V and for wavefunctions
! at each k-point

  ng=0
  do j3=-FFTgrid(3)/2,FFTgrid(3)/2
    if (mod(j3+FFTgrid(3)/2,peinf%npes).ne.peinf%inode) cycle
    kvec(3)=dble(j3)
    do j2=-FFTgrid(2)/2,FFTgrid(2)/2
      kvec(2)=dble(j2)
      do j1=-FFTgrid(1)/2,FFTgrid(1)/2
        kvec(1)=dble(j1)
        ekin=DOT_PRODUCT(kvec,MATMUL(bdot,kvec))
        if (ekin.lt.ecutrho) ng=ng+1
      enddo
    enddo
  enddo
#ifdef MPI
  i=ng
  call MPI_Allreduce(i,ng,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

  SAFE_ALLOCATE(ngk, (nk))
  
  ngk(:)=0
  do ik=1,nk
    do j3=-FFTgrid(3)/2,FFTgrid(3)/2
      if (mod(j3+FFTgrid(3)/2,peinf%npes).ne.peinf%inode) cycle
      kvec(3)=kpt(3,ik)+dble(j3)
      do j2=-FFTgrid(2)/2,FFTgrid(2)/2
        kvec(2)=kpt(2,ik)+dble(j2)
        do j1=-FFTgrid(1)/2,FFTgrid(1)/2
          kvec(1)=kpt(1,ik)+dble(j1)
          ekin=DOT_PRODUCT(kvec,MATMUL(bdot,kvec))
          if (ekin.lt.ecutwfn) ngk(ik)=ngk(ik)+1
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  SAFE_ALLOCATE(idummy, (nk))
  idummy(:)=ngk(:)
  call MPI_Allreduce(idummy,ngk,nk,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
  SAFE_DEALLOCATE(idummy)
#endif
  
  ngkmax=0
  do ik=1,nk
    if (ngk(ik).gt.ngkmax) ngkmax=ngk(ik)
  enddo
  
  call timacc(3,2)
  
  call timacc(2,1)
  
!------------------------
! Distribute G-vectors over processors, set internal variables,
! and print report to stdout

  if (mod(ng,peinf%npes).eq.0) then
    ng_l=ng/peinf%npes
  else
    ng_l=ng/peinf%npes+1
  endif
  ng_g=ng_l*peinf%npes
  
  if (mod(ngkmax,peinf%npes).eq.0) then
    ngk_l=ngkmax/peinf%npes
  else
    ngk_l=ngkmax/peinf%npes+1
  endif
  ngk_g=ngk_l*peinf%npes
  
  nFFTgridpts=FFTgrid(1)*FFTgrid(2)*FFTgrid(3)
  
  if (wfng_ref_flag.and.wfng_ref_band.ge.1.and.wfng_ref_band.le.nb) then
    i=wfng_ref_kpoint
    k=wfng_ref_band+(wfng_ref_spin-1)*nb
    if (i.ge.1.and.i.le.nk.and.k.ge.1.and.k.le.nb*ns) then
      de=wfng_ref_energy-enlo(k,i)
      do i1=1,nk
        do k1=1,nb*ns
          enlo(k1,i1)=enlo(k1,i1)+de
        enddo
      enddo
    endif
  endif
  
  nbmin=1
  nbmax=nb
  
  if (wfng_band_flag) then
    nbmin=wfng_band_min
    nbmax=wfng_band_max
    if (nbmin.eq.0) nbmin=1
    if (nbmax.eq.0) nbmax=nb
    if (nbmin.gt.nb) nbmin=nb
    if (nbmax.gt.nb) nbmax=nb
  endif
  
  if (wfng_energy_flag) then
    nbmin=nb+1
    do i=1,nk
      do j=1,ns
        do k=1,nb
          if (enlo(k+(j-1)*nb,i).gt.wfng_energy_min- &
            TOL_Small.and.nbmin.gt.k) nbmin=k
        enddo
      enddo
    enddo
    if (nbmin.gt.nb) nbmin=nb
    nbmax=0
    do i=1,nk
      do j=1,ns
        do k=1,nb
          if (enlo(k+(j-1)*nb,i).lt.wfng_energy_max+ &
            TOL_Small.and.nbmax.lt.k) nbmax=k
        enddo
      enddo
    enddo
    if (nbmax.lt.1) nbmax=1
    if (nbmax.lt.nbmin) then
      ib=(nbmin+nbmax)/2
      nbmin=ib
      nbmax=ib
    endif
  endif
  
  nbnum=nbmax-nbmin+1

  if (peinf%inode.eq.0) then
    
    write(s1,101)nsymc
    write(s2,101)nk
    write(s3,101)ns
    write(s4,101)nb
    write(s5,101)ng
    write(s6,101)ngkmax
    write(6,110)TRUNC(s1),TRUNC(s2),TRUNC(s3),TRUNC(s4),TRUNC(s5),TRUNC(s6)
    write(s1,101)ng_l
    write(s2,101)ng_g
    write(s3,101)ngk_l
    write(s4,101)ngk_g
    write(6,120)TRUNC(s1),TRUNC(s2),TRUNC(s3),TRUNC(s4)
    write(s1,101)FFTgrid(1)
    write(s2,101)FFTgrid(2)
    write(s3,101)FFTgrid(3)
    write(s4,101)nFFTgridpts
    write(6,130)TRUNC(s1),TRUNC(s2),TRUNC(s3),TRUNC(s4)
    write(s1,101)nbmin
    write(s2,101)nbmax
    write(s3,101)nbnum
    write(6,140)TRUNC(s1),TRUNC(s2),TRUNC(s3)

  endif

!------------------------
! Store SIESTA eigenvalues for output

  SAFE_ALLOCATE(en, (nbnum,nk,ns))
  SAFE_ALLOCATE(occ, (nbnum,nk,ns))
  SAFE_ALLOCATE(ifmin, (nk,ns))
  SAFE_ALLOCATE(ifmax, (nk,ns))
  
  do ik=1,nk
    do is=1,ns
      do ib=1,nbnum
        en(ib,ik,is)=enlo(nbmin+ib-1+(is-1)*nb,ik)/RYD
      enddo
    enddo
  enddo
  efermi=efermi/RYD
  
  do is=1,ns
    do ik=1,nk
      nbocc=0
      do ib=1,nbnum
        if (en(ib,ik,is).lt.efermi) then
          nbocc=nbocc+1
          occ(ib,ik,is)=1.0d0
        else
          occ(ib,ik,is)=0.0d0
        endif
      enddo
      if (nbmin.le.nbocc) then
        ifmin(ik,is)=1
        ifmax(ik,is)=nbocc
      else
        ifmin(ik,is)=0
        ifmax(ik,is)=0
      endif
    enddo
  enddo
  
  SAFE_DEALLOCATE(enlo)

!------------------------
! Reset output files

  if (peinf%inode.eq.0) then
    call open_file(8,file=wfng_output_file,status='replace',form='unformatted')
    if (ierr.eq.0) call close_file(8)
  endif
#ifdef MPI
  call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  if (ierr.ne.0) then
    call die('failed to write wavefunction file ' // TRUNC(wfng_output_file), only_root_writes = .true.)
  endif
  
  if (rhog_flag) then
    if (peinf%inode.eq.0) then
      call open_file(9,file=rhog_output_file,status='replace',form='unformatted')
      if (ierr.eq.0) call close_file(9)
    endif
#ifdef MPI
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
    if (ierr.ne.0) then
      call die('failed to write charge density file ' // TRUNC(rhog_output_file), only_root_writes = .true.)
    endif
  endif
  
  if (vxcg_flag) then
    if (peinf%inode.eq.0) then
      call open_file(10,file=vxcg_output_file,status='replace',form='unformatted')
      if (ierr.eq.0) call close_file(10)
    endif
#ifdef MPI
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
    if (ierr.ne.0) then
      call die('failed to write exchange-correlation potential file ' // &
        TRUNC(vxcg_output_file), only_root_writes = .true.)
    endif
  endif
  
  call timacc(2,2)
  
  call timacc(3,1)

!------------------------
! Construct list of G-vectors for rho and V and indices
! of G-vectors for wavefunctions at each k-point
  
  SAFE_ALLOCATE(gvec, (3,ng))

  ig=0
  do j3=-FFTgrid(3)/2,FFTgrid(3)/2
    kvec(3)=dble(j3)
    do j2=-FFTgrid(2)/2,FFTgrid(2)/2
      kvec(2)=dble(j2)
      do j1=-FFTgrid(1)/2,FFTgrid(1)/2
        kvec(1)=dble(j1)
        ekin=DOT_PRODUCT(kvec,MATMUL(bdot,kvec))
        if (ekin.lt.ecutrho) then
          ig=ig+1
          gvec(1,ig)=j1
          gvec(2,ig)=j2
          gvec(3,ig)=j3
        endif
      enddo
    enddo
  enddo
  
  SAFE_ALLOCATE(isort, (ngk_g))
  SAFE_ALLOCATE(isort_d, (ngk_l,nk))
  
  do ik=1,nk
    isort(:)=0
    igk=0
    do ig=1,ng
      do j=1,3
        kvec(j)=kpt(j,ik)+dble(gvec(j,ig))
      enddo
      ekin=DOT_PRODUCT(kvec,MATMUL(bdot,kvec))
      if (ekin.lt.ecutwfn) then
        igk=igk+1
        isort(igk)=ig
      endif
    enddo
#ifdef MPI
    call MPI_Barrier(MPI_COMM_WORLD,mpierr)
    call MPI_Scatter(isort,ngk_l,MPI_INTEGER,isort_d(:,ik),ngk_l, &
      MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#else
    isort_d(:,ik)=isort(:)
#endif
  enddo
  
  call timacc(3,2)
  
  call timacc(2,1)

!------------------------
! Allocate arrays for wavefunctions

  SAFE_ALLOCATE(wfn, (ngk_g,ns))
  SAFE_ALLOCATE(wfn_d, (ngk_l,nbnum,ns,nk))

!------------------------
! Open files

  if (peinf%inode.eq.0) then
    call open_file(8,file=wfng_output_file,status='old',form='unformatted')
  endif
  
  if (rhog_flag.and.peinf%inode.eq.0) then
    call open_file(9,file=rhog_output_file,status='old',form='unformatted')
  endif
  
  if (vxcg_flag.and.peinf%inode.eq.0) then
    call open_file(10,file=vxcg_output_file,status='old',form='unformatted')
  endif
  
!------------------------
! Initialize SIESTA wavefunction file

  if (nk.eq.1.and.wfng_gamma_real) then
    ! .cube files from siesta-2/denchar
    np=1
  else
    ! .REAL.cube & .IMAG.cube files from siesta-3/denchar
    np=2
  endif

!------------------------
! Initialize FFT

  if (fparafft) then
    call fft_set_p(FFTgrid,Nplane,Nrod)
  endif
  
  if (fparafft) then
    SAFE_ALLOCATE(crod, (peinf%npes))
    SAFE_ALLOCATE(csph, (peinf%npes))
    SAFE_ALLOCATE(drod, (peinf%npes))
    SAFE_ALLOCATE(dsph, (peinf%npes))
    SAFE_ALLOCATE(irod, (1,1))
    SAFE_ALLOCATE(isph, (1))
    SAFE_ALLOCATE(buffer_rod, (1))
    SAFE_ALLOCATE(buffer_sph, (1))
    SAFE_ALLOCATE(fftbox_2D, (FFTgrid(1),FFTgrid(2),Nplane))
    SAFE_ALLOCATE(fftbox_1D, (FFTgrid(3),Nrod))
    SAFE_ALLOCATE(buffer_2D, (Nrod,Nplane,peinf%npes))
    SAFE_ALLOCATE(buffer_1D, (Nrod,Nplane,peinf%npes))
  else
    SAFE_ALLOCATE(inv_indx, (FFTgrid(1),FFTgrid(2),FFTgrid(3)))
    SAFE_ALLOCATE(fftbox_3D, (FFTgrid(1),FFTgrid(2),FFTgrid(3)))
  endif

  scale=sqrt(celvol)/dble(FFTgrid(1)*FFTgrid(2)*FFTgrid(3))
  
  call timacc(2,2)

!------------------------
! Loop over SIESTA wavefunction files

  do ik=1,nk
    
    call timacc(3,1)

#ifdef MPI
    call MPI_Barrier(MPI_COMM_WORLD,mpierr)
    call MPI_Gather(isort_d(:,ik),ngk_l,MPI_INTEGER, &
      isort,ngk_l,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(isort,ngk_g,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#else
    isort(:)=isort_d(:,ik)
#endif

!------------------------
! Invert indices of G-vectors for current k-point

    if (fparafft) then
      call fft_map_p(0,ngk(ik),ngk_l,Nrod,FFTgrid,isort,gvec, &
        mrod,msph,crod,csph,drod,dsph,irod,isph)
      SAFE_DEALLOCATE(irod)
      SAFE_DEALLOCATE(isph)
      SAFE_DEALLOCATE(buffer_rod)
      SAFE_DEALLOCATE(buffer_sph)
      SAFE_ALLOCATE(irod, (2,MAX(1,mrod)))
      SAFE_ALLOCATE(isph, (MAX(1,msph)))
      SAFE_ALLOCATE(buffer_rod, (MAX(1,mrod)))
      SAFE_ALLOCATE(buffer_sph, (MAX(1,msph)))
      call fft_map_p(1,ngk(ik),ngk_l,Nrod,FFTgrid,isort,gvec, &
        mrod,msph,crod,csph,drod,dsph,irod,isph)
    else ! fparafft
      call fft_map_s(ngk(ik),FFTgrid,isort,gvec,inv_indx)
    endif ! fparafft
    
    call timacc(3,2)
    
    do is=1,ns
      do ib=1,nbnum
        
        call timacc(4,1)
        
        if (fparafft) then
          fftbox_2D(:,:,:)=(0.0d0,0.0d0)
        else
          fftbox_3D(:,:,:)=(0.0d0,0.0d0)
        endif
        
        do ip=1,np
          
          call wfnrfilename(ik,is,nbmin+ib-1,ip,ns,np,systemlabel,siesta_input_file)

!------------------------
! Read SIESTA wavefunction file

          if (peinf%inode.eq.0) write(6,803) TRUNC(siesta_input_file)
          
          call read_cube(fparafft,7,siesta_input_file,a,al,ip, &
            Nplane,FFTgrid,fftbox_3D,fftbox_2D,ierr)
          
          if (ierr.ne.0) then
            call die('failed to read SIESTA wavefunction file ' // TRUNC(siesta_input_file))
          endif
          
        enddo ! ip
        
#ifdef MPI
        if (.not.fparafft) then
          call MPI_Bcast(fftbox_3D(1,1,1),FFTgrid(1)*FFTgrid(2)*FFTgrid(3),MPI_COMPLEX_DPC,0,MPI_COMM_WORLD,mpierr)
        endif
#endif
        
        call timacc(4,2)
        
        call timacc(5,1)

!------------------------
! Perform FFT from R-space to G-space

        if (fparafft) then
          if (peinf%inode.eq.0) write(6,804)
          call fft_r2g_p(FFTgrid,Nplane,Nrod,mrod,msph,crod, &
            csph,drod,dsph,irod,isph,scale,fftbox_2D,fftbox_1D, &
            buffer_2D,buffer_1D,buffer_rod,buffer_sph,wfn_d(:,ib,is,ik))
        else
          if (peinf%inode.eq.0) write(6,805)
          call fft_r2g_s(FFTgrid,scale,inv_indx,fftbox_3D,wfn(:,is))
#ifdef MPI
          call MPI_Barrier(MPI_COMM_WORLD,mpierr)
          call MPI_Scatter(wfn(:,is),ngk_l,MPI_COMPLEX_DPC, &
            wfn_d(:,ib,is,ik),ngk_l,MPI_COMPLEX_DPC, &
            0,MPI_COMM_WORLD,mpierr)
#else
          wfn_d(:,ib,is,ik)=wfn(:,is)
#endif
        endif

        call timacc(5,2)
        
      enddo ! ib
    enddo ! is
    
  enddo ! ik
  
  if (peinf%inode.eq.0) write(6,*)
  
  call timacc(6,1)

!------------------------
! Normalize wavefunctions in G-space (complex version)
! Construct real wavefunctions in G-space (real version)

#ifdef CPLX
  if (peinf%inode.eq.0) write(6,806)
  nbstart=1
  nbend=nbnum
  call norm_wfng(ngk_l,nbstart,nbend,nbnum,ns,nk,wfn_d)
#else
  if (peinf%inode.eq.0) write(6,807)
  call real_wfng(ngk_l,nbnum,ns,nk,en,wfn_d)
#endif

  call timacc(6,2)
  
  call timacc(7,1)

!------------------------
! Write header of wavefunction file

  if (peinf%inode.eq.0) then
    write(6,808) TRUNC(wfng_output_file)
    sheader = 'WFN'
    iflavor = 0
    call write_binary_header(8, sheader, iflavor, ns, ng, nsymc, &
      cell_symmetry, na, nk, nbnum, ngkmax, ecutrho, ecutwfn, FFTgrid, &
      kgrid, kshift, celvol, al, a, adot, recvol, bl, b, bdot, &
      rotc, tauc, as, ap, ngk, kw, kpt, ifmin, ifmax, en, occ)
    call write_binary_gvectors(8, ng, ng, gvec)
  endif
  
!------------------------
! Loop over k-points

  SAFE_ALLOCATE(data_buf, (ngk_g, ns))
  
  do ik=1,nk

!------------------------
! Gather indices of g-vectors for current k-point
! Write list of g-vectors to wavefunction file

#ifdef MPI
    call MPI_Barrier(MPI_COMM_WORLD,mpierr)
    call MPI_Gather(isort_d(:,ik),ngk_l,MPI_INTEGER, &
      isort,ngk_l,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#else
    isort(:)=isort_d(:,ik)
#endif
    call write_binary_gvectors(8, ngk(ik), ng, gvec, gindex = isort)

!------------------------
! Gather band and spin wavefunctions for current k-point
! Write them to wavefunction file

    do ib=1,nbnum
      do is=1,ns
#ifdef MPI
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
        call MPI_Gather(wfn_d(:,ib,is,ik),ngk_l, &
          MPI_COMPLEX_DPC,wfn(:,is),ngk_l, &
          MPI_COMPLEX_DPC,0,MPI_COMM_WORLD,mpierr)
#else
        wfn(:,is)=wfn_d(:,ib,is,ik)
#endif
      enddo ! is
      if (peinf%inode.eq.0) then
        do is=1,ns
          do ig=1,ngk_g
            data_buf(ig,is) = SCALARIFY(wfn(ig,is))
          enddo
        enddo
      endif
      call write_binary_data(8, ngk(ik), ngk_g, ns, data_buf)
    enddo ! ib
    
  enddo ! ik
  
  SAFE_DEALLOCATE(data_buf)
  
  call timacc(7,2)
  
  call timacc(2,1)
  
!------------------------
! Deallocate unneeded arrays and reallocate isort

  SAFE_DEALLOCATE(wfn)
  SAFE_DEALLOCATE(wfn_d)
  SAFE_DEALLOCATE(isort_d)
  
  if (rhog_flag.or.vxcg_flag) then
    SAFE_DEALLOCATE(isort)
    SAFE_ALLOCATE(isort, (ng))
  endif

!------------------------
! Find index of G=(0,0,0) for rho(G=(0,0,0)) and Vxc(G=(0,0,0))

  if (rhog_flag.or.vxcg_flag) then
    ierr=0
    ig0=0
    do ig=1,ng
      if (gvec(1,ig).eq.0.and.gvec(2,ig).eq.0.and.gvec(3,ig).eq.0) ig0=ig
    enddo
    if (ig0.lt.1.or.ig0.gt.ng) ierr=1
    if (ierr.ne.0) then
      call die('failed to find zero G-vector')
    endif
  endif
  
  call timacc(2,2)
  
  call timacc(3,1)

!------------------------
! Invert indices of G-vectors for charge density and potential

  if (rhog_flag.or.vxcg_flag) then
    isort(:)=0
    do ig=1,ng
      isort(ig)=ig
    enddo ! ig
    if (fparafft) then
      call fft_map_p(0,ng,ng_l,Nrod,FFTgrid,isort,gvec, &
        mrod,msph,crod,csph,drod,dsph,irod,isph)
      SAFE_DEALLOCATE(irod)
      SAFE_DEALLOCATE(isph)
      SAFE_DEALLOCATE(buffer_rod)
      SAFE_DEALLOCATE(buffer_sph)
      SAFE_ALLOCATE(irod, (2,MAX(1,mrod)))
      SAFE_ALLOCATE(isph, (MAX(1,msph)))
      SAFE_ALLOCATE(buffer_rod, (MAX(1,mrod)))
      SAFE_ALLOCATE(buffer_sph, (MAX(1,msph)))
      call fft_map_p(1,ng,ng_l,Nrod,FFTgrid,isort,gvec, &
        mrod,msph,crod,csph,drod,dsph,irod,isph)
    else ! fparafft
      call fft_map_s(ng,FFTgrid,isort,gvec,inv_indx)
    endif ! fparafft
  endif ! rhog_flag.or.vxcg_flag
  
  call timacc(3,2)

!------------------------
! Loop over SIESTA charge density files

  if (rhog_flag) then
    
    call timacc(2,1)
    
    scale=celvol/dble(FFTgrid(1)*FFTgrid(2)*FFTgrid(3))
    
    SAFE_ALLOCATE(rho, (ng,ns))
    if (fparafft) then
      SAFE_ALLOCATE(rho_buf, (ng_g))
      SAFE_ALLOCATE(rho_d, (ng_l))
    endif
    
    call timacc(2,2)
    
    do is=1,ns
      
      call timacc(4,1)
      
      if (fparafft) then
        fftbox_2D(:,:,:)=(0.0d0,0.0d0)
      else
        fftbox_3D(:,:,:)=(0.0d0,0.0d0)
      endif
      
      call rhorfilename(is,ns,systemlabel,siesta_input_file)

!------------------------
! Read SIESTA charge density file
      
      if (peinf%inode.eq.0) write(6,809) TRUNC(siesta_input_file)
      
      call read_cube(fparafft,7,siesta_input_file,a,al,1, &
        Nplane,FFTgrid,fftbox_3D,fftbox_2D,ierr)

      if (ierr.ne.0) then
        call die('failed to read SIESTA charge-density file ' // TRUNC(siesta_input_file))
      endif
      
#ifdef MPI
      if (.not.fparafft) then
        call MPI_Bcast(fftbox_3D(1,1,1),FFTgrid(1)*FFTgrid(2)*FFTgrid(3),MPI_COMPLEX_DPC,0,MPI_COMM_WORLD,mpierr)
      endif
#endif
      
      call timacc(4,2)
      
      call timacc(5,1)

!------------------------
! Perform FFT from R-space to G-space

      if (fparafft) then
        if (peinf%inode.eq.0) write(6,804)
        call fft_r2g_p(FFTgrid,Nplane,Nrod,mrod,msph,crod, &
          csph,drod,dsph,irod,isph,scale,fftbox_2D,fftbox_1D, &
          buffer_2D,buffer_1D,buffer_rod,buffer_sph,rho_d)
#ifdef MPI
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
        call MPI_Gather(rho_d,ng_l,MPI_COMPLEX_DPC, &
          rho_buf,ng_l,MPI_COMPLEX_DPC,0,MPI_COMM_WORLD,mpierr)
#else
        rho_buf(:)=rho_d(:)
#endif
        if (peinf%inode.eq.0) then
          do ig=1,ng
            rho(ig,is)=rho_buf(ig)
          enddo
        endif
      else
        if (peinf%inode.eq.0) write(6,805)
        call fft_r2g_s(FFTgrid,scale,inv_indx,fftbox_3D,rho(:,is))
      endif
      
      call timacc(5,2)
      
    enddo ! is
    
    if (peinf%inode.eq.0) write(6,*)
    
    call timacc(7,1)

!------------------------
! Normalize charge density in G-space
! Print rho(G = 0) to stdout
    
    if (peinf%inode.eq.0) then
      relec=0.0d0
      do is=1,ns
        relec=relec+dble(rho(ig0,is))
      enddo
      nelec=nint(relec)
      rdummy=abs(relec-dble(nelec))
      if (nelec.ne.0.and.rdummy.gt.TOL_Zero) then
        rdummy=dble(nelec)/relec
        do is=1,ns
          do ig=1,ng
            rho(ig,is)=rho(ig,is)*rdummy
          enddo
        enddo
        write(s1,102)relec
        write(s2,102)dble(nelec)
        write(6,210)TRUNC(s1),TRUNC(s2)
      else
        write(s1,102)relec
        write(6,211)TRUNC(s1)
      endif
    endif

!------------------------
! Write charge density file

    SAFE_ALLOCATE(data_buf, (ng, ns))
    
    if (peinf%inode.eq.0) then
      write(6,810) TRUNC(rhog_output_file)
      sheader = 'RHO'
      iflavor = 0
      call write_binary_header(9, sheader, iflavor, ns, ng, nsymc, &
        cell_symmetry, na, nk, nbnum, ngkmax, ecutrho, ecutwfn, FFTgrid, &
        kgrid, kshift, celvol, al, a, adot, recvol, bl, b, bdot, &
        rotc, tauc, as, ap, ngk, kw, kpt, ifmin, ifmax, en, occ)
      call write_binary_gvectors(9, ng, ng, gvec)
      do is=1,ns
        do ig=1,ng
          data_buf(ig,is) = SCALARIFY(rho(ig,is))
        enddo
      enddo
      call write_binary_data(9, ng, ng, ns, data_buf)
    endif
    
    SAFE_DEALLOCATE(data_buf)
    SAFE_DEALLOCATE(rho)
    if (fparafft) then
      SAFE_DEALLOCATE(rho_buf)
      SAFE_DEALLOCATE(rho_d)
    endif
    
    call timacc(7,2)
    
  endif ! rhog_flag

!------------------------
! Loop over SIESTA total and electrostatic potential files

  if (vxcg_flag) then
    
    call timacc(2,1)
    
    scale=1.0d0/dble(FFTgrid(1)*FFTgrid(2)*FFTgrid(3))
    
    SAFE_ALLOCATE(vxc, (ng,ns))
    if (fparafft) then
      SAFE_ALLOCATE(vxc_buf, (ng_g))
      SAFE_ALLOCATE(vxc_d, (ng_l))
    endif
    
    call timacc(2,2)
    
    do is=1,ns
      
      call timacc(4,1)
      
      if (fparafft) then
        fftbox_2D(:,:,:)=(0.0d0,0.0d0)
      else
        fftbox_3D(:,:,:)=(0.0d0,0.0d0)
      endif
      
      call vtrfilename(is,ns,systemlabel,siesta_input_file)

!------------------------
! Read SIESTA total potential file

      if (peinf%inode.eq.0) write(6,811) TRUNC(siesta_input_file)
      
      call read_cube(fparafft,7,siesta_input_file,a,al,1, &
        Nplane,FFTgrid,fftbox_3D,fftbox_2D,ierr)
      
      if (ierr.ne.0) then
        call die('failed to read SIESTA total potential file ' // TRUNC(siesta_input_file))
      endif
      
      call vhrfilename(systemlabel,siesta_input_file)

!------------------------
! Read SIESTA electrostatic potential file

      if (peinf%inode.eq.0) write(6,812) TRUNC(siesta_input_file)
      
      call read_cube(fparafft,7,siesta_input_file,a,al,2, &
        Nplane,FFTgrid,fftbox_3D,fftbox_2D,ierr)
      
      if (ierr.ne.0) then
        call die('failed to read SIESTA electrostatic potential file ' // TRUNC(siesta_input_file))
      endif
      
#ifdef MPI
      if (.not.fparafft) then
        call MPI_Bcast(fftbox_3D(1,1,1),FFTgrid(1)*FFTgrid(2)*FFTgrid(3),MPI_COMPLEX_DPC,0,MPI_COMM_WORLD,mpierr)
      endif
#endif
      
      call timacc(4,2)
      
      call timacc(2,1)

!------------------------
! Compute exchange-correlation potential in R-space

      if (fparafft) then
        do k=1,Nplane
          do j=1,FFTgrid(2)
            do i=1,FFTgrid(1)
              fftbox_2D(i,j,k)=CMPLX(dble(fftbox_2D(i,j,k))-IMAG(fftbox_2D(i,j,k)),0.0d0)
            enddo
          enddo
        enddo
      else
        do k=1,FFTgrid(3)
          do j=1,FFTgrid(2)
            do i=1,FFTgrid(1)
              fftbox_3D(i,j,k)=CMPLX(dble(fftbox_3D(i,j,k))-IMAG(fftbox_3D(i,j,k)),0.0d0)
            enddo
          enddo
        enddo
      endif
      
      call timacc(2,2)
      
      call timacc(5,1)
      
!------------------------
! Perform FFT from R-space to G-space

      if (fparafft) then
        if (peinf%inode.eq.0) write(6,804)
        call fft_r2g_p(FFTgrid,Nplane,Nrod,mrod,msph,crod, &
          csph,drod,dsph,irod,isph,scale,fftbox_2D,fftbox_1D, &
          buffer_2D,buffer_1D,buffer_rod,buffer_sph,vxc_d)
#ifdef MPI
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
        call MPI_Gather(vxc_d,ng_l,MPI_COMPLEX_DPC, &
          vxc_buf,ng_l,MPI_COMPLEX_DPC,0,MPI_COMM_WORLD,mpierr)
#else
        vxc_buf(:)=vxc_d(:)
#endif
        if (peinf%inode.eq.0) then
          do ig=1,ng
            vxc(ig,is)=vxc_buf(ig)
          enddo
        endif
      else
        if (peinf%inode.eq.0) write(6,805)
        call fft_r2g_s(FFTgrid,scale,inv_indx,fftbox_3D,vxc(:,is))
      endif
      
      call timacc(5,2)
      
    enddo ! is
    
    if (peinf%inode.eq.0) write(6,*)
    
    call timacc(7,1)

!------------------------
! Print Vxc(G = 0) to stdout

    if (peinf%inode.eq.0) then
      rdummy=0.0d0
      do is=1,ns
        rdummy=rdummy+dble(vxc(ig0,is))
      enddo
      rdummy=rdummy*RYD/dble(ns)
      write(s1,102)rdummy
      write(6,220)TRUNC(s1)
    endif

!------------------------
! Write exchange-correlation potential file

    SAFE_ALLOCATE(data_buf, (ng, ns))
    
    if (peinf%inode.eq.0) then
      write(6,813) TRUNC(vxcg_output_file)
      sheader = 'VXC'
      iflavor = 0
      call write_binary_header(10, sheader, iflavor, ns, ng, nsymc, &
        cell_symmetry, na, nk, nbnum, ngkmax, ecutrho, ecutwfn, FFTgrid, &
        kgrid, kshift, celvol, al, a, adot, recvol, bl, b, bdot, &
        rotc, tauc, as, ap, ngk, kw, kpt, ifmin, ifmax, en, occ)
      call write_binary_gvectors(10, ng, ng, gvec)
      do is=1,ns
        do ig=1,ng
          data_buf(ig,is) = SCALARIFY(vxc(ig,is))
        enddo
      enddo
      call write_binary_data(10, ng, ng, ns, data_buf)
    endif
    
    SAFE_DEALLOCATE(data_buf)
    SAFE_DEALLOCATE(vxc)
    if (fparafft) then
      SAFE_DEALLOCATE(vxc_buf)
      SAFE_DEALLOCATE(vxc_d)
    endif
    
    call timacc(7,2)
    
  endif ! vxcg_flag
  
  call timacc(2,1)
  
!------------------------
! Done FFT

  if (fparafft) then
    SAFE_DEALLOCATE(crod)
    SAFE_DEALLOCATE(csph)
    SAFE_DEALLOCATE(drod)
    SAFE_DEALLOCATE(dsph)
    SAFE_DEALLOCATE(irod)
    SAFE_DEALLOCATE(isph)
    SAFE_DEALLOCATE(buffer_rod)
    SAFE_DEALLOCATE(buffer_sph)
    SAFE_DEALLOCATE_P(fftbox_2D)
    SAFE_DEALLOCATE(fftbox_1D)
    SAFE_DEALLOCATE(buffer_2D)
    SAFE_DEALLOCATE(buffer_1D)
  else
    SAFE_DEALLOCATE(inv_indx)
    SAFE_DEALLOCATE_P(fftbox_3D)
  endif

  call destroy_fftw_plans()

!------------------------
! Close files

  if (peinf%inode.eq.0) then
    call close_file(8)
  endif
  
  if (rhog_flag.and.peinf%inode.eq.0) then
    call close_file(9)
  endif
  
  if (vxcg_flag.and.peinf%inode.eq.0) then
    call close_file(10)
  endif

!------------------------
! Deallocate remaining arrays

  SAFE_DEALLOCATE(kw)
  SAFE_DEALLOCATE(kpt)
  SAFE_DEALLOCATE(en)
  SAFE_DEALLOCATE(occ)
  SAFE_DEALLOCATE(ifmin)
  SAFE_DEALLOCATE(ifmax)
  SAFE_DEALLOCATE(ngk)
  SAFE_DEALLOCATE_P(as)
  SAFE_DEALLOCATE_P(ap)
  SAFE_DEALLOCATE_P(ap_lat)
  SAFE_DEALLOCATE(gvec)
  SAFE_DEALLOCATE(isort)
  
  call timacc(2,2)
  
  call timacc(1,2)

!------------------------
! Time Accounting

  routnam(1)='TOTAL:'
  routnam(2)='INIT:'
  routnam(3)='G-SPACE:'
  routnam(4)='INPUT:'
  routnam(5)='FFT:'
#ifdef CPLX
  routnam(6)='WFN_NORM:'
#else
  routnam(6)='WFN_REAL:'
#endif
  routnam(7)='OUTPUT:'
  
  if(peinf%inode.eq.0) then
    write(6,9000)
    do i=2,7
      call timacc(i,3,tsec,ncount)
      write(6,9001) routnam(i),tsec(1),tsec(2),ncount
    enddo
    write(6,*)
    call timacc(1,3,tsec,ncount)
    write(6,9002) routnam(1),tsec(1),tsec(2)
    write(6,*)
9000 format(/,22x,"CPU (s)",8x,"WALL (s)",11x,"#",/)
9001 format(a16,f13.3,3x,f13.3,3x,i9)
9002 format(a16,f13.3,3x,f13.3)
  endif

!------------------------
! Finish

#ifdef MPI
  call MPI_Finalize(mpierr)
#endif

101 format(i16)
102 format(f16.6)
103 format(a)

110 format(1x,'nsym =',1x,a,2x,'nk =',1x,a,2x,'ns =',1x,a,2x,'nb =', &
      1x,a,2x,'ng =',1x,a,2x,'ngkmax =',1x,a)
120 format(1x,'ng_l =',1x,a,2x,'ng_g =',1x,a,2x,'ngk_l =',1x,a,2x, &
      'ngk_g =',1x,a)
130 format(1x,'FFTgrid =',1x,'(',1x,a,1x,a,1x,a,1x,')',2x,'nFFTgridpts =',1x,a)
140 format(1x,'nbmin =',1x,a,2x,'nbmax =',1x,a,2x,'nbnum =',1x,a,/)
  
210 format(1x,'rho(G = 0) =',1x,a,1x,'renormalized to',1x,a,/)
211 format(1x,'rho(G = 0) =',1x,a,/)
220 format(1x,'Vxc(G = 0) =',1x,a,1x,'eV'/)
  
701 format(/,1x,'WARNING: kgrid is set to zero in the wavefunction file', &
      /,10x,'this could lead to problems in Sigma and absorption',/)
702 format(/,1x,'WARNING: an energy window is incompatible with a range of bands',/)
  
801 format(1x,'Reading parameters from file',1x,a,/)
802 format(1x,'Reading SIESTA file',1x,a,/)
803 format(1x,'Reading SIESTA wavefunction from file',1x,a)
804 format(1x,'Performing parallel FFT from R-space to G-space')
805 format(1x,'Performing serial FFT from R-space to G-space')
806 format(1x,'Normalizing wavefunctions in G-space',/)
807 format(1x,'Constructing real wavefunctions in G-space',/)
808 format(1x,'Writing wavefunctions to file',1x,a,/)
809 format(1x,'Reading SIESTA charge density from file',1x,a)
810 format(1x,'Writing charge density to file',1x,a,/)
811 format(1x,'Reading SIESTA total potential from file',1x,a)
812 format(1x,'Reading SIESTA electrostatic potential from file',1x,a)
813 format(1x,'Writing exchange-correlation potential to file',1x,a,/)
  
end program siesta2bgw
