!===============================================================================
!
! Program:
!
! (1) sapo             Originally By gsm        Last Modified 12/24/2012 (gsm)
!
! The SAPO (Simple Approximate Physical Orbitals) code reads DFT
! wavefunctions from WFN file, generates additional wavefunctions on top
! of them, and writes both to another WFN file. It can generate plane waves,
! decompose them into irreducible representations of the space group of the
! Bravais lattice or the crystal, and orthonormalize them with respect to
! DFT wavefunctions. It can insert SIESTA wavefunctions read from auxiliary
! WFN file in between plane waves and orthonormalize them altogether.
! It can apply a random variation to the plane waves and SIESTA wavefunctions
! or generate completely random wavefunctions and orthonormalize them.
! It can correct eigenvalues during orthonormalization or perform subspace
! diagonalization. It can drop the plane waves and SIESTA wavefunctions which
! have a large overlap with DFT wavefunctions. It can keep wavefunctions with
! the eigenvalues in a given energy range. It can extract eigenvalues and
! plane wave coefficients from WFN file for plotting purposes.
!
! The SAPO method is described in Phys. Rev. Lett. 107, 186404 (2011)
!
! 2012-12-24 (gsm): added iterative Davidson diagonalization, see
! examples/DFT/Si2_sapo for details
!
! Note: wfng_input_file should contain all occupied orbitals (and may or
! may not contain some unoccupied orbitals), otherwise occupations written
! to wfng_output_file will be wrong.
!
! Note: Symmetrization of planewaves (sapo_symmetry .gt. 0) has no effect
! on a SAPO calculation since it is just a linear combination in a degenerate
! subspace. It is currently disabled because it requires symmetry subroutines
! from Quantum ESPRESSO. To enable it open MeanField/SAPO/pw.f90, comment out
! die, and uncomment s_axis_to_cart, find_group, set_irr_rap, divide_class.
!
! Input is read from file sapo.inp.
!
!===============================================================================

#include "f_defs.h"

program sapo

  use global_m
  use check_inversion_m
  use fft_parallel_m
  use groupk_m
  use random_m
  use sort_m
  use symmetries_m
  use wfn_rho_vxc_io_m
  use pw_m
  use ortho_m
  use hdiag_m
  use norm_m
  use write_program_header_m
  implicit none

!------------------------
! Allocatable arrays

  integer, pointer :: atyp(:)
  real(DP), pointer :: apos(:,:)
  real(DP), allocatable :: kw(:)
  real(DP), allocatable :: kpt(:,:)
  integer, allocatable :: ifmin(:,:)
  integer, allocatable :: ifmax(:,:)
  real(DP), allocatable :: occ(:,:,:)
  real(DP), allocatable :: en(:,:,:)
  real(DP), allocatable :: ensrt(:)
  integer, allocatable :: gvec(:,:)
  integer, allocatable :: gvectmp(:,:)
  integer, allocatable :: ngk(:)
  SCALAR, allocatable :: wfn(:,:)
  SCALAR, allocatable :: wfnaux(:,:)
  SCALAR, allocatable :: wfn_d(:,:,:,:)
  SCALAR, allocatable :: wfnsrt_d(:,:)
  integer, allocatable :: index_vec(:)
  integer, allocatable :: isort(:)
  integer, allocatable :: srtmap(:)
  integer, allocatable :: idummy(:)
  integer, allocatable :: isort_d(:,:)
  integer, allocatable :: pw_num(:,:,:)
  integer, allocatable :: pw_ind(:,:,:,:,:)
  integer, allocatable :: isrti(:)
  integer, allocatable :: isrtort(:,:,:)
  SCALAR, allocatable :: vsc(:,:)
  SCALAR, allocatable :: vsc_d(:,:)
  integer, allocatable :: idrop(:)
  integer, allocatable :: ityp(:)
  integer, allocatable :: nh(:)
  real(DP), allocatable :: deeq(:,:,:,:)
  complex(DPC), allocatable :: vkb(:,:)
  complex(DPC), allocatable :: vkb_d(:,:,:,:)
  SCALAR, allocatable :: wfnsp(:)
  SCALAR, allocatable :: wfnspdummy(:)
  real(DP), allocatable :: ekin(:)
  integer, allocatable :: isrt(:)
  integer, allocatable :: rdeg(:,:)
  integer, allocatable :: ideg(:)
  real(DP), allocatable :: proj(:,:)
  real(DP), allocatable :: projdummy(:,:)
  real(DP), allocatable :: amplx(:)
  real(DP), allocatable :: amply(:)
  real(DP), allocatable :: ampldummy(:)
  real(DP), allocatable :: comp(:)
  real(DP), allocatable :: compdummy(:)
  real(DP), allocatable :: en_i(:,:,:)
  integer, pointer :: atyp_x(:)
  real(DP), pointer :: apos_x(:,:)
  integer, allocatable :: ngk_x(:)
  real(DP), allocatable :: kw_x(:)
  real(DP), allocatable :: kpt_x(:,:)
  integer, allocatable :: ifmin_x(:,:)
  integer, allocatable :: ifmax_x(:,:)
  real(DP), allocatable :: en_x(:,:,:)
  real(DP), allocatable :: occ_x(:,:,:)
  real(DP), pointer :: apos_lat(:,:)

!------------------------
! Local variables

  logical :: sapo_energy_match,sapo_print_ir,aux_flag, &
    sapo_random_norm,sapo_overlap_flag,sapo_orthonormal, &
    sapo_ortho_energy,sapo_energy_sort,sapo_hamiltonian, &
    sapo_ham_resinc,sapo_do_all_bands,sapo_energy_range, &
    sapo_check_norm,sapo_eigenvalue,sapo_amplitude, &
    pw_flag,vkbg_flag,out_flag,pwsym,ffastsp,ffastsp_ib, &
    forder,fecor,fsym,fout,fshift,fvkbg,fparafft,fparala, &
    fblas,fblas_ort,fblas_ham
  integer :: sapo_band_number,sapo_planewave_min, &
    sapo_planewave_max,sapo_symmetry,aux_band_min, &
    aux_band_max,sapo_random,sapo_ortho_block, &
    sapo_ortho_order,sapo_plot_kpoint,sapo_plot_spin, &
    sapo_plot_bandmin,sapo_plot_bandmax,sapo_plot_pwmin, &
    sapo_plot_pwmax,sapo_projection,sapo_ampl_num, &
    sapo_ham_nrestart,sapo_ham_ndiag,sapo_ham_ndim, &
    cell_symmetry,iflavor,i,j,k,iadd,iout,ierr,info, &
    nFFTgridpts,is,ik,ib,ib1,ib2,jb,ig,jg,nsym,nsyml,nsymk, &
    npwdeg,ns,nk,nb,nbinp,nbaux,nbauxtot,nbpw,nbtot,nbmax, &
    nbstart,nbend,ng,ngkmax,ng_l,ng_g,ngk_l,ngk_g,b2g, &
    ndeg,mdeg,ipw,ipw1,ipw2,iat,isp,ih,jh,nat,nsp,ikb,nkb, &
    nhm,seed,ncount,iblock,iorder,idir,FFTgrid(3),kgrid(3), &
    rot(3,3,48),rotl(3,3,48),rotk(3,3,48),values(8), &
    ns_x,ng_x,nsym_x,cell_symmetry_x,nat_x,nk_x,ngkmax_x, &
    FFTgrid_x(3),kgrid_x(3),rot_x(3,3,48),real_or_complex, &
    spacegroup
  real(DP) :: sapo_energy_shift,aux_energy_shift,sapo_random_ampl, &
    sapo_overlap_max,sapo_energy_min,sapo_energy_max,sapo_ampl_del, &
    sapo_ampl_brd,sapo_ham_tol,ecutrho,ecutwfn,arnd,brnd,x, &
    celvol,recvol,al,bl,eshift,dotmax,dotcur,tsec(2),a(3,3), &
    b(3,3),adot(3,3),bdot(3,3),kshift(3),kvec(3),tau(3,48), &
    ecutrho_x,ecutwfn_x,kshift_x(3),celvol_x,al_x,a_x(3,3), &
    adot_x(3,3),recvol_x,bl_x,b_x(3,3),bdot_x(3,3),tau_x(3,48), &
    taul(3,48)
  SCALAR :: z
  character :: stitle*32,sdate*32,stime*32
  character(len=256) :: wfng_input_file,wfng_aux_file, &
    vscg_input_file,vkbg_input_file,wfng_output_file,tmpstr
  character(len=16) :: s1,s2,s3,s4,s5,s6
  character(len=16) :: routnam(26)
  character :: sheader*3
  character(len=21) :: symbol

!------------------------
! Parameters

  integer, parameter :: drnd = 500000
  integer, parameter :: nrnd = 5000000
  character(len=256), parameter :: sapo_input_file = 'sapo.inp'

!------------------------
! Input namelist

  namelist /input_sapo/ wfng_input_file, wfng_aux_file, &
    vscg_input_file, vkbg_input_file, wfng_output_file, &
    sapo_band_number, sapo_planewave_min, sapo_planewave_max, &
    sapo_energy_shift, sapo_energy_match, sapo_symmetry, &
    sapo_print_ir, aux_flag, aux_band_min, aux_band_max, &
    aux_energy_shift, sapo_random, sapo_random_ampl, &
    sapo_random_norm, sapo_overlap_flag, sapo_overlap_max, &
    sapo_orthonormal, sapo_ortho_block, sapo_ortho_order, &
    sapo_ortho_energy, sapo_energy_sort, sapo_hamiltonian, &
    sapo_ham_nrestart, sapo_ham_ndiag, sapo_ham_ndim, &
    sapo_ham_tol, sapo_ham_resinc, sapo_do_all_bands, &
    sapo_energy_range, sapo_energy_min, sapo_energy_max, &
    sapo_check_norm, sapo_plot_kpoint, sapo_plot_spin, &
    sapo_plot_bandmin, sapo_plot_bandmax, sapo_plot_pwmin, &
    sapo_plot_pwmax, sapo_eigenvalue, sapo_projection, &
    sapo_amplitude, sapo_ampl_num, sapo_ampl_del, sapo_ampl_brd

!------------------------
! Inirialize MPI

  call peinfo_init()

#ifdef MPI
  fblas_ort=.true.
  fblas_ham=.true.
  fparafft=.true.
#else
  fblas_ort=.true.
  fblas_ham=.true.
  fparafft=.false.
#endif

#ifdef USESCALAPACK
  fparala=.true.
#else
  fparala=.false.
#endif

#ifdef CPLX
  real_or_complex = 2
#else
  real_or_complex = 1
#endif

!------------------------
! Check libraries

  ierr=0
  if (fparala) then
#ifndef USESCALAPACK
    call die('ScaLAPACK is not available.')
#endif
  endif ! fparala

!------------------------
! Initialize timer

  call timacc(0,0)
  
  call timacc(1,1)
  
  call timacc(2,1)

!------------------------
! Write header

  call write_program_header('SAPO', .false.)

!------------------------
! Read input file

  if (peinf%inode.eq.0) then
    write(6,801) trim(sapo_input_file)
    wfng_input_file = ''
    wfng_aux_file = ''
    vscg_input_file = ''
    vkbg_input_file = ''
    wfng_output_file = ''
    sapo_band_number = -1
    sapo_planewave_min = 0
    sapo_planewave_max = 0
    sapo_energy_shift = 0.0d0
    sapo_energy_match = .false.
    sapo_symmetry = 0
    sapo_print_ir = .false.
    aux_flag = .false.
    aux_band_min = 0
    aux_band_max = 0
    aux_energy_shift = 0.0d0
    sapo_random = 0
    sapo_random_ampl = 0.0d0
    sapo_random_norm = .false.
    sapo_overlap_flag = .false.
    sapo_overlap_max = 0.0d0
    sapo_orthonormal = .false.
    sapo_ortho_block = 0
    sapo_ortho_order = 0
    sapo_ortho_energy = .false.
    sapo_energy_sort = .false.
    sapo_hamiltonian = .false.
    sapo_ham_nrestart = 0
    sapo_ham_ndiag = 0
    sapo_ham_ndim = 0
    sapo_ham_tol = 0.0d0
    sapo_ham_resinc = .false.
    sapo_do_all_bands = .false.
    sapo_energy_range = .false.
    sapo_energy_min = 0.0d0
    sapo_energy_max = 0.0d0
    sapo_check_norm = .false.
    sapo_plot_kpoint = 0
    sapo_plot_spin = 0
    sapo_plot_bandmin = 0
    sapo_plot_bandmax = 0
    sapo_plot_pwmin = 0
    sapo_plot_pwmax = 0
    sapo_eigenvalue = .false.
    sapo_projection = 0
    sapo_amplitude = .false.
    sapo_ampl_num = 0
    sapo_ampl_del = 0.0d0
    sapo_ampl_brd = 0.0d0
    call open_file(55,file=sapo_input_file,status='old',form='formatted')
    read(55,input_sapo,iostat=ierr)
    if (ierr.eq.0) call close_file(55)
    if (sapo_band_number.lt.-1.or. &
      sapo_planewave_min.lt.0.or. &
      sapo_planewave_max.lt.sapo_planewave_min.or. &
      sapo_symmetry.lt.0.or.sapo_symmetry.gt.2.or. &
      aux_band_min.lt.0.or. &
      aux_band_max.lt.aux_band_min.or. &
      sapo_random.lt.0.or.sapo_random.gt.2.or. &
      sapo_random_ampl.lt.-TOL_Zero.or. &
      sapo_overlap_max.lt.-TOL_Zero.or. &
      sapo_overlap_max.gt.1.0d0+TOL_Zero.or. &
      sapo_ortho_block.lt.0.or.sapo_ortho_block.gt.2.or. &
      sapo_ortho_order.lt.0.or.sapo_ortho_order.gt.1.or. &
      sapo_energy_max.lt.(sapo_energy_min-TOL_Zero).or. &
      sapo_plot_kpoint.lt.0.or. &
      sapo_plot_spin.lt.0.or. &
      sapo_plot_bandmin.lt.0.or. &
      sapo_plot_bandmax.lt.sapo_plot_bandmin.or. &
      sapo_plot_pwmin.lt.0.or. &
      sapo_plot_pwmax.lt.sapo_plot_pwmin.or. &
      sapo_projection.lt.0.or.sapo_projection.gt.7.or. &
      (sapo_amplitude.and.(sapo_ampl_num.le.0.or. &
      sapo_ampl_del.lt.TOL_Small.or.sapo_ampl_brd.lt.TOL_Small))) &
      ierr=1
    if (ierr.eq.0) then
      pw_flag=sapo_planewave_max.ne.0
      if (pw_flag.and.sapo_energy_match.and. &
        sapo_planewave_min.lt.2) then
        write(0,921)
        sapo_energy_match=.false.
      endif
      if (pw_flag.and.sapo_energy_match.and. &
        abs(sapo_energy_shift).gt.TOL_Zero) then
        write(0,922)
        sapo_energy_shift=0.0d0
      endif
      if (sapo_symmetry.gt.0.and..not.pw_flag) then
        write(0,923)
        sapo_symmetry=0
      endif
      if (aux_flag.and.len(trim(wfng_aux_file)).eq.0) then
        write(0,924)
        aux_flag=.false.
      endif
      if (sapo_random.ne.0.and..not.sapo_do_all_bands.and. &
        .not.pw_flag.and..not.aux_flag) then
        write(0,937)
        sapo_random=0
      endif
      if (sapo_random.eq.2) then
        if (sapo_symmetry.gt.0) write(0,925)
        if (aux_flag) write(0,926)
      endif
      if (sapo_overlap_flag.and..not.pw_flag.and..not.aux_flag) then
        write(0,927)
        sapo_overlap_flag=.false.
      endif
      if (sapo_overlap_flag.and.(sapo_planewave_min.gt.1.or. &
        aux_band_min.gt.1)) then
        write(0,928)
      endif
      if (sapo_orthonormal.and..not.sapo_do_all_bands.and. &
        .not.pw_flag.and..not.aux_flag) then
        write(0,929)
        sapo_orthonormal=.false.
      endif
      if (sapo_orthonormal.and.sapo_ortho_block.eq.2.and. &
        (.not.pw_flag.or..not.aux_flag)) then
        write(0,930)
        if (.not.aux_flag) sapo_ortho_block=0
        if (.not.pw_flag) sapo_ortho_block=1
      endif
      if (sapo_overlap_flag.and.sapo_orthonormal.and. &
        pw_flag.and.aux_flag.and.sapo_ortho_block.ne.2) then
        write(0,931)
        sapo_overlap_flag=.false.
      endif
      if (sapo_orthonormal.and.sapo_ortho_energy.and. &
        .not.pw_flag.and..not.aux_flag) then
        write(0,932)
        sapo_ortho_energy=.false.
      endif
      if (sapo_energy_sort.and..not.pw_flag.and..not.aux_flag) then
        write(0,933)
        sapo_energy_sort=.false.
      endif
      if (sapo_hamiltonian.and.len(trim(vscg_input_file)).eq.0) then
        write(0,934)
        sapo_hamiltonian=.false.
      endif
      if (sapo_hamiltonian.and..not.sapo_do_all_bands.and. &
        .not.pw_flag.and..not.aux_flag) then
        write(0,935)
        sapo_hamiltonian=.false.
      endif
      if (sapo_energy_range.and..not.(sapo_energy_sort.or. &
        sapo_hamiltonian)) then
        write(0,936)
        sapo_energy_range=.false.
      endif
      vkbg_flag=sapo_hamiltonian.and.len(trim(vkbg_input_file)).ne.0
      ffastsp=pw_flag.and.sapo_random.eq.0.and.sapo_orthonormal
      out_flag=len(trim(wfng_output_file)).ne.0
    endif
  endif
#ifdef MPI
  call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  if (ierr.ne.0) then
    call die('failed to read input file ' // trim(sapo_input_file), only_root_writes = .true.)
  endif
#ifdef MPI
  call MPI_Bcast(wfng_input_file,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_aux_file,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(vscg_input_file,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(vkbg_input_file,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(wfng_output_file,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_band_number,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_planewave_min,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_planewave_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_energy_shift,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_energy_match,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_symmetry,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_print_ir,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(aux_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(aux_band_min,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(aux_band_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(aux_energy_shift,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_random,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_random_ampl,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_random_norm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_overlap_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_overlap_max,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_orthonormal,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ortho_block,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ortho_order,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ortho_energy,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_energy_sort,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_hamiltonian,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ham_nrestart,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ham_ndiag,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ham_ndim,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ham_tol,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ham_resinc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_do_all_bands,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_energy_range,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_energy_min,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_energy_max,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_check_norm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_plot_kpoint,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_plot_spin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_plot_bandmin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_plot_bandmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_plot_pwmin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_plot_pwmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_eigenvalue,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_projection,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_amplitude,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ampl_num,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ampl_del,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(sapo_ampl_brd,1,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(pw_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(vkbg_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(ffastsp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(out_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif

!------------------------
! Read header of input wavefunction file

  if (peinf%inode.eq.0) then
    call open_file(7,file=wfng_input_file,status='old',form='unformatted')
  endif ! peinf%inode.eq.0
  sheader = 'WFN'
  iflavor = 0
  call read_binary_header(7, sheader, iflavor, ns, ng, nsym, &
    cell_symmetry, nat, nk, nb, ngkmax, ecutrho, ecutwfn, FFTgrid, &
    kgrid, kshift, celvol, al, a, adot, recvol, bl, b, bdot, &
    rot, tau, atyp, apos, ngk, kw, kpt, ifmin, ifmax, en_i, occ, &
    dont_warn_kgrid = .true.)

!------------------------
! Check if real/complex version is appropriate

  call check_inversion(real_or_complex, nsym, rot, ns, .true., .true., tnp = tau)

!------------------------
! Read header of auxiliary wavefunction file

  nbauxtot = 0
  if (aux_flag) then
    if (peinf%inode.eq.0) then
      call open_file(8,file=wfng_aux_file,status='old',form='unformatted')
    endif ! peinf%inode.eq.0
    sheader = 'WFN'
    iflavor = 0
    call read_binary_header(8, sheader, iflavor, ns_x, ng_x, nsym_x, &
      cell_symmetry_x, nat_x, nk_x, nbauxtot, ngkmax_x, ecutrho_x, ecutwfn_x, &
      FFTgrid_x, kgrid_x, kshift_x, celvol_x, al_x, a_x, adot_x, recvol_x, &
      bl_x, b_x, bdot_x, rot_x, tau_x, atyp_x, apos_x, ngk_x, kw_x, kpt_x, &
      ifmin_x, ifmax_x, en_x, occ_x, dont_warn_kgrid = .true.)
    
    ierr=0
    if (ns_x .ne. ns) ierr = 1
    if (ng_x .ne. ng) ierr = 2
    if (nsym_x .ne. nsym) ierr = 3
    if (cell_symmetry_x .ne. cell_symmetry) ierr = 4
    if (nat_x .ne. nat) ierr = 5
    if (nk_x .ne. nk) ierr = 6
    if (ngkmax_x .ne. ngkmax) ierr = 7
    if (abs(ecutrho_x - ecutrho) .gt. TOL_Small) ierr = 8
    if (abs(ecutwfn_x - ecutwfn) .gt. TOL_Small) ierr = 9
    if (any(FFTgrid_x(1:3) .ne. FFTgrid(1:3))) ierr = 10
    if (any(kgrid_x(1:3) .ne. kgrid(1:3))) ierr = 11
    if (any(abs(kshift_x(1:3) - kshift(1:3)) .gt. TOL_Small)) ierr = 12
    if (abs(celvol_x / celvol - 1.0d0) .gt. TOL_Small) ierr = 13
    if (abs(recvol_x / recvol - 1.0d0) .gt. TOL_Small) ierr = 14
    if (any(abs(al_x * a_x(1:3, 1:3) / al - a(1:3, 1:3)) .gt. TOL_Small)) ierr = 15
    if (any(abs(bl_x * b_x(1:3, 1:3) / bl - b(1:3, 1:3)) .gt. TOL_Small)) ierr = 16
    !if (any(rot_x(1:3, 1:3, 1:nsym_x) .ne. rot(1:3, 1:3, 1:nsym))) ierr = 17
    !if (any(abs(tau_x(1:3, 1:nsym_x) - tau(1:3, 1:nsym)) .gt. TOL_Small)) ierr = 18
    if (any(atyp_x(1:nat_x) .ne. atyp(1:nat))) ierr = 19
    if (any(abs(al_x * apos_x(1:3, 1:nat_x) / al - apos(1:3, 1:nat)) .gt. TOL_Small)) ierr = 20
    if (any(ngk_x(1:nk_x) .ne. ngk(1:nk))) ierr = 21
    if (any(abs(kw_x(1:nk_x) - kw(1:nk)) .gt. TOL_Small)) ierr = 22
    if (any(abs(kpt_x(1:3, 1:nk_x) - kpt(1:3, 1:nk)) .gt. TOL_Small)) ierr = 23
    
    if (ierr.ne.0) then
      call die('different content in wavefunction files ' // trim(wfng_input_file) // ' and ' // &
        trim(wfng_aux_file), only_root_writes = .true.)
    endif ! ierr
    
    SAFE_DEALLOCATE_P(atyp_x)
    SAFE_DEALLOCATE_P(apos_x)
    SAFE_DEALLOCATE(ngk_x)
    SAFE_DEALLOCATE(kw_x)
    SAFE_DEALLOCATE(kpt_x)
    SAFE_DEALLOCATE(ifmin_x)
    SAFE_DEALLOCATE(ifmax_x)
    SAFE_DEALLOCATE(occ_x)
  endif ! aux_flag

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
  
  nbinp=sapo_band_number
  if (nbinp.eq.-1.or.nbinp.gt.nb) nbinp=nb
  
  if (sapo_planewave_min.eq.0) sapo_planewave_min=1
  if (sapo_planewave_max.eq.0) sapo_planewave_max=1
  if (sapo_planewave_min.gt.ng) sapo_planewave_min=ng
  if (sapo_planewave_max.gt.ng) sapo_planewave_max=ng
  if (pw_flag) then
    nbpw=sapo_planewave_max-sapo_planewave_min+1
    b2g=nbinp+1-sapo_planewave_min
  else
    nbpw=0
    b2g=0
  endif
  
  if (aux_band_min.eq.0) aux_band_min=1
  if (aux_band_max.eq.0) aux_band_max=nbauxtot
  if (aux_band_min.gt.nbauxtot) aux_band_min=nbauxtot
  if (aux_band_max.gt.nbauxtot) aux_band_max=nbauxtot
  if (aux_flag) then
    nbaux=aux_band_max-aux_band_min+1
  else
    nbaux=0
  endif
  
  nbtot=nbinp+nbpw+nbaux
  nbmax=nbtot
  
  if (peinf%inode.eq.0) then
    
    write(s1,101)nsym
    write(s2,101)nk
    write(s3,101)ns
    write(s4,101)nb
    write(s5,101)ng
    write(s6,101)ngkmax
    write(6,110)TRUNC(s1),TRUNC(s2),TRUNC(s3),TRUNC(s4), &
      TRUNC(s5),TRUNC(s6)
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
    write(s1,101)nbinp
    write(s2,101)nbpw
    write(s3,101)b2g
    write(s4,101)nbaux
    write(s5,101)nbtot
    write(6,140)TRUNC(s1),TRUNC(s2),TRUNC(s3),TRUNC(s4), &
      TRUNC(s5)
    
  endif

!------------------------
! Distribute input and auxiliary eigenvalues

  SAFE_ALLOCATE(en, (nbmax,nk,ns))
  
  do is=1,ns
    do ik=1,nk
      do ib=1,nbinp
        en(ib,ik,is)=en_i(ib,ik,is)
      enddo
    enddo
  enddo
  do is=1,ns
    do ik=1,nk
      do ib=nbinp+1,nbinp+nbpw
        en(ib,ik,is)=0.0d0
      enddo
    enddo
  enddo
  if (aux_flag) then
    do is=1,ns
      do ik=1,nk
        do ib=aux_band_min,aux_band_max
          jb=ib-aux_band_max+nbtot
          en(jb,ik,is)=en_x(ib,ik,is)+aux_energy_shift/RYD
        enddo
      enddo
    enddo
  endif
  
  SAFE_DEALLOCATE(en_i)
  if (aux_flag) then
    SAFE_DEALLOCATE(en_x)
  endif
  
#ifdef MPI
  call MPI_Bcast(en(1,1,1),nbmax*nk*ns,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif

!------------------------
! Allocate distributed arrays for wavefunctions

  SAFE_ALLOCATE(wfn_d, (ngk_l,nbmax,ns,nk))
  SAFE_ALLOCATE(isort_d, (ngk_l,nk))

!------------------------
! Reset output wavefunction file

  if (out_flag) then
    if (peinf%inode.eq.0) then
      call open_file(9,file=wfng_output_file,status='replace',form='unformatted')
      if (ierr.eq.0) call close_file(9)
    endif
#ifdef MPI
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
    if (ierr.ne.0) then
      call die('failed to write wavefunction file ' // trim(wfng_output_file), only_root_writes = .true.)
    endif
    if (peinf%inode.eq.0) then
      call open_file(9,file=wfng_output_file,status='old',form='unformatted')
    endif
  endif ! out_flag
  
  call timacc(2,2)
  
  call timacc(3,1)

!------------------------
! Read master list of G-vectors

  if (peinf%inode.eq.0) then
    write(6,802) trim(wfng_input_file)
  endif
  
  SAFE_ALLOCATE(gvec, (3, ng))
  call read_binary_gvectors(7, ng, ng, gvec)

!------------------------
! Allocate and compute index_vec indices

  SAFE_ALLOCATE(index_vec, (nFFTgridpts))
  index_vec=0
  do ig=1,ng
    iadd=((gvec(1,ig)+FFTgrid(1)/2)*FFTgrid(2)+gvec(2,ig)+FFTgrid(2)/2)* &
      FFTgrid(3)+gvec(3,ig)+FFTgrid(3)/2+1
    index_vec(iadd)=ig
  enddo

!------------------------
! Loop over k-points

  SAFE_ALLOCATE(gvectmp, (3, ngkmax))
  SAFE_ALLOCATE(wfn, (ngk_g,ns))
  SAFE_ALLOCATE(isort, (ngk_g))
  
  do ik=1,nk

!------------------------
! Read g-vectors for current k-point
! Determine their indices in list of g-vectors
! Distribute index arrays over processors

    call read_binary_gvectors(7, ngk(ik), ngkmax, gvectmp, bcast = .false.)
    if (peinf%inode.eq.0) then
      ierr=0
      isort(:)=0
      do ig=1,ngk(ik)
        info=0
        iout=((gvectmp(1,ig)+FFTgrid(1)/2)*FFTgrid(2)+gvectmp(2,ig)+FFTgrid(2)/2)* &
          FFTgrid(3)+gvectmp(3,ig)+FFTgrid(3)/2+1
        if (iout.ge.1.and.iout.le.nFFTgridpts) then
          iout=index_vec(iout)
          if (iout.ge.1.and.iout.le.ng) then
            if (gvectmp(1,ig).ne.gvec(1,iout).or.gvectmp(2,ig).ne. &
              gvec(2,iout).or.gvectmp(3,ig).ne.gvec(3,iout)) info=-1
          else
            info=-1
          endif
        else
          info=-1
        endif
        isort(ig)=iout
        ierr=ierr+info
      enddo
    endif
#ifdef MPI
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
    if (ierr.ne.0) then
      write(s1,101)ik
      call die('failed to find G-vector for k-point ' // TRUNC(s1), only_root_writes = .true.)
    endif
#ifdef MPI
    call MPI_Barrier(MPI_COMM_WORLD,mpierr)
    call MPI_Scatter(isort,ngk_l,MPI_INTEGER,isort_d(:,ik),ngk_l, &
      MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#else
    isort_d(:,ik)=isort(:)
#endif

!------------------------
! Read band and spin wavefunctions for current k-point
! and distribute over processors

    do ib=1,nbinp
      wfn(:,:)=ZERO
      call read_binary_data(7, ngk(ik), ngk_g, ns, wfn, bcast = .false.)
      do is=1,ns
#ifdef MPI
        ! FHJ: The following barrier is likely superfluous, but can help
        ! mitigate buffer overflows. Uncomment it if you experience MPI errors.
        !call MPI_Barrier(MPI_COMM_WORLD,mpierr) ! Uncomment if you have MPI errors
        call MPI_Scatter(wfn(:,is),ngk_l,MPI_SCALAR, &
          wfn_d(:,ib,is,ik),ngk_l,MPI_SCALAR, &
          0,MPI_COMM_WORLD,mpierr)
#else
        wfn_d(:,ib,is,ik)=wfn(:,is)
#endif
      enddo
    enddo
    
    if (ik .lt. nk) then
      do ib=nbinp+1,nb
        call read_binary_data(7, ngk(ik), ngk_g, ns, wfn, dont_read = .true.)
      enddo
    endif
    
  enddo ! ik
  
  SAFE_DEALLOCATE(gvectmp)
  SAFE_DEALLOCATE(wfn)
  SAFE_DEALLOCATE(isort)

  call timacc(3,2)
  
  call timacc(4,1)

!------------------------
! Read auxiliary wavefunction file

  if (aux_flag) then
    
    if (peinf%inode.eq.0) then
      write(6,802) trim(wfng_aux_file)
    endif
    
    SAFE_ALLOCATE(gvectmp, (3, ng))
    call read_binary_gvectors(8, ng, ng, gvectmp, bcast = .false.)
    SAFE_DEALLOCATE(gvectmp)

!------------------------
! Loop over k-points

    SAFE_ALLOCATE(gvectmp, (3, ngkmax))
    SAFE_ALLOCATE(wfn, (ngk_g,ns))
    SAFE_ALLOCATE(isort, (ngk_g))
    SAFE_ALLOCATE(srtmap, (ngk_g))
    SAFE_ALLOCATE(idummy, (ngk_g))
    SAFE_ALLOCATE(wfnaux, (ngkmax,ns))
    
    do ik=1,nk

!------------------------
! Read g-vectors for current k-point
! Determine their indices in list of g-vectors

      call read_binary_gvectors(8, ngk(ik), ngkmax, gvectmp, bcast = .false.)
      if (peinf%inode.eq.0) then
        ierr=0
        isort(:)=0
        do ig=1,ngk(ik)
          info=0
          iout=((gvectmp(1,ig)+FFTgrid(1)/2)*FFTgrid(2)+gvectmp(2,ig)+FFTgrid(2)/2)* &
            FFTgrid(3)+gvectmp(3,ig)+FFTgrid(3)/2+1
          if (iout.ge.1.and.iout.le.nFFTgridpts) then
            iout=index_vec(iout)
            if (iout.ge.1.and.iout.le.ng) then
              if (gvectmp(1,ig).ne.gvec(1,iout).or.gvectmp(2,ig).ne. &
                gvec(2,iout).or.gvectmp(3,ig).ne.gvec(3,iout)) info=-1
            else
              info=-1
            endif
          else
            info=-1
          endif
          isort(ig)=iout
          ierr=ierr+info
        enddo
      endif
#ifdef MPI
      call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
      if (ierr.ne.0) then
        write(s1,101)ik
        call die('failed to find G-vector for k-point ' // TRUNC(s1), only_root_writes = .true.)
      endif
#ifdef MPI
      call MPI_Bcast(isort,ngk_g,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif

!------------------------
! Build a sorting map between aux isort indices
! and inp isort_d distributed indices

      srtmap(:)=0
      do ig=1,ngk(ik)
        do jg=1,ngk_l
          if (isort(ig).eq.isort_d(jg,ik)) &
            srtmap(ig)=peinf%inode*ngk_l+jg
        enddo
      enddo
#ifdef MPI
      idummy(:)=srtmap(:)
      call MPI_Reduce(idummy,srtmap,ngk_g, &
        MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
#endif

!------------------------
! Read band and spin wavefunctions for current k-point
! reorder them according to the sorting map
! and distribute over processors

      do ib=1,aux_band_min-1
        call read_binary_data(8, ngk(ik), ngkmax, ns, wfnaux, dont_read = .true.)
      enddo
      do ib=aux_band_min,aux_band_max
        jb=ib-aux_band_max+nbtot
        call read_binary_data(8, ngk(ik), ngkmax, ns, wfnaux, bcast = .false.)
        if (peinf%inode.eq.0) then
          wfn(:,:)=ZERO
          do is=1,ns
            do ig=1,ngk(ik)
              wfn(srtmap(ig),is)=wfnaux(ig,is)
            enddo
          enddo
        endif
        do is=1,ns
#ifdef MPI
          !call MPI_Barrier(MPI_COMM_WORLD,mpierr) ! Uncomment if you have MPI errors
          call MPI_Scatter(wfn(:,is),ngk_l,MPI_SCALAR, &
            wfn_d(:,jb,is,ik),ngk_l,MPI_SCALAR, &
            0,MPI_COMM_WORLD,mpierr)
#else
          wfn_d(:,jb,is,ik)=wfn(:,is)
#endif
        enddo
      enddo
      do ib=aux_band_max+1,nbauxtot
        call read_binary_data(8, ngk(ik), ngkmax, ns, wfnaux, dont_read = .true.)
      enddo
      
    enddo ! ik
    
    SAFE_DEALLOCATE(gvectmp)
    SAFE_DEALLOCATE(wfn)
    SAFE_DEALLOCATE(isort)
    SAFE_DEALLOCATE(srtmap)
    SAFE_DEALLOCATE(idummy)
    SAFE_DEALLOCATE(wfnaux)
    
  endif ! aux_flag
  
  call timacc(4,2)
  
  call timacc(11,1)

!------------------------
! Check self-consistent potential file

  if (sapo_hamiltonian) then
    if (peinf%inode.eq.0) then
      call open_file(10,file=vscg_input_file,status='old',form='unformatted')
    endif ! peinf%inode.eq.0
    sheader = 'VXC'
    iflavor = 0
    call read_binary_header(10, sheader, iflavor, ns_x, ng_x, nsym_x, &
      cell_symmetry_x, nat_x, nk_x, nbauxtot, ngkmax_x, ecutrho_x, ecutwfn_x, &
      FFTgrid_x, kgrid_x, kshift_x, celvol_x, al_x, a_x, adot_x, recvol_x, &
      bl_x, b_x, bdot_x, rot_x, tau_x, atyp_x, apos_x, ngk_x, kw_x, kpt_x, &
      ifmin_x, ifmax_x, en_x, occ_x, dont_warn_kgrid = .true.)

    ierr=0
    if (ns_x .ne. ns) ierr = 1
    if (ng_x .ne. ng) ierr = 2
    if (nsym_x .ne. nsym) ierr = 3
    if (cell_symmetry_x .ne. cell_symmetry) ierr = 4
    if (nat_x .ne. nat) ierr = 5
    if (abs(ecutrho_x - ecutrho) .gt. TOL_Small) ierr = 8
    if (any(FFTgrid_x(1:3) .ne. FFTgrid(1:3))) ierr = 10
    if (abs(celvol_x / celvol - 1.0d0) .gt. TOL_Small) ierr = 13
    if (abs(recvol_x / recvol - 1.0d0) .gt. TOL_Small) ierr = 14
    if (any(abs(al_x * a_x(1:3, 1:3) / al - a(1:3, 1:3)) .gt. TOL_Small)) ierr = 15
    if (any(abs(bl_x * b_x(1:3, 1:3) / bl - b(1:3, 1:3)) .gt. TOL_Small)) ierr = 16
    !if (any(rot_x(1:3, 1:3, 1:nsym_x) .ne. rot(1:3, 1:3, 1:nsym))) ierr = 17
    !if (any(abs(tau_x(1:3, 1:nsym_x) - tau(1:3, 1:nsym)) .gt. TOL_Small)) ierr = 18
    if (any(atyp_x(1:nat_x) .ne. atyp(1:nat))) ierr = 19
    if (any(abs(al_x * apos_x(1:3, 1:nat_x) / al - apos(1:3, 1:nat)) .gt. TOL_Small)) ierr = 20

    if (ierr.ne.0) then
      call die('failed to read potential file ' // trim(vscg_input_file), only_root_writes = .true.)
    endif ! ierr.ne.0

    SAFE_DEALLOCATE_P(atyp_x)
    SAFE_DEALLOCATE_P(apos_x)
  endif ! sapo_hamiltonian

!------------------------
! Check Kleinman-Bylander projectors file

  if (vkbg_flag) then
    if (peinf%inode.eq.0) then
      call open_file(11,file=vkbg_input_file,status='old',form='unformatted')
      
      read(11) stitle, sdate, stime
      read(11) ns_x, ng_x, nsym_x, cell_symmetry_x, nat_x, ecutrho_x, &
        nk_x, nsp, nkb, nhm, ngkmax_x, ecutwfn_x
      read(11) (FFTgrid_x(i), i = 1, 3), (kgrid_x(i), i = 1, 3), (kshift_x(i), i = 1, 3)
      read(11) celvol_x, al_x, ((a_x(j, i), j = 1, 3), i = 1, 3), &
        ((adot_x(j, i), j = 1, 3), i = 1, 3)
      read(11) recvol_x, bl_x, ((b_x(j, i), j = 1, 3), i = 1, 3), &
        ((bdot_x(j, i), j = 1, 3), i = 1, 3)
      read(11) (((rot_x(k, j, i), k = 1, 3), j = 1, 3), i = 1, nsym_x)
      read(11) ((tau_x(j, i), j = 1, 3), i = 1, nsym_x)
      
      SAFE_ALLOCATE(atyp_x, (nat_x))
      SAFE_ALLOCATE(apos_x, (3, nat_x))
      SAFE_ALLOCATE(ngk_x, (nk_x))
      SAFE_ALLOCATE(kw_x, (nk_x))
      SAFE_ALLOCATE(kpt_x, (3, nk_x))
      
      read(11) ((apos_x(j, i), j = 1, 3), atyp_x(i), i = 1, nat_x)
      read(11) (ngk_x(i), i = 1, nk_x)
      read(11) (kw_x(i), i = 1, nk_x)
      read(11) ((kpt_x(j, i), j = 1, 3), i = 1, nk_x)
      
      ierr=0
      if (stitle(1:11) /= 'VKB-Complex') ierr = -1
      if (ns_x .ne. ns) ierr = 1
      if (ng_x .ne. ng) ierr = 2
      if (nsym_x .ne. nsym) ierr = 3
      if (cell_symmetry_x .ne. cell_symmetry) ierr = 4
      if (nat_x .ne. nat) ierr = 5
      if (nk_x .ne. nk) ierr = 6
      if (ngkmax_x .ne. ngkmax) ierr = 7
      if (abs(ecutrho_x - ecutrho) .gt. TOL_Small) ierr = 8
      if (abs(ecutwfn_x - ecutwfn) .gt. TOL_Small) ierr = 9
      if (any(FFTgrid_x(1:3) .ne. FFTgrid(1:3))) ierr = 10
      if (any(kgrid_x(1:3) .ne. kgrid(1:3))) ierr = 11
      if (any(abs(kshift_x(1:3) - kshift(1:3)) .gt. TOL_Small)) ierr = 12
      if (abs(celvol_x / celvol - 1.0d0) .gt. TOL_Small) ierr = 13
      if (abs(recvol_x / recvol - 1.0d0) .gt. TOL_Small) ierr = 14
      if (any(abs(al_x * a_x(1:3, 1:3) / al - a(1:3, 1:3)) .gt. TOL_Small)) ierr = 15
      if (any(abs(bl_x * b_x(1:3, 1:3) / bl - b(1:3, 1:3)) .gt. TOL_Small)) ierr = 16
      !if (any(rot_x(1:3, 1:3, 1:nsym_x) .ne. rot(1:3, 1:3, 1:nsym))) ierr = 17
      !if (any(abs(tau_x(1:3, 1:nsym_x) - tau(1:3, 1:nsym)) .gt. TOL_Small)) ierr = 18
      if (any(atyp_x(1:nat_x) .ne. atyp(1:nat))) ierr = 19
      if (any(abs(al_x * apos_x(1:3, 1:nat_x) / al - apos(1:3, 1:nat)) .gt. TOL_Small)) ierr = 20
      if (any(ngk_x(1:nk_x) .ne. ngk(1:nk))) ierr = 21
      if (any(abs(kw_x(1:nk_x) - kw(1:nk)) .gt. TOL_Small)) ierr = 22
      if (any(abs(kpt_x(1:3, 1:nk_x) - kpt(1:3, 1:nk)) .gt. TOL_Small)) ierr = 23
      
      SAFE_DEALLOCATE_P(atyp_x)
      SAFE_DEALLOCATE_P(apos_x)
      SAFE_DEALLOCATE(ngk_x)
      SAFE_DEALLOCATE(kw_x)
      SAFE_DEALLOCATE(kpt_x)
    endif ! peinf%inode.eq.0
#ifdef MPI
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
    if (ierr.ne.0) then
      call die('failed to read potential file ' // trim(vkbg_input_file), only_root_writes = .true.)
    endif ! ierr.ne.0
    
#ifdef MPI
    call MPI_Bcast(nsp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nkb,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nhm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
    SAFE_ALLOCATE(ityp, (nat))
    SAFE_ALLOCATE(nh, (nsp))
    SAFE_ALLOCATE(deeq, (nhm,nhm,nat,ns))
    if (peinf%inode.eq.0) then
      read(11) (ityp(iat), iat = 1, nat)
      read(11) (nh(isp), isp = 1, nsp)
      read(11) ((((deeq(jh, ih, iat, is), jh = 1, nhm), ih = 1, nhm), iat = 1, nat), is = 1, ns)
      SAFE_ALLOCATE(gvectmp, (3, ng))
      call read_binary_gvectors(11, ng, ng, gvectmp, bcast = .false.)
      SAFE_DEALLOCATE(gvectmp)
    endif
#ifdef MPI
    call MPI_Bcast(ityp,nat,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nh,nsp,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(deeq(1,1,1,1),nhm*nhm*nat*ns,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif
  endif ! vkbg_flag
  
  call timacc(11,2)
  
  call timacc(2,1)

!------------------------
! Generate symmetries of the bravais lattice or the crystal

  if (sapo_symmetry.eq.2) then
    SAFE_ALLOCATE(apos_lat, (3, nat))
    apos_lat = matmul(transpose(b),apos) / al ! convert to crystal coords
    call get_symmetries(nat, atyp, apos_lat, a, FFTgrid, cell_symmetry, nsyml, rotl, taul, spacegroup, symbol)
    SAFE_DEALLOCATE_P(apos_lat)
  elseif (sapo_symmetry.eq.1) then
    nsyml=nsym
    rotl(:,:,:)=rot(:,:,:)
  else
    nsyml=1
    rotl(:,:,1)=0
    do i=1,3
      rotl(i,i,1)=1
    enddo
  endif

!------------------------
! Initialize random numbers

  if (sapo_random.gt.0) then
    call date_and_time(VALUES=values)
    seed=((values(5)*24+values(6))*60+values(7))*1000+values(8)
    seed=seed+peinf%inode*drnd
    call genrand_init(put=seed)
    do i=1,nrnd
      call genrand_real4(x)
    enddo
    if (sapo_random.eq.1) then
      ! add a small random variation to the wavefunctions
      arnd=1.0d0
      brnd=sapo_random_ampl
    elseif (sapo_random.eq.2) then
      ! generate completely random wavefunctions
      arnd=0.0d0
      brnd=1.0d0
    endif
  endif

!------------------------
! Allocate arrays for fast scalar products

  if (ffastsp) then
    if (sapo_symmetry.eq.0) then
      npwdeg=1
    else
      npwdeg=48
    endif
    SAFE_ALLOCATE(pw_num, (nbmax,ns,nk))
    SAFE_ALLOCATE(pw_ind, (2,npwdeg,nbmax,ns,nk))
    pw_num(:,:,:)=0
    pw_ind(:,:,:,:,:)=0
  else ! ffastsp
    SAFE_ALLOCATE(pw_num, (1,1,1))
    SAFE_ALLOCATE(pw_ind, (1,1,1,1,1))
  endif ! ffastsp
  
  call timacc(2,2)
  
  call timacc(5,1)

!------------------------
! Generate plane waves and optionally decompose them
! into irreducible representations of the space group
! of the Bravais lattice or the crystal

  if (pw_flag) then
    
    if (peinf%inode.eq.0) then
      if (sapo_symmetry.eq.0) write(6,804)
      if (sapo_symmetry.eq.1) write(6,805)
      if (sapo_symmetry.eq.2) write(6,806)
    endif
    
    fsym=sapo_symmetry.gt.0
    fout=sapo_print_ir
    fshift=sapo_energy_match
    eshift=sapo_energy_shift
    
    call pwgen(fsym,fout,ffastsp,fshift,nk,ns,nbinp, &
      nbpw,nbmax,nsyml,b2g,npwdeg,ngk_l,ng,rotl,gvec, &
      pw_num,pw_ind,isort_d,eshift,a,b,bdot,kpt,en,wfn_d)
    
  endif ! pw_flag
  
  call timacc(5,2)
  
  call timacc(6,1)

!------------------------
! Add random variation to wavefunctions
! Normalize randomized wavefunctions

  if (sapo_random.ne.0) then

    if (peinf%inode.eq.0) then
      if (sapo_random.eq.1) then
        write(s1,102)brnd
        write(6,807)TRUNC(s1)
      elseif (sapo_random.eq.2) then
        write(6,808)
      endif ! sapo_random
    endif

    if (sapo_do_all_bands) then
      nbstart=1
    else
      nbstart=nbinp+1
    endif
    nbend=nbtot

    do ik=1,nk
      do is=1,ns

        call random_wfng(nk,ns,nbstart,nbend,nbmax,ngk_l,ik,is,ngk, &
          arnd,brnd,wfn_d)

      enddo ! is
    enddo ! ik

    if (sapo_random_norm) then

      call norm_wfng(ngk_l,nbstart,nbend,nbmax,ns,nk,wfn_d)

    endif ! sapo_random_norm

  endif ! sapo_random.ne.0
  
  call timacc(6,2)
  
  call timacc(7,1)
  
!------------------------
! Drop PW & AUX states that have a large overlap with DFT states,
! drop PW states that have a large overlap with AUX states

  if (sapo_overlap_flag) then
    
    if (peinf%inode.eq.0) write(6,809)
    
    SAFE_ALLOCATE(idrop, (nbinp+1:nbtot))
    SAFE_ALLOCATE(wfnsp, (max(nbinp,nbaux)))
#ifdef MPI
    SAFE_ALLOCATE(wfnspdummy, (max(nbinp,nbaux)))
#endif
    
    ib1=0
    do ik=1,nk
      do is=1,ns
        ! check overlap of PW & AUX states with DFT states
        if (peinf%inode.eq.0) then
          write(s1,101)ik
          write(s2,101)is
          write(tmpstr,601)TRUNC(s1),TRUNC(s2)
          call open_file(12,file=tmpstr,status='replace',form='formatted')
        endif ! peinf%inode.eq.0
        idrop(:)=0
        do ib=nbinp+1,nbtot
          if (ffastsp) then
            ffastsp_ib=pw_num(ib,is,ik).ne.0
          else ! ffastsp
            ffastsp_ib=.false.
          endif ! ffastsp
          wfnsp(:)=ZERO
          if (ffastsp_ib) then
            do jb=1,nbinp
              z=ZERO
              do i=1,pw_num(ib,is,ik)
                if (peinf%inode.eq.pw_ind(1,i,ib,is,ik)) &
                  z=z+MYCONJG(wfn_d(pw_ind(2,i,ib,is,ik),jb,is,ik))* &
                  wfn_d(pw_ind(2,i,ib,is,ik),ib,is,ik)
              enddo ! i
              wfnsp(jb)=z
            enddo ! jb
          else ! ffastsp_ib
            do jb=1,nbinp
              z=ZERO
              do ig=1,ngk_l
                z=z+MYCONJG(wfn_d(ig,jb,is,ik))*wfn_d(ig,ib,is,ik)
              enddo ! ig
              wfnsp(jb)=z
            enddo ! jb
          endif ! ffastsp_ib
#ifdef MPI
          wfnspdummy(:)=wfnsp(:)
          call MPI_Allreduce(wfnspdummy,wfnsp,nbinp, &
            MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
          if (maxval(abs(wfnsp(1:nbinp))).gt.sapo_overlap_max) idrop(ib)=1
          if (peinf%inode.eq.0) then
            if (ierr.eq.0) write(12,603)ib,maxval(abs(wfnsp(1:nbinp)))
          endif ! peinf%inode.eq.0
        enddo ! ib
        if (peinf%inode.eq.0) then
          if (ierr.eq.0) call close_file(12)
        endif ! peinf%inode.eq.0
        ! check overlap of PW states with AUX states
        if (pw_flag.and.aux_flag) then
          if (peinf%inode.eq.0) then
            write(s1,101)ik
            write(s2,101)is
            write(tmpstr,602)TRUNC(s1),TRUNC(s2)
            call open_file(12,file=tmpstr,status='replace',form='formatted')
          endif ! peinf%inode.eq.0
          do ib=nbinp+1,nbinp+nbpw
            if (idrop(ib).eq.1) cycle
            if (ffastsp) then
              ffastsp_ib=pw_num(ib,is,ik).ne.0
            else ! ffastsp
              ffastsp_ib=.false.
            endif ! ffastsp
            wfnsp(:)=ZERO
            if (ffastsp_ib) then
              do jb=nbinp+nbpw+1,nbtot
                z=ZERO
                do i=1,pw_num(ib,is,ik)
                  if (peinf%inode.eq.pw_ind(1,i,ib,is,ik)) &
                    z=z+MYCONJG(wfn_d(pw_ind(2,i,ib,is,ik),jb,is,ik))* &
                    wfn_d(pw_ind(2,i,ib,is,ik),ib,is,ik)
                enddo ! i
                wfnsp(jb-nbinp-nbpw)=z
              enddo ! jb
            else ! ffastsp_ib
              do jb=nbinp+nbpw+1,nbtot
                z=ZERO
                do ig=1,ngk_l
                  z=z+MYCONJG(wfn_d(ig,jb,is,ik))*wfn_d(ig,ib,is,ik)
                enddo ! ig
                wfnsp(jb-nbinp-nbpw)=z
              enddo ! jb
            endif ! ffastsp_ib
#ifdef MPI
            wfnspdummy(:)=wfnsp(:)
            call MPI_Allreduce(wfnspdummy,wfnsp,nbaux, &
              MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
            if (maxval(abs(wfnsp(1:nbaux))).gt.sapo_overlap_max) idrop(ib)=1
            if (peinf%inode.eq.0) then
              if (ierr.eq.0) write(12,603)ib,maxval(abs(wfnsp(1:nbaux)))
            endif ! peinf%inode.eq.0
          enddo ! ib
          if (peinf%inode.eq.0) then
            if (ierr.eq.0) call close_file(12)
          endif ! peinf%inode.eq.0
        endif ! pw_flag.and.aux_flag
        ! keep record of maximum number of PW & AUX states we drop 
        if (ib1.lt.sum(idrop(nbinp+1:nbtot))) ib1=sum(idrop(nbinp+1:nbtot))
        ! drop PW & AUX states that have large overlap
        do ib=nbinp+1,nbtot
          ib2=sum(idrop(nbinp+1:ib))
          if (idrop(ib).eq.0.and.ib2.gt.0) then
            if (ffastsp) then
              pw_num(ib-ib2,is,ik)=pw_num(ib,is,ik)
              pw_ind(:,:,ib-ib2,is,ik)=pw_ind(:,:,ib,is,ik)
            endif ! ffastsp
            en(ib-ib2,ik,is)=en(ib,ik,is)
            wfn_d(:,ib-ib2,is,ik)=wfn_d(:,ib,is,ik)
          endif ! idrop
        enddo ! ib
      enddo ! is
    enddo ! ik
    
    if (ib1.gt.0) then
      nbtot=nbtot-ib1
      if (.not.aux_flag) then
        nbpw=nbpw-ib1
      elseif (.not.pw_flag) then
        nbaux=nbaux-ib1
      else
        nbpw=0
        nbaux=0
      endif
      if (peinf%inode.eq.0) then
        write(s1,101)ib1
        write(s2,101)nbtot
        write(6,810)TRUNC(s1),TRUNC(s2)
      endif ! peinf%inode.eq.0
    else ! ib1.gt.0
      if (peinf%inode.eq.0) write(6,*)
    endif ! ib1.gt.0
    
    SAFE_DEALLOCATE(idrop)
    SAFE_DEALLOCATE(wfnsp)
#ifdef MPI
    SAFE_DEALLOCATE(wfnspdummy)
#endif
    
  endif ! sapo_overlap_flag
  
  call timacc(7,2)
  
  call timacc(8,1)

!------------------------
! Orthonormalize plane wave and auxiliary wavefunctions
! with respect to input wavefunctions

  if (sapo_orthonormal) then
    
    if (peinf%inode.eq.0) write(6,811)
    
    forder=sapo_ortho_block.ne.0.or.sapo_ortho_order.ne.0
    
    if (forder) then
      ! build the index array that specifies in which
      ! order to orthonormalize the wavefunctions
      SAFE_ALLOCATE(isrtort, (nbtot,ns,nk))
      iblock=sapo_ortho_block
      iorder=sapo_ortho_order
      call orderindex(iblock,iorder,nk,ns,nbtot,nbinp,nbpw,nbaux, &
        nbmax,isrtort,en)
      idir=1
      call orderarray(ffastsp,idir,nk,ns,nbtot,nbmax,ngk_l,npwdeg, &
        isrtort,pw_num,pw_ind,en,wfn_d)
    endif
    
    fblas=fblas_ort
    fecor=sapo_ortho_energy
    if (sapo_do_all_bands) then
      nbstart=1
    else
      nbstart=nbinp+1
    endif
    nbend=nbtot

    call orthonormalize(ffastsp,fblas,fecor,nk,ns,nbstart, &
      nbend,nbmax,ngk_l,npwdeg,pw_num,pw_ind,en,wfn_d)
    
    if (forder) then
      idir=2
      call orderarray(ffastsp,idir,nk,ns,nbtot,nbmax,ngk_l,npwdeg, &
        isrtort,pw_num,pw_ind,en,wfn_d)
      SAFE_DEALLOCATE(isrtort)
    endif ! forder
    
  endif ! sapo_orthonormal
  
  call timacc(8,2)
  
  call timacc(2,1)

!------------------------
! Deallocate arrays for fast scalar products

  SAFE_DEALLOCATE(pw_num)
  SAFE_DEALLOCATE(pw_ind)
  
  call timacc(2,2)
  
  call timacc(9,1)

!------------------------
! Sort plane wave and auxiliary wavefunctions by eigenvalues

  if (sapo_energy_sort) then
    
    if (peinf%inode.eq.0) write(6,812)
    
    SAFE_ALLOCATE(ensrt, (nbtot-nbinp))
    SAFE_ALLOCATE(isrt, (nbtot-nbinp))
    SAFE_ALLOCATE(wfnsrt_d, (ngk_l,nbtot-nbinp))
    
    do ik=1,nk
      do is=1,ns
        do ib=nbinp+1,nbtot
          ensrt(ib-nbinp)=en(ib,ik,is)
        enddo ! ib
        call sortrx(nbtot-nbinp,ensrt,isrt)
        do ib=nbinp+1,nbtot
          en(ib,ik,is)=ensrt(isrt(ib-nbinp))
        enddo ! ib
        do ib=nbinp+1,nbtot
          wfnsrt_d(:,ib-nbinp)=wfn_d(:,ib,is,ik)
        enddo ! ib
        do ib=nbinp+1,nbtot
          wfn_d(:,ib,is,ik)=wfnsrt_d(:,isrt(ib-nbinp))
        enddo ! ib
      enddo ! is
    enddo ! ik
    
    SAFE_DEALLOCATE(ensrt)
    SAFE_DEALLOCATE(isrt)
    SAFE_DEALLOCATE(wfnsrt_d)
    
  endif ! sapo_energy_sort
  
  call timacc(9,2)
  
  call timacc(10,1)

!------------------------
! Construct and diagonalize Hamiltonian, update plane wave
! and auxiliary wavefunctions and eigenvalues

  if (sapo_hamiltonian) then
    
    SAFE_ALLOCATE(vsc_d, (ng_l,ns))
    if (vkbg_flag) then
      SAFE_ALLOCATE(vkb_d, (ngk_l,nkb,ns,nk))
    endif ! vkbg_flag
    
    call timacc(11,1)
    
    if (peinf%inode.eq.0) write(6,813) trim(vscg_input_file)

    SAFE_ALLOCATE(gvectmp, (3, ng))
    call read_binary_gvectors(10, ng, ng, gvectmp, bcast = .false.)
    if (peinf%inode .eq. 0) then
      info=0
      do ig=1,ng
        if (gvectmp(1,ig).ne.gvec(1,ig).or.gvectmp(2,ig).ne. &
          gvec(2,ig).or.gvectmp(3,ig).ne.gvec(3,ig)) info=-1
      enddo
      if (info .ne. 0) then
        call die('inconsistent G-vectors in potential file ' // trim(vscg_input_file), only_root_writes = .true.)
      endif
    endif ! peinf%inode
    SAFE_DEALLOCATE(gvectmp)

    SAFE_ALLOCATE(vsc, (ng_g, ns))
    call read_binary_data(10, ng, ng_g, ns, vsc, bcast = .false.)
    do is=1,ns
      do ig=ng+1,ng_g
        vsc(ig,is)=ZERO
      enddo
#ifdef MPI
      !call MPI_Barrier(MPI_COMM_WORLD,mpierr) ! Uncomment if you have MPI errors
      call MPI_Scatter(vsc(:,is),ng_l,MPI_SCALAR, &
        vsc_d(:,is),ng_l,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#else
      vsc_d(:,is)=vsc(:,is)
#endif
    enddo
    SAFE_DEALLOCATE(vsc)
    if (ierr.ne.0) then
      call die('failed to read potential file ' // trim(vscg_input_file), only_root_writes = .true.)
    endif

    if (peinf%inode.eq.0) then
      call close_file(10)
    endif ! peinf%inode.eq.0
    
    if (vkbg_flag) then
      
      if (peinf%inode.eq.0) write(6,813) trim(vkbg_input_file)
      
      SAFE_ALLOCATE(isort, (ngk_g))
      SAFE_ALLOCATE(gvectmp, (3, ngkmax))
      SAFE_ALLOCATE(vkb, (ngk_g, 1))
      do ik=1,nk
#ifdef MPI
        !call MPI_Barrier(MPI_COMM_WORLD,mpierr) ! Uncomment if you have MPI errors
        call MPI_Gather(isort_d(:,ik),ngk_l,MPI_INTEGER, &
          isort,ngk_l,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#else
        isort(:)=isort_d(:,ik)
#endif
        call read_binary_gvectors(11, ngk(ik), ngkmax, gvectmp, bcast = .false.)
        
        if (peinf%inode .eq. 0) then
          do ig=1,ngk(ik)
            info=0
            iout=((gvectmp(1,ig)+FFTgrid(1)/2)*FFTgrid(2)+gvectmp(2,ig)+FFTgrid(2)/2)* &
              FFTgrid(3)+gvectmp(3,ig)+FFTgrid(3)/2+1
            if (iout.ge.1.and.iout.le.nFFTgridpts) then
              iout=index_vec(iout)
              if (iout.ge.1.and.iout.le.ng) then
                if (gvectmp(1,ig).ne.gvec(1,iout).or.gvectmp(2,ig).ne. &
                  gvec(2,iout).or.gvectmp(3,ig).ne.gvec(3,iout)) info=-1
              else
                info=-1
              endif
            else
              info=-1
            endif
            if (iout.ne.isort(ig)) info=-1
            ierr=ierr+info
          enddo ! ig
          if (info .ne. 0) then
            call die('inconsistent G-vectors in potential file ' // trim(vkbg_input_file), only_root_writes = .true.)
          endif
        endif ! peinf%inode
        
        do is=1,ns
          do ikb=1,nkb
            vkb(:,:)=(0.0d0,0.0d0)
            call read_binary_complex_data(11, ngk(ik), ngk_g, 1, vkb, bcast = .false.)
#ifdef MPI
            !call MPI_Barrier(MPI_COMM_WORLD,mpierr) ! Uncomment if you have MPI errors
            call MPI_Scatter(vkb(:,1),ngk_l,MPI_COMPLEX_DPC, &
              vkb_d(:,ikb,is,ik),ngk_l,MPI_COMPLEX_DPC, &
              0,MPI_COMM_WORLD,mpierr)
#else
            vkb_d(:,ikb,is,ik)=vkb(:,1)
#endif
          enddo ! ikb
        enddo ! is
      enddo ! ik
      SAFE_DEALLOCATE(isort)
      SAFE_DEALLOCATE(gvectmp)
      SAFE_DEALLOCATE(vkb)
      
      if (peinf%inode.eq.0) then
        call close_file(11)
      endif ! peinf%inode.eq.0
      
    endif ! vkbg_flag
        
    call timacc(11,2)
    
    if (peinf%inode.eq.0) write(6,814)
    
    fblas=fblas_ham
    fvkbg=vkbg_flag
    if (sapo_do_all_bands) then
      nbstart=1
    else
      nbstart=nbinp+1
    endif
    nbend=nbtot

    call hdiag(fparafft,fparala,fblas,fvkbg,nk,ns,nbstart, &
      nbend,nbmax,ng,ng_l,ng_g,ngk_l,ngk_g,nat,nsp,nkb,nhm,FFTgrid, &
      ngk,ityp,nh,gvec,isort_d,celvol,bdot,kpt,vsc_d,deeq,vkb_d, &
      en,wfn_d,sapo_ham_nrestart,sapo_ham_ndiag,sapo_ham_ndim, &
      sapo_ham_tol,sapo_ham_resinc)
    
    SAFE_DEALLOCATE(vsc_d)
    if (vkbg_flag) then
      SAFE_DEALLOCATE(ityp)
      SAFE_DEALLOCATE(nh)
      SAFE_DEALLOCATE(deeq)
      SAFE_DEALLOCATE(vkb_d)
    endif ! vkbg_flag
    
  endif ! sapo_hamiltonian
  
  call timacc(10,2)

!------------------------
! Deallocate index_vec indices

  SAFE_DEALLOCATE(index_vec)
  
  call timacc(21,1)
  
!------------------------
! Keep plane wave and auxiliary eigenvalues in a given energy range

  if (sapo_energy_range) then
    
    if (peinf%inode.eq.0) write(6,815)
    
    ib1=0
    ib2=0
    do ib=nbinp+1,nbtot
      if (any(en(ib,1:nk,1:ns).lt.sapo_energy_min-TOL_Zero)) then
        ib1=ib1+1
      endif
      if (all(en(ib,1:nk,1:ns).gt.sapo_energy_max+TOL_Zero)) then
        ib2=ib2+1
      endif
    enddo ! ib
    
    if (ib1.gt.0) then
      if (peinf%inode.eq.0) then
        write(0,941)
      endif
      do ik=1,nk
        do is=1,ns
          do ib=nbinp+1,nbtot-ib1-ib2
            en(ib,ik,is)=en(ib+ib1,ik,is)
          enddo ! ib
          do ib=nbinp+1,nbtot-ib1-ib2
            wfn_d(:,ib,is,ik)=wfn_d(:,ib+ib1,is,ik)
          enddo ! ib
        enddo ! is
      enddo ! ik
    endif ! ib1.gt.0
    
    if (ib1.gt.0.or.ib2.gt.0) then
      nbtot=nbtot-ib1-ib2
      if (.not.aux_flag) then
        nbpw=nbpw-ib1-ib2
      elseif (.not.pw_flag) then
        nbaux=nbaux-ib1-ib2
      else
        nbpw=0
        nbaux=0
      endif
      if (peinf%inode.eq.0) then
        write(s1,101)ib1+ib2
        write(s2,103)sapo_energy_min
        write(s3,103)sapo_energy_max
        write(s4,101)nbtot
        write(6,816)TRUNC(s1),TRUNC(s2),TRUNC(s3),TRUNC(s4)
      endif ! peinf%inode.eq.0
    else ! ib1.gt.0.or.ib2.gt.0
      if (peinf%inode.eq.0) write(6,*)
    endif ! ib1.gt.0.or.ib2.gt.0
    
  endif ! sapo_energy_range

  call timacc(21,2)

  call timacc(2,1)

!------------------------
! Reallocate occupations

  SAFE_ALLOCATE(occ_x, (nb, nk, ns))
  occ_x(1:nb, 1:nk, 1:ns) = occ(1:nb, 1:nk, 1:ns)
  SAFE_DEALLOCATE(occ)
  SAFE_ALLOCATE(occ, (nbtot, nk, ns))
  occ(1:nb, 1:nk, 1:ns) = occ_x(1:nb, 1:nk, 1:ns)
  SAFE_DEALLOCATE(occ_x)
  occ(nb+1:nbtot, 1:nk, 1:ns) = 0.0d0
  
  call timacc(2,2)

  if (out_flag) then
  
    call timacc(22,1)

!------------------------
! Write header of wavefunction file

    if (peinf%inode.eq.0) then
      write(6,803) trim(wfng_output_file)
      sheader = 'WFN'
      iflavor = 0
      call write_binary_header(9, sheader, iflavor, ns, ng, nsym, &
        cell_symmetry, nat, nk, nbtot, ngkmax, ecutrho, ecutwfn, FFTgrid, &
        kgrid, kshift, celvol, al, a, adot, recvol, bl, b, bdot, &
        rot, tau, atyp, apos, ngk, kw, kpt, ifmin, ifmax, en, occ)
      call write_binary_gvectors(9, ng, ng, gvec)
    endif

!------------------------
! Loop over k-points

    SAFE_ALLOCATE(wfn, (ngk_g,ns))
    SAFE_ALLOCATE(isort, (ngk_g))
  
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
      call write_binary_gvectors(9, ngk(ik), ng, gvec, gindex = isort)
        
!------------------------
! Gather band and spin wavefunctions for current k-point
! Write them to wavefunction file

      do ib=1,nbtot
        do is=1,ns
#ifdef MPI
          !call MPI_Barrier(MPI_COMM_WORLD,mpierr) ! Uncomment if you have MPI errors
          call MPI_Gather(wfn_d(:,ib,is,ik),ngk_l, &
            MPI_SCALAR,wfn(:,is),ngk_l, &
            MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#else
          wfn(:,is)=wfn_d(:,ib,is,ik)
#endif
        enddo
        call write_binary_data(9, ngk(ik), ngk_g, ns, wfn)
      enddo
    
    enddo
  
    SAFE_DEALLOCATE(wfn)
    SAFE_DEALLOCATE(isort)
  
    call timacc(22,2)

  endif ! out_flag
  
  call timacc(2,1)

!------------------------
! Close files

  if (peinf%inode.eq.0) then
    call close_file(7)
  endif
  
  if (aux_flag.and.peinf%inode.eq.0) then
    call close_file(8)
  endif
  
  if (out_flag.and.peinf%inode.eq.0) then
    call close_file(9)
  endif
  
  call timacc(2,2)
  
  call timacc(23,1)

!------------------------
! Check the orthonormality of the generated wavefunctions

  if (sapo_check_norm) then
    
    SAFE_ALLOCATE(wfnsp, (nbtot))
#ifdef MPI
    SAFE_ALLOCATE(wfnspdummy, (nbtot))
#endif
    
    dotmax=0.0d0
    do ik=1,nk
      do is=1,ns
        if (peinf%inode.eq.0) then
          write(s1,101)ik
          write(s2,101)is
          write(6,301)TRUNC(s1),TRUNC(s2)
        endif
        do ib=1,nbtot
          wfnsp(:)=ZERO
          do jb=1,nbtot
            z=ZERO
            do ig=1,ngk_l
              z=z+MYCONJG(wfn_d(ig,jb,is,ik))*wfn_d(ig,ib,is,ik)
            enddo ! ig
            wfnsp(jb)=z
          enddo ! jb
#ifdef MPI
          wfnspdummy(:)=wfnsp(:)
          call MPI_Allreduce(wfnspdummy,wfnsp,nbtot, &
            MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
          do jb=1,nbtot
            if (jb.eq.ib) then
              dotcur = abs(wfnsp(jb)-ONE)
            else
              dotcur = abs(wfnsp(jb))
            endif
            if (dotmax.lt.dotcur) dotmax = dotcur
            if (dotcur.gt.TOL_Small.and.peinf%inode.eq.0) then
              write(s1,101)ib
              write(s2,101)jb
#ifdef CPLX
              write(s3,102)dble(wfnsp(jb))
              write(s4,102)IMAG(wfnsp(jb))
              write(0,303)TRUNC(s1),TRUNC(s2),TRUNC(s3),TRUNC(s4)
#else
              write(s3,102)wfnsp(jb)
              write(0,302)TRUNC(s1),TRUNC(s2),TRUNC(s3)
#endif
            endif
          enddo ! jb
        enddo ! ib
      enddo ! is
    enddo ! ik
    if (peinf%inode.eq.0) write(6,304) dotmax
    
    SAFE_DEALLOCATE(wfnsp)
#ifdef MPI
    SAFE_DEALLOCATE(wfnspdummy)
#endif
    
  endif ! sapo_check_norm
  
  call timacc(23,2)
  
  call timacc(24,1)

!------------------------
! Plot the energy eigenvalues of the generated wavefunctions
! with respect to band index

  if (sapo_eigenvalue) then
    ! determine range of bands
    if (sapo_plot_bandmin.eq.0.and.sapo_plot_bandmax.eq.0) then
      ib1=1
      ib2=nbtot
    else
      ib1=max(1,sapo_plot_bandmin)
      ib2=min(nbtot,sapo_plot_bandmax)
    endif
    !
    do ik=1,nk
      if (sapo_plot_kpoint.ne.0.and.ik.ne.sapo_plot_kpoint) cycle
      do is=1,ns
        if (sapo_plot_spin.ne.0.and.is.ne.sapo_plot_spin) cycle
        if (peinf%inode.eq.0) then
          write(s1,101)ik
          write(s2,101)is
          write(6,401)TRUNC(s1),TRUNC(s2)
        endif
        if (peinf%inode.eq.0) then
          write(s1,101)ik
          write(s2,101)is
          write(tmpstr,402)TRUNC(s1),TRUNC(s2)
          call open_file(9,file=tmpstr,status='replace',form='formatted')
          if (ierr.eq.0) then
            write(9,403)TRUNC(s1),TRUNC(s2)
            do ib=ib1,ib2
              write(9,404)ib,en(ib,ik,is)*RYD
            enddo ! ib
          endif
          if (ierr.eq.0) call close_file(9)
        endif
#ifdef MPI
        call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
        if (ierr.ne.0) then
          call die('failed to write eigenvalue file ' // trim(tmpstr), only_root_writes = .true.)
        endif
      enddo ! is
    enddo ! ik
    !
    if (peinf%inode.eq.0) write(6,*)
  endif ! sapo_eigenvalue
  
  call timacc(24,2)
  
  call timacc(25,1)

!------------------------
! Plot projections of the generated wavefunctions
! with respect to band index

  if (iand(sapo_projection,1).ne.0) then
    ! determine range of bands
    if (sapo_plot_bandmin.eq.0.and.sapo_plot_bandmax.eq.0) then
      ib1=1
      ib2=nbtot
    else
      ib1=max(1,sapo_plot_bandmin)
      ib2=min(nbtot,sapo_plot_bandmax)
    endif
    ! determine range of plane waves
    if (sapo_plot_pwmin.eq.0.and.sapo_plot_pwmax.eq.0) then
      ipw1=1
      ipw2=ngkmax
    else
      ipw1=max(1,sapo_plot_pwmin)
      ipw2=min(ngkmax,sapo_plot_pwmax)
    endif
    !
    SAFE_ALLOCATE(ekin, (ng))
    SAFE_ALLOCATE(isrt, (ng))
    SAFE_ALLOCATE(rdeg, (2,ng))
    SAFE_ALLOCATE(ideg, (ng))
    SAFE_ALLOCATE(isrti, (ng))
    SAFE_ALLOCATE(comp, (ib1:ib2))
    SAFE_ALLOCATE(compdummy, (ib1:ib2))
    !
    do ik=1,nk
      if (sapo_plot_kpoint.ne.0.and.ik.ne.sapo_plot_kpoint) cycle
      ! sort G-vectors by kinetic energies
      call groupk(kpt(:,ik),nsyml,rotl,nsymk,rotk)
      do ig=1,ng
        do i=1,3
          kvec(i)=kpt(i,ik)+dble(gvec(i,ig))
        enddo ! i
        ekin(ig)=DOT_PRODUCT(kvec,MATMUL(bdot,kvec))
      enddo ! ig
      call sortrx(ng,ekin,isrt,gvec=gvec)
      pwsym=.true.
      call pwdeg(pwsym,nsymk,rotk,ng,gvec,ekin,isrt,ndeg,mdeg,rdeg,ideg)
      ! invert isrt array
      isrti(:)=0
      do ig=1,ng
        isrti(isrt(ig))=ig
      enddo ! ig
      do is=1,ns
        if (sapo_plot_spin.ne.0.and.is.ne.sapo_plot_spin) cycle
        do ipw=ipw1,ipw2
          if (peinf%inode.eq.0) then
            write(s1,101)ik
            write(s2,101)is
            write(s3,101)ipw
            write(6,411)TRUNC(s1),TRUNC(s2),TRUNC(s3)
          endif
          comp(:)=0.0d0
          do ig=1,ngk_l
            jg=isort_d(ig,ik)
            if (jg.ne.0) then
              jg=isrti(jg)
              if (jg.eq.ipw) &
                comp(ib1:ib2)=abs(wfn_d(ig,ib1:ib2,is,ik))
            endif
          enddo ! ig
#ifdef MPI
          compdummy(:)=comp(:)
          call MPI_Reduce(compdummy,comp,ib2-ib1+1, &
            MPI_REAL_DP,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
#endif
          if (peinf%inode.eq.0) then
            write(s1,101)ik
            write(s2,101)is
            write(s3,101)ipw
            write(tmpstr,412)TRUNC(s1),TRUNC(s2),TRUNC(s3)
            call open_file(9,file=tmpstr,status='replace',form='formatted')
            if (ierr.eq.0) then
              write(9,413)TRUNC(s1),TRUNC(s2),TRUNC(s3)
              do ib=ib1,ib2
                write(9,414)ib,comp(ib)
              enddo ! ib
            endif
            if (ierr.eq.0) call close_file(9)
          endif
#ifdef MPI
          call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
          if (ierr.ne.0) then
            call die('failed to write projection file ' // trim(tmpstr), only_root_writes = .true.)
          endif
        enddo ! ipw
      enddo ! is
    enddo ! ik
    !
    !
    SAFE_DEALLOCATE(ekin)
    SAFE_DEALLOCATE(isrt)
    SAFE_DEALLOCATE(rdeg)
    SAFE_DEALLOCATE(ideg)
    SAFE_DEALLOCATE(isrti)
    SAFE_DEALLOCATE(comp)
    SAFE_DEALLOCATE(compdummy)
    !
    if (peinf%inode.eq.0) write(6,*)
  endif ! iand(sapo_projection,1).ne.0
  
!------------------------
! Plot projections of the generated wavefunctions
! with respect to plane wave

  if (iand(sapo_projection,2).ne.0) then
    ! determine range of bands
    if (sapo_plot_bandmin.eq.0.and.sapo_plot_bandmax.eq.0) then
      ib1=1
      ib2=nbtot
    else
      ib1=max(1,sapo_plot_bandmin)
      ib2=min(nbtot,sapo_plot_bandmax)
    endif
    ! determine range of plane waves
    if (sapo_plot_pwmin.eq.0.and.sapo_plot_pwmax.eq.0) then
      ipw1=1
      ipw2=ngkmax
    else
      ipw1=max(1,sapo_plot_pwmin)
      ipw2=min(ngkmax,sapo_plot_pwmax)
    endif
    !
    SAFE_ALLOCATE(ekin, (ng))
    SAFE_ALLOCATE(isrt, (ng))
    SAFE_ALLOCATE(rdeg, (2,ng))
    SAFE_ALLOCATE(ideg, (ng))
    SAFE_ALLOCATE(isrti, (ng))
    SAFE_ALLOCATE(comp, (ipw1:ipw2))
    SAFE_ALLOCATE(compdummy, (ipw1:ipw2))
    !
    do ik=1,nk
      if (sapo_plot_kpoint.ne.0.and.ik.ne.sapo_plot_kpoint) cycle
      ! sort G-vectors by kinetic energies
      call groupk(kpt(:,ik),nsyml,rotl,nsymk,rotk)
      do ig=1,ng
        do i=1,3
          kvec(i)=kpt(i,ik)+dble(gvec(i,ig))
        enddo ! i
        ekin(ig)=DOT_PRODUCT(kvec,MATMUL(bdot,kvec))
      enddo ! ig
      call sortrx(ng,ekin,isrt,gvec=gvec)
      pwsym=.true.
      call pwdeg(pwsym,nsymk,rotk,ng,gvec,ekin,isrt,ndeg,mdeg,rdeg,ideg)
      ! invert isrt array
      isrti(:)=0
      do ig=1,ng
        isrti(isrt(ig))=ig
      enddo ! ig
      do is=1,ns
        if (sapo_plot_spin.ne.0.and.is.ne.sapo_plot_spin) cycle
        do ib=ib1,ib2
          if (peinf%inode.eq.0) then
            write(s1,101)ik
            write(s2,101)is
            write(s3,101)ib
            write(6,421)TRUNC(s1),TRUNC(s2),TRUNC(s3)
          endif
          comp(:)=0.0d0
          do ig=1,ngk_l
            jg=isort_d(ig,ik)
            if (jg.ne.0) then
              jg=isrti(jg)
              if (jg.ge.ipw1.and.jg.le.ipw2) &
                comp(jg)=abs(wfn_d(ig,ib,is,ik))
            endif
          enddo ! ig
#ifdef MPI
          compdummy(:)=comp(:)
          call MPI_Reduce(compdummy,comp,ipw2-ipw1+1, &
            MPI_REAL_DP,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
#endif
          if (peinf%inode.eq.0) then
            write(s1,101)ik
            write(s2,101)is
            write(s3,101)ib
            write(tmpstr,422)TRUNC(s1),TRUNC(s2),TRUNC(s3)
            call open_file(9,file=tmpstr,status='replace',form='formatted')
            if (ierr.eq.0) then
              write(9,423)TRUNC(s1),TRUNC(s2),TRUNC(s3)
              do ipw=ipw1,ipw2
                write(9,424)ipw,comp(ipw)
              enddo ! ipw
            endif
            if (ierr.eq.0) call close_file(9)
          endif
#ifdef MPI
          call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
          if (ierr.ne.0) then
            call die('failed to write projection file ' // trim(tmpstr), only_root_writes = .true.)
          endif
        enddo ! ib
      enddo ! is
    enddo ! ik
    !
    SAFE_DEALLOCATE(ekin)
    SAFE_DEALLOCATE(isrt)
    SAFE_DEALLOCATE(rdeg)
    SAFE_DEALLOCATE(ideg)
    SAFE_DEALLOCATE(isrti)
    SAFE_DEALLOCATE(comp)
    SAFE_DEALLOCATE(compdummy)
    !
    if (peinf%inode.eq.0) write(6,*)
  endif ! iand(sapo_projection,2).ne.0

!------------------------
! Plot projections of the generated wavefunctions
! with respect to band index and plane wave

  if (iand(sapo_projection,4).ne.0) then
    ! determine range of bands
    if (sapo_plot_bandmin.eq.0.and.sapo_plot_bandmax.eq.0) then
      ib1=1
      ib2=nbtot
    else
      ib1=max(1,sapo_plot_bandmin)
      ib2=min(nbtot,sapo_plot_bandmax)
    endif
    ! determine range of plane waves
    if (sapo_plot_pwmin.eq.0.and.sapo_plot_pwmax.eq.0) then
      ipw1=1
      ipw2=ngkmax
    else
      ipw1=max(1,sapo_plot_pwmin)
      ipw2=min(ngkmax,sapo_plot_pwmax)
    endif
    !
    SAFE_ALLOCATE(ekin, (ng))
    SAFE_ALLOCATE(isrt, (ng))
    SAFE_ALLOCATE(rdeg, (2,ng))
    SAFE_ALLOCATE(ideg, (ng))
    SAFE_ALLOCATE(isrti, (ng))
    SAFE_ALLOCATE(proj, (ipw1:ipw2,ib1:ib2))
    SAFE_ALLOCATE(projdummy, (ipw1:ipw2,ib1:ib2))
    !
    do ik=1,nk
      if (sapo_plot_kpoint.ne.0.and.ik.ne.sapo_plot_kpoint) cycle
      ! sort G-vectors by kinetic energies
      call groupk(kpt(:,ik),nsyml,rotl,nsymk,rotk)
      do ig=1,ng
        do i=1,3
          kvec(i)=kpt(i,ik)+dble(gvec(i,ig))
        enddo ! i
        ekin(ig)=DOT_PRODUCT(kvec,MATMUL(bdot,kvec))
      enddo ! ig
      call sortrx(ng,ekin,isrt,gvec=gvec)
      pwsym=.true.
      call pwdeg(pwsym,nsymk,rotk,ng,gvec,ekin,isrt,ndeg,mdeg,rdeg,ideg)
      ! invert isrt array
      isrti(:)=0
      do ig=1,ng
        isrti(isrt(ig))=ig
      enddo ! ig
      ! calculate projections
      do is=1,ns
        if (sapo_plot_spin.ne.0.and.is.ne.sapo_plot_spin) cycle
        if (peinf%inode.eq.0) then
          write(s1,101)ik
          write(s2,101)is
          write(6,431)TRUNC(s1),TRUNC(s2)
        endif
        proj(:,:)=0.0d0
        do ib=ib1,ib2
          do ig=1,ngk_l
            jg=isort_d(ig,ik)
            if (jg.ne.0) then
              ipw=isrti(jg)
              if (ipw.ge.ipw1.and.ipw.le.ipw2) &
                proj(ipw,ib)=abs(wfn_d(ig,ib,is,ik))
            endif
          enddo ! ig
        enddo ! ib
#ifdef MPI
        projdummy(:,:)=proj(:,:)
        call MPI_Reduce(projdummy(1,1),proj(1,1),(ipw2-ipw1+1)*(ib2-ib1+1), &
          MPI_REAL_DP,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
#endif
        if (peinf%inode.eq.0) then
          write(s1,101)ik
          write(s2,101)is
          write(s3,101)ib1
          write(s4,101)ib2
          write(s5,101)ipw1
          write(s6,101)ipw2
          write(tmpstr,432)TRUNC(s1),TRUNC(s2)
          call open_file(9,file=tmpstr,status='replace',form='formatted')
          write(9,433)TRUNC(s1),TRUNC(s2),TRUNC(s3),TRUNC(s4), &
            TRUNC(s5),TRUNC(s6)
          if (ierr.eq.0) then
            do ipw=ipw1,ipw2
              write(9,434)(proj(ipw,ib),ib=ib1,ib2)
            enddo ! ipw
          endif
          if (ierr.eq.0) call close_file(9)
        endif
#ifdef MPI
        call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
        if (ierr.ne.0) then
          call die('failed to write projection file ' // trim(tmpstr), only_root_writes = .true.)
        endif
      enddo ! is
    enddo ! ik
    !
    SAFE_DEALLOCATE(ekin)
    SAFE_DEALLOCATE(isrt)
    SAFE_DEALLOCATE(rdeg)
    SAFE_DEALLOCATE(ideg)
    SAFE_DEALLOCATE(isrti)
    SAFE_DEALLOCATE(proj)
    SAFE_DEALLOCATE(projdummy)
    !
    if (peinf%inode.eq.0) write(6,*)
  endif ! iand(sapo_projection,4).ne.0
  
  call timacc(25,2)
  
  call timacc(26,1)

!------------------------
! Plot the squared absolute values of amplitudes of the generated
! wavefunctions with respect to kinetic energies of plane waves

  if (sapo_amplitude) then
    ! determine range of bands
    if (sapo_plot_bandmin.eq.0.and.sapo_plot_bandmax.eq.0) then
      ib1=1
      ib2=nbtot
    else
      ib1=max(1,sapo_plot_bandmin)
      ib2=min(nbtot,sapo_plot_bandmax)
    endif
    !
    SAFE_ALLOCATE(ekin, (ng))
    SAFE_ALLOCATE(amplx, (sapo_ampl_num))
    SAFE_ALLOCATE(amply, (sapo_ampl_num))
    SAFE_ALLOCATE(ampldummy, (sapo_ampl_num))
    !
    do i=1,sapo_ampl_num
      amplx(i)=sapo_ampl_del*dble(i-1)
    enddo ! i
    !
    do ik=1,nk
      if (sapo_plot_kpoint.ne.0.and.ik.ne.sapo_plot_kpoint) cycle
      ! compute kinetic energies of G-vectors
      do ig=1,ng
        do i=1,3
          kvec(i)=kpt(i,ik)+dble(gvec(i,ig))
        enddo ! i
        ekin(ig)=DOT_PRODUCT(kvec,MATMUL(bdot,kvec))
      enddo ! ig
      do is=1,ns
        if (sapo_plot_spin.ne.0.and.is.ne.sapo_plot_spin) cycle
        do ib=ib1,ib2
          if (peinf%inode.eq.0) then
            write(s1,101)ik
            write(s2,101)is
            write(s3,101)ib
            write(6,441)TRUNC(s1),TRUNC(s2),TRUNC(s3)
          endif
          amply(:)=0.0d0
          do ig=1,ngk_l
            jg=isort_d(ig,ik)
            if (jg.ne.0) then
              x=abs(wfn_d(ig,ib,is,ik))**2
              do i=1,sapo_ampl_num
                amply(i)=amply(i)+x/(sqrt(2.0d0*PI_D)*sapo_ampl_brd)* &
                  exp(-(amplx(i)-ekin(jg))**2/(2.0d0*sapo_ampl_brd**2))
              enddo ! i
            endif
          enddo ! ig
#ifdef MPI
          ampldummy(:)=amply(:)
          call MPI_Reduce(ampldummy,amply,sapo_ampl_num, &
            MPI_REAL_DP,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
#endif
          if (peinf%inode.eq.0) then
            write(s1,101)ik
            write(s2,101)is
            write(s3,101)ib
            write(tmpstr,442)TRUNC(s1),TRUNC(s2),TRUNC(s3)
            call open_file(9,file=tmpstr,status='replace',form='formatted')
            if (ierr.eq.0) then
              write(9,443)TRUNC(s1),TRUNC(s2),TRUNC(s3)
              do i=1,sapo_ampl_num
                write(9,444)amplx(i),amply(i)
              enddo ! i
            endif
            if (ierr.eq.0) call close_file(9)
          endif
#ifdef MPI
          call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
          if (ierr.ne.0) then
            call die('failed to write amplitude file ' // trim(tmpstr), only_root_writes = .true.)
          endif
        enddo ! ib
      enddo ! is
    enddo ! ik
    !
    SAFE_DEALLOCATE(ekin)
    SAFE_DEALLOCATE(amplx)
    SAFE_DEALLOCATE(amply)
    SAFE_DEALLOCATE(ampldummy)
    !
    if (peinf%inode.eq.0) write(6,*)
  endif ! sapo_amplitude
  
  call timacc(26,2)
  
  call timacc(2,1)

!------------------------
! Deallocate arrays

  SAFE_DEALLOCATE_P(atyp)
  SAFE_DEALLOCATE_P(apos)
  SAFE_DEALLOCATE(kw)
  SAFE_DEALLOCATE(kpt)
  SAFE_DEALLOCATE(ifmin)
  SAFE_DEALLOCATE(ifmax)
  SAFE_DEALLOCATE(occ)
  SAFE_DEALLOCATE(en)
  SAFE_DEALLOCATE(gvec)
  SAFE_DEALLOCATE(ngk)
  SAFE_DEALLOCATE(wfn_d)
  SAFE_DEALLOCATE(isort_d)
  
  call timacc(2,2)
  
  call timacc(1,2)

!------------------------
! Time Accounting

  routnam(1)='TOTAL:'
  routnam(2)='INIT:'
  routnam(3)='INPUT:'
  routnam(4)='INPUT_AUX:'
  routnam(5)='GEN_PW:'
  routnam(6)='RANDOM:'
  routnam(7)='DROP:'
  routnam(8)='ORTHONORM:'
  routnam(9)='SORT:'
  routnam(10)='HAMILTONIAN:'
  routnam(11)='HAM (INPUT):'
  routnam(12)='HAM (COMM):'
  routnam(13)='HAM (PLAN):'
  routnam(14)='HAM (FFT):'
  routnam(15)='HAM (KIN):'
  routnam(16)='HAM (VSC):'
  routnam(17)='HAM (VKB):'
  routnam(18)='HAM (GATHER):'
  routnam(19)='HAM (DIAG):'
  routnam(20)='HAM (WFN):'
  routnam(21)='RANGE:'
  routnam(22)='OUTPUT:'
  routnam(23)='CHECK_NORM:'
  routnam(24)='PLOT_ENER:'
  routnam(25)='PLOT_PROJ:'
  routnam(26)='PLOT_AMPL:'
  
  if(peinf%inode.eq.0) then
    write(6,9000)
    do i=2,26
      if (.not.sapo_hamiltonian.and.i.ge.10.and.i.le.19) cycle
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
103 format(f16.3)
  
110 format(1x,'nsym =',1x,a,2x,'nk =',1x,a,2x,'ns =',1x,a,2x,'nb =', &
      1x,a,2x,'ng =',1x,a,2x,'ngkmax =',1x,a)
120 format(1x,'ng_l =',1x,a,2x,'ng_g =',1x,a,2x,'ngk_l =',1x,a,2x, &
      'ngk_g =',1x,a)
130 format(1x,'FFTgrid =',1x,'(',1x,a,1x,a,1x,a,1x,')',2x,'nFFTgridpts =',1x,a)
140 format(1x,'nbinp =',1x,a,2x,'nbpw =',1x,a,2x,'b2g =',1x,a,2x, &
      'nbaux =',1x,a,2x,'nbtot =',1x,a,/)
  
301 format(1x,'Checking orthonormality for ik =',1x,a,1x,'is =', &
      1x,a)
302 format(4x,'WARNING: ib =',1x,a,1x,'jb =',1x,a,1x, &
      'scalar product =',1x,a)
303 format(4x,'WARNING: ib =',1x,a,1x,'jb =',1x,a,1x, &
      'scalar product =',1x,'(',1x,a,1x,a,1x,')')
304 format(/,1x,'Maximum deviation from orthonormality:',/, &
      4x,'max ( < nk | mk > - delta_nm ) =',f20.15,/)
  
401 format(1x,'Generating eigenvalues for ik =',1x,a,1x,'is =',1x,a)
402 format(1x,'ener_k',a,'_s',a,'.dat')
403 format('# E(kpt=',a,',spin=',a,') (eV) vs band')
404 format(i6,f12.6)
  
411 format(1x,'Generating projections for ik =',1x,a,1x,'is =', &
      1x,a,1x,'ipw =',1x,a)
412 format(1x,'proj_k',a,'_s',a,'_pw',a,'.dat')
413 format('# |psi(kpt=',a,',spin=',a,',pw=',a,')| vs band')
414 format(i6,f10.6)
  
421 format(1x,'Generating projections for ik =',1x,a,1x,'is =', &
      1x,a,1x,'ib =',1x,a)
422 format(1x,'proj_k',a,'_s',a,'_b',a,'.dat')
423 format('# |psi(kpt=',a,',spin=',a,',band=',a,')| vs pw')
424 format(i6,f10.6)
  
431 format(1x,'Generating projections for ik =',1x,a,1x,'is =',1x,a)
432 format(1x,'proj_k',a,'_s',a,'.dat')
433 format('# |psi(kpt=',a,',spin=',a,')| vs', &
      1x,'(row/band=',a,':',a,',col/pw=',a,':',a,')')
434 format(10000f10.6)
  
441 format(1x,'Generating amplitudes for ik =',1x,a,1x,'is =', &
      1x,a,1x,'ib =',1x,a)
442 format(1x,'ampl_k',a,'_s',a,'_b',a,'.dat')
443 format('# |psi(kpt=',a,',spin=',a,',band=',a,')|^2', &
      1x,'vs E_kin(pw) (Ry)')
444 format(2f10.6)
  
601 format(1x,'overlap_dft_k',a,'_s',a,'.dat')
602 format(1x,'overlap_aux_k',a,'_s',a,'.dat')
603 format(i6,f12.6)
  
801 format(1x,'Reading parameters from file',1x,a,/)
802 format(1x,'Reading wavefunctions from file',1x,a,/)
803 format(1x,'Writing wavefunctions to file',1x,a,/)
804 format(1x,'Generating plane waves',/)
805 format(1x,'Generating plane waves and decomposing them into irreducible representations',/,&
      4x,'of the space group of the crystal',/)
806 format(1x,'Generating plane waves and decomposing them into irreducible representations',/,&
      4x,'of the space group of the Bravais lattice',/)
807 format(1x,'Applying random variation with amplitude',1x,a,/)
808 format(1x,'Generating completely random wavefunctions',/)
809 format(1x,'Checking overlap of PW & AUX wavefunctions with DFT wavefunctions')
810 format(1x,'Dropping',1x,a,1x,'PW & AUX states that have a large overlap with DFT states',/,&
      3x,'Keeping',1x,a,1x,'DFT & PW & AUX states',/)
811 format(1x,'Orthonormalizing PW & AUX wavefunctions w.r.t. DFT',/)
812 format(1x,'Sorting PW & AUX wavefunctions by eigenvalues',/)
813 format(1x,'Reading potential from file',1x,a,/)
814 format(1x,'Constructing and diagonalizing Hamiltonian',/)
815 format(1x,'Checking if new PW & AUX eigenvalues fall outside the energy range')
816 format(1x,'Dropping',1x,a,1x,'PW & AUX states outside the energy range of',1x,a,1x,'to',1x,a,1x,'Ry',/,&
      3x,'Keeping',1x,a,1x,'DFT & PW & AUX states',/)
  
921 format(1x,'WARNING:',1x,'cannot match kinetic energies if starting from 1st plane wave',/)
922 format(1x,'WARNING:',1x,'cannot use energy shift if matching kinetic energies',/)
923 format(1x,'WARNING:',1x,'cannot use symmetries if there are no plane waves',/)
924 format(1x,'WARNING:',1x,'using AUX wavefunctions requires AUX input file',/)
925 format(1x,'WARNING:',1x,'random numbers overwrite irreducible representations',/)
926 format(1x,'WARNING:',1x,'random numbers overwrite auxiliary wavefunctions',/)
927 format(1x,'WARNING:',1x,'dropping PW & AUX states requires PW or AUX input',/)
928 format(1x,'WARNING:',1x,'set min PW & AUX to 0 since dropping states with a large overlap',/)
929 format(1x,'WARNING:',1x,'orthonormalizing wavefunctions requires PW or AUX input',/,3x,'or sapo_do_all_bands',/)
930 format(1x,'WARNING:',1x,'orthonormalizing PW & AUX together requires both PW & AUX input',/)
931 format(1x,'WARNING:',1x,'dropping PW & AUX states incompatible with not ordered',/,&
      3x,'in energy orthonormalization',/)
932 format(1x,'WARNING:',1x,'correcting PW & AUX eigenvalues requires PW or AUX input',/)
933 format(1x,'WARNING:',1x,'sorting PW & AUX eigenvalues requires PW or AUX input',/)
934 format(1x,'WARNING:',1x,'constructing Hamiltonian requires local potential input file',/)
935 format(1x,'WARNING:',1x,'constructing Hamiltonian requires PW or AUX input',/,3x,'or sapo_do_all_bands',/)
936 format(1x,'WARNING:',1x,'setting the range of PW & AUX eigenvalues requires sorting them or diagonalizing Hamiltonian',/)
937 format(1x,'WARNING:',1x,'randomizing wavefunctions requires PW or AUX input',/,3x,'or sapo_do_all_bands',/)
  
941 format(1x,'WARNING:',1x,'lowest PW & AUX state overlaps in energy with highest DFT state',/)

end program sapo
