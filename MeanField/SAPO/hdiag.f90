
!===============================================================================
!
! Routines:
!
! 1. hdiag()             Originally By gsm      Last Modified 9/3/2010 (gsm)
!
!    Computes the residual vectors, diagonalizes the subspace Hamiltonian,
!    or iteratively diagonalizes the Hamiltonian
!
! 2. hpsi()              Originally By gsm      Last Modified 9/3/2010 (gsm)
!
!    Applies the Hamiltonian to a wavefunction
!
! 3. hsub()              Originally By gsm      Last Modified 9/3/2010 (gsm)
!
!    Constructs and diagonalizes the subspace Hamiltonian
!
! 4. dsub()              Originally By gsm      Last Modified 9/3/2010 (gsm)
!
!    Constructs the diagonal subspace Hamiltonian and diagonal eigenvectors
!
! 5. rpsi()              Originally By gsm      Last Modified 9/3/2010 (gsm)
!
!    Rotates wavefunctions
!
! 6. npsi()              Originally By gsm      Last Modified 9/3/2010 (gsm)
!
!    Computes norms of wavefunctions and optionally normalizes wavefunctions
!
! 7. calc_precond()      Originally By gsm      Last Modified 9/3/2010 (gsm)
!
!    Computes a preconditioner for a component of a residual vector
!
!===============================================================================

#include "f_defs.h"

module hdiag_m

  use global_m
  use blas_m
  use lapack_m
  use scalapack_m
  use fftw_m
  use fft_parallel_m
  use misc_m
  use pw_m

  implicit none

  private

  public :: hdiag

contains

subroutine hdiag(fparafft,fparala,fblas,fvkbg,nk,ns,nbstart, &
  nbend,nbmax,ng,ng_l,ng_g,ngk_l,ngk_g,nat,nsp,nkb,nhm,FFTgrid, &
  ngk,ityp,nh,gvec,isort_d,celvol,bdot,kpt,vsc_d,deeq,vkb_d, &
  en,wfn_d,sapo_ham_nrestart,sapo_ham_ndiag,sapo_ham_ndim, &
  sapo_ham_tol,sapo_ham_resinc)

! This subroutine computes the residual vectors, diagonalizes the
! subspace Hamiltonian, or iteratively diagonalizes the Hamiltonian.
! On entrance takes initial eigenvalues and eigenvectors in arrays
! en(nbstart:nbend,nk,ns) and wfn_d(ngk_l,nbstart:nbend,ns,nk) set
! in SAPO. On exit returns final eigenvalues and eigenvectors in the
! same arrays
!
! For computing residual vectors set sapo_ham_ndiag = 0
!
! For subspace diagonalization set sapo_ham_ndiag = 1
!
! For iterative diagonalization set sapo_ham_ndiag > 1, sapo_ham_ndim > 1,
! sapo_ham_tol > 0.0d0
!
! If some k-points/spins didn`t converge during iterative diagonalization,
! try setting sapo_ham_nrestart > 1
!
! Some important variables for iterative diagonalization:
!   nvec - number of eigenvectors (wavefunctions)
!   nmat - leading dimension of the subspace Hamiltonian arrays
!          (nmat = nvec * sapo_ham_ndim)
!   nbas - number of the basis functions (nvec <= nbas <= nmat)
!   npre - number of the basis functions kept from previous iteration
!   nres - number of residual vectors added to the basis functions
!          (nres <= nvec)
!
! There is no point of setting sapo_ham_ndim > (sapo_ham_ndiag + 1).
! Since nbas is initialized to nvec and nbas is incremented by nres
! on each iteration, it requires at least (sapo_ham_ndim - 1)
! iterations for nbas to reach nmat
!
! There is no point of setting sapo_ham_ndiag > 1 if sapo_ham_ndim = 1.
! Subsequent iterations will be identical if the basis set is not
! allowed to expand

  logical, intent(in) :: fparafft,fparala,fblas,fvkbg
  integer, intent(in) :: nk,ns,nbstart,nbend,nbmax,ng,ng_l,ng_g, &
    ngk_l,ngk_g,nat,nsp,nkb,nhm
  integer, intent(in) :: FFTgrid(3)
  integer, intent(in) :: ngk(:) !< (nk)
  integer, intent(in) :: ityp(:) !< (nat)
  integer, intent(in) :: nh(:) !< (nsp)
  integer, intent(in) :: gvec(:,:) !< (3,ng)
  integer, intent(in) :: isort_d(:,:) !< (ngk_l,nk)
  real(DP), intent(in) :: celvol
  real(DP), intent(in) :: bdot(3,3)
  real(DP), intent(in) :: kpt(:,:) !< (3,nk)
  SCALAR, intent(in) :: vsc_d(:,:) !< (ng_l,ns)
  real(DP), intent(in) :: deeq(:,:,:,:) !< (nhm,nhm,nat,ns)
  complex(DPC), intent(in) :: vkb_d(:,:,:,:) !< (ngk_l,nkb,ns,nk)
  real(DP), intent(inout) :: en(:,:,:) !< (nbmax,nk,ns)
  SCALAR, intent(inout) :: wfn_d(:,:,:,:) !< (ngk_l,nbmax,ns,nk)
  integer, intent(in) :: sapo_ham_nrestart,sapo_ham_ndiag,sapo_ham_ndim
  real(DP), intent(in) :: sapo_ham_tol
  logical, intent(in) :: sapo_ham_resinc

  logical :: is_diag,is_precond,is_updatewfn,is_error
  integer :: i,j,k,ik,is,ib,im,jm,ig,irestart,idiag, &
    idiagmax,idiagtot,nvec,nbas,nmat,nres,npre,Nplane,Nrod, &
    m,nz,npr,npc,lwork,lrwork,liwork,liclustr,lgap,mrod,msph, &
    iat,isp,ih,ikb,ijkb0,info,nmpinode,dummyiwork(1), &
    dummyifail(1),dummyiclustr(1),desca(9)
  real(DP) :: minus_scale,plus_scale,plus_scale_v,abstol, &
    orfac,demax,detot,resmax,resold,restot,v_0,vdummy,precond, &
    amem,rmem,size_bool,size_int,size_real,size_cplx,size_scal, &
    impart_hpsi,impart_hsub,impart_ssub,impart_dummy, &
    dummyw(1),dummyrwork(1),dummygap(1),kvec(3)
  complex(DPC) :: ps1,ar
  SCALAR :: dummya(1,1),dummyb(1,1),dummyz(1,1),dummywork(3)
  character(len=16) :: s1,s2,s3
  
  real(DP), allocatable :: enbuf(:)
  real(DP), allocatable :: enold(:)
  SCALAR, allocatable :: wfnbuf_d(:,:)
  SCALAR, allocatable :: wfnnew_d(:,:)
  SCALAR, allocatable :: residual_d(:,:)
  complex(DPC), allocatable :: wfnc_d(:)
  complex(DPC), allocatable :: wfnc(:)
  complex(DPC), allocatable :: hpsic_d(:)
  SCALAR, allocatable :: hpsi_d(:,:)
  integer, allocatable :: isort(:)
  real(DP), allocatable :: ekin(:)
  real(DP), allocatable :: norm(:)
  real(DP), allocatable :: ndum(:)
  real(DP), allocatable :: h_gg(:)
  complex(DPC), allocatable :: betapsi(:)
  complex(DPC), allocatable :: dummykb(:)
  complex(DPC), allocatable :: ps(:)
  logical, allocatable :: is_converged(:,:)
  integer, pointer :: crod(:)
  integer, pointer :: csph(:)
  integer, pointer :: drod(:)
  integer, pointer :: dsph(:)
  integer, pointer :: irod(:,:)
  integer, pointer :: isph(:)
  complex(DPC), pointer :: buffer_rod(:)
  complex(DPC), pointer :: buffer_sph(:)
  complex(DPC), pointer :: fftbox_2D(:,:,:)
  complex(DPC), pointer :: fftbox_1D(:,:)
  complex(DPC), pointer :: buffer_2D(:,:,:)
  complex(DPC), pointer :: buffer_1D(:,:,:)
  complex(DPC), pointer :: vlr_d(:,:,:,:)
  integer, pointer :: inv_indx(:,:,:)
  complex(DPC), pointer :: fftbox_3D(:,:,:)
  complex(DPC), allocatable :: vlg(:)
  complex(DPC), pointer :: vlr(:,:,:,:)
  complex(DPC), allocatable :: vlg_d(:)
  SCALAR, allocatable :: hcol(:)
  SCALAR, allocatable :: scol(:)
  SCALAR, allocatable :: dummyc(:)
  SCALAR, allocatable :: work(:)
  real(DP), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  integer, allocatable :: iclustr(:)
  real(DP), allocatable :: gap(:)
  SCALAR, allocatable :: dummym(:,:)
  SCALAR, allocatable :: hmat(:,:)
  SCALAR, allocatable :: hbuf(:,:)
  SCALAR, allocatable :: smat(:,:)
  SCALAR, allocatable :: sbuf(:,:)
  SCALAR, allocatable :: zmat(:,:)
  SCALAR, allocatable :: zbuf(:,:)
  real(DP), allocatable :: w(:)
  integer, allocatable :: ifail(:)

  type (scalapack) :: scal

  PUSH_SUB(hdiag)

!------------------------
! Initialize

  minus_scale=sqrt(celvol)/dble(FFTgrid(1)*FFTgrid(2)*FFTgrid(3))
  plus_scale=1.0d0/sqrt(celvol)
  plus_scale_v=1.0d0
  if (fparafft) then
    call fft_set_p(FFTgrid,Nplane,Nrod)
  endif ! fparafft

  nvec=nbend-nbstart+1
  nmat=nvec*MAX(1,sapo_ham_ndim)
  is_diag=(sapo_ham_ndiag.gt.0)
  is_precond=(sapo_ham_ndiag.gt.1.and.sapo_ham_ndim.gt.1)

  abstol=0.0d0
  ! default sizes for lapack/scalapack arrays, as used in passing to subroutines, regardless of whether will be used
  lwork = 1
  lrwork = 1
  liwork = 1
  liclustr = 1
  lgap = 1
  if (fparala) then
    call blacs_setup(scal, nmat, .true.)
#ifdef USESCALAPACK
    orfac=1.0d-6
    call descinit(desca,nmat,nmat,scal%nbl,scal%nbl,0,0,scal%icntxt, &
      MAX(1,scal%npr),info)
#ifdef CPLX
    call PZHEGVX(1,'V','A','U',nmat,dummya,1,1,desca,dummyb,1,1,desca, &
      0.0d0,0.0d0,0,0,abstol,m,nz,dummyw,orfac,dummyz,1,1,desca, &
      dummywork,-1,dummyrwork,-1,dummyiwork,-1,dummyifail,dummyiclustr, &
      dummygap,info)
    lwork=nint(dble(dummywork(1)))
    lrwork=nint(dummyrwork(1))
#else
    call PDSYGVX(1,'V','A','U',nmat,dummya,1,1,desca,dummyb,1,1,desca, &
      0.0d0,0.0d0,0,0,abstol,m,nz,dummyw,orfac,dummyz,1,1,desca, &
      dummywork,-1,dummyiwork,-1,dummyifail,dummyiclustr,dummygap,info)
    lwork=nint(dummywork(1))
#endif
    liwork=dummyiwork(1)
    liclustr=2*scal%nprow*scal%npcol
    lgap=scal%nprow*scal%npcol
#endif
  else ! fparala
#ifdef CPLX
    call ZHEGVX(1,'V','A','U',nmat,dummya,nmat,dummyb,nmat,0.0d0,0.0d0, &
      0,0,abstol,m,dummyw,dummyz,nmat,dummywork,-1,dummyrwork,dummyiwork, &
      dummyifail,info)
    lwork=nint(dble(dummywork(1)))
    lrwork=7*nmat
#else
    call DSYGVX(1,'V','A','U',nmat,dummya,nmat,dummyb,nmat,0.0d0,0.0d0, &
      0,0,abstol,m,dummyw,dummyz,nmat,dummywork,-1,dummyiwork,dummyifail, &
      info)
    lwork=nint(dummywork(1))
#endif
    liwork=5*nmat
  endif ! fparala

!------------------------
! Memory report

  call procmem(amem,nmpinode)

#ifdef NOSIZEOF
  size_bool=4.0d0
  size_int=4.0d0
  size_real=8.0d0
  size_cplx=16.0d0
#else
  size_bool=dble(sizeof(is_diag))
  size_int=dble(sizeof(i))
  size_real=dble(sizeof(amem))
  size_cplx=dble(sizeof(ar))
#endif
  size_scal=dble(sizeof_scalar())

  rmem=0.0d0

  ! input arrays
  !
  rmem=rmem+3d0*dble(ng)*size_int
  rmem=rmem+dble(ngk_l)*dble(nk)*size_int
  rmem=rmem+dble(ng_l)*dble(ns)*size_scal
  if (fvkbg) then
    rmem=rmem+dble(nhm)*dble(nhm)*dble(nat)*dble(ns)*size_real
    rmem=rmem+dble(ngk_l)*dble(nkb)*dble(ns)*dble(nk)*size_cplx
  endif ! fvkbg
  rmem=rmem+dble(nbmax)*dble(nk)*dble(ns)*size_real
  rmem=rmem+dble(ngk_l)*dble(nbmax)*dble(ns)*dble(nk)*size_scal

  ! allocated arrays
  !
  rmem=rmem+dble(nvec)*size_real
  rmem=rmem+dble(nvec)*size_real
  rmem=rmem+dble(ngk_l)*dble(nmat)*size_scal
  rmem=rmem+dble(ngk_l)*dble(nvec)*size_scal
  rmem=rmem+dble(ngk_l)*dble(nvec)*size_scal
  rmem=rmem+dble(ngk_l)*size_cplx
  rmem=rmem+dble(ngk_g)*size_cplx
  rmem=rmem+dble(ngk_l)*size_cplx
  rmem=rmem+dble(ngk_l)*dble(nmat)*size_scal
  rmem=rmem+dble(ng_g)*size_int
  rmem=rmem+dble(ng)*size_real
  rmem=rmem+dble(nmat)*size_real
#ifdef MPI
  rmem=rmem+dble(nmat)*size_real
#endif
  !
  if (is_precond) then
    rmem=rmem+dble(ngk_l)*size_real
  endif ! is_precond
  !
  if (fvkbg) then
    rmem=rmem+dble(nkb)*size_cplx
#ifdef MPI
    rmem=rmem+dble(nkb)*size_cplx
#endif
    rmem=rmem+dble(nkb)*size_cplx
  endif ! fvkbg
  !
  rmem=rmem+dble(ns*nk)*size_bool
  !
  if (fparafft) then
    rmem=rmem+dble(peinf%npes)*size_int
    rmem=rmem+dble(peinf%npes)*size_int
    rmem=rmem+dble(peinf%npes)*size_int
    rmem=rmem+dble(peinf%npes)*size_int
    ! skip irod, isph, buffer_rod, and buffer_sph
    ! these arrays will be reallocated later
    rmem=rmem+product(dble(FFTgrid(1:2)))*dble(Nplane)*size_cplx
    rmem=rmem+dble(FFTgrid(3))*dble(Nrod)*size_cplx
    rmem=rmem+dble(Nrod)*dble(Nplane)*dble(peinf%npes)*size_cplx
    rmem=rmem+dble(Nrod)*dble(Nplane)*dble(peinf%npes)*size_cplx
    rmem=rmem+product(dble(FFTgrid(1:2)))*dble(Nplane)*ns*size_cplx
  else ! fparafft
    rmem=rmem+product(dble(FFTgrid(1:3)))*size_int
    rmem=rmem+product(dble(FFTgrid(1:3)))*size_cplx
    rmem=rmem+dble(ng_g)*size_cplx
    rmem=rmem+product(dble(FFTgrid(1:3)))*dble(ns)*size_cplx
  endif ! fparafft
  rmem=rmem+dble(ng_l)*size_cplx
  !
  if (fparala) then
    npr=MAX(1,scal%npr)
    npc=MAX(1,scal%npc)
    rmem=rmem+dble(nmat)*size_scal
    rmem=rmem+dble(nmat)*size_scal
#ifdef MPI
    rmem=rmem+dble(nmat)*size_scal
#endif
#ifdef USESCALAPACK
    rmem=rmem+dble(lwork)*size_scal
    rmem=rmem+dble(lrwork)*size_real
    rmem=rmem+dble(liwork)*size_int
    rmem=rmem+dble(liclustr)*size_int
    rmem=rmem+dble(lgap)*size_real
#endif
  else ! fparala
    npr=nmat
    npc=nmat
#ifdef MPI
    rmem=rmem+dble(npr)*dble(npc)*size_scal
#endif
    rmem=rmem+dble(lwork)*size_scal
    rmem=rmem+dble(lrwork)*size_real
    rmem=rmem+dble(liwork)*size_int
  endif ! fparala
  rmem=rmem+dble(npr)*dble(npc)*size_scal
  rmem=rmem+dble(npr)*dble(npc)*size_scal
  rmem=rmem+dble(npr)*dble(npc)*size_scal
  rmem=rmem+dble(npr)*dble(npc)*size_scal
  rmem=rmem+dble(npr)*dble(npc)*size_scal
  rmem=rmem+dble(npr)*dble(npc)*size_scal
  rmem=rmem+dble(nmat)*size_real
  rmem=rmem+dble(nmat)*size_int

  if (peinf%inode.eq.0) then
    write(s1,222) amem/1024.0d0**2
    write(s2,211) nmpinode
    write(s3,222) rmem/1024.0d0**2
    write(6,233) TRUNC(s1), TRUNC(s2), TRUNC(s3)
  endif ! peinf%inode.eq.0

!------------------------
! Allocate arrays

  SAFE_ALLOCATE(enbuf, (nvec))
  SAFE_ALLOCATE(enold, (nvec))
  SAFE_ALLOCATE(wfnbuf_d, (ngk_l,nmat))
  SAFE_ALLOCATE(wfnnew_d, (ngk_l,nvec))
  SAFE_ALLOCATE(residual_d, (ngk_l,nvec))
  SAFE_ALLOCATE(wfnc_d, (ngk_l))
  SAFE_ALLOCATE(wfnc, (ngk_g))
  SAFE_ALLOCATE(hpsic_d, (ngk_l))
  SAFE_ALLOCATE(hpsi_d, (ngk_l,nmat))
  SAFE_ALLOCATE(isort, (ng_g))
  SAFE_ALLOCATE(ekin, (ng))
  SAFE_ALLOCATE(norm, (nmat))
#ifdef MPI
  SAFE_ALLOCATE(ndum, (nmat))
#else
  SAFE_ALLOCATE(ndum, (1))
#endif
  !
  if (is_precond) then
    SAFE_ALLOCATE(h_gg, (ngk_l))
  endif ! is_precond
  !
  if (fvkbg) then
    SAFE_ALLOCATE(betapsi, (nkb))
#ifdef MPI
    SAFE_ALLOCATE(dummykb, (nkb))
#else
    SAFE_ALLOCATE(dummykb, (1))
#endif
    SAFE_ALLOCATE(ps, (nkb))
  endif ! fvkbg
  !
  SAFE_ALLOCATE(is_converged, (ns,nk))
  !
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
    SAFE_ALLOCATE(vlr_d, (FFTgrid(1),FFTgrid(2),Nplane,ns))
    SAFE_ALLOCATE(inv_indx, (1,1,1))
    SAFE_ALLOCATE(fftbox_3D, (1,1,1))
    SAFE_ALLOCATE(vlr, (1,1,1,1))
    ! vlg is not passed to subroutines in case of fparafft = .true., so it remains unallocated
  else ! fparafft
    SAFE_ALLOCATE(crod, (1))
    SAFE_ALLOCATE(csph, (1))
    SAFE_ALLOCATE(drod, (1))
    SAFE_ALLOCATE(dsph, (1))
    SAFE_ALLOCATE(irod, (1,1))
    SAFE_ALLOCATE(isph, (1))
    SAFE_ALLOCATE(buffer_rod, (1))
    SAFE_ALLOCATE(buffer_sph, (1))
    SAFE_ALLOCATE(fftbox_2D, (1,1,1))
    SAFE_ALLOCATE(fftbox_1D, (1,1))
    SAFE_ALLOCATE(buffer_2D, (1,1,1))
    SAFE_ALLOCATE(buffer_1D, (1,1,1))
    SAFE_ALLOCATE(vlr_d, (1,1,1,1))
    SAFE_ALLOCATE(inv_indx, (FFTgrid(1),FFTgrid(2),FFTgrid(3)))
    SAFE_ALLOCATE(fftbox_3D, (FFTgrid(1),FFTgrid(2),FFTgrid(3)))
    SAFE_ALLOCATE(vlr, (FFTgrid(1),FFTgrid(2),FFTgrid(3),ns))
    SAFE_ALLOCATE(vlg, (ng_g))
  endif ! fparafft
  SAFE_ALLOCATE(vlg_d, (ng_l))
  !
  if (fparala) then
    npr=MAX(1,scal%npr)
    npc=MAX(1,scal%npc)
    SAFE_ALLOCATE(hcol, (nmat))
    SAFE_ALLOCATE(scol, (nmat))
#ifdef MPI
    SAFE_ALLOCATE(dummyc, (nmat))
#else
    SAFE_ALLOCATE(dummyc, (1))
#endif
    SAFE_ALLOCATE(dummym, (1,1))
  else ! fparala
    npr=nmat
    npc=nmat
    SAFE_ALLOCATE(hcol, (1))
    SAFE_ALLOCATE(scol, (1))
    SAFE_ALLOCATE(dummyc, (1))
#ifdef MPI
    SAFE_ALLOCATE(dummym, (npr,npc))
#else
    SAFE_ALLOCATE(dummym, (1,1))
#endif
  endif ! fparala
  SAFE_ALLOCATE(work, (lwork))
  SAFE_ALLOCATE(rwork, (lrwork))
  SAFE_ALLOCATE(iwork, (liwork))
  SAFE_ALLOCATE(ifail, (nmat))
  SAFE_ALLOCATE(iclustr, (liclustr))
  SAFE_ALLOCATE(gap, (lgap))
  SAFE_ALLOCATE(hmat, (npr,npc))
  SAFE_ALLOCATE(hbuf, (npr,npc))
  SAFE_ALLOCATE(smat, (npr,npc))
  SAFE_ALLOCATE(sbuf, (npr,npc))
  SAFE_ALLOCATE(zmat, (npr,npc))
  SAFE_ALLOCATE(zbuf, (npr,npc))
  SAFE_ALLOCATE(w, (nmat))

!------------------------
! Invert indices of G-vectors for potential

  call timacc(14,1)
    
  isort(:)=0
  do ig=1,ng_g
    isort(ig)=ig
  enddo ! ig
  if (fparafft) then
    call fft_map_p(0,ng,ng_l,Nrod,FFTgrid,isort,gvec, &
      mrod,msph,crod,csph,drod,dsph,irod,isph)
    SAFE_DEALLOCATE_P(irod)
    SAFE_DEALLOCATE_P(isph)
    SAFE_DEALLOCATE_P(buffer_rod)
    SAFE_DEALLOCATE_P(buffer_sph)
    SAFE_ALLOCATE(irod, (2,MAX(1,mrod)))
    SAFE_ALLOCATE(isph, (MAX(1,msph)))
    SAFE_ALLOCATE(buffer_rod, (MAX(1,mrod)))
    SAFE_ALLOCATE(buffer_sph, (MAX(1,msph)))
    call fft_map_p(1,ng,ng_l,Nrod,FFTgrid,isort,gvec, &
      mrod,msph,crod,csph,drod,dsph,irod,isph)
  else ! fparafft
    call fft_map_s(ng,FFTgrid,isort,gvec,inv_indx)
  endif ! fparafft

!------------------------
! Bring potential from G-space to R-space

  call logit('Bringing pot to real space')
  do is=1,ns

    do ig=1,ng_l
      vlg_d(ig)=COMPLEXIFY(vsc_d(ig,is))
    enddo ! ig

    if (fparafft) then

      call fft_g2r_p(FFTgrid,Nplane,Nrod,mrod,msph,crod, &
        csph,drod,dsph,irod,isph,plus_scale_v,vlg_d,fftbox_2D, &
        fftbox_1D,buffer_2D,buffer_1D,buffer_rod,buffer_sph)
      vlr_d(:,:,:,is)=fftbox_2D(:,:,:)

    else ! fparafft

#ifdef MPI
      call MPI_Gather(vlg_d,ng_l,MPI_COMPLEX_DPC, &
        vlg,ng_l,MPI_COMPLEX_DPC,0,MPI_COMM_WORLD,mpierr)
#else
      vlg(:)=vlg_d(:)
#endif
      call fft_g2r_s(FFTgrid,plus_scale_v,inv_indx,vlg,fftbox_3D)
      vlr(:,:,:,is)=fftbox_3D(:,:,:)

    endif ! fparafft

  enddo ! is
  call logit('Ok')

  call timacc(14,2)

!------------------------
! Compute average potential

  if (is_precond) then

    call timacc(20,1)
    !
    v_0=0.0d0
    !
    if (fparafft) then
      !
      do is=1,ns
        do k=1,Nplane
          do j=1,FFTgrid(2)
            do i=1,FFTgrid(1)
              v_0=v_0+vlr_d(i,j,k,is)
            enddo ! i
          enddo ! j
        enddo ! k
      enddo ! is
      !
#ifdef MPI
      vdummy=v_0
      call MPI_Allreduce(vdummy,v_0,1, &
        MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
      !
    else ! fparafft
      !
      do is=1,ns
        do k=1,FFTgrid(3)
          do j=1,FFTgrid(2)
            do i=1,FFTgrid(1)
              v_0=v_0+vlr(i,j,k,is)
            enddo ! i
          enddo ! j
        enddo ! k
      enddo ! is
      !
    endif ! fparafft
    !
    v_0=v_0/dble(FFTgrid(1)*FFTgrid(2)*FFTgrid(3)*ns)
    !
    call timacc(20,2)

  endif ! is_precond

!------------------------
! Iterative diagonalization loop

  is_converged(:,:)=.false.

  call logit('Got to restart loop')
  restart_loop: do irestart = 1, MAX(1, sapo_ham_nrestart)

    idiagtot=0
    detot=0.0d0
    restot=0.0d0

    do ik=1,nk

      if (all(is_converged(:,ik))) cycle

!------------------------
! Compute kinetic energies of G-vectors

      call timacc(15,1)
      
      do ig=1,ng
        do i=1,3
          kvec(i)=kpt(i,ik)+dble(gvec(i,ig))
        enddo
        ekin(ig)=DOT_PRODUCT(kvec,MATMUL(bdot,kvec))
      enddo

      call timacc(15,2)

!------------------------
! Gather indices of G-vectors for current k-point

      call timacc(12,1)

#ifdef MPI
      call MPI_Gather(isort_d(:,ik),ngk_l,MPI_INTEGER, &
        isort,ngk_l,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(isort,ngk_g,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#else
      ! need explicit range here to avoid out of bounds error
      ! because in serial version isort is of size (ng) and
      ! isort_d is of size (ngkmax,nk)
      isort(1:ngk_l)=isort_d(1:ngk_l,ik)
#endif

      call timacc(12,2)

!------------------------
! Invert indices of G-vectors for current k-point

      call timacc(14,1)
      
      call logit('Inverting indices of G-vectors')
      if (fparafft) then
        call fft_map_p(0,ngk(ik),ngk_l,Nrod,FFTgrid,isort,gvec, &
          mrod,msph,crod,csph,drod,dsph,irod,isph)
        SAFE_DEALLOCATE_P(irod)
        SAFE_DEALLOCATE_P(isph)
        SAFE_DEALLOCATE_P(buffer_rod)
        SAFE_DEALLOCATE_P(buffer_sph)
        SAFE_ALLOCATE(irod, (2,MAX(1,mrod)))
        SAFE_ALLOCATE(isph, (MAX(1,msph)))
        SAFE_ALLOCATE(buffer_rod, (MAX(1,mrod)))
        SAFE_ALLOCATE(buffer_sph, (MAX(1,msph)))
        call fft_map_p(1,ngk(ik),ngk_l,Nrod,FFTgrid,isort,gvec, &
          mrod,msph,crod,csph,drod,dsph,irod,isph)
      else ! fparafft
        call fft_map_s(ngk(ik),FFTgrid,isort,gvec,inv_indx)
      endif ! fparafft
      call logit('Ok')
      
      call timacc(14,2)
      
      do is=1,ns

        if (is_converged(is,ik)) cycle

!------------------------
! Diagonal part of Hamiltonian

        if (is_precond) then

          call timacc(20,1)
          !
          do ig=1,ngk_l
            if (isort_d(ig,ik).ne.0) then
              h_gg(ig)=ekin(isort_d(ig,ik))+v_0
            else
              h_gg(ig)=0.0d0
            endif
          enddo
          !
          if (fvkbg) then
            ijkb0=0
            do isp=1,nsp
              do iat=1,nat
                if (ityp(iat).eq.isp) then
                  do ih=1,nh(isp)
                    ikb=ijkb0+ih
                    ps1=deeq(ih,ih,iat,is)
                    do ig=1,ngk_l
                      if (isort_d(ig,ik).ne.0) then
                        ar=vkb_d(ig,ikb,is,ik)*CONJG(vkb_d(ig,ikb,is,ik))
                        h_gg(ig)=h_gg(ig)+ps1*ar
                      endif
                    enddo
                  enddo
                  ijkb0=ijkb0+nh(isp)
                endif
              enddo
            enddo
          endif
          !
          call timacc(20,2)

        endif ! is_precond

        call timacc(20,1)

!------------------------
! Starting eigenvalues and wavefunctions

        do ib=nbstart,nbend
          im=ib-nbstart+1
          enbuf(im)=en(ib,ik,is)
        enddo ! ib

        do ib=nbstart,nbend
          im=ib-nbstart+1
          do ig=1,ngk_l
            wfnbuf_d(ig,im)=wfn_d(ig,ib,is,ik)
          enddo ! ig
        enddo ! ib

!------------------------
! Initialize subspace Hamiltonian and eigenvector matrices,
! only needed for computing residual vectors (sapo_ham_ndiag = 0)
! or if subspace diagonalization fails on the first iteration

        call logit('Calling main dsub')
        call dsub(fparala,nvec,nmat,npr,npc,scal, &
          enbuf,hcol,scol,hmat,smat,zmat)
        call logit('Ok')

!------------------------
! Initialize counters for the basis functions

        npre=0
        nbas=nvec
        resmax=INF

        call timacc(20,2)

!------------------------
! Iterative diagonalization loop

        diag_loop: do idiag = 1, MAX(1, sapo_ham_ndiag)

          call timacc(20,1)

!------------------------
! Normalize wavefunctions

          if (is_diag) then

            call npsi(.true.,fblas,npre,nbas,nmat,ngk_l,wfnbuf_d,norm,ndum)

          endif ! is_diag

!------------------------
! Store old eigenvectors

          zbuf(:,:)=zmat(:,:)

          call timacc(20,2)

!------------------------
! Apply Hamiltonian to wavefunctions

          impart_hpsi = 0.0d0
          call logit('Calculating H*Psi')
          do im=npre+1,nbas

            call hpsi(fparafft,fvkbg,nk,ns,nmat,ng,ngk_l,ngk_g, &
              ik,is,im,Nplane,Nrod,nat,nsp,nkb,nhm,FFTgrid,mrod,msph,crod, &
              csph,drod,dsph,irod,isph,ityp,nh,ekin,inv_indx,isort_d, &
              wfnbuf_d,wfnc_d,wfnc,vlr,vlr_d,vkb_d,deeq,betapsi,dummykb, &
              ps,fftbox_3D,fftbox_2D,fftbox_1D,buffer_2D,buffer_1D, &
              buffer_rod,buffer_sph,minus_scale,plus_scale,hpsic_d, &
              hpsi_d,impart_hpsi)

          enddo ! im
          call logit('Ok')
#ifdef MPI
          impart_dummy=impart_hpsi
          call MPI_Reduce(impart_dummy,impart_hpsi,1, &
            MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#endif
          if (peinf%inode.eq.0.and.impart_hpsi.gt.TOL_Small) then
            write(0,444) 'hpsi', impart_hpsi, idiag, is, ik, kpt(:,ik)
            write(0,455)
          endif

!------------------------
! Construct and diagonalize subspace Hamiltonian

          if (is_diag) then

            impart_hsub=0.0d0
            impart_ssub=0.0d0
            call logit('Diagonalizing subspace Hamiltonian')
            call hsub(fparala,fblas,npre,nvec,nbas,nmat,npr,npc, &
              ngk_l,lwork,lrwork,liwork,liclustr,lgap,scal,abstol, &
              orfac,desca,wfnbuf_d,hpsi_d,work,rwork,iwork,ifail, &
              iclustr,gap,hcol,scol,dummyc,dummym,w,hmat,hbuf,smat, &
              sbuf,zmat,impart_hsub,impart_ssub,is_error)
            call logit('Done diagonalizing subspace Hamiltonian')
#ifdef MPI
            impart_dummy=impart_hsub
            call MPI_Reduce(impart_dummy,impart_hsub,1, &
              MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
            impart_dummy=impart_ssub
            call MPI_Reduce(impart_dummy,impart_ssub,1, &
              MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#endif
            if (peinf%inode.eq.0.and.impart_hsub.gt.TOL_Small) &
              write(0,444) 'hsub', impart_hsub, idiag, is, ik, kpt(:,ik)
            if (peinf%inode.eq.0.and.impart_ssub.gt.TOL_Small) &
              write(0,444) 'ssub', impart_ssub, idiag, is, ik, kpt(:,ik)

          else ! is_diag

            is_error=.false.

          endif ! is_diag

          resold=resmax

          if (is_error) then

            call timacc(20,1)

!------------------------
! If subspace diagonalization failed restore the old
! eigenvectors and reset the variables for residual
! vectors and eigenvalues

            nbas=MAX(nvec,npre)
            zmat(:,:)=zbuf(:,:)
            resmax=INF
            nres=0
            demax=INF

            call timacc(20,2)

          else ! is_error

            call timacc(20,1)

!------------------------
! Rotate wavefunctions, wfnnew_d = zmat . wfnbuf_d

            call logit('Rotating WFNs')
            call rpsi(fparala,fblas,nvec,nbas,nmat,npr,npc,ngk_l, &
              scal,zmat,wfnbuf_d,wfnnew_d,scol,dummyc)
            call logit('Ok')

!------------------------
! Rotate hpsi vectors, residual_d = zmat . hpsi_d

            call logit('Rotating residues')
            call rpsi(fparala,fblas,nvec,nbas,nmat,npr,npc,ngk_l, &
              scal,zmat,hpsi_d,residual_d,scol,dummyc)
            call logit('Ok')

!------------------------
! Compute residual vectors, R = (H - E) psi

            call logit('Computing residues')
            do im=1,nvec
              do ig=1,ngk_l
                residual_d(ig,im)=residual_d(ig,im)-enbuf(im)* &
                  wfnnew_d(ig,im)
              enddo ! ig
            enddo ! im
            call logit('Ok')

!------------------------
! Compute norms of residual vectors

            call logit('Computing norm of residues')
            call npsi(.false.,fblas,0,nvec,nvec,ngk_l,residual_d,norm,ndum)
            call logit('Ok')

!------------------------
! Compute maximum norm of residual vectors

            resmax=MAXVAL(norm(1:nvec))

            call timacc(20,2)

!------------------------
! If computing the residual vectors only (is_diag = .false.)
! exit iterative diagonalization loop

            if (.not.is_diag) exit diag_loop

            call timacc(20,1)

!------------------------
! Check if convergence is reached

            is_converged(is,ik) = (resmax .le. sapo_ham_tol)

!------------------------
! Count unconverged wavefunctions

            jm=0
            do im=1,nvec
              if (norm(im).gt.sapo_ham_tol) jm=jm+1
            enddo ! im
            nres=jm

!------------------------
! Store old eigenvalues

            do im=1,nvec
              enold(im)=enbuf(im)
            enddo ! im

!------------------------
! Update eigenvalues

            do im=1,nvec
              enbuf(im)=w(im)
            enddo ! im

!------------------------
! Find the maximum change in eigenvalues

            demax=0.0d0
            do im=1,nvec
              if (demax.lt.abs(enold(im)-enbuf(im))) &
                demax=abs(enold(im)-enbuf(im))
            enddo ! im

            call timacc(20,2)

          endif ! is_error

          call timacc(20,1)

!------------------------
! Check whether to update the wavefunctions or to expand
! the basis set. If the size of the new basis set exceeds
! the maximum (nbas + nres > nmat) or we are on the last
! iteration (idiag = sapo_ham_ndiag) or the convergence is
! reached (is_converged(is,ik)) or subspace diagonalization
! failed (is_error) or the maximum norm of residual vectors
! is increasing (only checked if sapo_ham_resinc = .true.)
! then update the wavefunctions, otherwise expand the basis set

          is_updatewfn = ((nbas + nres .gt. nmat) .or. &
            (idiag .eq. MAX(1, sapo_ham_ndiag)) .or. &
            is_converged(is,ik) .or. &
            is_error .or. &
            (sapo_ham_resinc .and. resmax .gt. resold + TOL_Zero))

          if (is_updatewfn) then

!------------------------
! Rotate wavefunctions and hpsi vectors,
! update subspace Hamiltonian and eigenvector matrices,
! reset number of the basis functions

            call logit('Rotating WFNs (1)')
            call rpsi(fparala,fblas,nvec,nbas,nmat,npr,npc,ngk_l, &
              scal,zmat,wfnbuf_d,wfnnew_d,scol,dummyc)
            call logit('Ok')

            do im=1,nvec
              do ig=1,ngk_l
                wfnbuf_d(ig,im)=wfnnew_d(ig,im)
              enddo ! ig
            enddo ! im

            call logit('Rotating WFNs (2)')
            call rpsi(fparala,fblas,nvec,nbas,nmat,npr,npc,ngk_l, &
              scal,zmat,hpsi_d,wfnnew_d,scol,dummyc)
            call logit('Ok')

            do im=1,nvec
              do ig=1,ngk_l
                hpsi_d(ig,im)=wfnnew_d(ig,im)
              enddo ! ig
            enddo ! im

            call logit('Calling dsub')
            call dsub(fparala,nvec,nmat,npr,npc,scal, &
              enbuf,hcol,scol,hmat,smat,zmat)
            call logit('Ok')

            npre=nvec
            nbas=nvec

          else ! is_updatewfn

!------------------------
! Precondition residual vectors of unconverged wavefunctions,
! add them to the basis functions, increment number of the
! basis functions

            call logit('Preconditioning residues')
            jm=nbas
            do im=1,nvec
              if (norm(im).gt.sapo_ham_tol) then
                jm=jm+1
                do ig=1,ngk_l
                  if (isort_d(ig,ik).ne.0) then
                    call calc_precond(h_gg(ig),enbuf(im),precond)
                    wfnbuf_d(ig,jm)=precond*residual_d(ig,im)
                  else ! isort_d(ig,ik).ne.0
                    wfnbuf_d(ig,jm)=ZERO
                  endif ! isort_d(ig,ik).ne.0
                enddo ! ig
              endif ! norm(im).gt.sapo_ham_tol
            enddo ! im
            call logit('Ok')

            npre=nbas
            nbas=nbas+nres

          endif ! is_updatewfn

!------------------------
! Store the number of diagonalization steps

          idiagmax = idiag

!------------------------
! Print convergence information

          if (peinf%inode.eq.0) write(6,333) &
            ik, is, idiagmax, demax, resmax

          call timacc(20,2)

!------------------------
! If the convergence for this k-point/spin is reached
! exit iterative diagonalization loop

          if (is_converged(is,ik)) exit diag_loop

        enddo diag_loop ! idiag

!------------------------
! Final eigenvalues and wavefunctions

        call timacc(20,1)

        do ib=nbstart,nbend
          im=ib-nbstart+1
          en(ib,ik,is)=enbuf(im)
        enddo ! ib

        do ib=nbstart,nbend
          im=ib-nbstart+1
          do ig=1,ngk_l
            wfn_d(ig,ib,is,ik)=wfnbuf_d(ig,im)
          enddo ! ig
        enddo ! ib

        call timacc(20,2)

!------------------------
! Store total maximum change in eigenvalues,
! total maximum norm of residual vectors,
! total number of diagonalization steps

        if (is_diag) then
          if (idiagtot .lt. idiagmax) idiagtot = idiagmax
          if (detot .lt. demax) detot = demax
        endif ! is_diag
        if (restot .lt. resmax) restot = resmax

      enddo ! is
      
    enddo ! ik

!------------------------
! Print total convergence information

    if (peinf%inode.eq.0) then
      write(6,*)
      if (is_diag) then
        write(6,344) idiagtot, sapo_ham_ndiag
        write(6,355) detot
      endif ! is_diag
      write(6,366) restot
    endif ! peinf%inode.eq.0

!------------------------
! Set random initial wavefunctions
! for iterative diagonalization restart

    call timacc(20,1)

    if (irestart.lt.sapo_ham_nrestart) then
      do ik=1,nk
        do is=1,ns
          if (.not.is_converged(is,ik)) then
            call random_wfng(nk,ns,nbstart,nbend,nbmax,ngk_l,ik,is,ngk, &
              0.0d0,1.0d0,wfn_d)
          endif ! .not.is_converged(is,ik)
        enddo ! is
      enddo ! ik
    endif ! irestart.lt.sapo_ham_nrestart

    call timacc(20,2)

!------------------------
! If the convergence for all k-points/spins is reached
! exit restart loop

    if (all(is_converged(:,:))) exit restart_loop

  enddo restart_loop ! irestart

!------------------------
! Print warning if convergence has not been achieved

  if (peinf%inode.eq.0.and..not.all(is_converged(:,:))) write(0,466)

!------------------------
! Release BLACS

  if (fparala) then
#ifdef USESCALAPACK
    call blacs_exit(1)
#endif
  endif ! fparala
  
!------------------------
! Deallocate arrays

  SAFE_DEALLOCATE(enbuf)
  SAFE_DEALLOCATE(enold)
  SAFE_DEALLOCATE(wfnbuf_d)
  SAFE_DEALLOCATE(wfnnew_d)
  SAFE_DEALLOCATE(residual_d)
  SAFE_DEALLOCATE(wfnc_d)
  SAFE_DEALLOCATE(wfnc)
  SAFE_DEALLOCATE(hpsic_d)
  SAFE_DEALLOCATE(hpsi_d)
  SAFE_DEALLOCATE(isort)
  SAFE_DEALLOCATE(ekin)
  SAFE_DEALLOCATE(norm)
  SAFE_DEALLOCATE(ndum)
  !
  if (is_precond) then
    SAFE_DEALLOCATE(h_gg)
  endif ! is_precond
  !
  if (fvkbg) then
    SAFE_DEALLOCATE(betapsi)
    SAFE_DEALLOCATE(dummykb)
    SAFE_DEALLOCATE(ps)
  endif ! fvkbg
  !
  SAFE_DEALLOCATE(is_converged)
  !
  SAFE_DEALLOCATE_P(crod)
  SAFE_DEALLOCATE_P(csph)
  SAFE_DEALLOCATE_P(drod)
  SAFE_DEALLOCATE_P(dsph)
  SAFE_DEALLOCATE_P(irod)
  SAFE_DEALLOCATE_P(isph)
  SAFE_DEALLOCATE_P(buffer_rod)
  SAFE_DEALLOCATE_P(buffer_sph)
  SAFE_DEALLOCATE_P(fftbox_2D)
  SAFE_DEALLOCATE_P(fftbox_1D)
  SAFE_DEALLOCATE_P(buffer_2D)
  SAFE_DEALLOCATE_P(buffer_1D)
  SAFE_DEALLOCATE_P(vlr_d)
  SAFE_DEALLOCATE_P(inv_indx)
  SAFE_DEALLOCATE_P(fftbox_3D)
  SAFE_DEALLOCATE_P(vlr)
  if (.not.fparafft) then
    SAFE_DEALLOCATE(vlg)
  endif ! .not.fparafft
  SAFE_DEALLOCATE(vlg_d)
  !
  SAFE_DEALLOCATE(hcol)
  SAFE_DEALLOCATE(scol)
  SAFE_DEALLOCATE(dummyc)
  SAFE_DEALLOCATE(dummym)
  SAFE_DEALLOCATE(work)
  SAFE_DEALLOCATE(rwork)
  SAFE_DEALLOCATE(iwork)
  SAFE_DEALLOCATE(ifail)
  SAFE_DEALLOCATE(iclustr)
  SAFE_DEALLOCATE(gap)
  SAFE_DEALLOCATE(hmat)
  SAFE_DEALLOCATE(hbuf)
  SAFE_DEALLOCATE(smat)
  SAFE_DEALLOCATE(sbuf)
  SAFE_DEALLOCATE(zmat)
  SAFE_DEALLOCATE(zbuf)
  SAFE_DEALLOCATE(w)

  call destroy_fftw_plans()

  POP_SUB(hdiag)

  return

211 format(i16)
222 format(f16.1)
233 format(1x,'Memory available:',1x,a,1x,'MB per MPI task,',1x,a,1x, &
  'MPI tasks per node',/,1x,'Memory required for execution:',1x,a,1x, &
  'MB per MPI task',/)

333 format(3x,"ik =",i5,3x,"is =",i3,3x,"idiag =",i5,3x, &
  "dE(Ry) =",g20.15,3x,"residual =",g20.15)
344 format(1x,"Maximum number of iterations:",/, &
  4x,i6,2x,"of",2x,i6,/)
355 format(1x,"Maximum change in eigenvalues (Ry):",/, &
  4x,"max | delta E_nk | =",g20.15,/)
366 format(1x,"Maximum norm of residual vectors (Ry):",/, &
  4x,"max | R_nk | =",g20.15,/)

444 format(1x,"WARNING: max(Im(",a4,")) =",e13.6," for idiag =",i3, &
  " is =",i2," ik =",i3,/,3x,"k = (",3f10.6,")",/)
455 format(1x,"!!! MAKE SURE INVERSION SYMMETRY HAS NO FRACTIONAL", &
  " TRANSLATION !!!",/)
466 format(1x,"WARNING: Convergence has not been achieved.",/)

end subroutine hdiag

!===============================================================================

subroutine hpsi(fparafft,fvkbg,nk,ns,nmat,ng,ngk_l,ngk_g, &
  ik,is,im,Nplane,Nrod,nat,nsp,nkb,nhm,FFTgrid,mrod,msph,crod, &
  csph,drod,dsph,irod,isph,ityp,nh,ekin,inv_indx,isort_d, &
  wfnbuf_d,wfnc_d,wfnc,vlr,vlr_d,vkb_d,deeq,betapsi,dummykb, &
  ps,fftbox_3D,fftbox_2D,fftbox_1D,buffer_2D,buffer_1D, &
  buffer_rod,buffer_sph,minus_scale,plus_scale,hpsic_d, &
  hpsi_d,impart_hpsi)

! This subroutine applies the Hamiltonian to a wavefunction.
! On entrance takes wfnbuf_d(1:ngk_l,im) set in hdiag. On exit
! returns hpsi_d(1:ngk_l,im), in real flavor sets impart_hpsi
! to max imaginary part of hpsi_d(1:ngk_l,im)

  logical, intent(in) :: fparafft,fvkbg
  integer, intent(in) :: nk,ns,nmat,ng,ngk_l,ngk_g, &
  ik,is,im,Nplane,Nrod,nat,nsp,nkb,nhm
  integer, intent(in) :: FFTgrid(3)
  integer, intent(inout) :: mrod,msph
  integer, intent(inout) :: crod(:) !< (peinf%npes)
  integer, intent(inout) :: csph(:) !< (peinf%npes)
  integer, intent(inout) :: drod(:) !< (peinf%npes)
  integer, intent(inout) :: dsph(:) !< (peinf%npes)
  integer, intent(inout) :: irod(:,:) !< (2,MAX(1,mrod))
  integer, intent(inout) :: isph(:) !< (MAX(1,msph))
  integer, intent(in) :: ityp(:) !< (nat)
  integer, intent(in) :: nh(:) !< (nsp)
  real(DP), intent(in) :: ekin(:) !< (ng)
  integer, intent(in) :: inv_indx(:,:,:) !< (FFTgrid(1),FFTgrid(2),FFTgrid(3))
  integer, intent(in) :: isort_d(:,:) !< (ngk_l,nk)
  SCALAR, intent(in) :: wfnbuf_d(:,:) !< (ngk_l,nmat)
  complex(DPC), intent(out) :: wfnc_d(:) !< (ngk_l)
  complex(DPC), intent(out) :: wfnc(:) !< (ngk_g)
  complex(DPC), intent(in) :: vlr(:,:,:,:) !< (FFTgrid(1),FFTgrid(2),FFTgrid(3),ns)
  complex(DPC), intent(in) :: vlr_d(:,:,:,:) !< (FFTgrid(1),FFTgrid(2),Nplane,ns)
  complex(DPC), intent(in) :: vkb_d(:,:,:,:) !< (ngk_l,nkb,ns,nk)
  real(DP), intent(in) :: deeq(:,:,:,:) !< (nhm,nhm,nat,ns)
  complex(DPC), intent(out) :: betapsi(:) !< (nkb)
  complex(DPC), intent(inout) :: dummykb(:) !< (nkb)
  complex(DPC), intent(out) :: ps(:) !< (nkb)
  complex(DPC), intent(out) :: fftbox_3D(:,:,:) !< (FFTgrid(1),FFTgrid(2),FFTgrid(3))
  complex(DPC), intent(out) :: fftbox_2D(:,:,:) !< (FFTgrid(1),FFTgrid(2),Nplane)
  complex(DPC), intent(out) :: fftbox_1D(:,:) !< (FFTgrid(3),Nrod)
  complex(DPC), intent(out) :: buffer_2D(:,:,:) !< (Nrod,Nplane,peinf%npes)
  complex(DPC), intent(out) :: buffer_1D(:,:,:) !< (Nrod,Nplane,peinf%npes)
  complex(DPC), intent(out) :: buffer_rod(:) !< (MAX(1,mrod))
  complex(DPC), intent(out) :: buffer_sph(:) !< (MAX(1,msph))
  real(DP), intent(in) :: minus_scale,plus_scale
  complex(DPC), intent(out) :: hpsic_d(:) !< (ngk_l)
  SCALAR, intent(inout) :: hpsi_d(:,:) !< (ngk_l,nmat)
  real(DP), intent(inout) :: impart_hpsi

  integer :: i,j,k,ig,iat,isp,ih,jh,ikb,jkb,ijkb0

  PUSH_SUB(hpsi)

!------------------------
! Convert wavefunction from scalar to complex

  call timacc(20,1)

  do ig=1,ngk_l
    wfnc_d(ig)=COMPLEXIFY(wfnbuf_d(ig,im))
  enddo ! ig

  call timacc(20,2)

!------------------------
! Bring wavefunction to R-space

  call timacc(14,1)
  
  if (fparafft) then
    call fft_g2r_p(FFTgrid,Nplane,Nrod,mrod,msph,crod, &
      csph,drod,dsph,irod,isph,plus_scale,wfnc_d,fftbox_2D, &
      fftbox_1D,buffer_2D,buffer_1D,buffer_rod,buffer_sph)
  else
#ifdef MPI
    call MPI_Gather(wfnc_d,ngk_l,MPI_COMPLEX_DPC, &
      wfnc,ngk_l,MPI_COMPLEX_DPC,0,MPI_COMM_WORLD,mpierr)
#else
    wfnc(:)=wfnc_d(:)
#endif
    call fft_g2r_s(FFTgrid,plus_scale,inv_indx,wfnc,fftbox_3D)
  endif

  call timacc(14,2)

!------------------------
! Multiply wavefunction by local potential

  call timacc(16,1)
  
  if (fparafft) then
    do k=1,Nplane
      do j=1,FFTgrid(2)
        do i=1,FFTgrid(1)
          fftbox_2D(i,j,k)=vlr_d(i,j,k,is)*fftbox_2D(i,j,k)
        enddo ! i
      enddo ! j
    enddo ! k
  else ! fparafft
    do k=1,FFTgrid(3)
      do j=1,FFTgrid(2)
        do i=1,FFTgrid(1)
          fftbox_3D(i,j,k)=vlr(i,j,k,is)*fftbox_3D(i,j,k)
        enddo ! i
      enddo ! j
    enddo ! k
  endif ! fparafft
  
  call timacc(16,2)

!------------------------
! Bring the product back to G-space
  
  call timacc(14,1)
  
  if (fparafft) then
    call fft_r2g_p(FFTgrid,Nplane,Nrod,mrod,msph,crod, &
      csph,drod,dsph,irod,isph,minus_scale,fftbox_2D,fftbox_1D, &
      buffer_2D,buffer_1D,buffer_rod,buffer_sph,hpsic_d)
  else
    call fft_r2g_s(FFTgrid,minus_scale,inv_indx,fftbox_3D,wfnc)
#ifdef MPI
    call MPI_Scatter(wfnc,ngk_l,MPI_COMPLEX_DPC, &
      hpsic_d,ngk_l,MPI_COMPLEX_DPC,0,MPI_COMM_WORLD,mpierr)
#else
    hpsic_d(:)=wfnc(:)
#endif
  endif
  
  call timacc(14,2)

!------------------------
! Add non-local pseudopotential

  call timacc(17,1)
  
  if (fvkbg) then
    
    betapsi(:)=(0.0d0,0.0d0)
    do ikb=1,nkb
      betapsi(ikb) = betapsi(ikb) + sum(CONJG(vkb_d(:,ikb,is,ik))*wfnc_d(:))
      ! FHJ: The following line is equivalent to the previous, but should be
      ! faster. I`m keeping it commented out because we don`t have a test
      ! case for complex SAPO.
      !betapsi(ikb) = betapsi(ikb) + zdotc(ngk_l, vkb_d(:,ikb,is,ik), 1, wfnc_d(:), 1)
    enddo ! ikb
    
#ifdef MPI
    dummykb(:)=betapsi(:)
    call MPI_Allreduce(dummykb,betapsi,nkb, &
      MPI_COMPLEX_DPC,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
    
    ps(:)=(0.0d0,0.0d0)
    ijkb0=0
    do isp=1,nsp
      do iat=1,nat
        if (ityp(iat).eq.isp) then
          do ih=1,nh(isp)
            ikb=ijkb0+ih
            do jh=1,nh(isp)
              jkb=ijkb0+jh
              ps(jkb)=ps(jkb)+ &
                deeq(jh,ih,iat,is)*betapsi(ikb)
            enddo ! jh
          enddo ! ih
          ijkb0=ijkb0+nh(isp)
        endif ! ityp(iat).eq.isp
      enddo ! iat
    enddo ! isp
    
    do ikb=1,nkb
      do ig=1,ngk_l
        ! FHJ: One should technically enclose the following line with a
        ! "if (isort_d(ig,ik)/=0)". I`m not doing that because vkb_d is zero
        ! when isort_d(ig)==0, and conditionals penalize vectorization.
        hpsic_d(ig)=hpsic_d(ig)+ &
          vkb_d(ig,ikb,is,ik)*ps(ikb)
      enddo ! ig
    enddo ! ikb
    
  endif ! fvkbg
  
  call timacc(17,2)
  
!------------------------
! Add kinetic energy

  call timacc(15,1)
  
  do ig=1,ngk_l
    if (isort_d(ig,ik).ne.0) &
      hpsic_d(ig)=hpsic_d(ig)+ &
        ekin(isort_d(ig,ik))*wfnc_d(ig)
  enddo ! ig
  
  call timacc(15,2)

!------------------------
! Convert from complex to scalar

  call timacc(12,1)

  do ig=1,ngk_l
    if (isort_d(ig,ik).ne.0) then
      hpsi_d(ig,im)=SCALARIFY(hpsic_d(ig))
#ifndef CPLX
      impart_hpsi = MAX(impart_hpsi, abs(IMAG(hpsic_d(ig))))
#endif
    else
      hpsi_d(ig,im)=ZERO
    endif
  enddo ! ig

  call timacc(12,2)

  POP_SUB(hpsi)

  return

end subroutine hpsi

!===============================================================================

subroutine hsub(fparala,fblas,npre,nvec,nbas,nmat,npr,npc, &
  ngk_l,lwork,lrwork,liwork,liclustr,lgap,scal,abstol, &
  orfac,desca,wfnbuf_d,hpsi_d,work,rwork,iwork,ifail, &
  iclustr,gap,hcol,scol,dummyc,dummym,w,hmat,hbuf,smat, &
  sbuf,zmat,impart_hsub,impart_ssub,is_error)

! This subroutine constructs and diagonalizes the subspace
! Hamiltonian. On entrance takes wfnbuf_d(1:ngk_l,1:nbas)
! set in hdiag and hpsi_d(1:ngk_l,1:nbas) computed in hpsi.
! On exit returns eigenvalues in w(1:nvec) and eigenvectors in
! zmat(1:npr,1:npc)/(1:nbas,1:nvec) (if fparala = .true./.false.),
! in complex flavor sets impart_hsub/impart_ssub to max imaginary
! part of hmat(im,im)/smat(im,im)

  logical, intent(in) :: fparala,fblas
  integer, intent(in) :: npre,nvec,nbas,nmat,npr,npc,ngk_l, &
    lwork,lrwork,liwork,liclustr,lgap
  type (scalapack), intent(in) :: scal
  real(DP), intent(in) :: abstol,orfac
  integer, intent(in) :: desca(9)
  SCALAR, intent(in) :: wfnbuf_d(:,:) !< (ngk_l,nmat)
  SCALAR, intent(in) :: hpsi_d(:,:) !< (ngk_l,nmat)
  SCALAR, intent(inout) :: work(:) !< (lwork)
  real(DP), intent(inout) :: rwork(:) !< (lrwork)
  real(DP), intent(inout) :: gap(:) !< (lgap)
  integer, intent(inout) :: iwork(:) !< (liwork)
  integer, intent(inout) :: ifail(:) !< (nmat)
  integer, intent(inout) :: iclustr(:) !< (liclustr)
  SCALAR, intent(inout) :: hcol(:) !< (nmat)
  SCALAR, intent(inout) :: scol(:) !< (nmat)
  SCALAR, intent(inout) :: dummyc(:) !< (nmat)
  SCALAR, intent(inout) :: dummym(:,:) !< (nmat,nmat)
  real(DP), intent(out) :: w(:) !< (nmat)
  SCALAR, intent(inout) :: hmat(:,:) !< (npr,npc)
  SCALAR, intent(inout) :: hbuf(:,:) !< (npr,npc)
  SCALAR, intent(inout) :: smat(:,:) !< (npr,npc)
  SCALAR, intent(inout) :: sbuf(:,:) !< (npr,npc)
  SCALAR, intent(inout) :: zmat(:,:) !< (npr,npc)
  real(DP), intent(inout) :: impart_hsub
  real(DP), intent(inout) :: impart_ssub
  logical, intent(out) :: is_error

  integer :: im,jm,ig,icurr,icol,irow,icolm,irowm,m,nz,info
  character :: RANGE

  PUSH_SUB(hsub)

!------------------------
! Construct subspace Hamiltonian

  call timacc(18,1)

  hbuf(:,:)=ZERO
  sbuf(:,:)=ZERO

  call logit('Constructing Hamiltonian')
  if (fparala) then
    icurr=0
    do im=1,nbas
      hcol(:)=ZERO
      scol(:)=ZERO
      if (im.gt.npre) then
        if (fblas) then
          call X(gemv)('C',ngk_l,im,ONE,wfnbuf_d(:,:),ngk_l, &
            hpsi_d(:,im),1,ZERO,hcol(:),1)
          call X(gemv)('C',ngk_l,im,ONE,wfnbuf_d(:,:),ngk_l, &
            wfnbuf_d(:,im),1,ZERO,scol(:),1)
        else ! fblas
          do jm=1,im
            do ig=1,ngk_l
              hcol(jm)=hcol(jm)+MYCONJG(wfnbuf_d(ig,jm))* &
                hpsi_d(ig,im)
              scol(jm)=scol(jm)+MYCONJG(wfnbuf_d(ig,jm))* &
                wfnbuf_d(ig,im)
            enddo ! ig
          enddo ! jm
        endif ! fblas
#ifdef MPI
        dummyc(:)=hcol(:)
        call MPI_Allreduce(dummyc,hcol,im, &
          MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
        dummyc(:)=scol(:)
        call MPI_Allreduce(dummyc,scol,im, &
          MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
#ifdef CPLX
        impart_hsub = MAX(impart_hsub, abs(IMAG(hcol(im))))
        impart_ssub = MAX(impart_ssub, abs(IMAG(scol(im))))
        hcol(im)=CMPLX(dble(hcol(im)),0.0d0)
        scol(im)=CMPLX(dble(scol(im)),0.0d0)
#endif
      endif ! im.gt.npre
      icol=mod(int(((im-1)/scal%nbl)+TOL_Small),scal%npcol)
      ! jm up to nmat instead of nbas to correctly increment icurr
      do jm=1,nmat
        irow=mod(int(((jm-1)/scal%nbl)+TOL_Small),scal%nprow)
        if (irow.eq.scal%myprow.and.icol.eq.scal%mypcol) then
          icurr=icurr+1
          icolm=int((icurr-1)/scal%npr+TOL_Small)+1
          irowm=mod((icurr-1),scal%npr)+1
          if (im.le.npre) then
            hbuf(irowm,icolm)=hmat(irowm,icolm)
            sbuf(irowm,icolm)=smat(irowm,icolm)
          else ! im.le.npre
            hbuf(irowm,icolm)=hcol(jm)
            sbuf(irowm,icolm)=scol(jm)
          endif ! im.le.npre
        endif
      enddo ! jm
    enddo ! im
  else ! fparala
    do im=1,npre
      do jm=1,im
        hbuf(jm,im)=hmat(jm,im)
        sbuf(jm,im)=smat(jm,im)
      enddo ! jm
    enddo ! im
    if (fblas) then
      call X(gemm)('C','N',nbas,nbas-npre,ngk_l,ONE,wfnbuf_d(:,:), &
        ngk_l,hpsi_d(:,npre+1),ngk_l,ZERO,hbuf(:,npre+1),nmat)
      call X(gemm)('C','N',nbas,nbas-npre,ngk_l,ONE,wfnbuf_d(:,:), &
        ngk_l,wfnbuf_d(:,npre+1),ngk_l,ZERO,sbuf(:,npre+1),nmat)
    else ! fblas
      do im=npre+1,nbas
        do jm=1,im
          do ig=1,ngk_l
            hbuf(jm,im)=hbuf(jm,im)+ &
              MYCONJG(wfnbuf_d(ig,jm))*hpsi_d(ig,im)
            sbuf(jm,im)=sbuf(jm,im)+ &
              MYCONJG(wfnbuf_d(ig,jm))*wfnbuf_d(ig,im)
          enddo ! ig
        enddo ! jm
      enddo ! im
    endif ! fblas
#ifdef CPLX
    do im=npre+1,nbas
      impart_hsub = MAX(impart_hsub, abs(IMAG(hbuf(im,im))))
      impart_ssub = MAX(impart_ssub, abs(IMAG(sbuf(im,im))))
      hbuf(im,im)=CMPLX(dble(hbuf(im,im)),0.0d0)
      sbuf(im,im)=CMPLX(dble(sbuf(im,im)),0.0d0)
    enddo ! im
#endif
  endif ! fparala
  call logit('Ok')

  call timacc(18,2)

!------------------------
! Assemble subspace Hamiltonian

  call timacc(12,1)

  if (.not.fparala) then
#ifdef MPI
    dummym(:,:)=hbuf(:,:)
    call MPI_Allreduce(dummym(1,1),hbuf(1,1),nmat*nmat, &
      MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
    dummym(:,:)=sbuf(:,:)
    call MPI_Allreduce(dummym(1,1),sbuf(1,1),nmat*nmat, &
      MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
  endif ! .not.fparala

  call timacc(12,2)

!------------------------
! Store subspace Hamiltonian
! (because LAPACK/ScaLAPACK destroys input arrays)

  call timacc(18,1)

  hmat(:,:)=hbuf(:,:)
  smat(:,:)=sbuf(:,:)

  call timacc(18,2)

!------------------------
! Diagonalize subspace Hamiltonian

  call timacc(19,1)

  call logit('Doing diagonalization')
  if (nvec.eq.nbas) then
    RANGE='A'
  else ! nvec.eq.nbas
    RANGE='I'
  endif ! nvec.eq.nbas
  if (fparala) then
#ifdef USESCALAPACK
#ifdef CPLX
    call PZHEGVX(1,'V',RANGE,'U',nbas,hbuf,1,1,desca,sbuf,1,1,desca, &
      0.0d0,0.0d0,1,nvec,abstol,m,nz,w,orfac,zmat,1,1,desca,work, &
      lwork,rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)
#else
    call PDSYGVX(1,'V',RANGE,'U',nbas,hbuf,1,1,desca,sbuf,1,1,desca, &
      0.0d0,0.0d0,1,nvec,abstol,m,nz,w,orfac,zmat,1,1,desca,work, &
      lwork,iwork,liwork,ifail,iclustr,gap,info)
#endif
    is_error = (info.ne.0.or.m.ne.nvec.or.nz.ne.nvec)
#else
    is_error = .true.
#endif
  else ! fparala
#ifdef CPLX
    call ZHEGVX(1,'V',RANGE,'U',nbas,hbuf,nmat,sbuf,nmat,0.0d0,0.0d0, &
      1,nvec,abstol,m,w,zmat,nmat,work,lwork,rwork,iwork,ifail,info)
#else
    call DSYGVX(1,'V',RANGE,'U',nbas,hbuf,nmat,sbuf,nmat,0.0d0,0.0d0, &
      1,nvec,abstol,m,w,zmat,nmat,work,lwork,iwork,ifail,info)
#endif
    is_error = (info.ne.0.or.m.ne.nvec)
  endif ! fparala
  call logit('Ok')

  call timacc(19,2)

!  if (is_error) then
!    call die('hsub: failed to diagonalize subspace Hamiltonian')
!  endif ! is_error

  if (peinf%inode.eq.0.and.is_error) write(0,555)

  POP_SUB(hsub)

  return

555 format(1x, &
  "WARNING: subspace diagonalization failed, resetting the basis functions")

end subroutine hsub

!===============================================================================

subroutine dsub(fparala,nvec,nmat,npr,npc,scal, &
  enbuf,hcol,scol,hmat,smat,zmat)

! This subroutine constructs the diagonal subspace Hamiltonian
! and diagonal eigenvectors. On entrance takes eigenvalues in
! enbuf(1:nvec). On exit returns hmat(1:npr,1:npc)/(1:nvec,1:nvec),
! smat(1:npr,1:npc)/(1:nvec,1:nvec), zmat(1:npr,1:npc)/(1:nvec,1:nvec)
! (if fparala = .true./.false.)

  logical, intent(in) :: fparala
  integer, intent(in) :: nvec,nmat,npr,npc
  type (scalapack), intent(in) :: scal
  real(DP), intent(in) :: enbuf(:) !< (nvec)
  SCALAR, intent(inout) :: hcol(:) !< (nmat)
  SCALAR, intent(inout) :: scol(:) !< (nmat)
  SCALAR, intent(out) :: hmat(:,:) !< (npr,npc)
  SCALAR, intent(out) :: smat(:,:) !< (npr,npc)
  SCALAR, intent(out) :: zmat(:,:) !< (npr,npc)

  integer :: im,jm,icurr,icol,irow,icolm,irowm

  PUSH_SUB(dsub)

!------------------------
! Construct subspace Hamiltonian

  call timacc(18,1)

  hmat(:,:)=ZERO
  smat(:,:)=ZERO
  zmat(:,:)=ZERO

  if (fparala) then
    icurr=0
    do im=1,nvec
      hcol(:)=ZERO
      scol(:)=ZERO
      hcol(im)=SCALARIFY(enbuf(im))
      scol(im)=ONE
      icol=mod(int(((im-1)/scal%nbl)+TOL_Small),scal%npcol)
      ! jm up to nmat instead of nbas to correctly increment icurr
      do jm=1,nmat
        irow=mod(int(((jm-1)/scal%nbl)+TOL_Small),scal%nprow)
        if (irow.eq.scal%myprow.and.icol.eq.scal%mypcol) then
          icurr=icurr+1
          icolm=int((icurr-1)/scal%npr+TOL_Small)+1
          irowm=mod((icurr-1),scal%npr)+1
          hmat(irowm,icolm)=hcol(jm)
          smat(irowm,icolm)=scol(jm)
          zmat(irowm,icolm)=scol(jm)
        endif
      enddo ! jm
    enddo ! im
  else ! fparala
    do im=1,nvec
      hmat(im,im)=SCALARIFY(enbuf(im))
      smat(im,im)=ONE
      zmat(im,im)=ONE
    enddo ! im
  endif ! fparala

  call timacc(18,2)

  POP_SUB(dsub)

  return

end subroutine dsub

!===============================================================================

subroutine rpsi(fparala,fblas,nvec,nbas,nmat,npr,npc,ngk_l, &
  scal,zmat,wfnbuf_d,wfnnew_d,scol,dummyc)

! This subroutine rotates wavefunctions. On entrance takes
! wavefunctions in wfnbuf_d(1:ngk_l,1:nbas) set in hdiag
! and eigenvectors in zmat(1:npr,1:npc)/(1:nbas,1:nvec)
! (if fparala = .true./.false.) computed in hsub. On exit
! returns rotated wavefunctions in wfnnew_d(1:ngk_l,1:nvec)

  logical, intent(in) :: fparala,fblas
  integer, intent(in) :: nvec,nbas,nmat,npr,npc,ngk_l
  type (scalapack), intent(in) :: scal
  SCALAR, intent(in) :: zmat(:,:) !< (npr,npc)
  SCALAR, intent(in) :: wfnbuf_d(:,:) !< (ngk_l,nmat)
  SCALAR, intent(out) :: wfnnew_d(:,:) !< (ngk_l,nvec)
  SCALAR, intent(inout) :: scol(:) !< (nmat)
  SCALAR, intent(inout) :: dummyc(:) !< (nmat)

  integer :: im,jm,ig,icurr,icol,irow,icolm,irowm

  PUSH_SUB(rpsi)

  wfnnew_d(:,:)=ZERO

  if (fparala) then
    icurr=0
    do im=1,nvec
      scol(:)=ZERO
      icol=mod(int(((im-1)/scal%nbl)+TOL_Small),scal%npcol)
      ! jm up to nmat instead of nbas to correctly increment icurr
      do jm=1,nmat
        irow=mod(int(((jm-1)/scal%nbl)+TOL_Small),scal%nprow)
        if (irow.eq.scal%myprow.and.icol.eq.scal%mypcol) then
          icurr=icurr+1
          icolm=int((icurr-1)/scal%npr+TOL_Small)+1
          irowm=mod((icurr-1),scal%npr)+1
          scol(jm)=zmat(irowm,icolm)
        endif
      enddo ! jm
#ifdef MPI
      dummyc(:)=scol(:)
      call MPI_Allreduce(dummyc,scol,nbas, &
        MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
      if (fblas) then
        call X(gemv)('N',ngk_l,nbas,ONE,wfnbuf_d(:,:),ngk_l, &
          scol(:),1,ZERO,wfnnew_d(:,im),1)
      else ! fblas
        do jm=1,nbas
          do ig=1,ngk_l
            wfnnew_d(ig,im)=wfnnew_d(ig,im)+ &
              wfnbuf_d(ig,jm)*scol(jm)
          enddo ! ig
        enddo ! jm
      endif ! fblas
    enddo ! im
  else ! fparala
    if (fblas) then
      call X(gemm)('N','N',ngk_l,nvec,nbas,ONE,wfnbuf_d(:,:), &
        ngk_l,zmat(:,:),nmat,ZERO,wfnnew_d(:,:),ngk_l)
    else ! fblas
      do im=1,nvec
        do jm=1,nbas
          do ig=1,ngk_l
            wfnnew_d(ig,im)=wfnnew_d(ig,im)+ &
              wfnbuf_d(ig,jm)*zmat(jm,im)
          enddo ! ig
        enddo ! jm
      enddo ! im
    endif ! fblas
  endif ! fparala

  POP_SUB(rpsi)

  return

end subroutine rpsi

!===============================================================================

subroutine npsi(fnorm,fblas,npre,nvec,nmat,ngk_l,wfnbuf_d,norm,ndum)

! This subroutine computes norms of wavefunctions and optionally
! normalizes wavefunctions. On entrance takes wavefunctions in
! wfnbuf_d(1:ngk_l,1:nvec) set in hdiag. On exit returns norms
! in norm(1:nvec), if fnorm = .true. normalizes wavefunctions
! in wfnbuf_d(1:ngk_l,1:nvec)

  logical, intent(in) :: fnorm,fblas
  integer, intent(in) :: npre,nvec,nmat,ngk_l
  SCALAR, intent(inout) :: wfnbuf_d(:,:) !< (ngk_l,nmat)
  real(DP), intent(out) :: norm(:) !< (nmat)
  real(DP), intent(inout) :: ndum(:) !< (nmat)

  integer :: im,ig
  SCALAR :: inorm

  PUSH_SUB(npsi)

  norm(:)=0.0d0

  if (fblas) then
    do im=npre+1,nvec
#ifdef CPLX
      norm(im)=dble(ZDOTC(ngk_l,wfnbuf_d(:,im),1,wfnbuf_d(:,im),1))
#else
      norm(im)=DDOT(ngk_l,wfnbuf_d(:,im),1,wfnbuf_d(:,im),1)
#endif
    enddo ! im
  else ! fblas
    do im=npre+1,nvec
      do ig=1,ngk_l
        norm(im)=norm(im)+MYCONJG(wfnbuf_d(ig,im))* &
          wfnbuf_d(ig,im)
      enddo ! ig
    enddo ! im
  endif ! fblas

#ifdef MPI
  ndum(:)=norm(:)
  call MPI_Allreduce(ndum,norm,nvec, &
    MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

  do im=npre+1,nvec
    norm(im)=sqrt(norm(im))
  enddo ! im

  if (fnorm) then
    if (fblas) then
      do im=npre+1,nvec
        inorm=SCALARIFY(1.0d0/norm(im))
        call X(SCAL)(ngk_l,inorm,wfnbuf_d(:,im),1)
      enddo ! im
    else ! fblas
      do im=npre+1,nvec
        inorm=SCALARIFY(1.0d0/norm(im))
        do ig=1,ngk_l
          wfnbuf_d(ig,im)=wfnbuf_d(ig,im)*inorm
        enddo ! ig
      enddo ! im
    endif ! fblas
  endif ! fnorm

  POP_SUB(npsi)

  return

end subroutine npsi

!===============================================================================

subroutine calc_precond(h_gg, en, precond)

! This subroutine computes a preconditioner for a component of
! a residual vector. On entrance takes a diagonal part of the
! Hamiltonian in h_gg and an eigenvalue in en computed in hdiag.
! On exit returns a preconditioner in precond

  real(DP), intent(in) :: h_gg, en
  real(DP), intent(out) :: precond

  real(DP), parameter :: eps = 1.0d-4
  real(DP) :: del, denom

! do not use PUSH_SUB / POP_SUB as this subroutine is called
! too many times (for each k-point, spin, band, G-vector, iteration)

!  denom = h_gg - en
!  if (abs(denom) .lt. eps) denom = sign(eps, denom)

!  del = h_gg - en
!  denom = sqrt(1.0d0 + del**2)

  del = h_gg - en
  denom = 1.0d0 + del + sqrt(1.0d0 + (del - 1.0d0)**2)

  precond = 1.0d0 / denom

  return

end subroutine calc_precond

end module hdiag_m

