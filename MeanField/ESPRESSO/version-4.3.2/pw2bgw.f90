!!! this file contains post-espresso-5.1 commits: r11037, r11038 !!!

!-------------------------------------------------------------------------------
!
!   pw2bgw.f90
!   $MPIRUN pw2bgw.x -in pw2bgw.in > pw2bgw.out
!   converts Quantum ESPRESSO output files to BerkeleyGW input files
!   based on espresso-4.3.2/PP/pw_export.f90
!   written by G. Samsonidze (April 2008)
!
!   The format of the input file for pw2bgw is described in the sample file pw2bgw.inp.
!
!   NOTE: you CANNOT use USPP, PAW, or spinors for BerkeleyGW!
!
!   real wavefunctions are constructed by applying the Gram-Schmidt process
!   based on paratecSGL/src/para/gwreal.f90
!
!   make sure to obtain the whole sphere of G-vectors from espresso
!   do not use "K_POINTS gamma" as it produces only half the sphere
!   use "K_POINTS { tpiba | automatic | crystal }" even for the
!   Gamma-point calculation
!
!-------------------------------------------------------------------------------
!
!   Compilation: use install.sh, or manually follow these steps:
!
!   copy pw2bgw.f90 and make.bgw to espresso-4.3.2/PP
!   add this line to the end of espresso-4.3.2/PP/Makefile: "include make.bgw"
!   make as follows : cd espresso-4.3.2 ; make pw ph ; cd PP ; make pw2bgw.x
!
!-------------------------------------------------------------------------------
!
!   subroutines
!
!   write_wfng  - generates complex wavefunctions in G-space (normalized to 1)
!   real_wfng   - constructs real wavefunctions by applying the Gram-Schmidt
!                 process (called from write_wfng)
!   write_rhog  - generates complex charge density in G-space
!                 (units of the number of electronic states per unit cell)
!   write_vxcg  - generates complex exchange-correlation potential in G-space
!                 (units of Rydberg) [only local part of Vxc]
!   write_vxc0  - prints complex exchange-correlation potential at G=0
!                 (units of eV) [only local part of Vxc]
!   write_vxc - calculates matrix elements of exchange-correlation potential
!                 in R-space [only local part of Vxc] or G-space [supports
!                 non-local Vxc]. In units of eV.
!#BEGIN_INTERNAL_ONLY
!                 Also used for matrix elements
!                 from DFPT calculation (R is only dvscf, G includes bare too).
!   check_dfpt_modes - determines how many modes are available in dvscf file for DFPT
!   initialize_phonon - performs some initialization steps for DFPT for ionic perturbations
!   initialize_electric - performs some initialization steps for DFPT for electric perturbations
!   dfpt_dwf - alternate method of calculating DFPT matrix elements, for testing
!#END_INTERNAL_ONLY
!   write_vnlg  - generates Kleinman-Bylander projectors in G-space
!                 (units of Rydberg)
!   check_inversion - checks whether real/complex version is appropriate
!   write_header - writes header for VXC and RHO files
!   k_pools - determines start and end k-points for a pool
!   set_spin - sets parameters associated with LSDA and spinor cases
!
!-------------------------------------------------------------------------------
!
!   Quantum ESPRESSO stores the wavefunctions in is-ik-ib-ig order
!   BerkeleyGW stores the wavefunctions in ik-ib-is-ig order
!   the outer loop is over is(QE)/ik(BGW) and the inner loop is over ig
!   ik = k-point index, is = spin index, ib = band index, ig = G-vector index
!
!   write_wfng reverts the order of is and ik using smap and kmap arrays,
!   distributes wavefunctions over processors by ig (either real case or
!   spin-polarized case), calls real_wfng that applies the Gram-Schmidt
!   process (real case), reverts the order of is and ib (spin-polarized
!   case), and writes wavefunctions to disk
!
!-------------------------------------------------------------------------------

#include "f_defs.h"

module pw2bgw_m

  implicit none

  private

  public :: &
    write_wfng, &
    real_wfng, &
    write_rhog, &
    write_vxcg, &
    write_vxc0, &
    write_vxc, &
    write_vnlg, &
    check_inversion, &
    check_dfpt_modes, &
    initialize_phonon, &
    initialize_electric, &
    dfpt_dwf

contains

!-------------------------------------------------------------------------------

SUBROUTINE write_wfng ( output_file_name, real_or_complex, &
  wfng_kgrid, wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, &
  wfng_dk3, wfng_occupation, wfng_nvmin, wfng_nvmax )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, &
    ibrav, symm_type
  USE constants, ONLY : pi, tpi, eps6
  USE grid_dimensions, ONLY : nr1, nr2, nr3
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, g, mill, ecutrho
  USE io_files, ONLY : iunwfc, nwordwfc
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat, atm, ityp, tau
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, wk, ngk, nks, nkstot
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum, mp_max, mp_get, mp_bcast, mp_barrier
  USE mp_global, ONLY : mpime, nproc, world_comm, me_pool, &
    root_pool, npool, nproc_pool, intra_pool_comm
  USE mp_wave, ONLY : mergewf
  USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
  USE symm_base, ONLY : s, ftau, nsym
  USE wavefunctions_module, ONLY : evc
  USE wvfct, ONLY : npwx, nbnd, npw, et, wg, g2kin, ecutwfc
#ifdef __PARA
  USE parallel_include, ONLY : MPI_DOUBLE_COMPLEX
#endif

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  logical, intent (in) :: wfng_kgrid
  integer, intent (in) :: wfng_nk1
  integer, intent (in) :: wfng_nk2
  integer, intent (in) :: wfng_nk3
  real (DP), intent (in) :: wfng_dk1
  real (DP), intent (in) :: wfng_dk2
  real (DP), intent (in) :: wfng_dk3
  logical, intent (in) :: wfng_occupation
  integer, intent (in) :: wfng_nvmin
  integer, intent (in) :: wfng_nvmax

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  logical :: proc_wf, bad_kgrid
  integer :: unit, i, j, k, cell_symmetry, nrecord
  integer :: id, ib, ik, iks, ike, is, ig, ierr
  integer :: nd, ntran, nb, nk_l, nk_g, ns, nst, nsf, ng_l, ng_g
  integer :: ngg, npw_g, npwx_g
  integer :: local_pw, ipsour, igwx, ngkdist_g, ngkdist_l, npw_g2
  real (DP) :: alat2, recvol, dr1, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: kmap ( : )
  integer, allocatable :: smap ( : )
  integer, allocatable :: ifmin ( : )
  integer, allocatable :: ifmax ( : )
  integer, allocatable :: itmp ( : )
  integer, allocatable :: ngk_g ( : )
  integer, allocatable :: ipmask ( : )
  integer, allocatable :: igwk ( : )
  integer, allocatable :: igwf_l2g ( : )
  integer, allocatable :: g_g ( :, : )
  integer, allocatable :: igk_l2g ( :, : )
  real (DP), allocatable :: et_g ( :, : )
  real (DP), allocatable :: wg_g ( :, : )
  real (DP), allocatable :: energy ( :, : )
  complex (DP), allocatable :: wfng ( : )
  complex (DP), allocatable :: wfng_buf ( :, : )
  complex (DP), allocatable :: wfng_dist ( :, :, : )

  INTEGER, EXTERNAL :: atomic_number

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true. )

  IF ( real_or_complex .EQ. 1 .OR. nspin .GT. 1 ) THEN
    proc_wf = .TRUE.
  ELSE
    proc_wf = .FALSE.
  ENDIF

  bad_kgrid = .FALSE.
  IF ( wfng_kgrid ) THEN
    IF ( wfng_nk1 .LE. 0 .OR. wfng_nk2 .LE. 0 .OR. wfng_nk3 .LE. 0 ) &
      bad_kgrid = .TRUE.
  ELSE
    IF ( nk1 .LE. 0 .OR. nk2 .LE. 0 .OR. nk3 .LE. 0 ) &
      bad_kgrid = .TRUE.
  ENDIF
  IF ( bad_kgrid .AND. ionode ) THEN
    WRITE ( 0, 101 )
  ENDIF

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("WFN-Real",24X)' )
  ELSE
    WRITE ( stitle, '("WFN-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1
  nd = 3

  nb = nbnd
  nk_l = nks
  nk_g = nkstot
  call set_spin(ns,nst,nsf,nspin)

  ng_l = ngm
  ng_g = ngm_g

  CALL k_pools(iks, ike)

  ALLOCATE ( kmap ( nk_g ) )
  ALLOCATE ( smap ( nk_g ) )

  DO i = 1, nk_g
    j = ( i - 1 ) / ns
    k = i - 1 - j * ns
    kmap ( i ) = j + k * ( nk_g / ns ) + 1
    smap ( i ) = k + 1
  ENDDO
  ierr = 0
  DO i = 1, nk_g
    ik = kmap ( i )
    is = smap ( i )
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. is .NE. isk ( ik ) ) &
      ierr = ierr + 1
  ENDDO
  CALL mp_max ( ierr )
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_wfng', 'smap', ierr )

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( ibrav .GE. 1 .AND. ibrav .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( ibrav .GE. 4 .AND. ibrav .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( ibrav .GE. 6 .AND. ibrav .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_wfng', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2, dr1 )
    t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( nr1 )
    t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( nr2 )
    t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( nr3 )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  ALLOCATE ( et_g ( nb, nk_g ) )

  DO ik = 1, nk_l
    DO ib = 1, nb
      et_g ( ib, ik ) = et ( ib, ik )
    ENDDO
  ENDDO
#ifdef __PARA
  CALL poolrecover ( et_g, nb, nk_g, nk_l )
  CALL mp_bcast ( et_g, ionode_id )
#endif

  ALLOCATE ( wg_g ( nb, nk_g ) )
  ALLOCATE ( ifmin ( nk_g ) )
  ALLOCATE ( ifmax ( nk_g ) )

  IF ( wfng_occupation ) THEN

    DO ik = 1, nk_g
      DO ib = 1, nb
        IF ( ib .GE. wfng_nvmin .AND. ib .LE. wfng_nvmax ) THEN
          wg_g ( ib, ik ) = 1.0D0
        ELSE
          wg_g ( ib, ik ) = 0.0D0
        ENDIF
      ENDDO
    ENDDO
    DO ik = 1, nk_g
      ifmin ( ik ) = wfng_nvmin
    ENDDO
    DO ik = 1, nk_g
      ifmax ( ik ) = wfng_nvmax
    ENDDO

  ELSE

    DO ik = 1, nk_l
      DO ib = 1, nb
        wg_g ( ib, ik ) = wg ( ib, ik ) 
        IF ( abs ( wk ( ik ) ) .GT. eps6 ) THEN
          wg_g ( ib, ik ) = wg_g ( ib, ik ) / wk ( ik ) 
        ENDIF
      ENDDO
    ENDDO
#ifdef __PARA
    CALL poolrecover ( wg_g, nb, nk_g, nk_l )
#endif
    DO ik = 1, nk_g
      ifmin ( ik ) = 0
    ENDDO
    DO ik = 1, nk_g
      ifmax ( ik ) = 0
    ENDDO
    DO ik = 1, nk_g
      DO ib = 1, nb
        IF ( wg_g( ib, ik ) .GT. 0.5D0 ) THEN
          IF ( ifmin ( ik ) .EQ. 0 ) ifmin ( ik ) = ib
          ifmax ( ik ) = ib
        ENDIF
      ENDDO
    ENDDO

  ENDIF

  ALLOCATE ( g_g ( nd, ng_g ) )

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO
  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  CALL mp_sum ( g_g, intra_pool_comm )

  ALLOCATE ( igk_l2g ( npwx, nk_l ) )

  ALLOCATE ( itmp ( npwx ) )
  DO ik = 1, nk_l
    DO i = 1, npwx
      itmp ( i ) = 0
    ENDDO
    npw = npwx
    CALL gk_sort ( xk ( 1, ik + iks - 1 ), ng_l, g, ecutwfc / tpiba2, &
      npw, itmp ( 1 ), g2kin )
    DO ig = 1, npw
      igk_l2g ( ig, ik ) = ig_l2g ( itmp ( ig ) )
    ENDDO
    DO ig = npw + 1, npwx
      igk_l2g ( ig, ik ) = 0
    ENDDO
    ngk ( ik ) = npw
  ENDDO
  DEALLOCATE ( itmp )

  ALLOCATE ( ngk_g ( nk_g ) )

  DO ik = 1, nk_g
    ngk_g ( ik ) = 0
  ENDDO
  DO ik = 1, nk_l
    ngk_g ( ik + iks - 1 ) = ngk ( ik )
  ENDDO
  CALL mp_sum ( ngk_g )

  npw_g = MAXVAL ( igk_l2g ( :, : ) )
  CALL mp_max ( npw_g )

  npwx_g = MAXVAL ( ngk_g ( : ) )

  CALL cryst_to_cart ( nk_g / ns, xk, at, - 1 )

  IF ( ionode ) THEN
    OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )
    WRITE ( unit ) stitle, sdate, stime
    WRITE ( unit ) nspin, ng_g, ntran, cell_symmetry, nat, ecutrho, &
      nk_g / ns, nb, npwx_g, ecutwfc
    IF ( wfng_kgrid ) THEN
      WRITE ( unit ) nr1, nr2, nr3, wfng_nk1, wfng_nk2, wfng_nk3, &
        wfng_dk1, wfng_dk2, wfng_dk3
    ELSE
      WRITE ( unit ) nr1, nr2, nr3, nk1, nk2, nk3, &
        0.5D0 * dble ( k1 ), 0.5D0 * dble ( k2 ), 0.5D0 * dble ( k3 )
    ENDIF
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nk_g / ns )
    WRITE ( unit ) ( wk ( ik ) * dble ( nst ) / 2.0D0, ik = 1, nk_g / ns )
    WRITE ( unit ) ( ( xk ( id, ik ), id = 1, nd ), ik = 1, nk_g / ns )
    WRITE ( unit ) ( ifmin ( ik ), ik = 1, nk_g )
    WRITE ( unit ) ( ifmax ( ik ), ik = 1, nk_g )
    WRITE ( unit ) ( ( et_g ( ib, ik ), ib = 1, nb ), ik = 1, nk_g )
    WRITE ( unit ) ( ( wg_g ( ib, ik ), ib = 1, nb ), ik = 1, nk_g )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
  ENDIF

  CALL cryst_to_cart ( nk_g / ns, xk, bg, 1 )

  DEALLOCATE ( wg_g )
  DEALLOCATE ( ifmax )
  DEALLOCATE ( ifmin )

  ALLOCATE ( igwk ( npwx_g ) )

  IF ( proc_wf ) THEN
    IF ( nspin == 4) THEN
      IF ( MOD ( npwx_g, nproc ) .EQ. 0 ) THEN
        ngkdist_l = ( 2 * npwx_g ) / nproc
      ELSE
        ngkdist_l = ( 2 * npwx_g ) / nproc + 1
      ENDIF
    ELSE
      IF ( MOD ( npwx_g, nproc ) .EQ. 0 ) THEN
        ngkdist_l = npwx_g / nproc
      ELSE
        ngkdist_l = npwx_g / nproc + 1
      ENDIF
    ENDIF
    ngkdist_g = ngkdist_l * nproc
    IF ( real_or_complex .EQ. 1 ) &
    ALLOCATE ( energy ( nb, ns ) )
    ALLOCATE ( wfng_buf ( ngkdist_g, ns ) )
    ALLOCATE ( wfng_dist ( ngkdist_l, nb, ns ) )
  ENDIF

  DO i = 1, nk_g

    ik = kmap ( i )
    is = smap ( i )

    IF ( real_or_complex .EQ. 1 ) THEN
      DO ib = 1, nb
        energy ( ib, is ) = et_g ( ib, i )
      ENDDO
    ENDIF

    DO j = 1, npwx_g
      igwk ( j ) = 0
    ENDDO
    ALLOCATE ( itmp ( npw_g ) )
    DO j = 1, npw_g
      itmp ( j ) = 0
    ENDDO
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      DO ig = 1, ngk ( ik - iks + 1 )
        itmp ( igk_l2g ( ig, ik - iks + 1 ) ) = igk_l2g ( ig, ik - iks + 1 )
      ENDDO
    ENDIF
    CALL mp_sum ( itmp )
    ngg = 0
    DO ig = 1, npw_g
      IF ( itmp ( ig ) .EQ. ig ) THEN
        ngg = ngg + 1
        igwk ( ngg ) = ig
      ENDIF
    ENDDO
    DEALLOCATE ( itmp )

    IF ( ionode ) THEN
      IF ( is .EQ. 1 ) THEN
        WRITE ( unit ) nrecord
        WRITE ( unit ) ngk_g ( ik )
        WRITE ( unit ) ( ( g_g ( id, igwk ( ig ) ), id = 1, nd ), &
          ig = 1, ngk_g ( ik ) )
      ENDIF
    ENDIF

    local_pw = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      CALL davcio ( evc, nwordwfc, iunwfc, ik - iks + 1, - 1 )
      local_pw = ngk ( ik - iks + 1 )
    ENDIF

    ALLOCATE ( igwf_l2g ( local_pw ) )

    DO ig = 1, local_pw
      igwf_l2g ( ig ) = 0
    ENDDO
    DO ig = 1, local_pw
      ngg = igk_l2g ( ig, ik - iks + 1 )
      DO j = 1, ngk_g ( ik )
        IF ( ngg .EQ. igwk ( j ) ) THEN
          igwf_l2g ( ig ) = j
          EXIT
        ENDIF
      ENDDO
    ENDDO

    ALLOCATE ( ipmask ( nproc ) )
    DO j = 1, nproc
      ipmask ( j ) = 0
    ENDDO
    ipsour = ionode_id
    IF ( npool .GT. 1 ) THEN
      IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
        IF ( me_pool .EQ. root_pool ) ipmask ( mpime + 1 ) = 1
      ENDIF
      CALL mp_sum ( ipmask )
      DO j = 1, nproc
        IF ( ipmask ( j ) .EQ. 1 ) ipsour = j - 1
      ENDDO
    ENDIF
    DEALLOCATE ( ipmask )

    igwx = 0
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) &
      igwx = MAXVAL ( igwf_l2g ( 1 : local_pw ) )
    CALL mp_max ( igwx, intra_pool_comm )
    IF ( ipsour .NE. ionode_id ) &
      CALL mp_get ( igwx, igwx, mpime, ionode_id, ipsour, 1 )
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. igwx .NE. ngk_g ( ik ) ) &
      ierr = 1
    CALL mp_max ( ierr )
    IF ( ierr .GT. 0 ) &
      CALL errore ( 'write_wfng', 'igwx ngk_g', ierr )

    IF (nspin == 4) THEN
      igwx = 2 * igwx
    ENDIF
    ALLOCATE ( wfng ( MAX ( 1, igwx ) ) )

    IF (nspin == 4) THEN
      npw_g2 = MAXVAL ( ig_l2g ( 1: local_pw) )
      CALL mp_max (npw_g2)
    ENDIF

    DO ib = 1, nb

      DO j = 1, igwx
        wfng ( j ) = ( 0.0D0, 0.0D0 )
      ENDDO
      IF ( npool .GT. 1 ) THEN
        IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
          CALL mergewf ( evc ( :, ib ), wfng, local_pw, igwf_l2g, &
            me_pool, nproc_pool, root_pool, intra_pool_comm )
        ENDIF
        IF ( ipsour .NE. ionode_id ) THEN
          CALL mp_get ( wfng, wfng, mpime, ionode_id, ipsour, ib, &
            world_comm )
        ENDIF
      ELSE
        CALL mergewf ( evc ( :, ib ), wfng, local_pw, igwf_l2g, &
          mpime, nproc, ionode_id, world_comm )

        IF (nspin == 4) THEN
          CALL mp_barrier ( )
          CALL mergewf ( evc ( npwx_g + 1 : npwx_g + local_pw, ib), & 
            wfng (local_pw + 1 : 2 * local_pw ), &
           local_pw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
        ENDIF
      ENDIF

      IF ( proc_wf ) THEN
        DO ig = 1, igwx
          wfng_buf ( ig, is ) = wfng ( ig )
        ENDDO
        DO ig = igwx + 1, ngkdist_g
          wfng_buf ( ig, is ) = ( 0.0D0, 0.0D0 )
        ENDDO
#ifdef __PARA
        CALL mp_barrier ( )
        CALL MPI_Scatter ( wfng_buf ( :, is ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
        wfng_dist ( :, ib, is ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
        ionode_id, world_comm, ierr )
        IF ( ierr .GT. 0 ) &
          CALL errore ( 'write_wfng', 'mpi_scatter', ierr )
#else
        DO ig = 1, ngkdist_g
          wfng_dist ( ig, ib, is ) = wfng_buf ( ig, is )
        ENDDO
#endif
      ELSE
        IF ( ionode ) THEN
          WRITE ( unit ) nrecord
          WRITE ( unit ) ngk_g ( ik )
          WRITE ( unit ) ( wfng ( ig ), ig = 1, igwx )
        ENDIF
      ENDIF

    ENDDO

    DEALLOCATE ( wfng )
    DEALLOCATE ( igwf_l2g )

    IF ( proc_wf .AND. is .EQ. ns ) THEN
      IF ( real_or_complex .EQ. 1 ) THEN
        CALL start_clock ( 'real_wfng' )
        CALL real_wfng ( ik, ngkdist_l, nb, ns, energy, wfng_dist )
        CALL stop_clock ( 'real_wfng' )
      ENDIF
      DO ib = 1, nb
        DO is = 1, ns
#ifdef __PARA
          CALL mp_barrier ( )
          CALL MPI_Gather ( wfng_dist ( :, ib, is ), ngkdist_l, &
          MPI_DOUBLE_COMPLEX, wfng_buf ( :, is ), ngkdist_l, &
          MPI_DOUBLE_COMPLEX, ionode_id, world_comm, ierr )
          IF ( ierr .GT. 0 ) &
            CALL errore ( 'write_wfng', 'mpi_gather', ierr )
#else
          DO ig = 1, ngkdist_g
            wfng_buf ( ig, is ) = wfng_dist ( ig, ib, is )
          ENDDO
#endif
        ENDDO
        IF ( ionode ) THEN
          WRITE ( unit ) nrecord
          WRITE ( unit ) ngk_g ( ik )
          IF ( real_or_complex .EQ. 1 ) THEN
            WRITE ( unit ) ( ( dble ( wfng_buf ( ig, is ) ), &
              ig = 1, igwx ), is = 1, ns )
          ELSE
            WRITE ( unit ) ( ( wfng_buf ( ig, is ), &
              ig = 1, igwx ), is = 1, ns )
          ENDIF
        ENDIF
      ENDDO
    ENDIF

  ENDDO

  DEALLOCATE ( igwk )
  DEALLOCATE ( ngk_g )
  DEALLOCATE ( igk_l2g )
  DEALLOCATE ( et_g )

  IF ( proc_wf ) THEN
    IF ( real_or_complex .EQ. 1 ) &
    DEALLOCATE ( energy )
    DEALLOCATE ( wfng_buf )
    DEALLOCATE ( wfng_dist )
  ENDIF

  IF ( ionode ) THEN
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( g_g )
  DEALLOCATE ( smap )
  DEALLOCATE ( kmap )

  CALL mp_barrier ( )

  RETURN

101 FORMAT ( /, 5X, "WARNING: kgrid is set to zero in the wavefunction file.", &
             /, 14X, "The resulting file will only be usable as the fine grid in inteqp.", / )

END SUBROUTINE write_wfng

!-------------------------------------------------------------------------------

SUBROUTINE real_wfng ( ik, ngkdist_l, nb, ns, energy, wfng_dist )

  USE constants, ONLY : eps6
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode
  USE mp, ONLY : mp_sum

  IMPLICIT NONE

  integer, intent (in) :: ik, ngkdist_l, nb, ns
  real (DP), intent (in) :: energy ( nb, ns )
  complex (DP), intent (inout) :: wfng_dist ( ngkdist_l, nb, ns )

  real (DP), PARAMETER :: eps2 = 1.0D-2
  real (DP), PARAMETER :: eps5 = 1.0D-5

  character :: tmpstr*80
  integer :: i, j, k, is, ib, jb, ig, inum, deg, mdeg, inc
  integer :: dimension_span, reduced_span, ierr
  real (DP) :: x
  integer, allocatable :: imap ( :, : )
  integer, allocatable :: inums ( : )
  integer, allocatable :: inull ( : )
  integer, allocatable :: null_map ( :, : )
  real (DP), allocatable :: psi ( :, : )
  real (DP), allocatable :: phi ( :, : )
  real (DP), allocatable :: vec ( : )
  complex (DP), allocatable :: wfc ( : )

  mdeg = 1
  DO is = 1, ns
    DO ib = 1, nb - 1
      deg = 1
      DO jb = ib + 1, nb
        IF ( abs ( energy ( ib, is ) - energy ( jb, is ) ) &
          .LT. eps5 * dble ( jb - ib + 1 ) ) deg = deg + 1
      ENDDO
      IF ( deg .GT. mdeg ) mdeg = deg
    ENDDO
  ENDDO
  mdeg = mdeg * 2

  ALLOCATE ( imap ( nb, ns ) )
  ALLOCATE ( inums ( ns ) )
  ALLOCATE ( inull ( nb ) )
  ALLOCATE ( null_map ( mdeg, nb ) )

  DO is = 1, ns
    inum = 1
    DO ib = 1, nb
      IF ( ib .EQ. nb ) THEN
        imap ( inum, is ) = ib
        inum = inum + 1
      ELSEIF ( abs ( energy ( ib, is ) - &
        energy ( ib + 1, is ) ) .GT. eps5 ) THEN
        imap ( inum, is ) = ib
        inum = inum + 1
      ENDIF
    ENDDO
    inum = inum - 1
    inums ( is ) = inum
  ENDDO

  ALLOCATE ( wfc ( ngkdist_l ) )
  ALLOCATE ( psi ( ngkdist_l, mdeg ) )
  ALLOCATE ( phi ( ngkdist_l, mdeg ) )
  ALLOCATE ( vec ( ngkdist_l ) )

  DO is = 1, ns
    inc = 1
    inum = inums ( is )
    DO i = 1, inum
      inull ( i ) = 1
      DO ib = inc, imap ( i, is )
        DO ig = 1, ngkdist_l
          wfc ( ig ) = wfng_dist ( ig, ib, is )
        ENDDO
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + dble ( wfc ( ig ) ) **2
        ENDDO
        CALL mp_sum ( x )
        IF ( x .LT. eps2 ) null_map ( inull ( i ), i ) = 0
        IF ( x .GT. eps2 ) null_map ( inull ( i ), i ) = 1
        inull ( i ) = inull ( i ) + 1
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + aimag ( wfc ( ig ) ) **2
        ENDDO
        CALL mp_sum ( x )
        IF ( x .LT. eps2 ) null_map ( inull ( i ), i ) = 0
        IF ( x .GT. eps2 ) null_map ( inull ( i ), i ) = 1
        inull ( i ) = inull ( i ) + 1
      ENDDO
      inull ( i ) = inull ( i ) - 1
      inc = imap ( i, is ) + 1
    ENDDO
    inc = 1
    ib = 1
    DO i = 1, inum
      k = 1
      DO j = 1, 2 * ( imap ( i, is ) - inc ) + 1, 2
        IF ( null_map ( j, i ) .EQ. 1 .OR. &
          null_map ( j + 1, i ) .EQ. 1 ) THEN
          DO ig = 1, ngkdist_l
            wfc ( ig ) = wfng_dist ( ig, ib, is )
          ENDDO
          IF ( null_map ( j, i ) .EQ. 1 ) THEN
            DO ig = 1, ngkdist_l
              phi ( ig, k ) = dble ( wfc ( ig ) )
            ENDDO
            k = k + 1
          ENDIF
          IF ( null_map ( j + 1, i ) .EQ. 1 ) THEN
            DO ig = 1, ngkdist_l
              phi ( ig, k ) = aimag ( wfc ( ig ) )
            ENDDO
            k = k + 1
          ENDIF
          ib = ib + 1
        ENDIF
      ENDDO
      dimension_span = k - 1
      IF ( dimension_span .EQ. 0 ) THEN
        ierr = 201
        WRITE ( tmpstr, 201 ) ik, is, inc
        CALL errore ( 'real_wfng', tmpstr, ierr )
      ENDIF
      DO j = 1, dimension_span
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + phi ( ig, j ) **2
        ENDDO
        CALL mp_sum ( x )
        x = sqrt ( x )
        DO ig = 1, ngkdist_l
          phi ( ig, j ) = phi ( ig, j ) / x
        ENDDO
      ENDDO
!
! the Gram-Schmidt process begins
!
      reduced_span = 1
      DO ig = 1, ngkdist_l
        psi ( ig, 1 ) = phi ( ig, 1 )
      ENDDO
      DO j = 1, dimension_span - 1
        DO ig = 1, ngkdist_l
          vec ( ig ) = phi ( ig, j + 1 )
        ENDDO
        DO k = 1, reduced_span
          x = 0.0D0
          DO ig = 1, ngkdist_l
            x = x + phi ( ig, j + 1 ) * psi ( ig, k )
          ENDDO
          CALL mp_sum ( x )
          DO ig = 1, ngkdist_l
            vec ( ig ) = vec ( ig ) - psi ( ig, k ) * x
          ENDDO
        ENDDO
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + vec ( ig ) **2
        ENDDO
        CALL mp_sum ( x )
        x = sqrt ( x )
        IF ( x .GT. eps6 ) THEN
          reduced_span = reduced_span + 1
          DO ig = 1, ngkdist_l
            psi ( ig, reduced_span ) = vec ( ig ) / x
          ENDDO
        ENDIF
      ENDDO
!
! the Gram-Schmidt process ends
!
      IF ( reduced_span .LT. imap ( i, is ) - inc + 1 ) THEN
        ierr = 202
        WRITE ( tmpstr, 202 ) ik, is, inc
        CALL errore ( 'real_wfng', tmpstr, ierr )
      ENDIF
      DO ib = inc, imap ( i, is )
        DO ig = 1, ngkdist_l
          wfng_dist ( ig, ib, is ) = &
          CMPLX ( psi ( ig, ib - inc + 1 ), 0.0D0 )
        ENDDO
      ENDDO
      inc = imap ( i, is ) + 1
    ENDDO
  ENDDO

  DEALLOCATE ( vec )
  DEALLOCATE ( phi )
  DEALLOCATE ( psi )
  DEALLOCATE ( wfc )
  DEALLOCATE ( null_map )
  DEALLOCATE ( inull )
  DEALLOCATE ( inums )
  DEALLOCATE ( imap )

  RETURN

201 FORMAT("failed Gram-Schmidt dimension span for kpoint =",i6," spin =",i2," band =",i6)
202 FORMAT("failed Gram-Schmidt reduced span for kpoint =",i6," spin =",i2," band =",i6)

END SUBROUTINE real_wfng

!-------------------------------------------------------------------------------

SUBROUTINE write_rhog ( output_file_name, real_or_complex )

  USE cell_base, ONLY : omega
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, atm, ityp, tau 
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_global, ONLY : intra_pool_comm
  USE scf, ONLY : rho
  USE symm_base, ONLY : s, ftau, nsym

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex

  character :: stitle*32
  integer :: unit, is, ig
  integer :: ns, nst, nsf, ng_l, ng_g
  integer :: nrecord
  complex (DP), allocatable :: rhog_g ( :, : )

  INTEGER, EXTERNAL :: atomic_number

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true. )

  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("RHO-Real",24X)' )
  ELSE
    WRITE ( stitle, '("RHO-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1

  ng_l = ngm
  ng_g = ngm_g
  call set_spin(ns,nst,nsf,nspin)

  ALLOCATE ( rhog_g ( ng_g, ns ) )

  DO is = 1, ns
    DO ig = 1, ng_g
      rhog_g ( ig, is ) = ( 0.0D0, 0.0D0 )
    ENDDO
  ENDDO

  DO is = 1, ns
    DO ig = 1, ng_l
      rhog_g ( ig_l2g ( ig ), is ) = rho%of_g ( ig, is )
    ENDDO
  ENDDO

  CALL mp_sum ( rhog_g, intra_pool_comm )

  DO is = 1, ns
    DO ig = 1, ng_g
      rhog_g ( ig, is ) = rhog_g ( ig, is ) * CMPLX ( omega, 0.0D0 )
    ENDDO
  ENDDO

  IF ( ionode ) OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )

  CALL write_header(unit, stitle)

  IF ( ionode ) THEN
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    IF ( real_or_complex .EQ. 1 ) THEN
      WRITE ( unit ) ( ( dble ( rhog_g ( ig, is ) ), &
        ig = 1, ng_g ), is = 1, ns )
    ELSE
      WRITE ( unit ) ( ( rhog_g ( ig, is ), &
        ig = 1, ng_g ), is = 1, ns )
    ENDIF
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( rhog_g )

  RETURN

END SUBROUTINE write_rhog

!-------------------------------------------------------------------------------

SUBROUTINE write_vxcg ( output_file_name, real_or_complex, &
  vxc_zero_rho_core, input_dft, exx_flag )

  USE cell_base, ONLY : omega
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE funct, ONLY : enforce_input_dft
  USE grid_dimensions, ONLY : nrxx
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, nl, mill, ecutrho
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, atm, ityp, tau 
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_global, ONLY : intra_pool_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE symm_base, ONLY : s, ftau, nsym
  USE wavefunctions_module, ONLY : psic
#ifdef EXX
  USE funct, ONLY : start_exx, stop_exx
#endif

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  logical, intent (in) :: vxc_zero_rho_core
  character ( len = 20 ), intent (in) :: input_dft
  logical, intent (in) :: exx_flag

  character :: stitle*32
  integer :: unit, is, ir, ig
  integer :: ns, nst, nsf, nr, ng_l, ng_g
  integer :: nrecord
  real (DP), allocatable :: vxcr_g ( :, : )
  complex (DP), allocatable :: vxcg_g ( :, : )

  INTEGER, EXTERNAL :: atomic_number

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true. )

  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("VXC-Real",24X)' )
  ELSE
    WRITE ( stitle, '("VXC-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1

  call set_spin(ns,nst,nsf,nspin)
  nr = nrxx
  ng_l = ngm
  ng_g = ngm_g

  ALLOCATE ( vxcr_g ( nr, nsf ) )
  ALLOCATE ( vxcg_g ( ng_g, nsf ) )

  DO is = 1, nsf
    DO ig = 1, ng_g
      vxcg_g ( ig, is ) = ( 0.0D0, 0.0D0 )
    ENDDO
  ENDDO

  CALL enforce_input_dft ( input_dft )
#ifdef EXX
  IF ( exx_flag ) CALL start_exx ( )
#endif
  vxcr_g ( :, : ) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  CALL v_xc ( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g )
#ifdef EXX
  IF ( exx_flag ) CALL stop_exx ( )
#endif
  DO is = 1, ns
    DO ir = 1, nr
      psic ( ir ) = CMPLX ( vxcr_g ( ir, is ), 0.0D0 )
    ENDDO
    CALL fwfft ( 'Dense', psic, dfftp )
    DO ig = 1, ng_l
      vxcg_g ( ig_l2g ( ig ), is ) = psic ( nl ( ig ) )
    ENDDO
  ENDDO

  CALL mp_sum ( vxcg_g, intra_pool_comm )

  IF ( ionode ) OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )

  CALL write_header(unit, stitle)

  IF ( ionode ) THEN
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    IF ( real_or_complex .EQ. 1 ) THEN
      WRITE ( unit ) ( ( dble ( vxcg_g ( ig, is ) ), &
        ig = 1, ng_g ), is = 1, ns )
    ELSE
      WRITE ( unit ) ( ( vxcg_g ( ig, is ), ig = 1, ng_g ), is = 1, nsf )
    ENDIF
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( vxcg_g )
  DEALLOCATE ( vxcr_g )

  RETURN

END SUBROUTINE write_vxcg

!-------------------------------------------------------------------------------

SUBROUTINE write_vxc0 ( output_file_name, vxc_zero_rho_core, input_dft, &
  exx_flag )

  USE constants, ONLY : RYTOEV
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE funct, ONLY : enforce_input_dft
  USE grid_dimensions, ONLY : nrxx
  USE gvect, ONLY : ngm, nl, mill
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_global, ONLY : intra_pool_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions_module, ONLY : psic
#ifdef EXX
  USE funct, ONLY : start_exx, stop_exx
#endif

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  logical, intent (in) :: vxc_zero_rho_core
  character ( len = 20 ), intent (in) :: input_dft
  logical, intent (in) :: exx_flag

  integer :: unit
  integer :: is, ir, ig
  integer :: ns, nst, nsf, nr, ng_l
  real (DP), allocatable :: vxcr_g ( :, : )
  complex (DP), allocatable :: vxc0_g ( : )

  unit = 4

  nr = nrxx
  ng_l = ngm
  call set_spin(ns,nst,nsf,nspin)

  ALLOCATE ( vxcr_g ( nr, nsf ) )
  ALLOCATE ( vxc0_g ( nsf ) )

  DO is = 1, nsf
    vxc0_g ( is ) = ( 0.0D0, 0.0D0 )
  ENDDO

  CALL enforce_input_dft ( input_dft )
#ifdef EXX
  IF ( exx_flag ) CALL start_exx ( )
#endif
  vxcr_g ( :, : ) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  CALL v_xc ( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g )
#ifdef EXX
  IF ( exx_flag ) CALL stop_exx ( )
#endif
  DO is = 1, nsf
    DO ir = 1, nr
      psic ( ir ) = CMPLX ( vxcr_g ( ir, is ), 0.0D0 )
    ENDDO
    CALL fwfft ( 'Dense', psic, dfftp )
    DO ig = 1, ng_l
      IF ( mill ( 1, ig ) .EQ. 0 .AND. mill ( 2, ig ) .EQ. 0 .AND. &
        mill ( 3, ig ) .EQ. 0 ) vxc0_g ( is ) = psic ( nl ( ig ) )
    ENDDO
  ENDDO

  CALL mp_sum ( vxc0_g, intra_pool_comm )

  DO is = 1, nsf
    vxc0_g ( is ) = vxc0_g ( is ) * CMPLX ( RYTOEV, 0.0D0 )
  ENDDO

  IF ( ionode ) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    WRITE ( unit, 101 )
    DO is = 1, nsf
      WRITE ( unit, 102 ) is, vxc0_g ( is )
    ENDDO
    WRITE ( unit, 103 )
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  DEALLOCATE ( vxcr_g )
  DEALLOCATE ( vxc0_g )

  RETURN

101 FORMAT ( /, 5X, "--------------------------------------------", &
             /, 5X, "spin    Re Vxc(G=0) (eV)    Im Vxc(G=0) (eV)", &
             /, 5X, "--------------------------------------------" )
102 FORMAT ( 5X, I1, 3X, 2F20.15 )
103 FORMAT ( 5X, "--------------------------------------------", / )

END SUBROUTINE write_vxc0

!-------------------------------------------------------------------------------

SUBROUTINE write_vxc (vxc_integral, output_file_name, diag_nmin, diag_nmax, &
  offdiag_nmin, offdiag_nmax, vxc_zero_rho_core, input_dft, exx_flag, &
  dfpt_type, dfpt_mode)

  USE becmod, ONLY : bec_type, allocate_bec_type, deallocate_bec_type
  USE constants, ONLY : rytoev, eps6
  USE control_ph, ONLY : nbnd_occ ! initialize!
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc
  USE eqv, ONLY : dpsi, dvpsi
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE funct, ONLY : enforce_input_dft
  USE grid_dimensions, ONLY : nr1x, nr2x, nr3x, nrxx
  USE gvect, ONLY : ngm, g, nl
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nkstot
  USE lsda_mod, ONLY : nspin, isk
  USE modes, ONLY : u
  USE mp, ONLY : mp_sum
  USE mp_global, ONLY : intra_pool_comm, inter_pool_comm
  USE phus, ONLY : becp1
  USE units_ph, ONLY : iudvscf, lrdrho
  USE scf, ONLY : rho, rho_core, rhog_core
  USE uspp, ONLY : nkb
  USE wavefunctions_module, ONLY : evc, psic
  USE wvfct, ONLY : npwx, npw, nbnd, igk, g2kin, ecutwfc
#ifdef EXX
  USE funct, ONLY : start_exx, stop_exx, exx_is_active
  USE exx, ONLY : vexx
#endif

  IMPLICIT NONE

  character (len = 1), intent (in) :: vxc_integral ! 'g' = G-space, 'r' = R-space
  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin
  integer, intent (inout) :: offdiag_nmax
  logical, intent (in) :: vxc_zero_rho_core
  character (len = 20), intent (in) :: input_dft
  logical, intent (in) :: exx_flag
  integer, intent (in) :: dfpt_type ! 0 = none, 1 = ionic, 2 = electric
  integer, intent (in) :: dfpt_mode

  integer :: ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2
  integer :: ns, nst, nsf
  real(DP) :: dummyr, maxvaldvscf
  complex (DP) :: dummy
  complex (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  real (DP), allocatable :: vxcr (:, :)
  complex (DP), allocatable :: dvscf (:, :)
  complex (DP), allocatable :: psic2 (:)
  complex (DP), allocatable :: hpsi (:)
  type(bec_type) :: becp2 ! we do not care, but need to pass this

  IF ( vxc_integral /= 'r' .AND. vxc_integral /= 'g' ) &
    CALL errore ( 'write_vxc', 'vxc_integral', 1 )

  IF ( nspin .EQ. 4) CALL errore ( 'write_vxc', 'nspin', 4 )

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_vxc', 'diag_nmin > diag_nmax', diag_nmin )
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    if(ionode) write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF

  if(offdiag_nmin > offdiag_nmax) then
    call errore ( 'write_vxc', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
  endif
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    if(ionode) write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    offdiag_nmax = nbnd
  ENDIF

  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)
  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  unit = 4

  call set_spin(ns, nst, nsf, nspin)
  CALL k_pools(iks, ike)

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0)
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (nrxx, nsf))
  IF (noffdiag .GT. 0) ALLOCATE (psic2 (nrxx))
  if(vxc_integral == 'g') ALLOCATE (hpsi (nrxx))

  if(dfpt_type == 0) then

    CALL enforce_input_dft (input_dft)
#ifdef EXX
    IF (exx_flag) CALL start_exx ()
#endif

    vxcr (:, :) = 0.0D0
    IF ( vxc_zero_rho_core ) THEN
      rho_core ( : ) = 0.0D0
      rhog_core ( : ) = ( 0.0D0, 0.0D0 )
    ENDIF
    CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr)
  
  else

    IF(dfpt_type == 1) WRITE(6,'(a,i6)') 'Calculating ionic matrix elements for mode ', dfpt_mode
    IF(dfpt_type == 2) WRITE(6,'(a,i6)') 'Calculating electric matrix elements for polarization ', dfpt_mode
    IF(vxc_integral == 'r') WRITE(6,'(a)') 'Only the Hartree and XC parts of the perturbation are being included.'

    ! iudvscf, lrdrho set in check_dfpt_modes
    ALLOCATE (dvscf (nrxx, nsf))
    CALL davcio_drho(dvscf, lrdrho, iudvscf, dfpt_mode, - 1 )
    write(6,'(a)') 'Read dvscf from disk.'

    ! dvscf is real at q=0
    if(any(aimag(dvscf(:,:)) > eps6)) then
      write(6,'(a)') 'This mode has imaginary dvscf (max = ', maxval(abs(aimag(dvscf(:,:)))), ') so not calculating.'
      RETURN
    endif
    vxcr(:,:) = dble(dvscf(:,:))

    maxvaldvscf = maxval(abs(dvscf(:,:)))
    DEALLOCATE(dvscf)

    if(maxvaldvscf < 1e-9) then
      WRITE(6,'(a)') 'This mode has dvscf = 0 (max abs = ', maxvaldvscf, ') so not calculating.'
      RETURN
    endif

    if(dfpt_type == 2) call allocate_bec_type (nkb, nbnd, becp2)
  endif

  DO ik = iks, ike
    CALL gk_sort (xk (1, ik - iks + 1), ngm, g, ecutwfc &
      / tpiba2, npw, igk, g2kin)
    CALL davcio (evc, nwordwfc, iunwfc, ik - iks + 1, -1)

    if(vxc_integral == 'g') then
      if(dfpt_type == 1) then
        ! initialize dvpsi: bare ionic pert applied to all the wfns 'evc'
        ! gradient of local and non-local part of pseudopotential
        dvpsi(:,:) = (0d0, 0d0)
        call dvqpsi_us(ik, u (1, dfpt_mode), .false.)
      else if(dfpt_type == 2) then
        ! for electric field, from PH/dvpsi_e. check initialization of this stuff!
        ! calculate the commutator [H,x_ipol]  psi > and store it in dvpsi
        ! dpsi used as workspace
        dpsi(:,:) = (0d0, 0d0)
        dvpsi(:,:) = (0d0, 0d0)
        call commutator_Hx_psi (ik, nbnd_occ(ik), becp1(ik), becp2, dfpt_mode, dvpsi, dpsi )
      endif
    endif

    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (nl (igk (ig))) = evc (ig, ib)
        ENDDO
        CALL invfft ('Dense', psic, dfftp)

        if(vxc_integral == 'g') then
          DO ir = 1, nrxx
            psic (ir) = psic (ir) * vxcr (ir, isk (ik))
          ENDDO
          CALL fwfft ('Dense', psic, dfftp)

          if(dfpt_type == 1) then
            ! add bare part of ionic perturbation
            DO ig = 1, npw
              psic (nl (igk(ig))) = psic (nl (igk(ig))) + dvpsi(ig, ib)
            ENDDO
          endif

          hpsi (:) = (0.0D0, 0.0D0)
          DO ig = 1, npw
            hpsi (ig) = psic (nl (igk (ig)))
          ENDDO
          psic (:) = (0.0D0, 0.0D0)
          DO ig = 1, npw
            psic (ig) = evc (ig, ib)
          ENDDO
#ifdef EXX
          IF (exx_is_active () .and. dfpt_type == 0) &
            CALL vexx (npwx, npw, 1, psic, hpsi)
#endif
          dummy = (0.0D0, 0.0D0)
          DO ig = 1, npw
            dummy = dummy + conjg (psic (ig)) * hpsi (ig)
          ENDDO
          dummy = dummy * rytoev
          CALL mp_sum (dummy, intra_pool_comm)
          mtxeld (ib - diag_nmin + 1, ik) = dummy
        else
          dummyr = 0.0D0
          DO ir = 1, nrxx
            dummyr = dummyr + vxcr (ir, isk (ik)) &
              * (dble (psic (ir)) **2 + aimag (psic (ir)) **2)
          ENDDO
          dummyr = dummyr * rytoev / dble (nr1x * nr2x * nr3x)
          CALL mp_sum (dummyr, intra_pool_comm)
          mtxeld (ib - diag_nmin + 1, ik) = dummyr
        endif
      ENDDO
    ENDIF
    IF (noffdiag .GT. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (nl (igk (ig))) = evc (ig, ib)
        ENDDO
        CALL invfft ('Dense', psic, dfftp)

        if(vxc_integral == 'g') then
          DO ir = 1, nrxx
            psic (ir) = psic (ir) * vxcr (ir, isk (ik))
          ENDDO
          CALL fwfft ('Dense', psic, dfftp)

          if(dfpt_type == 1) then
            ! add bare part of ionic perturbation
            DO ig = 1, npw
              psic (nl (igk(ig))) = psic (nl (igk(ig))) + dvpsi(ig, ib)
            ENDDO
          endif

          hpsi (:) = (0.0D0, 0.0D0)
          DO ig = 1, npw
            hpsi (ig) = psic (nl (igk (ig)))
          ENDDO
          psic (:) = (0.0D0, 0.0D0)
          DO ig = 1, npw
            psic (ig) = evc (ig, ib)
          ENDDO
#ifdef EXX
          IF (exx_is_active () .and. dfpt_type == 0) &
            CALL vexx (npwx, npw, 1, psic, hpsi)
#endif
        endif
        DO ib2 = offdiag_nmin, offdiag_nmax
          psic2 (:) = (0.0D0, 0.0D0)
          dummy = (0.0D0, 0.0D0)
          if(vxc_integral == 'g') then
            DO ig = 1, npw
              psic2 (ig) = evc (ig, ib2)
            ENDDO
            DO ig = 1, npw
              dummy = dummy + conjg (psic2 (ig)) * hpsi (ig)
            ENDDO
          else
            DO ig = 1, npw
              psic2 (nl (igk (ig))) = evc (ig, ib2)
            ENDDO
            CALL invfft ('Dense', psic2, dfftp)
            DO ir = 1, nrxx
              dummy = dummy + CMPLX (vxcr (ir, isk (ik)), 0.0D0) &
                * conjg (psic2 (ir)) * psic (ir)
            ENDDO
            dummy = dummy / dble (nr1x * nr2x * nr3x)
          endif
          dummy = dummy * rytoev
          CALL mp_sum (dummy, intra_pool_comm)
          mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, ik) = dummy
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  ! from dvpsi_e
  if(dfpt_type == 2 .and. nkb > 0) call deallocate_bec_type (becp2)

#ifdef EXX
  IF (exx_flag) CALL stop_exx ()
#endif

  DEALLOCATE (vxcr)
  IF (noffdiag .GT. 0) DEALLOCATE (psic2)
  if(vxc_integral == 'g') DEALLOCATE (hpsi)

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    DO ik = 1, nkstot / ns
      WRITE (unit, 101) xk(:, ik), ns * ndiag, &
        ns * noffdiag **2
      DO is = 1, ns
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            WRITE (unit, 102) is, ib, mtxeld &
              (ib - diag_nmin + 1, ik + (is - 1) * nkstot / ns)
          ENDDO
        ENDIF
        IF (noffdiag .GT. 0) THEN
          DO ib = offdiag_nmin, offdiag_nmax
            DO ib2 = offdiag_nmin, offdiag_nmax
              WRITE (unit, 103) is, ib2, ib, mtxelo &
                (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                ik + (is - 1) * nkstot / ns)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

  101 FORMAT (3F13.9, 2I8)
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_vxc

!-------------------------------------------------------------------------------

SUBROUTINE write_vnlg (output_file_name, wfng_kgrid, wfng_nk1, wfng_nk2, &
  wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3)

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, &
    ibrav, symm_type
  USE constants, ONLY : pi, tpi, eps6
  USE grid_dimensions, ONLY : nr1, nr2, nr3
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, g, mill, ecutrho
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat, atm, ityp, tau, nsp
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, wk, ngk, nks, nkstot
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum, mp_max, mp_get, mp_barrier
  USE mp_global, ONLY : mpime, nproc, world_comm, me_pool, &
    root_pool, npool, nproc_pool, intra_pool_comm
  USE mp_wave, ONLY : mergewf
  USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
  USE symm_base, ONLY : s, ftau, nsym
  USE uspp, ONLY : nkb, vkb, deeq
  USE uspp_param, ONLY : nhm, nh
  USE wvfct, ONLY : npwx, npw, g2kin, ecutwfc

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  logical, intent (in) :: wfng_kgrid
  integer, intent (in) :: wfng_nk1
  integer, intent (in) :: wfng_nk2
  integer, intent (in) :: wfng_nk3
  real (DP), intent (in) :: wfng_dk1
  real (DP), intent (in) :: wfng_dk2
  real (DP), intent (in) :: wfng_dk3

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  integer :: i, j, k, ierr, ik, is, ig, ikb, iat, isp, ih, jh, &
    unit, iks, ike, npw_g, npwx_g, ngg, ipsour, &
    igwx, local_pw, id, nd, ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, dr1, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: kmap ( : )
  integer, allocatable :: smap ( : )
  integer, allocatable :: gvec ( :, : )
  integer, allocatable :: ngk_g ( : )
  integer, allocatable :: igk_l2g ( :, : )
  integer, allocatable :: itmp ( : )
  integer, allocatable :: igwk ( : )
  integer, allocatable :: igwf_l2g ( : )
  integer, allocatable :: ipmask ( : )
  complex (DP), allocatable :: vkb_g ( : )

  INTEGER, EXTERNAL :: atomic_number

  IF ( nkb == 0 ) RETURN

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  WRITE ( stitle, '("VKB-Complex",21X)' )

  IF ( nspin .EQ. 4) CALL errore ( 'write_vnlg', 'nspin', 4)

  unit = 4
  nrecord = 1
  nd = 3

  CALL k_pools(iks, ike)

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( ibrav .GE. 1 .AND. ibrav .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( ibrav .GE. 4 .AND. ibrav .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( ibrav .GE. 6 .AND. ibrav .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vnlg', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2, dr1 )
    t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( nr1 )
    t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( nr2 )
    t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( nr3 )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ALLOCATE ( kmap ( nkstot ) )
  ALLOCATE ( smap ( nkstot ) )

  DO i = 1, nkstot
    j = ( i - 1 ) / nspin
    k = i - 1 - j * nspin
    kmap ( i ) = j + k * ( nkstot / nspin ) + 1
    smap ( i ) = k + 1
  ENDDO
  ierr = 0
  DO i = 1, nkstot
    ik = kmap ( i )
    is = smap ( i )
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. is .NE. isk ( ik ) ) &
      ierr = ierr + 1
  ENDDO
  CALL mp_max ( ierr )
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vnlg', 'smap', ierr )

  ALLOCATE ( gvec ( 3, ngm_g ) )
  gvec = 0
  DO ig = 1, ngm
    gvec ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    gvec ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    gvec ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  CALL mp_sum ( gvec, intra_pool_comm )

  ALLOCATE ( ngk_g ( nkstot ) )
  ALLOCATE ( igk_l2g ( npwx, nks ) )
  ngk_g = 0
  igk_l2g = 0
  ALLOCATE ( itmp ( npwx ) )
  DO ik = 1, nks
    itmp = 0
    npw = npwx
    CALL gk_sort ( xk ( 1, ik + iks - 1 ), ngm, g, ecutwfc / tpiba2, &
      npw, itmp ( 1 ), g2kin )
    DO ig = 1, npw
      igk_l2g ( ig, ik ) = ig_l2g ( itmp ( ig ) )
    ENDDO
    ngk ( ik ) = npw
  ENDDO
  DEALLOCATE ( itmp )
  DO ik = 1, nks
    ngk_g ( ik + iks - 1 ) = ngk ( ik )
  ENDDO
  CALL mp_sum ( ngk_g )
  npw_g = MAXVAL ( igk_l2g ( :, : ) )
  CALL mp_max ( npw_g )
  npwx_g = MAXVAL ( ngk_g ( : ) )

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'unformatted', status = 'replace')
    WRITE ( unit ) stitle, sdate, stime
    WRITE ( unit ) nspin, ngm_g, ntran, cell_symmetry, nat, ecutrho, &
      nkstot / nspin, nsp, nkb, nhm, npwx_g, ecutwfc
    IF ( wfng_kgrid ) THEN
      WRITE ( unit ) nr1, nr2, nr3, wfng_nk1, wfng_nk2, wfng_nk3, &
        wfng_dk1, wfng_dk2, wfng_dk3
    ELSE
      WRITE ( unit ) nr1, nr2, nr3, nk1, nk2, nk3, &
        0.5D0 * dble ( k1 ), 0.5D0 * dble ( k2 ), 0.5D0 * dble ( k3 )
    ENDIF
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nkstot / nspin )
    WRITE ( unit ) ( wk ( ik ) * dble ( nspin ) / 2.0D0, ik = 1, nkstot / nspin )
    WRITE ( unit ) ( ( xk ( id, ik ), id = 1, nd ), ik = 1, nkstot / nspin )
    WRITE ( unit ) ( ityp ( iat ), iat = 1, nat )
    WRITE ( unit ) ( nh ( isp ), isp = 1, nsp )
    WRITE ( unit ) ( ( ( ( deeq ( jh, ih, iat, is ), &
      jh = 1, nhm ), ih = 1, nhm ), iat = 1, nat ), is = 1, nspin )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ngm_g
    WRITE ( unit ) ( ( gvec ( id, ig ), id = 1, nd ), ig = 1, ngm_g )
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  ALLOCATE ( igwk ( npwx_g ) )

  DO i = 1, nkstot

    ik = kmap ( i )
    is = smap ( i )

    igwk = 0

    ALLOCATE ( itmp ( npw_g ) )
    itmp = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      DO ig = 1, ngk ( ik - iks + 1 )
        itmp ( igk_l2g ( ig, ik - iks + 1 ) ) = igk_l2g ( ig, ik - iks + 1 )
      ENDDO
    ENDIF
    CALL mp_sum ( itmp )
    ngg = 0
    DO ig = 1, npw_g
      IF ( itmp ( ig ) .EQ. ig ) THEN
        ngg = ngg + 1
        igwk ( ngg ) = ig
      ENDIF
    ENDDO
    DEALLOCATE ( itmp )

    IF ( ionode ) THEN
      IF ( is .EQ. 1 ) THEN
        WRITE ( unit ) nrecord
        WRITE ( unit ) ngk_g ( ik )
        WRITE ( unit ) ( ( gvec ( j, igwk ( ig ) ), j = 1, 3 ), &
          ig = 1, ngk_g ( ik ) )
      ENDIF
    ENDIF

    local_pw = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      ALLOCATE ( itmp ( npwx ) )
      npw = npwx
      CALL gk_sort ( xk ( 1, ik ), ngm, g, ecutwfc / tpiba2, &
        npw, itmp ( 1 ), g2kin )
      CALL init_us_2 ( npw, itmp, xk ( 1, ik ), vkb )
      local_pw = ngk ( ik - iks + 1 )
      DEALLOCATE ( itmp )
    ENDIF

    ALLOCATE ( igwf_l2g ( local_pw ) )
    igwf_l2g = 0
    DO ig = 1, local_pw
      ngg = igk_l2g ( ig, ik - iks + 1 )
      DO j = 1, ngk_g ( ik )
        IF ( ngg .EQ. igwk ( j ) ) THEN
          igwf_l2g ( ig ) = j
          EXIT
        ENDIF
      ENDDO
    ENDDO

    ALLOCATE ( ipmask ( nproc ) )
    ipmask = 0
    ipsour = ionode_id
    IF ( npool .GT. 1 ) THEN
      IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
        IF ( me_pool .EQ. root_pool ) ipmask ( mpime + 1 ) = 1
      ENDIF
      CALL mp_sum ( ipmask )
      DO j = 1, nproc
        IF ( ipmask ( j ) .EQ. 1 ) ipsour = j - 1
      ENDDO
    ENDIF
    DEALLOCATE ( ipmask )

    igwx = 0
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) &
      igwx = MAXVAL ( igwf_l2g ( 1 : local_pw ) )
    CALL mp_max ( igwx, intra_pool_comm )
    IF ( ipsour .NE. ionode_id ) &
      CALL mp_get ( igwx, igwx, mpime, ionode_id, ipsour, 1 )
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. igwx .NE. ngk_g ( ik ) ) &
      ierr = 1
    CALL mp_max ( ierr )
    IF ( ierr .GT. 0 ) &
      CALL errore ( 'write_vnlg', 'igwx ngk_g', ierr )

    ALLOCATE ( vkb_g ( MAX ( 1, igwx ) ) )

    DO ikb = 1, nkb

      vkb_g = ( 0.0D0, 0.0D0 )
      IF ( npool .GT. 1 ) THEN
        IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
          CALL mergewf ( vkb ( :, ikb ), vkb_g, local_pw, igwf_l2g, &
            me_pool, nproc_pool, root_pool, intra_pool_comm )
        ENDIF
        IF ( ipsour .NE. ionode_id ) THEN
          CALL mp_get ( vkb_g, vkb_g, mpime, ionode_id, ipsour, ikb, &
            world_comm )
        ENDIF
      ELSE
        CALL mergewf ( vkb ( :, ikb ), vkb_g, local_pw, igwf_l2g, &
          mpime, nproc, ionode_id, world_comm )
      ENDIF

      IF ( ionode ) THEN
        WRITE ( unit ) nrecord
        WRITE ( unit ) igwx
        WRITE ( unit ) ( vkb_g ( ig ), ig = 1, igwx )
      ENDIF

    ENDDO

    DEALLOCATE ( vkb_g )
    DEALLOCATE ( igwf_l2g )

  ENDDO

  IF ( ionode ) THEN
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( igwk )
  DEALLOCATE ( igk_l2g )
  DEALLOCATE ( ngk_g )
  DEALLOCATE ( gvec )
  DEALLOCATE ( smap )
  DEALLOCATE ( kmap )

  RETURN

END SUBROUTINE write_vnlg

!-------------------------------------------------------------------------------

SUBROUTINE check_dfpt_modes (imode, dfpt_type, nmode)

  USE control_ph, ONLY : tmp_dir_ph
  USE grid_dimensions, ONLY : nr1x, nr2x, nr3x
  USE io_files, ONLY : diropn
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_bcast
  USE units_ph, ONLY : iudvscf, lrdrho

  IMPLICIT NONE

  integer, intent(in)  :: imode
  integer, intent(in)  :: dfpt_type
  integer, intent(out) :: nmode

  integer :: statb(13)
  logical :: exst
  integer*8 :: unf_recl
  character*10 :: filename

  iudvscf = 27
  lrdrho = 2 * nr1x * nr2x * nr3x * nspin ! from PH/openfilq.f90

  IF(ionode) THEN
    
    IF(dfpt_type == 1) THEN
      filename = "dvscf"
    ELSE
      filename = "dvscf.E"
    ENDIF

    CALL diropn (iudvscf, filename, lrdrho, exst, tmp_dir_ = tmp_dir_ph )
    IF(.NOT. exst) CALL errore ('check_dfpt_modes', 'could not open file', 7)

    ! from io_files.f90, diropn
#if defined(__SX6)
#  define DIRECT_IO_FACTOR 1
#else
#  define DIRECT_IO_FACTOR 8 
#endif

    ! NOTE: ifort needs -assume byterecl or else nrec will be 4 times too large
    unf_recl = DIRECT_IO_FACTOR * int(lrdrho, kind=kind(unf_recl))
    CALL fstat ( iudvscf, statb) ! see EPW, readdvscf.f90
    nmode = statb(8) / unf_recl
    IF(nmode * unf_recl /= statb(8) .and. ionode) &
      WRITE(0,*) 'WARNING: dvscf file does not contain an integral number of records.', &
      statb(8) * 1d0 / unf_recl
    IF(ionode) WRITE(6,*) 'Number of modes in dvscf file: ',  nmode
    IF(ionode) WRITE(6,*) 'Starting with mode number: ',  imode
    IF (imode > nmode) CALL errore('pw2bgw', 'starting mode is not present in dvscf file', iudvscf)

    IF(dfpt_type == 1) then
      IF(nmode < 3 * nat) WRITE(0,'(a)') 'WARNING: dvscf file does not contain all modes.'
      IF(nmode > 3 * nat) WRITE(0,'(a)') 'WARNING: dvscf file contains too many modes.'
    ENDIF
    IF(dfpt_type == 2) then
      IF(nmode < 3) WRITE(0,'(a)') 'WARNING: dvscf_e file does not contain all polarizations.'
      IF(nmode > 3) WRITE(0,'(a)') 'WARNING: dvscf_e file contains too many polarizations.'
    ENDIF
  ENDIF

  CALL mp_bcast ( nmode, ionode_id )

  RETURN

END SUBROUTINE check_dfpt_modes

!-------------------------------------------------------------------------------

subroutine check_inversion(real_or_complex, ntran, mtrx, nspin, warn, real_need_inv)

! check_inversion    Originally By D. Strubbe    Last Modified 10/14/2010
! Check whether our choice of real/complex version is appropriate given the
! presence or absence of inversion symmetry.

  USE io_global, ONLY : ionode

  implicit none

  integer, intent(in) :: real_or_complex
  integer, intent(in) :: ntran
  integer, intent(in) :: mtrx(3, 3, 48)
  integer, intent(in) :: nspin
  logical, intent(in) :: warn ! set to false to suppress warnings, for converters
  logical, intent(in) :: real_need_inv ! use for generating routines to block real without inversion

  integer :: invflag, isym, ii, jj, itest

  invflag = 0
  do isym = 1, ntran
    itest = 0
    do ii = 1, 3
      do jj = 1, 3
        if(ii .eq. jj) then
          itest = itest + (mtrx(ii, jj, isym) + 1)**2
        else
          itest = itest + mtrx(ii, jj, isym)**2
        endif
      enddo
    enddo
    if(itest .eq. 0) invflag = invflag + 1
    if(invflag .gt. 1) call errore('check_inversion', 'More than one inversion symmetry operation is present.', invflag)
  enddo

  if(real_or_complex .eq. 2) then
    if(invflag .ne. 0 .and. warn .and. nspin == 1) then
      if(ionode) write(0, '(a)') 'WARNING: Inversion symmetry is present. The real version would be faster.'
    endif
  else
    if(invflag .eq. 0) then
      if(real_need_inv) then
        call errore('check_inversion', 'The real version cannot be used without inversion symmetry.', -1)
      endif
      if(ionode) then
        write(0, '(a)') 'WARNING: Inversion symmetry is absent in symmetries used to reduce k-grid.'
        write(0, '(a)') 'Be sure inversion is still a spatial symmetry, or you must use complex version instead.'
      endif
    endif
    if(nspin > 1) then
      call errore('check_inversion', &
        'Real version may only be used for spin-unpolarized calculations.', nspin)
    endif
  endif

  return

end subroutine check_inversion

!-------------------------------------------------------------------------------
! Set up some things for DFPT electron-phonon matrix elements.
subroutine initialize_phonon()

  USE atom, ONLY : msh, rgrid
  USE becmod, ONLY : calbec, allocate_bec_type
  USE cell_base, ONLY : omega, tpiba, tpiba2
  USE control_ph, ONLY : current_iq, trans, tmp_dir_ph
  USE eqv, ONLY : dpsi, dvpsi, vlocq
  USE gvect, ONLY : g, ngm
  USE io_files, ONLY : tmp_dir
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, ntyp => nsp
  USE kinds, ONLY : DP
  USE klist, ONLY : nks, xk, nkstot
  USE modes, ONLY : u, npert, name_rap_mode
  USE noncollin_module, ONLY : npol
  USE ph_restart, ONLY : ph_readfile
  USE phus, ONLY : alphap, becp1
  USE qpoint, ONLY : npwq, igkq, xq, eigqts, ikks, ikqs, nksq  ! phcom
  USE uspp, ONLY : nkb, vkb
  USE uspp_param, ONLY : upf
  USE wavefunctions_module, ONLY : evc
  USE wvfct, ONLY : igk, nbnd, npw, npwx, g2kin, ecutwfc

  IMPLICIT NONE

  integer :: nt, ik, ipol, ibnd, ig, ierr, ii, jj, iks, ike
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:)

  if (ionode) then
    WRITE ( 6, '(a)') 'Performing initialization for ionic perturbation.'
  endif

  tmp_dir_ph = trim(tmp_dir) // '_ph0/'

  if ( ANY (upf(1:ntyp)%nlcc) ) then
    CALL errore ( 'initialize_phonon', 'NLCC not implemented.', ierr )
  endif

  ! if gamma_gamma tricks are used, the dvscf file and wfns do not mean what we think they do...
  if(nkstot == 1 .and. ionode) write(0,'(a)') &
    "WARNING: Be sure you are not using a phonon run with gamma_gamma tricks. Must be nogg = .true."

! allocate_phq
  allocate (u ( 3 * nat, 3 * nat))
  allocate (npert ( 3 * nat))
  allocate (name_rap_mode( 3 * nat))

  current_iq = 1
  trans = .true.
!  lgamma = .true.  ! using this causes crash from iotk in ph_readfile??
  CALL ph_readfile ('data_u', ierr) ! reads from XML, for uact
  IF ( ierr /= 0 ) CALL errore ( 'initialize_phonon', 'failed to open phonon XML', ierr )
  if(ionode) then
    OPEN (unit = 16, file = 'displacements.dat', form = 'formatted', status = 'replace')
    write(6,'(a)') 'displacements read from XML file'
    do ii = 1, 3 * nat
      do jj = 1, 3 * nat
        ! we must use (Re, Im) format for this to recognized as complex by a Fortran read
        write(16,'(a,f12.7,a,f12.7,a)',advance='no') '  (', dble(u(jj, ii)), ',', aimag(u(jj, ii)), ')'
      enddo
      write(16,*)
    enddo
    CLOSE(16)
  endif

! allocate_phq
  allocate (dvpsi ( npwx*npol , nbnd))

!  if (lgamma) then
  igkq => igk

! phq_init
!     IF ( lgamma ) THEN
  npwq = npw

  allocate(eigqts(nat))
  eigqts(:) = CMPLX(1d0, 0d0)

  xq(1:3) = 0d0

! allocate_phq.f90

  if(ionode) write(6,'(a,i6)') 'Number of species: ', ntyp
  allocate (vlocq ( ngm , ntyp))

! phq_init.f90:

  ! the fourier components of the local potential at q+G 
  vlocq(:,:) = 0.D0 
  ! 
  DO nt = 1, ntyp
     ! 
     IF (upf(nt)%tcoulombp) then
        CALL setlocq_coul ( xq, upf(nt)%zp, tpiba2, ngm, g, omega, vlocq(1,nt) ) 
     ELSE 
        CALL setlocq( xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r,& 
                   upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, g, omega, & 
                   vlocq(1,nt) ) 
     END IF 
    ! 
!     write(6,*) 'init setlocq nt = ', nt
!     write(6,*) vlocq(:, nt)
  END DO

! deallocate_phq.f90:53:  if(allocated(vlocq)) deallocate (vlocq)

! initialize_ph
!  IF ( lgamma ) THEN
  nksq = nks
  ALLOCATE(ikks(nksq), ikqs(nksq))
  DO ik=1,nksq
    ikks(ik) = ik
    ikqs(ik) = ik
  ENDDO

  ALLOCATE (becp1(nksq))
  ALLOCATE (alphap(3,nksq))

  DO ik=1,nksq
     call allocate_bec_type ( nkb, nbnd, becp1(ik) )
     DO ipol=1,3
        call allocate_bec_type ( nkb, nbnd, alphap(ipol,ik) )
     ENDDO
  END DO
!  CALL allocate_bec_type ( nkb, nbnd, becp )

! phq_init
     ! ... e) we compute the becp terms which are used in the rest of 
     ! ...    the code 
     ! 
  ALLOCATE( aux1( npwx*npol, nbnd ) )

  CALL k_pools(iks, ike)

  do ik = 1, nksq

    CALL gk_sort (xk (1, ik - iks + 1), ngm, g, ecutwfc &
      / tpiba2, npw, igk, g2kin)
    CALL init_us_2 ( npw, igk, xk ( 1, ik - iks + 1), vkb )

!      npw = npwx

!    CALL davcio (evc, nwordwfc, iunwfc, ik - iks + 1, -1) ?
     CALL calbec (npw, vkb, evc, becp1(ik) )
     !
     ! ... e') we compute the derivative of the becp term with respect to an 
     !         atomic displacement 
     ! 
     DO ipol = 1, 3 
        aux1=(0.d0,0.d0) 
        DO ibnd = 1, nbnd 
           DO ig = 1, npw 
              aux1(ig,ibnd) = evc(ig,ibnd) * tpiba * ( 0.D0, 1.D0 ) * & 
                              ( xk(ipol,ik) + g(ipol,igk(ig)) ) 
           END DO 
!           IF (noncolin) THEN 
!              DO ig = 1, npw 
!                 aux1(ig+npwx,ibnd)=evc(ig+npwx,ibnd)*tpiba*(0.D0,1.D0)*& 
!                           ( xk(ipol,ikk) + g(ipol,igk(ig)) ) 
!              END DO 
!           END IF 
        END DO 
        CALL calbec (npw, vkb, aux1, alphap(ipol,ik) ) 
     END DO
   enddo

   deallocate(aux1)

   return
end subroutine initialize_phonon

!-------------------------------------------------------------------------------
! Set up some things for DFPT electron-electric matrix elements.
subroutine initialize_electric()

  USE atom, ONLY : msh, rgrid
  USE becmod, ONLY : calbec, allocate_bec_type
  USE cell_base, ONLY : omega, tpiba, tpiba2
  USE control_ph, ONLY : current_iq, trans, tmp_dir_ph
  USE eqv, ONLY : dpsi, dvpsi, vlocq
  USE gvect, ONLY : g, ngm
  USE io_files, ONLY : tmp_dir
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, ntyp => nsp
  USE kinds, ONLY : DP
  USE klist, ONLY : nks, xk, nkstot
  USE modes, ONLY : u, npert, name_rap_mode
  USE noncollin_module, ONLY : npol
  USE ph_restart, ONLY : ph_readfile
  USE phus, ONLY : alphap, becp1
  USE qpoint, ONLY : npwq, igkq, xq, eigqts, ikks, ikqs, nksq  ! phcom
  USE uspp, ONLY : nkb, vkb
  USE uspp_param, ONLY : upf
  USE wavefunctions_module, ONLY : evc
  USE wvfct, ONLY : igk, nbnd, npw, npwx, g2kin, ecutwfc

  IMPLICIT NONE

  if (ionode) then
    WRITE ( 6, '(a)') 'Performing initialization for electric perturbation.'
  endif

  tmp_dir_ph = trim(tmp_dir) // '_ph0/'

  ! if gamma_gamma tricks are used, the dvscf file and wfns do not mean what we think they do...
  if(nkstot == 1 .and. ionode) write(0,'(a)') &
    "WARNING: Be sure you are not using a phonon run with gamma_gamma tricks. Must be nogg = .true."

  allocate ( dpsi ( npwx*npol , nbnd))

  ! more is needed!

  return

end subroutine initialize_electric

!-------------------------------------------------------------------------------

subroutine write_header(unit, stitle)

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, &
    ibrav, symm_type
  USE constants, ONLY : pi, tpi, eps6
  USE grid_dimensions, ONLY : nr1, nr2, nr3
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, atm, ityp, tau 
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_global, ONLY : intra_pool_comm
  USE scf, ONLY : rho
  USE symm_base, ONLY : s, ftau, nsym

  IMPLICIT NONE

  integer, intent(in) :: unit
  character(len=32), intent(in) :: stitle

  character :: cdate*9, ctime*9, sdate*32, stime*32
  integer :: id, ig, i, j, k, ierr
  integer :: nd, ns, nst, nsf, ng_l, ng_g
  integer :: ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, dr1, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: g_g ( :, : )

  INTEGER, EXTERNAL :: atomic_number

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)

  nrecord = 1
  nd = 3

  ng_l = ngm
  ng_g = ngm_g
  call set_spin(ns,nst,nsf,nspin)

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( ibrav .GE. 1 .AND. ibrav .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( ibrav .GE. 4 .AND. ibrav .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( ibrav .GE. 6 .AND. ibrav .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_header', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2, dr1 )
    t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( nr1 )
    t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( nr2 )
    t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( nr3 )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ALLOCATE ( g_g ( nd, ng_g ) )

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO

  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO

  CALL mp_sum ( g_g, intra_pool_comm )

  IF(ionode) then
    WRITE ( unit ) stitle, sdate, stime
    WRITE ( unit ) nspin, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) nr1, nr2, nr3
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
  ENDIF

  DEALLOCATE ( g_g )
  RETURN
  
end subroutine write_header

!-------------------------------------------------------------------------------

! Alternate method of calculating DFPT matrix elements, for testing.
! Project variation of occupied wavefunctions (dwf files) onto wavefunctions.
! Must hack ph code to preserve dwf's. Must be done one mode at a time.
! Only gives occ-unocc mtxels, which are actually not the ones needed in forces.
subroutine dfpt_dwf()

  USE cell_base, ONLY : tpiba2
  USE constants, ONLY : rytoev
  USE control_ph, ONLY : tmp_dir_ph
  USE eqv, ONLY : dpsi
  USE gvect, ONLY : g, ngm
  USE io_files, ONLY : diropn, prefix, nwordwfc, iunwfc, nd_nmbr
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nelec
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_global, ONLY : intra_pool_comm
  USE noncollin_module, ONLY : npol
  USE qpoint, ONLY : nksq
  USE units_ph, ONLY : iudwf, lrdwf
  USE wavefunctions_module, ONLY : evc
  USE wvfct, ONLY : nbnd, npw, npwx, igk, g2kin, ecutwfc, et

  IMPLICIT NONE

  integer :: nrec, ipert, ik, ib, ib2, nbnd_occ, iks, ike
  logical :: exst
  complex(DP), allocatable :: mtxel(:,:)

  nbnd_occ = nelec
  if(nspin == 1) nbnd_occ = nbnd_occ / 2d0
  allocate (dpsi (npwx*npol, nbnd_occ / nspin))
  allocate (mtxel (nbnd_occ, nbnd))

  ! from openfilq.f90
  iudwf = 22
  ! we do not actually know the number of bands in scf/ph. we guess...
  lrdwf = 2 * nbnd_occ * npwx * npol
  CALL diropn (iudwf, 'dwf', lrdwf, exst, tmp_dir_ = tmp_dir_ph)
  IF (.NOT.exst) CALL errore ('dfpt_dwf','file '//trim(prefix)//'.dwf not found', 1)

  ! this would be the index of the perturbation within the irrep.
  ! we do not have this info though, so we will just deal with the first one
  ipert = 1

  CALL k_pools(iks, ike)
  do ik = iks, ike
    nrec = (ipert - 1) * nksq + ik
    CALL gk_sort (xk (1, ik - iks + 1), ngm, g, ecutwfc &
      / tpiba2, npw, igk, g2kin)
    CALL davcio (evc, nwordwfc, iunwfc, ik - iks + 1, -1)
    call davcio (dpsi, lrdwf, iudwf, nrec, -1)
!    write(37,*) evc
!    write(38,*) dpsi
    do ib = 1, nbnd_occ
      do ib2 = 1, nbnd
        ! zgemm?
!        mtxel(ib, ib2) = sum(conjg(evc(1:npw, ib)) * evc(1:npw, ib2))
        mtxel(ib, ib2) = sum(conjg(dpsi(1:npw, ib)) * evc(1:npw, ib2))
        CALL mp_sum (mtxel(ib, ib2), intra_pool_comm)
        write(6,'(2i6,2f12.6)') ib, ib2, mtxel(ib, ib2) * (et(ib2, ik) - et(ib, ik)) * rytoev
      enddo
    enddo
  enddo

  deallocate(dpsi)
  deallocate(mtxel)

  RETURN
  
end subroutine dfpt_dwf

!-------------------------------------------------------------------------------

! determine start and end k-points for this pool
! taken from PW/pw_restart.f90
subroutine k_pools(iks, ike)

  USE klist, ONLY : nkstot
  USE mp_global, ONLY : kunit, my_pool_id, npool

  IMPLICIT NONE

  integer, intent(out) :: iks, ike

  integer :: nkbl, nkl, nkr

! ... find out the number of pools
!  npool = nproc / nproc_pool
! apparently this is already initialized. --DAS

! ... find out number of k-points blocks
  nkbl = nkstot / kunit
! ... k-points per pool  
  nkl = kunit * ( nkbl / npool )
! ... find out the remainder
  nkr = ( nkstot - nkl * npool ) / kunit
! ... Assign the remainder to the first nkr pools
  IF ( my_pool_id .LT. nkr ) nkl = nkl + kunit
! ... find out the index of the first k-point in this pool
  iks = nkl * my_pool_id + 1
  IF ( my_pool_id .GE. nkr ) iks = iks + nkr * kunit
! ... find out the index of the last k-point in this pool
  ike = iks + nkl - 1

  RETURN

end subroutine k_pools

!-------------------------------------------------------------------------------

! Set spin-related parameters ns, nst, nsf
!
subroutine set_spin(ns, nst, nsf, nspin)

  IMPLICIT NONE

  integer, intent(out) :: ns, nst, nsf
  integer, intent(in) :: nspin

  IF ( nspin == 4 ) THEN
    ns = 1
    nst = 2
    nsf = 4
  ELSE
    ns = nspin
    nst = nspin
    nsf = nspin
  ENDIF

  RETURN

end subroutine set_spin

!-------------------------------------------------------------------------------

end module pw2bgw_m

PROGRAM pw2bgw

  USE constants, ONLY : eps12
  USE control_flags, ONLY : gamma_only
  USE environment, ONLY : environment_start, environment_end
  USE io_files, ONLY : prefix, tmp_dir, outdir
  USE io_global, ONLY : ionode, ionode_id
  USE kinds, ONLY : DP
  USE klist, ONLY : nkstot
  USE mp, ONLY : mp_bcast
  USE mp_global, ONLY : mp_startup
  USE paw_variables, ONLY : okpaw
  use scf, ONLY: rho_core, rhog_core
  USE uspp, ONLY : okvan
  USE wvfct, ONLY : nbnd

  USE pw2bgw_m

  IMPLICIT NONE

  character(len=6) :: codename = 'PW2BGW'

  integer :: real_or_complex
  logical :: wfng_flag
  character ( len = 256 ) :: wfng_file
  logical :: wfng_kgrid
  integer :: wfng_nk1
  integer :: wfng_nk2
  integer :: wfng_nk3
  real (DP) :: wfng_dk1
  real (DP) :: wfng_dk2
  real (DP) :: wfng_dk3
  logical :: wfng_occupation
  integer :: wfng_nvmin
  integer :: wfng_nvmax
  logical :: rhog_flag
  character ( len = 256 ) :: rhog_file
  logical :: vxcg_flag
  character ( len = 256 ) :: vxcg_file
  logical :: vxc0_flag
  character ( len = 256 ) :: vxc0_file
  logical :: vxc_flag
  character ( len = 256 ) :: vxc_file
  integer :: dfpt_type  ! 0 = none, 1 = ionic, 2 = electric
  logical :: dfpt_from_dwf
  character ( len = 256 ) :: dfpt_file
  character :: vxc_integral
  integer :: vxc_diag_nmin
  integer :: vxc_diag_nmax
  integer :: vxc_offdiag_nmin
  integer :: vxc_offdiag_nmax
  integer :: dfpt_mode_start
  integer :: dfpt_mode_end
  logical :: vxc_zero_rho_core
  character ( len = 20 ) :: input_dft
  logical :: exx_flag
  logical :: vnlg_flag
  character ( len = 256 ) :: vnlg_file

  NAMELIST / input_pw2bgw / prefix, outdir, &
    real_or_complex, wfng_flag, wfng_file, wfng_kgrid, &
    wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3, &
    wfng_occupation, wfng_nvmin, wfng_nvmax, rhog_flag, rhog_file, &
    vxcg_flag, vxcg_file, vxc0_flag, vxc0_file, vxc_flag, &
    vxc_file, vxc_integral, vxc_diag_nmin, vxc_diag_nmax, &
    vxc_offdiag_nmin, vxc_offdiag_nmax, vxc_zero_rho_core, &
!#BEGIN_INTERNAL_ONLY
    dfpt_type, dfpt_from_dwf, dfpt_file, dfpt_mode_start, dfpt_mode_end, &
!#END_INTERNAL_ONLY
    input_dft, exx_flag, vnlg_flag, vnlg_file

  integer :: ii, ios, imode, nmode
  character ( len = 256 ) :: output_file_name, tmpstr
  logical :: die_spinors

  character (len=256), external :: trimcheck
  character (len=1), external :: lowercase

#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( codename )

  ! assign defaults
  prefix = 'prefix'
  CALL get_env ( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM ( outdir ) == ' ' ) outdir = './'
  real_or_complex = 2
  wfng_flag = .FALSE.
  wfng_file = 'WFN'
  wfng_kgrid = .FALSE.
  wfng_nk1 = 0
  wfng_nk2 = 0
  wfng_nk3 = 0
  wfng_dk1 = 0.0D0
  wfng_dk2 = 0.0D0
  wfng_dk3 = 0.0D0
  wfng_occupation = .FALSE.
  wfng_nvmin = 0
  wfng_nvmax = 0
  rhog_flag = .FALSE.
  rhog_file = 'RHO'
  vxcg_flag = .FALSE.
  vxcg_file = 'VXC'
  vxc0_flag = .FALSE.
  vxc0_file = 'vxc0.dat'
  vxc_flag = .FALSE.
  vxc_file = 'vxc.dat'
  vxc_integral = 'g'
  vxc_diag_nmin = 0
  vxc_diag_nmax = 0
  vxc_offdiag_nmin = 0
  vxc_offdiag_nmax = 0
  vxc_zero_rho_core = .TRUE.
  input_dft = 'sla+pz'
  exx_flag = .FALSE.
  vnlg_flag = .FALSE.
  vnlg_file = 'VNL'
  dfpt_type = 0
  dfpt_from_dwf = .FALSE.
  dfpt_file = 'dfpt.dat'
  dfpt_mode_start = 1
  dfpt_mode_end = 0 ! do all

  IF ( ionode ) THEN
    CALL input_from_file ( )
    READ ( 5, input_pw2bgw, iostat = ios )
    IF ( ios /= 0 ) CALL errore ( codename, 'input_pw2bgw', abs ( ios ) )

    DO ii = 1, LEN_TRIM (vxc_integral)
      vxc_integral(ii:ii) = lowercase (vxc_integral(ii:ii))
    END DO

    IF ( real_or_complex /= 1 .AND. real_or_complex /= 2 ) &
      CALL errore ( codename, 'real_or_complex', 1 )
    IF ( vxc_integral /= 'r' .AND. vxc_integral /= 'g' ) &
      CALL errore ( codename, 'vxc_integral', 1 )
  ENDIF

  tmp_dir = trimcheck ( outdir )
  CALL mp_bcast ( outdir, ionode_id )
  CALL mp_bcast ( tmp_dir, ionode_id )
  CALL mp_bcast ( prefix, ionode_id )
  CALL mp_bcast ( real_or_complex, ionode_id )
  CALL mp_bcast ( wfng_flag, ionode_id )
  CALL mp_bcast ( wfng_file, ionode_id )
  CALL mp_bcast ( wfng_kgrid, ionode_id )
  CALL mp_bcast ( wfng_nk1, ionode_id )
  CALL mp_bcast ( wfng_nk2, ionode_id )
  CALL mp_bcast ( wfng_nk3, ionode_id )
  CALL mp_bcast ( wfng_dk1, ionode_id )
  CALL mp_bcast ( wfng_dk2, ionode_id )
  CALL mp_bcast ( wfng_dk3, ionode_id )
  CALL mp_bcast ( wfng_occupation, ionode_id )
  CALL mp_bcast ( wfng_nvmin, ionode_id )
  CALL mp_bcast ( wfng_nvmax, ionode_id )
  CALL mp_bcast ( rhog_flag, ionode_id )
  CALL mp_bcast ( rhog_file, ionode_id )
  CALL mp_bcast ( vxcg_flag, ionode_id )
  CALL mp_bcast ( vxcg_file, ionode_id )
  CALL mp_bcast ( vxc0_flag, ionode_id )
  CALL mp_bcast ( vxc0_file, ionode_id )
  CALL mp_bcast ( vxc_flag, ionode_id )
  CALL mp_bcast ( vxc_integral, ionode_id )
  CALL mp_bcast ( vxc_file, ionode_id )
  CALL mp_bcast ( vxc_diag_nmin, ionode_id )
  CALL mp_bcast ( vxc_diag_nmax, ionode_id )
  CALL mp_bcast ( vxc_offdiag_nmin, ionode_id )
  CALL mp_bcast ( vxc_offdiag_nmax, ionode_id )
  CALL mp_bcast ( vxc_zero_rho_core, ionode_id )
  CALL mp_bcast ( input_dft, ionode_id )
  CALL mp_bcast ( exx_flag, ionode_id )
  CALL mp_bcast ( vnlg_flag, ionode_id )
  CALL mp_bcast ( vnlg_file, ionode_id )
  CALL mp_bcast ( dfpt_type, ionode_id )
  CALL mp_bcast ( dfpt_from_dwf, ionode_id )
  CALL mp_bcast ( dfpt_file, ionode_id )
  CALL mp_bcast ( dfpt_mode_start, ionode_id )
  CALL mp_bcast ( dfpt_mode_end, ionode_id )

  if(dfpt_type < 0 .or. dfpt_type > 2) &
    call errore ( 'pw2bgw', 'unknown dfpt_type', dfpt_type)
  if(dfpt_type /= 0 .and. dfpt_mode_start < 1) &
    call errore ( 'pw2bgw', 'dfpt_mode_start < 1', dfpt_mode_start)

  if (ionode) then
    WRITE ( 6, '(a)') 'Reading basic information from files.'
  endif
  CALL read_file ( )

  if (ionode) then
    if (MAX (MAXVAL (ABS (rho_core (:) ) ), MAXVAL (ABS (rhog_core (:) ) ) ) &
      .LT. eps12) then
      WRITE ( 6, '(/,5x,"NLCC is absent")' )
    else
      WRITE ( 6, '(/,5x,"NLCC is present")' )
    endif
  endif
  if (okvan) call errore ( 'pw2bgw', 'BGW cannot use USPP.', 3 )
  if (okpaw) call errore ( 'pw2bgw', 'BGW cannot use PAW.', 4 )
  if (gamma_only) call errore ( 'pw2bgw', 'BGW cannot use gamma-only run.', 5 )
  die_spinors = (nspin == 4)
  die_spinors = .false.
  if (die_spinors) call errore ( 'pw2bgw', 'BGW cannot use spinors.', 6 )
  if (real_or_complex == 1 .AND. vxc_flag .AND. vxc_offdiag_nmax > 0) &
    call errore ( 'pw2bgw', 'Off-diagonal matrix elements of Vxc ' // &
    'with real wavefunctions are not implemented, compute them in ' // &
    'Sigma using VXC.', 7)

  CALL openfil_pp ( )

  if ( ionode ) WRITE ( 6, '("")' )

  if ( ionode ) write ( 6, '(a, i8)' ) "Total number of bands    available: ", nbnd
  if ( ionode ) write ( 6, '(a, i8)' ) "Total number of k-points available: ", nkstot

  IF ( wfng_flag ) THEN
    output_file_name = TRIM ( outdir ) // '/' // TRIM ( wfng_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_wfng")' )
    CALL start_clock ( 'write_wfng' )
    CALL write_wfng ( output_file_name, real_or_complex, &
      wfng_kgrid, wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, &
      wfng_dk3, wfng_occupation, wfng_nvmin, wfng_nvmax )
    CALL stop_clock ( 'write_wfng' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_wfng",/)' )
  ENDIF

  IF ( rhog_flag ) THEN
    output_file_name = TRIM ( outdir ) // '/' // TRIM ( rhog_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_rhog")' )
    CALL start_clock ( 'write_rhog' )
    CALL write_rhog ( output_file_name, real_or_complex )
    CALL stop_clock ( 'write_rhog' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_rhog",/)' )
  ENDIF

  IF ( vxcg_flag ) THEN
    output_file_name = TRIM ( outdir ) // '/' // TRIM ( vxcg_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vxcg")' )
    CALL start_clock ( 'write_vxcg' )
    CALL write_vxcg ( output_file_name, real_or_complex, &
      vxc_zero_rho_core, input_dft, exx_flag )
    CALL stop_clock ( 'write_vxcg' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vxcg",/)' )
  ENDIF

  IF ( vxc0_flag ) THEN
    output_file_name = TRIM ( outdir ) // '/' // TRIM ( vxc0_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc0")' )
    CALL start_clock ( 'write_vxc0' )
    CALL write_vxc0 ( output_file_name, vxc_zero_rho_core, input_dft, &
      exx_flag )
    CALL stop_clock ( 'write_vxc0' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc0",/)' )
  ENDIF

  IF ( vxc_flag ) THEN
    output_file_name = TRIM ( outdir ) // '/' // TRIM ( vxc_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc")' )
    CALL start_clock ( 'write_vxc' )
    CALL write_vxc ( vxc_integral, output_file_name, &
      vxc_diag_nmin, vxc_diag_nmax, &
      vxc_offdiag_nmin, vxc_offdiag_nmax, &
      vxc_zero_rho_core, input_dft, exx_flag, 0, 0 )
    CALL stop_clock ( 'write_vxc' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc",/)' )
  ENDIF

  IF ( vnlg_flag ) THEN
    output_file_name = TRIM ( outdir ) // '/' // TRIM ( vnlg_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vnlg")' )
    CALL start_clock ( 'write_vnlg' )
    CALL write_vnlg ( output_file_name, wfng_kgrid, wfng_nk1, &
      wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3 )
    CALL stop_clock ( 'write_vnlg' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vnlg",/)' )
  ENDIF

  IF ( dfpt_type /= 0 ) THEN
    IF(dfpt_type == 1) THEN
      call initialize_phonon()
    ELSE
      call initialize_electric()
    ENDIF
    CALL check_dfpt_modes ( dfpt_mode_start, dfpt_type, nmode )
    IF(dfpt_mode_end == 0) dfpt_mode_end = nmode

    IF(dfpt_from_dwf) CALL dfpt_dwf()

    DO imode = dfpt_mode_start, dfpt_mode_end
      WRITE(tmpstr,'(i6.6)') imode
      IF ( ionode ) WRITE ( 6, '(5x,"call write_dfpt")' )
      CALL start_clock ( 'write_dfpt' )

!      output_file_name = TRIM ( outdir ) // '/' // TRIM ( dfpt_file ) // "_dvscf_mode" // TRIM ( tmpstr )
!      CALL write_vxc ( 'r', output_file_name, &
!        vxc_diag_nmin, vxc_diag_nmax, &
!        vxc_diag_nmin, vxc_diag_nmax, & ! we need the whole matrix here
!        input_dft, .false., dfpt_type, imode)

      output_file_name = TRIM ( outdir ) // '/' // TRIM ( dfpt_file ) // "_mode" // TRIM ( tmpstr )
      CALL write_vxc ( 'g', output_file_name, &
        vxc_diag_nmin, vxc_diag_nmax, &
        vxc_diag_nmin, vxc_diag_nmax, & ! we need the whole matrix here
        vxc_zero_rho_core, input_dft, .false., dfpt_type, imode)
      CALL stop_clock ( 'write_dfpt' )
      IF ( ionode ) WRITE ( 6, '(5x,"done write_dfpt",/)' )
    ENDDO
  ENDIF

  IF ( ionode ) WRITE ( 6, * )
  IF ( wfng_flag ) CALL print_clock ( 'write_wfng' )
  IF ( rhog_flag ) CALL print_clock ( 'write_rhog' )
  IF ( vxcg_flag ) CALL print_clock ( 'write_vxcg' )
  IF ( vxc0_flag ) CALL print_clock ( 'write_vxc0' )
  IF ( vxc_flag )  CALL print_clock ( 'write_vxc' )
  IF ( dfpt_type /= 0 ) CALL print_clock ( 'write_dfpt' )
  IF ( vnlg_flag ) CALL print_clock ( 'write_vnlg' )
  IF ( wfng_flag .AND. real_or_complex .EQ. 1 ) THEN
    IF ( ionode ) WRITE ( 6, '(/,5x,"Called by write_wfng:")' )
    CALL print_clock ( 'real_wfng' )
  ENDIF

  CALL environment_end ( codename )
  CALL stop_pp ( )

  STOP

END PROGRAM pw2bgw
