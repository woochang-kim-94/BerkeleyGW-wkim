!-------------------------------------------------------------------------------
!
!   bgw2pw.f90
!   $MPIRUN bgw2pw.x -in bgw2pw.in > bgw2pw.out
!   converts BerkeleyGW wavefunction file WFN to Quantum ESPRESSO format
!   based on espresso-4.3.2/PP/pw_export.f90
!   written by G. Samsonidze (February 2010)
!   The format of the input file is described in the sample file bgw2pw.inp.
!
!-------------------------------------------------------------------------------
!
!   Compilation: use install.sh, or manually follow these steps:
!
!   copy bgw2pw.f90 and make.bgw to espresso-4.3.2/PP
!   add this line to end of espresso-4.3.2/PP/Makefile: "include make.bgw"
!   make as follows : cd espresso-4.3.2 ; make pw ; cd PP ; make bgw2pw.x
!
!-------------------------------------------------------------------------------
!
!   subroutines
!
!   write_evc - generates eigenvalues and wavefunctions in espresso format
!   write_cd  - generates charge density in espresso format
!   check_inversion - checks whether real/complex version is appropriate
!
!-------------------------------------------------------------------------------
!
!   bgw2pw reads common parameters from file prefix.save/data-file.xml and
!   writes files prefix.save/charge-density.dat (charge density in R-space),
!   prefix.save/gvectors.dat (G-vectors for charge density and potential),
!   prefix.save/K$n/eigenval.xml (eigenvalues and occupations for nth k-point),
!   prefix.save/K$n/evc.dat (wavefunctions in G-space for nth k-point), and
!   prefix.save/K$n/gkvectors.dat (G-vectors for nth k-point)
!
!   bgw2pw doesn't modify file prefix.save/data-file.xml so make changes
!   to this file manually (for example you may change the number of bands)
!   copy data-file.xml from scf-dir to nscf-dir rather than linking to it
!   if you make a symbolic link from nscf-dir to scf-dir your original file
!   in scf-dir may be overwritten
!
!-------------------------------------------------------------------------------

#include "f_defs.h"

module bgw2pw_m

  implicit none

  private

  public :: &
    write_evc, &
    write_cd, &
    check_inversion

contains

!-------------------------------------------------------------------------------

SUBROUTINE write_evc ( input_file_name, real_or_complex, &
  wfng_nband, output_dir_name )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg
  USE constants, ONLY : eps6
  USE grid_dimensions, ONLY : nr1, nr2, nr3
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, g
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat
  USE iotk_module, ONLY : iotk_attlenx, iotk_free_unit, iotk_open_write, &
    iotk_write_begin, iotk_write_attr, iotk_write_empty, iotk_write_dat, &
    iotk_write_end, iotk_close_write, iotk_index
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nks, nkstot
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_bcast, mp_sum, mp_max, mp_barrier
  USE mp_global, ONLY : nproc, kunit, npool, my_pool_id, world_comm, &
    intra_pool_comm
  USE symm_base, ONLY : s, nsym
  USE xml_io_base, ONLY : create_directory, kpoint_dir, wfc_filename
#ifdef __PARA
  USE parallel_include, ONLY : MPI_INTEGER, MPI_DOUBLE_COMPLEX
#endif
  USE wvfct, ONLY : npwx, g2kin, ecutwfc

  IMPLICIT NONE

  character ( len = 256 ), intent ( in ) :: input_file_name
  integer, intent ( in ) :: real_or_complex
  integer, intent ( in ) :: wfng_nband
  character ( len = 256 ), intent ( in ) :: output_dir_name

  logical :: f1, f2
  integer :: ierr, i, j, iu, ik, is, ib, ig, jg, fg, ir, &
    na, nk, ns, nb, nbgw, ng, ngkmax, ntran, cell_symmetry, &
    nkbl, nkl, nkr, iks, ike, npw, npw_g, ngkdist_l, ngkdist_g, &
    igk_l2g, irecord, nrecord, ng_irecord, nr ( 3 )
  real ( DP ) :: ecutrho, ecutwfn, celvol, recvol, al, bl, xdel, &
    a ( 3, 3 ), b ( 3, 3 ), adot ( 3, 3 ), bdot ( 3, 3 )
  character :: sdate*32, stime*32, stitle*32
  character ( len = 256 ) :: filename
  character ( iotk_attlenx ) :: attr

  integer, allocatable :: itmp ( : )
  integer, allocatable :: ngk ( : )
  integer, allocatable :: gvec ( :, : )
  integer, allocatable :: igk_buf ( : )
  integer, allocatable :: igk_dist ( :, : )
  integer, allocatable :: gk_buf ( :, : )
  integer, allocatable :: gk_dist ( :, : )
  real ( DP ), allocatable :: k ( :, : )
  real ( DP ), allocatable :: en ( :, :, : )
  real ( DP ), allocatable :: oc ( :, :, : )
  real ( DP ), allocatable :: wfngr ( :, : )
  complex ( DP ), allocatable :: wfngc ( :, : )
  complex ( DP ), allocatable :: wfng_buf ( :, : )
  complex ( DP ), allocatable :: wfng_dist ( :, :, :, : )

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true. )

  IF ( ionode ) CALL iotk_free_unit ( iu )

  IF ( ionode ) THEN
    OPEN ( unit = iu, file = TRIM ( input_file_name ), &
      form = 'unformatted', status = 'old' )
    READ ( iu ) stitle, sdate, stime
  ENDIF

  CALL mp_bcast ( stitle, ionode_id )
  IF(real_or_complex == 1) THEN
    IF(stitle(1:8) /= 'WFN-Real') THEN
      WRITE(0,'(a)') 'Reading from file ' // TRIM(input_file_name) // ':'
      WRITE(0,'(a)') 'Expected header "WFN-Real" but found "' // TRIM(stitle) // '"'
      CALL errore ( 'write_evc', 'file header', 1 )
    ENDIF
  ELSE
    IF(stitle(1:11) /= 'WFN-Complex') THEN
      WRITE(0,'(a)') 'Reading from file ' // TRIM(input_file_name) // ':'
      WRITE(0,'(a)') 'Expected header "WFN-Complex" but found "' // TRIM(stitle) // '"'
      CALL errore ( 'write_evc', 'file header', 1 )
    ENDIF
  ENDIF

  IF ( ionode ) THEN
    READ ( iu ) ns, ng, ntran, cell_symmetry, na, ecutrho, nk, nb, ngkmax, ecutwfn
    READ ( iu ) ( nr ( ir ), ir = 1, 3 )
    READ ( iu ) celvol, al, ( ( a ( j, i ), j = 1, 3 ), i = 1, 3 ), &
      ( ( adot ( j, i ), j = 1, 3 ), i = 1, 3 )
    READ ( iu ) recvol, bl, ( ( b ( j, i ), j = 1, 3 ), i = 1, 3 ), &
      ( ( bdot ( j, i ), j = 1, 3 ), i = 1, 3 )
  ENDIF

  CALL mp_bcast ( ns, ionode_id )
  CALL mp_bcast ( ng, ionode_id )
  CALL mp_bcast ( ntran, ionode_id )
  CALL mp_bcast ( cell_symmetry, ionode_id )
  CALL mp_bcast ( na, ionode_id )
  CALL mp_bcast ( ecutrho, ionode_id )
  CALL mp_bcast ( nk, ionode_id )
  CALL mp_bcast ( nb, ionode_id )
  CALL mp_bcast ( ngkmax, ionode_id )
  CALL mp_bcast ( ecutwfn, ionode_id )
  CALL mp_bcast ( nr, ionode_id )
  CALL mp_bcast ( celvol, ionode_id )
  CALL mp_bcast ( al, ionode_id )
  CALL mp_bcast ( a, ionode_id )
  CALL mp_bcast ( adot, ionode_id )
  CALL mp_bcast ( recvol, ionode_id )
  CALL mp_bcast ( bl, ionode_id )
  CALL mp_bcast ( b, ionode_id )
  CALL mp_bcast ( bdot, ionode_id )

  IF ( ns .NE. nspin ) CALL errore ( 'write_evc', 'ns', 1 )
  IF ( ng .NE. ngm_g ) CALL errore ( 'write_evc', 'ng', 1 )
  IF ( na .NE. nat ) CALL errore ( 'write_evc', 'na', 1 )
  IF ( nk .NE. nkstot / nspin ) CALL errore ( 'write_evc', 'nk', 1 )
  IF ( nr ( 1 ) .NE. nr1 .OR. nr ( 2 ) .NE. nr2 .OR. nr ( 3 ) .NE. nr3 ) &
    CALL errore ( 'write_evc', 'nr', 1 )
  IF ( abs ( celvol - omega ) .GT. eps6 ) &
    CALL errore ( 'write_evc', 'unit cell volume', 1 )
  xdel = 0.0D0
  DO i = 1, 3
    DO j = 1, 3
      xdel = xdel + abs ( al * a ( j, i ) - alat * at ( j, i ) )
    ENDDO
  ENDDO
  IF ( xdel .GT. eps6 ) &
    CALL errore ( 'write_evc', 'direct lattice vectors', 1 )
  xdel = 0.0D0
  DO i = 1, 3
    DO j = 1, 3
      xdel = xdel + abs ( bl * b ( j, i ) - tpiba * bg ( j, i ) )
    ENDDO
  ENDDO
  IF ( xdel .GT. eps6 ) &
    CALL errore ( 'write_evc', 'reciprocal lattice vectors', 1 )

  nbgw = nb
  IF ( wfng_nband .GT. 0 .AND. wfng_nband .LT. nb ) nb = wfng_nband

  IF ( MOD ( ngkmax, nproc ) .EQ. 0 ) THEN
    ngkdist_l = ngkmax / nproc
  ELSE
    ngkdist_l = ngkmax / nproc + 1
  ENDIF
  ngkdist_g = ngkdist_l * nproc

  ALLOCATE ( ngk ( nk ) )
  ALLOCATE ( k ( 3, nk ) )
  ALLOCATE ( en ( nb, nk, ns ) )
  ALLOCATE ( oc ( nb, nk, ns ) )
  ALLOCATE ( gvec ( 3, ng ) )
  ALLOCATE ( igk_buf ( ngkdist_g ) )
  ALLOCATE ( igk_dist ( ngkdist_l, nk ) )
  ALLOCATE ( gk_buf ( 3, ngkdist_g ) )
  ALLOCATE ( gk_dist ( 3, ngkdist_l ) )
  IF ( real_or_complex .EQ. 1 ) THEN
    ALLOCATE ( wfngr ( ngkmax, ns ) )
  ELSE
    ALLOCATE ( wfngc ( ngkmax, ns ) )
  ENDIF
  ALLOCATE ( wfng_buf ( ngkdist_g, ns ) )
  ALLOCATE ( wfng_dist ( ngkdist_l, nb, ns, nk ) )

  IF ( ionode ) THEN
    READ ( iu )
    READ ( iu )
    READ ( iu )
    READ ( iu ) ( ngk ( ik ), ik = 1, nk )
    READ ( iu )
    READ ( iu ) ( ( k ( ir, ik ), ir = 1, 3 ), ik = 1, nk )
    READ ( iu )
    READ ( iu )
    READ ( iu ) ( ( ( en ( ib, ik, is ), ib = 1, nb ), ik = 1, nk ), is = 1, ns )
    READ ( iu ) ( ( ( oc ( ib, ik, is ), ib = 1, nb ), ik = 1, nk ), is = 1, ns )
    READ ( iu ) nrecord
    ig = 1
    DO irecord = 1, nrecord
      READ ( iu ) ng_irecord
      READ ( iu ) ( ( gvec ( ir, jg ), ir = 1, 3 ), jg = ig, ig + ng_irecord - 1 )
      ig = ig + ng_irecord
    ENDDO
  ENDIF

  IF ( ionode ) THEN
    CALL cryst_to_cart ( nk, k, bg, 1 )
    DO is = 1, ns
      DO ik = 1, nk
        DO ib = 1, nb
          en ( ib, ik, is ) = en ( ib, ik, is ) / 2.0D0
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  CALL mp_bcast ( ngk, ionode_id )
  CALL mp_bcast ( k, ionode_id )
  CALL mp_bcast ( en, ionode_id )
  CALL mp_bcast ( oc, ionode_id )
  CALL mp_bcast ( gvec, ionode_id )

  fg = 0
  DO ig = 1, ngm
    IF ( ( gvec ( 1, ig_l2g ( ig ) ) .NE. mill ( 1, ig ) ) .OR. &
         ( gvec ( 2, ig_l2g ( ig ) ) .NE. mill ( 2, ig ) ) .OR. &
         ( gvec ( 3, ig_l2g ( ig ) ) .NE. mill ( 3, ig ) ) ) &
      fg = fg + 1
  ENDDO
  CALL mp_sum ( fg, intra_pool_comm )
  IF ( fg .GT. 0 ) THEN
    IF ( ionode ) WRITE ( 6, 101 )
    gvec ( :, : ) = 0
    DO ig = 1, ngm
      gvec ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
      gvec ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
      gvec ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
    ENDDO
    CALL mp_sum ( gvec, intra_pool_comm )
  ENDIF

  DO ik = 1, nk

    IF ( ionode ) THEN
      WRITE(6, '(a,i8)') 'Working on k-point ', ik

      READ ( iu ) nrecord
      ig = 1
      DO irecord = 1, nrecord
        READ ( iu ) ng_irecord
        READ ( iu ) ( ( gk_buf ( ir, jg ), ir = 1, 3 ), jg = ig, ig + ng_irecord - 1 )
        ig = ig + ng_irecord
      ENDDO
      DO ig = ngk ( ik ) + 1, ngkdist_g
        DO ir = 1, 3
          gk_buf ( ir, ig ) = 0
        ENDDO
      ENDDO
    ENDIF
#ifdef __PARA
    CALL mp_barrier ( )
    CALL MPI_Scatter ( gk_buf, 3 * ngkdist_l, MPI_INTEGER, &
    gk_dist, 3 * ngkdist_l, MPI_INTEGER, &
    ionode_id, world_comm, ierr )
    IF ( ierr .GT. 0 ) CALL errore ( 'write_evc', 'mpi_scatter', ierr )
#else
    DO ig = 1, ngkdist_g
      DO ir = 1, 3
        gk_dist ( ir, ig ) = gk_buf ( ir, ig )
      ENDDO
    ENDDO
#endif
    DO ig = 1, ngkdist_l
      DO jg = 1, ng
        IF ( gk_dist ( 1, ig ) .EQ. gvec ( 1, jg ) .AND. &
             gk_dist ( 2, ig ) .EQ. gvec ( 2, jg ) .AND. &
             gk_dist ( 3, ig ) .EQ. gvec ( 3, jg ) ) THEN
          igk_dist ( ig, ik ) = jg
          EXIT
        ENDIF
      ENDDO
    ENDDO

    DO ib = 1, nb
      IF ( ionode ) THEN
        READ ( iu ) nrecord
        ig = 1
        DO irecord = 1, nrecord
          READ ( iu ) ng_irecord
          IF ( real_or_complex .EQ. 1 ) THEN
            READ ( iu ) ( ( wfngr ( jg, is ), jg = ig, ig + ng_irecord - 1 ), is = 1, ns )
          ELSE
            READ ( iu ) ( ( wfngc ( jg, is ), jg = ig, ig + ng_irecord - 1 ), is = 1, ns )
          ENDIF
          ig = ig + ng_irecord
        ENDDO
        DO is = 1, ns
          IF ( real_or_complex .EQ. 1 ) THEN
            DO ig = 1, ngk ( ik )
              wfng_buf ( ig, is ) = CMPLX ( wfngr ( ig, is ), 0.0D0 )
            ENDDO
          ELSE
            DO ig = 1, ngk ( ik )
              wfng_buf ( ig, is ) = wfngc ( ig, is )
            ENDDO
          ENDIF
          DO ig = ngk ( ik ) + 1, ngkdist_g
            wfng_buf ( ig, is ) = ( 0.0D0, 0.0D0 )
          ENDDO
        ENDDO
      ENDIF
#ifdef __PARA
      DO is = 1, ns
        CALL mp_barrier ( )
        CALL MPI_Scatter ( wfng_buf ( :, is ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
        wfng_dist ( :, ib, is, ik ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
        ionode_id, world_comm, ierr )
        IF ( ierr .GT. 0 ) CALL errore ( 'write_evc', 'mpi_scatter', ierr )
      ENDDO
#else
      DO is = 1, ns
        DO ig = 1, ngkdist_g
          wfng_dist ( ig, ib, is, ik ) = wfng_buf ( ig, is )
        ENDDO
      ENDDO
#endif
    ENDDO

    IF ( ik .LT. nk ) THEN
      DO ib = nb + 1, nbgw
        IF ( ionode ) THEN
          READ ( iu ) nrecord
          DO irecord = 1, nrecord
            READ ( iu )
            READ ( iu )
          ENDDO
        ENDIF
      ENDDO
    ENDIF

  ENDDO

  IF ( ionode ) THEN
    CLOSE ( unit = iu, status = 'keep' )
  ENDIF

  DEALLOCATE ( gk_buf )
  DEALLOCATE ( gk_dist )
  IF ( real_or_complex .EQ. 1 ) THEN
    DEALLOCATE ( wfngr )
  ELSE
    DEALLOCATE ( wfngc )
  ENDIF

  CALL mp_bcast ( ngk, ionode_id )

  nkbl = nkstot / kunit
  nkl = kunit * ( nkbl / npool )
  nkr = ( nkstot - nkl * npool ) / kunit
  IF ( my_pool_id .LT. nkr ) nkl = nkl + kunit
  iks = nkl * my_pool_id + 1
  IF ( my_pool_id .GE. nkr ) iks = iks + nkr * kunit
  ike = iks + nkl - 1

  npw_g = 0
  ALLOCATE ( itmp ( npwx ) )
  DO ik = 1, nks
    DO ig = 1, npwx
      itmp ( ig ) = 0
    ENDDO
    npw = npwx
    CALL gk_sort ( xk ( 1, ik + iks - 1 ), ngm, g, ecutwfc / tpiba2, &
      npw, itmp ( 1 ), g2kin )
    DO ig = 1, npw
      igk_l2g = ig_l2g ( itmp ( ig ) )
      IF ( igk_l2g .GT. npw_g ) npw_g = igk_l2g
    ENDDO
  ENDDO
  DEALLOCATE ( itmp )
  CALL mp_max ( npw_g )

  CALL create_directory ( output_dir_name )
  DO ik = 1, nk
    CALL create_directory ( kpoint_dir ( output_dir_name, ik ) )
  ENDDO

  filename = TRIM ( output_dir_name ) // '/gvectors.dat'

  IF ( ionode ) THEN
    CALL iotk_open_write ( iu, FILE = TRIM ( filename ), SKIP_ROOT = .TRUE., BINARY = .TRUE. )
    CALL iotk_write_begin ( iu, "G-VECTORS" )
    CALL iotk_write_attr ( attr, "nr1s", nr ( 1 ), FIRST = .TRUE. )
    CALL iotk_write_attr ( attr, "nr2s", nr ( 2 ) )
    CALL iotk_write_attr ( attr, "nr3s", nr ( 3 ) )
    CALL iotk_write_attr ( attr, "gvect_number", ng )
    CALL iotk_write_attr ( attr, "gamma_only", .FALSE. )
    CALL iotk_write_attr ( attr, "units", "crystal" )
    CALL iotk_write_empty ( iu, "INFO", ATTR = attr )
    CALL iotk_write_dat ( iu, "g", gvec ( 1 : 3, 1 : ng ), COLUMNS = 3 )
    CALL iotk_write_end ( iu, "G-VECTORS" )
    CALL iotk_close_write ( iu )
  ENDIF

  DO ik = 1, nk

    filename = TRIM ( wfc_filename ( output_dir_name, 'gkvectors', ik ) )

    IF ( ionode ) THEN
      CALL iotk_open_write ( iu, FILE = TRIM ( filename ), ROOT="GK-VECTORS", BINARY = .TRUE. )
      CALL iotk_write_dat ( iu, "NUMBER_OF_GK-VECTORS", ngk ( ik ) )
      CALL iotk_write_dat ( iu, "MAX_NUMBER_OF_GK-VECTORS", ngkmax )
      CALL iotk_write_dat ( iu, "GAMMA_ONLY", .FALSE. )
      CALL iotk_write_attr ( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_dat ( iu, "K-POINT_COORDS", k ( :, ik ), ATTR = attr )
    ENDIF
#ifdef __PARA
    CALL mp_barrier ( )
    CALL MPI_Gather ( igk_dist ( :, ik ) , ngkdist_l, MPI_INTEGER, &
    igk_buf, ngkdist_l, MPI_INTEGER, &
    ionode_id, world_comm, ierr )
    IF ( ierr .GT. 0 ) CALL errore ( 'write_evc', 'mpi_gather', ierr )
#else
    DO ig = 1, ngkdist_g
      igk_buf ( ig ) = igk_dist ( ig, ik )
    ENDDO
#endif
    IF ( ionode ) THEN
      CALL iotk_write_dat ( iu, "INDEX", igk_buf ( 1 : ngk ( ik ) ) )
      CALL iotk_write_dat ( iu, "GRID", gvec ( 1 : 3, igk_buf ( 1 : ngk ( ik ) ) ), COLUMNS = 3 )
      CALL iotk_close_write ( iu )
    ENDIF

    DO is = 1, ns

      IF ( ns .GT. 1 ) THEN
        filename = TRIM ( wfc_filename ( output_dir_name, 'eigenval', ik, is, EXTENSION = 'xml' ) )
      ELSE
        filename = TRIM ( wfc_filename ( output_dir_name, 'eigenval', ik, EXTENSION = 'xml' ) )
      ENDIF

      IF ( ionode ) THEN
        CALL iotk_open_write ( iu, FILE = TRIM ( filename ), BINARY = .FALSE. )
        CALL iotk_write_attr ( attr, "nbnd", nb, FIRST = .TRUE. )
        CALL iotk_write_attr ( attr, "ik", ik )
        IF ( ns .GT. 1 ) CALL iotk_write_attr ( attr, "ispin", is )
        CALL iotk_write_empty ( iu, "INFO", ATTR = attr )
        CALL iotk_write_attr ( attr, "UNITS", "Hartree", FIRST = .TRUE. )
        CALL iotk_write_empty ( iu, "UNITS_FOR_ENERGIES", ATTR = attr )
        CALL iotk_write_dat ( iu, "EIGENVALUES", en ( :, ik, is ) )
        CALL iotk_write_dat ( iu, "OCCUPATIONS", oc ( :, ik, is ) )
        CALL iotk_close_write ( iu )
      ENDIF

      IF ( ns .GT. 1 ) THEN
        filename = TRIM ( wfc_filename ( output_dir_name, 'evc', ik, is ) )
      ELSE
        filename = TRIM ( wfc_filename ( output_dir_name, 'evc', ik ) )
      ENDIF

      IF ( ionode ) THEN
        CALL iotk_open_write ( iu, FILE = TRIM ( filename ), ROOT = "WFC", BINARY = .TRUE. )
        CALL iotk_write_attr ( attr, "ngw", npw_g, FIRST = .TRUE. )
        CALL iotk_write_attr ( attr, "igwx", ngk ( ik ) )
        CALL iotk_write_attr ( attr, "gamma_only", .FALSE. )
        CALL iotk_write_attr ( attr, "nbnd", nb )
        CALL iotk_write_attr ( attr, "ik", ik )
        CALL iotk_write_attr ( attr, "nk", nk )
        CALL iotk_write_attr ( attr, "ispin", is )
        CALL iotk_write_attr ( attr, "nspin", ns )
        CALL iotk_write_attr ( attr, "scale_factor", 1.0D0 )
        CALL iotk_write_empty ( iu, "INFO", attr )
      ENDIF
      DO ib = 1, nb
#ifdef __PARA
        CALL mp_barrier ( )
        CALL MPI_Gather ( wfng_dist ( :, ib, is, ik ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
        wfng_buf ( :, is ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
        ionode_id, world_comm, ierr )
        IF ( ierr .GT. 0 ) CALL errore ( 'write_evc', 'mpi_gather', ierr )
#else
        DO ig = 1, ngkdist_g
          wfng_buf ( ig, is ) = wfng_dist ( ig, ib, is, ik )
        ENDDO
#endif
        IF ( ionode ) CALL iotk_write_dat ( iu, "evc" // iotk_index ( ib ), wfng_buf ( 1 : ngk ( ik ), is ) )
      ENDDO
      IF ( ionode ) CALL iotk_close_write ( iu )

    ENDDO
  ENDDO

  DEALLOCATE ( ngk )
  DEALLOCATE ( k )
  DEALLOCATE ( en )
  DEALLOCATE ( oc )
  DEALLOCATE ( gvec )
  DEALLOCATE ( igk_buf )
  DEALLOCATE ( igk_dist )
  DEALLOCATE ( wfng_buf )
  DEALLOCATE ( wfng_dist )

  CALL mp_barrier ( )

  RETURN

  101 FORMAT ( 5X, "WARNING: reordering G-vectors" )

END SUBROUTINE write_evc

!-------------------------------------------------------------------------------

SUBROUTINE write_cd ( input_file_name, real_or_complex, output_dir_name )

  USE cell_base, ONLY : omega, alat, tpiba, at, bg
  USE constants, ONLY : eps6
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : invfft
  USE grid_dimensions, ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, nl, mill
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat
  USE io_rho_xml, ONLY : write_rho
  USE iotk_module, ONLY : iotk_attlenx, iotk_free_unit, iotk_open_write, &
    iotk_write_begin, iotk_write_attr, iotk_write_empty, iotk_write_dat, &
    iotk_write_end, iotk_close_write
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_bcast, mp_sum
  USE mp_global, ONLY : intra_pool_comm
  USE scf, ONLY : rho
  USE symm_base, ONLY : s, nsym
  USE wavefunctions_module, ONLY : psic
  USE xml_io_base, ONLY : create_directory

  IMPLICIT NONE

  character ( len = 256 ), intent ( in ) :: input_file_name
  integer, intent ( in ) :: real_or_complex
  character ( len = 256 ), intent ( in ) :: output_dir_name

  logical :: f1, f2
  integer :: i, j, iu, is, ig, jg, fg, ir, ns, ng, &
    ntran, cell_symmetry, na, irecord, nrecord, ng_irecord, nr ( 3 )
  real ( DP ) :: ecutrho, celvol, recvol, al, bl, xdel, &
    a ( 3, 3 ), b ( 3, 3 ), adot ( 3, 3 ), bdot ( 3, 3 )
  character :: sdate*32, stime*32, stitle*32
  character ( len = 256 ) :: filename
  character ( iotk_attlenx ) :: attr

  integer, allocatable :: gvec ( :, : )
  real ( DP ), allocatable :: rhogr ( :, : )
  complex ( DP ), allocatable :: rhogc ( :, : )
  complex ( DP ), allocatable :: rhog ( :, : )
  complex ( DP ), allocatable :: rhog_unsrt ( :, : )

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true. )

  IF ( ionode ) CALL iotk_free_unit ( iu )

  IF ( ionode ) THEN
    OPEN ( unit = iu, file = TRIM ( input_file_name ), &
      form = 'unformatted', status = 'old' )
    READ ( iu ) stitle, sdate, stime
  ENDIF

  CALL mp_bcast ( stitle, ionode_id )
  f1 = real_or_complex == 1 .AND. stitle(1:8) == 'RHO-Real'
  f2 = real_or_complex == 2 .AND. stitle(1:11) == 'RHO-Complex'
  IF ( ( .NOT. f1 ) .AND. ( .NOT. f2 ) ) &
    CALL errore ( 'write_cd', 'file header', 1 )

  IF ( ionode ) THEN
    READ ( iu ) ns, ng, ntran, cell_symmetry, na, ecutrho
    READ ( iu ) ( nr ( ir ), ir = 1, 3 )
    READ ( iu ) celvol, al, ( ( a ( j, i ), j = 1, 3 ), i = 1, 3 ), &
      ( ( adot ( j, i ), j = 1, 3 ), i = 1, 3 )
    READ ( iu ) recvol, bl, ( ( b ( j, i ), j = 1, 3 ), i = 1, 3 ), &
      ( ( bdot ( j, i ), j = 1, 3 ), i = 1, 3 )
  ENDIF

  CALL mp_bcast ( ns, ionode_id )
  CALL mp_bcast ( ng, ionode_id )
  CALL mp_bcast ( ntran, ionode_id )
  CALL mp_bcast ( cell_symmetry, ionode_id )
  CALL mp_bcast ( na, ionode_id )
  CALL mp_bcast ( ecutrho, ionode_id )
  CALL mp_bcast ( nr, ionode_id )
  CALL mp_bcast ( celvol, ionode_id )
  CALL mp_bcast ( al, ionode_id )
  CALL mp_bcast ( a, ionode_id )
  CALL mp_bcast ( adot, ionode_id )
  CALL mp_bcast ( recvol, ionode_id )
  CALL mp_bcast ( bl, ionode_id )
  CALL mp_bcast ( b, ionode_id )
  CALL mp_bcast ( bdot, ionode_id )

  IF ( ns .NE. nspin ) CALL errore ( 'write_cd', 'ns', 1 )
  IF ( ng .NE. ngm_g ) CALL errore ( 'write_cd', 'ng', 1 )
  IF ( na .NE. nat ) CALL errore ( 'write_cd', 'na', 1 )
  IF ( nr ( 1 ) .NE. nr1 .OR. nr ( 2 ) .NE. nr2 .OR. nr ( 3 ) .NE. nr3 ) &
    CALL errore ( 'write_cd', 'nr', 1 )
  IF ( abs ( celvol - omega ) .GT. eps6 ) &
    CALL errore ( 'write_cd', 'unit cell volume', 1 )
  xdel = 0.0D0
  DO i = 1, 3
    DO j = 1, 3
      xdel = xdel + abs ( al * a ( j, i ) - alat * at ( j, i ) )
    ENDDO
  ENDDO
  IF ( xdel .GT. eps6 ) &
    CALL errore ( 'write_cd', 'direct lattice vectors', 1 )
  xdel = 0.0D0
  DO i = 1, 3
    DO j = 1, 3
      xdel = xdel + abs ( bl * b ( j, i ) - tpiba * bg ( j, i ) )
    ENDDO
  ENDDO
  IF ( xdel .GT. eps6 ) &
    CALL errore ( 'write_cd', 'reciprocal lattice vectors', 1 )

  ALLOCATE ( gvec ( 3, ng ) )
  IF ( real_or_complex .EQ. 1 ) THEN
    ALLOCATE ( rhogr ( ng, ns ) )
  ELSE
    ALLOCATE ( rhogc ( ng, ns ) )
  ENDIF
  ALLOCATE ( rhog ( ng, ns ) )

  IF ( ionode ) THEN
    READ ( iu )
    READ ( iu )
    READ ( iu )
    READ ( iu ) nrecord
    ig = 1
    DO irecord = 1, nrecord
      READ ( iu ) ng_irecord
      READ ( iu ) ( ( gvec ( ir, jg ), ir = 1, 3 ), jg = ig, ig + ng_irecord - 1 )
      ig = ig + ng_irecord
    ENDDO
    READ ( iu ) nrecord
    ig = 1
    DO irecord = 1, nrecord
      READ ( iu ) ng_irecord
      IF ( real_or_complex .EQ. 1 ) THEN
        READ ( iu ) ( ( rhogr ( jg, is ), jg = ig, ig + ng_irecord - 1 ), is = 1, ns )
      ELSE
        READ ( iu ) ( ( rhogc ( jg, is ), jg = ig, ig + ng_irecord - 1 ), is = 1, ns )
      ENDIF
      ig = ig + ng_irecord
    ENDDO
    DO is = 1, ns
      IF ( real_or_complex .EQ. 1 ) THEN
        DO ig = 1, ng
          rhog ( ig, is ) = CMPLX ( rhogr ( ig, is ), 0.0D0 )
        ENDDO
      ELSE
        DO ig = 1, ng
          rhog ( ig, is ) = rhogc ( ig, is )
        ENDDO
      ENDIF
    ENDDO
    CLOSE ( unit = iu, status = 'keep' )
  ENDIF

  IF ( real_or_complex .EQ. 1 ) THEN
    DEALLOCATE ( rhogr )
  ELSE
    DEALLOCATE ( rhogc )
  ENDIF

  IF ( ionode ) THEN
    DO is = 1, ns
      DO ig = 1, ng
        rhog ( ig, is ) = rhog ( ig, is ) / CMPLX ( omega, 0.0D0 )
      ENDDO
    ENDDO
  ENDIF

  CALL mp_bcast ( gvec, ionode_id )
  CALL mp_bcast ( rhog, ionode_id )

  fg = 0
  DO ig = 1, ngm
    IF ( ( gvec ( 1, ig_l2g ( ig ) ) .NE. mill ( 1, ig ) ) .OR. &
         ( gvec ( 2, ig_l2g ( ig ) ) .NE. mill ( 2, ig ) ) .OR. &
         ( gvec ( 3, ig_l2g ( ig ) ) .NE. mill ( 3, ig ) ) ) &
      fg = fg + 1
  ENDDO
  CALL mp_sum ( fg, intra_pool_comm )
  IF ( fg .GT. 0 ) THEN
    IF ( ionode ) WRITE ( 6, 101 )
    ALLOCATE ( rhog_unsrt ( ng, ns ) )
    rhog_unsrt ( :, : ) = rhog ( :, : )
    rhog ( :, : ) = ( 0.0D0, 0.0D0 )
    DO ig = 1, ng
      DO jg = 1, ngm
        IF ( ( mill ( 1, jg ) .EQ. gvec ( 1, ig ) ) .AND. &
             ( mill ( 2, jg ) .EQ. gvec ( 2, ig ) ) .AND. &
             ( mill ( 3, jg ) .EQ. gvec ( 3, ig ) ) ) THEN
          DO is = 1, ns
            rhog ( ig_l2g ( jg ), is ) = rhog_unsrt ( ig, is )
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE ( rhog_unsrt )
    CALL mp_sum ( rhog, intra_pool_comm )
    gvec ( :, : ) = 0
    DO ig = 1, ngm
      gvec ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
      gvec ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
      gvec ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
    ENDDO
    CALL mp_sum ( gvec, intra_pool_comm )
  ENDIF

  DO is = 1, ns
    DO ig = 1, ngm
      rho%of_g ( ig, is ) = rhog ( ig_l2g ( ig ), is )
    ENDDO
  ENDDO

  DEALLOCATE ( rhog )

  DO is = 1, ns
    DO ig = 1, ngm
      psic ( nl ( ig ) ) = rho%of_g ( ig, is )
    ENDDO
    CALL invfft ( 'Dense', psic, dfftp )
    DO ir = 1, nrxx
      rho%of_r ( ir, is ) = psic ( ir )
    ENDDO
  ENDDO

  CALL create_directory ( output_dir_name )

  filename = TRIM ( output_dir_name ) // '/gvectors.dat'

  IF ( ionode ) THEN
    CALL iotk_open_write ( iu, FILE = TRIM ( filename ), SKIP_ROOT = .TRUE., BINARY = .TRUE. )
    CALL iotk_write_begin ( iu, "G-VECTORS" )
    CALL iotk_write_attr ( attr, "nr1s", nr1, FIRST = .TRUE. )
    CALL iotk_write_attr ( attr, "nr2s", nr2 )
    CALL iotk_write_attr ( attr, "nr3s", nr3 )
    CALL iotk_write_attr ( attr, "gvect_number", ng )
    CALL iotk_write_attr ( attr, "gamma_only", .FALSE. )
    CALL iotk_write_attr ( attr, "units", "crystal" )
    CALL iotk_write_empty ( iu, "INFO", ATTR = attr )
    CALL iotk_write_dat ( iu, "g", gvec ( 1 : 3, 1 : ng ), COLUMNS = 3 )
    CALL iotk_write_end ( iu, "G-VECTORS" )
    CALL iotk_close_write ( iu )
  ENDIF

  DEALLOCATE ( gvec )

  CALL write_rho ( rho%of_r, nspin )

  RETURN

  101 FORMAT ( 5X, "WARNING: reordering G-vectors" )

END SUBROUTINE write_cd

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
    if(invflag .ne. 0 .and. warn) then
      if(ionode) write(6, '(a)') 'WARNING: Inversion symmetry is present. The real version would be faster.'
    endif
  else
    if(invflag .eq. 0) then
      if(real_need_inv) then
        call errore('check_inversion', 'The real version cannot be used without inversion symmetry.', -1)
      endif
      if(ionode) then
        write(6, '(a)') 'WARNING: Inversion symmetry is absent in symmetries used to reduce k-grid.'
        write(6, '(a)') 'Be sure inversion is still a spatial symmetry, or you must use complex version instead.'
      endif
    endif
    if(nspin .eq. 2) then
      call errore('check_inversion', &
        'Time-reversal symmetry is absent in spin-polarized calculation. Complex version must be used.', -2)
    endif
  endif

  return

end subroutine check_inversion

end module bgw2pw_m

!-------------------------------------------------------------------------------

PROGRAM bgw2pw

  USE environment, ONLY : environment_start, environment_end
  USE io_files, ONLY : prefix, tmp_dir, outdir
  USE io_global, ONLY : ionode, ionode_id
  USE kinds, ONLY : DP
  USE mp, ONLY : mp_bcast
  USE mp_global, ONLY : mp_startup
  use bgw2pw_m

  IMPLICIT NONE

  character(len=6) :: codename = 'BGW2PW'

  integer :: real_or_complex
  logical :: wfng_flag
  character ( len = 256 ) :: wfng_file
  integer :: wfng_nband
  logical :: rhog_flag
  character ( len = 256 ) :: rhog_file

  NAMELIST / input_bgw2pw / prefix, outdir, &
    real_or_complex, wfng_flag, wfng_file, wfng_nband, &
    rhog_flag, rhog_file

  integer :: ios
  character ( len = 256 ) :: input_file_name
  character ( len = 256 ) :: output_dir_name

  character (len=256), external :: trimcheck

#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( codename )

  prefix = 'prefix'
  CALL get_env ( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM ( outdir ) == ' ' ) outdir = './'
  real_or_complex = 2
  wfng_flag = .FALSE.
  wfng_file = 'WFN'
  wfng_nband = 0
  rhog_flag = .FALSE.
  rhog_file = 'RHO'

  IF ( ionode ) THEN
    CALL input_from_file ( )
    READ ( 5, input_bgw2pw, iostat = ios )
    IF ( ios /= 0 ) CALL errore ( codename, 'input_bgw2pw', abs ( ios ) )
    IF ( real_or_complex /= 1 .AND. real_or_complex /= 2 ) &
      CALL errore ( codename, 'real_or_complex', 1 )
    IF ( wfng_nband .LT. 0 ) CALL errore ( codename, 'wfng_nband', 1 )
  ENDIF

  tmp_dir = trimcheck ( outdir )
  CALL mp_bcast ( outdir, ionode_id )
  CALL mp_bcast ( tmp_dir, ionode_id )
  CALL mp_bcast ( prefix, ionode_id )
  CALL mp_bcast ( real_or_complex, ionode_id )
  CALL mp_bcast ( wfng_flag, ionode_id )
  CALL mp_bcast ( wfng_file, ionode_id )
  CALL mp_bcast ( wfng_nband, ionode_id )
  CALL mp_bcast ( rhog_flag, ionode_id )
  CALL mp_bcast ( rhog_file, ionode_id )

  CALL read_file ( )
  CALL openfil_pp ( )

  IF ( wfng_flag ) THEN
    input_file_name = TRIM ( outdir ) // '/' // TRIM ( wfng_file )
    output_dir_name = TRIM ( outdir ) // '/' // TRIM ( prefix ) // '.save'
    IF ( ionode ) WRITE ( 6, '(5x,"call write_evc")' )
    CALL start_clock ( 'write_evc' )
    CALL write_evc ( input_file_name, real_or_complex, wfng_nband, &
      output_dir_name )
    CALL stop_clock ( 'write_evc' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_evc",/)' )
  ENDIF

  IF ( rhog_flag ) THEN
    input_file_name = TRIM ( outdir ) // '/' // TRIM ( rhog_file )
    output_dir_name = TRIM ( outdir ) // '/' // TRIM ( prefix ) // '.save'
    IF ( ionode ) WRITE ( 6, '(5x,"call write_cd")' )
    CALL start_clock ( 'write_cd' )
    CALL write_cd ( input_file_name, real_or_complex, output_dir_name )
    CALL stop_clock ( 'write_cd' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_cd",/)' )
  ENDIF

  IF ( ionode ) WRITE ( 6, * )
  IF ( wfng_flag ) CALL print_clock ( 'write_evc' )
  IF ( rhog_flag ) CALL print_clock ( 'write_cd' )

  CALL environment_end ( codename )
  CALL stop_pp ( )

  STOP

END PROGRAM bgw2pw
