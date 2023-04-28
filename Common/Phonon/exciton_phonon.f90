!==============================================================================
!
! Program exciton_phonon
!
! Compute the exciton-phonon coupling matrix elements.
!
! Originally by GKA (2018)
!
!==============================================================================

#include "f_defs.h"

program exciton_phonon

#ifdef HDF5
  use hdf5
#endif
  use hdf5_io_m
  use global_m
  use gkq_m
  use evecs_m
  use xct_ph_maps_m
  use io_utils_m
  use compute_gss_m
  use write_program_header_m
  implicit none

  type(gkq_t) :: gkq, gss
  type(evecs_t) :: evecs_q1, evecs_q2
  character(len=30) :: gkq_fname, gss_fname, evecs_q1_fname, evecs_q2_fname
  character(len=30) :: prog

  integer :: stdout=6
  integer :: error
  integer :: nx1, nx2, ns, nq, nq_xct
  integer :: ix1, ix2, is, iq, iq_xct, ik
  integer :: neig_q1, neig_q2
  integer :: ii
  integer :: nocc
  real(DP) :: qpt(3)   !< Current phonon q-point
  integer, allocatable :: kq_map(:)  !< Mapping between k and k+q
  integer, allocatable :: v_map(:)  !< Mapping between valence bands in A and gkq
  integer, allocatable :: c_map(:)  !< Mapping between conduction bands in A and gkq 

  type(progress_info) :: prog_info !< a user-friendly progress report

  prog = 'exciton_phonon'

  gkq_fname = 'gkq.h5'
  gss_fname = 'gss.h5'

  ! TODO: Should enforce using hdf5
  evecs_q1_fname = 'eigenvectors.h5'
  evecs_q2_fname = 'eigenvectors_q.h5'

  ! Initialize environment

  call peinfo_init()

#ifdef HDF5
  call h5open_f(error)
#endif

  ! Write header
  call write_program_header('Phonon/exciton_phonon', .false.)

  ! ---------------------------------------------------------------------------
  ! Read the gkq matrix

  write(6, *) 'Reading the gkq.h5 file'

#ifdef HDF5
  call gkq%read_and_broadcast_gkq_header_hdf5(gkq_fname)
#endif
  call gkq%setup_paral(peinf, 1)  ! Using only a single processor atm
  !call gkq%setup_paral(peinf, np_qpt=peinf%npes)  ! Serial code atm
  call gkq%alloc()
#ifdef HDF5
  call gkq%read_and_broadcast_gkq_data_hdf5(gkq_fname)
#endif

  qpt = gkq%qpts(:,1)

  ! Print out some info
  !write(6, *) 'qpt=',qpt
  !write(6, *) 'Kpoints list: '
  !do ik = 1, gkq%nk
  !  write(6, '(3f12.8)') (gkq%kpts(ii,ik), ii=1,3)
  !end do

  ! ---------------------------------------------------------------------------
  ! Open and initialize xct eigenvectors (assume tda)

  write(6, *) 'Initializing exciton eigenvectors'

  call evecs_q1%read_file(fname=evecs_q1_fname)
  call evecs_q2%read_file(fname=evecs_q2_fname)

  neig_q1 = evecs_q1%neig
  neig_q2 = evecs_q2%neig

  ! We will be reading a single eigenvector at a time at q1
  ! but all the eigenvectors at q2.
  evecs_q1%meig = 1 ; evecs_q1%neig = 1
  !evecs_q2%meig = 1 ; evecs_q2%neig = 1
  evecs_q2%meig = evecs_q2%neig

  ! Allocate memory
  !call evecs_q1%alloc(with_Avc=.true.)
  !call evecs_q2%alloc(with_Avc=.true.)

  ! TODO: Should perform consistency checks for evecs_q1, evecs_q2, and gkq.

  ! ---------------------------------------------------------------------------

  ! FIXME : nocc should be the TOTAL number of occupied bands.
  ! For now, we assume that all the valence bands were used to solve the BSE
  ! Ideally, this information should be written in the header
  ! of the eigenvector file.
  ! An easier but less desirable solution is to specify it in the input file.
  nocc = evecs_q1%nv

  ! Compute the mapping between the k grid and the k+q grid
  call get_kq_map(kq_map, evecs_q1%kpts, qpt, evecs_q1%nk)  

  ! Compute the mapping between bands in gkq and bands in Avc
  call get_band_maps(v_map, c_map, evecs_q1%nv, evecs_q1%nc, nocc)

  ! ---------------------------------------------------------------------------
  ! Initialize the gss matrix

  ! dimensions
  !nx1 = neig_q1
  nx2 = neig_q2
  nx1 = 4   ! FIXME hard-code a small number of initial exciton states.
  ns = 1  ; is = 1    ! For now, consider non-spin-polarized case
  nq = 1  ; iq = 1    ! Consider a single qpoint
  nq_xct = 1  ; iq_xct = 1

  call gss%init(  &
      gkq%ns,     &  !< ns
      nq_xct,     &  !< nk
      nx1,        &  !< nband
      nx2,        &  !< mband
      nq,         &  !< nq
      gkq%nat     &  !< nat
      )

  call gss%setup_paral(peinf, 1)
  call gss%alloc()
  gss%qpts = gkq%qpts
  gss%kpts(:,:) = ZERO
  gss%g_nu = ZERO

  ! ---------------------------------------------------------------------------
  ! Iterate over pairs of exciton states and compute gss

  call progress_init(prog_info, 'calculation excito-phonon matrix elements', 'i_xct', nx1)

  do ix1=1,nx1

    call progress_step(prog_info)

    do ix2=1,nx2

      ! Main computation of the gss matrix elements
      call compute_gss(gss%g_nu(ix2,ix1,iq_xct,:,:,iq), gkq%g_nu(:,:,:,:,:,iq),&
                       evecs_q1%Avc(:,:,:,:,ix1), evecs_q2%Avc(:,:,:,:,ix2),&
                       kq_map, v_map, c_map, &
                       evecs_q1%ns, evecs_q1%nk, evecs_q1%nc, evecs_q1%nv, &
                       gkq%mband, gkq%nband, gkq%nmode)

    end do
  end do

  call progress_free(prog_info)

  ! ---------------------------------------------------------------------------
  ! Write the result

  write(6, *) 'Writing exciton-phonon matrix elements to disc'

#ifdef HDF5
  call gss%write_hdf5(gss_fname)
#endif

  ! ---------------------------------------------------------------------------

  write(6, *) 'All done'

  ! Free memory
  SAFE_DEALLOCATE(kq_map)
  SAFE_DEALLOCATE(v_map)
  SAFE_DEALLOCATE(c_map)
  call evecs_q1%free()
  call evecs_q2%free()
  call gkq%free()
  call gss%free()

  ! Terminate environment

#ifdef HDF5
  call h5close_f(error)
#endif

#ifdef MPI
  call MPI_Finalize(mpierr)
#endif

end program exciton_phonon
