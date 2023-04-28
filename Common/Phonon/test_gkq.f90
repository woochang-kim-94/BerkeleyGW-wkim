!==============================================================================
!
! Program test_gkq
!
! Test the gkq_t object by writing and reading a file,
! then comparing the data read against the original data.
!
! Originally by GKA (2018)
!
!==============================================================================


program test_gkq

#ifdef HDF5
  use hdf5
  use hdf5_io_m
#endif

  use gkq_m

  type(gkq_t) :: gkq, gkq_ref
  character(len=30) :: fname
  character(len=10) :: prog

  real(DP) :: checksum

  integer :: stdout=6
  integer :: error, ndiff, ndiff_tot

  prog = 'test_gkq:'
  fname = 'test_gkq.h5'

  call peinfo_init()

#ifdef HDF5
  call h5open_f(error)
#endif

  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, 'Start'

  ! Initialize gkq object with arbitrary dimensions
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, 'Initializing'
  call gkq%init( &
      1,  &   ! ns
      4,  &   ! nk
      6,  &   ! nband
      7,  &   ! mband
      3,  &   ! nq
      5   &   ! nat
      )

#ifdef MPI
  if (peinf%npes .ge. 1) then
    call gkq%setup_paral(peinf, np_qpt=2)
  end if
#endif

  ! Allocate memory
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, 'Allocating'
  call gkq%alloc()

  ! Set some arbirary values for arrays
  call gkq%set_amu((/12.011_dp, 12.011_dp, 12.011_dp, 12.011_dp, 12.011_dp/))
  call gkq%set_kpts((/0.0_dp, 0.0_dp,  0.0_dp,  &
                      0.0_dp, 0.0_dp,  0.25_dp, &
                      0.0_dp, 0.25_dp, 0.25_dp, &
                      0.0_dp, 0.0_dp,  0.5_dp/))

  ! Set some arbitrary value for the q-points
  call gkq%set_qpts_global((/0.0_dp, 0.0_dp, 0.0_dp, &
                             0.0_dp, 0.0_dp, 0.5_dp, &
                             0.0_dp, 0.0_dp, 0.25_dp/))

  ! Set some arbitrary value for the main g-matrix
  gkq%g_nu = float(peinf%inode)

  ! Copy the arrays onto a reference object
  call gkq%copy(gkq_ref)

#ifdef HDF5
  ! Write the file
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, 'Writing the files'
  call gkq%write_hdf5(fname)
#endif

  ! Flush the memory, for testing purposes
  call gkq%free()
  call gkq%alloc()

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif

#ifdef HDF5
  ! Read the file
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, 'Reading the files'
  call gkq%read_hdf5(fname)

  ! Compare data against reference
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, 'Comparing results'
  call gkq%compare(gkq_ref, ndiff, verbose=.True.)
#else
  ndiff = 0
#endif

#ifdef MPI
  call MPI_REDUCE(ndiff, ndiff_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
#else
  ndiff_tot = ndiff
#endif

  ! Report the result
  if (peinf%inode .eq. 0) then
    write(stdout, '(a,1x,a,1x,i3)') prog, 'Number of errors:', ndiff_tot
    if (ndiff_tot .eq. 0) then
        write(stdout, '(a,1x,a)') prog, 'Success'
    else
        write(stdout, '(a,1x,a)') prog, 'Failure'
    end if
  end if

  ! Free memory
  call gkq%free()
  call gkq_ref%free()

  ! Report
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, 'Done'

#ifdef HDF5
  call h5close_f(error)
#endif

#ifdef MPI
  call MPI_Finalize(mpierr)
#endif

end program test_gkq
