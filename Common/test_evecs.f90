!==============================================================================
!
! Program test_evecs
!
! Test the evecs_t object io.
!
! Originally by GKA (2019)
!
!==============================================================================

#include "f_defs.h"

program test_evecs

#ifdef HDF5
  use hdf5
  use hdf5_io_m
#endif

  use evecs_m
  use global_m

  type (evecs_t) :: evecs, evecs_ref
  character(len=30) :: fname
  character(len=12) :: prog
  character(len=512) :: msg

  real(DP) :: checksum

  integer :: ieig
  integer :: stdout=6
  integer :: error, ndiff, ndiff_tot

  prog = 'test_evecs: '
#ifdef HDF5
  fname = 'test_evecs.h5'
#else
  fname = 'test_evecs'
#endif

  ndiff_tot = 0

  ! MPI initialization
  call peinfo_init()


#ifdef HDF5
  ! HDF5 initialization
  call h5open_f(error)
#endif

  msg = 'Start'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)


  msg = 'Initializing object'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%init(ns=1, nk=2, nv=3, nc=4, neig=5, tda=.true.)
  call evecs%print_out_header_info(6)


  msg = 'Allocating'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%alloc(with_Avc=.true.)


  msg = 'Assigning'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  evecs%evals = 42.0D0
  evecs%u_r = 1.0D0
  evecs%Avc = 1.0D0
  evecs%has_Avc = .true.
  call evecs%print_out_dimensions(6)


  msg = 'Copying'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%copy(evecs_ref)
  call evecs%compare(evecs_ref, ndiff, verbose=.True.)
  ndiff_tot = ndiff_tot + ndiff


  msg = 'Writing files in serial'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%write_file(fname=fname)


  msg = 'Free memory'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%free()


  msg = 'Reading files'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  evecs%is_distributed = .false.
  call evecs%read_file(fname)
  call evecs%print_out_dimensions(6)
  call evecs%compare(evecs_ref, ndiff, verbose=.True.)
  ndiff_tot = ndiff_tot + ndiff


  msg = 'Free memory'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%free()


  msg = 'Reading files partially'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%read_file(fname, neig_read=2, ieig_offset=1)
  call evecs%print_out_dimensions(6)


  msg = 'Free memory'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%free()


  msg = 'Reading files sequentially'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%open_read(fname)
  call evecs%read_header(broadcast=.false.)
  call evecs%print_out_header_info(6)
  evecs%meig = 1
  call evecs%alloc(with_Avc=.true.)
  do ieig=1,evecs%neig
    call evecs%read_next_eigenvector()
    call evecs%reshape_Avc()
  end do
  call evecs%close_file()


  msg = 'Free memory'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%free()


  msg = 'Reading files'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  !evecs%is_distributed = .false.
  call evecs%read_file(fname)
  call evecs%print_out_dimensions(6)
  call evecs%compare(evecs_ref, ndiff, verbose=.True.)
  ndiff_tot = ndiff_tot + ndiff


  msg = 'Writing files in parallel'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  ! GKA: Note that we setup the parallel distribution scheme
  !      but we do not rearrange the data accordingly. It will not matter here
  !      because all eigenvectos are set to the same constant value.
  call evecs%setup_paral_distribution()
  call evecs%write_file(fname=fname)


  msg = 'Free memory'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%free()


  msg = 'Reading files'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%setup_paral_default()
  call evecs%read_file(fname, neig_read=5)
  call evecs%print_out_dimensions(6)
  call evecs%compare(evecs_ref, ndiff, verbose=.True.)
  ndiff_tot = ndiff_tot + ndiff


  msg = 'Free memory'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%free()


  msg = 'Reading files in parallel'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%read_file(fname, distribute=.true.)
  call evecs%print_out_dimensions(6)


  msg = 'Writing files in parallel'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%write_file(fname=fname)


  msg = 'Free memory'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%free()


  msg = 'Reading files'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  evecs%is_distributed = .false.
  call evecs%read_file(fname, neig_read=5)
  call evecs%print_out_dimensions(6)
  call evecs%compare(evecs_ref, ndiff, verbose=.True.)
  ndiff_tot = ndiff_tot + ndiff


  msg = 'Free memory'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%free()


  msg = 'Reading files on master only'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  evecs%is_distributed = .false.
  call evecs%read_file(fname, distribute=.false., master_only=.true.)
  call evecs%print_out_dimensions(6)
  if (evecs%is_master) then
    call evecs%compare(evecs_ref, ndiff, verbose=.True.)
  else
    ndiff = 0
  end if
  ndiff_tot = ndiff_tot + ndiff


  msg = 'Free memory'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call evecs%free()
  call evecs_ref%free()

  msg = 'Done'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)

  ! Sum up total number of errors
#ifdef MPI
  ndiff = 0
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


#ifdef HDF5
  msg = 'Closing HDF5'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call h5close_f(error)
#endif


#ifdef MPI
  msg = 'Finalizing MPI'
  if (peinf%inode .eq. 0) write(stdout, '(a,1x,a)') prog, trim(msg)
  call MPI_Finalize(mpierr)
#endif

end program test_evecs
