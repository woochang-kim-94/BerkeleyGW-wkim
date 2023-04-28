!==============================================================================
!
! Program print_gkq
!
! Read a gkq file and print out the content.
!
! Originally by GKA (2018)
!
!==============================================================================

#include "f_defs.h"

program print_gkq

  use gkq_m

  type(gkq_t) :: gkq
  character(len=60) :: fname
  character(len=10) :: prog
  character(len=60) :: usage

  integer :: error
  integer :: is, ik, iband, jband, imode, iq
  integer :: nargs

  prog = 'print_gkq: '
  fname = 'gkq.h5'

  usage = 'Usage: print_gkq.x gkq_file.h5'

  nargs = command_argument_count()
  if (nargs < 1) then
    call die(usage)
  endif
  call get_command_argument(1, fname)


  call peinfo_init()

#ifdef HDF5
  call h5open_f(error)
#endif

  ! Initialize gkq object with arbitrary dimensions
  !if (peinf%inode .eq. 0) write(6,'(a,x,a)') prog, 'Starting'

  ! Read the gkq file
  !call gkq%read_hdf5(fname)
  !call gkq%setup_paral(peinf, np_qpt=1)
#ifdef HDF5
  call gkq%read_and_broadcast_gkq_header_hdf5(fname)
#endif
  call gkq%setup_paral(peinf, np_qpt=1)
  call gkq%alloc()
#ifdef HDF5
  call gkq%read_and_broadcast_gkq_data_hdf5(fname)
#endif

  if (peinf%inode .eq. 0) then

    write(6,'(a12,i4)') 'ns = ', gkq%ns
    write(6,'(a12,i4)') 'nk = ', gkq%nk
    write(6,'(a12,i4)') 'nband = ', gkq%nband
    write(6,'(a12,i4)') 'mband = ', gkq%mband
    write(6,'(a12,i4)') 'nat = ', gkq%nat
    write(6,'(a12,i4)') 'ndir = ', gkq%ndir
    write(6,'(a12,i4)') 'nmode = ', gkq%nmode
    write(6,'(a12,i4)') 'nq = ', gkq%nq

    write(6,'(a12)') ' ik  kpt'
    do ik = 1, gkq%nk
      write(6,'(i6, 3f12.8)') ik, gkq%kpts(:,ik)
    end do

    write(6,'(a12)') ' iq  qpt'
    do iq = 1, gkq%nq_me
      write(6,'(i6, 3f12.8)') iq, gkq%qpts(:,iq)
    end do

#ifdef CPLX
    write(6,'(6a6,2a12)') 'iq','imode','is','ik','iband','jband','Re g', 'Im g'
#else
    write(6,'(6a6,a12)') 'iq','imode','is','ik','iband','jband','Re g'
#endif
    do iq = 1, gkq%nq
      do imode = 1, gkq%nmode
        do is = 1, gkq%ns
          do ik = 1, gkq%nk
            do iband = 1, gkq%nband
              do jband = 1, gkq%mband
#ifdef CPLX
                write(6,'(6i6,2f15.8)') iq,imode,is,ik,iband,jband, &
                  dble(gkq%g_nu(jband,iband,ik,is,imode,iq)), &
                  IMAG(gkq%g_nu(jband,iband,ik,is,imode,iq))
#else
                write(6,'(6i6,f15.8)') iq,imode,is,ik,iband,jband, &
                  dble(gkq%g_nu(jband,iband,ik,is,imode,iq))
#endif
              end do
            end do
          end do
        end do
      end do
    end do

  end if

  call gkq%free()

  !if (peinf%inode .eq. 0) write(6,'(a,x,a)') prog, 'Done'

#ifdef HDF5
  call h5close_f(error)
#endif

#ifdef MPI
  call MPI_Finalize(mpierr)
#endif

end program print_gkq
