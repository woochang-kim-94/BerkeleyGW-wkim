!=================================================================================
!
! Module read_matrix
!
! (1) read_matrix_d()           Originally by JRD       Last Modified 5/1/2008 (JRD)
!
! This program reads a distributed matrix like chimat or epsmat to file.
!
! (2) read_matrix_f()           Originally by JRD       Last Modified 9/10/2010 (gsm)
!
! Modification of read_matrix_d for full-frequency.
!
!=================================================================================

#include "f_defs.h"

module read_matrix_m

  use global_m
  use scalapack_m
  use epsread_hdf5_m
#ifdef HDF5
  use hdf5
#endif

  implicit none
  
  private

  public :: &
    read_matrix_d, &
    read_matrix_d_hdf5, &
    read_matrix_f, &
    read_matrix_f_hdf5


contains

subroutine read_matrix_d(scal,matrix,nmtx,iunit)
  type (scalapack), intent(in) :: scal
  SCALAR, intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iunit

  PUSH_SUB(read_matrix_d)

  call read_matrix_d_(scal,matrix,nmtx,iunit=iunit)

  POP_SUB(read_matrix_d)

end subroutine read_matrix_d


subroutine read_matrix_d_hdf5(scal,matrix,nmtx,fname,iq,is)
  type (scalapack), intent(in) :: scal
  SCALAR, intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  character(len=*), intent(in) :: fname
  integer, intent(in) :: iq
  integer, intent(in) :: is

  PUSH_SUB(read_matrix_d_hdf5)

  call read_matrix_d_(scal,matrix,nmtx,fname=fname,iq=iq,is=is)

  POP_SUB(read_matrix_d_hdf5)

end subroutine read_matrix_d_hdf5


subroutine read_matrix_d_(scal,matrix,nmtx,iunit,fname,iq,is)
  type (scalapack), intent(in) :: scal
  SCALAR, intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in), optional :: iunit
  character(len=*), intent(in), optional :: fname
  integer, intent(in), optional :: iq
  integer, intent(in), optional :: is
  
  integer :: ii, jj, jjj
  
#ifdef USESCALAPACK
  SCALAR, allocatable :: tempcol(:,:)
  integer :: irow, icol, irowm, icolm, ncol_read, col_blks
  integer :: icurr
#endif
  logical :: use_hdf5
  integer :: c1, c2, cr

  PUSH_SUB(read_matrix_d_)

  if (.not.present(iunit).and..not.(present(fname).and.present(iq))) then
    call die("Not enough arguments to read_matrix_d_", only_root_writes=.true.)
  endif
  if (present(iunit).and.(present(fname).or.present(iq))) then
    call die("Too many arguments to read_matrix_d_", only_root_writes=.true.)
  endif
  if ((present(fname).or.present(iq)).and..not.(present(fname).and.present(iq))) then
    call die("Inconsistent arguments to read_matrix_d_", only_root_writes=.true.)
  endif
  use_hdf5 = present(fname).and.present(iq)
#ifndef HDF5
  if (use_hdf5) then
    call die("read_matrix_d_ was not compiled with HDF5 support.", only_root_writes=.true.)
  endif
#endif

  if (peinf%verb_debug .and. peinf%inode==0) then
    if (use_hdf5) then
      write(6,*) 'Reading matrix: ', nmtx, fname
    else
      write(6,*) 'Reading matrix: ', nmtx, iunit
    endif
    write(6,*)
  endif

#ifdef USESCALAPACK
  
  ncol_read = 1

  if (use_hdf5) then
    ! read chuncks of contigous memory (512 Mb by default)
#ifdef HDF5
    ncol_read = INT( peinf%maxmem_read_h5_mat * 1.0D+06 / (dble(nmtx)*16.0D+00) )
    ncol_read = MIN(nmtx,ncol_read)
    ncol_read = MAX(1,ncol_read)
    if ( peinf%inode==0 ) then
      write(6,'(1x,A,F8.1)') 'Reading matrix using memory blocks of (Mb): ', peinf%maxmem_read_h5_mat
      write(6,'(1x,A,I8)')   'Reading matrix with column block of size:   ', ncol_read
    end if
#endif
  end if

  SAFE_ALLOCATE(tempcol, (nmtx,ncol_read))
  
  icurr=0
  
  do jjj = 1, nmtx, ncol_read

    col_blks = MIN(ncol_read,(nmtx-jjj+1))

!        if (peinf%inode .eq. 0) then
!          write(6,*) ' In loop: ', ii
!        endif

   !XDebugPerfX call SYSTEM_CLOCK( c1, cr)

    if (peinf%inode .eq. 0) then
      if (use_hdf5) then
#ifdef HDF5
        call read_eps_matrix_col_hdf5(tempcol,jjj,col_blks,nmtx,iq,is,fname)
#endif
      else
        read(iunit) (tempcol(ii,1),ii=1,nmtx)
      endif

    endif
    
    !XDebugPerfX call SYSTEM_CLOCK( c2, cr)
    !XDebugPerfX if (peinf%inode .eq. 0) write(*,*) 'READ', dble(c2-c1)/dble(cr)
    !XDebugPerfX call SYSTEM_CLOCK( c1, cr)

    call MPI_BCAST(tempcol,nmtx*ncol_read,MPI_SCALAR,0, &
      MPI_COMM_WORLD,mpierr)
    
    !XDebugPerfX call SYSTEM_CLOCK( c2, cr)
    !XDebugPerfX if (peinf%inode .eq. 0) write(*,*) 'BCast', dble(c2-c1)/dble(cr)
    !XDebugPerfX call SYSTEM_CLOCK( c1, cr)

    do jj = jjj, (jjj+col_blks-1)

      icol=MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
      if (icol .eq. scal%mypcol) then
        do ii = 1, nmtx
          irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
          if (irow .eq. scal%myprow) then
            icurr=icurr+1
            icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
            irowm=MOD((icurr-1),scal%npr)+1
            matrix(irowm,icolm)=tempcol(ii,jj-jjj+1)
          endif
        enddo
      endif

    end do
    
    !XDebugPerfX call SYSTEM_CLOCK( c2, cr)
    !XDebugPerfX if (peinf%inode .eq. 0) write(*,*) 'Redistr', dble(c2-c1)/dble(cr)

    !XXX MDB This barrier is not necessary, MPI_BCAST is an implicit barrier
    ! call MPI_barrier(MPI_COMM_WORLD,mpierr)
    
  enddo
  
  SAFE_DEALLOCATE(tempcol)
  
#else
  
  if(peinf%inode .eq. 0) then
    do jj = 1, nmtx
      if (use_hdf5) then
#ifdef HDF5
        call read_eps_matrix_col_hdf5(matrix(:,jj:jj),jj,1,nmtx,iq,is,fname)
#endif
      else
        read(iunit) (matrix(ii, jj), ii = 1, nmtx)
      endif
    enddo
  endif

#endif
 
  POP_SUB(read_matrix_d_)

  return
end subroutine read_matrix_d_


!=================================================================================


!> FHJ: Front end for read_matrix_f_ for Fortran binary files. See that routine for more info.
subroutine read_matrix_f(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, iunit, advanced)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq
  integer, intent(in) :: nfreq_in_group
  complex(DPC), intent(out) :: retarded(:,:,:) !< (nfreq_in_group,scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: nfreq_group
  integer, intent(in) :: iunit
  complex(DPC), optional, intent(out) :: advanced(:,:,:) !< (nfreq_in_group,scal%npr,scal%npc)

  PUSH_SUB(read_matrix_f)

  call read_matrix_f_(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, iunit=iunit, advanced=advanced)

  POP_SUB(read_matrix_f)

end subroutine read_matrix_f


!> FHJ: Front end for read_matrix_f_ for HDF5 files. See that routine for more info.
subroutine read_matrix_f_hdf5(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, fname, iq, is, advanced)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq
  integer, intent(in) :: nfreq_in_group
  complex(DPC), intent(out) :: retarded(:,:,:) !< (nfreq_in_group,scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: nfreq_group
  character(len=*), intent(in) :: fname
  integer, intent(in) :: iq
  integer, intent(in) :: is
  complex(DPC), optional, intent(out) :: advanced(:,:,:) !< (nfreq_in_group,scal%npr,scal%npc)

  PUSH_SUB(read_matrix_f_hdf5)

  call read_matrix_f_(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, &
    fname=fname, iq=iq, is=is, advanced=advanced)

  POP_SUB(read_matrix_f_hdf5)

end subroutine read_matrix_f_hdf5


!> FHJ: This routines the full-frequency chiR/epsR matrix from a file, and
!! optionally chiA/epsA (note: you shouldn`t really need chiA, ever...)
!! If using HDF5, we only read the retarded part. If legacy
!! Fortran binary, we read the retarded and skip the advanced. The final
!! matrix will be distributed in a ScaLAPACK layout given by scal. Note that
!! this routine is pretty innefficient, but this is not a core component
!! of BGW as it`s only used if you read_chi or use the eps*omega utility.
subroutine read_matrix_f_(scal, nfreq, nfreq_in_group, retarded, nmtx, &
  nfreq_group, iunit, fname, iq, is, advanced)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq
  integer, intent(in) :: nfreq_in_group
  complex(DPC), intent(out) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  integer, intent(in) :: nmtx
  integer, intent(in) :: nfreq_group
  integer, intent(in), optional :: iunit
  character(len=*), intent(in), optional :: fname
  integer, intent(in), optional :: iq
  integer, intent(in), optional :: is
  complex(DPC), intent(out), optional :: advanced(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)

  integer :: ig_glob, igp_glob, ifreq, ifreq_para, freq_grp_ind
#ifdef USESCALAPACK
  complex(DPC), allocatable :: tempcolR(:,:)
  complex(DPC), allocatable :: tempcolA(:,:)
#endif
  logical :: use_hdf5, want_advanced

  PUSH_SUB(read_matrix_f_)

  want_advanced = .false.
#ifdef CPLX
  want_advanced = present(advanced)
#endif
  if (.not.present(iunit).and..not.(present(fname).and.present(iq))) then
    call die("Not enough arguments to read_matrix_f_", only_root_writes=.true.)
  endif
  if (present(iunit).and.(present(fname).or.present(iq))) then
    call die("Too many arguments to read_matrix_f_", only_root_writes=.true.)
  endif
  if ((present(fname).or.present(iq)).and..not.(present(fname).and.present(iq))) then
    call die("Inconsistent arguments to read_matrix_f_", only_root_writes=.true.)
  endif
  use_hdf5 = present(fname).and.present(iq)
#ifndef HDF5
  if (use_hdf5) then
    call die("read_matrix_f_ was not compiled with HDF5 support.", only_root_writes=.true.)
  endif
#endif
  
  if (peinf%verb_debug .and. peinf%inode==0) then
    if (use_hdf5) then
      write(6,*) ' Reading matrix: ', nmtx, fname
    else
      write(6,*) ' Reading matrix: ', nmtx, iunit
    endif
    write(6,*)
  endif
  
  if (peinf%npes>1) then
#ifdef USESCALAPACK

    if (use_hdf5) then
#ifdef HDF5
      if (nfreq_group/=1) then
        call die('Reading polarizability matrix with frequency group not implemnted', &
          only_root_writes=.true.)
      endif
      if (want_advanced) then
        call die('Reading advanced matrix not supported', &
          only_root_writes=.true.)
      endif
      call read_eps_matrix_par_hdf5_general_f(scal, retarded, nmtx, nmtx, &
        nfreq, iq, is, fname, '/mats/matrix', MPI_COMM_WORLD)
#endif
    else !use_hdf5

      SAFE_ALLOCATE(tempcolR, (nmtx,nfreq))
      if (want_advanced) then
        SAFE_ALLOCATE(tempcolA, (nmtx,nfreq))
      endif

      do igp_glob = 1, nmtx
        if (peinf%inode .eq. 0) then
            do ig_glob = 1, nmtx
              read(iunit) (tempcolR(ig_glob,ifreq),ifreq=1,nfreq)
            enddo
#ifdef CPLX
            if (want_advanced) then
              do ig_glob = 1, nmtx
                read(iunit) (tempcolA(ig_glob,ifreq),ifreq=1,nfreq)
              enddo
            else
              do ig_glob = 1, nmtx
                read(iunit)
              enddo
            endif
#endif
        endif
        call MPI_BCAST(tempcolR, nfreq*nmtx, MPI_COMPLEX_DPC, 0, MPI_COMM_WORLD, mpierr)
        if (want_advanced) then
          call MPI_BCAST(tempcolA, nfreq*nmtx, MPI_COMPLEX_DPC, 0, MPI_COMM_WORLD, mpierr)
        endif

        do ifreq = 1, nfreq
          freq_grp_ind = mod(ifreq-1,nfreq_group)
          ifreq_para = (ifreq+nfreq_group-1)/nfreq_group
          if (want_advanced) then
            call save_local_matrix(tempcolR(:,ifreq), tempcolA(:,ifreq))
          else
            call save_local_matrix(tempcolR(:,ifreq))
          endif
        enddo
      enddo
    endif !use_hdf5
    SAFE_DEALLOCATE(tempcolR)
    if (want_advanced) then
      SAFE_DEALLOCATE(tempcolA)
    endif
#endif
  else !USESCALAPACK & peinf%npes>1

    if (use_hdf5) then
#ifdef HDF5
      call read_eps_matrix_ser_f_hdf5(retarded(1:nmtx,1:nmtx,1:nfreq), nmtx, nfreq, iq, is, fname)
#endif
    else !use_hdf5
      do igp_glob = 1, nmtx
        do ig_glob = 1, nmtx
          read(iunit) (retarded(ig_glob, igp_glob, ifreq), ifreq = 1, nfreq)
        enddo
#ifdef CPLX
        if (want_advanced) then
          do ig_glob = 1, nmtx
            read(iunit) (advanced(ig_glob, igp_glob, ifreq), ifreq = 1, nfreq)
          enddo
        else
          do ig_glob = 1, nmtx
            read(iunit)
          enddo
        endif
#endif
      enddo
    endif !use_hdf5
  endif !peinf%npes==1
  
  POP_SUB(read_matrix_f_)

contains

  !> Saves the global `tempcolR`, such as column of the dielectric matrix,
  !! to the local array corresponding to the distributed matrix `retarted`.
  subroutine save_local_matrix(colR, colA)
    complex(DPC), intent(in) :: colR(:)
    complex(DPC), intent(in), optional :: colA(:)

    integer :: igp_loc, ig_glob, ig_loc

    PUSH_SUB(read_matrix_f_.save_local_matrix)

    if ((freq_grp_ind==peinf%igroup_f) .and. (scal%mypcol==&
         INDXG2P(igp_glob, scal%nbl, scal%mypcol, 0, scal%npcol))) then
      igp_loc = INDXG2L(igp_glob, scal%nbl, scal%mypcol, 0, scal%npcol)
      do ig_loc = 1, scal%npr
        ig_glob = INDXL2G(ig_loc, scal%nbl, scal%myprow, 0, scal%nprow)
        retarded(ig_loc, igp_loc, ifreq_para) = colR(ig_glob)
        if (want_advanced) then
          advanced(ig_loc, igp_loc, ifreq_para) = colA(ig_glob)
        endif
      enddo
    endif

    POP_SUB(read_matrix_f_.save_local_matrix)

  end subroutine save_local_matrix

end subroutine read_matrix_f_

end module read_matrix_m
