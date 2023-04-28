!=================================================================================
!
! Routines:
!
! (1) epscopy()                                 Last Modified 2014 (FHJ)
!
!     Copy dielectric matrices from eps0mat(.h5)/epsmat(.h5) files to memory.
!
!     This routine reads the eps0mat and epsmat files (units 10 and 11)
!     and copies the relevant information about those dielectric matrices
!     to memory (epsmpi). Processor 0 will actually read data and redistribute
!     them to other processors. The eps0mat can have more than one q-point,
!     which is useful when performing voronoi averages for the q->0 point.
!
!==================================================================================

#include "f_defs.h"

module epscopy_m

  use global_m
  use misc_m
  use io_utils_m
  use epsread_hdf5_m
  use scalapack_m
  use timing_m, only: timing => sigma_timing

  implicit none

  private

  public :: epscopy_init, epscopy

  type :: comm_buffer
    complex(DPC), allocatable :: msg(:,:)
    integer, allocatable :: col_global_indx(:)
    integer :: nrow, ncol
    integer :: proc
  end type

contains

!> Determines the number of q-points, the maximum rank of the dielectric matrix,
!! number of frequencies
subroutine epscopy_init(sig, neps)
  type (siginfo), intent(inout) :: sig
  integer, intent(out) :: neps !< Max. rank of dielectric matrix

  integer :: ngarbage1, nmtx, itrash, ig, igp, iq
  integer :: nFreq, nfreq_imag, freq_dep_flag
  real(DP), allocatable :: dFreqGrid(:)
  complex(DPC), allocatable :: dFreqBrd(:)
  logical, allocatable :: is_imag(:)
  logical :: file_exists
  logical :: subtrunc_flags(6)
  integer :: neig_sub, ierr
  real(DP) :: chi_eigenvalue_cutoff
#ifdef HDF5
  integer :: ngarbage3(3), nmtxmax
  real(DP) :: dgarbage1
#endif

  PUSH_SUB(epscopy_init)

  sig%subspace_q0 = .FALSE.
  sig%matrix_in_subspace_basis_q0 = .FALSE.
  sig%keep_full_eps_static_q0 = .FALSE.
  sig%subspace = .FALSE.
  sig%matrix_in_subspace_basis = .FALSE.
  sig%keep_full_eps_static = .FALSE.
  subtrunc_flags(:) = .FALSE.
  sig%neig_sub_max = 0

  if (sig%freq_dep>=0.and.sig%freq_dep<=3) then
#ifdef HDF5
    if (sig%use_hdf5) then
      if (peinf%inode==0) then
        call read_eps_grid_sizes_hdf5(ngarbage1, sig%nq0, dgarbage1, sig%nFreq, &
          sig%nfreq_imag, nmtxmax, ngarbage3, freq_dep_flag, 'eps0mat.h5')
        neps = nmtxmax

        ! Consistency check, before continuing
        if (freq_dep_flag==2.and.sig%freq_dep/=2) &
          call die('eps0mat is frequency-dependent, but this Sigma calculation is not.')
        if (freq_dep_flag==0.and.(sig%freq_dep==2.or.sig%freq_dep==3)) then
          call die('This Sigma calculation is frequency-dependent, but eps0mat is not.')
        endif
        if (freq_dep_flag==3.and.sig%freq_dep/=3) then
          call die('This Sigma calculation is for GN GPP, but eps0mat is not.')
        endif

        ! Check if subspace method has been used
        neig_sub = 0
        call read_eps_subspace_info(sig%subspace_q0, &
                                    sig%matrix_in_subspace_basis_q0, &
                                    sig%keep_full_eps_static_q0, &
                                    neig_sub, chi_eigenvalue_cutoff, 'eps0mat.h5')
        sig%neig_sub_max = MAX(neig_sub, sig%neig_sub_max)
        ! Consistency check for subspace

        if ( sig%subspace_q0 ) then
          if ( sig%matrix_in_subspace_basis_q0 ) then
#if !defined MPI || !defined USESCALAPACK
            call die('Subspace method only works with MPI and SCALAPACK')
#endif
          end if
          if ( sig%freq_dep/=2 ) then
            call die('Subspace approximation implemented only for freq_dep = 2')
          end if
        end if

        inquire(file="epsmat.h5", exist=file_exists) 
        if (file_exists) then
          sig%igamma = 0
        else
          sig%igamma = 1
        endif

        if(sig%igamma/=0) then ! Gamma-only calculation
          sig%nq1 = 0
        else
          call read_eps_grid_sizes_hdf5(ngarbage1, sig%nq1, dgarbage1, nfreq, nfreq_imag, &
            nmtxmax, ngarbage3, freq_dep_flag, 'epsmat.h5')
          if (nFreq/=sig%nFreq) then
            call die('nFreq mismatch between eps0mat.h5 and epsmat.h5')
          endif
          if (nfreq_imag/=sig%nfreq_imag) then
            call die('nfreq_imag mismatch between eps0mat.h5 and epsmat.h5')
          endif
          if (nmtxmax>neps) neps = nmtxmax

          ! check subspace also for epsmat.h5
          neig_sub = 0
          call read_eps_subspace_info(sig%subspace, &
                                      sig%matrix_in_subspace_basis, &
                                      sig%keep_full_eps_static, &
                                      neig_sub, chi_eigenvalue_cutoff, 'epsmat.h5')
          sig%neig_sub_max = MAX(neig_sub, sig%neig_sub_max)
          ! Consistency check for subspace

          if ( sig%subspace ) then
            if ( sig%matrix_in_subspace_basis ) then
#if !defined MPI || !defined USESCALAPACK
              call die('Subspace method only works with MPI and SCALAPACK')
#endif
            end if
            if ( sig%freq_dep/=2 ) then
              call die('Subspace approximation implemented only for freq_dep = 2')
            end if
          end if ! sig%subspace

        endif ! .not.sig%igamma/=0
      endif !peinf%inode==0
    else
#endif
      if (peinf%inode==0) then
        call open_file(unit=10,file='eps0mat',form='unformatted',status='old')
        read(10)
        ierr = 0
        read(10,IOSTAT=ierr) freq_dep_flag, sig%nFreq, &
                 sig%subspace_q0, sig%matrix_in_subspace_basis_q0, sig%keep_full_eps_static_q0
        IF(ierr .NE. 0) THEN
          ! probably this is due because of the absence of the subspace flags
          sig%subspace_q0 = .FALSE.
          sig%matrix_in_subspace_basis_q0 = .FALSE.
          sig%keep_full_eps_static_q0 = .FALSE.
        END IF

! make sure that in the case of subspace MPI and scalapack are linked
#if !defined MPI || !defined USESCALAPACK
        IF(sig%subspace_q0) THEN
          call die('Subspace method only works with MPI and SCALAPACK')
        END IF
#endif

        ! Consistency check, before continuing
        if (freq_dep_flag==2.and.sig%freq_dep/=2) &
          call die('eps0mat is frequency-dependent, but this Sigma calculation is not.')
        if (freq_dep_flag==0.and.sig%freq_dep==2) then
          call die('This Sigma calculation is frequency-dependent, but eps0mat is not.')
        endif
        if (sig%freq_dep/=2 .and. sig%subspace_q0) then
          call die('Subspace approximation implemented only for freq_dep = 2')
        end if

        read(10)
        if (sig%freq_dep==2.or.sig%freq_dep==3) then
          SAFE_ALLOCATE(dFreqGrid, (sig%nFreq))
          SAFE_ALLOCATE(dFreqBrd, (sig%nFreq))
          SAFE_ALLOCATE(is_imag, (sig%nFreq))
          read(10) dFreqGrid, dFreqBrd
          is_imag = (dabs(dFreqGrid)<TOL_ZERO)
          is_imag(1) = .false.
          sig%nfreq_imag = count(is_imag)
          SAFE_DEALLOCATE(dFreqGrid)
          SAFE_DEALLOCATE(dFreqBrd)
          SAFE_DEALLOCATE(is_imag)
        else
          read(10)
        endif
        read(10)
        read(10)
        read(10)
        read(10) sig%nq0
        read(10)
        !XXXXXXXX
        IF(sig%subspace_q0 .AND. sig%matrix_in_subspace_basis_q0) THEN
          do iq=1,sig%nq0
            read(10) itrash, neps
            read(10) ! ekin
            read(10) ! q
            read(10) ! vcoul
            ! full epsilon omega zero
            IF(sig%keep_full_eps_static_q0) THEN
              do igp=1, neps
                read(10) ! For epsRDyn only
              end do 
            END IF
            ! eigenvectors
            read(10) neig_sub
            IF(sig%neig_sub_max < neig_sub) sig%neig_sub_max = neig_sub
            do igp=1, neps
              read(10) ! For epsRDyn only
            end do
            ! subspace matrices
            do igp=1, neig_sub
              do ig=1, neig_sub
                read(10)
              end do 
#ifdef CPLX
              do ig=1, neig_sub
                read(10) ! For epsADyn (empty line...)
              enddo
#endif
            end do
          end do ! iq
        ELSE
          read(10) itrash, neps
        END IF
        !XXXXXXXX
        call close_file(10)      
        call open_file(unit=10,file='eps0mat',form='unformatted',status='old')

        inquire(file="epsmat", exist=file_exists) 
        if (file_exists) then
          sig%igamma = 0
        else
          sig%igamma = 1
        endif
        if(sig%igamma/=0) then ! Gamma-only calculation
          sig%nq1 = 0
        else
          call open_file(unit=11,file='epsmat',form='unformatted',status='old')
          read(11)
          ierr = 0
          read(11,IOSTAT=ierr) ngarbage1, nFreq, & 
                   sig%subspace, sig%matrix_in_subspace_basis, sig%keep_full_eps_static
          IF(ierr .NE. 0) THEN
            sig%subspace = .FALSE.
            sig%matrix_in_subspace_basis = .FALSE.
            sig%keep_full_eps_static = .FALSE.
          END IF

! make sure that in the case of subspace MPI and scalapack are linked
#if !defined MPI || !defined USESCALAPACK
        IF(sig%subspace) THEN
          call die('Subspace method only works with MPI and SCALAPACK')
        END IF
#endif

          if (nFreq/=sig%nFreq) then
            call die('nFreq mismatch between eps0mat and epsmat')
          endif
          if (sig%freq_dep/=2 .and. sig%subspace) then
            call die('Subspace approximation implemented only for freq_dep = 2')
          end if

          read(11)
          if (sig%freq_dep==2.or.sig%freq_dep==3) then
            SAFE_ALLOCATE(dFreqGrid, (sig%nFreq))
            SAFE_ALLOCATE(dFreqBrd, (sig%nFreq))
            SAFE_ALLOCATE(is_imag, (sig%nFreq))
            read(11) dFreqGrid, dFreqBrd
            is_imag = (dabs(dFreqGrid)<TOL_ZERO)
            is_imag(1) = .false.
            if (sig%nfreq_imag/=count(is_imag)) then
              call die('nfreq_imag mismatch between eps0mat and epsmat')
            endif
            SAFE_DEALLOCATE(dFreqGrid)
            SAFE_DEALLOCATE(dFreqBrd)
            SAFE_DEALLOCATE(is_imag)
          else
            read(11)
          endif
          read(11)
          read(11)
          read(11)
          read(11) sig%nq1
          read(11)
        
          do iq=1,sig%nq1
            read(11) itrash, nmtx
            read(11)
            read(11)

            IF(sig%subspace .AND. sig%matrix_in_subspace_basis) THEN
              read(11) ! vcoul
              ! full epsilon omega zero
              IF(sig%keep_full_eps_static) THEN
                do igp=1, nmtx
                  read(11) ! For epsRDyn only
                end do
              END IF
              ! eigenvectors
              read(11) neig_sub
              do igp=1, nmtx
                read(11) ! For epsRDyn only
              end do
              ! subspace matrices
              do igp=1, neig_sub
                do ig=1, neig_sub
                  read(11)
                end do
#ifdef CPLX
                do ig=1, neig_sub
                  read(11) ! For epsADyn (empty line...)
                enddo
#endif
              end do
              ! write(*,*) neig_sub, nmtx
              if(sig%neig_sub_max < neig_sub) sig%neig_sub_max=neig_sub
            ELSE
              if (sig%freq_dep==0.or.sig%freq_dep==1) then
                do ig=1,nmtx
                  read(11)
                enddo
              endif
              if (sig%freq_dep==2.or.sig%freq_dep==3) then
                do igp=1,nmtx
                  do ig=1,nmtx
                    read(11) ! For epsRDyn
                  enddo
#ifdef CPLX
                  do ig=1,nmtx
                    read(11) ! For epsADyn
                  enddo
#endif
                enddo
              endif
            END IF
            if(neps<nmtx) neps = nmtx
          enddo
          call close_file(11)
          call open_file(unit=11,file='epsmat',form='unformatted',status='old')
        endif
      endif ! peinf%inode==0
#ifdef HDF5
    endif
#endif
#ifdef MPI
    call MPI_Bcast(sig%igamma, 1, MPI_INTEGER, 0 ,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(sig%nq0, 1, MPI_INTEGER, 0 ,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(sig%nq1, 1, MPI_INTEGER, 0 ,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(neps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(sig%nFreq, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(sig%nfreq_imag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
!XXXX
    subtrunc_flags(1) = sig%subspace_q0 
    subtrunc_flags(2) = sig%matrix_in_subspace_basis_q0 
    subtrunc_flags(3) = sig%keep_full_eps_static_q0 
    subtrunc_flags(4) = sig%subspace 
    subtrunc_flags(5) = sig%matrix_in_subspace_basis 
    subtrunc_flags(6) = sig%keep_full_eps_static 
    call MPI_Bcast(subtrunc_flags, 6, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
    sig%subspace_q0 = subtrunc_flags(1)
    sig%matrix_in_subspace_basis_q0 = subtrunc_flags(2)
    sig%keep_full_eps_static_q0 = subtrunc_flags(3)
    sig%subspace = subtrunc_flags(4)
    sig%matrix_in_subspace_basis = subtrunc_flags(5)
    sig%keep_full_eps_static = subtrunc_flags(6)
    call MPI_Bcast(sig%neig_sub_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
!XXXX
#endif
    if (sig%freq_dep==2.or.sig%freq_dep==3) then
      SAFE_ALLOCATE(sig%dFreqGrid,(sig%nFreq))
      SAFE_ALLOCATE(sig%dFreqBrd,(sig%nFreq))
    endif
  else ! sig%freq_dep>=0.and.sig_freq_dep<=2
    neps = 0
    sig%nFreq = 0
  endif
  if (sig%nq0/=1.and..not.sig%subsample) call die('epscopy_init: nq0/=1')
  if (peinf%inode==0) then
    write(6,'(/1x,a,i0)') 'Total number of frequencies in the dielectric matrix: ', sig%nfreq
    write(6,'(1x,a,i0/)') 'Number of imag. frequencies in the dielectric matrix: ', sig%nfreq_imag
    !XXXXX
    if(sig%subspace_q0) then
      write(6,'(1x,a)') 'Eps0mat calculated with subspace truncation method'
      IF(sig%matrix_in_subspace_basis_q0) THEN
        write(6,'(1x,a)') 'Eps0mat matrices written within the subspace of omega 0'
      END IF
      IF(sig%keep_full_eps_static_q0) THEN
        write(6,'(1x,a)') 'Eps0mat retains the full EpsInv for omega=0'
      END IF
      write(6,*)
    end if
    if(sig%subspace .AND. sig%igamma==0) then
      write(6,'(1x,a)') 'Epsmat calculated with subspace truncation method'
      IF(sig%matrix_in_subspace_basis) THEN
        write(6,'(1x,a)') 'Epsmat matrices written within the subspace of omega 0'
      END IF
      IF(sig%keep_full_eps_static) THEN
        write(6,'(1x,a)') 'Epsmat retains the full EpsInv for omega=0'
      END IF
      write(6,*)
    end if
    !XXXXX
  endif

  ! check if we can perform the calc in the subspace basis
  if(sig%do_sigma_subspace) then
    if(sig%igamma==0) then
      ! k-point, we need both
      if(sig%matrix_in_subspace_basis_q0 .AND. sig%matrix_in_subspace_basis) then
        if (peinf%inode==0) write(6,'(1x,a)') 'DIRECT evalution of sigma in subspace basis'
      else
        if (peinf%inode==0) write(6,'(1x,a)') 'Epsilon matrices not in subspace basis!'
        if (peinf%inode==0) write(6,'(1x,a)') 'Direct evalution of sigma in subspace basis turned off'
        sig%do_sigma_subspace = .false.
      end if
    else
      if(sig%matrix_in_subspace_basis_q0) then
        if (peinf%inode==0) write(6,'(1x,a)') 'DIRECT evalution of sigma in subspace basis'
      else
        if (peinf%inode==0) write(6,'(1x,a)') 'Epsilon matrix not in subspace basis!'
        if (peinf%inode==0) write(6,'(1x,a)') 'Direct evalution of sigma in subspace basis turned off'
        sig%do_sigma_subspace = .false.
      end if
    end if
  end if

  POP_SUB(epscopy_init)

end subroutine epscopy_init


subroutine epscopy(crys,gvec,sig,neps,epsmpi,epshead,iunit_eps,fne)
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (siginfo), intent(inout) :: sig
  integer, intent(in) :: neps !< number of G-vectors up to epsilon cutoff defined in sigma.inp
  type (epsmpiinfo), intent(inout) :: epsmpi
  SCALAR, intent(out) :: epshead !< for full frequency, this is the retarded static head.
  integer, intent(in) :: iunit_eps
  character*20, intent(in) :: fne

  SCALAR, allocatable :: eps(:)
  complex(DPC), allocatable :: epsRDyn(:,:)
  complex(DPC), allocatable :: epsADyn(:,:)
  logical :: is_static
#ifdef HDF5
  integer, allocatable :: sub_neig_of_q(:)
#endif

  PUSH_SUB(epscopy)

  epshead = ZERO
  is_static = sig%freq_dep/=2 .AND. sig%freq_dep/=3

  ! FHJ: Initialize buffers
  if (is_static) then
    SAFE_ALLOCATE(epsmpi%eps, (neps,epsmpi%ngpown,sig%nq))
  else
    if(sig%do_sigma_subspace) then
      ! if doing calcultion within subspace allocate stuff
      SAFE_ALLOCATE(sig%epssub%eps_sub, (sig%neig_sub_max, sig%epssub%Nbas_own, sig%nFreq, sig%nq))
      !XXXX SAFE_ALLOCATE(sig%epssub%eigenvec_sub, (sig%epssub%ngpown_sub, sig%neig_sub_max,  sig%nq))
      SAFE_ALLOCATE(sig%epssub%eigenvec_sub, (neps, sig%epssub%Nbas_own,  sig%nq))
      SAFE_ALLOCATE(sig%epssub%eps_wings_rows, (neps, sig%nFreq,  sig%nq))
      SAFE_ALLOCATE(sig%epssub%eps_wings_cols, (neps, sig%nFreq,  sig%nq))
      SAFE_ALLOCATE(sig%epssub%eps_wings_correction_rows, (neps, sig%nFreq))
      SAFE_ALLOCATE(sig%epssub%eps_wings_correction_cols, (neps, sig%nFreq))
      sig%epssub%eps_sub = (0.0d0,0.0d0)
      sig%epssub%eigenvec_sub = (0.0d0,0.0d0)
      sig%epssub%eps_wings_rows = (0.0d0,0.0d0)
      sig%epssub%eps_wings_cols = (0.0d0,0.0d0)
      sig%epssub%eps_wings_correction_rows = (0.0d0,0.0d0)
      sig%epssub%eps_wings_correction_cols = (0.0d0,0.0d0)
      SAFE_ALLOCATE(sig%epssub%eps_sub_info, (3, 2, sig%nq))
      SAFE_ALLOCATE(sig%epssub%wing_pos, (sig%nq))
      sig%epssub%wing_pos = 0
      SAFE_ALLOCATE(sig%epssub%vcoul_sub, (neps, sig%nq))
      sig%epssub%vcoul_sub = 0.0d0
      !XXX deallocate and reallocate epsmpi%epsR with a single frequency
      !XXX SAFE_DEALLOCATE_P(epsmpi%epsR)
      SAFE_ALLOCATE(epsmpi%epsR, (neps,epsmpi%ngpown,1,sig%nq))
    else
      ! standard case
      SAFE_ALLOCATE(epsmpi%epsR, (neps,epsmpi%ngpown,sig%nFreq,sig%nq))
      if (sig%need_advanced) then
        SAFE_ALLOCATE(epsmpi%epsA, (neps,epsmpi%ngpown,sig%nFreq,sig%nq))
      endif
    end if
  endif

  ! FHJ: Temp buffers. 
  if (is_static) then
    SAFE_ALLOCATE(eps, (neps))
  else
    SAFE_ALLOCATE(epsRDyn, (neps,sig%nFreq))
    if (sig%need_advanced) then
      SAFE_ALLOCATE(epsADyn, (neps,sig%nFreq))
    endif
  endif
  !------------------------------------------------
  ! FHJ: Read dielectric matrices from eps0mat(.h5)
  !------------------------------------------------
  call read_epsmat(.true.)
#ifdef MPI
  call MPI_Bcast(epshead, 1, MPI_SCALAR, 0, MPI_COMM_WORLD, mpierr)
#endif
  if (dble(epshead)<1d-3 .and. sig%iscreen==SCREEN_SEMICOND .and. peinf%inode==0) then
    write(0,'(a)') 'WARNING: You are using semiconductor screening, but the'
    write(0,'(a)') 'head of epsilon inverse is very small and seems metallic.'
  endif
  !------------------------------------------------
  ! FHJ: Read dielectric matrices from epsmat(.h5)
  !------------------------------------------------
  if (sig%igamma==0) call read_epsmat(.false.)

#ifdef MPI
  call MPI_Bcast(sig%qgrid,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
  epsmpi%qk(:,:) = sig%qpt(:,:)

  ! FHJ: Free temp. buffers
  if (is_static) then
    SAFE_DEALLOCATE(eps)
  else
    SAFE_DEALLOCATE(epsRDyn)
    if (sig%need_advanced) then
      SAFE_DEALLOCATE(epsADyn)
    endif
  endif

  POP_SUB(epscopy)

  return

contains


  subroutine read_epsmat(is_q0)
    logical, intent(in) :: is_q0

    character(len=16) :: fname
    integer :: qoffset, nq_file, iunit
    integer :: freq_dep_flag, nFreq, nfreq_imag, ng_old, ngq, nmtx
    integer :: ii, ig, igp, nq_tmp, iq, qgrid(3), gg(3), iout, iw, dest, igp_loc
    real(DP) :: ecuts, qk(3), ekin
    real(DP), allocatable :: ekold(:)
#ifdef HDF5
    integer :: nmtxmax
    integer, allocatable :: nmtx_of_q(:), isrtinvdummy(:)
#endif
    integer, allocatable :: isrtq(:), isrtqi(:), isrtold(:), gvecs_old(:,:)
    type(progress_info) :: prog_info !< a user-friendly progress report
    character(len=6) :: sname
    character(len=11) :: sdate
    logical :: matrix_in_subspace_basis, keep_full_eps_static
    real(DP), allocatable :: vcoul(:)

    PUSH_SUB(epscopy.read_epsmat)

    if (is_q0) then
      fname = 'eps0mat'
      qoffset = 0
      nq_file = sig%nq0
      iunit = 10
    else
      fname = 'epsmat'
      qoffset = sig%nq0
      nq_file = sig%nq1
      iunit = 11
    endif
    if (sig%use_hdf5) fname = TRUNC(fname) // '.h5'

    if (peinf%inode==0) then
      SAFE_ALLOCATE(ekold, (gvec%ng))
      SAFE_ALLOCATE(isrtold, (gvec%ng))
      SAFE_ALLOCATE(isrtq, (gvec%ng))
      SAFE_ALLOCATE(isrtqi, (gvec%ng))
    endif

    matrix_in_subspace_basis = .FALSE.
    keep_full_eps_static = .FALSE.
    if (is_q0) then
      if(sig%subspace_q0) then
        matrix_in_subspace_basis = sig%matrix_in_subspace_basis_q0
        keep_full_eps_static = sig%keep_full_eps_static_q0
      end if
    else
      if(sig%subspace) then
        matrix_in_subspace_basis = sig%matrix_in_subspace_basis
        keep_full_eps_static = sig%keep_full_eps_static
      end if
    end if

    call timing%start(timing%io_total)
!------------------------------------------------------------------------------
! Read q-grid, q-points, G-vectors and freq. grid
!------------------------------------------------------------------------------
    if (peinf%inode==0) then
#ifdef HDF5
      if (sig%use_hdf5) then
        sname = 'chiGG0'
        sdate = 'nodate'
        call read_eps_grid_sizes_hdf5(ng_old, nq_tmp, ecuts, nfreq, nfreq_imag, &
          nmtxmax, qgrid, freq_dep_flag, TRUNC(fname))
        if (.not.is_static.and.is_q0) then
          call read_eps_freqgrid_hdf5(nFreq, sig%dFreqGrid, sig%dFreqBrd, TRUNC(fname))
        endif

        SAFE_ALLOCATE(nmtx_of_q, (nq_file))
        call read_eps_qgrid_hdf5(nq_tmp, sig%qpt(:,qoffset+1:), nmtx_of_q, TRUNC(fname))
        ! Note: it seems like the old format stores isort up to ngq, and the
        ! HDF5 stores up to ng_old
        SAFE_ALLOCATE(gvecs_old, (3,ng_old))
        call read_eps_old_gvecs_hdf5(ng_old, gvecs_old, TRUNC(fname))
        ! subspace: read number of eigenvectors
        if(matrix_in_subspace_basis) then
          SAFE_ALLOCATE(sub_neig_of_q, (nq_file))
          sub_neig_of_q = 0
          call read_eps_sub_neig_of_q(nq_tmp, sub_neig_of_q, TRUNC(fname))
        end if
      else
#endif
        read(iunit) sname, sdate
        read(iunit) freq_dep_flag, nFreq
        read(iunit) qgrid(1:3)
        if (.not.is_static.and.is_q0) then
          read(iunit) sig%dFreqGrid, sig%dFreqBrd
        else
          read(iunit)
        endif
        read(iunit)
        read(iunit)
        read(iunit) ecuts
        read(iunit) nq_tmp, ((sig%qpt(ii,qoffset+iq), ii=1,3), iq=1,nq_tmp)
        read(iunit) ng_old
        backspace(iunit)
        SAFE_ALLOCATE(gvecs_old, (3, ng_old))
        read(iunit) ng_old, ((gvecs_old(ii,iq), ii=1,3), iq=1,ng_old)
#ifdef HDF5
      endif
#endif
      if (nq_tmp/=nq_file) call die('nq mismatch for '//TRUNC(fname))
      ! FHJ: only substitute sig%qgrid (used in voronoi averages) if sig%qgrid
      ! is empty and this is epsmat(.h5), if this file is available
      if (all(sig%qgrid(1:3)==0).and.((.not.is_q0).or.(sig%igamma==1))) then
        sig%qgrid(:) = qgrid(:)
      endif
    endif ! peinf%inode==0
    call timing%stop(timing%io_total)

!------------------------------------------------------------------------------
! Bcast q-points and allocate/bcast full-frequency and epsilon buffers
!------------------------------------------------------------------------------
#ifdef MPI
    call MPI_Bcast(sig%qpt(1,qoffset+1), 3*nq_file, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
#endif
    if (is_q0) then
      sig%q0vec(:) = sig%qpt(:,1)
      ! We need to manually set the q0 point in sig%qpt to zero
      if (.not.sig%subsample) sig%qpt(:,:sig%nq0) = 0d0
      if (.not.is_static) then
#ifdef MPI
        call MPI_Bcast(sig%dFreqGrid,sig%nFreq,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(sig%dFreqBrd,sig%nFreq,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
#endif
      endif
    endif !is_q0


    call timing%start(timing%io_total)
!------------------------------------------------------------------------------
! Loop over all q-points, map G-vectors and read/store dielectric matrices
!------------------------------------------------------------------------------
    call progress_init(prog_info, 'reading '//TRUNC(fname), 'q-point', nq_file)
    do iq=1,nq_file
      call progress_step(prog_info, iq)
      call logit('Storing eps to memory')

      ! FHJ: Map old isrt (from epsilon gspace) to new isrt (to gvec gspace)
      !--------------------------------------------------------=------------
      if (peinf%inode==0) then
#ifdef HDF5
        if (sig%use_hdf5) then
          SAFE_ALLOCATE(isrtinvdummy, (ng_old))
          call read_eps_gvecsofq_hdf5(ng_old,isrtold,isrtinvdummy,ekold,iq,TRUNC(fname))
          SAFE_DEALLOCATE(isrtinvdummy)
          nmtx = nmtx_of_q(iq)
          ngq = ng_old
          ! in case of subspace read the coulomb potential
          if(matrix_in_subspace_basis) then
            SAFE_ALLOCATE(vcoul, (nmtx))
            call read_vcoul_hdf5(vcoul, iq, TRUNC(fname))
          end if
        else
#endif
          read(iunit) ngq, nmtx, (isrtold(ig),ii,ig=1,ngq)
          read(iunit) (ekold(ig),ig=1,ngq)
          read(iunit) (qk(ii),ii=1,3)
          ! in case of subspace read the coulomb potential
          if(matrix_in_subspace_basis) then
            SAFE_ALLOCATE(vcoul, (nmtx))
            read(iunit) (vcoul(ig),ig=1,nmtx)
          end if
#ifdef HDF5
        endif
#endif
        isrtq = 0
        isrtqi = 0
        qk = sig%qpt(:,iq+qoffset)
        if (is_q0) qk(:) = 0d0
        ! FHJ: the following check seems to be important when the WFN file used in the
        ! epsilon calculation has a cutoff larger than that of the WFN file used here.
        ngq = min(ngq, gvec%ng)
        do ig=1,ngq
          ! FHJ: note that size(ekold) is gvec%ng
          if (isrtold(ig) <= gvec%ng) then
            if (ekold(isrtold(ig))<=sig%ecutb) then
              gg(:) = gvecs_old(:, isrtold(ig))
              call findvector(iout, gg, gvec)
              isrtq(ig) = iout
              isrtqi(iout) = ig
              if (ig==1) then
                ! just check the first so we do not waste time
                ekin = DOT_PRODUCT(gvec%components(:,iout)+qk(:), MATMUL(crys%bdot, gvec%components(:,iout)+qk(:)))
                if (abs(ekin-ekold(isrtold(ig))) > TOL_Zero) then
                  write(0,*) 'iq = ', iq, ' ig = ',ig, ' qk = ', qk
                  write(0,*) 'epsmat: ekold(isrtold(i)) = ', ekold(isrtold(ig)), ' ekin = ', ekin
                  call die("Incorrect kinetic energies in epsmat.")
                endif
              endif
            endif  ! if (ekold(isrtold(i))<=sig%ecutb)
          endif
        enddo  ! do ig=1,ngq
      endif  ! if (peinf%inode==0)
#ifdef MPI
      call MPI_Bcast(nmtx,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
      if (peinf%inode==0) then
        epsmpi%isrtq(:,iq+qoffset) = isrtq(:)
        epsmpi%isrtqi(:,iq+qoffset) = isrtqi(:)
      endif
      epsmpi%nmtx(iq+qoffset) = nmtx

      IF(matrix_in_subspace_basis) THEN
        ! broadcast coulomb potential
        IF(peinf%inode/=0) THEN
          SAFE_ALLOCATE(vcoul, (nmtx))
          vcoul = 0.0_DP
        END IF
#ifdef MPI
        call MPI_Bcast(vcoul,nmtx,MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
#endif

#if defined MPI && defined USESCALAPACK
        call timing%start(timing%epscopy_sub)
        call epscopy_subspace(nmtx,vcoul,keep_full_eps_static,iq,qoffset,iunit,fname)
        call timing%stop(timing%epscopy_sub)
#endif

        IF(sig%need_advanced) THEN
          write(6,'(1x,a)') &
          'Subspace truncation method only implemented for sig%freq_dep==2 and sig%freq_dep_method==2'
          call die('need_advanced is an invalid option for subspace + cd_integration_method.')
        END IF

        if (is_q0.and.iq==1.and.peinf%inode==0) then
          epshead = epsmpi%epsR(1,1,1,1)
        endif

        SAFE_DEALLOCATE(vcoul)
      ELSE ! subspace read
#ifdef HDF5
        if (sig%use_hdf5) then
          ! FHJ: Read & store epsinv === HDF5 ===
          !--------------------------------------
          if (is_static) then
            call timing%start(timing%epscopy_io)
            call read_eps_matrix_par_hdf5(epsmpi%eps(:,:,iq+qoffset), &
              epsmpi%nb, peinf%pool_rank, peinf%npes_pool, nmtx, &
              iq, 1, fname)
            call timing%stop(timing%epscopy_io)
          else
            call timing%start(timing%epscopy_io)
            if (sig%need_advanced) then
              call read_eps_matrix_par_f_hdf5(epsmpi%epsR(:,:,:,iq+qoffset), &
                epsmpi%nb, peinf%pool_comm, peinf%my_pool, peinf%npools, nmtx, &
                sig%nFreq, iq, 1, fname, advanced=epsmpi%epsA(:,:,:,iq+qoffset))
            else
              call read_eps_matrix_par_f_hdf5(epsmpi%epsR(:,:,:,iq+qoffset), &
                epsmpi%nb, peinf%pool_comm, peinf%my_pool, peinf%npools, nmtx, &
                sig%nFreq, iq, 1, fname)
            endif
            call timing%stop(timing%epscopy_io)
          endif
          if (is_q0.and.iq==1.and.peinf%inode==0) then
            if (is_static) then
              epshead = epsmpi%eps(1,1,1)
            else
              epshead = epsmpi%epsR(1,1,1,1)
            endif
          endif
        else !HDF5
#endif
          ! FHJ: Read & store epsinv === Fortran binary  ===
          !-------------------------------------------------
          do igp=1,nmtx
            if (is_static) then
              if (peinf%inode==0) then
                call timing%start(timing%epscopy_io)
                read(iunit) (eps(ig),ig=1,nmtx)
                call timing%stop(timing%epscopy_io)
                if (is_q0.and.iq==1.and.igp==1) epshead = eps(1)
                call timing%start(timing%epscopy_comm)
                call timing%stop(timing%epscopy_comm)
              endif
  
              call timing%start(timing%epscopy_comm)
#ifdef MPI
              call MPI_Bcast(eps, neps, MPI_SCALAR, 0, MPI_COMM_WORLD, mpierr)
#endif
              call timing%stop(timing%epscopy_comm)
            else ! is_static
              if (peinf%inode==0) then
                do ig=1,nmtx
                  call timing%start(timing%epscopy_io)
                  read(iunit) (epsRDyn(ig,iw),iw=1,sig%nFreq) ! Retarded part
                  call timing%stop(timing%epscopy_io)
                  if (is_q0.and.iq==1.and.igp==1) epshead = epsRDyn(1,1)
                  call timing%start(timing%epscopy_comm)
                  call timing%stop(timing%epscopy_comm)
                enddo
#ifdef CPLX
                do ig=1,nmtx
                  call timing%start(timing%epscopy_io)
                  if (sig%need_advanced) then
                    read(iunit) (epsADyn(ig,iw),iw=1,sig%nFreq) ! Advanced part
                  else
                    read(iunit) ! Advanced part
                  endif
                  call timing%stop(timing%epscopy_io)
                enddo
#endif
              endif
  
              call timing%start(timing%epscopy_comm)
#ifdef MPI
              call MPI_Bcast(epsRDyn, sig%nFreq*neps, MPI_COMPLEX_DPC, 0, MPI_COMM_WORLD, mpierr)
              if (sig%need_advanced) then
                call MPI_Bcast(epsADyn, sig%nFreq*neps, MPI_COMPLEX_DPC, 0, MPI_COMM_WORLD, mpierr)
              endif
#endif
              call timing%stop(timing%epscopy_comm)
            endif ! is_static
  
            dest = INDXG2P(igp, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
            igp_loc = INDXG2L(igp, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
            call timing%start(timing%epscopy_comm)
            if (peinf%pool_rank==dest) then
              if (is_static) then
                epsmpi%eps(:,igp_loc,iq+qoffset) = eps(:)
              else
                epsmpi%epsR(:,igp_loc,:,iq+qoffset) = epsRDyn(:,:)
                if (sig%need_advanced) then
                  epsmpi%epsA(:,igp_loc,:,iq+qoffset) = epsADyn(:,:)
                endif
              endif ! is_static
            endif ! if (peinf%pool_rank==dest)
            call timing%stop(timing%epscopy_comm)
          enddo  ! do ig=1,nmtx
#ifdef HDF5
        endif !HDF5
#endif
      END IF ! subspace read
    enddo ! iq
    call progress_free(prog_info)

    call timing%stop(timing%io_total)
    
!------------------------------------------------------------------------------
! Done!
!------------------------------------------------------------------------------

    if (peinf%inode==0) then
#ifdef HDF5
      if (sig%use_hdf5) then
        SAFE_DEALLOCATE(nmtx_of_q)
        if(matrix_in_subspace_basis) then
          SAFE_DEALLOCATE(sub_neig_of_q)
        end if
      else
#endif
        call close_file(iunit)
#ifdef HDF5
      endif
#endif
      SAFE_DEALLOCATE(ekold)
      SAFE_DEALLOCATE(isrtold)
      SAFE_DEALLOCATE(isrtq)
      SAFE_DEALLOCATE(isrtqi)
    endif

    POP_SUB(epscopy.read_epsmat)

  end subroutine read_epsmat

#if defined MPI && defined USESCALAPACK
  subroutine epscopy_subspace(nmtx,vcoul,keep_full_eps_static,iq,qoffset,iunit,fname)
    integer  :: nmtx
    real(DP) :: vcoul(nmtx)
    logical  :: keep_full_eps_static
    integer  :: iq, qoffset, iunit
    character(len=16), intent(in) :: fname

    complex(DPC), allocatable :: eigenvect(:,:), C_aux(:,:), buffer(:,:)
    complex(DPC), allocatable :: eps_sub_f(:,:,:), E_aux(:,:,:), buffer_f(:,:,:)
    complex(DPC), allocatable :: eps1Aux(:,:), C_Pgemm(:,:)
    complex(DPC) :: C_send(nmtx), C_recv(nmtx)
    real(DP) :: vc
    type (scalapack) :: scal, scal_sub, scal_1d, scal_aux
    integer, allocatable :: row_indices(:,:), col_indices(:,:)
    integer, allocatable :: row_indices_sub(:,:), col_indices_sub(:,:)
    integer :: desca(9), desc_sub(9), info
    integer :: max_npr, max_npc, tag
    integer :: nrow_loc_1d, ncol_loc_1d
    integer :: nrow_loc_sub, ncol_loc_sub, max_npr_sub, max_npc_sub
    integer :: igp, ig, ipe, neig_sub, Nfreq, ifreq, ig_l, ig_g, igp_l, igp_g
    integer :: i_local, i_global, j_local, j_global 
    integer :: icurr, i, irow, j, icol, irowm, icolm 
 
    integer :: Nmsg_to_send, Nmsg_to_recv
    type(comm_buffer), allocatable :: buf_send(:), buf_recv(:)
    integer, allocatable :: pe2grid(:,:), grid2pe(:,:), pool2pe(:,:) 
    integer, allocatable :: num_col_to_send(:), num_col_to_recv(:)
    integer, allocatable :: proc2msg_send(:), local_col_count(:)
    integer, allocatable :: req_send(:), req_recv(:)
    integer, allocatable :: my_status(:,:)
    integer :: iproc_dum, iproc_send, iproc_rec, pool_send, iproc_row_rec
    integer :: Nbuf_send, Nbuf_recv, imsg_recv, imsg_send, my_col
    ! variables for sigma-subspace
    complex(DPC), allocatable :: buffer_sub_eigen(:,:),  buffer_sub_eps(:,:,:)
    integer :: indx_start, indx_end, sub_size
    integer :: iproc_pool, iproc_row, iproc_col
    integer :: max_head, gg(3), iout
    integer :: ipe_wing, my_pos_loc_wing
    ! variables for fast vcoul scaling
    real(DP), allocatable :: row_vcoul_aux(:), col_vcoul_aux(:)
    integer, allocatable :: my_diag_indeces(:,:)
    integer :: idiag, my_num_diag_elem
    ! variables for cyclic redistribution of eigenvectors and subspace epsilon
    logical :: cyclic_redistr
    complex(DPC), allocatable :: rec_buffer_sub_eigen(:,:), send_buffer_sub_eigen(:,:) 
    complex(DPC), allocatable :: rec_buffer_sub_eps(:,:,:), send_buffer_sub_eps(:,:,:)
    integer :: max_npc_lt_Neig, npc_lt_Neig
    integer :: req_send_singl, tag_send, req_rec_singl, tag_rec
    integer :: irec_static, isend_static
    integer :: actual_send, actual_rec, ipe_real
    logical :: do_copy
    integer :: ncol_to_copy
    integer, allocatable :: copy_col_inx(:)

    PUSH_SUB(epscopy.epscopy_subspace)

    ! set default cyclic_redistr flag
    cyclic_redistr = .true.
    if( sig%sub_collective_redistr ) cyclic_redistr = .false.

    ! set up the blas env for the full matrices
    call blacs_setup(scal, nmtx, .true.)   

    ! exchange info
    SAFE_ALLOCATE(scal%nprd, (peinf%npes))
    SAFE_ALLOCATE(scal%npcd, (peinf%npes))
    scal%nprd = 0
    scal%npcd = 0
    scal%nprd(peinf%inode + 1) = scal%npr
    scal%npcd(peinf%inode + 1) = scal%npc
    !
    call MPI_ALLREDUCE(MPI_IN_PLACE,scal%nprd,peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,scal%npcd,peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)
    ! get max number of row/col
    max_npr = MAXVAL(scal%nprd)
    max_npc = MAXVAL(scal%npcd)

    ! precompute indices
    SAFE_ALLOCATE(row_indices, (max_npr, peinf%npes))
    row_indices = 0
    DO i_local = 1, scal%npr
      i_global = indxl2g(i_local,  scal%nbl, scal%myprow, 0, scal%nprow)
      row_indices(i_local,peinf%inode + 1) = i_global
    END DO
    SAFE_ALLOCATE(col_indices, (max_npc, peinf%npes))
    col_indices = 0
    DO i_local = 1, scal%npc
      i_global = indxl2g(i_local,  scal%nbl, scal%mypcol, 0, scal%npcol)
      col_indices(i_local,peinf%inode + 1) = i_global
    END DO
    call MPI_ALLREDUCE(MPI_IN_PLACE,row_indices,peinf%npes*max_npr, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,col_indices,peinf%npes*max_npc, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)

    !XXXXXXXXXX Here we allocate the full matrix on the first inode=0 (can be problematic in term of memory...)
    IF(keep_full_eps_static) THEN
      ! exchange info (low memory communication scheme)
      SAFE_ALLOCATE(pool2pe, (peinf%npools,peinf%npes_pool))
      pool2pe = 0
      IF(peinf%my_pool >= 0 .AND. peinf%pool_rank >= 0) THEN
        pool2pe(peinf%my_pool+1, peinf%pool_rank+1) = peinf%inode
      END IF
      call MPI_ALLREDUCE(MPI_IN_PLACE, pool2pe, peinf%npools*peinf%npes_pool, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
      !
      ! go with h5
#ifdef HDF5
      if (sig%use_hdf5) then
        ! read from .h5 
        ! here we try to read only for a single pool of process using a 1d-cyclic 
        ! scalapack layout
        scal_1d%nprow   = 1
        scal_1d%npcol   = peinf%npes_pool
        scal_1d%myprow  = 0
        scal_1d%mypcol  = peinf%pool_rank
        scal_1d%nbl     = epsmpi%nb
        scal_1d%icntxt  = scal%icntxt
        nrow_loc_1d = numroc(nmtx, scal_1d%nbl, scal_1d%myprow, 0, scal_1d%nprow)
        ncol_loc_1d = numroc(nmtx, scal_1d%nbl, scal_1d%mypcol, 0, scal_1d%npcol)
        scal_1d%npr     = nrow_loc_1d
        scal_1d%npc     = ncol_loc_1d
        ! write(*,*) peinf%inode, peinf%pool_rank, peinf%my_pool, nrow_loc_1d, ncol_loc_1d, epsmpi%ngpown

        ! temporary buffer
        SAFE_ALLOCATE(eps1Aux, (scal_1d%npr, scal_1d%npc))
        eps1Aux = (0.0d0, 0.0d0)
        if( peinf%my_pool == 0) then
          !
          call timing%start(timing%epscopy_io)
          call read_eps_matrix_par_hdf5_general(scal_1d, eps1Aux, nmtx, nmtx, iq, 1, fname, &
                                                'mats/matrix_fulleps0', peinf%pool_comm)
          call timing%stop(timing%epscopy_io)
          !
          call timing%start(timing%epscopy_comm)
          do iproc_dum = 1, peinf%npools - 1
            iproc_send = pool2pe(iproc_dum+1, peinf%pool_rank+1)
            tag = 0
            CALL MPI_SEND(eps1Aux, scal_1d%npr*scal_1d%npc, &
                          MPI_COMPLEX_DPC, iproc_send, tag, MPI_COMM_WORLD, mpierr)
          end do
          call timing%stop(timing%epscopy_comm)
          ! nrow_loc_1d should be equal to nmtx
          epsmpi%epsR(1:nrow_loc_1d,1:ncol_loc_1d,1,iq+qoffset) = eps1Aux(:,:)
        else
          iproc_dum = pool2pe(1, peinf%pool_rank+1)
          ! this condition should ensure that the processes left over will not be 
          ! involved in communication
          if(peinf%my_pool >= 0 .AND. peinf%pool_rank >= 0) then
            tag = 0
            CALL MPI_RECV(eps1Aux, scal_1d%npr*scal_1d%npc, &
                          MPI_COMPLEX_DPC, iproc_dum, tag, MPI_COMM_WORLD, mpistatus, mpierr)
          end if
          epsmpi%epsR(1:nrow_loc_1d,1:ncol_loc_1d,1,iq+qoffset) = eps1Aux(:,:)
        end if
        SAFE_DEALLOCATE(eps1Aux)

      else  ! sig%use_hdf5
#endif
        ! standard binary read
        !
        IF(peinf%inode == 0) THEN
          DO j_global = 1, nmtx
            ! read full column
            read(iunit) (C_send(ig), ig=1, nmtx)
            ! 
            ! global to pool_rank (in the new sigma distribution)
            pool_send = INDXG2P( j_global, epsmpi%nb, iproc_dum, 0, peinf%npes_pool)
            DO iproc_dum = 0, peinf%npools - 1
              ! select the process that I will have to send
              iproc_send = pool2pe(iproc_dum+1,pool_send+1)
              IF(iproc_send == 0) THEN
                j_local = INDXG2L( j_global, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
                epsmpi%epsR(1:nmtx,j_local,1,iq+qoffset) = C_send(1:nmtx)
              ELSE
                tag = 0
                CALL MPI_SEND(C_send, nmtx, &
                              MPI_COMPLEX_DPC, iproc_send, tag, MPI_COMM_WORLD, mpierr)
              END IF
            END DO
          END DO
        ELSE
          DO j_global = 1, nmtx
            ! global to pool_rank (in the new sigma distribution)
            pool_send = INDXG2P( j_global, epsmpi%nb, iproc_dum, 0, peinf%npes_pool)
            DO iproc_dum = 0, peinf%npools - 1
              ! select the process that I will have to send
              iproc_send = pool2pe(iproc_dum+1,pool_send+1)
              IF(iproc_send == peinf%inode) THEN
                j_local = INDXG2L( j_global, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
                tag = 0
                CALL MPI_RECV(C_recv, nmtx, &
                              MPI_COMPLEX_DPC, 0, tag, MPI_COMM_WORLD, mpistatus, mpierr)
                epsmpi%epsR(1:nmtx,j_local,1,iq+qoffset) = C_recv
              END IF
            END DO
          END DO
        END IF
        ! 
#ifdef HDF5
      end if ! h5
#endif
      !
      SAFE_DEALLOCATE(pool2pe)
    END IF ! keep_full_eps_static
 
    SAFE_ALLOCATE(eigenvect, (scal%npr, scal%npc))
    eigenvect = (0.0d0,0.0d0)

#ifdef HDF5
    if (sig%use_hdf5) then
      ! read from .h5
      ! get number of eigenvectors
      neig_sub = 0
      IF(peinf%inode==0) THEN
        neig_sub = sub_neig_of_q(iq)
      END IF
      call MPI_Bcast(neig_sub,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

      scal_aux%nprow   = scal%nprow
      scal_aux%npcol   = scal%npcol
      scal_aux%myprow  = scal%myprow
      scal_aux%mypcol  = scal%mypcol
      scal_aux%nbl     = scal%nbl
      scal_aux%icntxt  = scal%icntxt
      scal_aux%npr     = scal%npr
      scal_aux%npc     = numroc(neig_sub, scal%nbl, scal%mypcol, 0, scal%npcol)
      !
      call timing%start(timing%epscopy_io)
      if(.true.) then
        call read_eps_matrix_par_hdf5_general(scal_aux, eigenvect, nmtx, neig_sub, iq, 1, fname, &
                                              'mats/matrix_eigenvec', MPI_COMM_WORLD)
      else
        call read_eps_matrix_par_hdf5_general(scal, eigenvect, nmtx, nmtx, iq, 1, fname, &
                                              'mats/matrix_eigenvec', MPI_COMM_WORLD)
      end if
      call timing%stop(timing%epscopy_io)
      !
    else  ! sig%use_hdf5
#endif
      ! standard binary read
      !
      neig_sub = 0
      IF(peinf%inode==0) THEN
        read(iunit) neig_sub
      END IF
      call MPI_Bcast(neig_sub,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

      SAFE_ALLOCATE(buffer, (max_npr,max_npc))
      buffer = (0.0d0,0.0d0)

      ! read eigenvector and redistribute
      IF(peinf%inode==0) THEN
        call timing%start(timing%sub_io_vec)
        SAFE_ALLOCATE(C_aux, (nmtx,neig_sub))
        C_aux = (0.0d0,0.0d0)
        do igp=1, neig_sub
          read(iunit) (C_aux(ig,igp),ig=1,nmtx)
        end do
        do igp=neig_sub + 1, nmtx
          read(iunit) 
        end do
        call timing%stop(timing%sub_io_vec)
        ! prepare the messages and send them
        DO ipe = 1, peinf%npes - 1
          call timing%start(timing%sub_prep_vec)
          buffer = (0.0d0,0.0d0)
          DO j_local=1, scal%npcd(ipe + 1)
            j_global = col_indices(j_local, ipe + 1)
            IF(j_global > neig_sub) CYCLE
            DO i_local=1, scal%nprd(ipe + 1)
              i_global = row_indices(i_local, ipe + 1)
              buffer(i_local, j_local) = C_aux(i_global,j_global)
            END DO
          END DO
          call timing%stop(timing%sub_prep_vec)
          ! send
          call timing%start(timing%sub_comm_vec)
          tag = 0
          CALL MPI_SEND(buffer, max_npr*max_npc, MPI_COMPLEX_DPC, ipe, tag, MPI_COMM_WORLD,mpierr)
          call timing%stop(timing%sub_comm_vec)
        END DO
        ! myself
        DO j_local=1, scal%npc
          j_global = col_indices(j_local, peinf%inode + 1)
          IF(j_global > neig_sub) CYCLE
          DO i_local=1, scal%npr
            i_global = row_indices(i_local, peinf%inode + 1)
            eigenvect(i_local, j_local) = C_aux(i_global,j_global)
          END DO
        END DO
        ! done

        !XXX if(peinf%inode == 0) then
        !XXX   do ig = 1, neig_sub
        !XXX     do igp = 1, nmtx
        !XXX       write(2000,*) ig, igp, C_aux(igp,ig)
        !XXX     end do
        !XXX   end do
        !XXX   write(2000,*)
        !XXX end if

        SAFE_DEALLOCATE(C_aux)
      ELSE
        ! receive my stuff
        tag = 0
        CALL MPI_RECV(buffer, max_npr*max_npc, MPI_COMPLEX_DPC, 0, tag, MPI_COMM_WORLD, mpistatus, mpierr)
        eigenvect(1:scal%npr, 1:scal%npc) = buffer(1:scal%npr, 1:scal%npc)
      END IF
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      SAFE_DEALLOCATE(buffer)
      !  
#ifdef HDF5
    end if ! h5
#endif

      !XXX do j_local = 1, scal%npc
      !XXX   do i_local = 1, scal%npr
      !XXX     write(1000+peinf%inode,*) i_local, j_local, eigenvect(i_local,j_local)
      !XXX   end do
      !XXX end do
      !XXX flush(1000+peinf%inode) 
      !XXX write(*,*) 'AAAAA'

    ! All the same for the subspace epsilon (all frequencies in one)
    ! First create the scalapack env, retain the block/grid-layout as that 
    ! of the eigenvalue matrix
    ! get the row/col size for the subspace matrices
    ! call MPI_Bcast(neig_sub,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    nrow_loc_sub = numroc(neig_sub, scal%nbl, scal%myprow, 0, scal%nprow)
    ncol_loc_sub = numroc(neig_sub, scal%nbl, scal%mypcol, 0, scal%npcol)
    scal_sub%nprow   = scal%nprow
    scal_sub%npcol   = scal%npcol
    scal_sub%myprow  = scal%myprow
    scal_sub%mypcol  = scal%mypcol
    scal_sub%nbl     = scal%nbl
    scal_sub%icntxt  = scal%icntxt
    scal_sub%npr     = nrow_loc_sub
    scal_sub%npc     = ncol_loc_sub
    ! 
    SAFE_ALLOCATE(scal_sub%nprd, (peinf%npes))
    SAFE_ALLOCATE(scal_sub%npcd, (peinf%npes))
    scal_sub%nprd = 0
    scal_sub%npcd = 0
    scal_sub%nprd(peinf%inode + 1) = scal_sub%npr
    scal_sub%npcd(peinf%inode + 1) = scal_sub%npc
    !
    call MPI_ALLREDUCE(MPI_IN_PLACE,scal_sub%nprd,peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,scal_sub%npcd,peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)
    ! get max number of row/col
    max_npr_sub = MAXVAL(scal_sub%nprd)
    max_npc_sub = MAXVAL(scal_sub%npcd)

    ! precompute indices for subspace matrices
    SAFE_ALLOCATE(row_indices_sub, (max_npr_sub, peinf%npes))
    row_indices_sub = 0
    DO i_local = 1, scal_sub%npr
      i_global = indxl2g(i_local,  scal_sub%nbl, scal_sub%myprow, 0, scal_sub%nprow)
      row_indices_sub(i_local, peinf%inode + 1) = i_global
    END DO
    SAFE_ALLOCATE(col_indices_sub, (max_npc_sub, peinf%npes))
    col_indices_sub = 0
    DO i_local = 1, scal_sub%npc
      i_global = indxl2g(i_local,  scal_sub%nbl, scal_sub%mypcol, 0, scal_sub%npcol)
      col_indices_sub(i_local, peinf%inode + 1) = i_global
    END DO
    call MPI_ALLREDUCE(MPI_IN_PLACE, row_indices_sub,peinf%npes*max_npr_sub, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, col_indices_sub,peinf%npes*max_npc_sub, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,mpierr)

    ! here we go with copy/redistribution
    Nfreq = sig%Nfreq
    SAFE_ALLOCATE(eps_sub_f, (MAX(1,scal_sub%npr), MAX(1,scal_sub%npc), Nfreq))
    eps_sub_f = (0.0d0,0.0d0)
#ifdef HDF5
    if (sig%use_hdf5) then
      ! read from .h5
      call timing%start(timing%epscopy_io)
      call read_eps_matrix_par_hdf5_general_f(scal_sub, eps_sub_f, neig_sub, neig_sub, Nfreq, iq, 1, fname, &
                                              'mats/matrix_subspace', MPI_COMM_WORLD)
      call timing%stop(timing%epscopy_io)
      !
    else  ! sig%use_hdf5
#endif
      ! binary case
      SAFE_ALLOCATE(buffer_f, (Nfreq, max_npr_sub, max_npc_sub))
      buffer_f  = (0.0d0,0.0d0)
      IF(peinf%inode==0) THEN
        call timing%start(timing%sub_io_eps)
        SAFE_ALLOCATE(E_aux, (Nfreq,neig_sub,neig_sub))
        E_aux = (0.0d0,0.0d0)
        do igp= 1, neig_sub
          do ig = 1, neig_sub
            read(iunit) (E_aux(ifreq,ig,igp),ifreq=1,Nfreq)
          end do
#ifdef CPLX
          do ig=1, neig_sub
            read(iunit) ! For epsADyn (empty line...)
          enddo
#endif
        end do
        call timing%stop(timing%sub_io_eps)
        ! prepare the messages and send them
        DO ipe = 1, peinf%npes - 1
          call timing%start(timing%sub_prep_eps)
          buffer_f = (0.0d0,0.0d0)
          DO j_local=1, scal_sub%npcd(ipe + 1)
            j_global = col_indices_sub(j_local, ipe + 1)
            DO i_local=1, scal_sub%nprd(ipe + 1)
              i_global = row_indices_sub(i_local, ipe + 1)
              buffer_f(1:Nfreq, i_local, j_local) = E_aux(1:Nfreq, i_global, j_global)
            END DO
          END DO
          call timing%stop(timing%sub_prep_eps)
          ! send
          call timing%start(timing%sub_comm_eps)
          tag = 0
          CALL MPI_SEND(buffer_f, Nfreq*max_npr_sub*max_npc_sub, MPI_COMPLEX_DPC, ipe, tag, MPI_COMM_WORLD, mpierr)
          call timing%stop(timing%sub_comm_eps)
        END DO
        ! myself
        DO j_local=1, scal_sub%npc
          j_global = col_indices_sub(j_local, peinf%inode + 1)
          DO i_local=1, scal_sub%npr
            i_global = row_indices_sub(i_local, peinf%inode + 1)
            eps_sub_f(i_local, j_local, 1:Nfreq) = E_aux(1:Nfreq, i_global, j_global)
          END DO
        END DO
        ! done
        SAFE_DEALLOCATE(E_aux)
      ELSE
        ! receive my stuff
        tag = 0
        CALL MPI_RECV(buffer_f, max_npr_sub*max_npc_sub*Nfreq, MPI_COMPLEX_DPC, 0, tag, MPI_COMM_WORLD, mpistatus, mpierr)
        DO ifreq = 1, Nfreq
          eps_sub_f(1:scal_sub%npr, 1:scal_sub%npc, ifreq) = buffer_f(ifreq, 1:scal_sub%npr, 1:scal_sub%npc)
        END DO
      END IF
      SAFE_DEALLOCATE(buffer_f)
      !
#ifdef HDF5
    end if ! sig%use_hdf5
#endif
 
    ! initialize descriptors
    call descinit(desca, nmtx, nmtx, scal%nbl, scal%nbl, 0, 0, &
                  scal%icntxt, max(1,scal%npr), info)
    if(info < 0) then
      write(0,'(a,i3,a)') 'Argument number ', -info, ' had an illegal value on entry.'
      call die("descinit error for descaA in subspace epscopy")
    else if(info > 0) then
      write(0,*) 'info = ', info
      call die("descinit error for descaA in subspace epscopy")
    endif

    call descinit(desc_sub, neig_sub, neig_sub, scal_sub%nbl, scal_sub%nbl, 0, 0, &
                  scal_sub%icntxt, max(1,nrow_loc_sub), info)
    if(info < 0) then
      write(0,'(a,i3,a)') 'Argument number ', -info, ' had an illegal value on entry.'
      call die("descinit error for desca_sub in subspace epscopy")
    else if(info > 0) then
      write(0,*) 'info = ', info
      call die("descinit error for desca_sub in subspace epscopy")
    endif

    ! 0) initialize, map proc to rank_pool
    SAFE_ALLOCATE(grid2pe, (scal%nprow, scal%npcol))
    grid2pe = 0
    grid2pe(scal%myprow+1, scal%mypcol+1) = peinf%inode
    call MPI_ALLREDUCE(MPI_IN_PLACE, grid2pe, scal%nprow*scal%npcol , MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
    !
    SAFE_ALLOCATE(pe2grid, (2,peinf%npes))
    pe2grid = 0
    pe2grid(1, peinf%inode+1) = scal%myprow
    pe2grid(2, peinf%inode+1) = scal%mypcol
    call MPI_ALLREDUCE(MPI_IN_PLACE, pe2grid, 2*peinf%npes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
    !
    SAFE_ALLOCATE(pool2pe, (peinf%npools,peinf%npes_pool))
    pool2pe = 0
    IF(peinf%my_pool >= 0 .AND. peinf%pool_rank >= 0) THEN
      pool2pe(peinf%my_pool+1, peinf%pool_rank+1) = peinf%inode
    END IF
    call MPI_ALLREDUCE(MPI_IN_PLACE, pool2pe, peinf%npools*peinf%npes_pool, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
    ! loop over my local column and calculate the how many of them I have to send to all other proc
    SAFE_ALLOCATE(num_col_to_send, (peinf%npes))
    num_col_to_send = 0
    do j_local = 1, scal%npc
      ! local to global (in the epsilon distribution)
      j_global = col_indices(j_local, peinf%inode + 1) 
      IF(j_global > nmtx) CYCLE
      ! global to pool_rank (in the new sigma distribution)
      pool_send = INDXG2P( j_global, epsmpi%nb, iproc_dum, 0, peinf%npes_pool)
      DO iproc_dum = 0, peinf%npools - 1
        ! select the process that I will have to send
        iproc_send = pool2pe(iproc_dum+1,pool_send+1)
        num_col_to_send(iproc_send+1) = num_col_to_send(iproc_send+1) + 1
      END DO
    end do
    ! loop over my local column (in the new sigma distribution) and precompute how many columns I have to receive from
    ! all other processes, this will happen only if peinf%my_pool >= 0 (I hope...)
    SAFE_ALLOCATE(num_col_to_recv, (peinf%npes))
    num_col_to_recv = 0
    IF(peinf%my_pool >= 0) THEN
      do j_local = 1, epsmpi%ngpown
        ! local to global (in the sigma distribution)
        j_global = INDXL2G(j_local, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
        IF(j_global > nmtx) CYCLE
        ! get the process in the epsilon distribution that hold this column
        iproc_row_rec = INDXG2P( j_global, scal%nbl, iproc_dum, 0, scal%npcol)
        DO iproc_dum = 0, scal%nprow - 1
          iproc_rec = grid2pe(iproc_dum+1, iproc_row_rec+1)
          num_col_to_recv(iproc_rec+1) = num_col_to_recv(iproc_rec+1) + 1
        END DO
      end do
    END IF
    
    ! allocate receiving buffer
    Nbuf_recv = 0
    DO ipe = 1, peinf%npes
      IF(num_col_to_recv(ipe) > 0) Nbuf_recv = Nbuf_recv + 1 
    END DO
    SAFE_ALLOCATE(buf_recv, (Nbuf_recv))
    imsg_recv = 0
    DO ipe = 1, peinf%npes
      IF(num_col_to_recv(ipe) > 0) THEN
        imsg_recv = imsg_recv + 1
        buf_recv(imsg_recv)%proc = ipe - 1
        buf_recv(imsg_recv)%ncol = num_col_to_recv(ipe)
        buf_recv(imsg_recv)%nrow = scal%nprd(ipe)
        SAFE_ALLOCATE(buf_recv(imsg_recv)%msg, (buf_recv(imsg_recv)%nrow,buf_recv(imsg_recv)%ncol))
        SAFE_ALLOCATE(buf_recv(imsg_recv)%col_global_indx, (buf_recv(imsg_recv)%ncol))
        buf_recv(imsg_recv)%col_global_indx = -100
      END IF
    END DO

    ! allocate sending buffer
    Nbuf_send = 0
    DO ipe = 1, peinf%npes
      IF(num_col_to_send(ipe) > 0) Nbuf_send = Nbuf_send + 1
    END DO
    SAFE_ALLOCATE(proc2msg_send, (peinf%npes))
    proc2msg_send = 0
    SAFE_ALLOCATE(buf_send, (Nbuf_send))
    imsg_send = 0
    DO ipe = 1, peinf%npes
      IF(num_col_to_send(ipe) > 0) THEN
        imsg_send = imsg_send + 1
        buf_send(imsg_send)%proc = ipe - 1
        buf_send(imsg_send)%ncol = num_col_to_send(ipe)
        buf_send(imsg_send)%nrow = scal%npr
        SAFE_ALLOCATE(buf_send(imsg_send)%msg, (buf_send(imsg_send)%nrow, buf_send(imsg_send)%ncol))
        SAFE_ALLOCATE(buf_send(imsg_send)%col_global_indx, (buf_send(imsg_send)%ncol))
        buf_send(imsg_send)%col_global_indx = -100
        ! 
        proc2msg_send(ipe) = imsg_send
      END IF
    END DO
    SAFE_ALLOCATE(local_col_count, (Nbuf_send))
    SAFE_ALLOCATE(req_send, (2 * Nbuf_send))
    SAFE_ALLOCATE(my_status, (MPI_STATUS_SIZE, 2*Nbuf_send))
    SAFE_ALLOCATE(req_recv, (2 * Nbuf_recv))

    ! in the case of sigma subspace find which is the head/wings position for this 
    ! specific q-point
    if(sig%do_sigma_subspace) then
      max_head = 1
      if((peinf%inode==0)) then
        gg = 0
        call findvector(iout, gg, gvec)
        max_head = epsmpi%isrtqi(iout,iq+qoffset)
        call MPI_Bcast(max_head,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      else
        call MPI_Bcast(max_head,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      end if
      sig%epssub%wing_pos(iq+qoffset) = max_head
    end if

    ! prepare the coul_Aux to avoid recomputing many times the SQRT(vcoul)
    call timing%start(timing%epscopy_vcoul)
    SAFE_ALLOCATE(row_vcoul_aux, (scal%npr)) ! row -> * SQRT(vc)
    SAFE_ALLOCATE(col_vcoul_aux, (scal%npc)) ! col -> * 1/SQRT(vc)
    ! initialize to one since it will be a scaling factor
    row_vcoul_aux = 1.0d+00
    do ig_l = 1, scal%npr
      ig_g = indxl2g(ig_l, scal%nbl, scal%myprow, 0, scal%nprow)
      row_vcoul_aux(ig_l) = sqrt(vcoul(ig_g))
    enddo
    col_vcoul_aux = 1.0d+00
    do igp_l = 1, scal%npc
      igp_g = indxl2g(igp_l, scal%nbl, scal%mypcol, 0, scal%npcol)
      vc = 1.0d+00 / sqrt(vcoul(igp_g))
      col_vcoul_aux(igp_l) = vc
    end do
    ! diagonal elements
    my_num_diag_elem = 0
    do igp_l = 1, scal%npc
      igp_g = indxl2g(igp_l, scal%nbl, scal%mypcol, 0, scal%npcol)
      do ig_l = 1, scal%npr
        ig_g = indxl2g(ig_l, scal%nbl, scal%myprow, 0, scal%nprow)
        if (ig_g==igp_g) my_num_diag_elem = my_num_diag_elem + 1
      enddo
    enddo
    ! allocate index array for diagonal elements
    SAFE_ALLOCATE(my_diag_indeces, (2,my_num_diag_elem))
    my_diag_indeces = 0
    idiag = 0
    do igp_l = 1, scal%npc
      igp_g = indxl2g(igp_l, scal%nbl, scal%mypcol, 0, scal%npcol)
      vc = sqrt(vcoul(igp_g))
      do ig_l = 1, scal%npr
        ig_g = indxl2g(ig_l, scal%nbl, scal%myprow, 0, scal%nprow)
        if (ig_g==igp_g) then
          idiag = idiag + 1
          my_diag_indeces(1, idiag) = ig_l
          my_diag_indeces(2, idiag) = igp_l
        end if
      enddo
    enddo
    call timing%stop(timing%epscopy_vcoul)

    ! for the subspace epsinv ONE has been already scaled from the diagonal...
    ! expand the subspace inverse epsilon 
    ! allocate temporary matrices
    SAFE_ALLOCATE(eps1Aux, (scal%npr, scal%npc))
    SAFE_ALLOCATE(C_Pgemm, (scal%npr, scal%npc))
    ! loop over frequencies
    DO ifreq = 1, Nfreq
      call timing%start(timing%epscopy_pgemm)
      CALL pzgemm('N','N', nmtx, neig_sub, neig_sub, (1.0d0,0.0d0), eigenvect, 1, 1, desca, &
                  eps_sub_f(:,:,ifreq), 1, 1, desc_sub, (0.0d0,0.0d0), &
                  C_Pgemm, 1, 1, desca)
      if((.not.sig%do_sigma_subspace) .or. ifreq == 1) then
        ! here we always calculate the full epsilon^{-1} for the first frequency, this
        ! is apparently important if we need to calculate the static CH/SH
        CALL pzgemm('N','C', nmtx, nmtx, neig_sub, (1.0d0,0.0d0), C_Pgemm, 1, 1, desca, &
                    eigenvect, 1, 1, desca, (0.0d0,0.0d0), &
                    eps1Aux, 1, 1, desca)
      else
         CALL pzgemm('N','C', max_head, nmtx, neig_sub, (1.0d0,0.0d0), C_Pgemm, 1, 1, desca, &
                     eigenvect, 1, 1, desca, (0.0d0,0.0d0), &
                     eps1Aux, 1, 1, desca)
         do j_local = 1, scal%npc
           j_global = col_indices(j_local, peinf%inode + 1)
           if( j_global <= max_head) eps1Aux(:,j_local) = (0.0d0,0.0d0)
         end do
         CALL pzgemm('N','C', nmtx, max_head, neig_sub, (1.0d0,0.0d0), C_Pgemm, 1, 1, desca, &
                     eigenvect, 1, 1, desca, (1.0d0,0.0d0), &
                     eps1Aux, 1, 1, desca)
      end if
      call timing%stop(timing%epscopy_pgemm)
      ! here we have the full matrix for frequency ifreq, restore one on diagonal
      ! and "unsymmetrize" with vcoul
      ! FHJ: WARNING: never perform a nested loop over the global rows and
      ! columns, as these dimensions may be huge, and the code will spend
      ! a long time doing nothing. Instead, always loop over the local rows
      ! and columns and use indxl2g to get the corresponding global index.
      call timing%start(timing%epscopy_vcoul)
      if(.true.) then
        ! new fast scheme
        ! sum one on diagonal
        do idiag = 1, my_num_diag_elem
          i_local = my_diag_indeces(1, idiag)
          j_local = my_diag_indeces(2, idiag)
          eps1Aux(i_local, j_local) = eps1Aux(i_local, j_local) + ONE
        end do
        ! scale with coulomb operator
        do igp_l = 1, scal%npc
          vc = col_vcoul_aux(igp_l)
          do ig_l = 1, scal%npr
            eps1Aux(ig_l,igp_l) = eps1Aux(ig_l,igp_l) * row_vcoul_aux(ig_l) * vc
          end do
        end do
        !
      else
        ! keep old slow scheme for debug
        do igp_l = 1, scal%npc
          igp_g = indxl2g(igp_l, scal%nbl, scal%mypcol, 0, scal%npcol)
          vc = sqrt(vcoul(igp_g))
          do ig_l = 1, scal%npr
            ig_g = indxl2g(ig_l, scal%nbl, scal%myprow, 0, scal%nprow)
            if (ig_g==igp_g) eps1Aux(ig_l,igp_l) = eps1Aux(ig_l,igp_l) + ONE
            eps1Aux(ig_l,igp_l) = eps1Aux(ig_l,igp_l) * sqrt(vcoul(ig_g)) / vc
          enddo
        enddo
        !
      end if
      call timing%stop(timing%epscopy_vcoul)

      ! in case of subspace in sigma, we just save in the fully expanded 
      ! form the static case (first frequency, equal to 0)
      IF((.not.sig%do_sigma_subspace) .OR. ifreq == 1) THEN
        ! here we have the epsinv matrix as output from a standard epsilon calculation
        ! redistribute according to the new layout as required in sigma (for omega=0 there
        ! is the possibility to keep the full epsinv)
        IF((.NOT.(keep_full_eps_static)) .OR. ifreq > 1) THEN
          ! post all messages to be received
          DO imsg_recv = 1, Nbuf_recv
            ! WRITE(2000+peinf%inode,*) peinf%inode, imsg_recv, buf_recv(imsg_recv)%proc
            tag = 0
            CALL MPI_IRECV(buf_recv(imsg_recv)%msg, buf_recv(imsg_recv)%ncol*buf_recv(imsg_recv)%nrow, &
                           MPI_COMPLEX_DPC, buf_recv(imsg_recv)%proc, tag, MPI_COMM_WORLD, req_recv((imsg_recv-1)*2+1), &
                           mpierr)
            tag = 0
            CALL MPI_IRECV(buf_recv(imsg_recv)%col_global_indx, buf_recv(imsg_recv)%ncol, &
                           MPI_INTEGER, buf_recv(imsg_recv)%proc, tag, MPI_COMM_WORLD, req_recv(imsg_recv*2), &
                           mpierr)
          END DO
          ! fill send buffer
          ! loop over local column
          local_col_count = 0
          do j_local = 1, scal%npc
            j_global = col_indices(j_local, peinf%inode + 1)
            IF(j_global > nmtx) CYCLE
            pool_send = INDXG2P( j_global, epsmpi%nb, iproc_dum, 0, peinf%npes_pool)
            DO iproc_dum = 0, peinf%npools - 1
              iproc_send = pool2pe(iproc_dum+1,pool_send+1)
              imsg_send = proc2msg_send(iproc_send+1)
              local_col_count(imsg_send) = local_col_count(imsg_send) + 1
              buf_send(imsg_send)%msg(1:scal%npr,local_col_count(imsg_send)) = eps1Aux(1:scal%npr,j_local)
              buf_send(imsg_send)%col_global_indx(local_col_count(imsg_send)) = j_global
            END DO
          end do
          ! send messeges
          DO imsg_send = 1, Nbuf_send
            ! WRITE(3000+peinf%inode,*) peinf%inode, imsg_send, buf_send(imsg_send)%proc
            tag = 0
            CALL MPI_ISEND(buf_send(imsg_send)%msg, buf_send(imsg_send)%ncol*buf_send(imsg_send)%nrow, &
                           MPI_COMPLEX_DPC, buf_send(imsg_send)%proc, tag, MPI_COMM_WORLD, req_send((imsg_send-1)*2+1), &
                           mpierr)
            tag = 0
            CALL MPI_ISEND(buf_send(imsg_send)%col_global_indx, buf_send(imsg_send)%ncol, &
                           MPI_INTEGER, buf_send(imsg_send)%proc, tag, MPI_COMM_WORLD, req_send(imsg_send*2), &
                           mpierr)
          END DO
          ! collect messages
          DO imsg_recv = 1, Nbuf_recv
            CALL MPI_WAIT(req_recv((imsg_recv-1)*2+1), mpistatus, mpierr)
            CALL MPI_WAIT(req_recv(imsg_recv*2), mpistatus, mpierr)
            ! fill in epsmpi
            ! iproc_row_rec = pe2grid(1, buf_recv(imsg_recv)%proc+1)
            DO j_local = 1, buf_recv(imsg_recv)%ncol
              j_global = buf_recv(imsg_recv)%col_global_indx(j_local)
              my_col = INDXG2L(j_global, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
              DO i_local = 1, buf_recv(imsg_recv)%nrow
                i_global = row_indices(i_local, buf_recv(imsg_recv)%proc+1)
                epsmpi%epsR(i_global,my_col,ifreq,iq+qoffset) = buf_recv(imsg_recv)%msg(i_local,j_local)
              END DO
            END DO
          END DO

          ! wait for all
          CALL MPI_WAITALL(Nbuf_send*2,req_send,my_status,mpierr)
        END IF 
        call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      END IF ! .not.subspace and ifreq>1

      IF(sig%do_sigma_subspace) THEN
        ! in case of subspace here we need to save the wings of the full matrix
        ! for simplicity we always read from eps1Aux, even in the case (keep_full_eps_static and ifreq==1)

        sig%epssub%eps_wings_cols(:,ifreq,iq+qoffset) = (0.0d0,0.0d0)
        ipe_wing = INDXG2P(max_head, scal%nbl, iproc_dum, 0, scal%npcol)
        if (ipe_wing == scal%mypcol) then
          my_pos_loc_wing = INDXG2L(max_head, scal%nbl, scal%mypcol, 0, scal%npcol)
          do i_local = 1,  scal%npr
            i_global = row_indices(i_local, peinf%inode + 1)
            sig%epssub%eps_wings_cols(i_global,ifreq,iq+qoffset) = eps1Aux(i_local,my_pos_loc_wing)
          end do
        end if

        sig%epssub%eps_wings_rows(:,ifreq,iq+qoffset) = (0.0d0,0.0d0)
        ipe_wing = INDXG2P(max_head, scal%nbl, iproc_dum, 0, scal%nprow)
        if (ipe_wing == scal%myprow) then
          my_pos_loc_wing = INDXG2L(max_head, scal%nbl, scal%myprow, 0, scal%nprow)
          do j_local = 1, scal%npc
            j_global = col_indices(j_local, peinf%inode + 1)
            sig%epssub%eps_wings_rows(j_global,ifreq,iq+qoffset) = eps1Aux(my_pos_loc_wing,j_local)
            if (j_global == max_head) then
              sig%epssub%eps_wings_rows(max_head,ifreq,iq+qoffset) = &
              sig%epssub%eps_wings_rows(max_head,ifreq,iq+qoffset) - 1.0d0
            end if
          end do
        end if

      END IF

    END DO ! ifreq

    SAFE_DEALLOCATE(row_vcoul_aux)
    SAFE_DEALLOCATE(col_vcoul_aux)
    SAFE_DEALLOCATE(my_diag_indeces)

    ! check if we have to redistribute eigenvectors for sigma-subspace calc
    if(sig%do_sigma_subspace) then
      call timing%start(timing%epscopy_redstr)
      !
      ! sum-up wings from previous step
      call MPI_Allreduce(MPI_IN_PLACE, sig%epssub%eps_wings_rows(1:Neps,1:sig%nFreq,iq+qoffset), &
                         sig%epssub%neps * sig%nFreq, MPI_COMPLEX_DPC, MPI_SUM, MPI_COMM_WORLD, mpierr)
      call MPI_Allreduce(MPI_IN_PLACE, sig%epssub%eps_wings_cols(1:Neps,1:sig%nFreq,iq+qoffset), &
                         sig%epssub%neps * sig%nFreq, MPI_COMPLEX_DPC, MPI_SUM, MPI_COMM_WORLD, mpierr)

      if ( cyclic_redistr ) then
        ! cyclic style redistribution (not the most efficient but easier than other methods)
        ! start allocating buffers and initialize
        ! eigenvectors
        ! figure out how many col do this porc own with global index less than Neig
        ! (assuming that if j_local < j_local' then j_global < j_global')
        npc_lt_Neig = 0
        do j_local = 1, scal%npc
          j_global = col_indices(j_local, peinf%inode + 1)
          if ( j_global > neig_sub ) exit
          npc_lt_Neig = npc_lt_Neig + 1
        end do
        max_npc_lt_Neig = npc_lt_Neig
        call mpi_allreduce(MPI_IN_PLACE, max_npc_lt_Neig, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
        !
        SAFE_ALLOCATE(buffer_sub_eigen,      (max_npr, max_npc_lt_Neig))
        SAFE_ALLOCATE(rec_buffer_sub_eigen,  (max_npr, max_npc_lt_Neig))
        SAFE_ALLOCATE(send_buffer_sub_eigen, (max_npr, max_npc_lt_Neig))
!$OMP PARALLEL do private (j_local, i_local) collapse(2)
        do j_local = 1, max_npc_lt_Neig
          do i_local = 1, max_npr
            buffer_sub_eigen(i_local, j_local)      = (0.0d0,0.0d0)
            rec_buffer_sub_eigen(i_local, j_local)  = (0.0d0,0.0d0)
            send_buffer_sub_eigen(i_local, j_local) = (0.0d0,0.0d0)
          end do
        end do
!$OMP END PARALLEL DO
        !
!$OMP PARALLEL do private (j_local, i_local) collapse(2)
        do j_local = 1, npc_lt_Neig
          do i_local = 1, scal%npr
            buffer_sub_eigen(i_local, j_local)      = eigenvect(i_local, j_local)
            send_buffer_sub_eigen(i_local, j_local) = eigenvect(i_local, j_local)
          end do
        end do
!$OMP END PARALLEL DO

        ! calculate my index range and save info
        do_copy = .true.
        sig%epssub%eps_sub_info(:,2,iq+qoffset) = 0
        indx_start = peinf%pool_rank * sig%epssub%Nbas_own_max + 1
        if(indx_start <= neig_sub) then
          indx_end = MIN((peinf%pool_rank + 1) * sig%epssub%Nbas_own_max, neig_sub)
          sub_size = indx_end - indx_start + 1
          sig%epssub%eps_sub_info(1,2,iq+qoffset) = indx_start
          sig%epssub%eps_sub_info(2,2,iq+qoffset) = indx_end
        else
          indx_start = 0
          indx_end   = 0
          sub_size   = 0
          do_copy = .false.
        end if
        sig%epssub%eps_sub_info(3,2,iq+qoffset) = nmtx

        ! keep track of the indeces to copy
        SAFE_ALLOCATE(copy_col_inx, (max_npc_lt_Neig))
        ! get process for communication
        isend_static = MOD(peinf%inode + 1 + peinf%npes, peinf%npes)
        irec_static  = MOD(peinf%inode - 1 + peinf%npes, peinf%npes)
        ! start cyclic loop (starting with myself)
        ipe_real = peinf%inode + 1
        do ipe = 1, peinf%npes
          !
          actual_send = MOD(peinf%inode + ipe + peinf%npes, peinf%npes)
          actual_rec  = MOD(peinf%inode - ipe + peinf%npes, peinf%npes)
          !
          ! post msgs
          tag_rec  = 1
          tag_send = 1
          req_rec_singl  = MPI_REQUEST_NULL
          req_send_singl = MPI_REQUEST_NULL
          CALL MPI_Irecv(rec_buffer_sub_eigen, max_npr*max_npc_lt_Neig, &
                         MPI_COMPLEX_DPC, irec_static, &
                         tag_rec, MPI_COMM_WORLD, req_rec_singl, mpierr)
          CALL MPI_Isend(send_buffer_sub_eigen, max_npr*max_npc_lt_Neig, &
                         MPI_COMPLEX_DPC, isend_static, &
                         tag_send, MPI_COMM_WORLD, req_send_singl, mpierr)

          if ( do_copy ) then
            ! precompute number of column to copy
            ncol_to_copy = 0
            copy_col_inx = 0
            do j_local = 1, MIN(max_npc_lt_Neig, scal%npcd(ipe_real))
              j_global = col_indices(j_local, ipe_real)
              if ( j_global >= indx_start .and. j_global <= indx_end) then
                ncol_to_copy = ncol_to_copy + 1
                copy_col_inx( ncol_to_copy ) = j_local
              end if
              if (j_global > indx_end) exit
            end do
            ! copy message on local buffer
!$OMP PARALLEL do private (j, j_local, j_global, i_local, i_global) 
            do j = 1, ncol_to_copy
              j_local  = copy_col_inx( j )
              j_global = col_indices(j_local, ipe_real) - indx_start + 1
              do i_local = 1, scal%nprd(ipe_real)
                i_global = row_indices(i_local, ipe_real)
                sig%epssub%eigenvec_sub(i_global, j_global, iq+qoffset) =  buffer_sub_eigen(i_local,j_local)
              end do
            end do
!$OMP END PARALLEL DO
            !
          end if
          
          ! finalize communication
          CALL MPI_Wait(req_rec_singl,  mpistatus, mpierr)
          ! swap buffers 
!$OMP PARALLEL do private (j_local, i_local) collapse(2)
          do j_local = 1, max_npc_lt_Neig
            do i_local = 1, max_npr
              buffer_sub_eigen(i_local, j_local) = rec_buffer_sub_eigen(i_local, j_local) 
            end do
          end do
!$OMP END PARALLEL DO
          ! the same for the sendig buffer
          CALL MPI_Wait(req_send_singl, mpistatus, mpierr)
          ! swap buffers 
!$OMP PARALLEL do private (j_local, i_local) collapse(2)
          do j_local = 1, max_npc_lt_Neig
            do i_local = 1, max_npr
              send_buffer_sub_eigen(i_local, j_local) = rec_buffer_sub_eigen(i_local, j_local) 
            end do
          end do
!$OMP END PARALLEL DO

          ! be ready for the next cycle
          ipe_real = actual_rec + 1

        end do ! ipe
        ! deallocate stuff
        SAFE_DEALLOCATE(buffer_sub_eigen)
        SAFE_DEALLOCATE(rec_buffer_sub_eigen)
        SAFE_DEALLOCATE(send_buffer_sub_eigen)
        SAFE_DEALLOCATE(copy_col_inx)

        ! the same for the subspace epsilon
        SAFE_ALLOCATE(buffer_sub_eps,      (max_npr_sub, max_npc_sub, Nfreq))
        SAFE_ALLOCATE(rec_buffer_sub_eps,  (max_npr_sub, max_npc_sub, Nfreq))
        SAFE_ALLOCATE(send_buffer_sub_eps, (max_npr_sub, max_npc_sub, Nfreq))
!$OMP PARALLEL do private (j_local, i_local, ifreq) collapse(2)
        do ifreq = 1, Nfreq
          do j_local = 1, max_npc_sub
            do i_local = 1, max_npr_sub
              buffer_sub_eps(i_local, j_local, ifreq)      = (0.0d0,0.0d0)
              rec_buffer_sub_eps(i_local, j_local, ifreq)  = (0.0d0,0.0d0)
              send_buffer_sub_eps(i_local, j_local, ifreq) = (0.0d0,0.0d0)
            end do
          end do
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL do private (j_local, i_local, ifreq) collapse(2)
        do ifreq = 1, Nfreq
          do j_local = 1, scal_sub%npc
            do i_local = 1, scal_sub%npr
              buffer_sub_eps(i_local, j_local, ifreq)      = eps_sub_f(i_local, j_local, ifreq)
              send_buffer_sub_eps(i_local, j_local, ifreq) = eps_sub_f(i_local, j_local, ifreq)
            end do
          end do
        end do
!$OMP END PARALLEL DO

        ! calculate my index range and save info
        do_copy = .true.
        sig%epssub%eps_sub_info(:,1,iq+qoffset) = 0
        indx_start = peinf%pool_rank * sig%epssub%Nbas_own_max + 1
        if(indx_start <= neig_sub) then
          indx_end = MIN((peinf%pool_rank + 1) * sig%epssub%Nbas_own_max, neig_sub)
          sub_size = indx_end - indx_start + 1
          sig%epssub%eps_sub_info(1,1,iq+qoffset) = indx_start
          sig%epssub%eps_sub_info(2,1,iq+qoffset) = indx_end
        else
          indx_start = 0
          indx_end   = 0
          sub_size   = 0
          do_copy = .false.
        end if
        sig%epssub%eps_sub_info(3,1,iq+qoffset) = neig_sub

        ! keep track of the indeces to copy
        SAFE_ALLOCATE(copy_col_inx, (max_npc_sub))
        ! get process for communication
        isend_static = MOD(peinf%inode + 1 + peinf%npes, peinf%npes)
        irec_static  = MOD(peinf%inode - 1 + peinf%npes, peinf%npes)
        ! start cyclic loop (starting with myself)
        ipe_real = peinf%inode + 1
        do ipe = 1, peinf%npes
          !
          actual_send = MOD(peinf%inode + ipe + peinf%npes, peinf%npes)
          actual_rec  = MOD(peinf%inode - ipe + peinf%npes, peinf%npes)
          !
          ! post msgs
          tag_rec  = 1
          tag_send = 1
          req_rec_singl  = MPI_REQUEST_NULL
          req_send_singl = MPI_REQUEST_NULL
          CALL MPI_Irecv(rec_buffer_sub_eps, max_npr_sub*max_npc_sub*Nfreq, &
                         MPI_COMPLEX_DPC, irec_static, &
                         tag_rec, MPI_COMM_WORLD, req_rec_singl, mpierr)
          CALL MPI_Isend(send_buffer_sub_eps, max_npr_sub*max_npc_sub*Nfreq, &
                         MPI_COMPLEX_DPC, isend_static, &
                         tag_send, MPI_COMM_WORLD, req_send_singl, mpierr)

          if ( do_copy ) then
            ! precompute number of column to copy
            ncol_to_copy = 0
            copy_col_inx = 0
            do j_local = 1, scal_sub%npcd(ipe_real)
              j_global = col_indices_sub(j_local, ipe_real)
              if ( j_global >= indx_start .and. j_global <= indx_end) then
                ncol_to_copy = ncol_to_copy + 1
                copy_col_inx( ncol_to_copy ) = j_local
              end if
              if (j_global > indx_end) exit
            end do
            ! copy message on local buffer
!$OMP PARALLEL do private (j, j_local, j_global, i_local, i_global, ifreq) 
            do ifreq = 1, Nfreq
              do j = 1, ncol_to_copy
                j_local  = copy_col_inx( j )
                j_global = col_indices_sub(j_local, ipe_real) - indx_start + 1
                do i_local = 1, scal_sub%nprd(ipe_real)
                  i_global = row_indices_sub(i_local, ipe_real)
                  sig%epssub%eps_sub(i_global, j_global, ifreq, iq+qoffset) = &
                                    buffer_sub_eps(i_local, j_local, ifreq) 
                end do
              end do
            end do
!$OMP END PARALLEL DO
          end if ! do_copy

          ! finalize communication
          CALL MPI_Wait(req_rec_singl,  mpistatus, mpierr)
          ! swap buffers 
!$OMP PARALLEL do private (j_local, i_local, ifreq) collapse(2)
          do ifreq = 1, Nfreq
            do j_local = 1, max_npc_sub
              do i_local = 1, max_npr_sub
                buffer_sub_eps(i_local, j_local, ifreq) = rec_buffer_sub_eps(i_local, j_local, ifreq)
              end do
            end do
          end do
!$OMP END PARALLEL DO
          ! the same for the sendig buffer
          CALL MPI_Wait(req_send_singl, mpistatus, mpierr)
          ! swap buffers 
!$OMP PARALLEL do private (j_local, i_local, ifreq) collapse(2)
          do ifreq = 1, Nfreq
            do j_local = 1, max_npc_sub
              do i_local = 1, max_npr_sub
                send_buffer_sub_eps(i_local, j_local, ifreq) = rec_buffer_sub_eps(i_local, j_local, ifreq)
              end do
            end do
          end do
!$OMP END PARALLEL DO

          ! be ready for the next cycle
          ipe_real = actual_rec + 1

        end do ! ipe
        !
        SAFE_DEALLOCATE(copy_col_inx)
        SAFE_DEALLOCATE(buffer_sub_eps)
        SAFE_DEALLOCATE(rec_buffer_sub_eps)
        SAFE_DEALLOCATE(send_buffer_sub_eps)
        ! DONE
      else  ! cyclic_redistr
        ! redistribution based on collective operations
        if (peinf%inode==0) write(6,'(1x,a)') 'Redistribute using collective operations!'
        !
        SAFE_ALLOCATE(buffer_sub_eps,   (neig_sub, sig%epssub%Nbas_own_max, Nfreq))
        !XXX SAFE_ALLOCATE(buffer_sub_eigen, (sig%epssub%ngpown_sub_max, neig_sub))
        SAFE_ALLOCATE(buffer_sub_eigen, (nmtx, sig%epssub%Nbas_own_max))
        do iproc_pool = 1, peinf%npes_pool

          ! eigenvectors
          indx_start = (iproc_pool-1) * sig%epssub%Nbas_own_max + 1
          if(indx_start <= neig_sub) then
            indx_end = MIN(iproc_pool * sig%epssub%Nbas_own_max, neig_sub)
            sub_size = indx_end - indx_start + 1
            
            buffer_sub_eigen = (0.0d0,0.0d0)
            icurr = 0
            do j_global = indx_start, indx_end
              icurr = icurr + 1
              iproc_col = INDXG2P( j_global, scal%nbl, iproc_dum, 0, scal%npcol)
              if(iproc_col == scal%mypcol) then
                 j_local = INDXG2L(j_global, scal%nbl, scal%mypcol, 0, scal%npcol)
                 do i_local = 1, scal%npr
                   i_global = row_indices(i_local, peinf%inode + 1)
                   if(i_global > nmtx) cycle
                   buffer_sub_eigen(i_global, icurr) = eigenvect(i_local, j_local)
                 end do
              end if
            end do

            iproc_send = pool2pe(1,iproc_pool)
            if(iproc_send == peinf%inode) then
              call MPI_Reduce(MPI_IN_PLACE, buffer_sub_eigen, &
                              nmtx*sig%epssub%Nbas_own_max, MPI_COMPLEX_DPC,&
                              MPI_SUM, iproc_send, MPI_COMM_WORLD, mpierr)
              ! copy data 
              sig%epssub%eigenvec_sub(1:nmtx, 1:sig%epssub%Nbas_own, iq+qoffset) = &
                     buffer_sub_eigen(1:nmtx, 1:sig%epssub%Nbas_own)
            else
              call MPI_Reduce(buffer_sub_eigen, buffer_sub_eigen, &
                              nmtx*sig%epssub%Nbas_own_max, MPI_COMPLEX_DPC,&
                              MPI_SUM, iproc_send, MPI_COMM_WORLD, mpierr)
            end if
          end if

          ! save information
          if(iproc_pool-1 == peinf%pool_rank) then
            sig%epssub%eps_sub_info(:,2,iq+qoffset) = 0
            if(indx_start <= neig_sub) then
              sig%epssub%eps_sub_info(1,2,iq+qoffset) = indx_start
              sig%epssub%eps_sub_info(2,2,iq+qoffset) = indx_end
            end if
            sig%epssub%eps_sub_info(3,2,iq+qoffset) = nmtx
          end if

          ! subspace epsilon
          indx_start = (iproc_pool-1) * sig%epssub%Nbas_own_max + 1
          if(indx_start <= neig_sub) then
            indx_end = MIN(iproc_pool * sig%epssub%Nbas_own_max, neig_sub)
            sub_size = indx_end - indx_start + 1

            buffer_sub_eps = (0.0d0,0.0d0)
            icurr = 0
            do j_global = indx_start, indx_end
              icurr = icurr + 1
              iproc_col = INDXG2P( j_global, scal_sub%nbl, iproc_dum, 0, scal_sub%npcol)
              if(iproc_col == scal_sub%mypcol) then
                j_local = INDXG2L(j_global, scal_sub%nbl, scal_sub%mypcol, 0, scal_sub%npcol)
                do i_local = 1, scal_sub%npr
                  i_global = row_indices_sub(i_local, peinf%inode + 1)
                  if(i_global > neig_sub) cycle
                  buffer_sub_eps(i_global, icurr, :) = eps_sub_f(i_local, j_local, :)
                end do
              end if
            end do

            iproc_send = pool2pe(1,iproc_pool)
            if(iproc_send == peinf%inode) then
              call MPI_Reduce(MPI_IN_PLACE, buffer_sub_eps, &
                              neig_sub * sig%epssub%Nbas_own_max * Nfreq, MPI_COMPLEX_DPC,&
                              MPI_SUM, iproc_send, MPI_COMM_WORLD, mpierr)
              ! copy data 
              sig%epssub%eps_sub(1:neig_sub, 1:sig%epssub%Nbas_own, 1:Nfreq, iq+qoffset) = &
              buffer_sub_eps(1:neig_sub, 1:sig%epssub%Nbas_own, 1:Nfreq)
            else
              call MPI_Reduce(buffer_sub_eps, buffer_sub_eps, &
                              neig_sub * sig%epssub%Nbas_own_max * Nfreq, MPI_COMPLEX_DPC,&
                              MPI_SUM, iproc_send, MPI_COMM_WORLD, mpierr)
            end if
          end if

          ! save information
          if(iproc_pool-1 == peinf%pool_rank) then
            sig%epssub%eps_sub_info(:,1,iq+qoffset) = 0
            if(indx_start <= neig_sub) then
              sig%epssub%eps_sub_info(1,1,iq+qoffset) = indx_start
              sig%epssub%eps_sub_info(2,1,iq+qoffset) = indx_end
            end if
            sig%epssub%eps_sub_info(3,1,iq+qoffset) = neig_sub
          end if

        end do ! iproc_pool
        SAFE_DEALLOCATE(buffer_sub_eigen)
        SAFE_DEALLOCATE(buffer_sub_eps)

        ! replicate over pools
        if(peinf%inode == pool2pe(1, peinf%pool_rank + 1)) then
          do pool_send = 2, peinf%npools
            iproc_send = pool2pe(pool_send,  peinf%pool_rank + 1)
            tag = 0
            CALL MPI_SEND(sig%epssub%eps_sub(:,:,:,iq+qoffset), &
                          sig%neig_sub_max * sig%epssub%Nbas_own * sig%nFreq, &
                          MPI_COMPLEX_DPC, iproc_send, tag, MPI_COMM_WORLD, mpierr)
            tag = 0 
            CALL MPI_SEND(sig%epssub%eigenvec_sub(:,:,iq+qoffset), &
                          sig%epssub%neps * sig%epssub%Nbas_own, &
                          MPI_COMPLEX_DPC, iproc_send, tag, MPI_COMM_WORLD, mpierr)
          end do 
        else
          if(peinf%my_pool >= 0) then
            iproc_rec = pool2pe(1, peinf%pool_rank + 1)
            tag = 0
            CALL MPI_RECV(sig%epssub%eps_sub(:,:,:,iq+qoffset), &
                          sig%neig_sub_max * sig%epssub%Nbas_own * sig%nFreq, &
                          MPI_COMPLEX_DPC, iproc_rec, tag, MPI_COMM_WORLD, mpistatus, mpierr)
            tag = 0
            CALL MPI_RECV(sig%epssub%eigenvec_sub(:,:,iq+qoffset), &
                          sig%epssub%neps * sig%epssub%Nbas_own, &
                          MPI_COMPLEX_DPC, iproc_rec, tag, MPI_COMM_WORLD, mpistatus, mpierr)
          end if
        end if
        ! 
      end if ! cyclic_redistr

      ! copy coulomb operator
      sig%epssub%vcoul_sub(1:nmtx,iq+qoffset) = vcoul(1:nmtx)
      !
      call timing%stop(timing%epscopy_redstr)
    end if  ! do_sigma_subspace

    SAFE_DEALLOCATE(eps1Aux)
    SAFE_DEALLOCATE(C_Pgemm)
    SAFE_DEALLOCATE(eigenvect)
    SAFE_DEALLOCATE(eps_sub_f)

    ! deallocate buffer
    DO imsg_recv = 1, Nbuf_recv
      SAFE_DEALLOCATE( buf_recv(imsg_recv)%msg )
      SAFE_DEALLOCATE( buf_recv(imsg_recv)%col_global_indx )
    END DO
    SAFE_DEALLOCATE(buf_recv)
    DO imsg_send = 1, Nbuf_send
      SAFE_DEALLOCATE( buf_send(imsg_send)%msg )
      SAFE_DEALLOCATE( buf_send(imsg_send)%col_global_indx )
    END DO
    SAFE_DEALLOCATE(buf_send)
    SAFE_DEALLOCATE(proc2msg_send)
    SAFE_DEALLOCATE(local_col_count)
    SAFE_DEALLOCATE(req_send)
    SAFE_DEALLOCATE(my_status)
    SAFE_DEALLOCATE(req_recv)

    !WRITE(*,*) sig%need_advanced

    POP_SUB(epscopy.epscopy_subspace)

  end subroutine
#endif

end subroutine epscopy

end module epscopy_m
