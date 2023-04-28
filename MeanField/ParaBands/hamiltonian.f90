!===============================================================================
!
! Routines:
!
! 1. ham_solve()            Originally By FHJ    Last Modified  Apr/2015 (FHJ)
!
!    Builds the dense PW Hamiltonian and diagonalizes it, returning the
!    number of requested bands and energies.
!
! 2. ham_build_column()     Originally By FHJ    Last Modified  Apr/2015 (FHJ)
!
!    Build a column of the Hamiltonian in PW representation.
!
!===============================================================================

#include "f_defs.h"
#if defined MPI && !defined USESCALAPACK
  #error ScaLAPACK is required for MPI builds.
#endif

module hamiltonian_m

  use global_m
  use misc_m,          only: findvector
  use io_utils_m,      only: progress_info, progress_init, progress_step, progress_free
  use input_utils_m,   only: kinetic_energies
  use inread_m,        only: pb_params_t
  use pseudopot_m,     only: pseudopot_t
  use distribution_m,  only: distrib_mat_t
  use diag_driver_m,   only: diag_driver

  implicit none

  private

  public :: ham_solve

contains


!> Builds and diagonalizes the dense DFT Hamiltonian in PW representation
!> On exit, returns:
!> - band energies in en.
!> - band-distributed wavefunctions (eigenvectors) in wfn_d.
subroutine ham_solve(params, dm_wfn, gvec, kp, crys, isortk, vsc, pseudopot, ik, en, wfn_d)
  use, intrinsic :: iso_fortran_env

  type(pb_params_t), intent(in) :: params
  class(distrib_mat_t), intent(in) :: dm_wfn
  type(gspace), intent(in) :: gvec
  type(kpoints), intent(in) :: kp
  type(crystal), intent(in) :: crys
  integer, intent(in) :: isortk(kp%ngkmax)
  SCALAR, intent(in) :: vsc(gvec%ng, kp%nspin)
  class(pseudopot_t), allocatable, intent(in) :: pseudopot
  integer, intent(in) :: ik
  real(DP), intent(inout) :: en(params%nb, kp%nspin)
  SCALAR, intent(inout) :: wfn_d(dm_wfn%Ml, dm_wfn%Nl, kp%nspin)

  integer :: is, icol_l, icol_g, igp_g, ispinorp_g, ispinor_deeq
  integer :: ib, nprint
  integer :: ngk !< Number of gvectors for this k-points
  integer :: ngkx !< ngk * nspinor

  complex(DPC), allocatable :: vkb_atom(:,:)
  integer :: nh, iat

  SCALAR, allocatable :: ham_d(:,:)
  real(DP) :: ekin(gvec%ng)
  type(distrib_mat_t) :: dm_ham
  type(progress_info) :: prog_info !< a user-friendly progress report
  integer(int64) :: clock_count
  real(DP) :: clock_inv_rate, t0, t1

  PUSH_SUB(ham_solve)

  ngkx = dm_wfn%M
  ngk = ngkx / kp%nspinor
!------------------------
! Setup scalapack buffer and layout for 1d-distributed Hamiltonian
  call dm_ham%setup_mat_1d(ngkx, ngkx, params%kpp)
  SAFE_ALLOCATE(ham_d, (dm_ham%Ml, dm_ham%Nl))

!------------------------
! Compute kinetic energies of G-vectors
  call timacc(41,1)
  call kinetic_energies(gvec, crys%bdot, ekin, kp%rk(1:3,ik))
  call timacc(41,2)

  do is = 1, kp%nspin
!------------------------
! Build Hamiltonian for each spin component

    if (kp%nspin>1.and.params%kpp%inode==0) then
      write(params%kpp%iunit,'(1x,a,i0," / ",i0/)') 'Dealing with spin ', is, kp%nspin
    endif

    if (params%solver_alg >= -1) then ! FHJ: Don`t build H for dummy solver

    call timacc(42,1)
    call progress_init(prog_info, 'building local part of H', 'distributed column', &
      DIVUP(dm_ham%Nown,peinf%nthreads), should_print=params%kpp%inode==0, iunit=params%kpp%iunit)
    if (params%kpp%inode==0) FLUSH(params%kpp%iunit)
    ! FHJ: Loop over all local columns (which combine G-vector and spinor indices)
    !$OMP PARALLEL DO PRIVATE(icol_l, icol_g, igp_g, ispinorp_g) DEFAULT(shared) SCHEDULE(static)
    do icol_l = 1, dm_ham%Nl
      if (get_thread_id()==0) then
        call progress_step(prog_info)
        if (params%kpp%inode==0) FLUSH(params%kpp%iunit)
      endif
      icol_g = params%kpp%inode*dm_ham%Nb + icol_l
      if (icol_g>ngkx) cycle
      igp_g = mod(icol_g-1, ngk) + 1
      ispinorp_g = DIVUP(icol_g, ngk)
      ! Local potential is diagonal in the spinor component
      ham_d(:,icol_l) = ZERO
      call ham_build_column_loc(params, gvec%ng, kp%ngkmax, ngk, &
        igp_g, ekin, isortk, ham_d(1+(ispinorp_g-1)*ngk:(ispinorp_g)*ngk,icol_l), gvec, vsc(:,is))
    enddo ! igcol_l
    !$OMP END PARALLEL DO
    call progress_free(prog_info)
    call timacc(42,2)

    if (allocated(pseudopot)) then
      SAFE_ALLOCATE(vkb_atom, (kp%ngkmax,pseudopot%nhm))
      call timacc(43,1)
      call progress_init(prog_info, 'building non-local part of H', 'distributed column', &
        DIVUP(dm_ham%Nown,peinf%nthreads)*pseudopot%nat, &
        should_print=params%kpp%inode==0, iunit=params%kpp%iunit)
      if (params%kpp%inode==0) FLUSH(params%kpp%iunit)
      do iat = 1, pseudopot%nat
        nh = pseudopot%get_nh_for_atom(iat)
        if (nh==0) cycle
        call pseudopot%get_vkb_for_atom(iat, is, vkb_atom, params%kpp)

        ! FHJ: Loop over all local columns (which combine G-vector and spinor indices)
        !$OMP PARALLEL DO PRIVATE(icol_l, icol_g, igp_g, ispinorp_g, ispinor_deeq) DEFAULT(shared) SCHEDULE(static)
        do icol_l = 1, dm_ham%Nl
          if (get_thread_id()==0) then
            call progress_step(prog_info)
            if (params%kpp%inode==0) FLUSH(params%kpp%iunit)
          endif
          icol_g = params%kpp%inode*dm_ham%Nb + icol_l
          if (icol_g>ngkx) cycle
          igp_g = mod(icol_g-1, ngk) + 1
          ispinorp_g = DIVUP(icol_g, ngk)
          if (kp%nspinor==1) then
            call ham_build_column_nl_col(kp%ngkmax, ngk, igp_g, &
              pseudopot%nhm, nh, vkb_atom(:,1:nh), &
              pseudopot%deeq(:,:,iat,is), ham_d(:,icol_l))
          else
            ispinor_deeq = ispinorp_g
            call ham_build_column_nl_nc(kp%ngkmax, ngk, igp_g, &
              pseudopot%nhm, nh, vkb_atom(:,1:nh), &
              pseudopot%deeq_nc(:,:,iat,ispinor_deeq), ham_d(1:ngk,icol_l))
            ispinor_deeq = ispinorp_g + 2
            call ham_build_column_nl_nc(kp%ngkmax, ngk, igp_g, &
              pseudopot%nhm, nh, vkb_atom(:,1:nh), &
              pseudopot%deeq_nc(:,:,iat,ispinor_deeq), ham_d(ngk+1:2*ngk,icol_l))
          endif
        enddo
        !$OMP END PARALLEL DO
      enddo
      call timacc(43,2)
      SAFE_DEALLOCATE(vkb_atom)
      call progress_free(prog_info)
    endif

    endif ! solver_alg>=-1
    !if (ik==1) then
    !  call logit('Dumping Hamiltonian', params%kpp%inode==0, params%kpp%iunit)
    !  call ham_dump()
    !endif

    if (params%kpp%inode==0) then
      call system_clock(count_rate=clock_count)
      clock_inv_rate = 1d0/clock_count
      call system_clock(count=clock_count)
      t0 = clock_count * clock_inv_rate
      write(params%kpp%iunit,'(1x,a)') 'Diagonalizing H.'
      FLUSH(params%kpp%iunit)
    endif

    call logit('Calling diag_driver', params%kpp%inode==0, params%kpp%iunit)
    call timacc(45,1)
    call diag_driver(dm_ham, ham_d, dm_wfn, wfn_d(:,:,is), en(:,is), params)
    call timacc(45,2)
    call logit('Done calling diag_driver', params%kpp%inode==0, params%kpp%iunit)

    if (params%kpp%inode==0) then
      call system_clock(count=clock_count)
      t1 = clock_count * clock_inv_rate
      write(params%kpp%iunit,'(1x,a)') 'Done diagonalizing H.'
      write(params%kpp%iunit,'(1x,a,f0.3,a/)') 'Time elapsed: ', t1-t0, ' s.'
      if (params%has_wfn_in) then
        nprint = min(params%nb, maxval(kp%ifmax)+2)
      else
        nprint = count(en(:,is) < pseudopot%ef + TOL_SMALL)
        nprint = min(params%nb, nprint+2)
      endif
      write(params%kpp%iunit,'(1x,a,i0,a)') 'First ', nprint, ' eigenvalues:'
      write(params%kpp%iunit,'(1x,a4,2(1x,a11))') '#','en (Ry)', 'en (eV)'
      do ib = 1, nprint
        write(params%kpp%iunit,'(1x,i4,2(1x,f11.6))') ib, en(ib,is), en(ib,is)*ryd
      enddo
      write(params%kpp%iunit,'()')
      FLUSH(params%kpp%iunit)
    endif
  enddo !is
  SAFE_DEALLOCATE(ham_d)

  POP_SUB(ham_solve)

contains

  subroutine ham_dump
    integer :: igcol_g, igcol_l, ipe
    integer, external :: INDXG2P, INDXG2L
    SCALAR :: ham_buf(dm_ham%M)

    do igcol_g = 1, dm_ham%N
      ipe = INDXG2P(igcol_g, dm_ham%Nb, 0, 0, params%kpp%npes)
      igcol_l = INDXG2L(igcol_g, dm_ham%Nb, 0, 0, params%kpp%npes)
      if (ipe==params%kpp%inode) then
        ! I`ve got the band
        ham_buf(:) = ZERO
        ham_buf(1:ngk) = ham_d(1:ngk,igcol_l)
        if (params%kpp%inode/=0) then
          ! Need to send it to root
#ifdef MPI
          call MPI_Send(ham_buf(1), kp%ngkmax, MPI_SCALAR, &
            0, igcol_g, params%kpp%comm, mpierr)
#endif
        endif
      endif
      if (params%kpp%inode==0) then
        if (ipe/=params%kpp%inode) then
          ! Root already has ham_buf
#ifdef MPI
          call MPI_Recv(ham_buf(1), kp%ngkmax, MPI_SCALAR, &
            ipe, igcol_g, params%kpp%comm, MPI_STATUS_IGNORE, mpierr)
#endif
        endif
        write(666,'(10000(es22.15,1x))') MYCONJG(ham_buf)
      endif
    enddo

  end subroutine ham_dump

end subroutine ham_solve


!===============================================================================


!> Builds a column of the Hamiltonian in plane wave representation.
!! On exit, it returns the igcol_l`th local colum in ham_col(1:ngk_g,igcol_l)
!!
!! We only add the local pot. and kinetic energy here
subroutine ham_build_column_loc(params, ng, ngkmax, ngk, igp_g, ekin, isortk, &
  ham_col, gvec, vsc)

  type(pb_params_t), intent(in) :: params
  integer, intent(in) :: ng, ngkmax, ngk, igp_g
  real(DP), intent(in) :: ekin(ng)
  integer, intent(in) :: isortk(ngkmax)
  SCALAR, intent(out) :: ham_col(ngk)
  type(gspace), intent(in) :: gvec
  SCALAR, intent(in) :: vsc(ng)

  complex(DPC) :: ham_tmp(ngk)
  integer :: gg(3), idx, ig

  PUSH_SUB(ham_build_column_loc)

  ham_tmp(:) = (0d0, 0d0)

!------------------------
! Local pseudopotential
  if (params%has_vsc) then
    ! If psi(G) is a delta function, psi(G) = delta_{G,G_in}, then
    ! vsc(G) (*) psi(G) is simply a shift in G-space, vsc(G-G_in)
    do ig = 1, ngk
      gg = gvec%components(:,isortk(ig)) - gvec%components(:,isortk(igp_g))
      call findvector(idx, gg, gvec)
      if (idx>0) ham_tmp(ig) = vsc(idx)
    enddo
  endif

!------------------------
! Kinetic energy
  ham_tmp(igp_g) = ham_tmp(igp_g) + ekin(isortk(igp_g))

!------------------------
! Convert from complex to scalar
  do ig = 1, ngk
    ham_col(ig) = SCALARIFY(ham_tmp(ig))
  enddo

  POP_SUB(ham_build_column_loc)

end subroutine ham_build_column_loc


!> Builds a column of the Hamiltonian in plane wave representation.
!! On exit, it returns the igcol_l`th local colum in ham_col(1:ngk_g,igcol_l)
!!
!! We only add the non-local part of the PP (KB projectors) here. Some notes:
!!  - A somewhat similar code is found in QE (PW/src/add_vuspsi.f90)
!!  - What we really compute is the NL contribution from the KB projectors
!!    from a specific atom "a":
!!    V_NL^a = \sum_{h,h`}^{N_h(a)}
!!           \ket{\beta_{a h s k}} d_{h h` a s} \bra{\beta_{a h` s k}}
!!  - We first perform the multilpication over h` with matmul, since nh is small
!!  - We then perform the leftmost multiplication over h with BLAS
subroutine ham_build_column_nl_col(ngkmax, ngk, igp_g, nhm, nh, vkb_atom, &
  deeq, ham_col)

  integer, intent(in) :: ngkmax, ngk, igp_g, nhm, nh
  complex(DPC), intent(in) :: vkb_atom(ngkmax,nh)
  real(DP), intent(in) :: deeq(nhm,nhm)
  SCALAR, intent(inout) :: ham_col(ngk)

  complex(DPC) :: ps(nh), betapsi(nh), ham_tmp(ngk)

  PUSH_SUB(ham_build_column_nl_col)

  betapsi(:) = CONJG(vkb_atom(igp_g,:))
  ps(:) = matmul(deeq(1:nh,1:nh), betapsi)
  call zgemv('N', ngk, nh, (1d0,0d0), vkb_atom(1,1), ngkmax, &
    ps(1), 1, (0d0,0d0), ham_tmp(1), 1)
  ham_col = ham_col + SCALARIFY(ham_tmp)

  POP_SUB(ham_build_column_nl_col)

end subroutine ham_build_column_nl_col

!> Builds a column of the Hamiltonian in plane wave representation.
!! On exit, it returns the igcol_l`th local colum in ham_col(1:ngk_g,igcol_l)
!!
!! We only add the non-local part of the PP (KB projectors) here. Some notes:
!!  - A somewhat similar code is found in QE (PW/src/add_vuspsi.f90)
!!  - What we really compute is the NL contribution from the KB projectors
!!    from a specific atom "a":
!!    V_NL^a = \sum_{h,h`}^{N_h(a)}
!!           \ket{\beta_{a h s k}} d_{h h` a s} \bra{\beta_{a h` s k}}
!!  - We first perform the multilpication over h` with matmul, since nh is small
!!  - We then perform the leftmost multiplication over h with BLAS
subroutine ham_build_column_nl_nc(ngkmax, ngk, igp_g, nhm, nh, vkb_atom, &
  deeq, ham_col)

  integer, intent(in) :: ngkmax, ngk, igp_g, nhm, nh
  complex(DPC), intent(in) :: vkb_atom(ngkmax,nh)
  complex(DPC), intent(in) :: deeq(nhm,nhm)
  SCALAR, intent(inout) :: ham_col(ngk)

  complex(DPC) :: ps(nh), betapsi(nh), ham_tmp(ngk)

  PUSH_SUB(ham_build_column_nl_nc)

  betapsi(:) = CONJG(vkb_atom(igp_g,:))
  ps(:) = matmul(deeq(1:nh,1:nh), betapsi)
  call zgemv('N', ngk, nh, (1d0,0d0), vkb_atom(1,1), ngkmax, &
    ps(1), 1, (0d0,0d0), ham_tmp(1), 1)
  ham_col = ham_col + SCALARIFY(ham_tmp)

  POP_SUB(ham_build_column_nl_nc)

end subroutine ham_build_column_nl_nc


end module hamiltonian_m
