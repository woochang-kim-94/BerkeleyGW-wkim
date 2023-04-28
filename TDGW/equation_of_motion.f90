!===============================================================================
! 
! Types:
! (1) td_equation
!
! Routines:
!
! (1) init_green_function()
!     initialize lesser green's function
!
!===============================================================================

#include "f_defs.h"

module equation_of_motion_m

  use global_m
  use utility_m
  use timing_m, only: timing => tdgw_timing

  implicit none

  private

  public :: &
    td_equation, &
    init_parallelization, &
    init_green_function, &
    init_hamiltonian, &
    solve_bse_hamiltonian, &
    get_interaction, &
    init_td_equation_from_copy, &
    free_td_equation

  type td_equation
    integer :: iq          ! q-point label, derived type object contains data for iq
    integer :: nk          ! number of total k-points 
    integer :: nb          ! number of total bands
    integer :: nvb         ! number of valence bands
    integer :: ncb         ! number of conduction bands
    integer :: nspin       ! number of spin channels
    integer :: nspinor     ! number of spinor components

    real(DP) :: itlabel    ! time step label (changing during time evolution)
                           ! we make this variable real to allow flexibility
    real(DP) :: efield     ! efield at each time step (changing during time evolution)

    real(DP) :: qpt(3)     ! q-point
    real(DP), allocatable :: kpts(:,:)  ! k-point list
    real(DP), allocatable :: kqpts(:,:) ! k+q point list
    integer, allocatable  :: kqmap(:)   ! map index of k+q to k list

    real(DP) :: celvol         ! cell volume
    real(DP) :: blat           ! reciprocal lattice constant
    real(DP) :: bvec(3,3)      ! reciprocal lattice vectors

    real(DP), allocatable :: band_energy(:,:)         ! equilibrium band energy (nb, nk)
    complex(DPC), allocatable :: gf(:,:,:)            ! lesser green`s function (nb, nb, nk)
                                                      !    G_{mk, nk+q}, k ordered as kpts, m->c, n->v (in BSE-TDA limit)
                                                      !   this mn->cv convention is to follow formula translator
    complex(DPC), allocatable :: hamiltonian(:,:,:)   ! hamiltonian (nb, nb, nk) [(nc, nv, nk)]
                                                      !    k ordered as kpts 
    complex(DPC), allocatable :: interaction(:,:,:)   ! interaction term (v+W)G = V + Sigma, (nb, nb, nk) [(nc, nv, nk)]
                                                      !    k ordered as kpts
    complex(DPC), allocatable :: dipole(:,:,:)        ! dipole transitions (nb, nb, nk) [(nc, nv, nk)]
                                                      !    k ordered as kpts
    complex(DPC), allocatable :: vwint(:,:)           ! interaction matrix elements: bare v + screened W, distributed 
                                                      !    (kcv, k`c`v`)

    integer :: mpi_k_start    ! MPI loop bounds
    integer :: mpi_k_end
    integer :: mpi_kp_start
    integer :: mpi_kp_end
    integer :: mpi_qp_start
    integer :: mpi_qp_end

  end type td_equation

contains

!=========================================================================
!> initialize parallelization parameters
subroutine init_parallelization(nq, tdevol)

  integer :: nq
  type(td_equation), intent(inout) :: tdevol
  !
  integer :: kstart, kend
  integer :: inode

  PUSH_SUB(init_parallelization)

  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  ! find distribution of (k, k`) pairs
  call mpi_dispatcher(tdevol%nk, tdevol%nk, tdevol%mpi_k_start, tdevol%mpi_k_end, tdevol%mpi_kp_start, tdevol%mpi_kp_end)

  ! find distribution of (k, q`) pairs
  call mpi_dispatcher(tdevol%nk, nq, kstart, kend, tdevol%mpi_qp_start, tdevol%mpi_qp_end)

  if ((kstart .ne. tdevol%mpi_k_start) .or. (kend .ne. tdevol%mpi_k_end)) then
    call die('Error: MPI distribution inconsistent for (k, k`) and (k, q`) pairs.')
  endif

  ! print information if iq == 1
  if(tdevol%iq == 1) then
    do inode = 1, peinf%npes
      if ((peinf%inode == inode - 1) .and. (inode <= 10)) then
        write(*,*) 'Parallel indices of inode', peinf%inode
        write(*,*) '   k_start,  k_end', tdevol%mpi_k_start, tdevol%mpi_k_end
        write(*,*) '  kp_start, kp_end', tdevol%mpi_kp_start, tdevol%mpi_kp_end
        write(*,*) '  qp_start, qp_end', tdevol%mpi_qp_start, tdevol%mpi_qp_end
      endif
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    enddo

    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    if (peinf%inode == 0) then
      write(*,*) '...'
      write(*,*) '...'
      write(*,*) '...'
      write(*,*) 'Only printing parallelization information for the first q and first few k'
      write(*,*)
    endif
  endif

  POP_SUB(init_parallelization)

end subroutine init_parallelization


!=========================================================================
!> initialize green`s function
subroutine init_green_function(tdevol)

  type(td_equation), intent(inout) :: tdevol

  integer :: ib, ik
  real(DP) :: qnorm

  PUSH_SUB(init_green_function)

  call timing%start(timing%init_gf_hamil)

  SAFE_ALLOCATE(tdevol%gf, (tdevol%nb, tdevol%nb, tdevol%nk))

  tdevol%gf(:,:,:) = (0.0d0,0.0d0)

  tdevol%itlabel = 1.0d0   ! initialize time step
                           ! itlabel = 1 records equilibrium starting time axis
                           ! this means time = -inf, before applying any field
  tdevol%efield = 0.0d0  ! initially no field
                         ! since it represents time = -inf
                         ! note that although by setup, it = 1^+ can have finite efector
                         ! but here when (it == 1), we emphasize that it records equilibrium situation
                         ! the real time dynamics starts at it = 2

  ! iq = 0 must be gamma point
  ! q = 0 point is special in TD-aGW formalism since it initializes nonzero elements in the equation of motion
  !   even for finite-q case and even in the linear-response limit, q = 0 point is always needed

  qnorm = sqrt(dot_product(tdevol%qpt, tdevol%qpt))

  if (tdevol%iq == 1) then
    if (qnorm .gt. TOL_Small) then
      call die('The first finite_q point must be Gamma point, and is always needed in any time-dependent simulations.')
    endif

    do ik = 1, tdevol%nk
      do ib = 1, tdevol%nvb
        tdevol%gf(ib,ib,ik) = (0.0d0, 1.0d0)  ! diagonal (in bands and k (i.e. q=0)) elements are:
                                              ! 1.0j for valence and 0.0 for conduction in the initial condition
                                              ! all other matrix elements are zero
      enddo
    enddo
  else
    if (qnorm .lt. TOL_Small) then
      if (peinf%inode == 0) write(*,*) 'Warning: Gamma finite_q point must be listed as the first point.'
    endif
  endif

  call timing%stop(timing%init_gf_hamil)

  POP_SUB(init_green_function)

end subroutine init_green_function

!=========================================================================
!> initialize hamiltonian
subroutine init_hamiltonian(tdevol)

  type(td_equation), intent(inout) :: tdevol

  integer :: ib, ik

  PUSH_SUB(init_hamiltonian)

  call timing%start(timing%init_gf_hamil)

  SAFE_ALLOCATE(tdevol%hamiltonian, (tdevol%nb, tdevol%nb, tdevol%nk))
  SAFE_ALLOCATE(tdevol%interaction, (tdevol%nb, tdevol%nb, tdevol%nk))  ! we initialize interaction here

  ! tdevol%hamiltonian saves H^0 - V^0[G^0] - \Sigma^0[G^0]
  ! H^0: equilibrium hamiltonian contains band energies as diagonals
  ! V^0 + \Sigma^0 = interaction^0: equilibrium interaction term as function of equilibrium green`s function
  ! we subtract this interaction^0 term at the initialization stage to make time propagation straightforward

  do ik = 1, tdevol%nk
    do ib = 1, tdevol%nb
      tdevol%hamiltonian(ib, ib, ik) = tdevol%band_energy(ib, ik)  ! band energies ordered w.r.t. kpts
    enddo
  enddo

  ! get equilibrium interaction terms
  call get_interaction(tdevol)

  ! H^0 - interaction^0
  tdevol%hamiltonian(:,:,:) = tdevol%hamiltonian(:,:,:) - tdevol%interaction(:,:,:)

  call timing%stop(timing%init_gf_hamil)

  POP_SUB(init_hamiltonian)

end subroutine init_hamiltonian

!=========================================================================
!> solve BSE hamiltonian => linear-response limit of TD-aGW
subroutine solve_bse_hamiltonian(tdevol)

  type(td_equation), intent(in) :: tdevol   ! solving general finite-Q BSE for this tdevol

  complex(DPC), allocatable :: bse_hamiltonian(:,:)
  complex(DPC), allocatable :: evecs(:,:)
  real(DP), allocatable :: evals(:)
  integer :: ib, jb, mb, nb        ! label full band range
  integer :: icb, jvb, mcb, nvb    ! label valence/conduction band range
  integer :: ik, jk
  integer :: ik_loc, jk_loc ! local k indices for distributed data saved in vwint
  integer :: i, j, ih, jh
  integer :: ncount
  integer :: dimens
  integer :: nval
  complex(DPC), allocatable :: work(:)
  real(DP), allocatable :: rwork(:)
  integer, allocatable :: iwork(:), ifail(:)
  integer :: info

  PUSH_SUB(solve_bse_hamiltonian)

  call timing%start(timing%solve_bse)

  ! TODO: supposedly only root node compute

  ! matrix dimension
  ! within Tamm-Dancoff approximation (TDA), including vc<->v'c' (excitations) transitions
  ! extended kernel includes vc <-> v'c' (excitations) and cv <-> c'v' (deexcitations) transitions
  dimens = tdevol%nvb * tdevol%ncb * tdevol%nk  ! consider BSE-TDA

  SAFE_ALLOCATE(bse_hamiltonian, (dimens, dimens))
  SAFE_ALLOCATE(evecs, (dimens, dimens))
  SAFE_ALLOCATE(evals, (dimens))

  ! initialize
  bse_hamiltonian(:,:) = (0.0d0,0.0d0)

  ! diagonal elements: E_{ck} - E_{vk+Q} for excitations
  ! serial version:
  !do ik = 1, tdevol%nk
  !  do jk = 1, tdevol%nk
  do ik = tdevol%mpi_k_start, tdevol%mpi_k_end
    do jk = tdevol%mpi_kp_start, tdevol%mpi_kp_end

      ik_loc = ik - tdevol%mpi_k_start + 1
      jk_loc = jk - tdevol%mpi_kp_start + 1

      do icb = 1, tdevol%ncb
        ib = tdevol%nvb + icb
        do jvb = 1, tdevol%nvb
          jb = jvb

          do mcb = 1, tdevol%ncb
            mb = tdevol%nvb + mcb
            do nvb = 1, tdevol%nvb
              nb = nvb

              i = (ik_loc - 1) * tdevol%nb * tdevol%nb + (ib - 1) * tdevol%nb + jb  ! kcv
              j = (jk_loc - 1) * tdevol%nb * tdevol%nb + (mb - 1) * tdevol%nb + nb  ! k`c`v`

              ih = (ik - 1) * tdevol%ncb * tdevol%nvb + (icb - 1) * tdevol%nvb + jvb
              jh = (jk - 1) * tdevol%ncb * tdevol%nvb + (mcb - 1) * tdevol%nvb + nvb

              ! H_{kcv, k`c`v`} = <vck|H|v`c`k`>
              bse_hamiltonian(ih, jh) = -tdevol%vwint(i, j)

              if (ih == jh) then
                ! diagonal band energy difference (E_{ck} - E_{vk+Q})
                bse_hamiltonian(ih, ih) = bse_hamiltonian(ih, ih) &
                                         + tdevol%band_energy(ib, ik) - tdevol%band_energy(jb, tdevol%kqmap(ik))
              endif ! ih == jh
            enddo ! nvb and nb
          enddo ! mcb and mb
        enddo ! jvb and jb
      enddo  ! icb and ib
    enddo ! jk
  enddo ! ik

  ! collect data
  ncount = dimens * dimens
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  call timing%start(timing%communication)
  ! MPI_IN_PLACE requires intent of recvbuf to be inout, bse_hamiltonian is declared locally
  call MPI_Allreduce(MPI_IN_PLACE, bse_hamiltonian, ncount, MPI_COMPLEX_DPC, MPI_SUM, MPI_COMM_WORLD, mpierr)
  call timing%stop(timing%communication)

  SAFE_ALLOCATE(work, (2*dimens))
  SAFE_ALLOCATE(rwork, (7*dimens))
  SAFE_ALLOCATE(iwork, (5*dimens))
  SAFE_ALLOCATE(ifail, (dimens))

  ! we call lapack function ZHEEVX to diagonalize
  call ZHEEVX('V','A','U',dimens,bse_hamiltonian,dimens,0.0d0,0.0d0,0,0,0.0d0, &
              & nval,evals,evecs,dimens,work,2*dimens,rwork,iwork,ifail,info)

  ! TODO: write BSE-TDA eigenvalue file
  !if (peinf%inode == 0) then
  ! for the processor dealing with the first set of k and q`, print all q
  if ((tdevol%mpi_k_start == 1) .and. (tdevol%mpi_qp_start == 1)) then
    write(*,*)
    write(*,*) 'Solving BSE-TDA for iq =', tdevol%iq
    write(*,*) 'Finite Q shift = ', tdevol%qpt(:)
    write(*,*) 'LAPACK ZHEEVX info', info
    write(*,*)
    write(*,*) 'eigenvalues (eV)'
    write(*,*) evals(1:30) * ryd
    write(*,*)
  endif

  SAFE_DEALLOCATE(work)
  SAFE_DEALLOCATE(rwork)
  SAFE_DEALLOCATE(iwork)
  SAFE_DEALLOCATE(ifail)
  SAFE_DEALLOCATE(evals)
  SAFE_DEALLOCATE(evecs)
  SAFE_DEALLOCATE(bse_hamiltonian)

  call timing%stop(timing%solve_bse)

  POP_SUB(solve_bse_hamiltonian)

end subroutine solve_bse_hamiltonian


!=========================================================================
!> get interaction GW and GV
subroutine get_interaction(tdevol)

  type(td_equation), intent(inout) :: tdevol  ! with input vwint and gf, get output interaction
                                              ! for a specific finite-Q tdevol,
                                              ! interaction(Q) = gf(Q) * vwint(Q)
                                              ! where vwint(Q) are matrix elements so they do not contain finite Q phase
                                              ! that said, a finite-Q interaction term is constructed with its own gf and vwint
  integer :: ik, jk, ib, jb, mb, nb
  integer :: ik_loc, jk_loc  ! local k indices for distributed data saved in vwint
  integer :: i, j
  integer :: ncount
  complex(DPC) :: imag_one

  PUSH_SUB(get_interaction)

  call timing%start(timing%get_interaction)

  imag_one = (0.0d0,1.0d0)
  tdevol%interaction(:,:,:) = (0.0d0,0.0d0)

  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  ! Here is the one of the most numerical intensive loops
  ! we parallel over ik (k) and jk (k`)
  do ik = tdevol%mpi_k_start, tdevol%mpi_k_end   ! serial version: do ik = 1, tdevol%nk
    do jk = tdevol%mpi_kp_start, tdevol%mpi_kp_end   ! serial version: do jk = 1, tdevol%nk

      ik_loc = ik - tdevol%mpi_k_start + 1
      jk_loc = jk - tdevol%mpi_kp_start + 1

      do ib = 1, tdevol%nb
        do jb = 1, tdevol%nb  ! vb in TDA
          do mb = 1, tdevol%nb
            do nb = 1, tdevol%nb  ! vb in TDA
              i = (ik_loc - 1) * tdevol%nb * tdevol%nb + (ib - 1) * tdevol%nb + jb  ! kcv
              j = (jk_loc - 1) * tdevol%nb * tdevol%nb + (mb - 1) * tdevol%nb + nb  ! k`c`v`
  
              ! interaction(c,v,k) = imag_one * sum_{c`v`k`} G(c`,v`,k`) * vwint(kcv, k`c`v`)
              tdevol%interaction(ib, jb, ik) = tdevol%interaction(ib, jb, ik) & 
                                               + imag_one * tdevol%gf(mb, nb, jk) * tdevol%vwint(i, j)
  
            enddo ! nb
          enddo ! mb
        enddo ! jb
      enddo ! ib
    enddo ! jk
  enddo ! ik

  ! collect data
  ncount = tdevol%nk * tdevol%nb * tdevol%nb
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  call timing%start(timing%communication)
  ! MPI_IN_PLACE requires intent of recvbuf to be inout
  call MPI_Allreduce(MPI_IN_PLACE, tdevol%interaction, ncount, MPI_COMPLEX_DPC, MPI_SUM, MPI_COMM_WORLD, mpierr)
  call timing%stop(timing%communication)

  call timing%stop(timing%get_interaction)

  POP_SUB(get_interaction)

end subroutine get_interaction


!=========================================================================
!> initialize one object of tdevol by copying from another
subroutine init_td_equation_from_copy(tdevol_src, tdevol)

  type(td_equation), intent(in) :: tdevol_src
  type(td_equation), intent(inout) :: tdevol
  !
  integer :: my_nk, my_nkp

  PUSH_SUB(init_td_equation_from_copy)

  tdevol%iq = tdevol_src%iq
  tdevol%nk = tdevol_src%nk
  tdevol%nb = tdevol_src%nb
  tdevol%nvb = tdevol_src%nvb
  tdevol%ncb = tdevol_src%ncb
  tdevol%nspin = tdevol_src%nspin
  tdevol%nspinor = tdevol_src%nspinor

  tdevol%itlabel = tdevol_src%itlabel
  tdevol%efield = tdevol_src%efield

  tdevol%qpt(:) = tdevol_src%qpt(:)

  SAFE_ALLOCATE(tdevol%kpts, (3, tdevol%nk))
  SAFE_ALLOCATE(tdevol%kqpts, (3, tdevol%nk))
  SAFE_ALLOCATE(tdevol%kqmap, (tdevol%nk))

  tdevol%kpts(:,:) = tdevol_src%kpts(:,:)
  tdevol%kqpts(:,:) = tdevol_src%kqpts(:,:)
  tdevol%kqmap(:) = tdevol_src%kqmap(:)

  tdevol%celvol = tdevol_src%celvol
  tdevol%blat = tdevol_src%blat
  tdevol%bvec(:,:) = tdevol_src%bvec(:,:)

  SAFE_ALLOCATE(tdevol%band_energy, (tdevol%nb, tdevol%nk))
  SAFE_ALLOCATE(tdevol%gf, (tdevol%nb, tdevol%nb, tdevol%nk))
  SAFE_ALLOCATE(tdevol%hamiltonian, (tdevol%nb, tdevol%nb, tdevol%nk))
  SAFE_ALLOCATE(tdevol%interaction, (tdevol%nb, tdevol%nb, tdevol%nk))
  SAFE_ALLOCATE(tdevol%dipole, (tdevol%nb, tdevol%nb, tdevol%nk))
  ! vwint allocated at the end of this subroutine

  tdevol%band_energy(:,:) = tdevol_src%band_energy(:,:)
  tdevol%gf(:,:,:) = tdevol_src%gf(:,:,:)
  tdevol%hamiltonian(:,:,:) = tdevol_src%hamiltonian(:,:,:)
  tdevol%interaction(:,:,:) = tdevol_src%interaction(:,:,:)
  tdevol%dipole(:,:,:) = tdevol_src%dipole(:,:,:)
  ! vwint copied at the end of this subroutine

  tdevol%mpi_k_start = tdevol_src%mpi_k_start
  tdevol%mpi_k_end = tdevol_src%mpi_k_end
  tdevol%mpi_kp_start = tdevol_src%mpi_kp_start
  tdevol%mpi_kp_end = tdevol_src%mpi_kp_end
  tdevol%mpi_qp_start = tdevol_src%mpi_qp_start
  tdevol%mpi_qp_end = tdevol_src%mpi_qp_end

  ! vwint is distributed
  my_nk  = tdevol%mpi_k_end  - tdevol%mpi_k_start  + 1
  my_nkp = tdevol%mpi_kp_end - tdevol%mpi_kp_start + 1

  if ((my_nk .gt. 0) .and. (my_nkp .gt. 0)) then
    SAFE_ALLOCATE(tdevol%vwint, (tdevol%nb*tdevol%nb*my_nk, tdevol%nb*tdevol%nb*my_nkp))
    tdevol%vwint(:,:) = tdevol_src%vwint(:,:)
  endif

  POP_SUB(init_td_equation_from_copy)

end subroutine init_td_equation_from_copy


!=========================================================================
!> free td_equation object
subroutine free_td_equation(tdevol)

  type(td_equation), intent(inout) :: tdevol

  PUSH_SUB(init_td_equation_from_copy)

  SAFE_DEALLOCATE(tdevol%kpts)
  SAFE_DEALLOCATE(tdevol%kqpts)
  SAFE_DEALLOCATE(tdevol%kqmap)
  SAFE_DEALLOCATE(tdevol%band_energy)
  SAFE_DEALLOCATE(tdevol%gf)
  SAFE_DEALLOCATE(tdevol%hamiltonian)
  SAFE_DEALLOCATE(tdevol%interaction)
  SAFE_DEALLOCATE(tdevol%dipole)
  SAFE_DEALLOCATE(tdevol%vwint)

  POP_SUB(free_td_equation)

end subroutine free_td_equation

end module equation_of_motion_m
