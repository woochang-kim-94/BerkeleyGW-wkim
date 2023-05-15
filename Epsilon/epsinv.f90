!=================================================================================
!
! Routines:
!
! (1) epsinv()          Originally By (?)                Last Modified 5/1/2008 (JRD)
!
!     This routine:
!
!     1. Calculates epsilon based on chi.
!     2. Inverts epsilon.
!     3. Writes the result to unit = 12 if q0="zero" and unit 13 otherwise.
!
!=================================================================================

#include "f_defs.h"

module epsinv_m

  use global_m
  use inversion_m
  use io_utils_m
  use misc_m
  use scalapack_m
  use timing_m, only: timing => epsilon_timing
  use vcoul_generator_m
  use write_matrix_m
  implicit none

  private

  public :: &
#if defined MPI && defined USESCALAPACK
    invert_subspace, &
#endif
    epsinv

contains

subroutine epsinv(gvec, pol, ekin, q0, is_q0, crys, scal, kp, omega_plasma, iq, E_rpa)
  type (gspace), intent(in):: gvec
  type (polarizability), intent(in):: pol
  real(DP), intent(in):: ekin(gvec%ng)
  real(DP), intent(in):: q0(3)
  logical, intent(in):: is_q0
  type (crystal), intent(in):: crys
  type (scalapack), intent(in):: scal
  type (kpoints), intent(in):: kp
  real(DP), intent(in):: omega_plasma
  integer, intent(in):: iq
  real(DP), optional, intent(out):: E_rpa

  integer:: qgrid(3)
  real(DP):: q0vec(3)
  type (twork_scell):: work_scell
  integer:: i, j, iw, jj, ii, is, js, iunit, my_iq, ifreq
  integer:: irow, icol, icurr, irowm, icolm, icountr, icountc
  integer:: ig_l, ig_g, igp_l, igp_g
  integer, allocatable:: isorti(:)
  integer:: iscreen, nfq, iparallel, isize, ifreq_para, freq_grp_ind
  real(DP):: vc, oneoverq, avgcut
  real(DP):: epssum1R, epssum2R, epssum1R_rel, epssum2R_rel
#ifdef CPLX
  real(DP):: epssum1A, epssum2A, epssum1A_rel, epssum2A_rel
#endif
  SCALAR:: chitmp
  SCALAR, allocatable:: eps(:,:), ewng(:)
  real(DP), allocatable:: epsTemplate(:,:)
  real(DP), allocatable:: epsdiag(:,:,:), epsdiagt(:,:,:), vcoultemp(:)
  real(DP), allocatable:: vcoul(:)
#if defined(USESCALAPACK) && defined(CPLX)
  integer:: info, desc(9), jj_g
  complex(DPC), allocatable:: vcoul_dist(:), vcoul_inv_dist(:), wcoul_r(:,:)
#endif
  complex(DPC), allocatable:: chiRDyntmp(:)
  complex(DPC), allocatable:: epsRDyn(:,:,:), epsADyn(:,:,:)
  complex(DPC), allocatable:: epsRDyn_head(:), epsRDyn_head_temp(:)
  ! Auxiliary matrix for inversion
  complex(DPC), allocatable:: eps1Aux(:,:) 
  ! Auxiliary matrices for subspace truncation method
  complex(DPC), allocatable:: epsRDyn_sub(:,:,:), eps1Aux_sub(:,:)
  logical:: subspace, keep_full_eps_static
  integer:: nrow_loc_sub, ncol_loc_sub, neig_sub
  SCALAR:: epsheaddummy, wcoul0
  character*80:: filename
  type(progress_info):: prog_info
  type (scalapack):: scal_sub, scal_aux
  real(DP), allocatable:: integ_rpa_pol(:)
  integer:: freq_offset
  ! for timing (temporary)
  real(DP):: t1, t2
  ! to avoid large freq allocation in the subspace case
  integer:: nfreq_group_alloc

  PUSH_SUB(epsinv)

  SAFE_ALLOCATE(vcoul, (pol%nmtx))

  subspace = .FALSE.
#if defined MPI && defined USESCALAPACK
  if(pol%subspace .and. (.not. pol%need_full_chi)) subspace = .true.
#endif

  if(pol%freq_dep == 0) then
    SAFE_ALLOCATE(eps, (scal%npr, scal%npc))
    SAFE_ALLOCATE(ewng, (pol%nmtx))
  else
    SAFE_ALLOCATE(chiRDyntmp, (pol%nfreq_in_group))
    !XXXX
    nfreq_group_alloc = pol%nfreq_in_group
    if (subspace) then
      if ( pol%matrix_in_subspace_basis ) then
         nfreq_group_alloc = 1
      end if
    end if
    SAFE_ALLOCATE(epsRDyn, (scal%npr, scal%npc, nfreq_group_alloc))
#ifdef CPLX
    if (.not.pol%use_hdf5) then
      SAFE_ALLOCATE(epsADyn, (scal%npr, scal%npc, nfreq_group_alloc))
    endif
#endif
    !XXXX
    SAFE_ALLOCATE(eps1Aux, (scal%npr, scal%npc))
    ! subspace truncation specific stuff
    IF(subspace) THEN
      nrow_loc_sub = pol%nrow_local_sub
      ncol_loc_sub = pol%ncol_local_sub
      neig_sub       = pol%neig_sub
      SAFE_ALLOCATE(epsRDyn_sub, (MAX(1, nrow_loc_sub), MAX(1, ncol_loc_sub), pol%nfreq_in_group))
      SAFE_ALLOCATE(eps1Aux_sub, (MAX(1, nrow_loc_sub), MAX(1, ncol_loc_sub)))
      epsRDyn_sub = (0.d0, 0.d0)
    END IF
  endif

  
  SAFE_ALLOCATE(isorti, (gvec%ng))
      
!------------------------------
! Invert isrtx

!
! SIB: isorti is the inverse sort order for pol%isrtx.
! pol%isrtx has the sort indices for |q0+gvec%components|^2
!

! JRD XXX is the initialization necessary. probably very slow

  if (pol%freq_dep == 0) then
    eps(:,:)=ZERO
  else
    epsRDyn(:,:,:)=(0.d0, 0.d0)
#ifdef CPLX
    if (.not.pol%use_hdf5) then
      epsADyn(:,:,:)=(0.d0, 0.d0)
    endif
#endif
  endif

  
  vcoul(:)=0.0d0
  do i = 1, gvec%ng
    isorti(pol%isrtx(i)) = i
  end do


!-------------- Construct Dielectric Matrix---------------------------

!
! e(q+g, q+g`) = del(g, g`) - (8pi/(q+g)**2) chi(q+g, q+g`).  For spin-polarized
! calc., e(G+q, G`+q)=del(G, G`)- (8PI/(G+q)^2) SUM_spin chi(G+q, G`+q, ispin)
! Which is pol%chi(j, 1) as compiled in epsilon_main.f90.  If q--> 0, we have to treat
! the wings separately
!
! SIB:  using the Rydberg as our unit of energy
! if pol%icutv is on (coulomb cutoff) then we multiply the coulomb
! interaction by the appropriate factor (1-cos(vcut*|q0+g|))
!

!      if (peinf%inode .eq. 0) then
!        write(6, *) ' '
!        write(6, *) 'Calculating Coulomb Potential'
!      endif

  icurr = 0

! Generator Coulomb Interaction Array Vcoul

! For Epsilon, We want to treat all types of screening the same for vcoul.  Because
! we calculate it exactly

  avgcut = TOL_ZERO
  iscreen = 0
  ! FHJ: The following are not used when peinf%jobtypeeval == 1
  nfq = 0
  q0vec = 0d0
  iparallel = 1
  ! FHJ: FIXME-supercell truncation actually cares about qgrid.
  ! But is it doing the correct thing?
  qgrid(:) = 1

  call timing%start(timing%epsinv_vcoul)
  
  epsheaddummy = 0.0d0
  wcoul0 = 0.0d0
  
  IF(subspace) THEN
    ! we already have vcoul
    vcoul(:) = pol%vcoul_sub(:)
  ELSE

    if(pol%nfreq_group .eq. 1) then  
      call vcoul_generator(pol%icutv, pol%truncval, gvec, crys%bdot, crys%celvol, &
        nfq, pol%nmtx, pol%isrtx, iscreen, q0, q0vec, vcoul, &
        pol%iwritecoul, iparallel, avgcut, oneoverq, qgrid, epsheaddummy, &
        work_scell, .false.,wcoul0)
    else
      call vcoul_generator(pol%icutv, pol%truncval, gvec, crys%bdot, crys%celvol, &
        nfq, pol%nmtx, pol%isrtx, iscreen, q0, q0vec, vcoul, &
        pol%iwritecoul, iparallel, avgcut, oneoverq, qgrid, epsheaddummy, &
        work_scell, .false.,wcoul0, nfreq_group = pol%nfreq_group)
    endif

  END IF

  call timing%stop(timing%epsinv_vcoul)

!      write(6, *) 'Done VCoul'

  if (pol%freq_dep == 0) then
    do i = 1, pol%nmtx
      
      irow = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL), scal%nprow)
      if(irow .ne. scal%myprow) cycle
      
      vc = vcoul(i)

!-------------------------
! Construct non-interacting epsilon
! Static case

      do j = 1, pol%nmtx
        icol = MOD(INT(((j-1)/scal%nbl)+TOL_SMALL), scal%npcol)
        if(icol .eq. scal%mypcol) then
          icurr = icurr+1
          irowm = INT((icurr-1)/scal%npc+TOL_SMALL)+1
          icolm = MOD((icurr-1), scal%npc)+1

          eps(irowm, icolm) = ZERO
          if (i .eq. j) eps(irowm, icolm) = ONE
          chitmp = pol%chi(irowm, icolm, 1)
          eps(irowm, icolm) = eps(irowm, icolm) - vc*chitmp
        endif
      end do
    end do

  else  ! freq_dep > 0 below
  
    IF(subspace) THEN
      ! at this point we have the symmetrized polarizability 
      ! (already scaled by the SQRT of the coulomb potential) 
      ! for all frequency defined within the subspace of omega = 0, 
      ! additionally for omega = 0 we have the full (original basis) polarizability 
      ! we calculate epsilon just changing sing of chiRDyn and adding one on the diagonal
      epsRDyn_sub(:,:,:) = - pol%chiRDyn(:,:,:,1)   
      icountr = 0
      icountc = 0
      DO i = 1, neig_sub
        icol = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL), scal%npcol)
        irow = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL), scal%nprow)
        IF(irow == scal%myprow) icountr = icountr+1 
        IF(icol == scal%mypcol) icountc = icountc+1
        IF(irow /= scal%myprow) CYCLE
        IF(icol /= scal%mypcol) CYCLE
        epsRDyn_sub(icountr, icountc, :) = epsRDyn_sub(icountr, icountc, :) + 1.0D+00
      END DO
      ! same for omega = 0 full epsilon
      epsRDyn(:,:,1) = - pol%chiRDyn_sym_omega0(:,:)
      icountr = 0
      icountc = 0
      DO i = 1, pol%nmtx
        icol = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL), scal%npcol)
        irow = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL), scal%nprow)
        IF(irow == scal%myprow) icountr = icountr+1
        IF(icol == scal%mypcol) icountc = icountc+1
        IF(irow /= scal%myprow) CYCLE
        IF(icol /= scal%mypcol) CYCLE
        epsRDyn(icountr, icountc, 1) = epsRDyn(icountr, icountc, 1) + 1.0D+00
      END DO

    ELSE
      SAFE_ALLOCATE(epsTemplate, (scal%npr, scal%npc))
      SAFE_ALLOCATE(vcoultemp, (scal%npr))
      epsTemplate = 0D0
  
! JRD Set diagonal elements to 1 and define vcoultemp
  
      icountc = 0
      icountr = 0
      do i = 1, pol%nmtx
        icol = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL), scal%npcol)
        irow = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL), scal%nprow)
  
        if(icol .eq. scal%mypcol) icountc = icountc+1
        if(irow .eq. scal%myprow) icountr = icountr+1
  
        if(irow .ne. scal%myprow) cycle
        vcoultemp(icountr) = vcoul(i)
        if(icol .ne. scal%mypcol) cycle
  
        epsTemplate(icountr, icountc) = 1D0
      enddo
  
! JRD XXX if we block this we can probably hold epsTemplate in cache
      do ifreq = 1, pol%nfreq_in_group
        do i = 1, scal%npc
            epsRDyn(:,i, ifreq) = epsTemplate(:,i)- &
              vcoultemp(:)*pol%chiRDyn(:,i, ifreq, 1)
        end do
      end do
  
      SAFE_DEALLOCATE(epsTemplate)
      SAFE_DEALLOCATE(vcoultemp)
    END IF

  endif


3999 format(1x, a, i6, a, 2es25.15e3)
  if (peinf%inode == 0 .and. pol%freq_dep == 0) then
    write(6, 3999) 'q-pt ', iq, ': Head of Epsilon         = ', eps(1, 1)
    write(6, 3999) 'q-pt ', iq, ': Epsilon(2, 2)            = ', eps(2, 2)
  endif
  if (peinf%inode == 0 .and. pol%freq_dep > 0) then
    write(6, 3999) 'q-pt ', iq, ': Head of Epsilon         = ', epsRDyn(1, 1, 1)
    write(6, 3999) 'q-pt ', iq, ': Epsilon(2, 2)            = ', epsRDyn(2, 2, 1)
  endif

!-------------------------------------------------------------
! Print head versus frequency 


  if(pol%nfreq_group .eq. 1 .and. peinf%rank_f .eq. 0 .and. pol%freq_dep .eq. 2) then

    SAFE_ALLOCATE(epsRDyn_head, (pol%nfreq))
    epsRDyn_head = CMPLX(0d0, 0d0)
    IF(subspace) THEN
      epsRDyn_head(:)=epsRDyn_sub(1, 1, :)
    ELSE
      epsRDyn_head(:)=epsRDyn(1, 1, :)
    END IF
    

  elseif(pol%nfreq_group > 1 .and. peinf%inode .lt. peinf%npes .and. peinf%rank_f .eq. 0 .and. pol%freq_dep .eq. 2) then
    SAFE_ALLOCATE(epsRDyn_head, (pol%nfreq))
    SAFE_ALLOCATE(epsRDyn_head_temp, (pol%nfreq))
    epsRDyn_head = CMPLX(0d0, 0d0)
    epsRDyn_head_temp = CMPLX(0d0, 0d0)
! JRD XXX this nonsense may be slow
    do ifreq = 1, pol%nfreq
      freq_grp_ind = mod(ifreq-1, pol%nfreq_group)
      ifreq_para=(ifreq+pol%nfreq_group-1)/pol%nfreq_group
      if(freq_grp_ind .eq. peinf%rank_mtxel) then
        epsRDyn_head_temp(ifreq)=epsRDyn(1, 1, ifreq_para)
      endif
    enddo
#ifdef MPI
    call MPI_REDUCE(epsRDyn_head_temp(1), epsRDyn_head(1), pol%nfreq, MPI_COMPLEX_DPC, MPI_SUM, 0, &
     peinf%mtxel_comm, mpierr)
#endif
  endif

  if (peinf%rank_mtxel .eq. 0 .and. peinf%rank_f .eq. 0 .and. pol%freq_dep .eq. 2) then
    write(52, '("# q= ",3f12.6, " nmtx= ",i0)') q0(:), pol%nmtx
    write(52, *)
    do iw = 1, pol%nfreq
      write(52, '(f12.6, 4(1x, es13.6))') pol%dFreqGrid(iw), &
        dble(epsRDyn_head(iw)), aimag(epsRDyn_head(iw))
    enddo
  endif

!------------ Here we invert the epsilon matrix-----------------------------
!

  call timing%start(timing%epsinv_invert)


#ifdef USESCALAPACK
  if(pol%freq_dep .eq. 0) then
    call X(invert_with_scalapack)(pol%nmtx, scal, eps)
  endif

#ifdef CPLX
  ! FHJ: Prepares stuff to get the advanced epsilon from the retarded epsilon later on.
  ! We use the fact that W_A = (W_R)^H to write
  ! epsinv_A = (epsinv_R*v)^H * (1/v)
  ! We only do this if we use legacy Fortran binary Format. We never write the
  ! advanced matrix in the epsmat.h5 file.
  if (pol%freq_dep > 0.and..not.pol%use_hdf5) then
    SAFE_ALLOCATE(wcoul_r, (scal%npr, scal%npc))
    wcoul_r = (0d0, 0d0)
    call descinit(desc, pol%nmtx, pol%nmtx, scal%nbl, scal%nbl, 0, 0, &
      scal%icntxt, max(1, scal%npr), info)
    ! FHJ: vcoul_(inv_)dist contains the (inverse of the) Coulomb interaction
    ! for the distributed columns of the epsilon matrix that I own, and
    ! vcoul(jj_g) is never zero.
    SAFE_ALLOCATE(vcoul_dist, (scal%npc))
    SAFE_ALLOCATE(vcoul_inv_dist, (scal%npc))
    do jj = 1, scal%npc
      jj_g = indxl2g(jj, scal%nbl, scal%mypcol, 0, scal%npcol)
      vcoul_dist(jj) = vcoul(jj_g)
      vcoul_inv_dist(jj) = CMPLX(1d0/vcoul(jj_g), 0d0)
    enddo
  endif
#endif

  if (pol%freq_dep > 0) then
    call progress_init(prog_info, 'inversion of dielectric matrix', 'frequency', pol%nfreq_in_group)
    ! MDB: allocate stuff for RPA correlation energy integration
    ! DVF: we do this even for non-RPA calculations because gnu 5+does not allow optional arguments
    ! to subroutines that are not in modules, like invert_subspace
    SAFE_ALLOCATE(integ_rpa_pol, (pol%nfreq_imag))
    !
    do iw = 1, pol%nfreq_in_group
      call progress_step(prog_info)
! JRD XXX copy here no longer necessary
      IF(subspace) THEN
        !XXXXXXXXXX
        keep_full_eps_static = pol%keep_full_eps_static
        !XXXXXXXXXX
        IF(iw == 1 .AND. keep_full_eps_static) THEN
          ! for omega = 0 invert both epsilon matrices (subspace and not)
          ! to get an estimation of the error introduced by the subspace truncation
          eps1Aux(:,:) = epsRDyn(:,:,iw)
          call timing%start(timing%epsinv_omega_0)
          call zinvert_with_scalapack(pol%nmtx, scal, eps1Aux)
          call timing%stop(timing%epsinv_omega_0)
          epsRDyn(:,:,iw) = eps1Aux(:,:)
        END IF
        eps1Aux_sub(:,:) = epsRDyn_sub(:,:,iw)
        CALL invert_subspace(nrow_loc_sub, ncol_loc_sub, neig_sub, pol%nmtx, &
                               scal, pol, vcoul, eps1Aux_sub, eps1Aux, iw, integ_rpa_pol)
        if (iw /= 1 .or. (.not.keep_full_eps_static)) then
          if(nfreq_group_alloc == pol%nfreq_in_group) epsRDyn(:,:,iw) = eps1Aux(:,:)
          if(peinf%rank_f .eq. 0) epsRDyn_head(iw) =  eps1Aux(1, 1)
        endif
        ! save the subspace inverse epsilon retarded
        epsRDyn_sub(:,:,iw) = eps1Aux_sub(:,:)

        if ((pol%keep_full_eps_static .and. (iw == 1)) .or. .not.pol%matrix_in_subspace_basis) then
          ! Unsymmetrize dielectric matrix.
          ! FHJ: WARNING: never perform a nested loop over the global rows and
          ! columns, as these dimensions may be huge, and the code will spend
          ! a long time doing nothing. Instead, always loop over the local rows
          ! and columns and use indxl2g to get the corresponding global index.
          do igp_l = 1, scal%npc
            igp_g = indxl2g(igp_l, scal%nbl, scal%mypcol, 0, scal%npcol)
            vc = sqrt(vcoul(igp_g))
            do ig_l = 1, scal%npr
              ig_g = indxl2g(ig_l, scal%nbl, scal%myprow, 0, scal%nprow)
              epsRDyn(ig_l, igp_l, iw) = epsRDyn(ig_l, igp_l, iw) * sqrt(vcoul(ig_g)) / vc
            enddo
          enddo
        endif

      ELSE
        eps1Aux(:,:) = epsRDyn(:,:,iw)
        call zinvert_with_scalapack(pol%nmtx, scal, eps1Aux)
        epsRDyn(:,:,iw) = eps1Aux(:,:)
      END IF
#ifdef CPLX
      if (.not.pol%use_hdf5) then
        if ((.not.subspace) .or. (pol%keep_full_eps_static .and. iw == 1) .or. &
          (.not.pol%matrix_in_subspace_basis)) then
          ! FHJ: Scale columns of epsinv^R by vcoul to get W^R
          do jj = 1, scal%npc
            wcoul_r(:,jj) = eps1Aux(:,jj) * vcoul_dist(jj)
          enddo
          ! FHJ: Calculate adjoint of W^R and scale columns by 1/vcoul to get epsinv^A
          call pzgeadd('C', pol%nmtx, pol%nmtx, ONE, wcoul_r, 1, 1, desc, ZERO, eps1Aux, 1, 1, desc)
          do jj = 1, scal%npc
            eps1Aux(:,jj) = eps1Aux(:,jj) * vcoul_inv_dist(jj)
          enddo
          epsADyn(:,:,iw) = eps1Aux(:,:)
        else
          epsADyn(:,:,min(iw, nfreq_group_alloc)) = ZERO
        endif
      endif
#endif
    enddo
    call progress_free(prog_info)
    ! DVF: write out RPA information
    if(pol%do_rpa) then
      if(peinf%inode .eq. 0 .and. peinf%verb_max) then
        freq_offset = pol%nFreq-pol%nfreq_imag
        do iw = 1, pol%nfreq_imag
          write(*,*) IMAG(pol%dFreqBrd(iw+freq_offset)), integ_rpa_pol(iw) ! /(2.0D+00*PI_D)
        enddo
      endif
      E_rpa = 0.0D+00
      do iw = 1, pol%nfreq_imag
        E_rpa = E_rpa+pol%rpa_freq_grid(pol%nfreq_imag-iw+1) * integ_rpa_pol(iw) * 0.5D+00
      enddo
      E_rpa = E_rpa / (2.0D+00*PI_D)
      if(peinf%inode .eq. 0) write(*,*) "qpoint contribution to the total RPA energy = ", E_rpa
      SAFE_DEALLOCATE(integ_rpa_pol) 
    endif  ! do_rpa
    !
  endif

#else

! Serial Version

  if (pol%freq_dep == 0) then
    call X(invert_serial)(pol%nmtx, eps)
  else
    call progress_init(prog_info, 'inversion of dielectric matrix', 'frequency', pol%nfreq_in_group)
    do iw = 1, pol%nfreq_in_group
      call progress_step(prog_info)
! JRD XXX Copy no longer necessary
      eps1Aux(:,:) = epsRDyn(:,:,iw)
      call zinvert_serial(pol%nmtx, eps1Aux)
      epsRDyn(:,:,iw) = eps1Aux(:,:)
#ifdef CPLX
      if (.not.pol%use_hdf5) then
        ! FHJ: Scale columns of epsinv^R by vcoul to get W^R
        do jj = 1, pol%nmtx
          eps1Aux(:,jj) = eps1Aux(:,jj) * vcoul(jj)
        enddo
        ! FHJ: Calculate adjoint of W^R and scale columns by 1/vcoul to get epsinv^A
        eps1Aux = transpose(conjg(eps1Aux))
        do jj = 1, pol%nmtx
          eps1Aux(:,jj) = eps1Aux(:,jj) / vcoul(jj)
        enddo
        epsADyn(:,:,iw) = eps1Aux(:,:)
      endif
#endif
    enddo
    call progress_free(prog_info)
  endif

#endif

  call timing%stop(timing%epsinv_invert)

! Done inverting
!-----------------------------------------------------------------------------

  if (peinf%inode == 0 .and. pol%freq_dep == 0) then
    write(6, 3999) 'q-pt ', iq, ': Head of Epsilon Inverse = ', eps(1, 1)
    write(6, 3999) 'q-pt ', iq, ': Epsilon Inverse(2, 2)    = ', eps(2, 2)
  endif
  if (peinf%inode == 0 .and. pol%freq_dep > 0) then
    if(subspace) then
      if((.not. pol%keep_full_eps_static) .and. pol%matrix_in_subspace_basis) then
        write(6, '(1x, A)') 'Within the subspace basis:'
      end if
    end if
    write(6, 3999) 'q-pt ', iq, ': Head of Epsilon Inverse = ', epsRDyn(1, 1, 1)
    write(6, 3999) 'q-pt ', iq, ': Epsilon Inverse(2, 2)    = ', epsRDyn(2, 2, 1)
  endif

!----------- Print out independent matrix elements---------------------------

! JRD XXX this nonsense may be slow

  if(pol%nfreq_group .eq. 1 .and. peinf%rank_f .eq. 0 .and. pol%freq_dep .eq. 2) then

    if(nfreq_group_alloc == pol%nfreq_in_group) then
      epsRDyn_head(:)=epsRDyn(1, 1, :)
    endif ! (nfreq_group_alloc == pol%nfreq_in_group)

  elseif(pol%nfreq_group > 1 .and. peinf%inode .lt. peinf%npes .and. peinf%rank_f .eq. 0 .and. pol%freq_dep .eq. 2) then

    if(nfreq_group_alloc == pol%nfreq_in_group) then
      epsRDyn_head = CMPLX(0d0, 0d0)
      epsRDyn_head_temp = CMPLX(0d0, 0d0)
      do ifreq = 1, pol%nfreq
        freq_grp_ind = mod(ifreq-1, pol%nfreq_group)
        ifreq_para=(ifreq+pol%nfreq_group-1)/pol%nfreq_group
        if(freq_grp_ind .eq. peinf%rank_mtxel) then
          epsRDyn_head_temp(ifreq)=epsRDyn(1, 1, ifreq_para)
        endif
      enddo
#ifdef MPI
      call MPI_REDUCE(epsRDyn_head_temp(1), epsRDyn_head(1), pol%nfreq, MPI_COMPLEX_DPC, MPI_SUM, 0, &
       peinf%mtxel_comm, mpierr)
#endif
    endif ! (nfreq_group_alloc == pol%nfreq_in_group)

  endif

  if (peinf%inode == 0 .and. pol%freq_dep == 2) then
    write(51, '("# q= ",3f12.6, " nmtx= ",i0)') q0(:), pol%nmtx
    write(51, *)
    do iw = 1, pol%nfreq
      write(51, '(f12.6, 4(1x, es15.6e3))') pol%dFreqGrid(iw), &
        dble(epsRDyn_head(iw)), aimag(epsRDyn_head(iw))
    enddo
  endif

  if (peinf%inode .eq. 0 .and. pol%freq_dep .eq. 0) then
  
    ! JRD Warn User about possible lack of symmetry
    if (is_q0) then
      write(7, *)
      write(7, *) 'For q0 points, you should check the symmetry (eps(G, G'') = eps*(-G, -G'')) by'
      write(7, *) 'using the eps0sym code. Wavefunction convergence, as well as a finite q-shift'
      write(7, *) 'may cause this property of eps(G, G'') to be broken.'
      write(6, *)
      write(6, *) 'For q0 points, you should check the symmetry (eps(G, G'') = eps*(-G, -G'')) by'
      write(6, *) 'using the eps0sym code. Wavefunction convergence, as well as a finite q-shift'
      write(6, *) 'may cause this property of eps(G, G'') to be broken.'
    endif
    
    write(7, 4000) kp%nspin
    do i = 1, scal%npr
      is = scal%isrtxrow(i)
      do j = 1, scal%npc
        js = scal%isrtxcol(j)
        if (i .eq. j .or. i .eq. j+1) then
          write(7, 4200) gvec%components(1:3, is), gvec%components(1:3, js), eps(i, j)
        endif
      end do
    end do
  end if

! JRD the i and j loops are out of order below
  if (peinf%inode == 0 .and. pol%freq_dep > 0) then
    write(7, 4001) kp%nspin
    do iw = 1, min(pol%nfreq_in_group, nfreq_group_alloc)
      do i = 1, scal%npr
        is = scal%isrtxrow(i)
        do j = 1, scal%npc
          js = scal%isrtxcol(j)
          if (i .eq. j .or. i .eq. j+1) then
            write(7, 4300) gvec%components(1:3, is), gvec%components(1:3, js), epsRDyn(i, j, iw)
          endif
        end do
      end do
    end do
  end if

4000 format(/ /,13x, 'g',19x, 'gp',9x, &
       'inverse epsilon           nspin= ',1i1)
4001 format(/ /,13x, 'g',19x, 'gp',9x, &
#ifdef CPLX
       'inverse epsilon RDyn/ADyn nspin= ',1i1)
#else
  'inverse epsilon RDyn      nspin= ',1i1)
#endif

4200 format(5x, 3i5, 5x, 3i5, 5x, 2f13.8)
4300 format(5x, 3i5, 5x, 3i5, 5x, 4f13.8)

!---------- Full-Frequency Sum-Rule-----------------------------------------

  epssum1R = 0D0
  epssum2R = 0D0

#ifdef CPLX
  epssum1A = 0D0
  epssum2A = 0D0
#endif

  ! FHJ: These sum rules are*wrong*for systems without inversion symmetry.
  ! The sum rules apply to the Hermitian/Antihermitian components of the
  ! dielectric matrix. See D. J. Johnson, Phys. Rev. B 9, 4475 (1974), 
  ! eqns. 2.1 to 2.7. Note that the derivation of the complex GPP is right, 
  ! even if the notation is misleading (see comment in LITERATURE.html).
  ! Anyone care to fix this??
  if (peinf%inode == 0 .and. pol%freq_dep == 2 .and. pol%freq_dep_method /= 2) then
    do iw = 2, pol%nFreq-pol%nfreq_imag
      epssum1R = epssum1R+(1D0*Ryd/pol%dFreqGrid(iw))* &
        IMAG(epsRDyn_head(iw))*(pol%dFreqGrid(iw)-pol%dFreqGrid(iw-1))/Ryd
      epssum2R = epssum2R+(pol%dFreqGrid(iw)/Ryd)* &
        IMAG(epsRDyn_head(iw))*(pol%dFreqGrid(iw)-pol%dFreqGrid(iw-1))/Ryd
    enddo
    
    epssum1R=(2D0*epssum1R/Pi_D)+1D0
    epssum1R_rel=(epssum1R)/dble(epsRDyn(1, 1, 1))
    epssum2R_rel=(-1D0*epssum2R)/((Pi_D/2D0)*omega_plasma**2)
    
! Ref: Hybertsen & Louie PRB 34, 5390 (1986), eq. 29 and Appendix A
    write(6, *) ' '
    write(6, *) 'Full Frequency: Sum rules for head:'
    write(6, *) 'Int((1/w)*Im(eps^-1(w))) =', epssum1R_rel*100D0, ' % of exact'
    write(6, *) 'Int((w)*Im(eps^-1(w))) =', epssum2R_rel*100D0, ' % of exact'
  endif


!---------- Write inverse dielectric matrices to file-----------------------

  call timing%start(timing%io_total)
  call timing%start(timing%epsinv_i_o)
  if (peinf%inode == 0) write(6, '(/1x, a)') 'Writing dielectric matrix to file'
  !XXXXXXXX
  if (subspace) then
    if (keep_full_eps_static) then
      if (peinf%inode == 0) write(6, '(1x, a/)') 'Subspace: Full Epsilon will be retained for omega = 0' 
    else
      if (peinf%inode == 0) write(6, '(1x, a/)') 'Subspace: Full Epsilon will NOT be retained for omega = 0'
    end if
  end if
  !XXXXXXXX
  
  if (is_q0) then
    filename = 'eps0mat.h5'
  else
    filename = 'epsmat.h5'
  endif

#ifdef HDF5
  if (pol%use_hdf5) then
    my_iq = iq
    ! If this is not a q->0 point, write to (iq-nq0) q-point
    if (.not.is_q0) my_iq = iq-pol%nq0
    if (peinf%inode .eq. 0) then
      call write_gvec_indices_hdf(gvec%ng, pol%isrtx, isorti, ekin, my_iq, filename)
      call write_vcoul_hdf(vcoul, my_iq, filename)
    endif

! JRD: Write diagonal elements of Matrix

    if (pol%freq_dep .eq. 0) then 
      isize = SCALARSIZE
    else
      isize = 2
    endif

    SAFE_ALLOCATE(epsdiag, (isize, pol%nmtx, 1))
    epsdiag = 0D0

#ifdef USESCALAPACK
! DVF: Only the first frequency group contains omega = 0. With 
! no parallel frequencies, all processors are in the first frequency group.
    if(peinf%igroup_f .eq. 0) then
      SAFE_ALLOCATE(epsdiagt, (isize, pol%nmtx, 1))
      epsdiagt = 0D0
      do jj = 1, pol%nmtx
        icol = MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL), scal%npcol)
        if (icol .eq. scal%mypcol) then
          ii = jj
          irow = MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL), scal%nprow)
          if (irow .eq. scal%myprow) then
            if (pol%freq_dep .eq. 0) then
              epsdiagt(1, jj, 1)=dble(eps(scal%imyrowinv(jj), scal%imycolinv(jj)))
#ifdef CPLX
              epsdiagt(2, jj, 1)=IMAG(eps(scal%imyrowinv(jj), scal%imycolinv(jj)))
#endif
            else
              epsdiagt(1, jj, 1)=dble(epsRDyn(scal%imyrowinv(jj), scal%imycolinv(jj), 1))
              epsdiagt(2, jj, 1)=IMAG(epsRDyn(scal%imyrowinv(jj), scal%imycolinv(jj), 1))
            endif
          endif
        endif
      enddo
      if(pol%nfreq_group .gt. 1) then
        call MPI_reduce(epsdiagt, epsdiag, isize*pol%nmtx, MPI_DOUBLE_PRECISION, MPI_SUM, 0, peinf%freq_comm, mpierr)
      else
        call MPI_reduce(epsdiagt, epsdiag, isize*pol%nmtx, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      endif
      SAFE_DEALLOCATE(epsdiagt)
    endif
#else
    if (peinf%inode .eq. 0) then
      do jj = 1, pol%nmtx
        if (pol%freq_dep == 0) then
          epsdiag(1, jj, 1) = dble(eps(jj, jj))
#ifdef CPLX
          epsdiag(2, jj, 1) = IMAG(eps(jj, jj))
#endif
        else
          epsdiag(1, jj, 1) = dble(epsRDyn(jj, jj, 1))
          epsdiag(2, jj, 1) = IMAG(epsRDyn(jj, jj, 1))
        endif
      enddo
    endif    
#endif

    if (peinf%inode .eq. 0) then
      call write_matrix_diagonal_hdf(epsdiag, pol%nmtx, my_iq, isize, filename)
    endif
    SAFE_DEALLOCATE(epsdiag)

  else
#endif
    iunit = 13
    if (is_q0) iunit = 12
    if (peinf%inode .eq. 0) then
      write(iunit) gvec%ng, pol%nmtx, &
        (pol%isrtx(i), isorti(i), i = 1, gvec%ng)
      write(iunit) (ekin(i), i = 1, gvec%ng)
      write(iunit) (q0(i), i = 1, 3)
    endif
#ifdef HDF5
  endif
#endif
! WK add option for writing epsmat
  if (pol%freq_dep == 0) then
#ifdef HDF5
    if (pol%use_hdf5) then
      if (.not.pol%skip_epsmat_write) then
        !write full epsmat
        !> WK: Do serial io
        call write_matrix_d_hdf(scal, eps, pol%nmtx, my_iq, 1, filename)
        !> WK: Do serial io
!#ifdef USESCALAPACK
!        call write_matrix_d_par_hdf(scal, eps, pol%nmtx, my_iq, 1, filename)
!#else
!        call write_matrix_d_hdf(scal, eps, pol%nmtx, my_iq, 1, filename)
!#endif
      else
        !skip writing
      endif

    else !(.not.pol%use_hdf5) 
#endif
      if (.not.pol%skip_epsmat_write) then
        call write_matrix_d(scal, eps, pol%nmtx, iunit)
      else
        !skip writing
      
      endif
#ifdef HDF5
    endif
#endif

  else  ! freq_dep /= 0 below

    IF(subspace .AND. pol%matrix_in_subspace_basis) THEN
      ! here we write out the subspace epsinv matrices
#ifdef HDF5
      if (pol%use_hdf5) then
        ! here HDF5
        !XXX
        !XXX call die("WRITING SUBSPACE EPSINV in HDF5 not yet implemented")  
        !XXX
        ! for now subspace only work in combination with MPI/SCALAPACK 
        ! use existing routine to write to H5 format
        ! vcoul has already been written
#ifdef USESCALAPACK
        if ( keep_full_eps_static ) then
          call cpu_time(t1)
          if (peinf%inode == 0) write(6, '(1x, a)') 'Writing Full Inv Epsilon for omega = 0 (H5)'
          eps1Aux(:,:) = epsRDyn(:,:,1)
          call write_matrix_d_par_hdf_sub(scal, eps1Aux, pol%nmtx, pol%nmtx, my_iq, 1, filename, &
                                          'mats/matrix_fulleps0')
          call cpu_time(t2)
          if (peinf%inode == 0) write(6, '(1x, A, F12.2, A)') 'Done:', t2-t1, ' s'
        end if
        ! eigenvectors, here we use an auxiliary scalapack environment, such that 
        ! we`ll write only the needed eigenvectors (this maybe also useful if we want 
        ! to change the scalapack layout)
        call cpu_time(t1)
        if (peinf%inode == 0) write(6, '(1x, a)') 'Writing Eigenvectors of Epsilon omega = 0 (H5)'
        scal_aux%nprow   = scal%nprow
        scal_aux%npcol   = scal%npcol
        scal_aux%myprow  = scal%myprow
        scal_aux%mypcol  = scal%mypcol
        scal_aux%nbl     = scal%nbl
        scal_aux%icntxt  = scal%icntxt
        scal_aux%npr     = scal%npr
        scal_aux%npc     = numroc(neig_sub, scal%nbl, scal%mypcol, 0, scal%npcol)
        ! copy whole data in local buffer
        eps1Aux(:,:) =  pol%eigenvect_omega0(:,:)
        if(.true.) then
          call write_matrix_d_par_hdf_sub(scal_aux, eps1Aux, pol%nmtx, neig_sub, my_iq, 1, filename, &
                                          'mats/matrix_eigenvec')
        else
          call write_matrix_d_par_hdf_sub(scal, eps1Aux, pol%nmtx, pol%nmtx, my_iq, 1, filename, &
                                          'mats/matrix_eigenvec')
        end if
        call cpu_time(t2)
        if (peinf%inode == 0) write(6, '(1x, A, F12.2, A)') 'Done:', t2-t1, ' s'

        ! write subspace inverse matrix
        call cpu_time(t1)
        if (peinf%inode == 0) write(6, '(1x, a)') 'Writing Subspace Inv Epsilon for all other frequencies (H5)'
        if (peinf%inode == 0) write(6, '(1x, a, i5)') 'Number of frequencies:', pol%nFreq
        scal_sub%nprow   = scal%nprow
        scal_sub%npcol   = scal%npcol
        scal_sub%myprow  = scal%myprow
        scal_sub%mypcol  = scal%mypcol
        scal_sub%nbl     = scal%nbl
        scal_sub%icntxt  = scal%icntxt
        scal_sub%npr     = nrow_loc_sub
        scal_sub%npc     = ncol_loc_sub
        call write_matrix_f_par_hdf(scal_sub, pol%nFreq_in_group, epsRDyn_sub, &
                                    neig_sub, my_iq, 1, filename, pol%nfreq_group, &
                                    dset_name='mats/matrix_subspace')
        call cpu_time(t2)
        if (peinf%inode == 0) write(6, '(1x, A, F12.2, A)') 'Done:', t2-t1, ' s'
#else
        call die("WRITING SUBSPACE EPSINV in HDF5 only works in combination with SCALAPACK")
#endif
      else  ! pol%use_hdf5 within subspace
#endif
        ! here no HDF5
        ! write the coulomb potential
        if (peinf%inode .eq. 0) then
           ! already written: 
           !  gvec%ng, pol%nmtx, pol%isrtx, isorti, ekin, q0
           ! Write coulomb potential
           write(iunit) (vcoul(i), i = 1, pol%nmtx)
        end if
        ! if we want to use the full epsinv at omega zero write it out here
        IF(keep_full_eps_static) THEN
          if (peinf%inode == 0) write(6, '(1x, a)') 'Writing Full Inv Epsilon for omega = 0'
          eps1Aux(:,:) = epsRDyn(:,:,1)
          call write_matrix_d_sub(scal, eps1Aux, pol%nmtx, iunit)
        END IF
        ! write eigenvectors (just consider it as a square matrix)
        if (peinf%inode == 0) write(6, '(1x, a)') 'Writing Eigenvectors of Epsilon omega = 0'
        ! write the basis size here
        if (peinf%inode == 0) then
          write(iunit)  neig_sub
        end if
        eps1Aux(:,:) = pol%eigenvect_omega0(:,:)
        call write_matrix_d_sub(scal, eps1Aux, pol%nmtx, iunit, neig_sub)
        ! and all subspace matrices (including that for omega = zero), 
        ! we don`t write the advace, just recalculate it in sigma as done in the case of use_hdf5
        ! create the sub scalapack environment
        if (peinf%inode == 0) write(6, '(1x, a)') 'Writing Subspace Inv Epsilon for all other frequencies'
        if (peinf%inode == 0) write(6, '(1x, a, i5)') 'Number of frequencies:', pol%nFreq
        scal_sub%nprow   = scal%nprow
        scal_sub%npcol   = scal%npcol
        scal_sub%myprow  = scal%myprow
        scal_sub%mypcol  = scal%mypcol
        scal_sub%nbl     = scal%nbl
        scal_sub%icntxt  = scal%icntxt
        scal_sub%npr     = nrow_loc_sub
        scal_sub%npc     = ncol_loc_sub
        call write_matrix_f(scal_sub, pol%nFreq, epsRDyn_sub, neig_sub, iunit, pol%nfreq_group)
#ifdef HDF5
      endif
#endif
    ELSE  ! else (write subspace matrices)
#ifdef HDF5
      if (pol%use_hdf5) then
#ifdef USESCALAPACK
        call write_matrix_f_par_hdf(scal, pol%nFreq_in_group, epsRDyn, &
          pol%nmtx, my_iq, 1, filename, pol%nfreq_group)
#else
        call write_matrix_f_hdf(scal, pol%nFreq, epsRDyn, pol%nmtx, my_iq, 1, filename)
#endif
      else
#endif
        call write_matrix_f(scal, pol%nFreq, epsRDyn, pol%nmtx, iunit, pol%nfreq_group &
#ifdef CPLX
          ,advanced = epsADyn&
#endif
          )
#ifdef HDF5
      endif
#endif
    END IF  ! write subspace matrices
  endif

  if (peinf%inode == 0) write(6, '(1x, a/)') 'Ok'
  call timing%stop(timing%epsinv_i_o)
  call timing%stop(timing%io_total)

! Finished writing eps
!------------------------------------------------------------------------

  SAFE_DEALLOCATE(vcoul)
  SAFE_DEALLOCATE(isorti)

  if(pol%freq_dep == 0) then
    SAFE_DEALLOCATE(eps)
    SAFE_DEALLOCATE(ewng)
  else
    SAFE_DEALLOCATE(epsRDyn)
    SAFE_DEALLOCATE(eps1Aux)
    SAFE_DEALLOCATE(chiRDyntmp)
#ifdef CPLX
    if (.not.pol%use_hdf5) then
      SAFE_DEALLOCATE(epsADyn)
#ifdef USESCALAPACK
      SAFE_DEALLOCATE(wcoul_r)
      SAFE_DEALLOCATE(vcoul_dist)
      SAFE_DEALLOCATE(vcoul_inv_dist)
#endif
    endif
#endif
  endif
  if (pol%freq_dep == 2) then
    if (pol%nfreq_group == 1 .and. peinf%rank_f == 0) then
      SAFE_DEALLOCATE(epsRDyn_head)
    elseif (pol%nfreq_group > 1 .and. peinf%inode < peinf%npes .and. peinf%rank_f == 0) then
      SAFE_DEALLOCATE(epsRDyn_head)
      SAFE_DEALLOCATE(epsRDyn_head_temp)
    endif
  endif
  
  POP_SUB(epsinv)
  
  return
end subroutine epsinv

#if defined MPI && defined USESCALAPACK
subroutine invert_subspace(nrow_loc_sub, ncol_loc_sub, neig_sub, nmtx, &
                           scal, pol, vcoul, eps1Aux_sub, eps1Aux, iw, integ_rpa_pol)
  integer, intent(in):: nrow_loc_sub, ncol_loc_sub, neig_sub, nmtx
  type (scalapack), intent(in):: scal 
  type (polarizability), intent(in):: pol
  real(DP), intent(in):: vcoul(nmtx) 
  complex(DPC), intent(inout):: eps1Aux_sub(nrow_loc_sub, ncol_loc_sub)
  complex(DPC), intent(out):: eps1Aux(scal%npr, scal%npc)
  integer, intent(in):: iw
  real(DP), intent(out):: integ_rpa_pol(pol%nfreq_imag)

  type (scalapack):: scal_sub
  complex(DPC), allocatable:: C_Pgemm(:,:)
  integer:: desca(9), desc_sub(9), info
  integer:: i, j, irow, icol, icurr, irowm, icolm
  real(DP):: integ_rpa_val

  PUSH_SUB(invert_subspace)

  ! invert matrix 
  ! create the subspace scalapack environment
  scal_sub%nprow   = scal%nprow
  scal_sub%npcol   = scal%npcol
  scal_sub%myprow  = scal%myprow
  scal_sub%mypcol  = scal%mypcol
  scal_sub%nbl     = scal%nbl
  scal_sub%icntxt  = scal%icntxt
  scal_sub%npr     = nrow_loc_sub
  scal_sub%npc     = ncol_loc_sub
  call timing%start(timing%epsinv_omega_neq_0)
  if(pol%do_rpa) then
    call zinvert_with_scalapack(neig_sub, scal_sub, eps1Aux_sub, integ_rpa_val)
  else
    call zinvert_with_scalapack(neig_sub, scal_sub, eps1Aux_sub)
  endif
  call timing%stop(timing%epsinv_omega_neq_0)

  if(iw > pol%nFreq-pol%nfreq_imag .and. pol%do_rpa) then
    integ_rpa_pol(iw-(pol%nFreq-pol%nfreq_imag)) = integ_rpa_val
  endif

  ! subtract one from diagonal 
  icurr = 0
  do i = 1, neig_sub
    irow = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL), scal%nprow)
    if(irow .ne. scal%myprow) cycle
    do j = 1, neig_sub
      icol = MOD(INT(((j-1)/scal%nbl)+TOL_SMALL), scal%npcol)
      if(icol .eq. scal%mypcol) then
        icurr = icurr+1
        irowm = INT((icurr-1)/ncol_loc_sub+TOL_SMALL)+1
        icolm = MOD((icurr-1), ncol_loc_sub)+1

        IF(i == j) eps1Aux_sub(irowm, icolm) = eps1Aux_sub(irowm, icolm) - ONE

      endif
    end do
  end do

  if (pol%keep_full_eps_static .or. .not.pol%matrix_in_subspace_basis) then
    call descinit(desca, nmtx, nmtx, scal%nbl, scal%nbl, 0, 0, &
                  scal%icntxt, max(1, scal%npr), info)
    if(info < 0) then
      write(0, '(a, i3, a)') 'Argument number ', -info, ' had an illegal value on entry.'
      call die("descinit error for descaA in subspace inversion")
    else if(info > 0) then
      write(0, *) 'info = ', info
      call die("descinit error for descaA in subspace inversion")
    endif

    call descinit(desc_sub, neig_sub, neig_sub, scal%nbl, scal%nbl, 0, 0, &
                  scal%icntxt, max(1, nrow_loc_sub), info)
    if(info < 0) then
      write(0, '(a, i3, a)') 'Argument number ', -info, ' had an illegal value on entry.'
      call die("descinit error for desca_sub in subspace inversion")
    else if(info > 0) then
      write(0, *) 'info = ', info
      call die("descinit error for desca_sub in subspace inversion")
    endif

    ! go back to original basis
    SAFE_ALLOCATE(C_Pgemm, (scal%npr, scal%npc))
    call timing%start(timing%subspace_pgemm)
    CALL pzgemm('N','N', nmtx, neig_sub, neig_sub, (1.0d0, 0.0d0), pol%eigenvect_omega0, 1, 1, desca, &
                eps1Aux_sub(:,:), 1, 1, desc_sub, (0.0d0, 0.0d0), &
                C_Pgemm, 1, 1, desca)
    CALL pzgemm('N','C', nmtx, nmtx, neig_sub, (1.0d0, 0.0d0), C_Pgemm, 1, 1, desca, &
                pol%eigenvect_omega0, 1, 1, desca, (0.0d0, 0.0d0), &
                eps1Aux, 1, 1, desca)
    call timing%stop(timing%subspace_pgemm)
    SAFE_DEALLOCATE(C_Pgemm)

    ! restore one on diagonal
    icurr = 0
    do i = 1, nmtx

      irow = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL), scal%nprow)
      if(irow .ne. scal%myprow) cycle

      do j = 1, nmtx
        icol = MOD(INT(((j-1)/scal%nbl)+TOL_SMALL), scal%npcol)
        if(icol .eq. scal%mypcol) then
          icurr = icurr+1
          irowm = INT((icurr-1)/scal%npc+TOL_SMALL)+1
          icolm = MOD((icurr-1), scal%npc)+1

          IF(i == j) eps1Aux(irowm, icolm) = eps1Aux(irowm, icolm) + ONE
        endif
      end do
    end do
  else
    eps1Aux(:,:) = (0.0d0, 0.0d0)
    ! here we copy the head of the inverse eps in subspace basis, so it will not 
    ! appear zero on standard output (restore one on diagonal to be consistent)
    eps1Aux(1, 1) = eps1Aux_sub(1, 1) + ONE
    if(nrow_loc_sub > 1 .and. ncol_loc_sub > 1) then
      eps1Aux(2, 2) = eps1Aux_sub(2, 2) + ONE
    end if
  endif

  POP_SUB(invert_subspace)

end subroutine invert_subspace
#endif

end module epsinv_m
