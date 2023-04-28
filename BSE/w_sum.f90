!=================================================================================
!
! Routines:
!
! (1) w_sum()         Originally By JRD       Last Modified 4/1/2012 (JRD)
!
!  Multiply Valence-Valence matrix elements by W to create temporary arrays
!  for the head, wing and body.
!
!  This routine scales as N^3, but is nested within the the igp loop 
!  in mtxel_kernel. Thus, it represents an N^4 step. Doing the multiplication
!  here is cheaper than doing it in the N^5 g_sum subroutine.
!
!=================================================================================

#include "f_defs.h"

module w_sum_m

  use accel_linalg_m
  use algos_kernel_m
  use blas_m
  use global_m

  implicit none

  public :: w_sum_cpu, w_sum_openacc, w_sum_omp_target

  private

contains

  subroutine w_sum_cpu(xct,wptcol,ofs1,ofs1p,n1,n1p,temph,tempw,tempb,m11p,indinvigp,ng_eps)
    type (xctinfo), intent(in) :: xct
    SCALAR, intent(in) :: wptcol(:)
    !> offset (i.e., add this number to map a local index to the global band index)
    integer, intent(in) :: ofs1, ofs1p
    !> number of bands for each wfn
    integer, intent(in) :: n1, n1p
    SCALAR, intent(inout) :: tempw(:,:,:,:), tempb(:,:,:,:), temph(:,:,:), m11p(:,:,:,:)
    integer, intent(in) :: indinvigp
    integer, intent(in) :: ng_eps

    SCALAR, allocatable :: m11p_conj(:,:,:)
    integer :: isv, ig, i1, i1p, gi1, gi1p
    logical :: use_omp

    PUSH_SUB(w_sum_cpu)

    ! NAG gives me a runtime error with threading here. Did not have time
    ! yet to find out if the error is real or a compiler bug.
#ifdef NAG
    use_omp = .false.
#else
    use_omp = .true.
#endif

    ! JRD: We allocate a new temporary array in order to get better cache performance
    SAFE_ALLOCATE(m11p_conj,(n1,n1p,xct%nspin))
    m11p_conj(:,:,:) = MYCONJG(m11p(indinvigp,:,:,:))

    do isv=1,xct%nspin
      if (indinvigp .eq. 1) then
        temph(ofs1+1:ofs1+n1, ofs1p+1:ofs1p+n1p, isv) = wptcol(1)*m11p_conj(1:n1, 1:n1p, isv)

        !$OMP PARALLEL PRIVATE(i1p, gi1p, i1, gi1, ig) DEFAULT(SHARED) IF(use_omp)
        do i1p = 1, n1p
          gi1p = ofs1p+i1p
          do i1 = 1, n1
            gi1 = ofs1+i1
            !$OMP DO
            do ig=2,ng_eps
              tempw(ig, gi1, gi1p, isv) = wptcol(ig) * m11p_conj(i1, i1p, isv)
            enddo
            !$OMP END DO NOWAIT
          enddo
        enddo
        !$OMP END PARALLEL

      else

        tempw(1, ofs1+1:ofs1+n1, ofs1p+1:ofs1p+n1p, isv) = tempw(1, ofs1+1:ofs1+n1, ofs1p+1:ofs1p+n1p, isv) + &
          wptcol(1)*m11p_conj(1:n1, 1:n1p, isv)

        !$OMP PARALLEL PRIVATE(i1p, gi1p, i1, gi1, ig) DEFAULT(SHARED) IF(use_omp)
        do i1p = 1, n1p
          gi1p = ofs1p+i1p
          do i1 = 1, n1
            gi1 = ofs1+i1
            !$OMP DO
            do ig=2,ng_eps
              tempb(ig, gi1, gi1p, isv) = tempb(ig, gi1, gi1p, isv) + &
                wptcol(ig) * m11p_conj(i1, i1p, isv)
            enddo
            !$OMP END DO NOWAIT
          enddo
        enddo
        !$OMP END PARALLEL

      endif

    enddo

    SAFE_DEALLOCATE(m11p_conj)

    POP_SUB(w_sum_cpu)

  end subroutine w_sum_cpu

  subroutine w_sum_openacc(xct, wptcol_mat, m11p_all, temph, tempw, tempb, indinv, ng_eps, &
                       phinv, nepsmin, ngpown_ipe, invband, nbands, ipe)
   type (xctinfo), intent(in) :: xct
   SCALAR, intent(in) :: wptcol_mat(:,:)
   SCALAR, intent(in), target :: m11p_all(:,:,:,:,:)
   SCALAR, intent(inout) :: tempw(:,:,:,:), tempb(:,:,:,:), temph(:,:,:)
   integer, intent(in) :: indinv(:), nbands(:)
   SCALAR, intent(in)  :: phinv(:)
   integer, intent(in) :: ng_eps, ngpown_ipe, nepsmin, invband, ipe
#ifdef OPENACC

   integer :: t1, t1p, t1_max
   integer :: ofs1, ofs1p
   integer :: n1, n1p
   SCALAR, pointer :: m11p(:,:,:,:)
   integer :: wfn_l(2), wfn_r(2)
   integer :: my_igp_loc, my_igp, t1_all

   SCALAR, allocatable :: m11p_conj(:,:,:), m11p_conj_mat(:,:,:,:,:)
   SCALAR, allocatable :: wm_mat(:,:,:,:,:)
   integer :: isv, ig, i1, i1p, gi1, gi1p
   integer :: indinvigp
   integer :: s1, s2, s3, s4, s5

   integer, allocatable :: all_igp(:,:)
   logical :: indinvigp_is_1
   integer :: indinvigp_pos_1, my_igp_loc_pos_1

    PUSH_SUB(w_sum_openacc)

    s1 = SIZE(m11p_all, 1)
    s2 = SIZE(m11p_all, 2)
    s3 = SIZE(m11p_all, 3)
    s4 = SIZE(m11p_all, 4)
    s5 = SIZE(m11p_all, 5)
    SAFE_ALLOCATE( m11p_conj_mat, (ngpown_ipe,s2,s3,s4,s5) )
    m11p_conj_mat = ZERO

    indinvigp_is_1  = .false.
    indinvigp_pos_1 = 1
    do my_igp_loc = 1, ngpown_ipe
      my_igp = INDXL2G(my_igp_loc, xct%nb, ipe, 0, peinf%npes)

      if (phinv(my_igp)==0) cycle
      if (my_igp > nepsmin) cycle

      indinvigp = indinv(my_igp)

      if ( indinvigp == 1 ) then
        indinvigp_is_1   = .true.
        indinvigp_pos_1  = indinvigp
        my_igp_loc_pos_1 = my_igp_loc
        cycle
      end if

      m11p_conj_mat(my_igp_loc,:,:,:,:) = MYCONJG( m11p_all(indinvigp,:,:,:,:) )

    end do

    SAFE_ALLOCATE(wm_mat, (ng_eps, s2,s3,s4,s5) )

    !$acc data present(tempb_accel, temph_accel, tempw_accel) copyin(m11p_conj_mat, m11p_all, wptcol_mat) create(wm_mat)

    call accel_xgemm('n', 'n', &
                   ng_eps, s2*s3*s4*s5, ngpown_ipe, &
                   ONE, &
                   wptcol_mat,    ng_eps, &
                   m11p_conj_mat, ngpown_ipe, &
                   ZERO, &
                   wm_mat, ng_eps, &
                   w_sum_algo)

    t1_max = 1
    if (xct%extended_kernel) t1_max = 2
    ! Loop over "generalized valence WFNs"
    t1_all = 0
    do t1=1,t1_max
      ofs1 = invband*(t1-1)
      n1 = nbands(t1)
      do t1p=1,t1_max
        t1_all = t1_all + 1
        ofs1p = invband*(t1p-1)
        n1p = nbands(t1p)
        ! This would be Mvvp within TDA
        wfn_l = (/t1,1/)  ! <t1,k|
        wfn_r = (/t1p,2/) ! |t1p,kp>

        !$acc parallel
        do isv=1,xct%nspin
          if ( indinvigp_is_1 ) then
            !$acc loop gang vector collapse(2)
            do i1p = 1, n1p
              do i1 = 1, n1
                gi1p = ofs1p+i1p
                gi1 = ofs1+i1
                temph_accel( gi1, gi1p, isv) = wptcol_mat(1,my_igp_loc_pos_1) * MYCONJG( m11p_all(indinvigp_pos_1, i1, i1p, isv, t1_all) )
              end do
            end do
            !$acc end loop

            !$acc loop gang vector collapse(3)
            do i1p = 1, n1p
              do i1 = 1, n1
                do ig=2,ng_eps
                  gi1p = ofs1p+i1p
                  gi1 = ofs1+i1
                  tempw_accel(ig, gi1, gi1p, isv) = tempw_accel(ig, gi1, gi1p, isv) + &
                                         wptcol_mat(ig,my_igp_loc_pos_1) * MYCONJG( m11p_all(indinvigp_pos_1, i1, i1p, isv, t1_all) )
                enddo
              enddo
            enddo
            !$acc end loop
          end if

          !$acc loop gang vector collapse(2)
          do i1p = 1, n1p
            do i1 = 1, n1
              gi1p = ofs1p+i1p
              gi1 = ofs1+i1
              tempw_accel(1, gi1, gi1p, isv) = tempw_accel(1, gi1, gi1p, isv) + &
                                         wm_mat(1, i1, i1p, isv, t1_all)
            end do
          end do
          !$acc end loop

          !$acc loop gang vector collapse(3)
          do i1p = 1, n1p
            do i1 = 1, n1
              do ig=2,ng_eps
                gi1p = ofs1p+i1p
                gi1 = ofs1+i1
                tempb_accel(ig, gi1, gi1p, isv) = tempb_accel(ig, gi1, gi1p, isv) + &
                                            wm_mat(ig, i1, i1p, isv, t1_all)
              enddo
            enddo
          enddo
          !$acc end loop

        end do
        !$acc end parallel

      end do ! t1p
    end do ! t1

    !$acc end data

    SAFE_DEALLOCATE( m11p_conj_mat )
    SAFE_DEALLOCATE(wm_mat)

#else
    PUSH_SUB(w_sum_openacc)

    call die("OpenACC version of w_sum requested, but OpenACC not compiled&
            & into this executable", only_root_writes = .true.)
#endif
    POP_SUB(w_sum_openacc)
  end subroutine w_sum_openacc

  subroutine w_sum_omp_target(xct, wptcol_mat, m11p_all, temph, tempw, tempb, indinv, ng_eps, &
                       phinv, nepsmin, ngpown_ipe, invband, nbands, ipe)
   type (xctinfo), intent(in) :: xct
   SCALAR, intent(in) :: wptcol_mat(:,:)
   SCALAR, intent(in), target :: m11p_all(:,:,:,:,:)
   SCALAR, intent(inout) :: tempw(:,:,:,:), tempb(:,:,:,:), temph(:,:,:)
   integer, intent(in) :: indinv(:), nbands(:)
   SCALAR, intent(in)  :: phinv(:)
   integer, intent(in) :: ng_eps, ngpown_ipe, nepsmin, invband, ipe
#ifdef OMP_TARGET

   integer :: t1, t1p, t1_max
   integer :: ofs1, ofs1p
   integer :: n1, n1p
   SCALAR, pointer :: m11p(:,:,:,:)
   integer :: wfn_l(2), wfn_r(2)
   integer :: my_igp_loc, my_igp, t1_all

   SCALAR, allocatable :: m11p_conj(:,:,:), m11p_conj_mat(:,:,:,:,:)
   SCALAR, allocatable :: wm_mat(:,:,:,:,:)
   integer :: isv, ig, i1, i1p, gi1, gi1p
   integer :: indinvigp
   integer :: s1, s2, s3, s4, s5

   integer, allocatable :: all_igp(:,:)
   logical :: indinvigp_is_1
   integer :: indinvigp_pos_1, my_igp_loc_pos_1

    PUSH_SUB(w_sum_omp_target)

    s1 = SIZE(m11p_all, 1)
    s2 = SIZE(m11p_all, 2)
    s3 = SIZE(m11p_all, 3)
    s4 = SIZE(m11p_all, 4)
    s5 = SIZE(m11p_all, 5)
    SAFE_ALLOCATE( m11p_conj_mat, (ngpown_ipe,s2,s3,s4,s5) )
    m11p_conj_mat = ZERO

    indinvigp_is_1  = .false.
    indinvigp_pos_1 = 1
    do my_igp_loc = 1, ngpown_ipe
      my_igp = INDXL2G(my_igp_loc, xct%nb, ipe, 0, peinf%npes)

      if (phinv(my_igp)==0) cycle
      if (my_igp > nepsmin) cycle

      indinvigp = indinv(my_igp)

      if ( indinvigp == 1 ) then
        indinvigp_is_1   = .true.
        indinvigp_pos_1  = indinvigp
        my_igp_loc_pos_1 = my_igp_loc
        cycle
      end if

      m11p_conj_mat(my_igp_loc,:,:,:,:) = MYCONJG( m11p_all(indinvigp,:,:,:,:) )

    end do

    SAFE_ALLOCATE(wm_mat, (ng_eps, s2,s3,s4,s5) )

    !$omp target data map(to: m11p_conj_mat, m11p_all, wptcol_mat) map(alloc: wm_mat)

    call accel_xgemm('n', 'n', &
                   ng_eps, s2*s3*s4*s5, ngpown_ipe, &
                   ONE, &
                   wptcol_mat,    ng_eps, &
                   m11p_conj_mat, ngpown_ipe, &
                   ZERO, &
                   wm_mat, ng_eps, &
                   w_sum_algo)

    t1_max = 1
    if (xct%extended_kernel) t1_max = 2
    ! Loop over "generalized valence WFNs"
    t1_all = 0
    do t1=1,t1_max
      ofs1 = invband*(t1-1)
      n1 = nbands(t1)
      do t1p=1,t1_max
        t1_all = t1_all + 1
        ofs1p = invband*(t1p-1)
        n1p = nbands(t1p)
        ! This would be Mvvp within TDA
        wfn_l = (/t1,1/)  ! <t1,k|
        wfn_r = (/t1p,2/) ! |t1p,kp>

        do isv=1,xct%nspin
          if ( indinvigp_is_1 ) then
            !$omp target teams loop default(none) &
            !$omp&     shared(n1p, n1, ofs1p, ofs1, isv, my_igp_loc_pos_1, indinvigp_pos_1, t1_all) &
            !$omp&     shared(temph_accel, wptcol_mat, m11p_all) &
            !$omp&     private(i1p, i1, gi1p, gi1) &
            !$omp& collapse(2)
            do i1p = 1, n1p
              do i1 = 1, n1
                gi1p = ofs1p+i1p
                gi1 = ofs1+i1
                temph_accel( gi1, gi1p, isv) = wptcol_mat(1,my_igp_loc_pos_1) * MYCONJG( m11p_all(indinvigp_pos_1, i1, i1p, isv, t1_all) )
              end do
            end do
            !$omp end target teams loop

            !$omp target teams loop default(none) &
            !$omp&     shared(n1p, n1, ng_eps, ofs1p, ofs1, isv, my_igp_loc_pos_1, indinvigp_pos_1, t1_all) &
            !$omp&     shared(tempw_accel, wptcol_mat, m11p_all) &
            !$omp&     private(i1p, i1, ig, gi1p, gi1) &
            !$omp& collapse(3)
            do i1p = 1, n1p
              do i1 = 1, n1
                do ig=2,ng_eps
                  gi1p = ofs1p+i1p
                  gi1 = ofs1+i1
                  tempw_accel(ig, gi1, gi1p, isv) = tempw_accel(ig, gi1, gi1p, isv) + &
                                         wptcol_mat(ig,my_igp_loc_pos_1) * MYCONJG( m11p_all(indinvigp_pos_1, i1, i1p, isv, t1_all) )
                enddo
              enddo
            enddo
            !$omp end target teams loop
          end if

          !$omp target teams loop default(none) &
          !$omp&     shared(n1p, n1, ofs1p, ofs1, isv, t1_all) &
          !$omp&     shared(tempw_accel, wm_mat) &
          !$omp&     private(i1p, i1, gi1p, gi1) &
          !$omp& collapse(2)
          do i1p = 1, n1p
            do i1 = 1, n1
              gi1p = ofs1p+i1p
              gi1 = ofs1+i1
              tempw_accel(1, gi1, gi1p, isv) = tempw_accel(1, gi1, gi1p, isv) + &
                                         wm_mat(1, i1, i1p, isv, t1_all)
            end do
          end do
          !$omp end target teams loop

          !$omp target teams loop default(none) &
          !$omp&     shared(n1p, n1, ng_eps, ofs1p, ofs1, isv, t1_all) &
          !$omp&     shared(tempb_accel, wm_mat) &
          !$omp&     private(i1p, i1, ig, gi1p, gi1) &
          !$omp& collapse(3)
          do i1p = 1, n1p
            do i1 = 1, n1
              do ig=2,ng_eps
                gi1p = ofs1p+i1p
                gi1 = ofs1+i1
                tempb_accel(ig, gi1, gi1p, isv) = tempb_accel(ig, gi1, gi1p, isv) + &
                                            wm_mat(ig, i1, i1p, isv, t1_all)
              enddo
            enddo
          enddo
          !$omp end target teams loop

        end do

      end do ! t1p
    end do ! t1

    !$omp end target data

    SAFE_DEALLOCATE( m11p_conj_mat )
    SAFE_DEALLOCATE(wm_mat)

#else
    PUSH_SUB(w_sum_omp_target)

    call die("OpenMP Target version of w_sum requested, but OpenMP Target not compiled&
            & into this executable", only_root_writes = .true.)
#endif
    POP_SUB(w_sum_omp_target)
  end subroutine w_sum_omp_target

end module w_sum_m
