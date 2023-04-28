!===============================================================================
!
! Routines:
!
! (1)  inteqp         Originally by JRD       Last Edited: 10/11/2010 (JRD)
!
! Extrapolates Eqp corrections from the coarse grid to the fine grid by
! a wavefunction-based plus linear interpolation scheme that preserves band crossings/character. 
! The code reads eqp_co.dat and writes eqp.dat and
! eqp_q.dat. The DFT wavefunctions for the interpolation are read
! from WFN_co, WFN_fi, and WFNq_fi files.
!
! See absorption.inp for options and keywords.
!
!================================================================================

#include "f_defs.h"

program inteqp

  use bse_init_m
  use fullbz_m
  use global_m
  use input_fi_m
  use input_q_m
  use inread_m
  use intwfn_m
  use timing_m, only: timing => bse_timing
  use write_program_header_m
  implicit none

  type (crystal) :: crys
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (eqpinfo) :: eqp
  type (xctinfo) :: xct
  type (flags) :: flag
  type (grid) :: kg_fi,kgq_fi,kg_co,kgq_co
  type (kpoints) :: kp_fi,kpq_fi
  type (int_wavefunction) :: intwfnc
  type (int_wavefunction) :: intwfnv

  integer :: ii,ncount
  integer :: iunit_c,iunit_v
  real(DP) :: vol,omega_plasma
  real(DP) :: tsec(2),tmin(2),tmax(2)
  
  character*16, allocatable :: routnam(:)
  integer, allocatable :: routsrt(:)
  integer, allocatable :: fi2co_wfn(:,:),indexq_fi(:)
  real(DP), allocatable :: kco(:,:)
  character :: filename*20
  SCALAR, allocatable :: dcc(:,:,:,:,:),dvv(:,:,:,:,:)
  real(DP), allocatable :: intp_coefs(:,:)

  call peinfo_init()

!----------------------
! Initialize timer
  call timing%init()
  call timing%start(timing%total)

!---------------------------
! Write header

  call write_program_header('BSE/IntEqp', .false.)

!---------------------------
! Read inteqp.inp
  
  call logit('Calling inread_inteqp')
  call open_file(8,file='inteqp.inp',form='formatted',status='old')
  call inread(eqp,xct,flag)
  call close_file(8)

! FHJ: Initialize xct%nkpt_co and dimentionality of the problem
  call bse_init(xct,flag)

!--------------------------
! Read wavefunctions on the fine grid

  call logit('Calling input')
  call timing%start(timing%input)
  call input_fi(crys,gvec,kg_fi,kp_fi,syms,eqp,xct,flag,omega_plasma,.false.,intwfnc)

  vol = xct%nktotal*crys%celvol
  if (peinf%inode.eq.0) then
    write(6,*) ' '
    write(6,*) 'More Job Parameters: '
    write(6,'(a,f32.14,a)') ' Crystal volume = ',vol,' a.u.'
    write(6,*) 'Number of valence bands = ',xct%nvb_fi
    write(6,*) 'Number of cond. bands   = ',xct%ncb_fi
    write(6,*) 'Number of spins   = ',xct%nspin
    write(6,*) ' '
  endif
  call timing%stop(timing%input)

  SAFE_ALLOCATE(indexq_fi, (xct%nkpt_fi))
  SAFE_ALLOCATE(xct%indexq_fi, (xct%nkpt_fi))

!  if (flag%vm.ne.1.or.flag%dtm.ne.1) then ! both are always 0 in this code --DAS
  call timing%start(timing%input_q)
  call logit('Calling input_q')
  call input_q(kp_fi,crys,gvec,kg_fi,kgq_fi,kpq_fi,syms,xct,indexq_fi,eqp,flag,intwfnv)
  call timing%stop(timing%input_q)

!------------------------------
! Calculate the transformation matrices from coarse grid wavefunctions
! FHJ: These are the final transformation coefs that will be used to interpolate
! the kernel. However, we might use an unrestricted version of dvv/dcc to
! interpolate eqp if xct%unrestricted_transf==.true..
  SAFE_ALLOCATE(dvv, (xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
  SAFE_ALLOCATE(dcc, (xct%ncb_fi,xct%n2b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
  SAFE_ALLOCATE(kco, (3,xct%nkpt_co))
  SAFE_ALLOCATE(fi2co_wfn, (xct%npts_intp_kernel,xct%nkpt_fi))
  SAFE_ALLOCATE(intp_coefs, (xct%npts_intp_kernel, xct%nkpt_fi))

  call logit('Calling intwfn')
  call timing%start(timing%intwfn)
  call intwfn(kp_fi,crys,syms,xct,flag,gvec,kg_fi,kgq_fi,kg_co,kgq_co,dcc,dvv,&
    kco,fi2co_wfn,indexq_fi,eqp,intwfnv,intwfnc,intp_coefs)
  call timing%stop(timing%intwfn)

  call kp_fi%free()
  call kpq_fi%free()
  call dealloc_grid(kg_fi)
  call dealloc_grid(kgq_fi)

  SAFE_DEALLOCATE_P(xct%ifmax)
  SAFE_DEALLOCATE_P(xct%ifmaxq)

  SAFE_DEALLOCATE_P(intwfnc%cgk)
  SAFE_DEALLOCATE_P(intwfnv%cgk)
  SAFE_DEALLOCATE_P(intwfnc%isort)
  SAFE_DEALLOCATE_P(intwfnv%isort)

  SAFE_DEALLOCATE_P(eqp%ecqp)
  SAFE_DEALLOCATE_P(eqp%evqp)
  
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

!--------------------------------
! Time accounting

  call timing%stop(timing%total)
  call timing%print()

  call write_memory_usage()

#ifdef MPI
  call MPI_FINALIZE(mpierr)
#endif
  
end program inteqp
