!=================================================================================
!
! Routines:
!
! (1) nonlinearoptics(main)  Originally By JRD    Last Modified 4/24/2019 (GKA)
!
!     See the README file for more information and usage.
!
!=================================================================================

#include "f_defs.h"

program nonlinearoptics

#ifdef HDF5
  use hdf5
#endif
  use global_m
  use fullbz_m
  use genwf_m
  use mtxel_optical_m
  use readasvck_tp_m
  use readasvck_uf_m
  use timing_m, only: timing => extra_timing
  use absp3d_m
  use absp_tp_m
  use absp_uf_m
  use input_m
  use input_q_m
  use inread_m
  use write_program_header_m
  implicit none

  type (crystal) :: crys
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (eqpinfo) :: eqp
  type (xctinfo) :: xct
  type (flags) :: flag
  type (grid) :: kg_fi,kgq_fi
  type (wavefunction) :: wfnc_fi
  type (wavefunction) :: wfnvq_fi
  type (work_genwf) :: work, workq
  type (int_wavefunction) :: intwfn
  type (windowinfo) :: win
  type (otherinfo) :: other
  
  integer :: ii,ncount,ntim,nmat,nblock,neig,i
  integer :: ik,ikq,ikt,ikrq,ikcvs,ic,iv,is
  integer :: nmax,nk, unit
  integer :: error
  real(DP) :: vol
  real(DP) :: tsec(2),tmin(2),tmax(2)
  
  character*16, allocatable :: routnam(:)
  integer, allocatable :: indexq_fi(:)
  real(DP), allocatable :: k(:,:)
  integer, allocatable :: ikmax(:)
  
  SCALAR, allocatable :: &
    s1(:),s1k(:,:,:),dummy(:)
  
  real(DP), allocatable :: osctotal(:),enarray(:,:),enarraytp(:)
  real(DP), allocatable :: oscarray(:,:)
  character*20 :: filename
  
#ifdef HDF5
  call h5open_f(error)
#endif
  call peinfo_init()
  
  call timing%init()
  call timing%start(timing%total)
  
  call write_program_header('NonLinearOptics', .false.)

!-------------------------------
! Read nonlinearoptics.inp

  call logit('Calling inread')
  
  call open_file(8,file='nonlinearoptics.inp',form='formatted',status='old')
  if (peinf%inode .eq. 0) then
    write(6,*) ' '
  endif
  call inread(eqp,xct,flag,nmax,neig,win,other)
  call close_file(8)

! JRD: Dumb debugging

!      write(6,*) 'nkpt,ncband,nvband,nspin,neig',
!     > xct%nkpt_fi,xct%ncb_fi,xct%nvb_fi,xct%nspin,neig
!      write(6,*) 'Shift =', xct%shift
!      write(6,*) 'Pol =', xct%pol

  if(flag%vm.eq.2) then
    if (peinf%inode.eq.0) then
      write(0,*) 'WARNING: read_eps2_moments not supported in this code. Ignoring keyword.'
    endif
    flag%vm=0
  endif

!------------------------------------
! Read wavefunctions on the fine grid

  call logit('Calling input')
  call timing%start(timing%input)
  call input(crys,gvec,kg_fi,syms,eqp,xct,flag)
  
  if (peinf%inode.eq.0) then
    write(6,*) ' '
    write(6,*) 'Done input routine'
    write(6,*) ' '
  endif

! JRD: Dumb Debugging

!      write(6,*) 'nkpt,ncband,nvband,nspin,neig',
!     > xct%nkpt_fi,xct%ncb_fi,xct%nvb_fi,xct%nspin,neig

  nmat = xct%nkpt_fi*(xct%ncb_fi+xct%nvb_fi)*xct%nspin*xct%nspin*(xct%ncb_fi+xct%nvb_fi)
  neig = xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin

! JRD: Dumb Debugging

!      write(6,*) 'nmat=',nmat

  SAFE_ALLOCATE(k, (3,xct%nkpt_fi))
  SAFE_ALLOCATE(ikmax, (nmat))
  
  vol = xct%nktotal*crys%celvol
  if (peinf%inode.eq.0) then
    write(6,'(a,f32.14,a)') ' Crystal volume = ',vol,' a.u.'
    write(6,*) 'number of valence bands = ',xct%nvb_fi
    write(6,*) 'number of cond. bands   = ',xct%ncb_fi
    write(6,*) 'number of spins   = ',xct%nspin
    write(6,'(a,f7.4,a)') ' Broadening: ',xct%eta,' eV'
    write(6,*) ' '
  endif
  call timing%stop(timing%input)
  
  SAFE_ALLOCATE(indexq_fi, (xct%nkpt_fi))
  if (flag%vm.ne.1) then
    call timing%start(timing%input_q)
    call logit('Calling inputq')
    call input_q(crys,gvec,kg_fi,kgq_fi,syms,xct,indexq_fi,flag%bzq)
    call timing%stop(timing%input_q)
  endif

!-------------------------------
! Calculate the velocity (or momentum) matrix elements.
! Each PE calculate a small number of them. At the end, share data
!
! If flag%vm.eq.1, skip this part and just read the matrix elements
! from "vmtxel_nl".

  if (peinf%inode .eq.0) write(6,*) ' '
  call logit('Calculating v/p matels')
  
  call timing%start(timing%vmtxel)
  
  nblock=(xct%ncb_fi+xct%nvb_fi)*(xct%ncb_fi+xct%nvb_fi)

! JRD - I guess we define nmat again for safety?

  nmat= xct%nspin*xct%nkpt_fi*nblock
  
  SAFE_ALLOCATE(s1, (nmat))
  SAFE_ALLOCATE(s1k, (xct%ncb_fi+xct%nvb_fi,xct%ncb_fi+xct%nvb_fi,xct%nspin))
  s1= 0.d0
  s1k= 0.d0

  if (flag%vm.eq.0) then
    xct%iwriteint = 0    

    do ikt=1, peinf%ikt(peinf%inode+1)
      ik = peinf%ik(peinf%inode+1,ikt)
      ikq = indexq_fi(ik)
      ikrq = kg_fi%indr(ik)
      
      call genwf(crys,gvec,kg_fi,syms,wfnc_fi,ik,ik,xct%nspin,xct%ncb_fi,&
                 work,intwfn,xct%iwriteint,is_cond=.true.)
      
      call genwf(crys,gvec,kgq_fi,syms,wfnvq_fi,ik,ikq,xct%nspin,xct%nvb_fi,&
                 workq,intwfn,xct%iwriteint,is_cond=.false.)
      
      if (flag%opr.eq.0) then
        call mtxel_v(wfnc_fi,wfnvq_fi,gvec,xct%qshift,xct%nvb_fi+xct%ncb_fi,xct%nvb_fi+xct%ncb_fi,s1k)
        call die("WARNING!! We don't do velocity yet!")
      elseif (flag%opr.eq.1) then
        call mtxel_m(crys,wfnc_fi,wfnvq_fi,gvec,eqp,xct%pol,xct%ncb_fi+xct%nvb_fi,xct%ncb_fi+xct%nvb_fi,&
          s1k,ik,.false.,kg_fi%f(:,ik))
! JRD The below commented out code is BAD BAD for now.  It leads to NaN for matrix elements between the same band 
!            call mtxel_m(crys,wfnc_fi,wfnvq_fi,gvec,eqp,xct%pol,xct%ncb_fi+xct%nvb_fi,xct%ncb_fi+xct%nvb_fi,&
!              s1k,ik,.true.)
      endif
      do is=1,xct%nspin
        do ic=1,xct%ncb_fi+xct%nvb_fi
          do iv=1,xct%nvb_fi+xct%ncb_fi
            ikcvs= is + (iv - 1 + (ic - 1 + (ik - 1)*(xct%ncb_fi+xct%nvb_fi))*(xct%nvb_fi+xct%ncb_fi))*xct%nspin
!            if (peinf%inode.eq.0) then
!              write(6,*) s1k(1,1,1)
!            endif
            s1(ikcvs) = s1k(ic,iv,is)
          enddo
        enddo
      enddo
      SAFE_DEALLOCATE_P(wfnc_fi%cg)
      SAFE_DEALLOCATE_P(wfnc_fi%isort)
      SAFE_DEALLOCATE_P(wfnvq_fi%cg)
      SAFE_DEALLOCATE_P(wfnvq_fi%isort)
    enddo

    ! typedefs initializes all of these ikolds to 0
    if(work%ikold.ne.0) then
      SAFE_DEALLOCATE_P(work%cg)
      SAFE_DEALLOCATE_P(work%ph)
      SAFE_DEALLOCATE_P(work%ind)
      SAFE_DEALLOCATE_P(work%isort)
    endif
    if(workq%ikold.ne.0) then
      SAFE_DEALLOCATE_P(workq%cg)
      SAFE_DEALLOCATE_P(workq%ph)
      SAFE_DEALLOCATE_P(workq%ind)
      SAFE_DEALLOCATE_P(workq%isort)
    endif

! JRD: Debug

    if (peinf%inode.eq.0) then
      write(6,*) s1(1),s1(2)
    endif

!----------------------------
! Share matrix elements

#ifdef MPI
    SAFE_ALLOCATE(dummy, (nmat))
    dummy = s1
    call MPI_ALLREDUCE(dummy,s1,nmat,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
    SAFE_DEALLOCATE(dummy)
#endif

    if (flag%vm /= 1) then
      write(filename,'(a,i4.4)') 'INT_VWFNQ_', peinf%inode
      unit = 128+(2*peinf%inode)+2
      call open_file(unit, filename, status='old')
      call close_file(unit, delete = .true.) ! files INT_VWFNQ_*
    endif

    if (flag%vm == 0) then
      write(filename,'(a,i4.4)') 'INT_CWFN_', peinf%inode
      unit = 128+(2*peinf%inode)+1
      call open_file(unit, filename, status='old')
      call close_file(unit, delete = .true.) ! files INT_CWFN_*
    endif
    
    if (peinf%inode.eq.0) then
      write(6,*) ' '
      write(6,*) 'writing matrix elements into vmtxel_nl'
      call open_file(16,file='vmtxel_nl',form='unformatted',status='replace')
      call open_file(17,file='vmtxel_nl.dat',status='replace')
      write(16) xct%nkpt_fi,xct%ncb_fi+xct%nvb_fi,xct%nvb_fi+xct%ncb_fi,xct%nspin,flag%opr
      write(17,*) xct%nkpt_fi,xct%ncb_fi+xct%nvb_fi,xct%nvb_fi+xct%ncb_fi,xct%nspin,flag%opr
      write(16) (s1(ikcvs),ikcvs=1,nmat)
      write(17,*) (s1(ikcvs),ikcvs=1,nmat)
      call close_file(16)
      call close_file(17)
    endif
    
  else
    
    if (peinf%inode.eq.0) then
      write(6,*) 'reading matrix elements from vmtxel_nl'
      call open_file(16,file='vmtxel_nl',form='unformatted',status='replace')
      read(16) ik,ic,iv,is,ii
      if (ik.ne.xct%nkpt_fi.or.ic.ne.(xct%ncb_fi+xct%nvb_fi).or. &
        iv.ne.(xct%nvb_fi+xct%ncb_fi) .or.is.ne.xct%nspin.or.ii.ne.flag%opr) then
        write(0,*) ik, ic, iv, is, ii
        call die('parameter mismatch in vmtxel_nl')
      endif
      
      read(16) (s1(ikcvs),ikcvs=1,nmat)
      call close_file(16)
    endif
#ifdef MPI
    call MPI_BCAST(s1,nmat,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#endif

  endif

  call dealloc_grid(kg_fi)
  call dealloc_grid(kgq_fi)
  
  call timing%stop(timing%vmtxel)
  
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
  
  SAFE_ALLOCATE(enarray, (nmat,win%nstates))
  SAFE_ALLOCATE(oscarray, (nmat,win%nstates))

! JRD Dumb Debugging

!      if (peinf%inode.eq.0) then
!        write(6,*) ' '
!        write(6,*) 'Test Optical Matrix Elements:'
!        write(6,*) s1(1),s1(2)
!        write(6,*) ' '
!      endif

  call timing%start(timing%readasvck)
  
  if (flag%job .eq. 0) then
    
    call readasvck_uf(win,nmat,neig,s1,xct,oscarray, &
      enarray,nk,ikmax,xct%nkpt_fi,k)
    
  elseif (flag%job .eq. 1) then

    ! Two-particles calculation

    SAFE_ALLOCATE(osctotal, (neig))
    SAFE_ALLOCATE(enarraytp, (neig))
    
    !if (peinf%inode .eq. 0) write(6,*) "About to call readasvck_tp"  ! DEBUG
    call readasvck_tp(win,nmat,neig,s1,xct, &
      osctotal,enarraytp,nk,k,xct%nkpt_fi)
    
    if (peinf%inode .eq. 0) then
      write(6,*) ' '
      write(6,*) 'Finished calculation of matrix elements'
      write(6,*) ' '
    endif
    
  end if
  
  call timing%stop(timing%readasvck)
  
  if (flag%job .eq. 0) then
    
    if (peinf%inode.eq.0) then
      call absp_uf(win,xct%eta,xct%nspin,xct%nspinor,neig,oscarray,enarray,vol,nmat,flag)
    endif
    
    if (flag%job .eq. 0) then
      if (peinf%inode.eq.0.and.other%ithreeD.eq.1) then
        call absp3d(xct%nspin,xct%nspinor,xct%eta,neig,oscarray(:,1), &
          enarray(:,1),vol,nmat,nk,ikmax,other%keta, &
          other%knx,other%kny,other%knz,k,flag)
      endif
    endif
    
    SAFE_DEALLOCATE(oscarray)
    
  else ! Two-Photon

    if (peinf%inode .eq. 0) then
      write(6,*) 'Starting Absorption.'
      call absp_tp(xct%nspin,xct%nspinor,xct%eta,neig,osctotal, &
        enarraytp,vol,neig,flag)
    endif
    
    SAFE_DEALLOCATE(osctotal)
    
  endif


!-------------------------------
! Time accounting

  call timing%print()
  
  call write_memory_usage()
  
#ifdef HDF5
  call h5close_f(error)
#endif

#ifdef MPI
  call MPI_FINALIZE(mpierr)
#endif
  
end program nonlinearoptics
