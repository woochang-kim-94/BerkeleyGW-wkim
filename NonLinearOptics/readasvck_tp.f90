!=============================================================================
!
! Routines:
!
! (1) readasvck_tp()    By JRD, GKA        Last Modified 4/24/2019 (GKA)
!
!     This routine reads in eigenvectors and calculates (by calling oscstrength)
!     the two-photon oscillator strength of a system.
!
!=============================================================================

#include "f_defs.h"

module readasvck_tp_m

  use global_m
  use oscstrength_m
  use timing_m, only: timing => extra_timing
  use evecs_m
  implicit none

  private
  public :: readAsvck_tp

contains

! Compute two-photon oscillator strength
subroutine readasvck_tp(win,nmat,neig,s1,xct,osctotal,enarray,nk,k,nkpt)
  type (windowinfo), intent(in) :: win
  integer, intent(in) :: nmat
  integer, intent(in) :: neig
  SCALAR, intent(in) :: s1(:) !< (nmat)
  type (xctinfo), intent(inout) :: xct
  real(DP), intent(out) :: osctotal(:) !< (neig), enarray(neig)
  real(DP), intent(out) :: enarray(:) !< (neig)
  integer, intent(out) :: nk
  integer, intent(in) :: nkpt
  real(DP), intent(out) :: k(:,:) !< (3,nkpt)

  integer :: ns, nc, nv, m, ic, iv, is
  integer :: nmatold, meig
  integer :: ikcvs,iimax
  integer :: inode
  integer :: ik, ieig, ieig_local, jeig, jeig_local
  integer :: nloop, nwritefreq
  real(DP) :: nmax
  SCALAR :: osc,osc0
  integer, allocatable :: imax(:,:),imax2(:,:)
  real(DP), allocatable :: kosc(:), osc0arrayt(:)
  real(DP), allocatable :: dummy0t(:)
  real(DP), allocatable :: max1(:,:),max2(:,:)
  SCALAR, allocatable :: Af(:,:,:,:), A1(:,:,:,:)
  SCALAR, allocatable :: osc0array(:),sumosc(:),sumosc2(:)
  type(evecs_t) :: evecs
  
  PUSH_SUB(readasvck_tp)

! -----------------------------------------

! Initialize 
  evecs%use_hdf5 = xct%use_hdf5

! Read and divide states
  if (peinf%inode .eq. 0) then
    write(6,*) ' '
    write(6,*) "------------------------------------------"
    write(6,*) "Two-photon calculation"
    write(6,*) ' '
    write(6,*) "Reading eigenvectors"
  end if
  call evecs%read_file(distribute=.true., neig_read=neig)
  call evecs%print_out_dimensions(6)


  ns = evecs%ns
  nv = evecs%nv
  nc = evecs%nc
  nk = evecs%nk
  nmatold = evecs%nmat
  meig = evecs%meig
  k(:,:) = evecs%kpts(:,:)
  enarray(:) = evecs%evals(:)

  ! JRD: Debug
  !if (peinf%inode .eq. 0) then
  !  write(6,*) " "
  !  write(6,*) 'From   A: ns,nv,nc,nk = ',ns,nv,nc,nk
  !  write(6,*) 'From inp: ns,nv,nc,nk = ',xct%nspin,xct%nvb_fi,xct%ncb_fi,xct%nkpt_fi
  !  write(6,*) 'neig, nmatold:', neig, nmatold
  !  write(6,*) 'Proc, meig: ', peinf%inode,meig
  !  write(6,*) ' '
  !end if
  
  SAFE_ALLOCATE(osc0arrayt, (nmatold))
  SAFE_ALLOCATE(kosc, (nk))
  SAFE_ALLOCATE(Af, (ns,nv,nc,nk))
  SAFE_ALLOCATE(A1, (ns,nv,nc,nk))
  SAFE_ALLOCATE(osc0array, (meig))
  
  !-------------------------------
  ! Calculate osctrength from |0>

  osc0arrayt = 0D0
  do ieig_local = 1, evecs%meig

    osc0 = 0D0
    do ik=1,nk
      do ic=1,nc
        do iv=1,nv
          do is=1,ns
            ikcvs = is + ns * ( &
              iv - 1 + (ic + nv - 1 + (ik - 1) * (nv + nc)) * (nv+nc))
            osc0 = osc0 + MYCONJG(evecs%Avc(is,nv-iv+1,ic,ik,ieig_local)) * s1(ikcvs)
          enddo
        enddo
      enddo
    enddo
    
    ! Check whether correct phase
    ieig = evecs%global_ieigs(ieig_local)
    osc0array(ieig_local) = osc0
    osc0arrayt(ieig) = abs(osc0)**2
    
  enddo
  
! JRD: Share osc0 for debugging purposes.
#ifdef MPI
    SAFE_ALLOCATE(dummy0t, (nmatold))
    dummy0t = osc0arrayt
    call MPI_ALLREDUCE(dummy0t,osc0arrayt,nmatold,MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
    SAFE_DEALLOCATE(dummy0t)
#endif
  
  if (peinf%inode .eq. 0) then
    write(6,*) ' '
    write(6,*) 'Starting Calculation of Two-Photon Osc. Strength'
    write(6,*) ' '
  endif

!-------------------------------
! Loop over Final States
  
  SAFE_ALLOCATE(sumosc, (neig))
  SAFE_ALLOCATE(sumosc2, (neig))
  SAFE_ALLOCATE(max1, (peinf%npes,neig))
  SAFE_ALLOCATE(imax, (peinf%npes,neig))
  SAFE_ALLOCATE(max2, (peinf%npes,neig))
  SAFE_ALLOCATE(imax2, (peinf%npes,neig))
  
  max1=0d0
  imax=0
  
  sumosc = 0D0
  max1 = 0D0
  imax = 0D0
  
  if ( win%nstates > 0) then
    nloop = min(win%nstates,neig)
  else
    nloop=neig
  endif
  
  ! Frequency of progress info message
  nwritefreq = 1
  if (neig .ge. 100) nwritefreq = 10
  if (neig .ge. 1000) nwritefreq = 100
  
  do ieig = 1, nloop
    
    ! Write progress message
    if (mod(ieig,neig/nwritefreq) .eq. 0) then
      if (peinf%inode .eq. 0) then
        write(6,*) 'Starting Final State', ieig, neig
      endif
    endif

    ! Broadcast Final State
    if (evecs%who_owns(ieig) .eq. peinf%inode) then
      Af(:,:,:,:) = evecs%Avc(:,:,:,:,evecs%local_ieigs(ieig))
    endif

    call timing%start(timing%os_comm)

#ifdef MPI
    call MPI_BCAST(Af(1,1,1,1),neig,MPI_SCALAR,evecs%who_owns(ieig), &
      MPI_COMM_WORLD,mpierr)
#endif

    call timing%stop(timing%os_comm)
    call timing%start(timing%os_sums)

!-------------------------------
! Loop over states proc owns

    sumosc(ieig)=0d0

    do jeig_local = 1, meig

      jeig = evecs%global_ieigs(jeig_local)

      ! We need the term jeig=ieig for normalization including single-photon spectra.

      A1(:,:,:,:)=evecs%Avc(:,:,:,:,jeig_local)
      
      osc=0d0
      
      call oscstrength(nmat,nk,ns,xct,A1,Af,s1,osc)
      
      sumosc(ieig) = sumosc(ieig)+MYCONJG(osc)*osc0array(jeig_local) &
        / ((enarray(jeig) - (enarray(ieig)/2d0)) / ryd)

      nmax = abs(osc*osc0array(jeig_local) &
        / ((enarray(jeig) - (enarray(ieig)/2d0)) / ryd)) ** 2
      
      if (nmax .gt. max1(peinf%inode+1,ieig)) then
        max1(peinf%inode+1,ieig) = nmax
        imax(peinf%inode+1,ieig) = jeig
      endif
      
    enddo

    call timing%stop(timing%os_sums)

  enddo

  call evecs%free()

!-------------------------------
! Now collapse osctotal

#ifdef MPI
  sumosc2 = sumosc
  max2=max1
  imax2=imax
  call MPI_ALLREDUCE(sumosc2,sumosc,neig, &
    MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD,mpierr)
  call MPI_ALLREDUCE(max2(1,1),max1(1,1),neig*peinf%npes, &
    MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
  call MPI_ALLREDUCE(imax2(1,1),imax(1,1),neig*peinf%npes, &
    MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

! Write osc.dat file
  if (peinf%inode .eq. 0) then
    call open_file(11,file='osc.dat',form='formatted',status='replace')
  endif

  do ieig = 1, neig
    
    nmax = 0d0
    
! Find highest contributing middle state
! and square osctotal
    do inode = 1,peinf%npes
      if (max1(inode,ieig) .gt. nmax) then
        nmax = max1(inode,ieig)
        iimax = imax(inode,ieig)
      endif
    enddo
    
    osctotal(ieig)=abs(sumosc(ieig))**2
    
    if (peinf%inode .eq. 0) then
      write(11,'(i8,3e16.8,i8,e16.8)') ieig,enarray(ieig),osc0arrayt(ieig),osctotal(ieig),iimax,nmax
    endif
    
  enddo
  if (peinf%inode .eq. 0) then
    call close_file(11)
  endif
  
  SAFE_DEALLOCATE(sumosc2)
  SAFE_DEALLOCATE(sumosc)
  SAFE_DEALLOCATE(A1)
  SAFE_DEALLOCATE(Af)
  SAFE_DEALLOCATE(osc0array)
  SAFE_DEALLOCATE(osc0arrayt)
  SAFE_DEALLOCATE(max1)
  SAFE_DEALLOCATE(max2)
  SAFE_DEALLOCATE(imax)
  SAFE_DEALLOCATE(imax2)

  POP_SUB(readasvck_tp)

  return
end subroutine readasvck_tp

end module readasvck_tp_m
