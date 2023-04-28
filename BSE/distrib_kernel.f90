!==================================================================================
!
! Routines:
!
! (1) distrib_kernel()       Originally By MLT    Last Modified 01/10/2013 (FHJ)
!
!      nv = # of valence bands
!      nc = # of conduction bands
!      nk = # of k-points
!      np = # of MPI processes
!
! ! TODO: full_kernel calculations only works with parallelization over k^2
!
!==============================================================================

#include "f_defs.h"

module distrib_kernel_m

  use createpools_m
  use fftw_m
  use global_m
  use misc_m
  implicit none

  private

  public :: &
    distrib_kernel

contains


! FHJ: FIXME - pools is a terrible terminology here. We should refer to pools
! as the groups that are trivially distributed over (e.g., bands).
subroutine distrib_kernel(xct,ng,kg,kgq,gvec)
  type (xctinfo), intent(inout) :: xct
  integer, intent(in) :: ng
  type (grid), intent(in) :: kg
  type (grid), intent(in) :: kgq !< this is only for the not-working finite-q kernel
  type (gspace), intent(inout) :: gvec

  integer :: npools,ipool,ipoolrank,iholdperown
  integer :: nproc, nown, npown
  integer :: ik,ikp,kvp,jcp,iownmax,iownwfmax,iownkmax
  integer :: ii,kv,jc,ipe,ierr,npes_over_2nkpt2
  integer :: ijk,nv,nc,nq
  integer :: Nrod,Nplane,Nfft(3),dNfft(3),dkmax(3),nmpinode
  integer :: cmax,vmax
  real(DP) :: mem,rmem,rmem2,rmemtemp1,rmemtemp2
  real(DP) :: scale,dscale
  character*128 :: tmpstr
  character*16 :: tmpstr1,tmpstr2
  integer, allocatable :: iown(:)

  PUSH_SUB(distrib_kernel)
      
!-----------------------------
! Determine the Parallelization Scheme

  npes_over_2nkpt2 = peinf%npes/(2*(xct%nkpt_co)**2)

  if (peinf%inode==0) write(6,'()')
  if (npes_over_2nkpt2 .ge. xct%ncb_co**2 .or. xct%ilowmem .eq. 1) then

    xct%ivpar = 1
    xct%icpar = 1
    
    ! JIM: integer division can be done in two steps this prevents overflow of
    ! nkpt**2*ncband**2 and subsequent error in computing nproc
    nproc = peinf%npes/xct%nkpt_co**2
    nproc = nproc/xct%ncb_co**2
    
    if (nproc .eq. 0) then
      npools = 0
    else
      call createpools(xct%nvb_co,xct%nvb_co,nproc,npools,nown,npown)
    endif

    peinf%nck = (xct%nkpt_co*xct%ncb_co*xct%nvb_co)**2
    if (peinf%inode == 0) then
      write(6,'(1x,a)') 'Parallelization strategy: (k,c,v)^2'
      write(6,'(1x,a)') '- Distributing all (k,c,v)^2 transitions to all processors'
    endif

  else if (peinf%npes .ge. 2*(xct%nkpt_co)**2) then
    
    xct%ivpar = 0
    xct%icpar = 1
    nproc = peinf%npes/(xct%nkpt_co**2)
    
    call createpools(xct%ncb_co,xct%ncb_co,nproc,npools,nown,npown)
    
    peinf%nck = (xct%nkpt_co*xct%ncb_co)**2
    if (peinf%inode == 0) then
      write(6,'(1x,a)') 'Parallelization strategy: (k,c)^2'
      write(6,'(1x,a)') '- Breaking transitions into (v)^2 blocks'
      write(6,'(1x,a)') '- Distributing each block with (k,c)^2 transitions to blocks of processors'
    endif

  else

    xct%ivpar = 0
    xct%icpar = 0
    nproc = peinf%npes
    
    call createpools(xct%nkpt_co,xct%nkpt_co,nproc,npools,nown,npown)

    peinf%nck = (xct%nkpt_co)**2
    if (peinf%inode == 0) then
      write(6,'(1x,a)') 'Parallelization strategy: (k)^2'
      write(6,'(1x,a)') '- Breaking transitions into (v,c)^2 blocks'
      write(6,'(1x,a)') '- Distributing each block with (k)^2 transitions to blocks of processors'
    endif

  endif

  if (xct%extended_kernel.and.(xct%ivpar/=0.or.xct%icpar/=0)) then
    write(0,'(a)') 'Can only run extended kernel calculation with parallelization over (k)^2.'
    write(0,'(a)') 'For now, please make sure #proc <= nk^2.'
    call die('Extended kernel calculation requires parallelization over (k)^2.', &
      only_root_writes=.true.)
  endif
  
  if (npools == 0) then
    if (mod((xct%nkpt_co*xct%ncb_co*xct%nvb_co)**2,peinf%npes).eq.0) then
      peinf%nckpe = (xct%nkpt_co*xct%ncb_co*xct%nvb_co)**2 / peinf%npes
    else
      peinf%nckpe = (xct%nkpt_co*xct%ncb_co*xct%nvb_co)**2 / peinf%npes + 1
    endif
  else
    peinf%nckpe = nown * npown
  endif
  
  if (peinf%inode==0) then
    if (npools == 0) then
      write(6,'(1x,a)') 'Parallelizing using the low-memory option'
    else
      write(6,'(1x,a)') 'Parallelizing calculation with blocks and pools:'
      write(6,'(1x,a,i0)') '- Number of transitions per block: ', peinf%nck
      write(6,'(1x,a,i0)') '- Number of MPI ranks per block: ', nproc
      write(6,'(1x,a,i0)') '- Number of transitions per MPI rank: ', peinf%nckpe
    endif
    write(6,'()')
  endif

!---------------------------------
! Determine if the conditions for a decent load-balancing are met

  if (mod(peinf%nck,peinf%npes)/=0 .and. peinf%inode==0 .and. peinf%verb_high) then
    write(6,'(1x,a)') 'NOTE: We are not running with an ideal balance. Number of parallelized blocks'
    write(6,'(1x,2(a,i0)/)') 'nck=', peinf%nck,', not a multiple of the number of PEs=', peinf%npes
  endif

!---------------------------------

  SAFE_ALLOCATE(peinf%ik, (peinf%npes,peinf%nckpe))
  SAFE_ALLOCATE(peinf%ic, (peinf%npes,peinf%nckpe))
  SAFE_ALLOCATE(peinf%iv, (peinf%npes,peinf%nckpe))
  SAFE_ALLOCATE(peinf%ikp, (peinf%npes,peinf%nckpe))
  SAFE_ALLOCATE(peinf%icp, (peinf%npes,peinf%nckpe))
  SAFE_ALLOCATE(peinf%ivp, (peinf%npes,peinf%nckpe))
  SAFE_ALLOCATE(peinf%ipev, (peinf%npes,xct%nvb_co,xct%nkpt_co))
  SAFE_ALLOCATE(peinf%ipec, (peinf%npes,xct%ncb_co,xct%nkpt_co))
  SAFE_ALLOCATE(peinf%ipek, (peinf%npes,xct%nkpt_co))
  if (xct%qflag .ne.1) then
    SAFE_ALLOCATE(peinf%ipekq, (peinf%npes,xct%nkpt_co))
    SAFE_ALLOCATE(peinf%iownwfkq, (peinf%npes))
    peinf%ipekq=0
    peinf%iownwfkq=0
  endif
  SAFE_ALLOCATE(iown, (peinf%npes))
  SAFE_ALLOCATE(peinf%iownwfv, (peinf%npes))
  SAFE_ALLOCATE(peinf%iownwfc, (peinf%npes))
  SAFE_ALLOCATE(peinf%iownwfk, (peinf%npes))
  
  if (xct%ivpar .eq. 1) then
    SAFE_ALLOCATE(peinf%wown, (xct%nvb_co,xct%ncb_co,xct%nkpt_co,xct%nvb_co,xct%ncb_co,xct%nkpt_co))
  else if (xct%icpar .eq. 1) then
    SAFE_ALLOCATE(peinf%wown, (1,xct%ncb_co,xct%nkpt_co, 1,xct%ncb_co,xct%nkpt_co))
  else
    SAFE_ALLOCATE(peinf%wown, (1,1,xct%nkpt_co,1,1,xct%nkpt_co))
  endif

  peinf%ik=0
  peinf%ic=0
  peinf%iv=0
  peinf%ikp=0
  peinf%icp=0
  peinf%ivp=0
  peinf%ipev=0
  peinf%ipec=0
  peinf%ipek=0
  peinf%iownwfv=0
  peinf%iownwfc=0
  peinf%iownwfk=0
  peinf%wown=0
  
  ipe=0
  iown=0
  
  if (xct%ivpar .eq. 1) then
    iholdperown = 1
  else if (xct%icpar .eq. 1) then
    !FHJ: Note: this is not supported for full kernel calculations
    iholdperown = (xct%n1b_co)**2
  else
    iholdperown = (xct%n1b_co*xct%n2b_co)**2
  endif

  if (xct%ivpar .eq. 1 .or. npools .eq. 0) then
    vmax=xct%nvb_co
  else
    vmax=1
  endif

  if (xct%icpar .eq. 1 .or. npools .eq. 0) then
    cmax=xct%ncb_co
  else
    cmax=1
  endif 

  peinf%mypown = 0 

  !FHJ: This is the distribution of the WFNs. It is exactly the same for 
  ! full_kernel calculations *if* we parallelize over k^2.

  do ikp=1,xct%nkpt_co
    do ik=1,xct%nkpt_co
      do jcp=1,cmax
        do jc=1,cmax
          do kvp=1,vmax
            do kv=1,vmax
              
              if (npools .eq. 0) then
                ipe = mod((((ik-1)*xct%nkpt_co+(ikp-1))*xct%ncb_co**2 &
                  +(jc-1)*xct%ncb_co+(jcp-1))*xct%nvb_co**2 &
                  +(kv-1)*xct%nvb_co+(kvp-1),peinf%npes)
              elseif (xct%icpar .ne. 1) then
                ipoolrank = (ik-1)/(npown)
                ipool = (ikp-1)/(nown)
                ipe = ipool*(nproc/npools)+ipoolrank
              else if (xct%ivpar .ne. 1) then
                ipoolrank = (jc-1)/(npown)
                ipool = (jcp-1)/(nown)
                ipe = ipool*(nproc/npools)+ipoolrank &
                  +((ik-1)*xct%nkpt_co+(ikp-1))*nproc
              else
                ipoolrank = (kv-1)/(npown)
                ipool = (kvp-1)/(nown)
                ipe = ipool*(nproc/npools)+ipoolrank &
                  +(((ik-1)*xct%nkpt_co+(ikp-1))*xct%ncb_co**2 &
                  +(jc-1)*xct%ncb_co+(jcp-1))*nproc
              endif
              
              ipe = ipe + 1
              
              if ((xct%ivpar .eq. 1 .or. (kv .eq. 1 .and. kvp .eq. 1)) .and. &
                (xct%icpar .eq. 1 .or. (jc .eq. 1 .and. jcp .eq. 1))) then
                iown(ipe)=iown(ipe)+1
                if (npools>0 .and. ipe==peinf%inode+1) then
                  if (xct%ivpar .eq. 1) then
                    if ((mod(kvp,nown) .eq. 1 .or. nown .eq. 1)) then 
                      peinf%mypown = peinf%mypown + 1
                    endif
                  else if (xct%icpar .eq. 1) then
                    if ((mod(jcp,nown) .eq. 1 .or. nown .eq. 1)) then
                      peinf%mypown = peinf%mypown + 1
                    endif
                  else if ((mod(ikp,nown) .eq. 1 .or. nown .eq. 1)) then
                    peinf%mypown = peinf%mypown + 1
                  endif
                endif
                if (ipe==peinf%inode+1) then
                  if (xct%ivpar .eq. 1) then
                    peinf%wown(kvp,jcp,ikp,kv,jc,ik) = (iown(ipe)-1) * iholdperown + 1
                  else if (xct%icpar .eq. 1) then
                    peinf%wown(1,jcp,ikp,1,jc,ik) = (iown(ipe)-1) * iholdperown + 1
                  else
                    peinf%wown(1,1,ikp,1,1,ik) = (iown(ipe)-1) * iholdperown + 1
                  endif
                endif
                peinf%iv(ipe,iown(ipe)) = kv
                peinf%ivp(ipe,iown(ipe)) = kvp
                peinf%ic(ipe,iown(ipe)) = jc
                peinf%icp(ipe,iown(ipe)) = jcp
                peinf%ik(ipe,iown(ipe)) = ik
                peinf%ikp(ipe,iown(ipe)) = ikp
              endif
              
              if (peinf%ipec(ipe,jc,kg%indr(ik)).eq.0) then
                if (xct%icpar .eq. 1) then
                  peinf%iownwfc(ipe)=peinf%iownwfc(ipe)+1
                  peinf%ipec(ipe,jc,kg%indr(ik))=peinf%iownwfc(ipe)
                else if (jc .eq. 1) then
                  do ijk =1, xct%ncb_co
                    peinf%ipec(ipe,ijk,kg%indr(ik))=peinf%iownwfc(ipe)+ijk
                  enddo
                  peinf%iownwfc(ipe)=peinf%iownwfc(ipe)+xct%ncb_co
                endif
              endif
              
              if (peinf%ipec(ipe,jcp,kg%indr(ikp)).eq.0) then
                if (xct%icpar .eq. 1) then
                  peinf%iownwfc(ipe)=peinf%iownwfc(ipe)+1
                  peinf%ipec(ipe,jcp,kg%indr(ikp))=peinf%iownwfc(ipe)
                else if (jcp .eq. 1) then
                  do ijk =1, xct%ncb_co
                    peinf%ipec(ipe,ijk,kg%indr(ikp))=peinf%iownwfc(ipe)+ijk
                  enddo
                  peinf%iownwfc(ipe)=peinf%iownwfc(ipe)+xct%ncb_co
                endif
              endif
              
              if (xct%qflag.ne.0) then
                if (peinf%ipev(ipe,kv,kgq%indr(ik)).eq.0) then
                  if (xct%ivpar .eq. 1) then
                    peinf%iownwfv(ipe)=peinf%iownwfv(ipe)+1
                    peinf%ipev(ipe,kv,kgq%indr(ik))=peinf%iownwfv(ipe)
                  else if (kv .eq. 1) then
                    do ijk =1, xct%nvb_co
                      peinf%ipev(ipe,ijk,kgq%indr(ik))=peinf%iownwfv(ipe)+ijk
                    enddo
                    peinf%iownwfv(ipe)=peinf%iownwfv(ipe)+xct%nvb_co
                  endif
                endif
                
                if (peinf%ipev(ipe,kvp,kgq%indr(ikp)).eq.0) then
                  if (xct%ivpar .eq. 1) then
                    peinf%iownwfv(ipe)=peinf%iownwfv(ipe)+1
                    peinf%ipev(ipe,kvp,kgq%indr(ikp))=peinf%iownwfv(ipe)
                  else if (kvp .eq. 1) then
                    do ijk =1, xct%nvb_co
                      peinf%ipev(ipe,ijk,kgq%indr(ikp))=peinf%iownwfv(ipe)+ijk
                    enddo
                    peinf%iownwfv(ipe)=peinf%iownwfv(ipe)+xct%nvb_co
                  endif
                endif

                else if (xct%qflag.eq.0) then

                  if (peinf%ipev(ipe,kv,kgq%indr(xct%indexq(ik))).eq.0) then
                    if (xct%ivpar .eq. 1) then
                      peinf%iownwfv(ipe)=peinf%iownwfv(ipe)+1
                      peinf%ipev(ipe,kv,kgq%indr(xct%indexq(ik)))=peinf%iownwfv(ipe)
                    else if (kv .eq. 1) then
                      do ijk =1, xct%nvb_co
                        peinf%ipev(ipe,ijk,kgq%indr(xct%indexq(ik)))=peinf%iownwfv(ipe)+ijk
                      enddo
                      peinf%iownwfv(ipe)=peinf%iownwfv(ipe)+xct%nvb_co
                    endif
                  endif
                  
                  if (peinf%ipev(ipe,kvp,kgq%indr(xct%indexq(ikp))).eq.0) then
                    if (xct%ivpar .eq. 1) then
                      peinf%iownwfv(ipe)=peinf%iownwfv(ipe)+1
                      peinf%ipev(ipe,kvp,kgq%indr(xct%indexq(ikp)))=peinf%iownwfv(ipe)
                    else if (kvp .eq. 1) then
                      do ijk =1, xct%nvb_co
                        peinf%ipev(ipe,ijk,kgq%indr(xct%indexq(ikp)))=peinf%iownwfv(ipe)+ijk
                      enddo
                      peinf%iownwfv(ipe)=peinf%iownwfv(ipe)+xct%nvb_co
                    endif
                  endif
                endif
              
              if (peinf%ipek(ipe,kg%indr(ik)).eq.0) then
                peinf%iownwfk(ipe)=peinf%iownwfk(ipe)+1
                peinf%ipek(ipe,kg%indr(ik))=peinf%iownwfk(ipe)
              endif
              
              if (peinf%ipek(ipe,kg%indr(ikp)).eq.0) then
                peinf%iownwfk(ipe)=peinf%iownwfk(ipe)+1
                peinf%ipek(ipe,kg%indr(ikp))=peinf%iownwfk(ipe)
              endif
              
              if (xct%qflag.eq.2) then
                if (peinf%ipekq(ipe,kgq%indr(ik)).eq.0) then
                  peinf%iownwfkq(ipe)=peinf%iownwfkq(ipe)+1
                  peinf%ipekq(ipe,kgq%indr(ik))=peinf%iownwfkq(ipe)
                endif
              
                if (peinf%ipekq(ipe,kgq%indr(ikp)).eq.0) then
                  peinf%iownwfkq(ipe)=peinf%iownwfkq(ipe)+1
                  peinf%ipekq(ipe,kgq%indr(ikp))=peinf%iownwfkq(ipe)
                endif

              else if (xct%qflag.eq.0) then
                if (peinf%ipekq(ipe,kgq%indr(xct%indexq(ik))).eq.0) then
                  peinf%iownwfkq(ipe)=peinf%iownwfkq(ipe)+1
                  peinf%ipekq(ipe,kgq%indr(xct%indexq(ik)))=peinf%iownwfkq(ipe)
                endif
              
                if (peinf%ipekq(ipe,kgq%indr(xct%indexq(ikp))).eq.0) then
                  peinf%iownwfkq(ipe)=peinf%iownwfkq(ipe)+1
                  peinf%ipekq(ipe,kgq%indr(xct%indexq(ikp)))=peinf%iownwfkq(ipe)
                endif
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  if (npools .eq. 0) then
    !write(6,*) 'ds1',peinf%inode,peinf%mypown,peinf%npown,peinf%nckpe
    peinf%mypown = 1
    peinf%npown = 1
  else
    !write(6,*) 'ds2',peinf%inode,peinf%mypown,peinf%npown,peinf%nckpe
    peinf%npown = npown
  endif

!---------------------------------
! Calculate Max of iown and iownwf Their Memory Cost

  iownmax = 0
  iownwfmax = 0
  iownkmax = 0

  do ii=1,peinf%npes
    if (iown(ii) .gt. iownmax) iownmax = iown(ii)
    if ((peinf%iownwfv(ii) + peinf%iownwfc(ii)) .gt. iownwfmax) &
      iownwfmax = (peinf%iownwfc(ii)+peinf%iownwfv(ii))
  enddo
  
  do ii=1,peinf%npes
    if (peinf%iownwfk(ii) .gt. iownkmax) iownkmax = peinf%iownwfk(ii)
    if (xct%qflag.ne.1) then  
      if (peinf%iownwfkq(ii) .gt. iownkmax) iownkmax = peinf%iownwfkq(ii)
    endif
enddo
  
  peinf%myown = iown(peinf%inode + 1)
  
  if (iownmax .ne. peinf%nckpe) then
    write(tmpstr1,'(i10)') peinf%nckpe
    write(tmpstr2,'(i10)') iownmax
    tmpstr = 'nckpe estimate wrong, nckpe = ' // &
      TRUNC(tmpstr1) // ', iownmax = ' // TRUNC(tmpstr2)
    call die(tmpstr)
  endif

  SAFE_DEALLOCATE(iown)

!---------------------------------
! Determine the available memory

  call procmem(mem,nmpinode)
  mem = mem/1024.0d0**2
  if (peinf%inode==0) then
    write(6,'(1x,a,f0.1,a)') 'Memory available: ', mem, ' MB per PE'
  endif
  
!---------------------------------
! JRD: Report Memory

  if (peinf%inode .eq. 0) then

    call open_file(11,file='epsmat',form='unformatted',status='old',iostat=ierr)
    if (ierr.eq.0) then
      read(11)
      read(11)
      read(11)
      read(11)
      read(11)
      read(11)
      read(11)
      read(11) nq
      nq=nq+1
      call close_file(11)
    else
      nq=1
    endif

    rmem=0D0

! Storing epsilon
    if (xct%bLowComm) then
      rmem=rmem+dble(xct%neps*xct%neps*nq)
      !write(6,*) "Epsilon Low Comm", dble(xct%neps*xct%neps*nq)
    else
      rmem=rmem+dble(xct%neps*xct%neps*nq/peinf%npes)
      !write(6,*) "Epsilon Low Comm", dble(xct%neps*xct%neps*nq/peinf%npes)
    endif
! Storing intwfnv and intwfnc and wfnv wnfc wfnvp wfncp
    rmem=rmem+dble(iownwfmax+4)*dble(ng)
    !write(6,*) "intwfn etc...", dble(iownwfmax+4)*dble(ng)
! Storing bsemats
    rmem=rmem+4D0*dble(iownmax)*dble(iholdperown)
    !write(6,*) "bsemats", 4D0*dble(iownmax)*dble(iholdperown)
! Storing bsemat write temp arrays
    rmem=rmem+8D0*dble(xct%n1b_co*xct%n2b_co*xct%nkpt_co)
    ! FHJ: this is probablly an overshot for full_kernel
    nv=xct%n1b_co
    nc=xct%n2b_co

! JRD: Since mvv,mcc,tempb,tempw etc.. are not allocated
! at the same time as mvc etc.. we see which one is bigger

    rmemtemp1=0D0
    rmemtemp2=0D0

! Note this assumes that bare coulomb cutoff is same as wf cutoff
! Storing tempw(_old), tempb(_old), note mvv is never
! allocated at same time as tempw_old or tempb_old. So, it 
! is neglibible.
    if (xct%ivpar .eq. 1) then
      rmemtemp1=rmemtemp1+4D0*dble(ng)
      !write(6,*) "tempb(w) +old", 4D0*dble(ng) 
    else ! no 
      rmemtemp1=rmemtemp1+4D0*dble(ng*nv*nv)
      !write(6,*) "tempb(w) +old", 4D0*dble(ng*nv*nv) 
    endif
! Storing mcc, mccold
    if (xct%icpar .eq. 1) then
      rmemtemp1=rmemtemp1+2D0*dble(ng)
      !write(6,*) "mcc +old", 2D0*dble(ng) 
    else ! no mccold in this case
      rmemtemp1=rmemtemp1+1D0*dble(ng*nc*nc)
      !write(6,*) "mcc -old", 1D0*dble(ng*nc*nc) 
    endif

! Storing mvc, mvpcp, mvcold, mvpcpold
    if (xct%ivpar .eq. 1 .and. xct%icpar .eq. 1) then
      rmemtemp2=rmemtemp2+4D0*dble(ng)
      !write(6,*) "mvc +old", 4D0*dble(ng) 
    else 
      rmemtemp2=rmemtemp2+2D0*dble(ng*nv*nc)
      !write(6,*) "mvc -old", 2D0*dble(ng*nv*nc) 
    endif

! outtemp in gx_sum
    rmemtemp2=rmemtemp2+dble(iholdperown)    

    if (rmemtemp1 .gt. rmemtemp2) then
      rmem = rmem + rmemtemp1
    else
      rmem = rmem + rmemtemp2
    endif

! FFTBOXES
    call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)
    rmem=rmem+2D0*dble(Nfft(1)*Nfft(2)*Nfft(3))

    rmem=rmem * sizeof_scalar()

! Storing intwfnc%isort intwfnv%isort
    rmem=rmem+dble(2*iownkmax*gvec%ng)*4D0
! Array gvec%index_vec in input_kernel
    rmem=rmem+dble(gvec%nFFTgridpts)*4.0d0
    rmem = rmem/1024.0d0**2
    write(6,'(1x,a,f0.1,a)') 'Memory required for execution: ', rmem, ' MB per PE'
  endif

!---------------------------------------------------------
! (gsm) Determine the amount of memory required for Vcoul

! random numbers
  rmem=0.0D0
! (gsm) We don`t do random numbers in kernel anymore
!      if (xct%icutv/=TRUNC_BOX) then
! arrays ran, qran, and qran2
! (ran is deallocated before qran2 is allocated)
!        rmem=rmem+6.0D0*dble(nmc)*8.0D0
!      endif
! various truncation schemes
  call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)
  rmem2=0.0d0
! cell wire truncation
  if (xct%icutv==TRUNC_WIRE) then
    dkmax(1) = gvec%FFTgrid(1) * n_in_wire
    dkmax(2) = gvec%FFTgrid(2) * n_in_wire
    dkmax(3) = 1
    call setup_FFT_sizes(dkmax,dNfft,dscale)
! array fftbox_2D
    rmem2=rmem2+dble(dNfft(1))*dble(dNfft(2))*16.0d0
! array inv_indx
    rmem2=rmem2+dble(Nfft(1))*dble(Nfft(2))*dble(Nfft(3))*4.0d0
! array qran
    rmem2=rmem2+3.0D0*dble(nmc)*8.0D0
  endif
! cell box truncation (parallel version only)
  if (xct%icutv==TRUNC_BOX) then
    dkmax(1:3) = gvec%FFTgrid(1:3) * n_in_box
    call setup_FFT_sizes(dkmax,dNfft,dscale)
    if (mod(dNfft(3),peinf%npes) == 0) then
      Nplane = dNfft(3)/peinf%npes
    else
      Nplane = dNfft(3)/peinf%npes+1
    endif
    if (mod(dNfft(1)*dNfft(2),peinf%npes) == 0) then
      Nrod = (dNfft(1)*dNfft(2))/peinf%npes
    else
      Nrod = (dNfft(1)*dNfft(2))/peinf%npes+1
    endif
! array fftbox_2D
    rmem2=rmem2+dble(dNfft(1))*dble(dNfft(2))*dble(Nplane)*16.0d0
! array fftbox_1D
    rmem2=rmem2+dble(dNfft(3))*dble(Nrod)*16.0d0
! array dummy
!        rmem2=rmem2+dble(dNfft(1))*dble(dNfft(2))*16.0d0
! arrays dummy1 and dummy2
    rmem2=rmem2+dble(Nrod)*dble(peinf%npes+1)*16.0d0
! array inv_indx
    rmem2=rmem2+dble(Nfft(1))*dble(Nfft(2))*dble(Nfft(3))*4.0d0
  endif
  if (rmem2 .gt. rmem) rmem = rmem2
  rmem = rmem/1024.0d0**2
  if (peinf%inode==0) then
    write(6,'(1x,a,f0.1,a)') 'Extra memory required for vcoul: ', rmem, ' MB per PE'
  endif

!---------------------------------

  if (peinf%inode==0 .and. peinf%verb_high) then
    write(6,'()')
    write(6,*) iownmax*iholdperown, ' elements per PE'
    write(6,*) iownwfmax, ' wavefunctions stored per PE'
    write(6,*) ng, ' G-vectors per wavefunction'
    write(6,'()')
  endif

  POP_SUB(distrib_kernel)

  return
end subroutine distrib_kernel

end module distrib_kernel_m
