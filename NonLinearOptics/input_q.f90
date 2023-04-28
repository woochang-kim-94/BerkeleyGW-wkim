#include "f_defs.h"

module input_q_m

  use checkbz_m
  use fullbz_m
  use global_m
  use input_utils_m
  use misc_m
  use wfn_rho_vxc_io_m
  implicit none

  private

  public :: &
    input_q

contains

!-----------------------------------------------------------------------
subroutine input_q(crys,gvec,kg,kgq,syms,xct,indexq,flagbz)
!-----------------------------------------------------------------------
!
!     input: crys, gvec, kg,  syms, xct types
!
!     output: kgq type
!             indexq(1:xct%nkpt_fi_fi) : mapping of points in shifted grid
!             INT_VWFNQ_* files
!
!
!
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (grid), intent(in) :: kg
  type (grid), intent(out) :: kgq
  type (symmetry), intent(in) :: syms
  type (xctinfo), intent(inout) :: xct
  integer, intent(out) :: indexq(xct%nkpt_fi)
  integer, intent(in) :: flagbz

  type (crystal) :: crysq
  type (wavefunction) :: wfnv
  type (kpoints) :: kpq
  type (symmetry) :: symsq
  character :: filenamev*20
  character :: tmpfn*16
  integer :: iunit_v,iwrite
  integer :: ii,jj,kk,ik,ikq,irkq
  integer :: dest,tag,dummygvec(1,1)
  real(DP) :: delta,qq(3)
  real(DP) :: tol
  
  integer, allocatable :: isend(:)
  SCALAR, allocatable :: cg(:,:)
  
  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvecq, gvec_kpt
  logical :: irkq_match

  logical :: skip_checkbz
  
  PUSH_SUB(input_q)

!-----------------------------------------------------------------------
!     Read data in WFNq_fi

  if(peinf%inode.eq.0) call open_file(unit=26,file='WFNq_fi',form='unformatted',status='old')

  sheader = 'WFN'
  iflavor = 0
  call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, symsq, crysq, warn = .false.)

  call check_header('WFN_fi', kpq, gvec, syms, crys, 'WFNq_fi', kpq, gvecq, symsq, crysq, is_wfn = .true.)
  ! there is no kp object in this code??
  
  ! efermi is not currently used
  call find_efermi(xct%rfermi, xct%efermi, xct%efermi_input, kpq, kpq%mnband, 1, &
    "unshifted grid", should_search = .false., should_update = .false., write7 = .false.)
  
  call read_binary_gvectors(26, gvec%ng, gvec%ng, dummygvec, dont_read = .true.)
  
  kpq%nvband=minval(kpq%ifmax(:,:)-kpq%ifmin(:,:))+1
!
! Manual override of band numbers
!
  if (xct%vmax.ne.0) then
    kpq%nvband = xct%vmax-xct%vmin+1
    if (peinf%inode.eq.0) then
      write(0,*)
      write(0,*) '*** Overriding min/max occupied state counts'
      write(0,*) '*** for fine wavefunctions on shifted grid.'
      write(0,*) '*** New kpq%nvband=',kpq%nvband
      write(0,*)
    endif
  endif
  if(xct%nvb_fi.gt.kpq%nvband) then
    write(0,*) 'You requested ',xct%nvb_fi, &
      ' valence bands but WFNq_fi contains only ',kpq%nvband,'.'
    call die("not enough valence bands", only_root_writes = .true.)
  endif

!
!     Define shifted irreducible BZ, kgq%r, and define the shift vector
!     (if not done before). Make sure it is right!
!
  kgq%nr=kpq%nrk
  SAFE_ALLOCATE(kgq%r, (3,kgq%nr))
  kgq%r(1:3,1:kgq%nr)=kpq%rk(1:3,1:kpq%nrk)
  xct%qshift= sqrt( DOT_PRODUCT(xct%shift(:), &
    MATMUL(crys%bdot,xct%shift(:) )) )
!         write(6,*) kgq%r(:,1), kg%r(:,1)
!      if (xct%skpt.eq.0.and.xct%qshift.eq.0.d0) then
!         write(6,*) 'Bad Bad - qshift eq 0'
!         xct%shift(:)=kgq%r(:,1)-kg%r(:,1)
!         xct%qshift= sqrt( DOT_PRODUCT(xct%shift(:),
!     >                  MATMUL(crys%bdot,xct%shift(:) )) )
!      endif
  if(peinf%inode.eq.0) write(6,90) xct%shift(:),xct%qshift
90 format(/,2x,'Shift vector : ',3f9.5,2x,'Length =',f8.5,/)
!-----------------------------------------------------------------------
!     Define the polarization vector (used only with momentum operator)
!
  xct%lpol=sqrt(DOT_PRODUCT(xct%pol,MATMUL(crys%bdot,xct%pol)))
  if (xct%lpol.eq.0.d0) then
    write(0,*) 'WARNING: lpol = 0'
    xct%pol = xct%shift(:)
  endif
  xct%lpol=sqrt(DOT_PRODUCT(xct%pol,MATMUL(crys%bdot,xct%pol)))
  if(peinf%inode.eq.0) write(6,92) xct%pol(:),xct%lpol
  
92 format(/,2x,'Polarization vector : ',3f9.5,2x,'Length =',f8.5,/)
  
!-----------------------------------------------------------------------
!     Generate full brillouin zone from irreducible wedge, rk -> fk
!
  if (flagbz.eq.1) then
    call fullbz(crys,symsq,kgq,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  else
    call fullbz(crys,symsq,kgq,symsq%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  endif
  tmpfn='WFNq_fi'
  if (.not. skip_checkbz) then
    call checkbz(kgq%nf,kgq%f,kpq%kgrid,kpq%shift,crys%bdot, &
      tmpfn,'k',.true.,xct%freplacebz,xct%fwritebz)
  endif
!
  if (flagbz.eq.0.and.peinf%inode.eq.0) write(6,*) &
    'Using symmetries to expand the shifted grid sampling'
  if (flagbz.eq.1.and.peinf%inode.eq.0) write(6,*) &
    'No symmetries used in the shifted grid sampling'
!
!-----------------------------------------------------------------------
!     Find correspondence with fk from WFN_fi
!
!     indexq : correspondence between a k-point in the full BZ, kg%f, and
!       its shifted vector, in kgq%f
!     tol : tolerance
!
  tol = 1.d-6
  do ik=1,kg%nf
    ikq=0
    delta=0.1d0
    do while((delta.gt.tol).and.(ikq.lt.kgq%nf))
      ikq=ikq+1
      qq(:) = kg%f(:,ik)-(kgq%f(:,ikq)-xct%shift(:))
      do kk=1,3
        qq(kk) = qq(kk) - anint( qq(kk) )
      enddo
      delta=sqrt((qq(1))**2+(qq(2))**2+(qq(3))**2)
    enddo
    if(delta.gt.tol) then
      if(peinf%inode.eq.0) then
        write(0,*) ' Could not find point equivalent to ', (kg%f(ii,ik),ii=1,3)
      endif
      call die('k-point mismatch between WFN_fi and WFNq_fi.')
    else
!
!     make sure that kgq%f(:,ikq)-kg%f(:,ik) = shift vector
!     Near the zone edge, they may differ by a lattice vector
!
      do jj=1,3
        ii = nint( kgq%f(jj,ikq)-kg%f(jj,ik) )
        kgq%f(jj,ikq) = kgq%f(jj,ikq) - dble(ii)
        kgq%kg0(jj,ikq) = kgq%kg0(jj,ikq) - ii
      enddo
      qq(:) = kg%f(:,ik)-(kgq%f(:,ikq)-xct%shift(:))
      delta=sqrt((qq(1))**2+(qq(2))**2+(qq(3))**2)
      if (delta.gt.tol) then
        call die('k-point mismatch between WFN_fi and WFNq_fi. Wrong shift')
      endif
      indexq(ik)=ikq
    endif
  enddo
!-----------------------------------------------------------------------
!     Read the wavefunctions and create INT_VWFNQ_*
!
  wfnv%nband=xct%nvb_fi+xct%ncb_fi
  wfnv%nspin=kpq%nspin
  wfnv%nspinor=kpq%nspinor

  if(peinf%inode.lt.10000) then
    write(filenamev,'(a,i4.4)') 'INT_VWFNQ_', peinf%inode
  else
    call die("input: cannot use more than 10000 nodes.")
  endif
  iunit_v=128+(2*peinf%inode)+2
  call open_file(iunit_v,file=filenamev,form='unformatted',status='replace')
!
  SAFE_ALLOCATE(wfnv%isort, (gvec%ng))

  do irkq = 1, kpq%nrk

    wfnv%isort = 0

    irkq_match = .false.
    do ii=1,kg%nf
      if (irkq == kgq%indr(indexq(ii))) then
        irkq_match = .true.
        exit
      endif
    enddo
!
!    Skip this k-point if there is no k-point in kg%f that
!    corresponds to it
!
    SAFE_ALLOCATE(gvec_kpt%components, (3, kpq%ngk(irkq)))
    call read_binary_gvectors(26, kpq%ngk(irkq), kpq%ngk(irkq), gvec_kpt%components)
    
    SAFE_ALLOCATE(cg, (kpq%ngk(irkq), kpq%nspin*kpq%nspinor))
    if(irkq_match) then
      do ii = 1, kpq%ngk(irkq)
        call findvector(wfnv%isort(ii), gvec_kpt%components(:, ii), gvec)
        if(wfnv%isort(ii) == 0) call die('input_q: could not find gvec')
      enddo
!
      wfnv%ng = kpq%ngk(irkq)
      SAFE_ALLOCATE(wfnv%cg, (wfnv%ng,wfnv%nband,wfnv%nspin*wfnv%nspinor))

!       Determine which PEs will write the valence bands for this k-point
      iwrite=0
      do ii=1, peinf%ikt(peinf%inode+1)
        if(kgq%indr(indexq(peinf%ik(peinf%inode+1,ii))).eq. &
          irkq) then
          iwrite=1
          exit
        endif
      enddo
          
!       Determine to which PEs the valence bands for this k-point
!       need to be sent...
      SAFE_ALLOCATE(isend, (peinf%npes))
      isend=0
      if(peinf%inode.eq.0) then
        do jj=2,peinf%npes
          do ii=1, peinf%ikt(jj)
            if(kgq%indr(indexq(peinf%ik(jj,ii))).eq.irkq) then
              isend(jj)=1
              exit
            endif
          enddo
        enddo
      endif
    endif
!
!       Loop over the bands
!
    do ii=1,kpq%mnband
      call read_binary_data(26, kpq%ngk(irkq), kpq%ngk(irkq), kpq%nspin*kpq%nspinor, cg)
      
      if(.not. irkq_match) cycle
      if(peinf%inode.eq.0) then
!           Check normalization of this band
        do kk=1,kpq%nspin
          call checknorm('WFNq_fi',ii,irkq,kpq%ngk(irkq),kk,kpq%nspinor,cg(:,:))
        enddo
      endif

!         If ii is one of the selected valence band...
      if(ii.lt.(kpq%nvband+xct%ncb_fi+1).and.ii.gt.(kpq%nvband-xct%nvb_fi)) then
#ifdef MPI
        if(peinf%inode.eq.0) then
          do jj=2,peinf%npes
            if(isend(jj).eq.1) then
              dest=jj-1
              tag=1000+dest
              call MPI_SEND(cg,kpq%ngk(irkq)*kpq%nspin*kpq%nspinor,MPI_SCALAR, &
                dest,tag,MPI_COMM_WORLD,mpierr)
            endif
          enddo
        else
          if(iwrite.eq.1) then
            tag=1000+peinf%inode
            call MPI_RECV(cg,kpq%ngk(irkq)*kpq%nspin*kpq%nspinor,MPI_SCALAR,0,tag, &
              MPI_COMM_WORLD,mpistatus,mpierr)
          endif
        endif
#endif
        if(iwrite.eq.1) &
          wfnv%cg(1:wfnv%ng,ii-kpq%nvband+xct%nvb_fi,1:wfnv%nspin*wfnv%nspinor)=cg(1:wfnv%ng,1:wfnv%nspin*wfnv%nspinor)

      endif !ii is one of the selected valence band
      
    enddo

    SAFE_DEALLOCATE(cg)
    SAFE_DEALLOCATE_P(gvec_kpt%components)

    if(.not. irkq_match) cycle
    
    if(iwrite.eq.1) then      
      write(iunit_v) irkq,wfnv%ng,wfnv%nband,wfnv%nspin,wfnv%nspinor
      write(iunit_v) (wfnv%isort(ii),ii=1,gvec%ng), &
        (((wfnv%cg(ii,jj,kk),ii=1,wfnv%ng),jj=1,wfnv%nband),kk=1,wfnv%nspin*wfnv%nspinor)
    endif
    
    SAFE_DEALLOCATE(isend)
    SAFE_DEALLOCATE_P(wfnv%cg)
    
  enddo !end loop over k-points

  SAFE_DEALLOCATE_P(wfnv%isort)

  if (peinf%inode.eq.0) then
    call kpq%free()
  endif

!------------------------------
! Write out info about xtal

  if(peinf%inode.eq.0) then
    write(6,4004)
4004 format(/,2x,'crystal wavefunctions read from unit WFNq_fi')
    write(6,3007) kgq%nr
3007 format(/,6x,'nrk= ',i6,26x)
    write(6,'(12x,3f10.4)') ((kgq%r(ii,jj),ii=1,3),jj=1,kgq%nr)
    write(6,3070) kgq%nf,kgq%sz
3070 format(/,6x,'  fine grid     nfk= ',i6,4x,'ksz=',f10.5)
    
    call close_file(26)
  endif !end if(inode.eq.0)
  call close_file(iunit_v)

  ! only needed for comm_disk
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

  POP_SUB(input_q)

  return
end subroutine input_q

end module input_q_m
