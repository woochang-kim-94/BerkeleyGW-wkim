#include "f_defs.h"

module input_m

  use checkbz_m
  use distrib_m
  use fullbz_m
  use global_m
  use input_utils_m
  use misc_m
  use scissors_m
  use sort_m
  use wfn_rho_vxc_io_m
  implicit none

  private

  public :: &
    input

contains

!-----------------------------------------------------------------------
subroutine input(crys,gvec,kg,syms,eqp,xct,flag)
!-----------------------------------------------------------------------
!
!       Read data from file WFN_fi and initialize variables
!
!     input: xct types
!
!     output: crys,gvec,kg,syms,eqp types
!             peinf type (from distrib.f90)
!             INT_VWFN_* and INT_CWFN_* files (if flag%vm.ne.1)
!
!
!
  type (crystal), intent(out) :: crys
  type (gspace), intent(out) :: gvec
  type (grid), intent(out) :: kg
  type (symmetry), intent(out) :: syms
  type (eqpinfo), intent(inout) :: eqp
  type (xctinfo), intent(inout) :: xct
  type (flags), intent(in) :: flag

  type (wavefunction) :: wfnc
  type (kpoints) :: kp
  character :: filenamec*20
  character :: tmpfn*16
  integer :: iunit_c,iwrite
  integer :: ii,jj,kk,irk
  integer :: dest,tag
  integer :: irks
  real(DP) :: diffvol,vcell,kt(3),div
  real(DP) :: tol
  real(DP), allocatable :: ek_tmp(:)
  integer, allocatable :: indxk(:),k_tmp(:,:)
  
  integer, allocatable :: index(:),isend(:)
  SCALAR, allocatable :: cg(:,:)
  
  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvec_kpt

  logical :: skip_checkbz
  
  PUSH_SUB(input)

  if(peinf%inode == 0) call open_file(25,file='WFN_fi',form='unformatted',status='old')

  sheader = 'WFN'
  iflavor = 0
  call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys)
  
  call logit('input: reading gvec info')
  
  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  call read_binary_gvectors(25, gvec%ng, gvec%ng, gvec%components)
  
  call get_volume(vcell,crys%bdot)
  diffvol=abs(crys%celvol-vcell)
  if (diffvol.gt.1.0d-6) then
    call die('volume mismatch', only_root_writes = .true.)
  endif
  
  ! efermi is not currently used
  call find_efermi(xct%rfermi, xct%efermi, xct%efermi_input, kp, kp%mnband, 1, &
    "unshifted grid", should_search = .true., should_update = .true., write7 = .false.)
  !call find_efermi(xct%rfermi, xct%efermi, xct%efermi_input, kp, kp%mnband, &
  !  "unshifted grid", should_search = .false., should_update = .false., write7 = .false.)

  kp%nvband=minval(kp%ifmax(:,:)-kp%ifmin(:,:))+1
  kp%ncband=kp%mnband-maxval(kp%ifmax(:,:))

! Manual override of band numbers

  if (xct%vmax.ne.0) then
    kp%nvband = xct%vmax-xct%vmin+1
    kp%ncband = kp%mnband-xct%vmax
    if (peinf%inode.eq.0) then
      write(6,*)
      write(6,*) '*** Overriding min/max occupied state counts'
      write(6,*) '*** for fine wavefunctions.'
      write(6,*) '*** New kp%nvband=',kp%nvband, ' kp%ncband=',kp%ncband
      write(6,*)
    endif
  endif

! Manual override of regular grid
! it matters only if you want to perform minicell averages with epsrm2.f

  if (xct%rgrid(1).ne.0) kp%kgrid = xct%rgrid

  xct%nspin=kp%nspin

  if(xct%nvb_fi.gt.kp%nvband) then
    if(peinf%inode.eq.0) then
      write(0,*) 'You requested ',xct%nvb_fi, &
        ' valence bands but WFN_fi contains only ',kp%nvband,'.'
    endif
    call die("valence band mismatch", only_root_writes = .true.)
  endif
  if(xct%ncb_fi.gt.kp%ncband) then
    if(peinf%inode.eq.0) then
      write(0,*) 'You requested ',xct%ncb_fi, &
        ' conduction bands but WFN_fi contains only ',kp%ncband,'.'
    endif
    call die("valence band mismatch", only_root_writes = .true.)
  endif

!-----------------------------------------------------------------------
!     Read the k-point sampling from kpoints (if it exists) or from
!     WFN_fi. In either case, this sampling will define the irreducible
!     Brillouin zone.
!
  if (xct%read_kpoints) then
    write(0,*) 'WARNING: reading k-points from file'
    if (peinf%inode.eq.0) then
      call open_file(9,file='kpoints',form='formatted',status='old')
      read(9,*) kg%nr
      SAFE_ALLOCATE(kg%r, (3,kg%nr))
      do ii=1,kg%nr
        read(9,*) (kg%r(jj,ii),jj=1,3),div
        kg%r(:,ii) = kg%r(:,ii)/div
      enddo
      call close_file(9)
    endif
#ifdef MPI
    call MPI_BCAST(kg%nr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,mpierr)
    if(peinf%inode.ne.0) then
      SAFE_ALLOCATE(kg%r, (3,kg%nr))
    endif
    call MPI_BCAST(kg%r, 3 * kg%nr, MPI_REAL_DP, 0, MPI_COMM_WORLD,mpierr)
#endif
!
!     indxk : stores the correspondence between k-points kg%r and kp%rk
!     (it is used to select the set of wavefunctions to be stored)
!     tol : tolerance in the coordinates of k-points
!
    tol = 1.d-4
    SAFE_ALLOCATE(indxk, (kg%nr))
    indxk=0
    do jj=1,kg%nr
      do ii=1,kp%nrk
        kt(:) = kg%r(:,jj) - kp%rk(:,ii)
        if ((abs(kt(1)).lt.tol).and.(abs(kt(2)).lt.tol) &
          .and.(abs(kt(3)).lt.tol)) then
          if (indxk(jj).ne.0) write(0,*) 'WARNING: multiple ', &
            'definition of k-point',jj,indxk(jj),kg%r(:,jj)
          indxk(jj)=ii
        endif
      enddo
!
!     If some k-point listed in kg%r is not found in WFN_fi, indxk
!     will store zero. Later, the job will stop in genwf.
!
      if (indxk(jj).eq.0) write(0,*) 'WARNING: could not find ', &
        'vector ',kg%r(:,jj),' in WFN_fi'
    enddo
  else
    kg%nr=kp%nrk
    SAFE_ALLOCATE(kg%r, (3,kg%nr))
    kg%r(1:3,1:kg%nr)=kp%rk(1:3,1:kp%nrk)
    SAFE_ALLOCATE(indxk, (kg%nr))
    do ii=1,kg%nr
      indxk(ii) = ii
    enddo
  endif
!-----------------------------------------------------------------------
!     Generate full Brillouin zone from irreducible wedge, rk -> fk
!
!     If flag%bz0.eq.1, only Identity will be used as
!     symmetry operation. In this case, kg%r (irreducible BZ) and kg%f
!     (full BZ) will be identical.
!
  if (flag%bz0.eq.1) then
    call fullbz(crys,syms,kg,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  else
    call fullbz(crys,syms,kg,syms%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  endif
  tmpfn='WFN_fi'

! JRD Dumb Debugging
!  write(6,*) "1st kpoint in kp%r and kp%f"
!  write(6,*) kg%r(:,1)
!  write(6,*) kg%f(:,1)

  if (.not. skip_checkbz) then
    call checkbz(kg%nf,kg%f,kp%kgrid,kp%shift,crys%bdot, &
      tmpfn,'k',.true.,xct%freplacebz,xct%fwritebz)
  endif
!
  if (flag%bz0.eq.0.and.peinf%inode.eq.0) write(6,*) &
    'Using symmetries to expand the fine-grid sampling'
  if (flag%bz0.eq.1.and.peinf%inode.eq.0) write(6,*) &
    'No symmetries used in the fine-grid sampling'
!
  xct%nkpt_fi = kg%nf

  if (.not. xct%read_kpoints) xct%nktotal=xct%nkpt_fi

! JRD Dumb Debugging
  !write(6,*) "After checkbz in kp%r and kp%f"
  !write(6,*) kg%r(:,1)
  !write(6,*) kg%f(:,1)

!-----------------------------------------------------------------------
!     Work with energy arrays: ec,ev,eclda,evlda (all in Ryd!)
!     xct%ev            QP valence energies
!     xct%ec            QP conduction energies
!     xct%evlda         LDA valence energies (used only
!        with momentum operator)
!     xct%eclda         LDA conduction energies (idem)
!      The label of bands in xct%ev and xct%evlda is again reversed:
!      from higher to lower energies.
!
!     Start rescaling energies to eV

  kp%el(:,:,:)=RYD*kp%el(:,:,:)
  SAFE_ALLOCATE(eqp%evqp, (xct%nvb_fi+xct%ncb_fi,xct%nkpt_fi,xct%nspin))
  SAFE_ALLOCATE(eqp%ecqp, (xct%ncb_fi+xct%nvb_fi,xct%nkpt_fi,xct%nspin))
  SAFE_ALLOCATE(eqp%evlda, (xct%nvb_fi+xct%ncb_fi,xct%nkpt_fi,xct%nspin))
  SAFE_ALLOCATE(eqp%eclda, (xct%ncb_fi+xct%nvb_fi,xct%nkpt_fi,xct%nspin))

  do ii=1,kp%mnband
    do jj=1,xct%nkpt_fi
      do kk=1,xct%nspin

!       If ii is one of the selected valence band... Jack Changed
        if((ii.le.(kp%nvband+xct%ncb_fi)).and. &
          (ii.gt.kp%nvband-xct%nvb_fi)) then
          eqp%evqp(ii-kp%nvband+xct%nvb_fi,jj,kk)=scissors_function(kp%el(ii,indxk(kg%indr(jj)),kk), eqp%scis%val)
          eqp%evlda(ii-kp%nvband+xct%nvb_fi,jj,kk)= &
            kp%el(ii,indxk(kg%indr(jj)),kk)
        endif
        
!       If ii is one of the selected conduction band...
        if((ii.gt.kp%nvband-xct%nvb_fi).and. &
          (ii.le.kp%nvband+xct%ncb_fi)) then
          eqp%ecqp(ii-kp%nvband+xct%nvb_fi,jj,kk)=scissors_function(kp%el(ii,indxk(kg%indr(jj)),kk), eqp%scis%cond)
          eqp%eclda(ii-kp%nvband+xct%nvb_fi,jj,kk)= &
            kp%el(ii,indxk(kg%indr(jj)),kk)
        endif
        
      enddo
    enddo
  enddo
!
!     Change energy unit to Ryd
!
  eqp%evqp(:,:,:)=eqp%evqp(:,:,:)/RYD
  eqp%ecqp(:,:,:)=eqp%ecqp(:,:,:)/RYD
  eqp%evlda(:,:,:)=eqp%evlda(:,:,:)/RYD
  eqp%eclda(:,:,:)=eqp%eclda(:,:,:)/RYD

  if (peinf%inode.eq.0) call scissors_write(6, eqp%scis)
  
!-----------------------------------------------------------------------
!       Distribute kpoints among the PEs
!
  call logit('input:  calling distrib')
  call distrib(xct)
!
!-----------------------------------------------------------------------
!     Read info for g-vectors from WFN_fi
!
  if (flag%vm.ne.0) then
    if (peinf%inode.eq.0) write(6,*) ' Bypassing INT_CWFN_* and gvec'
    SAFE_DEALLOCATE_P(gvec%components)
    call write_transitions()
    POP_SUB(input)
    return
  endif

!-----------------------------------------------------------------------
!     Order g-vectors with respect to their kinetic energy
!
  call logit('input:  reordering gvecs')
  SAFE_ALLOCATE(index, (gvec%ng))
  SAFE_ALLOCATE(gvec%ekin, (gvec%ng))
  call kinetic_energies(gvec, crys%bdot, gvec%ekin)
  call sortrx(gvec%ng, gvec%ekin, index, gvec = gvec%components)
!
  SAFE_ALLOCATE(ek_tmp, (gvec%ng))
  ek_tmp = gvec%ekin
  SAFE_ALLOCATE(k_tmp, (3,gvec%ng))
  k_tmp = gvec%components
  do ii=1,gvec%ng
    gvec%ekin(ii) = ek_tmp(index(ii))
    gvec%components(:,ii) = k_tmp(:,index(ii))
  enddo
  SAFE_DEALLOCATE(ek_tmp)
  SAFE_DEALLOCATE(k_tmp)
  SAFE_DEALLOCATE(index)

  call gvec_index(gvec)

!-----------------------------------------------------------------------
!     Read the wavefunctions and create INT_CWFN_*
!
  if (flag%vm.ne.0) then
    if (peinf%inode.eq.0) write(6,*) ' Bypassing INT_CWFN_*'
    SAFE_DEALLOCATE_P(gvec%components)
    call write_transitions()
    POP_SUB(input)
    return
  endif
  
  wfnc%nband=xct%ncb_fi+xct%nvb_fi
  wfnc%nspin=kp%nspin
  wfnc%nspinor=kp%nspinor
  
  if(peinf%inode.lt.10000) then
    write(filenamec,'(a,i4.4)') 'INT_CWFN_', peinf%inode
  else
    call die("input: cannot use more than 10000 nodes.")
  endif
  iunit_c=128+(2*peinf%inode)+1
  call open_file(iunit_c,file=filenamec,form='unformatted',status='replace')

  SAFE_ALLOCATE(wfnc%isort, (gvec%ng))
  
  do irk=1,kp%nrk

    wfnc%isort = 0

    if (peinf%inode==0 .and. peinf%verb_debug) then
      write(6,*) 'proc 0 starting kpoint = ', irk
    endif

    irks = 0
    do ii=1,kg%nr
      if (irk == indxk(ii)) then
        irks=ii
        exit
      endif
    enddo
    
! JRD Dumb Debugging
    ! write(6,*) "irk,irks", irk,irks

    SAFE_ALLOCATE(gvec_kpt%components, (3, kp%ngk(irk)))
    call read_binary_gvectors(25, kp%ngk(irk), kp%ngk(irk), gvec_kpt%components)
    
    SAFE_ALLOCATE(cg, (kp%ngk(irk),kp%nspin*kp%nspinor))

    if(irks > 0) then
      do ii = 1, kp%ngk(irk)
        call findvector(wfnc%isort(ii), gvec_kpt%components(:, ii), gvec)
        if(wfnc%isort(ii) == 0) call die('input: could not find gvec')
      enddo

      wfnc%ng = kp%ngk(irk)
      SAFE_ALLOCATE(wfnc%cg, (wfnc%ng,wfnc%nband,wfnc%nspin*wfnc%nspinor))

! Determine which PEs will write the wavefunctions for this k-point

      iwrite=0
      do ii=1, peinf%ikt(peinf%inode+1)
        if(kg%indr(peinf%ik(peinf%inode+1,ii)).eq.irks) then
          iwrite=1
          exit
        endif
      enddo

!       Determine to which PEs the wavefunctions for this k-point
!       need to be sent...

      SAFE_ALLOCATE(isend, (peinf%npes))
      isend=0
      if(peinf%inode.eq.0) then
        do jj=2,peinf%npes
          do ii=1, peinf%ikt(jj)
            if(kg%indr(peinf%ik(jj,ii)).eq.irks) then
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
    do ii=1,kp%mnband
      call read_binary_data(25, kp%ngk(irk), kp%ngk(irk), kp%nspin*kp%nspinor, cg)
      
      if(irks == 0) cycle
      
      if(peinf%inode.eq.0) then
!           Check normalization of this band
        do kk=1,kp%nspin
          call checknorm('WFN_fi',ii,irks,kp%ngk(irk),kk,kp%nspinor,cg(:,:))
        enddo
      endif
      
!         If ii is one of the selected bands...
      if(ii.gt.(kp%nvband-xct%nvb_fi).and. &
        ii.lt.(kp%nvband+xct%ncb_fi+1)) then
#ifdef MPI
        if(peinf%inode.eq.0) then
          do jj=2,peinf%npes
            if(isend(jj).eq.1) then
              dest=jj-1
              tag=1000+dest
              call MPI_SEND(cg,kp%ngk(irk)*kp%nspin*kp%nspinor,MPI_SCALAR, &
                dest,tag,MPI_COMM_WORLD,mpierr)
            endif
          enddo
        else
          if(iwrite.eq.1) then
            tag=1000+peinf%inode
            call MPI_RECV(cg,kp%ngk(irk)*kp%nspin*kp%nspinor,MPI_SCALAR,0,tag, &
              MPI_COMM_WORLD,mpistatus,mpierr)
          endif
        endif
#endif
        
        if(iwrite.eq.1) &
          wfnc%cg(1:wfnc%ng,ii-kp%nvband+xct%nvb_fi,1:wfnc%nspin*wfnc%nspinor)= &
          cg(1:wfnc%ng,1:wfnc%nspin*wfnc%nspinor)
        
      endif !ii is one of the selected bands
!
    enddo

    SAFE_DEALLOCATE(cg)
    SAFE_DEALLOCATE_P(gvec_kpt%components)
    
    if(irks == 0) cycle
    
    if(iwrite.eq.1) then
      
      write(iunit_c) irks,wfnc%ng,wfnc%nband,wfnc%nspin,wfnc%nspinor
      write(iunit_c) (wfnc%isort(ii),ii=1,gvec%ng), &
        (((wfnc%cg(ii,jj,kk),ii=1,wfnc%ng),jj=1,wfnc%nband),kk=1,wfnc%nspin*wfnc%nspinor)
    endif
    
    SAFE_DEALLOCATE(isend)
    SAFE_DEALLOCATE_P(wfnc%cg)

  enddo !end loop over k-points


  SAFE_DEALLOCATE_P(wfnc%isort)
  SAFE_DEALLOCATE(indxk)
  call close_file(iunit_c)
  
  call write_transitions()

#ifdef MPI
  ! only for comm_disk
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif


  POP_SUB(input)
  return
  
contains
  
  subroutine write_transitions()
    
    PUSH_SUB(input.write_transitions)

    if(peinf%inode.eq.0) then
      write(6,3004)
3004  format(/,2x,'crystal wavefunctions read from unit WFN_fi')
      write(6,3007) kg%nr
3007  format(/,6x,'nrk= ',i6,26x)
      write(6,'(12x,3f10.4)') ((kg%r(ii,jj),ii=1,3),jj=1,kg%nr)
      write(6,3070) kg%nf,kg%sz
3070  format(/,6x,'  fine grid     nfk= ',i6,4x,'ksz=',f10.5)
      
      call close_file(25)
    endif !end if(inode.eq.0)
    
    call kp%free()

    POP_SUB(input.write_transitions)
  end subroutine write_transitions
  
end subroutine input

end module input_m
