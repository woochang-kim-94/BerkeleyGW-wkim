!================================================================================
!
! Module eqpcor_m
!
! 1. eqpcor()        Originally By gsm       Last Modified 8/4/2015 (FHJ)
!
!    Reads quasiparticle energy corrections from eqp.dat-type files.
!    Such files are made by bin/eqp.py.
!
!    There are some differences in how eqpcor is called from different codes.
!
!    In Epsilon, eqpcor is called from proc 0, so set inode = 0 and npes = 1.
!    In other places, set inode = peinf%inode and npes = peinf%npes.
!
!    Epsilon and BSE require quasiparticle energies in Rydbergs, so set
!    irydflag = 1. Sigma (outer) requires energies in eVs, so set irydflag = 0.
!
!    In Epsilon and Sigma, the quasiparticle energies are returned in
!    array eqp(nbmin:nbmax,1:kp%nrk,1:kp%nspin). In BSE, the valence and
!    conduction energies indexed with respect to the Fermi level are
!    returned in arrays eqpv(:,1:kp%nrk,1:kp%nspin) and eqpc(:,1:kp%nrk,1:kp%nspin).
!
!    In inteqp, set ivalflag = 2 to return the difference Eqp - Edft.
!    In other places, set ivalflag = 0 to return Eqp.
!
!    DO NOT try to use both styles, they will overwrite each other!
!
!================================================================================

#include "f_defs.h"

module eqpcor_m

  use global_m
  implicit none

  private

  public :: eqpcor

contains

subroutine eqpcor(fn,inode,npes,kp,nbmin,nbmax, &
  nbvnum,nbcnum,eqp,eqpv,eqpc,irydflag,ivalflag,dont_write,kpt_required,allow_trs)

  character(len=32), intent(in) :: fn
  integer, intent(in) :: inode,npes
  type(kpoints), intent(in) :: kp !< however, usually kp%el will be passed to eqp, eqpv, eqpc
  integer, intent(in) :: nbmin,nbmax,nbvnum,nbcnum
  real(DP), intent(inout) :: eqp(:,:,:)      !< (nb,  kp%nrk, kp%nspin)
  real(DP), intent(inout) :: eqpv(:,:,:)     !< (nbv, kp%nrk, kp%nspin)
  real(DP), intent(inout) :: eqpc(:,:,:)     !< (nbc, kp%nrk, kp%nspin)
  integer, intent(in) :: irydflag,ivalflag
  logical, optional, intent(in) :: dont_write !< silence output
  !> (kp%nrk) Which k-points are required in the eqp file? Defaults to all
  logical, optional, intent(in) :: kpt_required(:)
  !> allow use of time-reversal symmetry to find k-points? default is .true.
  logical, optional, intent(in) :: allow_trs

  integer :: eof, jk, kk, is, ib, isb, nsb, nk_found, &
    nbvmin, nbvmax, nbcmin, nbcmax, corsmin, corsmax, corbmin, corbmax
  integer :: itrs, itrs_max
  real(DP) :: dummydft, dummyqp, k(3), dk(3)
  character(100) :: errmsg
  logical :: dowrite, allow_trs_
  logical, allocatable :: kpt_present(:)

  ! eqp.py from sigma.log may not give better precision than this.
  real(DP), parameter :: TOL_eqp = 2d-6
  
  PUSH_SUB(eqpcor)

  dowrite = .true.
  if (present(dont_write)) dowrite = .not. dont_write
  dowrite = dowrite .and. peinf%verb_medium
  allow_trs_ = .true.
  if (present(allow_trs)) allow_trs_ = allow_trs
  itrs_max = 0
  if (allow_trs_) itrs_max = 1
  if (inode==0) then
    SAFE_ALLOCATE(kpt_present, (kp%nrk))
    kpt_present(:) = .false.
    write(6,'(/1x,a)') "Reading quasiparticle energy corrections from "//trim(fn)
    call open_file(9,file=trim(fn),form='formatted',status='old')

    if (ivalflag.eq.2) then
      if(dowrite) write(6,922)
922 format(/,6x,"band",3x,"k-point",6x,"spin",3x,"DeltaElda (eV)",2x,"DeltaEqp (eV)")
    else
      if(dowrite) write(6,902)
902 format(/,6x,"band",3x,"k-point",6x,"spin",6x,"Elda (eV)",7x,"Eqp (eV)")
    endif
    corsmin=kp%nspin+1
    corsmax=0
    corbmin=kp%mnband+1
    corbmax=0
    nk_found=0
    do while (nk_found<kp%nrk) !This is the k-point loop
      read(9,*,iostat=eof) k(:), nsb
      if (eof/=0) then
        if (present(kpt_required)) exit
        write(errmsg,'(a)') 'Missing k-points in file '//trim(fn)
        call die(errmsg, only_root_writes=.true.)
      endif
      kk=0
      trs_loop: do itrs = 0, itrs_max
        ! FHJ: This is either k or -k
        k = (1 - 2*itrs) * k
        do jk=1,kp%nrk
          dk(:) = k(:) - kp%rk(:,jk)
          ! this allows Umklapp, i.e. either 0.1 or -0.9 to be supplied
          dk = MIN_RANGE(dk)
          if (all(dabs(dk)<TOL_SMALL)) then
            kk = jk
            kpt_present(kk) = .true.
            nk_found = nk_found + 1 
            exit trs_loop
          endif
        enddo
      enddo trs_loop
      do isb=1,nsb
        read(9,*,iostat=eof)is,ib,dummydft,dummyqp
        if(kk == 0) cycle
        if (ivalflag.eq.2) then
          dummyqp=dummyqp-dummydft
        endif
        nbvmin=kp%ifmax(kk,is)-nbvnum+1
        nbvmax=kp%ifmax(kk,is)
        nbcmin=kp%ifmax(kk,is)+1
        nbcmax=kp%ifmax(kk,is)+nbcnum
        if (eof/=0) then
          write(errmsg,'(a)') 'Wrong contents of k-point blocks in file '//trim(fn)
          call die(errmsg, only_root_writes=.true.)
        endif
        if (is.lt.corsmin) corsmin=is
        if (is.gt.corsmax) corsmax=is
        if (ib.lt.corbmin) corbmin=ib
        if (ib.gt.corbmax) corbmax=ib
        if (ib.ge.nbmin.and.ib.le.nbmax) then
          if (irydflag.eq.0) then
            if(dowrite) write(6,903)ib,kk,is,eqp(ib,kk,is),dummyqp
            ! with ivalflag == 2, we do not have the DFT energy in eqp so cannot compare.
            if(ivalflag /= 2 .and. abs(dummydft - eqp(ib,kk,is)) > TOL_eqp) &
              call die("eqpcor mean-field energy mismatch")
            eqp(ib,kk,is)=dummyqp
          else
            if(dowrite) write(6,903)ib,kk,is,eqp(ib,kk,is)*RYD,dummyqp
            if(ivalflag /= 2 .and. abs(dummydft/RYD - eqp(ib,kk,is)) > TOL_eqp) &
              call die("eqpcor mean-field energy mismatch")
            eqp(ib,kk,is)=dummyqp/RYD
          endif
903           format(3i10,2f15.5)
        endif
        if (ib.ge.nbvmin.and.ib.le.nbvmax) then
          if (irydflag.eq.0) then
            if(dowrite) write(6,904)ib,kk,is,eqpv(nbvmax-ib+1,kk,is),dummyqp
            if(ivalflag /= 2 .and. abs(dummydft - eqpv(nbvmax-ib+1,kk,is)) > TOL_eqp) &
              call die("eqpcor mean-field energy mismatch")
            eqpv(nbvmax-ib+1,kk,is)=dummyqp
          else
            if(ivalflag /= 2 .and. abs(dummydft/RYD - eqpv(nbvmax-ib+1,kk,is)) > TOL_eqp) &
              call die("eqpcor mean-field energy mismatch")
            if(dowrite) write(6,904)ib,kk,is,eqpv(nbvmax-ib+1,kk,is)*RYD,dummyqp
            eqpv(nbvmax-ib+1,kk,is)=dummyqp/RYD
          endif
904           format("v",i9,2i10,2f15.5)
        endif
        if (ib.ge.nbcmin.and.ib.le.nbcmax) then
          if (irydflag.eq.0) then
            if(dowrite) write(6,905)ib,kk,is,eqpc(ib-nbcmin+1,kk,is),dummyqp
            if(ivalflag /= 2 .and. abs(dummydft - eqpc(ib-nbcmin+1,kk,is)) > TOL_eqp) &
              call die("eqpcor mean-field energy mismatch")
            eqpc(ib-nbcmin+1,kk,is)=dummyqp
          else
            if(dowrite) write(6,905)ib,kk,is,eqpc(ib-nbcmin+1,kk,is)*RYD,dummyqp
            if(ivalflag /= 2 .and. abs(dummydft/RYD - eqpc(ib-nbcmin+1,kk,is)) > TOL_eqp) &
              call die("eqpcor mean-field energy mismatch")
            eqpc(ib-nbcmin+1,kk,is)=dummyqp/RYD
          endif
905           format("c",i9,2i10,2f15.5)
        endif
      enddo
    enddo
    call close_file(9)
    if(any(kpt_present)) then
      nbvmin=minval(kp%ifmax(:,:))-nbvnum+1
      nbvmax=maxval(kp%ifmax(:,:))
      nbcmin=minval(kp%ifmax(:,:))+1
      nbcmax=maxval(kp%ifmax(:,:))+nbcnum
      if (corsmin>1.or.corsmax<kp%nspin) then
        write(errmsg,'(a)') 'Missing spins in file '//trim(fn)
        call die(errmsg, only_root_writes=.true.)
      endif
      if ((nbmin/=0.and.nbmax/=0).and.(corbmin>nbmin.or.corbmax<nbmax)) then
        write(errmsg,'(a)') 'Missing bands in file '//trim(fn)
        call die(errmsg, only_root_writes=.true.)
      endif
      if ((nbvnum/=0).and.(corbmin>nbvmin.or.corbmax<nbvmax)) then
        write(errmsg,'(a)') 'Missing valence bands in file '//trim(fn)
        call die(errmsg, only_root_writes=.true.)
      endif
      if ((nbcnum/=0).and.(corbmin>nbcmin.or.corbmax<nbcmax)) then
        write(errmsg,'(a)') 'Missing conduction bands in file '//trim(fn)
        call die(errmsg, only_root_writes=.true.)
      endif
    endif
    if (dowrite) write(6,*)
    ! FHJ: We might think there are missing kpts, so make sure all k-points
    ! we care about were present in the file.
    if (present(kpt_required)) then
      if (.not.all(kpt_present(:).or.(.not.kpt_required(:)))) then
        write(errmsg,'(a)') 'Missing k-points in file '//trim(fn)
        call die(errmsg, only_root_writes=.true.)
      endif
    endif
    SAFE_DEALLOCATE(kpt_present)
  endif

#ifdef MPI
  if (npes>1) then
    if (nbmin/=0.and.nbmax/=0) then
      call MPI_Bcast(eqp(1,1,1),size(eqp),MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
    endif
    if (nbvnum/=0) then
      call MPI_Bcast(eqpv(1,1,1),size(eqpv),MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
    endif
    if (nbcnum/=0) then
      call MPI_Bcast(eqpc(1,1,1),size(eqpc),MPI_REAL_DP,0,MPI_COMM_WORLD,mpierr)
    endif
  endif
#endif
  
  POP_SUB(eqpcor) 
  return
  
end subroutine eqpcor

end module eqpcor_m
