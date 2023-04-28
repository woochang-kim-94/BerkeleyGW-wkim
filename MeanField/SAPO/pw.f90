
!===============================================================================
!
! Routines:
!
! 1. gsetbas()        Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Generate basis functions for degenerate set of G-vectors.
!
! 2. gsetirr()        Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Decompose reducible representation of degenerate set of G-vectors
! into irreducible representations of the group of k-point.
!
! 3. gsetrep()        Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Construct reducible representation of degenerate set of G-vectors.
!
! 4. pwdeg()          Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Construct degenerate sets of G-vectors.
!
! 5. pwgen()          Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Generate plane waves and decompose them into irreducible representations
! of the space group of the Bravais lattice or the crystal.
!
! 6. random_wfng()    Originally By gsm        Last Modified 9/1/2010 (gsm)
!
! Generate random wave functions as opposed to plane waves.
!
!===============================================================================

#include "f_defs.h"

module pw_m

  use global_m
  use lapack_m
  use groupk_m
  use sort_m
  use random_m

  implicit none

  private

  public :: pwdeg, pwgen, random_wfng

contains

subroutine gsetbas(mset,gset,rotk,nclass,nelem,elem, &
  which_irr,char_mat,pwrap,ferr,pwmix)

  integer, intent(in) :: mset,gset(3,48),rotk(3,3,48),nclass, &
    nelem(12),elem(8,12),which_irr(12),pwrap(12)
  complex(DPC), intent(in) :: char_mat(12,12)
  logical, intent(out) :: ferr
  SCALAR, intent(out) :: pwmix(48,48)
  
  logical :: feq1,feq2,fbas(48)
  integer :: i,j,irap,iset,jset,kset,iclass,ielem,imix,ndeg,gvec(3), &
    grot(3)
  real(DP) :: xnorm
  complex(DPC) :: zchar,zphase,zdot,ibas(48,48)
  
  PUSH_SUB(gsetbas)
  
! Generate basis functions for degenerate set of G-vectors

  ferr=.false.
  
  imix=0
  do irap=1,nclass
    if (pwrap(irap).ne.0) then
      ndeg=pwrap(irap)*int(dble(char_mat(irap,1))+0.5d0)

! Construct the coefficients of the basis functions for irap

      ibas(:,:)=(0.0d0,0.0d0)
      do iset=1,mset
        do j=1,3
          gvec(j)=gset(j,iset)
        enddo
        do iclass=1,nclass
          zchar=char_mat(irap,which_irr(iclass))
          do ielem=1,nelem(iclass)
            do i=1,3
              grot(i)=0
              do j=1,3
                grot(i)=grot(i)+rotk(i,j,elem(ielem,iclass))*gvec(j)
              enddo
            enddo
            do jset=1,mset
              if (grot(1).eq.gset(1,jset).and.grot(2).eq. &
                gset(2,jset).and.grot(3).eq.gset(3,jset)) then
                kset=jset
                exit
              endif
            enddo
            ibas(kset,iset)=ibas(kset,iset)+zchar
          enddo
        enddo
      enddo

! Normalize the coefficients

      do iset=1,mset
        xnorm=0.0d0
        do jset=1,mset
          xnorm=xnorm+dble(ibas(jset,iset))**2+IMAG(ibas(jset,iset))**2
        enddo
        if (xnorm.gt.TOL_Zero) then
          xnorm=sqrt(xnorm)
        else
          xnorm=1.0d0
          ferr=.true.
        endif
        do jset=1,mset
          ibas(jset,iset)=ibas(jset,iset)/CMPLX(xnorm,0.0d0)
        enddo
      enddo

! Remove an arbitrary phase

      do iset=1,mset
        jset=0
        xnorm=0.0d0
        do kset=1,mset
          if (abs(ibas(kset,iset)).gt.xnorm+TOL_Small) then
            jset=kset
            xnorm=abs(ibas(kset,iset))
          endif
        enddo
        if (jset.eq.0) then
          jset=1
          ferr=.true.
        endif
        zphase=ibas(jset,iset)
        xnorm=abs(zphase)
        if (xnorm.le.TOL_Small) then
          xnorm=1.0d0
          ferr=.true.
        endif
        zphase=zphase/CMPLX(xnorm,0.0d0)
        do jset=1,mset
          ibas(jset,iset)=ibas(jset,iset)/zphase
        enddo
      enddo

! Find the inequivalent coefficients

      fbas(:)=.false.
      do iset=1,mset
        feq1=.true.
        do jset=1,mset
          if (fbas(jset)) then
            feq2=.true.
            do kset=1,mset
              if (abs(ibas(kset,jset)-ibas(kset,iset)).gt.TOL_Small) &
                feq2=.false.
            enddo
            if (feq2) feq1=.false.
          endif
        enddo
        if (feq1) fbas(iset)=.true.
      enddo

! Check the orthonormality of the coefficients

      do iset=1,mset
        if (fbas(iset)) then
          do jset=iset+1,mset
            if (fbas(jset)) then
              zdot=CMPLX(0.0d0,0.0d0)
              do kset=1,mset
                zdot=zdot+CONJG(ibas(kset,iset))*ibas(kset,jset)
              enddo
              if (abs(zdot).gt.TOL_Zero) then
                ! orthogonalize
                do kset=1,mset
                  ibas(kset,jset)=ibas(kset,jset)-zdot*ibas(kset,iset)
                enddo
                ! normalize
                xnorm=0.0d0
                do kset=1,mset
                  xnorm=xnorm+abs(ibas(kset,jset))**2
                enddo
                if (xnorm.gt.TOL_Zero) then
                  xnorm=sqrt(xnorm)
                else
                  ! remove from the list if norm is zero
                  fbas(jset)=.false.
                  xnorm=1.0d0
                endif
                do kset=1,mset
                  ibas(kset,jset)=ibas(kset,jset)/CMPLX(xnorm,0.0d0)
                enddo
              endif ! abs(zdot).gt.TOL_Small
            endif
          enddo
        endif
      enddo

! Store the coefficients of the basis functions for irap

      iset=0
      do jset=1,mset
        if (fbas(jset)) then
          iset=iset+1
          imix=imix+1
          do kset=1,mset
            pwmix(kset,imix)=SCALARIFY(ibas(kset,jset))
          enddo
          if (iset.eq.ndeg) exit
        endif
      enddo
      
    endif ! pwrap(irap).ne.0
  enddo ! irap

  POP_SUB(gsetbas)

  return

end subroutine gsetbas

!-------------------------------------------------------------------------------

subroutine gsetirr(nclass,which_irr,char_mat,pwclass,ferr,pwrap)

  integer, intent(in) :: nclass,which_irr(12),pwclass(12)
  complex(DPC), intent(in) :: char_mat(12,12)
  logical, intent(out) :: ferr
  integer, intent(out) :: pwrap(12)
  
  integer :: iclass,irap,info,ipiv(12)
  real(DP) :: xdummy
  complex(DPC) :: a(12,12),b(12)
  
  PUSH_SUB(gsetirr)
  
! Decompose reducible representation of degenerate set of G-vectors
! into irreducible representations of the group of k-point

  ferr=.false.
  
  do iclass=1,nclass
    do irap=1,nclass
      a(iclass,irap)=char_mat(irap,which_irr(iclass))
    enddo
  enddo
  do iclass=1,nclass
    b(iclass)=CMPLX(dble(pwclass(iclass)),0.0d0)
  enddo
  
  call ZGESV(nclass,1,a,12,ipiv,b,12,info)
  if (info.ne.0) ferr=.true.
  
  do irap=1,nclass
    xdummy=dble(b(irap))
    if (xdummy.ge.0.0d0) xdummy=xdummy+0.5d0
    if (xdummy.lt.0.0d0) xdummy=xdummy-0.5d0
    pwrap(irap)=int(xdummy)
    if (pwrap(irap).lt.0) ferr=.true.
  enddo
  
! Check if decomposition is correct

  do iclass=1,nclass
    b(iclass)=CMPLX(0.0d0,0.0d0)
  enddo
  do iclass=1,nclass
    do irap=1,nclass
      b(iclass)=b(iclass)+char_mat(irap,which_irr(iclass))* &
        CMPLX(dble(pwrap(irap)),0.0d0)
    enddo
  enddo
  do iclass=1,nclass
    if (abs(b(iclass)-CMPLX(dble(pwclass(iclass)),0.0d0)).gt.TOL_Small) &
      ferr=.true.
  enddo
  
  POP_SUB(gsetirr)
  
  return
  
end subroutine gsetirr

!-------------------------------------------------------------------------------

subroutine gsetrep(mset,gset,nsymk,rotk,nclass,nelem,elem,ferr,pwclass)
  
  integer, intent(in) :: mset,gset(3,48),nsymk,rotk(3,3,48),nclass, &
    nelem(12),elem(8,12)
  logical, intent(out) :: ferr
  integer, intent(out) :: pwclass(12)
  
  integer :: i,j,iset,isymk,iclass,ielem,gvec(3),grot(3),irep(48)
  
  PUSH_SUB(gsetrep)

! Construct reducible representation of degenerate set of G-vectors

  ferr=.false.
  
  do isymk=1,nsymk
    irep(isymk)=0
  enddo
  do iset=1,mset
    do j=1,3
      gvec(j)=gset(j,iset)
    enddo
    do isymk=1,nsymk
      do i=1,3
        grot(i)=0
        do j=1,3
          grot(i)=grot(i)+rotk(i,j,isymk)*gvec(j)
        enddo
      enddo
      if (grot(1).eq.gvec(1).and.grot(2).eq.gvec(2).and. &
        grot(3).eq.gvec(3)) irep(isymk)=irep(isymk)+1
    enddo
  enddo
  
! Check if different elements of the class have the same characters

  do iclass=1,nclass
    do ielem=2,nelem(iclass)
      if (irep(elem(ielem,iclass)).ne.irep(elem(1,iclass))) &
        ferr=.true.
    enddo
  enddo

! Output one character for each class

  do iclass=1,nclass
    pwclass(iclass)=irep(elem(1,iclass))
  enddo
  
  POP_SUB(gsetrep)

  return
  
end subroutine gsetrep

!-------------------------------------------------------------------------------

subroutine pwdeg(pwsym,nsymk,rotk,ng,gvec,ekin,isrt,ndeg,mdeg,rdeg,ideg)
  
  logical, intent(in) :: pwsym
  integer, intent(in) :: nsymk,ng,rotk(3,3,48)
  integer, intent(in) :: gvec(:,:) !< (3,ng)
  real(DP), intent(in) :: ekin(:) !< (ng)
  integer, intent(inout) :: isrt(:) !< (ng)
  integer, intent(out) :: ndeg,mdeg
  integer, intent(out) :: rdeg(:,:) !< (2,ng)
  integer, intent(out) :: ideg(:) !< (ng)
  
  integer :: i,j,p,q,n,m,ig,imap,isymk,ndeg1,mdeg1,grot(3)
  real(DP) :: edeg
  integer, allocatable :: rdeg1(:,:),gvec1(:,:),gmap(:),isrt1(:)
  
  PUSH_SUB(pwdeg)

! Construct degenerate sets of G-vectors. On output, ndeg is the number
! of degenerate sets, mdeg is the maximum number of G-vectors in all
! sets, rdeg(1,:) and rdeg(2,:) are the indices of the first and last
! G-vectors in each degenerate set, and ideg(:) is the index of the
! degenerate set for each g-vector.

  SAFE_ALLOCATE(rdeg1, (2,ng))

! Sort G-vectors according to degeneracies in kinetic energies.

  ndeg1=1
  rdeg1(1,ndeg1)=1
  edeg=ekin(isrt(1))
  do ig=2,ng
    if (abs(ekin(isrt(ig))-edeg).gt.TOL_Small) then
      rdeg1(2,ndeg1)=ig-1
      ndeg1=ndeg1+1
      rdeg1(1,ndeg1)=ig
      edeg=ekin(isrt(ig))
    endif
  enddo
  rdeg1(2,ndeg1)=ng
  mdeg1=1
  do i=1,ndeg1
    j=rdeg1(2,i)-rdeg1(1,i)+1
    if (mdeg1.lt.j) mdeg1=j
  enddo
  
  if (pwsym) then

! Divide each degenerate set of G-vectors into smaller subsets
! using the symmetries of the k-point. The isrt index array is
! rearranged into smaller subsets.

    ndeg=0
    mdeg=1
    SAFE_ALLOCATE(gvec1, (3,mdeg1))
    SAFE_ALLOCATE(gmap, (mdeg1))
    SAFE_ALLOCATE(isrt1, (mdeg1))
    do n=1,ndeg1
      m=rdeg1(2,n)-rdeg1(1,n)+1
      do i=1,m
        ig=isrt(rdeg1(1,n)+i-1)
        do j=1,3
          gvec1(j,i)=gvec(j,ig)
        enddo
      enddo
      do p=1,m
        gmap(p)=0
      enddo
      imap=0
      do p=1,m
        if (gmap(p).eq.0) then
          imap=imap+1
          gmap(p)=imap
          do isymk=1,nsymk
            do i=1,3
              grot(i)=0
              do j=1,3
                grot(i)=grot(i)+rotk(i,j,isymk)*gvec1(j,p)
              enddo
            enddo
            do q=1,m
              if (grot(1).eq.gvec1(1,q).and.grot(2).eq.gvec1(2,q) &
                .and.grot(3).eq.gvec1(3,q)) gmap(q)=imap
            enddo
          enddo
        endif
      enddo
      do i=1,m
        isrt1(i)=isrt(rdeg1(1,n)+i-1)
      enddo
      p=0
      do i=1,imap
        ndeg=ndeg+1
        rdeg(1,ndeg)=rdeg1(1,n)+p
        do j=1,m
          if (gmap(j).eq.i) then
            p=p+1
            isrt(rdeg1(1,n)+p-1)=isrt1(j)
          endif
        enddo
        rdeg(2,ndeg)=rdeg1(1,n)+p-1
        j=rdeg(2,ndeg)-rdeg(1,ndeg)+1
        if (mdeg.lt.j) mdeg=j
      enddo
    enddo
    SAFE_DEALLOCATE(gvec1)
    SAFE_DEALLOCATE(gmap)
    SAFE_DEALLOCATE(isrt1)

  else ! pwsym

! Return the degenerate sets derived from kinetic energies.

    ndeg=ndeg1
    mdeg=mdeg1
    do i=1,ndeg
      do j=1,2
        rdeg(j,i)=rdeg1(j,i)
      enddo
    enddo
    
  endif ! pwsym
  
  SAFE_DEALLOCATE(rdeg1)

! Invert rdeg index array and store it in ideg for output.

  i=1
  do ig=1,ng
    ideg(ig)=i
    if (ig.eq.rdeg(2,i)) i=i+1
  enddo
  
  POP_SUB(pwdeg)
  
  return
  
end subroutine pwdeg

!-------------------------------------------------------------------------------

subroutine pwgen(fsym,fout,ffastsp,fshift,nk,ns,nbinp, &
  nbpw,nbmax,nsyml,b2g,npwdeg,ngk_l,ng,rotl,gvec, &
  pw_num,pw_ind,isort_d,eshift,a,b,bdot,kpt,en,wfn_d)
  
  logical, intent(in) :: fsym,fout,ffastsp,fshift
  integer, intent(in) :: nk,ns,nbinp,nbpw,nbmax,nsyml, &
    b2g,npwdeg,ngk_l,ng
  integer, intent(in) :: rotl(3,3,48)
  integer, intent(in) :: gvec(:,:) !< (3,ng)
  integer, intent(inout) :: pw_num(:,:,:) !< (nbmax,ns,nk)
  integer, intent(inout) :: pw_ind(:,:,:,:,:) !< (2,npwdeg,nbmax,ns,nk)
  integer, intent(in) :: isort_d(:,:) !< (ngk_l,nk)
  real(DP), intent(in) :: eshift
  real(DP), intent(in) :: a(3,3)
  real(DP), intent(in) :: b(3,3)
  real(DP), intent(in) :: bdot(3,3)
  real(DP), intent(in) :: kpt(:,:) !< (3,nk)
  real(DP), intent(inout) :: en(:,:,:) !< (nbmax,nk,ns)
  SCALAR, intent(inout) :: wfn_d(:,:,:,:) !< (ngk_l,nbmax,ns,nk)
  
  logical :: pwsym,fstart,fend,ferr
  integer :: i,j,ik,is,ib,ig,ipw,isym,ndeg,mdeg,nsymk,nclass, &
    iset,mset,jset,inum,gset(3,48),rotk(3,3,48), &
    which_irr(12),nelem(12),elem(8,12),pwclass(12),pwrap(12), &
    pw_inode(48),pw_ig(48),pw_dummy(48)
!  integer :: code_group
  real(DP) :: kvec(3)
!  real(DP) :: rotc(3,3,48)
  complex(DPC) :: char_mat(12,12)
  SCALAR :: pwmix(48,48)
  character(len=64) :: s0
  character(len=16) :: s1,s2,s3,s4,s5
!  character(len=11) :: gname
!  character(len=15) :: name_rap(12)
!  character(len=5) :: name_class(12)
!  character(len=3) :: ir_ram(12)
  
  real(DP), allocatable :: de(:)
  real(DP), allocatable :: ekin(:)
  integer, allocatable :: isrt(:)
  integer, allocatable :: rdeg(:,:)
  integer, allocatable :: ideg(:)

  PUSH_SUB(pwgen)

! Generate plane waves and decompose them into
! irreducible representations of the space group
! of the Bravais lattice or the crystal.

  SAFE_ALLOCATE(de, (ns))
  SAFE_ALLOCATE(ekin, (ng))
  SAFE_ALLOCATE(isrt, (ng))
  SAFE_ALLOCATE(rdeg, (2,ng))
  SAFE_ALLOCATE(ideg, (ng))
  
  do ik=1,nk
    
    ! generate symmetries of the k-point
    call groupk(kpt(:,ik),nsyml,rotl,nsymk,rotk)
    ! transform symmetry matrices from crystal to cartesian axes
    do isym=1,nsymk
!      call s_axis_to_cart(rotk(:,:,isym),rotc(:,:,isym),a,b)
    enddo
    ! construct the character table
!    call find_group(nsymk,rotc,gname,code_group)
!    call set_irr_rap(code_group,nclass,char_mat,name_rap,name_class,ir_ram)
!    call divide_class(code_group,nsymk,rotc,nclass,nelem,elem,which_irr)

    ! to enable (sapo_symmetry .gt. 0) comment out die below and uncomment
    ! s_axis_to_cart, find_group, set_irr_rap, divide_class above
    if (fsym) call die("(sapo_symmetry .gt. 0) is disabled")
    
    do ig=1,ng
      do i=1,3
        kvec(i)=kpt(i,ik)+dble(gvec(i,ig))
      enddo
      ! compute kinetic energies of g-vectors
      ekin(ig)=DOT_PRODUCT(kvec,MATMUL(bdot,kvec))
    enddo
    ! sort kinetic energies in ascending order
    call sortrx(ng,ekin,isrt,gvec=gvec)
    ! construct degenerate sets of g-vectors
    ! first we employ degeneracies in kinetic energies
    ! second (if pwsym is true) we divide each degenerate set
    ! into smaller subsets using the symmetries of the k-point
    pwsym=.true.
    call pwdeg(pwsym,nsymk,rotk,ng,gvec,ekin,isrt,ndeg,mdeg,rdeg,ideg)
    
    ! check if plane waves start in the middle of degenerate set
    ipw=nbinp+1-b2g
    if (ipw.ne.rdeg(1,ideg(ipw))) then
      fstart=.true.
    else
      fstart=.false.
    endif
    ! check if plane waves end in the middle of degenerate set
    ipw=nbinp+nbpw-b2g
    if (ipw.ne.rdeg(2,ideg(ipw))) then
      fend=.true.
    else
      fend=.false.
    endif
    ! print warning if the above check is positive
    if ((fstart.or.fend).and.(peinf%inode.eq.0)) then
      write(s1,101)ik
      if (fstart.and..not.fend) write(0,961)TRUNC(s1)
      if (.not.fstart.and.fend) write(0,962)TRUNC(s1)
      if (fstart.and.fend) write(0,963)TRUNC(s1)
    endif
    
    ! optional energy shift for plane waves
    if (fshift) then
      ipw=nbinp-b2g
      do is=1,ns
        de(is)=en(nbinp,ik,is)-ekin(isrt(ipw))
      enddo
    else ! fshift
      do is=1,ns
        de(is)=eshift/RYD
      enddo
    endif ! fshift
    
    ! index of the previous degenerate set
    jset=0
    
    do ib=nbinp+1,nbinp+nbpw
      
      ! absolute index of the current plane wave
      ipw=ib-b2g
      ! index of the current degenerate set
      iset=ideg(ipw)
      ! relative index of plane wave in degenerate set
      inum=ipw-rdeg(1,iset)+1
      
      ! start a new degenerate set of G-vectors
      if (iset.ne.jset) then
        jset=iset
        pw_inode(:)=0
        pw_ig(:)=0
        i=0
        ! loop over all G-vectors in degenerate set
        do j=rdeg(1,iset),rdeg(2,iset)
          i=i+1
          ! store the current G-vector in array gset
          gset(:,i)=gvec(:,isrt(j))
          ! find the processor that holds the current
          ! G-vector and its local index on that processor
          do ig=1,ngk_l
            if (isort_d(ig,ik).eq.isrt(j)) then
              pw_inode(i)=peinf%inode
              pw_ig(i)=ig
            endif
          enddo
        enddo
        ! total number of G-vectors in the current set
        mset=i
#ifdef MPI
        pw_dummy(:)=pw_inode(:)
        call MPI_Allreduce(pw_dummy,pw_inode,mset, &
          MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
        pw_dummy(:)=pw_ig(:)
        call MPI_Allreduce(pw_dummy,pw_ig,mset, &
          MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
        
        if (fsym) then
          ! symmetrize G-vectors in the current set
          ! initialize variables for printing to stdout
          write(s1,101)ik
          write(s2,101)rdeg(1,iset)+b2g
          write(s3,101)rdeg(2,iset)+b2g
          write(s4,101)rdeg(1,iset)
          write(s5,101)rdeg(2,iset)
          if (rdeg(1,iset).eq.rdeg(2,iset)) then
            write(s0,201)TRUNC(s1),TRUNC(s2),TRUNC(s4)
          else
            write(s0,202)TRUNC(s1),TRUNC(s2),TRUNC(s3), &
              TRUNC(s4),TRUNC(s5)
          endif ! rdeg(1,iset).eq.rdeg(2,iset)
          ! construct reducible representation
          ! of degenerate set of G-vectors
          call gsetrep(mset,gset,nsymk,rotk,nclass,nelem,elem,ferr,pwclass)
          if (ferr.and.peinf%inode.eq.0) write(0,964)TRUNC(s0)
          ! decompose into irreducible representations
          call gsetirr(nclass,which_irr,char_mat,pwclass,ferr,pwrap)
          if (ferr.and.peinf%inode.eq.0) write(0,965)TRUNC(s0)
          ! generate basis functions
          call gsetbas(mset,gset,rotk,nclass,nelem,elem, &
            which_irr,char_mat,pwrap,ferr,pwmix)
          if (ferr.and.peinf%inode.eq.0) write(0,966)TRUNC(s0)
          ! print representations to stdout
          if (fout) then
            if (peinf%inode.eq.0) then
              write(6,701)trim(s0)
!              write(6,702)(name_class(i),i=1,nclass)
              write(6,703)(pwclass(which_irr(i)),i=1,nclass)
!              write(6,702)(name_rap(i),i=1,nclass)
              write(6,703)(pwrap(i),i=1,nclass)
              write(6,704)
              do i=1,mset
#ifdef CPLX
                write(6,705)(pwmix(j,i),j=1,mset)
#else
                write(6,706)(pwmix(j,i),j=1,mset)
#endif
              enddo
              write(6,*)
            endif
          endif ! fout
        else ! fsym
          ! use one plane wave for one band
          pwmix(:,:)=ZERO
          do i=1,mset
            pwmix(i,i)=ONE
          enddo
        endif ! fsym
        
      endif ! iset.ne.jset
      
      ! compute the energy of the plane wave
      do is=1,ns
        en(ib,ik,is)=ekin(isrt(ipw))+de(is)
      enddo
      
      ! generate the wavefunction of the plane wave
      do is=1,ns
        wfn_d(:,ib,is,ik)=ZERO
        do i=1,mset
          if (peinf%inode.eq.pw_inode(i)) &
            wfn_d(pw_ig(i),ib,is,ik)=pwmix(i,inum)
        enddo ! i
      enddo ! is
      
      ! store plane wave components of wavefunctions
      if (ffastsp) then
        do i=1,mset
          if (abs(pwmix(i,inum)).gt.TOL_Zero) then
            do is=1,ns
              pw_num(ib,is,ik)=pw_num(ib,is,ik)+1
              pw_ind(1,pw_num(ib,is,ik),ib,is,ik)=pw_inode(i)
              pw_ind(2,pw_num(ib,is,ik),ib,is,ik)=pw_ig(i)
            enddo ! is
          endif ! abs(pwmix(i,inum)).gt.TOL_Zero
        enddo ! i
      endif ! ffastsp
      
    enddo ! ib
  enddo ! ik
  
  SAFE_DEALLOCATE(de)
  SAFE_DEALLOCATE(ekin)
  SAFE_DEALLOCATE(isrt)
  SAFE_DEALLOCATE(rdeg)
  SAFE_DEALLOCATE(ideg)
  
  POP_SUB(pwgen)
  
  return
  
101 format(i16)
  
201 format('ik =',1x,a,1x,'ib =',1x,a,1x,'ig =',1x,a)
202 format('ik =',1x,a,1x,'ib =',1x,a,1x,'...',1x,a,1x,'ig =', &
      1x,a,1x,'...',1x,a)
  
701 format(1x,a)
702 format(4x,12a5)
703 format(12i5)
704 format(4x,'Basis-function matrix')
705 format(7x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3, &
      sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',26x,f6.3,sp,f6.3,ss, &
      '*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3, &
      sp,f6.3,ss,'*i',18x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss, &
      '*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',18x,f6.3, &
      sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss, &
      '*i',2x,f6.3,sp,f6.3,ss,'*i',18x,f6.3,sp,f6.3,ss,'*i',2x,f6.3, &
      sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss, &
      '*i',18x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3, &
      sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',18x,f6.3,sp,f6.3,ss, &
      '*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3, &
      sp,f6.3,ss,'*i',18x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss, &
      '*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',18x,f6.3, &
      sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss, &
      '*i',2x,f6.3,sp,f6.3,ss,'*i',18x,f6.3,sp,f6.3,ss,'*i',2x,f6.3, &
      sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss, &
      '*i',18x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3, &
      sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',18x,f6.3,sp,f6.3,ss, &
      '*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3,sp,f6.3,ss,'*i',2x,f6.3, &
      sp,f6.3,ss,'*i')
706 format(6x,8f8.4,24x,8f8.4,16x,8f8.4,16x,8f8.4,16x,8f8.4,16x,8f8.4)
  
961 format(1x,'WARNING: ik =',1x,a,1x,'plane waves start in',1x, &
      'the middle of degenerate set',/)
962 format(1x,'WARNING: ik =',1x,a,1x,'plane waves end in',1x, &
      'the middle of degenerate set',/)
963 format(1x,'WARNING: ik =',1x,a,1x,'plane waves start and',1x, &
      'end in the middle of degenerate set',/)
964 format(1x,'WARNING:',1x,a,1x,'reducible representation',/)
965 format(1x,'WARNING:',1x,a,1x,'irreducible representation',/)
966 format(1x,'WARNING:',1x,a,1x,'basis generation',/)
  
end subroutine pwgen

subroutine random_wfng(nk,ns,nbstart,nbend,nbmax,ngk_l,ik,is,ngk, &
  arnd,brnd,wfn_d)

  integer, intent(in) :: nk,ns,nbstart,nbend,nbmax,ngk_l,ik,is
  integer, intent(in) :: ngk(:) !< (nk)
  real(DP), intent(in) :: arnd,brnd
  SCALAR, intent(inout) :: wfn_d(:,:,:,:) !< (ngk_l,nbmax,ns,nk)

  integer :: ib,ig
  real(DP) :: x,y

  PUSH_SUB(random_wfng)

  do ib=nbstart,nbend
    do ig=1,ngk_l
      if (peinf%inode*ngk_l+ig.le.ngk(ik)) then
        call genrand_real4(x)
        x=x*2.0d0-1.0d0
#ifdef CPLX
        call genrand_real4(y)
        y=y*2.0d0-1.0d0
#endif
        wfn_d(ig,ib,is,ik)=arnd*wfn_d(ig,ib,is,ik)+brnd*SCALARIFY2(x,y)
      endif ! peinf%inode*ngk_l+ig.le.ngk(ik)
    enddo ! ig
  enddo ! ib

  POP_SUB(random_wfng)

  return
end subroutine random_wfng

end module pw_m

