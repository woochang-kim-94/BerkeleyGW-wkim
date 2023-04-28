!=========================================================================
!
! Program:
!
! (1) offdiag_QSGW           Originally By DVF        Last Modified 09/14/2010 (DVF)
!
! This routine reads in sigma_hp.log file, builds the Sigma matrix for each
! k-point and for each finite difference point (Ecor - dE, Ecor, Ecor +dE)
! according to Eq. (6) of Rohlfing & Louie PRB 62 4927, diagonalizes it with
! LAPACK and writes out the eigenvalues. This is a serial program, no MPI.
! If toff = -1/+1 the Hermitian matrix is constructed from the lower/upper
! triangle.
! 
! 06/26/10: modification of offdiag by D. Vigil 
! This routine builds the R matrix defined eq. (10) in Kotani, Van
! Schilfgaarde, and Faleev Phys. Rev. B 76 165106 (2007). The code is then
! carried out identically as before, with R replacing Sigma. See above. The
! two files read in are sig_hp_col.log and sig_hp_row.log, corresponding to
! the Sigma matrix being evaluated at the energies of its columns and rows, 
! respectively. See the Van Schilfgaarde PRB paper for details.
!
! For the spin-polarized case we need to postprocess vxc.dat_new file
! in order to combine spin-up and spin-down sections (similar to
! subroutine merge_vxc.dat in paratec). Spin-polarized usage is highly
! unlikely anyways. The complex version segfaults at the end during array 
! deallocation on some machines, but otherwise seems to work. This is being
! looked into (04/26/11). Appears to be due to an instability in the zheev 
! lapack routine. 
!=========================================================================

#include "f_defs.h"

program offdiag_QSGW

  use global_m
  use lapack_m
  use wfn_rho_vxc_io_m
  implicit none

  integer :: u_col,u_row,ierr_col,spin_col,spin_row,i_col,i_row,j_col,j_row
  integer :: l_col,l_row,ierr_row
  integer :: u_sig,u_eig,u_lda,ii
  
  integer :: freq_dep,bmin,bmax,loff,toff,fdf
  integer :: info,lda,ldvl,ldvr,lwork,nband,iw,nstart,nend,valmax
  integer :: valmaxcounter,jj,ll,mm

  real(DP) :: k_col(3),k_row(3),elda_col,elda_row,ecor_col,ecor_row
  real(DP) :: x_col,x_row,sx_col,sx_row,ch_col,ch_row,sig_row,sig_col
  real(DP) :: eqp0_col,eqp0_row,eqp1_col,eqp1_row,vxc_col,vxc_row,dw,z1,z2
#ifdef CPLX
  integer :: i2_col,i2_row,j2_col,j2_row,l2_col,l2_row
  real(DP) :: x2_col,sx2_col,ch2_col,vxc2_col,sig2_col
  real(DP) :: x2_row,sx2_row,ch2_row,sig2_row,vxc2_row
#endif
  character :: jobvl,jobvr,uplo
  character*4 :: numtype_col,numtype_row
  character*256 :: fl_col,fl_row,fl_eig,fl_sig,fl_lda,s_col,s_row,ss
  integer, allocatable :: isort(:)
  real(DP), allocatable :: ww(:)
#ifdef CPLX
  real(DP), allocatable :: rwork(:)
#endif
  SCALAR, allocatable :: ham(:,:),alda_col(:,:),alda_row(:,:),vl(:,:),vr(:,:), &
    work(:),V_qsgw(:,:),vxc_new(:,:),sig_matr(:,:),sig_matr_conj(:,:),vxc(:,:)

  type(kpoints) :: kp
  type(gspace) :: gvec
  type(symmetry) :: syms
  type(crystal) :: crys
  integer :: iflavor
  character(len=3) :: sheader

  u_col=21
  u_row=22
  u_eig=31
  u_sig=32
  u_lda=33
  uplo="U"
 
  write(jobvl,202)
  write(jobvr,802)
  write(fl_col,203)
  write(fl_row,803)
  write(fl_eig,806)
  write(fl_sig,807)
  write(fl_lda,808)


  freq_dep=-2
  bmin=0
  bmax=0
  loff=-3
  toff=-2
  fdf=-3
  dw=-1.0d0

  call open_file(unit=u_col,file=fl_col,status='old',form='formatted')
  ierr_col=0
  do while (ierr_col.eq.0)
    ! We read the frequency dependence, band index, etc. information
    ! from sig_hp_col.log only, since the column and row files
    ! should be the same
    read(u_col,206,iostat=ierr_col) s_col
    if (s_col(2:21).eq."frequency_dependence") read(s_col(22:),*) freq_dep
    if (s_col(2:11).eq."band_index") read(s_col(12:),*) bmin, bmax
    if (s_col(2:13).eq."sigma_matrix") read(s_col(14:),*) loff, toff
    if (s_col(2:23).eq."finite_difference_form") read(s_col(24:),*) fdf
    if (s_col(2:26).eq."finite_difference_spacing") read(s_col(27:),*) dw
  enddo
  call close_file(unit=u_col)
  if (freq_dep.lt.-1.or.freq_dep.gt.2.or. &
    bmin.lt.1.or.bmax.lt.bmin.or. &
    loff.lt.-2.or.(loff.gt.0.and.loff.lt.bmin).or.loff.gt.bmax.or. &
    toff.lt.-1.or.toff.gt.1.or. &
    fdf.lt.-2.or.fdf.gt.2.or. &
        dw.lt.TOL_Small) call die("Error reading file " // trim(fl_col))
  if (freq_dep.eq.2) call die(" Full frequency is not supported")
  
  if (fdf.eq.-1) then
    nstart = 1
    nend = 2
  elseif (fdf.eq.0) then
    nstart = 1
    nend = 3
  elseif (fdf.eq.1) then
    nstart = 2
    nend = 3
  else
    nstart = 2
    nend = 2
  endif
  
  nband=bmax-bmin+1
  lda=nband
  ldvl=1
  ldvr=nband
  lwork=2*nband

!DVF: What follows until the next comment is copied from wfnascbin.f90.
!I am reading in the wavefunction file in order to get the highest valence
!band.

  call open_file(unit=40,file='WFN_inner',form='unformatted',status='old')

  sheader = 'WFN'
  iflavor = 0
  call read_binary_header_type(40, sheader, iflavor, kp, gvec, syms, crys)  

  call close_file(40)

!Checking that all of the kpoints and spin polarizations
!have the same highest valence bands.Once this is checked, these 
!values are assigned to be the valence band minimum and 
!maximum for later band mixing

  valmaxcounter=0
  
  do ii = 1, kp%nrk - 1
    do jj = 1, kp%nspin
      if(kp%ifmax(ii, jj) == kp%ifmax(ii + 1, 1)) then
        valmaxcounter=valmaxcounter+1
      else
        call die("the maximum valence band is not the same for all k-points and/or polarizations")
      endif
    enddo
  enddo
  if(valmaxcounter == kp%nspin * (kp%nrk - 1)) then
    valmax = kp%ifmax(1,1)
  endif
  
  
  SAFE_ALLOCATE(isort, (nband))
  SAFE_ALLOCATE(ham,(lda,nband))
  SAFE_ALLOCATE(vxc,(lda,nband))
  SAFE_ALLOCATE(V_qsgw,(lda,nband))
  SAFE_ALLOCATE(alda_col, (lda,nband))
  SAFE_ALLOCATE(alda_row,(lda,nband))
  SAFE_ALLOCATE(sig_matr,(lda,nband))
  SAFE_ALLOCATE(vl, (ldvl,nband))
  SAFE_ALLOCATE(vr, (ldvr,nband))
  SAFE_ALLOCATE(work, (lwork))
  lwork=-1
  SAFE_ALLOCATE(ww, (nband))
  SAFE_ALLOCATE(sig_matr_conj,(lda,nband))
#ifdef CPLX
  SAFE_ALLOCATE(rwork, (2*nband))
  call zheev(jobvr,uplo,nband,ham,lda,ww,work,lwork,rwork,info)
  if (info.eq.0) lwork=int(dble(work(1)))
#else
  call dsyev(jobvr,uplo,nband,ham,lda,ww,work,lwork,info)
  if (info.eq.0) lwork=int(work(1))
#endif
  if (lwork.lt.1) then
    lwork=2*nband
  else
    SAFE_DEALLOCATE(work)
    SAFE_ALLOCATE(work, (lwork))
  endif
  
  V_qsgw(:,:) = ZERO
  alda_col(:,:) = ZERO
  alda_row(:,:) = ZERO
  vxc(:,:) = ZERO
  sig_matr_conj(:,:) = ZERO
  ham(:,:) = ZERO
  sig_matr(:,:) = ZERO
  
  call open_file(unit=u_col,file=fl_col,status='old',form='formatted')
  call open_file(unit=u_row,file=fl_row,status='old',form='formatted')
  call open_file(unit=u_sig,file=fl_sig,status='replace',form='formatted')
  call open_file(unit=u_eig,file=fl_eig,status='replace',form='formatted')
  call open_file(unit=u_lda,file=fl_lda,status='replace',form='formatted')
  call open_file(unit=44,file='vxc.dat_oldbasis',status='replace',form='formatted')
  call open_file(unit=45,file='Hamiltonian',status='replace',form='formatted')
  
  write(u_eig,'(i4,1x,i4)') bmin,bmax   ! Write the band indices to the start of eigenval_vec.log file - used in wfnmix_QSGW.f90
  
  ierr_col=0
  ierr_row=0
  do while (ierr_col.eq.0 .and. ierr_row.eq.0)
    read(u_col,206,iostat=ierr_col) s_col
    read(u_row,206,iostat=ierr_row) s_row
    if (s_col(8:10).eq."k =".and.s_row(8:10).eq."k =") then
      read(s_col(11:40),*,err=101) k_col(1:3)
      read(s_col(57:),*,err=101) spin_col
      read(s_row(11:40),*,err=602) k_row(1:3)
      read(s_row(57:),*,err=602) spin_row
      if (any(abs(k_row(1:3) - k_col(1:3)) .gt. TOL_Small)) &
        call die("k-point mismatch in sig_hp_row.log and sig_hp_col.log.")
      if (spin_row .ne. spin_col) & 
        call die("spin mismatch in sig_hp_row.log and sig_hp_col.log.")
      write(u_eig,205) k_col(1:3), spin_col
      read(u_col,*,err=101)
      read(u_col,*,err=101)
      read(u_row,*,err=602)
      read(u_row,*,err=602)
      write(*,*)
      do
        read(u_col,206,err=101) s_col
        read(u_row,206,err=602) s_row
        if (len(trim(s_col)) .eq. 0 .and. len(trim(s_row)) .eq. 0) exit
        read(s_col,207,err=101) i_col,elda_col,ecor_col,x_col,sx_col,ch_col,sig_col,vxc_col,eqp0_col,eqp1_col
        read(s_row,207,err=602) i_row,elda_row,ecor_row,x_row,sx_row,ch_row,sig_row,vxc_row,eqp0_row,eqp1_row
        write(u_lda,217,err=101) i_col,elda_col,eqp0_col,eqp1_col,sig_col,vxc_col
        alda_col(i_col-bmin+1,i_col-bmin+1)=elda_col
        alda_row(i_row-bmin+1,i_row-bmin+1)=elda_row
      enddo
      
      do iw = nstart, nend
        ham= 1.0d0/2.0d0*(alda_row+alda_col)
        if (iw.eq.1) write(31,501)  
        if (iw.eq.2) write(31,502)
        if (iw.eq.3) write(31,503)
        read(u_col,*,err=101)
        read(u_row,*,err=602)
        read(u_col,*,err=101)
        read(u_row,*,err=602)
        do
          read(u_col,206,err=101) s_col
          read(u_row,206,err=602) s_row
          if (len(trim(s_col)) .eq. 0 .and. len(trim(s_row)) .eq. 0) exit
          read(s_col,208) i_col,j_col,l_col,numtype_col,x_col,sx_col,ch_col,sig_col,vxc_col
          read(s_row,208) i_row,j_row,l_row,numtype_row,x_row,sx_row,ch_row,sig_row,vxc_row
          if(i_col .ne. i_row .or. j_col .ne. j_row) &
            call die("Matrix elements from sig_hp_row.log and sig_hp_col.log don't match!")
#ifdef CPLX
          read(u_col,206,err=101) s_col
          read(s_col,208) i2_col,j2_col,l2_col,numtype_col,x2_col,sx2_col,ch2_col,sig2_col,vxc2_col
          
          read(u_row,206,err=602) s_row
          read(s_row,208) i2_row,j2_row,l2_row,numtype_row,x2_row,sx2_row,ch2_row,sig2_row,vxc2_row
          
          if(i2_col .ne. i2_row .or. j2_col .ne. j2_row) &
            call die("Matrix elements from sig_hp_row.log and sig_hp_col.log don't match!")
#endif         
          vxc(i_col-bmin+1,j_col-bmin+1)= vxc(i_col-bmin+1,j_col-bmin+1) + &                  !Need this to construct
                                          SCALARIFY2(vxc_col,vxc2_col)                             !new Hamiltonian below
          V_qsgw(i_col-bmin+1,j_col-bmin+1) = V_qsgw(i_col-bmin+1,j_col-bmin+1)+ &                 !Adds up sigma evaluated at   
                                         1.0d0/2.0d0*SCALARIFY2(sig_col+sig_row,sig2_col+sig2_row) !row and column energies
        enddo


        !DVF:In the above do loop i_col,j_col,vxc_col should be equal to i_row,j_row,vcx_row. A test is done
        ! above to make sure the indices from the sig_hp_row.log and sig_hp_col.log files are the same. 
        ! After this check I go with the column values
        !Below the hamiltonian is constructed from the sum of sigma matrices evaluated at the energies of
        !the rows and columns and their hermitian conjugates. The van schilfgaarde potential is also
        !constructed

        sig_matr=V_qsgw
        sig_matr_conj=TRANSPOSE(MYCONJG(sig_matr))
        ham = ham + 1.0d0/2.0d0*(sig_matr + sig_matr_conj) - vxc
        V_qsgw = ZERO
        V_qsgw = 1.0d0/2.0d0*(sig_matr + sig_matr_conj)

!Write out the Hamiltonian in the old basis

        write(45,971) (k_col(ii),ii=1,3),nband,nband**2
        do jj=1,nband
          write(45,972) spin_col, jj, ham(jj,jj)
        enddo
        do jj=1,nband
          do ll=1,nband
            write(45,973) spin_col, jj, ll, ham(jj,ll)
          enddo
        enddo
        

        ! construct the Hermitian matrix from the lower triangle
        !
        if (toff.eq.-1) then
          do ii=1,nband
            do jj=ii+1,nband
              ham(ii,jj)=MYCONJG(ham(jj,ii))
            enddo
          enddo
        endif
        !
        ! construct the Hermitian matrix from the upper triangle
        !
        if (toff.eq.1) then
          do ii=1,nband
            do jj=1,ii-1
              ham(ii,jj)=MYCONJG(ham(jj,ii))
            enddo
          enddo
        endif
        !
        ! diagonalize with LAPACK
        !
#ifdef CPLX
        call zheev(jobvr,uplo,nband,ham,lda,ww,work,lwork,rwork,info)
#else
        call dsyev(jobvr,uplo,nband,ham,lda,ww,work,lwork,info)
#endif
        !
        ! sort and output eigenvalues and eigenvectors to file 
        !
        if (info.eq.0) then
          do ii=1,nband
            isort(ii)=ii
          enddo
          do ii=1,nband-1
            ll=0
            z1=1.0d6
            do jj=ii,nband
              z2=ww(isort(jj))
              if (z2.lt.z1) then
                ll=jj
                z1=z2
              endif
            enddo
            if (ll.gt.0) then
              jj=isort(ii)
              isort(ii)=isort(ll)
              isort(ll)=jj
            endif
          enddo
          
          do ii=1,nband
            write(u_eig,209) ii,ww(isort(ii))                 !Write eigenvalues to eigenval_vec.log
          enddo
          write(u_eig,*)
          do ii=1,nband
            do jj=1,nband
              write(u_eig,211) ham(jj,isort(ii))              !Write eigenvectors to eigenval_vec.log
            enddo
            write(u_eig,*)
          enddo
        else
          call die(" Error diagonalizing Sigma matrix")
        endif
      enddo
      
      SAFE_ALLOCATE(vxc_new,(nband,nband))
      vxc_new = 0.0

! Compute Vxc_new for use in next iteration
! vr is a matrix of right eigenvectors, with
! the eigenvectors as its columns

      do ii=1,nband
        do jj=1,nband
          do ll=1,nband
            do mm=1,nband
              vxc_new(ii,jj) = vxc_new(ii,jj) + MYCONJG(ham(ll,isort(ii)))*ham(mm,isort(jj))*V_qsgw(ll,mm)    
            enddo
          enddo
        enddo
      enddo
      
!The set of if statements below is need to give the proper precision
!for when sigma runs in the next iteration     


      do ll=1,3
        write(ss,815) k_col(ll)       
        if(ss(8:8)=="7") then
           ss(8:11)="6667"
        else
           if(ss(8:8)=="3") then
              ss(8:11)="3333"
           else
              if(ss(8:8)=="0") then
                 ss(8:11)="0000"
              endif
           endif
        endif
        read(ss,816) k_col(ll)
      enddo

!Write out the Van Schilfgaarde potential Vxc_new in the new basis

      write(u_sig,971) (k_col(ii),ii=1,3),nband,nband**2
      do jj=1,nband
        write(u_sig,972) spin_col, jj, COMPLEXIFY(vxc_new(jj,jj))
      enddo
      do jj=1,nband
        do ll=1,nband
          write(u_sig,973) spin_col, jj, ll, COMPLEXIFY(vxc_new(jj,ll))
        enddo
      enddo
      
!Write out the Van Schilfgaarde potential Vxc_new in the old basis

      write(44,971) (k_col(ii),ii=1,3),nband,nband**2
      do jj=1,nband
        write(44,972) spin_col, jj, COMPLEXIFY(V_qsgw(jj,jj))
      enddo
      do jj=1,nband
        do ll=1,nband
          write(44,973) spin_col, jj, ll, COMPLEXIFY(V_qsgw(jj,ll))
        enddo
      enddo


!Set potentials to zero so that you aren`t cumulatively adding the potentials from different k-points in each loop

      do ii=1,nband
        do jj=1,nband
          V_qsgw(ii,jj) = ZERO
          vxc(ii,jj) = ZERO
          sig_matr(ii,jj) = ZERO
          sig_matr_conj(ii,jj) = ZERO
        enddo
      enddo
      
      SAFE_DEALLOCATE(vxc_new) 

    endif
  enddo
  write(u_sig,*)
  write(44,*)
  
  call close_file(u_col)
  call close_file(u_row)
  call close_file(u_eig)
  call close_file(u_sig)
  call close_file(u_lda)
  call close_file(44)
  call close_file(45)

  call dealloc_header_type(sheader, crys, kp)
  
  SAFE_DEALLOCATE(isort)
  SAFE_DEALLOCATE(ham)
  SAFE_DEALLOCATE(vxc)
  SAFE_DEALLOCATE(sig_matr)
  SAFE_DEALLOCATE(V_qsgw)
  SAFE_DEALLOCATE(alda_col)
  SAFE_DEALLOCATE(alda_row)
  SAFE_DEALLOCATE(vl)
  SAFE_DEALLOCATE(vr)
  SAFE_DEALLOCATE(work)
  SAFE_DEALLOCATE(ww)
#ifdef CPLX
  SAFE_DEALLOCATE(rwork)
#endif
  SAFE_DEALLOCATE(sig_matr_conj)
  
  stop
  
101 call die("Error reading file " // trim(fl_col))
602 call die("Error reading file " // trim(fl_row))
  
202 format("N")
203 format("sig_hp_col.log")
205 format(/,1x,"k =",3f10.6,1x,"s =",i2)
206 format(a256)
207 format(i4,9f12.6)
208 format(3i4,3x,a4,3x,5f12.6)
209 format(1x,i4,2f12.6)
211 format(2f11.6)
217 format(i4,5f12.6)
501 format(/,1x,"Sig(Eo - dE)",/)
502 format(/,1x,"Sig(Eo)",/)
503 format(/,1x,"Sig(Eo + dE)",/)
802 format("V_qsgw")
803 format("sig_hp_row.log")
806 format("eigenval_vec.log") 
807 format("vxc.dat_new")
808 format("ldaoneshot.log")
815 format(f9.6)
816 format(f12.9)
971 format(3f13.9,2i8)
972 format(2i8,2f15.9)
973 format(3i8,2f15.9)
  
end program offdiag_QSGW
