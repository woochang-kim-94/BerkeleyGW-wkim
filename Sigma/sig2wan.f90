!==============================================================================================
!
! Utilities:
!
! (1) sig2wan           Originally By gsm               Last Modified 2/12/2010 (gsm)
! (1) sig2wan           Originally By gsm               Last Modified 4/29/2021 (MW)
!
! This program converts Sigma output file sigma_hp.log and Wannier90 input
! file prefix.nnkp to Wannier90 input file prefix.eig. The k-points in the
! irreducible wedge (read from sigma_hp.log) are mapped to the k-points in
! the full Brillouin zone (read from prefix.nnkp) using symmetries of the
! space group of the crystal (read from sigma_hp.log). The quasiparticle
! eigenvalues (read from sigma_hp.log) for all the k-points in the full
! Brillouin zone are written to prefix.eig. This program requires input
! file sig2wan.inp in the working directory (an example follows below).
! This utility ignores off-diagonal Sigma matrix elements.
!
! sig2wan.inp:
!
! sigma_hp.log ! Sigma output file to read k-points, eigenvalues and symmetries from
! spin         ! spin component to read from sigma_hp.log file
! eqpn         ! set to 0 or 1 to read eqp0 or eqp1 from sigma_hp.log file
! prefix.nnkp  ! Wannier90 input file to read k-points from
! prefix.eig   ! file where the output of sig2wan is written
! nbands       ! number of bands to write out
! ib_start     ! band index of the first band in prefix.eig
!
!==============================================================================================

#include "f_defs.h"

program sig2wan

  use global_m
  use misc_m
  implicit none

  character*256, parameter :: finp = "sig2wan.inp"
  character*256 :: fsig,fnnkp,feig,tmpstr1,tmpstr2
  integer :: spin,eqpn,eof,is,ikr,ikf,ib,isym,nkr,nkf,nb,nsym,nb_count,ib_min,ib_max
  integer :: i,j,g(3),sym(48,3,3), ib_start
  real(DP) :: elda,ecor,x,sx,ch,sig,vxc,eqp0,eqp1,k(3),znk
  logical :: reading_offdiags
  integer, allocatable :: f2r(:)
  real(DP), allocatable :: kr(:,:)
  real(DP), allocatable :: eig(:,:)
  real(DP), allocatable :: kf(:,:)

!-----------------------------
! read input file

  write(6,'(/,1x,"reading",1x,a,1x,"file",/)')TRUNC(finp)

  call open_file(55,file=TRUNC(finp),form='formatted',status='old')

  spin=0
  eqpn=-1
  read(55,*) fsig
  read(55,*) spin
  read(55,*) eqpn
  read(55,*) fnnkp
  read(55,*) feig
  read(55,*) nb
  read(55,*) ib_start
  
  call close_file(55)
  
  if (spin.ne.1.and.spin.ne.2) then
    call die('illegal spin value: only 1 or 2 allowed')
  endif
  if (eqpn.ne.0.and.eqpn.ne.1) then
    call die('illegal eqpn value: only 0 or 1 allowed')
  endif
  if (nb <= 0) then
    call die('illegal nbands value: must be positive')
  endif
  
  write(6,901) TRUNC(fsig), spin, eqpn, TRUNC(fnnkp), TRUNC(feig)

!-----------------------------
! read sig file

  write(6,'(1x,"reading",1x,a,1x,"file",/)')TRUNC(fsig)

  call open_file(55,file=TRUNC(fsig),form='formatted',status='old')

  nkr=0
  nsym=0
  eof=0
  do while(eof.eq.0)
    read(55,'(A)',iostat=eof)tmpstr1
    write(tmpstr2,'(A)')TRUNC(tmpstr1)
    if (len(trim(tmpstr2)).eq.0) cycle
    if (tmpstr2(1:14).eq."n = band index") exit
    if (tmpstr2(1:20).eq."frequency_dependence") cycle
    ! if (tmpstr2(1:10).eq."band_index") cycle
    if (tmpstr2(1:10).eq."band_index") then
       read(tmpstr2(11:),*) ib_min, ib_max
       cycle
    endif
    if (tmpstr2(1:12).eq."sigma_matrix") cycle
    if (tmpstr2(1:17).eq."finite_difference") cycle
    if (tmpstr2(1:10).eq."symmetries") cycle
    if (tmpstr2(1:5).eq."ntran") then
      read(tmpstr2(9:11),*)nsym
      cycle
    endif
    if (tmpstr2(1:4).eq."mtrx") cycle
    if (tmpstr2(1:3).eq."k =") then
      read(tmpstr2(51:51),*)is
      if (is.eq.spin) nkr=nkr+1
      cycle
    endif
    if (tmpstr2(1:1).eq."E") cycle
    if (tmpstr2(1:1).eq."n") cycle
  enddo

  write(*,'(1X,A,I5,A,I5,A)') "Band range in "//TRUNC(fsig)//": [", ib_min, ", ", ib_max, "]"
  write(*,'(1X,A,I5,A,I5,A)') "Band range in "//TRUNC(feig)//": [", ib_start, ", ", ib_start+nb-1, "]"  

  if ( nb > (ib_max - ib_min + 1)) then
    call die("Not enough bands in sigma_hp.log")
  endif

  if ( ib_start < ib_min ) then
    call die("ib_start < ib_min")
  endif

  if ( ib_start + nb - 1 > ib_max ) then
    call die("ib_start + nb - 1 > ib_max")
  endif  
  
  rewind(55)
  
  if (nkr.le.0) then
    call die('nkr <= 0')
  endif
  if (nb.le.0) then
    call die('nb <= 0')
  endif
  if (nsym.le.0) then
    call die('nsym <= 0')
  endif
  
  SAFE_ALLOCATE(kr, (3,nkr))
  SAFE_ALLOCATE(eig, (nb,nkr))
  kr(:,:)=0.0d0
  eig(:,:)=0.0d0
  sym(:,:,:)=0
  
  ikr=0
  eof=0
  reading_offdiags = .false.
  do while(eof.eq.0)
    read(55,'(A)',iostat=eof)tmpstr1
    write(tmpstr2,'(A)')TRUNC(tmpstr1)
    if (len(trim(tmpstr2)).eq.0) cycle
    if (tmpstr2(1:8).eq."========") exit
    if (tmpstr2(1:14).eq."n = band index") exit
    if (tmpstr2(1:20).eq."frequency_dependence") cycle
    if (tmpstr2(1:10).eq."band_index") cycle
    if (tmpstr2(1:12).eq."sigma_matrix") cycle
    if (tmpstr2(1:17).eq."finite_difference") cycle
    if (tmpstr2(1:10).eq."symmetries") cycle
    if (tmpstr2(1:5).eq."ntran") cycle
    if (tmpstr2(1:4).eq."mtrx") then
      read(tmpstr2(5:6),*)isym
      read(tmpstr2(9:),*)((sym(isym,i,j),i=1,3),j=1,3)
      cycle
    endif
    if (tmpstr2(1:3).eq."k =") then
      read(tmpstr2(51:51),*)is
      if (is.eq.spin) then
        ikr=ikr+1
        read(tmpstr2(4:33),*)(kr(i,ikr),i=1,3)
      endif
      reading_offdiags = .false.
      nb_count=0      
      cycle      
    endif
    if (tmpstr2(1:1).eq."E") cycle
    if (tmpstr2(1:16).eq."n            Emf") cycle
    if (tmpstr2(1:16).eq."n           elda") cycle
    if (tmpstr2(1:9).eq."n   m   l") then
      reading_offdiags = .true.
      cycle      
    endif
    if (is.eq.spin .and. .not. reading_offdiags) then
      read(tmpstr2,*) ib,elda,ecor,x,sx,ch,sig,vxc,eqp0,eqp1,znk
      if ( (ib >= ib_start) .and. nb_count <= nb ) then
        nb_count = nb_count + 1
        if (eqpn.eq.0) eig(nb_count,ikr)=eqp0
        if (eqpn.eq.1) eig(nb_count,ikr)=eqp1
      endif
    endif
  enddo
  
  call close_file(55)
  
  write(tmpstr1,*)nkr
  write(6,'(11x,"nkr =",1x,a)') TRUNC(tmpstr1)
  write(tmpstr1,*)nb
  write(6,'(12x,"nb =",1x,a)') TRUNC(tmpstr1)
  write(tmpstr1,*)nsym
  write(6,'(10x,"nsym =",1x,a,/)') TRUNC(tmpstr1)

!-----------------------------
! read nnkp file

  write(6,'(1x,"reading",1x,a,1x,"file",/)')TRUNC(fnnkp)
  
  call open_file(55,file=TRUNC(fnnkp),form='formatted',status='old')
  
  i=0
  ikf=0
  nkf=0
  eof=0
  do while(eof.eq.0)
    read(55,'(A)',iostat=eof)tmpstr1
    write(tmpstr2,'(A)')TRUNC(tmpstr1)
    if (len(trim(tmpstr2)).eq.0) cycle
    if (i.eq.2) then
      ikf=ikf+1
      read(tmpstr2,*)(kf(j,ikf),j=1,3)
      if (ikf.eq.nkf) i=0
      cycle
    endif
    if (i.eq.1) then
      read(tmpstr2,*)nkf
      SAFE_ALLOCATE(kf, (3,nkf))
      kf(:,:)=0.0d0
      i=2
      cycle
    endif
    if (tmpstr2(1:13).eq."begin kpoints") i=1
  enddo
  
  if (nkf.le.0) then
    call die('nkf <= 0')
  endif
  
  call close_file(55)
  
  write(tmpstr1,*)nkf
  write(6,'(11x,"nkf =",1x,a,/)') TRUNC(tmpstr1)

!-----------------------------
! unfold irreducible wedge

  write(6,'(1x,"unfolding irreducible wedge",/)')

  SAFE_ALLOCATE(f2r, (nkf))
  f2r(1:nkf) = 0
  
  do ikf = 1, nkf
    ikr_do: do ikr = 1, nkr
      do isym = 1, nsym
        k(1:3) = matmul(sym(isym, 1:3, 1:3), kr(1:3, ikr))
        !First check without folding and then fold
        if (all(abs(k(1:3) - kf(1:3, ikf)) .lt. TOL_Small)) then
          f2r(ikf) = ikr
          exit ikr_do
        endif
        call k_range(k, g, TOL_Small)
        if (all(abs(k(1:3) - kf(1:3, ikf)) .lt. TOL_Small)) then
          f2r(ikf) = ikr
          exit ikr_do
        endif
      enddo
    enddo ikr_do
    
    if (f2r(ikf) .eq. 0) then
      write(tmpstr1,'(a,i6)') "no match for ikf = ", ikf
      call die(tmpstr1)
    endif
  enddo
      
!-----------------------------
! write eig file

  write(6,'(1x,"writing",1x,a,1x,"file",/)')TRUNC(feig)
  
  call open_file(55,file=TRUNC(feig),form='formatted',status='replace')
  do ikf=1,nkf
    ikr=f2r(ikf)
    do ib=1,nb
      write(55,'(2i5,f18.12)') ib,ikf,eig(ib,ikr)
    enddo
  enddo
  
  call close_file(55)

!-----------------------------
! deallocate and finish

  SAFE_DEALLOCATE(kr)
  SAFE_DEALLOCATE(eig)
  SAFE_DEALLOCATE(kf)
  SAFE_DEALLOCATE(f2r)
  
901 FORMAT(6x,"sig file =",1x,a,/,10x,"spin =", &
      1x,i1,/,10x,"eqpn =",1x,i1,/,5x,"nnkp file =", &
      1x,a,/,6x,"eig file =",1x,a,/)
  
end program sig2wan
