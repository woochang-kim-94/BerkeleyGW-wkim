!=============================================================================
!
! Utilities:
!
! (1) wfnmix_spinor        Originally By BAB       Last Modified 5/14/2012 (BAB)
!
!     Based on DVF`s wfnmix_QSGW, this code takes an ordinary wavefunction
!     and mixes it into a spinor that would diagonalize either sigma_z 
!     or sigma_x, or is mixed randomly (for each degenerate subspace)
!
!     The purpose is to debug the spinor functionality of BGW.
!
!==============================================================================

#include "f_defs.h"

program wfnmix_spinor

  use blas_m
  use global_m
  use misc_m
  use wfn_rho_vxc_io_m
  use susymmetries_m
  use random_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp,kps
  type(gspace) :: gvec

  character(len=80)  :: infile, outfile, spin_vec
  integer :: nargs,ik,ib,jb,ispin,ispinor,spin_type
  
  complex(DPC) :: ovlp
  real(DP) :: ediff_back, ediff_forw
  
  SCALAR, allocatable :: wfn(:,:,:), wfn_new(:,:,:), wfn_temp(:,:)
  complex(DPC) :: umtrx(2,2)

  character(len=3) :: sheader
  integer :: iflavor
  
#ifndef CPLX
  call die('Can not use wfnmix_spinor in real case')
#endif

!---------------------------------
! Get file names from command-line arguments

  nargs = command_argument_count()
  
  if (nargs .ne. 3) then
    call die('Usage: wfnmix_spinor WFN WFN2 x_or_z_or_r')
  endif
  
  call get_command_argument(1,infile)
  call get_command_argument(2,outfile)
  call get_command_argument(3,spin_vec)
  if (spin_vec == 'x' .or. spin_vec == 'X') then
    spin_type = 0
  else if (spin_vec == 'z' .or. spin_vec == 'Z') then
    spin_type = 1
  else if (spin_vec == 'r' .or. spin_vec == 'R') then
    spin_type = 2
  else
    call die('Type "x" or "X" or "z" or "Z" or "r" or "R" only')
  endif

! Open units

  call open_file(unit=7,file=TRUNC(infile),form='unformatted',status='old')
  call open_file(unit=8,file=TRUNC(outfile),form='unformatted',status='replace')
  write(6,*) 'Converting file ',TRUNC(infile),' into spinor wavefunction file'

! Read/write data

  sheader = 'WFN'
  iflavor = 0
  call read_binary_header_type(7, sheader, iflavor, kp, gvec, syms, crys, warn = .false., dont_warn_kgrid = .true.)
  if (iflavor.eq.1) then
    call die('WFN must be type CPLX')
  endif

! Allocating arrays for kps type

  if (kp%nspinor.eq.2) then 
    call die('Input WFN file already has two spinor components')
  endif

  kps%nspinor=2
  kps%nspin=kp%nspin
  if (kps%nspin.ne.1) then
    call die('nspin not equal to unity')
  endif
  kps%mnband=2*kp%mnband
  kps%nvband=2*kp%nvband
  kps%ncband=2*kp%ncband
  kps%nrk=kp%nrk
  kps%kgrid(:)=kp%kgrid(:)
  kps%shift(:)=kp%shift(:)
  kps%ecutwfc=kp%ecutwfc

  SAFE_ALLOCATE( kps%ngk, (kp%nrk))
  kps%ngk(:)=kp%ngk(:)
  kps%ngkmax=kp%ngkmax        

  SAFE_ALLOCATE( kps%ifmin, (kp%nrk, kps%nspinor*kps%nspin))
  SAFE_ALLOCATE( kps%ifmax, (kp%nrk, kps%nspinor*kps%nspin))

  do ispinor=1,kps%nspinor
    kps%ifmin(:,ispinor) = 2*kp%ifmin(:,kp%nspin)-1
    kps%ifmax(:,ispinor) = 2*kp%ifmax(:,kp%nspin)
  enddo

  SAFE_ALLOCATE( kps%w, (kp%nrk))
  kps%w(:)=kp%w(:)

  SAFE_ALLOCATE( kps%rk, (UBOUND(kp%rk,1)-LBOUND(kp%rk,1)+1,kp%nrk))
  kps%rk(:,:)=kp%rk(:,:)

  SAFE_ALLOCATE( kps%el, (kps%mnband, UBOUND(kp%el,2)-LBOUND(kp%el,2)+1, UBOUND(kp%el,3)-LBOUND(kp%el,3)+1))
  do ib=1,kp%mnband
    kps%el(2*ib-1,:,:)=kp%el(ib,:,:)
    kps%el(2*ib,:,:)=kp%el(ib,:,:)
  enddo

  SAFE_ALLOCATE(kps%occ, (kps%mnband, kp%nrk, kps%nspin*kps%nspinor))
  do ispinor=1, kps%nspinor
    do ib=1,kp%mnband
      kps%occ(2*ib-1,:,ispinor)=kp%occ(ib,:,kp%nspinor)
      kps%occ(2*ib,:,ispinor)=kp%occ(ib,:,kp%nspinor)
    enddo
  enddo

  ! Now write new header info
  call write_binary_header_type(8, sheader, iflavor, kps, gvec, syms, crys)

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))

  call read_binary_gvectors(7, gvec%ng, gvec%ng, gvec%components)
  call write_binary_gvectors(8, gvec%ng, gvec%ng, gvec%components)

! Output info

  write(6,*) ' crystal volume: ',crys%celvol
  write(6,*) ' number of spins: ',kps%nspin
  write(6,*) ' number of spinors: ',kps%nspinor
  write(6,*) ' number of bands in file: ',kps%mnband
  
  SAFE_ALLOCATE(wfn, (kp%ngkmax, kp%mnband, kp%nspin))
  SAFE_ALLOCATE(wfn_new, (kps%ngkmax, kps%mnband, kps%nspin*kps%nspinor))

  do ik = 1, kp%nrk
    call read_binary_gvectors(7, kp%ngk(ik), kp%ngk(ik), gvec%components)
    call write_binary_gvectors(8, kps%ngk(ik), kps%ngk(ik), gvec%components)

    do ib = 1, kp%mnband
      call read_binary_data(7, kp%ngk(ik), kp%ngkmax, kp%nspin, wfn(:,ib,:))
      do ispin = 1, kp%nspin
        call checknorm('old wfn',ib,ik,kp%ngk(ik),ispin,kp%nspinor,wfn(1:kp%ngk(ik),ib,:))
      enddo
    enddo

    !
    ! We construct two-component spinor out of single-component wavefunctions here
    !
    do ib = 1, kps%mnband
      if (spin_type.eq.0) then
        if (mod(ib,2) == 1) then
          wfn_new(:,ib,1) = 1.0D0/SQRT(2.0d0)*wfn(:,ib/2+1,1) 
          wfn_new(:,ib,2) = 1.0D0/SQRT(2.0d0)*wfn(:,ib/2+1,1)
        else 
          wfn_new(:,ib,1) = 1.0D0/SQRT(2.0d0)*wfn(:,ib/2,1) 
          wfn_new(:,ib,2) = -1.0D0/SQRT(2.0d0)*wfn(:,ib/2,1)
        endif 
        call checknorm('new wfn',ib,ik,kps%ngk(ik),1,2,wfn_new(1:kps%ngk(ik),ib,:))
      else
        if (mod(ib,2) == 1) then
          wfn_new(:,ib,1) = wfn(:,ib/2+1,1)
          wfn_new(:,ib,2) = ZERO
        else 
          wfn_new(:,ib,1) = ZERO 
          wfn_new(:,ib,2) = wfn(:,ib/2,1) 
        endif
        call checknorm('new wfn',ib,ik,kps%ngk(ik),1,2,wfn_new(1:kps%ngk(ik),ib,:))
      endif
    enddo

    !! Random rotation part starts here

    if (spin_type.eq.2) then

      do ib = 1, kps%mnband

        ! this is the difference between previous band and this band, to see if we are in new subspace

        if (ib.ne.1) then
          do ispin=1, kps%nspin
            ediff_back = abs(kps%el(ib - 1, ik, ispin) - kps%el(ib, ik, ispin))
          enddo
        endif

        if (ib.eq.1 .or. ediff_back.gt.TOL_Degeneracy) then

          ! at minimum band value in degenerate subspace

          call random_rotation(umtrx)

          do jb = ib, kps%mnband

            ! this difference checks to see if we are still in new subspace
            ! note: last index in el is always 1, otherwise code crashes
            if (jb.ne.kps%mnband) then
              do ispin=1, kps%nspin
                ediff_forw = abs(kps%el(jb, ik, ispin) - kps%el(jb+1, ik, ispin))
              enddo
            endif

            SAFE_ALLOCATE(wfn_temp,(kps%ngkmax,kps%nspin*kps%nspinor))
            wfn_temp(:,:)=wfn_new(:,jb,:)
            wfn_new(:,jb,:) = MATMUL(wfn_temp,umtrx)
            SAFE_DEALLOCATE(wfn_temp)
 
            if (ediff_forw .gt. TOL_Degeneracy) exit

          enddo !jb

        endif 
 
        call checknorm('new wfn',ib,ik,kps%ngk(ik),1,2,wfn_new(1:kps%ngk(ik),ib,:))

      enddo !ib

    endif

    !! Random rotation part ends here

    do ib=1,kps%mnband
      do jb=1,ib-1
        ovlp = ZERO
        do ispinor=1,kps%nspinor*kps%nspin
          ovlp=ovlp+blas_dot(kp%ngk(ik),(wfn_new(:,ib,ispinor)),1,(wfn_new(:,jb,ispinor)),1)
        enddo
        if (abs(ovlp).gt.TOL_Small) then
          write(6,998) ik,ispinor,ib,jb,ovlp
998       format(1x,"dot_prod error: ik =",i4,1x,"is =",i2,1x,"ib =",i4,1x,"jb =",i4,1x,"prod =",2f12.6)
        endif

      enddo ! jb
    enddo ! ib

    do ib = 1,kps%mnband
      call write_binary_data(8, kps%ngk(ik), kps%ngkmax, kps%nspin*kps%nspinor, wfn_new(1:kps%ngkmax,ib,1:kps%nspin*kps%nspinor))
    enddo ! ib
 
  enddo ! ik

  call close_file(7)
  call close_file(8)

  SAFE_DEALLOCATE_P(gvec%components)
  call kps%free()

  call dealloc_header_type(sheader, crys, kp)
  call dealloc_header_type(sheader, crys, kps)
  
  write(6,*) 'Done '
  
end program wfnmix_spinor
