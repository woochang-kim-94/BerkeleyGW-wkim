!=============================================================================
!
! Utilities:
!
! (1) wfnreduce        On 3.x wfn files By PWD  Last Modified 10/10/2011 (PWD)
!                      BM had done this for the old format
!                      but this is based on mf_convert.f90 
!     Reduces the cutoff of WFN files, Real/Complex flavor
!
!==============================================================================

#include "f_defs.h"

program wfnreduce

  use global_m
  use wfn_rho_vxc_io_m
  use fftw_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp, kpr
  type(gspace) :: gvec, gvecr, gvecRho
  character*3 :: sheader, scutoff
  character*7 :: sflavor, sbandcut
  character*11 :: inform, outform
  character*256 :: infile, outfile, usage
  integer :: ik, ib, iflavor, ig, id, FFTgrid(3)
  integer :: nargs, bandcut, ngRho
  real(DP), pointer :: dwfn(:,:), dwfnr(:,:)
  complex(DPC), pointer :: zwfn(:,:), zwfnr(:,:)
  logical :: informat, outformat
  real(DP) :: ecutw, ecutrho, gcutm, kenergy
!  integer, allocatable :: k2(:,:),k2wfc(:,:),k2wfc_red(:,:)
  integer, allocatable :: ngWfc(:) ! (kpoint)
  integer, allocatable :: wfcindex(:)

  usage = 'Usage: wfnreduce newCutoff(Ry) lastBand infile outfile'


! Get file names from command-line arguments

  nargs = command_argument_count()

  if (nargs .ne. 4) then
    call die(usage)
  endif

  call get_command_argument(1, scutoff)
  call get_command_argument(2, sbandcut)
  call get_command_argument(3, infile)
  call get_command_argument(4, outfile)

  read(scutoff, *) ecutw
  read(sbandcut, *) bandcut

  ecutrho = ecutw * 4

  sheader = 'GET'
  iflavor = -1
  informat = .false.
  inform = 'unformatted'
  outformat = .false.
  outform = 'unformatted'

! First we need to read the header and entire file to build a new kp and gvec
  call open_file(unit=7,file=TRUNC(infile),form=inform,status='old')
  call read_header_type(7, informat, sheader, iflavor, kp, gvec, syms, crys, dont_warn_kgrid=.true.)

  write(6,'(3x,"File header:",17x,a,/)') sheader // '-' // TRUNC(sflavor)
  write(6,'(3x,"Crystal volume:",f32.14)') crys%celvol
  write(6,'(3x,"Number of G-vectors:",i12)') gvec%ng
  write(6,'(3x,"Number of spins:",i16)') kp%nspin
  if (sheader .eq. 'WFN') then
    write(6,'(3x,"Number of bands:",i16)') kp%mnband
    write(6,'(3x,"Number of k-points:",i13)') kp%nrk
    write(6,'(3x,"ngkmax:",i13)') kp%ngkmax
    write(6,'(3x,"ifmin(1,1):",i13)') kp%ifmin(1,1)
    write(6,'(3x,"ifmax(1,1):",i13)') kp%ifmax(1,1)
  endif
  write(6,*)

  if (iflavor .eq. 1) then
    sflavor = RFLAVOR
  else
    sflavor = CFLAVOR
  endif

  kpr = kp !this might not be safe do to later allocates on points that are cloned here. 
           !I`m too green at f90 to know.
  kpr%ngkmax = 0 !we`ll need to find out what this acutally is at ecutw

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))

  call read_gvectors(7, informat, gvec%ng, gvec%ng, gvec%components)

  write(6,'(3x,"Calculating GVec energies")')
  gvecr%ng=0
  do ig=1,gvec%ng
    kenergy=DOT_PRODUCT(gvec%components(:,ig),MATMUL(crys%bdot,gvec%components(:,ig)))
    if (kenergy.lt.ecutrho) then 
      gvecr%ng=gvecr%ng+1 
    end if
  end do

  ngRho = gvecr%ng
  SAFE_ALLOCATE(gvecr%components,(3,gvecr%ng))

  gvecr%components(1:3,1:gvecr%ng)=gvec%components(1:3,1:gvecr%ng)
  gcutm = ecutrho/crys%blat**2

  do id=1,3
    FFTgrid(id)=int(2.0d0*sqrt(gcutm)*sqrt(sum(crys%avec(1:3,id)**2)))+1
    do while (.not. check_FFT_size(FFTgrid(id),3))
      FFTgrid(id)=FFTgrid(id)+1
    end do
  end do

  gvecr%FFTgrid = FFTgrid

  if(iflavor == 1) then
    if (sheader .eq. 'WFN') then
      SAFE_ALLOCATE(dwfn, (kp%ngkmax, kp%nspin))
    else
      SAFE_ALLOCATE(dwfn, (gvec%ng, kp%nspin))
    endif
  else
    if (sheader .eq. 'WFN') then
      SAFE_ALLOCATE(zwfn, (kp%ngkmax, kp%nspin*kp%nspinor))
    else
      SAFE_ALLOCATE(zwfn, (gvec%ng, kp%nspin))
    endif
  endif

  !this will be reduced in the new file 
  SAFE_ALLOCATE(kpr%ngk, (kp%nrk))
  SAFE_ALLOCATE(ngWfc, (kp%nrk))
  !I`m going to try and see if we can get away without making new versions of the other pointer
  !values in the reduced kpoint structure.

  if (sheader .eq. 'WFN') then
    do ik = 1, kp%nrk
      call read_gvectors(7, informat, kp%ngk(ik), kp%ngkmax, gvec%components)
      write(6,*) 'kp%ngk(ik) ', ik, kp%ngk(ik)
      SAFE_ALLOCATE(wfcindex, (kp%ngk(ik)))
      gvecr%ng=0
      do ig=1,gvec%ng
        kenergy=DOT_PRODUCT(gvec%components(:,ig),MATMUL(crys%bdot,gvec%components(:,ig)))
        !    write(6,'(3x,"GVec energy", i12, f32.14)') ig, kenergy
        if (kenergy.lt.ecutw) then 
          gvecr%ng=gvecr%ng+1 
          wfcindex(gvecr%ng)=ig
        end if
      end do
      kpr%ngk(ik) = gvecr%ng
      if(kpr%ngkmax.lt.kpr%ngk(ik)) then
        write(6,*) 'updating kpr%ngkmax', ik, kpr%ngk(ik) 
        kpr%ngkmax = kpr%ngk(ik)
      end if
      ngWfc(ik) = gvecr%ng
      !skip this frist time through
!       do ig=1,gvecr%ng
!         gvecr%components(1,ig)=gvec%components(1,(wfcindex(ig)))
!         gvecr%components(2,ig)=gvec%components(2,(wfcindex(ig)))
!         gvecr%components(3,ig)=gvec%components(3,(wfcindex(ig)))
!       end do

!       call write_gvectors(8, outformat, gvecr%ng, gvecr%ng, gvecr%components)

      do ib = 1, kp%mnband
        if(iflavor == 1) then
          SAFE_ALLOCATE(dwfnr, (gvecr%ng,kp%nspin))
          call read_real_data(7, informat, kp%ngk(ik), kp%ngkmax, kp%nspin, dwfn)
!          call write_real_data(8, outformat, gvecr%ng, gvecr%ng, kp%nspin, dwfn)
          SAFE_DEALLOCATE_P(dwfnr)
        else
          SAFE_ALLOCATE(zwfnr, (gvecr%ng,kp%nspin*kp%nspinor))
          call read_complex_data(7, informat, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, zwfn)
!          write(6,*) 'managed to read wave function'
!          write(6,*) 'test access zwfn', zwfn(1,1)
!          call reduce_complex_wfn(gvec%ng, gvecr%ng, kp%nspin, zwfn, zwfnr)
!           call write_complex_data(8, outformat, gvecr%ng, gvecr%ng, kp%nspin*kp%nspinor, zwfnr)
          SAFE_DEALLOCATE_P(zwfnr)
        endif
      enddo
      SAFE_DEALLOCATE(wfcindex)
    enddo
  else
    if(iflavor == 1) then
      call read_real_data(7, informat, gvec%ng, gvec%ng, kp%nspin, dwfn)
      SAFE_ALLOCATE(dwfnr, (gvecr%ng,kp%nspin))
!      call write_real_data(8, outformat, kp%ngk(ik), gvecr%ng, kp%nspin, dwfn)
      SAFE_DEALLOCATE_P(dwfnr)
    else
      call read_complex_data(7, informat, gvec%ng, gvec%ng, kp%nspin*kp%nspinor, zwfn)
      SAFE_ALLOCATE(zwfnr, (gvecr%ng,kp%nspin*kp%nspinor))
!      call reduce_complex_wfn(gvec%ng, gvecr%ng, kp%nspin*kp%nspinor, zwfn, zwfnr)
!      call write_complex_data(8, outformat, gvecr%ng, gvecr%ng, kp%nspin*kp%nspinor, zwfnr)
      SAFE_DEALLOCATE_P(zwfnr)
    endif
  endif

  kpr%mnband = bandcut
  kpr%ecutwfc = ecutw

!okay now hopefully we know everything we need to for kpr
  call close_file(7)

  write(6,*) 'reopening file.'

! close file and open it and start over from scratch
  call open_file(unit=7,file=TRUNC(infile),form=inform,status='old')
  call open_file(unit=8,file=TRUNC(outfile),form=outform,status='replace')

  call read_header_type(7, informat, sheader, iflavor, kp, gvec, syms, crys, dont_warn_kgrid=.true.)

!  outformat = .true.
  SAFE_ALLOCATE(gvecRho%components, (3, ngRho))

  call read_gvectors(7, informat, gvec%ng, gvec%ng, gvec%components)

  gvecRho%components(1:3,1:gvecRho%ng)=gvec%components(1:3,1:gvecRho%ng)
  gcutm = ecutrho/crys%blat**2

  do id=1,3
    FFTgrid(id)=int(2.0d0*sqrt(gcutm)*sqrt(sum(crys%avec(1:3,id)**2)))+1
    do while (.not. check_FFT_size(FFTgrid(id),3))
      FFTgrid(id)=FFTgrid(id)+1
    end do
  end do

  gvecRho%FFTgrid = FFTgrid

  call write_header_type(8, outformat, sheader, iflavor, kpr, gvecr, syms, crys)
  write(6,*) 'writing base gvectors' 
  call write_gvectors(8, outformat, gvecr%ng, gvecr%ng, gvecr%components)
  write(6,*) 'wrote base gvectors'

  if(iflavor == 1) then
    if (sheader .eq. 'WFN') then
      SAFE_ALLOCATE(dwfn, (kp%ngkmax, kp%nspin))
    else
      SAFE_ALLOCATE(dwfn, (gvec%ng, kp%nspin))
    endif
  else
    if (sheader .eq. 'WFN') then
      SAFE_ALLOCATE(zwfn, (kp%ngkmax, kp%nspin*kp%nspinor))
    else
      SAFE_ALLOCATE(zwfn, (gvec%ng, kp%nspin*kp%nspinor))
    endif
  endif

  if (sheader .eq. 'WFN') then
    do ik = 1, kp%nrk
      call read_gvectors(7, informat, kp%ngk(ik), kp%ngkmax, gvec%components)
      SAFE_ALLOCATE(wfcindex, (kp%ngk(ik)))
      gvecr%ng=0
      do ig=1,gvec%ng
        kenergy=DOT_PRODUCT(gvec%components(:,ig),MATMUL(crys%bdot,gvec%components(:&
         &,ig)))
        if (kenergy.lt.ecutw) then 
          gvecr%ng=gvecr%ng+1 
          wfcindex(gvecr%ng)=ig
        end if
      end do
      gvecr%components(1:3,1:gvecr%ng)=gvec%components(1:3,(wfcindex(1:gvecr%ng)))
      call write_gvectors(8, outformat, gvecr%ng, gvecr%ng, gvecr%components)

      do ib = 1, kp%mnband
        if(ib.le.kpr%mnband) then
          if(iflavor == 1) then
            SAFE_ALLOCATE(dwfnr, (gvecr%ng,kp%nspin))
            call read_real_data(7, informat, kp%ngk(ik), kp%ngkmax, kp%nspin, dwfn)
            call write_real_data(8, outformat, gvecr%ng, gvecr%ng, kp%nspin, dwfn)
            SAFE_DEALLOCATE_P(dwfnr)
          else
            SAFE_ALLOCATE(zwfnr, (gvecr%ng,kp%nspin*kp%nspinor))
            call read_complex_data(7, informat, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, zwfn)
            !          write(6,*) 'managed to read wave function'
            !          write(6,*) 'test access zwfn', zwfn(1,1)
            call reduce_complex_wfn(gvec%ng, gvecr%ng, kp%nspin, kp%nspinor, zwfn, zwfnr)
            call write_complex_data(8, outformat, gvecr%ng, gvecr%ng, kp%nspin*kp%nspinor, zwfnr)
            SAFE_DEALLOCATE_P(zwfnr)
          endif
        else
          if(iflavor == 1) then
            SAFE_ALLOCATE(dwfnr, (gvecr%ng,kp%nspin))
            call read_real_data(7, informat, kp%ngk(ik), kp%ngkmax, kp%nspin, dwfn)
            SAFE_DEALLOCATE_P(dwfnr)
          else
            SAFE_ALLOCATE(zwfnr, (gvecr%ng,kp%nspin*kp%nspinor))
            call read_complex_data(7, informat, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, zwfn)
            SAFE_DEALLOCATE_P(zwfnr)
          endif
        endif
      enddo
    enddo
  else
    if(iflavor == 1) then
      call read_real_data(7, informat, gvec%ng, gvec%ng, kp%nspin, dwfn)
      SAFE_ALLOCATE(dwfnr, (gvecr%ng,kp%nspin))
      call write_real_data(8, outformat, kp%ngk(ik), gvecr%ng, kp%nspin, dwfn)
      SAFE_DEALLOCATE_P(dwfnr)
    else
      call read_complex_data(7, informat, gvec%ng, gvec%ng, kp%nspin*kp%nspinor, zwfn)
      SAFE_ALLOCATE(zwfnr, (gvecr%ng,kp%nspin))
      call reduce_complex_wfn(gvec%ng, gvecr%ng, kp%nspin, kp%nspinor, zwfn, zwfnr)
      call write_complex_data(8, outformat, gvecr%ng, gvecr%ng, kp%nspin*kp%nspinor, zwfnr)
      SAFE_DEALLOCATE_P(zwfnr)
    endif
  endif

!   if(iflavor == 1) then
!     SAFE_DEALLOCATE_P(dwfn)
!   else
!     SAFE_DEALLOCATE_P(zwfn)
!   endif

  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE_P(gvecr%components)
   

!   call dealloc_header_type(sheader, crys, kp)

!   call close_file(7)
!   call close_file(8)

!   write(6,'(3x,"Done",/)')

contains

subroutine reduce_complex_wfn(ngo, ngr, nspin, nspinor, zwfn, zwfnr)

  use blas_m
  implicit none

  integer, intent(in) :: ngo, ngr, nspin, nspinor
  complex(DPC), pointer, intent(in) :: zwfn(:,:)
  complex(DPC), pointer, intent(inout) :: zwfnr(:,:)

  real(DP) :: norm(2)
  integer :: ispin, ispinor

  PUSH_SUB(reduce_complex_wfn)

!  write(6,*) 'ngo, ngr, nspin ', ngo, ngr, nspin


! I`m pretending that in the new format the problem that caused wfcindex in the original
! wfnreduce doesn`t exist.  But it may. Must research

  do ispin = 1, nspin
    if (nspinor.eq.2) then
      norm(ispin) = 0.0D0
      do ispinor = 1, nspinor
        norm(ispin) = norm(ispin) +  blas_nrm2(ngr, zwfn(:,ispin*ispinor), 1)
      enddo
      write(6,*) 'Wave function has norm ', norm(ispin)
      do ispinor = 1, nspinor
        call zscal(ngr, COMPLEXIFY(ONE/norm(ispin)), zwfn(:,ispin*ispinor), 1)
        call zcopy(ngr, zwfn(:, ispin*ispinor), 1, zwfnr(:, ispin*ispinor), 1) 
      end do
    else
      norm(ispin) = blas_nrm2(ngr, zwfn(:,ispin), 1)
      write(6,*) 'Wave function has norm ', norm(ispin)
      call zscal(ngr, COMPLEXIFY(ONE/norm(ispin)), zwfn(:,ispin), 1)
      call zcopy(ngr, zwfn(:, ispin), 1, zwfnr(:, ispin), 1) 
    end if
  end do

  do ispin = 1, nspin
! I`m pretending that in the new format the problem that caused wfcindex in the original
! wfnreduce doesn`t exist.  But it may. Must research
    if (nspinor.eq.2) then
      norm(ispin) = 0.0D0
      do ispinor = 1, nspinor
        norm(ispin) = norm(ispin) + blas_nrm2(ngr, zwfnr(:,ispin*ispinor), 1)
      enddo    
      write(6,*) 'Wave function has norm ', norm(ispin)
    else
      norm(ispin) = blas_nrm2(ngr, zwfnr(:,ispin), 1)
      write(6,*) 'Wave function has norm ', norm(ispin)
    end if
  end do

  POP_SUB(reduce_complex_wfn)

end subroutine reduce_complex_wfn

    
end program wfnreduce

