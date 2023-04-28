!*************************************************************************

Program EPM2BGW
!*************************************************************************
! EPM calculates the band structure along various directions in
! k-space using a plane-wave basis and a fixed effective potential.
! Several choices of empirical pseudopotentials are provided.
!*************************************************************************
! Written by William Mattson and Richard M. Martin of the University
! of Illinois, based upon program written by K. Glassford and I. Souza.
! Modified by G. Samsonidze for use with BerkeleyGW (August 2008).
!*************************************************************************
! Note: atomic units (a.u.) are used throughout the program.
!*************************************************************************

#include "f_defs.h"

  use sysParams,           only : double, zero, one
  use TagHandlerMod,       only : tagHandlerT, TagHandlerInit, TagHandlerDestroy, &
                                  & FindTag
  use StructureMod,        only : StructureT, StructureInit, StructureDestroy
  use kPointsMod,          only : kPointsT, KPointsInit, KPointsDestroy
  use graphMod,            only : graphT, PlotBands
  use hamiltonianMod,      only : hamiltonianArrayT, HInit, HDestroy, HDiagonalize, & 
                                  & HPrint, HCheck
  use eigenStatesMod,      only : eigenStatesT, EigenStatesInit, EigenStatesDestroy, &
                                  & PrintEigenvalues
  use pwHamMod,            only : hamInfoT, HamInfoInit, HamInfoDestroy, &
                                  & HGenerate, FindHMaxSize, GridInit, &
                                  & PotentialRGenerate, KineticGenerate, &
                                  & CalculateDensity, GridDestroy
  use cgParamsMod,         only : cgParamsT, cgParamsInit
  use typeMod,             only : GridT
  use ConjGradMod,         only : cgEigenSystem,cgFFTDiagonalise
  use densityArrayMod,     only : densityArrayT, DensityArrayInit, &
                                  & DensityArrayDestroy

  use nrtype_m
  use message_m
  use wfn_rho_vxc_io_m,    only : write_binary_header, write_binary_gvectors, &
                                  & write_binary_real_data, write_binary_complex_data, write_matrix_elements
  use fftw_m,              only : check_FFT_size
  use check_inversion_m,   only : check_inversion
  use symmetries_m,        only : get_symmetries

  implicit none

  type(tagHandlerT),       pointer :: tagHandler
  type(hamInfoT),          pointer :: hamInfo
  type(StructureT),        pointer :: structure1
  type(kPointsT),          pointer :: kPoints
  type(eigenStatesT),      pointer :: eigenStates
  type(hamiltonianArrayT), pointer :: hamiltonian
  type(graphT),            pointer :: graphinfo
  type(cgParamsT)                  :: cgParams
  type(GridT),             pointer :: Grid
  type(densityArrayT),     pointer :: densityArray

  integer, allocatable :: gVectorsK(:,:)
  real(double), allocatable :: eigenVectorsReal(:,:,:,:)
  complex(double), allocatable :: eigenVectorsComplex(:,:,:,:)
  real(double), allocatable :: dataReal(:,:)
  complex(double), allocatable :: dataComplex(:,:), vxc_mtxel(:,:)

  integer, pointer :: atyp(:)
  real(double), pointer :: apos(:,:)
  integer, allocatable :: ngk(:)
  real(double), allocatable :: kw(:)
  real(double), allocatable :: kpt(:,:)
  integer, allocatable :: ifmin(:,:)
  integer, allocatable :: ifmax(:,:)
  real(double), allocatable :: energies(:,:,:)
  real(double), allocatable :: occupations(:,:,:)
  logical :: wfng_flag, rhog_flag, vxcg_flag, vxc_flag, v_of_q, disable_symmetries
  integer :: real_or_complex, ReportNo, error, i, j, k, &
   & ik, ispin, ib, ig, cell_symmetry, nat, nsym, nspin, &
   & vxc_diag_nmin, vxc_diag_nmax, vxc_offdiag_nmin, &
   & vxc_offdiag_nmax, ndiag, noffdiag, FFTgrid(3), &
   & rotation(3, 3, 48), spin_index(2), idiag, ioffdiag, spacegroup, &
   & identity(3,3)
  integer, allocatable :: diag(:), offdiag(:)
  real(double) :: celvol, recvol, al, bl, ecutwfn, ecutrho, &
   & gcutm, abstol, temperature, efermi, &
   & a(3, 3), b(3, 3), adot(3, 3), bdot(3, 3), translation(3, 48), ff(6)
  character(len=3) :: sheader
  character(len=256) :: wfng_file, rhog_file, vxcg_file, vxc_file, &
   & filename, tmpstr
  character*21 :: symbol

  integer, parameter :: Nfac = 3

! Initialization

  nullify(tagHandler)
  nullify(hamInfo)
  nullify(structure1)
  nullify(kPoints)
  nullify(eigenStates)
  nullify(hamiltonian)
  nullify(graphinfo)
  nullify(Grid)
  nullify(densityArray)

  nullify(atyp)
  nullify(apos)

  ReportNo = 110

  filename = 'input.tmp'
  call open_file(unit = 10, file = trim(filename), status = 'replace')
  i = 0
  do while (i .eq. 0)
     read(unit = 5, fmt = 999, iostat = i) tmpstr
     write(unit = 10, fmt = 999) tmpstr
  enddo
  call close_file(unit = 10)

  ! Necessary if using conjugate gradient method
  call Random_Seed

  ! Initialize "tag handler" to read the input from file 'filename'.
  call TagHandlerInit( tagHandler, 10, filename)

  ! Input parameter real_or_complex for BerkeleyGW files
  call FindTag( tagHandler, 'real_or_complex', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) real_or_complex
     if(error .eq. 0) then
        if(real_or_complex .ne. 1 .and. real_or_complex .ne. 2) error = 1
     endif
  endif
  if(error .ne. 0) real_or_complex = 2
  write(*,'(" real_or_complex ",i1)')real_or_complex

  ! Input parameter disable_symmetries for BerkeleyGW wavefunction file
  call FindTag( tagHandler, 'disable_symmetries', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) disable_symmetries
  endif
  if(error .ne. 0) disable_symmetries = .false.
  write(*,'(" disable_symmetries ",l1)')disable_symmetries

  ! Input parameter wfng_flag for BerkeleyGW wavefunction file
  call FindTag( tagHandler, 'wfng_flag', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) wfng_flag
  endif
  if(error .ne. 0) wfng_flag = .false.
  write(*,'(" wfng_flag ",l1)')wfng_flag

  ! Input parameter wfng_file for BerkeleyGW wavefunction file
  call FindTag( tagHandler, 'wfng_file', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) wfng_file
  endif
  if(error .ne. 0) write(wfng_file,'("WFN")')
  write(*,'(" wfng_file ",a)')trim(wfng_file)

  ! Input parameter rhog_flag for BerkeleyGW charge density file
  call FindTag( tagHandler, 'rhog_flag', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) rhog_flag
  endif
  if(error .ne. 0) rhog_flag = .false.
  write(*,'(" rhog_flag ",l1)')rhog_flag

  ! Input parameter rhog_file for BerkeleyGW charge density file
  call FindTag( tagHandler, 'rhog_file', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) rhog_file
  endif
  if(error .ne. 0) write(rhog_file,'("RHO")')
  write(*,'(" rhog_file ",a)')trim(rhog_file)

  ! Input parameter vxcg_flag for BerkeleyGW exchange-correlation potential file
  call FindTag( tagHandler, 'vxcg_flag', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) vxcg_flag
  endif
  if(error .ne. 0) vxcg_flag = .false.
  write(*,'(" vxcg_flag ",l1)')vxcg_flag

  ! Input parameter vxcg_file for BerkeleyGW exchange-correlation potential file
  call FindTag( tagHandler, 'vxcg_file', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) vxcg_file
  endif
  if(error .ne. 0) write(vxcg_file,'("VXC")')
  write(*,'(" vxcg_file ",a)')trim(vxcg_file)

  ! Input parameter vxc_flag for BerkeleyGW exchange-correlation matrix element file
  call FindTag( tagHandler, 'vxc_flag', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) vxc_flag
  endif
  if(error .ne. 0) vxc_flag = .false.
  write(*,'(" vxc_flag ",l1)')vxc_flag

  ! Input parameter vxc_file for BerkeleyGW exchange-correlation matrix element file
  call FindTag( tagHandler, 'vxc_file', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) vxc_file
  endif
  if(error .ne. 0) write(vxc_file,'("vxc.dat")')
  write(*,'(" vxc_file ",a)')trim(vxc_file)

  ! Input parameter vxc_diag_nmin for BerkeleyGW exchange-correlation matrix element file
  call FindTag( tagHandler, 'vxc_diag_nmin', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) vxc_diag_nmin
  endif
  if(error .ne. 0) vxc_diag_nmin = 0
  write(*,'(" vxc_diag_nmin ",i2)')vxc_diag_nmin

  ! Input parameter vxc_diag_nmax for BerkeleyGW exchange-correlation matrix element file
  call FindTag( tagHandler, 'vxc_diag_nmax', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) vxc_diag_nmax
  endif
  if(error .ne. 0) vxc_diag_nmax = 0
  write(*,'(" vxc_diag_nmax ",i2)')vxc_diag_nmax

  ! Input parameter vxc_offdiag_nmin for BerkeleyGW exchange-correlation matrix element file
  call FindTag( tagHandler, 'vxc_offdiag_nmin', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) vxc_offdiag_nmin
  endif
  if(error .ne. 0) vxc_offdiag_nmin = 0
  write(*,'(" vxc_offdiag_nmin ",i2)')vxc_offdiag_nmin

  ! Input parameter vxc_offdiag_nmax for BerkeleyGW exchange-correlation matrix element file
  call FindTag( tagHandler, 'vxc_offdiag_nmax', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) vxc_offdiag_nmax
  endif
  if(error .ne. 0) vxc_offdiag_nmax = 0
  write(*,'(" vxc_offdiag_nmax ",i2)')vxc_offdiag_nmax

  ! for nspin = 2, a second copy of the wavefunctions and eigenvalues will be written
  call FindTag( tagHandler, 'nspin', error)
  if(error .eq. 0) then
    read(unit = tagHandler%fileno, fmt = *, iostat = error) nspin
  endif
  if(error .ne. 0) nspin = 1
  write(*,'(" nspin ",i2)') nspin

  ! electronic temperature
  call FindTag( tagHandler, 'temperature', error)
  if(error .eq. 0) then
    read(unit = tagHandler%fileno, fmt = *, iostat = error) temperature
  endif
  if(error .ne. 0) temperature = 0d0
  write(*,'(" temperature ",f13.6," Ry")') temperature
  temperature = temperature * 0.5d0 ! convert to Ha

  ! Use form-factors from input file or hard-coded V(q) potentials
  call FindTag( tagHandler,'FormFactors', error)
  if(error .eq. 0) then
    read(unit = tagHandler%fileno, fmt = *, iostat = error) ff
    if(error /= 0) stop 'Error, could not read FormFactors.'
    v_of_q = .false.
    write(*,'("Using Form Factors V_3,8,11^S & V_3,4,11^A")')
  else
    ff(:) = 0.0d0
    v_of_q = .true.
    write(*,'("Using Potentials V(q) for Species")')
  endif

  ! Read LAPACK Absolute Tolerance
  call FindTag( tagHandler, 'AbsoluteTolerance', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) abstol
  endif
  if(error .ne. 0) abstol = -one
  write(*,'(" AbsoluteTolerance ",e13.6)')abstol

  ! Read the Conjugate Gradient Parameter
  call cgParamsInit(cgParams,tagHandler)

  ! This will call DimensionsInit, LatticeInit, and AtomsInit.
  call StructureInit( structure1, tagHandler )

  call KPointsInit( kPoints, graphInfo, structure1%ndim, structure1%lattice, tagHandler)

  ! Gvector list and structurefactor
  call HamInfoInit( hamInfo, structure1, tagHandler, v_of_q, ff )

  ! Find hMaxSize = max size of H at any k point
  call FindHMaxSize( hamInfo, structure1, kPoints )
  
  ! Allocates array for H and work arrays
  call HInit( hamiltonian, hamInfo%hMaxSize )
  hamiltonian%abstol = abstol

  ! Allocates array for eigenStates array
  call EigenStatesInit( eigenStates, kPoints%numBands, kPoints%numKPoints, &
                      & hamInfo%hMaxSize, -1.0d0, 1)
  eigenStates%eigenVectors = zero

  ! Initialise the FFT Grid and Plans
  Call GridInit(Grid, hamInfo, eigenStates, structure1%ndim)

! DAS -- do not write this useless file if not using CG
  if(cgParams%Switch.eq.1) then
    call open_file(unit=ReportNo,file='Report.cg',status='replace',position='rewind')
    write(ReportNo,*) 'Eigenvalues:'
    write(ReportNo,*) 'The FFT Grid Size +/-',Grid%Size(1)
    write(ReportNo,*) 'The number of GVectors:',Grid%NumgVec
    write(ReportNo,*) 'Eigenvalues:'
  endif

  ! Calculate the Potential and FFT to Vr
  Call PotentialRGenerate(hamInfo, Grid)

  if(cgParams%Switch .eq. 1) then 
     do k = 1, kPoints%numKPoints 
        call HGenerate( hamInfo, structure1, kPoints, eigenStates, &
             & k, hamiltonian%h, hamiltonian%sqH, hamiltonian%hSize )

        write(ReportNo, *) 'k point ', k

        call cgEigenSystem( hamiltonian, eigenStates, k, cgParams%tol, &
             & cgParams%period ,ReportNo )

        write(*, '(" k point ", i4, " out of ", i4)') k, kPoints%numKPoints

     end do ! k
  else
     do k = 1, kPoints%numKPoints

        call HGenerate( hamInfo, structure1, kPoints, eigenStates, &
             & k, hamiltonian%h, hamiltonian%sqH, hamiltonian%hSize )

! DAS -- do not write this useless file
!        write(ReportNo, *) 'k point ', k
     
        call HDiagonalize( hamiltonian, eigenStates, k )

! DAS -- make sure stupid LAPACK routines have genuinely solved the Hamiltonian, else die 
#ifdef DEBUG
        call HCheck( hamiltonian, eigenStates, k )
#endif

        write(*, '(" k point ", i4, " out of ", i4)') k, kPoints%numKPoints

     end do ! k
  end if

  if(cgParams%Switch.eq.1) then 
    call close_file(unit=ReportNo)
  endif

! Initialize data structures for BerkeleyGW output

  a(:,:) = structure1%lattice%aLatVec(:,:)
  al = sqrt(a(1,1)**2 + a(2,1)**2 + a(3,1)**2)
  a(:,:) = a(:,:) / al
  celvol = a(1,1) * (a(2,2) * a(3,3) - a(2,3) * a(3,2)) - &
           a(2,1) * (a(1,2) * a(3,3) - a(1,3) * a(3,2)) + &
           a(3,1) * (a(1,2) * a(2,3) - a(1,3) * a(2,2))
  bl = 2.0d0 * PI_D / al
  b(1,1) = (a(2,2) * a(3,3) - a(3,2) * a(2,3)) / celvol
  b(2,1) = (a(3,2) * a(1,3) - a(1,2) * a(3,3)) / celvol
  b(3,1) = (a(1,2) * a(2,3) - a(2,2) * a(1,3)) / celvol
  b(1,2) = (a(2,3) * a(3,1) - a(3,3) * a(2,1)) / celvol
  b(2,2) = (a(3,3) * a(1,1) - a(1,3) * a(3,1)) / celvol
  b(3,2) = (a(1,3) * a(2,1) - a(2,3) * a(1,1)) / celvol
  b(1,3) = (a(2,1) * a(3,2) - a(3,1) * a(2,2)) / celvol
  b(2,3) = (a(3,1) * a(1,2) - a(1,1) * a(3,2)) / celvol
  b(3,3) = (a(1,1) * a(2,2) - a(2,1) * a(1,2)) / celvol
  celvol = abs(celvol) * al**3
  recvol = (2.0d0 * PI_D)**3 / celvol
  adot(:,:) = 0.0d0
  do i=1,3
    do j=1,3
      do k=1,3
        adot(j,i) = adot(j,i) + a(k,j) * a(k,i)
      enddo
    enddo
  enddo
  adot(:,:) = adot(:,:) * al**2
  bdot(:,:) = 0.0d0
  do i=1,3
    do j=1,3
      do k=1,3
        bdot(j,i) = bdot(j,i) + b(k,j) * b(k,i)
      enddo
    enddo
  enddo
  bdot(:,:) = bdot(:,:) * bl**2

  ! kinetic energy cutoff for wave functions
  ! convert from Hartree to Rydberg
  ecutwfn=2.0d0*hamInfo%gVectors%energyCutOff
  ! kinetic energy cutoff for charge density and potential
  ecutrho=4.0d0*ecutwfn
  gcutm=ecutrho/bl**2

  call FindTag( tagHandler, 'FFTGrid', error)
  if(error .eq. 0) then
    read(unit = tagHandler%fileno, fmt = *, iostat = error) (FFTgrid(j),j=1,3)
    if(error /= 0) stop 'Error, could not read FFTgrid.'
    if(any(FFTgrid(:) <= 0)) stop 'Error, FFTGrid values must be positive.'
  else
    do i=1,3
      FFTgrid(i)=int(2.0d0*sqrt(gcutm)*sqrt(a(1,i)**2+a(2,i)**2+a(3,i)**2))+1
      do while (.not.check_FFT_size(FFTgrid(i), Nfac))
        FFTgrid(i)=FFTgrid(i)+1
      enddo
    enddo
  endif
  write(*,'(a,3i6)') 'FFT grid: ', FFTgrid(:)

  nat=structure1%atomBasis%totalNumAtoms
  SAFE_ALLOCATE(atyp, (nat))
  SAFE_ALLOCATE(apos, (3, nat))
  atyp(:) = structure1%atomBasis%atomicNumber(structure1%atomBasis%speciesOfAtom(:))
  apos(:,:) = structure1%atomBasis%atomCart(:,:)/al

! gsm: if using form-factors make sgam_at think we have different species
  if (.not. v_of_q .and. abs(ff(4)) + abs(ff(5)) + abs(ff(6)) .gt. TOL_Small) &
    atyp(1) = atyp(1) + 1000

  if(disable_symmetries) then
    cell_symmetry = 0
    nsym = 1
    write(*,'(a)') 'Symmetries have been disabled.'
    identity = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), shape(identity))
    rotation(:,:,nsym) = identity
  else
    call get_symmetries(nat, atyp, structure1%atomBasis%atomLat, a, FFTgrid, cell_symmetry, nsym, &
      rotation, translation, spacegroup, symbol)
    write(*,'(a,i3,a,a)') 'Space group ', spacegroup, ', symbol ', trim(symbol)
  endif

! gsm: restore atomic species back to normal
     if (.not. v_of_q .and. abs(ff(4)) + abs(ff(5)) + abs(ff(6)) .gt. TOL_Small) &
   atyp(1) = atyp(1) - 1000

  call check_inversion(real_or_complex, nsym, rotation(:,:,:), nspin, .true., .true., tnp = translation)

  if (wfng_flag) then

     write(*,*)

     SAFE_ALLOCATE(ngk, (kPoints%numKPoints))
     SAFE_ALLOCATE(kw, (kPoints%numKPoints))
     SAFE_ALLOCATE(kpt, (3, kPoints%numKPoints))
     ngk(1:kPoints%numKPoints)=eigenStates%numBasisVectors(1:kPoints%numKPoints)
     kw(1:kPoints%numKPoints)=kPoints%kPointWeights(1:kPoints%numKPoints)
     kpt(1:3,1:kPoints%numKPoints)=kPoints%kPointsLat(1:3,1:kPoints%numKPoints)

     SAFE_ALLOCATE(ifmin, (kPoints%numKPoints,nspin))
     SAFE_ALLOCATE(ifmax, (kPoints%numKPoints,nspin))
     ifmin(1:kPoints%numKPoints,1:nspin)=1
     ifmax(1:kPoints%numKPoints,1:nspin)=kPoints%numOccupiedBands
     

     ! Convert eigenvalues from Hartree to Rydberg for BerkeleyGW
     SAFE_ALLOCATE(energies, (kPoints%numBands, kPoints%numKPoints, nspin))
     do ispin=1,nspin
        do ik=1,kPoints%numKPoints
           do ib=1,kPoints%numBands
              energies(ib,ik,ispin)=eigenStates%eigenValues(ib,ik)*2.0d0
           enddo
        enddo
     enddo

     if(temperature > TOL_Small) then
       efermi = (minval(energies(kPoints%numOccupiedBands+1:,:,:)) + maxval(energies(1:kPoints%numOccupiedBands,:,:)))*0.5d0
       write(*,'(a,f13.6,a)') 'Fermi energy (mid-gap): ', efermi*2d0, ' Ry'
     endif
     
     SAFE_ALLOCATE(occupations, (kPoints%numBands, kPoints%numKPoints, nspin))
     do ispin=1,nspin
       do ik=1,kPoints%numKPoints
         do ib=1,kPoints%numBands
           if(temperature > TOL_Small) then
             occupations(ib, ik, ispin) = 1 / (exp((energies(ib,ik,ispin) - efermi)/temperature) + 1)
           else
             if(ib < ifmin(ik, ispin) .or. ib > ifmax(ik, ispin)) then
               occupations(ib, ik, ispin) = 0.0d0
             else
               occupations(ib, ik, ispin) = 1.0d0
             endif
           endif
         enddo
       enddo
     enddo

     ! Convert eigenvalues from Hartree to Rydberg for stdout
     do ik=1,kPoints%numKPoints
        do ib=1,kPoints%numBands
!           eigenStates%eigenValues(ib,ik)=eigenStates%eigenValues(ib,ik)*2.0d0*RYD
           eigenStates%eigenValues(ib,ik)=eigenStates%eigenValues(ib,ik)*2.0d0
        enddo
        call PrintEigenvalues(eigenStates,ik)
        write(*,*)
     enddo

     SAFE_ALLOCATE(gVectorsK, (3, eigenStates%maxBasisVectors))

     if(real_or_complex .eq. 1) then
        write(*,'(" Applying Gram-Schmidt process...")')
        SAFE_ALLOCATE(eigenVectorsReal, (eigenStates%maxBasisVectors, nspin, eigenStates%numStates, eigenStates%numKPoints))
        call real_wfng(eigenStates%numKPoints, eigenStates%numStates, &
           & nspin, eigenStates%maxBasisVectors, eigenStates%numBasisVectors, &
           & eigenStates%eigenValues, eigenStates%eigenVectors, eigenVectorsReal)
        write(*,'(" ...done!")')
        write(*,*)
     else
        SAFE_ALLOCATE(eigenVectorsComplex, (eigenStates%maxBasisVectors, nspin, eigenStates%numStates, eigenStates%numKPoints))
        do ik=1,kPoints%numKPoints
           do ib=1,kPoints%numBands
              do ispin=1,nspin
                 do ig=1,eigenStates%numBasisVectors(ik)
                    eigenVectorsComplex(ig, ispin, ib, ik) = &
                       eigenStates%eigenVectors(ig, ib, ik)
                 enddo
              enddo
           enddo
        enddo
     endif

     write(*,'(" Writing BerkeleyGW wavefunction file...")')

     call open_file(20,file=trim(wfng_file),status='replace',form='unformatted')

     sheader = 'WFN'
     call write_binary_header(20, sheader, real_or_complex, nspin, &
        hamInfo%gVectors%numGVectors, nsym, cell_symmetry, nat, &
        kPoints%numKPoints, kPoints%numBands, eigenStates%maxBasisVectors, &
        ecutrho, ecutwfn, FFTgrid, kPoints%kPointNumbers, kPoints%kPointOffsets, &
        celvol, al, a, adot, recvol, bl, b, bdot, rotation, translation, &
        atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations)

     call write_binary_gvectors(20, hamInfo%gVectors%numGVectors, &
        hamInfo%gVectors%numGVectors, hamInfo%gVectors%gVectors(:,:))

     do ik=1,kPoints%numKPoints
        do ig=1,eigenStates%numBasisVectors(ik)
           gVectorsK(:,ig)= &
              hamInfo%gVectors%gVectors(:,eigenStates%basisIndex(ig,ik))
        enddo
        call write_binary_gvectors(20, eigenStates%numBasisVectors(ik), &
           eigenStates%maxBasisVectors, gVectorsK(:,:))

        do ib=1,kPoints%numBands
           if(real_or_complex .eq. 1) then
              call write_binary_real_data(20, eigenStates%numBasisVectors(ik), &
                 eigenStates%maxBasisVectors, nspin, eigenVectorsReal(:, :, ib, ik))
           else
              call write_binary_complex_data(20, eigenStates%numBasisVectors(ik), &
                 eigenStates%maxBasisVectors, nspin, eigenVectorsComplex(:, :, ib, ik))
           endif
        enddo

     enddo

     call close_file(20)

     write(*,'(" ...done!")')

     SAFE_DEALLOCATE(gVectorsK)
     if (real_or_complex .eq. 1) then
        SAFE_DEALLOCATE(eigenVectorsReal)
     else
        SAFE_DEALLOCATE(eigenVectorsComplex)
     endif

  endif

  if (rhog_flag) then

     write(*,*)

     write(*,'(" Calculating charge density...")')

     call DensityArrayInit( DensityArray, structure1%ndim, &
        & FFTgrid, hamInfo%gVectors%numGVectors)

     call CalculateDensity( hamInfo, kPoints, eigenStates, densityArray )

     write(*,'(" ...done!")')

     write(*,*)

     if(real_or_complex .eq. 1) then
        SAFE_ALLOCATE(dataReal, (hamInfo%gVectors%numGVectors, nspin))
        do ispin=1,nspin
           do ig=1,hamInfo%gVectors%numGVectors
              dataReal(ig,ispin)= &
                 dble(densityArray%rhog(ig)*cmplx(2.0d0/dble(nspin),0.0d0,kind=double))
           enddo
        enddo
     else
        SAFE_ALLOCATE(dataComplex, (hamInfo%gVectors%numGVectors, nspin))
        do ispin=1,nspin
           do ig=1,hamInfo%gVectors%numGVectors
              dataComplex(ig,ispin)= &
                 densityArray%rhog(ig)*cmplx(2.0d0/dble(nspin),0.0d0,kind=double)
           enddo
        enddo
     endif

     write(*,'(" Writing BerkeleyGW charge density file...")')

     call open_file(20,file=trim(rhog_file),status='replace',form='unformatted')

     sheader = 'RHO'
     call write_binary_header(20, sheader, real_or_complex, nspin, &
        hamInfo%gVectors%numGVectors, nsym, cell_symmetry, nat, &
        kPoints%numKPoints, kPoints%numBands, eigenStates%maxBasisVectors, &
        ecutrho, ecutwfn, FFTgrid, kPoints%kPointNumbers, kPoints%kPointOffsets, &
        celvol, al, a, adot, recvol, bl, b, bdot, rotation, translation, &
        atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations)

     call write_binary_gvectors(20, hamInfo%gVectors%numGVectors, &
        hamInfo%gVectors%numGVectors, hamInfo%gVectors%gVectors(:,:))

     if(real_or_complex .eq. 1) then
        call write_binary_real_data(20, hamInfo%gVectors%numGVectors, &
           haminfo%gVectors%numGVectors, nspin, dataReal(:, :))
     else
        call write_binary_complex_data(20, hamInfo%gVectors%numGVectors, &
           haminfo%gVectors%numGVectors, nspin, dataComplex(:, :))
     endif

     call close_file(20)

     write(*,'(" ...done!")')

     call DensityArrayDestroy( densityArray )

     if(real_or_complex .eq. 1) then
        SAFE_DEALLOCATE(dataReal)
     else
        SAFE_DEALLOCATE(dataComplex)
     endif

  endif

  if (vxcg_flag) then

     write(*,*)

     if(real_or_complex .eq. 1) then
        SAFE_ALLOCATE(dataReal, (hamInfo%gVectors%numGVectors, nspin))
        do ispin=1,nspin
           do ig=1,hamInfo%gVectors%numGVectors
              dataReal(ig,ispin)=0.0d0
           enddo
        enddo
     else
        SAFE_ALLOCATE(dataComplex, (hamInfo%gVectors%numGVectors, nspin))
        do ispin=1,nspin
           do ig=1,hamInfo%gVectors%numGVectors
              dataComplex(ig,ispin)=(0.0d0,0.0d0)
           enddo
        enddo
     endif

     write(*,'(" Writing BerkeleyGW exchange-correlation potential file...")')

     call open_file(20,file=trim(vxcg_file),status='replace',form='unformatted')

     sheader = 'VXC'
     call write_binary_header(20, sheader, real_or_complex, nspin, &
        hamInfo%gVectors%numGVectors, nsym, cell_symmetry, nat, &
        kPoints%numKPoints, kPoints%numBands, eigenStates%maxBasisVectors, &
        ecutrho, ecutwfn, FFTgrid, kPoints%kPointNumbers, kPoints%kPointOffsets, &
        celvol, al, a, adot, recvol, bl, b, bdot, rotation, translation, &
        atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations)

     call write_binary_gvectors(20, hamInfo%gVectors%numGVectors, &
        hamInfo%gVectors%numGVectors, hamInfo%gVectors%gVectors(:,:))

     if(real_or_complex .eq. 1) then
        call write_binary_real_data(20, hamInfo%gVectors%numGVectors, &
           hamInfo%gVectors%numGVectors, nspin, dataReal(:, :))
     else
        call write_binary_complex_data(20, hamInfo%gVectors%numGVectors, &
           hamInfo%gVectors%numGVectors, nspin, dataComplex(:, :))
     endif

     call close_file(20)

     write(*,'(" ...done!")')

     if(real_or_complex .eq. 1) then
        SAFE_DEALLOCATE(dataReal)
     else
        SAFE_DEALLOCATE(dataComplex)
     endif

  endif

  SAFE_DEALLOCATE_P(atyp)
  SAFE_DEALLOCATE_P(apos)
  SAFE_DEALLOCATE(ngk)
  SAFE_DEALLOCATE(kw)
  SAFE_DEALLOCATE(kpt)
  SAFE_DEALLOCATE(ifmin)
  SAFE_DEALLOCATE(ifmax)
  SAFE_DEALLOCATE(energies)
  SAFE_DEALLOCATE(occupations)

  if (vxc_flag) then

     write(*,*)

     if (vxc_diag_nmin .lt. 1) vxc_diag_nmin = 1
     if (vxc_diag_nmax .gt. kPoints%numBands) vxc_diag_nmax = kPoints%numBands
     ndiag = MAX (vxc_diag_nmax - vxc_diag_nmin + 1, 0)
     if (vxc_offdiag_nmin .lt. 1) vxc_offdiag_nmin = 1
     if (vxc_offdiag_nmax .gt. kPoints%numBands) vxc_offdiag_nmax = kPoints%numBands
     noffdiag = MAX (vxc_offdiag_nmax - vxc_offdiag_nmin + 1, 0)
     noffdiag = noffdiag**2

     do ispin = 1, nspin
       spin_index(ispin) = ispin
     enddo
     SAFE_ALLOCATE(diag, (ndiag))
     do idiag = 1, ndiag
       diag(idiag) = vxc_diag_nmin + idiag - 1
     enddo
     SAFE_ALLOCATE(offdiag, (noffdiag))
     do ioffdiag = 1, noffdiag
       offdiag(ioffdiag) = vxc_offdiag_nmin + ioffdiag - 1
     enddo
     SAFE_ALLOCATE(vxc_mtxel, (ndiag + noffdiag, nspin))
     vxc_mtxel(:,:) = CMPLX(0d0, 0d0)

     write(*,'(" Writing BerkeleyGW exchange-correlation matrix element file...")')

     call open_file(20,file=trim(vxc_file),status='unknown',form='formatted')
     do ik=1,kPoints%numKPoints
       call write_matrix_elements(20, kPoints%kPointsLat(:,ik), nspin, ndiag, noffdiag, &
         spin_index, diag, offdiag, offdiag, vxc_mtxel(:,:))
     enddo
     call close_file(20)

     SAFE_DEALLOCATE(diag)
     SAFE_DEALLOCATE(offdiag)
     SAFE_DEALLOCATE(vxc_mtxel)

     write(*,'(" ...done!")')

  endif

  write(*,*)

  ! Clean up memory
  call EigenStatesDestroy( eigenStates )

  call HDestroy( hamiltonian )

  call GridDestroy( Grid)
                
  call HamInfoDestroy( hamInfo )

  call KPointsDestroy( kPoints )

  call StructureDestroy( structure1 )

  ! Clean up and close the input file
  call TagHandlerDestroy( tagHandler )

  call open_file(unit = 10, file = filename, status = 'old')
  call close_file(unit = 10, delete = .true.)

  999 format(a)

contains

!*************************************************************************
subroutine real_wfng (nk, nb, ns, ngkmax, ng, en, wfc, wfr)
!*************************************************************************
! Constructs real wavefunctions in G-space for systems
! with inversion symmetry by applying Gram-Schmidt process.
! Based on paratecSGL/src/para/gwreal.f90
!*************************************************************************

  use SysParams, only : double
  use message_m

  implicit none

  integer, intent(in) :: nk
  integer, intent(in) :: nb
  integer, intent(in) :: ns
  integer, intent(in) :: ngkmax
  integer, intent(in) :: ng(nk)
  real(double), intent(in) :: en(:, :) !< (nb, nk)
  complex(double), intent(in) :: wfc(:, :, :) !< (ngkmax, nb, nk)
  real(double), intent(out) :: wfr(:, :, :, :) !< (ngkmax, ns, nb, nk)

  real(double), parameter :: eps2 = 1.0d-2
  real(double), parameter :: eps5 = 1.0d-5
  real(double), parameter :: eps6 = 1.0d-6

  integer :: i, j, k, ik, ib, jb, is, ig, deg, mdeg, inc
  integer :: dimension_span, reduced_span
  real(double) :: x
  integer, allocatable :: inum (:)
  integer, allocatable :: imap (:, :)
  integer, allocatable :: inull(:)
  integer, allocatable :: null_map(:, :)
  real(double), allocatable :: psi(:, :)
  real(double), allocatable :: phi(:, :)
  real(double), allocatable :: vec(:)

  ! determine size of degenerate subspace
  mdeg = 1
  do ik = 1, nk
    do ib = 1, nb
      deg = 1
      do jb = ib + 1, nb
        if (abs(en(ib, ik) - en(jb, ik)) .lt. &
          eps5 * dble(jb - ib + 1)) deg = deg + 1
      enddo
      if (deg .gt. mdeg) mdeg = deg
    enddo
  enddo
  mdeg = mdeg * 2

  SAFE_ALLOCATE(imap, (nb, nk))
  SAFE_ALLOCATE(inum, (nk))
  SAFE_ALLOCATE(inull, (nb))
  SAFE_ALLOCATE(null_map, (mdeg, nb))

  do ik = 1, nk
    inum(ik) = 1
    do ib = 1, nb
      if (ib .eq. nb) then
        imap(inum(ik), ik) = ib
        inum(ik) = inum(ik) + 1
      elseif (abs(en(ib, ik) - en(ib + 1, ik)) .gt. eps5) then
        imap(inum(ik), ik) = ib
        inum(ik) = inum(ik) + 1
      endif
    enddo
    inum(ik) = inum(ik) - 1
  enddo

  SAFE_ALLOCATE(psi, (ngkmax, mdeg))
  SAFE_ALLOCATE(phi, (ngkmax, mdeg))
  SAFE_ALLOCATE(vec, (ngkmax))

  do ik = 1, nk
    inc = 1
    do i = 1, inum(ik)
      inull(i) = 1
      do ib = inc, imap(i, ik)
        x = 0.0d0
        do ig = 1, ng(ik)
          x = x + dble(wfc(ig, ib, ik))**2
        enddo
        if (x .lt. eps2) null_map(inull(i), i) = 0
        if (x .gt. eps2) null_map(inull(i), i) = 1
        inull(i) = inull(i) + 1
        x = 0.0d0
        do ig = 1, ng(ik)
          x = x + IMAG(wfc(ig, ib, ik))**2
        enddo
        if (x .lt. eps2) null_map(inull(i), i) = 0
        if (x .gt. eps2) null_map(inull(i), i) = 1
        inull(i) = inull(i) + 1
      enddo
      inull(i) = inull(i) - 1
      inc = imap(i, ik) + 1
    enddo
    inc = 1
    ib = 1
    do i = 1, inum(ik)
      k = 1
      do j = 1, 2 * (imap(i, ik) - inc) + 1, 2
        if (null_map(j, i) .eq. 1 .or. null_map(j + 1, i) .eq. 1) then
          if (null_map(j, i) .eq. 1) then
            do ig = 1, ng(ik)
              phi(ig, k) = dble(wfc(ig, ib, ik))
            enddo
            k = k + 1
          endif
          if (null_map(j + 1, i) .eq. 1) then
            do ig = 1, ng(ik)
              phi(ig, k) = IMAG(wfc(ig, ib, ik))
            enddo
            k = k + 1
          endif
          ib = ib + 1
        endif
      enddo
      dimension_span = k - 1
      if (dimension_span .eq. 0) then
        write(0,201)ik,inc
        stop
      endif
      do j = 1, dimension_span
        x = 0.0d0
        do ig = 1, ng(ik)
          x = x + phi(ig, j)**2
        enddo
        x = sqrt(x)
        do ig = 1, ng(ik)
          phi(ig, j) = phi(ig, j) / x
        enddo
      enddo
!
! the Gram-Schmidt process begins
!
      reduced_span = 1
      do ig = 1, ng(ik)
        psi(ig, 1) = phi(ig, 1)
      enddo
      do j = 1, dimension_span - 1
        do ig = 1, ng(ik)
          vec(ig) = phi(ig, j + 1)
        enddo
        do k = 1, reduced_span
          x = 0.0d0
          do ig = 1, ng(ik)
            x = x + phi(ig, j + 1) * psi(ig, k)
          enddo
          do ig = 1, ng(ik)
            vec(ig) = vec(ig) - psi(ig, k) * x
          enddo
        enddo
        x = 0.0d0
        do ig = 1, ng(ik)
          x = x + vec(ig)**2
        enddo
        x = sqrt(x)
        if (x .gt. eps6) then
          reduced_span = reduced_span + 1
          do ig = 1, ng(ik)
            psi(ig, reduced_span) = vec(ig) / x
          enddo
        endif
      enddo
!
! the Gram-Schmidt process ends
!
      if (reduced_span .lt. imap(i, ik) - inc + 1) then
        write(0,202) ik,inc
        stop
      endif
      do ib = inc, imap(i, ik)
        do is = 1, ns
          do ig = 1, ng(ik)
            wfr(ig, is, ib, ik) = psi(ig, ib - inc + 1)
          enddo
        enddo
      enddo
      inc = imap(i, ik) + 1
    enddo
  enddo

  SAFE_DEALLOCATE(inum)
  SAFE_DEALLOCATE(imap)
  SAFE_DEALLOCATE(inull)
  SAFE_DEALLOCATE(null_map)
  SAFE_DEALLOCATE(psi)
  SAFE_DEALLOCATE(phi)
  SAFE_DEALLOCATE(vec)

  return

  201 format(1x,"ERROR: failed Gram-Schmidt dimension span",/,14x, &
        & "k-point =",1x,i4.4,1x,"band =",1x,i4.4,/)
  202 format(1x,"ERROR: failed Gram-Schmidt reduced span",/,14x, &
        & "k-point =",1x,i4.4,1x,"band =",1x,i4.4,/)

end subroutine real_wfng

end Program EPM2BGW

