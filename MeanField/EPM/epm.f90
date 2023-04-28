!*************************************************************************
Program EPM
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

  use sysParams,           only : double, zero, one, hartree2eV, hbar, bohr, &
                                  echarge, emass
  use TagHandlerMod,       only : tagHandlerT, TagHandlerInit, TagHandlerDestroy, &
                                  & FindTag
  use StructureMod,        only : StructureT, StructureInit, StructureDestroy
  use kPointsMod,          only : kPointsT, KPointsInit, KPointsDestroy
  use graphMod,            only : graphT, PlotBands
  use hamiltonianMod,      only : hamiltonianArrayT, HInit, HDestroy, HDiagonalize, & 
                                  & HPrint
  use eigenStatesMod,      only : eigenStatesT, EigenStatesInit, EigenStatesDestroy, &
                                  & PrintEigenvalues
  use pwHamMod,            only : hamInfoT, HamInfoInit, HamInfoDestroy, &
                                  & HGenerate, FindHMaxSize, GridInit, &
                                  & PotentialRGenerate, KineticGenerate, GridDestroy
  use cgParamsMod,         only : cgParamsT, cgParamsInit
  use typeMod,             only : GridT
  use ConjGradMod,         only : cgEigenSystem,cgFFTDiagonalise

  use message_m, only : open_file, close_file
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

  logical :: gap_flag, mass_flag, v_of_q
  integer :: kpt, error, ReportNo, k, n, i, kv, kc
  real(double) :: homo, lumo, gap, e1, e2, e3, d12, d23
  real(double) :: kvPoint(3), kcPoint(3), dk, md(3,3)
  real(double) :: dv1, dv2, dv3, dc1, dc2, dc3, a, b, c
  real(double) :: mv1h, mv2h, mv3h, mv1l, mv2l, mv3l, mc1, mc2, mc3
  real(double) :: mvl, mvh, mct, mcl
  real(double) :: eps, abstol, ff(6)
  character(80) :: filename
  character(256) :: str
  character(256) :: gap_file
  character(256) :: mass_file

  nullify(tagHandler)
  nullify(hamInfo)
  nullify(structure1)
  nullify(kPoints)
  nullify(eigenStates)
  nullify(hamiltonian)
  nullify(graphinfo)
  nullify(Grid)

  ReportNo = 110

  filename = 'input.tmp'
  call open_file(unit = 10, file = filename, status = 'replace')
  i = 0
  do while (i .eq. 0)
     read(unit = 5, fmt = '(a)', iostat = i) str
     write(unit = 10, fmt = '(a)') str
  enddo
  call close_file(unit = 10)

  ! Necessary if using conjugate-gradients method
  call Random_Seed

  ! Initialize "tag handler" to read the input from file 'filename'.
  call TagHandlerInit( tagHandler, 10, filename)

  ! Input parameter gap_flag for BerkeleyGW band gap file
  call FindTag( tagHandler, 'gap_flag', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) gap_flag
  endif
  if(error .ne. 0) gap_flag = .false.
  write(*,'(" gap_flag ",l1)')gap_flag

  ! Input parameter gap_file for BerkeleyGW band gap file
  call FindTag( tagHandler, 'gap_file', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) gap_file
  endif
  if(error .ne. 0) write(gap_file,'("gap.dat")')
  write(*,'(" gap_file ",a)')trim(gap_file)

  ! Input parameter mass_flag for BerkeleyGW effective mass file
  call FindTag( tagHandler, 'mass_flag', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) mass_flag
  endif
  if(error .ne. 0) mass_flag = .false.
  write(*,'(" mass_flag ",l1)')mass_flag

  ! Input parameter mass_file for BerkeleyGW effective mass file
  call FindTag( tagHandler, 'mass_file', error)
  if(error .eq. 0) then
     read(unit = tagHandler%fileno, fmt = *, iostat = error) mass_file
  endif
  if(error .ne. 0) write(mass_file,'("mass.dat")')
  write(*,'(" mass_file ",a)')trim(mass_file)

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
  Call GridInit(Grid,hamInfo,eigenStates,structure1%ndim)

! DAS -- do not write this useless file if not using CG
  if(cgParams%Switch.eq.1) then
    call open_file(unit=ReportNo,file='Report.cg',status='replace', position='rewind')
    write(ReportNo,*) 'Eigenvalues:'
    write(ReportNo,*) 'The FFT Grid Size +/-',Grid%Size(1)
    write(ReportNo,*) 'The number of GVectors:',Grid%NumgVec
    write(ReportNo,*) 'Eigenvalues:'
  endif

  ! Calculate the Potential and FFT to Vr
  Call PotentialRGenerate(hamInfo,Grid)

  if(cgParams%Switch.eq.1) then 
     do kpt = 1, kPoints%numKPoints 
        call HGenerate( hamInfo, structure1, kPoints, eigenStates, &
             & kpt, hamiltonian%h, hamiltonian%sqH, hamiltonian%hSize )

        write(ReportNo,*) 'k:',kpt

        call cgEigenSystem( hamiltonian,eigenStates,kpt,cgParams%tol, &
             & cgParams%period,ReportNo)

        print *, kpt,'completed out of a total of ',kPoints%numKPoints

     end do ! kpt
  else
     do kpt = 1, kPoints%numKPoints

        call HGenerate( hamInfo, structure1, kPoints, eigenStates, &
             & kpt, hamiltonian%h, hamiltonian%sqH, hamiltonian%hSize )

! DAS -- do not write this useless file      
!        write(ReportNo,*) 'k:',kpt
     
        call HDiagonalize( hamiltonian, eigenStates, kpt )

        print *, kpt,'completed out of a total of ',kPoints%numKPoints
        
     end do ! kpt

  end if

  if(cgParams%Switch.eq.1) then 
    call close_file(unit=ReportNo)
  endif

  ! Convert eigenvalues from Hartree to eV
    do k = 1, kPoints%numKPoints
       do n = 1, kPoints%numBands
          eigenStates%eigenValues(n,k) = eigenStates%eigenValues(n,k) * hartree2eV
       enddo
    enddo

  ! Creates files bands.dat and bands.plt (Display with: gnuplot bands.plt)
  if( associated( graphinfo ) )    call PlotBands( graphinfo, eigenStates )

  if( .not. associated( graphinfo ) )  call PrintEigenvalues(eigenStates, kPoints%numKPoints)

  if (gap_flag .or. mass_flag) then

     homo=-1.0d6
     lumo=1.0d6

     do k=1,kPoints%numKPoints
        do n=1,kPoints%numBands
           if (n.eq.kPoints%numOccupiedBands &
           & .and.eigenStates%eigenValues(n,k).gt.homo) then
              kv=k
              homo=eigenStates%eigenValues(n,k)
           endif
           if (n.eq.kPoints%numOccupiedBands+1 &
           & .and.eigenStates%eigenValues(n,k).lt.lumo) then
              kc=k
              lumo=eigenStates%eigenValues(n,k)
           endif
        enddo
     enddo

  endif

  if (gap_flag) then

     write(*,*)

     if (kv.eq.1.or.kv.eq.kPoints%numKPoints) then
        homo=eigenStates%eigenValues(kPoints%numOccupiedBands,kv)
     else
        e1=eigenStates%eigenValues(kPoints%numOccupiedBands,kv-1)
        e2=eigenStates%eigenValues(kPoints%numOccupiedBands,kv)
        e3=eigenStates%eigenValues(kPoints%numOccupiedBands,kv+1)
        d12=0.0d0
        d23=0.0d0
        do i=1,3
           d12=d12+(kPoints%kPointsLat(1,kv)-kPoints%kPointsLat(1,kv-1))**2
           d23=d23+(kPoints%kPointsLat(1,kv+1)-kPoints%kPointsLat(1,kv))**2
        enddo
        d12=dsqrt(d12)
        d23=dsqrt(d23)
        c=e2
        a=((e1-c)*d23+(e3-c)*d12)/(d12*d23*(d12+d23))
        b=(e3-a*d23**2-c)/d23
        homo=c-b**2/(4.0d0*a)
     endif

     if (kc.eq.1.or.kc.eq.kPoints%numKPoints) then
        lumo=eigenStates%eigenValues(kPoints%numOccupiedBands+1,kc)
     else
        e1=eigenStates%eigenValues(kPoints%numOccupiedBands+1,kc-1)
        e2=eigenStates%eigenValues(kPoints%numOccupiedBands+1,kc)
        e3=eigenStates%eigenValues(kPoints%numOccupiedBands+1,kc+1)
        d12=0.0d0
        d23=0.0d0
        do i=1,3
           d12=d12+(kPoints%kPointsLat(1,kc)-kPoints%kPointsLat(1,kc-1))**2
           d23=d23+(kPoints%kPointsLat(1,kc+1)-kPoints%kPointsLat(1,kc))**2
        enddo
        d12=dsqrt(d12)
        d23=dsqrt(d23)
        c=e2
        a=((e1-c)*d23+(e3-c)*d12)/(d12*d23*(d12+d23))
        b=(e3-a*d23**2-c)/d23
        lumo=c-b**2/(4.0d0*a)
     endif

     gap=lumo-homo

     write(*,'(" Writing band gap...")')

     call open_file(20,file=trim(gap_file),status='replace',form='formatted')
     write(20,'(1x,"Eg =",1x,f9.6,1x,"eV")')gap
     call close_file(20)

     write(*,'(" ...done!")')

  endif

  if (mass_flag) then

     write(*,*)
     write(0,'(1x,a,/)') "WARNING: The effective-mass calculation assumes a bandstructure similar to bulk Si."

     kvPoint(:) = kPoints%kPointsLat(:,kv)
     kcPoint(:) = kPoints%kPointsLat(:,kc)

     kPoints%numKPoints = 14
     do k=1,7
       kPoints%kPointsLat(:,k)=kvPoint(:)
     enddo
     do k=8,14
       kPoints%kPointsLat(:,k)=kcPoint(:)
     enddo

     dk = 1.0d-3

     md(1,1) = 1.0d0/dsqrt(2.0d0)
     md(2,1) = 0.0d0
     md(3,1) = 1.0d0/dsqrt(2.0d0)
     md(1,2) = 0.0d0
     md(2,2) = 1.0d0/dsqrt(2.0d0)
     md(3,2) = 1.0d0/dsqrt(2.0d0)
     md(1,3) = 1.0d0/dsqrt(2.0d0)
     md(2,3) = 1.0d0/dsqrt(2.0d0)
     md(3,3) = 0.0d0

     kPoints%kPointsLat(:, 2)  = kPoints%kPointsLat(:, 2)  - dk * md(:,1)
     kPoints%kPointsLat(:, 3)  = kPoints%kPointsLat(:, 3)  + dk * md(:,1)
     kPoints%kPointsLat(:, 4)  = kPoints%kPointsLat(:, 4)  - dk * md(:,2)
     kPoints%kPointsLat(:, 5)  = kPoints%kPointsLat(:, 5)  + dk * md(:,2)
     kPoints%kPointsLat(:, 6)  = kPoints%kPointsLat(:, 6)  - dk * md(:,3)
     kPoints%kPointsLat(:, 7)  = kPoints%kPointsLat(:, 7)  + dk * md(:,3)

     kPoints%kPointsLat(:, 9)  = kPoints%kPointsLat(:, 9)  - dk * md(:,1)
     kPoints%kPointsLat(:, 10) = kPoints%kPointsLat(:, 10) + dk * md(:,1)
     kPoints%kPointsLat(:, 11) = kPoints%kPointsLat(:, 11) - dk * md(:,2)
     kPoints%kPointsLat(:, 12) = kPoints%kPointsLat(:, 12) + dk * md(:,2)
     kPoints%kPointsLat(:, 13) = kPoints%kPointsLat(:, 13) - dk * md(:,3)
     kPoints%kPointsLat(:, 14) = kPoints%kPointsLat(:, 14) + dk * md(:,3)

     do kpt = 1, kPoints%numKPoints
        call HGenerate( hamInfo, structure1, kPoints, eigenStates, &
             & kpt, hamiltonian%h, hamiltonian%sqH, hamiltonian%hSize )
        call HDiagonalize( hamiltonian, eigenStates, kpt )
     enddo

     n=kPoints%numOccupiedBands-2
     dv1=-(eigenStates%eigenValues(n, 2)+eigenStates%eigenValues(n, 3)-2.0d0*eigenStates%eigenValues(n, 1)) &
       /(dk*structure1%lattice%bLatVecLen(1))**2
     dv2=-(eigenStates%eigenValues(n, 4)+eigenStates%eigenValues(n, 5)-2.0d0*eigenStates%eigenValues(n, 1)) &
       /(dk*structure1%lattice%bLatVecLen(2))**2
     dv3=-(eigenStates%eigenValues(n, 6)+eigenStates%eigenValues(n, 7)-2.0d0*eigenStates%eigenValues(n, 1)) &
       /(dk*structure1%lattice%bLatVecLen(3))**2
     dv1=dv1*hartree2eV*echarge*bohr**2
     dv2=dv2*hartree2eV*echarge*bohr**2
     dv3=dv3*hartree2eV*echarge*bohr**2
     mv1l=hbar**2/dv1/emass
     mv2l=hbar**2/dv2/emass
     mv3l=hbar**2/dv3/emass

     n=kPoints%numOccupiedBands
     dv1=-(eigenStates%eigenValues(n, 2)+eigenStates%eigenValues(n, 3)-2.0d0*eigenStates%eigenValues(n, 1)) &
       /(dk*structure1%lattice%bLatVecLen(1))**2
     dv2=-(eigenStates%eigenValues(n, 4)+eigenStates%eigenValues(n, 5)-2.0d0*eigenStates%eigenValues(n, 1)) &
       /(dk*structure1%lattice%bLatVecLen(2))**2
     dv3=-(eigenStates%eigenValues(n, 6)+eigenStates%eigenValues(n, 7)-2.0d0*eigenStates%eigenValues(n, 1)) &
       /(dk*structure1%lattice%bLatVecLen(3))**2
     dv1=dv1*hartree2eV*echarge*bohr**2
     dv2=dv2*hartree2eV*echarge*bohr**2
     dv3=dv3*hartree2eV*echarge*bohr**2
     mv1h=hbar**2/dv1/emass
     mv2h=hbar**2/dv2/emass
     mv3h=hbar**2/dv3/emass

     n=kPoints%numOccupiedBands+1
     dc1=(eigenStates%eigenValues(n, 9)+eigenStates%eigenValues(n,10)-2.0d0*eigenStates%eigenValues(n, 8)) &
       /(dk*structure1%lattice%bLatVecLen(1))**2
     dc2=(eigenStates%eigenValues(n,11)+eigenStates%eigenValues(n,12)-2.0d0*eigenStates%eigenValues(n, 8)) &
       /(dk*structure1%lattice%bLatVecLen(2))**2
     dc3=(eigenStates%eigenValues(n,13)+eigenStates%eigenValues(n,14)-2.0d0*eigenStates%eigenValues(n, 8)) &
       /(dk*structure1%lattice%bLatVecLen(3))**2
     dc1=dc1*hartree2eV*echarge*bohr**2
     dc2=dc2*hartree2eV*echarge*bohr**2
     dc3=dc3*hartree2eV*echarge*bohr**2
     mc1=hbar**2/dc1/emass
     mc2=hbar**2/dc2/emass
     mc3=hbar**2/dc3/emass

! DAS -- This is where strong assumptions about the bandstructure were made by GSM.
!     mv=(dsqrt(mv1h*mv2h*mv3h)+dsqrt(mv1l*mv2l*mv3l))**(2.0d0/3.0d0)
!     mc=(6.0d0*dsqrt(mc1*mc2*mc3))**(2.0d0/3.0d0)

     eps=1.0d-6
     mvl=0.0d0
     mvh=0.0d0
     mct=0.0d0
     mcl=0.0d0
     if (abs(mv1l-mv2l).lt.eps.and.abs(mv2l-mv3l).lt.eps) mvl=mv1l
     if (abs(mv1h-mv2h).lt.eps.and.abs(mv2h-mv3h).lt.eps) mvh=mv1h
     if (abs(mc1-mc2).lt.eps) then
        mcl=mc3
        mct=mc1
     elseif (abs(mc2-mc3).lt.eps) then
        mcl=mc1
        mct=mc2
     elseif (abs(mc1-mc3).lt.eps) then
        mcl=mc2
        mct=mc1
     endif

     write(*,'(" Writing effective masses...")')

     call open_file(20,file=trim(mass_file),status='replace',form='formatted')
!     write(20,'(1x,"mh =",1x,f9.6,1x,"m0")')mv
!     write(20,'(1x,"me =",1x,f9.6,1x,"m0")')mc
     write(20,'("        light hole     - mlh = ",f9.6," m0")')mvl
     write(20,'("        heavy hole     - mhh = ",f9.6," m0")')mvh
     write(20,'(" longitudinal electron - mle = ",f9.6," m0")')mcl
     write(20,'("   transverse electron - mte = ",f9.6," m0")')mct
     call close_file(20)

     write(*,'(" ...done!")')

  endif

  write(*,*)

  ! Clean up memory
  call EigenStatesDestroy( eigenStates )

  call HDestroy( hamiltonian )

  call GridDestroy ( Grid )
                
  call HamInfoDestroy( hamInfo )

  call KPointsDestroy( kPoints )

  call StructureDestroy( structure1 )

  ! Clean up and close the input file
  call TagHandlerDestroy( tagHandler )

  call open_file(unit = 10, file = filename, status = 'old')
  call close_file(unit = 10, delete = .true.)

end Program EPM

