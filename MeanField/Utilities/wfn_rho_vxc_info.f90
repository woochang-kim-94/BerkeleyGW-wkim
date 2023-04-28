!=============================================================================
!
! Utilities:
!
! (1) wfn_rho_vxc_info    Originally By DAS      Last Modified 12/12/2011 (DAS)
!
!     Prints the contents of the header of a WFN, RHO, or VXC file in
!     a human-readable format.
!
!==============================================================================

#include "f_defs.h"

program wfn_rho_vxc_info

  use global_m
  use wfn_rho_vxc_io_m
  use wfn_io_hdf5_m
#ifdef HDF5
  use hdf5
#endif
  implicit none

  type(mf_header_t) :: mf
  character*256 :: infile, usage, unit_string, length_label, energy_label
  integer :: nargs, ii, iat, ik, isym, is, ierr
  logical :: is_atomic, use_hdf5
  real(DP) :: length_factor, energy_factor

  usage = 'Usage: wfn_rho_vxc_info.x wfn[.h5] [optional units: eVA/au]'

! Get file names from command-line arguments

  nargs = command_argument_count()

  if (nargs < 1 .or. nargs > 2) then
    call die(usage)
  endif

  call get_command_argument(1, infile)

  if(nargs == 2) then
    call get_command_argument(2, unit_string)
    if(trim(unit_string) == 'eVA') then
      is_atomic = .false.
    else if(trim(unit_string) == 'au') then
      is_atomic = .true.
    else
      call die("Unknown unit specifier " // trim(unit_string) // &
        ". Choose 'eVA' for eV/Angstrom, or 'au' for atomic units.")
    endif
  else
    is_atomic = .true.
  endif

  if(is_atomic) then
    length_label = "bohr"
    length_factor = 1
    energy_label = "Ry"
    energy_factor = 1
  else
    length_label = "Ang"
    length_factor = BOHR
    energy_label = "eV"
    energy_factor = RYD
  endif

  use_hdf5 = .false.
#ifdef HDF5
  use_hdf5 = index(infile, '.h5') == len(TRUNC(infile)) - 2
#endif
  if (.not.use_hdf5) then
    call open_file(unit=11, file=TRUNC(infile), form='unformatted', status='old')
    call read_mf_header(11, mf, warn=.false., dont_warn_kgrid=.true.)
    call close_file(11)
  else
#ifdef HDF5
    call h5open_f(ierr)
    call read_hdf5_mf_header(TRUNC(infile), mf)
    call h5close_f(ierr)
#endif
  endif

  write(6,'(a)') '====== GENERAL =====' 
  write(6,'(a,a)')  'Type: ', mf%sheader
  if(mf%iflavor == 1) then
    write(6,'(a)') 'Flavor: real'
  else
    write(6,'(a)') 'Flavor: complex'
  endif
  write(6,'(a,a)')  'Date created: ', mf%sdate
  write(6,'(a,a)')  'Time created: ', mf%stime
  write(6,'(a,i1)') 'Number of spins: ', mf%kp%nspin

  write(6,'(a)') '====== G-VECTORS =====' 
  write(6,'(a,i8)') 'Number of G-vectors: ', mf%gvec%ng
  write(6,'(a,f12.6,a)') 'Charge density cutoff: ', mf%gvec%ecutrho, ' Ry'
  write(6,'(a,3i8)') 'FFT grid: ', mf%gvec%FFTgrid(1:3)
  if(mf%sheader == 'WFN') then
    write(6,'(a,i8)') 'Max number of wfn G-vectors: ', mf%kp%ngkmax
    write(6,'(a,f12.6,a)') 'Wavefunction cutoff: ', mf%kp%ecutwfc, ' Ry'
  endif

  write(6,'(a)') '====== ATOMS =====' 
  write(6,'(a,i6)') 'Number of atoms: ', mf%crys%nat
  write(6,'(2a10,a30)') 'Index', 'Species', 'Coordinates (' // trim(length_label) // ')'
  do iat = 1, mf%crys%nat
    write(6,'(2i10,3f12.6)') iat, mf%crys%atyp(iat), mf%crys%apos(1:3, iat) * length_factor
  enddo

  write(6,'(a)') '====== LATTICE =====' 
  write(6,'(a,f12.6)') 'Cell volume (real space, ' // trim(length_label) // '^3): ', &
    mf%crys%celvol * length_factor**3
  write(6,'(a,f12.6)') 'Lattice constant (real space, ' // trim(length_label) // '): ', &
    mf%crys%alat  * length_factor
  write(6,'(a)') 'Lattice vectors (real space):'
  do ii = 1, 3
    write(6,'(10x,3f12.6)') mf%crys%avec(1:3, ii)
  enddo
  write(6,'(a)') 'Metric (real space, ' // trim(length_label) // '^2):'
  do ii = 1, 3
    write(6,'(10x,3f12.6)') mf%crys%adot(1:3, ii) * length_factor**2
  enddo

  write(6,'(a,f12.6)') 'Cell volume (reciprocal space, ' // trim(length_label) // '^-3): ', &
    mf%crys%recvol / length_factor**3
  write(6,'(a,f12.6)') 'Lattice constant (reciprocal space, ' // trim(length_label) // '^-1): ', &
    mf%crys%blat / length_factor
  write(6,'(a)') 'Lattice vectors (reciprocal space):'
  do ii = 1, 3
    write(6,'(10x,3f12.6)') mf%crys%bvec(1:3, ii)
  enddo
  write(6,'(a)') 'Metric (reciprocal space, ' // trim(length_label) // '^-2):'
  do ii = 1, 3
    write(6,'(10x,3f12.6)') mf%crys%bdot(1:3, ii) / length_factor**2
  enddo

  write(6,'(a)') '====== SYMMETRIES =====' 
  write(6,'(a,i2)') 'Number of symmetries: ', mf%syms%ntran
  if(mf%syms%cell_symmetry == 0) then
    write(6,'(a)') 'Symmetry type: cubic'
  else
    write(6,'(a)') 'Symmetry type: hexagonal'
  endif
  write(6,'(a7,a31,12x,a33)') 'Index', 'Rotation matrix', 'Fractional translations'
  do isym = 1, mf%syms%ntran
    write(6,'(i5,1x,a,2x,3(3i4,2x),3f12.6)') isym, ':', &
      mf%syms%mtrx(1:3, 1:3, isym), mf%syms%tnp(1:3, isym)
  enddo

  if(mf%sheader == 'WFN') then
    write(6,'(a)') '====== K-POINTS =====' 
    write(6,'(a,i8)') 'Number of k-points: ', mf%kp%nrk
    write(6,'(a,i8)') 'Number of bands: ', mf%kp%mnband
    write(6,'(a,3i4)') 'k-grid: ', mf%kp%kgrid(1:3)
    write(6,'(a,3f12.6)') 'k-shifts: ', mf%kp%shift(1:3)
    write(6,'(a)') '[ifmin = lowest occupied band, ifmax = highest occupied band, for each spin]'
    write(6,'(a8,a30,6x,a12,a22)',advance='no') 'Index', 'Coordinates (crystal)', 'Weight', 'Number of G-vectors'
    if(mf%kp%nspin == 1) then
      write(6,'(2a10)') 'ifmin', 'ifmax'
    else
      write(6,'(4a10)') 'ifmin1', 'ifmax1', 'ifmin2', 'ifmax2'
    endif
    do ik = 1, mf%kp%nrk
      write(6,'(i8,4f12.6,i22,4i10)') ik, mf%kp%rk(1:3, ik), mf%kp%w(ik), mf%kp%ngk(ik), &
        (mf%kp%ifmin(ik, is), mf%kp%ifmax(ik, is), is = 1, mf%kp%nspin)
    enddo

    write(6,'(a)') '====== ENERGIES (' // trim(energy_label) // ')/OCCUPATIONS ====='
    do is = 1, mf%kp%nspin
      if(mf%kp%nspin > 1) write(6,'(a,i2)') 'Spin ', is
      do ik = 1, mf%kp%nrk
        write(6,'(a,i6)') 'k-point ', ik
        write(6,'(9999999f14.6)') mf%kp%el(:, ik, is) * energy_factor
        write(6,'(9999999f14.6)') mf%kp%occ(:, ik, is)
      enddo
    enddo
  endif

  call dealloc_header_type(mf%sheader, mf%crys, mf%kp)

end program wfn_rho_vxc_info
