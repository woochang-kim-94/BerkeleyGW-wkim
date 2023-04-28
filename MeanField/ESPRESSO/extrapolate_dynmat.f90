! David Strubbe, Nov 2012
! This program should be considered highly experimental.
! Read dynamical matrix (k) and coordinates (x) from matdyn of QE ph.x (4.3.2),
! and forces output from BSE/forces.(real/cplx).x.
! Write out extrapolated coordinates and Stokes shift, based on
! assumption that matdyn is the same for the excited state.
! F = kx, E = kx^2 / 2
! x = k^-1 F, E = k^-1 F^2 / 2

#include "f_defs.h"

module extrapolate_dynmat_m

  use global_m
  implicit none

  public :: write_xsf

contains

  subroutine write_xsf(iunit, natoms, speclist, avec, coords)
    integer,  intent(in) :: iunit
    integer,  intent(in) :: natoms
    real(DP), intent(in) :: avec(:,:)   !< (3,3)
    integer,  intent(in) :: speclist(:) !< (natoms)
    real(DP), intent(in) :: coords(:)   !< (3*natoms)

    integer :: idir, iatom

    write(iunit, '(a)') 'CRYSTAL'
    write(iunit, '(a)') 'PRIMVEC'
    
    do idir = 1, 3
      write(iunit,'(3f15.9)') avec(1:3, idir) * BOHR
    enddo

    write(iunit, '(a)') 'PRIMCOORD'
    write(iunit, '(2i6)') natoms, 1

    do iatom = 1, natoms
      write(iunit, '(i6,3f15.9)') speclist(iatom), coords( (iatom - 1) * 3 + 1 : (iatom - 1) * 3 + 3) * BOHR
    enddo

  end subroutine write_xsf

end module extrapolate_dynmat_m

program extrapolate_dynmat

  use global_m
  use extrapolate_dynmat_m !-nodep
  use lapack_m
  implicit none

  integer :: iunit_dynmat, iunit_dynmat_r, iunit_forces, natoms, dim, lwork, info, ijunk, iatom, jatom, ii, &
    istart, jstart, ispec, nspec, idir, jdir
  real(DP) :: djunk(5), avec(3,3), alat, real_row(3), imag_row(3), coords_latt(3), real_force(3), dist(3)
  character :: linestr*256
  integer, allocatable :: ipiv(:), speclist(:)
  complex(DPC), allocatable :: dynmat(:,:), dynmat_inv(:,:), forces(:), temp(:), work(:), identity(:,:), identity2(:,:)
  real(DP), allocatable :: coords(:), coords_diff(:)
  complex(DPC) :: energy

  ! layout: dir 1, atom 1; dir 2, atom 1; dir 3, atom 1; dir 1, atom 2; dir 2, atom 2; etc.
  ! ii = (iatom - 1) * 3 + idir

  iunit_dynmat = 77
  call open_file(iunit_dynmat, 'matdyn', status='old')
  read(iunit_dynmat, '(a256)') linestr ! "Dynamical matrix file"
  !read(iunit_dynmat, *) linestr ! blank
  !write(6,*) linestr
  read(iunit_dynmat, *) nspec, natoms, ijunk, alat, djunk(1:5)
  read(iunit_dynmat, '(a256)') linestr ! cubic
  read(iunit_dynmat, *) avec(1:3, 1)
  read(iunit_dynmat, *) avec(1:3, 2)
  read(iunit_dynmat, *) avec(1:3, 3)
!  alat = 11.2054
  avec = avec * alat
  do ispec = 1, nspec
    read(iunit_dynmat, '(a256)') linestr ! atomic masses
  enddo

  dim = 3 * natoms
  SAFE_ALLOCATE(speclist, (natoms))
  SAFE_ALLOCATE(coords, (dim))
  do iatom = 1, natoms
    istart = (iatom - 1) * 3
    read(iunit_dynmat, *) ijunk, speclist(iatom), coords_latt(1:3)
    coords(istart + 1: istart + 3) = alat * coords_latt
  enddo

  call write_xsf(6, natoms, speclist, avec, coords)

  read(iunit_dynmat, *) linestr ! blank
  read(iunit_dynmat, '(a256)') linestr ! "Dynamical Matrix in cartesian axes"
  !read(iunit_dynmat, *) linestr ! blank
  read(iunit_dynmat, '(a256)') linestr ! q-point
  write(6,*) adjustl(trim(linestr))
  !read(iunit_dynmat, *) linestr ! blank

  SAFE_ALLOCATE(dynmat, (dim, dim))
  do ii = 1, natoms**2
    read(iunit_dynmat, *) iatom, jatom

    do jdir = 1, 3
      read(iunit_dynmat, *) real_row(1), imag_row(1), real_row(2), imag_row(2), real_row(3), imag_row(3)
      do idir = 1, 3
        dynmat( (iatom - 1) * 3 + idir, (jatom - 1) * 3 + jdir) = cmplx(real_row(idir), imag_row(idir), kind = DP)
      enddo
    enddo
  enddo
  
  read(iunit_dynmat, *) linestr ! blank
  read(iunit_dynmat, '(a256)') linestr ! "Diagonalizing the dynamical matrix"
  write(6,'(a)') "Parsing complete."
  call close_file(iunit_dynmat)

  if(any(aimag(dynmat(:,:)) > TOL_Zero)) then
    write(0,'(a)') 'WARNING: dynmat has imaginary parts.'
  endif

  write(6,'(a)') "Writing dynamical matrix as a function of interatomic distance to 'matdyn_r'."
  iunit_dynmat_r = 78
  call open_file(iunit_dynmat_r, 'matdyn_r', status='replace')
  write(iunit_dynmat_r, '(a)') '# atom 1, dir 1, atom2, dir 2, distance (Ang), |dynamical matrix element| (eV/Ang^2)'
  ! file 'matdyn' is in Ry atomic units, so bohr, Ry/bohr^2.

  do iatom = 1, natoms
    istart = (iatom - 1) * 3
    do jatom = 1, natoms
      jstart = (jatom - 1) * 3

      ! record the shortest distance between periodic images of the two atoms
      ! only implemented for tetragonal cells
      do idir = 1, 3
        dist(idir) = coords(istart+idir) - coords(jstart+idir)
        if(dist(idir) >  abs(avec(idir, idir))/2) dist(idir) = dist(idir) - abs(avec(idir, idir))
        if(dist(idir) < -abs(avec(idir, idir))/2) dist(idir) = dist(idir) + abs(avec(idir, idir))
      enddo

      do idir = 1, 3
        do jdir = 1, 3
          write(iunit_dynmat_r,*) iatom, idir, jatom, jdir, sqrt(sum((dist(1:3))**2)) * BOHR, &
            abs(dynmat( (iatom - 1) * 3 + idir, (jatom - 1) * 3 + jdir)) * RYD / BOHR**2
        enddo
      enddo
    enddo
  enddo

  call close_file(iunit_dynmat_r)

  SAFE_ALLOCATE(dynmat_inv, (dim, dim))
  SAFE_ALLOCATE(temp, (dim))
  SAFE_ALLOCATE(coords_diff, (dim))

  SAFE_ALLOCATE(ipiv, (dim))
  lwork = dim ! workspace query would be better
  SAFE_ALLOCATE(work, (lwork))
  dynmat_inv = dynmat
  call zgetrf(dim, dim, dynmat_inv, dim, ipiv, info)
  if(info /= 0) then
    write(0,'(a, i6)') 'zgetrf failed with info = ', info
    call die('zgetrf failure')
  endif
  call zgetri(dim, dynmat_inv, dim, ipiv, work, lwork, info)
  if(info /= 0) then
    write(0,'(a, i6)') 'zgetri failed with info = ', info
    call die('zgetri failure')
  endif
  SAFE_DEALLOCATE(ipiv)
  SAFE_DEALLOCATE(work)

  SAFE_ALLOCATE(identity,  (dim, dim))
  SAFE_ALLOCATE(identity2, (dim, dim))
  identity  = matmul(dynmat, dynmat_inv)
  identity2 = matmul(dynmat_inv, dynmat)

  do ii = 1, dim
    identity (ii, ii) = identity (ii, ii) - 1.0d0
    identity2(ii, ii) = identity2(ii, ii) - 1.0d0
  enddo

  write(6,*) 'Max discrepancies from identity: ', maxval(abs(identity)), maxval(abs(identity2))
  
  if(any(abs(identity(:,:)) > TOL_Small .or. abs(identity2(:,:)) > TOL_Small)) then
    call die("inversion check failed")
  endif
  SAFE_DEALLOCATE(identity)
  SAFE_DEALLOCATE(identity2)

  ! Read in forces, in eV/Ang
  iunit_forces = 78
  call open_file(iunit_forces, 'forces.dat', status = 'old')
  SAFE_ALLOCATE(forces, (dim))
  do iatom = 1, natoms
    istart = (iatom - 1) * 3
    read(iunit_forces, *) ijunk, real_force(1:3)
    forces(istart + 1 : istart + 3) = cmplx(real_force(1:3), 0d0, kind = DP)
  enddo
  call close_file(iunit_forces)

  ! x = k^-1 F
  coords_diff = dble(matmul(dynmat_inv, forces)) * BOHR / RYD

  ! E = x * k * x / 2
  temp = matmul(dynmat, coords_diff)
  energy = dble(dot_product(coords_diff, temp)) * 0.5d0

  write(6,*) 'Stokes shift = ', energy / RYD

  call write_xsf(6, natoms, speclist, avec, coords_diff)

  SAFE_DEALLOCATE(dynmat)
  SAFE_DEALLOCATE(coords)
  SAFE_DEALLOCATE(coords_diff)

end program extrapolate_dynmat
