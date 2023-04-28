!*************************************************************************
Module SysParams
!*************************************************************************
  implicit none

  save

  integer, parameter         :: double = kind(0.0D0) ! double precision
  real(double), parameter    :: zero   = 0.0D0
  real(double), parameter    :: half   = 0.5D0
  real(double), parameter    :: one    = 1.0D0      
  real(double), parameter    :: two    = 2.0D0
  real(double), parameter    :: four   = 4.0D0

  ! define some mathematical & physical constants

  real(double), parameter    :: epsilon = 1.0d-30
  real(double), parameter    :: twopi = 6.2831853071795864769252867665590D0
  real(double), parameter    :: log2 = 0.301030103010301030103010D0
  real(double), parameter    :: ryd2eV = 13.60569253D0    ! Rydberg atomic unit of energy
  real(double), parameter    :: hartree2eV = two * ryd2eV ! Hartree atomic unit of energy
  real(double), parameter    :: hBar2 = 1.0d0             ! in Hartree atomic units
  real(double), parameter    :: electronMass = 1.0d0      ! in Hartree atomic units
  real(double), parameter    :: hbar = 1.054571628d-34
  real(double), parameter    :: emass = 9.10938215D-31
  real(double), parameter    :: echarge = 1.602176487D-19
  real(double), parameter    :: bohr = 0.52917721092D-10
  real(double), parameter    :: bohr2ang = bohr * 1.0D10

  complex(double), parameter :: imaginary = (0.0d0, 1.0d0)
  integer, parameter         :: maxDimensions = 3

end Module SysParams

