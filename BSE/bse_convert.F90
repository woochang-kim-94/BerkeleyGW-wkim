!============================================================================
!
! Module bse_convert_m
!
!============================================================================

#include "f_defs.h"

module bse_convert_m

  use global_m
  use kernel_io_m

  implicit none

  private

  public ::              &
    bsemat_binasc,       &
    dtmat_binasc,        &
    vmtxel_binasc,       &
    eps2_moments_binasc, &
    eigenvectors_binasc, &
    bsemat_ascbin,       &
    dtmat_ascbin,        &
    vmtxel_ascbin,       &
    eps2_moments_ascbin, &
    eigenvectors_ascbin

contains

#define NAME(x) x ## _binasc
! this is a variadic macro
#define READWRITE(...) \
read(iunit) __VA_ARGS__; \
write(ounit, *) __VA_ARGS__

#include "bse_convert_inc.f90"

#undef NAME
#undef READWRITE

#define NAME(x) x ## _ascbin
#define READWRITE(...) \
read(iunit, *) __VA_ARGS__; \
write(ounit) __VA_ARGS__

#include "bse_convert_inc.f90"

end module bse_convert_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
