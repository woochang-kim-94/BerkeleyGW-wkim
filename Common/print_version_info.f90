!================================================================================
!
! Program:
!
! print_version_info      Originally By DAS        Last Modified by FHJ (Jan/19)
!
! Returns information on git repository date and commit hash.
!
!================================================================================

#include "f_defs.h"

program print_version_info

  use global_m
  use version_m
  implicit none

  character(len=256) :: string

  call get_version_info(string)

  write(6,'(a)') trim(string)

end program print_version_info
