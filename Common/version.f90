!================================================================================
!
! Modules:
!
! (1) version_m           Originally By FHJ
!
! Returns information on GIT repository date and commit hash.
!
!================================================================================

#include "f_defs.h"
#include "version_base.h"
#include "version_git.h"

module version_m

  implicit none

  private

  character(len=*), parameter :: version_base = VERSION_BASE
  logical, parameter :: is_git = IS_GIT
  character(len=*), parameter :: git_date = GIT_DATE
  character(len=*), parameter :: git_commit = GIT_COMMIT
  character(len=*), parameter :: git_author = GIT_AUTHOR

  public :: get_version_info

contains

  subroutine get_version_info(ainfo)
    character(len=256), intent(out) :: ainfo

    if (is_git) then
      ainfo = 'version ' // trim(version_base) // '+git ' // &
              trim(git_commit) // ' by ' // trim(git_author) // ' at ' // trim(git_date)
    else
      ainfo = 'version ' // trim(version_base)
    endif

  end subroutine get_version_info

end module version_m
