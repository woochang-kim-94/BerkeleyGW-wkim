!>===================================================================
!!
!! Module os_m      Originally by FHJ    Last Modified Jan/2018 (FHJ)
!!
!! Also, see some of the corresponding C routines in os_c.c
!!
!! Routines:
!!
!! - system
!! - sleep
!! - c_to_f_string
!! - get_host_name
!!
!!===================================================================

#include "f_defs.h"

module os_m

  use, intrinsic  :: iso_c_binding
  use nrtype_m
  implicit none

  private

  public ::        &
    system,        &
    sleep,         &
    c_to_f_string, &
    get_host_name

contains


  !> Call command `command` using call to POSIX `system` routine.
  subroutine system(command)
    character(len=*), intent(in) :: command

    character(len=1,kind=c_char) :: command_c(len(command)+1)
    integer(c_int) :: ret
    integer :: ii
    interface !int system(const char *command);
      integer(c_int) function system_c(command) bind(c, name='system')
        use, intrinsic :: iso_c_binding, only: c_int, c_char
        character(len=1,kind=c_char), intent(in) :: command(*)
      end function system_c
    end interface

    do ii = 1,  len(command)
      command_c(ii) = command(ii:ii)
    end do
    command_c(len(command)+1:len(command)+1) = c_null_char 

    ret = system_c(command_c)

  end subroutine system


  !> Sleep for `time` seconds. Internally calls POSIX `nanosecond`.
  subroutine sleep(time)
    real(DP), intent(in) :: time

    integer :: ret

    interface !int sleep_c(double s); // defined in os_c.c
      integer(c_int) function sleep_c(s) bind(c)
        use, intrinsic :: iso_c_binding, only: c_int, c_double
        real(c_double), value, intent(in) :: s ! pass by value, not by reference
      end function sleep_c
    end interface

    ret = sleep_c(time)

  end subroutine sleep


  !> Convert a C string to a Fotran string.
  function c_to_f_string(c_string) result(f_string)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    character(kind=c_char,len=1), intent(in) :: c_string(*)
    character(len=:), allocatable :: f_string
    integer(c_size_t) :: n

    interface !size_t strlen(const char *s);
      integer(c_size_t) function strlen(s) bind(c)
        use, intrinsic :: iso_c_binding, only: c_char, c_size_t
        character(len=1,kind=c_char), intent(in) :: s(*)
      end function strlen
    end interface

    n = strlen(c_string)
    allocate(character(len=n) :: f_string)
    f_string = transfer(c_string(1:n), f_string)
  end function c_to_f_string


  !> Get the host name for the current processor / MPI rank.
  function get_host_name() result(hostname)
    character(len=:), allocatable :: hostname
    integer, parameter :: max_hostname_len=1024
    character(len=max_hostname_len), target :: hostname_
    integer :: namelen, mpierr
#ifndef MPI
    integer(c_int) :: c_status
    character(kind=c_char,len=1) :: c_hostname(max_hostname_len)
    integer(c_size_t) :: c_namelen

    interface !int gethostname(char *name, size_t namelen);
      integer(c_int) function gethostname(name, namelen) bind(c)
        use, intrinsic :: iso_c_binding, only: c_char, c_int, c_size_t
        integer(c_size_t), value, intent(in) :: namelen
        character(len=1,kind=c_char), dimension(namelen), intent(inout) :: name
      end function gethostname
    end interface
#endif

#ifdef MPI
    ! FHJ: If we are using MPI, we call directly MPI_Get_processor_name
    namelen = max_hostname_len
    call MPI_Get_processor_name(hostname_, namelen, mpierr)
    namelen = len_trim(hostname_)
    allocate(character(len=namelen) :: hostname)
    hostname(1:namelen) = hostname_(1:namelen)
#else
    ! FHJ: Otherwise use system call to gethostname
    c_namelen = int(max_hostname_len, kind=c_size_t)
    c_status = gethostname(c_hostname, c_namelen)
    hostname = c_to_f_string(c_hostname)
#endif

  end function get_host_name

end module os_m
